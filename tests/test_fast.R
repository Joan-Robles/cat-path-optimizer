# Equivalencia del motor rápido (src/engine.cpp) con el motor en R.
# Ejecutar desde la raíz del proyecto:
#   Rscript -e 'testthat::test_file("tests/test_fast.R")'

library(testthat)

find_root <- function(path = getwd()) {
  while (!file.exists(file.path(path, "cats_game.Rproj"))) {
    parent <- dirname(path)
    if (parent == path) stop("No se encontró cats_game.Rproj")
    path <- parent
  }
  path
}
root <- find_root()
for (f in c("functions.R", "cards.R", "game.R", "policies.R", "fast.R")) source(file.path(root, f))
load_fast_engine(root)

# Ayudas -------------------------------------------------------------------

random_board <- function() {
  b <- array(0L, dim = c(4, 4, 4, 4))
  for (i in 1:4) for (j in 1:4)
    b[i, j, , ] <- rotate_card(card_types[[sample.int(6, 1)]], sample(c(0, 90, 180, 270), 1))
  b
}

# Estado al azar con varios gatos y hogares en el anillo, para forzar capturas
random_state <- function() {
  ring_type <- rep("", 16); ring_points <- rep(NA_real_, 16)
  pos <- sample.int(16, sample(4:10, 1))
  n_home <- sample(seq_len(length(pos) - 1), 1)
  ring_type[pos[seq_len(n_home)]] <- "home"
  cats <- pos[-seq_len(n_home)]
  ring_type[cats] <- "cat"
  ring_points[cats] <- sample(1:15, length(cats), replace = TRUE)
  list(board = random_board(), deck = build_deck(game_params()$deck_counts)[1:10],
       cat_deck = as.numeric(sample(1:15, sample(0:6, 1), replace = TRUE)),
       ring_type = ring_type, ring_points = ring_points, scores = 0, cats_won = 0L,
       n_players = 1L, player = 1L, turn = 1L, n_cats = 15L, params = game_params(),
       history = list())
}

fe_resolve_state <- function(s) {
  fe_resolve(board_codes(s$board), ring_codes(s$ring_type), s$ring_points,
             s$cat_deck, as.integer(s$n_cats - sum(s$cats_won)))
}

# Estados de partidas reales: jugadas al azar desde repartos distintos
game_states <- function(n_games, n_turns) {
  out <- list()
  for (g in seq_len(n_games)) {
    s <- new_game(game_params(), seed = 100 + g)
    pol <- policy_random(p_insert = 0.6)
    for (k in seq_len(n_turns)) {
      if (game_over(s)) break
      s <- play_turn(s, pol)
      out[[length(out) + 1L]] <- s
    }
  }
  out
}

# Pruebas ------------------------------------------------------------------

test_that("las tablas de códigos reconocen cada tipo y orientación", {
  for (t in seq_along(card_types)) for (r in 0:3) {
    code <- card_code(rotate_card(card_types[[t]], 90 * r))
    expect_equal(code %/% 4L, t - 1L)
    expect_equal(FE$mats[[code + 1L]], rotate_card(card_types[[t]], 90 * r))
  }
})

test_that("rotar e insertar sobre códigos coincide con el motor en R", {
  set.seed(1)
  for (k in 1:300) {
    b <- random_board()
    codes <- board_codes(b)
    cell <- sample.int(4, 2, replace = TRUE); q <- sample(c(1L, 3L), 1)
    rot <- codes; kk <- (cell[1] - 1) * 4 + cell[2]
    rot[kk] <- FE$canon[(rot[kk] %/% 4L) * 4L + (rot[kk] %% 4L + q) %% 4L + 1L]
    expect_equal(rot, board_codes(shift_one_new_card(b, cell, 90 * q)))

    card <- rotate_card(card_types[[sample.int(6, 1)]], sample(c(0, 90, 180, 270), 1))
    pos <- sample.int(16, 1)
    s <- list(board = b, deck = list())
    s2 <- apply_insertion(s, card, pos)
    ins <- fe_insert(codes, card_code(card), pos)
    expect_equal(ins$codes, board_codes(s2$board))
    expect_equal(ins$expelled, card_code(s2$deck[[1]]))
  }
})

test_that("resolve_captures coincide con el motor en R en estados al azar", {
  set.seed(2)
  n_capt <- 0
  for (k in 1:2000) {
    s <- random_state()
    ref <- resolve_captures(s)
    fe <- fe_resolve_state(s)
    expect_equal(fe$ring, ring_codes(ref$ring_type))
    expect_equal(fe$pts, ref$ring_points)
    expect_equal(fe$cats, ref$cats_won)
    expect_equal(fe$points, ref$scores)
    expect_equal(fe$cats_used, length(s$cat_deck) - length(ref$cat_deck))
    hist <- do.call(rbind, lapply(ref$history, function(h) data.frame(cat = h$cat, home = h$home)))
    if (!is.null(hist)) expect_equal(fe$history[, c("cat", "home")], hist)
    expect_equal(fe$reach, cat_reach(ref))
    n_capt <- n_capt + fe$cats
  }
  expect_gt(n_capt, 1000)  # los estados sí ejercitan capturas y cadenas
})

test_that("a una acción, los valores coinciden con los de policy_greedy", {
  set.seed(3)
  ev_rot <- function(s, act) {
    m <- masked_state(s)
    state_value(m, resolve_captures(apply_rotation(m, act$cell, act$angle)))
  }
  ev_ins <- function(s, card) {
    m <- masked_state(s)
    grid <- expand.grid(pos = 1:16, angle = distinct_rotations(card))
    max(mapply(function(p, a) state_value(m, resolve_captures(apply_insertion(m, card, p, a))),
               grid$pos, grid$angle))
  }
  states <- game_states(n_games = 6, n_turns = 5)
  for (s in states) {
    if (game_over(s)) next
    fs <- fe_state(s, belief = "mazo")
    ch <- fe_children(fs$codes, fs$ring, fs$pts, fs$cat_deck, fs$to_win, fs$unk, fs$known,
                      left = 1L, can_insert = TRUE, mandatory = TRUE,
                      use_points = FALSE, w_reach = 0.001)
    cands <- rotation_candidates(s)
    ref_rot <- vapply(cands, function(a) ev_rot(s, a), numeric(1))
    fe_cells <- (ch$rot[, 1] - 1)
    expect_equal(cbind(fe_cells %/% 4 + 1, fe_cells %% 4 + 1, 90 * ch$rot[, 2]),
                 t(vapply(cands, function(a) c(a$cell, a$angle), numeric(3))),
                 ignore_attr = TRUE)
    expect_equal(ch$rot[, 3], ref_rot)

    types <- table(vapply(s$deck, card_type_of, character(1)))
    ref_by_type <- vapply(names(types), function(t) ev_ins(s, card_types[[t]]), numeric(1))
    expect_equal(unname(ch$ins_by_type[match(names(types), names(card_types))]),
                 unname(ref_by_type))
    expect_equal(ch$ins_ev, sum(ref_by_type * types) / sum(types))
  }
})

test_that("la creencia exacta usa solo las cartas que nadie ha visto", {
  s <- new_game(game_params(), seed = 4)
  expect_equal(s$n_unknown, 10L)
  fs <- fe_state(s)
  expect_equal(sum(fs$unk), 10L)
  expect_length(fs$known, 0)
  d <- draw_card(s)
  s2 <- apply_insertion(d$state, d$card, pos = 5)
  fs2 <- fe_state(s2)
  expect_equal(s2$n_unknown, 9L)
  expect_equal(sum(fs2$unk), 9L)
  expect_equal(fs2$known, type_index(s2$deck[[10]]))
  expect_equal(sum(fe_state(s2, belief = "mazo")$unk), 10L)
})

test_that("la búsqueda sobre el turno completo juega partidas válidas", {
  for (rm in c(TRUE, FALSE)) {
    res <- play_game(game_params(cat_points = 1:15, rotation_mandatory = rm),
                     policy_search(horizon = "turno"), seed = 1)
    expect_true(res$completed)
    expect_equal(res$scores, sum(1:15))
  }
})
