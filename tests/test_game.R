# Pruebas del motor de reglas. Ejecutar desde la raíz del proyecto:
#   Rscript -e 'testthat::test_file("tests/test_game.R")'

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
for (f in c("functions.R", "cards.R", "game.R", "policies.R")) source(file.path(root, f))

# Ayudas -------------------------------------------------------------------

ct <- card_types
board_fill <- function(card) {
  b <- array(0L, dim = c(4, 4, 4, 4))
  for (i in 1:4) for (j in 1:4) b[i, j, , ] <- card
  b
}
test_state <- function(board, homes, cats, cat_points, cat_deck = numeric(0),
                       deck = rep(list(ct$puente), 10), params = game_params()) {
  ring_type <- rep("", 16); ring_points <- rep(NA_real_, 16)
  ring_type[homes] <- "home"
  ring_type[cats] <- "cat"
  ring_points[cats] <- cat_points
  list(board = board, deck = deck, cat_deck = cat_deck, ring_type = ring_type,
       ring_points = ring_points, scores = 0, cats_won = 0L, n_players = 1L,
       player = 1L, turn = 1L, n_cats = 15L, params = params, history = list())
}
types_of_board <- function(board) {
  as.vector(vapply(1:4, function(i) vapply(1:4, function(j) card_type_of(board[i, j, , ]),
                                           character(1)), character(4)))
}

# Cartas -------------------------------------------------------------------

test_that("cuatro rotaciones de 90 grados devuelven la carta original", {
  for (card in ct) {
    r <- card
    for (k in 1:4) r <- rotate_card(r, 90)
    expect_equal(r, card)
    expect_equal(rotate_card(rotate_card(card, 180), 180), card)
  }
})

test_that("la rotación de 90 grados es en sentido horario", {
  expect_equal(open_sides(ct$curva), c(TRUE, TRUE, FALSE, FALSE))                 # N-E
  expect_equal(open_sides(rotate_card(ct$curva, 90)), c(FALSE, TRUE, TRUE, FALSE)) # E-S
  expect_equal(open_sides(rotate_card(ct$curva, 270)), c(TRUE, FALSE, FALSE, TRUE)) # O-N
})

test_that("orientaciones distintas por tipo", {
  n_rot <- vapply(ct, function(m) length(distinct_rotations(m)), integer(1))
  expect_equal(unname(n_rot), c(4L, 2L, 4L, 2L, 1L, 1L))
})

test_that("una rotación gira solo un cuarto de vuelta", {
  s <- test_state(board_fill(ct$curva), homes = 9, cats = 5, cat_points = 1)
  expect_error(apply_rotation(s, c(1, 1), 180))
  expect_equal(apply_rotation(s, c(1, 1), 270)$board[1, 1, , ], rotate_card(ct$curva, 270))
  # Candidatos por tipo: curva y te 2; recta y doble curva 1 (90 = 270);
  # puente y cruz ninguno
  n_cand <- vapply(ct, function(card) length(rotation_candidates(
    test_state(board_fill(card), homes = 9, cats = 5, cat_points = 1))), integer(1))
  expect_equal(unname(n_cand), 16L * c(2L, 1L, 2L, 1L, 0L, 0L))
})

test_that("card_type_of reconoce cualquier rotación", {
  for (nm in names(ct)) for (a in c(0, 90, 180, 270))
    expect_equal(card_type_of(rotate_card(ct[[nm]], a)), nm)
})

# Anillo -------------------------------------------------------------------

test_that("las posiciones del anillo apuntan a la celda y el lado correctos", {
  expect_equal(ring_cell_side(1), c(4L, 1L, 4L))
  expect_equal(ring_cell_side(4), c(1L, 1L, 4L))
  expect_equal(ring_cell_side(5), c(1L, 1L, 1L))
  expect_equal(ring_cell_side(8), c(1L, 4L, 1L))
  expect_equal(ring_cell_side(9), c(1L, 4L, 2L))
  expect_equal(ring_cell_side(12), c(4L, 4L, 2L))
  expect_equal(ring_cell_side(13), c(4L, 4L, 3L))
  expect_equal(ring_cell_side(16), c(4L, 1L, 3L))
})

test_that("next_free_slot avanza en sentido horario y salta ocupadas", {
  ring <- rep("", 16)
  expect_equal(next_free_slot(ring, 16), 1L)
  ring[2:3] <- c("cat", "home")
  expect_equal(next_free_slot(ring, 1), 4L)
})

# Caminos ------------------------------------------------------------------

test_that("con puentes, las filas y columnas son caminos independientes", {
  b <- board_fill(ct$puente)
  expect_equal(nrow(connected_pairs(test_state(b, homes = 12, cats = 1, cat_points = 1))), 1L)
  expect_equal(nrow(connected_pairs(test_state(b, homes = 15, cats = 6, cat_points = 1))), 1L)
  expect_equal(nrow(connected_pairs(test_state(b, homes = 13, cats = 1, cat_points = 1))), 0L)
})

test_that("con rectas verticales solo conectan las columnas", {
  b <- board_fill(ct$recta)
  expect_equal(nrow(connected_pairs(test_state(b, homes = 15, cats = 6, cat_points = 1))), 1L)
  expect_equal(nrow(connected_pairs(test_state(b, homes = 12, cats = 1, cat_points = 1))), 0L)
})

# Inserción ----------------------------------------------------------------

test_that("insertar por arriba empuja la columna y la carta expulsada va al fondo del mazo", {
  b <- board_fill(ct$puente)
  b[1, 1, , ] <- ct$curva; b[2, 1, , ] <- ct$recta; b[3, 1, , ] <- ct$te; b[4, 1, , ] <- ct$cruz
  s <- test_state(b, homes = c(1, 9), cats = c(5, 13), cat_points = c(1, 1),
                  deck = list(ct$doble_curva, ct$recta))
  d <- draw_card(s)
  expect_equal(card_type_of(d$card), "doble_curva")
  s2 <- apply_insertion(d$state, d$card, pos = 5, angle = 0)
  col1 <- vapply(1:4, function(i) card_type_of(s2$board[i, 1, , ]), character(1))
  expect_equal(col1, c("doble_curva", "curva", "recta", "te"))
  expect_equal(vapply(s2$deck, card_type_of, character(1)), c("recta", "cruz"))
})

test_that("insertar por la izquierda y por abajo expulsa la carta del extremo opuesto", {
  b <- board_fill(ct$puente)
  b[4, 4, , ] <- ct$cruz; b[1, 4, , ] <- ct$te
  s <- test_state(b, homes = c(1, 9), cats = c(5, 13), cat_points = c(1, 1), deck = list())
  s_left <- apply_insertion(s, ct$recta, pos = 1, angle = 0)     # fila 4, entra por la izquierda
  expect_equal(card_type_of(s_left$board[4, 1, , ]), "recta")
  expect_equal(card_type_of(s_left$deck[[1]]), "cruz")
  s_down <- apply_insertion(s, ct$recta, pos = 13, angle = 0)    # columna 4, entra por abajo
  expect_equal(card_type_of(s_down$board[4, 4, , ]), "recta")
  expect_equal(card_type_of(s_down$board[3, 4, , ]), "cruz")
  expect_equal(card_type_of(s_down$deck[[1]]), "te")
})

# Capturas -----------------------------------------------------------------

# Columna 1 conectada de arriba (posición 5) a la izquierda de la fila 4 (posición 1)
column1_path_board <- function() {
  b <- board_fill(ct$puente)
  for (i in 1:3) b[i, 1, , ] <- ct$recta
  b[4, 1, , ] <- rotate_card(ct$curva, 270)  # O-N
  b
}

test_that("captura simple: el hogar avanza y entra un gato nuevo", {
  s <- test_state(column1_path_board(), homes = c(1, 9), cats = c(5, 13),
                  cat_points = c(3, 2), cat_deck = c(7, 8))
  s2 <- resolve_captures(s)
  expect_equal(s2$scores, 3)
  expect_equal(s2$cats_won, 1L)
  expect_equal(which(s2$ring_type == "home"), c(2L, 9L))
  expect_equal(which(s2$ring_type == "cat"), c(6L, 13L))
  expect_equal(s2$ring_points[6], 7)
  expect_equal(s2$cat_deck, 8)
})

test_that("dos gatos conectados al mismo hogar se ganan en la misma acción", {
  b <- column1_path_board()
  b[1, 1, , ] <- rotate_card(ct$te, 180)  # S-O-N
  s <- test_state(b, homes = c(1, 9), cats = c(4, 5), cat_points = c(3, 2), cat_deck = c(7, 8))
  s2 <- resolve_captures(s)
  expect_equal(s2$scores, 5)
  expect_equal(s2$cats_won, 2L)
  expect_equal(which(s2$ring_type == "home"), c(3L, 9L))  # el hogar avanzó dos veces
  expect_equal(which(s2$ring_type == "cat"), c(5L, 6L))
  expect_length(s2$cat_deck, 0)
})

test_that("reacción en cadena: el gato nuevo queda conectado al otro hogar", {
  b <- column1_path_board()
  b[1, 2, , ] <- ct$curva  # N-E: une la posición 6 con la fila 1 hacia el hogar 9
  s <- test_state(b, homes = c(1, 9), cats = c(5, 13), cat_points = c(3, 2), cat_deck = c(7, 8))
  expect_equal(nrow(connected_pairs(s)), 1L)
  s2 <- resolve_captures(s)
  expect_equal(s2$cats_won, 2L)
  expect_equal(s2$scores, 3 + 7)
  expect_equal(which(s2$ring_type == "home"), c(2L, 10L))
  expect_equal(which(s2$ring_type == "cat"), c(7L, 13L))
})

# Turnos -------------------------------------------------------------------

always_insert <- list(
  choose_action = function(state, allow_insert) list(type = "insert"),
  choose_insertion = function(state, card) list(pos = 5, angle = 0)
)

test_that("con rotación obligatoria no se puede insertar dos veces", {
  s <- new_game(game_params(rotation_mandatory = TRUE), seed = 3)
  expect_error(play_turn(s, always_insert), "debía rotar")
})

test_that("sin rotación obligatoria se puede insertar dos veces", {
  s <- new_game(game_params(rotation_mandatory = FALSE), seed = 3)
  s2 <- play_turn(s, always_insert)
  expect_length(s2$deck, 10)
  expect_equal(s2$turn, 1L)
})

# Preparación y partida completa -------------------------------------------

test_that("new_game respeta el reglamento y la composición del mazo", {
  params <- game_params()
  s <- new_game(params, seed = 11)
  expect_length(s$deck, 10)
  expect_length(s$cat_deck, 13)
  expect_equal(card_type_of(s$board[1, 4, , ]), "curva")
  expect_equal(card_type_of(s$board[4, 1, , ]), "curva")
  expect_equal(which(s$ring_type == "home"), c(1L, 9L))
  expect_equal(which(s$ring_type == "cat"), c(5L, 13L))
  expect_equal(nrow(connected_pairs(s)), 0L)
  all_types <- table(c(types_of_board(s$board), vapply(s$deck, card_type_of, character(1))))
  expect_equal(as.vector(all_types[names(params$deck_counts)]), unname(params$deck_counts))
})

test_that("las cartas se conservan durante una partida", {
  params <- game_params(rotation_mandatory = FALSE)
  s <- new_game(params, seed = 5)
  set.seed(5)
  pol <- policy_random(p_insert = 0.8)
  for (k in 1:30) if (!game_over(s)) s <- play_turn(s, pol)
  all_types <- table(c(types_of_board(s$board), vapply(s$deck, card_type_of, character(1))))
  expect_equal(as.vector(all_types[names(params$deck_counts)]), unname(params$deck_counts))
  expect_equal(sum(s$cats_won) + sum(s$ring_type == "cat") + length(s$cat_deck), 15)
})

test_that("game_params valida los parámetros", {
  expect_error(game_params(deck_counts = c(curva = 20, ovalo = 6)), "card_types")
  expect_error(game_params(deck_counts = c(curva = 1, recta = 25)), "2 curvas")
  expect_error(game_params(deck_counts = c(curva = 10)), "al menos 17")
  expect_silent(game_params(cat_points = 1:15, rotation_mandatory = FALSE))
})

test_that("una partida en solitario con la política codiciosa rescata los 15 gatos", {
  res <- play_game(game_params(cat_points = 1:15), policy_greedy(), seed = 1)
  expect_true(res$completed)
  expect_equal(res$cats_won, 15L)
  expect_equal(res$scores, sum(1:15))
})
