##########################################
# Motor rápido en C++ (src/engine.cpp) para evaluar jugadas
# Requiere functions.R, cards.R y game.R
#
# El tablero se representa como 16 códigos de carta (celda k = (i - 1) * 4 + j,
# por filas). Código 0..23 = (tipo - 1) * 4 + r, con la matriz de
# rotate_card(card_types[[tipo]], 90 * r).
##########################################

# Tabla de códigos: matrices, código canónico (primera matriz igual) y clave
fe_code_table <- function() {
  mats <- list()
  for (t in seq_along(card_types)) for (r in 0:3)
    mats[[length(mats) + 1L]] <- rotate_card(card_types[[t]], 90 * r)
  keys <- vapply(mats, function(m) paste(m, collapse = ""), character(1))
  canon <- match(keys, keys) - 1L
  conn <- array(0L, c(24, 4, 4))
  for (c in seq_along(mats)) conn[c, , ] <- mats[[c]]
  list(mats = mats, keys = keys, canon = canon, conn = conn,
       lookup = setNames(canon, keys))
}

# Si `make` no está en el PATH, usa el Rtools de RTOOLS45_HOME o RTOOLS44_HOME
# (R 4.6 compila con el Rtools 4.4 instalado en C:/RBuildTools/4.4)
use_rtools <- function() {
  if (.Platform$OS.type != "windows" || nzchar(Sys.which("make"))) return(invisible(FALSE))
  home <- Filter(function(h) nzchar(h) && dir.exists(h),
                 Sys.getenv(c("RTOOLS45_HOME", "RTOOLS44_HOME")))
  if (length(home) == 0) stop("No se encontró Rtools: instale Rtools o defina RTOOLS45_HOME.")
  home <- normalizePath(home[[1]], winslash = "/")
  Sys.setenv(RTOOLS45_HOME = home,
             PATH = paste(file.path(home, "usr", "bin"),
                          file.path(home, "x86_64-w64-mingw32.static.posix", "bin"),
                          Sys.getenv("PATH"), sep = ";"))
  invisible(TRUE)
}

# Compila (o reutiliza de la caché) el motor e inicializa sus tablas
load_fast_engine <- function(root = getwd(), cache_dir = file.path(root, ".rcpp_cache")) {
  use_rtools()
  Rcpp::sourceCpp(file.path(root, "src", "engine.cpp"), cacheDir = cache_dir,
                  env = globalenv())
  FE <<- fe_code_table()
  fe_init(FE$conn, FE$canon)
  invisible(TRUE)
}

card_code <- function(card) unname(FE$lookup[paste(card, collapse = "")])

board_codes <- function(board) {
  out <- integer(BOARD_N^2)
  for (i in seq_len(BOARD_N)) for (j in seq_len(BOARD_N))
    out[(i - 1L) * BOARD_N + j] <- card_code(board[i, j, , ])
  out
}

ring_codes <- function(ring_type) match(ring_type, c("", "cat", "home")) - 1L

type_index <- function(card) card_code(card) %/% 4L

# Estado para el motor. belief = "exacta": las primeras n_unknown cartas del
# mazo son desconocidas (solo se usa su multiconjunto de tipos) y el resto es
# la secuencia conocida de cartas expulsadas. belief = "mazo": todo el mazo se
# trata como desconocido, como en policy_greedy().
fe_state <- function(state, belief = c("exacta", "mazo")) {
  belief <- match.arg(belief)
  s <- masked_state(state)
  types <- vapply(s$deck, type_index, integer(1))
  n_unk <- if (belief == "mazo" || is.null(s$n_unknown)) length(types) else s$n_unknown
  unk_types <- types[seq_len(n_unk)]
  list(codes = board_codes(s$board), ring = ring_codes(s$ring_type),
       pts = s$ring_points, cat_deck = as.numeric(s$cat_deck),
       to_win = as.integer(s$n_cats - sum(s$cats_won)),
       unk = tabulate(unk_types + 1L, nbins = length(card_types)),
       known = types[setdiff(seq_along(types), seq_len(n_unk))],
       use_points = s$n_players > 1L)
}
