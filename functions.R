##########################################
# Random binary matrix generator
# Generates n x n matrix with 0s and 1s
##########################################

random_sym_matrix <- function(n, p = 0.2) {
  # Validaciones ultrarrápidas
  n <- as.integer(n)
  if (!is.finite(n) || n < 1L) stop("n must be a finite integer >= 1")
  p <- as.numeric(p)
  if (!is.finite(p)) stop("p must be finite")
  p <- if (p < 0) 0 else if (p > 1) 1 else p  # clamp
  
  M <- matrix(0L, n, n)
  k <- n * (n - 1L) / 2L
  if (k > 0L) {
    # Bernoulli mediante runif para evitar warnings de rbinom
    M[upper.tri(M)] <- as.integer(runif(k) < p)
    M <- M + t(M)
  }
  diag(M) <- 1L
  M
}

random_board <- function(game_size, p = 0.3, big = FALSE) {
  
  if(big){
    dims <- game_size + 2
  } else {
    dims <- game_size
  }
  
  arr <- array(0L, dim = c(dims, dims, game_size, game_size))
  
  for (i in 1:dims) {
    for (j in 1:dims) {
      arr[i, j, , ] <- random_sym_matrix(game_size, p)
    }
  }
  
  arr
}

##########################################
# Loop algebra over 1..n (vectorized)
# Suma con wrap-around: 1..n
# Comentarios en español
##########################################

# Normaliza a 1..n (p. ej., 0 -> n, n+1 -> 1)
.loop_norm <- function(x, n) {
  n <- as.integer(n)
  ((as.integer(x) - 1L) %% n) + 1L
}

# Suma en álgebra circular 1..n
loop_add <- function(a, b = 1, n = 4) {
  # a y b pueden ser vectores; se reciclan como en R base
  a2 <- .loop_norm(a, n)        # asegura dominio 1..n
  res <- .loop_norm(a2 + as.integer(b), n)
  res
}


shift_matrix <- function(mat, di = 0L, dj = 0L) {
  n <- nrow(mat)
  m <- ncol(mat)
  # Calcula nuevos índices vectorizados
  i_idx <- ((seq_len(n) - 1L + di) %% n) + 1L
  j_idx <- ((seq_len(m) - 1L + dj) %% m) + 1L
  # Reordena filas y columnas
  mat[i_idx, j_idx, drop = FALSE]
}



find_unique_matrices <- function(lst) {
  if (!is.list(lst)) stop("Input must be a list of matrices")
  
  # Convierte cada matriz en una cadena única (firma)
  keys <- vapply(lst, function(m) paste(c(m), collapse = ","), FUN.VALUE = character(1))
  
  # Encuentra índices únicos preservando el orden
  idx <- !duplicated(keys)
  
  lst[idx]
}

##########################################
# Lower-triangular off-diagonal indices (2 x k)
# Row 1 = i, Row 2 = j
# Comentarios en español
##########################################

lower_off_ij <- function(n) {
  n <- as.integer(n)
  if (!is.finite(n) || n < 1L) stop("n must be a finite integer >= 1")
  
  # Para cada columna j, hay (n - j) filas i en el triángulo inferior (sin diagonal)
  times <- n - seq_len(n)                 # (n-1, n-2, ..., 0)
  j <- rep.int(seq_len(n), times = times) # columnas repetidas
  if (length(j) == 0L) {
    return(matrix(integer(0), nrow = 2L, dimnames = list(c("i","j"), NULL)))
  }
  # i va de (j+1) a n en cada bloque; sequence() genera 1..times[j]
  i <- sequence(times) + rep.int(seq_len(n), times = times)
  
  res <- rbind(i = i, j = j)
  res
}

#############################################
# All 0/1 combinations as a matrix (2^n x n)
# Comments in Spanish
#############################################

all_binary_matrix <- function(n, order = c("lex", "gray"), exclude_zeros = TRUE) {
  # Validaciones rápidas
  n <- as.integer(n)
  if (!is.finite(n) || n < 0L) stop("n must be a finite integer >= 0")
  # Nota: 2^n filas -> crece exponencialmente; cuidado con memoria
  if (n > 30L) warning("n > 30: 2^n filas puede ser demasiado grande en memoria")
  
  order <- match.arg(order)
  rows <- bitwShiftL(1L, n)               # 2^n
  if (n == 0L) return(matrix(integer(0), nrow = 1L, ncol = 0L))
  
  # Valores base 0..2^n-1; opcionalmente en Gray code
  vals <- 0:(rows - 1L)
  if (order == "gray") vals <- bitwXor(vals, bitwShiftR(vals, 1L))
  
  # Máscara por columna (bit más significativo a la izquierda)
  masks <- bitwShiftL(1L, (n - 1L):0L)
  
  # Extrae bits por columna (vectorizado sobre filas; vapply es rápido)
  M <- vapply(
    masks,
    function(m) as.integer(bitwAnd(vals, m) != 0L),
    FUN.VALUE = integer(rows)
  )
  
  dimnames(M) <- NULL
  
  if (exclude_zeros) {
    M <- M[-1, ]
  }
  
  M
}

##########################################
# Insert matrix d and drop one slice in a 3D array
# where="begin": (d, a1, ..., a_{k-1})
# where="final": (a2, ..., ak, d)
# along: dimension (1..3)
# Comentarios en español
##########################################

array_insert_drop <- function(A, d, where = c("begin","final")) {
  stopifnot(is.array(A), length(dim(A)) == 3L)
  where <- match.arg(where)
  dA <- dim(A); k <- dA[1L]; X <- dA[2L]; Y <- dA[3L]
  
  # Normaliza d a 1×X×Y
  if (is.matrix(d)) {
    if (!all(dim(d) == c(X, Y))) stop("`d` debe ser matriz X×Y compatible.")
    d1 <- array(d, dim = c(1L, X, Y))
  } else if (is.array(d) && length(dim(d)) == 3L && all(dim(d) == c(1L, X, Y))) {
    d1 <- d
  } else {
    stop("`d` debe ser matriz X×Y o array 1×X×Y.")
  }
  
  out <- array(vector(typeof(c(A, d1)), length = length(A)), dim = dA)
  
  if (where == "begin") {
    out[1L, , ] <- d1[1L, , ]
    if (k > 1L) out[2L:k, , ] <- A[1L:(k - 1L), , ]
  } else { # "final"
    if (k > 1L) out[1L:(k - 1L), , ] <- A[2L:k, , ]
    out[k, , ] <- d1[1L, , ]
  }
  
  out
}

##########################################
# Insert d (4x4) en board (4x4x4x4) según pos 1..16
# where:
#  - 1:4   -> "left"  (insert al inicio de la dim 1 del bloque de fila i)
#  - 5:8   -> "up"    (insert al inicio de la dim 1 del bloque de columna j)
#  - 9:12  -> "right" (insert al final  de la dim 1 del bloque de fila i)
#  - 13:16 -> "down"  (insert al final  de la dim 1 del bloque de columna j)
# Comentarios en español
##########################################

array_insert_drop <- function(A, d, where = c("begin","final")) {
  stopifnot(is.array(A), length(dim(A)) == 3L)
  where <- match.arg(where)
  dA <- dim(A); k <- dA[1L]; X <- dA[2L]; Y <- dA[3L]
  
  # Normaliza d a 1×X×Y
  if (is.matrix(d)) {
    if (!all(dim(d) == c(X, Y))) stop("`d` debe ser matriz X×Y compatible.")
    d1 <- array(d, dim = c(1L, X, Y))
  } else if (is.array(d) && length(dim(d)) == 3L && all(dim(d) == c(1L, X, Y))) {
    d1 <- d
  } else {
    stop("`d` debe ser matriz X×Y o array 1×X×Y.")
  }
  
  out <- array(vector(typeof(c(A, d1)), length = length(A)), dim = dA)
  
  if (where == "begin") {
    out[1L, , ] <- d1[1L, , ]
    if (k > 1L) out[2L:k, , ] <- A[1L:(k - 1L), , ]
  } else { # "final"
    if (k > 1L) out[1L:(k - 1L), , ] <- A[2L:k, , ]
    out[k, , ] <- d1[1L, , ]
  }
  out
}

insert_new_card <- function(board, new_new_card, pos) {
  # Chequeos rápidos
  stopifnot(is.array(board), length(dim(board)) == 4L)
  if (!all(dim(board)[2:3] == dim(new_new_card)))
    stop("`new_new_card` debe ser matriz con dim(board)[2:3].")
  
  grp <- ((pos - 1L) %/% 4L) + 1L  # 1..4 por bloques de 4
  
  if (grp == 1L) {
    # LEFT: fila i = pos (1..4), inserta al inicio
    i <- 5L - pos
    new_block <- array_insert_drop(A = board[i, , , ],
                                   d = new_new_card,
                                   where = "begin")
    board[i, , , ] <- new_block
    
  } else if (grp == 2L) {
    # UP: columna j = pos-4 (1..4), inserta al inicio
    j <- pos - 4L
    new_block <- array_insert_drop(A = board[, j, , ], d = new_new_card, where = "begin")
    board[, j, , ] <- new_block
    
  } else if (grp == 3L) {
    # RIGHT: fila i = pos-8 (1..4), inserta al final
    i <- pos - 8L
    new_block <- array_insert_drop(A = board[i, , , ],
                                   d = new_new_card,
                                   where = "final")
    board[i, , , ] <- new_block
    
  } else if (grp == 4L) {
    # DOWN: columna j = pos-12 (1..4), inserta al final  <-- CASO 4
    j <- 17L - pos
    new_block <- array_insert_drop(A = board[, j, , ],
                                   d = new_new_card,
                                   where = "final")
    board[, j, , ] <- new_block
    
  } else {
    stop("`pos` debe estar en 1..16.")
  }
  
  board
}

# Mapea 90->-1, 180->2, 270->1 (vectorizado, ultra-rápido; NA para otros)

map_angle <- local({
  f1 <- function(a) switch(as.character(a), "90" = -1L, "180" = 2L, "270" = 1L, NA_integer_)
  function(x) vapply(x, f1, integer(1L))
})

shift_one_new_card <- function(board, cordenates = c(1, 1), amount = c(90, 180, 270)){
  
  # # Desplazamiento a la izquierda (todas las filas y columnas -1)
  # turn_270 <- function(mat) {
  #   shift_matrix(mat, di = 1L, dj = 1L)
  # }
  # 
  # # Desplazamiento a la derecha (todas las filas y columnas +1)
  # turn_90 <- function(mat) {
  #   shift_matrix(mat, di = -1L, dj = -1L)
  # }
  # 
  # turn_180 <- function(mat) {
  #   shift_matrix(mat, di = 2L, dj = 2L)
  # }
  
  sub_board <- board[cordenates[1], cordenates[2], , ]
  
  dir <- map_angle(amount)
  turned_sub_board <- shift_matrix(sub_board, di = dir, dj = dir)
  
  board[cordenates[1], cordenates[2], , ] <- turned_sub_board # change board
  
  board # return
  
}

handle_AB <- function(board, x = NULL, angle = NULL, pos = NULL, new_card = NULL) {
  # Predicados rápidos (enteros y rangos)
  is_int <- function(v) is.numeric(v) && all(!is.na(v)) && all(v == as.integer(v))
  
  is_A <- !is.null(x) && length(x) == 2L && is_int(x) && all(x >= 1L & x <= 4L) &&
    !is.null(angle) && is_int(angle) && length(angle) == 1L && angle %in% c(90L, 180L, 270L) &&
    is.null(pos) && is.null(new_card)
  
  is_B <- !is.null(pos) && is_int(pos) && length(pos) == 1L && pos >= 1L && pos <= 16L &&
    !is.null(new_card) && is.matrix(new_card) && all(dim(new_card) == c(4L, 4L)) &&
    is.null(x) && is.null(angle)
  
  if (is_A) {
    # --- CASE A: x (2D en 1:4) + angle (90/180/270) ---
    result <- shift_one_card(board = board, cordenates = x, amount = angle)
    
  } else if (is_B) {
    # --- CASE B: pos (1..16) + new_card (4x4) ---
    result <-  insert_card(board = board, new_card = new_card, pos = pos)
      
  } else {
    stop("Invalid input: expected A(x, angle) or B(pos, new_card).")
  }
  
  result # return
}

