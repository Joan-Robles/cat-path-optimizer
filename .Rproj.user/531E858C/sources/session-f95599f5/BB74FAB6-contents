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

random_board <- function(game_size, p = 0.2) {
  arr <- array(0L, dim = c(2, 2, game_size, game_size))
  
  for (i in 1:2) {
    for (j in 1:2) {
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

# Desplazamiento a la izquierda (todas las filas y columnas -1)
turn_270 <- function(mat) {
  shift_matrix(mat, di = 1L, dj = 1L)
}

# Desplazamiento a la derecha (todas las filas y columnas +1)
turn_90 <- function(mat) {
  shift_matrix(mat, di = -1L, dj = -1L)
}

turn_180 <- function(mat) {
  shift_matrix(mat, di = 2L, dj = 2L)
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


