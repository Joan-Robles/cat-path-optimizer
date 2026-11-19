##########################################
# Motor de reglas de "No sin mi gato" (Chabyrinthe, Cocktail Games)
# Requiere functions.R y cards.R
#
# Convenciones (comentarios en español):
#  - Lados de cada carta en sentido horario: 1 = N, 2 = E, 3 = S, 4 = O.
#  - board[i, j, , ] es la matriz 4x4 de conexiones entre lados de la carta
#    en la fila i, columna j.
#  - Anillo exterior de 16 posiciones, en sentido horario, con la misma
#    numeración que insert_new_card():
#      1:4   izquierda, filas 4..1
#      5:8   arriba,    columnas 1..4
#      9:12  derecha,   filas 1..4
#      13:16 abajo,     columnas 4..1
#    Insertar en la posición p mete la carta por esa casilla del anillo.
##########################################

N_SIDE <- 4L
BOARD_N <- 4L
RING_N <- 16L

##########################################
# Parámetros de la partida
# deck_counts: cartas de cada tipo de card_types (26 en el juego real)
# cat_points: puntos de cada carta de gato (15 en el juego real)
# rotation_mandatory: TRUE si una de las 2 acciones debe ser una rotación
#   (reglamento original); FALSE si se puede insertar dos veces
# corner_elbows: curvas obligatorias en las esquinas superior derecha e
#   inferior izquierda al preparar la partida
##########################################

game_params <- function(deck_counts = c(curva = 8, recta = 6, te = 4,
                                        doble_curva = 3, puente = 3, cruz = 2),
                        cat_points = rep(1, 15),
                        rotation_mandatory = TRUE,
                        corner_elbows = TRUE) {
  if (is.null(names(deck_counts)) || !all(names(deck_counts) %in% names(card_types)))
    stop("`deck_counts` debe tener nombres de card_types: ", paste(names(card_types), collapse = ", "))
  if (any(deck_counts < 0) || any(deck_counts != round(deck_counts)))
    stop("`deck_counts` debe contener enteros no negativos.")
  # Se necesitan 16 cartas en el tablero y al menos 1 en el mazo para insertar
  if (sum(deck_counts) < BOARD_N^2 + 1L)
    stop("El mazo necesita al menos ", BOARD_N^2 + 1L, " cartas.")
  if (corner_elbows && (is.na(deck_counts["curva"]) || deck_counts[["curva"]] < 2))
    stop("Con `corner_elbows = TRUE` el mazo necesita al menos 2 curvas.")
  if (length(cat_points) < 2 || any(!is.finite(cat_points)))
    stop("`cat_points` necesita al menos 2 valores finitos.")

  list(deck_counts = deck_counts, cat_points = cat_points,
       rotation_mandatory = isTRUE(rotation_mandatory),
       corner_elbows = isTRUE(corner_elbows))
}

##########################################
# Cartas: rotación, lados abiertos y tipo
##########################################

# Rota una carta 0, 90, 180 o 270 grados en sentido horario
rotate_card <- function(card, angle) {
  angle <- angle %% 360
  if (angle == 0) return(card)
  dir <- map_angle(angle)
  shift_matrix(card, di = dir, dj = dir)
}

# Un lado está abierto si conecta con al menos otro lado de la carta
open_sides <- function(card) {
  rowSums(card) - diag(card) > 0
}

# Orientaciones distintas de una carta (la cruz tiene 1, la recta 2, etc.)
distinct_rotations <- function(card) {
  angles <- c(0, 90, 180, 270)
  rots <- lapply(angles, function(a) rotate_card(card, a))
  keep <- !duplicated(vapply(rots, function(m) paste(m, collapse = ""), character(1)))
  angles[keep]
}

# Tipo de una carta según card_types (salvo rotación)
card_type_of <- function(card) {
  for (nm in names(card_types)) {
    for (a in c(0, 90, 180, 270)) {
      if (all(rotate_card(card_types[[nm]], a) == card)) return(nm)
    }
  }
  NA_character_
}

##########################################
# Anillo exterior
##########################################

# Celda y lado del tablero contiguos a la posición pos del anillo: c(i, j, lado)
ring_cell_side <- function(pos) {
  grp <- (pos - 1L) %/% 4L + 1L
  k <- (pos - 1L) %% 4L + 1L
  switch(grp,
         c(5L - k, 1L, 4L),  # izquierda, filas 4..1, lado O
         c(1L, k, 1L),       # arriba, columnas 1..4, lado N
         c(k, 4L, 2L),       # derecha, filas 1..4, lado E
         c(4L, 5L - k, 3L))  # abajo, columnas 4..1, lado S
}

# Siguiente posición libre en sentido horario, saltando las ocupadas
next_free_slot <- function(ring_type, from) {
  pos <- from
  for (step in seq_len(RING_N)) {
    pos <- loop_add(pos, 1L, RING_N)
    if (ring_type[pos] == "") return(pos)
  }
  stop("No hay posiciones libres en el anillo.")
}

##########################################
# Conectividad: grafo de lados
# Nodo = (celda, lado). Aristas: lados conectados dentro de una carta y
# lados abiertos enfrentados entre cartas vecinas.
##########################################

node_id <- function(i, j, side) ((i - 1L) * BOARD_N + (j - 1L)) * N_SIDE + side

board_graph <- function(board) {
  n_nodes <- BOARD_N^2 * N_SIDE
  A <- matrix(FALSE, n_nodes, n_nodes)
  open <- array(FALSE, c(BOARD_N, BOARD_N, N_SIDE))

  for (i in seq_len(BOARD_N)) for (j in seq_len(BOARD_N)) {
    M <- board[i, j, , ]
    ids <- node_id(i, j, 1:4)
    A[ids, ids] <- M == 1L
    open[i, j, ] <- open_sides(M)
  }
  # Vecinos horizontales: E (2) de (i, j) con O (4) de (i, j + 1)
  for (i in seq_len(BOARD_N)) for (j in seq_len(BOARD_N - 1L)) {
    if (open[i, j, 2] && open[i, j + 1L, 4]) {
      a <- node_id(i, j, 2L); b <- node_id(i, j + 1L, 4L)
      A[a, b] <- A[b, a] <- TRUE
    }
  }
  # Vecinos verticales: S (3) de (i, j) con N (1) de (i + 1, j)
  for (i in seq_len(BOARD_N - 1L)) for (j in seq_len(BOARD_N)) {
    if (open[i, j, 3] && open[i + 1L, j, 1]) {
      a <- node_id(i, j, 3L); b <- node_id(i + 1L, j, 1L)
      A[a, b] <- A[b, a] <- TRUE
    }
  }
  list(A = A, open = open)
}

# Nodos alcanzables desde start (vector lógico)
reachable_from <- function(A, start) {
  reach <- logical(nrow(A))
  reach[start] <- TRUE
  repeat {
    new_reach <- reach | colSums(A[reach, , drop = FALSE]) > 0
    if (all(new_reach == reach)) return(reach)
    reach <- new_reach
  }
}

# Nodo de entrada desde una posición del anillo (NA si ese lado está cerrado)
entry_node <- function(graph, pos) {
  cs <- ring_cell_side(pos)
  if (!graph$open[cs[1], cs[2], cs[3]]) return(NA_integer_)
  node_id(cs[1], cs[2], cs[3])
}

# Pares gato-hogar conectados, en orden del anillo
connected_pairs <- function(state, graph = board_graph(state$board)) {
  cats <- which(state$ring_type == "cat")
  homes <- which(state$ring_type == "home")
  home_nodes <- vapply(homes, function(h) entry_node(graph, h), integer(1))
  out <- data.frame(cat = integer(0), home = integer(0))
  for (c_pos in cats) {
    start <- entry_node(graph, c_pos)
    if (is.na(start)) next
    reach <- reachable_from(graph$A, start)
    hit <- homes[!is.na(home_nodes) & reach[ifelse(is.na(home_nodes), 1L, home_nodes)]]
    if (length(hit) > 0) out <- rbind(out, data.frame(cat = c_pos, home = hit[1]))
  }
  out
}

##########################################
# Preparación de la partida
##########################################

build_deck <- function(deck_counts) {
  types <- rep(names(deck_counts), times = deck_counts)
  lapply(types, function(t) card_types[[t]])
}

new_game <- function(params, n_players = 1L, seed = NULL, max_deals = 100L) {
  if (!is.null(seed)) set.seed(seed)

  for (deal in seq_len(max_deals)) {
    deck <- build_deck(params$deck_counts)
    deck <- deck[sample.int(length(deck))]
    board <- array(0L, dim = c(BOARD_N, BOARD_N, N_SIDE, N_SIDE))
    random_rot <- function(card) rotate_card(card, sample(c(0, 90, 180, 270), 1))

    cells <- as.matrix(expand.grid(i = seq_len(BOARD_N), j = seq_len(BOARD_N)))
    if (params$corner_elbows) {
      # Curvas obligatorias en (1, 4) y (4, 1); orientación al azar
      curves <- which(vapply(deck, card_type_of, character(1)) == "curva")[1:2]
      board[1, BOARD_N, , ] <- random_rot(deck[[curves[1]]])
      board[BOARD_N, 1, , ] <- random_rot(deck[[curves[2]]])
      deck <- deck[-curves]
      cells <- cells[!((cells[, 1] == 1 & cells[, 2] == BOARD_N) |
                         (cells[, 1] == BOARD_N & cells[, 2] == 1)), , drop = FALSE]
    }
    for (r in seq_len(nrow(cells))) {
      board[cells[r, 1], cells[r, 2], , ] <- random_rot(deck[[1]])
      deck <- deck[-1]
    }

    cat_deck <- sample(params$cat_points)
    ring_type <- rep("", RING_N)
    ring_points <- rep(NA_real_, RING_N)
    # Esquema 1 del reglamento: gato arriba de la columna 1 (5), hogar a la
    # derecha de la fila 1 (9), gato abajo de la columna 4 (13), hogar a la
    # izquierda de la fila 4 (1)
    ring_type[c(1, 9)] <- "home"
    ring_type[c(5, 13)] <- "cat"
    ring_points[c(5, 13)] <- cat_deck[1:2]
    cat_deck <- cat_deck[-(1:2)]

    # n_unknown: cartas de arriba del mazo que nadie ha visto. Las expulsadas
    # van al fondo y su tipo es conocido.
    state <- list(board = board, deck = deck, n_unknown = length(deck), cat_deck = cat_deck,
                  ring_type = ring_type, ring_points = ring_points,
                  scores = numeric(n_players), cats_won = integer(n_players),
                  n_players = as.integer(n_players), player = 1L, turn = 0L,
                  n_cats = length(params$cat_points), params = params,
                  history = list())
    # El reglamento no cubre un camino ya formado al inicio: se reparte de nuevo
    if (nrow(connected_pairs(state)) == 0) return(state)
  }
  stop("No se encontró un reparto inicial sin caminos tras ", max_deals, " intentos.")
}

##########################################
# Acciones
##########################################

# Una rotación gira la carta un cuarto de vuelta: 90 (derecha) o 270
# (izquierda). Girarla 180 grados cuesta las dos acciones del turno.
ROTATION_ANGLES <- c(90, 270)

apply_rotation <- function(state, cell, angle) {
  stopifnot(length(cell) == 2, all(cell %in% seq_len(BOARD_N)), angle %in% ROTATION_ANGLES)
  state$board <- shift_one_new_card(state$board, cordenates = cell, amount = angle)
  state
}

# Carta que sale del tablero al insertar por la posición pos
expelled_card <- function(board, pos) {
  grp <- (pos - 1L) %/% 4L + 1L
  k <- (pos - 1L) %% 4L + 1L
  switch(grp,
         board[5L - k, BOARD_N, , ],  # entra por la izquierda, sale por la derecha
         board[BOARD_N, k, , ],       # entra por arriba, sale por abajo
         board[k, 1L, , ],            # entra por la derecha, sale por la izquierda
         board[1L, 5L - k, , ])       # entra por abajo, sale por arriba
}

draw_card <- function(state) {
  card <- state$deck[[1]]
  state$deck <- state$deck[-1]
  if (!is.null(state$n_unknown) && state$n_unknown > 0) state$n_unknown <- state$n_unknown - 1L
  list(state = state, card = card)
}

# Inserta la carta robada (rotada angle grados) por la posición pos.
# La carta expulsada va al fondo del mazo.
apply_insertion <- function(state, card, pos, angle = 0) {
  stopifnot(pos %in% seq_len(RING_N), angle %in% c(0, 90, 180, 270))
  expelled <- expelled_card(state$board, pos)
  state$board <- insert_new_card(state$board, rotate_card(card, angle), pos)
  state$deck <- c(state$deck, list(expelled))
  state
}

##########################################
# Gatos ganados y reacciones en cadena
##########################################

resolve_captures <- function(state) {
  repeat {
    pairs <- connected_pairs(state)
    if (nrow(pairs) == 0) return(state)
    p <- state$player

    # 1) Todos los gatos conectados tras la acción se ganan a la vez
    for (r in seq_len(nrow(pairs))) {
      cat_pos <- pairs$cat[r]
      state$scores[p] <- state$scores[p] + state$ring_points[cat_pos]
      state$cats_won[p] <- state$cats_won[p] + 1L
      state$history[[length(state$history) + 1L]] <-
        list(turn = state$turn, player = p, cat = cat_pos, home = pairs$home[r],
             points = state$ring_points[cat_pos])
      state$ring_type[cat_pos] <- ""
      state$ring_points[cat_pos] <- NA_real_
    }

    # 2) Por cada gato ganado, su hogar avanza un lugar en sentido horario y
    #    entra un gato nuevo un lugar más allá del gato ganado. Si dos gatos
    #    usaron el mismo hogar, este avanza dos veces.
    home_now <- pairs$home
    for (r in seq_len(nrow(pairs))) {
      h <- home_now[r]
      state$ring_type[h] <- ""
      new_h <- next_free_slot(state$ring_type, h)
      state$ring_type[new_h] <- "home"
      home_now[home_now == h] <- new_h

      if (length(state$cat_deck) > 0) {
        new_cat <- next_free_slot(state$ring_type, pairs$cat[r])
        state$ring_type[new_cat] <- "cat"
        state$ring_points[new_cat] <- state$cat_deck[1]
        state$cat_deck <- state$cat_deck[-1]
      }
    }
    # 3) Se vuelve a revisar: los hogares y gatos desplazados pueden formar
    #    caminos nuevos (reacciones en cadena)
  }
}

game_over <- function(state) sum(state$cats_won) >= state$n_cats

##########################################
# Turnos y partida
# Una política es una lista con dos funciones:
#  - choose_action(state, allow_insert): devuelve
#      list(type = "rotate", cell = c(i, j), angle = 90 | 270)  o
#      list(type = "insert")
#  - choose_insertion(state, card): devuelve list(pos = 1..16, angle = 0..270)
#    (se llama después de robar la carta, como en el reglamento)
# state$action indica si es la primera o la segunda acción del turno.
##########################################

play_turn <- function(state, policy) {
  state$turn <- state$turn + 1L
  used_insert <- FALSE
  for (a in 1:2) {
    allow_insert <- !(state$params$rotation_mandatory && used_insert)
    state$action <- a  # acción del turno (1 o 2), para políticas que planean el turno
    choice <- policy$choose_action(state, allow_insert)
    if (identical(choice$type, "insert")) {
      if (!allow_insert) stop("La política insertó cuando debía rotar.")
      d <- draw_card(state)
      ins <- policy$choose_insertion(d$state, d$card)
      state <- apply_insertion(d$state, d$card, ins$pos, ins$angle)
      used_insert <- TRUE
    } else {
      state <- apply_rotation(state, choice$cell, choice$angle)
    }
    state <- resolve_captures(state)
    if (game_over(state)) break
  }
  state$player <- loop_add(state$player, 1L, state$n_players)
  state
}

# policies: una política o una lista con una por jugador
play_game <- function(params, policies, n_players = 1L, seed = NULL, max_turns = 200L) {
  state <- new_game(params, n_players = n_players, seed = seed)
  if (!is.null(policies$choose_action)) policies <- rep(list(policies), n_players)
  while (!game_over(state) && state$turn < max_turns) {
    state <- play_turn(state, policies[[state$player]])
  }
  list(state = state, turns = state$turn, completed = game_over(state),
       scores = state$scores, cats_won = state$cats_won)
}

# Escala del modo solitario según el número de turnos
solo_tier <- function(turns) {
  cut(turns, breaks = c(-Inf, 7, 9, 11, 13, 15, 17, Inf),
      labels = c("<=7", "8-9", "10-11", "12-13", "14-15", "16-17", ">=18"))
}
