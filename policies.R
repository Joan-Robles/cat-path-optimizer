##########################################
# Políticas de juego y simulación
# Requiere game.R
##########################################

# Índice del máximo con desempate al azar
pick_best <- function(values) {
  idx <- which(values >= max(values) - 1e-9)
  idx[sample.int(length(idx), 1L)]
}

##########################################
# Política aleatoria: referencia para comparar
# p_insert: probabilidad de insertar cuando está permitido
##########################################

policy_random <- function(p_insert = 0.5) {
  list(
    choose_action = function(state, allow_insert) {
      if (allow_insert && runif(1) < p_insert) return(list(type = "insert"))
      list(type = "rotate", cell = sample(seq_len(BOARD_N), 2, replace = TRUE),
           angle = sample(c(90, 180, 270), 1))
    },
    choose_insertion = function(state, card) {
      list(pos = sample(seq_len(RING_N), 1), angle = sample(distinct_rotations(card), 1))
    }
  )
}

##########################################
# Política codiciosa (miope, una acción a la vez)
# - Evalúa todas las rotaciones distintas y elige la que más gatos (solitario)
#   o puntos (varios jugadores) gana de inmediato.
# - Desempata con el alcance de la red conectada a los gatos.
# - Para decidir si inserta, estima el valor esperado de insertar como el
#   promedio, sobre las cartas que quedan en el mazo, del mejor valor de
#   inserción de cada tipo. No mira el orden del mazo (está boca abajo).
##########################################

# Copia para evaluar jugadas sin información oculta: los puntos de los gatos
# no revelados se reemplazan por su promedio
masked_state <- function(state) {
  if (length(state$cat_deck) > 0)
    state$cat_deck <- rep(mean(state$cat_deck), length(state$cat_deck))
  state$history <- list()
  state
}

# Número de nodos alcanzables desde los gatos (heurística de desempate)
cat_reach <- function(state, graph = board_graph(state$board)) {
  total <- 0
  for (c_pos in which(state$ring_type == "cat")) {
    start <- entry_node(graph, c_pos)
    if (!is.na(start)) total <- total + sum(reachable_from(graph$A, start))
  }
  total
}

state_value <- function(before, after) {
  p <- before$player
  gain <- if (before$n_players == 1L) {
    sum(after$cats_won) - sum(before$cats_won)
  } else {
    after$scores[p] - before$scores[p]
  }
  gain + 0.001 * cat_reach(after)
}

rotation_candidates <- function(state) {
  out <- list()
  for (i in seq_len(BOARD_N)) for (j in seq_len(BOARD_N)) {
    card <- state$board[i, j, , ]
    seen <- paste(card, collapse = "")
    for (a in c(90, 180, 270)) {
      key <- paste(rotate_card(card, a), collapse = "")
      if (key %in% seen) next
      seen <- c(seen, key)
      out[[length(out) + 1L]] <- list(type = "rotate", cell = c(i, j), angle = a)
    }
  }
  out
}

policy_greedy <- function(consider_insert = TRUE) {
  eval_rotation <- function(state, act) {
    s <- masked_state(state)
    state_value(s, resolve_captures(apply_rotation(s, act$cell, act$angle)))
  }
  eval_insertion <- function(state, card, pos, angle) {
    s <- masked_state(state)
    state_value(s, resolve_captures(apply_insertion(s, card, pos, angle)))
  }
  best_insertion <- function(state, card) {
    grid <- expand.grid(pos = seq_len(RING_N), angle = distinct_rotations(card))
    vals <- mapply(function(p, a) eval_insertion(state, card, p, a), grid$pos, grid$angle)
    k <- pick_best(vals)
    list(pos = grid$pos[k], angle = grid$angle[k], value = vals[k])
  }

  list(
    choose_action = function(state, allow_insert) {
      cands <- rotation_candidates(state)
      vals <- vapply(cands, function(act) eval_rotation(state, act), numeric(1))
      k <- pick_best(vals)
      if (allow_insert && consider_insert) {
        types <- table(vapply(state$deck, card_type_of, character(1)))
        best_by_type <- vapply(names(types),
                               function(t) best_insertion(state, card_types[[t]])$value,
                               numeric(1))
        expected <- sum(best_by_type * types) / sum(types)
        if (expected > vals[k] + 1e-9) return(list(type = "insert"))
      }
      cands[[k]]
    },
    choose_insertion = function(state, card) {
      b <- best_insertion(state, card)
      list(pos = b$pos, angle = b$angle)
    }
  )
}

##########################################
# Simulación de partidas en solitario
# Cada partida g usa la misma semilla para todas las políticas, de modo que
# todas parten del mismo reparto y del mismo orden de mazos
##########################################

simulate_solo <- function(params, policies, n_games = 20L, seed = 1L, max_turns = 200L) {
  rows <- list()
  for (g in seq_len(n_games)) {
    for (nm in names(policies)) {
      res <- play_game(params, policies[[nm]], n_players = 1L,
                       seed = seed + g, max_turns = max_turns)
      rows[[length(rows) + 1L]] <- data.frame(policy = nm, game = g, turns = res$turns,
                                              completed = res$completed)
    }
  }
  do.call(rbind, rows)
}
