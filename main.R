# Parameters --------------------------------------------------------------
rm(list = ls()) # limpiar
game_size <- 4 # Columns and rows

library("dplyr")
source("./functions.R")
source("./cards.R")
source("./game.R")
source("./policies.R")

# Parámetros del juego ----------------------------------------------------
# La composición del mazo y los puntos de los gatos son supuestos: el
# reglamento no los trae. Reemplazar por el conteo de las cartas reales.
params <- game_params(
  deck_counts = c(curva = 8, recta = 6, te = 4, doble_curva = 3, puente = 3, cruz = 2),
  cat_points = rep(1, 15),
  rotation_mandatory = TRUE # reglamento original; FALSE = se puede insertar dos veces
)

# random objs -------------------------------------------------------------

# Random card
card <- random_sym_matrix(game_size)

# Random board
small_board <- random_board(game_size = game_size)

# With all borders
board <- random_board(game_size = game_size, big = TRUE)

# Una partida en solitario ------------------------------------------------

partida <- play_game(params, policy_greedy(), n_players = 1, seed = 1)
partida$turns
solo_tier(partida$turns)

# Comparación de políticas en solitario -----------------------------------

resumen <- function(sim) {
  sim %>%
    group_by(policy) %>%
    summarise(partidas = n(), completas = mean(completed),
              turnos_media = mean(turns), turnos_mediana = median(turns),
              mejor = min(turns), peor = max(turns), .groups = "drop")
}

sim <- simulate_solo(params,
                     list(aleatoria = policy_random(), codiciosa = policy_greedy()),
                     n_games = 30, seed = 2016, max_turns = 400)
resumen(sim)
table(sim$policy, solo_tier(sim$turns))

# Efecto de la rotación obligatoria ---------------------------------------

params_libre <- params
params_libre$rotation_mandatory <- FALSE
sim_libre <- simulate_solo(params_libre, list(codiciosa_libre = policy_greedy()),
                           n_games = 30, seed = 2016)
resumen(bind_rows(sim %>% filter(policy == "codiciosa"), sim_libre))
