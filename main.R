
# Parameters --------------------------------------------------------------
rm(list = ls()) # limpiar
game_size <- 4 # Columns and rows

library("dplyr")
source("./functions.R")
source("./cards.R")

# random objs -------------------------------------------------------------


# Random card 
card <- random_sym_matrix(game_size)


# Random board
small_board <- random_board(game_size = game_size)

# With all borders
board <- random_board(game_size = game_size)


