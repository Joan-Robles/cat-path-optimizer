
# Parameters --------------------------------------------------------------
rm(list = ls()) # limpiar
game_size <- 4 # Columns and rows

library("dplyr")
source("./functions.R")
source("./cards.R")

# random objs -------------------------------------------------------------


# Random card 
card <- random_sym_matrix(game_size)

# Turn it

card90 <- turn_90(card)
card180 <- turn_180(card)
card270 <- turn_270(card)


# Random board
board <- random_board(game_size)


