
# List the six kind of cards

# Rotationally unique()

card_1 <- rbind(c(1, 1, 0, 0),
                c(1, 1, 0, 0),
                c(0, 0, 1, 0),
                c(0, 0, 0, 1))

card_2 <- rbind(c(1, 0, 1, 0),
                c(0, 1, 0, 0),
                c(1, 0, 1, 0),
                c(0, 0, 0, 1))

card_3 <- rbind(c(1, 1, 1, 0),
                c(1, 1, 1, 0),
                c(1, 1, 1, 0),
                c(0, 0, 0, 1))

card_4 <- rbind(c(1, 1, 0, 0),
                c(1, 1, 0, 0),
                c(0, 0, 1, 1),
                c(0, 0, 1, 1))

card_5 <- rbind(c(1, 0, 1, 0),
                c(0, 1, 0, 1),
                c(1, 0, 1, 0),
                c(0, 1, 0, 1))

card_6 <- rbind(c(1, 1, 1, 1),
                c(1, 1, 1, 1),
                c(1, 1, 1, 1),
                c(1, 1, 1, 1))

# Nombres de los tipos, en el orden en que se dan las cantidades del mazo.
# Lados en sentido horario: 1 = N, 2 = E, 3 = S, 4 = O
card_types <- list(
  curva       = card_1, # N-E
  recta       = card_2, # N-S
  te          = card_3, # N-E-S
  doble_curva = card_4, # N-E y S-O
  puente      = card_5, # N-S y E-O sin conexión entre sí
  cruz        = card_6  # todos conectados
)
card_types <- lapply(card_types, function(m) { storage.mode(m) <- "integer"; m })

