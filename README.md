# cat-path-optimizer

Motor de reglas y simulación en R para *No sin mi gato* (*Chabyrinthe*, Antoine Bauza, Cocktail Games, 2007), un juego de cartas en el que se rotan e insertan cartas de canaletas para llevar gatos callejeros hasta un hogar. El objetivo del repositorio es explorar heurísticas y algoritmos que mejoren las decisiones de juego.

## Reglas implementadas

Según el reglamento oficial:

- **Tablero:** 4×4 cartas de canaleta al azar, con curvas obligatorias en las esquinas superior derecha e inferior izquierda. Las cartas sobrantes forman el mazo, boca abajo.
- **Anillo exterior:** 16 posiciones alrededor del tablero (4 por lado). Al empezar hay un gato arriba de la columna 1, un hogar a la derecha de la fila 1, un gato abajo de la columna 4 y un hogar a la izquierda de la fila 4.
- **Turno: 2 acciones.**
  - *Rotar* una carta 90°, 180° o 270°.
  - *Insertar:* se roba la carta de arriba del mazo y se mete por una de las 16 posiciones del borde, empujando la fila o la columna. La carta que sale por el extremo opuesto va al fondo del mazo.
- **Rotación obligatoria (parámetro):** el reglamento original exige que una de las dos acciones sea una rotación. Algunas ediciones permiten repetir la misma acción.
- **Captura:** después de cada acción se ganan todos los gatos conectados a un hogar por un camino continuo de canaletas. Por cada gato ganado, su hogar avanza un lugar en sentido horario y entra un gato nuevo un lugar más allá del gato ganado; si la posición está ocupada, se salta. Después se vuelve a revisar el tablero, porque puede haber reacciones en cadena.
- **Fin:** cuando se ganan todos los gatos. En solitario se cuentan los turnos (escala del reglamento: 7 o menos es la mejor marca).

## Convenciones

- **Cartas:** matrices 4×4 de conexiones entre lados, numerados en sentido horario: 1 = N, 2 = E, 3 = S, 4 = O. Un lado está abierto si conecta con otro lado de la carta.
- **Rotación:** `rotate_card(card, 90)` rota en sentido horario.
- **Anillo:** numerado en sentido horario, con la misma numeración que `insert_new_card()`. Las posiciones 1–4 son la izquierda (filas 4 a 1), 5–8 arriba (columnas 1 a 4), 9–12 la derecha (filas 1 a 4) y 13–16 abajo (columnas 4 a 1). Insertar en la posición *p* mete la carta por esa casilla del anillo.

## Parámetros

```r
params <- game_params(
  deck_counts = c(curva = 8, recta = 6, te = 4, doble_curva = 3, puente = 3, cruz = 2),
  cat_points = rep(1, 15),
  rotation_mandatory = TRUE
)
```

| Parámetro | Valor por defecto | Estado |
|---|---|---|
| `deck_counts` | 26 cartas en los 6 tipos de `cards.R` | **Supuesto.** El juego trae 26 cartas, pero el reglamento no dice cuántas hay de cada tipo. |
| `cat_points` | 15 gatos de 1 punto | **Supuesto.** Cada gato tiene un valor distinto, que el reglamento no detalla. |
| `rotation_mandatory` | `TRUE` | Reglamento original. |
| `corner_elbows` | `TRUE` | Reglamento original. |

## Decisiones donde el reglamento no es explícito

- La carta insertada puede orientarse libremente.
- Si un gato queda conectado a los dos hogares, se usa el primero en el orden del anillo.
- Si el reparto inicial ya tiene un camino gato–hogar, se reparte de nuevo.
- Las curvas obligatorias de las esquinas tienen orientación al azar.

## Políticas

- `policy_random()`: acciones al azar. Sirve de referencia.
- `policy_greedy()`: evalúa todas las jugadas de una acción y elige la que más gatos gana de inmediato (puntos, si hay varios jugadores). Desempata por el tamaño de la red conectada a los gatos. Para decidir si inserta, estima el valor esperado sobre las cartas que quedan en el mazo, sin mirar su orden. No usa información oculta.

Una política es una lista con `choose_action(state, allow_insert)` y `choose_insertion(state, card)`. Para probar una estrategia nueva basta con escribir esas dos funciones.

## Resultados con los parámetros por defecto

30 partidas en solitario (`main.R`, semilla 2016):

| Política | Partidas completas | Turnos: media | Mediana | Mejor | Peor |
|---|---|---|---|---|---|
| Aleatoria (tope de 400 turnos) | 90 % | 266 | 270 | 116 | 400 |
| Codiciosa | 100 % | 10,6 | 11 | 5 | 15 |
| Codiciosa sin rotación obligatoria | 100 % | 8,6 | 9 | 6 | 12 |

Permitir dos inserciones por turno ahorra unos dos turnos por partida. Estos números dependen de la composición supuesta del mazo.

## Uso

```r
source("functions.R"); source("cards.R"); source("game.R"); source("policies.R")
res <- play_game(game_params(), policy_greedy(), n_players = 1, seed = 1)
res$turns
```

`main.R` ejecuta una partida de ejemplo y las comparaciones de arriba; tarda unos 4 minutos.

## Pruebas

```bash
Rscript -e 'testthat::test_file("tests/test_game.R")'
```

## Archivos

| Archivo | Contenido |
|---|---|
| `cards.R` | Los seis tipos de carta de canaleta. |
| `functions.R` | Utilidades: matrices, álgebra circular, inserción y rotación sobre el tablero. |
| `game.R` | Motor de reglas: preparación, caminos, acciones, capturas, turnos y partida. |
| `policies.R` | Políticas de juego y simulación en solitario. |
| `main.R` | Ejemplo de uso y comparación de políticas. |
| `tests/test_game.R` | Pruebas del motor. |
