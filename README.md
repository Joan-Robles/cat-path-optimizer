# cat-path-optimizer

Motor de reglas y simulación en R para *No sin mi gato* (*Chabyrinthe*, Antoine Bauza, Cocktail Games, 2007), un juego de cartas en el que se rotan e insertan cartas de canaletas para llevar gatos callejeros hasta un hogar. El objetivo del repositorio es explorar heurísticas y algoritmos que mejoren las decisiones de juego.

## Reglas implementadas

Según el reglamento oficial:

- **Tablero:** 4×4 cartas de canaleta al azar, con curvas obligatorias en las esquinas superior derecha e inferior izquierda. Las cartas sobrantes forman el mazo, boca abajo.
- **Anillo exterior:** 16 posiciones alrededor del tablero (4 por lado). Al empezar hay un gato arriba de la columna 1, un hogar a la derecha de la fila 1, un gato abajo de la columna 4 y un hogar a la izquierda de la fila 4.
- **Turno: 2 acciones.**
  - *Rotar* una carta un cuarto de vuelta: 90° (derecha) o 270° (izquierda). Girarla 180° cuesta las dos acciones del turno. En la recta y la doble curva 90° y 270° dan la misma orientación; el puente y la cruz no cambian al girar.
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

- `policy_search(horizon, belief, w_reach)`: usa el motor rápido en C++.
  - `horizon = "turno"`: expectimax sobre las acciones que quedan en el turno, con un nodo de azar por cada carta robada (unas 13 000 jugadas evaluadas por turno). `horizon = "accion"` decide una acción a la vez, como la codiciosa.
  - `belief = "exacta"`: la carta robada sale de las cartas que nadie ha visto. El mazo siempre tiene 10 cartas y las expulsadas van al fondo, así que después de 10 inserciones se conoce el orden de tipos del mazo. `belief = "mazo"` promedia sobre todo el mazo, como la codiciosa.
  - Con `horizon = "accion"` y `belief = "mazo"` da los mismos valores que `policy_greedy()`; solo cambia el desempate al azar.

Una política es una lista con `choose_action(state, allow_insert)` y `choose_insertion(state, card)`. Para probar una estrategia nueva basta con escribir esas dos funciones. `state$action` dice si es la primera o la segunda acción del turno, y `state$n_unknown` cuántas cartas de arriba del mazo no se han visto.

## Motor rápido

`src/engine.cpp` replica en C++ la conectividad, las capturas con reacciones en cadena, la inserción y la búsqueda. `fast.R` lo compila con Rcpp (`load_fast_engine()`) y traduce el estado de R. El motor en R sigue siendo la referencia: `tests/test_fast.R` compara los dos en 2000 estados al azar con capturas y cadenas, en 300 rotaciones e inserciones y en los valores de todas las jugadas de la codiciosa.

- El tablero se pasa como 16 códigos de carta: código 0–23 = (tipo − 1) × 4 + r, con la matriz de `rotate_card(card_types[[tipo]], 90 * r)`.
- La conectividad usa union-find sobre los 64 nodos (celda, lado).
- La primera compilación tarda unos 30 s; queda en `.rcpp_cache/` y los procesos siguientes la reutilizan.
- En Windows hace falta Rtools. Si `make` no está en el PATH, `load_fast_engine()` usa el de `RTOOLS45_HOME` o `RTOOLS44_HOME` (R 4.6 compila con el Rtools 4.4 de `C:/RBuildTools/4.4`).

## Resultados con los parámetros por defecto

30 partidas en solitario (`main.R`, semilla 2016), medidas cuando un giro de 180° costaba una sola acción; las cifras con la regla actual están en el banco de pruebas de abajo:

| Política | Partidas completas | Turnos: media | Mediana | Mejor | Peor |
|---|---|---|---|---|---|
| Aleatoria (tope de 400 turnos) | 90 % | 266 | 270 | 116 | 400 |
| Codiciosa | 100 % | 10,6 | 11 | 5 | 15 |
| Codiciosa sin rotación obligatoria | 100 % | 8,6 | 9 | 6 | 12 |

Permitir dos inserciones por turno ahorra unos dos turnos por partida. Estos números dependen de la composición supuesta del mazo.

### Banco de pruebas

`bench.R` juega cada política en los mismos repartos y compara por semilla. Hay dos conjuntos fijos: desarrollo (100 repartos, semillas 1001–1100) para ajustar pesos y evaluación (200 repartos, semillas 5001–5200) para reportar. Corre en paralelo y guarda un CSV en `results/` con el hash del commit.

```bash
Rscript bench.R eval oficial
```

Conjunto de evaluación, 200 repartos, reglas oficiales, rotaciones de 90° o 270° por acción:

| Política | Turnos: media | Mediana | Mejor | Peor | ≤ 7 turnos | Diferencia pareada con la codiciosa (IC 95 %) |
|---|---|---|---|---|---|---|
| Codiciosa (`policy_greedy`) | 11,3 | 11 | 6 | 18 | 3,5 % | — |
| Búsqueda por acción, mazo completo | 11,1 | 11 | 6 | 18 | 4 % | −0,15 (−0,46 a 0,17) |
| Búsqueda por acción, creencia exacta | 11,2 | 11 | 6 | 19 | 4 % | −0,12 (−0,42 a 0,23) |
| **Búsqueda por turno, creencia exacta** | **8,9** | **9** | 5 | 13 | **20 %** | **−2,38 (−2,69 a −2,08)** |

Sin rotación obligatoria (`Rscript bench.R eval libre`) la búsqueda por turno baja de 8,7 a 6,8 turnos de media (−1,86; IC −2,12 a −1,62) y termina en 7 turnos o menos en el 78 % de los repartos.

- La búsqueda sobre el turno completo es la mejora grande: gana en el 84 % de los repartos y pierde en el 8,5 %.
- La creencia exacta del mazo no cambia el resultado a una acción: la diferencia con el promedio sobre todo el mazo está dentro del ruido.
- Con 15 partidas en paralelo, la codiciosa tarda de 13 a 25 s por partida y la búsqueda por turno alrededor de 1 s.
- Antes del cambio de regla, cuando un giro de 180° costaba una acción, las cifras eran 11,1 y 8,5 turnos (CSV con sufijo `_giro180` en `results/`).

## Uso

```r
source("functions.R"); source("cards.R"); source("game.R"); source("policies.R")
res <- play_game(game_params(), policy_greedy(), n_players = 1, seed = 1)
res$turns

source("fast.R"); load_fast_engine()
res <- play_game(game_params(), policy_search(horizon = "turno"), n_players = 1, seed = 1)
```

`main.R` ejecuta una partida de ejemplo y las comparaciones de arriba; tarda unos 4 minutos.

## Pruebas

```bash
Rscript -e 'testthat::test_file("tests/test_game.R")'
```

```bash
Rscript -e 'testthat::test_file("tests/test_fast.R")'
```

`test_game.R` prueba las reglas (20 tests). `test_fast.R` prueba la equivalencia del motor rápido con el de R (6 tests, unas 14 000 expectativas, menos de 1 minuto).

## Archivos

| Archivo | Contenido |
|---|---|
| `cards.R` | Los seis tipos de carta de canaleta. |
| `functions.R` | Utilidades: matrices, álgebra circular, inserción y rotación sobre el tablero. |
| `game.R` | Motor de reglas: preparación, caminos, acciones, capturas, turnos y partida. |
| `policies.R` | Políticas de juego y simulación en solitario. |
| `fast.R` | Carga del motor rápido y traducción del estado de R. |
| `src/engine.cpp` | Motor rápido en C++: capturas, inserción y búsqueda. |
| `bench.R` | Banco de pruebas pareado de políticas en solitario. |
| `main.R` | Ejemplo de uso y comparación de políticas. |
| `tests/test_game.R` | Pruebas del motor. |
| `tests/test_fast.R` | Equivalencia del motor rápido con el motor en R. |
