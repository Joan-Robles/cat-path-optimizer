// Motor rápido de "No sin mi gato" para evaluar jugadas.
// Replica las reglas de game.R; el motor en R queda como referencia.
//
// Convenciones:
//  - Código de carta 0..23 = tipo * 4 + r, con la matriz de
//    rotate_card(card_types[[tipo]], 90 * r). Lo carga fe_init() desde R.
//  - Celda k = (i - 1) * 4 + (j - 1); nodo = k * 4 + (lado - 1).
//  - Anillo: 0 vacío, 1 gato, 2 hogar. Posiciones 0-based (p = pos - 1).
//  - Creencia del mazo: unk[t] cartas desconocidas de cada tipo arriba del
//    mazo y, debajo, la secuencia conocida de tipos (cartas expulsadas).

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
using namespace Rcpp;

static const int N_CODES = 24;
static const int N_TYPES = 6;
static int CONN[N_CODES][4][4];
static bool OPEN[N_CODES][4];
static int CANON[N_CODES];
static bool INITED = false;

struct St {
  int codes[16];
  int ring[16];
  double pts[16];
  int cat_next;            // siguiente gato del mazo de gatos
  int to_win;              // gatos que faltan para terminar
  int unk[N_TYPES];        // cartas desconocidas por tipo
  std::vector<int> known;  // tipos conocidos, en orden, debajo de las desconocidas
};

struct Ctx {
  std::vector<double> cat_deck;
  bool mandatory;
  bool use_points;
  double w_reach;
};

static inline int rot_code(int c, int q) { return (c / 4) * 4 + ((c % 4) + q) % 4; }

// Celda y lado (0-based) contiguos a la posición p del anillo
static inline void ring_cell_side(int p, int &cell, int &side) {
  int grp = p / 4, k = p % 4 + 1;
  int i, j;
  switch (grp) {
    case 0: i = 5 - k; j = 1; side = 3; break;  // izquierda, lado O
    case 1: i = 1; j = k; side = 0; break;      // arriba, lado N
    case 2: i = k; j = 4; side = 1; break;      // derecha, lado E
    default: i = 4; j = 5 - k; side = 2; break; // abajo, lado S
  }
  cell = (i - 1) * 4 + (j - 1);
}

static inline int uf_find(int *par, int x) {
  while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; }
  return x;
}
static inline void uf_unite(int *par, int a, int b) {
  a = uf_find(par, a); b = uf_find(par, b);
  if (a != b) par[a] = b;
}

// Etiqueta de componente y tamaño por nodo
static void components(const int *codes, int *comp, int *size) {
  int par[64];
  for (int x = 0; x < 64; x++) par[x] = x;
  for (int k = 0; k < 16; k++) {
    int c = codes[k];
    for (int a = 0; a < 4; a++) for (int b = a + 1; b < 4; b++)
      if (CONN[c][a][b]) uf_unite(par, k * 4 + a, k * 4 + b);
  }
  for (int i = 0; i < 4; i++) for (int j = 0; j < 3; j++) {
    int k = i * 4 + j, k2 = k + 1;
    if (OPEN[codes[k]][1] && OPEN[codes[k2]][3]) uf_unite(par, k * 4 + 1, k2 * 4 + 3);
  }
  for (int i = 0; i < 3; i++) for (int j = 0; j < 4; j++) {
    int k = i * 4 + j, k2 = k + 4;
    if (OPEN[codes[k]][2] && OPEN[codes[k2]][0]) uf_unite(par, k * 4 + 2, k2 * 4 + 0);
  }
  for (int x = 0; x < 64; x++) size[x] = 0;
  for (int x = 0; x < 64; x++) { comp[x] = uf_find(par, x); size[comp[x]]++; }
}

// Nodo de entrada desde la posición p (-1 si el lado está cerrado)
static inline int entry_node(const int *codes, int p) {
  int cell, side;
  ring_cell_side(p, cell, side);
  if (!OPEN[codes[cell]][side]) return -1;
  return cell * 4 + side;
}

static int next_free(const int *ring, int from) {
  int p = from;
  for (int step = 0; step < 16; step++) {
    p = (p + 1) % 16;
    if (ring[p] == 0) return p;
  }
  Rcpp::stop("No hay posiciones libres en el anillo.");
}

struct Capture { int cat, home; double points; };

// Gatos ganados y reacciones en cadena (como resolve_captures() en game.R).
// Deja en comp/size las componentes del tablero, que no cambia.
static void resolve(St &s, const Ctx &ctx, int &cats, double &points,
                    int *comp, int *size, std::vector<Capture> *hist = nullptr) {
  cats = 0; points = 0;
  components(s.codes, comp, size);
  int pc[16], ph[16];
  for (;;) {
    int n = 0;
    for (int c = 0; c < 16; c++) {
      if (s.ring[c] != 1) continue;
      int ec = entry_node(s.codes, c);
      if (ec < 0) continue;
      for (int h = 0; h < 16; h++) {
        if (s.ring[h] != 2) continue;
        int eh = entry_node(s.codes, h);
        if (eh >= 0 && comp[eh] == comp[ec]) { pc[n] = c; ph[n] = h; n++; break; }
      }
    }
    if (n == 0) return;
    for (int r = 0; r < n; r++) {
      cats++; points += s.pts[pc[r]];
      if (hist) hist->push_back({pc[r], ph[r], s.pts[pc[r]]});
      s.ring[pc[r]] = 0; s.pts[pc[r]] = NA_REAL;
    }
    s.to_win -= n;
    for (int r = 0; r < n; r++) {
      int h = ph[r];
      s.ring[h] = 0;
      int nh = next_free(s.ring, h);
      s.ring[nh] = 2;
      for (int r2 = 0; r2 < n; r2++) if (ph[r2] == h) ph[r2] = nh;
      if (s.cat_next < (int) ctx.cat_deck.size()) {
        int nc = next_free(s.ring, pc[r]);
        s.ring[nc] = 1;
        s.pts[nc] = ctx.cat_deck[s.cat_next++];
      }
    }
  }
}

static double cat_reach(const St &s, const int *comp, const int *size) {
  double total = 0;
  for (int c = 0; c < 16; c++) {
    if (s.ring[c] != 1) continue;
    int e = entry_node(s.codes, c);
    if (e >= 0) total += size[comp[e]];
  }
  return total;
}

// Inserta el código c por la posición p; devuelve el código expulsado
static int insert_code(int *codes, int c, int p) {
  int grp = p / 4, out;
  if (grp == 0) {        // fila 3 - p % 4, entra por la izquierda
    int i = 3 - p % 4;
    out = codes[i * 4 + 3];
    for (int j = 3; j > 0; j--) codes[i * 4 + j] = codes[i * 4 + j - 1];
    codes[i * 4] = c;
  } else if (grp == 1) { // columna p % 4, entra por arriba
    int j = p % 4;
    out = codes[12 + j];
    for (int i = 3; i > 0; i--) codes[i * 4 + j] = codes[(i - 1) * 4 + j];
    codes[j] = c;
  } else if (grp == 2) { // fila p % 4, entra por la derecha
    int i = p % 4;
    out = codes[i * 4];
    for (int j = 0; j < 3; j++) codes[i * 4 + j] = codes[i * 4 + j + 1];
    codes[i * 4 + 3] = c;
  } else {               // columna 15 - p, entra por abajo
    int j = 15 - p;
    out = codes[j];
    for (int i = 0; i < 3; i++) codes[i * 4 + j] = codes[(i + 1) * 4 + j];
    codes[12 + j] = c;
  }
  return out;
}

// Giros que dan orientaciones distintas de la carta c.
// Para una rotación (rotation = true) solo cuentan los cuartos de vuelta
// q = 1 (90) y q = 3 (270) que cambian la carta: girar 180 grados cuesta dos
// acciones. Para una inserción, q = 0..3.
static int distinct_turns(int c, bool rotation, int *qs) {
  static const int ROT_Q[2] = {1, 3}, INS_Q[4] = {0, 1, 2, 3};
  const int *cand = rotation ? ROT_Q : INS_Q;
  int n_cand = rotation ? 2 : 4;
  int seen[5], n_seen = 0, n = 0;
  if (rotation) seen[n_seen++] = CANON[c];
  for (int t = 0; t < n_cand; t++) {
    int q = cand[t];
    int k = CANON[rot_code(c, q)];
    bool dup = false;
    for (int t = 0; t < n_seen; t++) if (seen[t] == k) dup = true;
    if (dup) continue;
    seen[n_seen++] = k;
    qs[n++] = q;
  }
  return n;
}

struct Draw { int type; double prob; };

static int draws(const St &s, Draw *out) {
  int total = 0, n = 0;
  for (int t = 0; t < N_TYPES; t++) total += s.unk[t];
  if (total > 0) {
    for (int t = 0; t < N_TYPES; t++)
      if (s.unk[t] > 0) out[n++] = {t, (double) s.unk[t] / total};
  } else if (!s.known.empty()) {
    out[n++] = {s.known[0], 1.0};
  }
  return n;
}

static void take_draw(St &s, int type) {
  int total = 0;
  for (int t = 0; t < N_TYPES; t++) total += s.unk[t];
  if (total > 0) s.unk[type]--; else s.known.erase(s.known.begin());
}

static double value(const St &s, int left, bool can_insert, const Ctx &ctx);

static inline double gain_of(int cats, double points, const Ctx &ctx) {
  return ctx.use_points ? points : (double) cats;
}

// Valor de rotar la celda k q cuartos de vuelta y seguir
static double rotation_value(const St &s, int k, int q, int left, bool can_insert, const Ctx &ctx) {
  St s2 = s;
  s2.codes[k] = rot_code(s.codes[k], q);
  int cats, comp[64], size[64]; double points;
  resolve(s2, ctx, cats, points, comp, size);
  double g = gain_of(cats, points, ctx);
  if (left <= 1 || s2.to_win <= 0) return g + (s2.to_win > 0 ? ctx.w_reach * cat_reach(s2, comp, size) : 0);
  return g + value(s2, left - 1, can_insert, ctx);
}

// Valor de insertar el código c (ya robado) por p y seguir
static double insertion_value(const St &s, int c, int p, int left, bool can_insert_after, const Ctx &ctx) {
  St s2 = s;
  int out = insert_code(s2.codes, c, p);
  s2.known.push_back(out / 4);
  int cats, comp[64], size[64]; double points;
  resolve(s2, ctx, cats, points, comp, size);
  double g = gain_of(cats, points, ctx);
  if (left <= 1 || s2.to_win <= 0) return g + (s2.to_win > 0 ? ctx.w_reach * cat_reach(s2, comp, size) : 0);
  return g + value(s2, left - 1, can_insert_after, ctx);
}

static double best_insertion_value(const St &s_after_draw, int c, int left, bool can_after, const Ctx &ctx) {
  int qs[4]; int nq = distinct_turns(c, false, qs);
  double best = -std::numeric_limits<double>::infinity();
  for (int a = 0; a < nq; a++) for (int p = 0; p < 16; p++) {
    double v = insertion_value(s_after_draw, rot_code(c, qs[a]), p, left, can_after, ctx);
    if (v > best) best = v;
  }
  return best;
}

// Valor esperado de insertar sin conocer la carta: E sobre la carta robada
static double insert_ev(const St &s, int left, const Ctx &ctx, NumericVector *by_type = nullptr) {
  Draw d[N_TYPES]; int nd = draws(s, d);
  if (nd == 0) return -std::numeric_limits<double>::infinity();
  double ev = 0;
  bool can_after = !ctx.mandatory;
  for (int r = 0; r < nd; r++) {
    St s2 = s;
    take_draw(s2, d[r].type);
    double b = best_insertion_value(s2, d[r].type * 4, left, can_after, ctx);
    if (by_type) (*by_type)[d[r].type] = b;
    ev += d[r].prob * b;
  }
  return ev;
}

// Mejor valor esperado con `left` acciones por jugar
static double value(const St &s, int left, bool can_insert, const Ctx &ctx) {
  double best = -std::numeric_limits<double>::infinity();
  for (int k = 0; k < 16; k++) {
    int qs[2]; int nq = distinct_turns(s.codes[k], true, qs);
    for (int a = 0; a < nq; a++) {
      double v = rotation_value(s, k, qs[a], left, can_insert, ctx);
      if (v > best) best = v;
    }
  }
  if (can_insert) {
    double ev = insert_ev(s, left, ctx);
    if (ev > best) best = ev;
  }
  return best;
}

// Interfaz con R ---------------------------------------------------------

static St make_state(IntegerVector codes, IntegerVector ring, NumericVector pts,
                     int to_win, IntegerVector unk, IntegerVector known) {
  if (!INITED) stop("Llamar fe_init() antes de usar el motor.");
  St s;
  for (int k = 0; k < 16; k++) { s.codes[k] = codes[k]; s.ring[k] = ring[k]; s.pts[k] = pts[k]; }
  s.cat_next = 0; s.to_win = to_win;
  for (int t = 0; t < N_TYPES; t++) s.unk[t] = unk[t];
  s.known.assign(known.begin(), known.end());
  return s;
}

static Ctx make_ctx(NumericVector cat_deck, bool mandatory, bool use_points, double w_reach) {
  Ctx ctx;
  ctx.cat_deck.assign(cat_deck.begin(), cat_deck.end());
  ctx.mandatory = mandatory; ctx.use_points = use_points; ctx.w_reach = w_reach;
  return ctx;
}

// [[Rcpp::export]]
void fe_init(IntegerVector conn, IntegerVector canon) {
  // conn: arreglo 24 x 4 x 4 de R (orden por columnas)
  for (int c = 0; c < N_CODES; c++) {
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++)
      CONN[c][a][b] = conn[c + N_CODES * a + N_CODES * 4 * b];
    for (int a = 0; a < 4; a++) {
      int deg = 0;
      for (int b = 0; b < 4; b++) if (b != a) deg += CONN[c][a][b];
      OPEN[c][a] = deg > 0;
    }
    CANON[c] = canon[c];
  }
  INITED = true;
}

// Aplica resolve() y devuelve el estado resultante (para pruebas y el motor)
// [[Rcpp::export]]
List fe_resolve(IntegerVector codes, IntegerVector ring, NumericVector pts,
                NumericVector cat_deck, int to_win) {
  IntegerVector none(N_TYPES), known(0);
  St s = make_state(codes, ring, pts, to_win, none, known);
  Ctx ctx = make_ctx(cat_deck, true, false, 0);
  int cats, comp[64], size[64]; double points;
  std::vector<Capture> hist;
  resolve(s, ctx, cats, points, comp, size, &hist);
  IntegerVector ring_out(16); NumericVector pts_out(16);
  for (int k = 0; k < 16; k++) { ring_out[k] = s.ring[k]; pts_out[k] = s.pts[k]; }
  IntegerVector hc(hist.size()), hh(hist.size()); NumericVector hp(hist.size());
  for (size_t r = 0; r < hist.size(); r++) { hc[r] = hist[r].cat + 1; hh[r] = hist[r].home + 1; hp[r] = hist[r].points; }
  return List::create(_["ring"] = ring_out, _["pts"] = pts_out, _["cats"] = cats,
                      _["points"] = points, _["cats_used"] = s.cat_next,
                      _["reach"] = cat_reach(s, comp, size),
                      _["history"] = DataFrame::create(_["cat"] = hc, _["home"] = hh, _["points"] = hp));
}

// [[Rcpp::export]]
List fe_insert(IntegerVector codes, int code, int pos) {
  IntegerVector out = clone(codes);
  int expelled = insert_code(out.begin(), code, pos - 1);
  return List::create(_["codes"] = out, _["expelled"] = expelled);
}

// Valores de todas las jugadas de la primera acción.
// rot: matriz (celda k 1-based, q, valor); ins_ev: valor esperado de
// insertar (-Inf si no se puede); ins_by_type: mejor inserción por tipo.
// [[Rcpp::export]]
List fe_children(IntegerVector codes, IntegerVector ring, NumericVector pts,
                 NumericVector cat_deck, int to_win, IntegerVector unk, IntegerVector known,
                 int left, bool can_insert, bool mandatory, bool use_points, double w_reach) {
  St s = make_state(codes, ring, pts, to_win, unk, known);
  Ctx ctx = make_ctx(cat_deck, mandatory, use_points, w_reach);
  std::vector<double> rk, rq, rv;
  for (int k = 0; k < 16; k++) {
    int qs[2]; int nq = distinct_turns(s.codes[k], true, qs);
    for (int a = 0; a < nq; a++) {
      rk.push_back(k + 1); rq.push_back(qs[a]);
      rv.push_back(rotation_value(s, k, qs[a], left, can_insert, ctx));
    }
  }
  NumericMatrix rot(rk.size(), 3);
  for (size_t r = 0; r < rk.size(); r++) { rot(r, 0) = rk[r]; rot(r, 1) = rq[r]; rot(r, 2) = rv[r]; }
  NumericVector by_type(N_TYPES, NA_REAL);
  double ev = can_insert ? insert_ev(s, left, ctx, &by_type) : -std::numeric_limits<double>::infinity();
  return List::create(_["rot"] = rot, _["ins_ev"] = ev, _["ins_by_type"] = by_type);
}

// Valores de insertar la carta robada `code` (estado después de robar).
// Matriz (pos 1-based, q, valor), con pos variando más rápido.
// [[Rcpp::export]]
NumericMatrix fe_insertions(IntegerVector codes, IntegerVector ring, NumericVector pts,
                            NumericVector cat_deck, int to_win, IntegerVector unk, IntegerVector known,
                            int code, int left, bool mandatory, bool use_points, double w_reach) {
  St s = make_state(codes, ring, pts, to_win, unk, known);
  Ctx ctx = make_ctx(cat_deck, mandatory, use_points, w_reach);
  int qs[4]; int nq = distinct_turns(code, false, qs);
  NumericMatrix out(nq * 16, 3);
  int r = 0;
  for (int a = 0; a < nq; a++) for (int p = 0; p < 16; p++, r++) {
    out(r, 0) = p + 1; out(r, 1) = qs[a];
    out(r, 2) = insertion_value(s, rot_code(code, qs[a]), p, left, !mandatory, ctx);
  }
  return out;
}
