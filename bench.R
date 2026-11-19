##########################################
# Banco de pruebas de políticas en solitario
#
# Uso desde la raíz del proyecto:
#   Rscript bench.R [eval|dev] [oficial|libre] [núcleos]
# - eval: 200 repartos (semillas 5001:5200), para reportar.
# - dev: 100 repartos (semillas 1001:1100), para ajustar pesos.
# - oficial: rotation_mandatory = TRUE; libre: FALSE.
#
# Cada semilla fija el reparto, el orden del mazo y el de los gatos; después
# de new_game() el motor es determinista, así que todas las políticas juegan
# los mismos repartos y la comparación es pareada por semilla.
##########################################

SEEDS_DEV <- 1001:1100
SEEDS_EVAL <- 5001:5200

# Políticas como llamadas, para construirlas dentro de cada proceso
BENCH_POLICIES <- list(
  codiciosa = quote(policy_greedy()),
  accion_mazo = quote(policy_search(horizon = "accion", belief = "mazo")),
  accion_exacta = quote(policy_search(horizon = "accion", belief = "exacta")),
  turno_exacta = quote(policy_search(horizon = "turno", belief = "exacta"))
)

bench_setup <- function(root) {
  setwd(root)
  for (f in c("functions.R", "cards.R", "game.R", "policies.R", "fast.R")) source(f)
  load_fast_engine(root)
}

commit_id <- function(root) {
  h <- tryCatch(system2("git", c("-C", shQuote(root), "rev-parse", "--short", "HEAD"),
                        stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
  dirty <- tryCatch(length(system2("git", c("-C", shQuote(root), "status", "--porcelain"),
                                   stdout = TRUE, stderr = FALSE)) > 0, error = function(e) FALSE)
  paste0(h, if (dirty) "-dirty" else "")
}

# Juega todas las semillas con todas las políticas, en paralelo
run_bench <- function(params, policies = BENCH_POLICIES, seeds = SEEDS_EVAL,
                      cores = parallel::detectCores() - 1L, max_turns = 200L,
                      root = getwd()) {
  root <- normalizePath(root, winslash = "/")
  # Compilar en el proceso principal deja el motor en la caché para los demás
  load_fast_engine(root)
  jobs <- expand.grid(seed = seeds, policy = names(policies), stringsAsFactors = FALSE)
  run_one <- function(r) {
    pol <- eval(policies[[jobs$policy[r]]])
    t <- system.time(res <- play_game(params, pol, n_players = 1L,
                                      seed = jobs$seed[r], max_turns = max_turns))
    data.frame(policy = jobs$policy[r], seed = jobs$seed[r], turns = res$turns,
               completed = res$completed, secs = t[["elapsed"]])
  }
  environment(run_one) <- globalenv()
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl))
  # Uno por uno: sourceCpp reescribe archivos de la caché al cargarla y, en
  # paralelo, un proceso puede leerlos vacíos
  for (k in seq_along(cl)) parallel::clusterCall(cl[k], bench_setup, root)
  parallel::clusterExport(cl, c("jobs", "policies", "params", "max_turns"), envir = environment())
  rows <- parallel::parLapplyLB(cl, seq_len(nrow(jobs)), run_one)
  out <- do.call(rbind, rows)
  out$policy <- factor(out$policy, levels = names(policies))
  out
}

bench_summary <- function(res) {
  do.call(rbind, lapply(split(res, res$policy), function(d) {
    data.frame(policy = d$policy[1], partidas = nrow(d), censuradas = sum(!d$completed),
               media = mean(d$turns), mediana = median(d$turns), de = sd(d$turns),
               mejor = min(d$turns), peor = max(d$turns),
               pct_7 = 100 * mean(d$turns <= 7), s_partida = mean(d$secs))
  }))
}

# Diferencia pareada por semilla contra la política base (negativa = mejor)
bench_paired <- function(res, base = levels(res$policy)[1], n_boot = 5000L, seed = 1L) {
  wide <- reshape(res[, c("policy", "seed", "turns")], idvar = "seed",
                  timevar = "policy", direction = "wide")
  set.seed(seed)
  others <- setdiff(levels(res$policy), base)
  do.call(rbind, lapply(others, function(p) {
    d <- wide[[paste0("turns.", p)]] - wide[[paste0("turns.", base)]]
    boot <- replicate(n_boot, mean(sample(d, replace = TRUE)))
    pval <- if (all(d == 0)) 1 else suppressWarnings(wilcox.test(d, exact = FALSE)$p.value)
    data.frame(policy = p, base = base, dif_media = mean(d),
               ic95_inf = unname(quantile(boot, 0.025)), ic95_sup = unname(quantile(boot, 0.975)),
               mejora = mean(d < 0), empate = mean(d == 0), empeora = mean(d > 0),
               wilcoxon_p = pval)
  }))
}

bench_tiers <- function(res) table(res$policy, solo_tier(res$turns))

save_bench <- function(res, label, root = getwd()) {
  dir.create(file.path(root, "results"), showWarnings = FALSE)
  id <- commit_id(root)
  res$commit <- id
  path <- file.path(root, "results", sprintf("bench_%s_%s_%s.csv", label, format(Sys.Date()), id))
  write.csv(res, path, row.names = FALSE)
  path
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  set <- if (length(args) >= 1) args[1] else "eval"
  rules <- if (length(args) >= 2) args[2] else "oficial"
  cores <- if (length(args) >= 3) as.integer(args[3]) else parallel::detectCores() - 1L
  bench_setup(getwd())
  params <- game_params(rotation_mandatory = rules == "oficial")
  seeds <- if (set == "dev") SEEDS_DEV else SEEDS_EVAL
  t <- system.time(res <- run_bench(params, seeds = seeds, cores = cores))
  cat(sprintf("Conjunto %s, reglas %s: %d partidas en %.0f s\n\n", set, rules, nrow(res), t[["elapsed"]]))
  print(bench_summary(res), row.names = FALSE, digits = 3)
  cat("\nDiferencia pareada en turnos contra la codiciosa (negativa = mejor):\n")
  print(bench_paired(res), row.names = FALSE, digits = 3)
  cat("\nEscala del modo solitario:\n")
  print(bench_tiers(res))
  cat("\nResultados en", save_bench(res, paste(set, rules, sep = "_")), "\n")
}
