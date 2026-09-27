options(stringsAsFactors = FALSE)

fair_fssgl_fail <- function(...) stop(sprintf(...), call. = FALSE)
table_dir <- "results/tables/simulation/v3_fair"
data_dir <- "data/processed/simulation/v3_fair"
result_file <- file.path(table_dir, "fssgl_v3_fair_replicates.csv")
path_file <- file.path(table_dir, "fssgl_v3_fair_cv_path.csv")
summary_file <- file.path(table_dir, "fssgl_v3_fair_summary.csv")
paired_file <- file.path(table_dir, "fssgl_v3_fair_paired_differences.csv")
manifest_file <- file.path(data_dir, "fssgl_v3_fair_manifest.rds")
baseline_file <- file.path(table_dir, "baseline_v3_fair_replicates.csv")
for (path in c(result_file, path_file, summary_file, paired_file, manifest_file, baseline_file)) {
  if (!file.exists(path)) fair_fssgl_fail("Missing fair-comparison artifact: %s", path)
}

results <- read.csv(result_file, check.names = FALSE)
path <- read.csv(path_file, check.names = FALSE)
summary_table <- read.csv(summary_file, check.names = FALSE)
paired_table <- read.csv(paired_file, check.names = FALSE)
baseline_replicates <- read.csv(baseline_file, check.names = FALSE)
baseline <- unique(baseline_replicates[
  c("p", "rep", "seed", "fold_seed")
])
expected_repetitions <- c(`10` = 100L, `20` = 100L, `50` = 50L, `60` = 50L, `100` = 50L)
if (nrow(results) != sum(expected_repetitions) ||
    any(duplicated(results[c("p", "rep")]))) {
  fair_fssgl_fail("Fair FSSGL must have 350 unique p/rep results")
}
for (p_value in as.integer(names(expected_repetitions))) {
  reps <- sort(results$rep[results$p == p_value])
  if (!identical(reps, seq_len(expected_repetitions[[as.character(p_value)]]))) {
    fair_fssgl_fail("Unexpected FSSGL repetitions at p=%d", p_value)
  }
}
paired <- merge(
  results[c("p", "rep", "seed", "fold_seed")], baseline,
  by = c("p", "rep"), suffixes = c("_fssgl", "_baseline")
)
if (nrow(paired) != nrow(results) ||
    any(paired$seed_fssgl != paired$seed_baseline) ||
    any(paired$fold_seed_fssgl != paired$fold_seed_baseline)) {
  fair_fssgl_fail("FSSGL and baselines are not paired on data and CV folds")
}
allowed_status <- c("ok", "nonconverged", "grid_boundary", "tuning_error", "fit_error")
if (anyNA(results$status) || any(!results$status %in% allowed_status)) {
  fair_fssgl_fail("Unknown or missing FSSGL result status")
}
ok <- results[results$status == "ok", , drop = FALSE]
if (nrow(ok) && (any(!ok$strict_converged) ||
    any(ok$selected_at_lower_boundary |
        (ok$selected_at_upper_boundary & !ok$upper_null_plateau)) ||
    any(ok$final_path_fits < 1L) ||
    any(!is.finite(as.matrix(ok[c(
      "covariate_tpr", "covariate_fdr", "coefficient_relative_error",
      "test_coeff_rmse", "tuning_runtime_sec", "final_runtime_sec",
      "total_runtime_sec"
    )]))))) {
  fair_fssgl_fail("Reportable FSSGL fits have invalid convergence, grid, or metrics")
}
if (any(duplicated(path[c("p", "rep", "fold_id", "multiplier")]))) {
  fair_fssgl_fail("Duplicate FSSGL CV-path row")
}
if (nrow(path) && (any(!path$fold_id %in% seq_len(5L)) ||
    any(!path$p %in% as.integer(names(expected_repetitions))))) {
  fair_fssgl_fail("Invalid FSSGL CV-path dimension or fold")
}

metric_columns <- c(
  "covariate_tpr", "covariate_fdr", "coefficient_relative_error",
  "test_coeff_rmse", "tuning_runtime_sec", "final_runtime_sec", "total_runtime_sec"
)
recomputed <- if (nrow(ok)) do.call(rbind, lapply(split(ok, ok$p), function(block) {
  do.call(rbind, lapply(metric_columns, function(metric) {
    values <- block[[metric]]
    data.frame(
      p = block$p[[1L]], metric = metric, n = length(values),
      mean = mean(values), sd = stats::sd(values),
      mcse = stats::sd(values) / sqrt(length(values))
    )
  }))
})) else data.frame()
key <- function(data) paste(data$p, data$metric, sep = "|")
recomputed <- recomputed[order(key(recomputed)), , drop = FALSE]
summary_table <- summary_table[order(key(summary_table)), , drop = FALSE]
if (!identical(key(recomputed), key(summary_table)) ||
    any(recomputed$n != summary_table$n) ||
    !isTRUE(all.equal(
      recomputed[c("mean", "sd", "mcse")],
      summary_table[c("mean", "sd", "mcse")],
      tolerance = 1e-10, check.attributes = FALSE
    ))) {
  fair_fssgl_fail("Fair FSSGL summary does not reproduce reportable replicates")
}

direct_methods <- c("basis_group_scad", "fpca_group_scad", "structured_group_lasso")
direct <- baseline_replicates[
  baseline_replicates$method_id %in% direct_methods, , drop = FALSE
]
matched <- merge(direct, results, by = c("p", "rep"), suffixes = c("_baseline", "_fssgl"))
if (nrow(matched) != 3L * nrow(results) || nrow(paired_table) != 60L) {
  fair_fssgl_fail("Direct paired-comparison evidence is incomplete")
}
paired_key <- function(data) paste(data$p, data$method_id, data$metric, sep = "|")
for (row_id in seq_len(nrow(paired_table))) {
  row <- paired_table[row_id, , drop = FALSE]
  block <- matched[matched$p == row$p & matched$method_id == row$method_id, , drop = FALSE]
  difference <- block[[paste0(row$metric, "_baseline")]] -
    block[[paste0(row$metric, "_fssgl")]]
  if (length(difference) != row$n ||
      abs(mean(difference) - row$mean_baseline_minus_fssgl) > 1e-10 ||
      abs(stats::sd(difference) - row$sd) > 1e-10 ||
      abs(row$mcse - row$sd / sqrt(row$n)) > 1e-10 ||
      abs(row$ci_low - (row$mean_baseline_minus_fssgl - 1.96 * row$mcse)) > 1e-10 ||
      abs(row$ci_high - (row$mean_baseline_minus_fssgl + 1.96 * row$mcse)) > 1e-10) {
    fair_fssgl_fail("Paired comparison does not reproduce row %d", row_id)
  }
}
if (any(duplicated(paired_key(paired_table)))) {
  fair_fssgl_fail("Duplicate paired-comparison summary key")
}

manifest <- readRDS(manifest_file)
if (!identical(manifest$experiment_id, "fssgl_v3_fair_common_fold_comparison") ||
    !identical(as.integer(manifest$design$p), as.integer(names(expected_repetitions))) ||
    !identical(as.integer(manifest$design$repetitions), unname(expected_repetitions)) ||
    !identical(manifest$design$n_folds, 5L) ||
    !identical(manifest$design$test_set_used_for_tuning, FALSE)) {
  fair_fssgl_fail("Fair FSSGL manifest does not match the registered design")
}
cat(sprintf(
  "Fair FSSGL validation passed: %d/350 reportable fits; statuses: %s.\n",
  nrow(ok), paste(names(table(results$status)), table(results$status), collapse = ", ")
))
