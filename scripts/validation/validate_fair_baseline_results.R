options(stringsAsFactors = FALSE)

fair_fail <- function(...) stop(sprintf(...), call. = FALSE)
table_dir <- "results/tables/simulation/v3_fair"
data_dir <- "data/processed/simulation/v3_fair"
replicate_file <- file.path(table_dir, "baseline_v3_fair_replicates.csv")
summary_file <- file.path(table_dir, "baseline_v3_fair_summary.csv")
manifest_file <- file.path(data_dir, "baseline_v3_fair_manifest.rds")

for (path in c(replicate_file, summary_file, manifest_file)) {
  if (!file.exists(path)) fair_fail("Missing fair-baseline artifact: %s", path)
}

replicates <- read.csv(replicate_file, check.names = FALSE)
summary_table <- read.csv(summary_file, check.names = FALSE)
required_columns <- c(
  "tuning_protocol_version", "method_id", "p", "rep", "seed", "fold_seed",
  "validation_loss", "status", "covariate_tpr", "covariate_fdr",
  "coefficient_relative_error", "test_coeff_rmse", "tuning_runtime_sec",
  "final_runtime_sec", "runtime_sec"
)
missing_columns <- setdiff(required_columns, names(replicates))
if (length(missing_columns)) {
  fair_fail(
    "Fair-baseline replicates are missing columns: %s",
    paste(missing_columns, collapse = ", ")
  )
}
if (nrow(replicates) != 2050L) {
  fair_fail("Fair-baseline experiment must contain 2050 rows, found %d", nrow(replicates))
}
if (any(duplicated(replicates[c("p", "rep", "method_id")]))) {
  fair_fail("Fair-baseline experiment contains duplicate p/rep/method rows")
}
if (any(is.na(replicates$status)) || any(replicates$status != "ok")) {
  fair_fail("Fair-baseline experiment contains non-ok rows")
}

expected_repetitions <- c(`10` = 100L, `20` = 100L, `50` = 50L, `60` = 50L, `100` = 50L)
low_methods <- c(
  "basis_ridge", "basis_group_scad", "fpca_ridge", "fpca_group_lasso",
  "fpca_group_scad", "fpca_aenet", "structured_group_lasso", "kernel_ridge"
)
high_methods <- c("basis_group_scad", "fpca_group_scad", "structured_group_lasso")
for (p_value in as.integer(names(expected_repetitions))) {
  subset <- replicates[replicates$p == p_value, , drop = FALSE]
  expected_methods <- if (p_value <= 20L) low_methods else high_methods
  counts <- table(subset$method_id)
  if (!setequal(names(counts), expected_methods) ||
      any(counts != expected_repetitions[[as.character(p_value)]])) {
    fair_fail("Unexpected method or repetition count at p=%d", p_value)
  }
}

expected_seed <- ifelse(
  replicates$p <= 20L,
  2026101000L + replicates$p * 1000L + replicates$rep,
  2026110000L + replicates$p * 1000L + replicates$rep
)
if (any(replicates$seed != expected_seed) ||
    any(replicates$fold_seed != expected_seed + 900000L)) {
  fair_fail("Fair-baseline DGP or fold seeds do not follow the registered rule")
}
if (any(replicates$validation_loss != "response-coefficient RMSE")) {
  fair_fail("Fair-baseline validation loss is inconsistent")
}

selection_methods <- c(
  "basis_group_scad", "fpca_group_lasso", "fpca_group_scad",
  "fpca_aenet", "structured_group_lasso"
)
selection_rows <- replicates[replicates$method_id %in% selection_methods, , drop = FALSE]
probability_columns <- c("covariate_tpr", "covariate_fdr")
if (any(!is.finite(as.matrix(selection_rows[probability_columns]))) ||
    any(as.matrix(selection_rows[probability_columns]) < 0) ||
    any(as.matrix(selection_rows[probability_columns]) > 1)) {
  fair_fail("Fair-baseline selection metrics are missing or outside [0,1]")
}
common_finite_columns <- c(
  "test_coeff_rmse", "tuning_runtime_sec", "final_runtime_sec", "runtime_sec"
)
if (any(!is.finite(as.matrix(replicates[common_finite_columns]))) ||
    any(as.matrix(replicates[common_finite_columns]) < 0)) {
  fair_fail("Fair-baseline prediction or timing metrics are missing or negative")
}
coefficient_rows <- replicates[replicates$method_id != "kernel_ridge", , drop = FALSE]
if (any(!is.finite(coefficient_rows$coefficient_relative_error)) ||
    any(coefficient_rows$coefficient_relative_error < 0)) {
  fair_fail("Coefficient-capable baselines have missing or negative coefficient error")
}
if (any(is.finite(replicates$coefficient_relative_error[
  replicates$method_id == "kernel_ridge"
]))) {
  fair_fail("Kernel ridge unexpectedly reports an explicit coefficient-surface error")
}

metric_columns <- c(
  "covariate_tpr", "covariate_fdr", "coefficient_relative_error",
  "test_coeff_rmse", "tuning_runtime_sec", "final_runtime_sec", "runtime_sec"
)
recomputed <- do.call(rbind, lapply(split(replicates, list(replicates$p, replicates$method_id)), function(block) {
  do.call(rbind, lapply(metric_columns, function(metric) {
    values <- block[[metric]]
    values <- values[is.finite(values)]
    if (!length(values)) return(NULL)
    data.frame(
      p = block$p[[1L]], method_id = block$method_id[[1L]], metric = metric,
      n = length(values), mean = mean(values), sd = stats::sd(values),
      mcse = stats::sd(values) / sqrt(length(values))
    )
  }))
}))
summary_key <- function(data) paste(data$p, data$method_id, data$metric, sep = "|")
recomputed <- recomputed[order(summary_key(recomputed)), , drop = FALSE]
summary_table <- summary_table[order(summary_key(summary_table)), , drop = FALSE]
if (!identical(summary_key(recomputed), summary_key(summary_table)) ||
    any(recomputed$n != summary_table$n) ||
    max(abs(recomputed$mean - summary_table$mean), na.rm = TRUE) > 1e-10 ||
    max(abs(recomputed$sd - summary_table$sd), na.rm = TRUE) > 1e-10 ||
    max(abs(recomputed$mcse - summary_table$mcse), na.rm = TRUE) > 1e-10) {
  fair_fail("Fair-baseline summary does not reproduce the replicate evidence")
}

manifest <- readRDS(manifest_file)
if (!identical(manifest$experiment_id, "baseline_v3_fair_common_fold_comparison") ||
    !identical(as.integer(manifest$design$p), as.integer(names(expected_repetitions))) ||
    !identical(as.integer(manifest$design$repetitions), unname(expected_repetitions)) ||
    !identical(manifest$design$n_folds, 5L) ||
    !identical(manifest$design$validation_loss, "response-coefficient RMSE per held-out curve") ||
    !identical(manifest$design$test_set_used_for_tuning, FALSE)) {
  fair_fail("Fair-baseline manifest does not match the registered design")
}

cat("Fair baseline validation passed: 2050/2050 fits are complete and reproducible.\n")
