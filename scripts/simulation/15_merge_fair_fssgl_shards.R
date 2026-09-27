# Merge completed deterministic shards without rerunning or changing any fit.
library(data.table)

shard_count <- as.integer(Sys.getenv("FSSGL_FAIR_SHARD_COUNT", "5"))
if (length(shard_count) != 1L || is.na(shard_count) || shard_count < 2L) {
  stop("FSSGL_FAIR_SHARD_COUNT must be at least two for merging.")
}
table_dir <- "results/tables/simulation/v3_fair"
data_dir <- "data/processed/simulation/v3_fair"
suffixes <- paste0("_shard", seq_len(shard_count), "of", shard_count)
artifact_paths <- function(stem, extension, directory) {
  file.path(directory, paste0(stem, suffixes, extension))
}
result_paths <- artifact_paths("fssgl_v3_fair_replicates", ".csv", table_dir)
cv_paths <- artifact_paths("fssgl_v3_fair_cv_path", ".csv", table_dir)
manifest_paths <- artifact_paths("fssgl_v3_fair_manifest", ".rds", data_dir)
for (path in c(result_paths, cv_paths, manifest_paths)) {
  if (!file.exists(path)) stop("Missing completed shard artifact: ", path)
}

manifests <- lapply(manifest_paths, readRDS)
for (manifest in manifests[-1L]) {
  if (!identical(manifest$experiment_id, manifests[[1L]]$experiment_id) ||
      !identical(manifest$parameters, manifests[[1L]]$parameters) ||
      !identical(manifest$design, manifests[[1L]]$design)) {
    stop("Shard manifests disagree on the experiment design.")
  }
}
results <- rbindlist(lapply(result_paths, fread), use.names = TRUE, fill = TRUE)
expected_repetitions <- c(`10` = 100L, `20` = 100L, `50` = 50L, `60` = 50L, `100` = 50L)
if (nrow(results) != sum(expected_repetitions) ||
    anyDuplicated(results, by = c("p", "rep"))) {
  stop("Shards do not contain 350 unique p/rep results.")
}
for (p_value in as.integer(names(expected_repetitions))) {
  reps <- sort(results[p == p_value, rep])
  if (!identical(reps, seq_len(expected_repetitions[[as.character(p_value)]]))) {
    stop("Incomplete repetitions at p=", p_value)
  }
}
cv_path <- rbindlist(lapply(cv_paths, fread), use.names = TRUE, fill = TRUE)
if (anyDuplicated(cv_path, by = c("p", "rep", "fold_id", "multiplier"))) {
  stop("Duplicate p/rep/fold/multiplier CV rows across shards.")
}
setorder(results, p, rep)
if (nrow(cv_path)) setorder(cv_path, p, rep, fold_id, multiplier)
fwrite(results, file.path(table_dir, "fssgl_v3_fair_replicates.csv"))
fwrite(cv_path, file.path(table_dir, "fssgl_v3_fair_cv_path.csv"))

reportable <- results[status == "ok" & strict_converged == TRUE]
metric_columns <- c(
  "covariate_tpr", "covariate_fdr", "coefficient_relative_error", "test_coeff_rmse",
  "tuning_runtime_sec", "final_runtime_sec", "total_runtime_sec"
)
summary <- if (nrow(reportable)) {
  melt(
    reportable,
    id.vars = "p",
    measure.vars = metric_columns,
    variable.name = "metric",
    value.name = "value"
  )[, .(
    n = .N,
    mean = mean(value),
    sd = stats::sd(value),
    mcse = stats::sd(value) / sqrt(.N)
  ), by = .(p, metric)]
} else {
  data.table()
}
fwrite(summary, file.path(table_dir, "fssgl_v3_fair_summary.csv"))

baseline <- fread(file.path(table_dir, "baseline_v3_fair_replicates.csv"))[
  method_id %in% c("basis_group_scad", "fpca_group_scad", "structured_group_lasso")
]
paired <- merge(
  baseline[, c("p", "rep", "method_id", metric_columns[1L:4L]), with = FALSE],
  results[, c("p", "rep", metric_columns[1L:4L]), with = FALSE],
  by = c("p", "rep"),
  suffixes = c("_baseline", "_fssgl")
)
if (nrow(paired) != 3L * nrow(results)) {
  stop("The direct baselines and FSSGL are not paired on every repetition.")
}
paired_summary <- rbindlist(lapply(metric_columns[1L:4L], function(metric) {
  differences <- paired[[paste0(metric, "_baseline")]] -
    paired[[paste0(metric, "_fssgl")]]
  paired[, value := differences][, .(
    metric = metric,
    n = .N,
    mean_baseline_minus_fssgl = mean(value),
    sd = stats::sd(value),
    mcse = stats::sd(value) / sqrt(.N)
  ), by = .(p, method_id)]
}))
paired_summary[, `:=`(
  ci_low = mean_baseline_minus_fssgl - 1.96 * mcse,
  ci_high = mean_baseline_minus_fssgl + 1.96 * mcse
)]
setorder(paired_summary, p, method_id, metric)
fwrite(
  paired_summary,
  file.path(table_dir, "fssgl_v3_fair_paired_differences.csv")
)

manifest <- manifests[[1L]]
manifest$design$execution_shards <- shard_count
saveRDS(manifest, file.path(data_dir, "fssgl_v3_fair_manifest.rds"))
print(results[, .N, by = .(p, status)])
cat("Merged ", shard_count, " fair-FSSGL shards.\n", sep = "")
