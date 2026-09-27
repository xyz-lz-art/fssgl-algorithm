# Merge the main experiment and its common-seed baseline refits.

library(data.table)

table_dir <- "results/tables/simulation/v2_formal"
fssgl <- fread(file.path(table_dir, "fssgl_v2_main_replicates.csv"))
baselines <- fread(file.path(table_dir, "functional_baselines_v2_replicates.csv"))
fssgl[, `:=`(
  method_id = "fssgl_v2",
  supports_selection = TRUE,
  supports_coefficients = TRUE,
  solver_converged = strict_converged
)]

metric_columns <- c(
  "covariate_tpr", "covariate_fpr", "covariate_fdr",
  "selected_covariate_count", "group_tpr", "group_fpr", "group_fdr",
  "selected_group_count", "coefficient_relative_error",
  "active_coefficient_relative_error", "inactive_coefficient_norm",
  "test_coeff_rmse", "runtime_sec", "final_iter"
)
common_columns <- unique(c(
  "method_id", "p", "rep", "seed", "status", "supports_selection",
  "supports_coefficients", "solver_converged", metric_columns
))
comparison <- rbindlist(list(
  fssgl[, intersect(common_columns, names(fssgl)), with = FALSE],
  baselines[, intersect(common_columns, names(baselines)), with = FALSE]
), fill = TRUE)
comparison[, (metric_columns) := lapply(.SD, as.numeric), .SDcols = metric_columns]

expected_methods <- c(
  "fssgl_v2", "basis_ridge", "fpca_ridge", "fpca_group_lasso",
  "fpca_aenet", "fpca_group_scad", "structured_group_lasso", "kernel_ridge"
)
completion <- comparison[, .(
  n_attempted = .N,
  n_ok = sum(status == "ok"),
  error_count = sum(status != "ok"),
  solver_convergence_rate = mean(solver_converged[!is.na(solver_converged)], na.rm = TRUE)
), by = .(method_id, p)]
if (
  nrow(completion) != length(expected_methods) * 2L ||
    any(completion$n_ok != 100L) ||
    !setequal(completion$method_id, expected_methods)
) {
  print(completion)
  stop("The formal baseline comparison is incomplete.")
}

long <- melt(
  comparison[status == "ok"],
  id.vars = c("method_id", "p", "rep", "seed"),
  measure.vars = metric_columns,
  variable.name = "metric",
  value.name = "value"
)
summary_long <- long[, .(
  n_available = sum(!is.na(value)),
  mean = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE),
  sd = if (sum(!is.na(value)) <= 1L) NA_real_ else sd(value, na.rm = TRUE)
), by = .(method_id, p, metric)]
summary_wide <- dcast(
  summary_long[, mean_sd := fifelse(
    is.na(mean), NA_character_, sprintf("%.4f (%.4f)", mean, sd)
  )],
  method_id + p ~ metric,
  value.var = "mean_sd"
)

reference <- long[method_id == "fssgl_v2", .(
  p, rep, seed, metric, reference_value = value
)]
paired <- merge(
  long[method_id != "fssgl_v2"],
  reference,
  by = c("p", "rep", "seed", "metric")
)[!is.na(value) & !is.na(reference_value)]
paired[, difference := value - reference_value]
paired_summary <- paired[, {
  difference_sd <- sd(difference)
  difference_se <- difference_sd / sqrt(.N)
  .(
    n_pairs = .N,
    mean_difference = mean(difference),
    sd_difference = difference_sd,
    ci95_lower = mean(difference) - qt(0.975, .N - 1L) * difference_se,
    ci95_upper = mean(difference) + qt(0.975, .N - 1L) * difference_se
  )
}, by = .(method_id, p, metric)]

fwrite(comparison, file.path(table_dir, "functional_method_comparison_v2_replicates.csv"))
fwrite(summary_long, file.path(table_dir, "functional_method_comparison_v2_summary_long.csv"))
fwrite(summary_wide, file.path(table_dir, "functional_method_comparison_v2_summary_mean_sd.csv"))
fwrite(paired_summary, file.path(table_dir, "functional_method_comparison_v2_paired_differences.csv"))
fwrite(completion, file.path(table_dir, "functional_method_comparison_v2_completion.csv"))

print(summary_wide[, .(
  method_id, p, covariate_tpr, covariate_fpr, covariate_fdr,
  coefficient_relative_error, test_coeff_rmse
)])
print(completion)
