# Formal robustness experiment with a paired reference and three difficult regimes.

library(data.table)

source("R/parameters.R")
source("R/fssgl/basis_design.R")
source("R/fssgl/penalty_weights.R")
source("R/fssgl/weighted_membership_solver.R")
source("R/fssgl/solver.R")
source("R/fssgl/simulation_design.R")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "results/tables/simulation/v2_formal")
processed_dir <- file.path(root, "data/processed/simulation/v2_formal")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

parameters <- fssgl_parameters()
scenarios <- data.table(
  scenario_id = c("baseline", "snr_low", "correlation_high", "within_group_sparse"),
  snr = c(2, 1, 2, 2),
  rho_common = c(0.10, 0.10, 0.20, 0.10),
  rho_group = c(0.30, 0.30, 0.55, 0.30),
  active_covariates_per_group = c(2L, 2L, 2L, 1L)
)
scenario_grid <- scenarios[, .(rep = seq_len(100L)), by = names(scenarios)]
checkpoint_file <- file.path(processed_dir, "fssgl_v2_robustness_running.rds")
completed <- if (file.exists(checkpoint_file)) readRDS(checkpoint_file) else data.table()
results <- list()
result_index <- 1L

for (i in seq_len(nrow(scenario_grid))) {
  scenario <- scenario_grid[i]
  if (
    nrow(completed) > 0L &&
      nrow(completed[scenario_id == scenario$scenario_id & rep == scenario$rep]) > 0L
  ) next
  seed <- 2026103000L + scenario$rep
  cat("FSSGL robustness: ", scenario$scenario_id, ", rep=", scenario$rep, "/100\n", sep = "")
  dgp <- generate_fssgl_simulation_dgp(
    n_train = 40L,
    n_test = 24L,
    n_covariates = 20L,
    kx = 4L,
    ky = 4L,
    n_groups = 4L,
    n_active_groups = 2L,
    n_active_covariates_per_group = scenario$active_covariates_per_group,
    rho_common = scenario$rho_common,
    rho_group = scenario$rho_group,
    snr = scenario$snr,
    surface_smoothness = "moderate",
    signal_scale = 1,
    seed = seed
  )
  started <- proc.time()[["elapsed"]]
  fit_result <- tryCatch(
    fit_fssgl(
      x_coef = dgp$x_train,
      y_coef = dgp$y_train,
      structural_membership = dgp$structural_membership,
      parameters = parameters,
      verbose = FALSE
    ),
    error = function(e) e
  )
  runtime_sec <- proc.time()[["elapsed"]] - started
  identification <- data.table(
    algorithm_version = FSSGL_ALGORITHM_VERSION,
    scenario_id = scenario$scenario_id,
    rep = scenario$rep,
    seed = seed,
    p = 20L,
    n_train = 40L,
    n_test = 24L,
    n_groups = 4L,
    active_group_count = 2L,
    active_covariate_count = length(dgp$active_covariates),
    active_covariates_per_group = scenario$active_covariates_per_group,
    snr = scenario$snr,
    rho_common = scenario$rho_common,
    rho_group = scenario$rho_group,
    noise_variance_true = dgp$noise_sd^2
  )
  if (inherits(fit_result, "error")) {
    result <- cbind(
      identification,
      data.table(status = "error", error_message = conditionMessage(fit_result))
    )
  } else {
    metrics <- evaluate_fssgl_simulation_fit(
      fit_result,
      dgp,
      posterior_cutoff = parameters$posterior_cutoff,
      runtime_sec = runtime_sec
    )
    result <- cbind(
      identification,
      data.table(
        status = "ok",
        error_message = NA_character_,
        theta_covariate_final = fit_result$fit$hyperparameters$theta_covariate_final,
        theta_group_final = fit_result$fit$hyperparameters$theta_group_final,
        sigma2_final = fit_result$fit$hyperparameters$sigma2_final
      ),
      metrics
    )
  }
  results[[result_index]] <- result
  result_index <- result_index + 1L
  saveRDS(rbindlist(c(list(completed), results), fill = TRUE), checkpoint_file)
}

replicates <- rbindlist(c(list(completed), results), fill = TRUE)
setorder(replicates, scenario_id, rep)
metric_columns <- c(
  "covariate_tpr", "covariate_fpr", "covariate_fdr", "selected_covariate_count",
  "group_tpr", "group_fpr", "group_fdr", "selected_group_count",
  "coefficient_relative_error", "active_coefficient_relative_error",
  "inactive_coefficient_norm", "test_coeff_rmse", "runtime_sec", "final_iter",
  "final_beta_change", "objective_tail_rel_change", "final_residual_variance",
  "theta_covariate_final", "theta_group_final", "sigma2_final",
  "noise_variance_true"
)
summary_long <- melt(
  replicates[status == "ok"],
  id.vars = "scenario_id",
  measure.vars = metric_columns,
  variable.name = "metric",
  value.name = "value"
)[, .(
  n = .N,
  mean = mean(value, na.rm = TRUE),
  sd = sd(value, na.rm = TRUE),
  median = median(value, na.rm = TRUE)
), by = .(scenario_id, metric)]
summary_wide <- dcast(
  summary_long[, mean_sd := sprintf("%.4f (%.4f)", mean, sd)],
  scenario_id ~ metric,
  value.var = "mean_sd"
)

paired_metrics <- c(
  "covariate_tpr", "covariate_fpr", "covariate_fdr", "selected_covariate_count",
  "group_tpr", "group_fpr", "coefficient_relative_error", "test_coeff_rmse"
)
paired <- merge(
  replicates[scenario_id != "baseline"],
  replicates[scenario_id == "baseline", c("rep", "seed", paired_metrics), with = FALSE],
  by = c("rep", "seed"),
  suffixes = c("", "_baseline")
)
paired_summary <- rbindlist(lapply(paired_metrics, function(metric) {
  paired[, {
    difference <- get(metric) - get(paste0(metric, "_baseline"))
    standard_error <- sd(difference) / sqrt(.N)
    list(
      n = .N,
      mean_difference = mean(difference),
      sd_difference = sd(difference),
      ci95_lower = mean(difference) - qt(0.975, .N - 1L) * standard_error,
      ci95_upper = mean(difference) + qt(0.975, .N - 1L) * standard_error
    )
  }, by = scenario_id][, metric := metric]
}))
convergence <- replicates[, .(
  n_attempted = .N,
  n_ok = sum(status == "ok"),
  strict_rate = mean(strict_converged[status == "ok"], na.rm = TRUE),
  error_count = sum(status != "ok")
), by = scenario_id]

fwrite(replicates, file.path(table_dir, "fssgl_v2_robustness_replicates.csv"))
fwrite(summary_long, file.path(table_dir, "fssgl_v2_robustness_summary_long.csv"))
fwrite(summary_wide, file.path(table_dir, "fssgl_v2_robustness_summary_mean_sd.csv"))
fwrite(paired_summary, file.path(table_dir, "fssgl_v2_robustness_paired_vs_baseline.csv"))
fwrite(convergence, file.path(table_dir, "fssgl_v2_robustness_convergence.csv"))
saveRDS(
  new_experiment_manifest(
    experiment_id = "fssgl_v2_robustness_formal",
    parameters = parameters,
    design = list(scenarios = scenarios, repetitions = 100L, paired_seeds = TRUE),
    algorithm_version = FSSGL_ALGORITHM_VERSION
  ),
  file.path(processed_dir, "fssgl_v2_robustness_manifest.rds")
)
if (nrow(replicates) == nrow(scenario_grid)) unlink(checkpoint_file)

print(summary_wide[, .(
  scenario_id, covariate_tpr, covariate_fpr, covariate_fdr,
  group_tpr, coefficient_relative_error, test_coeff_rmse
)])
print(convergence)
