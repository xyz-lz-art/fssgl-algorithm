# Formal frozen-parameter main simulation: 100 repetitions at p = 10 and 20.

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
dgp_defaults <- fssgl_main_dgp_defaults()
scenario_grid <- CJ(p = c(10L, 20L), rep = seq_len(100L))
checkpoint_file <- file.path(processed_dir, "fssgl_v2_main_running.rds")
completed <- if (file.exists(checkpoint_file)) readRDS(checkpoint_file) else data.table()
results <- list()
result_index <- 1L

for (i in seq_len(nrow(scenario_grid))) {
  scenario <- scenario_grid[i]
  if (
    nrow(completed) > 0L &&
      nrow(completed[p == scenario$p & rep == scenario$rep]) > 0L
  ) next
  n_groups <- if (scenario$p == 10L) 2L else 4L
  n_active_groups <- if (scenario$p == 10L) 1L else 2L
  seed <- 2026101000L + scenario$p * 1000L + scenario$rep
  cat("FSSGL main: p=", scenario$p, ", rep=", scenario$rep, "/100\n", sep = "")

  dgp <- do.call(generate_fssgl_simulation_dgp, c(dgp_defaults, list(
    n_covariates = scenario$p,
    n_groups = n_groups,
    n_active_groups = n_active_groups,
    n_active_covariates_per_group = 2L,
    seed = seed
  )))
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
    p = scenario$p,
    rep = scenario$rep,
    seed = seed,
    n_train = dgp$dimensions$n_train,
    n_test = dgp$dimensions$n_test,
    n_groups = n_groups,
    active_group_count = n_active_groups,
    active_covariate_count = length(dgp$active_covariates),
    snr = dgp$settings$snr,
    rho_common = dgp$settings$rho_common,
    rho_group = dgp$settings$rho_group,
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
setorder(replicates, p, rep)
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
  id.vars = "p",
  measure.vars = metric_columns,
  variable.name = "metric",
  value.name = "value"
)[, .(
  n = .N,
  mean = mean(value, na.rm = TRUE),
  sd = sd(value, na.rm = TRUE),
  median = median(value, na.rm = TRUE),
  q05 = quantile(value, 0.05, na.rm = TRUE),
  q95 = quantile(value, 0.95, na.rm = TRUE)
), by = .(p, metric)]
summary_wide <- dcast(
  summary_long[, mean_sd := sprintf("%.4f (%.4f)", mean, sd)],
  p ~ metric,
  value.var = "mean_sd"
)
convergence <- replicates[, .(
  n_attempted = .N,
  n_ok = sum(status == "ok"),
  n_strict = sum(strict_converged %in% TRUE, na.rm = TRUE),
  strict_rate = mean(strict_converged[status == "ok"], na.rm = TRUE),
  relaxed_rate = mean(relaxed_converged[status == "ok"], na.rm = TRUE),
  error_count = sum(status != "ok")
), by = p]

fwrite(replicates, file.path(table_dir, "fssgl_v2_main_replicates.csv"))
fwrite(summary_long, file.path(table_dir, "fssgl_v2_main_summary_long.csv"))
fwrite(summary_wide, file.path(table_dir, "fssgl_v2_main_summary_mean_sd.csv"))
fwrite(convergence, file.path(table_dir, "fssgl_v2_main_convergence.csv"))
saveRDS(
  new_experiment_manifest(
    experiment_id = "fssgl_v2_main_formal",
    parameters = parameters,
    design = list(dgp_defaults = dgp_defaults, scenario_grid = scenario_grid),
    algorithm_version = FSSGL_ALGORITHM_VERSION
  ),
  file.path(processed_dir, "fssgl_v2_main_manifest.rds")
)
if (nrow(replicates) == nrow(scenario_grid)) unlink(checkpoint_file)

print(summary_wide[, .(
  p, covariate_tpr, covariate_fpr, covariate_fdr,
  coefficient_relative_error, test_coeff_rmse, runtime_sec
)])
print(convergence)
