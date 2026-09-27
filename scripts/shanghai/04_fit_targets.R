# Fit the current FSSGL model on the 12-target Shanghai set.
# Summary outputs are written for all targets. A complete fit object is retained
# for the prespecified surface target so that manuscript figures are reproducible.

library(data.table)

source("R/parameters.R")
source("R/fssgl/basis_design.R")
source("R/fssgl/penalty_weights.R")
source("R/fssgl/weighted_membership_solver.R")
source("R/application/shanghai_workflow.R")

root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
processed_dir <- file.path(root, "data/processed/shanghai_metroflow")
metadata_dir <- file.path(root, "data/metadata/shanghai_metroflow")
geometry_dir <- file.path(processed_dir, "orthonormal")
out_dir <- file.path(geometry_dir, "multi_targets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

detail_targets <- "Pudong International Airport"

targets <- shanghai_target_names()
params <- shanghai_main_parameters()

summaries <- list()

for (target in targets) {
  cat("Fitting multi-target", target, "\n")
  result <- tryCatch(
    fit_shanghai_target_fssgl(
      target_name = target,
      params = params,
      processed_dir = processed_dir,
      metadata_dir = metadata_dir,
      kx = 5,
      ky = 5
    ),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    summaries[[target]] <- data.table(
      config_id = params$config_id,
      target = target,
      status = "error",
      error_message = conditionMessage(result)
    )
    next
  }

  fit <- result$fit
  hist <- fit$history
  final <- hist[.N]
  diagnostics <- evaluate_shanghai_fit(result)

  if (target %in% detail_targets) {
    stem <- clean_shanghai_name(target)
    saveRDS(
      list(
        target = result$design_obj$target,
        kx = result$design_obj$kx,
        ky = result$design_obj$ky,
        basis_x = result$design_obj$basis_x,
        basis_y = result$design_obj$basis_y,
        basis_geometry = result$design_obj$basis_geometry,
        params = params,
        fit = fit
      ),
      file.path(out_dir, paste0("fssgl_", stem, ".rds"))
    )

    fwrite(fit$station_posterior, file.path(out_dir, paste0("station_posterior_", stem, ".csv")))
    fwrite(fit$line_posterior, file.path(out_dir, paste0("line_posterior_", stem, ".csv")))
    fwrite(hist, file.path(out_dir, paste0("fit_history_", stem, ".csv")))
    fwrite(fit$station_posterior[selected == TRUE], file.path(out_dir, paste0("selected_stations_", stem, ".csv")))
    fwrite(fit$line_posterior[selected == TRUE], file.path(out_dir, paste0("selected_lines_", stem, ".csv")))
  }

  summaries[[target]] <- data.table(
    config_id = params$config_id,
    target = target,
    status = "ok",
    final_iter = final$iter,
    final_objective = final$objective,
    final_beta_change = final$beta_change,
    objective_tail_rel_change = fit$convergence$objective_tail_rel_change,
    strict_converged = fit$convergence$strict_converged,
    relaxed_converged = fit$convergence$relaxed_converged,
    selected_line_count = final$selected_line_count,
    selected_station_count = final$selected_station_count,
    selected_lines = paste(fit$line_posterior[selected == TRUE, line], collapse = ";"),
    top_line = fit$line_posterior[1, line],
    top_line_prob = fit$line_posterior[1, posterior_slab_prob],
    top_station = fit$station_posterior[1, name],
    top_station_prob = fit$station_posterior[1, posterior_slab_prob],
    final_residual_variance = fit$convergence$final_residual_variance,
    test_coeff_rmse = diagnostics$test_coeff_rmse,
    test_curve_rmse = diagnostics$test_curve_rmse,
    test_curve_mae = diagnostics$test_curve_mae
  )
}

summary <- rbindlist(summaries, fill = TRUE)
summary[, converged_beta := final_beta_change < params$relaxed_tol]
summary[, stable_objective_tail := objective_tail_rel_change < params$objective_tail_tol]
fwrite(summary, file.path(out_dir, "multi_targets_summary.csv"))
fwrite(as.data.table(params), file.path(out_dir, "multi_targets_params.csv"))
saveRDS(
  new_experiment_manifest(
    experiment_id = "shanghai_multi_target",
    parameters = params,
    design = list(
      targets = targets,
      kx = 5L,
      ky = 5L,
      basis_geometry = FSSGL_BASIS_GEOMETRY_VERSION,
      preprocessing_version = SHANGHAI_PREPROCESSING_VERSION,
      split = "chronological_80_20"
    )
  ),
  file.path(out_dir, "multi_targets_manifest.rds")
)

cat("\nMulti-target summary\n")
print(summary)
