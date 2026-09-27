# Ablation experiments for the Shanghai real-data study.
# Variants:
#   1. covariate-only
#   2. covariate+group
#
# Default output is limited to summary CSVs under
# data/processed/shanghai_metroflow/orthonormal/ablation/. Set write_detail_outputs <- TRUE
# only when per-target fit objects/posteriors are needed for diagnostics.

library(data.table)

source("R/parameters.R")
source("R/fssgl/basis_design.R")
source("R/fssgl/penalty_weights.R")
source("R/fssgl/weighted_membership_solver.R")
source("R/application/shanghai_workflow.R")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
processed_dir <- file.path(root, "data/processed/shanghai_metroflow")
metadata_dir <- file.path(root, "data/metadata/shanghai_metroflow")
geometry_dir <- file.path(processed_dir, "orthonormal")
out_dir <- file.path(geometry_dir, "ablation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
main_summary_file <- file.path(geometry_dir, "multi_targets", "multi_targets_summary.csv")
main_summary <- if (file.exists(main_summary_file)) fread(main_summary_file) else data.table()
checkpoint_file <- file.path(out_dir, "ablation_running.rds")
completed_file <- file.path(out_dir, "ablation_completed.rds")
checkpoint_state <- if (file.exists(checkpoint_file)) readRDS(checkpoint_file) else NULL
completed <- if (
  is.list(checkpoint_state) &&
    identical(checkpoint_state$preprocessing_version, SHANGHAI_PREPROCESSING_VERSION)
) checkpoint_state$summary else data.table()

write_detail_outputs <- FALSE

targets <- shanghai_target_names()
base_params <- shanghai_main_parameters()

ablation_grid <- data.table(
  ablation_id = c(
    "A_covariate_only",
    "B_covariate_group"
  ),
  use_group_penalty = c(FALSE, TRUE),
  line_penalty_scale = c(0, base_params$line_penalty_scale)
)

summaries <- if (nrow(completed)) list(completed) else list()
counter <- length(summaries) + 1L

for (a_i in seq_len(nrow(ablation_grid))) {
  ablation <- as.list(ablation_grid[a_i])
  if (write_detail_outputs) {
    variant_dir <- file.path(out_dir, ablation$ablation_id)
    dir.create(variant_dir, recursive = TRUE, showWarnings = FALSE)
  }
  params <- base_params
  params$line_penalty_scale <- ablation$line_penalty_scale
  params$use_group_penalty <- ablation$use_group_penalty
  params$config_id <- ablation$ablation_id

  for (target in targets) {
    target_name <- target
    if (nrow(completed) &&
        nrow(completed[ablation_id == ablation$ablation_id & target == target_name])) {
      cat("Skipping completed", ablation$ablation_id, target, "\n")
      next
    }
    cat("Ablation", ablation$ablation_id, "target", target, "\n")
    if (ablation$use_group_penalty && nrow(main_summary)) {
      main_row <- main_summary[target == target_name]
      if (nrow(main_row) != 1L || main_row$status != "ok") {
        stop("Current main result is missing or invalid for ", target)
      }
      summaries[[counter]] <- data.table(
        ablation_id = ablation$ablation_id,
        target = target,
        status = "ok",
        use_group_penalty = TRUE,
        final_iter = main_row$final_iter,
        final_objective = main_row$final_objective,
        final_beta_change = main_row$final_beta_change,
        objective_tail_rel_change = main_row$objective_tail_rel_change,
        strict_converged = main_row$strict_converged,
        relaxed_converged = main_row$relaxed_converged,
        selected_station_count = main_row$selected_station_count,
        diagnostic_line_count = main_row$selected_line_count,
        selected_line_count = main_row$selected_line_count,
        selected_lines = main_row$selected_lines,
        top_line = main_row$top_line,
        top_line_prob = main_row$top_line_prob,
        top_station = main_row$top_station,
        top_station_prob = main_row$top_station_prob,
        final_residual_variance = main_row$final_residual_variance,
        test_coeff_rmse = main_row$test_coeff_rmse,
        test_curve_rmse = main_row$test_curve_rmse,
        test_curve_mae = main_row$test_curve_mae
      )
      counter <- counter + 1L
      saveRDS(list(
        preprocessing_version = SHANGHAI_PREPROCESSING_VERSION,
        summary = rbindlist(summaries, fill = TRUE)
      ), checkpoint_file)
      next
    }
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
      summaries[[counter]] <- data.table(
        ablation_id = ablation$ablation_id,
        target = target,
        status = "error",
        error_message = conditionMessage(result)
      )
      counter <- counter + 1L
      saveRDS(list(
        preprocessing_version = SHANGHAI_PREPROCESSING_VERSION,
        summary = rbindlist(summaries, fill = TRUE)
      ), checkpoint_file)
      next
    }

    fit <- result$fit
    hist <- fit$history
    final <- hist[.N]
    diagnostics <- evaluate_shanghai_fit(result)
    selected_lines <- fit$line_posterior[selected == TRUE, line]

    if (write_detail_outputs) {
      stem <- clean_shanghai_name(target)
      saveRDS(
        list(
          target = result$design_obj$target,
          kx = result$design_obj$kx,
          ky = result$design_obj$ky,
          params = params,
          fit = fit
        ),
        file.path(variant_dir, paste0("fssgl_", stem, ".rds"))
      )
      fwrite(hist, file.path(variant_dir, paste0("fit_history_", stem, ".csv")))
      fwrite(fit$station_posterior, file.path(variant_dir, paste0("station_posterior_", stem, ".csv")))
      fwrite(fit$line_posterior, file.path(variant_dir, paste0("line_posterior_", stem, ".csv")))
    }

    summaries[[counter]] <- data.table(
      ablation_id = ablation$ablation_id,
      target = target,
      status = "ok",
      use_group_penalty = ablation$use_group_penalty,
      final_iter = final$iter,
      final_objective = final$objective,
      final_beta_change = final$beta_change,
      objective_tail_rel_change = fit$convergence$objective_tail_rel_change,
      strict_converged = fit$convergence$strict_converged,
      relaxed_converged = fit$convergence$relaxed_converged,
      selected_station_count = final$selected_station_count,
      diagnostic_line_count = final$selected_line_count,
      selected_line_count = if (ablation$use_group_penalty) final$selected_line_count else NA_integer_,
      selected_lines = if (ablation$use_group_penalty) paste(selected_lines, collapse = ";") else NA_character_,
      top_line = fit$line_posterior[1, line],
      top_line_prob = fit$line_posterior[1, posterior_slab_prob],
      top_station = fit$station_posterior[1, name],
      top_station_prob = fit$station_posterior[1, posterior_slab_prob],
      final_residual_variance = fit$convergence$final_residual_variance,
      test_coeff_rmse = diagnostics$test_coeff_rmse,
      test_curve_rmse = diagnostics$test_curve_rmse,
      test_curve_mae = diagnostics$test_curve_mae
    )
    counter <- counter + 1L
    saveRDS(list(
      preprocessing_version = SHANGHAI_PREPROCESSING_VERSION,
      summary = rbindlist(summaries, fill = TRUE)
    ), checkpoint_file)
  }
}

summary <- rbindlist(summaries, fill = TRUE)
summary[, converged_beta := final_beta_change < base_params$relaxed_tol]
summary[, stable_objective_tail := objective_tail_rel_change < base_params$objective_tail_tol]
fwrite(summary, file.path(out_dir, "ablation_summary.csv"))

variant_summary <- summary[
  status == "ok",
  .(
    n_targets = .N,
    all_relaxed_converged = all(relaxed_converged),
    mean_beta_change = mean(final_beta_change),
    max_beta_change = max(final_beta_change),
    mean_tail_change = mean(objective_tail_rel_change),
    max_tail_change = max(objective_tail_rel_change),
    mean_station_count = mean(selected_station_count),
    mean_line_count = if (all(is.na(selected_line_count))) NA_real_ else mean(selected_line_count, na.rm = TRUE),
    mean_diagnostic_line_count = mean(diagnostic_line_count),
    mean_residual_variance = mean(final_residual_variance),
    mean_test_curve_rmse = mean(test_curve_rmse),
    max_test_curve_rmse = max(test_curve_rmse),
    mean_test_curve_mae = mean(test_curve_mae)
  ),
  by = .(ablation_id, use_group_penalty)
]
setorder(variant_summary, use_group_penalty)
fwrite(variant_summary, file.path(out_dir, "ablation_variant_summary.csv"))
saveRDS(
  list(
    preprocessing_version = SHANGHAI_PREPROCESSING_VERSION,
    summary = summary
  ),
  completed_file
)
if (file.exists(checkpoint_file)) unlink(checkpoint_file)
saveRDS(
  new_experiment_manifest(
    experiment_id = "shanghai_structural_ablation",
    parameters = base_params,
    design = list(
      targets = targets,
      variants = ablation_grid,
      kx = 5L,
      ky = 5L,
      basis_geometry = FSSGL_BASIS_GEOMETRY_VERSION,
      preprocessing_version = SHANGHAI_PREPROCESSING_VERSION
    )
  ),
  file.path(out_dir, "ablation_manifest.rds")
)

cat("\nAblation variant summary\n")
print(variant_summary)
