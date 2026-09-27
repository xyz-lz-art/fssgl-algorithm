# Expanding-window rolling evaluation for all 12 Shanghai targets.
# The first fit uses 60 training days; each origin predicts the next seven days,
# and the origin advances by seven days until the end of the observed period.
# Fixed parameters were previously calibrated on the 98-day development set;
# early origins are retrospective diagnostics, not nested out-of-time tuning.
# Evaluation conditions on the complete same-day predictor curves.

library(data.table)

source("R/parameters.R")
source("R/fssgl/basis_design.R")
source("R/fssgl/penalty_weights.R")
source("R/fssgl/weighted_membership_solver.R")
source("R/application/shanghai_workflow.R")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
processed_dir <- file.path(root, "data/processed/shanghai_metroflow")
metadata_dir <- file.path(root, "data/metadata/shanghai_metroflow")
out_dir <- file.path(processed_dir, "orthonormal/rolling_evaluation")
table_dir <- file.path(root, "results/tables/shanghai_metroflow/orthonormal/rolling_evaluation")
figure_dir <- file.path(root, "results/figures/shanghai_metroflow/orthonormal/rolling_evaluation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

targets <- shanghai_target_names()
params <- shanghai_main_parameters()
initial_train_days <- 60L
horizon_days <- 7L
step_days <- 7L
origins <- seq(initial_train_days, 123L - horizon_days, by = step_days)

rolling_design <- rbindlist(lapply(seq_along(origins), function(window_id) {
  origin <- origins[window_id]
  data.table(
    window_id = window_id,
    train_start_day = 1L,
    train_end_day = origin,
    n_train_days = origin,
    test_start_day = origin + 1L,
    test_end_day = origin + horizon_days,
    n_test_days = horizon_days
  )
}))
fwrite(rolling_design, file.path(out_dir, "rolling_window_design.csv"))

checkpoint_file <- file.path(out_dir, "rolling_running.rds")
completed_file <- file.path(out_dir, "rolling_completed.rds")
empty_state <- list(
  preprocessing_version = SHANGHAI_PREPROCESSING_VERSION,
  fit_summary = data.table(),
  selected_stations = data.table(),
  selected_lines = data.table()
)
state_candidate <- if (file.exists(checkpoint_file)) {
  readRDS(checkpoint_file)
} else if (file.exists(completed_file)) {
  readRDS(completed_file)
} else {
  empty_state
}
state <- if (
  is.list(state_candidate) &&
    identical(state_candidate$preprocessing_version, SHANGHAI_PREPROCESSING_VERSION)
) state_candidate else empty_state

shanghai_data <- load_shanghai_data(processed_dir, metadata_dir)
date_values <- as.character(shanghai_data$date_index$date)
save_checkpoint <- function() saveRDS(state, checkpoint_file)

for (target_name in targets) {
  for (window_id in rolling_design$window_id) {
    current_window <- window_id
    already_done <- nrow(state$fit_summary) > 0L && nrow(state$fit_summary[
      target == target_name & window_id == current_window
    ]) > 0L
    if (already_done) {
      cat("Skipping completed rolling fit", target_name, window_id, "\n")
      next
    }

    window <- rolling_design[window_id == current_window]
    train_days <- seq.int(window$train_start_day, window$train_end_day)
    test_days <- seq.int(window$test_start_day, window$test_end_day)
    if (max(train_days) >= min(test_days)) stop("Rolling window leaks future days into training.")

    cat("Rolling fit", target_name, "window", window_id, "\n")
    result <- tryCatch(
      fit_shanghai_target_fssgl(
        target_name = target_name,
        params = params,
        processed_dir = processed_dir,
        metadata_dir = metadata_dir,
        kx = 5,
        ky = 5,
        train_days = train_days,
        test_days = test_days,
        shanghai_data = shanghai_data
      ),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      state$fit_summary <- rbind(
        state$fit_summary,
        data.table(
          target = target_name,
          window_id = window_id,
          status = "error",
          error_message = conditionMessage(result)
        ),
        fill = TRUE
      )
      save_checkpoint()
      next
    }

    fit <- result$fit
    final <- fit$history[.N]
    diagnostics <- evaluate_shanghai_fit(result)
    y_sd <- result$target_obj$scaling$y_sd
    state$fit_summary <- rbind(
      state$fit_summary,
      data.table(
        target = target_name,
        window_id = window_id,
        status = "ok",
        train_start_date = date_values[min(train_days)],
        train_end_date = date_values[max(train_days)],
        test_start_date = date_values[min(test_days)],
        test_end_date = date_values[max(test_days)],
        n_train_days = length(train_days),
        n_test_days = length(test_days),
        final_iter = final$iter,
        final_beta_change = final$beta_change,
        objective_tail_rel_change = fit$convergence$objective_tail_rel_change,
        strict_converged = fit$convergence$strict_converged,
        relaxed_converged = fit$convergence$relaxed_converged,
        selected_station_count = final$selected_station_count,
        selected_line_count = final$selected_line_count,
        test_coeff_rmse = diagnostics$test_coeff_rmse,
        test_curve_rmse_standardized = diagnostics$test_curve_rmse,
        test_curve_mae_standardized = diagnostics$test_curve_mae,
        test_curve_rmse_original = diagnostics$test_curve_rmse * y_sd,
        test_curve_mae_original = diagnostics$test_curve_mae * y_sd,
        training_response_sd = y_sd,
        final_residual_variance = fit$convergence$final_residual_variance
      ),
      fill = TRUE
    )

    selected_stations <- fit$station_posterior[selected == TRUE, .(
      target = target_name,
      window_id = window_id,
      station,
      name,
      lines,
      posterior_slab_prob,
      beta_norm
    )]
    selected_lines <- fit$line_posterior[selected == TRUE, .(
      target = target_name,
      window_id = window_id,
      line,
      posterior_slab_prob,
      line_norm
    )]
    state$selected_stations <- rbind(state$selected_stations, selected_stations, fill = TRUE)
    state$selected_lines <- rbind(state$selected_lines, selected_lines, fill = TRUE)
    save_checkpoint()
  }
}

fit_summary <- state$fit_summary
selected_stations <- state$selected_stations
selected_lines <- state$selected_lines
ok <- fit_summary[status == "ok"]

target_summary <- ok[, .(
  n_windows = .N,
  strict_convergence_rate = mean(strict_converged),
  relaxed_convergence_rate = mean(relaxed_converged),
  mean_test_curve_rmse_standardized = mean(test_curve_rmse_standardized),
  sd_test_curve_rmse_standardized = sd(test_curve_rmse_standardized),
  max_test_curve_rmse_standardized = max(test_curve_rmse_standardized),
  mean_test_curve_mae_standardized = mean(test_curve_mae_standardized),
  mean_test_curve_rmse_original = mean(test_curve_rmse_original),
  mean_test_curve_mae_original = mean(test_curve_mae_original),
  mean_selected_stations = mean(selected_station_count),
  sd_selected_stations = sd(selected_station_count),
  mean_selected_lines = mean(selected_line_count),
  mean_iterations = mean(final_iter)
), by = target]
setorder(target_summary, mean_test_curve_rmse_standardized)

window_summary <- ok[, .(
  n_targets = .N,
  all_strict_converged = all(strict_converged),
  mean_test_curve_rmse_standardized = mean(test_curve_rmse_standardized),
  mean_test_curve_mae_standardized = mean(test_curve_mae_standardized),
  mean_selected_stations = mean(selected_station_count),
  mean_selected_lines = mean(selected_line_count)
), by = .(window_id, train_end_date, test_start_date, test_end_date, n_train_days)]
setorder(window_summary, window_id)

overall_summary <- ok[, .(
  n_fits = .N,
  n_targets = uniqueN(target),
  n_windows = uniqueN(window_id),
  strict_convergence_rate = mean(strict_converged),
  relaxed_convergence_rate = mean(relaxed_converged),
  mean_test_curve_rmse_standardized = mean(test_curve_rmse_standardized),
  median_test_curve_rmse_standardized = median(test_curve_rmse_standardized),
  max_test_curve_rmse_standardized = max(test_curve_rmse_standardized),
  mean_test_curve_mae_standardized = mean(test_curve_mae_standardized),
  mean_selected_stations = mean(selected_station_count),
  mean_selected_lines = mean(selected_line_count)
)]

station_frequency <- selected_stations[, .(
  selected_windows = uniqueN(window_id),
  mean_slab_probability_when_selected = mean(posterior_slab_prob),
  mean_beta_norm_when_selected = mean(beta_norm)
), by = .(target, station, name, lines)]
station_frequency[, selection_frequency := selected_windows / length(origins)]
setorder(station_frequency, target, -selection_frequency)

line_frequency <- selected_lines[, .(
  selected_windows = uniqueN(window_id),
  mean_slab_probability_when_selected = mean(posterior_slab_prob),
  mean_line_norm_when_selected = mean(line_norm)
), by = .(target, line)]
line_frequency[, selection_frequency := selected_windows / length(origins)]
setorder(line_frequency, target, -selection_frequency)

fwrite(fit_summary, file.path(out_dir, "rolling_fit_summary.csv"))
fwrite(selected_stations, file.path(out_dir, "rolling_selected_stations.csv"))
fwrite(selected_lines, file.path(out_dir, "rolling_selected_lines.csv"))
fwrite(target_summary, file.path(table_dir, "rolling_target_summary.csv"))
fwrite(window_summary, file.path(table_dir, "rolling_window_summary.csv"))
fwrite(overall_summary, file.path(table_dir, "rolling_overall_summary.csv"))
fwrite(station_frequency, file.path(table_dir, "rolling_station_frequency.csv"))
fwrite(line_frequency, file.path(table_dir, "rolling_line_frequency.csv"))

png(file.path(figure_dir, "rolling_error_summary.png"), width = 2000, height = 1100, res = 160)
old_par <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(5, 14, 4, 2))
boxplot(test_curve_rmse_standardized ~ target, data = ok, las = 1,
        horizontal = TRUE, cex.axis = 0.8, ylab = "",
        col = "#1B6CA8", xlab = "Standardized curve RMSE", main = "Rolling Error by Target")
par(mar = c(5, 5, 4, 2))
plot(window_summary$window_id, window_summary$mean_test_curve_rmse_standardized,
     type = "b", pch = 19, col = "#C73E1D", xlab = "Rolling window",
     ylab = "Mean standardized curve RMSE", main = "Error across Time")
grid(col = "#DDDDDD")
par(old_par)
dev.off()

saveRDS(state, completed_file)
if (file.exists(checkpoint_file)) unlink(checkpoint_file)
saveRDS(
  new_experiment_manifest(
    experiment_id = "shanghai_expanding_window_evaluation",
    parameters = params,
    design = list(
      targets = targets,
      initial_train_days = initial_train_days,
      horizon_days = horizon_days,
      step_days = step_days,
      origins = origins,
      kx = 5L,
      ky = 5L,
      basis_geometry = FSSGL_BASIS_GEOMETRY_VERSION,
      preprocessing_version = SHANGHAI_PREPROCESSING_VERSION,
      leakage_rule = "max_training_day_strictly_before_minimum_test_day"
    )
  ),
  file.path(out_dir, "rolling_manifest.rds")
)

cat("\nRolling overall summary\n")
print(overall_summary)
cat("\nRolling target summary\n")
print(target_summary)
