# Frozen independent-pilot calibration used to reproduce the current v2 table.
# New submission comparisons use 11_fair_fssgl_comparison.R, which shares folds
# and curve-level loss with the baselines and rejects boundary/nonconverged fits.

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
figure_dir <- file.path(root, "results/figures/simulation/v2_formal")
p_values <- c(50L, 60L, 100L)
n_repetitions <- 50L
n_pilot <- 10L
multiplier_grid <- c(0.001953125, 0.00390625, 0.0078125, 0.015625, 0.03125, 0.0625)
dgp_defaults <- fssgl_main_dgp_defaults()
base_parameters <- fssgl_parameters()
experiment_version <- "independent_pilot_holdout_v1"

result_file <- file.path(table_dir, "fssgl_v2_high_dimensional_recalibrated_replicates.csv")
pilot_file <- file.path(table_dir, "fssgl_v2_high_dimensional_pilot_calibration.csv")
checkpoint_file <- file.path(processed_dir, "fssgl_v2_high_dimensional_recalibration_running.rds")
state <- if (file.exists(checkpoint_file)) readRDS(checkpoint_file) else NULL
if (
  is.null(state) && file.exists(result_file) && file.exists(pilot_file) &&
    file.exists(file.path(table_dir, "fssgl_v2_high_dimensional_selected_multipliers.csv"))
) {
  state <- list(
    experiment_version = experiment_version,
    pilot = fread(pilot_file),
    selected = fread(file.path(table_dir, "fssgl_v2_high_dimensional_selected_multipliers.csv")),
    results = fread(result_file)
  )
}
if (is.null(state) || !identical(state$experiment_version, experiment_version)) {
  state <- list(
    experiment_version = experiment_version,
    pilot = data.table(),
    selected = data.table(),
    results = data.table()
  )
  saveRDS(state, checkpoint_file)
}

scaled_parameters <- function(multiplier) {
  parameters <- base_parameters
  parameters$covariate_penalty_scale <- base_parameters$covariate_penalty_scale * multiplier
  parameters$group_penalty_scale <- base_parameters$group_penalty_scale * multiplier
  parameters
}

# Ten independent pilot samples per dimension.  Days 1:32 are used for fitting
# and 33:40 for validation; none of these seeds occur in the reported experiment.
if (!nrow(state$selected)) {
  pilot_rows <- list()
  pilot_index <- 1L
  for (p in p_values) {
    for (pilot_rep in seq_len(n_pilot)) {
      seed <- 2026990000L + p * 100L + pilot_rep
      dgp <- do.call(generate_fssgl_simulation_dgp, c(dgp_defaults, list(
        n_covariates = p, n_groups = p %/% 5L, n_active_groups = 2L,
        n_active_covariates_per_group = 2L, seed = seed
      )))
      training <- 1:32
      validation <- 33:40
      validation_design <- build_fof_design(
        dgp$x_train[, validation, , drop = FALSE],
        n_response_basis = dgp$dimensions$ky
      )$design
      for (multiplier in multiplier_grid) {
        cat("Pilot calibration: p=", p, ", pilot=", pilot_rep, "/", n_pilot,
            ", multiplier=", multiplier, "\n", sep = "")
        started <- proc.time()[["elapsed"]]
        fit <- tryCatch(fit_fssgl(
          x_coef = dgp$x_train[, training, , drop = FALSE],
          y_coef = dgp$y_train[training, , drop = FALSE],
          structural_membership = dgp$structural_membership,
          parameters = scaled_parameters(multiplier),
          verbose = FALSE
        ), error = function(e) e)
        runtime_sec <- proc.time()[["elapsed"]] - started
        if (inherits(fit, "error")) {
          rmse <- Inf
          status <- "error"
        } else {
          prediction <- matrix(
            drop(validation_design %*% fit$fit$beta),
            nrow = length(validation), ncol = dgp$dimensions$ky
          )
          rmse <- sqrt(mean((prediction - dgp$y_train[validation, , drop = FALSE])^2))
          status <- "ok"
        }
        pilot_rows[[pilot_index]] <- data.table(
          p = p, pilot_rep = pilot_rep, seed = seed, multiplier = multiplier,
          validation_rmse = rmse, status = status, runtime_sec = runtime_sec
        )
        pilot_index <- pilot_index + 1L
      }
    }
  }
  state$pilot <- rbindlist(pilot_rows)
  state$selected <- state$pilot[status == "ok", .(
    mean_validation_rmse = mean(validation_rmse),
    sd_validation_rmse = sd(validation_rmse)
  ), by = .(p, multiplier)][order(p, mean_validation_rmse)][, .SD[1L], by = p]
  saveRDS(state, checkpoint_file)
}

for (p in p_values) {
  multiplier <- state$selected[["multiplier"]][state$selected[["p"]] == p]
  for (rep in seq_len(n_repetitions)) {
    if (nrow(state$results) && any(state$results[["p"]] == p & state$results[["rep"]] == rep)) next
    seed <- 2026110000L + p * 1000L + rep
    dgp <- do.call(generate_fssgl_simulation_dgp, c(dgp_defaults, list(
      n_covariates = p, n_groups = p %/% 5L, n_active_groups = 2L,
      n_active_covariates_per_group = 2L, seed = seed
    )))
    cat("Recalibrated FSSGL: p=", p, ", rep=", rep, "/", n_repetitions,
        ", multiplier=", multiplier, "\n", sep = "")
    started <- proc.time()[["elapsed"]]
    fit <- tryCatch(fit_fssgl(
      x_coef = dgp$x_train, y_coef = dgp$y_train,
      structural_membership = dgp$structural_membership,
      parameters = scaled_parameters(multiplier), verbose = FALSE
    ), error = function(e) e)
    runtime_sec <- proc.time()[["elapsed"]] - started
    if (inherits(fit, "error")) {
      result <- data.table(
        algorithm_version = FSSGL_ALGORITHM_VERSION, method_id = "fssgl_v2_recalibrated",
        p = p, rep = rep, seed = seed, status = "error",
        error_message = conditionMessage(fit), selected_multiplier = multiplier,
        runtime_sec = runtime_sec
      )
    } else {
      result <- cbind(data.table(
        algorithm_version = FSSGL_ALGORITHM_VERSION, method_id = "fssgl_v2_recalibrated",
        p = p, rep = rep, seed = seed, status = "ok", error_message = NA_character_,
        selected_multiplier = multiplier
      ), evaluate_fssgl_simulation_fit(
        fit, dgp, posterior_cutoff = base_parameters$posterior_cutoff,
        runtime_sec = runtime_sec
      ))
    }
    state$results <- rbindlist(list(state$results, result), fill = TRUE)
    saveRDS(state, checkpoint_file)
  }
}

setorder(state$results, p, rep)
fwrite(state$results, result_file)
fwrite(state$pilot, pilot_file)
fwrite(state$selected, file.path(table_dir, "fssgl_v2_high_dimensional_selected_multipliers.csv"))

other <- fread(file.path(table_dir, "fssgl_v2_high_dimensional_replicates.csv"))
pilot_cost <- state$pilot[, .(
  amortized_pilot_runtime_sec = sum(runtime_sec) / n_repetitions
), by = p]
state$results <- merge(state$results, pilot_cost, by = "p", all.x = TRUE)
state$results[, comparison_runtime_sec := runtime_sec + amortized_pilot_runtime_sec]
other[, comparison_runtime_sec := runtime_sec]
combined <- rbindlist(list(
  state$results,
  other[method_id %in% c("fpca_group_scad", "structured_group_lasso")]
), fill = TRUE)
metrics <- c("covariate_tpr", "covariate_fdr", "coefficient_relative_error", "test_coeff_rmse", "comparison_runtime_sec")
long <- melt(combined[status == "ok"], id.vars = c("method_id", "p", "rep"), measure.vars = metrics)
long[variable == "comparison_runtime_sec", variable := "runtime_sec"]
summary <- long[, .(
  n = .N, mean = mean(value), sd = sd(value), mcse = sd(value) / sqrt(.N), median = median(value)
), by = .(method_id, p, metric = variable)]
fwrite(summary, file.path(table_dir, "fssgl_v2_high_dimensional_comparison_summary.csv"))

saveRDS(new_experiment_manifest(
  experiment_id = "fssgl_v2_high_dimensional_recalibration",
  parameters = list(base = base_parameters, multiplier_grid = multiplier_grid),
  design = list(
    version = experiment_version, p = p_values, n_repetitions = n_repetitions,
    pilot_repetitions = n_pilot, pilot_split = "32 training / 8 validation",
    pilot_seed_rule = "2026990000 + p * 100 + pilot_rep",
    evaluation_seed_rule = "2026110000 + p * 1000 + rep"
  ),
  algorithm_version = FSSGL_ALGORITHM_VERSION
), file.path(processed_dir, "fssgl_v2_high_dimensional_recalibration_manifest.rds"))
if (nrow(state$results) == length(p_values) * n_repetitions) unlink(checkpoint_file)

method_ids <- c("fssgl_v2_recalibrated", "fpca_group_scad", "structured_group_lasso")
method_labels <- c(fssgl_v2_recalibrated = "FSSGL (recalibrated)", fpca_group_scad = "FPCA group SCAD", structured_group_lasso = "Structured group lasso")
method_colors <- c(fssgl_v2_recalibrated = "#176B67", fpca_group_scad = "#8B3A62", structured_group_lasso = "#2F855A")
pdf(file.path(figure_dir, "fssgl_v2_high_dimensional_scaling.pdf"), width = 8.2, height = 6.5, family = "Helvetica")
par(mfrow = c(2, 2), mar = c(4, 4.2, 2.2, 0.8), las = 1)
for (metric_name in c("covariate_tpr", "covariate_fdr", "test_coeff_rmse", "runtime_sec")) {
  values <- summary[metric == metric_name]
  ylim <- range(c(values$mean - 1.96 * values$mcse, values$mean + 1.96 * values$mcse), finite = TRUE)
  if (metric_name %in% c("covariate_tpr", "covariate_fdr")) ylim <- c(0, 1.02)
  plot(NA, xlim = range(p_values), ylim = ylim, xlab = "Number of functional predictors p",
       ylab = switch(metric_name, covariate_tpr = "Predictor TPR", covariate_fdr = "Predictor FDR",
                     test_coeff_rmse = "Test coefficient RMSE", runtime_sec = "Runtime including tuning (seconds)"),
       main = switch(metric_name, covariate_tpr = "(a) Support recovery", covariate_fdr = "(b) False discoveries",
                     test_coeff_rmse = "(c) Held-out prediction", runtime_sec = "(d) Computation including tuning"))
  grid(nx = NA, ny = NULL, col = "#E5E5E5")
  for (method_key in method_ids) {
    rows <- values[values[["method_id"]] == method_key][order(p)]
    lines(rows$p, rows$mean, type = "b", pch = 19, col = method_colors[[method_key]], lwd = 1.4)
    lower <- rows$mean - 1.96 * rows$mcse
    upper <- rows$mean + 1.96 * rows$mcse
    nonzero <- is.finite(lower) & is.finite(upper) & upper - lower > 1e-12
    if (any(nonzero)) arrows(rows$p[nonzero], lower[nonzero], rows$p[nonzero], upper[nonzero],
                             angle = 90, code = 3, length = 0.025, col = method_colors[[method_key]])
  }
  if (metric_name == "covariate_tpr") legend("bottomleft", method_labels, col = method_colors,
                                               lty = 1, pch = 19, bty = "n", cex = 0.72)
}
dev.off()

print(state$selected)
print(dcast(summary, method_id + p ~ metric, value.var = "mean"))
