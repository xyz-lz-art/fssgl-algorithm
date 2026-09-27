# Training-only penalty recalibration for each fitted basis dimension.
#
# Within every repetition, the same observed curves are projected to K = 3, 4,
# and 5. For each K, five-fold curve-level cross-validation selects a common
# multiplier for the covariate and group penalty scales. The untouched test set
# is used only after refitting on all 40 training curves.

library(data.table)

source("R/parameters.R")
source("R/fssgl/basis_design.R")
source("R/fssgl/penalty_weights.R")
source("R/fssgl/weighted_membership_solver.R")
source("R/fssgl/solver.R")
source("R/fssgl/simulation_design.R")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "results/tables/simulation/v2_formal")
figure_dir <- file.path(root, "results/figures/simulation/v2_formal")
processed_dir <- file.path(root, "data/processed/simulation/v2_formal")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

base_parameters <- fssgl_parameters()
dgp_defaults <- fssgl_main_dgp_defaults()
basis_dimensions <- 3:5
reference_dimension <- 5L
n_repetitions <- 100L
n_folds <- 5L
penalty_multipliers <- c(0.125, 0.25, 0.5, 1, 2, 4)
grid_x <- seq(0, 1, length.out = 51L)
grid_y <- seq(0, 1, length.out = 51L)

basis_x_reference <- make_bspline_basis(grid_x, df = reference_dimension)
basis_y_reference <- make_bspline_basis(grid_y, df = reference_dimension)
basis_x <- setNames(lapply(basis_dimensions, function(k) {
  make_bspline_basis(grid_x, df = k)
}), basis_dimensions)
basis_y <- setNames(lapply(basis_dimensions, function(k) {
  make_bspline_basis(grid_y, df = k)
}), basis_dimensions)

make_observed_curves <- function(dgp) {
  x_train <- array(
    NA_real_,
    dim = c(dgp$dimensions$n_covariates, dgp$dimensions$n_train, length(grid_x))
  )
  x_test <- array(
    NA_real_,
    dim = c(dgp$dimensions$n_covariates, dgp$dimensions$n_test, length(grid_x))
  )
  for (j in seq_len(dgp$dimensions$n_covariates)) {
    x_train[j, , ] <- dgp$x_train[j, , ] %*% t(basis_x_reference)
    x_test[j, , ] <- dgp$x_test[j, , ] %*% t(basis_x_reference)
  }
  list(
    x_train = x_train,
    x_test = x_test,
    y_train = dgp$y_train %*% t(basis_y_reference),
    y_test = dgp$y_test %*% t(basis_y_reference)
  )
}

surface_relative_error <- function(beta_hat, dgp, current_basis_x, current_basis_y) {
  block_size <- ncol(current_basis_x) * ncol(current_basis_y)
  squared_error <- 0
  squared_truth <- 0
  weights <- outer(
    attr(current_basis_x, "quadrature_weights"),
    attr(current_basis_y, "quadrature_weights")
  )
  for (j in seq_len(dgp$dimensions$n_covariates)) {
    columns <- ((j - 1L) * block_size + 1L):(j * block_size)
    coefficient_hat <- matrix(
      beta_hat[columns],
      nrow = ncol(current_basis_x),
      ncol = ncol(current_basis_y)
    )
    surface_hat <- evaluate_coefficient_surface(
      coefficient_hat,
      current_basis_x,
      current_basis_y
    )
    surface_true <- evaluate_coefficient_surface(
      dgp$surfaces[[j]],
      basis_x_reference,
      basis_y_reference
    )
    squared_error <- squared_error + sum((surface_hat - surface_true)^2 * weights)
    squared_truth <- squared_truth + sum(surface_true^2 * weights)
  }
  sqrt(squared_error) / (sqrt(squared_truth) + 1e-8)
}

selection_metrics <- function(fit, dgp) {
  selected_covariates <- fit$covariate_posterior[
    posterior_slab_prob >= base_parameters$posterior_cutoff,
    covariate_id
  ]
  selected_groups <- unique(dgp$structural_membership[
    covariate_id %in% selected_covariates,
    group_id
  ])
  inactive_covariates <- setdiff(
    seq_len(dgp$dimensions$n_covariates),
    dgp$active_covariates
  )
  inactive_groups <- setdiff(dgp$truth_groups$group_id, dgp$active_groups)
  data.table(
    covariate_tpr = mean(dgp$active_covariates %in% selected_covariates),
    covariate_fpr = mean(inactive_covariates %in% selected_covariates),
    covariate_fdr = if (length(selected_covariates) == 0L) 0 else {
      mean(!selected_covariates %in% dgp$active_covariates)
    },
    selected_covariate_count = length(selected_covariates),
    group_tpr = mean(dgp$active_groups %in% selected_groups),
    group_fpr = mean(inactive_groups %in% selected_groups),
    group_fdr = if (length(selected_groups) == 0L) 0 else {
      mean(!selected_groups %in% dgp$active_groups)
    },
    selected_group_count = length(selected_groups)
  )
}

curve_l2_squared_error <- function(prediction_coef, response_curves, response_basis) {
  prediction_curves <- prediction_coef %*% t(response_basis)
  weights <- attr(response_basis, "quadrature_weights")
  rowSums(sweep((prediction_curves - response_curves)^2, 2L, weights, `*`))
}

checkpoint_file <- file.path(processed_dir, "fssgl_v2_basis_recalibration_running.rds")
empty_state <- list(replicates = data.table(), tuning = data.table())
state <- if (file.exists(checkpoint_file)) readRDS(checkpoint_file) else empty_state

for (replicate_id in seq_len(n_repetitions)) {
  seed <- 2026107000L + replicate_id
  dgp <- do.call(generate_fssgl_simulation_dgp, modifyList(dgp_defaults, list(
    n_covariates = 20L,
    kx = reference_dimension,
    ky = reference_dimension,
    n_groups = 4L,
    n_active_groups = 2L,
    n_active_covariates_per_group = 2L,
    seed = seed
  )))
  observed <- make_observed_curves(dgp)
  set.seed(seed + 100000L)
  fold_id <- sample(rep(seq_len(n_folds), length.out = dgp$dimensions$n_train))

  for (k in basis_dimensions) {
    already_done <- nrow(state$replicates) > 0L && nrow(state$replicates[
      rep == replicate_id & basis_dimension == k
    ]) > 0L
    if (already_done) next
    cat("Basis recalibration: K=", k, ", rep=", replicate_id, "/", n_repetitions, "\n", sep = "")

    current_basis_x <- basis_x[[as.character(k)]]
    current_basis_y <- basis_y[[as.character(k)]]
    x_train <- project_predictor_array(observed$x_train, current_basis_x)
    x_test <- project_predictor_array(observed$x_test, current_basis_x)
    y_train <- project_curves_l2(observed$y_train, current_basis_y)

    tuning_rows <- vector("list", length(penalty_multipliers))
    tuning_started <- proc.time()[["elapsed"]]
    for (candidate_id in seq_along(penalty_multipliers)) {
      multiplier <- penalty_multipliers[candidate_id]
      candidate_parameters <- base_parameters
      candidate_parameters$covariate_penalty_scale <-
        base_parameters$covariate_penalty_scale * multiplier
      candidate_parameters$group_penalty_scale <-
        base_parameters$group_penalty_scale * multiplier
      squared_error <- numeric(dgp$dimensions$n_train)
      candidate_errors <- character()

      for (fold in seq_len(n_folds)) {
        validation_index <- which(fold_id == fold)
        fitting_index <- which(fold_id != fold)
        fold_fit <- tryCatch(
          fit_fssgl(
            x_coef = x_train[, fitting_index, , drop = FALSE],
            y_coef = y_train[fitting_index, , drop = FALSE],
            structural_membership = dgp$structural_membership,
            parameters = candidate_parameters,
            verbose = FALSE
          ),
          error = function(e) e
        )
        if (inherits(fold_fit, "error")) {
          candidate_errors <- c(candidate_errors, conditionMessage(fold_fit))
          squared_error[validation_index] <- Inf
          next
        }
        validation_design <- build_fof_design(
          x_train[, validation_index, , drop = FALSE],
          n_response_basis = k
        )
        prediction_coef <- matrix(
          drop(validation_design$design %*% fold_fit$fit$beta),
          nrow = length(validation_index),
          ncol = k
        )
        squared_error[validation_index] <- curve_l2_squared_error(
          prediction_coef,
          observed$y_train[validation_index, , drop = FALSE],
          current_basis_y
        )
      }

      tuning_rows[[candidate_id]] <- data.table(
        rep = replicate_id,
        seed = seed,
        basis_dimension = k,
        multiplier = multiplier,
        covariate_penalty_scale = candidate_parameters$covariate_penalty_scale,
        group_penalty_scale = candidate_parameters$group_penalty_scale,
        cv_curve_l2_rmse = sqrt(mean(squared_error)),
        cv_error_count = length(candidate_errors),
        cv_error_message = if (length(candidate_errors)) {
          paste(unique(candidate_errors), collapse = " | ")
        } else {
          NA_character_
        }
      )
    }
    tuning_elapsed <- proc.time()[["elapsed"]] - tuning_started
    tuning_table <- rbindlist(tuning_rows)
    finite_candidates <- tuning_table[is.finite(cv_curve_l2_rmse)]
    if (!nrow(finite_candidates)) stop("All calibration candidates failed.")
    setorder(finite_candidates, cv_curve_l2_rmse, -multiplier)
    chosen <- finite_candidates[1L]

    selected_parameters <- base_parameters
    selected_parameters$covariate_penalty_scale <- chosen$covariate_penalty_scale
    selected_parameters$group_penalty_scale <- chosen$group_penalty_scale
    final_started <- proc.time()[["elapsed"]]
    rescue_refit <- FALSE
    fit_result <- tryCatch(
      fit_fssgl(
        x_coef = x_train,
        y_coef = y_train,
        structural_membership = dgp$structural_membership,
        parameters = selected_parameters,
        verbose = FALSE
      ),
      error = function(e) e
    )
    if (!inherits(fit_result, "error") &&
        !isTRUE(fit_result$fit$convergence$strict_converged)) {
      rescue_parameters <- selected_parameters
      rescue_parameters$max_outer_iter <- 2000L
      fit_result <- tryCatch(
        fit_fssgl(
          x_coef = x_train,
          y_coef = y_train,
          structural_membership = dgp$structural_membership,
          parameters = rescue_parameters,
          verbose = FALSE
        ),
        error = function(e) e
      )
      rescue_refit <- TRUE
    }
    final_runtime_sec <- proc.time()[["elapsed"]] - final_started

    identification <- data.table(
      algorithm_version = FSSGL_ALGORITHM_VERSION,
      rep = replicate_id,
      seed = seed,
      basis_dimension = k,
      kx = k,
      ky = k,
      reference_dimension = reference_dimension,
      n_train = dgp$dimensions$n_train,
      n_test = dgp$dimensions$n_test,
      p = dgp$dimensions$n_covariates,
      n_groups = dgp$dimensions$n_groups,
      active_covariate_count = length(dgp$active_covariates),
      active_group_count = length(dgp$active_groups),
      n_folds = n_folds,
      selected_multiplier = chosen$multiplier,
      selected_covariate_penalty_scale = chosen$covariate_penalty_scale,
      selected_group_penalty_scale = chosen$group_penalty_scale,
      selected_cv_curve_l2_rmse = chosen$cv_curve_l2_rmse,
      tuning_runtime_sec = tuning_elapsed,
      final_runtime_sec = final_runtime_sec,
      rescue_refit = rescue_refit
    )

    if (inherits(fit_result, "error")) {
      result <- cbind(
        identification,
        data.table(status = "error", error_message = conditionMessage(fit_result))
      )
    } else {
      test_design <- build_fof_design(x_test, n_response_basis = k)
      prediction_coef <- matrix(
        drop(test_design$design %*% fit_result$fit$beta),
        nrow = dgp$dimensions$n_test,
        ncol = k
      )
      result <- cbind(
        identification,
        data.table(status = "ok", error_message = NA_character_),
        selection_metrics(fit_result$fit, dgp),
        data.table(
          coefficient_surface_relative_error = surface_relative_error(
            fit_result$fit$beta,
            dgp,
            current_basis_x,
            current_basis_y
          ),
          test_curve_l2_rmse = sqrt(mean(curve_l2_squared_error(
            prediction_coef,
            observed$y_test,
            current_basis_y
          ))),
          final_iter = fit_result$fit$convergence$final_iter,
          strict_converged = fit_result$fit$convergence$strict_converged,
          relaxed_converged = fit_result$fit$convergence$relaxed_converged,
          final_beta_change = fit_result$fit$convergence$final_beta_change,
          objective_tail_rel_change = fit_result$fit$convergence$objective_tail_rel_change
        )
      )
    }

    tuning_table[, selected := multiplier == chosen$multiplier]
    state$tuning <- rbind(state$tuning, tuning_table, fill = TRUE)
    state$replicates <- rbind(state$replicates, result, fill = TRUE)
    saveRDS(state, checkpoint_file)
  }
}

replicates <- state$replicates
tuning <- state$tuning
if (!"rescue_refit" %in% names(replicates)) replicates[, rescue_refit := FALSE]
replicates[is.na(rescue_refit), rescue_refit := FALSE]
setorder(replicates, basis_dimension, rep)
setorder(tuning, basis_dimension, rep, multiplier)
expected_rows <- n_repetitions * length(basis_dimensions)
if (nrow(replicates) != expected_rows || anyDuplicated(replicates[, .(rep, basis_dimension)])) {
  stop("Recalibration output is incomplete or contains duplicate scenarios.")
}
if (any(replicates$status != "ok")) stop("At least one final recalibrated fit failed.")
if (nrow(tuning) != expected_rows * length(penalty_multipliers)) {
  stop("Calibration-path output is incomplete.")
}
if (any(tuning$cv_error_count != 0L) || any(!is.finite(tuning$cv_curve_l2_rmse))) {
  stop("At least one cross-validation fit failed.")
}

metric_columns <- c(
  "covariate_tpr", "covariate_fpr", "covariate_fdr", "selected_covariate_count",
  "group_tpr", "group_fpr", "group_fdr", "selected_group_count",
  "coefficient_surface_relative_error", "test_curve_l2_rmse",
  "selected_multiplier", "selected_cv_curve_l2_rmse",
  "tuning_runtime_sec", "final_runtime_sec",
  "final_iter", "final_beta_change", "objective_tail_rel_change"
)
replicates[, (metric_columns) := lapply(.SD, as.numeric), .SDcols = metric_columns]
summary_long <- melt(
  replicates,
  id.vars = "basis_dimension",
  measure.vars = metric_columns,
  variable.name = "metric",
  value.name = "value"
)[, .(
  n = .N,
  mean = mean(value),
  sd = sd(value),
  mcse = sd(value) / sqrt(.N),
  median = median(value),
  q05 = quantile(value, 0.05),
  q95 = quantile(value, 0.95)
), by = .(basis_dimension, metric)]
summary_wide <- dcast(
  summary_long[, mean_sd := sprintf("%.4f (%.4f)", mean, sd)],
  basis_dimension ~ metric,
  value.var = "mean_sd"
)
multiplier_frequency <- replicates[, .(
  selection_count = .N,
  selection_frequency = .N / n_repetitions
), by = .(basis_dimension, selected_multiplier)]
setorder(multiplier_frequency, basis_dimension, selected_multiplier)
convergence <- replicates[, .(
  n_attempted = .N,
  n_ok = sum(status == "ok"),
  strict_rate = mean(strict_converged),
  relaxed_rate = mean(relaxed_converged),
  error_count = sum(status != "ok")
), by = basis_dimension]

fixed_file <- file.path(table_dir, "fssgl_v2_basis_sensitivity_replicates.csv")
if (!file.exists(fixed_file)) stop("Run 07_basis_dimension_sensitivity.R first.")
fixed <- fread(fixed_file)
comparison_metrics <- c(
  "covariate_tpr", "covariate_fpr", "covariate_fdr", "group_tpr",
  "coefficient_surface_relative_error", "test_curve_l2_rmse"
)
comparison <- rbindlist(list(
  fixed[, c(list(calibration = "Fixed"), lapply(.SD, as.numeric)),
        .SDcols = c("basis_dimension", comparison_metrics)],
  replicates[, c(list(calibration = "Recalibrated"), lapply(.SD, as.numeric)),
             .SDcols = c("basis_dimension", comparison_metrics)]
), use.names = TRUE)
comparison_long <- melt(
  comparison,
  id.vars = c("calibration", "basis_dimension"),
  measure.vars = comparison_metrics,
  variable.name = "metric",
  value.name = "value"
)[, .(
  n = .N,
  mean = mean(value),
  sd = sd(value),
  mcse = sd(value) / sqrt(.N)
), by = .(calibration, basis_dimension, metric)]

fwrite(replicates, file.path(table_dir, "fssgl_v2_basis_recalibration_replicates.csv"))
fwrite(tuning, file.path(table_dir, "fssgl_v2_basis_recalibration_cv_path.csv"))
fwrite(summary_long, file.path(table_dir, "fssgl_v2_basis_recalibration_summary_long.csv"))
fwrite(summary_wide, file.path(table_dir, "fssgl_v2_basis_recalibration_summary_mean_sd.csv"))
fwrite(multiplier_frequency, file.path(table_dir, "fssgl_v2_basis_recalibration_multiplier_frequency.csv"))
fwrite(convergence, file.path(table_dir, "fssgl_v2_basis_recalibration_convergence.csv"))
fwrite(comparison_long, file.path(table_dir, "fssgl_v2_basis_recalibration_comparison.csv"))

pdf(file.path(figure_dir, "fssgl_v2_basis_recalibration.pdf"), width = 8.8, height = 6.6)
old_par <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), mar = c(4, 4.2, 2.5, 1), las = 1)
colors <- c(Fixed = "#999999", Recalibrated = "#1B6CA8")
plot_comparison <- function(metric_name, ylab, ylim = NULL) {
  values <- comparison_long[metric == metric_name]
  if (is.null(ylim)) {
    bounds <- range(values$mean - 1.96 * values$mcse, values$mean + 1.96 * values$mcse)
    padding <- max(diff(bounds) * 0.12, 0.01)
    ylim <- bounds + c(-padding, padding)
  }
  plot(NA, xlim = range(basis_dimensions), ylim = ylim, xaxt = "n",
       xlab = "Basis dimension (KX = KY)", ylab = ylab)
  axis(1, at = basis_dimensions)
  grid(col = "#DDDDDD")
  for (label in names(colors)) {
    current <- values[calibration == label][order(basis_dimension)]
    lines(current$basis_dimension, current$mean, type = "b", pch = 19,
          lwd = 1.5, col = colors[label])
    arrows(current$basis_dimension, current$mean - 1.96 * current$mcse,
           current$basis_dimension, current$mean + 1.96 * current$mcse,
           angle = 90, code = 3, length = 0.05, col = colors[label])
  }
  legend("topleft", legend = names(colors), col = colors, lty = 1, pch = 19,
         bty = "n", cex = 0.85)
}
plot_comparison("covariate_tpr", "Predictor TPR", ylim = c(0, 1.02))
plot_comparison("covariate_fdr", "Predictor FDR", ylim = c(0, 0.35))
plot_comparison("coefficient_surface_relative_error", "Surface relative error")
plot_comparison("test_curve_l2_rmse", "Test curve L2 RMSE")
par(old_par)
dev.off()

saveRDS(
  new_experiment_manifest(
    experiment_id = "fssgl_v2_basis_dimension_recalibration",
    parameters = list(
      base = base_parameters,
      penalty_multipliers = penalty_multipliers,
      n_folds = n_folds,
      selection_rule = "minimum pooled out-of-fold curve L2 RMSE"
    ),
    design = list(
      reference_dimension = reference_dimension,
      fitted_dimensions = basis_dimensions,
      common_curve_grids = list(x = grid_x, y = grid_y),
      n_repetitions = n_repetitions,
      p = 20L,
      dgp_defaults = dgp_defaults,
      paired_curve_generation = TRUE,
      common_folds_across_dimensions = TRUE,
      test_set_used_for_calibration = FALSE
    ),
    algorithm_version = FSSGL_ALGORITHM_VERSION
  ),
  file.path(processed_dir, "fssgl_v2_basis_recalibration_manifest.rds")
)
if (file.exists(checkpoint_file)) unlink(checkpoint_file)

print(summary_wide[, .(
  basis_dimension,
  covariate_tpr,
  covariate_fdr,
  coefficient_surface_relative_error,
  test_curve_l2_rmse,
  selected_multiplier
)])
print(multiplier_frequency)
print(convergence)
