# Paired basis-dimension sensitivity experiment for FSSGL.
#
# Each repetition is generated once in a K = 5 reference basis. The resulting
# predictor and response curves are evaluated on common dense grids and then
# projected into K = 3, 4, and 5 orthonormal B-spline bases. Thus all three fits
# within a repetition use exactly the same observed curves.

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

parameters <- fssgl_parameters()
dgp_defaults <- fssgl_main_dgp_defaults()
basis_dimensions <- 3:5
reference_dimension <- 5L
n_repetitions <- 100L
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
    posterior_slab_prob >= parameters$posterior_cutoff,
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

scenario_grid <- CJ(rep = seq_len(n_repetitions), basis_dimension = basis_dimensions)
checkpoint_file <- file.path(processed_dir, "fssgl_v2_basis_sensitivity_running.rds")
completed <- if (file.exists(checkpoint_file)) readRDS(checkpoint_file) else data.table()
results <- list()
result_index <- 1L

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

  predictor_train_curves <- array(
    NA_real_,
    dim = c(dgp$dimensions$n_covariates, dgp$dimensions$n_train, length(grid_x))
  )
  predictor_test_curves <- array(
    NA_real_,
    dim = c(dgp$dimensions$n_covariates, dgp$dimensions$n_test, length(grid_x))
  )
  for (j in seq_len(dgp$dimensions$n_covariates)) {
    predictor_train_curves[j, , ] <- dgp$x_train[j, , ] %*% t(basis_x_reference)
    predictor_test_curves[j, , ] <- dgp$x_test[j, , ] %*% t(basis_x_reference)
  }
  response_train_curves <- dgp$y_train %*% t(basis_y_reference)
  response_test_curves <- dgp$y_test %*% t(basis_y_reference)

  for (k in basis_dimensions) {
    if (
      nrow(completed) > 0L &&
        nrow(completed[rep == replicate_id & basis_dimension == k]) > 0L
    ) next
    cat("Basis sensitivity: K=", k, ", rep=", replicate_id, "/", n_repetitions, "\n", sep = "")

    current_basis_x <- basis_x[[as.character(k)]]
    current_basis_y <- basis_y[[as.character(k)]]
    x_train <- project_predictor_array(predictor_train_curves, current_basis_x)
    x_test <- project_predictor_array(predictor_test_curves, current_basis_x)
    y_train <- project_curves_l2(response_train_curves, current_basis_y)

    started <- proc.time()[["elapsed"]]
    fit_result <- tryCatch(
      fit_fssgl(
        x_coef = x_train,
        y_coef = y_train,
        structural_membership = dgp$structural_membership,
        parameters = parameters,
        verbose = FALSE
      ),
      error = function(e) e
    )
    runtime_sec <- proc.time()[["elapsed"]] - started
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
      active_group_count = length(dgp$active_groups)
    )

    if (inherits(fit_result, "error")) {
      result <- cbind(
        identification,
        data.table(status = "error", error_message = conditionMessage(fit_result))
      )
    } else {
      design_test <- build_fof_design(x_test, n_response_basis = k)
      prediction_coef <- matrix(
        drop(design_test$design %*% fit_result$fit$beta),
        nrow = dgp$dimensions$n_test,
        ncol = k
      )
      prediction_curves <- prediction_coef %*% t(current_basis_y)
      response_weights <- attr(current_basis_y, "quadrature_weights")
      curve_squared_error <- rowSums(
        sweep((prediction_curves - response_test_curves)^2, 2L, response_weights, `*`)
      )
      metrics <- cbind(
        selection_metrics(fit_result$fit, dgp),
        data.table(
          coefficient_surface_relative_error = surface_relative_error(
            fit_result$fit$beta,
            dgp,
            current_basis_x,
            current_basis_y
          ),
          test_curve_l2_rmse = sqrt(mean(curve_squared_error)),
          runtime_sec = runtime_sec,
          final_iter = fit_result$fit$convergence$final_iter,
          strict_converged = fit_result$fit$convergence$strict_converged,
          relaxed_converged = fit_result$fit$convergence$relaxed_converged,
          final_beta_change = fit_result$fit$convergence$final_beta_change,
          objective_tail_rel_change = fit_result$fit$convergence$objective_tail_rel_change
        )
      )
      result <- cbind(
        identification,
        data.table(status = "ok", error_message = NA_character_),
        metrics
      )
    }
    results[[result_index]] <- result
    result_index <- result_index + 1L
    saveRDS(rbindlist(c(list(completed), results), fill = TRUE), checkpoint_file)
  }
}

replicates <- rbindlist(c(list(completed), results), fill = TRUE)
setorder(replicates, basis_dimension, rep)
metric_columns <- c(
  "covariate_tpr", "covariate_fpr", "covariate_fdr", "selected_covariate_count",
  "group_tpr", "group_fpr", "group_fdr", "selected_group_count",
  "coefficient_surface_relative_error", "test_curve_l2_rmse", "runtime_sec",
  "final_iter", "final_beta_change", "objective_tail_rel_change"
)
replicates[, (metric_columns) := lapply(.SD, as.numeric), .SDcols = metric_columns]
summary_long <- melt(
  replicates[status == "ok"],
  id.vars = "basis_dimension",
  measure.vars = metric_columns,
  variable.name = "metric",
  value.name = "value"
)[, .(
  n = .N,
  mean = mean(value, na.rm = TRUE),
  sd = sd(value, na.rm = TRUE),
  mcse = sd(value, na.rm = TRUE) / sqrt(.N),
  median = median(value, na.rm = TRUE),
  q05 = quantile(value, 0.05, na.rm = TRUE),
  q95 = quantile(value, 0.95, na.rm = TRUE)
), by = .(basis_dimension, metric)]
summary_wide <- dcast(
  summary_long[, mean_sd := sprintf("%.4f (%.4f)", mean, sd)],
  basis_dimension ~ metric,
  value.var = "mean_sd"
)
convergence <- replicates[, .(
  n_attempted = .N,
  n_ok = sum(status == "ok"),
  strict_rate = mean(strict_converged[status == "ok"], na.rm = TRUE),
  relaxed_rate = mean(relaxed_converged[status == "ok"], na.rm = TRUE),
  error_count = sum(status != "ok")
), by = basis_dimension]

expected_rows <- nrow(scenario_grid)
if (nrow(replicates) != expected_rows || anyDuplicated(replicates[, .(rep, basis_dimension)])) {
  stop("Basis-sensitivity output is incomplete or contains duplicate scenarios.")
}
if (any(replicates$status != "ok")) {
  stop("At least one basis-sensitivity fit returned an error.")
}

fwrite(replicates, file.path(table_dir, "fssgl_v2_basis_sensitivity_replicates.csv"))
fwrite(summary_long, file.path(table_dir, "fssgl_v2_basis_sensitivity_summary_long.csv"))
fwrite(summary_wide, file.path(table_dir, "fssgl_v2_basis_sensitivity_summary_mean_sd.csv"))
fwrite(convergence, file.path(table_dir, "fssgl_v2_basis_sensitivity_convergence.csv"))

pdf(file.path(figure_dir, "fssgl_v2_basis_sensitivity.pdf"), width = 8.8, height = 6.6)
old_par <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), mar = c(4, 4.2, 2.5, 1), las = 1)
plot_metric <- function(metric_name, ylab, ylim = NULL) {
  values <- summary_long[summary_long$metric == metric_name]
  if (is.null(ylim)) {
    bounds <- range(
      values$mean - 1.96 * values$mcse,
      values$mean + 1.96 * values$mcse
    )
    padding <- max(diff(bounds) * 0.15, 0.01)
    ylim <- bounds + c(-padding, padding)
  }
  plot(
    values$basis_dimension,
    values$mean,
    type = "b",
    pch = 19,
    lwd = 1.5,
    col = "#1B6CA8",
    xlab = "Basis dimension (KX = KY)",
    ylab = ylab,
    xaxt = "n",
    ylim = ylim
  )
  axis(1, at = basis_dimensions)
  arrows(
    values$basis_dimension,
    values$mean - 1.96 * values$mcse,
    values$basis_dimension,
    values$mean + 1.96 * values$mcse,
    angle = 90,
    code = 3,
    length = 0.05,
    col = "#1B6CA8"
  )
  grid(col = "#DDDDDD")
}
plot_metric("covariate_tpr", "Predictor TPR", ylim = c(0, 1.02))
plot_metric("covariate_fdr", "Predictor FDR", ylim = c(0, max(summary_long[metric == "covariate_fdr", mean + 1.96 * mcse]) * 1.15))
plot_metric("coefficient_surface_relative_error", "Surface relative error")
plot_metric("test_curve_l2_rmse", "Test curve L2 RMSE")
par(old_par)
dev.off()

saveRDS(
  new_experiment_manifest(
    experiment_id = "fssgl_v2_basis_dimension_sensitivity",
    parameters = parameters,
    design = list(
      reference_dimension = reference_dimension,
      fitted_dimensions = basis_dimensions,
      common_curve_grids = list(x = grid_x, y = grid_y),
      n_repetitions = n_repetitions,
      p = 20L,
      dgp_defaults = dgp_defaults,
      paired_curve_generation = TRUE
    ),
    algorithm_version = FSSGL_ALGORITHM_VERSION
  ),
  file.path(processed_dir, "fssgl_v2_basis_sensitivity_manifest.rds")
)
if (file.exists(checkpoint_file)) unlink(checkpoint_file)

print(summary_wide[, .(
  basis_dimension,
  covariate_tpr,
  covariate_fdr,
  coefficient_surface_relative_error,
  test_curve_l2_rmse,
  runtime_sec
)])
print(convergence)
