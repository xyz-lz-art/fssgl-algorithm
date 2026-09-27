# High-dimensional comparison at p > n with common data and folds.

library(data.table)

source("R/parameters.R")
source("R/fssgl/basis_design.R")
source("R/fssgl/penalty_weights.R")
source("R/fssgl/weighted_membership_solver.R")
source("R/fssgl/solver.R")
source("R/fssgl/tuning.R")
source("R/fssgl/simulation_design.R")
source("R/baselines/functional_methods.R")
source("R/diagnostics/identifiability.R")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "results/tables/simulation/v2_formal")
processed_dir <- file.path(root, "data/processed/simulation/v2_formal")
figure_dir <- file.path(root, "results/figures/simulation/v2_formal")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

p_values <- c(50L, 60L, 100L)
n_repetitions <- as.integer(Sys.getenv("FSSGL_HIGH_DIM_REPS", "50"))
method_ids <- c("fssgl_v2", "fpca_group_scad", "structured_group_lasso")
parameters <- fssgl_parameters()
dgp_defaults <- fssgl_main_dgp_defaults()
baseline_config <- fssgl_baseline_parameters()

result_file <- file.path(table_dir, "fssgl_v2_high_dimensional_replicates.csv")
diagnostic_file <- file.path(table_dir, "fssgl_v2_identifiability_diagnostics.csv")
checkpoint_file <- file.path(processed_dir, "fssgl_v2_high_dimensional_running.rds")
completed_sources <- list()
if (file.exists(result_file)) completed_sources[[length(completed_sources) + 1L]] <- fread(result_file)
if (file.exists(checkpoint_file)) completed_sources[[length(completed_sources) + 1L]] <- readRDS(checkpoint_file)
completed <- if (length(completed_sources)) {
  unique(rbindlist(completed_sources, fill = TRUE), by = c("p", "rep", "method_id"), fromLast = TRUE)
} else data.table()
results <- list()
result_index <- 1L

for (p in p_values) {
  for (rep in seq_len(n_repetitions)) {
    seed <- 2026110000L + p * 1000L + rep
    dgp <- do.call(generate_fssgl_simulation_dgp, c(dgp_defaults, list(
      n_covariates = p,
      n_groups = p %/% 5L,
      n_active_groups = 2L,
      n_active_covariates_per_group = 2L,
      seed = seed
    )))
    folds <- make_baseline_cv_folds(dgp$dimensions$n_train, 5L, seed + 900000L)
    for (method_id in method_ids) {
      current_p <- p
      current_rep <- rep
      current_method <- method_id
      if (
        nrow(completed) > 0L &&
          nrow(completed[
            completed[["p"]] == current_p & completed[["rep"]] == current_rep &
              completed[["method_id"]] == current_method
          ]) > 0L
      ) next
      cat("High-dimensional: p=", p, ", rep=", rep, "/", n_repetitions,
          ", method=", method_id, "\n", sep = "")
      started <- proc.time()[["elapsed"]]
      fit_result <- tryCatch({
        if (method_id == "fssgl_v2") {
          fit_fssgl(
            x_coef = dgp$x_train,
            y_coef = dgp$y_train,
            structural_membership = dgp$structural_membership,
            parameters = parameters,
            verbose = FALSE
          )
        } else {
          fit_functional_baseline(
            method_id = method_id,
            x_coef = dgp$x_train,
            y_coef = dgp$y_train,
            folds = folds,
            config = baseline_config,
            structural_membership = dgp$structural_membership
          )
        }
      }, error = function(e) e)
      runtime_sec <- proc.time()[["elapsed"]] - started
      identification <- data.table(
        algorithm_version = FSSGL_ALGORITHM_VERSION,
        method_id = method_id,
        p = p,
        rep = rep,
        seed = seed,
        n_train = dgp$dimensions$n_train,
        n_test = dgp$dimensions$n_test,
        n_groups = dgp$dimensions$n_groups,
        active_covariate_count = length(dgp$active_covariates)
      )
      if (inherits(fit_result, "error")) {
        result <- cbind(identification, data.table(
          status = "error", error_message = conditionMessage(fit_result),
          runtime_sec = runtime_sec
        ))
      } else if (method_id == "fssgl_v2") {
        result <- cbind(
          identification,
          data.table(status = "ok", error_message = NA_character_, solver_converged = fit_result$fit$convergence$strict_converged),
          evaluate_fssgl_simulation_fit(fit_result, dgp, parameters$posterior_cutoff, runtime_sec)
        )
      } else {
        result <- cbind(
          identification,
          data.table(status = "ok", error_message = NA_character_, solver_converged = fit_result$solver_converged),
          evaluate_functional_baseline(fit_result, dgp, runtime_sec)
        )
      }
      results[[result_index]] <- result
      result_index <- result_index + 1L
      saveRDS(rbindlist(c(list(completed), results), fill = TRUE), checkpoint_file)
    }
  }
}

replicates <- unique(
  rbindlist(c(list(completed), results), fill = TRUE),
  by = c("p", "rep", "method_id"), fromLast = TRUE
)
setorder(replicates, p, rep, method_id)
fwrite(replicates, result_file)
if (nrow(replicates) == length(p_values) * n_repetitions * length(method_ids)) unlink(checkpoint_file)

metrics <- c("covariate_tpr", "covariate_fdr", "coefficient_relative_error", "test_coeff_rmse", "runtime_sec")
long <- melt(replicates[status == "ok"], id.vars = c("method_id", "p", "rep"), measure.vars = metrics)
summary <- long[, .(
  n = .N,
  mean = mean(value),
  sd = sd(value),
  mcse = sd(value) / sqrt(.N),
  median = median(value)
), by = .(method_id, p, metric = variable)]
fwrite(summary, file.path(table_dir, "fssgl_v2_high_dimensional_summary.csv"))

# Diagnostics use 50 common DGPs at every dimension, including the two original settings.
diagnostics <- list()
diagnostic_index <- 1L
for (p in c(10L, 20L, p_values)) {
  for (rep in seq_len(50L)) {
    seed <- if (p <= 20L) 2026101000L + p * 1000L + rep else 2026110000L + p * 1000L + rep
    n_groups <- if (p == 10L) 2L else if (p == 20L) 4L else p %/% 5L
    dgp <- do.call(generate_fssgl_simulation_dgp, c(dgp_defaults, list(
      n_covariates = p, n_groups = n_groups, n_active_groups = if (p == 10L) 1L else 2L,
      n_active_covariates_per_group = 2L, seed = seed
    )))
    diag <- functional_predictor_identifiability(dgp$x_train)
    diagnostics[[diagnostic_index]] <- as.data.table(c(list(p = p, rep = rep, seed = seed), diag))
    diagnostic_index <- diagnostic_index + 1L
  }
}
diagnostic_rows <- rbindlist(diagnostics, fill = TRUE)
fwrite(diagnostic_rows, diagnostic_file)
diagnostic_summary <- diagnostic_rows[, .(
  joint_columns = unique(n_columns),
  mean_joint_rank = mean(numerical_rank),
  mean_effective_rank = mean(entropy_effective_rank),
  mean_rank_99pct = mean(rank_99pct),
  median_nonzero_gram_condition = median(nonzero_gram_condition),
  mean_minimum_block_effective_rank = mean(minimum_block_effective_rank),
  mean_max_abs_profile_correlation = mean(maximum_absolute_profile_correlation),
  design_rank_deficient_fraction = mean(design_rank_deficient),
  severe_conditioning_fraction = mean(severely_ill_conditioned)
), by = p]
fwrite(diagnostic_summary, file.path(table_dir, "fssgl_v2_identifiability_summary.csv"))

saveRDS(new_experiment_manifest(
  experiment_id = "fssgl_v2_high_dimensional",
  parameters = list(fssgl = parameters, baselines = baseline_config),
  design = list(p = p_values, n_repetitions = n_repetitions, n_train = 40L, n_test = 24L,
                methods = method_ids, group_size = 5L, active_covariates = 4L),
  algorithm_version = FSSGL_ALGORITHM_VERSION
), file.path(processed_dir, "fssgl_v2_high_dimensional_manifest.rds"))

method_labels <- c(fssgl_v2 = "FSSGL", fpca_group_scad = "FPCA group SCAD", structured_group_lasso = "Structured group lasso")
method_colors <- c(fssgl_v2 = "#176B67", fpca_group_scad = "#8B3A62", structured_group_lasso = "#2F855A")
figure_file <- file.path(figure_dir, "fssgl_v2_high_dimensional_scaling.pdf")
pdf(figure_file, width = 8.2, height = 6.5, family = "Helvetica")
par(mfrow = c(2, 2), mar = c(4, 4.2, 2.2, 0.8), las = 1)
for (metric_name in c("covariate_tpr", "covariate_fdr", "test_coeff_rmse", "runtime_sec")) {
  values <- summary[metric == metric_name]
  ylim <- range(c(values$mean - 1.96 * values$mcse, values$mean + 1.96 * values$mcse), finite = TRUE)
  if (metric_name %in% c("covariate_tpr", "covariate_fdr")) ylim <- c(0, 1.02)
  plot(NA, xlim = range(p_values), ylim = ylim, xlab = "Number of functional predictors p",
       ylab = switch(metric_name, covariate_tpr = "Predictor TPR", covariate_fdr = "Predictor FDR",
                     test_coeff_rmse = "Test coefficient RMSE", runtime_sec = "Runtime (seconds)"),
       main = switch(metric_name, covariate_tpr = "(a) Support recovery", covariate_fdr = "(b) False discoveries",
                     test_coeff_rmse = "(c) Held-out prediction", runtime_sec = "(d) End-to-end computation"))
  grid(nx = NA, ny = NULL, col = "#E5E5E5")
  for (method_id in method_ids) {
    method_key <- method_id
    rows <- values[values[["method_id"]] == method_key][order(p)]
    lines(rows$p, rows$mean, type = "b", pch = 19, col = method_colors[[method_id]], lwd = 1.4)
    arrows(rows$p, rows$mean - 1.96 * rows$mcse, rows$p, rows$mean + 1.96 * rows$mcse,
           angle = 90, code = 3, length = 0.025, col = method_colors[[method_id]])
  }
  if (metric_name == "covariate_tpr") legend("bottomleft", method_labels, col = method_colors, lty = 1, pch = 19, bty = "n", cex = 0.72)
}
dev.off()

print(dcast(summary, method_id + p ~ metric, value.var = "mean"))
print(diagnostic_summary)
