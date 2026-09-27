# Formal structural-information ablation on 100 common-seed data sets.

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
dgp_defaults <- modifyList(fssgl_main_dgp_defaults(), list(
  n_covariates = 20L,
  n_groups = 4L,
  n_active_groups = 2L,
  n_active_covariates_per_group = 2L
))
variants <- c("correct_groups", "covariate_only", "permuted_groups")
scenario_grid <- CJ(rep = seq_len(100L), variant = variants)
checkpoint_file <- file.path(processed_dir, "fssgl_v2_structure_ablation_running.rds")
completed <- if (file.exists(checkpoint_file)) readRDS(checkpoint_file) else data.table()
results <- list()
result_index <- 1L

make_permuted_membership <- function(membership, seed) {
  out <- copy(membership)[order(covariate_id)]
  original <- out$group_id
  original_pairs <- outer(original, original, `==`)
  set.seed(seed)
  for (attempt in seq_len(100L)) {
    candidate <- sample(original, replace = FALSE)
    candidate_pairs <- outer(candidate, candidate, `==`)
    if (!identical(candidate_pairs, original_pairs)) break
  }
  if (identical(candidate_pairs, original_pairs)) {
    stop("Could not construct a genuinely different balanced grouping.")
  }
  off_diagonal <- upper.tri(original_pairs)
  disagreement <- mean(original_pairs[off_diagonal] != candidate_pairs[off_diagonal])
  out[, group_id := candidate]
  list(membership = out, pair_disagreement = disagreement)
}

for (i in seq_len(nrow(scenario_grid))) {
  scenario <- scenario_grid[i]
  if (
    nrow(completed) > 0L &&
      nrow(completed[rep == scenario$rep & variant == scenario$variant]) > 0L
  ) next
  seed <- 2026105000L + scenario$rep
  cat("FSSGL ablation: ", scenario$variant, ", rep=", scenario$rep, "/100\n", sep = "")
  dgp <- do.call(generate_fssgl_simulation_dgp, c(dgp_defaults, list(seed = seed)))
  fit_membership <- dgp$structural_membership
  fit_parameters <- parameters
  pair_disagreement <- 0
  if (scenario$variant == "covariate_only") {
    fit_parameters$group_penalty_scale <- 0
  } else if (scenario$variant == "permuted_groups") {
    permuted <- make_permuted_membership(dgp$structural_membership, seed + 700000L)
    fit_membership <- permuted$membership
    pair_disagreement <- permuted$pair_disagreement
  }
  started <- proc.time()[["elapsed"]]
  fit_result <- tryCatch(
    fit_fssgl(
      x_coef = dgp$x_train,
      y_coef = dgp$y_train,
      structural_membership = fit_membership,
      parameters = fit_parameters,
      verbose = FALSE
    ),
    error = function(e) e
  )
  runtime_sec <- proc.time()[["elapsed"]] - started
  identification <- data.table(
    algorithm_version = FSSGL_ALGORITHM_VERSION,
    variant = scenario$variant,
    rep = scenario$rep,
    seed = seed,
    p = 20L,
    n_train = 40L,
    n_test = 24L,
    n_groups = 4L,
    active_group_count = 2L,
    active_covariate_count = length(dgp$active_covariates),
    group_penalty_scale = fit_parameters$group_penalty_scale,
    grouping_pair_disagreement = pair_disagreement,
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
      posterior_cutoff = fit_parameters$posterior_cutoff,
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
setorder(replicates, variant, rep)
metric_columns <- c(
  "covariate_tpr", "covariate_fpr", "covariate_fdr", "selected_covariate_count",
  "group_tpr", "group_fpr", "group_fdr", "selected_group_count",
  "coefficient_relative_error", "active_coefficient_relative_error",
  "inactive_coefficient_norm", "test_coeff_rmse", "runtime_sec", "final_iter",
  "theta_covariate_final", "theta_group_final", "sigma2_final"
)
summary_long <- melt(
  replicates[status == "ok"],
  id.vars = "variant",
  measure.vars = metric_columns,
  variable.name = "metric",
  value.name = "value"
)[, .(
  n = .N,
  mean = mean(value, na.rm = TRUE),
  sd = sd(value, na.rm = TRUE)
), by = .(variant, metric)]
summary_wide <- dcast(
  summary_long[, mean_sd := sprintf("%.4f (%.4f)", mean, sd)],
  variant ~ metric,
  value.var = "mean_sd"
)

paired_metrics <- c(
  "covariate_tpr", "covariate_fpr", "covariate_fdr", "selected_covariate_count",
  "coefficient_relative_error", "test_coeff_rmse"
)
paired <- merge(
  replicates[variant != "correct_groups"],
  replicates[variant == "correct_groups", c("rep", "seed", paired_metrics), with = FALSE],
  by = c("rep", "seed"),
  suffixes = c("", "_correct")
)
paired_summary <- rbindlist(lapply(paired_metrics, function(metric) {
  paired[, {
    difference <- get(metric) - get(paste0(metric, "_correct"))
    standard_error <- sd(difference) / sqrt(.N)
    list(
      n = .N,
      mean_difference = mean(difference),
      sd_difference = sd(difference),
      ci95_lower = mean(difference) - qt(0.975, .N - 1L) * standard_error,
      ci95_upper = mean(difference) + qt(0.975, .N - 1L) * standard_error
    )
  }, by = variant][, metric := metric]
}))
convergence <- replicates[, .(
  n_attempted = .N,
  n_ok = sum(status == "ok"),
  strict_rate = mean(strict_converged[status == "ok"], na.rm = TRUE),
  error_count = sum(status != "ok")
), by = variant]

fwrite(replicates, file.path(table_dir, "fssgl_v2_structure_ablation_replicates.csv"))
fwrite(summary_long, file.path(table_dir, "fssgl_v2_structure_ablation_summary_long.csv"))
fwrite(summary_wide, file.path(table_dir, "fssgl_v2_structure_ablation_summary_mean_sd.csv"))
fwrite(paired_summary, file.path(table_dir, "fssgl_v2_structure_ablation_paired_vs_correct.csv"))
fwrite(convergence, file.path(table_dir, "fssgl_v2_structure_ablation_convergence.csv"))
saveRDS(
  new_experiment_manifest(
    experiment_id = "fssgl_v2_structure_ablation_formal",
    parameters = parameters,
    design = list(dgp_defaults = dgp_defaults, variants = variants, repetitions = 100L),
    algorithm_version = FSSGL_ALGORITHM_VERSION
  ),
  file.path(processed_dir, "fssgl_v2_structure_ablation_manifest.rds")
)
if (nrow(replicates) == nrow(scenario_grid)) unlink(checkpoint_file)

print(summary_wide[, .(
  variant, covariate_tpr, covariate_fpr, covariate_fdr,
  selected_covariate_count, coefficient_relative_error, test_coeff_rmse
)])
print(convergence)
