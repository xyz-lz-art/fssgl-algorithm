# Baselines paired with 11_fair_fssgl_comparison.R. All methods receive the
# identical generated data, curve-level folds, and response-coefficient RMSE.
# Outputs are isolated from the frozen v2 evidence.

library(data.table)

source("R/parameters.R")
source("R/fssgl/basis_design.R")
source("R/fssgl/tuning.R")
source("R/fssgl/simulation_design.R")
source("R/baselines/functional_methods.R")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "results/tables/simulation/v3_fair")
processed_dir <- file.path(root, "data/processed/simulation/v3_fair")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

parse_integer_list <- function(value, default) {
  if (!nzchar(value)) return(default)
  out <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
  if (!length(out) || anyNA(out) || any(out < 1L)) {
    stop("Environment list must contain comma-separated positive integers.")
  }
  unique(out)
}

p_values <- parse_integer_list(
  Sys.getenv("FSSGL_FAIR_P", ""),
  c(10L, 20L, 50L, 60L, 100L)
)
forced_repetitions <- as.integer(Sys.getenv("FSSGL_FAIR_REPS", "0"))
repetition_count <- function(p) {
  if (forced_repetitions > 0L) forced_repetitions else if (p <= 20L) 100L else 50L
}
methods_for_p <- function(p) {
  if (p <= 20L) {
    c(
      "basis_ridge", "basis_group_scad", "fpca_ridge", "fpca_group_lasso",
      "fpca_group_scad", "fpca_aenet", "structured_group_lasso", "kernel_ridge"
    )
  } else {
    c("basis_group_scad", "fpca_group_scad", "structured_group_lasso")
  }
}
if (!requireNamespace("grpreg", quietly = TRUE)) {
  stop(
    "The fair baseline experiment requires grpreg; run the project setup before fitting."
  )
}

n_folds <- 5L
dgp_defaults <- fssgl_main_dgp_defaults()
baseline_config <- fssgl_baseline_parameters()
scenario_grid <- rbindlist(lapply(p_values, function(p) {
  CJ(p = p, rep = seq_len(repetition_count(p)), method_id = methods_for_p(p))
}))
result_file <- file.path(table_dir, "baseline_v3_fair_replicates.csv")
checkpoint_file <- file.path(processed_dir, "baseline_v3_fair_running.rds")
results <- if (file.exists(checkpoint_file)) {
  readRDS(checkpoint_file)
} else if (file.exists(result_file)) {
  fread(result_file)
} else {
  data.table()
}
# Completed rows survive an interrupted or partially failed run. Error rows are
# retried after the dependency or implementation problem has been corrected.
if (nrow(results) && "status" %in% names(results)) {
  results <- results[status != "error"]
}

for (scenario_id in seq_len(nrow(scenario_grid))) {
  scenario <- scenario_grid[scenario_id]
  if (nrow(results) && any(
    results$p == scenario$p & results$rep == scenario$rep &
      results$method_id == scenario$method_id
  )) next

  p <- scenario$p
  replicate_id <- scenario$rep
  method_id <- scenario$method_id
  seed <- if (p <= 20L) {
    2026101000L + p * 1000L + replicate_id
  } else {
    2026110000L + p * 1000L + replicate_id
  }
  n_groups <- if (p == 10L) 2L else if (p == 20L) 4L else p %/% 5L
  n_active_groups <- if (p == 10L) 1L else 2L
  dgp <- do.call(generate_fssgl_simulation_dgp, c(dgp_defaults, list(
    n_covariates = p,
    n_groups = n_groups,
    n_active_groups = n_active_groups,
    n_active_covariates_per_group = 2L,
    seed = seed
  )))
  fold_seed <- seed + 900000L
  folds <- make_curve_cv_folds(dgp$dimensions$n_train, n_folds, fold_seed)
  cat(
    "Fair baseline: p=", p, ", rep=", replicate_id, ", method=", method_id,
    "\n", sep = ""
  )

  total_started <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    fit_functional_baseline(
      method_id = method_id,
      x_coef = dgp$x_train,
      y_coef = dgp$y_train,
      folds = folds,
      config = baseline_config,
      structural_membership = dgp$structural_membership
    ),
    error = function(e) e
  )
  total_runtime_sec <- proc.time()[["elapsed"]] - total_started
  identification <- data.table(
    tuning_protocol_version = FSSGL_TUNING_PROTOCOL_VERSION,
    method_id = method_id,
    p = p,
    rep = replicate_id,
    seed = seed,
    fold_seed = fold_seed,
    validation_loss = "response-coefficient RMSE"
  )
  if (inherits(fit, "error")) {
    result <- cbind(identification, data.table(
      status = "error",
      error_message = conditionMessage(fit),
      total_runtime_sec = total_runtime_sec
    ))
  } else {
    metrics <- evaluate_functional_baseline(fit, dgp, total_runtime_sec)
    result <- cbind(identification, data.table(
      status = if (isFALSE(fit$solver_converged)) "nonconverged" else "ok",
      error_message = NA_character_,
      selected_lambda = if (!is.null(fit$tuning$selected_lambda)) {
        fit$tuning$selected_lambda
      } else {
        NA_real_
      },
      selected_fraction = if (!is.null(fit$tuning$selected_fraction)) {
        fit$tuning$selected_fraction
      } else {
        NA_real_
      }
    ), metrics)
  }
  results <- rbindlist(list(results, result), use.names = TRUE, fill = TRUE)
  saveRDS(results, checkpoint_file)
}

setorder(results, p, rep, method_id)
fwrite(results, result_file)
reportable <- results[status == "ok"]
metric_columns <- c(
  "covariate_tpr", "covariate_fdr", "coefficient_relative_error",
  "test_coeff_rmse", "tuning_runtime_sec", "final_runtime_sec", "runtime_sec"
)
summary <- melt(
  reportable,
  id.vars = c("p", "method_id"),
  measure.vars = metric_columns,
  variable.name = "metric",
  value.name = "value"
)[is.finite(value), .(
  n = .N,
  mean = mean(value),
  sd = stats::sd(value),
  mcse = stats::sd(value) / sqrt(.N)
), by = .(p, method_id, metric)]
fwrite(summary, file.path(table_dir, "baseline_v3_fair_summary.csv"))

saveRDS(
  new_experiment_manifest(
    experiment_id = "baseline_v3_fair_common_fold_comparison",
    parameters = list(
      selection_rules = paste(
        "minimum validation RMSE for ridge and kernel methods;",
        "one-standard-error sparse rule for group lasso, group SCAD,",
        "adaptive elastic net, and structured group lasso"
      ),
      configuration = baseline_config
    ),
    design = list(
      p = p_values,
      repetitions = vapply(p_values, repetition_count, integer(1L)),
      n_folds = n_folds,
      dgp_seed_rule = "shared with the paired FSSGL v3 experiment",
      fold_seed_rule = "DGP seed + 900000",
      validation_loss = "response-coefficient RMSE per held-out curve",
      test_set_used_for_tuning = FALSE,
      tuning_and_final_fit_timed_separately = TRUE
    ),
    algorithm_version = FSSGL_ALGORITHM_VERSION
  ),
  file.path(processed_dir, "baseline_v3_fair_manifest.rds")
)
if (nrow(results) == nrow(scenario_grid)) unlink(checkpoint_file)

print(results[, .N, by = .(p, method_id, status)])
print(summary)
