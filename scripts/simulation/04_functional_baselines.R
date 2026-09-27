# Refit all functional baselines on the exact 200 DGPs used by the main experiment.

library(data.table)

source("R/parameters.R")
source("R/fssgl/basis_design.R")
source("R/fssgl/tuning.R")
source("R/fssgl/simulation_design.R")
source("R/baselines/functional_methods.R")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "results/tables/simulation/v2_formal")
processed_dir <- file.path(root, "data/processed/simulation/v2_formal")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

method_ids <- c(
  "basis_ridge", "fpca_ridge", "fpca_group_lasso", "fpca_aenet",
  "fpca_group_scad", "structured_group_lasso", "kernel_ridge"
)
baseline_config <- fssgl_baseline_parameters()
dgp_defaults <- fssgl_main_dgp_defaults()
scenario_grid <- CJ(p = c(10L, 20L), rep = seq_len(100L), method_id = method_ids)
checkpoint_file <- file.path(processed_dir, "functional_baselines_v2_running.rds")
final_file <- file.path(table_dir, "functional_baselines_v2_replicates.csv")
completed_sources <- list()
if (file.exists(final_file)) completed_sources[[length(completed_sources) + 1L]] <- fread(final_file)
if (file.exists(checkpoint_file)) completed_sources[[length(completed_sources) + 1L]] <- readRDS(checkpoint_file)
completed <- if (length(completed_sources)) {
  unique(
    rbindlist(completed_sources, fill = TRUE),
    by = c("p", "rep", "method_id"),
    fromLast = TRUE
  )
} else {
  data.table()
}
results <- list()
result_index <- 1L

extract_tuning <- function(fit) {
  data.table(
    selected_lambda = if (!is.null(fit$tuning$selected_lambda)) fit$tuning$selected_lambda else NA_real_,
    selected_gamma_multiplier = if (!is.null(fit$tuning$selected_gamma_multiplier)) fit$tuning$selected_gamma_multiplier else NA_real_,
    selected_lambda_fraction = if (!is.null(fit$tuning$selected_fraction)) fit$tuning$selected_fraction else NA_real_,
    selected_alpha = if (!is.null(fit$tuning$selected_alpha)) fit$tuning$selected_alpha else NA_real_,
    fpca_x_components_mean = if (!is.null(fit$representation)) {
      mean(vapply(fit$representation$x_models, function(model) ncol(model$rotation), integer(1)))
    } else {
      NA_real_
    },
    fpca_y_components = if (!is.null(fit$representation)) {
      ncol(fit$representation$y_model$rotation)
    } else {
      NA_real_
    }
  )
}

for (i in seq_len(nrow(scenario_grid))) {
  scenario <- scenario_grid[i]
  if (
    nrow(completed) > 0L &&
      nrow(completed[p == scenario$p & rep == scenario$rep & method_id == scenario$method_id]) > 0L
  ) next
  seed <- 2026101000L + scenario$p * 1000L + scenario$rep
  n_groups <- if (scenario$p == 10L) 2L else 4L
  n_active_groups <- if (scenario$p == 10L) 1L else 2L
  dgp <- do.call(generate_fssgl_simulation_dgp, c(dgp_defaults, list(
    n_covariates = scenario$p,
    n_groups = n_groups,
    n_active_groups = n_active_groups,
    n_active_covariates_per_group = 2L,
    seed = seed
  )))
  folds <- make_baseline_cv_folds(
    n = dgp$dimensions$n_train,
    n_folds = baseline_config$cv_folds,
    seed = seed + 900000L
  )
  cat("functional baseline: p=", scenario$p, ", rep=", scenario$rep,
      "/100, method=", scenario$method_id, "\n", sep = "")
  started <- proc.time()[["elapsed"]]
  fit_result <- tryCatch(
    fit_functional_baseline(
      method_id = scenario$method_id,
      x_coef = dgp$x_train,
      y_coef = dgp$y_train,
      folds = folds,
      config = baseline_config,
      structural_membership = dgp$structural_membership
    ),
    error = function(e) e
  )
  runtime_sec <- proc.time()[["elapsed"]] - started
  identification <- data.table(
    reference_algorithm_version = FSSGL_ALGORITHM_VERSION,
    method_id = scenario$method_id,
    p = scenario$p,
    rep = scenario$rep,
    seed = seed,
    n_train = dgp$dimensions$n_train,
    n_test = dgp$dimensions$n_test,
    n_groups = n_groups,
    active_group_count = n_active_groups,
    active_covariate_count = length(dgp$active_covariates),
    cv_folds = baseline_config$cv_folds
  )
  if (inherits(fit_result, "error")) {
    result <- cbind(
      identification,
      data.table(status = "error", error_message = conditionMessage(fit_result))
    )
  } else {
    result <- cbind(
      identification,
      data.table(
        status = "ok",
        error_message = NA_character_,
        supports_selection = fit_result$supports_selection,
        supports_coefficients = fit_result$supports_coefficients
      ),
      extract_tuning(fit_result),
      evaluate_functional_baseline(fit_result, dgp, runtime_sec)
    )
  }
  results[[result_index]] <- result
  result_index <- result_index + 1L
  saveRDS(rbindlist(c(list(completed), results), fill = TRUE), checkpoint_file)
}

replicates <- rbindlist(c(list(completed), results), fill = TRUE)
setorder(replicates, p, rep, method_id)
fwrite(replicates, final_file)
saveRDS(
  new_experiment_manifest(
    experiment_id = "functional_baselines_v2_formal",
    parameters = baseline_config,
    design = list(
      dgp_defaults = dgp_defaults,
      scenario_grid = scenario_grid,
      main_seed_rule = "2026101000 + p * 1000 + rep"
    ),
    algorithm_version = FSSGL_ALGORITHM_VERSION
  ),
  file.path(processed_dir, "functional_baselines_v2_manifest.rds")
)
if (nrow(replicates) == nrow(scenario_grid)) unlink(checkpoint_file)

print(replicates[, .(
  n = .N,
  n_ok = sum(status == "ok"),
  mean_runtime_sec = mean(runtime_sec, na.rm = TRUE)
), by = .(method_id, p)])
