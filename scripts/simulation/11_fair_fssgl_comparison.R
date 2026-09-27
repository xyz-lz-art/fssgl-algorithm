# Submission-quality FSSGL tuning on the same curve-level folds and loss used by
# the comparison methods. This script writes v3 outputs and never overwrites the
# frozen v2 results reported in the current draft.

library(data.table)

source("R/parameters.R")
source("R/fssgl/basis_design.R")
source("R/fssgl/penalty_weights.R")
source("R/fssgl/weighted_membership_solver.R")
source("R/fssgl/solver.R")
source("R/fssgl/tuning.R")
source("R/fssgl/simulation_design.R")

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
n_folds <- 5L
selection_rule <- "one_se_sparsest"
multiplier_grid <- fssgl_penalty_multiplier_grid(
  core_grid = 2^seq(-11, 5, by = 2),
  lower_expansions = 3L
)
dgp_defaults <- fssgl_main_dgp_defaults()
scenario_grid <- rbindlist(lapply(p_values, function(p) {
  data.table(p = p, rep = seq_len(repetition_count(p)))
}))
shard_count <- as.integer(Sys.getenv("FSSGL_FAIR_SHARD_COUNT", "1"))
shard_id <- as.integer(Sys.getenv("FSSGL_FAIR_SHARD_ID", "1"))
if (length(shard_count) != 1L || is.na(shard_count) || shard_count < 1L ||
    length(shard_id) != 1L || is.na(shard_id) ||
    shard_id < 1L || shard_id > shard_count) {
  stop("FSSGL_FAIR_SHARD_ID must be between 1 and FSSGL_FAIR_SHARD_COUNT.")
}
if (shard_count > 1L) {
  scenario_grid <- scenario_grid[(rep - 1L) %% shard_count == shard_id - 1L]
}
suffix <- if (shard_count == 1L) "" else paste0("_shard", shard_id, "of", shard_count)

result_file <- file.path(table_dir, paste0("fssgl_v3_fair_replicates", suffix, ".csv"))
tuning_file <- file.path(table_dir, paste0("fssgl_v3_fair_cv_path", suffix, ".csv"))
checkpoint_file <- file.path(processed_dir, paste0("fssgl_v3_fair_running", suffix, ".rds"))
state <- if (file.exists(checkpoint_file)) {
  readRDS(checkpoint_file)
} else {
  list(results = data.table(), tuning = data.table())
}

for (scenario_id in seq_len(nrow(scenario_grid))) {
  scenario <- scenario_grid[scenario_id]
  if (nrow(state$results) &&
      any(state$results$p == scenario$p & state$results$rep == scenario$rep)) {
    next
  }
  p <- scenario$p
  replicate_id <- scenario$rep
  # Match the already specified baseline DGPs exactly: scripts 04 (p <= 20)
  # and 09 (p >= 50) use different experiment prefixes.
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
  folds <- make_curve_cv_folds(
    n = dgp$dimensions$n_train,
    n_folds = n_folds,
    seed = seed + 900000L
  )
  base_parameters <- fssgl_submission_parameters(p, dgp$dimensions$n_train)
  extension_factors <- if (p >= dgp$dimensions$n_train) c(1L, 2L) else 1L
  cat(
    "Fair FSSGL: p=", p, ", rep=", replicate_id, "/", repetition_count(p),
    "\n", sep = ""
  )

  tuning <- tryCatch(
    tune_fssgl_cv(
      x_coef = dgp$x_train,
      y_coef = dgp$y_train,
      structural_membership = dgp$structural_membership,
      folds = folds,
      multiplier_grid = multiplier_grid,
      parameters = base_parameters,
      selection_rule = selection_rule,
      extension_factors = extension_factors,
      require_interior_minimum = FALSE,
      verbose = FALSE
    ),
    error = function(e) e
  )

  identification <- data.table(
    algorithm_version = FSSGL_ALGORITHM_VERSION,
    tuning_protocol_version = FSSGL_TUNING_PROTOCOL_VERSION,
    p = p,
    rep = replicate_id,
    seed = seed,
    fold_seed = seed + 900000L,
    selection_rule = selection_rule
  )
  if (inherits(tuning, "error")) {
    result <- cbind(identification, data.table(
      status = "tuning_error",
      error_message = conditionMessage(tuning)
    ))
  } else {
    tuning_rows <- copy(tuning$path)
    tuning_rows[, `:=`(
      p = p,
      rep = replicate_id,
      seed = seed,
      fold_seed = seed + 900000L
    )]
    state$tuning <- rbindlist(
      list(state$tuning, tuning_rows), use.names = TRUE, fill = TRUE
    )

    boundary_selected <- tuning$selected_at_lower_boundary ||
      (tuning$selected_at_upper_boundary && !tuning$upper_null_plateau)
    if (boundary_selected) {
      result <- cbind(identification, data.table(
        status = "grid_boundary",
        error_message = paste(
          "Selected value lies at the",
          if (tuning$selected_at_lower_boundary) "lower" else "upper",
          "grid boundary; expand the grid before fitting or reporting."
        ),
        selected_multiplier = tuning$selected_multiplier,
        upper_null_plateau = tuning$upper_null_plateau,
        tuning_runtime_sec = tuning$tuning_runtime_sec,
        final_runtime_sec = NA_real_,
        total_runtime_sec = tuning$tuning_runtime_sec
      ))
    } else {
      final_started <- proc.time()[["elapsed"]]
      final_path <- tryCatch(
        fit_fssgl_selected_path(
          x_coef = dgp$x_train,
          y_coef = dgp$y_train,
          structural_membership = dgp$structural_membership,
          multiplier_grid = multiplier_grid,
          selected_multiplier = tuning$selected_multiplier,
          parameters = base_parameters,
          extension_factors = extension_factors,
          verbose = FALSE
        ),
        error = function(e) e
      )
      final_runtime_sec <- proc.time()[["elapsed"]] - final_started
      if (inherits(final_path, "error")) {
        result <- cbind(identification, data.table(
          status = "fit_error",
          error_message = conditionMessage(final_path),
          selected_multiplier = tuning$selected_multiplier,
          tuning_runtime_sec = tuning$tuning_runtime_sec,
          final_runtime_sec = final_runtime_sec,
          total_runtime_sec = tuning$tuning_runtime_sec + final_runtime_sec
        ))
      } else {
        fit <- final_path$fit
        status <- if (isTRUE(fit$fit$convergence$strict_converged)) {
          "ok"
        } else {
          "nonconverged"
        }
        metrics <- evaluate_fssgl_simulation_fit(
          fit,
          dgp,
          posterior_cutoff = tuning$selected_parameters$posterior_cutoff,
          runtime_sec = final_runtime_sec
        )
        result <- cbind(identification, data.table(
          status = status,
          error_message = NA_character_,
          selected_multiplier = tuning$selected_multiplier,
          selected_at_lower_boundary = tuning$selected_at_lower_boundary,
          selected_at_upper_boundary = tuning$selected_at_upper_boundary,
          upper_null_plateau = tuning$upper_null_plateau,
          final_path_fits = nrow(final_path$path),
          tuning_runtime_sec = tuning$tuning_runtime_sec,
          final_runtime_sec = final_runtime_sec,
          total_runtime_sec = tuning$tuning_runtime_sec + final_runtime_sec
        ), metrics)
        if (p == 20L && replicate_id == 1L && status == "ok") {
          saveRDS(
            list(dgp = dgp, fit = fit, tuning = tuning),
            file.path(processed_dir, "fssgl_v3_surface_example_p20_rep1.rds")
          )
        }
      }
    }
  }
  state$results <- rbindlist(
    list(state$results, result), use.names = TRUE, fill = TRUE
  )
  saveRDS(state, checkpoint_file)
}

setorder(state$results, p, rep)
if (nrow(state$tuning)) setorder(state$tuning, p, rep, fold_id, multiplier)
fwrite(state$results, result_file)
fwrite(state$tuning, tuning_file)

reportable <- if (all(c("status", "strict_converged") %in% names(state$results))) {
  state$results[status == "ok" & strict_converged == TRUE]
} else {
  data.table()
}
metric_columns <- c(
  "covariate_tpr", "covariate_fdr", "coefficient_relative_error", "test_coeff_rmse",
  "tuning_runtime_sec", "final_runtime_sec", "total_runtime_sec"
)
summary <- if (nrow(reportable)) {
  melt(
    reportable,
    id.vars = "p",
    measure.vars = metric_columns,
    variable.name = "metric",
    value.name = "value"
  )[, .(
    n = .N,
    mean = mean(value),
    sd = stats::sd(value),
    mcse = stats::sd(value) / sqrt(.N)
  ), by = .(p, metric)]
} else {
  data.table()
}
fwrite(summary, file.path(table_dir, paste0("fssgl_v3_fair_summary", suffix, ".csv")))

saveRDS(
  new_experiment_manifest(
    experiment_id = "fssgl_v3_fair_common_fold_comparison",
    parameters = list(
      parameter_rule = "fssgl_submission_parameters(p, n_train)",
      multiplier_grid = multiplier_grid,
      selection_rule = selection_rule,
      extension_factors = "1 for p<n; c(1,2) for p>=n",
      final_refit = "full-training strong-to-selected path, matching CV continuation"
    ),
    design = list(
      p = p_values,
      repetitions = vapply(p_values, repetition_count, integer(1L)),
      n_folds = n_folds,
      dgp_seed_rule = "2026101000 + 1000*p + rep for p<=20; 2026110000 + 1000*p + rep otherwise",
      fold_seed_rule = "DGP seed + 900000",
      validation_loss = "response-coefficient RMSE per held-out curve",
      test_set_used_for_tuning = FALSE,
      boundary_result_reportable = "only a verified upper all-zero plateau is reportable",
      strict_convergence_required = TRUE
    ),
    algorithm_version = FSSGL_ALGORITHM_VERSION
  ),
  file.path(processed_dir, paste0("fssgl_v3_fair_manifest", suffix, ".rds"))
)
if (nrow(state$results) == nrow(scenario_grid)) unlink(checkpoint_file)

print(state$results[, .N, by = .(p, status)])
print(summary)
