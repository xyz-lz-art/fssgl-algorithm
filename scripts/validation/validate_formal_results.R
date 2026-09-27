# Validate the complete frozen-parameter simulation suite.

library(data.table)

source("R/parameters.R")

table_dir <- "results/tables/simulation/v2_formal"
processed_dir <- "data/processed/simulation/v2_formal"

fail <- function(...) stop(sprintf(...), call. = FALSE)

read_required <- function(filename, columns) {
  path <- file.path(table_dir, filename)
  if (!file.exists(path)) fail("Missing result file: %s", path)
  out <- fread(path)
  missing_columns <- setdiff(columns, names(out))
  if (length(missing_columns) > 0L) {
    fail("%s is missing columns: %s", path, paste(missing_columns, collapse = ", "))
  }
  out
}

require_unique <- function(data, keys, label) {
  if (anyDuplicated(data, by = keys)) fail("%s has duplicate rows", label)
}

require_success <- function(data, label, require_strict = FALSE) {
  if (any(is.na(data$status)) || any(data$status != "ok")) {
    fail("%s contains failed fits", label)
  }
  if (require_strict && any(!data$strict_converged)) {
    fail("%s contains non-strictly-converged fits", label)
  }
}

main <- read_required(
  "fssgl_v2_main_replicates.csv",
  c("algorithm_version", "p", "rep", "seed", "status", "strict_converged")
)
if (nrow(main) != 200L || !setequal(main$p, c(10L, 20L))) {
  fail("Main experiment must have 200 rows at p = 10 and 20")
}
require_unique(main, c("p", "rep"), "main experiment")
require_success(main, "main experiment", require_strict = TRUE)
if (any(main$algorithm_version != FSSGL_ALGORITHM_VERSION) ||
    any(main[, .N, by = p]$N != 100L)) {
  fail("Main experiment version or cell counts are invalid")
}

robustness <- read_required(
  "fssgl_v2_robustness_replicates.csv",
  c("algorithm_version", "scenario_id", "rep", "seed", "status", "strict_converged")
)
expected_scenarios <- c("baseline", "snr_low", "correlation_high", "within_group_sparse")
if (nrow(robustness) != 400L || !setequal(robustness$scenario_id, expected_scenarios)) {
  fail("Robustness experiment must have four 100-row scenarios")
}
require_unique(robustness, c("scenario_id", "rep"), "robustness experiment")
require_success(robustness, "robustness experiment", require_strict = TRUE)
if (any(robustness[, .N, by = scenario_id]$N != 100L)) {
  fail("Robustness scenario counts are invalid")
}

responsibility <- read_required(
  "fssgl_v2_responsibility_diagnostic_replicates.csv",
  c(
    "scenario_id", "rep", "seed", "strict_converged", "final_iter",
    "covariate_tpr", "group_tpr", "theta_covariate_final", "theta_group_final"
  )
)
expected_responsibility_scenarios <- c("baseline", "within_group_sparse")
if (nrow(responsibility) != 200L ||
    !setequal(responsibility$scenario_id, expected_responsibility_scenarios) ||
    any(!responsibility$strict_converged) ||
    any(responsibility[, .N, by = scenario_id]$N != 100L)) {
  fail("Responsibility diagnostic must contain two strictly converged 100-fit scenarios")
}
require_unique(responsibility, c("scenario_id", "rep"), "responsibility diagnostic")
responsibility_reference <- robustness[
  scenario_id %in% expected_responsibility_scenarios,
  .(scenario_id, rep, seed)
]
if (!fsetequal(
  responsibility[, .(scenario_id, rep, seed)],
  responsibility_reference
)) {
  fail("Responsibility diagnostic does not match the paired robustness seeds")
}
responsibility_trace_file <- file.path(
  processed_dir, "fssgl_v2_responsibility_trajectories.rds"
)
if (!file.exists(responsibility_trace_file)) {
  fail("Missing responsibility trajectory evidence: %s", responsibility_trace_file)
}
responsibility_trace <- as.data.table(readRDS(responsibility_trace_file))
required_trace_columns <- c(
  "scenario_id", "rep", "iter", "final_iter", "level", "unit_id", "role",
  "beta_norm", "posterior_slab_prob", "transition_norm", "norm_to_transition"
)
if (length(setdiff(required_trace_columns, names(responsibility_trace))) ||
    any(!is.finite(responsibility_trace$beta_norm)) ||
    any(!is.finite(responsibility_trace$posterior_slab_prob)) ||
    any(!is.finite(responsibility_trace$transition_norm)) ||
    any(responsibility_trace$posterior_slab_prob < 0) ||
    any(responsibility_trace$posterior_slab_prob > 1) ||
    any(responsibility_trace$transition_norm <= 0)) {
  fail("Responsibility trajectory values or schema are invalid")
}
trace_counts <- responsibility_trace[, .(
  trace_rows = .N,
  trace_final_iter = max(iter),
  unit_count = uniqueN(paste(level, unit_id, sep = ":"))
), by = .(scenario_id, rep)]
trace_counts <- merge(
  trace_counts,
  responsibility[, .(scenario_id, rep, final_iter)],
  by = c("scenario_id", "rep")
)
if (nrow(trace_counts) != 200L ||
    any(trace_counts$trace_final_iter != trace_counts$final_iter) ||
    any(trace_counts$unit_count != 24L) ||
    any(trace_counts$trace_rows != 24L * trace_counts$final_iter)) {
  fail("Responsibility trajectory does not cover every unit and iteration")
}
responsibility_trajectory_summary <- read_required(
  "fssgl_v2_responsibility_trajectory_summary.csv",
  c(
    "scenario_id", "level", "progress_bin", "n_repetitions",
    "posterior_median", "norm_ratio_median"
  )
)
if (nrow(responsibility_trajectory_summary) != 84L ||
    any(responsibility_trajectory_summary$n_repetitions != 100L)) {
  fail("Responsibility trajectory summary must contain 21 points for four paths")
}

ablation <- read_required(
  "fssgl_v2_structure_ablation_replicates.csv",
  c("algorithm_version", "variant", "rep", "seed", "status", "strict_converged")
)
expected_variants <- c("correct_groups", "covariate_only", "permuted_groups")
if (nrow(ablation) != 300L || !setequal(ablation$variant, expected_variants)) {
  fail("Structure ablation must have three 100-row variants")
}
require_unique(ablation, c("variant", "rep"), "structure ablation")
require_success(ablation, "structure ablation", require_strict = TRUE)
if (any(ablation[, .N, by = variant]$N != 100L)) {
  fail("Structure-ablation variant counts are invalid")
}

baselines <- read_required(
  "functional_baselines_v2_replicates.csv",
  c("reference_algorithm_version", "method_id", "p", "rep", "seed", "status")
)
expected_baselines <- c(
  "basis_ridge", "fpca_ridge", "fpca_group_lasso", "fpca_aenet",
  "fpca_group_scad", "structured_group_lasso", "kernel_ridge"
)
if (nrow(baselines) != 1400L || !setequal(baselines$method_id, expected_baselines)) {
  fail("Baselines must have 1400 rows across seven methods")
}
require_unique(baselines, c("method_id", "p", "rep"), "baselines")
require_success(baselines, "baselines")
if (any(baselines[, .N, by = .(method_id, p)]$N != 100L) ||
    any(baselines$reference_algorithm_version != FSSGL_ALGORITHM_VERSION)) {
  fail("Baseline version or cell counts are invalid")
}

main_keys <- unique(main[, .(p, rep, seed)])
baseline_keys <- unique(baselines[, .(p, rep, seed)])
if (!fsetequal(main_keys, baseline_keys)) {
  fail("Baseline DGP keys do not exactly match the main experiment")
}

comparison <- read_required(
  "functional_method_comparison_v2_replicates.csv",
  c("method_id", "p", "rep", "seed", "status")
)
if (nrow(comparison) != 1600L ||
    !setequal(comparison$method_id, c("fssgl_v2", expected_baselines))) {
  fail("Combined method comparison must have 1600 rows across eight methods")
}
require_unique(comparison, c("method_id", "p", "rep"), "method comparison")
require_success(comparison, "method comparison")

highdim <- read_required(
  "fssgl_v2_high_dimensional_replicates.csv",
  c("algorithm_version", "method_id", "p", "rep", "seed", "status")
)
expected_highdim_methods <- c(
  "fssgl_v2", "fpca_group_scad", "structured_group_lasso"
)
if (nrow(highdim) != 450L ||
    !setequal(highdim$p, c(50L, 60L, 100L)) ||
    !setequal(highdim$method_id, expected_highdim_methods)) {
  fail("High-dimensional comparison must have 450 rows in nine 50-row cells")
}
require_unique(highdim, c("method_id", "p", "rep"), "high-dimensional comparison")
require_success(highdim, "high-dimensional comparison")
if (any(highdim[, .N, by = .(method_id, p)]$N != 50L) ||
    any(highdim$algorithm_version != FSSGL_ALGORITHM_VERSION)) {
  fail("High-dimensional comparison version or cell counts are invalid")
}

highdim_recalibrated <- read_required(
  "fssgl_v2_high_dimensional_recalibrated_replicates.csv",
  c(
    "algorithm_version", "method_id", "p", "rep", "seed", "status",
    "selected_multiplier", "strict_converged", "relaxed_converged", "final_iter"
  )
)
if (nrow(highdim_recalibrated) != 150L ||
    !setequal(highdim_recalibrated$p, c(50L, 60L, 100L))) {
  fail("Recalibrated high-dimensional FSSGL must have three 50-row cells")
}
require_unique(
  highdim_recalibrated,
  c("p", "rep"),
  "recalibrated high-dimensional FSSGL"
)
require_success(highdim_recalibrated, "recalibrated high-dimensional FSSGL")
if (any(highdim_recalibrated$algorithm_version != FSSGL_ALGORITHM_VERSION) ||
    any(highdim_recalibrated[, .N, by = p]$N != 50L) ||
    any(highdim_recalibrated$strict_converged) ||
    any(!highdim_recalibrated$relaxed_converged) ||
    any(highdim_recalibrated$final_iter != 500L)) {
  fail("Recalibrated high-dimensional version, counts, or disclosed convergence status are invalid")
}

highdim_pilot <- read_required(
  "fssgl_v2_high_dimensional_pilot_calibration.csv",
  c("p", "pilot_rep", "seed", "multiplier", "validation_rmse", "status")
)
if (nrow(highdim_pilot) != 180L ||
    anyDuplicated(highdim_pilot, by = c("p", "pilot_rep", "multiplier")) ||
    any(highdim_pilot$status != "ok") ||
    any(!is.finite(highdim_pilot$validation_rmse))) {
  fail("High-dimensional pilot calibration is incomplete, duplicated, or failed")
}

highdim_scales <- read_required(
  "fssgl_v2_high_dimensional_selected_multipliers.csv",
  c("p", "multiplier", "mean_validation_rmse")
)
expected_highdim_scales <- data.table(
  p = c(50L, 60L, 100L),
  multiplier = c(0.001953125, 0.00390625, 0.015625)
)
if (nrow(highdim_scales) != 3L ||
    !isTRUE(all.equal(
      highdim_scales[order(p), .(p, multiplier)],
      expected_highdim_scales,
      check.attributes = FALSE
    ))) {
  fail("High-dimensional selected multipliers do not match the frozen pilot result")
}
pilot_minima <- highdim_pilot[, .(
  mean_validation_rmse = mean(validation_rmse)
), by = .(p, multiplier)][order(p, mean_validation_rmse, multiplier), .SD[1L], by = p]
replicate_scales <- unique(
  highdim_recalibrated[, .(p, multiplier = selected_multiplier)]
)[order(p)]
if (!isTRUE(all.equal(
      pilot_minima[order(p), .(p, multiplier)],
      highdim_scales[order(p), .(p, multiplier)],
      check.attributes = FALSE
    )) ||
    !isTRUE(all.equal(
      replicate_scales,
      highdim_scales[order(p), .(p, multiplier)],
      check.attributes = FALSE
    ))) {
  fail(paste(
    "Frozen high-dimensional pilot minima, selected-multiplier table, and",
    "final replicate file disagree"
  ))
}

identifiability <- read_required(
  "fssgl_v2_identifiability_diagnostics.csv",
  c(
    "p", "rep", "seed", "n_rows", "n_columns", "numerical_rank",
    "entropy_effective_rank", "gram_condition", "nonzero_gram_condition",
    "design_rank_deficient", "severely_ill_conditioned"
  )
)
if (nrow(identifiability) != 250L ||
    !setequal(identifiability$p, c(10L, 20L, 50L, 60L, 100L)) ||
    anyDuplicated(identifiability, by = c("p", "rep")) ||
    any(identifiability$numerical_rank != 39L) ||
    any(!identifiability$design_rank_deficient)) {
  fail("Simulation identifiability diagnostics are incomplete or inconsistent")
}

basis_fixed <- read_required(
  "fssgl_v2_basis_sensitivity_replicates.csv",
  c("algorithm_version", "basis_dimension", "rep", "seed", "status", "strict_converged")
)
if (nrow(basis_fixed) != 300L || !setequal(basis_fixed$basis_dimension, 3:5)) {
  fail("Fixed-penalty basis experiment must have three 100-row dimensions")
}
require_unique(basis_fixed, c("basis_dimension", "rep"), "fixed-penalty basis experiment")
require_success(basis_fixed, "fixed-penalty basis experiment", require_strict = TRUE)
if (any(basis_fixed$algorithm_version != FSSGL_ALGORITHM_VERSION) ||
    any(basis_fixed[, .N, by = basis_dimension]$N != 100L)) {
  fail("Fixed-penalty basis experiment version or cell counts are invalid")
}

basis_recalibrated <- read_required(
  "fssgl_v2_basis_recalibration_replicates.csv",
  c(
    "algorithm_version", "basis_dimension", "rep", "seed", "status",
    "strict_converged", "selected_multiplier"
  )
)
if (nrow(basis_recalibrated) != 300L ||
    !setequal(basis_recalibrated$basis_dimension, 3:5)) {
  fail("Recalibrated basis experiment must have three 100-row dimensions")
}
require_unique(
  basis_recalibrated,
  c("basis_dimension", "rep"),
  "recalibrated basis experiment"
)
require_success(
  basis_recalibrated,
  "recalibrated basis experiment",
  require_strict = TRUE
)
if (any(basis_recalibrated$algorithm_version != FSSGL_ALGORITHM_VERSION) ||
    any(basis_recalibrated[, .N, by = basis_dimension]$N != 100L) ||
    !all(basis_recalibrated$selected_multiplier %in% c(0.125, 0.25, 0.5, 1, 2, 4))) {
  fail("Recalibrated basis experiment version, cell counts, or selected multipliers are invalid")
}

basis_cv <- read_required(
  "fssgl_v2_basis_recalibration_cv_path.csv",
  c("basis_dimension", "rep", "multiplier", "cv_curve_l2_rmse", "cv_error_count", "selected")
)
if (nrow(basis_cv) != 1800L ||
    anyDuplicated(basis_cv, by = c("basis_dimension", "rep", "multiplier")) ||
    any(!is.finite(basis_cv$cv_curve_l2_rmse)) ||
    any(basis_cv$cv_error_count != 0L) ||
    any(basis_cv[, sum(selected), by = .(basis_dimension, rep)]$V1 != 1L)) {
  fail("Basis recalibration path is incomplete, duplicated, failed, or lacks one selection per fit")
}

manifest_files <- c(
  "fssgl_v2_main_manifest.rds",
  "fssgl_v2_robustness_manifest.rds",
  "fssgl_v2_structure_ablation_manifest.rds",
  "functional_baselines_v2_manifest.rds",
  "fssgl_v2_basis_sensitivity_manifest.rds",
  "fssgl_v2_basis_recalibration_manifest.rds",
  "fssgl_v2_high_dimensional_manifest.rds",
  "fssgl_v2_high_dimensional_recalibration_manifest.rds",
  "fssgl_v2_responsibility_diagnostic_manifest.rds"
)
for (filename in manifest_files) {
  path <- file.path(processed_dir, filename)
  if (!file.exists(path)) fail("Missing formal manifest: %s", path)
  manifest <- readRDS(path)
  if (!identical(manifest$algorithm_version, FSSGL_ALGORITHM_VERSION) ||
      !identical(manifest$result_schema_version, FSSGL_RESULT_SCHEMA_VERSION)) {
    fail("Wrong version in manifest: %s", path)
  }
}
main_manifest <- readRDS(file.path(processed_dir, "fssgl_v2_main_manifest.rds"))
if (!isTRUE(all.equal(main_manifest$parameters, fssgl_parameters()))) {
  fail("Main manifest parameters differ from the frozen configuration")
}

running <- list.files(processed_dir, pattern = "_running[.]rds$", full.names = TRUE)
if (length(running) > 0L) {
  fail("Unfinished formal checkpoints remain: %s", paste(running, collapse = ", "))
}

cat("Formal FSSGL output validation passed.\n")
cat("  Main: 200 strictly converged fits\n")
cat("  Robustness: 400 strictly converged fits\n")
cat("  Responsibility trajectories: 200 strictly converged traced fits\n")
cat("  Structure ablation: 300 strictly converged fits\n")
cat("  Baselines: 1400 completed fits\n")
cat("  Combined method comparison: 1600 complete rows\n")
cat("  High-dimensional comparison: 450 completed fits\n")
cat("  High-dimensional pilot calibration: 180 candidates and 150 final fits\n")
cat("  Identifiability diagnostics: 250 simulated designs\n")
cat("  Basis sensitivity: 300 fixed and 300 recalibrated strictly converged fits\n")
cat("  Basis calibration: 1800 complete cross-validation candidates\n")
