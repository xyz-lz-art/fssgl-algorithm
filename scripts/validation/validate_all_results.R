options(stringsAsFactors = FALSE)

expected_shanghai_algorithm_version <- "normalized-adaptive-block-orthonormal-v2.0"
expected_basis_geometry_version <- "orthonormal-trapezoid-v1.0"
expected_shanghai_preprocessing_version <- "training_pointwise_function_centering_v2"
expected_schema_version <- "2.0"

fail <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

read_required_csv <- function(path, columns) {
  if (!file.exists(path)) {
    fail("Required result file is missing: %s", path)
  }
  result <- read.csv(path, check.names = FALSE)
  missing_columns <- setdiff(columns, names(result))
  if (length(missing_columns) > 0L) {
    fail("%s is missing columns: %s", path, paste(missing_columns, collapse = ", "))
  }
  result
}

require_unique <- function(data, keys, label) {
  duplicated_rows <- duplicated(data[keys])
  if (any(duplicated_rows)) {
    fail("%s contains %d duplicate rows by %s", label, sum(duplicated_rows), paste(keys, collapse = ", "))
  }
}

require_complete <- function(data, label, require_strict = TRUE) {
  if (any(is.na(data$status)) || any(data$status != "ok")) {
    fail("%s contains failed or missing status values", label)
  }
  if (require_strict && (!"strict_converged" %in% names(data) || any(!data$strict_converged))) {
    fail("%s contains non-strictly-converged fits", label)
  }
}

v2_validation_environment <- new.env(parent = globalenv())
source(
  "scripts/validation/validate_formal_results.R",
  local = v2_validation_environment
)

fair_baseline_validation_environment <- new.env(parent = globalenv())
source(
  "scripts/validation/validate_fair_baseline_results.R",
  local = fair_baseline_validation_environment
)

fair_fssgl_validation_environment <- new.env(parent = globalenv())
source(
  "scripts/validation/validate_fair_fssgl_results.R",
  local = fair_fssgl_validation_environment
)

shanghai_table_dir <- "results/tables/shanghai_metroflow/orthonormal"
shanghai_data_dir <- "data/processed/shanghai_metroflow/orthonormal"

real_main <- read_required_csv(
  file.path(shanghai_data_dir, "multi_targets", "multi_targets_summary.csv"),
  c("target", "status", "strict_converged", "selected_station_count", "test_curve_rmse")
)
if (nrow(real_main) != 12L) {
  fail("Shanghai main experiment must contain 12 targets")
}
require_unique(real_main, "target", "Shanghai main experiment")
require_complete(real_main, "Shanghai main experiment")

real_ablation <- read_required_csv(
  file.path(shanghai_data_dir, "ablation", "ablation_summary.csv"),
  c("ablation_id", "target", "status", "strict_converged", "selected_station_count", "test_curve_rmse")
)
if (nrow(real_ablation) != 24L || length(unique(real_ablation$ablation_id)) != 2L) {
  fail("Shanghai ablation must contain two variants over 12 targets")
}
require_unique(real_ablation, c("ablation_id", "target"), "Shanghai ablation")
require_complete(real_ablation, "Shanghai ablation")

compact_real <- read_required_csv(
  file.path(shanghai_table_dir, "multi_target_summary_compact.csv"),
  c("target", "selected_station_count", "test_curve_rmse", "final_residual_variance")
)
if (!identical(sort(compact_real$target), sort(real_main$target))) {
  fail("Compact Shanghai table does not match the detailed target set")
}

manifest_paths <- c(
  file.path(shanghai_data_dir, "multi_targets", "multi_targets_manifest.rds"),
  file.path(shanghai_data_dir, "ablation", "ablation_manifest.rds")
)
for (path in manifest_paths) {
  if (!file.exists(path)) {
    fail("Required experiment manifest is missing: %s", path)
  }
  manifest <- readRDS(path)
  if (!identical(manifest$algorithm_version, expected_shanghai_algorithm_version) ||
      !identical(manifest$result_schema_version, expected_schema_version) ||
      !identical(manifest$design$basis_geometry, expected_basis_geometry_version) ||
      !identical(
        manifest$design$preprocessing_version,
        expected_shanghai_preprocessing_version
      )) {
    fail("Manifest has an unexpected algorithm, schema, basis geometry, or preprocessing version: %s", path)
  }
}

running_files <- list.files(
  "data/processed",
  pattern = "_running(_shard[0-9]+of[0-9]+)?[.]rds$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(running_files) > 0L) {
  fail("Unfinished checkpoints remain in the active data tree: %s", paste(running_files, collapse = ", "))
}

cat("Current output validation passed.\n")
cat("  Formal simulation suite: validated by validate_formal_results.R\n")
cat("  Shanghai: 12 main fits and 24 ablation fits, all strictly converged\n")
source("scripts/validation/validate_shanghai_resampling.R", local = new.env(parent = globalenv()))
