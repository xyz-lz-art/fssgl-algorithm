# Finite-sample identifiability diagnostics for the Shanghai coefficient design.

library(data.table)

source("R/fssgl/basis_design.R")
source("R/application/shanghai_workflow.R")
source("R/parameters.R")
source("R/diagnostics/identifiability.R")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
processed_dir <- file.path(root, "data/processed/shanghai_metroflow")
metadata_dir <- file.path(root, "data/metadata/shanghai_metroflow")
table_dir <- file.path(root, "results/tables/shanghai_metroflow/orthonormal")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

shanghai_data <- load_shanghai_data(processed_dir, metadata_dir)
diagnostic_rows <- list()
pair_rows <- list()
for (target in shanghai_target_names()) {
  cat("Identifiability diagnostic:", target, "\n")
  prepared <- prepare_shanghai_target(
    target, processed_dir, metadata_dir, shanghai_data = shanghai_data
  )
  design <- build_shanghai_basis_design(prepared, kx = 5L, ky = 5L)
  diag <- functional_predictor_identifiability(design$x_coef)
  diagnostic_rows[[target]] <- as.data.table(c(list(target = target), diag))
  pairs <- as.data.table(top_functional_profile_correlations(
    design$x_coef, labels = prepared$predictor_station$name, top_n = 10L
  ))
  pairs[, target := ..target]
  pair_rows[[target]] <- pairs
}

diagnostics <- rbindlist(diagnostic_rows, fill = TRUE)
pairs <- rbindlist(pair_rows, fill = TRUE)
summary <- diagnostics[, .(
  targets = .N,
  predictor_count = unique(predictor_count),
  coefficient_columns = unique(n_columns),
  joint_rank_min = min(numerical_rank),
  joint_rank_max = max(numerical_rank),
  effective_rank_mean = mean(entropy_effective_rank),
  effective_rank_range = sprintf("%.2f--%.2f", min(entropy_effective_rank), max(entropy_effective_rank)),
  nonzero_gram_condition_median = median(nonzero_gram_condition),
  minimum_block_rank = min(minimum_block_rank),
  minimum_block_effective_rank = min(minimum_block_effective_rank),
  maximum_block_gram_condition = max(maximum_block_gram_condition),
  maximum_absolute_profile_correlation = max(maximum_absolute_profile_correlation),
  q95_absolute_profile_correlation_mean = mean(q95_absolute_profile_correlation),
  pair_fraction_above_0_9_mean = mean(pair_fraction_above_0_9),
  design_rank_deficient_targets = sum(design_rank_deficient),
  severely_ill_conditioned_targets = sum(severely_ill_conditioned)
)]

fwrite(diagnostics, file.path(table_dir, "shanghai_identifiability_by_target.csv"))
fwrite(pairs, file.path(table_dir, "shanghai_high_correlation_pairs.csv"))
fwrite(summary, file.path(table_dir, "shanghai_identifiability_summary.csv"))
print(summary)
print(pairs[order(-absolute_correlation)][1:20])
