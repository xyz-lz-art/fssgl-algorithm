# Recompute resampling summaries from replicate evidence without refitting.
library(data.table)
source("R/parameters.R")
source("R/application/shanghai_workflow.R")

data_root <- "data/processed/shanghai_metroflow/orthonormal"
table_root <- "results/tables/shanghai_metroflow/orthonormal"
read_evidence <- function(directory, filename, report = FALSE) {
  path <- file.path(if (report) table_root else data_root, directory, filename)
  if (!file.exists(path)) stop("Missing resampling evidence: ", path)
  fread(path)
}
same_numbers <- function(a, b) {
  stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-10)))
}
unique_rows <- function(x, columns) stopifnot(!anyDuplicated(x[, ..columns]))

s <- read_evidence("stability_selection", "stability_fit_summary.csv")
r <- read_evidence("rolling_evaluation", "rolling_fit_summary.csv")
for (x in list(s, r)) {
  stopifnot(setequal(x$target, shanghai_target_names()),
            all(x$status == "ok"), all(x$strict_converged),
            all(is.finite(x$final_beta_change)))
}
unique_rows(s, c("target", "replicate_id"))
unique_rows(r, c("target", "window_id"))
stopifnot(nrow(s) == 240L, nrow(r) == 108L,
          all(s[, .N, by = target]$N == 20L),
          all(r[, .N, by = target]$N == 9L))

sd <- read_evidence("stability_selection", "stability_subsample_design.csv")
unique_rows(sd, c("replicate_id", "day_pos"))
stopifnot(nrow(sd) == 980L, all(sd$day_pos %in% 1:98),
          all(sd[, .N, by = replicate_id]$N == 49L),
          all(sd[, .N, by = .(replicate_id, block_id)]$N == 7L),
          all(sd$block_id == (sd$day_pos - 1L) %/% 7L + 1L),
          all(s$n_train_days == 49L))

st <- read_evidence("stability_selection", "stability_target_summary.csv", TRUE)
for (unit in c("station", "line")) {
  rep <- read_evidence("stability_selection", paste0("stability_", unit, "_replicates.csv"))
  freq <- read_evidence("stability_selection", paste0("stability_", unit, "_frequency.csv"), TRUE)
  keys <- c("target", unit)
  unique_rows(rep, c(keys, "replicate_id"))
  unique_rows(freq, keys)
  stopifnot(all(rep$posterior_slab_prob >= 0 & rep$posterior_slab_prob <= 1))
  if (unit == "station") {
    stopifnot(nrow(rep) == 240L * 301L,
              all(rep$selected == (rep$posterior_slab_prob >= 0.5)))
  }
  calculated <- rep[, .(frequency = mean(selected), n = .N), by = keys]
  paired <- merge(calculated, freq, by = keys)
  stopifnot(nrow(paired) == nrow(freq), all(paired$n == 20L),
            all(paired$stable_selected == (paired$selection_frequency >= 0.8)))
  same_numbers(paired$frequency, paired$selection_frequency)
  counts <- rep[, .(selected_count = sum(selected)), by = .(target, replicate_id)]
  count_check <- merge(counts, s, by = c("target", "replicate_id"))
  same_numbers(count_check$selected_count, count_check[[paste0("selected_", unit, "_count")]])
  stable <- freq[, .(stable_count = sum(stable_selected)), by = target]
  stable_check <- merge(stable, st, by = "target")
  same_numbers(stable_check$stable_count, stable_check[[paste0("stable_", unit, "_count")]])
}

rd <- read_evidence("rolling_evaluation", "rolling_window_design.csv")
stopifnot(identical(rd$train_end_day, seq(60L, 116L, by = 7L)),
          all(rd$train_end_day < rd$test_start_day),
          all(rd$n_test_days == 7L),
          all(r$train_end_date < r$test_start_date),
          all(r$n_test_days == 7L))
stopifnot(identical(unlist(Map(seq.int, rd$test_start_day, rd$test_end_day)), 61:123))
dates <- fread("data/processed/shanghai_metroflow/date_index.csv")$date
joined <- merge(r, rd, by = "window_id", suffixes = c("", "_design"))
stopifnot(all(joined$n_train_days == joined$n_train_days_design),
          all(joined$train_end_date == dates[joined$train_end_day]),
          all(joined$test_start_date == dates[joined$test_start_day]),
          all(joined$test_end_date == dates[joined$test_end_day]))
same_numbers(r$test_curve_rmse_original,
             r$test_curve_rmse_standardized * r$training_response_sd)
same_numbers(r$test_curve_mae_original,
             r$test_curve_mae_standardized * r$training_response_sd)

for (key in c("target", "window_id")) {
  report <- read_evidence("rolling_evaluation",
                         if (key == "target") "rolling_target_summary.csv" else "rolling_window_summary.csv", TRUE)
  calculated <- r[, .(rmse = mean(test_curve_rmse_standardized),
                      mae = mean(test_curve_mae_standardized)), by = key]
  paired <- merge(calculated, report, by = key)
  stopifnot(nrow(paired) == nrow(report))
  same_numbers(paired$rmse, paired$mean_test_curve_rmse_standardized)
  same_numbers(paired$mae, paired$mean_test_curve_mae_standardized)
}
overall <- read_evidence("rolling_evaluation", "rolling_overall_summary.csv", TRUE)
same_numbers(mean(r$test_curve_rmse_standardized), overall$mean_test_curve_rmse_standardized)

for (unit in c("station", "line")) {
  selected <- read_evidence("rolling_evaluation", paste0("rolling_selected_", unit, "s.csv"))
  unique_rows(selected, c("target", "window_id", unit))
  counts <- selected[, .(count = .N), by = .(target, window_id)]
  paired <- merge(r, counts, by = c("target", "window_id"), all.x = TRUE)
  paired[is.na(count), count := 0L]
  same_numbers(paired$count, paired[[paste0("selected_", unit, "_count")]])
}

for (entry in list(c("stability_selection", "stability"), c("rolling_evaluation", "rolling"))) {
  m <- readRDS(file.path(data_root, entry[1], paste0(entry[2], "_manifest.rds")))
  stopifnot(identical(m$algorithm_version, FSSGL_WEIGHTED_MEMBERSHIP_VERSION),
            identical(m$design$basis_geometry, FSSGL_BASIS_GEOMETRY_VERSION),
            identical(m$design$preprocessing_version, SHANGHAI_PREPROCESSING_VERSION),
            identical(m$parameters, shanghai_main_parameters()),
            !file.exists(file.path(data_root, entry[1], paste0(entry[2], "_running.rds"))))
}

identifiability <- fread(file.path(table_root, "shanghai_identifiability_by_target.csv"))
identifiability_summary <- fread(file.path(table_root, "shanghai_identifiability_summary.csv"))
correlation_pairs <- fread(file.path(table_root, "shanghai_high_correlation_pairs.csv"))
stopifnot(
  nrow(identifiability) == 12L,
  setequal(identifiability$target, shanghai_target_names()),
  all(identifiability$n_rows == 98L),
  all(identifiability$n_columns == 1505L),
  all(identifiability$numerical_rank == 97L),
  all(is.infinite(identifiability$gram_condition)),
  all(identifiability$minimum_block_rank == 5L),
  all(identifiability$design_rank_deficient),
  all(identifiability$severely_ill_conditioned),
  nrow(identifiability_summary) == 1L,
  identifiability_summary$targets == 12L,
  identifiability_summary$joint_rank_min == 97L,
  identifiability_summary$joint_rank_max == 97L,
  identifiability_summary$design_rank_deficient_targets == 12L,
  identifiability_summary$severely_ill_conditioned_targets == 12L,
  nrow(correlation_pairs) > 0L,
  max(correlation_pairs$absolute_correlation) > 0.99
)
cat(paste0(
  "Shanghai validation passed: 240 stability fits, 108 rolling fits, and ",
  "12 target-specific identifiability diagnostics; summaries recomputed from evidence.\n"
))
