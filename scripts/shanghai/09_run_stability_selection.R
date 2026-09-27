# Block-subsample stability selection for the 12 Shanghai targets.
# Only the original 98-day development period is used. Each replicate samples
# seven of fourteen non-overlapping seven-day blocks (49 training days).

library(data.table)

source("R/parameters.R")
source("R/fssgl/basis_design.R")
source("R/fssgl/penalty_weights.R")
source("R/fssgl/weighted_membership_solver.R")
source("R/application/shanghai_workflow.R")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
processed_dir <- file.path(root, "data/processed/shanghai_metroflow")
metadata_dir <- file.path(root, "data/metadata/shanghai_metroflow")
out_dir <- file.path(processed_dir, "orthonormal/stability_selection")
table_dir <- file.path(root, "results/tables/shanghai_metroflow/orthonormal/stability_selection")
figure_dir <- file.path(root, "results/figures/shanghai_metroflow/orthonormal/stability_selection")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

targets <- shanghai_target_names()
params <- shanghai_main_parameters()
stability_seed <- 2026090201L
n_replicates <- 20L
block_length <- 7L
development_days <- 1:98
heldout_days <- 99:123
selection_threshold <- 0.80

blocks <- split(development_days, rep(seq_len(14L), each = block_length))
set.seed(stability_seed)
sampled_block_keys <- character()
sampled_blocks <- vector("list", n_replicates)
for (replicate_id in seq_len(n_replicates)) {
  repeat {
    candidate <- sort(sample(seq_along(blocks), length(blocks) / 2L, replace = FALSE))
    key <- paste(candidate, collapse = ";")
    if (!key %in% sampled_block_keys) break
  }
  sampled_block_keys <- c(sampled_block_keys, key)
  sampled_blocks[[replicate_id]] <- candidate
}

stability_design <- rbindlist(lapply(seq_len(n_replicates), function(replicate_id) {
  block_ids <- sampled_blocks[[replicate_id]]
  days <- sort(unlist(blocks[block_ids], use.names = FALSE))
  data.table(
    replicate_id = replicate_id,
    block_id = rep(block_ids, each = block_length),
    day_pos = days
  )
}))
fwrite(stability_design, file.path(out_dir, "stability_subsample_design.csv"))

checkpoint_file <- file.path(out_dir, "stability_running.rds")
completed_file <- file.path(out_dir, "stability_completed.rds")
empty_state <- list(
  preprocessing_version = SHANGHAI_PREPROCESSING_VERSION,
  fit_summary = data.table(),
  station_replicates = data.table(),
  line_replicates = data.table()
)
state_candidate <- if (file.exists(checkpoint_file)) {
  readRDS(checkpoint_file)
} else if (file.exists(completed_file)) {
  readRDS(completed_file)
} else {
  empty_state
}
state <- if (
  is.list(state_candidate) &&
    identical(state_candidate$preprocessing_version, SHANGHAI_PREPROCESSING_VERSION)
) state_candidate else empty_state

shanghai_data <- load_shanghai_data(processed_dir, metadata_dir)

save_checkpoint <- function() saveRDS(state, checkpoint_file)

for (target_name in targets) {
  for (replicate_id in seq_len(n_replicates)) {
    current_replicate <- replicate_id
    already_done <- nrow(state$fit_summary) > 0L && nrow(state$fit_summary[
      target == target_name & replicate_id == current_replicate
    ]) > 0L
    if (already_done) {
      cat("Skipping completed stability fit", target_name, replicate_id, "\n")
      next
    }

    train_days <- stability_design[replicate_id == current_replicate, day_pos]
    cat("Stability fit", target_name, "replicate", replicate_id, "\n")
    result <- tryCatch(
      fit_shanghai_target_fssgl(
        target_name = target_name,
        params = params,
        processed_dir = processed_dir,
        metadata_dir = metadata_dir,
        kx = 5,
        ky = 5,
        train_days = train_days,
        test_days = heldout_days,
        shanghai_data = shanghai_data
      ),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      state$fit_summary <- rbind(
        state$fit_summary,
        data.table(
          target = target_name,
          replicate_id = replicate_id,
          status = "error",
          error_message = conditionMessage(result)
        ),
        fill = TRUE
      )
      save_checkpoint()
      next
    }

    fit <- result$fit
    final <- fit$history[.N]
    state$fit_summary <- rbind(
      state$fit_summary,
      data.table(
        target = target_name,
        replicate_id = replicate_id,
        status = "ok",
        n_train_days = length(train_days),
        first_train_day = min(train_days),
        last_train_day = max(train_days),
        final_iter = final$iter,
        final_beta_change = final$beta_change,
        objective_tail_rel_change = fit$convergence$objective_tail_rel_change,
        strict_converged = fit$convergence$strict_converged,
        relaxed_converged = fit$convergence$relaxed_converged,
        selected_station_count = final$selected_station_count,
        selected_line_count = final$selected_line_count,
        final_residual_variance = fit$convergence$final_residual_variance
      ),
      fill = TRUE
    )

    station_rows <- fit$station_posterior[, .(
      target = target_name,
      replicate_id = replicate_id,
      station,
      name,
      lines,
      posterior_slab_prob,
      beta_norm,
      selected
    )]
    line_rows <- fit$line_posterior[, .(
      target = target_name,
      replicate_id = replicate_id,
      line,
      posterior_slab_prob,
      line_norm,
      selected
    )]
    state$station_replicates <- rbind(state$station_replicates, station_rows, fill = TRUE)
    state$line_replicates <- rbind(state$line_replicates, line_rows, fill = TRUE)
    save_checkpoint()
  }
}

fit_summary <- state$fit_summary
station_replicates <- state$station_replicates
line_replicates <- state$line_replicates

station_frequency <- station_replicates[, .(
  n_successful_fits = .N,
  selection_frequency = mean(selected),
  mean_slab_probability = mean(posterior_slab_prob),
  median_slab_probability = median(posterior_slab_prob),
  mean_beta_norm = mean(beta_norm)
), by = .(target, station, name, lines)]
station_frequency[, stable_selected := selection_frequency >= selection_threshold]
station_frequency[, frequency_rank := frank(-selection_frequency, ties.method = "min"), by = target]
setorder(station_frequency, target, -selection_frequency, -mean_slab_probability)

line_frequency <- line_replicates[, .(
  n_successful_fits = .N,
  selection_frequency = mean(selected),
  mean_slab_probability = mean(posterior_slab_prob),
  median_slab_probability = median(posterior_slab_prob),
  mean_line_norm = mean(line_norm)
), by = .(target, line)]
line_frequency[, stable_selected := selection_frequency >= selection_threshold]
line_frequency[, frequency_rank := frank(-selection_frequency, ties.method = "min"), by = target]
setorder(line_frequency, target, -selection_frequency, -mean_slab_probability)

target_summary <- fit_summary[, .(
  n_requested_fits = n_replicates,
  n_successful_fits = sum(status == "ok"),
  strict_convergence_rate = mean(strict_converged[status == "ok"]),
  relaxed_convergence_rate = mean(relaxed_converged[status == "ok"]),
  mean_selected_stations = mean(selected_station_count[status == "ok"]),
  sd_selected_stations = sd(selected_station_count[status == "ok"]),
  mean_selected_lines = mean(selected_line_count[status == "ok"]),
  mean_iterations = mean(final_iter[status == "ok"])
), by = target]
stable_station_counts <- station_frequency[, .(
  stable_station_count = sum(stable_selected),
  stations_frequency_at_least_half = sum(selection_frequency >= 0.50),
  top_station = name[1L],
  top_station_frequency = selection_frequency[1L]
), by = target]
stable_line_counts <- line_frequency[, .(
  stable_line_count = sum(stable_selected),
  lines_frequency_at_least_half = sum(selection_frequency >= 0.50),
  top_line = line[1L],
  top_line_frequency = selection_frequency[1L]
), by = target]
target_summary <- merge(target_summary, stable_station_counts, by = "target")
target_summary <- merge(target_summary, stable_line_counts, by = "target")
setorder(target_summary, target)

fwrite(fit_summary, file.path(out_dir, "stability_fit_summary.csv"))
fwrite(station_replicates, file.path(out_dir, "stability_station_replicates.csv"))
fwrite(line_replicates, file.path(out_dir, "stability_line_replicates.csv"))
fwrite(station_frequency, file.path(table_dir, "stability_station_frequency.csv"))
fwrite(line_frequency, file.path(table_dir, "stability_line_frequency.csv"))
fwrite(target_summary, file.path(table_dir, "stability_target_summary.csv"))

png(file.path(figure_dir, "stability_target_summary.png"), width = 2000, height = 1100, res = 160)
old_par <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(5, 14, 4, 2))
ordered <- target_summary[order(stable_station_count)]
barplot(ordered$stable_station_count, names.arg = ordered$target, las = 1,
        horiz = TRUE, cex.names = 0.8, col = "#2A9D8F",
        xlab = "Stations with frequency >= 0.80", main = "Stable Stations (20 subsamples)")
barplot(ordered$stable_line_count, names.arg = ordered$target, las = 1,
        horiz = TRUE, cex.names = 0.8, col = "#7B2CBF",
        xlab = "Lines with frequency >= 0.80", main = "Stable Lines (20 subsamples)")
par(old_par)
dev.off()

saveRDS(state, completed_file)
if (file.exists(checkpoint_file)) unlink(checkpoint_file)
saveRDS(
  new_experiment_manifest(
    experiment_id = "shanghai_block_stability_selection",
    parameters = params,
    design = list(
      targets = targets,
      n_replicates = n_replicates,
      seed = stability_seed,
      development_days = development_days,
      heldout_days_excluded_from_selection = heldout_days,
      block_length = block_length,
      blocks_per_subsample = length(blocks) / 2L,
      days_per_subsample = length(development_days) / 2L,
      selection_threshold = selection_threshold,
      kx = 5L,
      ky = 5L,
      basis_geometry = FSSGL_BASIS_GEOMETRY_VERSION,
      preprocessing_version = SHANGHAI_PREPROCESSING_VERSION
    )
  ),
  file.path(out_dir, "stability_manifest.rds")
)

cat("\nStability-selection target summary\n")
print(target_summary)
