# Create figures and compact tables for the Shanghai ablation experiment.

library(data.table)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
ablation_dir <- file.path(root, "data/processed/shanghai_metroflow/orthonormal/ablation")
figure_dir <- file.path(root, "results/figures/shanghai_metroflow/orthonormal")
table_dir <- file.path(root, "results/tables/shanghai_metroflow/orthonormal")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(ablation_dir, "ablation_summary.csv")
variant_file <- file.path(ablation_dir, "ablation_variant_summary.csv")
if (!file.exists(summary_file) || !file.exists(variant_file)) {
  stop("Missing ablation outputs. Run scripts/shanghai/06_run_structure_ablation.R first.")
}

dt <- fread(summary_file)
variant <- fread(variant_file)
dt <- dt[status == "ok"]

variant[, label := fifelse(
  use_group_penalty,
  "Covariate + group",
  "Covariate only"
)]
variant[, label := factor(label, levels = c(
  "Covariate only",
  "Covariate + group"
))]
setorder(variant, label)

plot_variant_summary <- function() {
  out_file <- file.path(figure_dir, "ablation_variant_summary.png")
  png(out_file, width = 1600, height = 1100, res = 160)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  par(mfrow = c(1, 3), mar = c(6, 4.5, 3, 1))
  cols <- c("#1B6CA8", "#C73E1D")

  barplot(
    variant$mean_test_curve_rmse,
    names.arg = variant$label,
    las = 2,
    col = cols,
    ylab = "Mean test curve RMSE",
    main = "Prediction Error"
  )
  grid(nx = NA, ny = NULL, col = "#DDDDDD")

  barplot(
    variant$mean_station_count,
    names.arg = variant$label,
    las = 2,
    col = cols,
    ylab = "Mean selected stations",
    main = "Covariate Sparsity"
  )
  grid(nx = NA, ny = NULL, col = "#DDDDDD")

  line_count <- variant$mean_line_count
  line_count[!is.finite(line_count)] <- 0
  barplot(
    line_count,
    names.arg = variant$label,
    las = 2,
    col = cols,
    ylab = "Mean selected lines",
    main = "Group Selection"
  )
  grid(nx = NA, ny = NULL, col = "#DDDDDD")

  out_file
}

plot_target_rmse_tradeoff <- function() {
  wide <- dcast(
    dt,
    target ~ ablation_id,
    value.var = c("test_curve_rmse", "selected_station_count")
  )
  wide[, rmse_delta_group_vs_covariate :=
    test_curve_rmse_B_covariate_group - test_curve_rmse_A_covariate_only]
  wide[, station_delta_group_vs_covariate :=
    selected_station_count_B_covariate_group - selected_station_count_A_covariate_only]
  setorder(wide, rmse_delta_group_vs_covariate)

  fwrite(wide, file.path(table_dir, "ablation_target_level_wide.csv"))

  out_file <- file.path(figure_dir, "ablation_target_tradeoff.png")
  png(out_file, width = 1500, height = 900, res = 160)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  par(mfrow = c(1, 2), mar = c(7, 4.5, 3, 1))
  barplot(
    wide$rmse_delta_group_vs_covariate,
    names.arg = wide$target,
    las = 2,
    col = ifelse(wide$rmse_delta_group_vs_covariate <= 0, "#2A9D8F", "#C73E1D"),
    ylab = "RMSE(covariate+group) - RMSE(covariate-only)",
    main = "Prediction Cost of Group Structure"
  )
  abline(h = 0, col = "#222222", lwd = 1)
  grid(nx = NA, ny = NULL, col = "#DDDDDD")

  barplot(
    wide$station_delta_group_vs_covariate,
    names.arg = wide$target,
    las = 2,
    col = ifelse(wide$station_delta_group_vs_covariate <= 0, "#1B6CA8", "#C73E1D"),
    ylab = "Station count(covariate+group) - count(covariate-only)",
    main = "Sparsity Gain of Group Structure"
  )
  abline(h = 0, col = "#222222", lwd = 1)
  grid(nx = NA, ny = NULL, col = "#DDDDDD")

  out_file
}

compact <- variant[, .(
  ablation_id,
  use_group_penalty,
  all_relaxed_converged,
  mean_station_count,
  mean_line_count,
  mean_diagnostic_line_count,
  mean_test_curve_rmse,
  max_test_curve_rmse,
  mean_test_curve_mae
)]
fwrite(compact, file.path(table_dir, "ablation_variant_summary_compact.csv"))

variant_png <- plot_variant_summary()
tradeoff_png <- plot_target_rmse_tradeoff()

cat("Saved ablation variant summary:", variant_png, "\n")
cat("Saved ablation target tradeoff:", tradeoff_png, "\n")
