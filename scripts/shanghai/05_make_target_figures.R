# Summarize the broader multi-target Shanghai experiment.

library(data.table)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
processed_dir <- file.path(root, "data/processed/shanghai_metroflow")
multi_dir <- file.path(processed_dir, "orthonormal/multi_targets")
figure_dir <- file.path(root, "results/figures/shanghai_metroflow/orthonormal")
table_dir <- file.path(root, "results/tables/shanghai_metroflow/orthonormal")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(multi_dir, "multi_targets_summary.csv")
if (!file.exists(summary_file)) {
  stop("Missing multi-target summary. Run scripts/shanghai/04_fit_targets.R first.")
}

summary <- fread(summary_file)
summary <- summary[status == "ok"]
summary[, target_short := target]
summary[nchar(target_short) > 22, target_short := paste0(substr(target_short, 1, 20), "...")]
summary[, target_short := factor(target_short, levels = target_short[order(test_curve_rmse)])]

out_file <- file.path(figure_dir, "multi_target_summary.png")
png(out_file, width = 1900, height = 1500, res = 160)
old_par <- par(no.readonly = TRUE)

par(mfrow = c(2, 2), mar = c(10, 4.5, 3, 1))

ordered <- summary[order(test_curve_rmse)]
barplot(
  ordered$test_curve_rmse,
  names.arg = ordered$target_short,
  cex.names = 0.72,
  las = 2,
  col = "#1B6CA8",
  ylab = "Test curve RMSE",
  main = "Held-out Curve Error"
)
grid(nx = NA, ny = NULL, col = "#DDDDDD")

ordered <- summary[order(selected_line_count, selected_station_count)]
barplot(
  ordered$selected_line_count,
  names.arg = ordered$target_short,
  cex.names = 0.72,
  las = 2,
  col = "#C73E1D",
  ylab = "Selected line count",
  main = "Group Selection"
)
grid(nx = NA, ny = NULL, col = "#DDDDDD")

barplot(
  ordered$selected_station_count,
  names.arg = ordered$target_short,
  cex.names = 0.72,
  las = 2,
  col = "#2A9D8F",
  ylab = "Selected station count",
  main = "Covariate Selection"
)
grid(nx = NA, ny = NULL, col = "#DDDDDD")

line_tokens <- unlist(strsplit(summary[nchar(selected_lines) > 0, selected_lines], ";", fixed = TRUE))
line_freq <- sort(table(line_tokens), decreasing = TRUE)
barplot(
  as.numeric(line_freq),
  names.arg = names(line_freq),
  cex.names = 0.78,
  las = 2,
  col = "#7B2CBF",
  ylab = "Frequency across targets",
  main = "Selected Line Frequency"
)
grid(nx = NA, ny = NULL, col = "#DDDDDD")

par(old_par)
dev.off()

fwrite(
  summary[, .(
    target,
    relaxed_converged,
    selected_line_count,
    selected_station_count,
    selected_lines,
    top_line,
    top_line_prob,
    top_station,
    top_station_prob,
    test_curve_rmse,
    final_residual_variance
  )],
  file.path(table_dir, "multi_target_summary_compact.csv")
)

cat("Saved multi-target summary:", out_file, "\n")
