# Build manuscript tables and figures from the frozen replicate-level data.

library(data.table)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "results/tables/simulation/v2_formal")
figure_dir <- file.path(root, "results/figures/simulation/v2_formal")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

main <- fread(file.path(table_dir, "fssgl_v2_main_replicates.csv"))
robust <- fread(file.path(table_dir, "fssgl_v2_robustness_replicates.csv"))
ablation <- fread(file.path(table_dir, "fssgl_v2_structure_ablation_replicates.csv"))
comparison <- fread(file.path(table_dir, "functional_method_comparison_v2_replicates.csv"))

metrics <- c(
  "covariate_tpr", "covariate_fpr", "covariate_fdr",
  "selected_covariate_count", "group_tpr", "group_fpr", "group_fdr",
  "selected_group_count", "coefficient_relative_error",
  "active_coefficient_relative_error", "inactive_coefficient_norm",
  "test_coeff_rmse", "runtime_sec", "final_iter"
)

summarize_metrics <- function(data, by, experiment) {
  data <- copy(data)
  measure_columns <- intersect(metrics, names(data))
  for (column in measure_columns) set(data, j = column, value = as.numeric(data[[column]]))
  long <- melt(
    data,
    id.vars = by,
    measure.vars = measure_columns,
    variable.name = "metric",
    value.name = "value"
  )
  out <- long[!is.na(value), .(
    n = .N,
    mean = mean(value),
    sd = sd(value),
    mcse = sd(value) / sqrt(.N)
  ), by = c(by, "metric")]
  out[, experiment := experiment]
  setcolorder(out, c("experiment", by, "metric", "n", "mean", "sd", "mcse"))
  out[]
}

main_summary <- summarize_metrics(main, "p", "main")
robust_summary <- summarize_metrics(robust, "scenario_id", "robustness")
ablation_summary <- summarize_metrics(ablation, "variant", "structure_ablation")
method_summary <- summarize_metrics(comparison, c("method_id", "p"), "method_comparison")

fwrite(main_summary, file.path(table_dir, "fssgl_v2_main_summary_mcse.csv"))
fwrite(robust_summary, file.path(table_dir, "fssgl_v2_robustness_summary_mcse.csv"))
fwrite(ablation_summary, file.path(table_dir, "fssgl_v2_structure_ablation_summary_mcse.csv"))
fwrite(method_summary, file.path(table_dir, "functional_method_comparison_v2_summary_mcse.csv"))
fwrite(
  rbindlist(list(main_summary, robust_summary, ablation_summary, method_summary), fill = TRUE),
  file.path(table_dir, "fssgl_v2_manuscript_summary_mcse.csv")
)

method_labels <- c(
  fssgl_v2 = "FSSGL",
  basis_ridge = "Basis ridge",
  fpca_ridge = "FPCA ridge",
  fpca_group_lasso = "FPCA group lasso",
  fpca_group_scad = "FPCA group SCAD",
  fpca_aenet = "FPCA-AEnet",
  structured_group_lasso = "Structured group lasso",
  kernel_ridge = "Kernel ridge"
)
sparse_methods <- c(
  "fssgl_v2", "fpca_group_lasso", "fpca_group_scad", "fpca_aenet",
  "structured_group_lasso"
)
all_methods <- names(method_labels)

variant_order <- c("correct_groups", "covariate_only", "permuted_groups")

scenario_order <- c("baseline", "snr_low", "correlation_high", "within_group_sparse")
scenario_labels <- c(
  baseline = "Baseline",
  snr_low = "Low SNR",
  correlation_high = "High correlation",
  within_group_sparse = "One active per group"
)
palette <- c(
  fssgl_v2 = "#176B67", basis_ridge = "#777777", fpca_ridge = "#3567A8",
  fpca_group_lasso = "#B64E3C", fpca_aenet = "#A26D13",
  fpca_group_scad = "#8B3A62",
  structured_group_lasso = "#2F855A", kernel_ridge = "#7A5A9E"
)

draw_grouped_metric <- function(summary, metric_name, methods, title, ylab, ylim = NULL) {
  values <- summary[metric == metric_name & method_id %in% methods]
  values[, method_id := factor(method_id, levels = methods)]
  setorder(values, method_id, p)
  x_base <- match(values$p, c(10L, 20L))
  method_index <- match(as.character(values$method_id), methods)
  offsets <- seq(-0.27, 0.27, length.out = length(methods))
  x <- x_base + offsets[method_index]
  half_width <- 1.96 * values$mcse
  if (is.null(ylim)) {
    limits <- range(c(values$mean - half_width, values$mean + half_width), na.rm = TRUE)
    padding <- max(diff(limits) * 0.18, 0.015)
    ylim <- limits + c(-padding, padding)
  }
  plot(x, values$mean, type = "n", xaxt = "n", xlab = "", ylab = ylab,
       xlim = c(0.55, 2.45), ylim = ylim, main = title)
  axis(1, at = 1:2, labels = c("p = 10", "p = 20"))
  grid(nx = NA, ny = NULL, col = "#E5E5E5")
  for (i in seq_len(nrow(values))) {
    method <- as.character(values$method_id[i])
    lower <- max(values$mean[i] - half_width[i], ylim[1L])
    upper <- min(values$mean[i] + half_width[i], ylim[2L])
    if (upper - lower > 1e-12) {
      arrows(x[i], lower, x[i], upper, angle = 90, code = 3, length = 0.025,
             col = palette[[method]], lwd = 1.1)
    }
    points(x[i], values$mean[i], pch = 14 + method_index[i], col = palette[[method]], cex = 1.05)
  }
  legend("topleft", method_labels[methods], col = palette[methods],
         pch = 14 + seq_along(methods), bty = "n", cex = 0.64)
}

main_figure <- file.path(figure_dir, "fssgl_v2_main_baseline.pdf")
pdf(main_figure, width = 8.2, height = 6.5, family = "Helvetica")
par(mfrow = c(2, 2), mar = c(3.8, 4.2, 2.3, 0.7), las = 1)
draw_grouped_metric(method_summary, "covariate_tpr", sparse_methods,
                    "(a) Support recovery", "Covariate TPR", c(0, 1.05))
draw_grouped_metric(method_summary, "covariate_fdr", sparse_methods,
                    "(b) False discoveries", "Covariate FDR", c(0, 0.5))
draw_grouped_metric(
                    method_summary, "coefficient_relative_error",
                    setdiff(all_methods, "kernel_ridge"),
                    "(c) Coefficient estimation", "Relative error")
draw_grouped_metric(method_summary, "test_coeff_rmse", all_methods,
                    "(d) Held-out prediction", "Test coefficient RMSE")
dev.off()

point_ci <- function(x, mean, mcse, color, pch, offset = 0, ylim = NULL) {
  xpos <- x + offset
  lower <- mean - 1.96 * mcse
  upper <- mean + 1.96 * mcse
  if (!is.null(ylim)) {
    lower <- pmax(lower, ylim[1L])
    upper <- pmin(upper, ylim[2L])
  }
  nonzero <- upper - lower > 1e-12
  if (any(nonzero)) {
    arrows(xpos[nonzero], lower[nonzero], xpos[nonzero], upper[nonzero],
           angle = 90, code = 3, length = 0.035, col = color, lwd = 1.2)
  }
  points(xpos, mean, pch = pch, col = color, bg = color, cex = 1.15)
}

metric_rows <- function(summary, metric_name, id_name, order_values) {
  out <- summary[metric == metric_name]
  out[match(order_values, get(id_name))]
}

ablation_figure <- file.path(figure_dir, "fssgl_v2_structure_ablation.pdf")
pdf(ablation_figure, width = 7.4, height = 6.2, family = "Helvetica")
par(mfrow = c(2, 2), mar = c(4.3, 4.1, 2.4, 0.8), las = 1)
x <- seq_along(variant_order)
short_variant_labels <- c("Correct", "No group", "Permuted")
teal <- "#176B67"; red <- "#B64E3C"; blue <- "#3567A8"
tpr <- metric_rows(ablation_summary, "covariate_tpr", "variant", variant_order)
plot(x, tpr$mean, type = "n", xaxt = "n", xlab = "", ylab = "Covariate TPR",
     ylim = c(0, 1.05), xlim = c(0.65, 3.35), main = "(a) Support recovery")
axis(1, at = x, labels = short_variant_labels, cex.axis = 0.82); grid(nx = NA, ny = NULL, col = "#E5E5E5")
point_ci(x, tpr$mean, tpr$mcse, teal, 21, ylim = c(0, 1.05))
fpr <- metric_rows(ablation_summary, "covariate_fpr", "variant", variant_order)
fdr <- metric_rows(ablation_summary, "covariate_fdr", "variant", variant_order)
plot(x, fpr$mean, type = "n", xaxt = "n", xlab = "", ylab = "Rate", ylim = c(0, 0.25),
     xlim = c(0.65, 3.35), main = "(b) False selections")
axis(1, at = x, labels = short_variant_labels, cex.axis = 0.82); grid(nx = NA, ny = NULL, col = "#E5E5E5")
point_ci(x, fpr$mean, fpr$mcse, blue, 21, -0.07, c(0, 0.25)); point_ci(x, fdr$mean, fdr$mcse, red, 22, 0.07, c(0, 0.25))
legend("topleft", c("FPR", "FDR"), col = c(blue, red), pch = c(21, 22), pt.bg = c(blue, red), bty = "n", cex = 0.80)
coef <- metric_rows(ablation_summary, "coefficient_relative_error", "variant", variant_order)
plot(x, coef$mean, type = "n", xaxt = "n", xlab = "", ylab = "Relative error", ylim = c(0.5, 0.76),
     xlim = c(0.65, 3.35), main = "(c) Coefficient estimation")
axis(1, at = x, labels = short_variant_labels, cex.axis = 0.82); grid(nx = NA, ny = NULL, col = "#E5E5E5")
point_ci(x, coef$mean, coef$mcse, teal, 21, ylim = c(0.5, 0.76))
rmse <- metric_rows(ablation_summary, "test_coeff_rmse", "variant", variant_order)
plot(x, rmse$mean, type = "n", xaxt = "n", xlab = "", ylab = "Test coefficient RMSE", ylim = c(0.95, 1.16),
     xlim = c(0.65, 3.35), main = "(d) Held-out prediction")
axis(1, at = x, labels = short_variant_labels, cex.axis = 0.82); grid(nx = NA, ny = NULL, col = "#E5E5E5")
point_ci(x, rmse$mean, rmse$mcse, blue, 21, ylim = c(0.95, 1.16))
dev.off()

draw_robust <- function(first, second = NULL, xlim, xlab, title, legend_labels = NULL) {
  y <- rev(seq_along(scenario_order))
  plot(first$mean, y, type = "n", yaxt = "n", ylab = "", xlab = xlab,
       xlim = xlim, ylim = c(0.5, length(y) + 0.5), main = title)
  axis(2, at = y, labels = unname(scenario_labels[scenario_order]), las = 1, cex.axis = 0.78)
  grid(nx = NULL, ny = NA, col = "#E5E5E5")
  first_lower <- pmax(first$mean - 1.96 * first$mcse, xlim[1L])
  first_upper <- pmin(first$mean + 1.96 * first$mcse, xlim[2L])
  first_nonzero <- first_upper - first_lower > 1e-12
  if (any(first_nonzero)) {
    arrows(first_lower[first_nonzero], y[first_nonzero] + 0.09,
           first_upper[first_nonzero], y[first_nonzero] + 0.09, angle = 90,
           code = 3, length = 0.03, col = teal, lwd = 1.1)
  }
  points(first$mean, y + 0.09, pch = 21, col = teal, bg = teal, cex = 1.05)
  if (!is.null(second)) {
    second_lower <- pmax(second$mean - 1.96 * second$mcse, xlim[1L])
    second_upper <- pmin(second$mean + 1.96 * second$mcse, xlim[2L])
    second_nonzero <- second_upper - second_lower > 1e-12
    if (any(second_nonzero)) {
      arrows(second_lower[second_nonzero], y[second_nonzero] - 0.09,
             second_upper[second_nonzero], y[second_nonzero] - 0.09, angle = 90,
             code = 3, length = 0.03, col = red, lwd = 1.1)
    }
    points(second$mean, y - 0.09, pch = 22, col = red, bg = red, cex = 1.0)
  }
  if (!is.null(legend_labels)) legend("bottomright", legend_labels, col = c(teal, red),
    pch = c(21, 22), pt.bg = c(teal, red), bty = "n", cex = 0.72)
}

robust_metric <- function(metric_name) metric_rows(robust_summary, metric_name, "scenario_id", scenario_order)
robust_figure <- file.path(figure_dir, "fssgl_v2_robustness.pdf")
pdf(robust_figure, width = 8.2, height = 8.2, family = "Helvetica")
par(mfrow = c(2, 2), mar = c(4.8, 7.1, 2.4, 0.7), las = 1)
draw_robust(robust_metric("covariate_tpr"), robust_metric("group_tpr"), c(0, 1.05),
            "True-positive rate", "(a) Support recovery", c("Covariate", "Group"))
draw_robust(robust_metric("covariate_fpr"), robust_metric("covariate_fdr"), c(0, 0.5),
            "Rate", "(b) False selections", c("Covariate FPR", "Covariate FDR"))
draw_robust(robust_metric("coefficient_relative_error"), NULL, c(0.45, 1.0),
            "Coefficient relative error", "(c) Coefficient estimation")
draw_robust(robust_metric("test_coeff_rmse"), NULL, c(0.75, 1.6),
            "Test coefficient RMSE", "(d) Held-out prediction")
dev.off()

cat("Wrote MCSE summaries and three manuscript figures.\n")
