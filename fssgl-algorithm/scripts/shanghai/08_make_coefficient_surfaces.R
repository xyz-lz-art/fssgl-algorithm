# Reconstruct current coefficient surfaces from the prespecified detailed fit.

library(data.table)

source("R/fssgl/basis_design.R")
source("R/application/shanghai_workflow.R")
source("R/visualization/coefficient_surfaces.R")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
processed_dir <- file.path(root, "data/processed/shanghai_metroflow")
metadata_dir <- file.path(root, "data/metadata/shanghai_metroflow")
multi_dir <- file.path(processed_dir, "orthonormal/multi_targets")
figure_dir <- file.path(root, "results/figures/shanghai_metroflow/orthonormal")
table_dir <- file.path(root, "results/tables/shanghai_metroflow/orthonormal")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

target_name <- "Pudong International Airport"
stem <- clean_shanghai_name(target_name)
fit_file <- file.path(multi_dir, paste0("fssgl_", stem, ".rds"))
station_file <- file.path(multi_dir, paste0("selected_stations_", stem, ".csv"))
if (!file.exists(fit_file) || !file.exists(station_file)) {
  stop("Missing current detailed fit. Run scripts/shanghai/04_fit_targets.R first.")
}

fit_obj <- readRDS(fit_file)
target_obj <- prepare_shanghai_target(target_name, processed_dir, metadata_dir)
design_obj <- build_shanghai_basis_design(target_obj, kx = fit_obj$kx, ky = fit_obj$ky)
if (!identical(fit_obj$basis_geometry, "orthonormal_trapezoid")) {
  stop("Detailed fit does not use the current orthonormal basis geometry.")
}
design_obj$basis_x <- fit_obj$basis_x
design_obj$basis_y <- fit_obj$basis_y
selected <- fread(station_file)[order(-posterior_slab_prob)]
selected <- selected[seq_len(min(4L, nrow(selected)))]
selected <- merge(
  selected[, .(station, name, lines, posterior_slab_prob)],
  design_obj$station_groups,
  by = c("station", "name", "lines"),
  all.x = TRUE
)
if (nrow(selected) == 0L || anyNA(selected$block_start)) {
  stop("No reconstructable selected stations were found.")
}

surfaces <- lapply(seq_len(nrow(selected)), function(i) {
  block_beta <- fit_obj$fit$beta[selected$block_start[i]:selected$block_end[i]]
  beta_matrix <- matrix(block_beta, nrow = fit_obj$kx, ncol = fit_obj$ky)
  evaluate_coefficient_surface(beta_matrix, design_obj$basis_x, design_obj$basis_y)
})
shared_limit <- symmetric_surface_limit(surfaces)
out_file <- file.path(figure_dir, paste0("coefficient_surface_", stem, "_current.png"))
plot_estimated_surface_grid(
  surfaces = surfaces,
  labels = sprintf(
    "%s\n%s, slab score %.3f",
    selected$name,
    selected$lines,
    selected$posterior_slab_prob
  ),
  output_file = out_file,
  x = surface_axis(design_obj$basis_x),
  y = surface_axis(design_obj$basis_y),
  columns = 2L,
  note = paste(
    sprintf(
      "Shared symmetric coefficient scale: [%.3f, %.3f].",
      shared_limit[1L], shared_limit[2L]
    ),
    "Local peaks should be interpreted cautiously because the joint design is weakly identified."
  )
)

fwrite(
  selected[, .(station, name, lines, posterior_slab_prob, block_start, block_end)],
  file.path(table_dir, paste0("coefficient_surface_", stem, "_current_stations.csv"))
)
cat("Saved current coefficient surfaces:", out_file, "\n")
