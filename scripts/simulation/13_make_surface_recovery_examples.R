# Build a truth-estimate-difference figure from the designated strictly
# converged fit saved by 11_fair_fssgl_comparison.R.

source("R/fssgl/basis_design.R")
source("R/visualization/coefficient_surfaces.R")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_file <- file.path(
  root,
  "data/processed/simulation/v3_fair/fssgl_v3_surface_example_p20_rep1.rds"
)
if (!file.exists(input_file)) {
  stop(
    "Missing the designated fair-tuning fit. Run ",
    "scripts/simulation/11_fair_fssgl_comparison.R first."
  )
}
bundle <- readRDS(input_file)
dgp <- bundle$dgp
fit <- bundle$fit
if (!isTRUE(fit$fit$convergence$strict_converged)) {
  stop("The designated surface example did not satisfy strict convergence.")
}

block_norms <- vapply(fit$design$predictor_blocks, function(columns) {
  sqrt(sum(fit$fit$beta[columns]^2))
}, numeric(1L))
inactive <- setdiff(
  seq_len(dgp$dimensions$n_covariates),
  dgp$active_covariates
)
inactive_example <- inactive[which.max(block_norms[inactive])]
predictor_ids <- c(head(dgp$active_covariates, 2L), inactive_example)
predictor_labels <- c(
  paste0("Active predictor ", predictor_ids[1:2]),
  paste0("Inactive predictor ", inactive_example)
)

grid <- seq(0, 1, length.out = 61L)
basis_x <- make_bspline_basis(grid, df = dgp$dimensions$kx)
basis_y <- make_bspline_basis(grid, df = dgp$dimensions$ky)
output_file <- file.path(
  root,
  "results/figures/simulation/v3_fair/fssgl_v3_surface_recovery_p20_rep1.png"
)
plot_surface_recovery(
  beta_true = dgp$beta_true,
  beta_hat = fit$fit$beta,
  predictor_blocks = fit$design$predictor_blocks,
  basis_x = basis_x,
  basis_y = basis_y,
  predictor_ids = predictor_ids,
  predictor_labels = predictor_labels,
  output_file = output_file,
  note = paste(
    "Blue denotes negative and red positive values.",
    "Truth and estimate share a scale; differences use a separate scale."
  )
)
cat("Saved surface-recovery figure:", output_file, "\n")
