# Run from the repository root: Rscript examples/quickstart.R
source("R/parameters.R")
source("R/fssgl/basis_design.R")
source("R/fssgl/penalty_weights.R")
source("R/fssgl/weighted_membership_solver.R")
source("R/fssgl/solver.R")
source("R/fssgl/simulation_design.R")

sample <- generate_fssgl_simulation_dgp(
  n_train = 24L,
  n_test = 12L,
  n_covariates = 6L,
  kx = 4L,
  ky = 4L,
  n_groups = 2L,
  n_active_groups = 1L,
  n_active_covariates_per_group = 2L,
  seed = 20260927L
)
fit <- fit_fssgl(
  x_coef = sample$x_train,
  y_coef = sample$y_train,
  structural_membership = sample$structural_membership,
  parameters = fssgl_parameters(),
  verbose = FALSE
)
summary <- evaluate_fssgl_simulation_fit(fit, sample)
print(summary[, .(
  selected_covariate_count,
  selected_group_count,
  covariate_tpr,
  covariate_fdr,
  coefficient_relative_error,
  test_coeff_rmse,
  strict_converged
)])
