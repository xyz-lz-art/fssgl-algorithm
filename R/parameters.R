# Canonical parameters for the formal simulations and Shanghai extension.

# A change to the objective, update equations, thresholds, convergence rule, or
# tuning protocol requires a new method-release identifier.
FSSGL_METHOD_RELEASE <- "fssgl-orthonormal-2026.09"
FSSGL_WEIGHTED_MEMBERSHIP_VERSION <- "normalized-adaptive-block-orthonormal-v2.0"
FSSGL_ALGORITHM_VERSION <- "composite-nested-gem-v2.0"
FSSGL_RESULT_SCHEMA_VERSION <- "2.0"
FSSGL_BASIS_GEOMETRY_VERSION <- "orthonormal-trapezoid-v1.0"
FSSGL_TUNING_PROTOCOL_VERSION <- "common-curve-folds-one-se-v1.0"

fssgl_parameters <- function(max_outer_iter = 500L) {
  list(
    lambda_covariate_spike = 40,
    lambda_covariate_slab = 1,
    theta_covariate = 0.35,
    lambda_group_spike = 60,
    lambda_group_slab = 1,
    theta_group = 0.35,
    covariate_penalty_scale = 3,
    group_penalty_scale = 0.5,
    posterior_cutoff = 0.5,
    ridge_lambda = 1,
    score_mode = "dimension_correct",
    use_variance_update = TRUE,
    update_theta = TRUE,
    theta_prior_covariate = c(1, 1),
    theta_prior_group = c(1, 1),
    theta_damping = 0.5,
    theta_bounds = c(0.02, 0.98),
    sigma2_init = 1,
    sigma2_damping = 0.25,
    variance_prior_shape = 100,
    variance_prior_scale = 101,
    max_outer_iter = as.integer(max_outer_iter),
    inner_max_iter = 1L,
    inner_tol = 1e-6,
    tol = 1e-4,
    relaxed_tol = 0.05,
    objective_tail_tol = 0.12,
    lipschitz_multiplier = 1.05
  )
}

# Settings for common-fold comparisons under the current method release.
fssgl_submission_parameters <- function(p, n_train = 40L) {
  p <- as.integer(p)
  n_train <- as.integer(n_train)
  if (length(p) != 1L || !is.finite(p) || p < 1L) {
    stop("p must be one positive integer.")
  }
  if (length(n_train) != 1L || !is.finite(n_train) || n_train < 4L) {
    stop("n_train must be one integer of at least four.")
  }
  high_dimensional <- p >= n_train
  parameters <- fssgl_parameters(
    max_outer_iter = if (high_dimensional) 750L else 500L
  )
  parameters$inner_max_iter <- if (high_dimensional) 20L else 5L
  parameters$inner_tol <- if (high_dimensional) 1e-5 else 1e-6
  parameters
}

fssgl_penalty_multiplier_grid <- function(
  core_grid = 2^seq(-11, 1, by = 2),
  lower_expansions = 3L
) {
  core_grid <- sort(unique(as.numeric(core_grid)))
  lower_expansions <- as.integer(lower_expansions)
  if (!length(core_grid) || any(!is.finite(core_grid)) || any(core_grid <= 0)) {
    stop("core_grid must contain finite positive multipliers.")
  }
  if (length(lower_expansions) != 1L || lower_expansions < 0L) {
    stop("lower_expansions must be one nonnegative integer.")
  }
  lower <- if (lower_expansions == 0L) {
    numeric()
  } else {
    min(core_grid) / 2^seq_len(lower_expansions)
  }
  sort(unique(c(lower, core_grid)))
}

fssgl_baseline_parameters <- function() {
  list(
    cv_folds = 5L,
    ridge_lambda_grid = 10^seq(-4, 2, length.out = 13L),
    fpca_variance_threshold = 0.95,
    group_lasso_lambda_fractions = exp(seq(log(1), log(0.02), length.out = 15L)),
    group_lasso_max_iter = 2000L,
    group_lasso_tol = 1e-6,
    group_scad_nlambda = 50L,
    group_scad_lambda_min = 0.01,
    group_scad_gamma = 4,
    aenet_alpha_grid = c(0.3, 0.5, 0.7, 0.9),
    aenet_lambda_fractions = exp(seq(log(1), log(0.02), length.out = 15L)),
    aenet_adaptive_gamma = 1,
    aenet_initial_ridge = 1,
    aenet_max_iter = 2000L,
    aenet_tol = 1e-6,
    structured_alpha_grid = c(0.25, 0.5, 0.75),
    structured_lambda_fractions = exp(seq(log(1), log(0.02), length.out = 15L)),
    structured_max_iter = 2000L,
    structured_tol = 1e-6,
    kernel_gamma_multipliers = c(0.25, 0.5, 1, 2, 4),
    kernel_lambda_grid = 10^seq(-4, 1, length.out = 8L)
  )
}

fssgl_main_dgp_defaults <- function() {
  list(
    n_train = 40L,
    n_test = 24L,
    kx = 4L,
    ky = 4L,
    rho_common = 0.10,
    rho_group = 0.30,
    snr = 2,
    surface_smoothness = "moderate",
    signal_scale = 1
  )
}

shanghai_target_names <- function() {
  c(
    "People's Square",
    "Xujiahui",
    "Pudong International Airport",
    "East Nanjing Road",
    "Hongqiao Railway Station",
    "Shanghai Railway Station",
    "Dishui Lake",
    "Century Avenue",
    "Longyang Road",
    "Shanghai South Railway Station",
    "Lujiazui",
    "Hanzhong Road"
  )
}

shanghai_main_parameters <- function(max_iter = 1500L) {
  list(
    config_id = "orthonormal_cov1000_group1400_scales0002_00004_v1",
    algorithm_version = FSSGL_WEIGHTED_MEMBERSHIP_VERSION,
    lambda_station_spike = 1000,
    lambda_station_slab = 2.4,
    theta_station = 0.45,
    lambda_line_spike = 1400,
    lambda_line_slab = 2.4,
    theta_line = 0.45,
    step_multiplier = 90,
    station_penalty_scale = 0.002,
    line_penalty_scale = 0.0004,
    posterior_cutoff = 0.5,
    ridge_lambda = 10,
    max_iter = as.integer(max_iter),
    tol = 1e-3,
    relaxed_tol = 0.012,
    objective_tail_tol = 0.08
  )
}

new_experiment_manifest <- function(
  experiment_id,
  parameters,
  design = list(),
  algorithm_version = FSSGL_WEIGHTED_MEMBERSHIP_VERSION
) {
  list(
    experiment_id = experiment_id,
    method_release = FSSGL_METHOD_RELEASE,
    algorithm_version = algorithm_version,
    result_schema_version = FSSGL_RESULT_SCHEMA_VERSION,
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    r_version = R.version.string,
    parameters = parameters,
    design = design
  )
}
