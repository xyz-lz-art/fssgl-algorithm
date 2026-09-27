root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

source(file.path(root, "R/parameters.R"))
source(file.path(root, "R/fssgl/basis_design.R"))
source(file.path(root, "R/fssgl/penalty_weights.R"))
source(file.path(root, "R/fssgl/weighted_membership_solver.R"))
source(file.path(root, "R/fssgl/solver.R"))
source(file.path(root, "R/fssgl/tuning.R"))
source(file.path(root, "R/fssgl/simulation_design.R"))
source(file.path(root, "R/application/shanghai_workflow.R"))
source(file.path(root, "R/baselines/functional_methods.R"))
source(file.path(root, "R/visualization/coefficient_surfaces.R"))

expect_true <- function(value, label) {
  if (!isTRUE(value)) stop("FAILED: ", label, call. = FALSE)
}

expect_equal <- function(value, expected, label, tolerance = 1e-8) {
  ok <- isTRUE(all.equal(value, expected, tolerance = tolerance, check.attributes = FALSE))
  if (!ok) {
    stop(
      "FAILED: ", label, "\nExpected: ", paste(expected, collapse = ", "),
      "\nObserved: ", paste(value, collapse = ", "),
      call. = FALSE
    )
  }
}

expect_error <- function(expr, pattern, label) {
  message <- tryCatch({
    force(expr)
    NA_character_
  }, error = conditionMessage)
  if (is.na(message) || !grepl(pattern, message, fixed = TRUE)) {
    stop("FAILED: ", label, call. = FALSE)
  }
}

cat("Running fixed method-release tests...\n")
expect_equal(
  FSSGL_METHOD_RELEASE,
  "fssgl-orthonormal-2026.09",
  "current method release is fixed"
)
release_manifest <- new_experiment_manifest("release-test", list())
expect_equal(
  release_manifest$method_release,
  FSSGL_METHOD_RELEASE,
  "experiment manifests record the method release"
)

cat("Running basis/design tests...\n")
day_split <- validate_shanghai_day_split(20L, 1:12, 13:20)
expect_equal(day_split$train_days, 1:12, "Shanghai custom training days")
expect_equal(day_split$test_days, 13:20, "Shanghai custom test days")
expect_error(
  validate_shanghai_day_split(20L, 1:12, 12:20),
  "must not overlap",
  "Shanghai train/test overlap is rejected"
)
fixture <- list(
  arr_in = array(seq_len(3L * 10L * 8L), c(3L, 10L, 8L)),
  arr_out = array(seq_len(3L * 10L * 8L) / 2, c(3L, 10L, 8L)),
  station_metadata = data.table::data.table(
    station_pos = 1:3, station = 1:3, name = c("Target", "A", "B"),
    lines = c("L1", "L1", "L2"), n_lines = 1L
  ),
  date_index = data.table::data.table(date = 1:10),
  time_index = data.table::data.table(time = 1:8),
  membership = data.table::data.table(stationID = 1:3, line = c("L1", "L1", "L2"))
)
split_before <- prepare_shanghai_target("Target", "", "", 1:6, 7:10, fixture)
fixture$arr_in[, 7:10, ] <- fixture$arr_in[, 7:10, ] + 1e6
fixture$arr_out[, 7:10, ] <- fixture$arr_out[, 7:10, ] + 1e6
split_after <- prepare_shanghai_target("Target", "", "", 1:6, 7:10, fixture)
expect_equal(split_before$X_train, split_after$X_train, "future inflows do not change training data")
expect_equal(split_before$Y_train, split_after$Y_train, "future outflows do not change training data")
expect_equal(split_before$scaling, split_after$scaling, "scaling uses training days only")
expect_equal(
  apply(split_before$X_train, c(1L, 3L), mean),
  matrix(0, nrow = dim(split_before$X_train)[1L], ncol = dim(split_before$X_train)[3L]),
  "Shanghai predictor mean curves are removed pointwise",
  tolerance = 1e-12
)
expect_equal(
  colMeans(split_before$Y_train),
  numeric(dim(split_before$Y_train)[2L]),
  "Shanghai response mean curve is removed pointwise",
  tolerance = 1e-12
)
expect_equal(
  vapply(seq_len(dim(split_before$X_train)[1L]), function(j) {
    stats::sd(as.vector(split_before$X_train[j, , ]))
  }, numeric(1L)),
  rep(1, dim(split_before$X_train)[1L]),
  "Shanghai centered predictor curves use one residual scale per station"
)
expect_equal(
  stats::sd(as.vector(split_before$Y_train)),
  1,
  "Shanghai centered response curves use one residual scale"
)
expect_equal(
  split_before$scaling$preprocessing_version,
  SHANGHAI_PREPROCESSING_VERSION,
  "Shanghai preprocessing version is recorded"
)
basis_x <- make_bspline_basis(seq(0, 1, length.out = 151L), df = 3L)
basis_y <- make_bspline_basis(seq(0, 1, length.out = 161L), df = 5L)
expect_equal(
  attr(basis_x, "gram"),
  diag(3L),
  "predictor basis is orthonormal under quadrature",
  tolerance = 1e-10
)
expect_equal(
  attr(basis_y, "gram"),
  diag(5L),
  "response basis is orthonormal under quadrature",
  tolerance = 1e-10
)

geometry_matrix <- matrix(seq_len(15L) / 15, nrow = 3L, ncol = 5L)
geometry_surface <- evaluate_coefficient_surface(geometry_matrix, basis_x, basis_y)
expect_equal(
  surface_l2_norm(geometry_surface, basis_x, basis_y),
  sqrt(sum(geometry_matrix^2)),
  "coefficient Frobenius norm equals surface L2 norm",
  tolerance = 1e-9
)

dimension_state <- fssgl_v2_slab_state(
  norm_value = 0,
  dimension = 16L,
  lambda_spike = 40,
  lambda_slab = 1,
  theta = 0.35,
  score_mode = "dimension_correct"
)
expect_equal(dimension_state$spike_rate, 160, "dimension-corrected spike rate at q = 16")
expect_equal(dimension_state$slab_rate, 4, "dimension-corrected slab rate at q = 16")

structured_membership_fixture <- data.table::data.table(
  covariate_id = 1:4,
  group_id = c(1L, 1L, 2L, 2L),
  membership_weight = 1
)
structured_groups <- build_structured_baseline_groups(
  array(0, dim = c(4L, 8L, 2L)),
  structured_membership_fixture
)
expect_equal(
  structured_groups$parent[[1L]],
  1:4,
  "structured baseline parent group contains its two predictor blocks"
)
structured_path <- structured_group_lasso_path(
  x = matrix(seq_len(64L) / 64, nrow = 8L),
  y = matrix(seq_len(16L) / 16, nrow = 8L),
  predictor_groups = structured_groups$predictor,
  parent_groups = structured_groups$parent,
  alpha = 0.5,
  lambda_fractions = c(1, 0.2),
  max_iter = 200L,
  tol = 1e-5
)

expect_equal(length(structured_path$fits), 2L, "structured baseline returns its penalty path")

predictor_coef <- c(0.5, -1, 0.75)
predictor_curve <- drop(basis_x %*% predictor_coef)
integrated_effect <- colSums(
  geometry_surface * predictor_curve * attr(basis_x, "quadrature_weights")
)
expected_effect <- drop(basis_y %*% drop(crossprod(geometry_matrix, predictor_coef)))
expect_equal(
  integrated_effect,
  expected_effect,
  "continuous integral equals coefficient-space prediction for K_X != K_Y",
  tolerance = 1e-9
)

x_coef <- array(seq_len(2 * 3 * 2), dim = c(2, 3, 2))
design <- build_fof_design(x_coef, n_response_basis = 2)
expect_equal(dim(design$design), c(6, 8), "stacked design dimensions")
expect_equal(design$predictor_blocks[[1]], 1:4, "first block ordering")
expect_equal(design$predictor_blocks[[2]], 5:8, "second block ordering")

# Unequal marginal dimensions make transposition and vectorization errors visible.
x_nonsquare <- array(c(1, 2, 3, 4, 5, 6), dim = c(1, 2, 3))
a_nonsquare <- matrix(c(0.5, -1, 2, 1.5, 0, -0.5), nrow = 3, ncol = 2)
nonsquare_design <- build_fof_design(x_nonsquare, n_response_basis = 2)$design
expect_equal(
  as.numeric(nonsquare_design %*% as.vector(a_nonsquare)),
  as.vector(x_nonsquare[1, , ] %*% a_nonsquare),
  "coefficient-space design orientation for K_X != K_Y"
)

cat("Running membership validation tests...\n")
blocks <- build_functional_groups(design$predictor_blocks)
valid_membership <- data.table::data.table(
  covariate_id = c(1L, 2L),
  group_id = c("g1", "g2"),
  membership_weight = c(1, 1)
)
validate_fssgl_inputs(design$design, rep(0, 6), blocks, valid_membership, 4L)

bad_membership <- data.table::copy(valid_membership)
bad_membership[1, membership_weight := 0.5]
expect_error(
  validate_fssgl_inputs(design$design, rep(0, 6), blocks, bad_membership, 4L),
  "sum to one",
  "membership weights must sum to one"
)

duplicate_membership <- rbind(valid_membership, valid_membership[1])
expect_error(
  validate_fssgl_inputs(design$design, rep(0, 6), blocks, duplicate_membership, 4L),
  "must be unique",
  "duplicate membership pairs are rejected"
)

cat("Running distributed membership-weight tests...\n")
overlap_membership <- data.table::data.table(
  covariate_id = c(1L, 1L, 2L),
  group_id = c("g1", "g2", "g2"),
  membership_weight = c(0.5, 0.5, 1)
)
validate_fssgl_inputs(design$design, rep(0, 6), blocks, overlap_membership, 4L)
overlap_norms <- compute_group_norms(c(`1` = 3, `2` = 4), overlap_membership)
expect_equal(
  overlap_norms[group_id == "g1", group_norm],
  sqrt(0.5 * 3^2),
  "weighted norm for a shared membership"
)
expect_equal(
  overlap_norms[group_id == "g2", group_norm],
  sqrt(0.5 * 3^2 + 4^2),
  "weighted norm combines shared and exclusive members"
)
overlap_groups <- overlap_membership[, .(
  group_penalty_weight = sqrt(4 * sum(membership_weight))
), by = group_id]
overlap_add <- group_penalty_contribution(
  overlap_groups,
  overlap_membership,
  c(g1 = 2, g2 = 3)
)
expect_equal(
  overlap_add[covariate_id == 1L, group_penalty_add],
  0.5 * sqrt(2) * 2 + 0.5 * sqrt(6) * 3,
  "shared covariate receives distributed group penalties"
)
expect_equal(
  overlap_add[covariate_id == 2L, group_penalty_add],
  sqrt(6) * 3,
  "exclusive covariate receives its full group share"
)

cat("Running weighted-membership solver smoke test...\n")
weighted_fit <- fit_fssgl_solver(
  x = design$design,
  y = rep(0, nrow(design$design)),
  covariate_blocks = blocks,
  structural_membership = overlap_membership,
  block_size = 4L,
  max_iter = 2L,
  verbose = FALSE
)
expect_true(all(is.finite(weighted_fit$beta)), "weighted solver coefficients are finite")
expect_equal(
  sort(unique(weighted_fit$structural_membership$covariate_id)),
  1:2,
  "weighted solver retains both covariates"
)

cat("Preparing deterministic simulation fixture...\n")
dgp <- generate_fssgl_simulation_dgp(
  n_train = 24,
  n_test = 12,
  n_covariates = 6,
  kx = 4,
  ky = 4,
  n_groups = 2,
  n_active_groups = 1,
  n_active_covariates_per_group = 2,
  rho_common = 0.05,
  rho_group = 0.20,
  snr = 3,
  surface_smoothness = "smooth",
  seed = 20260710
)

if (!requireNamespace("grpreg", quietly = TRUE)) {
  stop("Required test dependency grpreg is unavailable; SCAD tests cannot be skipped.")
}
scad_folds <- make_baseline_cv_folds(dim(dgp$x_train)[2L], 3L, 42L)
stopifnot(identical(
  match_grpreg_fraction_path(
    list(lambda = c(10, 5, 2)), c(1, 0.6, 0.1)
  ),
  c(1L, 2L, 3L)
))
scad_fit <- fit_fpca_group_scad_baseline(
  x_coef = dgp$x_train,
  y_coef = dgp$y_train,
  folds = scad_folds,
  variance_threshold = 0.95,
  nlambda = 50L
)
stopifnot(
  identical(scad_fit$method_id, "fpca_group_scad"),
  length(scad_fit$beta) == length(dgp$beta_true),
  all(scad_fit$selected_covariates %in% seq_len(dgp$dimensions$n_covariates)),
  all(is.finite(unlist(scad_fit$timing))),
  all(dim(scad_fit$tuning$losses) == c(3L, 50L)),
  scad_fit$tuning$returned_nlambda >= 1L,
  scad_fit$tuning$returned_nlambda <= 50L
)
basis_scad_config <- fssgl_baseline_parameters()
basis_scad_config$group_scad_nlambda <- 50L
basis_scad_fit <- fit_functional_baseline(
  method_id = "basis_group_scad",
  x_coef = dgp$x_train,
  y_coef = dgp$y_train,
  folds = scad_folds,
  config = basis_scad_config
)
basis_scad_prediction <- predict_functional_baseline(
  basis_scad_fit, dgp$x_test
)
stopifnot(
  identical(basis_scad_fit$method_id, "basis_group_scad"),
  identical(basis_scad_fit$basis_geometry, "orthonormal_bspline_coefficients"),
  length(basis_scad_fit$beta) == length(dgp$beta_true),
  all(basis_scad_fit$selected_covariates %in% seq_len(dgp$dimensions$n_covariates)),
  all(is.finite(unlist(basis_scad_fit$timing))),
  all(dim(basis_scad_fit$tuning$losses) == c(3L, 50L)),
  basis_scad_fit$tuning$returned_nlambda >= 1L,
  basis_scad_fit$tuning$returned_nlambda <= 50L,
  all(dim(basis_scad_prediction) == dim(dgp$y_test)),
  all(is.finite(basis_scad_prediction))
)

cat("Running formal slab-score and nested-prox tests...\n")
test_norm <- 0.2
test_dimension <- 4
test_theta <- 0.35
test_spike <- 10
test_slab <- 1
test_state <- fssgl_v2_slab_state(
  test_norm,
  test_dimension,
  test_spike,
  test_slab,
  test_theta,
  score_mode = "dimension_correct"
)
test_spike_rate <- sqrt(test_dimension) * test_spike
test_slab_rate <- sqrt(test_dimension) * test_slab
expected_log_odds <- log(test_theta / (1 - test_theta)) +
  test_dimension * log(test_slab_rate / test_spike_rate) +
  (test_spike_rate - test_slab_rate) * test_norm
expect_equal(
  qlogis(test_state$slab_probability),
  expected_log_odds,
  "dimension-correct slab log-odds"
)
transition_norm <- fssgl_v2_responsibility_transition(
  test_dimension,
  test_spike,
  test_slab,
  test_theta,
  score_mode = "dimension_correct"
)
transition_state <- fssgl_v2_slab_state(
  transition_norm,
  test_dimension,
  test_spike,
  test_slab,
  test_theta,
  score_mode = "dimension_correct"
)
expect_equal(
  transition_state$slab_probability,
  0.5,
  "responsibility transition norm yields equal spike and slab weights"
)

prox_blocks <- data.table::data.table(
  covariate_id = 1:2,
  block_start = c(1L, 3L),
  block_end = c(2L, 4L)
)
prox_membership <- data.table::data.table(
  covariate_id = 1:2,
  group_id = "g1",
  membership_weight = 1
)
prox_groups <- build_v2_group_structure(prox_blocks, prox_membership, block_size = 2L)
prox_value <- c(3, 4, 0, 5)
prox_observed <- prox_nested_fssgl(
  value = prox_value,
  step_size = 1,
  sample_size = 1,
  variance_scale = 1,
  covariate_penalty = c(1, 1),
  group_penalty = 2,
  covariate_blocks = prox_blocks,
  group_structure = prox_groups
)
after_block <- c(2.4, 3.2, 0, 4)
prox_expected <- (1 - 2 / sqrt(sum(after_block^2))) * after_block
expect_equal(prox_observed, prox_expected, "nested block-then-group proximal map")
expect_equal(
  update_fssgl_variance(c(1, -1), sample_size = 2, prior_shape = 0, prior_scale = 0),
  0.5,
  "Jeffreys residual-variance update"
)
expect_equal(
  update_fssgl_variance(c(1, -1), sample_size = 2, prior_shape = 2, prior_scale = 3),
  1,
  "proper inverse-gamma residual-variance update"
)

cat("Running formal FSSGL smoke test...\n")
v2_smoke_arguments <- list(
  x_coef = dgp$x_train,
  y_coef = dgp$y_train,
  structural_membership = dgp$structural_membership,
  lambda_covariate_spike = 80,
  lambda_covariate_slab = 1,
  theta_covariate = 0.30,
  lambda_group_spike = 80,
  lambda_group_slab = 1,
  theta_group = 0.30,
  covariate_penalty_scale = 0.02,
  group_penalty_scale = 0.01,
  score_mode = "dimension_correct",
  use_variance_update = TRUE,
  update_theta = TRUE,
  max_outer_iter = 5L,
  inner_max_iter = 3L,
  beta_init = ridge_dual_fit(
    build_fof_design(dgp$x_train, n_response_basis = 4)$design,
    as.vector(dgp$y_train),
    lambda = 1
  ),
  verbose = FALSE
)
v2_fit <- do.call(fit_fssgl_coefficients, v2_smoke_arguments)
v2_traced_fit <- do.call(
  fit_fssgl_coefficients,
  c(v2_smoke_arguments, list(trace_responsibilities = TRUE))
)
expect_true(all(is.finite(v2_fit$fit$beta)), "FSSGL coefficients are finite")
expect_equal(
  v2_traced_fit$fit$beta,
  v2_fit$fit$beta,
  "responsibility tracing does not change coefficient updates",
  tolerance = 0
)
expect_equal(
  v2_traced_fit$fit$history,
  v2_fit$fit$history,
  "responsibility tracing does not change iteration history",
  tolerance = 0
)
expected_trace_rows <- nrow(v2_traced_fit$fit$history) *
  (dgp$dimensions$n_covariates + dgp$dimensions$n_groups)
expect_equal(
  nrow(v2_traced_fit$fit$responsibility_trace),
  expected_trace_rows,
  "responsibility trace contains every covariate and group at every iteration"
)
expect_true(
  all(is.finite(v2_traced_fit$fit$responsibility_trace$transition_norm)) &&
    all(v2_traced_fit$fit$responsibility_trace$transition_norm >= 0),
  "responsibility transition norms are finite and nonnegative"
)
expect_true(
  v2_fit$fit$hyperparameters$theta_covariate_final > 0 &&
    v2_fit$fit$hyperparameters$theta_covariate_final < 1,
  "FSSGL theta update remains in the unit interval"
)
expect_true(
  v2_fit$fit$hyperparameters$sigma2_final > 0,
  "FSSGL variance update remains positive"
)

selected_parameters <- fssgl_parameters(max_outer_iter = 5L)
selected_fit <- fit_fssgl(
  x_coef = dgp$x_train,
  y_coef = dgp$y_train,
  structural_membership = dgp$structural_membership,
  parameters = selected_parameters,
  verbose = FALSE
)
expect_equal(
  selected_fit$fit$hyperparameters$variance_prior_shape,
  100,
  "formal FSSGL uses the frozen variance prior"
)
expect_equal(
  selected_fit$fit$hyperparameters$inner_max_iter,
  1L,
  "formal FSSGL uses one generalized-EM M-step"
)
expect_true(
  selected_fit$fit$hyperparameters$update_theta,
  "formal FSSGL updates theta"
)

cat("Running common-fold tuning and convergence-extension tests...\n")
shared_folds <- make_curve_cv_folds(n = dgp$dimensions$n_train, n_folds = 3L, seed = 42L)
expect_equal(
  shared_folds,
  make_baseline_cv_folds(dgp$dimensions$n_train, 3L, 42L),
  "FSSGL and baselines use identical curve folds"
)
expect_true(
  min(fssgl_penalty_multiplier_grid(core_grid = c(1, 2), lower_expansions = 2L)) == 0.25,
  "penalty grid expands below its initial lower boundary"
)
selection_fixture <- data.table::data.table(
  multiplier = c(0.5, 1, 2),
  mean_validation_rmse = c(1.00, 1.01, 1.04),
  se_validation_rmse = c(0.05, 0.05, 0.05),
  strict_convergence_rate = 1,
  all_folds_complete = TRUE
)
expect_equal(
  select_fssgl_cv_candidate(selection_fixture)$multiplier,
  2,
  "one-standard-error rule breaks ties toward the sparser larger penalty"
)
plateau_fixture <- data.table::data.table(
  multiplier = c(0.5, 2, 8),
  mean_validation_rmse = c(1.1, 1.2, 1.2),
  mean_selected_covariates = c(2, 0, 0),
  strict_convergence_rate = 1,
  all_folds_complete = TRUE
)
expect_true(
  is_stable_upper_null_plateau(plateau_fixture),
  "identical upper all-zero fits establish a resolved penalty plateau"
)
plateau_fixture$mean_selected_covariates[2L] <- 1
expect_true(
  !is_stable_upper_null_plateau(plateau_fixture),
  "a changing upper support remains an unresolved grid boundary"
)

extension_parameters <- fssgl_parameters(max_outer_iter = 2L)
extension_parameters$tol <- 1e-12
extended_fit <- fit_fssgl_until_converged(
  x_coef = dgp$x_train,
  y_coef = dgp$y_train,
  structural_membership = dgp$structural_membership,
  parameters = extension_parameters,
  extension_factors = c(1L, 2L),
  verbose = FALSE
)
expect_equal(
  extended_fit$fit$convergence$extension_stages,
  2L,
  "nonconverged fit receives a warm continuation stage"
)
expect_equal(
  extended_fit$fit$convergence$total_outer_iter,
  nrow(extended_fit$fit$history),
  "continuation history records every outer iteration"
)

tuning_parameters <- fssgl_parameters(max_outer_iter = 5L)
tuning_parameters$tol <- 10
tuning_folds <- make_curve_cv_folds(
  n = dgp$dimensions$n_train,
  n_folds = 2L,
  seed = 43L
)
tuning_fit <- tune_fssgl_cv(
  x_coef = dgp$x_train,
  y_coef = dgp$y_train,
  structural_membership = dgp$structural_membership,
  folds = tuning_folds,
  multiplier_grid = c(0.5, 1),
  parameters = tuning_parameters,
  selection_rule = "one_se_sparsest",
  extension_factors = 1L,
  require_interior_minimum = FALSE
)
expect_true(
  tuning_fit$selected_multiplier %in% c(0.5, 1) &&
    all(tuning_fit$summary$strict_convergence_rate == 1),
  "curve-level FSSGL tuning returns only strictly converged candidates"
)
selected_path_fit <- fit_fssgl_selected_path(
  x_coef = dgp$x_train,
  y_coef = dgp$y_train,
  structural_membership = dgp$structural_membership,
  multiplier_grid = c(0.5, 1),
  selected_multiplier = 0.5,
  parameters = tuning_parameters,
  extension_factors = 1L
)
expect_equal(
  selected_path_fit$path$multiplier,
  c(1, 0.5),
  "full-training refit follows the same strong-to-weak path as CV"
)
expect_true(
  all(is.finite(selected_path_fit$fit$fit$beta)),
  "selected-path refit returns finite coefficients"
)

cat("Running coefficient-surface visualization tests...\n")
surface_grid <- seq(0, 1, length.out = 31L)
surface_basis_x <- make_bspline_basis(surface_grid, df = dgp$dimensions$kx)
surface_basis_y <- make_bspline_basis(surface_grid, df = dgp$dimensions$ky)
surface_blocks <- build_fof_design(
  dgp$x_train,
  n_response_basis = dgp$dimensions$ky
)$predictor_blocks
surface_ids <- c(dgp$active_covariates[1L], setdiff(
  seq_len(dgp$dimensions$n_covariates), dgp$active_covariates
)[1L])
surface_values <- reconstruct_coefficient_surfaces(
  dgp$beta_true,
  surface_blocks,
  surface_basis_x,
  surface_basis_y,
  surface_ids
)
expect_true(
  max(abs(surface_values[[2L]])) < 1e-12,
  "inactive predictor reconstructs to an exactly zero surface"
)
surface_test_file <- tempfile(fileext = ".png")
surface_plot <- plot_surface_recovery(
  beta_true = dgp$beta_true,
  beta_hat = 0.8 * dgp$beta_true,
  predictor_blocks = surface_blocks,
  basis_x = surface_basis_x,
  basis_y = surface_basis_y,
  predictor_ids = surface_ids,
  predictor_labels = c("Active", "Inactive"),
  output_file = surface_test_file,
  note = "Automated rendering test",
  width = 900,
  row_height = 280,
  resolution = 120
)
expect_true(
  file.exists(surface_test_file) && file.info(surface_test_file)$size > 1000,
  "surface-recovery plot is rendered to a nonempty image"
)
expect_true(
  all(surface_plot$effect_limit == -rev(surface_plot$effect_limit)) &&
    all(surface_plot$difference_limit == -rev(surface_plot$difference_limit)),
  "truth, estimate, and difference use symmetric color limits"
)
unlink(surface_test_file)

cat("All tests passed.\n")
