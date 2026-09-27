# Generic simulation data generators for FSSGL experiments.
# Public simulation objects use covariate/group terminology. The only
# conversion to solver-internal column names happens immediately before fitting.

make_simulation_group_membership <- function(
  n_covariates,
  n_groups
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is required.")
  }
  if (n_covariates < 1 || n_groups < 1) {
    stop("n_covariates and n_groups must be positive.")
  }
  if (n_covariates %% n_groups != 0) {
    stop("The current simulation design requires equal-size non-overlapping groups.")
  }

  group_index <- ((seq_len(n_covariates) - 1L) %% n_groups) + 1L
  data.table::data.table(
    covariate_id = seq_len(n_covariates),
    group_id = paste0("group_", group_index),
    membership_weight = 1
  )
}

make_simulation_surface <- function(
  kx,
  ky,
  amplitude = 1,
  smoothness = c("smooth", "moderate", "rough"),
  phase = 0
) {
  smoothness <- match.arg(smoothness)
  sx <- seq(0, 1, length.out = kx)
  sy <- seq(0, 1, length.out = ky)

  base <- outer(
    sx,
    sy,
    function(x, y) {
      sin(pi * (x + phase)) * cos(pi * y) +
        0.5 * cos(2 * pi * x) * sin(pi * (y + phase))
    }
  )

  if (smoothness %in% c("moderate", "rough")) {
    base <- base + 0.35 * outer(
      sx,
      sy,
      function(x, y) sin(3 * pi * x + phase) * sin(2 * pi * y)
    )
  }
  if (smoothness == "rough") {
    base <- base + 0.25 * outer(
      sx,
      sy,
      function(x, y) cos(5 * pi * x) * cos(4 * pi * y + phase)
    )
  }

  base_norm <- sqrt(sum(base^2))
  if (!is.finite(base_norm) || base_norm == 0) {
    return(matrix(0, nrow = kx, ncol = ky))
  }
  amplitude * base / base_norm
}

generate_fssgl_simulation_dgp <- function(
  n_train = 60,
  n_test = 40,
  n_covariates = 20,
  kx = 5,
  ky = 5,
  n_groups = 4,
  n_active_groups = 2,
  n_active_covariates_per_group = 2,
  rho_common = 0.10,
  rho_group = 0.30,
  snr = 2,
  surface_smoothness = c("smooth", "moderate", "rough"),
  signal_scale = 1,
  seed = NULL
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is required.")
  }
  surface_smoothness <- match.arg(surface_smoothness)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (n_active_groups > n_groups) {
    stop("n_active_groups cannot exceed n_groups.")
  }
  if (rho_common < 0 || rho_group < 0 || rho_common + rho_group >= 1) {
    stop("rho_common and rho_group must be nonnegative and sum to less than 1.")
  }
  if (snr <= 0) {
    stop("snr must be positive.")
  }

  structural_membership <- make_simulation_group_membership(
    n_covariates = n_covariates,
    n_groups = n_groups
  )
  group_lookup <- structural_membership[order(covariate_id)]

  active_group_ids <- paste0("group_", seq_len(n_active_groups))
  active_covariates <- unlist(lapply(active_group_ids, function(gid) {
    ids <- group_lookup[group_id == gid, covariate_id]
    head(ids, n_active_covariates_per_group)
  }), use.names = FALSE)
  active_covariates <- sort(unique(active_covariates))
  active_groups <- sort(unique(structural_membership[covariate_id %in% active_covariates, group_id]))

  make_x <- function(n_sample) {
    common_factor <- matrix(rnorm(n_sample * kx), nrow = n_sample, ncol = kx)
    group_factor <- array(rnorm(n_groups * n_sample * kx), dim = c(n_groups, n_sample, kx))
    x <- array(0, dim = c(n_covariates, n_sample, kx))

    for (j in seq_len(n_covariates)) {
      g <- as.integer(sub("group_", "", group_lookup$group_id[j]))
      idiosyncratic <- matrix(rnorm(n_sample * kx), nrow = n_sample, ncol = kx)
      x[j, , ] <-
        sqrt(rho_common) * common_factor +
        sqrt(rho_group) * group_factor[g, , ] +
        sqrt(1 - rho_common - rho_group) * idiosyncratic
    }
    x
  }

  x_train <- make_x(n_train)
  x_test <- make_x(n_test)
  design_train <- build_fof_design(x_train, n_response_basis = ky)
  design_test <- build_fof_design(x_test, n_response_basis = ky)

  beta_true <- numeric(n_covariates * kx * ky)
  surfaces <- vector("list", n_covariates)
  for (j in seq_len(n_covariates)) {
    surface <- matrix(0, nrow = kx, ncol = ky)
    if (j %in% active_covariates) {
      surface <- make_simulation_surface(
        kx = kx,
        ky = ky,
        amplitude = signal_scale,
        smoothness = surface_smoothness,
        phase = j / max(n_covariates, 1)
      )
    }
    surfaces[[j]] <- surface
    beta_true[design_train$predictor_blocks[[j]]] <- as.vector(surface)
  }

  signal_train <- drop(design_train$design %*% beta_true)
  signal_sd <- stats::sd(signal_train)
  noise_sd <- if (is.finite(signal_sd) && signal_sd > 0) signal_sd / sqrt(snr) else 1
  y_train_vec <- signal_train + rnorm(length(signal_train), sd = noise_sd)
  y_test_vec <- drop(design_test$design %*% beta_true) + rnorm(n_test * ky, sd = noise_sd)

  truth_covariates <- data.table::data.table(
    covariate_id = seq_len(n_covariates),
    active = seq_len(n_covariates) %in% active_covariates
  )
  truth_groups <- data.table::data.table(
    group_id = paste0("group_", seq_len(n_groups)),
    active = paste0("group_", seq_len(n_groups)) %in% active_groups
  )

  list(
    x_train = x_train,
    y_train = matrix(y_train_vec, nrow = n_train, ncol = ky),
    x_test = x_test,
    y_test = matrix(y_test_vec, nrow = n_test, ncol = ky),
    beta_true = beta_true,
    surfaces = surfaces,
    structural_membership = structural_membership,
    truth_covariates = truth_covariates,
    truth_groups = truth_groups,
    active_covariates = active_covariates,
    active_groups = active_groups,
    noise_sd = noise_sd,
    dimensions = list(
      n_train = n_train,
      n_test = n_test,
      n_covariates = n_covariates,
      kx = kx,
      ky = ky,
      n_groups = n_groups
    ),
    settings = list(
      basis_geometry = "orthonormal_coefficient_coordinates",
      rho_common = rho_common,
      rho_group = rho_group,
      snr = snr,
      surface_smoothness = surface_smoothness,
      signal_scale = signal_scale
    )
  )
}

evaluate_fssgl_simulation_fit <- function(fit_obj, dgp, posterior_cutoff = 0.5, runtime_sec = NA_real_) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is required.")
  }
  fit <- fit_obj$fit
  covariate_posterior <- data.table::copy(fit$covariate_posterior)

  selected_covariates <- covariate_posterior[
    posterior_slab_prob >= posterior_cutoff,
    covariate_id
  ]
  selected_groups <- unique(dgp$structural_membership[
    covariate_id %in% selected_covariates,
    group_id
  ])

  covariate_eval <- merge(
    dgp$truth_covariates,
    data.table::data.table(
      covariate_id = dgp$truth_covariates$covariate_id,
      selected = dgp$truth_covariates$covariate_id %in% selected_covariates
    ),
    by = "covariate_id"
  )
  group_eval <- merge(
    dgp$truth_groups,
    data.table::data.table(
      group_id = dgp$truth_groups$group_id,
      selected = dgp$truth_groups$group_id %in% selected_groups
    ),
    by = "group_id"
  )

  rate <- function(active, selected, positive) {
    denom <- sum(active == positive)
    if (denom == 0) return(NA_real_)
    sum(selected & active == positive) / denom
  }
  fdr <- function(active, selected) {
    denom <- sum(selected)
    if (denom == 0) return(0)
    sum(selected & !active) / denom
  }

  design_test <- build_fof_design(dgp$x_test, n_response_basis = dgp$dimensions$ky)
  pred_test <- matrix(
    drop(design_test$design %*% fit$beta),
    nrow = dgp$dimensions$n_test,
    ncol = dgp$dimensions$ky
  )

  block_size <- dgp$dimensions$kx * dgp$dimensions$ky
  active_cols <- unlist(lapply(dgp$active_covariates, function(j) {
    start <- (j - 1L) * block_size + 1L
    start:(start + block_size - 1L)
  }), use.names = FALSE)
  inactive_covariates <- setdiff(seq_len(dgp$dimensions$n_covariates), dgp$active_covariates)
  inactive_norm <- if (length(inactive_covariates) == 0) {
    NA_real_
  } else {
    mean(vapply(inactive_covariates, function(j) {
      start <- (j - 1L) * block_size + 1L
      sqrt(sum(fit$beta[start:(start + block_size - 1L)]^2))
    }, numeric(1)))
  }

  data.table::data.table(
    covariate_tpr = rate(covariate_eval$active, covariate_eval$selected, TRUE),
    covariate_fpr = rate(covariate_eval$active, covariate_eval$selected, FALSE),
    covariate_fdr = fdr(covariate_eval$active, covariate_eval$selected),
    selected_covariate_count = length(selected_covariates),
    group_tpr = rate(group_eval$active, group_eval$selected, TRUE),
    group_fpr = rate(group_eval$active, group_eval$selected, FALSE),
    group_fdr = fdr(group_eval$active, group_eval$selected),
    selected_group_count = length(selected_groups),
    coefficient_relative_error = sqrt(sum((fit$beta - dgp$beta_true)^2)) /
      (sqrt(sum(dgp$beta_true^2)) + 1e-8),
    active_coefficient_relative_error = sqrt(sum((fit$beta[active_cols] - dgp$beta_true[active_cols])^2)) /
      (sqrt(sum(dgp$beta_true[active_cols]^2)) + 1e-8),
    inactive_coefficient_norm = inactive_norm,
    test_coeff_rmse = sqrt(mean((pred_test - dgp$y_test)^2)),
    runtime_sec = runtime_sec,
    final_iter = fit$convergence$final_iter,
    strict_converged = fit$convergence$strict_converged,
    relaxed_converged = fit$convergence$relaxed_converged,
    final_beta_change = fit$convergence$final_beta_change,
    objective_tail_rel_change = fit$convergence$objective_tail_rel_change,
    final_residual_variance = fit$convergence$final_residual_variance
  )
}
