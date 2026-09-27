# Composite FSSGL solver for the formal non-overlapping estimator.
# The Shanghai weighted-membership application uses a separate solver in
# weighted_membership_solver.R.

fssgl_v2_slab_state <- function(
  norm_value,
  dimension,
  lambda_spike,
  lambda_slab,
  theta,
  score_mode = c("dimension_correct", "normalized_scalar")
) {
  score_mode <- match.arg(score_mode)
  if (lambda_spike <= lambda_slab) {
    stop("lambda_spike must be larger than lambda_slab.")
  }
  if (any(!is.finite(norm_value)) || any(norm_value < 0)) {
    stop("norm_value must be finite and nonnegative.")
  }
  if (length(dimension) == 1L) {
    dimension <- rep(dimension, length(norm_value))
  }
  if (length(dimension) != length(norm_value) || any(dimension < 1)) {
    stop("dimension must be positive and conformable with norm_value.")
  }
  theta <- min(max(theta, 1e-8), 1 - 1e-8)

  if (score_mode == "normalized_scalar") {
    normalized_norm <- norm_value / sqrt(dimension)
    slab_probability <- ssgl_slab_probability(
      normalized_norm,
      1,
      lambda_spike,
      lambda_slab,
      theta
    )
    effective_rate <- sqrt(dimension) * ssgl_effective_lambda(
      slab_probability,
      lambda_spike,
      lambda_slab
    )
    return(list(
      slab_probability = slab_probability,
      effective_rate = effective_rate,
      spike_rate = sqrt(dimension) * lambda_spike,
      slab_rate = sqrt(dimension) * lambda_slab
    ))
  }

  spike_rate <- sqrt(dimension) * lambda_spike
  slab_rate <- sqrt(dimension) * lambda_slab
  log_slab <- log(theta) + dimension * log(slab_rate) - slab_rate * norm_value
  log_spike <- log1p(-theta) + dimension * log(spike_rate) - spike_rate * norm_value
  slab_probability <- stable_logistic_from_logs(log_slab, log_spike)

  list(
    slab_probability = slab_probability,
    effective_rate = slab_probability * slab_rate + (1 - slab_probability) * spike_rate,
    spike_rate = spike_rate,
    slab_rate = slab_rate
  )
}

# Norm at which the spike and slab responsibilities are equal. This is a
# diagnostic quantity only; it is not used by the coefficient updates.
fssgl_v2_responsibility_transition <- function(
  dimension,
  lambda_spike,
  lambda_slab,
  theta,
  score_mode = c("dimension_correct", "normalized_scalar")
) {
  score_mode <- match.arg(score_mode)
  dimension <- as.numeric(dimension)
  theta <- pmin(pmax(as.numeric(theta), 1e-8), 1 - 1e-8)
  if (length(theta) == 1L) theta <- rep(theta, length(dimension))
  if (length(theta) != length(dimension) || any(dimension < 1)) {
    stop("dimension and theta must be positive and conformable.")
  }
  if (lambda_spike <= lambda_slab) {
    stop("lambda_spike must be larger than lambda_slab.")
  }
  if (score_mode == "normalized_scalar") {
    crossing <- sqrt(dimension) * (
      log((1 - theta) / theta) + log(lambda_spike / lambda_slab)
    ) / (lambda_spike - lambda_slab)
  } else {
    spike_rate <- sqrt(dimension) * lambda_spike
    slab_rate <- sqrt(dimension) * lambda_slab
    crossing <- (
      log((1 - theta) / theta) + dimension * log(spike_rate / slab_rate)
    ) / (spike_rate - slab_rate)
  }
  pmax(crossing, 0)
}

fssgl_v2_soft_threshold <- function(value, threshold) {
  value_norm <- sqrt(sum(value^2))
  if (!is.finite(value_norm) || value_norm <= threshold) {
    return(numeric(length(value)))
  }
  (1 - threshold / value_norm) * value
}

update_fssgl_variance <- function(
  residual,
  sample_size = length(residual),
  prior_shape = 0,
  prior_scale = 0
) {
  max(
    (sum(residual^2) + 2 * prior_scale) /
      (sample_size + 2 * prior_shape + 2),
    1e-8
  )
}

validate_nonoverlapping_partition <- function(covariate_blocks, structural_membership) {
  membership_count <- structural_membership[, .N, by = covariate_id]
  if (any(membership_count$N != 1L)) {
    stop("The formal solver requires one structural group per covariate.")
  }
  if (any(abs(structural_membership$membership_weight - 1) > 1e-8)) {
    stop("Non-overlapping memberships must have weight one.")
  }
  if (!setequal(covariate_blocks$covariate_id, structural_membership$covariate_id)) {
    stop("The structural partition must cover every covariate exactly once.")
  }
  invisible(TRUE)
}

build_v2_group_structure <- function(covariate_blocks, structural_membership, block_size) {
  group_structure <- merge(
    structural_membership[, .(covariate_id, group_id)],
    covariate_blocks[, .(covariate_id, block_start, block_end)],
    by = "covariate_id",
    all.x = TRUE
  )
  group_structure[
    ,
    .(
      member_count = .N,
      group_dimension = block_size * .N,
      group_columns = list(unlist(Map(seq.int, block_start, block_end), use.names = FALSE))
    ),
    by = group_id
  ][order(group_id)]
}

compute_v2_group_norms <- function(beta, group_structure) {
  vapply(
    group_structure$group_columns,
    function(columns) sqrt(sum(beta[columns]^2)),
    numeric(1)
  )
}

prox_nested_fssgl <- function(
  value,
  step_size,
  sample_size,
  variance_scale,
  covariate_penalty,
  group_penalty,
  covariate_blocks,
  group_structure
) {
  out <- value
  for (j in seq_len(nrow(covariate_blocks))) {
    columns <- covariate_blocks$block_start[j]:covariate_blocks$block_end[j]
    threshold <- step_size * variance_scale * covariate_penalty[j] / sample_size
    out[columns] <- fssgl_v2_soft_threshold(out[columns], threshold)
  }
  for (g in seq_len(nrow(group_structure))) {
    columns <- group_structure$group_columns[[g]]
    threshold <- step_size * variance_scale * group_penalty[g] / sample_size
    out[columns] <- fssgl_v2_soft_threshold(out[columns], threshold)
  }
  out
}

weighted_nested_objective <- function(
  x,
  y,
  beta,
  variance_scale,
  covariate_penalty,
  group_penalty,
  covariate_blocks,
  group_structure
) {
  sample_size <- nrow(x)
  residual <- y - drop(x %*% beta)
  covariate_norms <- compute_covariate_beta_norm(beta, covariate_blocks)
  group_norms <- compute_v2_group_norms(beta, group_structure)
  sum(residual^2) / (2 * sample_size) +
    variance_scale * (
      sum(covariate_penalty * covariate_norms) +
        sum(group_penalty * group_norms)
    ) / sample_size
}

solve_weighted_nested_mstep <- function(
  x,
  y,
  beta_init,
  variance_scale,
  covariate_penalty,
  group_penalty,
  covariate_blocks,
  group_structure,
  lipschitz_init,
  max_iter = 1L,
  tol = 1e-6,
  backtracking_factor = 2,
  max_backtracking = 20L
) {
  sample_size <- nrow(x)
  beta <- beta_init
  lipschitz <- lipschitz_init
  objective <- weighted_nested_objective(
    x, y, beta, variance_scale, covariate_penalty, group_penalty,
    covariate_blocks, group_structure
  )
  final_change <- Inf
  total_backtracking <- 0L

  for (inner_iter in seq_len(max_iter)) {
    beta_old <- beta
    residual <- drop(x %*% beta - y)
    gradient <- drop(crossprod(x, residual)) / sample_size
    accepted <- FALSE

    for (bt in 0:max_backtracking) {
      step_size <- 1 / lipschitz
      proposal <- prox_nested_fssgl(
        beta - step_size * gradient,
        step_size,
        sample_size,
        variance_scale,
        covariate_penalty,
        group_penalty,
        covariate_blocks,
        group_structure
      )
      proposal_objective <- weighted_nested_objective(
        x, y, proposal, variance_scale, covariate_penalty, group_penalty,
        covariate_blocks, group_structure
      )
      if (proposal_objective <= objective + 1e-12) {
        beta <- proposal
        objective <- proposal_objective
        total_backtracking <- total_backtracking + bt
        accepted <- TRUE
        break
      }
      lipschitz <- lipschitz * backtracking_factor
    }
    if (!accepted) {
      stop("Nested M-step backtracking failed to find a descent update.")
    }

    final_change <- sqrt(sum((beta - beta_old)^2)) /
      (sqrt(sum(beta_old^2)) + 1e-8)
    if (final_change < tol) {
      break
    }
  }

  list(
    beta = beta,
    objective = objective,
    final_change = final_change,
    iterations = inner_iter,
    lipschitz = lipschitz,
    backtracking_steps = total_backtracking
  )
}

fit_fssgl_core <- function(
  x,
  y,
  covariate_blocks,
  structural_membership,
  block_size,
  lambda_covariate_spike = 160,
  lambda_covariate_slab = 1,
  theta_covariate = 0.35,
  lambda_group_spike = 160,
  lambda_group_slab = 1,
  theta_group = 0.35,
  covariate_penalty_scale = 3,
  group_penalty_scale = 0.5,
  posterior_cutoff = 0.5,
  score_mode = c("dimension_correct", "normalized_scalar"),
  use_variance_update = FALSE,
  update_theta = FALSE,
  theta_prior_covariate = c(1, 1),
  theta_prior_group = c(1, 1),
  theta_damping = 0.5,
  theta_bounds = c(0.02, 0.98),
  sigma2_damping = 0.5,
  variance_prior_shape = 0,
  variance_prior_scale = 0,
  max_outer_iter = 500L,
  inner_max_iter = 1L,
  inner_tol = 1e-6,
  tol = 1e-4,
  relaxed_tol = 0.05,
  objective_tail_tol = 0.12,
  beta_init = NULL,
  sigma2_init = NULL,
  lipschitz_multiplier = 1.05,
  trace_responsibilities = FALSE,
  verbose = FALSE
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is required.")
  }
  score_mode <- match.arg(score_mode)
  covariate_blocks <- data.table::as.data.table(covariate_blocks)
  structural_membership <- data.table::as.data.table(structural_membership)
  if (!"membership_weight" %in% names(structural_membership)) {
    structural_membership[, membership_weight := 1]
  }
  validate_fssgl_inputs(x, y, covariate_blocks, structural_membership, block_size)
  validate_nonoverlapping_partition(covariate_blocks, structural_membership)
  if (length(theta_prior_covariate) != 2L || any(theta_prior_covariate <= 0) ||
      length(theta_prior_group) != 2L || any(theta_prior_group <= 0)) {
    stop("Theta priors must contain two positive beta-shape parameters.")
  }
  if (length(theta_bounds) != 2L || theta_bounds[1] <= 0 ||
      theta_bounds[2] >= 1 || theta_bounds[1] >= theta_bounds[2]) {
    stop("theta_bounds must lie strictly inside (0, 1).")
  }
  if (theta_damping <= 0 || theta_damping > 1) {
    stop("theta_damping must lie in (0, 1].")
  }
  if (sigma2_damping <= 0 || sigma2_damping > 1) {
    stop("sigma2_damping must lie in (0, 1].")
  }
  if (variance_prior_shape < 0 || variance_prior_scale < 0) {
    stop("Inverse-gamma variance-prior parameters must be nonnegative.")
  }
  if (!is.logical(trace_responsibilities) ||
      length(trace_responsibilities) != 1L || is.na(trace_responsibilities)) {
    stop("trace_responsibilities must be TRUE or FALSE.")
  }
  trace_responsibilities <- isTRUE(trace_responsibilities)

  sample_size <- nrow(x)
  coefficient_count <- ncol(x)
  y <- as.numeric(y)
  beta <- if (is.null(beta_init)) numeric(coefficient_count) else as.numeric(beta_init)
  if (length(beta) != coefficient_count) {
    stop("beta_init has wrong length.")
  }
  group_structure <- build_v2_group_structure(
    covariate_blocks, structural_membership, block_size
  )
  covariate_dimensions <- rep(block_size, nrow(covariate_blocks))
  group_dimensions <- group_structure$group_dimension
  theta_covariate_current <- theta_covariate
  theta_group_current <- theta_group
  initial_residual_variance <- mean((y - mean(y))^2)
  sigma2 <- if (use_variance_update) {
    if (is.null(sigma2_init)) initial_residual_variance else sigma2_init
  } else {
    1
  }
  sigma2 <- max(as.numeric(sigma2), 1e-8)

  largest_eigenvalue <- max(eigen(crossprod(x), symmetric = TRUE, only.values = TRUE)$values)
  lipschitz <- max(lipschitz_multiplier * largest_eigenvalue / sample_size, 1e-8)
  history <- data.table::data.table()
  responsibility_trace_rows <- list()
  responsibility_trace_index <- 1L

  for (outer_iter in seq_len(max_outer_iter)) {
    beta_old <- beta
    sigma2_old <- sigma2
    theta_covariate_old <- theta_covariate_current
    theta_group_old <- theta_group_current

    covariate_norms <- compute_covariate_beta_norm(beta, covariate_blocks)
    group_norms <- compute_v2_group_norms(beta, group_structure)
    covariate_state <- fssgl_v2_slab_state(
      covariate_norms, covariate_dimensions,
      lambda_covariate_spike, lambda_covariate_slab,
      theta_covariate_current, score_mode
    )
    group_state <- fssgl_v2_slab_state(
      group_norms, group_dimensions,
      lambda_group_spike, lambda_group_slab,
      theta_group_current, score_mode
    )
    covariate_penalty <- covariate_penalty_scale * covariate_state$effective_rate
    group_penalty <- group_penalty_scale * group_state$effective_rate

    mstep <- solve_weighted_nested_mstep(
      x = x,
      y = y,
      beta_init = beta,
      variance_scale = sigma2,
      covariate_penalty = covariate_penalty,
      group_penalty = group_penalty,
      covariate_blocks = covariate_blocks,
      group_structure = group_structure,
      lipschitz_init = lipschitz,
      max_iter = inner_max_iter,
      tol = inner_tol
    )
    beta <- mstep$beta
    lipschitz <- mstep$lipschitz
    residual <- y - drop(x %*% beta)

    covariate_norms <- compute_covariate_beta_norm(beta, covariate_blocks)
    group_norms <- compute_v2_group_norms(beta, group_structure)
    covariate_state <- fssgl_v2_slab_state(
      covariate_norms, covariate_dimensions,
      lambda_covariate_spike, lambda_covariate_slab,
      theta_covariate_current, score_mode
    )
    group_state <- fssgl_v2_slab_state(
      group_norms, group_dimensions,
      lambda_group_spike, lambda_group_slab,
      theta_group_current, score_mode
    )

    if (update_theta) {
      target_theta_covariate <-
        (theta_prior_covariate[1] + sum(covariate_state$slab_probability)) /
        (sum(theta_prior_covariate) + length(covariate_state$slab_probability))
      target_theta_group <-
        (theta_prior_group[1] + sum(group_state$slab_probability)) /
        (sum(theta_prior_group) + length(group_state$slab_probability))
      theta_covariate_current <-
        (1 - theta_damping) * theta_covariate_current + theta_damping * target_theta_covariate
      theta_group_current <-
        (1 - theta_damping) * theta_group_current + theta_damping * target_theta_group
      theta_covariate_current <- min(max(theta_covariate_current, theta_bounds[1]), theta_bounds[2])
      theta_group_current <- min(max(theta_group_current, theta_bounds[1]), theta_bounds[2])
    }
    if (use_variance_update) {
      target_sigma2 <- update_fssgl_variance(
        residual,
        sample_size,
        variance_prior_shape,
        variance_prior_scale
      )
      sigma2 <- (1 - sigma2_damping) * sigma2 + sigma2_damping * target_sigma2
    }

    if (trace_responsibilities) {
      trace_covariate_state <- fssgl_v2_slab_state(
        covariate_norms, covariate_dimensions,
        lambda_covariate_spike, lambda_covariate_slab,
        theta_covariate_current, score_mode
      )
      trace_group_state <- fssgl_v2_slab_state(
        group_norms, group_dimensions,
        lambda_group_spike, lambda_group_slab,
        theta_group_current, score_mode
      )
      covariate_transition <- fssgl_v2_responsibility_transition(
        covariate_dimensions, lambda_covariate_spike, lambda_covariate_slab,
        theta_covariate_current, score_mode
      )
      group_transition <- fssgl_v2_responsibility_transition(
        group_dimensions, lambda_group_spike, lambda_group_slab,
        theta_group_current, score_mode
      )
      responsibility_trace_rows[[responsibility_trace_index]] <- data.table::rbindlist(
        list(
          data.table::data.table(
            iter = outer_iter,
            level = "covariate",
            unit_id = as.character(covariate_blocks$covariate_id),
            dimension = covariate_dimensions,
            beta_norm = covariate_norms,
            posterior_slab_prob = trace_covariate_state$slab_probability,
            effective_rate = trace_covariate_state$effective_rate,
            penalty_scale = covariate_penalty_scale,
            penalty_weight = covariate_penalty_scale * trace_covariate_state$effective_rate,
            transition_norm = covariate_transition,
            theta = theta_covariate_current
          ),
          data.table::data.table(
            iter = outer_iter,
            level = "group",
            unit_id = as.character(group_structure$group_id),
            dimension = group_dimensions,
            beta_norm = group_norms,
            posterior_slab_prob = trace_group_state$slab_probability,
            effective_rate = trace_group_state$effective_rate,
            penalty_scale = group_penalty_scale,
            penalty_weight = group_penalty_scale * trace_group_state$effective_rate,
            transition_norm = group_transition,
            theta = theta_group_current
          )
        ),
        use.names = TRUE
      )
      responsibility_trace_index <- responsibility_trace_index + 1L
    }

    beta_change <- sqrt(sum((beta - beta_old)^2)) /
      (sqrt(sum(beta_old^2)) + 1e-8)
    theta_change <- max(
      abs(theta_covariate_current - theta_covariate_old),
      abs(theta_group_current - theta_group_old)
    )
    sigma2_change <- abs(sigma2 - sigma2_old) / (abs(sigma2_old) + 1e-8)
    history <- rbind(
      history,
      data.table::data.table(
        iter = outer_iter,
        working_objective = mstep$objective,
        residual_variance = sum(residual^2) / sample_size,
        sigma2 = sigma2,
        theta_covariate = theta_covariate_current,
        theta_group = theta_group_current,
        beta_change = beta_change,
        theta_change = theta_change,
        sigma2_change = sigma2_change,
        inner_iterations = mstep$iterations,
        inner_change = mstep$final_change,
        backtracking_steps = mstep$backtracking_steps
      )
    )

    if (verbose && (outer_iter == 1L || outer_iter %% 25L == 0L)) {
      message(
        "outer=", outer_iter,
        " beta_change=", signif(beta_change, 4),
        " theta_F=", signif(theta_covariate_current, 4),
        " theta_G=", signif(theta_group_current, 4),
        " sigma2=", signif(sigma2, 4)
      )
    }
    nuisance_converged <- (!update_theta || theta_change < tol) &&
      (!use_variance_update || sigma2_change < tol)
    if (outer_iter > 2L && beta_change < tol && nuisance_converged) {
      break
    }
  }

  covariate_norms <- compute_covariate_beta_norm(beta, covariate_blocks)
  group_norms <- compute_v2_group_norms(beta, group_structure)
  covariate_state <- fssgl_v2_slab_state(
    covariate_norms, covariate_dimensions,
    lambda_covariate_spike, lambda_covariate_slab,
    theta_covariate_current, score_mode
  )
  group_state <- fssgl_v2_slab_state(
    group_norms, group_dimensions,
    lambda_group_spike, lambda_group_slab,
    theta_group_current, score_mode
  )
  covariate_table <- data.table::copy(covariate_blocks)
  covariate_table[, `:=`(
    beta_norm = covariate_norms,
    posterior_slab_prob = covariate_state$slab_probability,
    effective_rate = covariate_state$effective_rate,
    selected = covariate_state$slab_probability >= posterior_cutoff
  )]
  group_table <- data.table::data.table(
    group_id = group_structure$group_id,
    member_count = group_structure$member_count,
    group_dimension = group_dimensions,
    beta_norm = group_norms,
    posterior_slab_prob = group_state$slab_probability,
    effective_rate = group_state$effective_rate
  )
  selected_covariates <- covariate_table[selected == TRUE, covariate_id]
  supported_groups <- unique(structural_membership[
    covariate_id %in% selected_covariates,
    group_id
  ])
  group_table[, `:=`(
    score_selected = posterior_slab_prob >= posterior_cutoff,
    selected = group_id %in% supported_groups
  )]
  data.table::setorder(covariate_table, -posterior_slab_prob, -beta_norm)
  data.table::setorder(group_table, -selected, -posterior_slab_prob, -beta_norm)

  final_iter <- nrow(history)
  objective_tail <- tail(history$working_objective, min(10L, final_iter))
  objective_tail_rel_change <- if (length(objective_tail) > 1L) {
    abs(tail(objective_tail, 1L) - objective_tail[1L]) /
      (abs(objective_tail[1L]) + 1e-8)
  } else {
    NA_real_
  }
  final_beta_change <- history$beta_change[final_iter]
  final_theta_change <- history$theta_change[final_iter]
  final_sigma2_change <- history$sigma2_change[final_iter]
  strict_converged <- final_beta_change < tol &&
    (!update_theta || final_theta_change < tol) &&
    (!use_variance_update || final_sigma2_change < tol)
  relaxed_converged <- final_beta_change < relaxed_tol &&
    objective_tail_rel_change < objective_tail_tol
  responsibility_trace <- if (length(responsibility_trace_rows)) {
    data.table::rbindlist(responsibility_trace_rows, use.names = TRUE)
  } else {
    data.table::data.table(
      iter = integer(), level = character(), unit_id = character(),
      dimension = numeric(), beta_norm = numeric(),
      posterior_slab_prob = numeric(), effective_rate = numeric(),
      penalty_scale = numeric(), penalty_weight = numeric(),
      transition_norm = numeric(), theta = numeric()
    )
  }
  responsibility_trace[, norm_to_transition := beta_norm / pmax(transition_norm, 1e-12)]
  responsibility_trace[, score_selected := posterior_slab_prob >= posterior_cutoff]

  list(
    beta = beta,
    fitted = y - residual,
    residual = residual,
    coefficient_estimates = beta,
    covariate_posterior = covariate_table,
    structure_posterior = group_table,
    structural_membership = structural_membership,
    history = history,
    responsibility_trace = responsibility_trace,
    convergence = list(
      strict_converged = strict_converged,
      relaxed_converged = relaxed_converged,
      final_iter = final_iter,
      final_beta_change = final_beta_change,
      final_theta_change = final_theta_change,
      final_sigma2_change = final_sigma2_change,
      objective_tail_rel_change = objective_tail_rel_change,
      final_residual_variance = sum(residual^2) / sample_size
    ),
    hyperparameters = list(
      score_mode = score_mode,
      use_variance_update = use_variance_update,
      update_theta = update_theta,
      theta_covariate_initial = theta_covariate,
      theta_group_initial = theta_group,
      theta_covariate_final = theta_covariate_current,
      theta_group_final = theta_group_current,
      sigma2_final = sigma2,
      sigma2_damping = sigma2_damping,
      variance_prior_shape = variance_prior_shape,
      variance_prior_scale = variance_prior_scale,
      inner_max_iter = inner_max_iter,
      inner_tol = inner_tol,
      max_outer_iter = max_outer_iter,
      covariate_penalty_scale = covariate_penalty_scale,
      group_penalty_scale = group_penalty_scale,
      posterior_cutoff = posterior_cutoff
    )
  )
}

fit_fssgl_coefficients <- function(
  x_coef,
  y_coef,
  structural_membership,
  ...
) {
  design <- build_fof_design(x_coef, n_response_basis = ncol(y_coef))
  covariate_blocks <- build_functional_groups(design$predictor_blocks)
  fit <- fit_fssgl_core(
    x = design$design,
    y = as.vector(y_coef),
    covariate_blocks = covariate_blocks,
    structural_membership = structural_membership,
    block_size = dim(x_coef)[3] * ncol(y_coef),
    ...
  )
  list(
    fit = fit,
    design = design,
    covariate_blocks = covariate_blocks,
    structural_membership = structural_membership,
    coefficient_estimates = fit$coefficient_estimates,
    covariate_posterior = fit$covariate_posterior,
    structure_posterior = fit$structure_posterior
  )
}

fit_fssgl <- function(
  x_coef,
  y_coef,
  structural_membership,
  beta_init = NULL,
  verbose = FALSE,
  parameters = fssgl_parameters(),
  trace_responsibilities = FALSE
) {
  required_parameters <- names(fssgl_parameters(max_outer_iter = 1L))
  missing_parameters <- setdiff(required_parameters, names(parameters))
  if (length(missing_parameters) > 0L) {
    stop(
      "FSSGL parameters are missing: ",
      paste(missing_parameters, collapse = ", ")
    )
  }
  if (is.null(beta_init)) {
    warm_design <- build_fof_design(x_coef, n_response_basis = ncol(y_coef))
    beta_init <- ridge_dual_fit(
      warm_design$design,
      as.vector(y_coef),
      lambda = parameters$ridge_lambda
    )
  }
  solver_parameters <- setdiff(required_parameters, "ridge_lambda")
  do.call(
    fit_fssgl_coefficients,
    c(
      list(
        x_coef = x_coef,
        y_coef = y_coef,
        structural_membership = structural_membership
      ),
      parameters[solver_parameters],
      list(
        beta_init = beta_init,
        trace_responsibilities = trace_responsibilities,
        verbose = verbose
      )
    )
  )
}

# Continue a nonconverged fit from its complete terminal state. This preserves
# the objective and update equations; it only increases numerical effort. The
# returned history records all stages so tuning can reject unresolved fits.
fit_fssgl_until_converged <- function(
  x_coef,
  y_coef,
  structural_membership,
  beta_init = NULL,
  verbose = FALSE,
  parameters = fssgl_parameters(),
  extension_factors = c(1L, 2L),
  trace_responsibilities = FALSE
) {
  extension_factors <- as.integer(extension_factors)
  if (!length(extension_factors) || extension_factors[1L] != 1L ||
      any(extension_factors < 1L) || is.unsorted(extension_factors)) {
    stop("extension_factors must be nondecreasing positive integers beginning with one.")
  }

  base_max_outer <- as.integer(parameters$max_outer_iter)
  stage_histories <- list()
  stage_responsibility_traces <- list()
  warm_beta <- beta_init
  fit <- NULL

  for (stage in seq_along(extension_factors)) {
    stage_parameters <- parameters
    stage_parameters$max_outer_iter <- base_max_outer * extension_factors[stage]
    fit <- fit_fssgl(
      x_coef = x_coef,
      y_coef = y_coef,
      structural_membership = structural_membership,
      beta_init = warm_beta,
      verbose = verbose,
      parameters = stage_parameters,
      trace_responsibilities = trace_responsibilities
    )
    stage_history <- data.table::copy(fit$fit$history)
    stage_history[, `:=`(stage = stage, stage_iter = iter)]
    stage_histories[[stage]] <- stage_history
    if (trace_responsibilities) {
      stage_trace <- data.table::copy(fit$fit$responsibility_trace)
      stage_trace[, `:=`(stage = stage, stage_iter = iter)]
      stage_responsibility_traces[[stage]] <- stage_trace
    }

    if (isTRUE(fit$fit$convergence$strict_converged)) break

    warm_beta <- fit$fit$beta
    parameters$theta_covariate <- fit$fit$hyperparameters$theta_covariate_final
    parameters$theta_group <- fit$fit$hyperparameters$theta_group_final
    parameters$sigma2_init <- fit$fit$hyperparameters$sigma2_final
  }

  combined_history <- data.table::rbindlist(stage_histories, use.names = TRUE, fill = TRUE)
  combined_history[, iter := seq_len(.N)]
  fit$fit$history <- combined_history
  if (trace_responsibilities) {
    combined_trace <- data.table::rbindlist(
      stage_responsibility_traces, use.names = TRUE, fill = TRUE
    )
    stage_offsets <- c(0L, cumsum(vapply(
      stage_histories, nrow, integer(1L)
    )))[seq_along(stage_histories)]
    combined_trace[, iter := stage_iter + stage_offsets[stage]]
    fit$fit$responsibility_trace <- combined_trace
  }
  fit$fit$convergence$stage_final_iter <- fit$fit$convergence$final_iter
  fit$fit$convergence$total_outer_iter <- nrow(combined_history)
  fit$fit$convergence$extension_stages <- length(stage_histories)
  fit$fit$convergence$max_outer_iter_schedule <-
    base_max_outer * extension_factors[seq_along(stage_histories)]
  fit$fit$convergence$final_iter <- nrow(combined_history)
  fit$fit$hyperparameters$max_outer_iter_schedule <-
    fit$fit$convergence$max_outer_iter_schedule
  fit
}
