# Baseline estimators for the function-on-function simulation study.
# Inputs use the same basis-coefficient representation as FSSGL.

flatten_predictor_coefficients <- function(x_coef) {
  if (length(dim(x_coef)) != 3L) {
    stop("x_coef must be a covariate x sample x basis array.")
  }
  do.call(cbind, lapply(seq_len(dim(x_coef)[1L]), function(j) {
    x_coef[j, , , drop = FALSE][1L, , ]
  }))
}

coefficient_matrix_to_block_vector <- function(coef_matrix, n_covariates, kx) {
  unlist(lapply(seq_len(n_covariates), function(j) {
    rows <- ((j - 1L) * kx + 1L):(j * kx)
    as.vector(coef_matrix[rows, , drop = FALSE])
  }), use.names = FALSE)
}

make_baseline_cv_folds <- function(n, n_folds = 5L, seed = 1L) {
  if (!exists("make_curve_cv_folds", mode = "function")) {
    stop("Source R/fssgl/tuning.R before constructing comparison folds.")
  }
  make_curve_cv_folds(n = n, n_folds = n_folds, seed = seed)
}

standardize_matrix <- function(x) {
  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  list(
    values = sweep(sweep(x, 2L, center, "-"), 2L, scale, "/"),
    center = center,
    scale = scale
  )
}

fit_multivariate_ridge <- function(x, y, lambda) {
  prep <- standardize_matrix(x)
  y_center <- colMeans(y)
  yc <- sweep(y, 2L, y_center, "-")
  q <- ncol(prep$values)
  coef_scaled <- solve(
    crossprod(prep$values) + as.numeric(lambda) * diag(q),
    crossprod(prep$values, yc)
  )
  coef <- sweep(coef_scaled, 1L, prep$scale, "/")
  intercept <- y_center - drop(prep$center %*% coef)
  list(coef = coef, intercept = intercept, lambda = lambda)
}

predict_multivariate_linear <- function(model, x) {
  sweep(x %*% model$coef, 2L, model$intercept, "+")
}

cv_multivariate_ridge <- function(x, y, folds, lambda_grid) {
  losses <- matrix(NA_real_, nrow = length(folds), ncol = length(lambda_grid))
  for (fold_id in seq_along(folds)) {
    valid <- folds[[fold_id]]
    train <- setdiff(seq_len(nrow(x)), valid)
    for (lambda_id in seq_along(lambda_grid)) {
      model <- fit_multivariate_ridge(x[train, , drop = FALSE], y[train, , drop = FALSE], lambda_grid[lambda_id])
      prediction <- predict_multivariate_linear(model, x[valid, , drop = FALSE])
      losses[fold_id, lambda_id] <- sqrt(mean((prediction - y[valid, , drop = FALSE])^2))
    }
  }
  mean_loss <- colMeans(losses)
  best <- which.min(mean_loss)
  list(
    selected_lambda = lambda_grid[best],
    cv_mean = mean_loss,
    cv_sd = apply(losses, 2L, stats::sd),
    losses = losses
  )
}

fit_basis_ridge_baseline <- function(x_coef, y_coef, folds, lambda_grid) {
  x <- flatten_predictor_coefficients(x_coef)
  tuning_started <- proc.time()[["elapsed"]]
  tuning <- cv_multivariate_ridge(x, y_coef, folds, lambda_grid)
  tuning_runtime_sec <- proc.time()[["elapsed"]] - tuning_started
  final_started <- proc.time()[["elapsed"]]
  model <- fit_multivariate_ridge(x, y_coef, tuning$selected_lambda)
  final_runtime_sec <- proc.time()[["elapsed"]] - final_started
  p <- dim(x_coef)[1L]
  kx <- dim(x_coef)[3L]
  list(
    method_id = "basis_ridge",
    supports_selection = FALSE,
    supports_coefficients = TRUE,
    model = model,
    beta = coefficient_matrix_to_block_vector(model$coef, p, kx),
    selected_covariates = NULL,
    tuning = tuning,
    timing = list(
      tuning_runtime_sec = tuning_runtime_sec,
      final_runtime_sec = final_runtime_sec
    )
  )
}

fit_pca_projection <- function(x, variance_threshold = 0.95, max_components = NULL) {
  center <- colMeans(x)
  xc <- sweep(x, 2L, center, "-")
  decomposition <- svd(xc, nu = 0L)
  variance <- decomposition$d^2
  if (!any(is.finite(variance)) || sum(variance) <= 1e-12) {
    rotation <- diag(ncol(x))[, 1L, drop = FALSE]
  } else {
    cumulative <- cumsum(variance) / sum(variance)
    component_count <- which(cumulative >= variance_threshold)[1L]
    if (is.null(max_components)) max_components <- ncol(x)
    component_count <- max(1L, min(component_count, max_components, ncol(decomposition$v)))
    rotation <- decomposition$v[, seq_len(component_count), drop = FALSE]
  }
  list(center = center, rotation = rotation)
}

fit_fpca_representation <- function(
  x_coef,
  y_coef,
  variance_threshold = 0.95,
  max_x_components = NULL,
  max_y_components = NULL
) {
  p <- dim(x_coef)[1L]
  x_models <- vector("list", p)
  x_scores <- vector("list", p)
  groups <- vector("list", p)
  column_start <- 1L

  for (j in seq_len(p)) {
    xj <- x_coef[j, , , drop = FALSE][1L, , ]
    x_models[[j]] <- fit_pca_projection(xj, variance_threshold, max_x_components)
    x_scores[[j]] <- sweep(xj, 2L, x_models[[j]]$center, "-") %*% x_models[[j]]$rotation
    groups[[j]] <- column_start:(column_start + ncol(x_scores[[j]]) - 1L)
    column_start <- max(groups[[j]]) + 1L
  }

  y_model <- fit_pca_projection(y_coef, variance_threshold, max_y_components)
  y_scores <- sweep(y_coef, 2L, y_model$center, "-") %*% y_model$rotation
  list(
    x_models = x_models,
    y_model = y_model,
    x_scores = do.call(cbind, x_scores),
    y_scores = y_scores,
    groups = groups
  )
}

transform_fpca_predictors <- function(x_coef, representation) {
  do.call(cbind, lapply(seq_along(representation$x_models), function(j) {
    xj <- x_coef[j, , , drop = FALSE][1L, , ]
    sweep(xj, 2L, representation$x_models[[j]]$center, "-") %*%
      representation$x_models[[j]]$rotation
  }))
}

reconstruct_fpca_coefficients <- function(score_coef, representation, ky) {
  p <- length(representation$x_models)
  kx <- length(representation$x_models[[1L]]$center)
  beta <- numeric(p * kx * ky)
  surfaces <- vector("list", p)
  for (j in seq_len(p)) {
    score_rows <- representation$groups[[j]]
    surface <- representation$x_models[[j]]$rotation %*%
      score_coef[score_rows, , drop = FALSE] %*%
      t(representation$y_model$rotation)
    surfaces[[j]] <- surface
    block <- ((j - 1L) * kx * ky + 1L):(j * kx * ky)
    beta[block] <- as.vector(surface)
  }
  list(beta = beta, surfaces = surfaces)
}

fit_fpca_ridge_baseline <- function(
  x_coef,
  y_coef,
  folds,
  lambda_grid,
  variance_threshold = 0.95
) {
  tuning_started <- proc.time()[["elapsed"]]
  losses <- matrix(NA_real_, nrow = length(folds), ncol = length(lambda_grid))
  n <- dim(x_coef)[2L]
  for (fold_id in seq_along(folds)) {
    valid <- folds[[fold_id]]
    train <- setdiff(seq_len(n), valid)
    representation <- fit_fpca_representation(
      x_coef[, train, , drop = FALSE],
      y_coef[train, , drop = FALSE],
      variance_threshold = variance_threshold
    )
    valid_scores <- transform_fpca_predictors(x_coef[, valid, , drop = FALSE], representation)
    for (lambda_id in seq_along(lambda_grid)) {
      model <- fit_multivariate_ridge(
        representation$x_scores,
        representation$y_scores,
        lambda_grid[lambda_id]
      )
      predicted_scores <- predict_multivariate_linear(model, valid_scores)
      prediction <- sweep(
        predicted_scores %*% t(representation$y_model$rotation),
        2L,
        representation$y_model$center,
        "+"
      )
      losses[fold_id, lambda_id] <- sqrt(mean((prediction - y_coef[valid, , drop = FALSE])^2))
    }
  }
  mean_loss <- colMeans(losses)
  best <- which.min(mean_loss)
  tuning <- list(
    selected_lambda = lambda_grid[best],
    cv_mean = mean_loss,
    cv_sd = apply(losses, 2L, stats::sd),
    losses = losses
  )
  tuning_runtime_sec <- proc.time()[["elapsed"]] - tuning_started

  final_started <- proc.time()[["elapsed"]]
  representation <- fit_fpca_representation(x_coef, y_coef, variance_threshold)
  model <- fit_multivariate_ridge(
    representation$x_scores,
    representation$y_scores,
    tuning$selected_lambda
  )
  reconstructed <- reconstruct_fpca_coefficients(model$coef, representation, ncol(y_coef))
  final_runtime_sec <- proc.time()[["elapsed"]] - final_started
  list(
    method_id = "fpca_ridge",
    supports_selection = FALSE,
    supports_coefficients = TRUE,
    representation = representation,
    model = model,
    beta = reconstructed$beta,
    selected_covariates = NULL,
    tuning = tuning,
    timing = list(
      tuning_runtime_sec = tuning_runtime_sec,
      final_runtime_sec = final_runtime_sec
    )
  )
}

group_lasso_path <- function(
  x,
  y,
  groups,
  lambda_fractions,
  max_iter = 2000L,
  tol = 1e-6
) {
  prep <- standardize_matrix(x)
  y_center <- colMeans(y)
  yc <- sweep(y, 2L, y_center, "-")
  n <- nrow(x)
  response_dim <- ncol(y)
  group_weights <- vapply(groups, function(index) sqrt(length(index) * response_dim), numeric(1))
  gradient_norms <- vapply(seq_along(groups), function(g) {
    norm(crossprod(prep$values[, groups[[g]], drop = FALSE], yc) / n, type = "F") /
      group_weights[g]
  }, numeric(1))
  lambda_max <- max(gradient_norms)
  if (!is.finite(lambda_max) || lambda_max <= 0) lambda_max <- 1

  eigenvalues <- eigen(crossprod(prep$values) / n, symmetric = TRUE, only.values = TRUE)$values
  lipschitz <- max(eigenvalues, na.rm = TRUE)
  if (!is.finite(lipschitz) || lipschitz <= 0) lipschitz <- 1

  fractions <- sort(unique(lambda_fractions), decreasing = TRUE)
  coefficient_scaled <- matrix(0, nrow = ncol(x), ncol = response_dim)
  fits <- vector("list", length(fractions))

  for (lambda_id in seq_along(fractions)) {
    lambda <- fractions[lambda_id] * lambda_max
    accelerated <- coefficient_scaled
    momentum <- 1
    converged <- FALSE

    for (iter in seq_len(max_iter)) {
      old <- coefficient_scaled
      gradient <- crossprod(prep$values, prep$values %*% accelerated - yc) / n
      candidate <- accelerated - gradient / lipschitz
      for (g in seq_along(groups)) {
        rows <- groups[[g]]
        block <- candidate[rows, , drop = FALSE]
        block_norm <- norm(block, type = "F")
        threshold <- lambda * group_weights[g] / lipschitz
        candidate[rows, ] <- if (!is.finite(block_norm) || block_norm <= threshold) {
          0
        } else {
          (1 - threshold / block_norm) * block
        }
      }

      new_momentum <- (1 + sqrt(1 + 4 * momentum^2)) / 2
      accelerated <- candidate + ((momentum - 1) / new_momentum) * (candidate - old)
      coefficient_scaled <- candidate
      momentum <- new_momentum
      relative_change <- norm(coefficient_scaled - old, type = "F") /
        (norm(old, type = "F") + 1e-8)
      if (iter > 2L && relative_change < tol) {
        converged <- TRUE
        break
      }
    }

    coefficient <- sweep(coefficient_scaled, 1L, prep$scale, "/")
    intercept <- y_center - drop(prep$center %*% coefficient)
    selected <- which(vapply(groups, function(rows) {
      norm(coefficient[rows, , drop = FALSE], type = "F") > 1e-8
    }, logical(1)))
    fits[[lambda_id]] <- list(
      coef = coefficient,
      intercept = intercept,
      lambda = lambda,
      lambda_fraction = fractions[lambda_id],
      selected_groups = selected,
      converged = converged,
      final_iter = iter
    )
  }
  names(fits) <- format(fractions, scientific = FALSE, trim = TRUE)
  list(fractions = fractions, lambda_max = lambda_max, fits = fits)
}

select_one_se_fraction <- function(losses, fractions) {
  means <- colMeans(losses)
  standard_errors <- apply(losses, 2L, stats::sd) / sqrt(nrow(losses))
  minimum <- which.min(means)
  eligible <- which(means <= means[minimum] + standard_errors[minimum])
  selected <- eligible[which.max(fractions[eligible])]
  list(
    selected_index = selected,
    selected_fraction = fractions[selected],
    cv_mean = means,
    cv_se = standard_errors,
    losses = losses
  )
}

structured_group_lasso_path <- function(
  x,
  y,
  predictor_groups,
  parent_groups,
  alpha,
  lambda_fractions,
  max_iter = 2000L,
  tol = 1e-6
) {
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between zero and one.")
  }
  x_center <- colMeans(x)
  xc <- sweep(x, 2L, x_center, "-")
  y_center <- colMeans(y)
  yc <- sweep(y, 2L, y_center, "-")
  n <- nrow(xc)
  response_dim <- ncol(yc)
  predictor_weights <- vapply(
    predictor_groups,
    function(rows) sqrt(length(rows) * response_dim),
    numeric(1L)
  )
  parent_weights <- vapply(
    parent_groups,
    function(rows) sqrt(length(rows) * response_dim),
    numeric(1L)
  )
  gradient_zero <- crossprod(xc, yc) / n
  predictor_bound <- max(vapply(seq_along(predictor_groups), function(j) {
    norm(gradient_zero[predictor_groups[[j]], , drop = FALSE], type = "F") /
      (alpha * predictor_weights[j])
  }, numeric(1L)))
  parent_bound <- max(vapply(seq_along(parent_groups), function(g) {
    norm(gradient_zero[parent_groups[[g]], , drop = FALSE], type = "F") /
      ((1 - alpha) * parent_weights[g])
  }, numeric(1L)))
  lambda_max <- max(predictor_bound, parent_bound)
  if (!is.finite(lambda_max) || lambda_max <= 0) lambda_max <- 1

  largest_eigenvalue <- max(
    eigen(crossprod(xc) / n, symmetric = TRUE, only.values = TRUE)$values,
    na.rm = TRUE
  )
  lipschitz <- if (is.finite(largest_eigenvalue) && largest_eigenvalue > 0) {
    1.01 * largest_eigenvalue
  } else {
    1
  }
  fractions <- sort(unique(lambda_fractions), decreasing = TRUE)
  coefficient <- matrix(0, nrow = ncol(xc), ncol = response_dim)
  fits <- vector("list", length(fractions))

  for (lambda_id in seq_along(fractions)) {
    lambda <- fractions[lambda_id] * lambda_max
    accelerated <- coefficient
    momentum <- 1
    converged <- FALSE
    for (iter in seq_len(max_iter)) {
      old <- coefficient
      gradient <- crossprod(xc, xc %*% accelerated - yc) / n
      candidate <- accelerated - gradient / lipschitz

      # Predictor blocks are nested within disjoint parent groups, so the
      # proximal map of the tree penalty is the child-before-parent composition.
      for (j in seq_along(predictor_groups)) {
        rows <- predictor_groups[[j]]
        block <- candidate[rows, , drop = FALSE]
        block_norm <- norm(block, type = "F")
        threshold <- lambda * alpha * predictor_weights[j] / lipschitz
        candidate[rows, ] <- if (!is.finite(block_norm) || block_norm <= threshold) {
          0
        } else {
          (1 - threshold / block_norm) * block
        }
      }
      for (g in seq_along(parent_groups)) {
        rows <- parent_groups[[g]]
        block <- candidate[rows, , drop = FALSE]
        block_norm <- norm(block, type = "F")
        threshold <- lambda * (1 - alpha) * parent_weights[g] / lipschitz
        candidate[rows, ] <- if (!is.finite(block_norm) || block_norm <= threshold) {
          0
        } else {
          (1 - threshold / block_norm) * block
        }
      }

      new_momentum <- (1 + sqrt(1 + 4 * momentum^2)) / 2
      accelerated <- candidate + ((momentum - 1) / new_momentum) * (candidate - old)
      coefficient <- candidate
      momentum <- new_momentum
      relative_change <- norm(coefficient - old, type = "F") /
        (norm(old, type = "F") + 1e-8)
      if (iter > 2L && relative_change < tol) {
        converged <- TRUE
        break
      }
    }
    selected <- which(vapply(predictor_groups, function(rows) {
      norm(coefficient[rows, , drop = FALSE], type = "F") > 1e-8
    }, logical(1L)))
    fits[[lambda_id]] <- list(
      coef = coefficient,
      intercept = y_center - drop(x_center %*% coefficient),
      lambda = lambda,
      lambda_fraction = fractions[lambda_id],
      alpha = alpha,
      selected_predictors = selected,
      converged = converged,
      final_iter = iter
    )
  }
  list(fractions = fractions, lambda_max = lambda_max, fits = fits)
}

build_structured_baseline_groups <- function(x_coef, structural_membership) {
  p <- dim(x_coef)[1L]
  kx <- dim(x_coef)[3L]
  membership <- data.table::as.data.table(structural_membership)
  counts <- membership[, .N, by = covariate_id]
  if (!setequal(counts$covariate_id, seq_len(p)) || any(counts$N != 1L)) {
    stop("The structured group-lasso baseline requires one parent group per predictor.")
  }
  predictor_groups <- lapply(seq_len(p), function(j) {
    ((j - 1L) * kx + 1L):(j * kx)
  })
  parent_ids <- sort(unique(membership$group_id))
  parent_groups <- lapply(parent_ids, function(current_group_id) {
    members <- membership[group_id == current_group_id, covariate_id]
    sort(unlist(predictor_groups[members], use.names = FALSE))
  })
  list(predictor = predictor_groups, parent = parent_groups, parent_ids = parent_ids)
}

fit_structured_group_lasso_baseline <- function(
  x_coef,
  y_coef,
  structural_membership,
  folds,
  alpha_grid,
  lambda_fractions,
  max_iter = 2000L,
  tol = 1e-6
) {
  tuning_started <- proc.time()[["elapsed"]]
  x <- flatten_predictor_coefficients(x_coef)
  groups <- build_structured_baseline_groups(x_coef, structural_membership)
  alphas <- sort(unique(alpha_grid))
  fractions <- sort(unique(lambda_fractions), decreasing = TRUE)
  candidates <- expand.grid(
    alpha = alphas,
    lambda_fraction = fractions,
    KEEP.OUT.ATTRS = FALSE
  )
  losses <- matrix(NA_real_, nrow = length(folds), ncol = nrow(candidates))
  selected_counts <- matrix(NA_real_, nrow = length(folds), ncol = nrow(candidates))

  for (fold_id in seq_along(folds)) {
    valid <- folds[[fold_id]]
    train <- setdiff(seq_len(nrow(x)), valid)
    for (alpha_id in seq_along(alphas)) {
      path <- structured_group_lasso_path(
        x[train, , drop = FALSE], y_coef[train, , drop = FALSE],
        groups$predictor, groups$parent, alphas[alpha_id], fractions,
        max_iter, tol
      )
      for (fraction_id in seq_along(fractions)) {
        candidate_id <- which(
          candidates$alpha == alphas[alpha_id] &
            candidates$lambda_fraction == fractions[fraction_id]
        )
        model <- path$fits[[fraction_id]]
        prediction <- predict_multivariate_linear(model, x[valid, , drop = FALSE])
        losses[fold_id, candidate_id] <- sqrt(mean(
          (prediction - y_coef[valid, , drop = FALSE])^2
        ))
        selected_counts[fold_id, candidate_id] <- length(model$selected_predictors)
      }
    }
  }

  cv_mean <- colMeans(losses)
  cv_se <- apply(losses, 2L, stats::sd) / sqrt(nrow(losses))
  minimum <- which.min(cv_mean)
  eligible <- which(cv_mean <= cv_mean[minimum] + cv_se[minimum])
  mean_selected <- colMeans(selected_counts)
  sparsest <- eligible[mean_selected[eligible] == min(mean_selected[eligible])]
  selected_id <- sparsest[order(
    -candidates$lambda_fraction[sparsest],
    -candidates$alpha[sparsest],
    cv_mean[sparsest]
  )[1L]]
  tuning_runtime_sec <- proc.time()[["elapsed"]] - tuning_started

  final_started <- proc.time()[["elapsed"]]
  path <- structured_group_lasso_path(
    x, y_coef, groups$predictor, groups$parent,
    candidates$alpha[selected_id], candidates$lambda_fraction[selected_id],
    max_iter, tol
  )
  model <- path$fits[[1L]]
  final_runtime_sec <- proc.time()[["elapsed"]] - final_started
  list(
    method_id = "structured_group_lasso",
    supports_selection = TRUE,
    supports_coefficients = TRUE,
    model = model,
    beta = coefficient_matrix_to_block_vector(
      model$coef, dim(x_coef)[1L], dim(x_coef)[3L]
    ),
    selected_covariates = model$selected_predictors,
    tuning = list(
      selected_alpha = candidates$alpha[selected_id],
      selected_lambda = model$lambda,
      selected_fraction = candidates$lambda_fraction[selected_id],
      candidates = candidates,
      cv_mean = cv_mean,
      cv_se = cv_se,
      mean_selected_covariates = mean_selected,
      losses = losses
    ),
    timing = list(
      tuning_runtime_sec = tuning_runtime_sec,
      final_runtime_sec = final_runtime_sec
    ),
    solver_converged = model$converged,
    final_iter = model$final_iter
  )
}

fit_fpca_group_lasso_baseline <- function(
  x_coef,
  y_coef,
  folds,
  lambda_fractions,
  variance_threshold = 0.95,
  max_iter = 2000L,
  tol = 1e-6
) {
  tuning_started <- proc.time()[["elapsed"]]
  fractions <- sort(unique(lambda_fractions), decreasing = TRUE)
  losses <- matrix(NA_real_, nrow = length(folds), ncol = length(fractions))
  n <- dim(x_coef)[2L]

  for (fold_id in seq_along(folds)) {
    valid <- folds[[fold_id]]
    train <- setdiff(seq_len(n), valid)
    representation <- fit_fpca_representation(
      x_coef[, train, , drop = FALSE],
      y_coef[train, , drop = FALSE],
      variance_threshold
    )
    path <- group_lasso_path(
      representation$x_scores,
      representation$y_scores,
      representation$groups,
      fractions,
      max_iter,
      tol
    )
    valid_scores <- transform_fpca_predictors(x_coef[, valid, , drop = FALSE], representation)
    for (fraction_id in seq_along(fractions)) {
      predicted_scores <- predict_multivariate_linear(path$fits[[fraction_id]], valid_scores)
      prediction <- sweep(
        predicted_scores %*% t(representation$y_model$rotation),
        2L,
        representation$y_model$center,
        "+"
      )
      losses[fold_id, fraction_id] <- sqrt(mean((prediction - y_coef[valid, , drop = FALSE])^2))
    }
  }

  tuning <- select_one_se_fraction(losses, fractions)
  tuning_runtime_sec <- proc.time()[["elapsed"]] - tuning_started
  final_started <- proc.time()[["elapsed"]]
  representation <- fit_fpca_representation(x_coef, y_coef, variance_threshold)
  path <- group_lasso_path(
    representation$x_scores,
    representation$y_scores,
    representation$groups,
    tuning$selected_fraction,
    max_iter,
    tol
  )
  model <- path$fits[[1L]]
  reconstructed <- reconstruct_fpca_coefficients(model$coef, representation, ncol(y_coef))
  final_runtime_sec <- proc.time()[["elapsed"]] - final_started
  list(
    method_id = "fpca_group_lasso",
    supports_selection = TRUE,
    supports_coefficients = TRUE,
    representation = representation,
    model = model,
    beta = reconstructed$beta,
    selected_covariates = model$selected_groups,
    tuning = tuning,
    timing = list(
      tuning_runtime_sec = tuning_runtime_sec,
      final_runtime_sec = final_runtime_sec
    ),
    solver_converged = model$converged,
    final_iter = model$final_iter
  )
}

# Map a common lambda/lambda_max grid to the path that grpreg actually returns.
# grpreg may stop a path early for a saturated fit, so folds cannot be required
# to have identical path lengths or addressed by the same raw column index.
match_grpreg_fraction_path <- function(fit, target_fractions) {
  lambda <- as.numeric(fit$lambda)
  if (!length(lambda) || any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("grpreg returned an empty or invalid penalty path.")
  }
  actual_fractions <- lambda / lambda[[1L]]
  vapply(target_fractions, function(target) {
    which.min(abs(log(actual_fractions) - log(target)))
  }, integer(1L))
}

# Group SCAD on the same orthonormal B-spline coefficient representation used by
# FSSGL. Each functional predictor is one group of K_X rows across all K_Y
# response coefficients. grpreg applies its standard group preprocessing after
# receiving this common representation.
fit_basis_group_scad_baseline <- function(
  x_coef,
  y_coef,
  folds,
  nlambda = 50L,
  lambda_min = 0.01,
  gamma = 4
) {
  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop("Package grpreg is required for the B-spline group-SCAD baseline.")
  }
  x <- flatten_predictor_coefficients(x_coef)
  p <- dim(x_coef)[1L]
  kx <- dim(x_coef)[3L]
  groups <- split(seq_len(ncol(x)), rep(seq_len(p), each = kx))
  group <- rep(seq_len(p), each = kx)
  nlambda <- as.integer(nlambda)
  fractions <- exp(seq(0, log(lambda_min), length.out = nlambda))
  losses <- matrix(NA_real_, nrow = length(folds), ncol = nlambda)
  selected_counts <- matrix(NA_real_, nrow = length(folds), ncol = nlambda)

  tuning_started <- proc.time()[["elapsed"]]
  for (fold_id in seq_along(folds)) {
    valid <- folds[[fold_id]]
    train <- setdiff(seq_len(nrow(x)), valid)
    fold_path <- grpreg::grpreg(
      X = x[train, , drop = FALSE],
      y = y_coef[train, , drop = FALSE],
      group = group,
      penalty = "grSCAD",
      nlambda = nlambda,
      lambda.min = lambda_min,
      gamma = gamma,
      max.iter = 10000L,
      eps = 1e-5
    )
    fold_path_ids <- match_grpreg_fraction_path(fold_path, fractions)
    predicted <- stats::predict(fold_path, X = x[valid, , drop = FALSE])
    coefficient_path <- stats::coef(fold_path)
    for (lambda_id in seq_len(nlambda)) {
      path_id <- fold_path_ids[[lambda_id]]
      prediction <- predicted[, , path_id, drop = FALSE][, , 1L]
      losses[fold_id, lambda_id] <- sqrt(mean(
        (prediction - y_coef[valid, , drop = FALSE])^2
      ))
      coefficient <- t(coefficient_path[, -1L, path_id, drop = FALSE][, , 1L])
      selected_counts[fold_id, lambda_id] <- sum(vapply(
        groups,
        function(rows) sqrt(sum(coefficient[rows, , drop = FALSE]^2)) > 1e-8,
        logical(1L)
      ))
    }
  }

  cv_mean <- colMeans(losses)
  cv_se <- apply(losses, 2L, stats::sd) / sqrt(nrow(losses))
  minimum <- which.min(cv_mean)
  eligible <- which(cv_mean <= cv_mean[minimum] + cv_se[minimum])
  mean_selected <- colMeans(selected_counts)
  sparsest <- eligible[mean_selected[eligible] == min(mean_selected[eligible])]
  selected_id <- sparsest[which.max(fractions[sparsest])]
  tuning_runtime_sec <- proc.time()[["elapsed"]] - tuning_started

  final_started <- proc.time()[["elapsed"]]
  final_fit <- grpreg::grpreg(
    X = x,
    y = y_coef,
    group = group,
    penalty = "grSCAD",
    nlambda = nlambda,
    lambda.min = lambda_min,
    gamma = gamma,
    max.iter = 10000L,
    eps = 1e-5
  )
  final_path_ids <- match_grpreg_fraction_path(final_fit, fractions)
  selected_path_id <- final_path_ids[[selected_id]]
  coefficient_with_intercept <- stats::coef(final_fit)[
    , , selected_path_id, drop = FALSE
  ][, , 1L]
  if (!is.matrix(coefficient_with_intercept)) {
    coefficient_with_intercept <- matrix(coefficient_with_intercept, nrow = 1L)
  }
  coefficient <- t(coefficient_with_intercept[, -1L, drop = FALSE])
  intercept <- coefficient_with_intercept[, 1L]
  selected_covariates <- which(vapply(groups, function(rows) {
    sqrt(sum(coefficient[rows, , drop = FALSE]^2)) > 1e-8
  }, logical(1L)))
  final_runtime_sec <- proc.time()[["elapsed"]] - final_started

  list(
    method_id = "basis_group_scad",
    supports_selection = TRUE,
    supports_coefficients = TRUE,
    basis_geometry = "orthonormal_bspline_coefficients",
    model = list(coef = coefficient, intercept = intercept),
    beta = coefficient_matrix_to_block_vector(coefficient, p, kx),
    selected_covariates = selected_covariates,
    tuning = list(
      selected_lambda = final_fit$lambda[selected_path_id],
      selected_fraction = fractions[selected_id],
      selected_path_fraction = final_fit$lambda[selected_path_id] /
        final_fit$lambda[[1L]],
      returned_nlambda = length(final_fit$lambda),
      gamma = gamma,
      nlambda = nlambda,
      lambda_min_fraction = lambda_min,
      cv_mean = cv_mean,
      cv_se = cv_se,
      mean_selected_covariates = mean_selected,
      losses = losses,
      lambda = final_fit$lambda
    ),
    timing = list(
      tuning_runtime_sec = tuning_runtime_sec,
      final_runtime_sec = final_runtime_sec
    ),
    solver_converged = TRUE,
    final_iter = final_fit$iter[selected_path_id]
  )
}

# FPCA plus grouped SCAD following Cai, Xue and Cao (2022). Each functional
# predictor is one SCAD group and all response-score equations are fitted
# jointly by grpreg. FPCA is re-estimated within every supplied fold, while
# validation is evaluated in the common response-coefficient geometry.
fit_fpca_group_scad_baseline <- function(
  x_coef,
  y_coef,
  folds,
  variance_threshold = 0.95,
  nlambda = 50L,
  lambda_min = 0.01,
  gamma = 4
) {
  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop("Package grpreg is required for the FPCA group-SCAD baseline.")
  }
  tuning_started <- proc.time()[["elapsed"]]
  n <- dim(x_coef)[2L]
  nlambda <- as.integer(nlambda)
  fractions <- exp(seq(0, log(lambda_min), length.out = nlambda))
  losses <- matrix(NA_real_, nrow = length(folds), ncol = nlambda)
  selected_counts <- matrix(NA_real_, nrow = length(folds), ncol = nlambda)

  for (fold_id in seq_along(folds)) {
    valid <- folds[[fold_id]]
    train <- setdiff(seq_len(n), valid)
    fold_representation <- fit_fpca_representation(
      x_coef[, train, , drop = FALSE],
      y_coef[train, , drop = FALSE],
      variance_threshold = variance_threshold
    )
    fold_group <- integer(ncol(fold_representation$x_scores))
    for (j in seq_along(fold_representation$groups)) {
      fold_group[fold_representation$groups[[j]]] <- j
    }
    fold_path <- grpreg::grpreg(
      X = fold_representation$x_scores,
      y = fold_representation$y_scores,
      group = fold_group,
      penalty = "grSCAD",
      nlambda = nlambda,
      lambda.min = lambda_min,
      gamma = gamma,
      max.iter = 10000L,
      eps = 1e-5
    )
    fold_path_ids <- match_grpreg_fraction_path(fold_path, fractions)
    valid_scores <- transform_fpca_predictors(
      x_coef[, valid, , drop = FALSE], fold_representation
    )
    predicted_scores <- stats::predict(fold_path, X = valid_scores)
    coefficient_path <- stats::coef(fold_path)
    for (lambda_id in seq_len(nlambda)) {
      path_id <- fold_path_ids[[lambda_id]]
      prediction <- sweep(
        predicted_scores[, , path_id, drop = FALSE][, , 1L] %*%
          t(fold_representation$y_model$rotation),
        2L,
        fold_representation$y_model$center,
        "+"
      )
      losses[fold_id, lambda_id] <- sqrt(mean(
        (prediction - y_coef[valid, , drop = FALSE])^2
      ))
      coefficient <- t(coefficient_path[, -1L, path_id, drop = FALSE][, , 1L])
      selected_counts[fold_id, lambda_id] <- sum(vapply(
        fold_representation$groups,
        function(rows) sqrt(sum(coefficient[rows, , drop = FALSE]^2)) > 1e-8,
        logical(1L)
      ))
    }
  }

  cv_mean <- colMeans(losses)
  cv_se <- apply(losses, 2L, stats::sd) / sqrt(nrow(losses))
  minimum <- which.min(cv_mean)
  eligible <- which(cv_mean <= cv_mean[minimum] + cv_se[minimum])
  mean_selected <- colMeans(selected_counts)
  sparsest <- eligible[mean_selected[eligible] == min(mean_selected[eligible])]
  selected_id <- sparsest[which.max(fractions[sparsest])]
  tuning_runtime_sec <- proc.time()[["elapsed"]] - tuning_started

  final_started <- proc.time()[["elapsed"]]
  representation <- fit_fpca_representation(
    x_coef,
    y_coef,
    variance_threshold = variance_threshold
  )
  group <- integer(ncol(representation$x_scores))
  for (j in seq_along(representation$groups)) {
    group[representation$groups[[j]]] <- j
  }
  final_fit <- grpreg::grpreg(
    X = representation$x_scores,
    y = representation$y_scores,
    group = group,
    penalty = "grSCAD",
    nlambda = nlambda,
    lambda.min = lambda_min,
    gamma = gamma,
    max.iter = 10000L,
    eps = 1e-5
  )
  final_path_ids <- match_grpreg_fraction_path(final_fit, fractions)
  selected_path_id <- final_path_ids[[selected_id]]
  coefficient_with_intercept <- stats::coef(final_fit)[
    , , selected_path_id, drop = FALSE
  ][, , 1L]
  if (!is.matrix(coefficient_with_intercept)) {
    coefficient_with_intercept <- matrix(coefficient_with_intercept, nrow = 1L)
  }
  coefficient <- t(coefficient_with_intercept[, -1L, drop = FALSE])
  intercept <- coefficient_with_intercept[, 1L]
  selected_covariates <- which(vapply(representation$groups, function(rows) {
    sqrt(sum(coefficient[rows, , drop = FALSE]^2)) > 1e-8
  }, logical(1)))
  reconstructed <- reconstruct_fpca_coefficients(
    coefficient,
    representation,
    ncol(y_coef)
  )
  final_runtime_sec <- proc.time()[["elapsed"]] - final_started

  list(
    method_id = "fpca_group_scad",
    supports_selection = TRUE,
    supports_coefficients = TRUE,
    representation = representation,
    model = list(coef = coefficient, intercept = intercept),
    beta = reconstructed$beta,
    selected_covariates = selected_covariates,
    tuning = list(
      selected_lambda = final_fit$lambda[selected_path_id],
      selected_fraction = fractions[selected_id],
      selected_path_fraction = final_fit$lambda[selected_path_id] /
        final_fit$lambda[[1L]],
      returned_nlambda = length(final_fit$lambda),
      gamma = gamma,
      nlambda = nlambda,
      lambda_min_fraction = lambda_min,
      cv_mean = cv_mean,
      cv_se = cv_se,
      mean_selected_covariates = mean_selected,
      losses = losses,
      lambda = final_fit$lambda
    ),
    timing = list(
      tuning_runtime_sec = tuning_runtime_sec,
      final_runtime_sec = final_runtime_sec
    ),
    solver_converged = TRUE,
    final_iter = NA_integer_
  )
}

adaptive_enet_path <- function(
  x,
  y,
  alpha,
  lambda_fractions,
  adaptive_gamma = 1,
  initial_ridge = 1,
  weight_epsilon = 1e-4,
  max_iter = 2000L,
  tol = 1e-6
) {
  if (!is.finite(alpha) || alpha <= 0 || alpha > 1) {
    stop("alpha must be in (0, 1].")
  }
  prep <- standardize_matrix(x)
  y_center <- colMeans(y)
  yc <- sweep(y, 2L, y_center, "-")
  n <- nrow(x)
  q <- ncol(x)

  initial_coef <- solve(
    crossprod(prep$values) + initial_ridge * diag(q),
    crossprod(prep$values, yc)
  )
  initial_norms <- sqrt(rowSums(initial_coef^2))
  adaptive_weights <- (initial_norms + weight_epsilon)^(-adaptive_gamma)

  gradient_at_zero <- crossprod(prep$values, yc) / n
  gradient_norms <- sqrt(rowSums(gradient_at_zero^2))
  lambda_max <- max(gradient_norms / (alpha * adaptive_weights))
  if (!is.finite(lambda_max) || lambda_max <= 0) lambda_max <- 1

  eigenvalues <- eigen(
    crossprod(prep$values) / n,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  lipschitz <- max(eigenvalues, na.rm = TRUE)
  if (!is.finite(lipschitz) || lipschitz <= 0) lipschitz <- 1

  fractions <- sort(unique(lambda_fractions), decreasing = TRUE)
  coefficient_scaled <- matrix(0, nrow = q, ncol = ncol(y))
  fits <- vector("list", length(fractions))

  for (lambda_id in seq_along(fractions)) {
    lambda <- fractions[lambda_id] * lambda_max
    accelerated <- coefficient_scaled
    momentum <- 1
    converged <- FALSE

    for (iter in seq_len(max_iter)) {
      old <- coefficient_scaled
      gradient <- crossprod(
        prep$values,
        prep$values %*% accelerated - yc
      ) / n
      candidate <- accelerated - gradient / lipschitz
      ridge_denominator <- 1 + lambda * (1 - alpha) / lipschitz

      for (row_id in seq_len(q)) {
        row_norm <- sqrt(sum(candidate[row_id, ]^2))
        threshold <- lambda * alpha * adaptive_weights[row_id] / lipschitz
        candidate[row_id, ] <- if (!is.finite(row_norm) || row_norm <= threshold) {
          0
        } else {
          ((1 - threshold / row_norm) / ridge_denominator) * candidate[row_id, ]
        }
      }

      new_momentum <- (1 + sqrt(1 + 4 * momentum^2)) / 2
      accelerated <- candidate + ((momentum - 1) / new_momentum) * (candidate - old)
      coefficient_scaled <- candidate
      momentum <- new_momentum
      relative_change <- norm(coefficient_scaled - old, type = "F") /
        (norm(old, type = "F") + 1e-8)
      if (iter > 2L && relative_change < tol) {
        converged <- TRUE
        break
      }
    }

    coefficient <- sweep(coefficient_scaled, 1L, prep$scale, "/")
    intercept <- y_center - drop(prep$center %*% coefficient)
    fits[[lambda_id]] <- list(
      coef = coefficient,
      intercept = intercept,
      alpha = alpha,
      lambda = lambda,
      lambda_fraction = fractions[lambda_id],
      selected_rows = which(sqrt(rowSums(coefficient^2)) > 1e-8),
      converged = converged,
      final_iter = iter
    )
  }

  list(
    fractions = fractions,
    lambda_max = lambda_max,
    adaptive_weights = adaptive_weights,
    fits = fits
  )
}

selected_covariates_from_score_rows <- function(selected_rows, groups) {
  which(vapply(groups, function(rows) any(rows %in% selected_rows), logical(1)))
}

fit_fpca_aenet_baseline <- function(
  x_coef,
  y_coef,
  folds,
  alpha_grid,
  lambda_fractions,
  variance_threshold = 0.95,
  adaptive_gamma = 1,
  initial_ridge = 1,
  max_iter = 2000L,
  tol = 1e-6
) {
  tuning_started <- proc.time()[["elapsed"]]
  alphas <- sort(unique(alpha_grid))
  fractions <- sort(unique(lambda_fractions), decreasing = TRUE)
  candidates <- expand.grid(
    alpha = alphas,
    lambda_fraction = fractions,
    KEEP.OUT.ATTRS = FALSE
  )
  losses <- matrix(NA_real_, nrow = length(folds), ncol = nrow(candidates))
  selected_counts <- matrix(NA_real_, nrow = length(folds), ncol = nrow(candidates))
  n <- dim(x_coef)[2L]

  for (fold_id in seq_along(folds)) {
    valid <- folds[[fold_id]]
    train <- setdiff(seq_len(n), valid)
    representation <- fit_fpca_representation(
      x_coef[, train, , drop = FALSE],
      y_coef[train, , drop = FALSE],
      variance_threshold
    )
    valid_scores <- transform_fpca_predictors(
      x_coef[, valid, , drop = FALSE],
      representation
    )

    for (alpha_id in seq_along(alphas)) {
      path <- adaptive_enet_path(
        representation$x_scores,
        representation$y_scores,
        alphas[alpha_id],
        fractions,
        adaptive_gamma,
        initial_ridge,
        max_iter = max_iter,
        tol = tol
      )
      for (fraction_id in seq_along(fractions)) {
        candidate_id <- which(
          candidates$alpha == alphas[alpha_id] &
            candidates$lambda_fraction == fractions[fraction_id]
        )
        model <- path$fits[[fraction_id]]
        predicted_scores <- predict_multivariate_linear(model, valid_scores)
        prediction <- sweep(
          predicted_scores %*% t(representation$y_model$rotation),
          2L,
          representation$y_model$center,
          "+"
        )
        losses[fold_id, candidate_id] <- sqrt(mean(
          (prediction - y_coef[valid, , drop = FALSE])^2
        ))
        selected_counts[fold_id, candidate_id] <- length(
          selected_covariates_from_score_rows(
            model$selected_rows,
            representation$groups
          )
        )
      }
    }
  }

  cv_mean <- colMeans(losses)
  cv_se <- apply(losses, 2L, stats::sd) / sqrt(nrow(losses))
  minimum <- which.min(cv_mean)
  eligible <- which(cv_mean <= cv_mean[minimum] + cv_se[minimum])
  mean_selected <- colMeans(selected_counts)
  sparsest <- eligible[mean_selected[eligible] == min(mean_selected[eligible])]
  ordering <- order(
    -candidates$lambda_fraction[sparsest],
    -candidates$alpha[sparsest],
    cv_mean[sparsest]
  )
  selected_id <- sparsest[ordering[1L]]
  tuning_runtime_sec <- proc.time()[["elapsed"]] - tuning_started

  final_started <- proc.time()[["elapsed"]]
  representation <- fit_fpca_representation(x_coef, y_coef, variance_threshold)
  path <- adaptive_enet_path(
    representation$x_scores,
    representation$y_scores,
    candidates$alpha[selected_id],
    candidates$lambda_fraction[selected_id],
    adaptive_gamma,
    initial_ridge,
    max_iter = max_iter,
    tol = tol
  )
  model <- path$fits[[1L]]
  selected_covariates <- selected_covariates_from_score_rows(
    model$selected_rows,
    representation$groups
  )
  reconstructed <- reconstruct_fpca_coefficients(
    model$coef,
    representation,
    ncol(y_coef)
  )
  final_runtime_sec <- proc.time()[["elapsed"]] - final_started

  list(
    method_id = "fpca_aenet",
    supports_selection = TRUE,
    supports_coefficients = TRUE,
    representation = representation,
    model = model,
    beta = reconstructed$beta,
    selected_covariates = selected_covariates,
    tuning = list(
      selected_alpha = candidates$alpha[selected_id],
      selected_lambda = model$lambda,
      selected_fraction = candidates$lambda_fraction[selected_id],
      adaptive_gamma = adaptive_gamma,
      candidates = candidates,
      cv_mean = cv_mean,
      cv_se = cv_se,
      mean_selected_covariates = mean_selected,
      losses = losses
    ),
    timing = list(
      tuning_runtime_sec = tuning_runtime_sec,
      final_runtime_sec = final_runtime_sec
    ),
    solver_converged = model$converged,
    final_iter = model$final_iter
  )
}

squared_distance_matrix <- function(x, y = x) {
  distances <- outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y)
  pmax(distances, 0)
}

fit_kernel_ridge_once <- function(x, y, gamma_multiplier, lambda) {
  prep <- standardize_matrix(x)
  distance <- squared_distance_matrix(prep$values)
  positive <- distance[upper.tri(distance) & distance > 1e-12]
  median_distance <- if (length(positive) > 0L) stats::median(positive) else 1
  gamma <- gamma_multiplier / median_distance
  kernel <- exp(-gamma * distance)
  y_center <- colMeans(y)
  yc <- sweep(y, 2L, y_center, "-")
  alpha <- solve(kernel + nrow(x) * lambda * diag(nrow(x)), yc)
  list(
    train_x = prep$values,
    x_center = prep$center,
    x_scale = prep$scale,
    y_center = y_center,
    alpha = alpha,
    gamma = gamma,
    gamma_multiplier = gamma_multiplier,
    lambda = lambda
  )
}

predict_kernel_ridge <- function(model, x) {
  standardized <- sweep(sweep(x, 2L, model$x_center, "-"), 2L, model$x_scale, "/")
  kernel <- exp(-model$gamma * squared_distance_matrix(standardized, model$train_x))
  sweep(kernel %*% model$alpha, 2L, model$y_center, "+")
}

fit_kernel_ridge_baseline <- function(
  x_coef,
  y_coef,
  folds,
  gamma_multipliers,
  lambda_grid
) {
  tuning_started <- proc.time()[["elapsed"]]
  x <- flatten_predictor_coefficients(x_coef)
  candidates <- expand.grid(
    gamma_multiplier = gamma_multipliers,
    lambda = lambda_grid,
    KEEP.OUT.ATTRS = FALSE
  )
  losses <- matrix(NA_real_, nrow = length(folds), ncol = nrow(candidates))
  for (fold_id in seq_along(folds)) {
    valid <- folds[[fold_id]]
    train <- setdiff(seq_len(nrow(x)), valid)
    for (candidate_id in seq_len(nrow(candidates))) {
      model <- fit_kernel_ridge_once(
        x[train, , drop = FALSE],
        y_coef[train, , drop = FALSE],
        candidates$gamma_multiplier[candidate_id],
        candidates$lambda[candidate_id]
      )
      prediction <- predict_kernel_ridge(model, x[valid, , drop = FALSE])
      losses[fold_id, candidate_id] <- sqrt(mean((prediction - y_coef[valid, , drop = FALSE])^2))
    }
  }
  mean_loss <- colMeans(losses)
  best <- which.min(mean_loss)
  tuning_runtime_sec <- proc.time()[["elapsed"]] - tuning_started
  final_started <- proc.time()[["elapsed"]]
  model <- fit_kernel_ridge_once(
    x,
    y_coef,
    candidates$gamma_multiplier[best],
    candidates$lambda[best]
  )
  final_runtime_sec <- proc.time()[["elapsed"]] - final_started
  list(
    method_id = "kernel_ridge",
    supports_selection = FALSE,
    supports_coefficients = FALSE,
    model = model,
    beta = NULL,
    selected_covariates = NULL,
    tuning = list(
      selected_gamma_multiplier = candidates$gamma_multiplier[best],
      selected_lambda = candidates$lambda[best],
      candidates = candidates,
      cv_mean = mean_loss,
      cv_sd = apply(losses, 2L, stats::sd),
      losses = losses
    ),
    timing = list(
      tuning_runtime_sec = tuning_runtime_sec,
      final_runtime_sec = final_runtime_sec
    )
  )
}

fit_functional_baseline <- function(
  method_id,
  x_coef,
  y_coef,
  folds,
  config,
  structural_membership = NULL
) {
  switch(
    method_id,
    basis_ridge = fit_basis_ridge_baseline(
      x_coef, y_coef, folds, config$ridge_lambda_grid
    ),
    fpca_ridge = fit_fpca_ridge_baseline(
      x_coef, y_coef, folds, config$ridge_lambda_grid, config$fpca_variance_threshold
    ),
    fpca_group_lasso = fit_fpca_group_lasso_baseline(
      x_coef, y_coef, folds, config$group_lasso_lambda_fractions,
      config$fpca_variance_threshold, config$group_lasso_max_iter,
      config$group_lasso_tol
    ),
    fpca_group_scad = fit_fpca_group_scad_baseline(
      x_coef, y_coef, folds, config$fpca_variance_threshold,
      config$group_scad_nlambda, config$group_scad_lambda_min,
      config$group_scad_gamma
    ),
    basis_group_scad = fit_basis_group_scad_baseline(
      x_coef, y_coef, folds, config$group_scad_nlambda,
      config$group_scad_lambda_min, config$group_scad_gamma
    ),
    fpca_aenet = fit_fpca_aenet_baseline(
      x_coef, y_coef, folds, config$aenet_alpha_grid,
      config$aenet_lambda_fractions, config$fpca_variance_threshold,
      config$aenet_adaptive_gamma, config$aenet_initial_ridge,
      config$aenet_max_iter, config$aenet_tol
    ),
    structured_group_lasso = fit_structured_group_lasso_baseline(
      x_coef, y_coef, structural_membership, folds,
      config$structured_alpha_grid, config$structured_lambda_fractions,
      config$structured_max_iter, config$structured_tol
    ),
    kernel_ridge = fit_kernel_ridge_baseline(
      x_coef, y_coef, folds, config$kernel_gamma_multipliers,
      config$kernel_lambda_grid
    ),
    stop("Unknown baseline method: ", method_id)
  )
}

predict_functional_baseline <- function(fit, x_coef) {
  if (fit$method_id %in% c(
    "basis_ridge", "basis_group_scad", "structured_group_lasso"
  )) {
    return(predict_multivariate_linear(fit$model, flatten_predictor_coefficients(x_coef)))
  }
  if (fit$method_id %in% c(
    "fpca_ridge", "fpca_group_lasso", "fpca_group_scad", "fpca_aenet"
  )) {
    scores <- transform_fpca_predictors(x_coef, fit$representation)
    predicted_scores <- predict_multivariate_linear(fit$model, scores)
    return(sweep(
      predicted_scores %*% t(fit$representation$y_model$rotation),
      2L,
      fit$representation$y_model$center,
      "+"
    ))
  }
  if (fit$method_id == "kernel_ridge") {
    return(predict_kernel_ridge(fit$model, flatten_predictor_coefficients(x_coef)))
  }
  stop("No prediction method for: ", fit$method_id)
}

evaluate_functional_baseline <- function(fit, dgp, runtime_sec = NA_real_) {
  prediction <- predict_functional_baseline(fit, dgp$x_test)
  test_rmse <- sqrt(mean((prediction - dgp$y_test)^2))
  selected <- fit$selected_covariates

  selection_metrics <- c(
    covariate_tpr = NA_real_, covariate_fpr = NA_real_, covariate_fdr = NA_real_,
    selected_covariate_count = NA_real_, group_tpr = NA_real_, group_fpr = NA_real_,
    group_fdr = NA_real_, selected_group_count = NA_real_
  )
  if (isTRUE(fit$supports_selection)) {
    truth <- dgp$truth_covariates
    selected_flag <- truth$covariate_id %in% selected
    active <- truth$active
    selected_groups <- unique(dgp$structural_membership[
      covariate_id %in% selected,
      group_id
    ])
    group_truth <- dgp$truth_groups
    group_selected <- group_truth$group_id %in% selected_groups
    safe_rate <- function(flag, selected_flag, positive) {
      denominator <- sum(flag == positive)
      if (denominator == 0L) return(NA_real_)
      sum(selected_flag & flag == positive) / denominator
    }
    safe_fdr <- function(flag, selected_flag) {
      if (sum(selected_flag) == 0L) return(0)
      sum(selected_flag & !flag) / sum(selected_flag)
    }
    selection_metrics <- c(
      covariate_tpr = safe_rate(active, selected_flag, TRUE),
      covariate_fpr = safe_rate(active, selected_flag, FALSE),
      covariate_fdr = safe_fdr(active, selected_flag),
      selected_covariate_count = length(selected),
      group_tpr = safe_rate(group_truth$active, group_selected, TRUE),
      group_fpr = safe_rate(group_truth$active, group_selected, FALSE),
      group_fdr = safe_fdr(group_truth$active, group_selected),
      selected_group_count = length(selected_groups)
    )
  }

  coefficient_metrics <- c(
    coefficient_relative_error = NA_real_,
    active_coefficient_relative_error = NA_real_,
    inactive_coefficient_norm = NA_real_
  )
  if (isTRUE(fit$supports_coefficients)) {
    block_size <- dgp$dimensions$kx * dgp$dimensions$ky
    active_columns <- unlist(lapply(dgp$active_covariates, function(j) {
      ((j - 1L) * block_size + 1L):(j * block_size)
    }), use.names = FALSE)
    inactive <- setdiff(seq_len(dgp$dimensions$n_covariates), dgp$active_covariates)
    inactive_norm <- mean(vapply(inactive, function(j) {
      columns <- ((j - 1L) * block_size + 1L):(j * block_size)
      sqrt(sum(fit$beta[columns]^2))
    }, numeric(1)))
    coefficient_metrics <- c(
      coefficient_relative_error = sqrt(sum((fit$beta - dgp$beta_true)^2)) /
        (sqrt(sum(dgp$beta_true^2)) + 1e-8),
      active_coefficient_relative_error = sqrt(sum((fit$beta[active_columns] - dgp$beta_true[active_columns])^2)) /
        (sqrt(sum(dgp$beta_true[active_columns]^2)) + 1e-8),
      inactive_coefficient_norm = inactive_norm
    )
  }

  data.table::as.data.table(as.list(c(
    selection_metrics,
    coefficient_metrics,
    test_coeff_rmse = test_rmse,
    tuning_runtime_sec = if (!is.null(fit$timing$tuning_runtime_sec)) {
      fit$timing$tuning_runtime_sec
    } else {
      NA_real_
    },
    final_runtime_sec = if (!is.null(fit$timing$final_runtime_sec)) {
      fit$timing$final_runtime_sec
    } else {
      NA_real_
    },
    runtime_sec = runtime_sec,
    solver_converged = if (!is.null(fit$solver_converged)) fit$solver_converged else NA,
    final_iter = if (!is.null(fit$final_iter)) fit$final_iter else NA
  )))
}
