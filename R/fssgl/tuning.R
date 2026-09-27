# Curve-level cross-validation for FSSGL. The same folds and response-scale
# RMSE can be passed to every comparison method.

make_curve_cv_folds <- function(n, n_folds = 5L, seed = 1L) {
  n <- as.integer(n)
  n_folds <- max(2L, min(as.integer(n_folds), n))
  if (length(n) != 1L || !is.finite(n) || n < 2L) {
    stop("n must be one integer of at least two.")
  }
  set.seed(seed)
  assignment <- sample(rep(seq_len(n_folds), length.out = n))
  lapply(seq_len(n_folds), function(k) which(assignment == k))
}

validate_curve_cv_folds <- function(folds, n) {
  if (!is.list(folds) || length(folds) < 2L || any(lengths(folds) == 0L)) {
    stop("folds must be a list containing at least two nonempty validation sets.")
  }
  indices <- unlist(folds, use.names = FALSE)
  if (!identical(sort(as.integer(indices)), seq_len(as.integer(n)))) {
    stop("folds must partition every curve exactly once.")
  }
  invisible(TRUE)
}

scale_fssgl_penalties <- function(parameters, multiplier) {
  multiplier <- as.numeric(multiplier)
  if (length(multiplier) != 1L || !is.finite(multiplier) || multiplier <= 0) {
    stop("multiplier must be one finite positive number.")
  }
  out <- parameters
  out$covariate_penalty_scale <- parameters$covariate_penalty_scale * multiplier
  out$group_penalty_scale <- parameters$group_penalty_scale * multiplier
  out
}

select_fssgl_cv_candidate <- function(
  summary,
  selection_rule = c("one_se_sparsest", "minimum")
) {
  selection_rule <- match.arg(selection_rule)
  summary <- data.table::as.data.table(summary)
  eligible <- summary[
    all_folds_complete & strict_convergence_rate == 1 & is.finite(mean_validation_rmse)
  ]
  if (!nrow(eligible)) {
    stop("No multiplier completed every fold with strict convergence.")
  }
  minimum <- eligible[order(mean_validation_rmse, -multiplier)][1L]
  if (selection_rule == "minimum") return(minimum)

  threshold <- minimum$mean_validation_rmse + minimum$se_validation_rmse
  eligible[mean_validation_rmse <= threshold + 1e-12][
    order(-multiplier, mean_validation_rmse)
  ][1L]
}

is_stable_upper_null_plateau <- function(summary, tolerance = 1e-10) {
  summary <- data.table::as.data.table(summary)[order(multiplier)]
  if (nrow(summary) < 2L) return(FALSE)
  upper_pair <- tail(summary, 2L)
  all(upper_pair$all_folds_complete) &&
    all(upper_pair$strict_convergence_rate == 1) &&
    all(is.finite(upper_pair$mean_validation_rmse)) &&
    all(abs(upper_pair$mean_selected_covariates) < tolerance) &&
    abs(diff(upper_pair$mean_validation_rmse)) < tolerance
}

fit_fssgl_selected_path <- function(
  x_coef,
  y_coef,
  structural_membership,
  multiplier_grid,
  selected_multiplier,
  parameters,
  extension_factors = c(1L, 2L),
  verbose = FALSE
) {
  multiplier_grid <- sort(unique(as.numeric(multiplier_grid)), decreasing = TRUE)
  if (length(selected_multiplier) != 1L ||
      !selected_multiplier %in% multiplier_grid) {
    stop("selected_multiplier must belong to multiplier_grid.")
  }
  path_multipliers <- multiplier_grid[multiplier_grid >= selected_multiplier]
  warm_beta <- NULL
  warm_theta_covariate <- NULL
  warm_theta_group <- NULL
  warm_sigma2 <- NULL
  path <- vector("list", length(path_multipliers))
  fit <- NULL
  for (index in seq_along(path_multipliers)) {
    multiplier <- path_multipliers[index]
    candidate_parameters <- scale_fssgl_penalties(parameters, multiplier)
    if (!is.null(warm_theta_covariate)) {
      candidate_parameters$theta_covariate <- warm_theta_covariate
      candidate_parameters$theta_group <- warm_theta_group
      candidate_parameters$sigma2_init <- warm_sigma2
    }
    started <- proc.time()[["elapsed"]]
    fit <- fit_fssgl_until_converged(
      x_coef = x_coef,
      y_coef = y_coef,
      structural_membership = structural_membership,
      beta_init = warm_beta,
      parameters = candidate_parameters,
      extension_factors = extension_factors,
      verbose = verbose
    )
    runtime_sec <- proc.time()[["elapsed"]] - started
    path[[index]] <- data.table::data.table(
      multiplier = multiplier,
      strict_converged = isTRUE(fit$fit$convergence$strict_converged),
      selected_covariate_count = sum(fit$fit$covariate_posterior$selected),
      runtime_sec = runtime_sec
    )
    warm_beta <- fit$fit$beta
    warm_theta_covariate <- fit$fit$hyperparameters$theta_covariate_final
    warm_theta_group <- fit$fit$hyperparameters$theta_group_final
    warm_sigma2 <- fit$fit$hyperparameters$sigma2_final
  }
  list(fit = fit, path = data.table::rbindlist(path))
}

tune_fssgl_cv <- function(
  x_coef,
  y_coef,
  structural_membership,
  folds,
  multiplier_grid = fssgl_penalty_multiplier_grid(),
  parameters = fssgl_submission_parameters(dim(x_coef)[1L], dim(x_coef)[2L]),
  selection_rule = c("one_se_sparsest", "minimum"),
  extension_factors = c(1L, 2L),
  require_interior_minimum = TRUE,
  verbose = FALSE
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is required.")
  }
  selection_rule <- match.arg(selection_rule)
  n <- dim(x_coef)[2L]
  validate_curve_cv_folds(folds, n)
  multiplier_grid <- sort(unique(as.numeric(multiplier_grid)))
  if (!length(multiplier_grid) || any(!is.finite(multiplier_grid)) ||
      any(multiplier_grid <= 0)) {
    stop("multiplier_grid must contain finite positive values.")
  }

  rows <- list()
  row_index <- 1L
  for (fold_id in seq_along(folds)) {
    validation <- folds[[fold_id]]
    training <- setdiff(seq_len(n), validation)
    validation_design <- build_fof_design(
      x_coef[, validation, , drop = FALSE],
      n_response_basis = ncol(y_coef)
    )$design
    warm_beta <- NULL
    warm_theta_covariate <- NULL
    warm_theta_group <- NULL
    warm_sigma2 <- NULL

    # Strong-to-weak continuation reduces path instability for the nonconvex fit.
    for (multiplier in sort(multiplier_grid, decreasing = TRUE)) {
      candidate_parameters <- scale_fssgl_penalties(parameters, multiplier)
      if (!is.null(warm_theta_covariate)) {
        candidate_parameters$theta_covariate <- warm_theta_covariate
        candidate_parameters$theta_group <- warm_theta_group
        candidate_parameters$sigma2_init <- warm_sigma2
      }
      started <- proc.time()[["elapsed"]]
      fit <- tryCatch(
        fit_fssgl_until_converged(
          x_coef = x_coef[, training, , drop = FALSE],
          y_coef = y_coef[training, , drop = FALSE],
          structural_membership = structural_membership,
          beta_init = warm_beta,
          parameters = candidate_parameters,
          extension_factors = extension_factors,
          verbose = verbose
        ),
        error = function(e) e
      )
      runtime_sec <- proc.time()[["elapsed"]] - started
      if (inherits(fit, "error")) {
        status <- "error"
        error_message <- conditionMessage(fit)
        strict_converged <- FALSE
        total_outer_iter <- NA_integer_
        selected_count <- NA_integer_
        validation_rmse <- Inf
      } else {
        warm_beta <- fit$fit$beta
        warm_theta_covariate <- fit$fit$hyperparameters$theta_covariate_final
        warm_theta_group <- fit$fit$hyperparameters$theta_group_final
        warm_sigma2 <- fit$fit$hyperparameters$sigma2_final
        status <- "ok"
        error_message <- NA_character_
        strict_converged <- isTRUE(fit$fit$convergence$strict_converged)
        total_outer_iter <- fit$fit$convergence$total_outer_iter
        selected_count <- sum(fit$fit$covariate_posterior$selected)
        prediction <- matrix(
          drop(validation_design %*% fit$fit$beta),
          nrow = length(validation),
          ncol = ncol(y_coef)
        )
        validation_rmse <- if (strict_converged) {
          sqrt(mean((prediction - y_coef[validation, , drop = FALSE])^2))
        } else {
          Inf
        }
      }
      rows[[row_index]] <- data.table::data.table(
        fold_id = fold_id,
        multiplier = multiplier,
        validation_rmse = validation_rmse,
        selected_covariate_count = selected_count,
        strict_converged = strict_converged,
        total_outer_iter = total_outer_iter,
        runtime_sec = runtime_sec,
        status = status,
        error_message = error_message
      )
      row_index <- row_index + 1L
    }
  }

  path <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  summary <- path[, .(
    mean_validation_rmse = if (all(is.finite(validation_rmse))) mean(validation_rmse) else Inf,
    sd_validation_rmse = if (all(is.finite(validation_rmse))) stats::sd(validation_rmse) else Inf,
    se_validation_rmse = if (all(is.finite(validation_rmse))) {
      stats::sd(validation_rmse) / sqrt(.N)
    } else {
      Inf
    },
    mean_selected_covariates = if (all(is.finite(selected_covariate_count))) {
      mean(selected_covariate_count)
    } else {
      NA_real_
    },
    strict_convergence_rate = mean(strict_converged),
    all_folds_complete = all(status == "ok") && all(is.finite(validation_rmse)),
    tuning_runtime_sec = sum(runtime_sec)
  ), by = multiplier][order(multiplier)]
  selected <- select_fssgl_cv_candidate(summary, selection_rule)
  selected_at_lower_boundary <- selected$multiplier == min(multiplier_grid)
  selected_at_upper_boundary <- selected$multiplier == max(multiplier_grid)
  upper_null_plateau <- selected_at_upper_boundary &&
    is_stable_upper_null_plateau(summary)
  if (isTRUE(require_interior_minimum) &&
      (selected_at_lower_boundary ||
       (selected_at_upper_boundary && !upper_null_plateau))) {
    boundary <- if (selected_at_lower_boundary) "lower" else "upper"
    stop(
      "Selected multiplier lies on the ", boundary,
      " grid boundary; expand the grid before reporting results."
    )
  }

  list(
    protocol_version = FSSGL_TUNING_PROTOCOL_VERSION,
    selection_rule = selection_rule,
    folds = folds,
    path = path,
    summary = summary,
    selected = selected,
    selected_multiplier = selected$multiplier,
    selected_parameters = scale_fssgl_penalties(parameters, selected$multiplier),
    selected_at_lower_boundary = selected_at_lower_boundary,
    selected_at_upper_boundary = selected_at_upper_boundary,
    upper_null_plateau = upper_null_plateau,
    tuning_runtime_sec = sum(path$runtime_sec)
  )
}
