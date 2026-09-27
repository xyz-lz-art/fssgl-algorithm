# Generic structured FSSGL fitting routine.
# Required inputs:
#   x, y: stacked design matrix and response vector.
#   covariate_blocks: one row per functional covariate block, with
#     covariate_id, block_start, and block_end columns.
#   structural_membership: mapping from covariate_id to a structural group, with optional
#     membership_weight for overlapping structural groups.
# The function returns coefficient estimates plus covariate- and structure-level
# slab-score summaries. Application-specific names should be added by the
# calling workflow, not inside this generic method file.

stable_logistic_from_logs <- function(log_yes, log_no) {
  m <- pmax(log_yes, log_no)
  exp(log_yes - m) / (exp(log_yes - m) + exp(log_no - m))
}

validate_fssgl_inputs <- function(
  x,
  y,
  covariate_blocks,
  structural_membership,
  block_size
) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x))) {
    stop("x must be a finite numeric matrix.")
  }
  if (!is.numeric(y) || length(y) != nrow(x) || any(!is.finite(y))) {
    stop("y must be finite and have one value per row of x.")
  }
  if (length(block_size) != 1L || !is.finite(block_size) || block_size < 1L) {
    stop("block_size must be a positive scalar.")
  }

  required_blocks <- c("covariate_id", "block_start", "block_end")
  if (!all(required_blocks %in% names(covariate_blocks))) {
    stop("covariate_blocks must contain covariate_id, block_start, and block_end.")
  }
  if (nrow(covariate_blocks) < 1L || anyNA(covariate_blocks[, ..required_blocks])) {
    stop("covariate_blocks must be nonempty and cannot contain missing identifiers or indices.")
  }
  if (anyDuplicated(covariate_blocks$covariate_id)) {
    stop("covariate_id must be unique in covariate_blocks.")
  }

  block_columns <- unlist(Map(seq.int, covariate_blocks$block_start, covariate_blocks$block_end))
  if (!identical(sort(as.integer(block_columns)), seq_len(ncol(x)))) {
    stop("Covariate blocks must partition all columns of x exactly once.")
  }
  block_lengths <- covariate_blocks$block_end - covariate_blocks$block_start + 1L
  if (any(block_lengths != block_size)) {
    stop("Every covariate block must have length block_size.")
  }

  required_membership <- c("covariate_id", "group_id", "membership_weight")
  if (!all(required_membership %in% names(structural_membership))) {
    stop("structural_membership must contain covariate_id, group_id, and membership_weight.")
  }
  if (nrow(structural_membership) < 1L || anyNA(structural_membership[, ..required_membership])) {
    stop("structural_membership must be nonempty and cannot contain missing values.")
  }
  if (anyDuplicated(structural_membership[, .(covariate_id, group_id)])) {
    stop("Each covariate_id/group_id membership pair must be unique.")
  }
  if (any(!is.finite(structural_membership$membership_weight)) ||
      any(structural_membership$membership_weight <= 0)) {
    stop("membership_weight values must be finite and positive.")
  }
  if (!setequal(structural_membership$covariate_id, covariate_blocks$covariate_id)) {
    stop("Every covariate block must appear in structural_membership, with no unknown covariates.")
  }
  weight_sums <- structural_membership[, .(weight_sum = sum(membership_weight)), by = covariate_id]
  if (any(abs(weight_sums$weight_sum - 1) > 1e-8)) {
    stop("Membership weights must sum to one for each covariate.")
  }

  invisible(TRUE)
}

ssgl_slab_probability <- function(group_norm, group_dim, lambda_spike, lambda_slab, theta) {
  if (lambda_spike <= lambda_slab) {
    stop("lambda_spike must be larger than lambda_slab.")
  }
  theta <- min(max(theta, 1e-8), 1 - 1e-8)

  log_slab <- log(theta) + group_dim * log(lambda_slab) - lambda_slab * group_norm
  log_spike <- log1p(-theta) + group_dim * log(lambda_spike) - lambda_spike * group_norm
  stable_logistic_from_logs(log_slab, log_spike)
}

ssgl_effective_lambda <- function(slab_prob, lambda_spike, lambda_slab) {
  (1 - slab_prob) * lambda_spike + slab_prob * lambda_slab
}

estimate_block_lipschitz <- function(x, covariate_blocks, block_size = NULL) {
  n <- nrow(x)
  if (is.null(block_size)) {
    block_size <- covariate_blocks$block_end[1] - covariate_blocks$block_start[1] + 1L
  }
  vapply(
    seq_len(nrow(covariate_blocks)),
    function(i) {
      cols <- covariate_blocks$block_start[i]:covariate_blocks$block_end[i]
      xg <- x[, cols, drop = FALSE]
      values <- eigen(crossprod(xg), symmetric = TRUE, only.values = TRUE)$values
      max(values, na.rm = TRUE) / n
    },
    numeric(1)
  )
}

compute_covariate_beta_norm <- function(beta, covariate_blocks) {
  vapply(
    seq_len(nrow(covariate_blocks)),
    function(i) {
      cols <- covariate_blocks$block_start[i]:covariate_blocks$block_end[i]
      sqrt(sum(beta[cols]^2))
    },
    numeric(1)
  )
}

compute_group_norms <- function(covariate_norms, structural_membership) {
  tmp <- data.table::copy(structural_membership)
  tmp[, covariate_norm := covariate_norms[covariate_id]]
  tmp[, .(group_norm = sqrt(sum(membership_weight * covariate_norm^2))), by = group_id]
}

group_penalty_contribution <- function(group_table, structural_membership, group_lambda_eff) {
  tmp <- merge(
    structural_membership[, .(group_id, covariate_id, membership_weight)],
    group_table[, .(group_id, group_penalty_weight)],
    by = "group_id",
    all.x = TRUE
  )
  tmp[, group_lambda_eff := group_lambda_eff[group_id]]
  tmp[
    ,
    .(group_penalty_add = sum(membership_weight * group_penalty_weight * group_lambda_eff)),
    by = covariate_id
  ]
}

drop_existing_columns <- function(dt, cols) {
  data.table::copy(dt)[, setdiff(names(dt), cols), with = FALSE]
}

apply_supported_group_selection <- function(
  group_table,
  structural_membership,
  covariate_table,
  posterior_cutoff
) {
  selected_covariate <- covariate_table[selected == TRUE, .(covariate_id)]
  supported_group <- merge(
    structural_membership[, .(group_id, covariate_id)],
    selected_covariate,
    by = "covariate_id"
  )[, .(has_selected_covariate = TRUE), by = group_id]

  out <- merge(
    drop_existing_columns(group_table, c("score_selected", "has_selected_covariate", "selected")),
    supported_group,
    by = "group_id",
    all.x = TRUE
  )
  out[is.na(has_selected_covariate), has_selected_covariate := FALSE]
  out[, score_selected := posterior_slab_prob >= posterior_cutoff]
  out[, selected := has_selected_covariate]
  out
}

block_surrogate_objective <- function(
  x_block,
  residual_without_block,
  beta_block,
  penalty,
  n
) {
  residual_block <- residual_without_block - drop(x_block %*% beta_block)
  sum(residual_block^2) / (2 * n) + penalty * sqrt(sum(beta_block^2)) / n
}

fit_fssgl_solver <- function(
  x,
  y,
  covariate_blocks,
  structural_membership,
  block_size,
  lambda_covariate_spike = 40,
  lambda_covariate_slab = 2,
  theta_covariate = 0.05,
  lambda_group_spike = 25,
  lambda_group_slab = 1.5,
  theta_group = 0.20,
  max_iter = 100,
  tol = 1e-5,
  relaxed_tol = 0.012,
  objective_tail_tol = 0.08,
  step_multiplier = 10,
  use_backtracking = TRUE,
  backtracking_factor = 2,
  max_backtracking = 8,
  descent_tol = 1e-10,
  covariate_penalty_scale = 1,
  group_penalty_scale = 1,
  posterior_cutoff = 0.5,
  beta_init = NULL,
  verbose = TRUE
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is required.")
  }

  covariate_blocks <- data.table::as.data.table(covariate_blocks)
  structural_membership <- data.table::as.data.table(structural_membership)
  if (!"membership_weight" %in% names(structural_membership)) {
    structural_membership[, membership_weight := 1]
  }
  validate_fssgl_inputs(x, y, covariate_blocks, structural_membership, block_size)
  if (length(posterior_cutoff) != 1L || !is.finite(posterior_cutoff) ||
      posterior_cutoff <= 0 || posterior_cutoff >= 1) {
    stop("posterior_cutoff must lie strictly between zero and one.")
  }

  group_table <- structural_membership[
    ,
    .(
      member_count = data.table::uniqueN(covariate_id),
      effective_member_count = sum(membership_weight),
      effective_block_size = block_size * sum(membership_weight),
      group_penalty_weight = sqrt(block_size * sum(membership_weight))
    ),
    by = group_id
  ]
  data.table::setorder(group_table, group_id)

  n <- nrow(x)
  p <- ncol(x)
  beta <- if (is.null(beta_init)) numeric(p) else as.numeric(beta_init)
  if (length(beta) != p) {
    stop("beta_init has wrong length.")
  }

  y <- as.numeric(y)
  residual <- y - drop(x %*% beta)
  lipschitz <- step_multiplier * estimate_block_lipschitz(
    x,
    covariate_blocks,
    block_size = block_size
  )
  lipschitz[!is.finite(lipschitz) | lipschitz <= 0] <- 1
  residual_variance <- max(mean(residual^2), 1e-8)

  objective <- numeric(max_iter)
  history <- data.table::data.table(
    iter = integer(),
    objective = numeric(),
    objective_rel_change = numeric(),
    rss_value = numeric(),
    residual_variance = numeric(),
    beta_change = numeric(),
    selected_covariate_count = integer(),
    selected_group_count = integer(),
    backtracking_steps = integer()
  )

  covariate_penalty_weight <- sqrt(block_size)

  for (iter in seq_len(max_iter)) {
    beta_old <- beta
    iter_backtracking_steps <- 0L

    covariate_norms <- compute_covariate_beta_norm(beta, covariate_blocks)
    group_norms <- compute_group_norms(covariate_norms, structural_membership)
    group_table <- merge(
      drop_existing_columns(group_table, c("group_norm", "group_slab_prob", "group_lambda_eff")),
      group_norms,
      by = "group_id",
      all.x = TRUE
    )
    group_table[is.na(group_norm), group_norm := 0]
    group_table[, group_score_norm := group_norm / group_penalty_weight]
    group_table[
      ,
      group_slab_prob := ssgl_slab_probability(
        group_score_norm,
        1,
        lambda_group_spike,
        lambda_group_slab,
        theta_group
      )
    ]
    group_table[, group_lambda_eff := ssgl_effective_lambda(group_slab_prob, lambda_group_spike, lambda_group_slab)]

    group_lambda_eff <- group_table$group_lambda_eff
    names(group_lambda_eff) <- group_table$group_id

    group_add <- group_penalty_contribution(group_table, structural_membership, group_lambda_eff)
    covariate_table <- data.table::copy(covariate_blocks)
    covariate_table[, covariate_norm := covariate_norms]
    covariate_table[, covariate_score_norm := covariate_norm / covariate_penalty_weight]
    covariate_table[
      ,
      covariate_slab_prob := ssgl_slab_probability(
        covariate_score_norm,
        1,
        lambda_covariate_spike,
        lambda_covariate_slab,
        theta_covariate
      )
    ]
    covariate_table[, covariate_lambda_eff := ssgl_effective_lambda(covariate_slab_prob, lambda_covariate_spike, lambda_covariate_slab)]
    covariate_table <- merge(covariate_table, group_add, by = "covariate_id", all.x = TRUE)
    covariate_table[is.na(group_penalty_add), group_penalty_add := 0]
    data.table::setorder(covariate_table, covariate_id)
    covariate_table[
      ,
      total_penalty :=
        covariate_penalty_scale * covariate_penalty_weight * covariate_lambda_eff +
        group_penalty_scale * group_penalty_add
    ]

    for (g in seq_len(nrow(covariate_table))) {
      cols <- covariate_table$block_start[g]:covariate_table$block_end[g]
      old <- beta[cols]
      x_block <- x[, cols, drop = FALSE]
      residual_without_block <- residual + drop(x_block %*% old)
      penalty_g <- covariate_table$total_penalty[g]
      old_obj <- block_surrogate_objective(
        x_block,
        residual_without_block,
        old,
        penalty_g,
        n
      )

      local_lipschitz <- lipschitz[g]
      new <- old
      for (bt in seq_len(max_backtracking + 1L)) {
        current_residual <- residual_without_block - drop(x_block %*% old)
        smooth_grad <- -drop(crossprod(x_block, current_residual))
        z <- old - smooth_grad / (n * local_lipschitz)
        z_norm <- sqrt(sum(z^2))
        threshold <- penalty_g / (n * local_lipschitz)

        candidate <- if (z_norm <= threshold || !is.finite(z_norm)) {
          numeric(length(cols))
        } else {
          (1 - threshold / z_norm) * z
        }

        new_obj <- block_surrogate_objective(
          x_block,
          residual_without_block,
          candidate,
          penalty_g,
          n
        )
        if (!use_backtracking || new_obj <= old_obj + descent_tol || bt > max_backtracking) {
          new <- candidate
          iter_backtracking_steps <- iter_backtracking_steps + bt - 1L
          lipschitz[g] <- local_lipschitz
          break
        }
        local_lipschitz <- local_lipschitz * backtracking_factor
      }

      beta[cols] <- new
      residual <- residual_without_block - drop(x_block %*% new)
    }

    covariate_norms <- compute_covariate_beta_norm(beta, covariate_blocks)
    group_norms <- compute_group_norms(covariate_norms, structural_membership)

    covariate_table[, beta_norm := covariate_norms]
    covariate_table[, posterior_score_norm := beta_norm / covariate_penalty_weight]
    covariate_table[
      ,
      posterior_slab_prob := ssgl_slab_probability(
        posterior_score_norm,
        1,
        lambda_covariate_spike,
        lambda_covariate_slab,
        theta_covariate
      )
    ]
    covariate_table[, selected := posterior_slab_prob >= posterior_cutoff]

    group_table <- merge(
      drop_existing_columns(group_table, c("group_norm", "group_slab_prob", "group_lambda_eff")),
      group_norms,
      by = "group_id",
      all.x = TRUE
    )
    group_table[is.na(group_norm), group_norm := 0]
    group_table[, posterior_score_norm := group_norm / group_penalty_weight]
    group_table[
      ,
      posterior_slab_prob := ssgl_slab_probability(
        posterior_score_norm,
        1,
        lambda_group_spike,
        lambda_group_slab,
        theta_group
      )
    ]
    group_table <- apply_supported_group_selection(
      group_table,
      structural_membership,
      covariate_table,
      posterior_cutoff
    )

    raw_rss <- sum(residual^2)
    residual_variance <- max(raw_rss / n, 1e-8)
    rss <- raw_rss / (2 * n)
    penalty_value <- sum(covariate_table$total_penalty * covariate_table$beta_norm) / n
    objective[iter] <- rss + penalty_value
    objective_rel_change <- if (iter == 1) {
      NA_real_
    } else {
      abs(objective[iter] - objective[iter - 1]) / (1 + abs(objective[iter - 1]))
    }
    beta_change <- sqrt(sum((beta - beta_old)^2)) / (sqrt(sum(beta_old^2)) + 1e-8)

    history <- rbind(
      history,
      data.table::data.table(
        iter = iter,
        objective = objective[iter],
        objective_rel_change = objective_rel_change,
        rss_value = rss,
        residual_variance = residual_variance,
        beta_change = beta_change,
        selected_covariate_count = sum(covariate_table$selected),
        selected_group_count = sum(group_table$selected),
        backtracking_steps = iter_backtracking_steps
      )
    )

    if (verbose && (iter == 1 || iter %% 10 == 0)) {
      message(
        "iter=", iter,
        " objective=", signif(objective[iter], 5),
        " beta_change=", signif(beta_change, 4),
        " selected covariates=", sum(covariate_table$selected),
        " selected groups=", sum(group_table$selected)
      )
    }

    if (iter > 2 && beta_change < tol) {
      break
    }
  }

  final_iter <- nrow(history)
  tail_n <- min(10L, final_iter)
  objective_tail <- tail(history$objective, tail_n)
  objective_tail_rel_change <- if (length(objective_tail) > 1) {
    abs(tail(objective_tail, 1) - objective_tail[1]) / (abs(objective_tail[1]) + 1e-8)
  } else {
    NA_real_
  }
  final_beta_change <- history$beta_change[final_iter]
  strict_converged <- final_beta_change < tol
  relaxed_converged <- isTRUE(final_beta_change < relaxed_tol) &&
    isTRUE(objective_tail_rel_change < objective_tail_tol)
  history[, strict_converged := FALSE]
  history[, relaxed_converged := FALSE]
  history[final_iter, strict_converged := strict_converged]
  history[final_iter, relaxed_converged := relaxed_converged]

  covariate_table <- covariate_table[order(-posterior_slab_prob, -beta_norm)]
  group_table <- group_table[order(-selected, -posterior_slab_prob, -group_norm)]

  list(
    beta = beta,
    fitted = y - residual,
    residual = residual,
    coefficient_estimates = beta,
    covariate_posterior = covariate_table,
    structure_posterior = group_table,
    structural_membership = structural_membership,
    history = history,
    convergence = list(
      strict_converged = strict_converged,
      relaxed_converged = relaxed_converged,
      final_iter = final_iter,
      final_beta_change = final_beta_change,
      objective_tail_rel_change = objective_tail_rel_change,
      final_residual_variance = residual_variance,
      strict_tol = tol,
      relaxed_tol = relaxed_tol,
      objective_tail_tol = objective_tail_tol
    ),
    hyperparameters = list(
      lambda_covariate_spike = lambda_covariate_spike,
      lambda_covariate_slab = lambda_covariate_slab,
      theta_covariate = theta_covariate,
      lambda_group_spike = lambda_group_spike,
      lambda_group_slab = lambda_group_slab,
      theta_group = theta_group,
      block_size = block_size,
      max_iter = max_iter,
      tol = tol,
      relaxed_tol = relaxed_tol,
      objective_tail_tol = objective_tail_tol,
      step_multiplier = step_multiplier,
      use_backtracking = use_backtracking,
      backtracking_factor = backtracking_factor,
      max_backtracking = max_backtracking,
      descent_tol = descent_tol,
      covariate_penalty_scale = covariate_penalty_scale,
      group_penalty_scale = group_penalty_scale,
      posterior_cutoff = posterior_cutoff
    )
  )
}
