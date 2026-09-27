# Finite-sample identifiability diagnostics for coefficient-score designs.

spectral_design_diagnostics <- function(x, relative_tolerance = 1e-8) {
  x <- as.matrix(x)
  xc <- sweep(x, 2L, colMeans(x), "-")
  singular_values <- svd(xc, nu = 0L, nv = 0L)$d
  largest <- if (length(singular_values)) max(singular_values) else 0
  threshold <- largest * relative_tolerance
  positive <- singular_values[singular_values > threshold]
  numerical_rank <- length(positive)
  full_column_rank <- numerical_rank == ncol(xc)
  gram_condition <- if (!full_column_rank || !length(positive)) {
    Inf
  } else {
    (max(positive) / min(positive))^2
  }
  nonzero_gram_condition <- if (length(positive) <= 1L) {
    NA_real_
  } else {
    (max(positive) / min(positive))^2
  }
  eigen_mass <- singular_values^2
  if (sum(eigen_mass) <= 0) {
    entropy_effective_rank <- 0
    rank_99pct <- 0L
  } else {
    probability <- eigen_mass / sum(eigen_mass)
    entropy_effective_rank <- exp(-sum(probability[probability > 0] * log(probability[probability > 0])))
    rank_99pct <- which(cumsum(probability) >= 0.99)[1L]
  }
  list(
    n_rows = nrow(xc),
    n_columns = ncol(xc),
    numerical_rank = numerical_rank,
    rank_ratio = numerical_rank / max(1L, ncol(xc)),
    entropy_effective_rank = entropy_effective_rank,
    rank_99pct = rank_99pct,
    gram_condition = gram_condition,
    nonzero_gram_condition = nonzero_gram_condition,
    design_rank_deficient = !full_column_rank,
    severely_ill_conditioned = is.finite(nonzero_gram_condition) &&
      nonzero_gram_condition >= 1e6
  )
}

functional_predictor_identifiability <- function(x_coef, relative_tolerance = 1e-8) {
  if (length(dim(x_coef)) != 3L) {
    stop("x_coef must be a predictor x sample x coefficient array.")
  }
  p <- dim(x_coef)[1L]
  n <- dim(x_coef)[2L]
  kx <- dim(x_coef)[3L]
  flattened <- do.call(cbind, lapply(seq_len(p), function(j) x_coef[j, , ]))
  joint <- spectral_design_diagnostics(flattened, relative_tolerance)
  block <- lapply(seq_len(p), function(j) {
    spectral_design_diagnostics(x_coef[j, , ], relative_tolerance)
  })
  block_rank <- vapply(block, `[[`, numeric(1), "numerical_rank")
  block_effective_rank <- vapply(block, `[[`, numeric(1), "entropy_effective_rank")
  block_condition <- vapply(block, `[[`, numeric(1), "gram_condition")

  profiles <- vapply(seq_len(p), function(j) as.vector(x_coef[j, , ]), numeric(n * kx))
  profile_correlation <- suppressWarnings(stats::cor(profiles))
  pairwise <- abs(profile_correlation[upper.tri(profile_correlation)])
  pairwise <- pairwise[is.finite(pairwise)]

  c(
    joint,
    list(
      predictor_count = p,
      basis_dimension = kx,
      minimum_block_rank = min(block_rank),
      median_block_effective_rank = stats::median(block_effective_rank),
      minimum_block_effective_rank = min(block_effective_rank),
      maximum_block_gram_condition = max(block_condition),
      rank_deficient_block_fraction = mean(block_rank < kx),
      maximum_absolute_profile_correlation = if (length(pairwise)) max(pairwise) else NA_real_,
      q95_absolute_profile_correlation = if (length(pairwise)) unname(stats::quantile(pairwise, 0.95)) else NA_real_,
      pair_fraction_above_0_9 = if (length(pairwise)) mean(pairwise >= 0.9) else NA_real_
    )
  )
}

top_functional_profile_correlations <- function(x_coef, labels = NULL, top_n = 10L) {
  p <- dim(x_coef)[1L]
  n <- dim(x_coef)[2L]
  kx <- dim(x_coef)[3L]
  if (is.null(labels)) labels <- as.character(seq_len(p))
  if (length(labels) != p) stop("labels must have one entry per predictor.")
  profiles <- vapply(seq_len(p), function(j) as.vector(x_coef[j, , ]), numeric(n * kx))
  correlation <- suppressWarnings(stats::cor(profiles))
  index <- which(upper.tri(correlation), arr.ind = TRUE)
  out <- data.frame(
    predictor_1 = labels[index[, 1L]],
    predictor_2 = labels[index[, 2L]],
    correlation = correlation[index],
    absolute_correlation = abs(correlation[index]),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$absolute_correlation), , drop = FALSE]
  out <- out[order(out$absolute_correlation, decreasing = TRUE), , drop = FALSE]
  utils::head(out, as.integer(top_n))
}
