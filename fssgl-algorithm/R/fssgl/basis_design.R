# Generic basis-expansion helpers for function-on-function regression.
# Use these functions before fitting FSSGL: build marginal bases, project
# observed curves to basis coefficients, then construct the stacked design.

trapezoid_weights <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 2L || any(!is.finite(x)) || is.unsorted(x, strictly = TRUE)) {
    stop("x must contain at least two finite, strictly increasing grid points.")
  }
  gaps <- diff(x)
  weights <- c(gaps[1L] / 2, (gaps[-1L] + gaps[-length(gaps)]) / 2, gaps[length(gaps)] / 2)
  weights / sum(weights)
}

make_bspline_basis <- function(x, df = 6, degree = 3) {
  if (df < 2L) stop("df must be at least two.")
  x <- as.numeric(x)
  x_range <- range(x)
  if (!all(is.finite(x_range)) || diff(x_range) <= 0) {
    stop("x must span a finite interval with positive length.")
  }
  x_scaled <- (x - x_range[1L]) / diff(x_range)
  degree <- min(as.integer(degree), as.integer(df) - 1L)
  raw_basis <- splines::bs(
    x_scaled,
    df = df,
    degree = degree,
    intercept = TRUE,
    Boundary.knots = c(0, 1)
  )
  weights <- trapezoid_weights(x_scaled)
  raw_gram <- crossprod(raw_basis, sweep(raw_basis, 1L, weights, `*`))
  chol_gram <- chol(raw_gram)
  basis <- raw_basis %*% backsolve(chol_gram, diag(ncol(chol_gram)))
  gram <- crossprod(basis, sweep(basis, 1L, weights, `*`))

  attr(basis, "x_scaled") <- x_scaled
  attr(basis, "quadrature_weights") <- weights
  attr(basis, "raw_basis") <- raw_basis
  attr(basis, "raw_gram") <- raw_gram
  attr(basis, "orthonormalization") <- chol_gram
  attr(basis, "gram") <- gram
  attr(basis, "geometry") <- "orthonormal_trapezoid"
  basis
}

# Weighted L2 projection onto a supplied marginal basis. Rows are independent
# curves and columns are grid points. Orthonormal bases have an identity Gram
# matrix up to quadrature error, but the solve is retained as a numerical guard.
project_curves_l2 <- function(curves, basis) {
  curves <- as.matrix(curves)
  if (ncol(curves) != nrow(basis)) {
    stop("curves and basis must use the same observation grid.")
  }
  weights <- attr(basis, "quadrature_weights")
  if (is.null(weights)) {
    weights <- rep(1 / nrow(basis), nrow(basis))
  }
  gram <- crossprod(basis, sweep(basis, 1L, weights, `*`))
  weighted_curves <- sweep(curves, 2L, weights, `*`)
  weighted_curves %*% basis %*% solve(gram)
}

project_curves_lm <- project_curves_l2

# Project a 3D predictor array with dimensions covariate x sample x grid.
project_predictor_array <- function(x_array, basis) {
  n_predictor <- dim(x_array)[1]
  n_sample <- dim(x_array)[2]
  n_basis <- ncol(basis)

  out <- array(
    NA_real_,
    dim = c(n_predictor, n_sample, n_basis),
    dimnames = list(dimnames(x_array)[[1]], dimnames(x_array)[[2]], paste0("kx", seq_len(n_basis)))
  )

  for (j in seq_len(n_predictor)) {
    out[j, , ] <- project_curves_l2(x_array[j, , , drop = FALSE][1, , ], basis)
  }

  out
}

evaluate_coefficient_surface <- function(coefficient_matrix, basis_x, basis_y) {
  coefficient_matrix <- as.matrix(coefficient_matrix)
  if (!identical(dim(coefficient_matrix), c(ncol(basis_x), ncol(basis_y)))) {
    stop("coefficient_matrix dimensions must match the two marginal bases.")
  }
  basis_x %*% coefficient_matrix %*% t(basis_y)
}

surface_l2_norm <- function(surface, basis_x, basis_y) {
  surface <- as.matrix(surface)
  if (!identical(dim(surface), c(nrow(basis_x), nrow(basis_y)))) {
    stop("surface dimensions must match the two marginal grids.")
  }
  weights_x <- attr(basis_x, "quadrature_weights")
  weights_y <- attr(basis_y, "quadrature_weights")
  sqrt(sum(surface^2 * outer(weights_x, weights_y)))
}

# Build the stacked coefficient-space design. If X_j is n by K_X and A_j is
# K_X by K_Y, this uses vec(X_j A_j) = (I_{K_Y} %x% X_j) vec(A_j).
# Each covariate block contains vec(A_j) in column-major order.
build_fof_design <- function(x_coef, n_response_basis) {
  n_predictor <- dim(x_coef)[1]
  n_sample <- dim(x_coef)[2]
  n_predictor_basis <- dim(x_coef)[3]

  design <- matrix(
    0,
    nrow = n_sample * n_response_basis,
    ncol = n_predictor * n_predictor_basis * n_response_basis
  )

  predictor_blocks <- vector("list", n_predictor)
  col_start <- 1

  for (j in seq_len(n_predictor)) {
    xj <- x_coef[j, , , drop = FALSE][1, , ]
    block <- kronecker(diag(n_response_basis), xj)
    cols <- col_start:(col_start + ncol(block) - 1)
    design[, cols] <- block
    predictor_blocks[[j]] <- cols
    col_start <- max(cols) + 1
  }

  list(design = design, predictor_blocks = predictor_blocks)
}

build_functional_groups <- function(predictor_blocks) {
  data.table::data.table(
    covariate_id = seq_along(predictor_blocks),
    block_start = vapply(predictor_blocks, min, integer(1)),
    block_end = vapply(predictor_blocks, max, integer(1))
  )
}
