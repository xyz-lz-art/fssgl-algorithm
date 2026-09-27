# Generic warm-start utility used by the experiment scripts.

ridge_dual_fit <- function(x, y, lambda = 1) {
  n <- nrow(x)
  gram <- tcrossprod(x)
  alpha <- solve(gram + lambda * diag(n), y)
  drop(crossprod(x, alpha))
}
