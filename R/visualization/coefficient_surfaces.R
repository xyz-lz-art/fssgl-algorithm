# Reusable coefficient-surface reconstruction and plotting helpers.

surface_axis <- function(basis) {
  axis <- attr(basis, "x_scaled")
  if (is.null(axis)) axis <- seq(0, 1, length.out = nrow(basis))
  as.numeric(axis)
}

coefficient_matrix_from_block <- function(beta, columns, kx, ky) {
  columns <- as.integer(columns)
  if (length(columns) != kx * ky || any(columns < 1L) || any(columns > length(beta))) {
    stop("Each coefficient block must contain exactly kx * ky valid columns.")
  }
  matrix(beta[columns], nrow = kx, ncol = ky)
}

reconstruct_coefficient_surfaces <- function(
  beta,
  predictor_blocks,
  basis_x,
  basis_y,
  predictor_ids = seq_along(predictor_blocks)
) {
  predictor_ids <- as.integer(predictor_ids)
  if (!length(predictor_ids) || any(!predictor_ids %in% seq_along(predictor_blocks))) {
    stop("predictor_ids must index predictor_blocks.")
  }
  lapply(predictor_ids, function(j) {
    coefficient_matrix <- coefficient_matrix_from_block(
      beta,
      predictor_blocks[[j]],
      ncol(basis_x),
      ncol(basis_y)
    )
    evaluate_coefficient_surface(coefficient_matrix, basis_x, basis_y)
  })
}

symmetric_surface_limit <- function(surfaces, fallback = 1) {
  values <- unlist(surfaces, use.names = FALSE)
  limit <- max(abs(values[is.finite(values)]), 0)
  if (!is.finite(limit) || limit <= 0) limit <- fallback
  c(-limit, limit)
}

open_surface_device <- function(output_file, width, height, resolution) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  extension <- tolower(tools::file_ext(output_file))
  if (extension == "png") {
    grDevices::png(output_file, width = width, height = height, res = resolution)
  } else if (extension == "pdf") {
    grDevices::pdf(
      output_file,
      width = width / resolution,
      height = height / resolution,
      family = "Helvetica"
    )
  } else {
    stop("output_file must have extension .png or .pdf.")
  }
}

draw_coefficient_surface <- function(
  surface,
  x,
  y,
  title,
  zlim,
  palette,
  xlab = "Predictor domain",
  ylab = "Response domain"
) {
  graphics::image(
    x, y, surface,
    col = palette,
    zlim = zlim,
    xlab = xlab,
    ylab = ylab,
    main = title,
    cex.main = 0.90,
    cex.lab = 0.82,
    cex.axis = 0.78,
    useRaster = TRUE
  )
  if (diff(range(surface, finite = TRUE)) > 1e-12) {
    graphics::contour(
      x, y, surface,
      add = TRUE,
      drawlabels = FALSE,
      col = "#22222270",
      lwd = 0.55
    )
  }
}

plot_surface_recovery <- function(
  beta_true,
  beta_hat,
  predictor_blocks,
  basis_x,
  basis_y,
  predictor_ids,
  predictor_labels = paste0("Predictor ", predictor_ids),
  output_file,
  note = NULL,
  width = 1800,
  row_height = 600,
  resolution = 180
) {
  if (length(beta_true) != length(beta_hat)) {
    stop("beta_true and beta_hat must have the same length.")
  }
  if (length(predictor_labels) != length(predictor_ids)) {
    stop("predictor_labels must match predictor_ids.")
  }
  truth <- reconstruct_coefficient_surfaces(
    beta_true, predictor_blocks, basis_x, basis_y, predictor_ids
  )
  estimate <- reconstruct_coefficient_surfaces(
    beta_hat, predictor_blocks, basis_x, basis_y, predictor_ids
  )
  difference <- Map(`-`, estimate, truth)
  effect_limit <- symmetric_surface_limit(c(truth, estimate))
  difference_limit <- symmetric_surface_limit(difference)
  x <- surface_axis(basis_x)
  y <- surface_axis(basis_y)
  palette <- grDevices::colorRampPalette(c("#2454A6", "#F7F7F7", "#B83227"))(101)

  open_surface_device(
    output_file,
    width = width,
    height = max(row_height * length(predictor_ids), row_height),
    resolution = resolution
  )
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(
    mfrow = c(length(predictor_ids), 3L),
    mar = c(4.2, 4.0, 3.2, 1.2),
    oma = c(if (is.null(note)) 0 else 2.8, 0, 0, 0)
  )
  for (i in seq_along(predictor_ids)) {
    draw_coefficient_surface(
      truth[[i]], x, y,
      paste0(predictor_labels[i], ": truth"),
      effect_limit, palette
    )
    draw_coefficient_surface(
      estimate[[i]], x, y,
      paste0(predictor_labels[i], ": estimate"),
      effect_limit, palette
    )
    draw_coefficient_surface(
      difference[[i]], x, y,
      paste0(predictor_labels[i], ": estimate - truth"),
      difference_limit, palette
    )
  }
  if (!is.null(note)) {
    graphics::mtext(note, side = 1, outer = TRUE, line = 0.6, cex = 0.78)
  }
  invisible(list(
    truth = truth,
    estimate = estimate,
    difference = difference,
    effect_limit = effect_limit,
    difference_limit = difference_limit,
    output_file = normalizePath(output_file, winslash = "/", mustWork = FALSE)
  ))
}

plot_estimated_surface_grid <- function(
  surfaces,
  labels,
  output_file,
  x = seq(0, 1, length.out = nrow(surfaces[[1L]])),
  y = seq(0, 1, length.out = ncol(surfaces[[1L]])),
  columns = 2L,
  note = NULL,
  width = 1600,
  panel_height = 600,
  resolution = 180
) {
  if (!length(surfaces) || length(labels) != length(surfaces)) {
    stop("surfaces and labels must be nonempty and have the same length.")
  }
  reference_dimension <- dim(as.matrix(surfaces[[1L]]))
  if (any(vapply(surfaces, function(z) {
    !identical(dim(as.matrix(z)), reference_dimension) || any(!is.finite(z))
  }, logical(1L)))) {
    stop("All surfaces must be finite matrices with common dimensions.")
  }
  columns <- max(1L, min(as.integer(columns), length(surfaces)))
  rows <- ceiling(length(surfaces) / columns)
  zlim <- symmetric_surface_limit(surfaces)
  palette <- grDevices::colorRampPalette(c("#2454A6", "#F7F7F7", "#B83227"))(101)
  open_surface_device(
    output_file,
    width = width,
    height = rows * panel_height,
    resolution = resolution
  )
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(
    mfrow = c(rows, columns),
    mar = c(4.2, 4.0, 3.3, 1.2),
    oma = c(if (is.null(note)) 0 else 2.8, 0, 0, 0)
  )
  for (i in seq_along(surfaces)) {
    draw_coefficient_surface(
      surfaces[[i]], x, y, labels[i], zlim, palette
    )
  }
  if (!is.null(note)) {
    graphics::mtext(note, side = 1, outer = TRUE, line = 0.6, cex = 0.78)
  }
  invisible(list(
    zlim = zlim,
    output_file = normalizePath(output_file, winslash = "/", mustWork = FALSE)
  ))
}
