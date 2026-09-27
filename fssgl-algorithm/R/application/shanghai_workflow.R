# Shanghai MetroFlow adapter for the generic FSSGL routines.
# This file is the only place where station/line semantics are attached to
# the generic covariate-block and structural-group structure used by R/fssgl/.

SHANGHAI_PREPROCESSING_VERSION <- "training_pointwise_function_centering_v2"

clean_shanghai_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", tolower(x))
  gsub("_+$", "", gsub("^_+", "", x))
}

load_shanghai_data <- function(processed_dir, metadata_dir) {
  clean_file <- file.path(processed_dir, "shanghai_metroflow_clean.RData")
  membership_file <- file.path(metadata_dir, "station_line_membership_long.csv")

  if (!file.exists(clean_file)) {
    stop("Missing shanghai_metroflow_clean.RData.")
  }
  if (!file.exists(membership_file)) {
    stop("Missing station_line_membership_long.csv.")
  }

  env <- new.env(parent = emptyenv())
  load(clean_file, envir = env)
  list(
    arr_in = env$arr_in,
    arr_out = env$arr_out,
    station_metadata = env$station_metadata,
    date_index = env$date_index,
    time_index = env$time_index,
    membership = data.table::fread(membership_file)
  )
}

validate_shanghai_day_split <- function(n_day, train_days, test_days) {
  normalize_days <- function(days, label) {
    days <- as.integer(days)
    if (!length(days) || anyNA(days) || any(days < 1L | days > n_day)) {
      stop(label, " must contain valid day positions in 1:", n_day, ".")
    }
    if (anyDuplicated(days)) stop(label, " contains duplicated day positions.")
    sort(days)
  }

  train_days <- normalize_days(train_days, "train_days")
  test_days <- normalize_days(test_days, "test_days")
  if (length(intersect(train_days, test_days))) {
    stop("train_days and test_days must not overlap.")
  }
  list(train_days = train_days, test_days = test_days)
}

prepare_shanghai_target <- function(
  target_name,
  processed_dir,
  metadata_dir,
  train_days = NULL,
  test_days = NULL,
  shanghai_data = NULL
) {
  if (is.null(shanghai_data)) {
    shanghai_data <- load_shanghai_data(processed_dir, metadata_dir)
  }

  target_pos <- shanghai_data$station_metadata[name == target_name, station_pos]
  if (length(target_pos) != 1) {
    stop("Target station not found or ambiguous: ", target_name)
  }

  predictor_pos <- setdiff(seq_len(dim(shanghai_data$arr_in)[1]), target_pos)
  predictor_station <- shanghai_data$station_metadata[
    predictor_pos,
    .(station_pos, station, name, lines, n_lines)
  ]

  n_day <- dim(shanghai_data$arr_in)[2]
  if (is.null(train_days) && is.null(test_days)) {
    train_days <- seq_len(floor(0.8 * n_day))
    test_days <- setdiff(seq_len(n_day), train_days)
  } else if (is.null(train_days) || is.null(test_days)) {
    stop("train_days and test_days must either both be supplied or both be omitted.")
  }
  split <- validate_shanghai_day_split(n_day, train_days, test_days)
  train_days <- split$train_days
  test_days <- split$test_days

  y_train <- shanghai_data$arr_out[target_pos, train_days, , drop = FALSE][1, , ]
  y_test <- shanghai_data$arr_out[target_pos, test_days, , drop = FALSE][1, , ]
  x_train <- shanghai_data$arr_in[predictor_pos, train_days, , drop = FALSE]
  x_test <- shanghai_data$arr_in[predictor_pos, test_days, , drop = FALSE]

  # The coefficient model has no functional intercept. Remove the training
  # mean curve at every grid point before applying one residual scale per
  # station. Using a single scalar mean would leave the common diurnal shape in
  # every curve and allow it to masquerade as cross-station association.
  station_mean_curve <- apply(x_train, c(1L, 3L), mean)
  station_sd <- numeric(nrow(station_mean_curve))
  for (j in seq_len(nrow(station_mean_curve))) {
    x_train[j, , ] <- sweep(x_train[j, , ], 2L, station_mean_curve[j, ], "-")
    x_test[j, , ] <- sweep(x_test[j, , ], 2L, station_mean_curve[j, ], "-")
    station_sd[j] <- stats::sd(as.vector(x_train[j, , ]))
    if (!is.finite(station_sd[j]) || station_sd[j] == 0) station_sd[j] <- 1
    x_train[j, , ] <- x_train[j, , ] / station_sd[j]
    x_test[j, , ] <- x_test[j, , ] / station_sd[j]
  }

  target_mean_curve <- colMeans(y_train)
  y_train <- sweep(y_train, 2L, target_mean_curve, "-")
  y_test <- sweep(y_test, 2L, target_mean_curve, "-")
  target_sd <- stats::sd(as.vector(y_train))
  if (target_sd == 0 || is.na(target_sd)) target_sd <- 1
  y_train <- y_train / target_sd
  y_test <- y_test / target_sd

  line_membership <- merge(
    predictor_station[, .(predictor_id = seq_len(.N), station, station_pos, name, n_lines)],
    shanghai_data$membership[, .(stationID, line)],
    by.x = "station",
    by.y = "stationID",
    all.x = TRUE
  )
  line_membership[, line_share := 1 / pmax(n_lines, 1)]
  line_membership[, covariate_id := predictor_id]
  line_membership[, group_id := line]
  line_membership[, membership_weight := line_share]

  list(
    target = shanghai_data$station_metadata[target_pos],
    predictor_station = predictor_station,
    line_membership = line_membership,
    train_days = train_days,
    test_days = test_days,
    date_index = shanghai_data$date_index,
    time_index = shanghai_data$time_index,
    X_train = x_train,
    X_test = x_test,
    Y_train = y_train,
    Y_test = y_test,
    scaling = list(
      preprocessing_version = SHANGHAI_PREPROCESSING_VERSION,
      x_mean_curve = station_mean_curve,
      x_sd = station_sd,
      y_mean_curve = target_mean_curve,
      y_sd = target_sd
    )
  )
}

build_shanghai_basis_design <- function(target_obj, kx = 5, ky = 5) {
  # Convert Shanghai inflow/outflow arrays to the generic FSSGL design.
  time_numeric <- seq_len(nrow(target_obj$time_index))
  basis_x <- make_bspline_basis(time_numeric, df = kx)
  basis_y <- make_bspline_basis(time_numeric, df = ky)

  x_coef <- project_predictor_array(target_obj$X_train, basis_x)
  y_coef <- project_curves_l2(target_obj$Y_train, basis_y)

  design_obj <- build_fof_design(x_coef, n_response_basis = ky)
  covariate_blocks <- build_functional_groups(design_obj$predictor_blocks)

  station_groups <- data.table::data.table(
    predictor_id = seq_along(target_obj$predictor_station$station),
    station = target_obj$predictor_station$station,
    name = target_obj$predictor_station$name,
    lines = target_obj$predictor_station$lines,
    n_lines = target_obj$predictor_station$n_lines,
    is_transfer = target_obj$predictor_station$n_lines > 1
  )
  station_groups <- cbind(covariate_blocks, station_groups)

  list(
    target = target_obj$target,
    kx = kx,
    ky = ky,
    basis_x = basis_x,
    basis_y = basis_y,
    basis_geometry = "orthonormal_trapezoid",
    x_coef = x_coef,
    y_coef = y_coef,
    y_vec = as.vector(y_coef),
    design = design_obj$design,
    predictor_blocks = design_obj$predictor_blocks,
    station_groups = station_groups,
    line_membership = target_obj$line_membership
  )
}

evaluate_shanghai_fit <- function(result) {
  # Test-set diagnostic for one Shanghai target. This is intentionally kept
  # outside R/fssgl/ because it depends on the MetroFlow train/test arrays.
  target_obj <- result$target_obj
  design_obj <- result$design_obj
  fit <- result$fit

  x_coef_test <- project_predictor_array(target_obj$X_test, design_obj$basis_x)
  y_coef_test <- project_curves_l2(target_obj$Y_test, design_obj$basis_y)
  design_test <- build_fof_design(x_coef_test, n_response_basis = design_obj$ky)$design

  pred_coef_vec <- drop(design_test %*% fit$beta)
  pred_coef <- matrix(pred_coef_vec, nrow = nrow(y_coef_test), ncol = design_obj$ky)
  pred_curves <- pred_coef %*% t(design_obj$basis_y)

  list(
    test_coeff_rmse = sqrt(mean((pred_coef - y_coef_test)^2)),
    test_curve_rmse = sqrt(mean((pred_curves - target_obj$Y_test)^2)),
    test_curve_mae = mean(abs(pred_curves - target_obj$Y_test))
  )
}

fit_shanghai_target_fssgl <- function(
  target_name,
  params,
  processed_dir,
  metadata_dir,
  kx = 5,
  ky = 5,
  train_days = NULL,
  test_days = NULL,
  shanghai_data = NULL
) {
  # Fit one Shanghai target station and attach station/line labels to the
  # generic FSSGL output for downstream reporting scripts.
  target_obj <- prepare_shanghai_target(
    target_name,
    processed_dir,
    metadata_dir,
    train_days = train_days,
    test_days = test_days,
    shanghai_data = shanghai_data
  )
  design_obj <- build_shanghai_basis_design(target_obj, kx = kx, ky = ky)
  block_size <- design_obj$kx * design_obj$ky
  if (is.null(params$relaxed_tol)) params$relaxed_tol <- 0.012
  if (is.null(params$objective_tail_tol)) params$objective_tail_tol <- 0.08
  if (is.null(params$posterior_cutoff)) params$posterior_cutoff <- 0.5

  beta_warm <- ridge_dual_fit(design_obj$design, design_obj$y_vec, lambda = params$ridge_lambda)

  fit <- fit_fssgl_solver(
    x = design_obj$design,
    y = design_obj$y_vec,
    covariate_blocks = design_obj$station_groups,
    structural_membership = design_obj$line_membership,
    block_size = block_size,
    lambda_covariate_spike = params$lambda_station_spike,
    lambda_covariate_slab = params$lambda_station_slab,
    theta_covariate = params$theta_station,
    lambda_group_spike = params$lambda_line_spike,
    lambda_group_slab = params$lambda_line_slab,
    theta_group = params$theta_line,
    max_iter = params$max_iter,
    tol = params$tol,
    relaxed_tol = params$relaxed_tol,
    objective_tail_tol = params$objective_tail_tol,
    step_multiplier = params$step_multiplier,
    covariate_penalty_scale = params$station_penalty_scale,
    group_penalty_scale = params$line_penalty_scale,
    posterior_cutoff = params$posterior_cutoff,
    beta_init = beta_warm,
    verbose = FALSE
  )
  fit$history <- data.table::copy(fit$history)
  data.table::setnames(
    fit$history,
    c("selected_covariate_count", "selected_group_count"),
    c("selected_station_count", "selected_line_count"),
    skip_absent = TRUE
  )

  fit$station_posterior <- data.table::copy(fit$covariate_posterior)
  fit$line_posterior <- data.table::copy(fit$structure_posterior)
  data.table::setnames(
    fit$line_posterior,
    c("group_id", "group_norm"),
    c("line", "line_norm"),
    skip_absent = TRUE
  )
  fit$line_membership <- data.table::copy(fit$structural_membership)

  list(target_obj = target_obj, design_obj = design_obj, fit = fit)
}
