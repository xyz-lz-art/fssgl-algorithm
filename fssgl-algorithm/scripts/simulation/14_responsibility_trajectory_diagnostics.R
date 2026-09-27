# Diagnose how isolated within-group signals move through the FSSGL spike/slab
# transitions. Tracing is observational and leaves every solver update unchanged.

library(data.table)

source("R/parameters.R")
source("R/fssgl/basis_design.R")
source("R/fssgl/penalty_weights.R")
source("R/fssgl/weighted_membership_solver.R")
source("R/fssgl/solver.R")
source("R/fssgl/simulation_design.R")

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "results/tables/simulation/v2_formal")
figure_dir <- file.path(root, "results/figures/simulation/v2_formal")
processed_dir <- file.path(root, "data/processed/simulation/v2_formal")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

n_repetitions <- as.integer(Sys.getenv("FSSGL_RESPONSIBILITY_REPS", "100"))
if (!is.finite(n_repetitions) || n_repetitions < 1L) {
  stop("FSSGL_RESPONSIBILITY_REPS must be a positive integer.")
}
parameters <- fssgl_parameters()
scenarios <- data.table(
  scenario_id = c("baseline", "within_group_sparse"),
  scenario_label = c("Two active per group", "One active per group"),
  active_covariates_per_group = c(2L, 1L)
)
scenario_grid <- scenarios[, .(rep = seq_len(n_repetitions)), by = names(scenarios)]
checkpoint_file <- file.path(processed_dir, "fssgl_v2_responsibility_trace_running.rds")
fit_file <- file.path(table_dir, "fssgl_v2_responsibility_diagnostic_replicates.csv")
trace_file <- file.path(processed_dir, "fssgl_v2_responsibility_trajectories.rds")
state <- if (file.exists(checkpoint_file)) {
  readRDS(checkpoint_file)
} else if (file.exists(fit_file) && file.exists(trace_file)) {
  list(fits = fread(fit_file), traces = readRDS(trace_file))
} else {
  list(fits = data.table(), traces = data.table())
}
state$fits <- state$fits[
  scenario_id %in% scenarios$scenario_id & rep <= n_repetitions
]
state$traces <- state$traces[
  scenario_id %in% scenarios$scenario_id & rep <= n_repetitions
]

for (scenario_index in seq_len(nrow(scenario_grid))) {
  scenario <- scenario_grid[scenario_index]
  if (nrow(state$fits) && any(
    state$fits$scenario_id == scenario$scenario_id & state$fits$rep == scenario$rep
  )) next

  seed <- 2026103000L + scenario$rep
  cat(
    "Responsibility diagnostic: ", scenario$scenario_id,
    ", rep=", scenario$rep, "/", n_repetitions, "\n", sep = ""
  )
  dgp <- generate_fssgl_simulation_dgp(
    n_train = 40L,
    n_test = 24L,
    n_covariates = 20L,
    kx = 4L,
    ky = 4L,
    n_groups = 4L,
    n_active_groups = 2L,
    n_active_covariates_per_group = scenario$active_covariates_per_group,
    rho_common = 0.10,
    rho_group = 0.30,
    snr = 2,
    surface_smoothness = "moderate",
    signal_scale = 1,
    seed = seed
  )
  started <- proc.time()[["elapsed"]]
  fit <- fit_fssgl(
    x_coef = dgp$x_train,
    y_coef = dgp$y_train,
    structural_membership = dgp$structural_membership,
    parameters = parameters,
    trace_responsibilities = TRUE,
    verbose = FALSE
  )
  runtime_sec <- proc.time()[["elapsed"]] - started
  metrics <- evaluate_fssgl_simulation_fit(
    fit, dgp, parameters$posterior_cutoff, runtime_sec
  )
  fit_row <- cbind(
    scenario[, .(
      scenario_id, scenario_label, active_covariates_per_group, rep
    )],
    data.table(
      seed = seed,
      strict_converged = fit$fit$convergence$strict_converged,
      final_iter = fit$fit$convergence$final_iter,
      theta_covariate_final = fit$fit$hyperparameters$theta_covariate_final,
      theta_group_final = fit$fit$hyperparameters$theta_group_final
    ),
    as.data.table(as.list(metrics))
  )

  covariate_roles <- merge(
    dgp$structural_membership[, .(covariate_id, group_id)],
    dgp$truth_covariates[, .(covariate_id, active)],
    by = "covariate_id"
  )
  active_groups <- dgp$truth_groups[active == TRUE, group_id]
  covariate_roles[, role := fifelse(
    active,
    "active",
    fifelse(group_id %in% active_groups, "inactive_neighbor", "inactive_group")
  )]
  covariate_roles[, `:=`(
    unit_id = as.character(covariate_id),
    unit_number = covariate_id
  )]
  group_roles <- dgp$truth_groups[, .(
    group_id,
    role = fifelse(active, "active", "inactive")
  )]
  group_roles[, `:=`(
    unit_id = as.character(group_id),
    unit_number = as.integer(sub("group_", "", group_id))
  )]

  trace <- copy(fit$fit$responsibility_trace)
  covariate_trace <- merge(
    trace[level == "covariate"],
    covariate_roles[, .(unit_id, unit_number, role)],
    by = "unit_id"
  )
  group_trace <- merge(
    trace[level == "group"],
    group_roles[, .(unit_id, unit_number, role)],
    by = "unit_id"
  )
  trace <- rbindlist(list(covariate_trace, group_trace), use.names = TRUE, fill = TRUE)
  trace[, `:=`(
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    rep = scenario$rep,
    seed = seed,
    final_iter = fit$fit$convergence$final_iter
  )]
  trace[, progress := fifelse(
    final_iter <= 1L, 1, (iter - 1) / (final_iter - 1)
  )]
  trace[, progress_bin := round(progress * 20) / 20]
  setcolorder(
    trace,
    c(
      "scenario_id", "scenario_label", "rep", "seed", "iter", "final_iter",
      "progress", "progress_bin", "level", "unit_id", "unit_number", "role",
      setdiff(names(trace), c(
        "scenario_id", "scenario_label", "rep", "seed", "iter", "final_iter",
        "progress", "progress_bin", "level", "unit_id", "unit_number", "role"
      ))
    )
  )
  state$fits <- rbindlist(list(state$fits, fit_row), use.names = TRUE, fill = TRUE)
  state$traces <- rbindlist(list(state$traces, trace), use.names = TRUE, fill = TRUE)
  saveRDS(state, checkpoint_file)
}

setorder(state$fits, scenario_id, rep)
setorder(state$traces, scenario_id, rep, iter, level, unit_number)
fwrite(
  state$fits,
  fit_file
)
saveRDS(
  state$traces,
  trace_file
)

active_iteration <- state$traces[role == "active", .(
  posterior_slab_prob = mean(posterior_slab_prob),
  beta_norm = mean(beta_norm),
  transition_norm = mean(transition_norm),
  norm_to_transition = mean(norm_to_transition),
  theta = mean(theta)
), by = .(scenario_id, scenario_label, rep, level, iter, progress)]
progress_grid <- seq(0, 1, by = 0.05)
active_path <- active_iteration[order(progress), {
  carry_forward <- function(value) {
    stats::approx(
      x = progress,
      y = value,
      xout = progress_grid,
      method = "constant",
      f = 0,
      rule = 2,
      ties = mean
    )$y
  }
  .(
    progress_bin = progress_grid,
    posterior_slab_prob = carry_forward(posterior_slab_prob),
    beta_norm = carry_forward(beta_norm),
    transition_norm = carry_forward(transition_norm),
    norm_to_transition = carry_forward(norm_to_transition),
    theta = carry_forward(theta)
  )
}, by = .(scenario_id, scenario_label, rep, level)]
trajectory_summary <- active_path[, .(
  n_repetitions = uniqueN(rep),
  posterior_median = median(posterior_slab_prob),
  posterior_q25 = quantile(posterior_slab_prob, 0.25),
  posterior_q75 = quantile(posterior_slab_prob, 0.75),
  norm_ratio_median = median(norm_to_transition),
  norm_ratio_q25 = quantile(norm_to_transition, 0.25),
  norm_ratio_q75 = quantile(norm_to_transition, 0.75),
  theta_median = median(theta)
), by = .(scenario_id, scenario_label, level, progress_bin)]
setorder(trajectory_summary, level, scenario_id, progress_bin)
fwrite(
  trajectory_summary,
  file.path(table_dir, "fssgl_v2_responsibility_trajectory_summary.csv")
)

final_trace <- state$traces[iter == final_iter]
final_summary <- final_trace[, .(
  n_repetitions = uniqueN(rep),
  n_units = .N,
  posterior_mean = mean(posterior_slab_prob),
  posterior_sd = sd(posterior_slab_prob),
  posterior_median = median(posterior_slab_prob),
  posterior_q25 = quantile(posterior_slab_prob, 0.25),
  posterior_q75 = quantile(posterior_slab_prob, 0.75),
  beta_norm_mean = mean(beta_norm),
  transition_norm_mean = mean(transition_norm),
  norm_ratio_mean = mean(norm_to_transition),
  norm_ratio_median = median(norm_to_transition)
), by = .(scenario_id, scenario_label, level, role)]
setorder(final_summary, level, role, scenario_id)
fwrite(
  final_summary,
  file.path(table_dir, "fssgl_v2_responsibility_final_summary.csv")
)

plot_trajectory_panel <- function(level_value, metric_prefix, ylab, reference) {
  panel <- trajectory_summary[level == level_value]
  median_name <- paste0(metric_prefix, "_median")
  lower_name <- paste0(metric_prefix, "_q25")
  upper_name <- paste0(metric_prefix, "_q75")
  ylim <- range(c(panel[[lower_name]], panel[[upper_name]], reference), finite = TRUE)
  padding <- 0.05 * diff(ylim)
  if (!is.finite(padding) || padding == 0) padding <- 0.05
  plot(
    NA, xlim = c(0, 1), ylim = ylim + c(-padding, padding),
    xlab = "Normalized iteration progress", ylab = ylab, las = 1
  )
  abline(h = reference, col = "grey55", lty = 2)
  colors <- c(baseline = "#176B67", within_group_sparse = "#C65D32")
  for (scenario_value in names(colors)) {
    values <- panel[scenario_id == scenario_value][order(progress_bin)]
    polygon(
      c(values$progress_bin, rev(values$progress_bin)),
      c(values[[lower_name]], rev(values[[upper_name]])),
      border = NA,
      col = grDevices::adjustcolor(colors[[scenario_value]], alpha.f = 0.18)
    )
    lines(
      values$progress_bin, values[[median_name]],
      col = colors[[scenario_value]], lwd = 2
    )
  }
}

figure_file <- file.path(figure_dir, "fssgl_v2_responsibility_trajectories.pdf")
pdf(figure_file, width = 8.2, height = 6.5, family = "Helvetica")
par(mfrow = c(2, 2), mar = c(4.2, 4.3, 2.5, 1), las = 1)
plot_trajectory_panel("covariate", "posterior", "Active predictor responsibility", 0.5)
title("(a) Predictor slab responsibility")
legend(
  "topleft",
  legend = scenarios$scenario_label,
  col = c("#176B67", "#C65D32"), lwd = 2, bty = "n", cex = 0.85
)
plot_trajectory_panel("group", "posterior", "Active group responsibility", 0.5)
title("(b) Group slab responsibility")
plot_trajectory_panel("covariate", "norm_ratio", "Predictor norm / transition norm", 1)
title("(c) Predictor transition ratio")
plot_trajectory_panel("group", "norm_ratio", "Group norm / transition norm", 1)
title("(d) Group transition ratio")
dev.off()

saveRDS(
  new_experiment_manifest(
    experiment_id = "fssgl_v2_responsibility_trajectory_diagnostic",
    parameters = parameters,
    design = list(
      scenarios = scenarios,
      n_repetitions = n_repetitions,
      paired_seed_rule = "2026103000 + rep",
      trace_phase = "post-update state entering the next outer iteration",
      active_unit_aggregation = paste(
        "within-replicate mean, last-state interpolation on 21 normalized",
        "progress points, then median and IQR"
      )
    ),
    algorithm_version = FSSGL_ALGORITHM_VERSION
  ),
  file.path(processed_dir, "fssgl_v2_responsibility_diagnostic_manifest.rds")
)
if (nrow(state$fits) == nrow(scenario_grid)) unlink(checkpoint_file)

print(state$fits[, .(
  n = .N,
  strict_convergence_rate = mean(strict_converged),
  predictor_tpr = mean(covariate_tpr),
  group_tpr = mean(group_tpr),
  theta_covariate = mean(theta_covariate_final),
  theta_group = mean(theta_group_final)
), by = .(scenario_id, scenario_label)])
print(final_summary[role == "active"])
