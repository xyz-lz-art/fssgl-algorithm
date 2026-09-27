# Portable task runner for project checks and explicit workflow entry points.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("Unable to determine scripts/run.R location.", call. = FALSE)
}

script_file <- normalizePath(
  sub("^--file=", "", script_argument[[1L]]),
  winslash = "/",
  mustWork = TRUE
)
root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
setwd(root)

tasks <- list(
  dependencies = "scripts/setup/check_dependencies.R",
  install = "scripts/setup/install_dependencies.R",
  syntax = "scripts/setup/check_syntax.R",
  doctor = "scripts/setup/doctor.R",
  test = "tests/run_tests.R",
  validate_outputs = "scripts/validation/validate_all_results.R"
)

run_script <- function(path) {
  cat("\n==>", path, "\n")
  task_environment <- new.env(parent = globalenv())
  sys.source(file.path(root, path), envir = task_environment)
  invisible(TRUE)
}

print_help <- function() {
  cat(
    "Usage: Rscript scripts/run.R <task>\n\n",
    "Quality tasks:\n",
    "  doctor       Inspect runtime, paths, dependencies, and optional tools\n",
    "  dependencies Check minimum R and package versions\n",
    "  install      Install missing/outdated direct dependencies\n",
    "  syntax       Parse every R source and entry-point file\n",
    "  test         Run deterministic unit and smoke tests\n",
    "  check        Run doctor, syntax checks, and tests\n",
    "  validate     Validate frozen outputs and manuscript values\n",
    "  all          Run check and validate\n\n",
    "Formal workflows remain explicit under scripts/simulation/ and scripts/shanghai/.\n",
    sep = ""
  )
}

arguments <- commandArgs(trailingOnly = TRUE)
task <- if (length(arguments) == 0L) "help" else arguments[[1L]]

task_sequence <- switch(
  task,
  help = character(),
  doctor = tasks$doctor,
  dependencies = tasks$dependencies,
  install = tasks$install,
  syntax = tasks$syntax,
  test = tasks$test,
  check = c(tasks$doctor, tasks$syntax, tasks$test),
  validate = tasks$validate_outputs,
  all = c(
    tasks$doctor,
    tasks$syntax,
    tasks$test,
    tasks$validate_outputs
  ),
  stop("Unknown task: ", task, ". Run with 'help' to list tasks.", call. = FALSE)
)

if (task == "help") {
  print_help()
} else {
  for (path in task_sequence) run_script(path)
  cat("\nTask '", task, "' completed successfully.\n", sep = "")
}
