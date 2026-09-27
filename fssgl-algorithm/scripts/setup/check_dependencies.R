root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(root, "R/dependencies.R"))

dependencies <- read_project_dependencies(root)
r_dependency <- dependencies[dependencies$package == "R", , drop = FALSE]
if (nrow(r_dependency) != 1L || is.na(r_dependency$minimum_version[[1L]])) {
  stop("DESCRIPTION must declare one minimum R version.", call. = FALSE)
}

problems <- character()
r_minimum <- r_dependency$minimum_version[[1L]]
if (utils::compareVersion(as.character(getRversion()), r_minimum) < 0) {
  problems <- c(
    problems,
    paste0("R ", getRversion(), " is older than required ", r_minimum)
  )
}

packages <- project_installable_dependencies(root)
for (index in seq_len(nrow(packages))) {
  package <- packages$package[[index]]
  minimum_version <- packages$minimum_version[[index]]
  if (!requireNamespace(package, quietly = TRUE)) {
    problems <- c(problems, paste0(package, " is not installed"))
    next
  }
  installed <- as.character(utils::packageVersion(package))
  if (
    !is.na(minimum_version) &&
      utils::compareVersion(installed, minimum_version) < 0
  ) {
    problems <- c(
      problems,
      paste0(package, " ", installed, " is older than required ", minimum_version)
    )
  }
}

if (length(problems) > 0L) {
  stop(paste(problems, collapse = "\n"), call. = FALSE)
}

cat("Dependency check passed under", R.version.string, "\n")
cat("R minimum:", r_minimum, "\n")
for (index in seq_len(nrow(packages))) {
  package <- packages$package[[index]]
  minimum <- packages$minimum_version[[index]]
  suffix <- if (is.na(minimum)) "" else paste0(" (minimum ", minimum, ")")
  cat(package, " ", as.character(utils::packageVersion(package)), suffix, "\n", sep = "")
}
