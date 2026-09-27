root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(root, "R/dependencies.R"))

dependencies <- read_project_dependencies(root)
r_dependency <- dependencies[dependencies$package == "R", , drop = FALSE]
if (nrow(r_dependency) != 1L || is.na(r_dependency$minimum_version[[1L]])) {
  stop("DESCRIPTION must declare one minimum R version.", call. = FALSE)
}
if (
  utils::compareVersion(
    as.character(getRversion()),
    r_dependency$minimum_version[[1L]]
  ) < 0
) {
  stop(
    "R ", getRversion(), " is older than required ",
    r_dependency$minimum_version[[1L]],
    call. = FALSE
  )
}

packages <- project_installable_dependencies(root)
needs_install <- vapply(seq_len(nrow(packages)), function(index) {
  package <- packages$package[[index]]
  minimum_version <- packages$minimum_version[[index]]
  !requireNamespace(package, quietly = TRUE) ||
    (!is.na(minimum_version) &&
      utils::compareVersion(
        as.character(utils::packageVersion(package)),
        minimum_version
      ) < 0)
}, logical(1L))

if (!any(needs_install)) {
  cat("All required packages already satisfy DESCRIPTION.\n")
} else {
  repository <- Sys.getenv("FSSGL_CRAN_REPO", unset = "https://cloud.r-project.org")
  to_install <- packages$package[needs_install]
  cat("Installing:", paste(to_install, collapse = ", "), "\n")
  utils::install.packages(to_install, repos = repository)
  source(file.path(root, "scripts/setup/check_dependencies.R"), local = TRUE)
}
