root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required_markers <- c("DESCRIPTION", "SSGL.Rproj", "R", "scripts", "tests")
missing_markers <- required_markers[!file.exists(file.path(root, required_markers))]
if (length(missing_markers) > 0L) {
  stop(
    "Run the doctor from the project root; missing: ",
    paste(missing_markers, collapse = ", "),
    call. = FALSE
  )
}

cat("Project root:", root, "\n")
cat("Runtime:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("Library paths:\n", paste0("  - ", .libPaths(), collapse = "\n"), "\n", sep = "")

locale_names <- c("LANG", "LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME")
locale_variables <- stats::setNames(Sys.getenv(locale_names), locale_names)
active_locale_variables <- locale_variables[nzchar(locale_variables)]
if (length(active_locale_variables) > 0L) {
  cat(
    "Locale environment:",
    paste(names(active_locale_variables), active_locale_variables, sep = "=", collapse = ", "),
    "\n"
  )
}
if (.Platform$OS.type == "windows" && any(active_locale_variables == "C.UTF-8")) {
  cat(
    "NOTE: C.UTF-8 is a Unix locale name and can warn on Windows. ",
    "Use scripts/run.ps1 to clear it for the R child process.\n",
    sep = ""
  )
}

cat("Project writable:", file.access(root, mode = 2L) == 0L, "\n")

tool_status <- function(tool) {
  path <- Sys.which(tool)
  if (nzchar(path)) path else "not found (optional)"
}
cat("pdflatex:", tool_status("pdflatex"), "\n")
cat("git:", tool_status("git"), "\n")

source(file.path(root, "scripts/setup/check_dependencies.R"), local = TRUE)
cat("Environment doctor passed.\n")
