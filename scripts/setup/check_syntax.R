root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
r_files <- list.files(
  root,
  pattern = "[.]R$",
  recursive = TRUE,
  full.names = TRUE
)

failures <- character()
for (file in r_files) {
  tryCatch(
    parse(file = file),
    error = function(error) {
      relative_file <- substring(file, nchar(root) + 2L)
      failures <<- c(failures, paste0(relative_file, ": ", conditionMessage(error)))
    }
  )
}

if (length(failures) > 0L) {
  stop("R syntax check failed:\n", paste(failures, collapse = "\n"), call. = FALSE)
}

cat("R syntax check passed for", length(r_files), "files.\n")
