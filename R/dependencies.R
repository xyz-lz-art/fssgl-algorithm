read_project_dependencies <- function(root = getwd()) {
  description_file <- file.path(root, "DESCRIPTION")
  if (!file.exists(description_file)) {
    stop("DESCRIPTION is missing: ", description_file, call. = FALSE)
  }

  description <- read.dcf(description_file)
  parse_field <- function(field) {
    if (!field %in% colnames(description)) {
      return(data.frame(
        package = character(),
        minimum_version = character(),
        field = character()
      ))
    }
    entries <- trimws(strsplit(description[1L, field], ",", fixed = TRUE)[[1L]])
    packages <- trimws(sub("[(].*$", "", entries))
    has_version <- grepl("[(]>=[[:space:]]*", entries)
    minimum_versions <- rep(NA_character_, length(entries))
    minimum_versions[has_version] <- sub(
      ".*[(]>=[[:space:]]*([^)]*)[)].*",
      "\\1",
      entries[has_version]
    )
    data.frame(
      package = packages,
      minimum_version = minimum_versions,
      field = field,
      stringsAsFactors = FALSE
    )
  }

  dependencies <- rbind(parse_field("Depends"), parse_field("Imports"))
  dependencies[nzchar(dependencies$package), , drop = FALSE]
}

project_installable_dependencies <- function(root = getwd()) {
  dependencies <- read_project_dependencies(root)
  base_packages <- c(
    "base", "compiler", "datasets", "graphics", "grDevices", "grid",
    "methods", "parallel", "splines", "stats", "stats4", "tcltk", "tools", "utils"
  )
  dependencies[
    dependencies$package != "R" & !dependencies$package %in% base_packages,
    ,
    drop = FALSE
  ]
}
