library(data.table)

root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(root, "data/raw/shanghai_metroflow/MetroFlow")
metadata_dir <- file.path(root, "data/metadata/shanghai_metroflow")
processed_dir <- file.path(root, "data/processed/shanghai_metroflow")

required_files <- c(
  file.path(raw_dir, "stationInfo.csv"),
  file.path(raw_dir, "metroData_InOutFlow.csv"),
  file.path(metadata_dir, "station_line_map.csv"),
  file.path(metadata_dir, "station_line_membership_long.csv"),
  file.path(processed_dir, "shanghai_metroflow_clean.RData")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

station_line_map <- fread(file.path(metadata_dir, "station_line_map.csv"))
membership <- fread(file.path(metadata_dir, "station_line_membership_long.csv"))
load(file.path(processed_dir, "shanghai_metroflow_clean.RData"))
flow <- fread(file.path(raw_dir, "metroData_InOutFlow.csv"))
setnames(flow, trimws(names(flow)))

flow_key <- c("station", "date", "startTime", "endTime")
duplicate_flow_keys <- flow[, .N, by = flow_key][N > 1L, .N]
expected_flow_rows <- uniqueN(flow$station) * uniqueN(flow$date) *
  uniqueN(flow[, .(startTime, endTime)])

checks <- list(
  n_station_metadata = nrow(station_line_map),
  n_station_array = dim(arr_in)[1],
  n_day = dim(arr_in)[2],
  n_time = dim(arr_in)[3],
  stations_without_line = sum(station_line_map$n_lines == 0 | station_line_map$lines == ""),
  transfer_stations = sum(station_line_map$is_transfer),
  membership_rows = nrow(membership),
  line_count = uniqueN(membership$line),
  raw_flow_rows = nrow(flow),
  expected_flow_rows = expected_flow_rows,
  duplicate_flow_keys = duplicate_flow_keys,
  arr_in_na = sum(is.na(arr_in)),
  arr_out_na = sum(is.na(arr_out))
)

stopifnot(checks$n_station_metadata == 302)
stopifnot(checks$n_station_array == 302)
stopifnot(checks$stations_without_line == 0)
stopifnot(checks$raw_flow_rows == checks$expected_flow_rows)
stopifnot(checks$duplicate_flow_keys == 0)
stopifnot(checks$arr_in_na == 0)
stopifnot(checks$arr_out_na == 0)

cat("Shanghai MetroFlow validation passed.\n")
for (nm in names(checks)) {
  cat(nm, ":", checks[[nm]], "\n")
}
