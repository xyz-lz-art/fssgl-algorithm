library(data.table)

root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(root, "data/raw/shanghai_metroflow/MetroFlow")
metadata_dir <- file.path(root, "data/metadata/shanghai_metroflow")
processed_dir <- file.path(root, "data/processed/shanghai_metroflow")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

flow_file <- file.path(raw_dir, "metroData_InOutFlow.csv")
station_file <- file.path(metadata_dir, "station_line_map.csv")

if (!file.exists(flow_file)) {
  stop("Missing metroData_InOutFlow.csv. Download/extract MetroFlow first.")
}
if (!file.exists(station_file)) {
  stop("Missing station_line_map.csv. Run 01_build_line_metadata.R first.")
}

flow <- fread(flow_file)
setnames(flow, trimws(names(flow)))

required_cols <- c("date", "timeslot", "startTime", "endTime", "station", "inFlow", "outFlow")
missing_cols <- setdiff(required_cols, names(flow))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

station_map <- fread(station_file)

if (anyNA(flow[, ..required_cols])) {
  stop("Raw flow data contain missing values in required columns.")
}
if (any(!is.finite(flow$inFlow)) || any(!is.finite(flow$outFlow))) {
  stop("Raw flow values must be finite.")
}
if (any(flow$inFlow < 0) || any(flow$outFlow < 0)) {
  stop("Raw flow values cannot be negative.")
}

dates <- sort(unique(flow$date))
stations <- sort(unique(flow$station))

station_index <- data.table(station = stations, station_pos = seq_along(stations))
date_index <- data.table(date = dates, day_pos = seq_along(dates))
time_index <- unique(flow[, .(startTime, endTime)])
setorder(time_index, startTime, endTime)
time_index[, time_pos := seq_len(.N)]

flow <- merge(flow, station_index, by = "station", all.x = TRUE)
flow <- merge(flow, date_index, by = "date", all.x = TRUE)
flow <- merge(flow, time_index[, .(startTime, endTime, time_pos)], by = c("startTime", "endTime"), all.x = TRUE)

key_cols <- c("station_pos", "day_pos", "time_pos")
duplicate_keys <- flow[, .N, by = key_cols][N > 1L]
if (nrow(duplicate_keys) > 0L) {
  stop("Raw flow data contain duplicate station/date/time records.")
}

n_station <- length(stations)
n_day <- length(dates)
n_time <- nrow(time_index)

expected_rows <- n_station * n_day * n_time
if (nrow(flow) != expected_rows) {
  stop(
    "Raw flow grid is incomplete: expected ", expected_rows,
    " station/date/time records but found ", nrow(flow), "."
  )
}

arr_in <- array(
  NA_real_,
  dim = c(n_station, n_day, n_time),
  dimnames = list(as.character(stations), as.character(dates), as.character(time_index$startTime))
)
arr_out <- arr_in

idx <- cbind(flow$station_pos, flow$day_pos, flow$time_pos)
arr_in[idx] <- flow$inFlow
arr_out[idx] <- flow$outFlow
if (anyNA(arr_in) || anyNA(arr_out)) {
  stop("Raw flow grid has missing station/date/time combinations.")
}

station_metadata <- merge(
  station_index,
  station_map,
  by.x = "station",
  by.y = "stationID",
  all.x = TRUE
)
setorder(station_metadata, station_pos)
if (anyNA(station_metadata$name) || anyNA(station_metadata$lines)) {
  stop("Station metadata do not cover every station in the flow arrays.")
}

save(
  arr_in,
  arr_out,
  station_metadata,
  date_index,
  time_index,
  file = file.path(processed_dir, "shanghai_metroflow_clean.RData")
)

fwrite(station_metadata, file.path(processed_dir, "station_metadata.csv"))
fwrite(date_index, file.path(processed_dir, "date_index.csv"))
fwrite(time_index, file.path(processed_dir, "time_index.csv"))

cat("arr_in dim:", paste(dim(arr_in), collapse = " x "), "\n")
cat("arr_out dim:", paste(dim(arr_out), collapse = " x "), "\n")
cat("stations:", n_station, "\n")
cat("days:", n_day, "\n")
cat("timeslots:", n_time, "\n")
cat("processed_dir:", processed_dir, "\n")
