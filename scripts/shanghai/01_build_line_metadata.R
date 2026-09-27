library(data.table)
if (!requireNamespace("rvest", quietly = TRUE)) {
  stop("Station-line web extraction requires rvest; install it before running this preprocessing script.")
}
library(rvest)

root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(root, "data/raw/shanghai_metroflow/MetroFlow")
metadata_dir <- file.path(root, "data/metadata/shanghai_metroflow")
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

station_file <- file.path(raw_dir, "stationInfo.csv")
if (!file.exists(station_file)) {
  stop("Missing stationInfo.csv. Download/extract MetroFlow first.")
}

station_info <- fread(station_file)
setnames(station_info, old = names(station_info)[1], new = "row_id")

normalize_name <- function(x) {
  x <- gsub("\\[[^]]*\\]", "", x)
  x <- gsub("\\([^)]*\\)", "", x)
  x <- gsub("’", "'", x)
  x <- gsub("–", "-", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

station_info[, name_norm := normalize_name(name)]

station_alias <- data.table(
  metroflow_name = c(
    "Jinjiang Amusement Park",
    "Zhongshan North Road",
    "Shanghai Circus City",
    "East Xujing",
    "Zhangjiang Hi-Tech Park",
    "Chuangxin Central Road",
    "Haitian 3rd Road",
    "Pudong International Airport",
    "Lingping Road",
    "Jingping Road",
    "South Waigaoqiao Free Trade Zone Station",
    "North Waigaoqiao Free Trade Zone Station",
    "Liu Hang",
    "Dahua San Road",
    "Great World",
    "South Songjiang Railway Station",
    "Zuibaichi",
    "Songjiang New City",
    "Hongqiao Terminal 1",
    "Hongqiao Terminal 2",
    "Jiaotong University",
    "South Huangpi Road",
    "Xintiandi",
    "Longhua Middle Road",
    "Jiading New City",
    "Hesha Hangcheng",
    "Hangtou East",
    "Shanghai Wild Animal Park",
    "Huinan East"
  ),
  wiki_name = c(
    "Jinjiang Park",
    "North Zhongshan Road",
    "Shanghai Circus World",
    "East Xujing",
    "Zhangjiang High Technology Park",
    "Middle Chuangxin Road",
    "Haitiansan Road",
    "Pudong Airport Terminal 1&2",
    "Linping Road",
    "Jinping Road",
    "South Waigaoqiao Free Trade Zone",
    "North Waigaoqiao Free Trade Zone",
    "Liuhang",
    "Dahuasan Road",
    "Dashijie",
    "Shanghai Songjiang Railway Station",
    "Zuibaichi Park",
    "Songjiang Xincheng",
    "Hongqiao Airport Terminal 1",
    "Hongqiao Airport Terminal 2",
    "Jiao Tong University",
    "Site of the First CPC National Congress · South Huangpi Road",
    "Site of the First CPC National Congress · Xintiandi",
    "Middle Longhua Road",
    "Jiading Xincheng",
    "Heshahangcheng",
    "East Hangtou",
    "Wild Animal Park",
    "East Huinan"
  )
)
station_alias[, metroflow_norm := normalize_name(metroflow_name)]
station_alias[, wiki_norm := normalize_name(wiki_name)]

station_info_key <- merge(
  station_info,
  station_alias[, .(metroflow_norm, wiki_norm)],
  by.x = "name_norm",
  by.y = "metroflow_norm",
  all.x = TRUE
)
station_info_key[, match_norm := fifelse(is.na(wiki_norm), name_norm, wiki_norm)]

wiki_url <- "https://en.wikipedia.org/wiki/List_of_Shanghai_Metro_stations"
doc <- read_html(wiki_url)
tables <- html_table(doc, fill = TRUE)

# The page tables after the transport navbox correspond to Line 1, Line 2, ...,
# Line 18, and Pujiang line. The MetroFlow period is May-August 2017, so we keep
# only lines operating during that period: Lines 1-13 and 16.
line_table_index <- data.table(
  line = paste0("Line ", c(1:13, 16)),
  table_index = c(2:14, 17)
)

extract_line_stations <- function(tbl, line_name) {
  cols <- names(tbl)
  station_col <- which(grepl("Station Name", cols))[1]
  opened_col <- which(cols == "Opened")[1]

  if (is.na(station_col)) {
    return(data.table())
  }

  out <- data.table(
    line = line_name,
    station_name_wiki = normalize_name(tbl[[station_col]]),
    opened = if (!is.na(opened_col)) as.character(tbl[[opened_col]]) else NA_character_
  )
  out <- out[
    station_name_wiki != "" &
      !station_name_wiki %in% c("English", "Station Name", "Opened", "Service Routes")
  ]
  unique(out)
}

line_membership <- rbindlist(
  lapply(seq_len(nrow(line_table_index)), function(i) {
    extract_line_stations(
      tables[[line_table_index$table_index[i]]],
      line_table_index$line[i]
    )
  }),
  fill = TRUE
)

line_membership[, name_norm := normalize_name(station_name_wiki)]

manual_membership <- data.table(
  name = c("East Xujing", "Dongchang Road", "New Jiangwan City", "Expo Avenue"),
  line = c("Line 2", "Line 2", "Line 10", "Line 13"),
  station_name_wiki = c("East Xujing", "Dongchang Road", "New Jiangwan City", "Expo Avenue"),
  opened = NA_character_
)
manual_membership[, name_norm := normalize_name(name)]

membership <- merge(
  station_info_key[, .(stationID, name, lon, lat, neighbour, match_norm)],
  line_membership[, .(line, station_name_wiki, opened, name_norm)],
  by.x = "match_norm",
  by.y = "name_norm",
  all.x = FALSE,
  all.y = FALSE
)

manual_join <- merge(
  station_info_key[, .(stationID, name, lon, lat, neighbour, match_norm)],
  manual_membership[, .(line, station_name_wiki, opened, name_norm)],
  by.x = "match_norm",
  by.y = "name_norm",
  all.x = FALSE,
  all.y = FALSE
)

membership <- rbindlist(list(membership, manual_join), use.names = TRUE, fill = TRUE)

membership <- unique(membership[
  order(line, stationID),
  .(stationID, name, line, station_name_wiki, opened, lon, lat, neighbour)
])

station_lines <- membership[
  ,
  .(
    lines = paste(sort(unique(line)), collapse = ";"),
    n_lines = uniqueN(line),
    is_transfer = uniqueN(line) > 1
  ),
  by = .(stationID, name)
]

station_line_map <- merge(
  station_info[, .(stationID, name, lon, lat, neighbour)],
  station_lines,
  by = c("stationID", "name"),
  all.x = TRUE
)

station_line_map[is.na(lines), `:=`(
  lines = "",
  n_lines = 0L,
  is_transfer = FALSE
)]

unmatched_station_info <- station_line_map[n_lines == 0]
unmatched_wiki_stations <- line_membership[
  !name_norm %in% station_info_key$match_norm,
  .(line, station_name_wiki, opened)
]

fwrite(station_info[, !"name_norm"], file.path(metadata_dir, "station_info.csv"))
fwrite(station_alias[, .(metroflow_name, wiki_name)], file.path(metadata_dir, "station_name_aliases.csv"))
fwrite(manual_membership[, .(name, line, note = "manual verified fallback for unmatched public table name")],
       file.path(metadata_dir, "station_line_manual_fallback.csv"))
fwrite(membership, file.path(metadata_dir, "station_line_membership_long.csv"))
fwrite(station_line_map, file.path(metadata_dir, "station_line_map.csv"))
fwrite(unmatched_station_info, file.path(metadata_dir, "station_unmatched_to_wikipedia.csv"))
fwrite(unmatched_wiki_stations, file.path(metadata_dir, "wikipedia_stations_unmatched_to_metroflow.csv"))

cat("stationInfo rows:", nrow(station_info), "\n")
cat("matched station-line memberships:", nrow(membership), "\n")
cat("stations with at least one line:", sum(station_line_map$n_lines > 0), "\n")
cat("stations without matched line:", nrow(unmatched_station_info), "\n")
cat("transfer stations:", sum(station_line_map$is_transfer), "\n")
cat("metadata_dir:", metadata_dir, "\n")
