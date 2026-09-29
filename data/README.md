# Frozen Shanghai analysis input

This directory contains a processed snapshot used by the Shanghai MetroFlow
analysis. It includes the station/line metadata and the cleaned array input
(`processed/shanghai_metroflow/shanghai_metroflow_clean.RData`) together with
its station, date, and time indices.

The raw MetroFlow passenger-flow files are deliberately not included. The
directory also excludes fitted objects, result tables, and figures. The
processed snapshot is intended to make the analysis input deterministic; the
Shanghai scripts can use it directly for fitting, stability selection, rolling
evaluation, and diagnostics.

The metadata files are retained because the weighted-membership application
uses the station-to-line mapping, including transfer-station membership. The
snapshot should be used with the method release recorded in `R/parameters.R`.

The fitting and downstream diagnostic scripts load this snapshot directly.
`scripts/shanghai/03_validate_processed_data.R` additionally checks the raw
flow grid and therefore requires the original raw files; it is not a
processed-only validation command.
