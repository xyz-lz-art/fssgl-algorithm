# Functional Spike-and-Slab Group Lasso (FSSGL)

Research code for structured function-on-function regression. The formal
estimator selects functional predictors and their predefined, disjoint groups
while estimating bivariate coefficient surfaces. A separate weighted-membership
solver handles predictors that belong to multiple groups in the Shanghai metro
application. The current method identifier is `fssgl-orthonormal-2026.09`.

## Contents

- `R/fssgl/`: basis construction, penalties, formal solver, cross-validation,
  and synthetic data generation.
- `R/application/`: curve preprocessing and application helpers.
- `R/baselines/`: ridge, FPCA, group SCAD, and structured group-lasso comparisons.
- `R/diagnostics/` and `R/visualization/`: design and surface diagnostics.
- `scripts/simulation/`: simulation and figure-generating workflows.
- `scripts/shanghai/`: Shanghai preprocessing and analysis workflows.
- `scripts/setup/`, `scripts/run.R`, `tests/`: environment checks and deterministic
  tests.

The repository contains code only. Raw MetroFlow data, processed observations,
fitted objects, result tables, manuscripts, and third-party article PDFs are not
included. Scripts that validate completed experiments require their
corresponding saved outputs; they are not part of the first-run check.
Manuscript-specific value checks and appendix tables remain in the writing
project.

## Run locally

Use R 4.5 or newer. From the repository root:

```text
Rscript scripts/run.R install
Rscript scripts/run.R check
Rscript examples/quickstart.R
```

On Windows, `scripts/run.ps1 check` locates R and clears inherited Unix locale
settings for that process. `DESCRIPTION` declares the core R packages. `rvest`
is optional and needed only when rebuilding station--line metadata from web
pages. The quickstart generates a small synthetic data set, fits FSSGL once,
and prints selection and error diagnostics.

`R/parameters.R` holds the documented experiment configurations. The formal
solver uses non-overlapping groups; the Shanghai weighted-membership extension
is an engineering implementation for overlapping line memberships, not an
equivalent solver for the formal objective. The simulation scripts write to
`data/processed/` and `results/` when run. These paths are ignored by Git.

No open-source redistribution license has been granted. Third-party data and
papers remain under their own terms.
