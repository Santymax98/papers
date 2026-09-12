# Conditional Dimension Reduction for Vine Copula Probabilities

Code and numerical results accompanying the article
*Conditional Dimension Reduction for Vine Copula Probabilities*.

The repository contains the research implementation used in the experiments,
frozen numerical outputs, data provenance, and external validation. The
journal Supplementary Material is submitted separately and is not duplicated
here; the computational results underlying it are included.

## Requirements

- Julia 1.12.6
- R 4.5.3
- `rvinecopulib` 0.7.3.1.0
- `mvtnorm` 1.3-3

The main Julia environment is defined by `code/Project.toml` and
`code/Manifest.toml`. The research version of `VineCopulas.jl` used in the
experiments is included under `code/vendor/VineCopulas/`.

## Contents

- `code/internal_benchmark/` — principal benchmark and frozen results
- `code/reproduction/` — integration-rule, high-dimensional, financial,
  Gaussian-specific, and figure workflows
- `code/benchmarks/jmva_extension/` — ordering and permutation-entropy scripts
- `code/reports/jmva_extension/` — frozen ordering and ordinal results
- `code/validation/rvinecopulib/` — independent `rvinecopulib` validation
- `data/` — data provenance and download helper
- `output/` — locally regenerated files (ignored by Git)

## Setup

From the repository directory:

```bash
julia --project=code -e 'using Pkg; Pkg.instantiate()'
julia --project=code/vendor/VineCopulas/benchmarks -e 'using Pkg; Pkg.instantiate()'
julia --project=code/reproduction/gaussian_specialized -e 'using Pkg; Pkg.instantiate()'
Rscript R/check_versions.R
```

Exact R-package installation commands are provided in
`R/install_exact_versions.R`.

## Data

The financial and permutation-entropy applications use the ten daily industry
portfolios from the Kenneth R. French Data Library. The source, sample period,
processing steps, and archive checksum are documented in `data/README.md`.

Download and verify the input with:

```bash
./data/download_kenneth_french.sh
```

## Check

Run the portable smoke tests with:

```bash
./check.sh
```

The check verifies the Julia environment and vendored `VineCopulas.jl`, runs a
small probability smoke test, checks the `rvinecopulib` density preflight, and
reconstructs the Gaussian comparison from the frozen raw CSV files. It writes
only to temporary directories and `output/`.

The full experiment scripts are kept in their corresponding directories. Some
benchmark campaigns use many randomized replicates and are intentionally not
part of the quick check.
