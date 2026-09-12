# Gaussian-specific comparison

This directory contains the four-case Gaussian comparison reported with the
article. It compares Reduced RQMC, the full-dimensional Indicator estimator,
`AdditionalDistributions.jl` 0.3.0, and `mvtnorm::GenzBretz`.

The retained CSV files are the numerical results used in the reported table.
The experiment uses 30 independent randomizations and nominal budget or
`maxpts=4096`.

Restore the Julia environment with:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

To rerun the Julia comparison without overwriting the frozen files:

```bash
OUT_DIR=/path/to/output julia --project=. run_gaussian_specialized_comparison.jl
```

To rerun the `mvtnorm` comparison:

```bash
OUT_DIR=/path/to/output ./run_mvtnorm_comparison.sh
```

To reconstruct the reported table from the frozen raw CSV files:

```bash
OUT_TABLE=/path/to/table.tex julia --project=. rebuild_table5_from_raw.jl
```
