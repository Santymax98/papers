# External rvinecopulib benchmark report

## Design

- CDF cases: 16
- External method: `rvinecopulib::pvinecop`
- Budgets: 256, 1024, 4096, 16384, 65536
- Replicates per external cell: 30
- Cores: 1
- General rectangles are excluded; no 2^d inclusion-exclusion baseline is manufactured.
- R and Julia models are audited by deterministic log-density equality before timing.
- Cross-language randomized replicates are independent; equal-budget comparisons below are cell-level, not replicate-paired.

## Model parity

- Density audit rows passed: 32/32

## Reference resolution

- External aggregate cells: 80
- Resolved external cells: 80/80

## Equal QMC budget: reduced vs rvinecopulib

- Resolved three-method cells: 80/80
- Reduced lower RMSE: 80/80
- Reduced lower median runtime: 11/80
- Reduced wins both: 11/80
- Median reduced/rvinecopulib RMSE ratio: 0.0108309946070393
- rvinecopulib lower RMSE than current indicator: 12/80

## Matched accuracy: rvinecopulib targets

Every resolved rvinecopulib target cell is retained; the fastest reduced cell with RMSE no larger is selected.

- Target cells: 80
- Targets matched by reduced: 80/80
- Median speedup: 6.215193666437292
- 10th percentile speedup: 0.6254340857301136
- 90th percentile speedup: 93.18349277736435
- Median speedup at the strictest rvinecopulib target per case: 90.50057783157824

Raw external replicates, model-parity audit, bootstrap summaries, matched-accuracy data, Pareto frontiers, and figure/table-ready CSVs are retained under this directory.
