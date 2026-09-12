# Publication benchmark report

## Design

- Mode: `final`
- Cases: 20
- Budgets: 256, 1024, 4096, 16384, 65536
- Replicates per cell: 30
- References meeting the aspirational absolute target: 15/20
- Aggregate-cell reference rule: split-half z <= 4.0 and uncertainty/RMSE <= 0.2
- The absolute reference target is reported as an audit field; it is not itself required for a method cell to be resolved.

## Equal-budget paired results

- Valid paired replicate comparisons: 3000
- Reduced error wins: 2902/3000 (96.73333333333333%)
- Reduced runtime wins: 1902/3000 (63.4%)
- Reduced wins both: 1834/3000 (61.13333333333333%)
- Median reduced/current absolute-error ratio: 0.009534882064673462

## Matched accuracy

Matched accuracy is evaluated for EVERY resolved current-method target cell; it does not select only the maximum speedup per case.

- Matched target cells: 100
- Median speedup over all matched targets: 14.311198909574312
- 10th percentile speedup: 1.0576002969352312
- 90th percentile speedup: 297.7597000123735
- Median speedup at the strictest resolved current target per case: 308.87522163501194

## Complexity check

For the compiled D-vine traversal, exact counts with s=d-2 are checked:
CDF: hfunc=s(3s+1)/2, hinv=s(s-1)/2, total=2s^2.
Interior rectangle: hfunc=2s(s+1), hinv=s(s-1)/2, total=s(5s+3)/2.
Genuine R-vine counts are reported as structure-dependent; no D-vine exact formula is imposed on them.
For genuine R-vines, operation counts are recorded empirically; a universal exact complexity formula is not claimed here without a separate proof.

## Reference-resolution audit

- Aggregate method cells: 200
- Cells meeting reference-resolution rule: 200
- Cells flagged unresolved: 0

Raw replicate data, references, bootstrap summaries, Pareto frontiers, and figure/table-ready CSV files are retained under this run directory.
