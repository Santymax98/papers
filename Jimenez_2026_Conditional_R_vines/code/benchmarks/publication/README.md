# Publication benchmark for reduced R-vine probabilities

This directory contains the reproducible benchmark battery intended for the
paper on numerical computation of probabilities for simplified absolutely
continuous R-vines.

It deliberately does **not** change the public `cdf` API.

## What is recorded

Every run has its own output directory and retains enough information to
reconstruct tables/figures without rerunning the expensive integrations:

- `run_metadata.toml`: Git commit/branch, Julia version, machine, threads,
  seeds, RQMC budgets, reference tolerances, and experiment mode.
- `case_manifest.csv`: case definitions, dimensions, lower/upper bounds,
  terminal pair type, and default sampling order.
- `references.csv`: final reference estimate and its uncertainty for every
  case.
- `references/<case>.toml`: per-case reference checkpoint, so reference
  generation is resumable.
- `references/<case>_history.csv`: adaptive high-dimensional reference history.
- `references/<case>_replicates.csv`: independent reference scramble estimates.
- `raw/<case>.csv`: replicate-level benchmark data written incrementally.
- `raw_results.csv`: consolidated replicate-level data.
- `plan_costs.csv`: algorithmic plan-construction costs, recorded separately
  from steady integration.
- `summary_by_case.csv`: RMSE, bootstrap CIs, runtime CIs, allocations and
  reference-resolution diagnostics.
- `paired_equal_N.csv`: current/reduced paired comparisons at equal `N`.
- `matched_accuracy.csv`: honest matched-accuracy comparisons. For every
  current-method target cell it chooses the fastest reduced cell with no larger
  RMSE; it does not cherry-pick only the largest speedup per case.
- `matched_accuracy_strictest.csv`: matched accuracy at the strictest resolved
  current target for each case.
- `empirical_convergence_rates.csv`: descriptive log-RMSE/log-N slopes. These
  are empirical diagnostics, not theoretical RQMC rate claims.
- `complexity_check.csv`: operation-count identities checked against the code.
- `pareto_frontier.csv`: error/runtime Pareto points.
- `figure_data/*.csv`: figure-ready data.
- `table_data/*.csv`: table-ready data.
- `REPORT.md`: compact human-readable summary.
- `environment/`: copies of benchmark/package project metadata.

Optional diagnostics produce:

- `regularity_raw.csv` and `regularity_by_stage.csv`;
- `sampling_order_study.csv`;
- `external_rvinecopulib/`: external `rvinecopulib` raw data, model-parity audits, three-method summaries, and figure/table-ready CSVs.

## Mathematical operation counts

Let `s = d - 2` be the conditioning dimension. The following exact identities
are for the compiled **D-vine traversal** used in the scaling benchmark.
Genuine R-vines can reuse conditional states differently and can require fewer
calls; their exact counts are therefore reported empirically rather than forced
into the D-vine formula.

For a D-vine CDF (`lower = 0`):

```text
hfunc / point = s(3s + 1)/2
hinv  / point = s(s - 1)/2
total h+hinv  = 2s^2 = 2(d-2)^2
terminal pair CDF calls = 1
```

For the fully interior D-vine rectangular cases in this battery:

```text
hfunc / point = 2s(s + 1)
hinv  / point = s(s - 1)/2
total h+hinv  = s(5s + 3)/2
terminal pair CDF calls = 4
```

Both D-vine counts are quadratic in dimension. For genuine R-vines the
battery records the actual structure-dependent counts. A universal exact
complexity formula for the generic R-vine traversal should not be claimed until
it has a separate proof. Do not infer algorithmic complexity from a small-sample
log-log regression against `d`.

## Reference strategy

Low dimensions (`d <= 5`) use deterministic HCubature on the reduced
`d-2`-dimensional integrand.

Higher dimensions use an **independent** reduced-RQMC reference run with:

- independent reference seed range;
- multiple Owen-scrambled Sobol replicates;
- adaptive `N` over `REFERENCE_POWERS`;
- an explicit target standard error;
- a split-half consistency diagnostic.

The analysis additionally flags a benchmark cell as reference-resolved only if

```text
reference uncertainty / observed RMSE <= MAX_REFERENCE_RMSE_RATIO
```

(default `0.20`). This prevents apparent high-precision method differences from
being reported below the resolution of the numerical reference.

## Pilot run

Run this first. It checks the machinery without committing to the expensive
publication experiment.

```bash
cd /path/to/VineCopulas
MODE=pilot \
OUT=benchmarks/results/publication/pilot_01 \
bash benchmarks/publication/run_all.sh
```

Defaults in pilot mode are intentionally small.

## Final run

The recommended starting configuration is:

```bash
cd /path/to/VineCopulas
MODE=final \
CASE_SET=core \
DIMS=4,5,10,20 \
POWERS=8,10,12,14,16 \
REPS=30 \
REFERENCE_R=32 \
REFERENCE_POWERS=16,18,20 \
REFERENCE_TARGET_SE=1e-7 \
BOOTSTRAP_B=2000 \
OUT=benchmarks/results/publication/final_01 \
bash benchmarks/publication/run_all.sh
```

Run with **one Julia thread** for timing comparability. The scripts record the
thread count and warn otherwise.

The final reference stage can be expensive in dimensions 10 and 20. It is
resumable: completed target-met references are stored per case and reused when
`RESUME=true` (the default).

## Run stages separately

```bash
OUT=benchmarks/results/publication/final_01 \
MODE=final \
julia -t 1 --project=benchmarks benchmarks/publication/preflight.jl

OUT=benchmarks/results/publication/final_01 \
MODE=final \
julia -t 1 --project=benchmarks benchmarks/publication/generate_references.jl

OUT=benchmarks/results/publication/final_01 \
MODE=final \
julia -t 1 --project=benchmarks benchmarks/publication/run_publication_benchmark.jl

OUT=benchmarks/results/publication/final_01 \
MODE=final \
julia -t 1 --project=benchmarks benchmarks/publication/analyze_publication_results.jl
```

## Optional regularity diagnostic

This is a numerical diagnostic of the diagonal transport stretching

```text
kappa_j = |partial v_{pi_j} / partial z_j|
```

estimated by centered finite differences. It is **not** a proof of mixed
Sobolev regularity or an RQMC convergence theorem.

```bash
OUT=benchmarks/results/publication/final_01 \
MODE=final \
REGULARITY_SAMPLES=512 \
julia -t 1 --project=benchmarks benchmarks/publication/run_regularity_diagnostics.jl
```

## Optional sampling-order study

Only cases for which `:smallest_width_first` actually changes the admissible
order are retained.

```bash
OUT=benchmarks/results/publication/final_01 \
MODE=final \
ORDER_N=1024 \
ORDER_REPS=30 \
julia -t 1 --project=benchmarks benchmarks/publication/run_order_study.jl
```

This is diagnostic. `smallest_width_first` is bound-motivated, not claimed to
be RQMC-optimal.

## Before a final paper run

1. Freeze the code commit/tag locally.
2. Ensure `git status` is clean.
3. Keep the Julia version fixed for the complete run.
4. Run `Pkg.test()` and `git diff --check`.
5. Run the pilot battery.
6. Run final references.
7. Inspect `references.csv` and ensure target references are resolved.
8. Run the final benchmark.
9. Run analysis once; figures and tables should subsequently be generated from
   retained CSVs, not by rerunning the estimators.
10. Preserve the complete output directory unchanged as the paper artifact.


## Reference-resolution policy

`REFERENCE_TARGET_SE` is an aspirational stopping target for reference generation.
A difficult reference that does not reach this absolute target is not automatically
discarded. In the publication analysis, a method cell is considered resolved when
(i) the independent-reference split-half diagnostic is stable and (ii) the reference
uncertainty is at most `MAX_REFERENCE_RMSE_RATIO` times the observed cell RMSE.
This makes the adequacy criterion relative to the numerical error actually being
measured while retaining the absolute target as an audit field.

## Final external comparison with `rvinecopulib`

Run this only after the internal final benchmark and `analyze_publication_results.jl`
have been frozen.  It does **not** modify or regenerate the publication references.

The external comparison is intentionally CDF-only.  It includes all publication
CDF cases, including the genuine 5D R-vine.  General rectangles are not converted
into a `2^d` inclusion-exclusion workload because that would be an artificial and
unfair use of `pvinecop()`.

Before any timing, Julia exports two deterministic interior log-density values per
case and R reconstructs the same vine.  `rvinecopulib_preflight.R` aborts unless
all model-parity checks pass.  It also verifies same-seed reproducibility of the
randomized `pvinecop()` path.

The final design uses the same budgets `2^8, ..., 2^16`, 30 independent external
replicates per cell, and `cores=1`.  R and Julia randomized designs are not forced
to share point sets, so cross-language results are compared at the aggregate-cell
level (RMSE and median wall time), never as falsely paired replicates.

Run:

```bash
MODE=final \
CASE_SET=core \
OUT=benchmarks/results/publication/final_01 \
caffeinate -dimsu bash benchmarks/publication/run_external_rvinecopulib.sh
```

The output is isolated under:

```text
benchmarks/results/publication/final_01/external_rvinecopulib/
```

Important outputs include `raw_results.csv`, `model_parity_audit.csv`,
`summary.csv`, `three_method_summary.csv`, `equal_budget.csv`,
`matched_accuracy.csv`, `matched_accuracy_strictest.csv`, `pareto_frontier.csv`,
`figure_data/`, `table_data/`, and `REPORT.md`.
