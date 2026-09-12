# Phase 4 adaptive high-dimensional scaling

Status: COMPLETE for the adaptive pilot design.

## Mathematical setup

An equicorrelated Gaussian copula with common correlation rho is represented as a Gaussian D-vine whose tree-t pair correlations are partial correlations `rho / (1 + (t-1)rho)`. This was validated numerically by comparing the D-vine log density against a direct Gaussian copula for d=5 and d=10.

- Maximum low-dimensional log-density discrepancy: 1.326583287664107e-11

Reference probabilities use the independent one-dimensional Gaussian-mixture representation, not the reduced estimator.

## Adaptive execution

- Primary dimensions: d=20,50,100.
- Stress dimensions: d=300,500.
- Primary rule: Sobol/Owen.
- Stress tests and high-dimensional Halton/Sobol pilots are checkpointed separately.
- Stress section enabled in this run: false.
- Halton/Sobol high-dimensional pilot enabled in this run: false.
- Projection threshold per method cell: 180.0 seconds.
- Full-dimensional stress indicator skipped by default: true.
  An attempted d=300, N=256, one-rep pilot was manually interrupted after more than two minutes inside the full inverse Rosenblatt path, so the stress section records Reduced results and explicit skipped indicator cells instead of fabricating a matched campaign.

Removing two dimensions is only a small relative dimension reduction at d=300 or d=500; any advantage observed there must come from indicator removal, density-Jacobian cancellation, and integrand regularity, not from dimension reduction alone.

## Outputs

- `highdim_reference_values.csv`
- `highdim_density_parity.csv`
- `highdim_raw.csv`
- `highdim_stress_raw.csv`
- `highdim_halton_sobol.csv`
- `highdim_summary.csv`
- `highdim_skipped_cells.csv` if any cells were capped by projected runtime
