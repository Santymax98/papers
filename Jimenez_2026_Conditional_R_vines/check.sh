#!/usr/bin/env bash
set -euo pipefail

ROOT="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$ROOT/output"

echo "Checking Julia environment and vendored VineCopulas..."
REPRO_ROOT="$ROOT" julia --project="$ROOT/code" -e '
    using Pkg
    Pkg.instantiate()
    using VineCopulas
    expected = realpath(joinpath(ENV["REPRO_ROOT"], "code", "vendor", "VineCopulas"))
    loaded = realpath(dirname(dirname(pathof(VineCopulas))))
    loaded == expected || error("VineCopulas loaded from $loaded instead of $expected")
'

julia --project="$ROOT/code" "$ROOT/code/run_tutorial.jl" > "$TMP/tutorial.txt"
grep -q "Smoke test completed successfully" "$TMP/tutorial.txt"

echo "Checking Gaussian-specific environment..."
julia --project="$ROOT/code/reproduction/gaussian_specialized" -e '
    using Pkg
    Pkg.instantiate()
    using AdditionalDistributions
'

echo "Checking R environment and rvinecopulib density parity..."
Rscript "$ROOT/R/check_versions.R" > "$TMP/R_versions.txt"
cp -R "$ROOT/code/internal_benchmark/frozen_results" "$TMP/rvine_internal"
OUT="$TMP/rvine_internal" julia -t 1 \
  --project="$ROOT/code/vendor/VineCopulas/benchmarks" \
  "$ROOT/code/validation/rvinecopulib/prepare_rvinecopulib_external.jl" \
  > "$TMP/rvinecopulib_prepare.txt"
OUT="$TMP/rvine_internal" Rscript \
  "$ROOT/code/validation/rvinecopulib/rvinecopulib_preflight.R" \
  > "$TMP/rvinecopulib_preflight.txt"
grep -q "preflight passed" "$TMP/rvinecopulib_preflight.txt"

echo "Checking frozen Gaussian comparison..."
OUT_TABLE="$TMP/gaussian_specialized_comparison.tex" julia \
  --project="$ROOT/code/reproduction/gaussian_specialized" \
  "$ROOT/code/reproduction/gaussian_specialized/rebuild_table5_from_raw.jl" \
  > "$TMP/gaussian_table.txt"
grep -q "2.501e-06" "$TMP/gaussian_specialized_comparison.tex"
grep -q "3.792e-04" "$TMP/gaussian_specialized_comparison.tex"

echo "PASS"
