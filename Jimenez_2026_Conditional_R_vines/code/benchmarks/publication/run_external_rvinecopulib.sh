#!/usr/bin/env bash
set -euo pipefail

: "${MODE:=final}"
: "${CASE_SET:=core}"
: "${OUT:=benchmarks/results/publication/final_01}"
: "${POWERS:=8,10,12,14,16}"
: "${EXTERNAL_REPS:=30}"
: "${EXTERNAL_CORES:=1}"
: "${EXTERNAL_RESUME:=true}"

export MODE CASE_SET OUT POWERS EXTERNAL_REPS EXTERNAL_CORES EXTERNAL_RESUME
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1

echo "External rvinecopulib publication campaign"
echo "  MODE=${MODE}"
echo "  CASE_SET=${CASE_SET}"
echo "  OUT=${OUT}"
echo "  POWERS=${POWERS}"
echo "  EXTERNAL_REPS=${EXTERNAL_REPS}"
echo "  EXTERNAL_CORES=${EXTERNAL_CORES}"

command -v Rscript >/dev/null 2>&1 || { echo "Rscript is not available" >&2; exit 1; }
Rscript -e 'if (!requireNamespace("rvinecopulib", quietly=TRUE)) stop("Install rvinecopulib with install.packages(\"rvinecopulib\")"); cat("rvinecopulib ", as.character(packageVersion("rvinecopulib")), "\n", sep="")'

julia -t 1 --project=benchmarks benchmarks/publication/prepare_rvinecopulib_external.jl
Rscript benchmarks/publication/rvinecopulib_preflight.R
Rscript benchmarks/publication/rvinecopulib_external.R
julia -t 1 --project=benchmarks benchmarks/publication/analyze_external_rvinecopulib.jl

echo "External campaign complete: ${OUT}/external_rvinecopulib"
