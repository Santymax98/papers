#!/usr/bin/env bash
set -euo pipefail

MODE="${MODE:-pilot}"
CASE_SET="${CASE_SET:-core}"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
OUT="${OUT:-benchmarks/results/publication/${RUN_ID}}"

export MODE CASE_SET OUT

if [[ "$MODE" == "final" ]]; then
  export DIMS="${DIMS:-4,5,10,20}"
  export POWERS="${POWERS:-8,10,12,14,16}"
  export REPS="${REPS:-30}"
  export REFERENCE_R="${REFERENCE_R:-32}"
  export REFERENCE_POWERS="${REFERENCE_POWERS:-16,18,20}"
  export REFERENCE_TARGET_SE="${REFERENCE_TARGET_SE:-1e-7}"
  export BOOTSTRAP_B="${BOOTSTRAP_B:-2000}"
else
  export DIMS="${DIMS:-4,5,10}"
  export POWERS="${POWERS:-8,12,16}"
  export REPS="${REPS:-4}"
  export REFERENCE_R="${REFERENCE_R:-8}"
  export REFERENCE_POWERS="${REFERENCE_POWERS:-14,16}"
  export REFERENCE_TARGET_SE="${REFERENCE_TARGET_SE:-5e-6}"
  export BOOTSTRAP_B="${BOOTSTRAP_B:-300}"
fi

mkdir -p "$OUT"
printf 'Publication run\n  MODE=%s\n  CASE_SET=%s\n  OUT=%s\n' "$MODE" "$CASE_SET" "$OUT"

julia -t 1 --project=benchmarks benchmarks/publication/preflight.jl
julia -t 1 --project=benchmarks benchmarks/publication/generate_references.jl
julia -t 1 --project=benchmarks benchmarks/publication/run_publication_benchmark.jl
julia -t 1 --project=benchmarks benchmarks/publication/analyze_publication_results.jl

echo "Publication benchmark complete: $OUT"
