#!/bin/sh
set -eu

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../../.." && pwd)
OUT_DIR=${OUT_DIR:-"$ROOT/output/gaussian_specialized"}
export OUT_DIR
mkdir -p "$OUT_DIR"

julia --project="$HERE" "$HERE/prepare_mvtnorm_inputs.jl"
Rscript "$HERE/run_mvtnorm_comparison.R"
