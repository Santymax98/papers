#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)}"
URL="https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/10_Industry_Portfolios_daily_CSV.zip"
EXPECTED="b1a23edba335ddd37223c46d75a6aa07827170dd6fd04d9a27249602e1054e4f"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

curl --fail --location --silent --show-error "$URL" --output "$TMP/data.zip"
ACTUAL="$(shasum -a 256 "$TMP/data.zip" | awk '{print $1}')"
if [[ "$ACTUAL" != "$EXPECTED" ]]; then
  echo "Checksum mismatch: expected $EXPECTED, received $ACTUAL" >&2
  exit 1
fi

FINANCIAL="$ROOT/code/data/10_Industry_Portfolios_daily_CSV.zip"
ORDINAL="$ROOT/code/benchmarks/jmva_extension/permutation_entropy/data/10_Industry_Portfolios_daily_CSV.zip"
mkdir -p "$(dirname "$FINANCIAL")" "$(dirname "$ORDINAL")"
cp "$TMP/data.zip" "$FINANCIAL"
cp "$TMP/data.zip" "$ORDINAL"
echo "Verified Kenneth French archive installed."

