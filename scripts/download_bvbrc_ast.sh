#!/usr/bin/env bash
# Download BV-BRC laboratory AST records for K. pneumoniae (taxon 573). Run from the
# repository root:  bash scripts/download_bvbrc_ast.sh
# Pages of 25,000 with a fixed sort order, so no record is skipped or repeated.
set -euo pipefail
mkdir -p data/bvbrc
OUT=data/bvbrc/kpneumoniae_amr_all.csv
Q='and(eq(taxon_id,573),eq(evidence,Laboratory%20Method))'
off=0
while :; do
  curl -fsSL -o page.csv "https://www.bv-brc.org/api/genome_amr/?${Q}&sort(+id)&limit(25000,${off})&http_accept=text/csv"
  n=$(( $(wc -l < page.csv) - 1 ))
  echo "offset $off: $n rows"
  [ "$off" -eq 0 ] && head -1 page.csv > "$OUT"
  tail -n +2 page.csv >> "$OUT"
  [ "$n" -lt 25000 ] && break
  off=$((off + 25000))
done
rm -f page.csv
wc -l "$OUT"
shasum -a 256 "$OUT"
