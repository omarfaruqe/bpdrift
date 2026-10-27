#!/usr/bin/env bash
# Replication cohort: K. pneumoniae BioSamples carrying an antibiogram (AST table),
# fetched from NCBI BioSample through E-utilities. Raw XML is saved in batches;
# analysis/ncbi_ast.py turns it into BV-BRC-style records.
#
#   bash scripts/download_ncbi_ast.sh            # optional: NCBI_API_KEY=... for faster requests
#
# Output: data/ncbi/ast/biosample_<start>.xml, data/ncbi/ast/query.txt (query, date, count)

set -euo pipefail
OUT=data/ncbi/ast
mkdir -p "$OUT"
rm -f "$OUT"/biosample_*.xml     # one consistent snapshot per run
EU=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
TERM='antibiogram[filter] AND "Klebsiella pneumoniae"[orgn]'
KEY=${NCBI_API_KEY:+&api_key=$NCBI_API_KEY}
PAUSE=${NCBI_API_KEY:+0.12}; PAUSE=${PAUSE:-0.4}
BATCH=500

ENC=$(python3 -c 'import sys, urllib.parse; print(urllib.parse.quote(sys.argv[1]))' "$TERM")
SEARCH=$(curl -sS --retry 3 "$EU/esearch.fcgi?db=biosample&term=$ENC&usehistory=y&retmax=0$KEY")
COUNT=$(printf '%s' "$SEARCH" | sed -n 's:.*<Count>\([0-9]*\)</Count>.*:\1:p' | head -1)
WEBENV=$(printf '%s' "$SEARCH" | sed -n 's:.*<WebEnv>\(.*\)</WebEnv>.*:\1:p')
QKEY=$(printf '%s' "$SEARCH" | sed -n 's:.*<QueryKey>\([0-9]*\)</QueryKey>.*:\1:p')
if [ -z "$COUNT" ] || [ "$COUNT" = 0 ]; then
    echo "no BioSamples matched; esearch said:" >&2; printf '%s\n' "$SEARCH" >&2; exit 1
fi
printf 'query\t%s\ndate\t%s\ncount\t%s\n' "$TERM" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$COUNT" > "$OUT/query.txt"
echo "BioSamples with an antibiogram: $COUNT"

for (( start=0; start<COUNT; start+=BATCH )); do
    f=$(printf '%s/biosample_%06d.xml' "$OUT" "$start")
    curl -sS --retry 5 --retry-delay 3 -o "$f" \
        "$EU/efetch.fcgi?db=biosample&query_key=$QKEY&WebEnv=$WEBENV&retstart=$start&retmax=$BATCH&retmode=xml$KEY"
    tail -c 200 "$f" | grep -q '</BioSampleSet>' || { echo "incomplete batch $f; rerun the script" >&2; exit 1; }
    printf '\r%d / %d' "$(( start + BATCH < COUNT ? start + BATCH : COUNT ))" "$COUNT"
    sleep "$PAUSE"
done
echo
ls "$OUT"/biosample_*.xml | wc -l | xargs echo "batches:"
shasum -a 256 "$OUT"/biosample_*.xml | shasum -a 256 | cut -c1-16 | xargs echo "checksum of all batches:"
