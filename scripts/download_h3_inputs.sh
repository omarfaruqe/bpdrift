#!/usr/bin/env bash
# Download the inputs for the genotype arm (H3). Run from the repository root:
#     bash scripts/download_h3_inputs.sh
# Everything lands in data/ (git-ignored). Re-running overwrites with the latest files,
# so record the printed checksums with the analysis.
set -euo pipefail

mkdir -p data/ncbi data/bvbrc
NCBI=https://ftp.ncbi.nlm.nih.gov/pathogen

echo "1/4 NCBI Pathogen Detection metadata for Klebsiella (AMR_genotypes per isolate)"
META=$(curl -fsSL "$NCBI/Results/Klebsiella/latest_snps/Metadata/" \
       | grep -oE 'PDG[0-9]+\.[0-9]+\.metadata\.tsv' | sort -u | tail -1)
echo "    $META"
curl -fL --retry 3 -o "data/ncbi/$META" "$NCBI/Results/Klebsiella/latest_snps/Metadata/$META"

echo "2/4 AMRFinderPlus reference gene catalogue (gene -> class / subclass)"
DB="$NCBI/Antimicrobial_resistance/AMRFinderPlus/database/latest"
curl -fL --retry 3 -o data/ncbi/ReferenceGeneCatalog.txt "$DB/ReferenceGeneCatalog.txt"
curl -fL --retry 3 -o data/ncbi/amrfinder_db_version.txt "$DB/version.txt"

echo "3/4 AMRFinderPlus point-mutation catalogue"
curl -fL --retry 3 -o data/ncbi/AMRProt-mutation.tsv "$DB/AMRProt-mutation.tsv"

echo "4/4 BV-BRC K. pneumoniae genomes -> BioSample / assembly accessions"
Q='eq(taxon_id,573)&select(genome_id,genome_name,biosample_accession,assembly_accession)'
OUT=data/bvbrc/kpneumoniae_genomes.csv
off=0
while :; do
  curl -fsSL -o page.csv "https://www.bv-brc.org/api/genome/?${Q}&sort(+genome_id)&limit(25000,${off})&http_accept=text/csv"
  n=$(( $(wc -l < page.csv) - 1 ))
  echo "    offset $off: $n rows"
  [ "$off" -eq 0 ] && head -1 page.csv > "$OUT"
  tail -n +2 page.csv >> "$OUT"
  [ "$n" -lt 25000 ] && break
  off=$((off + 25000))
done
rm -f page.csv

echo
echo "AMRFinderPlus database: $(cat data/ncbi/amrfinder_db_version.txt)"
shasum -a 256 data/ncbi/* "$OUT"
wc -l "data/ncbi/$META" "$OUT"
