#!/usr/bin/env bash
# Regenerate every committed result from the raw inputs and check it is identical.
#
#   bash scripts/reproduce.sh            # needs the inputs listed in results/table_s1_inputs.csv
#
# Steps: verify every input against the committed manifest (Table S1: path, size,
# SHA-256; no filesystem dates); test suite; re-interpretation of both cohorts
# under every EUCAST version; every analysis script; figures. Then compares
# results/ with the committed files: exit 0 and "REPRODUCED" when every file is
# byte-identical, otherwise lists the files that differ and exits 1. Exit 2 when
# the inputs do not match the manifest. Inputs are not downloaded here (see README).

set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PYTHON:-python}
export PYTHONPATH=src PYTHONDONTWRITEBYTECODE=1

MANIFEST=results/table_s1_inputs.csv
if [ -n "$(git status --porcelain -- results/)" ]; then
    echo "results/ has uncommitted changes; commit or discard them first" >&2; exit 2
fi

step() { printf '\n== %s\n' "$*"; }
step "inputs against $MANIFEST"
$PY analysis/inputs_manifest.py data --check "$MANIFEST" || exit 2

# every input path below comes from the verified manifest
listed() { awk -F, -v pat="$1" 'NR > 1 && $2 ~ pat {print $2}' "$MANIFEST"; }
by_version() { $PY -c 'import re, sys
key = lambda p: [int(x) for x in re.search(r"v_([0-9.]+)_", p).group(1).split(".")]
print("\n".join(sorted(sys.stdin.read().split(), key=key)))'; }
T=""; for f in $(listed '^data/breakpoints/.*[.]xlsx$' | by_version); do T="$T --table $f"; done
PD=$(listed '^data/ncbi/PDG.*[.]metadata[.]tsv$')
[ "$(printf '%s\n' "$PD" | grep -c .)" = 1 ] || { echo "manifest must list exactly one PDG metadata file" >&2; exit 2; }

step "tests"
$PY -m pytest -q -p no:cacheprovider

LINKS=mappings/bvbrc_to_eucast_enterobacterales.csv
GRP="--group Enterobacterales --group Enterobacteriaceae"
GENO="--pd-metadata $PD --genes data/ncbi/ReferenceGeneCatalog.txt --mutations data/ncbi/AMRProt-mutation.tsv
      --map mappings/amrfinder_subclass_to_antibiotic.csv --exclude mappings/amrfinder_excluded_families.csv"

step "BV-BRC cohort"
R=data/derived/kpneumoniae_reinterpreted.csv
$PY -m bpdrift reinterpret --bvbrc data/bvbrc/kpneumoniae_amr_all.csv --links $LINKS $GRP $T -o $R
$PY analysis/flip_rates.py $R -o results/flip_rates.csv
$PY analysis/as_reported.py $R -o results/as_reported.csv
$PY analysis/robustness.py $R -o results/robustness.csv --transitions results/transitions.csv
$PY analysis/by_source.py $R --bvbrc data/bvbrc/kpneumoniae_amr_all.csv -o results/by_source.csv
G="--genomes data/bvbrc/kpneumoniae_genomes.csv $GENO"
$PY analysis/genotype_errors.py $R $G --matched -o results/genotype_errors.csv \
    --summary results/genotype_join.txt --pairs-out data/derived/genotype_pairs.csv
$PY analysis/genotype_errors.py $R $G -o results/genotype_errors_unmatched.csv \
    --summary results/genotype_join_unmatched.txt
P=data/derived/genotype_pairs.csv
$PY analysis/genotype_diff.py $P -o results/genotype_diff.csv
$PY analysis/genotype_diff.py $P --i-as drop -o results/genotype_diff_s_r.csv
$PY analysis/genotype_diff.py $P --i-as r -o results/genotype_diff_i_as_r.csv
$PY analysis/genotype_diff.py $P --family all -o results/genotype_diff_family_all.csv

step "NCBI replication cohort"
$PY analysis/ncbi_ast.py data/ncbi/ast \
    --exclude-bvbrc data/bvbrc/kpneumoniae_amr_all.csv data/bvbrc/kpneumoniae_genomes.csv \
    -o data/derived/ncbi_ast_bvbrc_format.csv --genomes-out data/derived/ncbi_ast_genomes.csv \
    --summary results/ncbi/cohort.txt
R=data/derived/ncbi_reinterpreted.csv
$PY -m bpdrift reinterpret --bvbrc data/derived/ncbi_ast_bvbrc_format.csv --links $LINKS $GRP $T -o $R
$PY analysis/flip_rates.py $R -o results/ncbi/flip_rates.csv
$PY analysis/as_reported.py $R -o results/ncbi/as_reported.csv
$PY analysis/robustness.py $R -o results/ncbi/robustness.csv --transitions results/ncbi/transitions.csv
G="--genomes data/derived/ncbi_ast_genomes.csv $GENO"
$PY analysis/genotype_errors.py $R $G --matched -o results/ncbi/genotype_errors.csv \
    --summary results/ncbi/genotype_join.txt --pairs-out data/derived/ncbi_genotype_pairs.csv
$PY analysis/genotype_errors.py $R $G -o results/ncbi/genotype_errors_unmatched.csv \
    --summary results/ncbi/genotype_join_unmatched.txt
P=data/derived/ncbi_genotype_pairs.csv
$PY analysis/genotype_diff.py $P -o results/ncbi/genotype_diff.csv
$PY analysis/genotype_diff.py $P --i-as drop -o results/ncbi/genotype_diff_s_r.csv
$PY analysis/genotype_diff.py $P --i-as r -o results/ncbi/genotype_diff_i_as_r.csv
$PY analysis/genotype_diff.py $P --family all -o results/ncbi/genotype_diff_family_all.csv

step "figures"
$PY analysis/figures.py results -o results/figures

step "comparison with the committed results"
CHANGED=$(git status --porcelain -- results/)
N=$(git ls-files results/ | wc -l | tr -d ' ')
if [ -z "$CHANGED" ]; then
    echo "REPRODUCED: all $N committed result files are byte-identical ($(git rev-parse --short HEAD))"
else
    echo "DIFFERENT: these result files changed:"; echo "$CHANGED"; exit 1
fi
