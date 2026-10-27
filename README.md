# bpdrift

[![tests](https://github.com/omarfaruqe/bpdrift/actions/workflows/tests.yml/badge.svg)](https://github.com/omarfaruqe/bpdrift/actions/workflows/tests.yml)

Interpret minimum inhibitory concentration (MIC) values as S / I / R under a
specific, pinned version of the EUCAST clinical breakpoint tables, and report
where two versions disagree.

## Design principles

- **No breakpoint data is bundled.** EUCAST tables are copyright EUCAST. You
  download the official file yourself; bpdrift parses it and pins it by
  version string, retrieval date and SHA-256.
- **Censored MICs are handled explicitly.** A value such as `<=2` or `>8` that
  could fall on either side of a breakpoint returns `UNDECIDABLE`, never a
  guess.
- **Fail loudly on ambiguity.** Unparseable input raises an error.

## Install

```bash
pip install -e ".[dev]"          # library, CLI and tests
pip install -e ".[dev,analysis]" # plus matplotlib for the paper figures
pytest -q
```

## Usage

Download the EUCAST clinical breakpoint tables (.xlsx, "file for screen")
from eucast.org into `data/breakpoints/`, a folder git ignores.

```bash
# record exactly which file you used
bpdrift pin data/breakpoints/v_16.1_Breakpoint_Tables.xlsx \
    --version 16.1 --retrieved 2026-10-26 -o data/breakpoints/v16.1.pin.json

# interpret one MIC
bpdrift interpret "<=4" --table data/breakpoints/v_16.1_Breakpoint_Tables.xlsx \
    --group Enterobacterales --agent "Meropenem (indications other than meningitis)"

# compare two versions; one CSV row per breakpoint, with the dilutions that change category
bpdrift diff data/breakpoints/v_10.0_Breakpoint_Tables.xlsx \
    data/breakpoints/v_16.1_Breakpoint_Tables.xlsx \
    --group Enterobacterales --changed-only -o diff.csv
```

Rows are matched across versions by organism group and agent name only.
Rows EUCAST renamed or split are reported as removed and added, with
likely partners listed as candidates. To pair them, list them in a CSV
passed with `--map`:

```text
old_group,old_agent,new_group,new_agent
Enterobacterales,Ciprofloxacin,Enterobacterales,Ciprofloxacin (meningitis)
Enterobacterales,Ciprofloxacin,Enterobacterales,Ciprofloxacin (indications other than meningitis)
P.multocida,,Pasteurella,
```

A line with no agents renames a whole organism group. A mapping line
naming a row that does not exist is an error. An optional fifth column,
`note`, records why two rows correspond; lines starting with `#` are
comments.

## Re-interpreting public AST data

`bpdrift reinterpret` reads a BV-BRC `genome_amr` CSV export, links each
antibiotic name to the EUCAST row used in every table version through a
reviewed CSV (`mappings/bvbrc_to_eucast_enterobacterales.csv`), and writes
one row per MIC per version. Records it cannot use are counted, never
dropped silently.

```bash
bpdrift reinterpret --bvbrc data/bvbrc/kpneumoniae_amr_all.csv \
    --links mappings/bvbrc_to_eucast_enterobacterales.csv \
    --group Enterobacterales --group Enterobacteriaceae \
    --table data/breakpoints/v_8.1_Breakpoint_Tables.xlsx \
    --table data/breakpoints/v_16.1_Breakpoint_Tables.xlsx \
    -o data/derived/kpneumoniae_reinterpreted.csv
```

## Genotype predictor error rates

`bpdrift.genotype` scores a fixed AMRFinderPlus-based predictor against
labels from each breakpoint version (very major and major error rates with
Wilson 95% intervals). The drug-class mapping and the excluded intrinsic
genes are reviewable CSVs in `mappings/`.

## Reproducing the K. pneumoniae analysis

Inputs are downloaded, not committed: EUCAST tables by hand into
`data/breakpoints/`, the rest with the two download steps below.

```bash
bash scripts/download_bvbrc_ast.sh    # BV-BRC laboratory AST records, K. pneumoniae
bash scripts/download_h3_inputs.sh    # NCBI genotypes, AMRFinderPlus catalogues, BV-BRC accessions

B=data/breakpoints; T=""; for v in 8.1 9.0 10.0 13.1 14.0 15.0 16.1; do T="$T --table $B/v_${v}_Breakpoint_Tables.xlsx"; done
bpdrift reinterpret --bvbrc data/bvbrc/kpneumoniae_amr_all.csv \
    --links mappings/bvbrc_to_eucast_enterobacterales.csv \
    --group Enterobacterales --group Enterobacteriaceae $T \
    -o data/derived/kpneumoniae_reinterpreted.csv
python analysis/flip_rates.py data/derived/kpneumoniae_reinterpreted.csv -o results/flip_rates.csv
python analysis/as_reported.py data/derived/kpneumoniae_reinterpreted.csv -o results/as_reported.csv
python analysis/robustness.py data/derived/kpneumoniae_reinterpreted.csv \
    -o results/robustness.csv --transitions results/transitions.csv
G="--genomes data/bvbrc/kpneumoniae_genomes.csv --pd-metadata data/ncbi/PDG*.metadata.tsv \
   --genes data/ncbi/ReferenceGeneCatalog.txt --mutations data/ncbi/AMRProt-mutation.tsv \
   --map mappings/amrfinder_subclass_to_antibiotic.csv \
   --exclude mappings/amrfinder_excluded_families.csv"
# primary: identical genome-drug pairs under every version
python analysis/genotype_errors.py data/derived/kpneumoniae_reinterpreted.csv $G --matched \
    -o results/genotype_errors.csv --summary results/genotype_join.txt
# sensitivity: every pair labelled under each version
python analysis/genotype_errors.py data/derived/kpneumoniae_reinterpreted.csv $G \
    -o results/genotype_errors_unmatched.csv --summary results/genotype_join_unmatched.txt
python analysis/figures.py results -o results/figures
```

## Status

Version 0.1.0. A research tool, not for clinical use: it reports what a
breakpoint table says about an MIC, not what to treat with.

## Citation

See `CITATION.cff`, or use GitHub's "Cite this repository" button.

## Licence

Code: MIT. Breakpoint tables remain the property of their publishers.
