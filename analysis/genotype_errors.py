"""Genotype arm (H3): error rates of a fixed genotype predictor under each EUCAST version.

    python analysis/genotype_errors.py data/derived/kpneumoniae_reinterpreted.csv \
        --genomes data/bvbrc/kpneumoniae_genomes.csv \
        --pd-metadata data/ncbi/PDG*.metadata.tsv \
        --genes data/ncbi/ReferenceGeneCatalog.txt --mutations data/ncbi/AMRProt-mutation.tsv \
        --map mappings/amrfinder_subclass_to_antibiotic.csv \
        --exclude mappings/amrfinder_excluded_families.csv \
        --matched -o results/genotype_errors.csv --summary results/genotype_join.txt

Labels: main-analysis records only, one label per (genome, antibiotic, version)
when that genome's records agree. Genotypes: NCBI Pathogen Detection's
AMRFinderPlus calls (AMR_genotypes), joined to BV-BRC genomes by BioSample.
An isolate with an empty AMR_genotypes field is treated as not analysed and
left out: every K. pneumoniae genome carries at least its chromosomal
beta-lactamase and fosA, so an empty field means no calls were made.

With --matched, only genome-drug pairs labelled S, I or R under every version
are compared, so every version is scored on identical pairs and differences
between versions come from the labels alone.

Output: one row per (antibiotic or ALL, EUCAST version): genome-drug pairs
compared, then very major and
major error counts and rates (Wilson 95% intervals). Only phenotype labels
S, I and R are compared; undecidable labels are counted separately.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from bpdrift.genotype import (
    collapse_labels,
    error_rates,
    load_catalogue,
    Exclusions,
    load_exclusions,
    load_subclass_map,
    parse_amr_genotypes,
    predict,
)

COLUMNS = [
    "antibiotic", "eucast_version", "compared", "undecidable", "phenotype_r", "phenotype_not_r",
    "vme", "vme_rate", "vme_lo", "vme_hi", "me", "me_rate", "me_lo", "me_hi", "agreement_rate",
]


def _vkey(v):
    return tuple(int(x) for x in v.split("."))


def _column(fieldnames, *candidates):
    """First field whose normalised name matches a candidate."""
    canon = {f.lower().replace(" ", "_").lstrip("#"): f for f in fieldnames}
    for c in candidates:
        if c in canon:
            return canon[c]
    raise SystemExit(f"none of {candidates} in columns {list(fieldnames)[:12]}...")


def read_biosamples(path) -> dict[str, str]:
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh)
        gid = _column(r.fieldnames, "genome_id")
        bs = _column(r.fieldnames, "biosample_accession")
        return {row[gid].strip(): row[bs].strip() for row in r if row[bs].strip()}


def read_pd_genotypes(path) -> tuple[dict[str, str], int]:
    """biosample -> AMR_genotypes; biosamples with differing entries are dropped."""
    seen = defaultdict(set)
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        bs = _column(r.fieldnames, "biosample_acc")
        amr = _column(r.fieldnames, "amr_genotypes")
        for row in r:
            b = (row[bs] or "").strip()
            if b and b != "NULL":
                seen[b].add((row[amr] or "").strip())
    out = {b: next(iter(v)) for b, v in seen.items() if len(v) == 1}
    return out, len(seen) - len(out)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("reinterpreted")
    p.add_argument("--genomes", required=True)
    p.add_argument("--pd-metadata", required=True)
    p.add_argument("--genes", required=True)
    p.add_argument("--mutations")
    p.add_argument("--map", required=True)
    p.add_argument("--exclude", help="CSV of gene families or symbols the predictor ignores")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--pairs-out", help="also write one row per compared genome-drug-version pair "
                   "(genome, antibiotic, version, predicted resistant, phenotype) for paired analyses")
    p.add_argument("--matched", action="store_true",
                   help="compare only genome-drug pairs labelled S, I or R under every version")
    a = p.parse_args(argv)

    with open(a.reinterpreted, newline="", encoding="utf-8") as fh:
        rows = [(r["genome_id"], r["antibiotic"], r["eucast_version"], r["category"])
                for r in csv.DictReader(fh) if r["exclusion"] == ""]
    labels, label_conflicts = collapse_labels(rows)
    all_versions = {v for _, _, v in labels}
    unmatched = 0
    if a.matched:
        decided = defaultdict(set)
        for (g, ab, v), cat in labels.items():
            if cat in ("S", "I", "R"):
                decided[(g, ab)].add(v)
        keep = {k for k, vs in decided.items() if vs == all_versions}
        before = len({(g, ab) for g, ab, _ in labels})
        labels = {k: c for k, c in labels.items() if (k[0], k[1]) in keep}
        unmatched = before - len(keep)
    genomes = {g for g, _, _ in labels}

    biosample = read_biosamples(a.genomes)
    pd, pd_conflicts = read_pd_genotypes(a.pd_metadata)
    catalogue = load_catalogue(a.genes, a.mutations)
    smap = load_subclass_map(a.map)
    excluded = load_exclusions(a.exclude) if a.exclude else Exclusions()

    predictions, unknown, ignored = {}, Counter(), Counter()
    lost = Counter()
    for g in genomes:
        b = biosample.get(g)
        if not b:
            lost["no BioSample in BV-BRC"] += 1
            continue
        if b not in pd:
            lost["BioSample not in NCBI Pathogen Detection"] += 1
            continue
        calls = parse_amr_genotypes(pd[b])
        if not calls:
            lost["empty AMR_genotypes"] += 1
            continue
        pr = predict(calls, catalogue, smap, excluded_families=excluded.families,
                     excluded_symbols=excluded.symbols)
        unknown.update(pr.unknown)
        ignored.update(x.split("=")[-1] if "=" in x else catalogue.family_of(x) for x in pr.ignored)
        predictions[g] = pr.resistant

    versions = sorted({v for _, _, v in labels}, key=_vkey)
    targets = smap.targets
    pairs = defaultdict(list)          # (antibiotic, version) -> [(predicted, category)]
    undecidable = Counter()
    for (g, ab, v), cat in labels.items():
        if g not in predictions or ab not in targets:
            continue
        if cat in ("S", "I", "R"):
            pairs[(ab, v)].append((ab in predictions[g], cat))
            pairs[("ALL", v)].append((ab in predictions[g], cat))
        else:
            undecidable[(ab, v)] += 1
            undecidable[("ALL", v)] += 1

    if a.pairs_out:
        Path(a.pairs_out).parent.mkdir(parents=True, exist_ok=True)
        with open(a.pairs_out, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["genome_id", "antibiotic", "eucast_version", "predicted_resistant", "phenotype"])
            for (g, ab, v), cat in sorted(labels.items()):
                if g in predictions and ab in targets and cat in ("S", "I", "R"):
                    w.writerow([g, ab, v, int(ab in predictions[g]), cat])

    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    with open(a.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(COLUMNS)
        abs_ = ["ALL"] + sorted({ab for ab, _ in pairs if ab != "ALL"})
        for ab in abs_:
            for v in versions:
                if (ab, v) not in pairs:
                    continue
                e = error_rates(pairs[(ab, v)])
                w.writerow([ab, v, e.n, undecidable[(ab, v)], e.phenotype_r, e.phenotype_not_r,
                            e.vme, *(f"{x:.4f}" for x in e.vme_rate),
                            e.me, *(f"{x:.4f}" for x in e.me_rate), f"{e.agreement[0]:.4f}"])

    with open(a.summary, "w", encoding="utf-8") as fh:
        fh.write(f"genomes with main-analysis labels: {len(genomes)}\n")
        fh.write(f"genome-drug-version labels dropped for disagreeing records: {label_conflicts}\n")
        if a.matched:
            fh.write(f"matched: genome-drug pairs dropped for lacking an S/I/R label under some version: {unmatched}\n")
        for reason, n in sorted(lost.items(), key=lambda kv: (-kv[1], kv[0])):
            fh.write(f"genomes lost, {reason}: {n}\n")
        fh.write(f"genomes with genotype predictions: {len(predictions)}\n")
        fh.write(f"Pathogen Detection BioSamples dropped for differing entries: {pd_conflicts}\n")
        fh.write(f"antibiotics predicted: {len(targets)}\n")
        fh.write(f"excluded families: {', '.join(sorted(excluded.families)) or 'none'}\n")
        fh.write(f"excluded symbols: {', '.join(sorted(excluded.symbols)) or 'none'}\n")
        fh.write("element calls ignored (by method or excluded family):\n")
        for key, n in sorted(ignored.items(), key=lambda kv: (-kv[1], kv[0])):
            fh.write(f"    {key}: {n}\n")
        fh.write(f"element symbols not in the AMRFinderPlus catalogue: {sum(unknown.values())} "
                 f"({len(unknown)} distinct)\n")
        for sym, n in sorted(unknown.items(), key=lambda kv: (-kv[1], kv[0]))[:25]:
            fh.write(f"    {sym}: {n}\n")
    print(f"wrote {a.output} and {a.summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
