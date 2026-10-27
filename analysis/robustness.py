"""Robustness of the label-change estimates (reviewer requests).

    python analysis/robustness.py data/derived/kpneumoniae_reinterpreted.csv \
        -o results/robustness.csv --transitions results/transitions.csv

Main-analysis records only. Scenarios, each compared on units decidable under
both versions of a pair (paired comparison, so a pair's denominator is fixed):

    record         one unit per MIC record (the main analysis)
    genome         one unit per genome and antibiotic; a genome whose records
                   for that drug disagree under any version is dropped
    uncensored     records whose MIC is an exact value (no <, <=, >, >= and no
                   sign conflict), so no censoring is involved
    all_decidable  records decidable under every version: one fixed cohort, so
                   every version pair compares the same records

Intervals: Wilson score 95% intervals, and percentile 95% intervals from a
bootstrap that resamples genomes (2,000 replicates, fixed seed), which allows
for several records from the same genome. One joint resample of genomes is
used for every antibiotic and version pair.

Outputs counts and rates only; no breakpoint values are written.
"""

from __future__ import annotations

import argparse
import csv
import itertools
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from bpdrift.mic import parse_mic
from bpdrift.stats import wilson

SCENARIOS = ("record", "genome", "uncensored", "all_decidable")
DECIDABLE = {"S", "I", "R"}
REPLICATES = 2000
SEED = 20261027

COLUMNS = [
    "scenario", "antibiotic", "version_a", "version_b", "units", "genomes", "flips", "flip_rate",
    "wilson_lo", "wilson_hi", "boot_lo", "boot_hi", "r_line_flips", "r_line_rate",
    "r_line_boot_lo", "r_line_boot_hi",
]
TRANSITION_COLUMNS = ["scenario", "antibiotic", "version_a", "version_b", "from", "to", "units"]


def _vkey(v: str) -> tuple[int, ...]:
    return tuple(int(x) for x in v.split("."))


def load(path):
    """Main-analysis records: source_line -> dict(genome, antibiotic, exact, labels{version: category})."""
    recs: dict[str, dict] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["exclusion"] != "":
                continue
            r = recs.setdefault(row["source_line"], {
                "genome": row["genome_id"], "antibiotic": row["antibiotic"],
                "exact": parse_mic(row["mic"]).sign == "=" and row["sign_conflict"] == "0",
                "labels": {},
            })
            r["labels"][row["eucast_version"]] = row["category"]
    return recs


def units_for(scenario: str, recs: dict, versions: list[str]):
    """Units for a scenario as (genome, antibiotic, labels) and the count of dropped genome-drug pairs."""
    if scenario == "record":
        return [(r["genome"], r["antibiotic"], r["labels"]) for r in recs.values()], 0
    if scenario == "uncensored":
        return [(r["genome"], r["antibiotic"], r["labels"]) for r in recs.values() if r["exact"]], 0
    if scenario == "all_decidable":
        return [(r["genome"], r["antibiotic"], r["labels"]) for r in recs.values()
                if all(r["labels"].get(v) in DECIDABLE for v in versions)], 0
    if scenario == "genome":
        groups = defaultdict(list)
        for r in recs.values():
            groups[(r["genome"], r["antibiotic"])].append(r["labels"])
        out, dropped = [], 0
        for (g, ab), labs in groups.items():
            merged = {}
            for v in versions:
                cats = {lab[v] for lab in labs if v in lab}
                if len(cats) > 1:
                    merged = None
                    break
                if cats:
                    merged[v] = cats.pop()
            if merged is None:
                dropped += 1
            else:
                out.append((g, ab, merged))
        return out, dropped
    raise ValueError(scenario)


def pair_counts(units, a: str, b: str):
    """Per-genome arrays of (decidable units, flips, R-line flips) for one version pair."""
    n, f, rl = Counter(), Counter(), Counter()
    trans = Counter()
    for g, _, lab in units:
        x, y = lab.get(a), lab.get(b)
        if x in DECIDABLE and y in DECIDABLE:
            n[g] += 1
            trans[(x, y)] += 1
            if x != y:
                f[g] += 1
                if (x == "R") != (y == "R"):
                    rl[g] += 1
    return n, f, rl, trans


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("reinterpreted")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--transitions", required=True)
    p.add_argument("--replicates", type=int, default=REPLICATES)
    a = p.parse_args(argv)

    recs = load(a.reinterpreted)
    versions = sorted({v for r in recs.values() for v in r["labels"]}, key=_vkey)
    genomes = sorted({r["genome"] for r in recs.values()})
    gidx = {g: i for i, g in enumerate(genomes)}
    rng = np.random.default_rng(SEED)
    # joint genome resample: weight = how many times each genome is drawn
    weights = rng.multinomial(len(genomes), np.full(len(genomes), 1 / len(genomes)),
                              size=a.replicates).astype(np.float32)
    pairs = list(zip(versions, versions[1:]))
    if (versions[0], versions[-1]) not in pairs:
        pairs.append((versions[0], versions[-1]))
    all_pairs = list(itertools.combinations(versions, 2))

    out_rows, trans_rows, dropped_note = [], [], {}
    for scenario in SCENARIOS:
        units, dropped = units_for(scenario, recs, versions)
        dropped_note[scenario] = dropped
        by_ab = defaultdict(list)
        for u in units:
            by_ab[u[1]].append(u)
        strata = [("ALL", units)] + [(ab, by_ab[ab]) for ab in sorted(by_ab)]
        vecs, meta = [], []
        for va, vb in all_pairs:
            for name, us in strata:
                n, f, rl, trans = pair_counts(us, va, vb)
                tot, flips, rls = sum(n.values()), sum(f.values()), sum(rl.values())
                if tot == 0:
                    continue
                if (va, vb) in pairs:
                    for (x, y), k in sorted(trans.items()):
                        trans_rows.append([scenario, name, va, vb, x, y, k])
                    for c in (n, f, rl):
                        v = np.zeros(len(genomes), dtype=np.float32)
                        for g, k in c.items():
                            v[gidx[g]] = k
                        vecs.append(v)
                    boot = True
                else:
                    boot = False
                meta.append((name, va, vb, tot, len(n), flips, rls, boot))
        boot_ci = {}
        if vecs:
            sums = weights @ np.stack(vecs, axis=1)          # replicates x (3 per stratum-pair)
            j = 0
            for name, va, vb, *_rest, boot in meta:
                if not boot:
                    continue
                n_b, f_b, r_b = sums[:, j], sums[:, j + 1], sums[:, j + 2]
                j += 3
                with np.errstate(invalid="ignore", divide="ignore"):
                    fr, rr = f_b / n_b, r_b / n_b
                boot_ci[(name, va, vb)] = (*np.nanpercentile(fr, [2.5, 97.5]), *np.nanpercentile(rr, [2.5, 97.5]))
        for name, va, vb, tot, ng, flips, rls, boot in meta:
            fr, lo, hi = wilson(flips, tot)
            b = boot_ci.get((name, va, vb))
            fmt = (lambda x: f"{x:.4f}")
            bf = [fmt(x) for x in b] if b else ["", "", "", ""]
            out_rows.append([scenario, name, va, vb, tot, ng, flips, fmt(fr), fmt(lo), fmt(hi),
                             bf[0], bf[1], rls, fmt(rls / tot), bf[2], bf[3]])

    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    with open(a.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(COLUMNS)
        w.writerows(out_rows)
    with open(a.transitions, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(TRANSITION_COLUMNS)
        w.writerows(trans_rows)
    print(f"wrote {a.output} and {a.transitions}; genome-drug pairs dropped for disagreeing "
          f"records: {dropped_note['genome']}; bootstrap: {a.replicates} genome resamples, seed {SEED}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
