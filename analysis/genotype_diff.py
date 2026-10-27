"""Paired change in genotype-predictor error rates between two EUCAST versions.

    python analysis/genotype_diff.py data/derived/genotype_pairs.csv -o results/genotype_diff.csv

Input: the --pairs-out file of analysis/genotype_errors.py run with --matched,
so every genome-drug pair is labelled under both versions and the predictions
are identical; only the phenotype label differs.

For each antibiotic (and ALL), between the first and last version:

    very major error rate  VME / phenotype R          (I counted as not resistant)
    major error rate       ME  / phenotype S or I
    difference (b - a)     with a 95% percentile interval from a bootstrap that
                           resamples genomes (2,000 replicates, fixed seed), so
                           pairs from the same genome stay together
    discordant pairs       pairs that are an error (VME or ME) under one version
                           and correct under the other, with an exact two-sided
                           McNemar test of whether errors are as likely to appear
                           as to disappear

Counts and rates only.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

REPLICATES = 2000
SEED = 20261027
COLUMNS = [
    "antibiotic", "version_a", "version_b", "pairs", "genomes",
    "vme_a", "r_a", "vme_rate_a", "vme_b", "r_b", "vme_rate_b", "vme_diff", "vme_diff_lo", "vme_diff_hi",
    "me_a", "notr_a", "me_rate_a", "me_b", "notr_b", "me_rate_b", "me_diff", "me_diff_lo", "me_diff_hi",
    "error_to_correct", "correct_to_error", "mcnemar_p",
]


def _vkey(v):
    return tuple(int(x) for x in v.split("."))


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact McNemar p-value: binomial test of b vs c with p = 1/2."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    logs = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) - n * math.log(2)
            for i in range(k + 1)]
    m = max(logs)
    tail = math.exp(m) * sum(math.exp(x - m) for x in logs)
    return min(1.0, 2 * tail)


def error(pred: int, cat: str) -> str:
    if cat == "R":
        return "vme" if not pred else ""
    return "me" if pred else ""


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("pairs")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--replicates", type=int, default=REPLICATES)
    a = p.parse_args(argv)

    lab = defaultdict(dict)       # (genome, antibiotic) -> version -> (pred, cat)
    with open(a.pairs, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            lab[(r["genome_id"], r["antibiotic"])][r["eucast_version"]] = (int(r["predicted_resistant"]), r["phenotype"])
    versions = sorted({v for d in lab.values() for v in d}, key=_vkey)
    va, vb = versions[0], versions[-1]
    pairs = {k: d for k, d in lab.items() if va in d and vb in d}
    unpaired = len(lab) - len(pairs)
    if any(d[va][0] != d[vb][0] for d in pairs.values()):
        raise SystemExit("predictions differ between versions; expected a fixed predictor")

    genomes = sorted({g for g, _ in pairs})
    gi = {g: i for i, g in enumerate(genomes)}
    rng = np.random.default_rng(SEED)
    w = rng.multinomial(len(genomes), np.full(len(genomes), 1 / len(genomes)),
                        size=a.replicates).astype(np.float64)

    by_ab = defaultdict(list)
    for (g, ab), d in pairs.items():
        by_ab[ab].append((g, d))
        by_ab["ALL"].append((g, d))

    rows = []
    for ab in ["ALL"] + sorted(k for k in by_ab if k != "ALL"):
        # per-genome counts: vme_a, r_a, vme_b, r_b, me_a, notr_a, me_b, notr_b
        m = np.zeros((len(genomes), 8))
        b = c = 0
        for g, d in by_ab[ab]:
            (pa, ca), (pb, cb) = d[va], d[vb]
            i = gi[g]
            ea, eb = error(pa, ca), error(pb, cb)
            m[i, 0] += ea == "vme"; m[i, 1] += ca == "R"; m[i, 2] += eb == "vme"; m[i, 3] += cb == "R"
            m[i, 4] += ea == "me"; m[i, 5] += ca != "R"; m[i, 6] += eb == "me"; m[i, 7] += cb != "R"
            b += bool(ea) and not eb
            c += bool(eb) and not ea
        tot = m.sum(axis=0)
        boot = w @ m
        with np.errstate(invalid="ignore", divide="ignore"):
            vme_d = boot[:, 2] / boot[:, 3] - boot[:, 0] / boot[:, 1]
            me_d = boot[:, 6] / boot[:, 7] - boot[:, 4] / boot[:, 5]

        def rate(k, n):
            return k / n if n else float("nan")

        def ci(x):
            x = x[np.isfinite(x)]
            return (np.percentile(x, 2.5), np.percentile(x, 97.5)) if len(x) else (float("nan"),) * 2

        va_r, vb_r = rate(tot[0], tot[1]), rate(tot[2], tot[3])
        ma_r, mb_r = rate(tot[4], tot[5]), rate(tot[6], tot[7])
        f = lambda x: f"{x:.4f}"  # noqa: E731
        rows.append([ab, va, vb, len(by_ab[ab]), len({g for g, _ in by_ab[ab]}),
                     int(tot[0]), int(tot[1]), f(va_r), int(tot[2]), int(tot[3]), f(vb_r), f(vb_r - va_r),
                     *map(f, ci(vme_d)),
                     int(tot[4]), int(tot[5]), f(ma_r), int(tot[6]), int(tot[7]), f(mb_r), f(mb_r - ma_r),
                     *map(f, ci(me_d)),
                     b, c, f"{mcnemar_exact(b, c):.3g}"])

    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    with open(a.output, "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh, lineterminator="\n")
        wr.writerow(COLUMNS)
        wr.writerows(rows)
    print(f"wrote {a.output}: v{va} vs v{vb}, {len(pairs)} pairs ({unpaired} unpaired left out), "
          f"{a.replicates} genome resamples, seed {SEED}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
