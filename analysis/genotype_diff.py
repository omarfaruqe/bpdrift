"""Paired change in genotype-predictor error rates between two EUCAST versions.

    python analysis/genotype_diff.py data/derived/genotype_pairs.csv -o results/genotype_diff.csv
    python analysis/genotype_diff.py data/derived/genotype_pairs.csv --i-as drop \
        -o results/genotype_diff_s_r.csv

Input: the --pairs-out file of analysis/genotype_errors.py run with --matched,
so every genome-drug pair is labelled under both versions and the prediction is
identical; only the phenotype label differs.

Error definitions are operational, for a binary (resistant / not resistant)
predictor. By default EUCAST I ("susceptible, increased exposure") is grouped with
S as not resistant. Sensitivity analyses: --i-as r groups I with R (non-susceptible
versus susceptible); --i-as drop leaves out pairs labelled I under either version,
so every remaining pair is S or R under both. The testing family is the antibiotics
whose labels change between the versions (after dropping, under drop).

Three estimands per antibiotic (and ALL), between the first and last version,
each a difference (b - a):

    very major error rate  VME / phenotype R           denominator changes with the labels
    major error rate       ME  / phenotype not R       denominator changes with the labels
    overall error rate     (VME + ME) / all pairs      fixed denominator; the estimand of
                                                       the exact McNemar test on pairs that
                                                       change error status

Intervals: 95% percentile intervals from a bootstrap that resamples genomes
(2,000 replicates, fixed seed). Bootstrap sums are exact integer arithmetic.
<measure>_boot_valid is the number of replicates in which both rates are defined
(a replicate that draws no isolates in a denominator gives no estimate). ALL is the primary, pooled analysis.
Drug-level results are exploratory. Two multiplicity families are used, both over
the same set of antibiotics:

    Holm        one McNemar test per antibiotic (overall error rate): k tests
    Bonferroni  three interval estimates per antibiotic (VME, ME, overall): 3k
                intervals, each at level 1 - 0.05 / 3k, computed as
                estimate +/- z * bootstrap standard error

The set of antibiotics is, by default, those with at least one pair whose label
changes between the two versions (for the others every difference is zero by
construction). That set is chosen after seeing the labels, so it is an exploratory
choice; --family all uses every antibiotic compared instead (pre-specified, more
conservative).

Counts and rates only.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

from statistics import NormalDist

import numpy as np

REPLICATES = 2000
SEED = 20261027
MEASURES = ("vme", "me", "err")
COLUMNS = (
    ["antibiotic", "version_a", "version_b", "pairs", "genomes"]
    + [f"{m}_{x}" for m in MEASURES for x in ("num_a", "den_a", "rate_a", "num_b", "den_b", "rate_b",
                                               "diff", "lo", "hi", "lo_adj", "hi_adj", "boot_valid")]
    + ["error_to_correct", "correct_to_error", "mcnemar_p", "mcnemar_p_holm"]
)


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
    return min(1.0, 2 * math.exp(m) * sum(math.exp(x - m) for x in logs))


def holm(pvalues: list[float]) -> list[float]:
    """Holm step-down adjusted p-values, in the input order."""
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    adj, running = [0.0] * len(pvalues), 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(pvalues) - rank) * pvalues[i]))
        adj[i] = running
    return adj


def error(pred: int, cat: str) -> str:
    """Operational errors for a binary predictor; I counts as not resistant."""
    if cat == "R":
        return "vme" if not pred else ""
    return "me" if pred else ""


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("pairs")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--i-as", choices=("not_r", "r", "drop"), default="not_r",
                   help="how phenotype I enters the binary outcome: not_r (default; grouped with S), "
                        "r (grouped with R, i.e. non-susceptible), or drop (pairs labelled I under "
                        "either version are left out)")
    p.add_argument("--family", choices=("changed", "all"), default="changed",
                   help="drug-level testing family: antibiotics whose labels change between the two "
                        "versions (default; chosen after seeing the labels, so exploratory) or every "
                        "antibiotic compared (pre-specified, more conservative)")
    p.add_argument("--replicates", type=int, default=REPLICATES)
    a = p.parse_args(argv)

    lab = defaultdict(dict)       # (genome, antibiotic) -> version -> (pred, cat)
    with open(a.pairs, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            lab[(r["genome_id"], r["antibiotic"])][r["eucast_version"]] = (int(r["predicted_resistant"]),
                                                                          r["phenotype"])
    versions = sorted({v for d in lab.values() for v in d}, key=_vkey)
    va, vb = versions[0], versions[-1]
    pairs = {k: d for k, d in lab.items() if va in d and vb in d}
    if any(d[va][0] != d[vb][0] for d in pairs.values()):
        raise SystemExit("predictions differ between versions; expected a fixed predictor")
    changed = {ab for (g, ab), d in pairs.items() if d[va][1] != d[vb][1]}
    dropped_i = 0
    if a.i_as == "drop":
        before = len(pairs)
        pairs = {k: d for k, d in pairs.items() if d[va][1] != "I" and d[vb][1] != "I"}
        dropped_i = before - len(pairs)
        changed = {ab for (g, ab), d in pairs.items() if d[va][1] != d[vb][1]}
    elif a.i_as == "r":
        pairs = {k: {v: (pr, "R" if c == "I" else c) for v, (pr, c) in d.items()} for k, d in pairs.items()}

    genomes = sorted({g for g, _ in pairs})
    gi = {g: i for i, g in enumerate(genomes)}
    rng = np.random.default_rng(SEED)
    # integer resampling weights; all bootstrap sums are computed in exact int64
    # arithmetic (no floating-point BLAS), so they are identical on every platform
    w = rng.multinomial(len(genomes), np.full(len(genomes), 1 / len(genomes)), size=a.replicates)

    by_ab = defaultdict(list)
    for (g, ab), d in pairs.items():
        by_ab[ab].append((g, d))
        by_ab["ALL"].append((g, d))
    names = ["ALL"] + sorted(k for k in by_ab if k != "ALL")

    results = {}
    for ab in names:
        # per-genome counts, columns: (num_a, den_a, num_b, den_b) for vme, me, err
        m = np.zeros((len(genomes), 12), dtype=np.int64)
        b = c = 0
        for g, d in by_ab[ab]:
            (pa, ca), (pb, cb) = d[va], d[vb]
            i, ea, eb = gi[g], error(pa, ca), error(pb, cb)
            m[i] += [ea == "vme", ca == "R", eb == "vme", cb == "R",
                     ea == "me", ca != "R", eb == "me", cb != "R",
                     bool(ea), 1, bool(eb), 1]
            b += bool(ea) and not eb
            c += bool(eb) and not ea
        tot, boot = m.sum(axis=0), np.einsum("rg,gk->rk", w, m)
        if boot.min() < 0:
            raise SystemExit(f"{ab}: negative bootstrap count; integer overflow?")
        est = {}
        for j, name in enumerate(MEASURES):
            na, da, nb, db = tot[4 * j: 4 * j + 4]
            with np.errstate(invalid="ignore", divide="ignore"):
                diffs = boot[:, 4 * j + 2] / boot[:, 4 * j + 3] - boot[:, 4 * j] / boot[:, 4 * j + 1]
            ra = na / da if da else math.nan
            rb = nb / db if db else math.nan
            # a replicate is invalid only when it draws no pairs in a denominator
            # (rate undefined); the count of valid replicates is reported
            est[name] = (na, da, ra, nb, db, rb, rb - ra, diffs[np.isfinite(diffs)])
        results[ab] = (len(by_ab[ab]), len({g for g, _ in by_ab[ab]}), est, b, c, mcnemar_exact(b, c))

    # testing family: antibiotics with at least one pair whose label changes; for the
    # others every difference is zero by construction and nothing is tested
    drugs = [ab for ab in names if ab != "ALL" and (a.family == "all" or ab in changed)]
    n_tests = 3 * len(drugs)
    alpha_adj = 0.05 / n_tests if n_tests else 0.05
    z_adj = NormalDist().inv_cdf(1 - alpha_adj / 2)
    holm_p = dict(zip(drugs, holm([results[ab][5] for ab in drugs])))

    def ci(x):
        if len(x) == 0:
            return math.nan, math.nan
        return float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))

    def ci_adj(d, x):
        """Bonferroni interval: estimate +/- z * bootstrap SE (stable at small tail areas)."""
        if len(x) < 2 or math.isnan(d):
            return math.nan, math.nan
        se = float(np.std(x, ddof=1))
        return d - z_adj * se, d + z_adj * se

    f = lambda x: "" if isinstance(x, float) and math.isnan(x) else f"{x:.4f}"  # noqa: E731
    rows = []
    for ab in names:
        n, ng, est, b, c, pval = results[ab]
        row = [ab, va, vb, n, ng]
        for name in MEASURES:
            na, da, ra, nb, db, rb, d, diffs = est[name]
            lo, hi = ci(diffs)
            lo2, hi2 = ci_adj(d, diffs) if ab in holm_p else (math.nan, math.nan)
            row += [int(na), int(da), f(ra), int(nb), int(db), f(rb), f(d), f(lo), f(hi), f(lo2), f(hi2),
                    len(diffs)]
        row += [b, c, f"{pval:.3g}", f"{holm_p[ab]:.3g}" if ab in holm_p else ""]
        rows.append(row)

    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    with open(a.output, "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh, lineterminator="\n")
        wr.writerow(COLUMNS)
        wr.writerows(rows)
    print(f"wrote {a.output}: v{va} vs v{vb}, {len(pairs)} pairs"
          + f"; I as {a.i_as}" + (f" ({dropped_i} pairs with I left out)" if a.i_as == "drop" else "")
          + f"; testing family ({a.family}): {len(drugs)} antibiotics, Holm over {len(drugs)} McNemar tests, "
          f"Bonferroni over {n_tests} interval estimates "
          f"(Bonferroni intervals at {100 * (1 - alpha_adj):.2f}%); {a.replicates} genome resamples, seed {SEED}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
