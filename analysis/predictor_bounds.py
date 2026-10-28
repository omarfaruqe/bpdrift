"""How far could ANY fixed predictor's measured error move between two breakpoint versions?

    python analysis/predictor_bounds.py data/derived/genotype_pairs.csv -o results/predictor_bounds.csv

Input: the --pairs-out file of analysis/genotype_errors.py run with --matched
(identical genome-drug pairs under every version).

For a binary outcome (resistant or not), a pair whose binary label is the same
under both versions keeps its error status, whatever the predictor says. A pair
whose label changes always switches status: correct under one version, an error
under the other. So for any predictor held fixed,

    change in overall error rate = (c - b) / n,   with b + c = B,

where n is the number of pairs, B the number whose binary label changes, b the
pairs that stop being errors and c those that become errors. The change can
therefore never exceed B / n in size, for any predictor, and its extremes are
reached by a predictor that reproduces the old labels exactly (+B/n: every
changed pair becomes an error) or the new ones (-B/n). B depends only on the
labels. The output gives B / n per antibiotic, cohort and treatment of I, with
the error rates of the "old-label oracle" under the new version, and the observed
predictor's position in that range, (c - b) / B, from -1 to +1.

Counts and rates only.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from bpdrift.stats import wilson

TREATMENTS = {"not_r": {"S": 0, "I": 0, "R": 1}, "r": {"S": 0, "I": 1, "R": 1}, "drop": {"S": 0, "R": 1}}
COLUMNS = ["i_as", "antibiotic", "version_a", "version_b", "pairs", "label_changes", "bound", "bound_lo", "bound_hi",
           "oracle_vme_rate", "oracle_me_rate", "observed_err_change", "errors_lost", "errors_gained",
           "observed_position"]


def _vkey(v):
    return tuple(int(x) for x in v.split("."))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("pairs")
    p.add_argument("-o", "--output", required=True)
    a = p.parse_args(argv)

    lab = defaultdict(dict)
    with open(a.pairs, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            lab[(r["genome_id"], r["antibiotic"])][r["eucast_version"]] = (int(r["predicted_resistant"]), r["phenotype"])
    versions = sorted({v for d in lab.values() for v in d}, key=_vkey)
    va, vb = versions[0], versions[-1]

    rows = []
    for name, code in TREATMENTS.items():
        acc = defaultdict(lambda: defaultdict(int))
        for (g, ab), d in lab.items():
            if va not in d or vb not in d:
                continue
            (pred, ca), (_, cb) = d[va], d[vb]
            if ca not in code or cb not in code:          # 'drop' leaves out pairs labelled I
                continue
            ya, yb = code[ca], code[cb]
            for key in ("ALL", ab):
                s = acc[key]
                s["n"] += 1
                if ya != yb:
                    s["B"] += 1
                    # old-label oracle predicts ya; under the new labels it is wrong on every changed pair
                    s["oracle_vme"] += yb == 1        # truly R now, oracle said not R
                    s["oracle_me"] += yb == 0
                    err_a, err_b = pred != ya, pred != yb
                    s["b"] += err_a and not err_b
                    s["c"] += err_b and not err_a
                s["r_b"] += yb == 1
                s["notr_b"] += yb == 0
        for ab in ["ALL"] + sorted(k for k in acc if k != "ALL"):
            s = acc[ab]
            if s["n"] == 0:
                continue
            if s["b"] + s["c"] != s["B"]:
                raise SystemExit(f"{name} {ab}: changed pairs {s['B']} != errors lost + gained {s['b'] + s['c']}")
            bound, lo, hi = wilson(s["B"], s["n"])
            f = lambda x: f"{x:.4f}"  # noqa: E731
            rows.append([name, ab, va, vb, s["n"], s["B"], f(bound), f(lo), f(hi),
                         f(s["oracle_vme"] / s["r_b"]) if s["r_b"] else "",
                         f(s["oracle_me"] / s["notr_b"]) if s["notr_b"] else "",
                         f((s["c"] - s["b"]) / s["n"]), s["b"], s["c"],
                         f((s["c"] - s["b"]) / s["B"]) if s["B"] else ""])

    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    with open(a.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(COLUMNS)
        w.writerows(rows)
    print(f"wrote {a.output}: v{va} vs v{vb}, {len(rows)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
