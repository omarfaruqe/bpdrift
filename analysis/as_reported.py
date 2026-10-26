"""Agreement between BV-BRC's reported S/I/R calls and the re-interpreted labels.

    python analysis/as_reported.py data/derived/kpneumoniae_reinterpreted.csv \
        -o results/as_reported.csv

Main-analysis records only. A record is compared when it has a reported
call (Susceptible, Intermediate, Resistant) and its re-interpreted label
under that version is decidable (S, I, R).

    agreement         reported call == re-interpreted category
    r_line_agreement  both resistant, or both not resistant

Strata: every stated standard (CLSI, EUCAST, none, both), all years
together and each year separately. A stated EUCAST year compared with the
table version that was in force that year is the pipeline's validation:
disagreement there cannot come from breakpoint drift.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from bpdrift.stats import wilson

REPORTED = {"Susceptible": "S", "Intermediate": "I", "Resistant": "R"}
DECIDABLE = {"S", "I", "R"}

# Year each loaded EUCAST version took effect (EUCAST publishes on 1 January).
IN_FORCE = {"10.0": 2020, "13.1": 2023, "15.0": 2025, "16.1": 2026}

COLUMNS = [
    "standard", "year", "eucast_version", "in_force_that_year", "reported", "compared",
    "agreement", "agreement_rate", "agreement_lo", "agreement_hi",
    "r_line_agreement", "r_line_rate", "r_line_lo", "r_line_hi", "disagreements",
]


def _stratum_counts(rows):
    """(standard, year, version) -> [reported, compared, agree, r_agree, Counter of disagreements]."""
    acc = defaultdict(lambda: [0, 0, 0, 0, Counter()])
    for r in rows:
        if r["exclusion"] != "":
            continue
        reported = REPORTED.get(r["reported_phenotype"])
        if reported is None:
            continue
        std = r["reported_standard"] or "none"
        for year in ("all", r["reported_standard_year"] or "none"):
            a = acc[(std, year, r["eucast_version"])]
            a[0] += 1
            cat = r["category"]
            if cat not in DECIDABLE:
                continue
            a[1] += 1
            if cat == reported:
                a[2] += 1
            else:
                a[4][f"{reported}>{cat}"] += 1
            if (cat == "R") == (reported == "R"):
                a[3] += 1
    return acc


def _version_key(v):
    return tuple(int(p) for p in v.split("."))


def in_force(standard: str, year: str, version: str) -> str:
    """'yes' when a stated EUCAST year is the year this version took effect.

    EUCAST issues one main version a year, so for the loaded versions this is
    exact for 10.0, 15.0 and 16.1; 13.1 replaced 13.0 during 2023.
    """
    if standard != "EUCAST" or not year.isdigit() or version not in IN_FORCE:
        return ""
    return "yes" if int(year) == IN_FORCE[version] else "no"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("reinterpreted")
    p.add_argument("-o", "--output", required=True)
    args = p.parse_args(argv)

    with open(args.reinterpreted, newline="", encoding="utf-8") as fh:
        acc = _stratum_counts(csv.DictReader(fh))

    def order(k):
        std, year, v = k
        return (std, year != "all", year, _version_key(v))

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(COLUMNS)
        for key in sorted(acc, key=order):
            std, year, v = key
            reported, compared, agree, r_agree, dis = acc[key]
            a, rl = wilson(agree, compared), wilson(r_agree, compared)
            w.writerow([
                std, year, v, in_force(std, year, v), reported, compared,
                agree, *(f"{x:.4f}" for x in a), r_agree, *(f"{x:.4f}" for x in rl),
                ";".join(f"{k}:{n}" for k, n in dis.most_common()),
            ])
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
