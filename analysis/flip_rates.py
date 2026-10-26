"""Label changes between EUCAST versions for the re-interpreted BV-BRC data.

    python analysis/flip_rates.py data/derived/kpneumoniae_reinterpreted.csv \
        -o results/flip_rates.csv

Scenarios (see mappings/README.md for the decisions behind them):

    main          main-analysis records, signs as recorded
    sign_exact    main-analysis records, sign-conflict MICs taken as exact
    all_records   every linked record, including urinary-only breakpoints,
                  ratio-method combination MICs and conflicting repeats

One output row per (scenario, antibiotic or ALL, version pair). Counts and
rates only: no breakpoint values are written.
"""

from __future__ import annotations

import argparse
import csv
import itertools
from collections import defaultdict
from pathlib import Path

from bpdrift.interpret import Breakpoint, Category, interpret
from bpdrift.mic import MIC, parse_mic
from bpdrift.stats import compare

SCENARIOS = ("main", "sign_exact", "all_records")

COLUMNS = [
    "scenario", "antibiotic", "version_a", "version_b", "records", "both_decidable",
    "flips", "flip_rate", "flip_lo", "flip_hi",
    "r_line_flips", "r_line_rate", "r_line_lo", "r_line_hi",
    "undecidable_a", "undecidable_rate_a", "undecidable_b", "undecidable_rate_b",
    "no_breakpoint_a", "no_breakpoint_b", "transitions",
]


def _bp_value(raw: str) -> float | None:
    t = raw.strip()
    if t.startswith("(") and t.endswith(")"):
        t = t[1:-1]
    try:
        v = float(t)
    except ValueError:
        return None
    return v if v > 0 else None


def _label(row: dict, exact: bool) -> Category:
    category = Category(row["category"])
    if not exact or category is Category.NO_BREAKPOINT:
        return category
    s, r = _bp_value(row["eucast_S"]), _bp_value(row["eucast_R"])
    if s is None or r is None:
        return category
    mic = parse_mic(row["mic"])
    return interpret(MIC(mic.value), Breakpoint(s, r)).category


def load(path: str | Path) -> dict[str, dict[str, dict[str, dict[str, Category]]]]:
    """scenario -> antibiotic -> record key -> version -> category."""
    out = {s: defaultdict(lambda: defaultdict(dict)) for s in SCENARIOS}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            key, ab, v = row["source_line"], row["antibiotic"], row["eucast_version"]
            main = row["exclusion"] == ""
            out["all_records"][ab][key][v] = Category(row["category"])
            if main:
                out["main"][ab][key][v] = Category(row["category"])
                out["sign_exact"][ab][key][v] = _label(row, exact=row["sign_conflict"] == "1")
    return out


def _version_key(v: str) -> tuple[int, ...]:
    return tuple(int(p) for p in v.split("."))


def rows_for(scenario: str, antibiotic: str, labels: list[dict], versions: list[str]):
    for a, b in itertools.combinations(versions, 2):
        c = compare(labels, a, b)
        f, r, ua, ub = c.flip_rate, c.r_line_rate, c.undecidable_rate_a, c.undecidable_rate_b
        top = ";".join(f"{x}>{y}:{n}" for (x, y), n in c.transitions.most_common() if x != y)
        yield [
            scenario, antibiotic, a, b, c.records, c.both_decidable,
            c.flips, *(f"{x:.4f}" for x in f),
            c.r_line_flips, *(f"{x:.4f}" for x in r),
            c.undecidable_a, f"{ua[0]:.4f}", c.undecidable_b, f"{ub[0]:.4f}",
            c.no_breakpoint_a, c.no_breakpoint_b, top,
        ]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("reinterpreted")
    p.add_argument("-o", "--output", required=True)
    args = p.parse_args(argv)

    data = load(args.reinterpreted)
    versions = sorted({v for recs in data["all_records"].values() for r in recs.values() for v in r},
                      key=_version_key)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(COLUMNS)
        for scenario in SCENARIOS:
            by_ab = data[scenario]
            everything = [lab for recs in by_ab.values() for lab in recs.values()]
            w.writerows(rows_for(scenario, "ALL", everything, versions))
            for ab in sorted(by_ab):
                w.writerows(rows_for(scenario, ab, list(by_ab[ab].values()), versions))
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
