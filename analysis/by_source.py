"""Leave-one-source-out: does one depositing study drive the label-change estimates?

    python analysis/by_source.py data/derived/kpneumoniae_reinterpreted.csv \
        --bvbrc data/bvbrc/kpneumoniae_amr_all.csv -o results/by_source.csv

A record's source is the set of PubMed IDs BV-BRC gives for it (order ignored;
"none" when empty), which identifies the study that deposited it. Main-analysis
records only, compared between the first and last EUCAST version on records
decidable under both. For every source with records the output gives:

    own        the source's own change rate
    leave_out  the rate in all other sources together

for all antibiotics and for each antibiotic with at least one change. Counts and
rates only; no breakpoint values are written.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from bpdrift.stats import wilson

DECIDABLE = {"S", "I", "R"}
COLUMNS = [
    "source", "antibiotic", "version_a", "version_b", "genomes", "testing_standard", "platform",
    "own_units", "own_flips", "own_rate", "own_lo", "own_hi",
    "leave_out_units", "leave_out_flips", "leave_out_rate", "leave_out_lo", "leave_out_hi",
]


def _vkey(v):
    return tuple(int(x) for x in v.split("."))


def source_of(pubmed: str) -> str:
    ids = sorted({p.strip() for p in pubmed.replace(",", ";").split(";") if p.strip()})
    return ";".join(ids) or "none"


def read_sources(bvbrc_csv) -> dict[str, dict]:
    """source_line (as written by reinterpret) -> source, standard and platform of that line."""
    out = {}
    with open(bvbrc_csv, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            out[str(reader.line_num)] = {
                "source": source_of(row.get("PubMed", "")),
                "standard": (row.get("Testing Standard") or "").strip() or "not stated",
                "platform": (row.get("Laboratory Typing Platform") or "").strip() or "not stated",
            }
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("reinterpreted")
    p.add_argument("--bvbrc", required=True, help="the BV-BRC export the reinterpreted table was made from")
    p.add_argument("-o", "--output", required=True)
    a = p.parse_args(argv)

    sources = read_sources(a.bvbrc)
    recs: dict[str, dict] = {}
    with open(a.reinterpreted, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["exclusion"] != "":
                continue
            r = recs.setdefault(row["source_line"], {"genome": row["genome_id"], "ab": row["antibiotic"], "lab": {}})
            r["lab"][row["eucast_version"]] = row["category"]
    missing = [k for k in recs if k not in sources]
    if missing:
        raise SystemExit(f"{len(missing)} records have no line in {a.bvbrc} (first: {missing[0]}); "
                         "is it the export the table was made from?")
    versions = sorted({v for r in recs.values() for v in r["lab"]}, key=_vkey)
    va, vb = versions[0], versions[-1]

    units = defaultdict(Counter)    # (source, antibiotic) -> {"n", "f"}
    genomes, standards, platforms = defaultdict(set), defaultdict(Counter), defaultdict(Counter)
    for line, r in recs.items():
        s = sources[line]
        src = s["source"]
        genomes[src].add(r["genome"])
        standards[src][s["standard"]] += 1
        platforms[src][s["platform"]] += 1
        x, y = r["lab"].get(va), r["lab"].get(vb)
        if x in DECIDABLE and y in DECIDABLE:
            for ab in ("ALL", r["ab"]):
                units[(src, ab)]["n"] += 1
                units[(src, ab)]["f"] += x != y
    total = defaultdict(Counter)
    for (src, ab), c in units.items():
        total[ab].update(c)
    changing = sorted(ab for ab, c in total.items() if ab != "ALL" and c["f"] > 0)

    srcs = sorted(genomes, key=lambda s: -len(genomes[s]))
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    with open(a.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(COLUMNS)
        for src in srcs:
            for ab in ["ALL"] + changing:
                own = units.get((src, ab), Counter())
                rest_n, rest_f = total[ab]["n"] - own["n"], total[ab]["f"] - own["f"]
                if own["n"] == 0 and ab != "ALL":
                    continue
                o, l = wilson(own["f"], own["n"]), wilson(rest_f, rest_n)
                w.writerow([src, ab, va, vb, len(genomes[src]),
                            standards[src].most_common(1)[0][0], platforms[src].most_common(1)[0][0],
                            own["n"], own["f"], *(f"{x:.4f}" for x in o),
                            rest_n, rest_f, *(f"{x:.4f}" for x in l)])
    print(f"wrote {a.output}: {len(srcs)} sources, {len(changing)} antibiotics with changes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
