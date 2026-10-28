"""Checks requested in manuscript review (EUCAST conventions, validation sources, cohort scope).

    python analysis/review_checks.py --bvbrc-re data/derived/kpneumoniae_reinterpreted.csv \
        --bvbrc-export data/bvbrc/kpneumoniae_amr_all.csv \
        --ncbi-re data/derived/ncbi_reinterpreted.csv --ncbi-genomes data/derived/ncbi_ast_genomes.csv \
        --ncbi-export data/derived/ncbi_ast_bvbrc_format.csv \
        --links mappings/bvbrc_to_eucast_enterobacterales.csv --tables data/breakpoints/*.xlsx \
        -o results/review_checks.csv

All comparisons are main-analysis records, first against last EUCAST version,
decidable under both. Checks:

    dilution         for each antibiotic with changes: how many distinct MIC values the
                     changed records have, and the share at the most common one
    bracketed        changes where either version's breakpoint is bracketed (EUCAST's mark
                     for restricted use, e.g. aminoglycosides in systemic infection), and
                     the change rate with those records left out
    atu              changes whose exact MIC lies in either version's area of technical
                     uncertainty (from v9.0 on), and the change rate with them left out
    censored_bound   labels a rule that reads a censored MIC at its bound (<=x as x) would
                     categorise, although their true category is not determined
    validation       sources (PubMed IDs) of the calls stated as EUCAST 2020, and how many
                     carry a computational-method entry
    ncbi_host        NCBI cohort change rate for human-host BioSamples only
    ncbi_phenotype   NCBI reported phenotypes outside S/I/R (not used in agreement analyses)

Counts and rates only; no breakpoint or ATU values are written.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

from bpdrift.eucast import read_eucast_xlsx
from bpdrift.mic import parse_mic
from bpdrift.stats import wilson

DECIDABLE = {"S", "I", "R"}
COLUMNS = ["check", "cohort", "antibiotic", "measure", "value", "n", "rate", "lo", "hi"]


def _vkey(v):
    return tuple(int(x) for x in v.split("."))


def bracketed(raw: str) -> bool:
    t = raw.strip()
    return t.startswith("(") and t.endswith(")")


def atu_range(atu: str) -> tuple[float, float] | None:
    """'16' -> (16, 16); '0.5-1' -> (0.5, 1); anything else (blank, '-', notes) -> None."""
    nums = re.findall(r"\d+(?:\.\d+)?", atu or "")
    if not nums or len(nums) > 2:
        return None
    lo, hi = float(nums[0]), float(nums[-1])
    return (min(lo, hi), max(lo, hi))


def in_atu(mic: str, rng: tuple[float, float] | None) -> bool:
    if rng is None:
        return False
    m = parse_mic(mic)
    return m.sign == "=" and rng[0] <= m.value <= rng[1]


def load_records(path):
    """source_line -> {genome, ab, mic, standard, year, labels{v: cat}, S{v}, R{v}}"""
    recs = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["exclusion"] != "":
                continue
            r = recs.setdefault(row["source_line"], {
                "genome": row["genome_id"], "ab": row["antibiotic"], "mic": row["mic"],
                "standard": row["reported_standard"], "year": row["reported_standard_year"],
                "labels": {}, "S": {}, "R": {}})
            v = row["eucast_version"]
            r["labels"][v], r["S"][v], r["R"][v] = row["category"], row["eucast_S"], row["eucast_R"]
    return recs


def atu_table(links_csv, tables) -> dict[tuple[str, str], tuple[float, float] | None]:
    with open(links_csv, encoding="utf-8") as fh:
        lines = [t for t in fh if not t.startswith("#")]
    links = {(r["bvbrc_antibiotic"], r["eucast_version"]): r["eucast_agent"] for r in csv.DictReader(lines)}
    out = {}
    for path in tables:
        t = read_eucast_xlsx(path)
        sheet = [s for s in t.sheets if s.startswith("Enterobacter")][0]
        rows = {r.agent: r for r in t.for_sheet(sheet)}
        for (ab, v), agent in links.items():
            if v == t.version and agent in rows:
                out[(ab, v)] = atu_range(rows[agent].atu)
    return out


def changed_pairs(recs, va, vb):
    for key, r in recs.items():
        a, b = r["labels"].get(va), r["labels"].get(vb)
        if a in DECIDABLE and b in DECIDABLE:
            yield key, r, a != b


def rate_row(check, cohort, ab, measure, k, n):
    p, lo, hi = wilson(k, n)
    return [check, cohort, ab, measure, k, n, f"{p:.4f}", f"{lo:.4f}", f"{hi:.4f}"]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--bvbrc-re", required=True)
    p.add_argument("--bvbrc-export", required=True)
    p.add_argument("--ncbi-re", required=True)
    p.add_argument("--ncbi-genomes", required=True)
    p.add_argument("--ncbi-export", required=True)
    p.add_argument("--links", required=True)
    p.add_argument("--tables", nargs="+", required=True)
    p.add_argument("-o", "--output", required=True)
    a = p.parse_args(argv)

    atus = atu_table(a.links, a.tables)
    rows = []
    cohorts = {"bvbrc": load_records(a.bvbrc_re), "ncbi": load_records(a.ncbi_re)}
    for cohort, recs in cohorts.items():
        versions = sorted({v for r in recs.values() for v in r["labels"]}, key=_vkey)
        va, vb = versions[0], versions[-1]
        pairs = list(changed_pairs(recs, va, vb))
        n_all, k_all = len(pairs), sum(ch for _, _, ch in pairs)

        # dilution concentration
        mics = defaultdict(Counter)
        for _, r, ch in pairs:
            if ch:
                mics[r["ab"]][r["mic"]] += 1
        for ab in sorted(mics):
            c = mics[ab]
            rows.append(["dilution", cohort, ab, "distinct_mic_values_among_changes", len(c), sum(c.values()), "", "", ""])
            rows.append(rate_row("dilution", cohort, ab, "share_at_most_common_mic",
                                 c.most_common(1)[0][1], sum(c.values())))

        # bracketed and ATU
        br = [(r, ch) for _, r, ch in pairs
              if any(bracketed(r[x].get(v, "")) for x in ("S", "R") for v in (va, vb))]
        k_br = sum(ch for _, ch in br)
        rows.append(["bracketed", cohort, "ALL", "changes_with_bracketed_breakpoint", k_br, k_all, "", "", ""])
        rows.append(rate_row("bracketed", cohort, "ALL", "rate_without_bracketed_records", k_all - k_br, n_all - len(br)))
        rows.append(["bracketed", cohort, "ALL", "antibiotics_with_bracketed_rows",
                     len({r["ab"] for r, _ in br}), "", "", "", ""])
        at = [(r, ch) for _, r, ch in pairs
              if any(in_atu(r["mic"], atus.get((r["ab"], v))) for v in (va, vb))]
        k_at = sum(ch for _, ch in at)
        rows.append(["atu", cohort, "ALL", "changes_with_mic_in_atu", k_at, k_all, "", "", ""])
        rows.append(rate_row("atu", cohort, "ALL", "rate_without_atu_records", k_all - k_at, n_all - len(at)))
        for ab, c in sorted(Counter(r["ab"] for r, ch in at if ch).items()):
            rows.append(["atu", cohort, ab, "changes_with_mic_in_atu", c, sum(mics[ab].values()), "", "", ""])
        both = [(r, ch) for _, r, ch in pairs
                if not any(bracketed(r[x].get(v, "")) for x in ("S", "R") for v in (va, vb))
                and not any(in_atu(r["mic"], atus.get((r["ab"], v))) for v in (va, vb))]
        rows.append(rate_row("bracketed_and_atu", cohort, "ALL", "rate_without_either",
                             sum(ch for _, ch in both), len(both)))

        # censored values a read-at-the-bound rule would categorise
        for v in versions:
            und = sum(1 for r in recs.values() if r["labels"].get(v) == "UNDECIDABLE")
            lab = sum(1 for r in recs.values() if r["labels"].get(v) in DECIDABLE | {"UNDECIDABLE"})
            rows.append(rate_row("censored_bound", cohort, "ALL", f"undecidable_labels_v{v}", und, lab))

    # validation sources: BV-BRC calls stated as EUCAST 2020
    src = {}
    with open(a.bvbrc_export, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            src[str(reader.line_num)] = (";".join(sorted(x.strip() for x in row.get("PubMed", "").split(";") if x.strip())) or "none",
                                         bool((row.get("Computational Method") or "").strip()))
    eu20 = [k for k, r in cohorts["bvbrc"].items()
            if r["standard"].upper() == "EUCAST" and r["year"] == "2020"]
    by_src = Counter(src[k][0] for k in eu20)
    for s, n in by_src.most_common():
        rows.append(["validation", "bvbrc", "ALL", f"eucast2020_records_source_{s}", n, len(eu20), "", "", ""])
    rows.append(["validation", "bvbrc", "ALL", "eucast2020_records_with_computational_method",
                 sum(src[k][1] for k in eu20), len(eu20), "", "", ""])
    rows.append(["validation", "bvbrc", "ALL", "eucast2020_genomes",
                 len({cohorts["bvbrc"][k]["genome"] for k in eu20}), len(eu20), "", "", ""])

    # NCBI: human hosts only
    with open(a.ncbi_genomes, newline="", encoding="utf-8") as fh:
        host = {r["Genome ID"]: r["Host"] for r in csv.DictReader(fh)}
    human = lambda h: h.lower().startswith("homo sapiens")  # noqa: E731
    recs = cohorts["ncbi"]
    versions = sorted({v for r in recs.values() for v in r["labels"]}, key=_vkey)
    hp = [(r, ch) for _, r, ch in changed_pairs(recs, versions[0], versions[-1]) if human(host.get(r["genome"], ""))]
    rows.append(rate_row("ncbi_host", "ncbi", "ALL", "rate_human_hosts_only", sum(ch for _, ch in hp), len(hp)))
    rows.append(["ncbi_host", "ncbi", "ALL", "biosamples_human", sum(human(h) for h in host.values()), len(host), "", "", ""])
    rows.append(["ncbi_host", "ncbi", "ALL", "biosamples_host_not_stated", sum(not h for h in host.values()), len(host), "", "", ""])

    # NCBI phenotypes outside S/I/R
    with open(a.ncbi_export, newline="", encoding="utf-8") as fh:
        ph = Counter((r["Resistant Phenotype"] or "not stated") for r in csv.DictReader(fh)
                     if (r["Measurement Unit"] or "").lower() == "mg/l")
    total = sum(ph.values())
    for k, n in ph.most_common():
        if k not in ("Susceptible", "Intermediate", "Resistant"):
            rows.append(["ncbi_phenotype", "ncbi", "ALL", f"reported_{k.replace(' ', '_')}", n, total, "", "", ""])

    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    with open(a.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(COLUMNS)
        w.writerows(rows)
    print(f"wrote {a.output}: {len(rows)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
