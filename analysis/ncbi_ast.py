"""Replication cohort: NCBI BioSample antibiograms -> BV-BRC-format AST records.

    python analysis/ncbi_ast.py data/ncbi/ast \
        --exclude-bvbrc data/bvbrc/kpneumoniae_amr_all.csv data/bvbrc/kpneumoniae_genomes.csv \
        -o data/derived/ncbi_ast_bvbrc_format.csv --genomes-out data/derived/ncbi_ast_genomes.csv \
        --summary results/ncbi_cohort.txt

Reads the BioSample XML batches written by scripts/download_ncbi_ast.sh and
writes every antibiogram row in the column layout of a BV-BRC genome_amr
export, so the replication runs through exactly the same reader, antibiotic
links and analysis scripts as the main cohort. The BioSample accession takes
the place of the BV-BRC genome ID.

Independence: BioSamples that already have AST records in BV-BRC (matched by
BioSample accession through the BV-BRC genome table) are left out, so no
isolate contributes to both cohorts.

Mapping: antibiotic names use BV-BRC spelling for combinations
("piperacillin-tazobactam" -> "piperacillin/tazobactam"); measurement sign
"==" becomes "="; phenotypes are capitalised; standard "missing" becomes
empty. NCBI antibiograms state no testing-standard year, so that column is
empty.
"""

from __future__ import annotations

import argparse
import csv
import glob
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

BVBRC_COLUMNS = [
    "Taxon ID", "Genome ID", "Genome Name", "Antibiotic", "Resistant Phenotype", "Measurement",
    "Measurement Sign", "Measurement Value", "Measurement Unit", "Laboratory Typing Method",
    "Laboratory Typing Method Version", "Laboratory Typing Platform", "Vendor", "Testing Standard",
    "Testing Standard Year", "Computational Method", "Computational Method Version",
    "Computational Method Performance", "Evidence", "Source", "PubMed",
]
MISSING = {"", "missing", "not applicable", "not collected", "not provided", "na", "n/a"}
SIGNS = {"==": "=", "=": "=", "<": "<", "<=": "<=", ">": ">", ">=": ">="}


def _clean(text: str | None) -> str:
    t = (text or "").strip()
    return "" if t.lower() in MISSING else t


def antibiotic_name(ncbi: str) -> str:
    """NCBI 'piperacillin-tazobactam' -> BV-BRC 'piperacillin/tazobactam'."""
    return "/".join(p.strip() for p in ncbi.strip().lower().split("-"))


def parse_biosamples(paths) -> tuple[list[dict], dict[str, dict]]:
    """Antibiogram rows (one dict per row, with BioSample fields) and per-sample attributes."""
    rows, samples = [], {}
    for path in paths:
        for bs in ET.parse(path).getroot().iter("BioSample"):
            acc = bs.get("accession")
            org = bs.find("Description/Organism")
            attrs = {(a.get("harmonized_name") or a.get("attribute_name")): (a.text or "").strip()
                     for a in bs.iter("Attribute")}
            samples[acc] = {
                "organism": (org.get("taxonomy_name") if org is not None else "") or "",
                "taxonomy_id": (org.get("taxonomy_id") if org is not None else "") or "",
                "host": _clean(attrs.get("host")),
            }
            for table in bs.iter("Table"):
                if not (table.get("class") or "").startswith("Antibiogram"):
                    continue
                header = [(c.text or "").strip() for c in table.find("Header")]
                for r in table.find("Body"):
                    values = [(c.text or "").strip() for c in r]
                    rows.append({"biosample": acc, **dict(zip(header, values))})
    return rows, samples


def bvbrc_biosamples(amr_csv, genomes_csv) -> set[str]:
    """BioSamples of genomes that have AST records in a BV-BRC export."""
    with open(genomes_csv, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh)
        canon = {f.lower().replace(" ", "_"): f for f in r.fieldnames}
        gid, bs = canon["genome_id"], canon["biosample_accession"]
        to_bs = {row[gid].strip(): row[bs].strip() for row in r}
    with open(amr_csv, newline="", encoding="utf-8") as fh:
        genomes = {row["Genome ID"].strip() for row in csv.DictReader(fh)}
    return {to_bs[g] for g in genomes if to_bs.get(g)}


def to_bvbrc(row: dict, sample: dict) -> list[str]:
    sign = SIGNS.get(row.get("Measurement sign", "").strip(), "")
    value = _clean(row.get("Measurement"))
    shown = value if sign in ("", "=") else f"{sign}{value}"
    phen = _clean(row.get("Resistance phenotype"))
    std = _clean(row.get("Testing standard"))
    d = dict.fromkeys(BVBRC_COLUMNS, "")
    d.update({
        "Taxon ID": sample["taxonomy_id"], "Genome ID": row["biosample"], "Genome Name": sample["organism"],
        "Antibiotic": antibiotic_name(row.get("Antibiotic", "")),
        "Resistant Phenotype": phen[:1].upper() + phen[1:] if phen and phen != "not defined" else "",
        "Measurement": shown if value else "", "Measurement Sign": sign if value else "",
        "Measurement Value": value, "Measurement Unit": _clean(row.get("Measurement units")),
        "Laboratory Typing Method": _clean(row.get("Laboratory typing method")),
        "Laboratory Typing Method Version": _clean(row.get("Laboratory typing method version or reagent")),
        "Laboratory Typing Platform": _clean(row.get("Laboratory typing platform")),
        "Vendor": _clean(row.get("Vendor")), "Testing Standard": std,
        "Evidence": "Laboratory Method", "Source": "NCBI BioSample antibiogram",
    })
    return [d[c] for c in BVBRC_COLUMNS]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("xml_dir")
    p.add_argument("--exclude-bvbrc", nargs=2, metavar=("AMR_CSV", "GENOMES_CSV"))
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--genomes-out", required=True)
    p.add_argument("--summary", required=True)
    a = p.parse_args(argv)

    paths = sorted(glob.glob(str(Path(a.xml_dir) / "biosample_*.xml")))
    if not paths:
        raise SystemExit(f"no biosample_*.xml in {a.xml_dir}")
    rows, samples = parse_biosamples(paths)
    excluded = bvbrc_biosamples(*a.exclude_bvbrc) if a.exclude_bvbrc else set()
    kept = [r for r in rows if r["biosample"] not in excluded]
    kept_samples = sorted({r["biosample"] for r in kept})
    mic = [r for r in kept if (r.get("Measurement units") or "").strip() == "mg/L"]

    for path in (a.output, a.genomes_out, a.summary):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(a.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(BVBRC_COLUMNS)
        w.writerows(to_bvbrc(r, samples[r["biosample"]]) for r in kept)
    with open(a.genomes_out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["Genome ID", "BioSample Accession", "Host"])
        w.writerows([s, s, samples[s]["host"]] for s in kept_samples)
    hosts = Counter(samples[s]["host"] or "not stated" for s in kept_samples)
    with open(a.summary, "w", encoding="utf-8") as fh:
        fh.write(f"BioSamples downloaded: {len(samples)}\n")
        fh.write(f"antibiogram rows: {len(rows)}\n")
        fh.write(f"BioSamples left out, already with AST in BV-BRC: {len(set(samples) & excluded)}\n")
        fh.write(f"BioSamples kept: {len(kept_samples)}\n")
        fh.write(f"rows kept: {len(kept)} (MIC in mg/L: {len(mic)}; BioSamples with an MIC: "
                 f"{len({r['biosample'] for r in mic})})\n")
        fh.write("testing standard of kept MIC rows: " + ", ".join(
            f"{k or 'not stated'} {n}" for k, n in Counter(_clean(r.get('Testing standard')) for r in mic).most_common()) + "\n")
        fh.write("host of kept BioSamples: " + ", ".join(f"{k} {n}" for k, n in hosts.most_common(6)) + "\n")
    print(f"wrote {a.output}, {a.genomes_out} and {a.summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
