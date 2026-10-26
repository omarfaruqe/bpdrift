"""Read laboratory AST records exported from BV-BRC (genome_amr, CSV).

Export, for example for Klebsiella pneumoniae (taxon 573):

    https://www.bv-brc.org/api/genome_amr/?and(eq(taxon_id,573),
        eq(evidence,Laboratory%20Method))&sort(+id)&limit(25000,OFFSET)
        &http_accept=text/csv

Every record ends up either as an ASTRecord or as a counted exclusion
with a reason, so the import can be reported in full:

    no_value         no measurement, only an S/I/R call
    not_mic          a zone diameter (mm) or another non-MIC unit
    not_laboratory   a computational prediction
    unparseable      a value that is not a number
    duplicate        an exact repeat of an earlier line

Rules, each recorded on the record so it can be audited:

* The sign comes from the Measurement Sign field. When that field is
  empty, from the sign written in Measurement. When the two disagree
  (Measurement "32", sign ">"), the sign field wins and sign_conflict is
  set: the sign decides whether the value is censored, so a sensitivity
  analysis may want to treat these as exact.
* A combination value such as "16/304" keeps its first component, the
  concentration EUCAST breakpoints are stated in. ratio_format marks
  values written as a pair, which for amoxicillin-clavulanic acid means a
  2:1 ratio method (CLSI) rather than EUCAST's fixed clavulanate.
* Several values for the same genome and antibiotic are all kept;
  conflicting_values marks them.
"""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path

from bpdrift.mic import MIC, MICParseError

REQUIRED = [
    "Genome ID", "Genome Name", "Antibiotic", "Resistant Phenotype", "Measurement",
    "Measurement Sign", "Measurement Value", "Measurement Unit",
    "Laboratory Typing Method", "Testing Standard", "Testing Standard Year", "Evidence",
]
MIC_UNITS = {"mg/l", ""}  # an empty unit with a value is taken as mg/L; counted
_SIGN = re.compile(r"^\s*(<=|>=|=<|=>|<|>|=)")
_SIGN_ALIASES = {"=<": "<=", "=>": ">="}

EXCLUSION_REASONS = ("no_value", "not_mic", "not_laboratory", "unparseable", "duplicate")


class BVBRCFormatError(ValueError):
    """Raised when the file is not a BV-BRC genome_amr CSV export."""


@dataclass(frozen=True)
class ASTRecord:
    line: int                 # line number in the export
    genome_id: str
    genome_name: str
    antibiotic: str           # as BV-BRC writes it, e.g. "piperacillin/tazobactam"
    mic: MIC
    reported_phenotype: str   # "Susceptible", "Resistant", "Intermediate", ...
    method: str
    platform: str
    vendor: str
    standard: str             # "CLSI", "EUCAST", "EUCAST,CLSI" or ""
    standard_year: str
    pmid: str
    raw_measurement: str
    raw_sign: str
    raw_value: str
    sign_conflict: bool = False
    ratio_format: bool = False
    unit_assumed: bool = False
    conflicting_values: bool = False


@dataclass
class Import:
    records: list[ASTRecord] = field(default_factory=list)
    excluded: Counter = field(default_factory=Counter)
    total_lines: int = 0

    def summary(self) -> dict[str, int]:
        return {
            "lines": self.total_lines,
            "records": len(self.records),
            **{f"excluded_{k}": self.excluded.get(k, 0) for k in EXCLUSION_REASONS},
            "sign_conflict": sum(r.sign_conflict for r in self.records),
            "ratio_format": sum(r.ratio_format for r in self.records),
            "unit_assumed": sum(r.unit_assumed for r in self.records),
            "conflicting_values": sum(r.conflicting_values for r in self.records),
            "genomes": len({r.genome_id for r in self.records}),
        }


def _sign_from_text(text: str) -> str:
    m = _SIGN.match(text or "")
    if not m:
        return ""
    return _SIGN_ALIASES.get(m.group(1), m.group(1))


def _parse_value(text: str) -> tuple[float, bool]:
    """Number in a Measurement Value; first component of 'a/b'."""
    t = text.strip()
    t = _SIGN.sub("", t).strip()
    ratio = "/" in t
    first = t.split("/", 1)[0].strip()
    if not re.fullmatch(r"\d*\.?\d+", first):
        raise MICParseError(f"not a number: {text!r}")
    return float(first), ratio


def read_bvbrc_csv(path: str | Path) -> Import:
    result = Import()
    seen: set[tuple] = set()
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            raise BVBRCFormatError(f"{path}: missing columns {missing}")
        for line, row in enumerate(reader, start=2):
            result.total_lines += 1
            key = tuple(row.values())
            if key in seen:
                result.excluded["duplicate"] += 1
                continue
            seen.add(key)

            if row["Evidence"].strip() != "Laboratory Method":
                result.excluded["not_laboratory"] += 1
                continue
            raw_value = row["Measurement Value"].strip()
            if not raw_value:
                result.excluded["no_value"] += 1
                continue
            unit = row["Measurement Unit"].strip()
            if unit.lower() not in MIC_UNITS:
                result.excluded["not_mic"] += 1
                continue
            try:
                value, ratio = _parse_value(raw_value)
            except MICParseError:
                result.excluded["unparseable"] += 1
                continue
            if value <= 0:
                result.excluded["unparseable"] += 1
                continue

            field_sign = _SIGN_ALIASES.get(row["Measurement Sign"].strip(), row["Measurement Sign"].strip())
            text_sign = _sign_from_text(row["Measurement"]) or _sign_from_text(raw_value)
            sign = field_sign or text_sign or "="
            # The sign field is trusted; flag it whenever the reading says otherwise.
            conflict = bool(field_sign) and field_sign != (text_sign or "=")

            result.records.append(ASTRecord(
                line=line,
                genome_id=row["Genome ID"].strip(),
                genome_name=row["Genome Name"].strip(),
                antibiotic=row["Antibiotic"].strip(),
                mic=MIC(value, sign),
                reported_phenotype=row["Resistant Phenotype"].strip(),
                method=row["Laboratory Typing Method"].strip(),
                platform=(row.get("Laboratory Typing Platform") or "").strip(),
                vendor=(row.get("Vendor") or "").strip(),
                standard=row["Testing Standard"].strip(),
                standard_year=row["Testing Standard Year"].strip(),
                pmid=(row.get("PubMed") or "").strip(),
                raw_measurement=row["Measurement"].strip(),
                raw_sign=row["Measurement Sign"].strip(),
                raw_value=raw_value,
                sign_conflict=conflict,
                ratio_format=ratio,
                unit_assumed=(unit == ""),
            ))

    # several distinct values for one genome and antibiotic
    values = defaultdict(set)
    for r in result.records:
        values[(r.genome_id, r.antibiotic.casefold())].add((r.mic.sign, r.mic.value))
    clash = {k for k, v in values.items() if len(v) > 1}
    if clash:
        result.records = [
            replace(r, conflicting_values=True) if (r.genome_id, r.antibiotic.casefold()) in clash else r
            for r in result.records
        ]
    return result
