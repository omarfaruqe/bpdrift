"""Re-interpret source MIC records under several pinned breakpoint tables.

One output row per (record, table version). Records whose antibiotic has
no entry in the links file are counted, not interpreted. A link with no
EUCAST agent gives NO_BREAKPOINT for that version, with the link's note.

Each row also states whether it belongs to the main analysis and, if not,
why. The reasons, checked in this order:

    conflicting_values   the genome has more than one value for the drug
    ratio_method         a ratio-pair MIC for an inhibitor combination whose
                         EUCAST breakpoint assumes a fixed inhibitor dose
    urinary_only         the only applicable breakpoint is urinary-specific

Sign conflicts stay in the main analysis (the sign field is trusted) and
are flagged for a sensitivity analysis.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable

from bpdrift.bvbrc import ASTRecord
from bpdrift.diff import norm
from bpdrift.eucast import EucastTable
from bpdrift.interpret import Category, interpret
from bpdrift.links import LinkError, Links, resolve

# EUCAST tests these with a fixed inhibitor concentration; a value written as
# a ratio pair ("32/16") comes from a 2:1 ratio method and is not comparable.
FIXED_INHIBITOR = frozenset({"amoxicillin/clavulanic acid", "ampicillin/sulbactam"})

ORDER = (Category.S, Category.I, Category.R)


@dataclass(frozen=True)
class Reinterpretation:
    record: ASTRecord
    version: str
    agent: str
    cell: str
    s: str
    r: str
    restriction: str
    category: Category
    possible: str          # "S|I" for an undecidable value, the category otherwise
    reason: str
    exclusion: str         # "" when in the main analysis

    @property
    def in_main(self) -> bool:
        return self.exclusion == ""


def exclusion_reason(record: ASTRecord, restriction: str) -> str:
    if record.conflicting_values:
        return "conflicting_values"
    if record.ratio_format and norm(record.antibiotic) in FIXED_INHIBITOR:
        return "ratio_method"
    if restriction == "urinary":
        return "urinary_only"
    return ""


def reinterpret(
    records: Iterable[ASTRecord],
    tables: dict[str, EucastTable],
    links: Links,
    group: str,
) -> tuple[list[Reinterpretation], Counter]:
    """Interpret every linked record under every table.

    Returns (rows, unlinked) where unlinked counts records per antibiotic
    that the links file does not cover.
    """
    rows_by_link = resolve(links, tables, group)
    versions = sorted(tables, key=lambda v: tuple(int(p) for p in v.split(".")))
    out: list[Reinterpretation] = []
    unlinked: Counter = Counter()

    for rec in records:
        if not any(links.get(rec.antibiotic, v) for v in versions):
            unlinked[rec.antibiotic] += 1
            continue
        for v in versions:
            link = links.get(rec.antibiotic, v)
            if link is None:
                raise LinkError(f"{rec.antibiotic!r} is linked for some versions but not v{v}")
            row = rows_by_link[(norm(rec.antibiotic), v)]
            if row is None:
                out.append(Reinterpretation(
                    rec, v, "", "", "", "", link.restriction, Category.NO_BREAKPOINT,
                    "", link.note or "no applicable row", exclusion_reason(rec, link.restriction),
                ))
                continue
            result = interpret(rec.mic, row.breakpoint())
            possible = "|".join(c.value for c in ORDER if c in result.possible)
            out.append(Reinterpretation(
                rec, v, row.agent, row.cell, row.s.raw or row.s.kind, row.r.raw or row.r.kind,
                link.restriction, result.category, possible, result.reason,
                exclusion_reason(rec, link.restriction),
            ))
    return out, unlinked


COLUMNS = [
    "genome_id", "antibiotic", "mic", "source_line",
    "reported_phenotype", "reported_standard", "reported_standard_year", "method",
    "eucast_version", "eucast_agent", "eucast_cell", "eucast_S", "eucast_R",
    "category", "possible", "reason",
    "restriction", "exclusion", "sign_conflict", "ratio_format",
]


def as_row(x: Reinterpretation) -> list:
    r = x.record
    return [
        r.genome_id, r.antibiotic, str(r.mic), r.line,
        r.reported_phenotype, r.standard, r.standard_year, r.method,
        x.version, x.agent, x.cell, x.s, x.r,
        x.category.value, x.possible, x.reason,
        x.restriction, x.exclusion, int(r.sign_conflict), int(r.ratio_format),
    ]
