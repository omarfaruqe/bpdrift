"""Read MIC breakpoints from an official EUCAST breakpoint table (.xlsx).

The tables carry no public licence, so bpdrift never ships them. The user
downloads the file from eucast.org; this module reads it.

Layout, as found in v16.1 (2026), one worksheet per organism group:

    row n    <drug class> | "MIC breakpoints (mg/L)" | ...
    row n+1               | "S ≤" | "R >" | "ATU"
    rows     <agent>      |  S    |  R    |  ATU     (until column A is blank)

Some sheets (the anaerobes) hold one block per organism subgroup, each
under its own large bold heading ("Bacteroides spp.", "Prevotella spp."),
and repeat the same agents in every block. The heading is kept as the
row's subgroup, so those rows stay distinct.

Footnote markers are superscript runs inside the cell's rich text, so
"8" + superscript "2" is the value 8 with note 2, never 82. Reading the
plain cell value would silently corrupt such breakpoints; this module
splits the runs instead.

Every S and R cell must be a number, a bracketed number, or one of the
codes in CODES. Anything else raises EucastFormatError naming the cell.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from bpdrift.interpret import Breakpoint
from bpdrift.pin import Pin, verify

# Codes EUCAST uses in place of a breakpoint.
CODES = {
    "-": "none",   # susceptibility testing not recommended
    "IE": "IE",    # insufficient evidence
    "Note": "note",  # see the numbered or lettered note
    "NA": "NA",    # not applicable
    "IP": "IP",    # in preparation
}

MIC_HEADER = "MIC breakpoint"
# Sheets that hold no organism breakpoint table.
NON_ORGANISM_SHEETS = {
    "Content", "Changes", "Notes", "Guidance", "Dosages",
    "Technical uncertainty", "PK PD cutoffs",
}
PKPD_SHEET = "PK PD breakpoints"

_NUMBER = re.compile(r"\d*\.?\d+")
_VERSION = re.compile(r"EUCAST Clinical Breakpoint Tables v\.\s*(\d+\.\d+)")


class EucastFormatError(ValueError):
    """Raised when the workbook does not look like a EUCAST breakpoint table."""


@dataclass(frozen=True)
class BPValue:
    """One S or R cell."""

    kind: str                  # "value", "none", "IE", "note", "NA", "IP", "empty"
    value: float | None = None
    bracketed: bool = False    # EUCAST brackets: use with the caveat in the notes
    notes: tuple[str, ...] = ()
    raw: str = ""


@dataclass(frozen=True)
class MICRow:
    sheet: str
    drug_class: str
    agent: str                 # agent text as written, footnote markers removed
    agent_notes: tuple[str, ...]
    s: BPValue
    r: BPValue
    atu: str                   # area of technical uncertainty, as written
    row: int                   # 1-based worksheet row
    subgroup: str = ""         # organism subgroup heading, "" when the sheet has none

    @property
    def cell(self) -> str:
        return f"{self.sheet}!A{self.row}"

    @property
    def group(self) -> str:
        """The organism group the breakpoint applies to."""
        return self.subgroup or self.sheet

    def breakpoint(self) -> Breakpoint:
        """The pair for interpret(); incomplete unless both cells hold values."""
        if self.s.kind == "value" and self.r.kind == "value":
            return Breakpoint(self.s.value, self.r.value)
        return Breakpoint(None, None)


@dataclass
class EucastTable:
    version: str
    pin: Pin | None
    rows: list[MICRow] = field(default_factory=list)

    @property
    def sheets(self) -> list[str]:
        return list(dict.fromkeys(r.sheet for r in self.rows))

    def for_sheet(self, sheet: str) -> list[MICRow]:
        return [r for r in self.rows if r.sheet == sheet]


# ----------------------------------------------------------------------
# rich text helpers


def _split(value) -> tuple[str, tuple[str, ...]]:
    """Return (plain text, superscript notes) for a cell value."""
    from openpyxl.cell.rich_text import CellRichText, TextBlock

    if value is None:
        return "", ()
    if not isinstance(value, CellRichText):
        if isinstance(value, float) and value.is_integer():
            value = int(value)
        return str(value), ()

    text, notes = [], []
    for block in value:
        if isinstance(block, TextBlock):
            font = block.font
            if font is not None and font.vertAlign == "superscript":
                notes.extend(n for n in re.split(r"[,\s]+", block.text) if n)
                continue
            text.append(block.text)
        else:
            text.append(str(block))
    return "".join(text), tuple(notes)


def _parse_value(cell_value, ref: str) -> BPValue:
    text, notes = _split(cell_value)
    raw = text.strip()
    if raw == "":
        return BPValue("empty", notes=notes, raw=raw)
    if raw in CODES:
        return BPValue(CODES[raw], notes=notes, raw=raw)

    inner, bracketed = raw, False
    if raw.startswith("(") and raw.endswith(")"):
        inner, bracketed = raw[1:-1].strip(), True
    if _NUMBER.fullmatch(inner):
        value = float(inner)
        if value <= 0:
            raise EucastFormatError(f"{ref}: breakpoint must be positive, got {raw!r}")
        return BPValue("value", value, bracketed, notes, raw)

    raise EucastFormatError(f"{ref}: unrecognised breakpoint cell {raw!r}")


# ----------------------------------------------------------------------
# reading


def _detect_version(workbook) -> str:
    """Most common version string in the organism sheets' title cells."""
    found = Counter()
    for ws in workbook.worksheets:
        if ws.title in NON_ORGANISM_SHEETS or ws.title == PKPD_SHEET:
            continue
        for row in ws.iter_rows(min_row=1, max_row=3):
            for c in row:
                m = _VERSION.search(_split(c.value)[0])
                if m:
                    found[m.group(1)] += 1
    if not found:
        raise EucastFormatError("no 'EUCAST Clinical Breakpoint Tables v. N' title found")
    return found.most_common(1)[0][0]


SUBGROUP_MIN_SIZE = 12  # subgroup headings are 14 pt bold; notes are 8 pt


def _subgroup(ws, hdr_row: int, floor: int) -> str:
    """Nearest large bold heading in column A between floor and hdr_row.

    Row 1 is the sheet title and is never a subgroup.
    """
    for r in range(hdr_row - 1, max(floor, 1), -1):
        cell = ws.cell(r, 1)
        text = " ".join(_split(cell.value)[0].split())
        if not text or text.startswith("http"):
            continue
        font = cell.font
        if font is not None and font.b and (font.sz or 0) >= SUBGROUP_MIN_SIZE:
            return text
    return ""


def _read_sheet(ws) -> list[MICRow]:
    rows: list[MICRow] = []
    floor = 1  # last row of the previous block
    for header in ws["B"]:
        if not _split(header.value)[0].startswith(MIC_HEADER):
            continue
        hdr_row = header.row
        subgroup = _subgroup(ws, hdr_row, floor)
        s_label = _split(ws.cell(hdr_row + 1, 2).value)[0].strip()
        r_label = _split(ws.cell(hdr_row + 1, 3).value)[0].strip()
        if not s_label.endswith("≤") or r_label != "R >":
            raise EucastFormatError(
                f"{ws.title}!B{hdr_row + 1}: expected 'S ≤' / 'R >' under the MIC header, "
                f"got {s_label!r} / {r_label!r}"
            )
        drug_class = _split(ws.cell(hdr_row, 1).value)[0].strip()

        r = hdr_row + 2
        while r <= ws.max_row:
            agent, agent_notes = _split(ws.cell(r, 1).value)
            agent = " ".join(agent.split())
            if not agent:
                break
            if _split(ws.cell(r, 2).value)[0].startswith(MIC_HEADER):
                break
            rows.append(
                MICRow(
                    sheet=ws.title,
                    drug_class=drug_class,
                    agent=agent,
                    agent_notes=agent_notes,
                    s=_parse_value(ws.cell(r, 2).value, f"{ws.title}!B{r}"),
                    r=_parse_value(ws.cell(r, 3).value, f"{ws.title}!C{r}"),
                    atu=_split(ws.cell(r, 4).value)[0].strip(),
                    row=r,
                    subgroup=subgroup,
                )
            )
            r += 1
        floor = r
    return rows


def read_eucast_xlsx(
    path: str | Path,
    pin: Pin | None = None,
    include_pkpd: bool = False,
) -> EucastTable:
    """Read every organism sheet's MIC breakpoints.

    If a pin is given, the file must match its checksum and the version
    written in the workbook must match the pinned version.
    """
    import warnings

    import openpyxl

    if pin is not None:
        verify(path, pin)

    with warnings.catch_warnings():
        # openpyxl warns about drawing shapes it cannot keep; irrelevant here.
        warnings.simplefilter("ignore", UserWarning)
        wb = openpyxl.load_workbook(path, rich_text=True)

    version = _detect_version(wb)
    if pin is not None and pin.version != version:
        raise EucastFormatError(
            f"pinned as v{pin.version} but the workbook says v{version}"
        )

    table = EucastTable(version=version, pin=pin)
    for ws in wb.worksheets:
        if ws.title in NON_ORGANISM_SHEETS:
            continue
        if ws.title == PKPD_SHEET and not include_pkpd:
            continue
        table.rows.extend(_read_sheet(ws))

    if not table.rows:
        raise EucastFormatError("no MIC breakpoint blocks found")
    return table
