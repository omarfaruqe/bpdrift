"""Link source antibiotic names to the EUCAST row used in each table version.

A links file is a CSV, '#' lines being comments:

    bvbrc_antibiotic,eucast_version,eucast_agent,restriction,note
    meropenem,10.0,Meropenem,,
    meropenem,16.1,Meropenem (indications other than meningitis),,indication split
    doripenem,10.0,,,agent not in this version

An empty eucast_agent records, on purpose, that no row applies in that
version. restriction is "" or "urinary" (the only applicable breakpoint is
urinary-tract specific), so an analysis can include or exclude those.

Nothing is guessed: an antibiotic or version missing from the file is
reported as unlinked, never matched by name.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from bpdrift.diff import norm
from bpdrift.eucast import EucastTable, MICRow

LINK_COLUMNS = ["bvbrc_antibiotic", "eucast_version", "eucast_agent", "restriction", "note"]
RESTRICTIONS = {"", "urinary"}


class LinkError(ValueError):
    """Raised for a malformed links file or a link to a row that does not exist."""


@dataclass(frozen=True)
class Link:
    antibiotic: str
    version: str
    agent: str          # "" when no row applies
    restriction: str
    note: str
    line: int


@dataclass
class Links:
    entries: dict[tuple[str, str], Link] = field(default_factory=dict)

    def get(self, antibiotic: str, version: str) -> Link | None:
        return self.entries.get((norm(antibiotic), version))

    @property
    def antibiotics(self) -> set[str]:
        return {link.antibiotic for link in self.entries.values()}

    @property
    def versions(self) -> set[str]:
        return {link.version for link in self.entries.values()}


def load_links(path: str | Path) -> Links:
    with open(path, newline="", encoding="utf-8") as fh:
        numbered = [(n, t) for n, t in enumerate(fh, start=1)
                    if t.strip() and not t.lstrip().startswith("#")]
    if not numbered:
        raise LinkError(f"{path}: empty links file")
    header = [c.strip() for c in next(csv.reader([numbered[0][1]]))]
    if header != LINK_COLUMNS:
        raise LinkError(f"{path}: header must be {','.join(LINK_COLUMNS)}")

    links = Links()
    for line, text in numbered[1:]:
        cells = next(csv.reader([text]))
        if len(cells) != len(header):
            raise LinkError(f"{path} line {line}: {len(cells)} fields, expected {len(header)}")
        ab, version, agent, restriction, note = (" ".join(c.split()) for c in cells)
        if not ab or not version:
            raise LinkError(f"{path} line {line}: antibiotic and version are required")
        if restriction not in RESTRICTIONS:
            raise LinkError(f"{path} line {line}: restriction must be one of {sorted(RESTRICTIONS)}")
        if not agent and restriction:
            raise LinkError(f"{path} line {line}: a restriction needs an agent")
        key = (norm(ab), version)
        if key in links.entries:
            raise LinkError(
                f"{path} line {line}: {ab!r} v{version} already linked on line {links.entries[key].line}"
            )
        links.entries[key] = Link(ab, version, agent, restriction, note, line)
    return links


def resolve(
    links: Links, tables: dict[str, EucastTable], group: str | Iterable[str]
) -> dict[tuple[str, str], MICRow | None]:
    """Find the table row for every link; fail on any link that does not resolve.

    group may be several names for the same organism group, for a sheet
    EUCAST renamed between versions (Enterobacteriaceae -> Enterobacterales).

    Returns {(normalised antibiotic, version): row or None}.
    """
    given = [group] if isinstance(group, str) else list(group)
    names = {norm(g) for g in given}
    label = " / ".join(given)
    found: dict[tuple[str, str], MICRow | None] = {}
    for key, link in links.entries.items():
        table = tables.get(link.version)
        if table is None:
            raise LinkError(f"line {link.line}: no table loaded for v{link.version}")
        if not link.agent:
            found[key] = None
            continue
        rows = [r for r in table.rows
                if names & {norm(r.group), norm(r.sheet)} and norm(r.agent) == norm(link.agent)]
        if not rows:
            raise LinkError(f"line {link.line}: no {link.agent!r} in {label} of v{link.version}")
        if len(rows) > 1:
            raise LinkError(f"line {link.line}: {link.agent!r} is ambiguous in {label} of v{link.version}")
        found[key] = rows[0]
    return found
