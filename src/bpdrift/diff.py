"""Compare the MIC breakpoints of two EUCAST tables.

Rows are matched on (organism group, agent). The organism group is the
subgroup heading when a sheet has one (the anaerobes), otherwise the sheet.
Agent names are compared ignoring case, spacing and dash style, and
nothing else: there is no fuzzy matching.

What exact matching cannot pair is reported as REMOVED or ADDED. Rows
that probably correspond (same drug, different qualifier, such as
"Ciprofloxacin" splitting into "(meningitis)" and "(indications other
than meningitis)") are listed as candidates but never paired
automatically. To pair them, write them in a mapping file (see
load_mapping); every mapping line must name rows that exist, or the diff
stops with an error.

For each pair whose breakpoints differ, the diff lists the MIC dilutions
whose category changes between the two versions.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Iterable

from bpdrift.eucast import BPValue, EucastTable, MICRow
from bpdrift.interpret import Category, interpret

# Doubling dilutions in mg/L, as reported on MIC panels.
DILUTIONS = (
    0.002, 0.004, 0.008, 0.016, 0.03, 0.06, 0.125, 0.25, 0.5,
    1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
)

MAPPING_COLUMNS = ["old_group", "old_agent", "new_group", "new_agent"]


class DiffError(ValueError):
    """Raised for an ambiguous table or an invalid mapping."""


class ChangeKind(str, Enum):
    UNCHANGED = "unchanged"
    CHANGED = "changed"
    ADDED = "added"
    REMOVED = "removed"


@dataclass(frozen=True)
class Flip:
    """One MIC dilution whose category differs between the versions."""

    mic: float
    old: Category
    new: Category


@dataclass(frozen=True)
class Change:
    kind: ChangeKind
    old: MICRow | None
    new: MICRow | None
    match: str = ""                     # "exact", "mapped" or "" when unpaired
    flips: tuple[Flip, ...] = ()
    candidates: tuple[str, ...] = ()    # likely partners on the other side, unpaired

    @property
    def group(self) -> str:
        return (self.new or self.old).group

    @property
    def agent(self) -> str:
        return (self.new or self.old).agent


@dataclass
class Mapping:
    group_renames: dict[str, str] = field(default_factory=dict)
    pairs: list[tuple[str, str, str, str, int]] = field(default_factory=list)


@dataclass
class Diff:
    old_version: str
    new_version: str
    changes: list[Change]
    old_pin: object = None
    new_pin: object = None

    def of_kind(self, kind: ChangeKind) -> list[Change]:
        return [c for c in self.changes if c.kind is kind]

    def summary(self) -> dict[str, int]:
        counts = Counter(c.kind.value for c in self.changes)
        return {k.value: counts.get(k.value, 0) for k in ChangeKind}


# ----------------------------------------------------------------------
# names


def norm(text: str) -> str:
    """Comparison form of an agent or group name."""
    text = text.replace("–", "-").replace("—", "-")
    return " ".join(text.split()).casefold()


def base_name(agent: str) -> str:
    """Drug name without qualifiers: the text before the first '(' or ','."""
    n = norm(agent)
    for sep in ("(", ","):
        n = n.split(sep, 1)[0]
    return n.strip()


def related(a: str, b: str) -> bool:
    """True when two agent names probably refer to the same drug.

    Used only to suggest candidates, never to pair rows. The base names
    must be equal, or one must extend the other by a whole word
    ("Amoxicillin" and "Amoxicillin iv"). "Amoxicillin" and
    "Amoxicillin-clavulanic acid" are different drugs and are not related.
    """
    x, y = base_name(a), base_name(b)
    if x == y:
        return True
    short, long_ = sorted((x, y), key=len)
    return long_.startswith(short + " ")


# ----------------------------------------------------------------------
# mapping file


def load_mapping(path: str | Path) -> Mapping:
    """Read a CSV with columns old_group, old_agent, new_group, new_agent.

    A line with both agent columns empty renames an organism group, for
    example a sheet EUCAST renamed between versions. Any other line pairs
    one old row with one new row; repeat the old row on several lines to
    record a split. Group names are as they appear in each table.
    """
    mapping = Mapping()
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or [f.strip() for f in reader.fieldnames] != MAPPING_COLUMNS:
            raise DiffError(f"{path}: header must be {','.join(MAPPING_COLUMNS)}")
        for line, raw in enumerate(reader, start=2):
            og, oa, ng, na = (" ".join((raw.get(c) or "").split()) for c in MAPPING_COLUMNS)
            if not og or not ng:
                raise DiffError(f"{path} line {line}: old_group and new_group are required")
            if not oa and not na:
                if norm(og) in {norm(k) for k in mapping.group_renames}:
                    raise DiffError(f"{path} line {line}: group {og!r} renamed twice")
                mapping.group_renames[og] = ng
            elif oa and na:
                mapping.pairs.append((og, oa, ng, na, line))
            else:
                raise DiffError(
                    f"{path} line {line}: give both agents, or neither for a group rename"
                )
    return mapping


# ----------------------------------------------------------------------
# comparison


def _same(a: BPValue, b: BPValue) -> bool:
    return (a.kind, a.value, a.bracketed) == (b.kind, b.value, b.bracketed)


def flips(old: MICRow, new: MICRow) -> tuple[Flip, ...]:
    """Dilutions whose category differs; empty unless both pairs are complete."""
    a, b = old.breakpoint(), new.breakpoint()
    if not (a.is_complete and b.is_complete):
        return ()
    out = []
    for mic in DILUTIONS:
        ca, cb = interpret(mic, a).category, interpret(mic, b).category
        if ca is not cb:
            out.append(Flip(mic, ca, cb))
    return tuple(out)


def _pair(old: MICRow, new: MICRow, match: str) -> Change:
    if _same(old.s, new.s) and _same(old.r, new.r):
        return Change(ChangeKind.UNCHANGED, old, new, match)
    return Change(ChangeKind.CHANGED, old, new, match, flips(old, new))


def _index(rows: Iterable[MICRow], label: str, renames: dict[str, str]) -> dict:
    lowered = {norm(k): v for k, v in renames.items()}
    index = {}
    for r in rows:
        group = lowered.get(norm(r.group), r.group)
        key = (norm(group), norm(r.agent))
        if key in index:
            raise DiffError(
                f"{label} table has {r.agent!r} twice in group {r.group!r}: "
                f"{index[key].cell} and {r.cell}"
            )
        index[key] = r
    return index


def diff_tables(
    old: EucastTable,
    new: EucastTable,
    mapping: Mapping | None = None,
    groups: Iterable[str] | None = None,
) -> Diff:
    """Compare two tables, optionally limited to some organism groups or sheets."""
    mapping = mapping or Mapping()
    wanted = {norm(g) for g in groups} if groups else None

    def keep(r: MICRow) -> bool:
        return wanted is None or norm(r.group) in wanted or norm(r.sheet) in wanted

    old_rows = [r for r in old.rows if keep(r)]
    new_rows = [r for r in new.rows if keep(r)]
    old_idx = _index(old_rows, f"old (v{old.version})", mapping.group_renames)
    new_idx = _index(new_rows, f"new (v{new.version})", {})

    # The mapping is checked against the full tables, so a filter cannot hide a typo.
    full_old = _index(old.rows, f"old (v{old.version})", {})
    full_new = _index(new.rows, f"new (v{new.version})", {})
    for og in mapping.group_renames:
        if not any(k[0] == norm(og) for k in full_old):
            raise DiffError(f"mapping renames group {og!r}, which is not in the old table")

    changes: list[Change] = []
    used_old, used_new = set(), set()

    for og, oa, ng, na, line in mapping.pairs:
        ko, kn = (norm(og), norm(oa)), (norm(ng), norm(na))
        if ko not in full_old:
            raise DiffError(f"mapping line {line}: no agent {oa!r} in old group {og!r}")
        if kn not in full_new:
            raise DiffError(f"mapping line {line}: no agent {na!r} in new group {ng!r}")
        o, n = full_old[ko], full_new[kn]
        if not (keep(o) or keep(n)):
            continue
        changes.append(_pair(o, n, "mapped"))
        used_old.add(id(o))
        used_new.add(id(n))

    for key, o in old_idx.items():
        if id(o) in used_old:
            continue
        n = new_idx.get(key)
        if n is not None and id(n) not in used_new:
            changes.append(_pair(o, n, "exact"))
            used_old.add(id(o))
            used_new.add(id(n))

    lone_old = [(k, r) for k, r in old_idx.items() if id(r) not in used_old]
    lone_new = [(k, r) for k, r in new_idx.items() if id(r) not in used_new]

    def candidates(group_key: str, agent: str, pool) -> tuple[str, ...]:
        return tuple(r.agent for k, r in pool if k[0] == group_key and related(agent, r.agent))

    for k, o in lone_old:
        changes.append(Change(ChangeKind.REMOVED, o, None, candidates=candidates(k[0], o.agent, lone_new)))
    for k, n in lone_new:
        changes.append(Change(ChangeKind.ADDED, None, n, candidates=candidates(k[0], n.agent, lone_old)))

    order = {k: i for i, k in enumerate(ChangeKind)}
    changes.sort(key=lambda c: (norm(c.group), norm(c.agent), order[c.kind]))
    return Diff(old.version, new.version, changes, old.pin, new.pin)
