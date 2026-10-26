"""Interpret an MIC against one breakpoint pair, EUCAST style.

EUCAST breakpoints are written as two numbers in mg/L:

    S <= s     susceptible
    R >  r     resistant
    between    I ("susceptible, increased exposure")

When s == r there is no I category.

A censored MIC is an interval of possible true values. The result is a
category only when the whole interval falls in that category; otherwise it
is UNDECIDABLE and the categories it could belong to are reported.

Semantics are continuous: "<1" means any value below 1, not "the next
dilution down". This is the conservative reading; it never turns a
censored value into a category it might not belong to.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

from bpdrift.mic import MIC, parse_mic


class Category(str, Enum):
    S = "S"
    I = "I"  # noqa: E741  (EUCAST's own letter)
    R = "R"
    UNDECIDABLE = "UNDECIDABLE"
    NO_BREAKPOINT = "NO_BREAKPOINT"


class BreakpointError(ValueError):
    """Raised for an inconsistent breakpoint pair."""


@dataclass(frozen=True)
class Breakpoint:
    """One EUCAST breakpoint pair. Either value may be None when the table
    gives no breakpoint (a dash, "IE", or a note)."""

    s_le: float | None
    r_gt: float | None

    def __post_init__(self) -> None:
        for name in ("s_le", "r_gt"):
            v = getattr(self, name)
            if v is not None and (not math.isfinite(v) or v <= 0):
                raise BreakpointError(f"{name} must be positive and finite, got {v!r}")
        if self.s_le is not None and self.r_gt is not None and self.s_le > self.r_gt:
            raise BreakpointError(
                f"S breakpoint {self.s_le} is above R breakpoint {self.r_gt}"
            )

    @property
    def is_complete(self) -> bool:
        return self.s_le is not None and self.r_gt is not None

    @property
    def has_i_zone(self) -> bool:
        return self.is_complete and self.s_le < self.r_gt


@dataclass(frozen=True)
class Interpretation:
    category: Category
    possible: frozenset[Category] = field(default_factory=frozenset)
    reason: str = ""


def _interval(mic: MIC) -> tuple[float, bool, float, bool]:
    """Possible true values as (lo, lo_inclusive, hi, hi_inclusive)."""
    v = mic.value
    if mic.sign == "=":
        return v, True, v, True
    if mic.sign == "<=":
        return 0.0, False, v, True
    if mic.sign == "<":
        return 0.0, False, v, False
    if mic.sign == ">=":
        return v, True, math.inf, False
    return v, False, math.inf, False  # ">"


def interpret(mic: MIC | str | float, bp: Breakpoint) -> Interpretation:
    """Interpret one MIC under one breakpoint pair."""
    if not isinstance(mic, MIC):
        mic = parse_mic(mic)

    if not bp.is_complete:
        return Interpretation(Category.NO_BREAKPOINT, frozenset(), "breakpoint pair incomplete")

    s, r = bp.s_le, bp.r_gt
    lo, lo_incl, hi, _hi_incl = _interval(mic)

    possible = set()
    # S region: (0, s]
    if lo < s or (lo == s and lo_incl):
        possible.add(Category.S)
    # I region: (s, r], only when s < r
    if s < r and hi > s and (lo < r or (lo == r and lo_incl)):
        possible.add(Category.I)
    # R region: (r, inf)
    if hi > r:
        possible.add(Category.R)

    frozen = frozenset(possible)
    if len(frozen) == 1:
        (only,) = frozen
        return Interpretation(only, frozen, f"{mic} lies entirely in {only.value}")

    names = "/".join(c.value for c in (Category.S, Category.I, Category.R) if c in frozen)
    return Interpretation(
        Category.UNDECIDABLE,
        frozen,
        f"{mic} could be {names} under S<={s:g}, R>{r:g}",
    )
