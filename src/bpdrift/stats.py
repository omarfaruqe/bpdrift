"""Counting label changes between breakpoint versions.

Definitions used throughout (they are the paper's definitions):

* decidable: S, I or R. UNDECIDABLE and NO_BREAKPOINT are not.
* flip: a record decidable under both versions whose category differs.
* R-line flip: a flip that crosses the resistant line, i.e. S or I under
  one version and R under the other. This is the change that alters a
  binary resistant / not-resistant label.

Rates come with Wilson score 95% intervals.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Mapping

from bpdrift.interpret import Category

DECIDABLE = frozenset({Category.S, Category.I, Category.R})


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float, float]:
    """Proportion k/n with its Wilson score interval; (nan, nan, nan) when n is 0."""
    if n == 0:
        return math.nan, math.nan, math.nan
    if not 0 <= k <= n:
        raise ValueError(f"need 0 <= k <= n, got k={k}, n={n}")
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    lo = 0.0 if k == 0 else max(0.0, centre - half)
    hi = 1.0 if k == n else min(1.0, centre + half)
    return p, lo, hi


@dataclass(frozen=True)
class Comparison:
    version_a: str
    version_b: str
    records: int
    both_decidable: int
    flips: int
    r_line_flips: int
    undecidable_a: int
    undecidable_b: int
    no_breakpoint_a: int
    no_breakpoint_b: int
    transitions: Counter  # (category_a, category_b) -> count, decidable pairs only

    @property
    def flip_rate(self) -> tuple[float, float, float]:
        return wilson(self.flips, self.both_decidable)

    @property
    def r_line_rate(self) -> tuple[float, float, float]:
        return wilson(self.r_line_flips, self.both_decidable)

    @property
    def undecidable_rate_a(self) -> tuple[float, float, float]:
        return wilson(self.undecidable_a, self.records - self.no_breakpoint_a)

    @property
    def undecidable_rate_b(self) -> tuple[float, float, float]:
        return wilson(self.undecidable_b, self.records - self.no_breakpoint_b)


def compare(
    labels: Iterable[Mapping[str, Category]], version_a: str, version_b: str
) -> Comparison:
    """Compare per-record labels between two versions.

    labels: one mapping per record, version -> category. A record missing
    either version raises KeyError rather than being skipped silently.
    """
    n = both = flips = r_line = ua = ub = na = nb = 0
    transitions: Counter = Counter()
    for rec in labels:
        a, b = rec[version_a], rec[version_b]
        n += 1
        ua += a is Category.UNDECIDABLE
        ub += b is Category.UNDECIDABLE
        na += a is Category.NO_BREAKPOINT
        nb += b is Category.NO_BREAKPOINT
        if a in DECIDABLE and b in DECIDABLE:
            both += 1
            transitions[(a.value, b.value)] += 1
            if a is not b:
                flips += 1
                if (a is Category.R) != (b is Category.R):
                    r_line += 1
    return Comparison(version_a, version_b, n, both, flips, r_line, ua, ub, na, nb, transitions)
