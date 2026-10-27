"""Exhaustive check of interpret() against an independent brute-force oracle.

For every breakpoint pair (s <= r) on a doubling-dilution ladder and every MIC
on that ladder with every censoring sign, the set of possible categories is
computed a second, naive way: enumerate candidate true values inside the MIC's
interval and categorise each one directly. Category boundaries lie on the
ladder, so candidates at every ladder value, every geometric midpoint, and just
inside each open end hit every category region the interval touches.

Breakpoints here are synthetic, not EUCAST data.
"""

import itertools
import math

import pytest

from bpdrift.interpret import Breakpoint, Category, interpret
from bpdrift.mic import MIC

LADDER = [2.0 ** k for k in range(-6, 8)]          # 0.015625 .. 128 mg/L
SIGNS = ["=", "<=", "<", ">=", ">"]
EPS = 1e-9


def category_of(x: float, s: float, r: float) -> Category:
    if x <= s:
        return Category.S
    if x > r:
        return Category.R
    return Category.I


def inside(x: float, v: float, sign: str) -> bool:
    return {"=": x == v, "<=": x <= v, "<": x < v, ">=": x >= v, ">": x > v}[sign]


def candidates(v: float) -> list[float]:
    pts = set(LADDER)
    pts |= {math.sqrt(a * b) for a, b in zip(LADDER, LADDER[1:])}
    pts |= {LADDER[0] / 2, LADDER[-1] * 2, LADDER[-1] * 1e6}
    pts |= {v * (1 - EPS), v * (1 + EPS), v}
    return sorted(pts)


def oracle(v: float, sign: str, s: float, r: float) -> frozenset:
    return frozenset(category_of(x, s, r) for x in candidates(v) if inside(x, v, sign))


PAIRS = [(s, r) for s, r in itertools.product(LADDER, LADDER) if s <= r]


@pytest.mark.parametrize("s, r", PAIRS)
def test_interpret_matches_oracle(s, r):
    bp = Breakpoint(s, r)
    for v, sign in itertools.product(LADDER, SIGNS):
        got = interpret(MIC(v, sign), bp)
        want = oracle(v, sign, s, r)
        assert got.possible == want, (v, sign, s, r)
        if len(want) == 1:
            assert got.category is next(iter(want)), (v, sign, s, r)
        else:
            assert got.category is Category.UNDECIDABLE, (v, sign, s, r)


def test_exact_mic_is_always_decidable():
    for (s, r), v in itertools.product(PAIRS, LADDER):
        assert interpret(MIC(v, "="), Breakpoint(s, r)).category in (Category.S, Category.I, Category.R)


def test_case_count():
    # 105 breakpoint pairs x 14 values x 5 signs
    assert len(PAIRS) * len(LADDER) * len(SIGNS) == 7350
