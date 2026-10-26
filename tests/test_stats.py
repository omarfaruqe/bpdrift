import math

import pytest

from bpdrift.interpret import Category
from bpdrift.stats import compare, wilson

S, I, R, U, NB = Category.S, Category.I, Category.R, Category.UNDECIDABLE, Category.NO_BREAKPOINT


def test_wilson_known_values():
    # reference values: 10/100 -> 0.0552 to 0.1744 (Wilson, z = 1.96)
    p, lo, hi = wilson(10, 100)
    assert p == 0.1
    assert lo == pytest.approx(0.0552, abs=5e-4)
    assert hi == pytest.approx(0.1744, abs=5e-4)


def test_wilson_edges():
    p, lo, hi = wilson(0, 50)
    assert (p, lo) == (0.0, 0.0) and 0 < hi < 0.08
    p, lo, hi = wilson(50, 50)
    assert (p, hi) == (1.0, 1.0) and 0.92 < lo < 1
    assert all(math.isnan(x) for x in wilson(0, 0))
    with pytest.raises(ValueError):
        wilson(5, 4)


def test_compare_counts_every_case():
    labels = [
        {"a": S, "b": S},     # same
        {"a": I, "b": R},     # flip across the R line
        {"a": S, "b": I},     # flip, not across the R line
        {"a": R, "b": S},     # flip across the R line
        {"a": U, "b": R},     # undecidable under a: not compared
        {"a": S, "b": NB},    # no breakpoint under b: not compared
    ]
    c = compare(labels, "a", "b")
    assert (c.records, c.both_decidable, c.flips, c.r_line_flips) == (6, 4, 3, 2)
    assert (c.undecidable_a, c.undecidable_b, c.no_breakpoint_a, c.no_breakpoint_b) == (1, 0, 0, 1)
    assert c.flip_rate[0] == 0.75
    assert c.transitions[("I", "R")] == 1 and sum(c.transitions.values()) == 4
    # undecidable share excludes records with no breakpoint under that version
    assert c.undecidable_rate_a[0] == pytest.approx(1 / 6)
    assert c.undecidable_rate_b[0] == 0.0


def test_compare_missing_version_fails_loudly():
    with pytest.raises(KeyError):
        compare([{"a": S}], "a", "b")
