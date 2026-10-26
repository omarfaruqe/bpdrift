"""Interpretation tests.

All breakpoint values here are synthetic test fixtures, not EUCAST data.
"""

import pytest

from bpdrift.interpret import Breakpoint, BreakpointError, Category, interpret

S, I, R, U = Category.S, Category.I, Category.R, Category.UNDECIDABLE

WITH_I = Breakpoint(s_le=0.5, r_gt=2)   # S <= 0.5, I 1-2, R > 2
NO_I = Breakpoint(s_le=1, r_gt=1)       # S <= 1, R > 1


@pytest.mark.parametrize(
    "mic, expected",
    [
        # exact values, including both boundaries
        ("0.25", S), ("0.5", S), ("1", I), ("2", I), ("4", R),
        # censored low
        ("<=0.25", S), ("<=0.5", S), ("<0.5", S), ("<=1", U), ("<1", U),
        # censored high
        (">2", R), (">=4", R), (">=2", U), (">1", U), (">0.5", U),
    ],
)
def test_with_i_zone(mic, expected):
    assert interpret(mic, WITH_I).category is expected


@pytest.mark.parametrize(
    "mic, expected",
    [
        ("1", S), ("2", R), ("<=1", S), (">1", R),
        ("<=2", U), (">=1", U), ("<2", U), (">0.5", U),
    ],
)
def test_without_i_zone(mic, expected):
    assert interpret(mic, NO_I).category is expected


def test_undecidable_reports_possible_categories():
    result = interpret("<=1", WITH_I)
    assert result.category is U
    assert result.possible == {S, I}
    assert "S/I" in result.reason

    wide = interpret(">0.25", WITH_I)
    assert wide.possible == {S, I, R}


def test_breakpoint_lowered_makes_old_value_undecidable():
    """The censoring case the study measures: a panel floor that was below
    the old S breakpoint is above the new one."""
    old = Breakpoint(s_le=2, r_gt=2)
    new = Breakpoint(s_le=0.5, r_gt=0.5)
    assert interpret("<=1", old).category is S
    assert interpret("<=1", new).category is U


def test_incomplete_breakpoint():
    result = interpret("1", Breakpoint(s_le=None, r_gt=None))
    assert result.category is Category.NO_BREAKPOINT
    assert interpret("1", Breakpoint(s_le=1, r_gt=None)).category is Category.NO_BREAKPOINT


def test_accepts_numbers():
    assert interpret(4, WITH_I).category is R
    assert interpret(0.5, WITH_I).category is S


@pytest.mark.parametrize("s, r", [(4, 2), (0, 1), (-1, 1), (1, float("inf"))])
def test_bad_breakpoints_rejected(s, r):
    with pytest.raises(BreakpointError):
        Breakpoint(s_le=s, r_gt=r)


def test_i_zone_flag():
    assert WITH_I.has_i_zone
    assert not NO_I.has_i_zone
