"""Helpers of the review-checks script."""

import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "analysis" / "review_checks.py"
spec = importlib.util.spec_from_file_location("review_checks", SCRIPT)
rc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rc)


def test_bracketed():
    assert rc.bracketed("(8)") and rc.bracketed(" (0.5) ")
    assert not rc.bracketed("8") and not rc.bracketed("") and not rc.bracketed("(8")


def test_atu_range():
    assert rc.atu_range("16") == (16.0, 16.0)
    assert rc.atu_range("0.5-1") == (0.5, 1.0)
    assert rc.atu_range("") is None and rc.atu_range("-") is None and rc.atu_range("Note") is None
    assert rc.atu_range("1-2-4") is None


def test_in_atu_only_for_exact_mics():
    assert rc.in_atu("16", (16.0, 16.0))
    assert rc.in_atu("1", (0.5, 1.0)) and not rc.in_atu("2", (0.5, 1.0))
    assert not rc.in_atu(">16", (16.0, 16.0))       # censored: not known to lie in the ATU
    assert not rc.in_atu("16", None)


def test_rate_row():
    row = rc.rate_row("x", "c", "ALL", "m", 1, 4)
    assert row[:6] == ["x", "c", "ALL", "m", 1, 4] and row[6] == "0.2500"
