"""Run the parser on the user's own EUCAST tables, if any are present.

Place downloaded tables in data/breakpoints/ (git-ignored). The tests skip
when the folder is empty, so the public test suite needs no EUCAST data.
They check structure only and assert no breakpoint values.
"""

import re
from pathlib import Path

import pytest

from bpdrift.eucast import read_eucast_xlsx
from bpdrift.interpret import Category, interpret

DATA = Path(__file__).resolve().parents[1] / "data" / "breakpoints"
TABLES = sorted(DATA.glob("*.xlsx")) if DATA.is_dir() else []
KINDS = {"value", "none", "IE", "note", "NA", "IP", "empty"}

pytestmark = pytest.mark.skipif(not TABLES, reason="no EUCAST tables in data/breakpoints/")


def filename_version(path):
    m = re.search(r"v_?(\d+\.\d+)", path.name)
    return m.group(1) if m else None


@pytest.fixture(scope="module", params=TABLES, ids=lambda p: p.name)
def table(request):
    return request.param, read_eucast_xlsx(request.param)


def test_version_matches_filename(table):
    path, t = table
    expected = filename_version(path)
    if expected is None:
        pytest.skip(f"no version in file name {path.name}")
    assert t.version == expected


def test_has_many_rows_and_sheets(table):
    _, t = table
    assert len(t.rows) > 500
    assert len(t.sheets) > 15


def test_every_row_is_well_formed(table):
    _, t = table
    for r in t.rows:
        assert r.agent, r.cell
        assert r.s.kind in KINDS and r.r.kind in KINDS, r.cell
        # breakpoint() rejects S above R, so this checks every pair
        bp = r.breakpoint()
        if r.s.kind == "value" and r.r.kind == "value":
            assert bp.is_complete, r.cell


def test_a_good_share_have_numeric_breakpoints(table):
    _, t = table
    complete = [r for r in t.rows if r.breakpoint().is_complete]
    assert len(complete) > 0.2 * len(t.rows)


def test_every_complete_breakpoint_interprets(table):
    _, t = table
    for r in t.rows:
        bp = r.breakpoint()
        if bp.is_complete:
            low = interpret(f"<={bp.s_le:g}", bp).category
            high = interpret(f">{bp.r_gt:g}", bp).category
            assert low is Category.S and high is Category.R, r.cell
