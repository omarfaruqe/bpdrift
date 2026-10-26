"""CLI tests on synthetic workbooks (see test_eucast.build)."""

import csv
import json

import pytest

pytest.importorskip("openpyxl")
from test_eucast import build  # noqa: E402

from bpdrift.cli import main  # noqa: E402


@pytest.fixture
def tables(tmp_path):
    old = build(tmp_path / "old.xlsx", version="13.1", extra_value="1")
    new = build(tmp_path / "new.xlsx", version="16.1", extra_value="0.5")
    return old, new


def test_pin_prints_and_writes(tmp_path, tables, capsys):
    old, _ = tables
    assert main(["pin", str(old), "--version", "13.1", "--retrieved", "2026-10-26"]) == 0
    assert capsys.readouterr().out.startswith("EUCAST v13.1 (2026-10-26, sha256:")
    out = tmp_path / "old.pin.json"
    assert main(["pin", str(old), "--version", "13.1", "--retrieved", "2026-10-26", "-o", str(out)]) == 0
    assert json.loads(out.read_text())["version"] == "13.1"


def test_interpret(tables, capsys):
    _, new = tables
    assert main(["interpret", "<=0.25", "--table", str(new), "--group", "Testobacterales",
                 "--agent", "betacillin"]) == 0
    category, reason, where = capsys.readouterr().out.strip().split("\t")
    assert category == "S"
    assert where.startswith("EUCAST v16.1 Testobacterales!A")


def test_interpret_by_subgroup_and_ambiguity(tables, capsys):
    _, new = tables
    assert main(["interpret", "1", "--table", str(new), "--group", "Betoides spp.",
                 "--agent", "Alphacillin"]) == 0
    assert capsys.readouterr().out.startswith("R")
    # Alphacillin is in two subgroups of the Anaerobes sheet
    assert main(["interpret", "1", "--table", str(new), "--group", "Anaerobes",
                 "--agent", "Alphacillin"]) == 1
    assert "ambiguous" in capsys.readouterr().err


def test_interpret_unknown_agent(tables, capsys):
    _, new = tables
    assert main(["interpret", "1", "--table", str(new), "--group", "Testobacterales",
                 "--agent", "Nope"]) == 1
    assert "no agent 'Nope'" in capsys.readouterr().err


def test_diff_writes_csv_with_flips(tmp_path, tables, capsys):
    old, new = tables
    out = tmp_path / "diff.csv"
    assert main(["diff", str(old), str(new), "--changed-only", "-o", str(out)]) == 0
    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    (row,) = rows
    assert row["kind"] == "changed" and row["new_agent"] == "Omegamycin"
    assert (row["old_S"], row["new_S"]) == ("1", "0.5")
    assert row["flips"] == "1:S>I"
    assert (row["old_version"], row["new_version"]) == ("13.1", "16.1")
    assert "v13.1 -> v16.1" in capsys.readouterr().err


def test_diff_group_filter_to_stdout(tables, capsys):
    old, new = tables
    assert main(["diff", str(old), str(new), "--group", "Alphoides spp."]) == 0
    rows = list(csv.DictReader(capsys.readouterr().out.splitlines()))
    assert [(r["kind"], r["new_agent"]) for r in rows] == [("unchanged", "Alphacillin")]


def test_diff_checks_pins(tmp_path, tables, capsys):
    old, new = tables
    pin = tmp_path / "wrong.json"
    main(["pin", str(new), "--version", "16.1", "--retrieved", "2026-10-26", "-o", str(pin)])
    capsys.readouterr()
    assert main(["diff", str(old), str(new), "--old-pin", str(pin)]) == 1
    assert "does not match pin" in capsys.readouterr().err


def test_bad_mic_is_an_error_not_a_crash(tables, capsys):
    _, new = tables
    assert main(["interpret", "eight", "--table", str(new), "--group", "Testobacterales",
                 "--agent", "Betacillin"]) == 1
    assert "cannot parse MIC" in capsys.readouterr().err
