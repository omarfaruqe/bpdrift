"""Robustness script on a hand-made re-interpreted table."""

import csv
import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("numpy")

from bpdrift.reinterpret import COLUMNS as RE_COLUMNS

SCRIPT = Path(__file__).resolve().parents[1] / "analysis" / "robustness.py"
spec = importlib.util.spec_from_file_location("robustness", SCRIPT)
rob = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rob)


def rows_for(line, genome, ab, mic, cats, exclusion="", sign_conflict="0"):
    out = []
    for v, c in zip(("8.1", "10.0", "16.1"), cats):
        d = dict.fromkeys(RE_COLUMNS, "")
        d.update(source_line=line, genome_id=genome, antibiotic=ab, mic=mic, eucast_version=v,
                 category=c, exclusion=exclusion, sign_conflict=sign_conflict, ratio_format="0")
        out.append(d)
    return out


@pytest.fixture
def result(tmp_path):
    rows = []
    rows += rows_for(1, "g1", "alpha", "16", ["I", "R", "R"])            # flips 8.1->10.0
    rows += rows_for(2, "g1", "alpha", "16", ["I", "R", "R"])            # same genome, same result
    rows += rows_for(3, "g2", "alpha", ">8", ["UNDECIDABLE", "R", "R"])  # censored, not decidable at 8.1
    rows += rows_for(4, "g3", "alpha", "2", ["S", "S", "S"])
    rows += rows_for(5, "g4", "alpha", "4", ["S", "S", "R"])
    rows += rows_for(6, "g4", "alpha", "4", ["S", "I", "R"])             # g4 disagrees at 10.0 -> dropped at genome level
    rows += rows_for(7, "g5", "alpha", "4", ["S", "R", "R"], exclusion="urinary_only")
    src = tmp_path / "re.csv"
    with src.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=RE_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    out, trans = tmp_path / "rob.csv", tmp_path / "trans.csv"
    assert rob.main([str(src), "-o", str(out), "--transitions", str(trans), "--replicates", "200"]) == 0
    res = {(r["scenario"], r["antibiotic"], r["version_a"], r["version_b"]): r
           for r in csv.DictReader(out.open(encoding="utf-8"))}
    tr = list(csv.DictReader(trans.open(encoding="utf-8")))
    return res, tr


def test_record_level(result):
    res, _ = result
    r = res[("record", "ALL", "8.1", "16.1")]
    # records 1, 2, 4, 5, 6 decidable under both; 1, 2, 5, 6 flip
    assert (r["units"], r["genomes"], r["flips"], r["r_line_flips"]) == ("5", "3", "4", "4")
    assert r["boot_lo"] != "" and float(r["boot_lo"]) <= 0.8 <= float(r["boot_hi"])


def test_genome_level_collapses_and_drops_conflicts(result):
    res, _ = result
    r = res[("genome", "ALL", "8.1", "16.1")]
    # g1 counted once; g4 dropped (records disagree at 10.0); g2 undecidable at 8.1
    assert (r["units"], r["flips"]) == ("2", "1")


def test_uncensored_and_fixed_cohort(result):
    res, _ = result
    assert res[("uncensored", "ALL", "8.1", "10.0")]["units"] == "5"       # record 3 (>8) out
    fixed = [r for k, r in res.items() if k[0] == "all_decidable" and k[1] == "ALL"]
    assert {r["units"] for r in fixed} == {"5"}                            # same cohort for every pair


def test_transitions_and_no_breakpoint_columns(result):
    res, tr = result
    rec = {(t["from"], t["to"]): int(t["units"]) for t in tr
           if t["scenario"] == "record" and t["antibiotic"] == "ALL" and t["version_a"] == "8.1"
           and t["version_b"] == "16.1"}
    assert rec == {("I", "R"): 2, ("S", "R"): 2, ("S", "S"): 1}
    assert not set(rob.COLUMNS) & {"eucast_S", "eucast_R", "old_S", "new_S"}
