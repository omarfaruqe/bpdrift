"""Predictor-bound script on a hand-made pairs file."""

import csv
import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "analysis" / "predictor_bounds.py"
spec = importlib.util.spec_from_file_location("predictor_bounds", SCRIPT)
pb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pb)

ROWS = [
    ("g1", "alpha", "8.1", 1, "I"), ("g1", "alpha", "16.1", 1, "R"),   # I -> R, predicted R: error lost
    ("g2", "alpha", "8.1", 0, "I"), ("g2", "alpha", "16.1", 0, "R"),   # I -> R, predicted not R: error gained
    ("g3", "alpha", "8.1", 0, "S"), ("g3", "alpha", "16.1", 0, "R"),   # S -> R, predicted not R: error gained
    ("g4", "alpha", "8.1", 1, "R"), ("g4", "alpha", "16.1", 1, "R"),   # unchanged
]


@pytest.fixture
def res(tmp_path):
    src = tmp_path / "pairs.csv"
    with src.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["genome_id", "antibiotic", "eucast_version", "predicted_resistant", "phenotype"])
        w.writerows(ROWS)
    out = tmp_path / "b.csv"
    assert pb.main([str(src), "-o", str(out)]) == 0
    return {(r["i_as"], r["antibiotic"]): r for r in csv.DictReader(out.open(encoding="utf-8"))}


def test_i_not_resistant(res):
    r = res[("not_r", "alpha")]
    assert (r["pairs"], r["label_changes"], r["errors_lost"], r["errors_gained"]) == ("4", "3", "1", "2")
    assert float(r["bound"]) == pytest.approx(0.75)
    assert float(r["observed_err_change"]) == pytest.approx(0.25)          # (2 - 1) / 4
    assert float(r["observed_position"]) == pytest.approx(1 / 3, abs=1e-4)  # (2 - 1) / 3
    assert float(r["oracle_vme_rate"]) == pytest.approx(3 / 4)             # 3 changed pairs now R, of 4 R


def test_i_as_r_and_drop(res):
    r = res[("r", "alpha")]                    # I counts as R: only g3 (S -> R) changes
    assert (r["label_changes"], r["errors_gained"]) == ("1", "1")
    d = res[("drop", "alpha")]                 # pairs with I left out: g3 and g4
    assert (d["pairs"], d["label_changes"]) == ("2", "1")
