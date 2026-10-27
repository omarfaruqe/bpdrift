"""Paired genotype error-rate differences on a hand-made pairs file."""

import csv
import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("numpy")

SCRIPT = Path(__file__).resolve().parents[1] / "analysis" / "genotype_diff.py"
spec = importlib.util.spec_from_file_location("genotype_diff", SCRIPT)
gd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gd)

HEADER = ["genome_id", "antibiotic", "eucast_version", "predicted_resistant", "phenotype"]
ROWS = [
    # g1: predicted R; I under 8.1 (ME) -> R under 16.1 (correct): error disappears
    ("g1", "alpha", "8.1", 1, "I"), ("g1", "alpha", "16.1", 1, "R"),
    # g2: predicted not R; I -> R: an error appears (VME)
    ("g2", "alpha", "8.1", 0, "I"), ("g2", "alpha", "16.1", 0, "R"),
    # g3: predicted R; R in both: correct both times
    ("g3", "alpha", "8.1", 1, "R"), ("g3", "alpha", "16.1", 1, "R"),
    # g4: predicted not R; S -> R: error appears, survives --i-as drop
    ("g4", "alpha", "8.1", 0, "S"), ("g4", "alpha", "16.1", 0, "R"),
    # g5: only one version -> unpaired, left out
    ("g5", "alpha", "8.1", 0, "S"),
    # g6, g7: a second antibiotic, for the multiple-testing columns
    ("g6", "beta", "8.1", 1, "S"), ("g6", "beta", "16.1", 1, "S"),
    ("g7", "beta", "8.1", 0, "R"), ("g7", "beta", "16.1", 0, "R"),
]


def run(tmp_path, *extra):
    src = tmp_path / "pairs.csv"
    with src.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(HEADER)
        w.writerows(ROWS)
    out = tmp_path / "diff.csv"
    assert gd.main([str(src), "-o", str(out), "--replicates", "200", *extra]) == 0
    return {x["antibiotic"]: x for x in csv.DictReader(out.open(encoding="utf-8"))}


def test_mcnemar_and_holm():
    assert gd.mcnemar_exact(0, 0) == 1.0
    assert gd.mcnemar_exact(0, 5) == pytest.approx(2 / 32)
    assert gd.mcnemar_exact(3, 3) == 1.0
    assert gd.mcnemar_exact(1, 9) == pytest.approx(22 / 1024)
    assert gd.holm([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.06, 0.06])


def test_error_kinds_count_i_as_not_resistant():
    assert gd.error(0, "R") == "vme"
    assert gd.error(1, "I") == "me"
    assert gd.error(1, "S") == "me"
    assert gd.error(0, "I") == "" and gd.error(1, "R") == ""


def test_three_estimands(tmp_path):
    r = run(tmp_path)["alpha"]
    assert r["pairs"] == "4"
    assert (r["vme_num_a"], r["vme_den_a"], r["vme_num_b"], r["vme_den_b"]) == ("0", "1", "2", "4")
    assert (r["me_num_a"], r["me_den_a"], r["me_num_b"], r["me_den_b"]) == ("1", "3", "0", "0")
    assert (r["err_num_a"], r["err_den_a"], r["err_num_b"], r["err_den_b"]) == ("1", "4", "2", "4")
    assert float(r["err_diff"]) == pytest.approx(0.25, abs=1e-4)
    assert (r["error_to_correct"], r["correct_to_error"]) == ("1", "2")
    assert r["mcnemar_p_holm"] != "" and r["err_lo_adj"] != ""


def test_pooled_row_is_not_adjusted(tmp_path):
    res = run(tmp_path)
    assert res["ALL"]["mcnemar_p_holm"] == "" and res["ALL"]["vme_lo_adj"] == ""
    assert res["ALL"]["pairs"] == "6"
    # beta's labels never change: reported, but outside the testing family
    assert res["beta"]["mcnemar_p_holm"] == "" and res["beta"]["err_lo_adj"] == ""
    assert res["beta"]["err_diff"] == "0.0000"


def test_i_dropped(tmp_path):
    r = run(tmp_path, "--i-as", "drop")["alpha"]
    assert r["pairs"] == "2"                                   # g1, g2 dropped (I under 8.1)
    assert (r["vme_num_a"], r["vme_den_a"], r["vme_num_b"], r["vme_den_b"]) == ("0", "1", "1", "2")


def test_rejects_changing_predictions(tmp_path):
    src = tmp_path / "pairs.csv"
    src.write_text(",".join(HEADER) + "\ng1,alpha,8.1,1,R\ng1,alpha,16.1,0,R\n", encoding="utf-8")
    with pytest.raises(SystemExit, match="predictions differ"):
        gd.main([str(src), "-o", str(tmp_path / "o.csv"), "--replicates", "10"])


def test_i_as_r(tmp_path):
    r = run(tmp_path, "--i-as", "r")["alpha"]
    # I becomes R: under 8.1, g1 (predicted R) and g2 (predicted not R) count as R, so g2 is a VME at both versions
    assert (r["vme_num_a"], r["vme_den_a"], r["vme_num_b"], r["vme_den_b"]) == ("1", "3", "2", "4")
    assert (r["me_num_a"], r["me_den_a"]) == ("0", "1")
    assert (r["error_to_correct"], r["correct_to_error"]) == ("0", "1")   # only g4 (S -> R) changes
