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


def test_mcnemar_exact_known_values():
    assert gd.mcnemar_exact(0, 0) == 1.0
    assert gd.mcnemar_exact(0, 5) == pytest.approx(2 / 32)          # 2 * 0.5**5
    assert gd.mcnemar_exact(3, 3) == 1.0
    # binomial(10, 0.5): P(X <= 1) = 11/1024, two-sided 22/1024
    assert gd.mcnemar_exact(1, 9) == pytest.approx(22 / 1024)


def test_error_kinds_count_i_as_not_resistant():
    assert gd.error(0, "R") == "vme"
    assert gd.error(1, "I") == "me"       # predicted resistant, phenotype I -> major error
    assert gd.error(1, "S") == "me"
    assert gd.error(0, "I") == "" and gd.error(1, "R") == ""


def test_paired_differences(tmp_path):
    rows = [
        # g1: predicted R; I under 8.1 (ME) -> R under 16.1 (correct): error disappears
        ("g1", "alpha", "8.1", 1, "I"), ("g1", "alpha", "16.1", 1, "R"),
        # g2: predicted not R; I -> R: an error appears (VME)
        ("g2", "alpha", "8.1", 0, "I"), ("g2", "alpha", "16.1", 0, "R"),
        # g3: predicted R; R in both: correct both times
        ("g3", "alpha", "8.1", 1, "R"), ("g3", "alpha", "16.1", 1, "R"),
        # g4: only one version -> unpaired, left out
        ("g4", "alpha", "8.1", 0, "S"),
    ]
    src = tmp_path / "pairs.csv"
    with src.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["genome_id", "antibiotic", "eucast_version", "predicted_resistant", "phenotype"])
        w.writerows(rows)
    out = tmp_path / "diff.csv"
    assert gd.main([str(src), "-o", str(out), "--replicates", "200"]) == 0
    r = {x["antibiotic"]: x for x in csv.DictReader(out.open(encoding="utf-8"))}["alpha"]
    assert (r["pairs"], r["vme_a"], r["r_a"], r["vme_b"], r["r_b"]) == ("3", "0", "1", "1", "3")
    assert (r["me_a"], r["notr_a"], r["me_b"], r["notr_b"]) == ("1", "2", "0", "0")
    assert (r["error_to_correct"], r["correct_to_error"]) == ("1", "1")
    assert float(r["vme_diff"]) == pytest.approx(1 / 3, abs=1e-4)


def test_rejects_changing_predictions(tmp_path):
    src = tmp_path / "pairs.csv"
    src.write_text("genome_id,antibiotic,eucast_version,predicted_resistant,phenotype\n"
                   "g1,alpha,8.1,1,R\ng1,alpha,16.1,0,R\n", encoding="utf-8")
    with pytest.raises(SystemExit, match="predictions differ"):
        gd.main([str(src), "-o", str(tmp_path / "o.csv"), "--replicates", "10"])
