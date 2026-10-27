"""The analysis script on a hand-made re-interpreted table."""

import csv
import importlib.util
from pathlib import Path

from bpdrift.reinterpret import COLUMNS as RE_COLUMNS

SCRIPT = Path(__file__).resolve().parents[1] / "analysis" / "flip_rates.py"
spec = importlib.util.spec_from_file_location("flip_rates", SCRIPT)
flip_rates = importlib.util.module_from_spec(spec)
spec.loader.exec_module(flip_rates)


def row(line, ab, mic, version, s, r, category, exclusion="", sign_conflict=0, std="CLSI"):
    d = dict.fromkeys(RE_COLUMNS, "")
    d.update(source_line=line, antibiotic=ab, mic=mic, eucast_version=version, eucast_S=s, eucast_R=r,
             category=category, exclusion=exclusion, sign_conflict=sign_conflict, ratio_format=0,
             reported_standard=std)
    return d


def test_scenarios(tmp_path):
    rows = [
        # record 1: 16 mg/L, I under 8/16, R under 8/8 -> one flip in every scenario
        row(1, "alpha", "16", "10.0", "8", "16", "I", std="EUCAST"),
        row(1, "alpha", "16", "16.1", "8", "8", "R", std="EUCAST"),
        # record 2: sign conflict ">8": UNDECIDABLE under 8/16, R under 8/8; exact 8 -> S, S
        row(2, "alpha", ">8", "10.0", "8", "16", "UNDECIDABLE", sign_conflict=1),
        row(2, "alpha", ">8", "16.1", "8", "8", "R", sign_conflict=1),
        # record 3: excluded from the main analysis, flips S -> R
        row(3, "alpha", "4", "10.0", "(4)", "(4)", "S", exclusion="urinary_only"),
        row(3, "alpha", "4", "16.1", "2", "2", "R", exclusion="urinary_only"),
    ]
    src = tmp_path / "re.csv"
    with src.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=RE_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    out = tmp_path / "flips.csv"
    assert flip_rates.main([str(src), "-o", str(out)]) == 0
    res = {(r["scenario"], r["antibiotic"]): r for r in csv.DictReader(out.open(encoding="utf-8"))}

    main = res[("main", "alpha")]
    assert (main["records"], main["both_decidable"], main["flips"], main["undecidable_a"]) == ("2", "1", "1", "1")
    assert main["transitions"] == "I>R:1"

    exact = res[("sign_exact", "alpha")]
    assert (exact["both_decidable"], exact["flips"], exact["undecidable_a"]) == ("2", "1", "0")

    everything = res[("all_records", "ALL")]
    assert (everything["records"], everything["flips"], everything["r_line_flips"]) == ("3", "2", "2")

    eu = res[("eucast_stated", "alpha")]
    assert (eu["records"], eu["flips"]) == ("1", "1")
    cl = res[("clsi_stated", "alpha")]
    assert (cl["records"], cl["both_decidable"], cl["flips"]) == ("1", "0", "0")
