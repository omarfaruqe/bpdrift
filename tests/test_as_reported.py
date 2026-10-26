import csv
import importlib.util
from pathlib import Path

from bpdrift.reinterpret import COLUMNS as RE_COLUMNS

SCRIPT = Path(__file__).resolve().parents[1] / "analysis" / "as_reported.py"
spec = importlib.util.spec_from_file_location("as_reported", SCRIPT)
as_reported = importlib.util.module_from_spec(spec)
spec.loader.exec_module(as_reported)


def row(line, reported, version, category, std="EUCAST", year="2020", exclusion=""):
    d = dict.fromkeys(RE_COLUMNS, "")
    d.update(source_line=line, reported_phenotype=reported, eucast_version=version, category=category,
             reported_standard=std, reported_standard_year=year, exclusion=exclusion)
    return d


def test_agreement_strata(tmp_path):
    rows = [
        row(1, "Resistant", "10.0", "R"),
        row(2, "Intermediate", "10.0", "R"),          # disagrees, also across the R line
        row(3, "Susceptible", "10.0", "I"),           # disagrees, same side of the R line
        row(4, "Susceptible", "10.0", "UNDECIDABLE"), # reported but not compared
        row(5, "", "10.0", "S"),                      # no reported call: ignored
        row(6, "Resistant", "10.0", "R", exclusion="urinary_only"),  # not main: ignored
        row(7, "Resistant", "16.1", "R", std="CLSI", year="2017"),
    ]
    src = tmp_path / "re.csv"
    with src.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=RE_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    out = tmp_path / "ar.csv"
    assert as_reported.main([str(src), "-o", str(out)]) == 0
    res = {(r["standard"], r["year"], r["eucast_version"]): r for r in csv.DictReader(out.open(encoding="utf-8"))}

    e = res[("EUCAST", "2020", "10.0")]
    assert (e["reported"], e["compared"], e["agreement"], e["r_line_agreement"]) == ("4", "3", "1", "2")
    assert e["in_force_that_year"] == "yes"
    assert e["disagreements"] == "I>R:1;S>I:1"
    assert res[("EUCAST", "all", "10.0")]["compared"] == "3"
    assert res[("CLSI", "2017", "16.1")]["in_force_that_year"] == ""


def test_in_force():
    assert as_reported.in_force("EUCAST", "2020", "10.0") == "yes"
    assert as_reported.in_force("EUCAST", "2021", "10.0") == "no"
    assert as_reported.in_force("CLSI", "2020", "10.0") == ""
    assert as_reported.in_force("EUCAST", "none", "10.0") == ""
