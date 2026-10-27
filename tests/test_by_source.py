"""Leave-one-source-out script on a hand-made export and re-interpreted table."""

import csv
import importlib.util
from pathlib import Path

from bpdrift.reinterpret import COLUMNS as RE_COLUMNS

SCRIPT = Path(__file__).resolve().parents[1] / "analysis" / "by_source.py"
spec = importlib.util.spec_from_file_location("by_source", SCRIPT)
by_source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(by_source)


def test_source_key():
    assert by_source.source_of("28776045;28512093") == by_source.source_of("28512093; 28776045")
    assert by_source.source_of("") == "none"


def test_own_and_leave_out(tmp_path):
    export = tmp_path / "bv.csv"
    with export.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Genome ID", "PubMed", "Testing Standard", "Laboratory Typing Platform"])
        w.writerow(["g1", "111", "CLSI", "Vitek"])         # line 2
        w.writerow(["g2", "111", "CLSI", "Vitek"])         # line 3
        w.writerow(["g3", "222;333", "EUCAST", ""])        # line 4
        w.writerow(["g4", "333;222", "EUCAST", ""])        # line 5
        w.writerow(["g5", "", "", ""])                     # line 6

    def rows(line, genome, ab, a, b, exclusion=""):
        out = []
        for v, c in (("8.1", a), ("16.1", b)):
            d = dict.fromkeys(RE_COLUMNS, "")
            d.update(source_line=line, genome_id=genome, antibiotic=ab, eucast_version=v,
                     category=c, exclusion=exclusion)
            out.append(d)
        return out

    re_rows = rows(2, "g1", "alpha", "I", "R") + rows(3, "g2", "alpha", "S", "S") \
        + rows(4, "g3", "alpha", "S", "S") + rows(5, "g4", "beta", "I", "R") \
        + rows(6, "g5", "alpha", "UNDECIDABLE", "R") + rows(6, "g5", "alpha", "S", "R", exclusion="urinary_only")
    re = tmp_path / "re.csv"
    with re.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=RE_COLUMNS)
        w.writeheader()
        w.writerows(re_rows)
    out = tmp_path / "out.csv"
    assert by_source.main([str(re), "--bvbrc", str(export), "-o", str(out)]) == 0
    res = {(r["source"], r["antibiotic"]): r for r in csv.DictReader(out.open(encoding="utf-8"))}

    a = res[("111", "ALL")]
    assert (a["genomes"], a["own_units"], a["own_flips"], a["leave_out_units"], a["leave_out_flips"]) == \
        ("2", "2", "1", "2", "1")
    b = res[("222;333", "ALL")]          # g3 and g4: one source despite the reversed PubMed order
    assert (b["genomes"], b["own_units"], b["own_flips"], b["testing_standard"]) == ("2", "2", "1", "EUCAST")
    assert res[("none", "ALL")]["own_units"] == "0"     # undecidable under v8.1
    assert ("111", "beta") not in res                   # no beta records in that source
