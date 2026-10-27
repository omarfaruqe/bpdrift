import csv
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "analysis" / "inputs_manifest.py"
spec = importlib.util.spec_from_file_location("inputs_manifest", SCRIPT)
inputs_manifest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inputs_manifest)


def test_manifest(tmp_path):
    from test_eucast import build

    data = tmp_path / "data"
    (data / "breakpoints").mkdir(parents=True)
    (data / "bvbrc").mkdir()
    (data / "ncbi").mkdir()
    build(data / "breakpoints" / "v_16.1_Breakpoint_Tables.xlsx")
    (data / "bvbrc" / "amr.csv").write_text("a,b\n1,2\n")
    (data / "ncbi" / "PDG000000012.1.metadata.tsv").write_text("x\n")
    (data / "ncbi" / "amrfinder_db_version.txt").write_text("2026-08-07.1\n")
    (data / "ncbi" / ".DS_Store").write_text("ignored")
    out = tmp_path / "s1.csv"
    assert inputs_manifest.main([str(data), "-o", str(out)]) == 0
    rows = {r["file"]: r for r in csv.DictReader(out.open(encoding="utf-8"))}
    assert set(rows) == {"v_16.1_Breakpoint_Tables.xlsx", "amr.csv", "PDG000000012.1.metadata.tsv",
                         "amrfinder_db_version.txt"}
    assert rows["v_16.1_Breakpoint_Tables.xlsx"]["version"] == "v16.1"
    assert rows["PDG000000012.1.metadata.tsv"]["version"] == "PDG000000012.1"
    assert rows["amrfinder_db_version.txt"]["version"] == "database 2026-08-07.1"
    assert len(rows["amr.csv"]["sha256"]) == 64
