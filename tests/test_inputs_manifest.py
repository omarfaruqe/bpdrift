import csv
import importlib.util
import os
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "analysis" / "inputs_manifest.py"
spec = importlib.util.spec_from_file_location("inputs_manifest", SCRIPT)
inputs_manifest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inputs_manifest)


def make_data(tmp_path):
    from test_eucast import build

    data = tmp_path / "data"
    (data / "breakpoints").mkdir(parents=True)
    (data / "bvbrc").mkdir()
    (data / "ncbi" / "ast").mkdir(parents=True)
    build(data / "breakpoints" / "v_16.1_Breakpoint_Tables.xlsx")
    (data / "bvbrc" / "amr.csv").write_text("a,b\n1,2\n")
    (data / "ncbi" / "PDG000000012.1.metadata.tsv").write_text("x\n")
    (data / "ncbi" / "amrfinder_db_version.txt").write_text("2026-08-07.1\n")
    (data / "ncbi" / "ast" / "biosample_000000.xml").write_text("<BioSampleSet/>\n")
    (data / "ncbi" / ".DS_Store").write_text("ignored")
    return data


def test_manifest(tmp_path):
    data = make_data(tmp_path)
    out = tmp_path / "s1.csv"
    assert inputs_manifest.main([str(data), "-o", str(out)]) == 0
    rows = {r["path"]: r for r in csv.DictReader(out.open(encoding="utf-8"))}
    assert set(rows) == {"data/breakpoints/v_16.1_Breakpoint_Tables.xlsx", "data/bvbrc/amr.csv",
                         "data/ncbi/PDG000000012.1.metadata.tsv", "data/ncbi/amrfinder_db_version.txt",
                         "data/ncbi/ast/biosample_000000.xml"}
    assert rows["data/breakpoints/v_16.1_Breakpoint_Tables.xlsx"]["version"] == "v16.1"
    assert rows["data/ncbi/PDG000000012.1.metadata.tsv"]["version"] == "PDG000000012.1"
    assert rows["data/ncbi/amrfinder_db_version.txt"]["version"] == "database 2026-08-07.1"
    assert rows["data/ncbi/ast/biosample_000000.xml"]["source"].startswith("NCBI BioSample")
    assert len(rows["data/bvbrc/amr.csv"]["sha256"]) == 64
    assert "file_date" not in next(iter(rows.values()))


def test_manifest_ignores_modification_times(tmp_path):
    data = make_data(tmp_path)
    first, second = tmp_path / "a.csv", tmp_path / "b.csv"
    inputs_manifest.main([str(data), "-o", str(first)])
    f = data / "bvbrc" / "amr.csv"
    os.utime(f, (1_000_000_000, 1_000_000_000))              # same bytes, "downloaded" in 2001
    inputs_manifest.main([str(data), "-o", str(second)])
    assert first.read_bytes() == second.read_bytes()


def test_check(tmp_path, capsys):
    data = make_data(tmp_path)
    manifest = tmp_path / "s1.csv"
    inputs_manifest.main([str(data), "-o", str(manifest)])
    assert inputs_manifest.main([str(data), "--check", str(manifest)]) == 0

    (data / "bvbrc" / "amr.csv").write_text("a,b\n1,3\n")                     # same size, other bytes
    (data / "ncbi" / "PDG000000012.1.metadata.tsv").unlink()                   # missing
    (data / "ncbi" / "PDG000000012.2.metadata.tsv").write_text("y\n")          # unlisted extra
    assert inputs_manifest.main([str(data), "--check", str(manifest)]) == 1
    out = capsys.readouterr().out
    assert "checksum differs: data/bvbrc/amr.csv" in out
    assert "missing: data/ncbi/PDG000000012.1.metadata.tsv" in out
    assert "not in the manifest: data/ncbi/PDG000000012.2.metadata.tsv" in out
