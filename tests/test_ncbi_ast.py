"""NCBI antibiogram converter on a synthetic BioSample batch."""

import csv
import importlib.util
from pathlib import Path

from bpdrift.bvbrc import read_bvbrc_csv

SCRIPT = Path(__file__).resolve().parents[1] / "analysis" / "ncbi_ast.py"
spec = importlib.util.spec_from_file_location("ncbi_ast", SCRIPT)
ncbi_ast = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ncbi_ast)

HEADER = ("Antibiotic", "Resistance phenotype", "Measurement sign", "Measurement", "Measurement units",
          "Laboratory typing method", "Laboratory typing platform", "Vendor",
          "Laboratory typing method version or reagent", "Testing standard")


def sample(acc, rows, host="Homo sapiens"):
    head = "".join(f"<Cell>{h}</Cell>" for h in HEADER)
    body = "".join("<Row>" + "".join(f"<Cell>{v}</Cell>" if v else "<Cell/>" for v in r) + "</Row>" for r in rows)
    return (f'<BioSample accession="{acc}"><Description><Organism taxonomy_id="573" '
            f'taxonomy_name="Klebsiella pneumoniae"/><Comment><Table class="Antibiogram.1.0">'
            f"<Header>{head}</Header><Body>{body}</Body></Table></Comment></Description>"
            f'<Attributes><Attribute attribute_name="host" harmonized_name="host">{host}</Attribute>'
            f"</Attributes></BioSample>")


def test_convert_and_exclude(tmp_path):
    xml_dir = tmp_path / "ast"
    xml_dir.mkdir()
    (xml_dir / "biosample_000000.xml").write_text(
        '<?xml version="1.0"?><BioSampleSet>'
        + sample("SAMN1", [
            ("piperacillin-tazobactam", "resistant", "&gt;", "64", "mg/L", "MIC", "Vitek", "", "", "EUCAST"),
            ("gentamicin", "susceptible", "==", "1", "mg/L", "MIC", "Vitek", "", "", "CLSI"),
            ("aztreonam", "resistant", "==", "15", "mm", "disk diffusion", "", "", "", "missing"),
            ("imipenem", "not defined", "&lt;=", "0.25", "mg/L", "MIC", "", "", "", "missing"),
        ])
        + sample("SAMN2", [("amikacin", "susceptible", "&lt;=", "8", "mg/L", "MIC", "", "", "", "CLSI")])
        + "</BioSampleSet>", encoding="utf-8")
    amr = tmp_path / "bv_amr.csv"
    amr.write_text("Genome ID,Antibiotic\n573.1,amikacin\n", encoding="utf-8")
    genomes = tmp_path / "bv_genomes.csv"
    genomes.write_text("genome_id,biosample_accession\n573.1,SAMN2\n573.2,SAMN9\n", encoding="utf-8")
    out, gout, summ = tmp_path / "out.csv", tmp_path / "g.csv", tmp_path / "s.txt"
    assert ncbi_ast.main([str(xml_dir), "--exclude-bvbrc", str(amr), str(genomes), "-o", str(out),
                          "--genomes-out", str(gout), "--summary", str(summ)]) == 0

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert {r["Genome ID"] for r in rows} == {"SAMN1"}              # SAMN2 already in BV-BRC
    pt = next(r for r in rows if r["Antibiotic"] == "piperacillin/tazobactam")
    assert (pt["Measurement"], pt["Measurement Sign"], pt["Resistant Phenotype"], pt["Testing Standard"]) == \
        (">64", ">", "Resistant", "EUCAST")
    gm = next(r for r in rows if r["Antibiotic"] == "gentamicin")
    assert (gm["Measurement"], gm["Measurement Sign"]) == ("1", "=")
    imi = next(r for r in rows if r["Antibiotic"] == "imipenem")
    assert (imi["Resistant Phenotype"], imi["Testing Standard"]) == ("", "")

    # the converted file goes through the unchanged BV-BRC reader
    imp = read_bvbrc_csv(out)
    assert sorted(r.antibiotic for r in imp.records) == ["gentamicin", "imipenem", "piperacillin/tazobactam"]
    assert sum(imp.excluded.values()) == 1                           # the zone diameter

    assert list(csv.reader(gout.open(encoding="utf-8")))[1:] == [["SAMN1", "SAMN1", "Homo sapiens"]]
    text = summ.read_text()
    assert "BioSamples left out, already with AST in BV-BRC: 1" in text
    assert "BioSamples kept: 1" in text


def test_antibiotic_names():
    assert ncbi_ast.antibiotic_name("Trimethoprim-sulfamethoxazole") == "trimethoprim/sulfamethoxazole"
    assert ncbi_ast.antibiotic_name("amoxicillin-clavulanic acid") == "amoxicillin/clavulanic acid"
    assert ncbi_ast.antibiotic_name("meropenem") == "meropenem"
