"""End-to-end test of the H3 script on synthetic inputs."""

import csv
import importlib.util
from pathlib import Path

from bpdrift.reinterpret import COLUMNS as RE_COLUMNS

SCRIPT = Path(__file__).resolve().parents[1] / "analysis" / "genotype_errors.py"
spec = importlib.util.spec_from_file_location("genotype_errors", SCRIPT)
genotype_errors = importlib.util.module_from_spec(spec)
spec.loader.exec_module(genotype_errors)


def write_csv(path, header, rows, delimiter=","):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter=delimiter)
        w.writerow(header)
        w.writerows(rows)
    return path


def test_end_to_end(tmp_path):
    def re_row(genome, line, version, category, exclusion=""):
        d = dict.fromkeys(RE_COLUMNS, "")
        d.update(genome_id=genome, source_line=line, antibiotic="meropenem",
                 eucast_version=version, category=category, exclusion=exclusion)
        return [d[c] for c in RE_COLUMNS]

    re = write_csv(tmp_path / "re.csv", RE_COLUMNS, [
        # g1: KPC, R in both versions -> correct
        re_row("g1", 1, "10.0", "R"), re_row("g1", 1, "16.1", "R"),
        # g2: no carbapenemase; I under 10.0, R under 16.1 -> VME appears only in 16.1
        re_row("g2", 2, "10.0", "I"), re_row("g2", 2, "16.1", "R"),
        # g3: no BioSample match -> lost
        re_row("g3", 3, "10.0", "S"), re_row("g3", 3, "16.1", "S"),
        # g4: undecidable under 10.0
        re_row("g4", 4, "10.0", "UNDECIDABLE"), re_row("g4", 4, "16.1", "S"),
        # g5: excluded from main analysis -> ignored
        re_row("g5", 5, "10.0", "R", "urinary_only"), re_row("g5", 5, "16.1", "R", "urinary_only"),
    ])
    genomes = write_csv(tmp_path / "genomes.csv", ["Genome ID", "Genome Name", "BioSample Accession"],
                        [["g1", "K", "SAMN1"], ["g2", "K", "SAMN2"], ["g4", "K", "SAMN4"], ["g3", "K", ""]])
    pd = write_csv(tmp_path / "pd.tsv", ["#label", "biosample_acc", "AMR_genotypes"], [
        ["a", "SAMN1", "blaKPC-2=COMPLETE,blaSHV-1=COMPLETE"],
        ["b", "SAMN2", "blaSHV-1=COMPLETE,weird=COMPLETE"],
        ["c", "SAMN4", "blaSHV-1=COMPLETE"],
    ], delimiter="\t")
    genes = write_csv(tmp_path / "genes.txt",
                      ["allele", "gene_family", "class", "subclass"],
                      [["blaKPC-2", "blaKPC", "BETA-LACTAM", "CARBAPENEM"],
                       ["blaSHV-1", "blaSHV", "BETA-LACTAM", "BETA-LACTAM"]], delimiter="\t")
    smap = tmp_path / "map.csv"
    smap.write_text("subclass_token,antibiotic,note\nCARBAPENEM,meropenem,\nBETA-LACTAM,ampicillin,\n",
                    encoding="utf-8")
    out, summary = tmp_path / "out.csv", tmp_path / "summary.txt"
    assert genotype_errors.main([str(re), "--genomes", str(genomes), "--pd-metadata", str(pd),
                                 "--genes", str(genes), "--map", str(smap),
                                 "-o", str(out), "--summary", str(summary)]) == 0

    res = {(r["antibiotic"], r["eucast_version"]): r for r in csv.DictReader(out.open(encoding="utf-8"))}
    v10, v16 = res[("meropenem", "10.0")], res[("meropenem", "16.1")]
    assert (v10["compared"], v10["undecidable"], v10["phenotype_r"], v10["vme"], v10["me"]) == ("2", "1", "1", "0", "0")
    assert (v16["compared"], v16["phenotype_r"], v16["vme"], v16["me"]) == ("3", "2", "1", "0")
    assert res[("ALL", "16.1")]["vme"] == "1"

    text = summary.read_text()
    assert "genomes with main-analysis labels: 4" in text
    assert "genomes lost, no BioSample in BV-BRC: 1" in text
    assert "genomes with genotype predictions: 3" in text
    assert "weird: 1" in text

    # matched: g4 (undecidable under 10.0) leaves both versions, so they compare the same pairs
    out_m, summary_m = tmp_path / "out_m.csv", tmp_path / "summary_m.txt"
    assert genotype_errors.main([str(re), "--genomes", str(genomes), "--pd-metadata", str(pd),
                                 "--genes", str(genes), "--map", str(smap), "--matched",
                                 "-o", str(out_m), "--summary", str(summary_m)]) == 0
    m = {(r["antibiotic"], r["eucast_version"]): r for r in csv.DictReader(out_m.open(encoding="utf-8"))}
    assert m[("meropenem", "10.0")]["compared"] == m[("meropenem", "16.1")]["compared"] == "2"
    assert (m[("meropenem", "10.0")]["vme"], m[("meropenem", "16.1")]["vme"]) == ("0", "1")
    assert "lacking an S/I/R label under some version: 1" in summary_m.read_text()
