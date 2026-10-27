"""Genotype predictor and error-rate tests on synthetic catalogues."""

import pytest

from bpdrift.genotype import (
    GenotypeError,
    load_exclusions,
    load_excluded_families,
    collapse_labels,
    error_rates,
    load_catalogue,
    load_subclass_map,
    parse_amr_genotypes,
    predict,
)


@pytest.fixture
def catalogue(tmp_path):
    genes = tmp_path / "ReferenceGeneCatalog.txt"
    genes.write_text(
        "allele\tgene_family\twhitelisted_taxa\tproduct_name\tscope\ttype\tsubtype\tclass\tsubclass\n"
        "blaKPC-2\tblaKPC\t\tcarbapenemase\tcore\tAMR\tAMR\tBETA-LACTAM\tCARBAPENEM\n"
        "\tblaCTX-M\t\tESBL\tcore\tAMR\tAMR\tBETA-LACTAM\tCEPHALOSPORIN\n"
        "aac(6')-Ib\taac(6')-Ib\t\taminoglycoside acetyltransferase\tcore\tAMR\tAMR\tAMINOGLYCOSIDE\tAMIKACIN/KANAMYCIN/TOBRAMYCIN\n"
        "dfrA12\tdfrA\t\tDHFR\tcore\tAMR\tAMR\tTRIMETHOPRIM\tTRIMETHOPRIM\n"
        "sul1\tsul1\t\tDHPS\tcore\tAMR\tAMR\tSULFONAMIDE\tSULFONAMIDE\n"
        "qacE\tqacE\t\tefflux\tplus\tSTRESS\tBIOCIDE\tQUATERNARY AMMONIUM\t\n"
        "\toqxA\t\tefflux\tcore\tAMR\tAMR\tPHENICOL\tNITROFURANTOIN/PHENICOL/QUINOLONE/TIGECYCLINE\n",
        encoding="utf-8",
    )
    muts = tmp_path / "AMRProt-mutation.tsv"
    muts.write_text(
        "#taxgroup\taccession_version\tmutation_position\tstandard_mutation_symbol\treported_mutation_symbol\tclass\tsubclass\tmutated_protein_name\n"
        "Klebsiella_pneumoniae\tX\t83\tgyrA_S83I\tgyrA_S83I\tQUINOLONE\tQUINOLONE\tgyrase\n"
        "Klebsiella_pneumoniae\tZ\t135\tompK36_D135DGD\tompK36_D135DGD\tBETA-LACTAM\tCARBAPENEM\tporin\n"
        "Escherichia\tY\t83\tgyrA_S83L\tgyrA_S83L\tQUINOLONE\tQUINOLONE\tgyrase\n",
        encoding="utf-8",
    )
    return load_catalogue(genes, muts, version="test")


@pytest.fixture
def smap(tmp_path):
    p = tmp_path / "map.csv"
    p.write_text(
        "# synthetic\nsubclass_token,antibiotic,note\n"
        "CARBAPENEM,meropenem,\nCARBAPENEM,cefotaxime,\nCEPHALOSPORIN,cefotaxime,\n"
        "AMIKACIN,amikacin,\nTOBRAMYCIN,tobramycin,\nQUINOLONE,ciprofloxacin,\n"
        "TRIMETHOPRIM,trimethoprim/sulfamethoxazole,\n",
        encoding="utf-8",
    )
    return load_subclass_map(p)


@pytest.mark.parametrize(
    "text, expected",
    [
        ("blaKPC-2=COMPLETE,gyrA_S83I=POINT", [("blaKPC-2", "COMPLETE"), ("gyrA_S83I", "POINT")]),
        ("blaCTX-M=PARTIAL", [("blaCTX-M", "PARTIAL")]),
        ("sul1", [("sul1", "")]),
        ("", []), ("NULL", []), ("  ", []),
        ("aac(6')-Ib=COMPLETE, ,dfrA12=COMPLETE", [("aac(6')-Ib", "COMPLETE"), ("dfrA12", "COMPLETE")]),
    ],
)
def test_parse(text, expected):
    assert parse_amr_genotypes(text) == expected


def test_parse_rejects_empty_symbol():
    with pytest.raises(GenotypeError):
        parse_amr_genotypes("=COMPLETE")


def test_catalogue_keys_alleles_families_and_organism_mutations(catalogue):
    assert catalogue.lookup("blaKPC-2") == "CARBAPENEM"
    assert catalogue.lookup("blaKPC") == "CARBAPENEM"
    assert catalogue.lookup("blaCTX-M") == "CEPHALOSPORIN"
    assert catalogue.lookup("gyrA_S83I") == "QUINOLONE"
    assert catalogue.lookup("gyrA_S83L") is None      # E. coli mutation not loaded
    assert catalogue.lookup("qacE") is None           # no subclass: not an antibiotic determinant


def test_compound_subclass_and_hierarchy(catalogue, smap):
    p = predict(parse_amr_genotypes("aac(6')-Ib=COMPLETE,blaKPC-2=COMPLETE,mystery=COMPLETE"), catalogue, smap)
    assert p.resistant == {"amikacin", "tobramycin", "meropenem", "cefotaxime"}
    assert p.unknown == ["mystery"]


def test_sulfonamide_alone_predicts_nothing(catalogue, smap):
    assert predict([("sul1", "COMPLETE")], catalogue, smap).resistant == set()
    assert "trimethoprim/sulfamethoxazole" in predict([("dfrA12", "COMPLETE")], catalogue, smap).resistant


def test_error_rates():
    pairs = [(True, "R"), (False, "R"), (False, "R"), (True, "S"), (False, "I"), (False, "S")]
    e = error_rates(pairs)
    assert (e.n, e.phenotype_r, e.phenotype_not_r, e.vme, e.me) == (6, 3, 3, 2, 1)
    assert e.vme_rate[0] == pytest.approx(2 / 3)
    assert e.me_rate[0] == pytest.approx(1 / 3)
    assert e.agreement[0] == pytest.approx(3 / 6)


def test_error_rates_reject_undecidable():
    with pytest.raises(GenotypeError):
        error_rates([(True, "UNDECIDABLE")])


def test_collapse_labels_keeps_only_consistent():
    rows = [("g1", "Meropenem", "16.1", "R"), ("g1", "meropenem", "16.1", "R"),
            ("g2", "meropenem", "16.1", "S"), ("g2", "meropenem", "16.1", "R"),
            ("g3", "meropenem", "10.0", "I")]
    labels, dropped = collapse_labels(rows)
    assert labels == {("g1", "meropenem", "16.1"): "R", ("g3", "meropenem", "10.0"): "I"}
    assert dropped == 1


def test_bad_subclass_map(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("token,drug\nA,b\n", encoding="utf-8")
    with pytest.raises(GenotypeError, match="header"):
        load_subclass_map(p)


def test_family_of(catalogue):
    assert catalogue.family_of("blaKPC-2") == "blaKPC"
    assert catalogue.family_of("ompK36_D135DGD") == "ompK36"
    assert catalogue.family_of("oqxA") == "oqxA"


def test_excluded_families_and_methods(catalogue, smap, tmp_path):
    ex = tmp_path / "ex.csv"
    ex.write_text("# intrinsic\nfamily,reason\noqxA,intrinsic efflux\nompK36,porin\n", encoding="utf-8")
    fams = load_excluded_families(ex)
    assert fams == {"oqxA": "intrinsic efflux", "ompK36": "porin"}
    calls = parse_amr_genotypes(
        "oqxA,ompK36_D135DGD=POINT,gyrA_S83I=POINT,blaKPC-2=MISTRANSLATION,dfrA12=PARTIAL_END_OF_CONTIG")
    p = predict(calls, catalogue, smap, excluded_families=fams)
    assert p.resistant == {"ciprofloxacin", "trimethoprim/sulfamethoxazole"}
    assert sorted(p.ignored) == ["blaKPC-2=MISTRANSLATION", "ompK36_D135DGD", "oqxA"]
    # without the exclusions the porin mutation and efflux gene would predict resistance
    loose = predict(calls, catalogue, smap)
    assert {"meropenem", "ciprofloxacin"} <= loose.resistant


def test_mutation_catalogue_must_load_something(tmp_path):
    genes = tmp_path / "g.txt"
    genes.write_text("allele\tgene_family\tclass\tsubclass\nx\tx\tA\tB\n", encoding="utf-8")
    muts = tmp_path / "m.tsv"
    muts.write_text("#taxgroup\treported_mutation_symbol\tclass\tsubclass\nEscherichia\tgyrA_S83L\tQ\tQUINOLONE\n",
                    encoding="utf-8")
    with pytest.raises(GenotypeError, match="no point mutations"):
        load_catalogue(genes, muts)


def test_disruption_fallback(tmp_path):
    genes = tmp_path / "g.txt"
    genes.write_text(
        "allele\tgene_family\twhitelisted_taxa\tsubtype\tclass\tsubclass\n"
        "\tmgrB\tKlebsiella_pneumoniae\tPOINT_DISRUPT\tCOLISTIN\tCOLISTIN\n"
        "\tacrR\tEscherichia\tPOINT_DISRUPT\tEFFLUX\tFLUOROQUINOLONE\n"
        "mcr-1.1\tmcr-1\t\tAMR\tCOLISTIN\tCOLISTIN\n",
        encoding="utf-8",
    )
    cat = load_catalogue(genes)
    assert cat.lookup("mgrB_K43Ter") == "COLISTIN"
    assert cat.lookup("mgrB_G37YfsTer1") == "COLISTIN"
    assert cat.lookup("mgrB") is None              # the gene itself is not an element
    assert cat.lookup("mgrB_K43N") is None         # a missense change is not a disruption
    assert cat.lookup("acrR_W45Ter") is None       # disruption entry is for E. coli only
    assert cat.lookup("mcr-1.1") == "COLISTIN"


def test_symbol_level_exclusion(tmp_path):
    genes = tmp_path / "g.txt"
    genes.write_text(
        "allele\tgene_family\tclass\tsubclass\n"
        "\tfosA\tFOSFOMYCIN\tFOSFOMYCIN\n"
        "fosA3\tfosA\tFOSFOMYCIN\tFOSFOMYCIN\n"
        "fosA5\tfosA\tFOSFOMYCIN\tFOSFOMYCIN\n",
        encoding="utf-8",
    )
    cat = load_catalogue(genes)
    m = tmp_path / "m.csv"
    m.write_text("subclass_token,antibiotic\nFOSFOMYCIN,fosfomycin\n", encoding="utf-8")
    smap = load_subclass_map(m)
    ex_csv = tmp_path / "ex.csv"
    ex_csv.write_text("family,reason,match\nfosA,chromosomal,symbol\nfosA5,chromosomal,symbol\n"
                      "oqxA,efflux,family\noqxB,efflux,\n", encoding="utf-8")
    ex = load_exclusions(ex_csv)
    assert set(ex.symbols) == {"fosA", "fosA5"} and set(ex.families) == {"oqxA", "oqxB"}
    assert load_excluded_families(ex_csv) == ex.families

    def pred(text):
        return predict(parse_amr_genotypes(text), cat, smap, excluded_families=ex.families,
                       excluded_symbols=ex.symbols).resistant

    assert pred("fosA=EXACT") == set()                       # chromosomal copy ignored
    assert pred("fosA5=EXACT") == set()
    assert pred("fosA=EXACT,fosA3=EXACT") == {"fosfomycin"}  # acquired gene of the same family counts
    # excluding the whole family would also lose the acquired gene
    assert predict(parse_amr_genotypes("fosA3=EXACT"), cat, smap, excluded_families={"fosA"}).resistant == set()


def test_exclusion_match_must_be_known(tmp_path):
    p = tmp_path / "ex.csv"
    p.write_text("family,reason,match\nfosA,x,allele\n", encoding="utf-8")
    with pytest.raises(GenotypeError, match="match must be"):
        load_exclusions(p)


def test_intermediate_phenotype_counts_as_not_resistant():
    """EUCAST I is a susceptible category: predicted R with phenotype I is a major error,
    predicted not R with phenotype I is correct, and I never enters the VME denominator."""
    e = error_rates([(True, "I"), (False, "I"), (True, "R")])
    assert (e.phenotype_r, e.phenotype_not_r, e.vme, e.me) == (1, 2, 0, 1)
