"""BV-BRC import tests on a small synthetic export. Genome IDs are invented."""

import csv

import pytest

from bpdrift.bvbrc import BVBRCFormatError, read_bvbrc_csv

HEADER = [
    "Taxon ID", "Genome ID", "Genome Name", "Antibiotic", "Resistant Phenotype", "Measurement",
    "Measurement Sign", "Measurement Value", "Measurement Unit", "Laboratory Typing Method",
    "Laboratory Typing Method Version", "Laboratory Typing Platform", "Vendor", "Testing Standard",
    "Testing Standard Year", "Computational Method", "Computational Method Version",
    "Computational Method Performance", "Evidence", "Source", "PubMed",
]


def rec(genome="573.1", ab="meropenem", pheno="Susceptible", meas="0.5", sign="", value="0.5",
        unit="mg/L", method="Broth dilution", std="EUCAST", year="", evidence="Laboratory Method"):
    r = dict.fromkeys(HEADER, "")
    r.update({
        "Taxon ID": "573", "Genome ID": genome, "Genome Name": "Klebsiella pneumoniae test",
        "Antibiotic": ab, "Resistant Phenotype": pheno, "Measurement": meas,
        "Measurement Sign": sign, "Measurement Value": value, "Measurement Unit": unit,
        "Laboratory Typing Method": method, "Testing Standard": std,
        "Testing Standard Year": year, "Evidence": evidence,
    })
    return r


def write(tmp_path, rows, header=HEADER):
    p = tmp_path / "amr.csv"
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return p


def test_exclusions_are_counted(tmp_path):
    rows = [
        rec(),
        rec(),                                            # exact duplicate
        rec(genome="573.2", meas="", value=""),           # S/R call only
        rec(genome="573.3", unit="mm", value="22", meas="22", method="Disk diffusion"),
        rec(genome="573.4", evidence="Computational Prediction"),
        rec(genome="573.5", value="see note", meas="see note"),
        rec(genome="573.6", value="0", meas="0"),
    ]
    imp = read_bvbrc_csv(write(tmp_path, rows))
    s = imp.summary()
    assert s["lines"] == 7 and s["records"] == 1
    assert (s["excluded_duplicate"], s["excluded_no_value"], s["excluded_not_mic"],
            s["excluded_not_laboratory"], s["excluded_unparseable"]) == (1, 1, 1, 1, 2)


@pytest.mark.parametrize(
    "meas, sign, value, expected_sign, conflict",
    [
        ("<=1", "<=", "1", "<=", False),
        (">8", ">", "8", ">", False),
        ("32", ">", "32", ">", True),        # sign only in the sign field
        (">=16", "", "16", ">=", False),     # sign only in the measurement
        ("2", "=", "2", "=", False),
        ("2", "", "2", "=", False),
        ("<4", "<=", "4", "<=", True),       # both present and different
    ],
)
def test_sign_rules(tmp_path, meas, sign, value, expected_sign, conflict):
    (r,) = read_bvbrc_csv(write(tmp_path, [rec(meas=meas, sign=sign, value=value)])).records
    assert r.mic.sign == expected_sign
    assert r.mic.value == float(value)
    assert r.sign_conflict is conflict
    assert (r.raw_measurement, r.raw_sign, r.raw_value) == (meas, sign, value)


def test_combination_values_keep_first_component(tmp_path):
    rows = [
        rec(ab="trimethoprim/sulfamethoxazole", meas="16/304", sign=">=", value="16/304"),
        rec(genome="573.2", ab="piperacillin/tazobactam", meas="8", value="8"),
    ]
    a, b = read_bvbrc_csv(write(tmp_path, rows)).records
    assert (a.mic.value, a.mic.sign, a.ratio_format) == (16.0, ">=", True)
    assert (b.mic.value, b.ratio_format) == (8.0, False)


def test_conflicting_values_are_kept_and_flagged(tmp_path):
    rows = [
        rec(genome="573.1", value="1", meas="1"),
        rec(genome="573.1", value="2", meas="2"),
        rec(genome="573.1", ab="amikacin", value="4", meas="4"),
        rec(genome="573.9", value="1", meas="1"),
    ]
    imp = read_bvbrc_csv(write(tmp_path, rows))
    flags = {(r.genome_id, r.antibiotic, r.mic.value): r.conflicting_values for r in imp.records}
    assert flags == {
        ("573.1", "meropenem", 1.0): True,
        ("573.1", "meropenem", 2.0): True,
        ("573.1", "amikacin", 4.0): False,
        ("573.9", "meropenem", 1.0): False,
    }


def test_empty_unit_is_assumed_mg_per_l_and_flagged(tmp_path):
    (r,) = read_bvbrc_csv(write(tmp_path, [rec(unit="")])).records
    assert r.unit_assumed


def test_record_keeps_provenance(tmp_path):
    (r,) = read_bvbrc_csv(write(tmp_path, [rec(std="CLSI", year="2017")])).records
    assert (r.line, r.standard, r.standard_year, r.method) == (2, "CLSI", "2017", "Broth dilution")


def test_wrong_file_fails_loudly(tmp_path):
    p = tmp_path / "x.csv"
    p.write_text("a,b\n1,2\n", encoding="utf-8")
    with pytest.raises(BVBRCFormatError, match="missing columns"):
        read_bvbrc_csv(p)
