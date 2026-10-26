"""Re-interpretation tests on synthetic tables, links and records."""

import pytest

from bpdrift.bvbrc import ASTRecord
from bpdrift.eucast import BPValue, EucastTable, MICRow
from bpdrift.interpret import Category
from bpdrift.links import load_links
from bpdrift.mic import parse_mic
from bpdrift.reinterpret import COLUMNS, as_row, reinterpret

S, I, R, U, NB = Category.S, Category.I, Category.R, Category.UNDECIDABLE, Category.NO_BREAKPOINT


def bp(agent, s, r):
    return MICRow("Testobacterales", "Class", agent, (),
                  BPValue("value", float(s), raw=f"{s:g}"), BPValue("value", float(r), raw=f"{r:g}"), "", 10)


TABLES = {
    "10.0": EucastTable("10.0", None, [bp("Alphacillin", 8, 16), bp("Uricillin (UTI)", 2, 2)]),
    "16.1": EucastTable("16.1", None, [bp("Alphacillin", 8, 8), bp("Uricillin (UTI)", 2, 2),
                                         bp("Newcillin", 1, 1)]),
}


@pytest.fixture
def links(tmp_path):
    p = tmp_path / "links.csv"
    p.write_text(
        "bvbrc_antibiotic,eucast_version,eucast_agent,restriction,note\n"
        "alphacillin,10.0,Alphacillin,,\n"
        "alphacillin,16.1,Alphacillin,,\n"
        "uricillin,10.0,Uricillin (UTI),urinary,\n"
        "uricillin,16.1,Uricillin (UTI),urinary,\n"
        "newcillin,10.0,,,agent not in this version\n"
        "newcillin,16.1,Newcillin,,\n"
        "amoxicillin/clavulanic acid,10.0,Alphacillin,,\n"
        "amoxicillin/clavulanic acid,16.1,Alphacillin,,\n",
        encoding="utf-8",
    )
    return load_links(p)


def rec(ab="alphacillin", mic="16", genome="573.1", line=2, **flags):
    return ASTRecord(line, genome, "K. pneumoniae", ab, parse_mic(mic), "Intermediate",
                     "Broth dilution", "", "", "CLSI", "2017", "", mic, "", mic, **flags)


def cats(rows):
    return [(x.version, x.category) for x in rows]


def test_one_row_per_version_with_the_drift(links):
    rows, unlinked = reinterpret([rec(mic="16")], TABLES, links, "Testobacterales")
    assert cats(rows) == [("10.0", I), ("16.1", R)]
    assert not unlinked
    assert rows[1].agent == "Alphacillin" and rows[1].s == "8" and rows[1].r == "8"
    assert all(x.in_main for x in rows)


def test_undecidable_lists_possible(links):
    rows, _ = reinterpret([rec(mic=">=8")], TABLES, links, "Testobacterales")
    assert [(x.category, x.possible) for x in rows] == [(U, "S|I|R"), (U, "S|R")]


def test_deliberate_gap_is_no_breakpoint(links):
    rows, _ = reinterpret([rec(ab="newcillin", mic="0.5")], TABLES, links, "Testobacterales")
    assert cats(rows) == [("10.0", NB), ("16.1", S)]
    assert rows[0].reason == "agent not in this version"


def test_unlinked_antibiotics_are_counted_not_interpreted(links):
    rows, unlinked = reinterpret([rec(ab="polymyxin B"), rec(ab="polymyxin B", genome="573.2")],
                                 TABLES, links, "Testobacterales")
    assert rows == [] and unlinked == {"polymyxin B": 2}


@pytest.mark.parametrize(
    "record, reason",
    [
        (rec(ab="uricillin", mic="1"), "urinary_only"),
        (rec(ab="amoxicillin/clavulanic acid", ratio_format=True), "ratio_method"),
        (rec(ab="alphacillin", ratio_format=True), ""),           # ratio is fine elsewhere
        (rec(conflicting_values=True), "conflicting_values"),
        (rec(sign_conflict=True), ""),                             # kept, flagged
    ],
)
def test_main_analysis_exclusions(links, record, reason):
    rows, _ = reinterpret([record], TABLES, links, "Testobacterales")
    assert {x.exclusion for x in rows} == {reason}


def test_partially_linked_antibiotic_fails_loudly(tmp_path):
    p = tmp_path / "l.csv"
    p.write_text("bvbrc_antibiotic,eucast_version,eucast_agent,restriction,note\n"
                 "alphacillin,10.0,Alphacillin,,\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not v16.1"):
        reinterpret([rec()], TABLES, load_links(p), "Testobacterales")


def test_csv_row_shape(links):
    rows, _ = reinterpret([rec(mic="<=0.5", sign_conflict=True)], TABLES, links, "Testobacterales")
    out = as_row(rows[0])
    assert len(out) == len(COLUMNS)
    d = dict(zip(COLUMNS, out))
    assert (d["mic"], d["category"], d["eucast_version"], d["sign_conflict"]) == ("<=0.5", "S", "10.0", 1)
