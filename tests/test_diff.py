"""Diff tests on hand-built tables. All values are synthetic fixtures."""

import pytest

from bpdrift.diff import (
    ChangeKind,
    DiffError,
    base_name,
    diff_tables,
    load_mapping,
    norm,
    related,
)
from bpdrift.eucast import BPValue, EucastTable, MICRow
from bpdrift.interpret import Category

S, I, R = Category.S, Category.I, Category.R


def val(x):
    if isinstance(x, str):
        return BPValue({"-": "none", "IE": "IE"}[x], raw=x)
    return BPValue("value", float(x), raw=f"{x:g}")


def mk(agent, s, r, sheet="Testobacterales", subgroup="", row=10):
    return MICRow(sheet, "Class", agent, (), val(s), val(r), "", row, subgroup)


def table(version, *rows):
    return EucastTable(version=version, pin=None, rows=list(rows))


def one(diff, agent):
    (c,) = [c for c in diff.changes if c.agent == agent]
    return c


def test_names():
    assert norm("  Piperacillin–tazobactam ") == "piperacillin-tazobactam"
    assert base_name("Ciprofloxacin (meningitis)") == "ciprofloxacin"
    assert base_name("Imipenem, Morganellaceae") == "imipenem"


@pytest.mark.parametrize(
    "a, b, expected",
    [
        ("Amoxicillin", "Amoxicillin iv", True),
        ("Amoxicillin", "Amoxicillin oral (other indications)", True),
        ("Ciprofloxacin", "Ciprofloxacin (meningitis)", True),
        ("Imipenem", "Imipenem, Morganellaceae", True),
        ("Amoxicillin", "Amoxicillin-clavulanic acid", False),
        ("Imipenem", "Imipenem-relebactam", False),
        ("Cefotaxime", "Ceftriaxone", False),
    ],
)
def test_related(a, b, expected):
    assert related(a, b) is expected
    assert related(b, a) is expected


def test_unchanged_changed_added_removed():
    old = table("10.0", mk("Alphacillin", 8, 16), mk("Betacillin", 1, 2), mk("Gone", 1, 1))
    new = table("13.1", mk("ALPHACILLIN", 8, 8), mk("Betacillin", 1, 2), mk("Fresh", 2, 2))
    d = diff_tables(old, new)
    assert d.summary() == {"unchanged": 1, "changed": 1, "added": 1, "removed": 1}
    assert one(d, "ALPHACILLIN").match == "exact"
    assert one(d, "Gone").kind is ChangeKind.REMOVED
    assert one(d, "Fresh").kind is ChangeKind.ADDED


def test_flips_list_exact_dilutions():
    # R breakpoint lowered from >16 to >8: an MIC of 16 moves from I to R
    d = diff_tables(table("10.0", mk("Alphacillin", 8, 16)), table("13.1", mk("Alphacillin", 8, 8)))
    c = one(d, "Alphacillin")
    assert c.kind is ChangeKind.CHANGED
    assert [(f.mic, f.old, f.new) for f in c.flips] == [(16, I, R)]


def test_lowered_s_breakpoint_flips_several_dilutions():
    d = diff_tables(table("a", mk("X", 2, 8)), table("b", mk("X", 0.5, 8)))
    assert [(f.mic, f.old, f.new) for f in one(d, "X").flips] == [(1, S, I), (2, S, I)]


def test_breakpoint_appearing_is_a_change_without_flips():
    d = diff_tables(table("a", mk("X", "IE", "IE")), table("b", mk("X", 1, 1)))
    c = one(d, "X")
    assert c.kind is ChangeKind.CHANGED and c.flips == ()


def test_bracketing_alone_counts_as_a_change():
    old = table("a", mk("X", 2, 2))
    b = MICRow("Testobacterales", "Class", "X", (), BPValue("value", 2.0, True), BPValue("value", 2.0, True), "", 10)
    c = one(diff_tables(old, table("b", b)), "X")
    assert c.kind is ChangeKind.CHANGED and c.flips == ()


def test_candidates_suggested_but_never_paired():
    old = table("10.0", mk("Ciprofloxacin", 0.25, 0.5))
    new = table("13.1", mk("Ciprofloxacin (meningitis)", 0.125, 0.125, row=11),
                mk("Ciprofloxacin (indications other than meningitis)", 0.25, 0.5, row=12))
    d = diff_tables(old, new)
    removed = one(d, "Ciprofloxacin")
    assert removed.kind is ChangeKind.REMOVED
    assert set(removed.candidates) == {"Ciprofloxacin (meningitis)",
                                       "Ciprofloxacin (indications other than meningitis)"}
    assert one(d, "Ciprofloxacin (meningitis)").candidates == ("Ciprofloxacin",)


def write_map(tmp_path, *lines):
    p = tmp_path / "map.csv"
    p.write_text("old_group,old_agent,new_group,new_agent\n" + "\n".join(lines) + "\n", encoding="utf-8")
    return p


def test_mapping_records_a_split(tmp_path):
    old = table("10.0", mk("Ciprofloxacin", 0.25, 0.5))
    new = table("13.1", mk("Ciprofloxacin (meningitis)", 0.125, 0.125, row=11),
                mk("Ciprofloxacin (indications other than meningitis)", 0.25, 0.5, row=12))
    m = load_mapping(write_map(
        tmp_path,
        "Testobacterales,Ciprofloxacin,Testobacterales,Ciprofloxacin (meningitis)",
        "Testobacterales,Ciprofloxacin,Testobacterales,Ciprofloxacin (indications other than meningitis)",
    ))
    d = diff_tables(old, new, m)
    assert d.summary() == {"unchanged": 1, "changed": 1, "added": 0, "removed": 0}
    assert {c.match for c in d.changes} == {"mapped"}
    assert one(d, "Ciprofloxacin (meningitis)").flips  # 0.25 was S, now R


def test_group_rename(tmp_path):
    old = table("10.0", mk("Alphacillin", 1, 1, sheet="P.testida"))
    new = table("13.1", mk("Alphacillin", 1, 1, sheet="Testurella"))
    assert diff_tables(old, new).summary()["removed"] == 1
    m = load_mapping(write_map(tmp_path, "P.testida,,Testurella,"))
    d = diff_tables(old, new, m)
    assert d.summary() == {"unchanged": 1, "changed": 0, "added": 0, "removed": 0}


def test_subgroups_match_across_reorganised_sheets():
    # an organism subgroup moving to a new sheet still matches by group
    old = table("15.0", mk("Alphacillin", 2, 2, sheet="Anaerobes", subgroup="Alphoides spp."))
    new = table("16.1", mk("Alphacillin", 1, 1, sheet="Anaerobes G-neg", subgroup="Alphoides spp."))
    c = one(diff_tables(old, new), "Alphacillin")
    assert c.kind is ChangeKind.CHANGED and c.match == "exact"


def test_same_agent_in_two_subgroups_is_not_ambiguous():
    rows = [mk("Alphacillin", 2, 2, subgroup="A spp.", row=10), mk("Alphacillin", 1, 1, subgroup="B spp.", row=20)]
    d = diff_tables(table("a", *rows), table("b", *rows))
    assert d.summary()["unchanged"] == 2


def test_duplicate_rows_fail_loudly():
    dup = table("a", mk("X", 1, 1, row=10), mk("x", 2, 2, row=11))
    with pytest.raises(DiffError, match="twice"):
        diff_tables(dup, table("b", mk("X", 1, 1)))


def test_group_filter_by_sheet_or_subgroup():
    old = table("a", mk("X", 1, 1), mk("Y", 1, 1, sheet="Other"))
    new = table("b", mk("X", 2, 2), mk("Y", 1, 1, sheet="Other"))
    d = diff_tables(old, new, groups=["testobacterales"])
    assert [c.agent for c in d.changes] == ["X"]


@pytest.mark.parametrize(
    "line, message",
    [
        ("Testobacterales,Nope,Testobacterales,X", "no agent 'Nope' in old group"),
        ("Testobacterales,X,Testobacterales,Nope", "no agent 'Nope' in new group"),
        ("Missing,,Testobacterales,", "not in the old table"),
    ],
)
def test_mapping_errors_name_the_line(tmp_path, line, message):
    m = load_mapping(write_map(tmp_path, line))
    with pytest.raises(DiffError, match=message):
        diff_tables(table("a", mk("X", 1, 1)), table("b", mk("X", 1, 1)), m)


def test_mapping_note_column_and_comments(tmp_path):
    p = tmp_path / "noted.csv"
    p.write_text(
        "# reviewed 2026-10-26\n"
        "old_group,old_agent,new_group,new_agent,note\n"
        "\n"
        "# carbapenems\n"
        'Testobacterales,X,Testobacterales,"X, subgroup",indication split; see v13 note 2\n'
        "Testobacterales,Y,Testobacterales,Y2\n",
        encoding="utf-8",
    )
    m = load_mapping(p)
    assert [(oa, na, line) for _, oa, _, na, line in m.pairs] == [
        ("X", "X, subgroup", 5), ("Y", "Y2", 6),
    ]


def test_mapping_error_reports_file_line_despite_comments(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("# c\nold_group,old_agent,new_group,new_agent\n# c\nG,X,G,\n", encoding="utf-8")
    with pytest.raises(DiffError, match="line 4"):
        load_mapping(p)


@pytest.mark.parametrize(
    "content, message",
    [
        ("a,b,c\n", "header must be"),
        ("old_group,old_agent,new_group,new_agent,extra\n", "header must be"),
        ("# only a comment\n", "empty mapping"),
        ("old_group,old_agent,new_group,new_agent\nG,X,G,Y,too,many\n", "6 fields"),
        ("old_group,old_agent,new_group,new_agent\nG,X,G,\n", "line 2: give both agents"),
        ("old_group,old_agent,new_group,new_agent\n,X,G,Y\n", "line 2: old_group and new_group"),
        ("old_group,old_agent,new_group,new_agent\nG,,H,\ng,,K,\n", "renamed twice"),
    ],
)
def test_bad_mapping_files(tmp_path, content, message):
    p = tmp_path / "bad.csv"
    p.write_text(content, encoding="utf-8")
    with pytest.raises(DiffError, match=message):
        load_mapping(p)
