"""Parser tests on a synthetic workbook laid out like a EUCAST table.

Every agent name and value here is an invented test fixture, not EUCAST data.
"""

import pytest

openpyxl = pytest.importorskip("openpyxl")
from openpyxl.cell.rich_text import CellRichText, TextBlock  # noqa: E402
from openpyxl.cell.text import InlineFont  # noqa: E402

from bpdrift.eucast import EucastFormatError, read_eucast_xlsx  # noqa: E402
from bpdrift.interpret import Category, interpret  # noqa: E402
from bpdrift.pin import PinError, pin_file  # noqa: E402

SUP = InlineFont(vertAlign="superscript")


def sup(text, note):
    """Cell text followed by a superscript footnote marker."""
    return CellRichText(text, TextBlock(SUP, note))


def mic_block(ws, row, drug_class, agents):
    ws.cell(row, 1, drug_class)
    ws.cell(row, 2, "MIC breakpoints\n(mg/L)")
    ws.cell(row + 1, 2, "S ≤")
    ws.cell(row + 1, 3, "R >")
    ws.cell(row + 1, 4, "ATU")
    for i, (agent, s, r, atu) in enumerate(agents):
        ws.cell(row + 2 + i, 1, agent)
        ws.cell(row + 2 + i, 2, s)
        ws.cell(row + 2 + i, 3, r)
        if atu is not None:
            ws.cell(row + 2 + i, 4, atu)
    return row + 2 + len(agents)


def build(path, version="16.1", extra_value=None, header_r="R >"):
    wb = openpyxl.Workbook()
    wb.active.title = "Content"
    ws = wb.create_sheet("Testobacterales")
    ws.cell(1, 1, "Testobacterales")
    ws.cell(1, 9, f"EUCAST Clinical Breakpoint Tables v. {version}, valid from 2026-01-01")
    end = mic_block(ws, 8, sup("Alphacillins", "1"), [
        (sup("Alphacillin iv", "1"), "8", "8", None),
        ("Alphacillin-inhibitor", sup("8", "2"), sup("8", "2"), "16"),
        ("Alphacillin oral (other indications)", sup("(0.001)", "3"), sup("(8)", "3"), None),
        ("Betacillin", 0.5, 2, None),
        ("Gammacillin", "-", "-", None),
        ("Deltacillin", "IE", "IE", None),
        ("Epsicillin (screen only)", sup("Note", "5"), sup("Note", "5"), None),
        ("Zetacillin", "NA", "NA", None),
        ("Etacillin", "IP", "IP", None),
        (CellRichText("Thetacillin", TextBlock(SUP, "2"), ", subgroup X"), "2", "4", "4"),
    ])
    block2 = mic_block(ws, end + 2, "Omegamycins", [
        ("Omegamycin", extra_value if extra_value is not None else "1", "2", None),
    ])
    ws.cell(block2 - 2, 3, header_r)  # the R header of the second block

    pk = wb.create_sheet("PK PD breakpoints")
    pk.cell(1, 4, "EUCAST Clinical Breakpoint Tables v. 14.0")
    mic_block(pk, 6, "Penicillins", [("Alphacillin", "2", "8", None)])
    wb.save(path)
    return path


@pytest.fixture
def table(tmp_path):
    return read_eucast_xlsx(build(tmp_path / "synthetic.xlsx"))


def row(table, agent):
    (match,) = [r for r in table.rows if r.agent == agent]
    return match


def test_version_and_sheets(table):
    assert table.version == "16.1"
    assert table.sheets == ["Testobacterales"]  # PK-PD excluded by default
    assert len(table.rows) == 11


def test_superscript_notes_are_split_not_glued(table):
    r = row(table, "Alphacillin-inhibitor")
    assert r.s.kind == "value" and r.s.value == 8.0 and r.s.notes == ("2",)
    assert r.r.value == 8.0
    assert r.atu == "16"
    assert row(table, "Alphacillin iv").agent_notes == ("1",)
    assert row(table, "Alphacillin iv").drug_class == "Alphacillins"


def test_superscript_inside_agent_name(table):
    r = row(table, "Thetacillin, subgroup X")
    assert r.agent_notes == ("2",)
    assert (r.s.value, r.r.value) == (2.0, 4.0)


def test_bracketed_values(table):
    r = row(table, "Alphacillin oral (other indications)")
    assert r.s.bracketed and r.s.value == 0.001 and r.s.notes == ("3",)
    assert r.r.bracketed and r.r.value == 8.0


def test_numeric_cells(table):
    r = row(table, "Betacillin")
    assert (r.s.value, r.r.value) == (0.5, 2.0)


@pytest.mark.parametrize(
    "agent, kind",
    [("Gammacillin", "none"), ("Deltacillin", "IE"), ("Epsicillin (screen only)", "note"),
     ("Zetacillin", "NA"), ("Etacillin", "IP")],
)
def test_codes(table, agent, kind):
    r = row(table, agent)
    assert r.s.kind == kind and r.r.kind == kind
    assert interpret("1", r.breakpoint()).category is Category.NO_BREAKPOINT


def test_second_block_and_cell_reference(table):
    r = row(table, "Omegamycin")
    assert r.drug_class == "Omegamycins"
    assert r.cell == f"Testobacterales!A{r.row}"


def test_breakpoint_feeds_interpret(table):
    bp = row(table, "Betacillin").breakpoint()
    assert interpret("<=0.25", bp).category is Category.S
    assert interpret("1", bp).category is Category.I
    assert interpret(">=1", bp).category is Category.UNDECIDABLE


def test_include_pkpd(tmp_path):
    t = read_eucast_xlsx(build(tmp_path / "s.xlsx"), include_pkpd=True)
    assert "PK PD breakpoints" in t.sheets
    assert t.version == "16.1"  # the PK-PD sheet's own stale title is ignored


@pytest.mark.parametrize("bad", ["8 mg", "S", "0", "1/2"])
def test_unrecognised_value_fails_loudly(tmp_path, bad):
    path = build(tmp_path / "bad.xlsx", extra_value=bad)
    with pytest.raises(EucastFormatError, match=r"Testobacterales!B\d+"):
        read_eucast_xlsx(path)


def test_unexpected_header_fails_loudly(tmp_path):
    path = build(tmp_path / "hdr.xlsx", header_r="R ≥")
    with pytest.raises(EucastFormatError, match="expected 'S ≤' / 'R >'"):
        read_eucast_xlsx(path)


def test_pin_checks_checksum_and_version(tmp_path):
    path = build(tmp_path / "p.xlsx")
    good = pin_file(path, "16.1", "2026-10-26")
    assert read_eucast_xlsx(path, pin=good).pin == good

    wrong_version = pin_file(path, "15.0", "2026-10-26")
    with pytest.raises(EucastFormatError, match="pinned as v15.0"):
        read_eucast_xlsx(path, pin=wrong_version)

    other = build(tmp_path / "other.xlsx", version="13.1")
    with pytest.raises(PinError):
        read_eucast_xlsx(other, pin=good)


def test_not_a_eucast_table(tmp_path):
    path = tmp_path / "empty.xlsx"
    openpyxl.Workbook().save(path)
    with pytest.raises(EucastFormatError):
        read_eucast_xlsx(path)
