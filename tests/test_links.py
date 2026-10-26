"""Links tests; agents and values are synthetic fixtures."""

import pytest

from bpdrift.eucast import BPValue, EucastTable, MICRow
from bpdrift.links import LinkError, load_links, resolve

HEADER = "bvbrc_antibiotic,eucast_version,eucast_agent,restriction,note\n"


def write(tmp_path, body, header=HEADER):
    p = tmp_path / "links.csv"
    p.write_text("# comment\n" + header + body, encoding="utf-8")
    return p


def mk(agent, sheet="Testobacterales", subgroup=""):
    v = BPValue("value", 1.0, raw="1")
    return MICRow(sheet, "Class", agent, (), v, v, "", 10, subgroup)


TABLES = {
    "10.0": EucastTable("10.0", None, [mk("Alphacillin"), mk("Betacillin (systemic)")]),
    "16.1": EucastTable("16.1", None, [mk("Alphacillin (other)"), mk("Betacillin (systemic)")]),
}


def test_load_and_get(tmp_path):
    links = load_links(write(tmp_path,
        "alphacillin,10.0,Alphacillin,,\n"
        "alphacillin,16.1,Alphacillin (other),,indication split\n"
        "beta/cillin,10.0,,,agent not in this version\n"
        "beta/cillin,16.1,Betacillin (systemic),urinary,only row\n"))
    assert links.get("ALPHACILLIN", "16.1").note == "indication split"
    assert links.get("beta/cillin", "10.0").agent == ""
    assert links.get("beta/cillin", "16.1").restriction == "urinary"
    assert links.get("gammacillin", "10.0") is None
    assert links.versions == {"10.0", "16.1"}


def test_resolve_finds_rows_and_keeps_deliberate_gaps(tmp_path):
    links = load_links(write(tmp_path,
        "alphacillin,10.0,Alphacillin,,\nalphacillin,16.1,alphacillin (OTHER),,\n"
        "beta,10.0,,,none\n"))
    found = resolve(links, TABLES, "Testobacterales")
    assert found[("alphacillin", "16.1")].agent == "Alphacillin (other)"
    assert found[("beta", "10.0")] is None


@pytest.mark.parametrize(
    "body, message",
    [
        ("alphacillin,16.1,Alphacillin,,\n", "no 'Alphacillin' in Testobacterales of v16.1"),
        ("alphacillin,13.1,Alphacillin,,\n", "no table loaded for v13.1"),
    ],
)
def test_resolve_errors(tmp_path, body, message):
    with pytest.raises(LinkError, match=message):
        resolve(load_links(write(tmp_path, body)), TABLES, "Testobacterales")


def test_ambiguous_row(tmp_path):
    tables = {"10.0": EucastTable("10.0", None, [mk("X", subgroup="A"), mk("X", subgroup="B")])}
    with pytest.raises(LinkError, match="ambiguous"):
        resolve(load_links(write(tmp_path, "x,10.0,X,,\n")), tables, "Testobacterales")


@pytest.mark.parametrize(
    "body, message",
    [
        ("a,10.0,A,,\na,10.0,B,,\n", "already linked on line 3"),
        ("a,10.0,A,oral,\n", "restriction must be"),
        ("a,10.0,,urinary,\n", "a restriction needs an agent"),
        (",10.0,A,,\n", "antibiotic and version are required"),
        ("a,10.0,A\n", "3 fields"),
    ],
)
def test_bad_links(tmp_path, body, message):
    with pytest.raises(LinkError, match=message):
        load_links(write(tmp_path, body))


def test_bad_header(tmp_path):
    with pytest.raises(LinkError, match="header must be"):
        load_links(write(tmp_path, "", header="antibiotic,version\n"))
