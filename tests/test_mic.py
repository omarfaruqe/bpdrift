import pytest

from bpdrift import MIC, MICParseError, parse_mic


@pytest.mark.parametrize(
    "text, value, sign",
    [
        ("0.5", 0.5, "="),
        (".25", 0.25, "="),
        ("<=0.25", 0.25, "<="),
        ("≤0.25", 0.25, "<="),
        ("=<0.25", 0.25, "<="),
        ("< 0.25", 0.25, "<"),
        (">16", 16.0, ">"),
        (">=32", 32.0, ">="),
        ("≥ 32 mg/L", 32.0, ">="),
        ("8 ug/mL", 8.0, "="),
        ("8µg/mL", 8.0, "="),
        ("  2  ", 2.0, "="),
        ("==4", 4.0, "="),
    ],
)
def test_parse_valid(text, value, sign):
    mic = parse_mic(text)
    assert mic.value == value
    assert mic.sign == sign


def test_parse_numbers():
    assert parse_mic(4) == MIC(4.0)
    assert parse_mic(0.125) == MIC(0.125)


@pytest.mark.parametrize(
    "text",
    ["", "   ", "abc", "<", "<=-1", "0", "1/2", "4 mm", "S", "8 g/L", "1e3", "<= 2 <= 4"],
)
def test_parse_invalid_raises(text):
    with pytest.raises(MICParseError):
        parse_mic(text)


@pytest.mark.parametrize("bad", [True, None, [1], float("nan"), float("inf"), -2])
def test_parse_bad_types_raise(bad):
    with pytest.raises(MICParseError):
        parse_mic(bad)


def test_censoring_flag_and_str():
    assert not parse_mic("2").is_censored
    assert parse_mic(">2").is_censored
    assert str(parse_mic("<=0.25")) == "<=0.25"
    assert str(parse_mic("4")) == "4"


def test_unknown_sign_rejected():
    with pytest.raises(MICParseError):
        MIC(1.0, "~")
