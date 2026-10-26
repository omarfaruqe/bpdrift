import datetime as dt
import hashlib

import pytest

from bpdrift.pin import Pin, PinError, load_pin, pin_file, save_pin, verify


@pytest.fixture
def table(tmp_path):
    p = tmp_path / "fake_breakpoints.xlsx"
    p.write_bytes(b"not a real table, test bytes only")
    return p


def test_pin_records_checksum_and_metadata(table):
    pin = pin_file(table, "16.1", dt.date(2026, 10, 26))
    assert pin.sha256 == hashlib.sha256(table.read_bytes()).hexdigest()
    assert pin.size_bytes == table.stat().st_size
    assert pin.filename == "fake_breakpoints.xlsx"
    assert pin.retrieved == "2026-10-26"
    assert pin.label.startswith("EUCAST v16.1 (2026-10-26, sha256:")


def test_verify_passes_and_detects_change(table):
    pin = pin_file(table, "16.1", "2026-10-26")
    verify(table, pin)
    table.write_bytes(b"edited")
    with pytest.raises(PinError, match="does not match pin"):
        verify(table, pin)


def test_round_trip(table, tmp_path):
    pin = pin_file(table, "15.0", "2025-01-01")
    out = tmp_path / "pin.json"
    save_pin(pin, out)
    assert load_pin(out) == pin


def test_missing_file(tmp_path):
    with pytest.raises(PinError, match="no such file"):
        pin_file(tmp_path / "absent.xlsx", "16.1", "2026-10-26")


@pytest.mark.parametrize(
    "field, value",
    [("retrieved", "26/10/2026"), ("sha256", "abc"), ("version", " "), ("sha256", "A" * 64)],
)
def test_bad_pin_fields(field, value):
    good = dict(
        guideline="EUCAST", version="16.1", retrieved="2026-10-26",
        sha256="0" * 64, size_bytes=1, filename="x.xlsx",
    )
    good[field] = value
    with pytest.raises(PinError):
        Pin(**good)


def test_malformed_pin_file(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")
    with pytest.raises(PinError):
        load_pin(bad)
