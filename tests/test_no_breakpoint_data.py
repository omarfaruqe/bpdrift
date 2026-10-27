"""Guard: nothing in the repository or the built package may carry EUCAST breakpoint values.

EUCAST tables are copyrighted and not redistributed. Outputs of `bpdrift diff`
and `bpdrift reinterpret` repeat breakpoint values row by row, so they must
never be committed either. This test fails if a tracked file is an Excel table
or a CSV with the breakpoint columns those commands write.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BREAKPOINT_COLUMNS = {"old_S", "old_R", "new_S", "new_R", "eucast_S", "eucast_R"}


def tracked_files():
    if shutil.which("git") is None or not (ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True)
    return [ROOT / line for line in out.stdout.splitlines() if line]


def test_no_excel_tables_tracked():
    bad = [p.name for p in tracked_files() if p.suffix.lower() in {".xlsx", ".xls", ".xlsm"}]
    assert not bad, f"breakpoint tables must not be committed: {bad}"


def test_no_tracked_csv_has_breakpoint_columns():
    bad = []
    for p in tracked_files():
        if p.suffix.lower() != ".csv" or not p.exists():
            continue
        lines = [t for t in p.read_text(encoding="utf-8", errors="replace").splitlines()
                 if t.strip() and not t.lstrip().startswith("#")]
        if lines and BREAKPOINT_COLUMNS & {c.strip() for c in lines[0].split(",")}:
            bad.append(str(p.relative_to(ROOT)))
    assert not bad, f"these CSVs repeat breakpoint values and must not be committed: {bad}"
