"""Supplementary Table S1: every input file with its version and checksum.

    python analysis/inputs_manifest.py data -o results/table_s1_inputs.csv

Lists the EUCAST tables (version read from each workbook), the BV-BRC exports
and the NCBI files under data/, with size and SHA-256, so a reader can check
that a re-run uses byte-identical inputs. File contents are not copied.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
from pathlib import Path

from bpdrift.eucast import read_eucast_xlsx
from bpdrift.pin import sha256_of

SOURCES = {
    "breakpoints": "EUCAST clinical breakpoint table",
    "bvbrc": "BV-BRC export",
    "ncbi": "NCBI Pathogen Detection / AMRFinderPlus",
}


def describe(path: Path, folder: str) -> tuple[str, str]:
    if folder == "breakpoints" and path.suffix == ".xlsx":
        return SOURCES[folder], f"v{read_eucast_xlsx(path).version}"
    if path.name.endswith(".metadata.tsv"):
        return SOURCES[folder], path.name.split(".metadata")[0]
    if path.name == "amrfinder_db_version.txt":
        return SOURCES[folder], f"database {path.read_text().strip()}"
    return SOURCES.get(folder, folder), ""


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("data", help="the data/ folder")
    p.add_argument("-o", "--output", required=True)
    a = p.parse_args(argv)
    root = Path(a.data)
    rows = []
    for folder in SOURCES:
        for f in sorted((root / folder).glob("*")):
            if not f.is_file() or f.name.startswith("."):
                continue
            source, version = describe(f, folder)
            modified = dt.datetime.fromtimestamp(f.stat().st_mtime).date().isoformat()
            rows.append([source, f.name, version, f.stat().st_size, modified, sha256_of(f)])
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    with open(a.output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["source", "file", "version", "bytes", "file_date", "sha256"])
        w.writerows(rows)
    print(f"wrote {a.output}: {len(rows)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
