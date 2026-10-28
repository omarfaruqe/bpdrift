"""Supplementary Table S1: every input file with its version and checksum.

    python analysis/inputs_manifest.py data -o results/table_s1_inputs.csv     # write
    python analysis/inputs_manifest.py data --check results/table_s1_inputs.csv  # verify

Lists the EUCAST tables (version read from each workbook), the BV-BRC exports,
the NCBI files and the NCBI BioSample antibiogram batches under data/, with
path, size and SHA-256. File contents are not copied.

The table records only what identifies the bytes (path, size, checksum) and a
version read from the content, never filesystem dates, so identical inputs give
an identical table wherever and whenever they were downloaded. Retrieval dates
are stated in the paper's Methods.

--check verifies the inputs present under data/ against a committed table:
every listed file must exist with the same size and checksum, its source and
version columns must match what is derived from the file now, and no unlisted
file may sit in a listed folder (it could change what a glob in the pipeline
picks up). Exit status 0 when everything matches, 1 otherwise, with each
problem listed.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from bpdrift.eucast import read_eucast_xlsx
from bpdrift.pin import sha256_of

SOURCES = {
    "breakpoints": "EUCAST clinical breakpoint table",
    "bvbrc": "BV-BRC export",
    "ncbi": "NCBI Pathogen Detection / AMRFinderPlus",
    "ncbi/ast": "NCBI BioSample antibiograms (E-utilities)",
}
COLUMNS = ["source", "path", "version", "bytes", "sha256"]


def describe(path: Path, folder: str) -> tuple[str, str]:
    if folder == "breakpoints" and path.suffix == ".xlsx":
        return SOURCES[folder], f"v{read_eucast_xlsx(path).version}"
    if path.name.endswith(".metadata.tsv"):
        return SOURCES[folder], path.name.split(".metadata")[0]
    if path.name == "amrfinder_db_version.txt":
        return SOURCES[folder], f"database {path.read_text().strip()}"
    return SOURCES.get(folder, folder), ""


def inputs(root: Path) -> list[Path]:
    """Input files directly inside each source folder (subfolders are their own sources)."""
    out = []
    for folder in SOURCES:
        for f in sorted((root / folder).glob("*")):
            if f.is_file() and not f.name.startswith("."):
                out.append(f)
    return out


def rel(path: Path, root: Path) -> str:
    return (Path(root.name) / path.relative_to(root)).as_posix()


def write(root: Path, output: str) -> int:
    rows = []
    for f in inputs(root):
        folder = f.parent.relative_to(root).as_posix()
        source, version = describe(f, folder)
        rows.append([source, rel(f, root), version, f.stat().st_size, sha256_of(f)])
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(COLUMNS)
        w.writerows(rows)
    print(f"wrote {output}: {len(rows)} files")
    return 0


def check(root: Path, manifest: str) -> int:
    with open(manifest, newline="", encoding="utf-8") as fh:
        listed = {r["path"]: r for r in csv.DictReader(fh)}
    problems = []
    base = root.parent
    for path, r in listed.items():
        f = base / path
        if not f.is_file():
            problems.append(f"missing: {path}")
        elif f.stat().st_size != int(r["bytes"]):
            problems.append(f"size differs: {path} ({f.stat().st_size} bytes, expected {r['bytes']})")
        elif sha256_of(f) != r["sha256"]:
            problems.append(f"checksum differs: {path}")
        else:
            source, version = describe(f, f.parent.relative_to(root).as_posix())
            if (source, version) != (r["source"], r["version"]):
                problems.append(f"metadata differs: {path} (manifest: {r['source']!r}, {r['version']!r}; "
                                f"file: {source!r}, {version!r})")
    for f in inputs(root):
        if rel(f, root) not in listed:
            problems.append(f"not in the manifest: {rel(f, root)}")
    if problems:
        print(f"inputs do NOT match {manifest}:")
        for p in problems:
            print(f"    {p}")
        return 1
    print(f"inputs match {manifest}: {len(listed)} files, sizes, SHA-256, sources and versions identical")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("data", help="the data/ folder")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("-o", "--output", help="write the manifest here")
    g.add_argument("--check", metavar="MANIFEST", help="verify data/ against this manifest")
    a = p.parse_args(argv)
    root = Path(a.data)
    return write(root, a.output) if a.output else check(root, a.check)


if __name__ == "__main__":
    raise SystemExit(main())
