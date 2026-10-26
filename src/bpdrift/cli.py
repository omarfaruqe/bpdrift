"""Command-line interface.

    bpdrift pin FILE --version 16.1 --retrieved 2026-10-26 [-o FILE.pin.json]
    bpdrift interpret MIC --table FILE --group G --agent A [--pin P]
    bpdrift diff OLD NEW [--map map.csv] [--group G ...] [--changed-only]
                         [--old-pin P] [--new-pin P] [-o diff.csv]
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from bpdrift import __version__
from bpdrift.diff import ChangeKind, DiffError, diff_tables, load_mapping, norm
from bpdrift.eucast import EucastFormatError, read_eucast_xlsx
from bpdrift.interpret import interpret
from bpdrift.mic import MICParseError
from bpdrift.pin import PinError, load_pin, pin_file, save_pin

DIFF_COLUMNS = [
    "kind", "match",
    "old_version", "old_group", "old_agent", "old_cell", "old_S", "old_R",
    "new_version", "new_group", "new_agent", "new_cell", "new_S", "new_R",
    "flips", "candidates",
]


class CLIError(Exception):
    pass


def _pin_arg(path: str | None):
    return load_pin(path) if path else None


def _cmd_pin(args) -> int:
    pin = pin_file(args.file, args.version, args.retrieved, args.guideline)
    if args.output:
        save_pin(pin, args.output)
        print(f"wrote {args.output}: {pin.label}", file=sys.stderr)
    else:
        print(pin.label)
        print(f"sha256 {pin.sha256}")
    return 0


def _cmd_interpret(args) -> int:
    table = read_eucast_xlsx(args.table, pin=_pin_arg(args.pin))
    matches = [
        r for r in table.rows
        if norm(args.group) in (norm(r.group), norm(r.sheet)) and norm(r.agent) == norm(args.agent)
    ]
    if not matches:
        raise CLIError(f"no agent {args.agent!r} in group {args.group!r} (v{table.version})")
    if len(matches) > 1:
        where = ", ".join(f"{m.cell} ({m.group})" for m in matches)
        raise CLIError(f"{args.agent!r} is ambiguous in {args.group!r}: {where}; name the subgroup")
    row = matches[0]
    result = interpret(args.mic, row.breakpoint())
    print(f"{result.category.value}\t{result.reason}\tEUCAST v{table.version} {row.cell}")
    return 0


def _fmt_value(v) -> str:
    return v.raw if v.raw else v.kind


def _cmd_diff(args) -> int:
    old = read_eucast_xlsx(args.old, pin=_pin_arg(args.old_pin))
    new = read_eucast_xlsx(args.new, pin=_pin_arg(args.new_pin))
    mapping = load_mapping(args.map) if args.map else None
    diff = diff_tables(old, new, mapping, args.group or None)

    out = open(args.output, "w", newline="", encoding="utf-8") if args.output else sys.stdout
    try:
        writer = csv.writer(out)
        writer.writerow(DIFF_COLUMNS)
        for c in diff.changes:
            if args.changed_only and c.kind is ChangeKind.UNCHANGED:
                continue
            o, n = c.old, c.new
            writer.writerow([
                c.kind.value, c.match,
                old.version, o.group if o else "", o.agent if o else "", o.cell if o else "",
                _fmt_value(o.s) if o else "", _fmt_value(o.r) if o else "",
                new.version, n.group if n else "", n.agent if n else "", n.cell if n else "",
                _fmt_value(n.s) if n else "", _fmt_value(n.r) if n else "",
                ";".join(f"{f.mic:g}:{f.old.value}>{f.new.value}" for f in c.flips),
                " | ".join(c.candidates),
            ])
    finally:
        if out is not sys.stdout:
            out.close()

    counts = ", ".join(f"{k} {v}" for k, v in diff.summary().items())
    print(f"EUCAST v{old.version} -> v{new.version}: {counts}", file=sys.stderr)
    if args.output:
        print(f"wrote {args.output}", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="bpdrift", description="Version-pinned EUCAST MIC breakpoints.")
    p.add_argument("--version", action="version", version=f"bpdrift {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    sp = sub.add_parser("pin", help="record a breakpoint file's version, date and SHA-256")
    sp.add_argument("file")
    sp.add_argument("--version", required=True, dest="version", help="table version, e.g. 16.1")
    sp.add_argument("--retrieved", required=True, help="download date, YYYY-MM-DD")
    sp.add_argument("--guideline", default="EUCAST")
    sp.add_argument("-o", "--output", help="write the pin as JSON")
    sp.set_defaults(func=_cmd_pin)

    sp = sub.add_parser("interpret", help="interpret one MIC")
    sp.add_argument("mic", help="MIC in mg/L, e.g. '<=0.25' or 16")
    sp.add_argument("--table", required=True)
    sp.add_argument("--group", required=True, help="sheet or organism subgroup")
    sp.add_argument("--agent", required=True, help="agent exactly as written in the table")
    sp.add_argument("--pin")
    sp.set_defaults(func=_cmd_interpret)

    sp = sub.add_parser("diff", help="compare two breakpoint tables")
    sp.add_argument("old")
    sp.add_argument("new")
    sp.add_argument("--map", help="CSV mapping renamed or split rows")
    sp.add_argument("--group", action="append", help="limit to a sheet or subgroup (repeatable)")
    sp.add_argument("--changed-only", action="store_true", help="omit unchanged rows")
    sp.add_argument("--old-pin")
    sp.add_argument("--new-pin")
    sp.add_argument("-o", "--output", help="write CSV here instead of stdout")
    sp.set_defaults(func=_cmd_diff)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (CLIError, DiffError, EucastFormatError, MICParseError, PinError, OSError) as exc:
        print(f"bpdrift: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
