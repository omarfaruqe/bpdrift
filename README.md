# bpdrift

Interpret minimum inhibitory concentration (MIC) values as S / I / R under a
specific, pinned version of the EUCAST clinical breakpoint tables, and report
where two versions disagree.

## Design principles

- **No breakpoint data is bundled.** EUCAST tables are copyright EUCAST. You
  download the official file yourself; bpdrift parses it and pins it by
  version string, retrieval date and SHA-256.
- **Censored MICs are handled explicitly.** A value such as `<=2` or `>8` that
  could fall on either side of a breakpoint returns `UNDECIDABLE`, never a
  guess.
- **Fail loudly on ambiguity.** Unparseable input raises an error.

## Usage

Download the EUCAST clinical breakpoint tables (.xlsx, "file for screen")
from eucast.org into `data/breakpoints/`, a folder git ignores.

```bash
# record exactly which file you used
bpdrift pin data/breakpoints/v_16.1_Breakpoint_Tables.xlsx \
    --version 16.1 --retrieved 2026-10-26 -o data/breakpoints/v16.1.pin.json

# interpret one MIC
bpdrift interpret "<=4" --table data/breakpoints/v_16.1_Breakpoint_Tables.xlsx \
    --group Enterobacterales --agent "Meropenem (indications other than meningitis)"

# compare two versions; one CSV row per breakpoint, with the dilutions that change category
bpdrift diff data/breakpoints/v_10.0_Breakpoint_Tables.xlsx \
    data/breakpoints/v_16.1_Breakpoint_Tables.xlsx \
    --group Enterobacterales --changed-only -o diff.csv
```

Rows are matched across versions by organism group and agent name only.
Rows EUCAST renamed or split are reported as removed and added, with
likely partners listed as candidates. To pair them, list them in a CSV
passed with `--map`:

```text
old_group,old_agent,new_group,new_agent
Enterobacterales,Ciprofloxacin,Enterobacterales,Ciprofloxacin (meningitis)
Enterobacterales,Ciprofloxacin,Enterobacterales,Ciprofloxacin (indications other than meningitis)
P.multocida,,Pasteurella,
```

A line with no agents renames a whole organism group. A mapping line
naming a row that does not exist is an error.

## Status

Early development (v0.0.x). Not for clinical use.

## Licence

Code: MIT. Breakpoint tables remain the property of their publishers.
