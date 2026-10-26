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

## Status

Early development (v0.0.x). Not for clinical use.

## Licence

Code: MIT. Breakpoint tables remain the property of their publishers.
