# Changelog

## 0.1.0 (2026-10-27)

First release.

- Parse MIC values with their censoring sign; refuse anything ambiguous.
- Interpret an MIC against an EUCAST breakpoint pair, returning UNDECIDABLE
  when a censored value could fall in more than one category.
- Pin breakpoint files by version, retrieval date and SHA-256.
- Read EUCAST clinical breakpoint tables (.xlsx) from v8.1 (2018) to v16.1
  (2026): superscript footnotes separated from values, organism subgroups
  kept, spacer rows tolerated, and any unrecognised cell or empty table
  block reported with its cell reference.
- Compare two tables row by row, with an explicit mapping file for renamed
  or split rows and the MIC dilutions whose category changes.
- Import BV-BRC laboratory AST exports with every excluded record counted.
- Re-interpret imported MICs under several tables through reviewed
  antibiotic links.
- Flip-rate statistics with Wilson intervals; genotype predictor error rates.
- Command-line interface: pin, interpret, diff, reinterpret.
