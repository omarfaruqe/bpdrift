---
title: 'bpdrift: version-pinned interpretation of antimicrobial susceptibility data under changing EUCAST breakpoints'
tags:
  - Python
  - antimicrobial resistance
  - minimum inhibitory concentration
  - clinical breakpoints
  - EUCAST
  - genotype-phenotype prediction
authors:
  - name: Md Omar Faruqe
    affiliation: 1
affiliations:
  - name: Department of Computer Science and Engineering, University of Rajshahi, Bangladesh
    index: 1
date: 27 October 2026
bibliography: paper.bib
---

# Summary

Antimicrobial susceptibility is measured as a minimum inhibitory concentration
(MIC) and reported as susceptible (S), susceptible at increased exposure (I) or
resistant (R) by comparing the MIC with clinical breakpoints. The European
Committee on Antimicrobial Susceptibility Testing (EUCAST) revises its
breakpoint tables every year, so the same MIC can carry different labels
depending on the table used. `bpdrift` is a Python package and command-line
tool that interprets MICs under an explicitly pinned version of the EUCAST
tables, compares versions, and measures how much the labels of a dataset, and
the measured accuracy of a genotype-based resistance predictor, depend on the
version chosen.

# Statement of need

Public genotype-phenotype resources such as BV-BRC [@olson2023] pool
susceptibility records made under different standards and years, and machine
learning benchmarks for resistance prediction treat their labels as fixed.
Breakpoint revisions can change reported resistance with no change in the
organisms [@sallam2025], but tools rarely make the version an explicit,
checkable input. `bpdrift` was written to make it one, for researchers who
build or evaluate resistance-prediction benchmarks and for surveillance
analyses that compare data across years.

Four design decisions follow from that aim.

1. **No breakpoint data is bundled.** The EUCAST tables carry no licence
   permitting redistribution. Users download the official files; `bpdrift`
   parses them and pins each by version, retrieval date and SHA-256, so every
   label can be traced to the exact file that produced it.
2. **Censored values are never guessed.** An MIC such as "<=2" or ">8" is
   treated as the range of values it allows. It receives a category only if
   the whole range lies within it; otherwise the result is `UNDECIDABLE`,
   with the categories it could take.
3. **Matching is explicit.** EUCAST renames, splits and reorganises rows
   between versions (for example by indication, route or organism). `bpdrift`
   matches rows only on exact names and records every other correspondence
   in a reviewable mapping file; a mapping line naming a row that does not
   exist is an error.
4. **Nothing is dropped silently.** Unrecognised table cells, empty table
   blocks, malformed mappings and unusable source records stop the run or are
   counted and reported.

# Functionality

The package reads EUCAST clinical breakpoint tables in Excel format from
v8.1 (2018) to v16.1 (2026), separating superscript footnote markers from
values (a cell reading "82" is the value 8 with note 2) and keeping the
organism subgroup of each row. It provides MIC parsing and interpretation, a
row-by-row comparison of two tables that lists the MIC dilutions whose
category changes, an importer for BV-BRC laboratory AST exports, and
re-interpretation of every imported MIC under several tables. A genotype
module scores a fixed predictor built from AMRFinderPlus calls
[@feldgarden2021] against labels from each version, reporting very major and
major error rates. Proportions are reported with Wilson score intervals
[@wilson1927]. All functions are available from Python and from the
`bpdrift` command line (`pin`, `interpret`, `diff`, `reinterpret`).

Applied to 47,804 *Klebsiella pneumoniae* MICs from BV-BRC, re-interpretation
under seven EUCAST versions changed 2.9% of labels between 2018 and 2026,
concentrated in a few antibiotics (for example 13.7% for
trimethoprim-sulfamethoxazole), and moved the measured error rates of a fixed
genotype predictor (imipenem major errors fell from 15.9% to 7.8%). Reported
calls stated as EUCAST 2020 were reproduced exactly under the version then in
force. The analysis scripts and reviewed mapping files are in the
repository.

# Relation to existing software

The R package AMR [@berends2022] provides comprehensive tools for
antimicrobial resistance data analysis, including interpretation of MICs
under bundled EUCAST and CLSI breakpoints whose public use EUCAST and CLSI
endorse. `bpdrift` is narrower and complementary. It is pure Python, works
from the user's own copy of any EUCAST table version rather than a bundled
dataset, pins that copy by checksum, and focuses on differences between
versions and their effect on labelled datasets. It does not implement EUCAST
expert rules, intrinsic resistance or CLSI breakpoints.

# AI usage disclosure

The software, tests, analysis scripts and a first draft of this paper were
developed with the assistance of a large language model (Claude, Anthropic),
working under the author's direction. The author specified the design,
reviewed the mapping decisions, checked the results against their sources
and takes responsibility for the content.

# References
