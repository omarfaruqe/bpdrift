# Mapping files

Each file pairs EUCAST breakpoint rows that were renamed or split between
two consecutive table versions, for use with `bpdrift diff --map`.
Rows that keep the same name need no mapping.

Status: **reviewed** by Md Omar Faruqe, 2026-10-26. Scope:
Enterobacterales, agents relevant to *Klebsiella pneumoniae*.

## Rules

1. Renames (wording, punctuation, taxonomy) map one to one.
2. Indication or route splits (meningitis / other, iv / oral) map to every
   resulting row. The mapping records the split; the analysis chooses
   which row applies.
3. Organism splits map only to the row that covers *Klebsiella* spp.
   A drug with no such row is left unmapped and is reported as removed.

## Decisions for the *K. pneumoniae* analysis

| Split | Row used in the analysis | Reason |
|---|---|---|
| Meningitis / indications other than meningitis (meropenem, ciprofloxacin, cefotaxime, ceftriaxone) | indications other than meningitis | public AST isolates are not meningitis cases |
| iv / oral (amoxicillin-clavulanic acid, ampicillin) | iv | systemic-infection breakpoint, comparable with the pre-split row |
| Organism-specific rows with no *Klebsiella* row (oral fosfomycin from v13.1; *Proteus* and *Serratia* rows) | none | drug drops out of the *Klebsiella* analysis for that version |

These files hold agent names only, no breakpoint values.

## BV-BRC antibiotic links

`bvbrc_to_eucast_enterobacterales.csv` links each BV-BRC antibiotic name to
the EUCAST Enterobacterales row used for *K. pneumoniae* in each version.
Reviewed by Md Omar Faruqe, 2026-10-26. Beyond the decisions above:

| Decision | Choice |
|---|---|
| Aminoglycosides: systemic / urinary-tract rows | systemic |
| Agents whose only *Klebsiella* breakpoint is urinary-tract specific (cefazolin, cefpodoxime, trimethoprim) | kept, tagged `restriction=urinary`; excluded from the main analysis |
| Agents with no applicable numeric breakpoint in any version (e.g. tigecycline, nitrofurantoin: *E. coli* only; tetracycline, minocycline: "-"; cefoxitin: screen only; polymyxin B: no row) | not linked |
| Combination MICs reported as a ratio pair for amoxicillin-clavulanic acid and ampicillin-sulbactam (CLSI 2:1 method; EUCAST uses a fixed inhibitor concentration) | excluded from the main analysis |
| MICs whose sign field disagrees with the reading | sign field trusted; sensitivity analysis treats them as exact |
