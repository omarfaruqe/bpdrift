# Mapping files

Each file pairs EUCAST breakpoint rows that were renamed or split between
two consecutive table versions, for use with `bpdrift diff --map`.
Rows that keep the same name need no mapping.

Status: **draft, pending review.** Scope: Enterobacterales, agents
relevant to *Klebsiella pneumoniae*.

Rules applied:

1. Renames (wording, punctuation, taxonomy) map one to one.
2. Indication or route splits (meningitis / other, iv / oral) map to every
   resulting row; the analysis chooses the indication.
3. Organism splits map only to the row that covers *Klebsiella* spp.
   A drug with no such row is left unmapped and is reported as removed.

These files hold agent names only, no breakpoint values.
