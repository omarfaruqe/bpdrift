"""Genotype-based resistance predictions and their error rates.

The predictor is deliberately simple and fixed, because the question is how
its measured error changes when only the phenotype labels change (breakpoint
version), not how good it is:

    an isolate is predicted resistant to an antibiotic when it carries at
    least one AMR element whose AMRFinderPlus subclass maps to that
    antibiotic in a reviewable mapping file.

Inputs
------
* AMRFinderPlus element calls per isolate, as NCBI Pathogen Detection writes
  them in its AMR_genotypes column: "blaKPC-2=COMPLETE,gyrA_S83I=POINT".
* The AMRFinderPlus reference gene catalogue (ReferenceGeneCatalog.txt) and
  point-mutation catalogue (AMRProt-mutation.tsv), giving each element's
  subclass, for example CARBAPENEM or AMIKACIN/KANAMYCIN/TOBRAMYCIN.
* A mapping CSV from subclass token to antibiotic:

      subclass_token,antibiotic,note
      CARBAPENEM,meropenem,carbapenemase
      TRIMETHOPRIM,trimethoprim/sulfamethoxazole,dfrA

  A compound subclass ("AMIKACIN/KANAMYCIN/TOBRAMYCIN") is split on "/",
  and each token is looked up separately.

Error rates (CLSI M52 / FDA definitions, on resistant vs not resistant):

    very major error (VME) = predicted not resistant, phenotype R  / phenotype R
    major error (ME)       = predicted resistant, phenotype S or I / phenotype S or I
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from bpdrift.diff import norm
from bpdrift.stats import wilson

NO_CALLS = {"", "NULL", "null", "NA", "-"}

# AMRFinderPlus call methods that do not count as evidence of a working element:
# a truncated gene (PARTIAL), a premature stop (MISTRANSLATION, INTERNAL_STOP) or a
# distant homologue (HMM). PARTIAL_END_OF_CONTIG (a gene split by assembly) counts.
EXCLUDED_METHODS = frozenset({"PARTIAL", "MISTRANSLATION", "INTERNAL_STOP", "HMM"})


class GenotypeError(ValueError):
    """Raised for malformed genotype inputs."""


# ----------------------------------------------------------------------
# element calls


def parse_amr_genotypes(text: str) -> list[tuple[str, str]]:
    """'blaKPC-2=COMPLETE,gyrA_S83I=POINT' -> [('blaKPC-2', 'COMPLETE'), ('gyrA_S83I', 'POINT')]."""
    text = (text or "").strip()
    if text in NO_CALLS:
        return []
    out = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        symbol, _, method = item.partition("=")
        symbol = symbol.strip()
        if not symbol:
            raise GenotypeError(f"empty element symbol in {text!r}")
        out.append((symbol, method.strip()))
    return out


# ----------------------------------------------------------------------
# catalogue: element symbol -> subclass


@dataclass
class Catalogue:
    subclass: dict[str, str] = field(default_factory=dict)   # symbol -> subclass
    family: dict[str, str] = field(default_factory=dict)     # symbol -> gene family
    disrupt: dict[str, str] = field(default_factory=dict)    # gene -> subclass when disrupted
    version: str = ""

    def lookup(self, symbol: str) -> str | None:
        """Subclass of an element; a frameshift or premature stop (``mgrB_K43Ter``,
        ``mgrB_G37YfsTer1``) in a gene with a POINT_DISRUPT entry takes that entry."""
        sub = self.subclass.get(symbol)
        if sub is None and "_" in symbol:
            gene, change = symbol.split("_", 1)
            if "Ter" in change:
                sub = self.disrupt.get(gene)
        return sub

    def family_of(self, symbol: str) -> str:
        """Gene family for a gene, the gene name for a mutation (gyrA_S83I -> gyrA)."""
        return self.family.get(symbol) or symbol.split("_", 1)[0]


def _read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames:
            reader.fieldnames = [f.lstrip("#").strip() for f in reader.fieldnames]
        return list(reader)


def load_catalogue(
    gene_catalog: str | Path,
    mutation_catalog: str | Path | None = None,
    organism_prefix: str = "Klebsiella",
    version: str = "",
) -> Catalogue:
    """Read AMRFinderPlus ReferenceGeneCatalog.txt (and AMRProt-mutation.tsv).

    Genes are keyed by allele (e.g. blaKPC-2) and by gene family (e.g. blaKPC),
    since NCBI reports whichever level AMRFinderPlus could resolve. Point
    mutations are keyed by mutation symbol (e.g. gyrA_S83I) for organisms
    whose taxgroup starts with organism_prefix (e.g. Klebsiella_pneumoniae).
    """
    cat = Catalogue(version=version)
    rows = _read_tsv(gene_catalog)
    if not rows or "subclass" not in rows[0]:
        raise GenotypeError(f"{gene_catalog}: no 'subclass' column")
    for r in rows:
        sub = (r.get("subclass") or "").strip()
        if not sub:
            continue
        fam = (r.get("gene_family") or "").strip()
        if (r.get("subtype") or "").strip() == "POINT_DISRUPT":
            taxa = (r.get("whitelisted_taxa") or "").strip()
            if fam and (not taxa or taxa.startswith(organism_prefix)):
                cat.disrupt[fam] = sub
            continue
        for key in ("allele", "gene_family"):
            sym = (r.get(key) or "").strip()
            if sym:
                cat.subclass.setdefault(sym, sub)
                cat.family.setdefault(sym, fam or sym)
    if mutation_catalog is not None:
        mrows = _read_tsv(mutation_catalog)
        if not mrows or "subclass" not in mrows[0]:
            raise GenotypeError(f"{mutation_catalog}: no 'subclass' column")
        loaded = 0
        for r in mrows:
            org = (r.get("taxgroup") or r.get("Organism") or "").strip()
            sym = (r.get("reported_mutation_symbol") or r.get("mutation_symbol")
                   or r.get("standard_mutation_symbol") or "").strip()
            sub = (r.get("subclass") or "").strip()
            if sym and sub and org.startswith(organism_prefix):
                cat.subclass[sym] = sub
                loaded += 1
        if loaded == 0:
            raise GenotypeError(
                f"{mutation_catalog}: no point mutations for taxgroup {organism_prefix!r}"
            )
    return cat


# ----------------------------------------------------------------------
# subclass token -> antibiotics


@dataclass
class SubclassMap:
    antibiotics: dict[str, set[str]] = field(default_factory=dict)  # TOKEN -> {antibiotic}

    def predicts(self, subclass: str) -> set[str]:
        out: set[str] = set()
        for token in subclass.split("/"):
            out |= self.antibiotics.get(token.strip().upper(), set())
        return out

    @property
    def targets(self) -> set[str]:
        return set().union(*self.antibiotics.values()) if self.antibiotics else set()


def load_subclass_map(path: str | Path) -> SubclassMap:
    with open(path, newline="", encoding="utf-8") as fh:
        lines = [t for t in fh if t.strip() and not t.lstrip().startswith("#")]
    reader = csv.DictReader(lines)
    need = {"subclass_token", "antibiotic"}
    if not reader.fieldnames or not need <= set(reader.fieldnames):
        raise GenotypeError(f"{path}: header must include subclass_token,antibiotic")
    m = SubclassMap()
    for r in reader:
        token = r["subclass_token"].strip().upper()
        ab = norm(r["antibiotic"])
        if not token or not ab:
            raise GenotypeError(f"{path}: empty token or antibiotic in {r}")
        m.antibiotics.setdefault(token, set()).add(ab)
    return m


# ----------------------------------------------------------------------
# predictions


@dataclass
class Exclusions:
    """Elements the predictor ignores: whole gene families, or single gene symbols."""
    families: dict[str, str] = field(default_factory=dict)   # family -> reason
    symbols: dict[str, str] = field(default_factory=dict)    # exact symbol -> reason


def load_exclusions(path: str | Path) -> Exclusions:
    """CSV family,reason[,match]: match is 'family' (default) or 'symbol'.

    A family entry ignores every member (oqxA ignores oqxA, oqxA10, ...); a
    symbol entry ignores that exact call only, so the chromosomal fosA of
    K. pneumoniae can be ignored while acquired fosA3 still counts.
    """
    with open(path, newline="", encoding="utf-8") as fh:
        lines = [t for t in fh if t.strip() and not t.lstrip().startswith("#")]
    reader = csv.DictReader(lines)
    if not reader.fieldnames or "family" not in reader.fieldnames:
        raise GenotypeError(f"{path}: header must include family")
    ex = Exclusions()
    for r in reader:
        name = r["family"].strip()
        if not name:
            continue
        match = (r.get("match") or "family").strip().lower()
        if match not in ("family", "symbol"):
            raise GenotypeError(f"{path}: match must be 'family' or 'symbol', got {match!r} for {name}")
        (ex.families if match == "family" else ex.symbols)[name] = (r.get("reason") or "").strip()
    return ex


def load_excluded_families(path: str | Path) -> dict[str, str]:
    """Family-level entries of an exclusions CSV (see load_exclusions)."""
    return load_exclusions(path).families


@dataclass
class Prediction:
    resistant: set[str]                       # normalised antibiotic names
    elements: list[tuple[str, str]]           # (symbol, method) as called
    unknown: list[str]                        # symbols missing from the catalogue
    ignored: list[str]                        # excluded family or non-counting method


def predict(
    calls: Iterable[tuple[str, str]],
    catalogue: Catalogue,
    smap: SubclassMap,
    excluded_families: Iterable[str] = (),
    excluded_methods: Iterable[str] = EXCLUDED_METHODS,
    excluded_symbols: Iterable[str] = (),
) -> Prediction:
    calls = list(calls)
    fams, methods, syms = set(excluded_families), set(excluded_methods), set(excluded_symbols)
    resistant: set[str] = set()
    unknown: list[str] = []
    ignored: list[str] = []
    for symbol, method in calls:
        if method in methods:
            ignored.append(f"{symbol}={method}")
            continue
        sub = catalogue.lookup(symbol)
        if sub is None:
            unknown.append(symbol)
            continue
        if symbol in syms or catalogue.family_of(symbol) in fams:
            ignored.append(symbol)
            continue
        resistant |= smap.predicts(sub)
    return Prediction(resistant, calls, unknown, ignored)


# ----------------------------------------------------------------------
# error rates


@dataclass(frozen=True)
class ErrorRates:
    n: int
    phenotype_r: int
    phenotype_not_r: int
    vme: int
    me: int

    @property
    def vme_rate(self) -> tuple[float, float, float]:
        return wilson(self.vme, self.phenotype_r)

    @property
    def me_rate(self) -> tuple[float, float, float]:
        return wilson(self.me, self.phenotype_not_r)

    @property
    def agreement(self) -> tuple[float, float, float]:
        return wilson(self.n - self.vme - self.me, self.n)


def error_rates(pairs: Iterable[tuple[bool, str]]) -> ErrorRates:
    """pairs: (predicted resistant, phenotype category 'S' | 'I' | 'R').

    Other categories (UNDECIDABLE, NO_BREAKPOINT) must be filtered out first;
    passing one raises, so nothing is dropped silently.
    """
    c = Counter()
    for predicted_r, category in pairs:
        if category not in ("S", "I", "R"):
            raise GenotypeError(f"phenotype must be S, I or R, got {category!r}")
        c["n"] += 1
        if category == "R":
            c["r"] += 1
            if not predicted_r:
                c["vme"] += 1
        else:
            c["not_r"] += 1
            if predicted_r:
                c["me"] += 1
    return ErrorRates(c["n"], c["r"], c["not_r"], c["vme"], c["me"])


def collapse_labels(
    rows: Iterable[tuple[str, str, str, str]]
) -> tuple[dict[tuple[str, str, str], str], int]:
    """One label per (genome, antibiotic, version) from (genome, antibiotic, version, category).

    A genome with several records for the same drug and version keeps the
    label only when they all agree. Returns (labels, number of keys dropped
    because their records disagreed).
    """
    seen: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for genome, ab, version, category in rows:
        seen[(genome, norm(ab), version)].add(category)
    labels = {k: next(iter(v)) for k, v in seen.items() if len(v) == 1}
    return labels, len(seen) - len(labels)
