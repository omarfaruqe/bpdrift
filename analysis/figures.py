"""Paper figures from the results tables.

    python analysis/figures.py results/ -o results/figures

Figure 1  heatmap: share of records changing category at each EUCAST
          revision and overall, per antibiotic (from robustness.csv;
          drawn when it exists).
Figure S1 label flip rate between the first and last EUCAST version, per
          antibiotic, with Wilson 95% intervals; main analysis and the
          sign-exact sensitivity analysis.
Figure 2  agreement of BV-BRC's reported calls with the re-interpreted
          labels (resistant vs not resistant) under each EUCAST version,
          by the testing standard the record states.
Figure 3  very major and major error rates of the genotype predictor under
          each EUCAST version (drawn when results/genotype_errors.csv exists).

Colour: categorical slots in fixed order from a palette validated for
colour-vision deficiency; identity is also carried by marker shape and a
direct label, never by colour alone.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]  # validated, fixed order
INK, INK_2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
MIN_RECORDS = 200  # antibiotics with fewer decidable pairs are left out of Fig. 1

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.edgecolor": INK_2,
    "axes.labelcolor": INK, "xtick.color": INK_2, "ytick.color": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
    "svg.hashsalt": "bpdrift",  # deterministic element ids, so reruns give identical SVGs
})
METADATA = {"svg": {"Date": None}, "png": {}}  # no timestamps in the files


def _read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _pretty(ab: str) -> str:
    return ab.replace("/", "-")


def _vkey(v):
    return tuple(int(x) for x in v.split("."))


def figure_1(flips, out: Path, a=None, b=None):
    versions = sorted({r["version_a"] for r in flips} | {r["version_b"] for r in flips}, key=_vkey)
    a, b = a or versions[0], b or versions[-1]
    pick = {(r["scenario"], r["antibiotic"]): r for r in flips
            if r["version_a"] == a and r["version_b"] == b and r["antibiotic"] != "ALL"}
    drugs = [ab for (sc, ab), r in pick.items() if sc == "main" and int(r["both_decidable"]) >= MIN_RECORDS]
    drugs.sort(key=lambda ab: (float(pick[("main", ab)]["flip_rate"]), ab))  # largest at the top

    fig, ax = plt.subplots(figsize=(7.0, 0.22 * len(drugs) + 1.3))
    for i, ab in enumerate(drugs):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color="#f4f3f0", zorder=0, lw=0)
    for k, (scenario, label, marker, dy) in enumerate([
        ("main", "Main analysis", "o", 0.12),
        ("sign_exact", "Sign conflicts taken as exact", "D", -0.12),
    ]):
        for i, ab in enumerate(drugs):
            r = pick[(scenario, ab)]
            p, lo, hi = (100 * float(r[x]) for x in ("flip_rate", "flip_lo", "flip_hi"))
            ax.plot([lo, hi], [i + dy] * 2, color=SERIES[k], lw=2, solid_capstyle="round", zorder=2)
            ax.plot(p, i + dy, marker=marker, ms=6 if marker == "o" else 5, color=SERIES[k],
                    mec="white", mew=1, zorder=3, label=label if i == 0 else None)
            if scenario == "main" and int(r["flips"]) > 0:
                far = max(hi, 100 * float(pick[("sign_exact", ab)]["flip_hi"]))
                ax.annotate(f"{p:.1f}%  ({r['flips']}/{r['both_decidable']})", (far, i),
                            xytext=(8, 0), textcoords="offset points", va="center", fontsize=8, color=INK)
    ax.set_yticks(range(len(drugs)), [_pretty(d) for d in drugs])
    ax.set_ylim(-0.6, len(drugs) - 0.4)
    widest = max(100 * float(pick[(sc, ab)]["flip_hi"]) for sc in ("main", "sign_exact") for ab in drugs)
    ax.set_xlim(-0.6, widest * 1.45 + 2)
    ax.set_xlabel(f"Records whose S/I/R category differs between EUCAST v{a} and v{b} (%)")
    ax.xaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.legend(loc="lower right", frameon=False, fontsize=8)
    changed = sum(int(pick[("main", ab)]["flips"]) > 0 for ab in drugs)
    fig.text(0.01, 0.995, f"Breakpoint revisions from v{a} to v{b} changed labels for {changed} of "
             f"{len(drugs)} Klebsiella drugs", fontsize=10, fontweight="bold", va="top")
    fig.text(0.01, 0.005, "K. pneumoniae, BV-BRC laboratory MICs; records decidable under both versions; "
             "bars are Wilson 95% intervals.", fontsize=7, color=INK_2)
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    for ext in ("png", "svg"):
        fig.savefig(out / f"figureS1_flip_rates.{ext}", dpi=300, metadata=METADATA[ext])
    plt.close(fig)


def figure_2(agree, out: Path):
    versions = sorted({r["eucast_version"] for r in agree}, key=_vkey)
    strata = [  # (standard, year, label, marker)
        ("EUCAST", "2019", "EUCAST, stated year 2019", "v"),
        ("EUCAST", "2020", "EUCAST, stated year 2020", "o"),
        ("EUCAST", "all", "EUCAST, all years", "D"),
        ("none", "all", "No standard stated", "s"),
        ("CLSI", "all", "CLSI, all years", "^"),
    ]
    idx = {(r["standard"], r["year"], r["eucast_version"]): r for r in agree}
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    x = range(len(versions))
    # a stratum is drawn only when it has reported calls under every version
    strata = [s for s in strata if all((s[0], s[1], v) in idx for v in versions)]
    for k, (std, year, label, marker) in enumerate(strata):
        rows = [idx[(std, year, v)] for v in versions]
        p = [100 * float(r["r_line_rate"]) for r in rows]
        lo = [100 * float(r["r_line_lo"]) for r in rows]
        hi = [100 * float(r["r_line_hi"]) for r in rows]
        ax.plot(x, p, color=SERIES[k], lw=2, marker=marker, ms=7, mec="white", mew=1, label=label, zorder=3)
        ax.vlines(x, lo, hi, color=SERIES[k], lw=2, zorder=2)
        n = rows[-1]["compared"]
        ax.annotate(f"{label}  (n = {n})", (x[-1], p[-1]), xytext=(10, 0), textcoords="offset points",
                    va="center", fontsize=8, color=INK)
    valid = [idx[k] for k in (("EUCAST", "2019", "9.0"), ("EUCAST", "2020", "10.0")) if k in idx]
    note = "; ".join(f"stated {r['year']} vs v{r['eucast_version']} (then in force): "
                     f"{r['agreement']}/{r['compared']} reproduced" for r in valid)
    ax.text(-0.08, 100.45, note, fontsize=7.5, color=INK_2, va="bottom")
    ax.set_xticks(list(x), [f"v{v}" for v in versions])
    ax.set_xlim(-0.3, len(versions) - 1 + 1.9)
    ax.set_ylim(91.5, 101.2)
    ax.set_yticks(range(92, 101))
    ax.set_ylabel("Reported call agrees with\nre-interpreted label (%)")
    ax.set_xlabel("EUCAST version used to re-interpret the MIC")
    ax.yaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    fig.text(0.01, 0.995, "Reported labels match the EUCAST version they were made under, then drift",
             fontsize=10, fontweight="bold", va="top")
    fig.text(0.01, 0.01, "Agreement on resistant vs not resistant; K. pneumoniae, BV-BRC, main analysis.\n"
             f"n = records compared under v{versions[-1]}; bars are Wilson 95% intervals.", fontsize=7, color=INK_2)
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    for ext in ("png", "svg"):
        fig.savefig(out / f"figure2_reported_agreement.{ext}", dpi=300, metadata=METADATA[ext])
    plt.close(fig)


def figure_3(errors, out: Path, min_compared=200, min_shift=2.0):
    """Very major and major error rates of the genotype predictor across versions."""
    versions = sorted({r["eucast_version"] for r in errors}, key=_vkey)
    idx = {(r["antibiotic"], r["eucast_version"]): r for r in errors}
    shifts = []
    for ab in {r["antibiotic"] for r in errors} - {"ALL"}:
        rows = [idx.get((ab, v)) for v in (versions[0], versions[-1])]
        if not all(rows) or int(rows[1]["compared"]) < min_compared:
            continue
        d = max(abs(100 * (float(rows[1][k]) - float(rows[0][k]))) for k in ("vme_rate", "me_rate")
                if rows[0][k] != "nan" and rows[1][k] != "nan")
        if d >= min_shift:
            shifts.append((d, ab))
    drugs = [ab for _, ab in sorted(shifts, reverse=True)]
    if not drugs:
        return
    cols = min(4, len(drugs))
    nrows = -(-len(drugs) // cols)
    fig, axes = plt.subplots(nrows, cols, figsize=(7.0, 2.25 * nrows + 0.9), sharex=True, squeeze=False)
    x = range(len(versions))
    for i, ab in enumerate(drugs):
        ax = axes[i // cols][i % cols]
        for k, (key, label, marker) in enumerate([("vme", "Very major errors", "o"), ("me", "Major errors", "D")]):
            pts = [(j, idx.get((ab, v))) for j, v in enumerate(versions)]
            pts = [(j, r) for j, r in pts if r and r[f"{key}_rate"] != "nan"]
            xs = [j for j, _ in pts]
            p = [100 * float(r[f"{key}_rate"]) for _, r in pts]
            lo = [100 * float(r[f"{key}_lo"]) for _, r in pts]
            hi = [100 * float(r[f"{key}_hi"]) for _, r in pts]
            ax.vlines(xs, lo, hi, color=SERIES[k], lw=1.5, alpha=0.6, zorder=2)
            ax.plot(xs, p, color=SERIES[k], lw=2, marker=marker, ms=4.5, mec="white", mew=0.8,
                    label=label if i == 0 else None, zorder=3)
            ax.annotate(f"{p[-1]:.0f}%", (xs[-1], p[-1]), xytext=(4, 0), textcoords="offset points",
                        va="center", fontsize=7, color=INK)
        short = {"trimethoprim/sulfamethoxazole": "co-trimoxazole"}.get(ab, _pretty(ab))
        ax.set_title(short, fontsize=8.5, loc="left", color=INK)
        ax.set_ylim(bottom=0)
        ax.yaxis.grid(True, color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        ax.set_xticks(list(x), [f"v{v}" for v in versions], rotation=90, fontsize=6.5)
        ax.tick_params(axis="y", labelsize=7)
        ax.set_xlim(-0.4, len(versions) - 0.2)
    for j in range(len(drugs), nrows * cols):
        axes[j // cols][j % cols].set_visible(False)
    for row in axes:
        row[0].set_ylabel("Error rate (%)", fontsize=8)
    fig.legend(loc="upper right", frameon=False, fontsize=8, ncol=2, bbox_to_anchor=(0.995, 0.965))
    fig.text(0.01, 0.995, "The same genotype predictor scores differently under each breakpoint version",
             fontsize=10, fontweight="bold", va="top")
    fig.text(0.01, 0.005, "AMRFinderPlus calls (NCBI Pathogen Detection) vs re-interpreted MICs; K. pneumoniae, "
             "identical genome-drug pairs under every version.\nVery major = missed resistance; major = false resistance; bars are Wilson 95% "
             "intervals.\nAntibiotics whose error rate moves by at least 2 percentage points; y-axis scales "
             "differ by panel.", fontsize=6.5, color=INK_2)
    fig.tight_layout(rect=(0, 0.07, 1, 0.92))
    for ext in ("png", "svg"):
        fig.savefig(out / f"figure3_genotype_errors.{ext}", dpi=300, metadata=METADATA[ext])
    plt.close(fig)


# single-hue sequential ramp (light -> dark), from the validated palette
RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
        "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
ZERO = "#f0efec"


def figure_heatmap(rob, out: Path, scenario="record"):
    from matplotlib.colors import LinearSegmentedColormap

    rows = [r for r in rob if r["scenario"] == scenario and r["boot_lo"] != ""]
    versions = sorted({r["version_a"] for r in rows} | {r["version_b"] for r in rows}, key=_vkey)
    steps = list(zip(versions, versions[1:]))
    cols = steps + [(versions[0], versions[-1])]
    cell = {(r["antibiotic"], r["version_a"], r["version_b"]): r for r in rows}
    drugs = sorted({r["antibiotic"] for r in rows if r["antibiotic"] != "ALL" and int(r["flips"]) > 0},
                   key=lambda ab: -float(cell[(ab, versions[0], versions[-1])]["flip_rate"]))
    n_all = len({r["antibiotic"] for r in rows if r["antibiotic"] != "ALL"})
    names = [f"All {n_all} antibiotics"] + [_pretty(d) for d in drugs]
    keys = ["ALL"] + drugs
    vmax = max(float(cell[(k, a, b)]["flip_rate"]) for k in keys for a, b in cols if (k, a, b) in cell)
    cmap = LinearSegmentedColormap.from_list("blue", RAMP)

    fig, ax = plt.subplots(figsize=(7.0, 0.32 * len(keys) + 1.6))
    for i, k in enumerate(keys):
        y = i + (0 if i == 0 else 0.35)                     # gap under the ALL row
        for j, (a, b) in enumerate(cols):
            x = j + (0.35 if j == len(cols) - 1 else 0)    # gap before the overall column
            r = cell.get((k, a, b))
            if r is None:
                face, label, ink = "white", "", INK_2
            else:
                v, flips = float(r["flip_rate"]), int(r["flips"])
                face = ZERO if flips == 0 else cmap(min(1.0, 0.08 + 0.92 * v / vmax))
                dark = flips > 0 and v / vmax > 0.45
                ink = "white" if dark else (INK_2 if flips == 0 else INK)
                label = "0" if flips == 0 else (f"{100 * v:.1f}" if v >= 0.001 else "<0.1")
            ax.add_patch(plt.Rectangle((x + 0.03, y + 0.04), 0.94, 0.92, facecolor=face, edgecolor="none"))
            ax.text(x + 0.5, y + 0.5, label, ha="center", va="center", fontsize=7.5, color=ink)
    ax.set_xlim(0, len(cols) + 0.35)
    ax.set_ylim(len(keys) + 0.35, 0)
    ax.set_xticks([j + 0.5 + (0.35 if j == len(cols) - 1 else 0) for j in range(len(cols))],
                  [f"v{a}\n\u2192 v{b}" for a, b in steps] + [f"Overall\nv{versions[0]} \u2192 v{versions[-1]}"],
                  fontsize=7.5)
    ax.xaxis.tick_top()
    ax.set_yticks([i + 0.5 + (0 if i == 0 else 0.35) for i in range(len(keys))], names)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.text(0.01, 0.985, "Share of records whose S/I/R category changes, by EUCAST revision",
             fontsize=10, fontweight="bold", va="top")
    fig.text(0.01, 0.01, f"K. pneumoniae MICs from BV-BRC, main analysis; % of records decidable under both "
             f"versions. Darker = larger share; 0 = no change.\nRows: antibiotics with any change; "
             f"the other {n_all - len(drugs)} antibiotics never changed. v11.0 and v12.0 were not available.",
             fontsize=6.5, color=INK_2)
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    for ext in ("png", "svg"):
        fig.savefig(out / f"figure1_change_heatmap.{ext}", dpi=300, metadata=METADATA[ext])
    plt.close(fig)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("results", help="folder holding flip_rates.csv and as_reported.csv")
    p.add_argument("-o", "--output", required=True)
    args = p.parse_args(argv)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    figure_1(_read(Path(args.results) / "flip_rates.csv"), out)
    figure_2(_read(Path(args.results) / "as_reported.csv"), out)
    rob = Path(args.results) / "robustness.csv"
    if rob.exists():
        figure_heatmap(_read(rob), out)
    genotype = Path(args.results) / "genotype_errors.csv"
    if genotype.exists():
        figure_3(_read(genotype), out)
    print(f"wrote figures to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
