#!/usr/bin/env python3
"""
Compare several PFES runs against each other.

visual_pfes.py describes one run. Nothing in the repository looks at more than
one at a time, so every cross-run statement has to be assembled by hand. This
script closes that gap: point it at any number of run directories and it groups
them by objective, plots each metric as a mean with an across-repeat band, and
writes the per-run numbers out as tables.

    python analysis/compare_runs.py results/final/*/ -o results/comparison

Works with a single run per arm (it just draws no band and says so), so it can
be developed against two runs and pointed at ten later.

Arms are identified from the data, not from the directory name: the two
objectives are reconstructed from the logged columns and whichever reproduces
the logged score is the arm. A name that disagrees with the data is reported.

Author: Apostolos Fysekidis
"""

import argparse
import os
import sys
from collections import defaultdict

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Okabe-Ito, colour-vision safe. Structured blue, fold-only orange, as in the
# thesis figures so the two sets can sit side by side.
C = {"structured": "#0072B2", "foldonly": "#E69F00"}
C_INK, C_GRID, C_MUTE = "#22262B", "#E5E8EC", "#9AA1AE"
LABEL = {"structured": "Structured", "foldonly": "Fold-only"}


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def load_log(path):
    """Read a progress.log, skipping the hash-prefixed preamble."""
    rows, hdr = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = parts
                continue
            if len(parts) == len(hdr):
                rows.append(dict(zip(hdr, parts)))
    return rows


def fnum(row, key, default=np.nan):
    try:
        return float(row[key])
    except (KeyError, ValueError, TypeError):
        return default


def gen_of(row):
    g = row.get("gndx", "")
    return int(g[4:]) if g.startswith("gndx") and g[4:].isdigit() else -1


def detect_arm(rows, name_hint):
    """Identify the objective by reconstructing it from the logged columns.

    The structured objective multiplies in the three penalties and the contact
    booster; the fold-only one does not. Whichever reproduces `score` is the arm.
    """
    errs = {"structured": [], "foldonly": []}
    for r in rows[:3000]:
        L, N = fnum(r, "seq_len"), fnum(r, "num_conts")
        if not np.isfinite(L) or L <= 0:
            continue
        base = fnum(r, "ptm") * fnum(r, "mean_plddt") * fnum(r, "amp_prob")
        struct = (base * fnum(r, "prot_len_penalty") * fnum(r, "max_alpha_penalty")
                  * fnum(r, "max_beta_penalty") * ((N + L) / L))
        s = fnum(r, "score")
        errs["structured"].append(abs(struct - s))
        errs["foldonly"].append(abs(base - s))
    if not errs["structured"]:
        return name_hint or "structured", float("nan")
    means = {k: float(np.nanmean(v)) for k, v in errs.items()}
    arm = min(means, key=means.get)
    if name_hint and name_hint != arm:
        sys.stderr.write(
            f"  note: directory name says '{name_hint}' but the logged scores "
            f"reconstruct as '{arm}' (residuals {means}); using the data.\n")
    return arm, means[arm]


def load_run(path):
    log = os.path.join(path, "progress.log")
    if not os.path.isfile(log):
        return None
    rows = load_log(log)
    if not rows:
        return None
    base = os.path.basename(path.rstrip("/"))
    hint = ("foldonly" if "foldonly" in base or "fold-only" in base
            else "structured" if "structured" in base else None)
    arm, resid = detect_arm(rows, hint)
    return {"name": base, "path": path, "arm": arm, "rows": rows, "resid": resid}


# --------------------------------------------------------------------------- #
# Aggregation
# --------------------------------------------------------------------------- #
def per_generation(rows, col, agg="mean"):
    """Collapse one logged column to one value per generation."""
    buckets = defaultdict(list)
    for r in rows:
        g = gen_of(r)
        if g >= 0:
            buckets[g].append(fnum(r, col))
    gens = np.array(sorted(buckets))
    f = {"mean": np.nanmean, "max": np.nanmax, "min": np.nanmin}[agg]
    with np.errstate(all="ignore"):
        vals = np.array([f(buckets[g]) for g in gens])
    return gens, vals


def stack_runs(runs, col, agg="mean"):
    """Align repeats on a common generation axis and return (gens, matrix)."""
    series = [per_generation(r["rows"], col, agg) for r in runs]
    if not series:
        return np.array([]), np.zeros((0, 0))
    n = min(len(g) for g, _ in series)
    gens = series[0][0][:n]
    return gens, np.vstack([v[:n] for _, v in series])


def longest_helix(ss):
    """Longest run of H. PFES scores G/F/T/P as coil, so only H counts."""
    best = run = 0
    for ch in ss or "":
        run = run + 1 if ch == "H" else 0
        best = max(best, run)
    return best


def best_of_final_generation(rows):
    last = max(gen_of(r) for r in rows)
    return max((r for r in rows if gen_of(r) == last),
               key=lambda x: fnum(x, "score", -1))


def global_best(rows):
    return max(rows, key=lambda x: fnum(x, "score", -1))


def lineage(rows, start=None):
    """An ancestral line, oldest first, walked back through prev_id.

    Defaults to the best individual of the FINAL generation, which is what
    visual_pfes.extract_lineage does and therefore what lineage.tsv contains.
    That is deliberately not the same as the globally best individual: a run
    can peak mid-trajectory and drift down afterwards, as the structured arm
    does (best at generation 241, final generation scores lower). The global
    best still lies on this line, it is simply not its endpoint.

    Note this is not the number of distinct parents in the log either. Most
    individuals leave no descendant in the winning line at all.

    Counting caveat: visual_pfes.extract_lineage writes its starting node twice,
    so lineage.tsv has one more row than the line has members. This function
    returns unique members, so it reports one fewer than `wc -l lineage.tsv`.
    """
    by_id = {r.get("id"): r for r in rows}
    node = start or best_of_final_generation(rows)
    chain, seen = [], set()
    while node is not None and node.get("id") not in seen:
        seen.add(node.get("id"))
        chain.append(node)
        node = by_id.get(node.get("prev_id"))
    return list(reversed(chain))


# --------------------------------------------------------------------------- #
# Plot helpers
# --------------------------------------------------------------------------- #
def style(ax, xlabel=None, ylabel=None):
    ax.grid(color=C_GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(C_MUTE)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors="#5C6472", labelsize=8, length=3)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=8.5, color="#5C6472")
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=8.5, color="#5C6472")


def band(ax, groups, col, agg="mean"):
    """Mean across repeats with a min-max band, one colour per arm."""
    for arm, runs in sorted(groups.items()):
        gens, M = stack_runs(runs, col, agg)
        if not len(gens):
            continue
        mean = np.nanmean(M, axis=0)
        if M.shape[0] > 1:
            ax.fill_between(gens, np.nanmin(M, axis=0), np.nanmax(M, axis=0),
                            color=C[arm], alpha=0.16, lw=0, zorder=2)
            for row in M:                      # individual repeats, faint
                ax.plot(gens, row, color=C[arm], lw=0.5, alpha=0.35, zorder=3)
        ax.plot(gens, mean, color=C[arm], lw=1.8, zorder=4)


def legend(fig, groups, y=0.01):
    handles = [Line2D([], [], color=C[a], lw=1.8,
                      label=f"{LABEL[a]} (n={len(r)})")
               for a, r in sorted(groups.items())]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, y))


def save(fig, outdir, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(outdir, f"{stem}.{ext}"),
                    dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  {stem}")


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def fig_trajectories(groups, outdir):
    panels = [("seq_len", "Chain length (residues)"),
              ("score", "Fitness score")]
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.8))
    for ax, (col, lab) in zip(axes, panels):
        band(ax, groups, col)
        style(ax, "Generation", lab)
    axes[0].axhline(30, color="#5C6472", lw=0.9, ls=(0, (4, 3)))
    axes[0].annotate(r"$L_0$", xy=(1.0, 30), xycoords=("axes fraction", "data"),
                     fontsize=7.5, color="#5C6472", va="center", ha="left")
    legend(fig, groups, y=-0.04)
    fig.subplots_adjust(bottom=0.30, wspace=0.28)
    save(fig, outdir, "cmp_trajectories")


def fig_inputs(groups, outdir):
    panels = [("mean_plddt", "pLDDT"), ("ptm", "pTM"),
              ("num_conts", "Contacts"), ("amp_prob", "MACREL probability")]
    fig, axes = plt.subplots(1, 4, figsize=(7.4, 2.0))
    for ax, (col, lab) in zip(axes, panels):
        band(ax, groups, col)
        style(ax, "Generation", lab)
    legend(fig, groups, y=-0.06)
    fig.subplots_adjust(bottom=0.36, wspace=0.42)
    save(fig, outdir, "cmp_inputs")


def fig_helix(groups, outdir):
    fig, ax = plt.subplots(figsize=(4.8, 2.8))
    for arm, runs in sorted(groups.items()):
        series = []
        for r in runs:
            buckets = defaultdict(list)
            for row in r["rows"]:
                g = gen_of(row)
                if g >= 0:
                    buckets[g].append(longest_helix(row.get("ss", "")))
            gens = np.array(sorted(buckets))
            series.append((gens, np.array([np.mean(buckets[g]) for g in gens])))
        n = min(len(g) for g, _ in series)
        gens = series[0][0][:n]
        M = np.vstack([v[:n] for _, v in series])
        if M.shape[0] > 1:
            ax.fill_between(gens, M.min(0), M.max(0), color=C[arm],
                            alpha=0.16, lw=0)
        ax.plot(gens, M.mean(0), color=C[arm], lw=1.8)
    ax.axhline(20, color="#5C6472", lw=0.9, ls=(0, (4, 3)))
    style(ax, "Generation", "Longest α-helix (residues)")
    legend(fig, groups, y=-0.02)
    fig.subplots_adjust(bottom=0.32)
    save(fig, outdir, "cmp_helix")


def fig_hemolysis(groups, outdir):
    """Only meaningful if HemoPI2 actually ran (PFES_SKIP_HEMO unset)."""
    allzero = all(abs(fnum(row, "hemo_prob", 0.0)) < 1e-9
                  for runs in groups.values() for r in runs
                  for row in r["rows"][::37])
    if allzero:
        print("  cmp_hemolysis: skipped, hemo_prob is 0.0 everywhere "
              "(PFES_SKIP_HEMO was set for these runs)")
        return
    fig, ax = plt.subplots(figsize=(4.8, 2.8))
    band(ax, groups, "hemo_prob")
    style(ax, "Generation", "HemoPI2 probability")
    legend(fig, groups, y=-0.02)
    fig.subplots_adjust(bottom=0.32)
    save(fig, outdir, "cmp_hemolysis")


def endpoint_table(groups):
    """Per-run endpoint metrics: the input to the variance figure and the TSV."""
    out = []
    for arm, runs in sorted(groups.items()):
        for r in runs:
            last = max(gen_of(x) for x in r["rows"])
            final = [x for x in r["rows"] if gen_of(x) == last]
            best = global_best(r["rows"])   # the peptide actually reported
            seq = best.get("sequence", "")
            out.append({
                "run": r["name"], "arm": arm,
                "generations": last + 1,
                "final_len_mean": float(np.mean([fnum(x, "seq_len") for x in final])),
                "best_score": fnum(best, "score"),
                "winner_len": fnum(best, "seq_len"),
                "winner_helix": longest_helix(best.get("ss", "")),
                "winner_amp": fnum(best, "amp_prob"),
                "winner_hemo": fnum(best, "hemo_prob"),
                "winner_charge_KR_DE": (sum(seq.count(a) for a in "KR")
                                        - sum(seq.count(a) for a in "DE")),
                "lineage_nodes": len(lineage(r["rows"])),
                "winner": seq,
            })
    return out


def fig_variance(table, groups, outdir):
    """Between-arm separation read against within-arm scatter."""
    metrics = [("final_len_mean", "Final chain length"),
               ("winner_helix", "Winner longest α-helix"),
               ("best_score", "Best fitness")]
    arms = sorted(groups)
    fig, axes = plt.subplots(1, len(metrics), figsize=(7.4, 2.6))
    for ax, (key, lab) in zip(axes, metrics):
        for i, arm in enumerate(arms):
            vals = [t[key] for t in table if t["arm"] == arm]
            x = np.full(len(vals), i, dtype=float)
            x += np.linspace(-0.09, 0.09, len(vals)) if len(vals) > 1 else 0
            ax.scatter(x, vals, s=26, color=C[arm], zorder=4,
                       edgecolor="white", linewidth=0.7)
            if len(vals) > 1:
                m, sd = float(np.mean(vals)), float(np.std(vals, ddof=1))
                ax.errorbar(i, m, yerr=sd, color=C[arm], lw=1.4, capsize=5,
                            zorder=3, marker="_", markersize=16)
        ax.set_xticks(range(len(arms)))
        ax.set_xticklabels([LABEL[a] for a in arms], fontsize=8)
        ax.set_xlim(-0.5, len(arms) - 0.5)
        style(ax, None, lab)
    fig.subplots_adjust(wspace=0.42, bottom=0.16)
    save(fig, outdir, "cmp_variance")


def write_tables(table, groups, outdir):
    cols = ["run", "arm", "generations", "final_len_mean", "best_score",
            "winner_len", "winner_helix", "winner_amp", "winner_hemo",
            "winner_charge_KR_DE", "lineage_nodes", "winner"]
    p = os.path.join(outdir, "per_run.tsv")
    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for t in table:
            fh.write("\t".join(f"{t[c]:.3f}" if isinstance(t[c], float)
                               else str(t[c]) for c in cols) + "\n")
    print("  per_run.tsv")

    p = os.path.join(outdir, "summary.tsv")
    keys = [c for c in cols if c not in ("run", "arm", "winner")]
    with open(p, "w") as fh:
        fh.write("arm\tn\t" + "\t".join(f"{k}_mean\t{k}_sd" for k in keys) + "\n")
        for arm in sorted(groups):
            vals = [t for t in table if t["arm"] == arm]
            cells = []
            for k in keys:
                v = np.array([x[k] for x in vals], dtype=float)
                sd = float(np.std(v, ddof=1)) if len(v) > 1 else float("nan")
                cells += [f"{np.mean(v):.3f}", f"{sd:.3f}"]
            fh.write(f"{arm}\t{len(vals)}\t" + "\t".join(cells) + "\n")
    print("  summary.tsv")


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="run directories (each with progress.log)")
    ap.add_argument("-o", "--outdir", default="comparison")
    args = ap.parse_args()

    runs = [r for r in (load_run(p) for p in args.runs) if r]
    if not runs:
        sys.exit("no readable runs (each needs a progress.log)")

    groups = defaultdict(list)
    for r in runs:
        groups[r["arm"]].append(r)

    print(f"\n{len(runs)} run(s):")
    for arm in sorted(groups):
        for r in groups[arm]:
            g = max(gen_of(x) for x in r["rows"]) + 1
            print(f"  {r['name']:<34} {LABEL[arm]:<11} {g:>4} gens  "
                  f"{len(r['rows']):>6,} rows  (objective residual {r['resid']:.4f})")
    if any(len(v) < 2 for v in groups.values()):
        print("\n  note: at least one arm has a single run, so no band is drawn "
              "for it and its standard deviation is undefined.")

    os.makedirs(args.outdir, exist_ok=True)
    print(f"\nwriting to {args.outdir}/")
    fig_trajectories(groups, args.outdir)
    fig_inputs(groups, args.outdir)
    fig_helix(groups, args.outdir)
    fig_hemolysis(groups, args.outdir)
    table = endpoint_table(groups)
    fig_variance(table, groups, args.outdir)
    write_tables(table, groups, args.outdir)
    print()


if __name__ == "__main__":
    main()
