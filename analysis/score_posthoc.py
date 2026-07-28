#!/usr/bin/env python3
"""
Score evolved candidates after evolution, with a classifier that did not drive it.

MACREL drives selection, so a high MACREL probability on an evolved peptide is
partly a restatement of the objective rather than independent evidence. AMPlify
is architecturally unrelated (attentive BiLSTM against a random forest on
descriptors) and never sees the search, so its agreement or disagreement is
informative in a way MACREL's own score cannot be.

    python analysis/score_posthoc.py results/final/*/ -o results/comparison

Writes posthoc_scores.tsv and two figures:

  posthoc_classifiers   MACREL against AMPlify for every candidate. Points far
                        below the diagonal are peptides the objective rewarded
                        and an independent model does not recognise, which is
                        what a single-classifier objective conceals.
  posthoc_charge_length chain length against net charge, with natural AMPs from
                        the DRAMP-derived database as background, showing
                        whether candidates land inside the natural region.

AMPlify needs its own environment (python 3.6, old TensorFlow) and cannot share
the ESM3 one:

    export AMPLIFY_CMD="conda run -n amplify AMPlify"

Without it the script still runs, reports MACREL and HemoPI2, and says which
figure it could not draw.

Author: Apostolos Fysekidis
"""

import argparse
import os
import shlex
import subprocess
import sys
import tempfile
import glob as _glob
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                      # compare_runs.py, alongside
sys.path.insert(0, os.path.dirname(_HERE))     # score.py, at the repo root
from compare_runs import (load_run, gen_of, fnum, global_best,   # noqa: E402
                          longest_helix, C, LABEL, style, save)

# Net charge at pH 7.4. Histidine is counted at +0.1, its approximate
# protonated fraction at that pH; this reproduces the +5.2 reported for the
# structured winner (3 K, 3 R, 1 E, 2 H).
CHARGE = {**{a: 1.0 for a in "KR"}, **{a: -1.0 for a in "DE"}, "H": 0.1}


def net_charge(seq):
    return sum(CHARGE.get(a, 0.0) for a in seq)


def hydrophobic_fraction(seq):
    return sum(a in "AVLIMFWYC" for a in seq) / max(len(seq), 1)


# --------------------------------------------------------------------------- #
# AMPlify
# --------------------------------------------------------------------------- #
def amplify_scores(sequences, model="balanced"):
    """{sequence: probability}, or {} if AMPlify is unavailable.

    Mirrors amplify_score_batch on the fitness-amplify-* branches, including the
    header-drift tolerance: the probability column has been named differently
    across AMPlify releases.
    """
    if not sequences:
        return {}
    import pandas as pd
    cmd = shlex.split(os.environ.get("AMPLIFY_CMD", "AMPlify"))
    try:
        with tempfile.TemporaryDirectory() as tmp:
            fa = os.path.join(tmp, "in.fa")
            with open(fa, "w") as fh:
                for i, s in enumerate(sequences):
                    fh.write(f">seq{i}\n{s}\n")
            r = subprocess.run(cmd + ["-s", fa, "-od", tmp, "-of", "tsv",
                                      "-m", model],
                               capture_output=True, text=True, timeout=1800)
            if r.returncode != 0:
                raise RuntimeError(r.stderr.strip() or r.stdout.strip())
            outs = sorted(_glob.glob(os.path.join(tmp, "*.tsv")),
                          key=os.path.getmtime)
            if not outs:
                raise FileNotFoundError("no AMPlify .tsv produced")
            df = pd.read_csv(outs[-1], sep="\t")
            col = next((c for c in df.columns
                        if "prob" in str(c).lower() and "score" in str(c).lower()), None)
            col = col or next((c for c in df.columns
                               if "prob" in str(c).lower()), None)
            if col is None:
                raise ValueError(f"no probability column in {list(df.columns)}")
            seqcol = next((c for c in df.columns
                           if str(c).strip().lower() == "sequence"), None)
            if seqcol is None:
                return {s: float(v) for s, v in zip(sequences, df[col])}
            return {str(s): float(v) for s, v in zip(df[seqcol], df[col])}
    except FileNotFoundError:
        sys.stderr.write(
            '  AMPlify not found. Install it and set\n'
            '    export AMPLIFY_CMD="conda run -n amplify AMPlify"\n'
            '  Continuing without it.\n')
    except Exception as e:
        sys.stderr.write(f"  AMPlify failed ({type(e).__name__}: {e}).\n")
    return {}


def macrel_hemo(sequences):
    """MACREL probability and HemoPI2 risk, via the repository's own score.py."""
    try:
        import score
    except Exception as e:
        sys.stderr.write(f"  could not import score.py ({e})\n")
        return {}, {}
    out = score.macrel_score_batch(list(sequences))
    return ({s: v[0] for s, v in out.items()},
            {s: v[1] for s, v in out.items()})


# --------------------------------------------------------------------------- #
# Candidates
# --------------------------------------------------------------------------- #
def candidates(run, top_n):
    """The winner plus the best of the final generation, deduplicated."""
    rows = run["rows"]
    last = max(gen_of(r) for r in rows)
    final = sorted((r for r in rows if gen_of(r) == last),
                   key=lambda x: fnum(x, "score", -1), reverse=True)
    picked, seen = [], set()
    for r in [global_best(rows)] + final[:top_n]:
        s = r.get("sequence", "")
        if s and s not in seen:
            seen.add(s)
            picked.append(r)
    return picked


def natural_reference(root):
    """Natural AMPs for the background of the charge-length figure."""
    for p in (os.path.join(root, "data", "dbaasp.faa"),
              os.environ.get("PFES_AMP_DB", "")):
        if p and os.path.isfile(p):
            seqs = []
            for block in open(p).read().split(">")[1:]:
                s = "".join(block.split("\n")[1:]).strip().upper()
                if 8 <= len(s) <= 60 and set(s) <= set("ACDEFGHIKLMNPQRSTVWY"):
                    seqs.append(s)
            if seqs:
                return seqs, os.path.basename(p)
    # small curated fallback so the figure is never empty
    return ([
        "GIGKFLHSAKKFGKAFVGEIMNS",            # magainin 2
        "GIGAVLKVLTTGLPALISWIKRKRQQ",         # melittin
        "KWKLFKKIGAVLKVL",                    # part of cecropin-melittin
        "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES",  # LL-37
        "ILPWKWPWWPWRR",                      # indolicidin
        "RGGRLCYCRRRFCVCVGR",                 # protegrin-like
        "ACYCRIPACIAGERRYGTCIYQGRLWAFCC",     # human beta-defensin 1
        "GLLSKLWEEVEKVAGGIWEAFKR",            # unnamed helical AMP
    ], "curated fallback")


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def fig_classifiers(recs, outdir):
    if not any(r["amplify"] == r["amplify"] for r in recs):   # all NaN
        print("  posthoc_classifiers: skipped, no AMPlify scores")
        return
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    ax.plot([0, 1], [0, 1], color="#9AA1AE", lw=0.9, ls=(0, (4, 3)), zorder=1)
    for arm in sorted({r["arm"] for r in recs}):
        pts = [r for r in recs if r["arm"] == arm and r["amplify"] == r["amplify"]]
        ax.scatter([p["macrel"] for p in pts], [p["amplify"] for p in pts],
                   s=34, color=C[arm], alpha=0.85, edgecolor="white",
                   linewidth=0.6, zorder=3)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    style(ax, "MACREL probability (drove selection)", "AMPlify probability (independent)")
    fig.legend(handles=[Line2D([], [], marker="o", ls="none", ms=7,
                               mfc=C[a], mec="white", label=LABEL[a])
                        for a in sorted({r["arm"] for r in recs})],
               loc="lower center", ncol=2, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.01))
    fig.subplots_adjust(bottom=0.24, left=0.17)
    save(fig, outdir, "posthoc_classifiers")


def fig_charge_length(recs, natural, source, outdir):
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    nx = [len(s) for s in natural]
    ny = [net_charge(s) for s in natural]
    ax.scatter(nx, ny, s=9, color="#C7CCD4", alpha=0.55, lw=0, zorder=1)
    for arm in sorted({r["arm"] for r in recs}):
        pts = [r for r in recs if r["arm"] == arm]
        ax.scatter([p["length"] for p in pts], [p["charge"] for p in pts],
                   s=40, color=C[arm], edgecolor="white", linewidth=0.7, zorder=3)
    ax.axvline(30, color="#5C6472", lw=0.9, ls=(0, (4, 3)))
    style(ax, "Chain length (residues)", "Net charge at pH 7.4")
    fig.legend(handles=[Line2D([], [], marker="o", ls="none", ms=5.5,
                               mfc="#C7CCD4", mec="none",
                               label=f"natural AMPs (n={len(natural)})")] +
               [Line2D([], [], marker="o", ls="none", ms=7, mfc=C[a],
                       mec="white", label=LABEL[a])
                for a in sorted({r["arm"] for r in recs})],
               loc="lower center", ncol=3, frameon=False, fontsize=8,
               bbox_to_anchor=(0.5, -0.02))
    fig.subplots_adjust(bottom=0.30, left=0.13)
    save(fig, outdir, "posthoc_charge_length")
    print(f"  natural reference: {source}, {len(natural)} peptides")


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("-o", "--outdir", default="comparison")
    ap.add_argument("-n", "--top-n", type=int, default=10,
                    help="candidates from the final generation per run (default 10)")
    ap.add_argument("--model", default="balanced",
                    choices=["balanced", "imbalanced"])
    args = ap.parse_args()

    runs = [r for r in (load_run(p) for p in args.runs) if r]
    if not runs:
        sys.exit("no readable runs")

    recs = []
    for run in runs:
        for r in candidates(run, args.top_n):
            seq = r["sequence"]
            recs.append({"run": run["name"], "arm": run["arm"], "sequence": seq,
                         "length": len(seq), "charge": net_charge(seq),
                         "hydrophobic": hydrophobic_fraction(seq),
                         "helix": longest_helix(r.get("ss", "")),
                         "score": fnum(r, "score"),
                         "macrel_logged": fnum(r, "amp_prob")})
    seqs = sorted({r["sequence"] for r in recs})
    print(f"\n{len(runs)} run(s), {len(recs)} candidates, {len(seqs)} unique sequences")

    print("scoring...")
    mac, hemo = macrel_hemo(seqs)
    amp = amplify_scores(seqs, args.model)
    for r in recs:
        s = r["sequence"]
        r["macrel"] = mac.get(s, r["macrel_logged"])
        r["hemopi2"] = hemo.get(s, float("nan"))
        r["amplify"] = amp.get(s, float("nan"))

    os.makedirs(args.outdir, exist_ok=True)
    cols = ["run", "arm", "length", "charge", "hydrophobic", "helix", "score",
            "macrel", "amplify", "hemopi2", "sequence"]
    with open(os.path.join(args.outdir, "posthoc_scores.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(recs, key=lambda x: (-x["score"],)):
            fh.write("\t".join(f"{r[c]:.3f}" if isinstance(r[c], float)
                               else str(r[c]) for c in cols) + "\n")
    print("  posthoc_scores.tsv")

    fig_classifiers(recs, args.outdir)
    natural, source = natural_reference(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    fig_charge_length(recs, natural, source, args.outdir)

    ok = [r for r in recs if r["amplify"] == r["amplify"]]
    if ok:
        d = np.array([r["macrel"] - r["amplify"] for r in ok])
        print(f"\n  MACREL - AMPlify: mean {d.mean():+.3f}, "
              f"{(d > 0.3).sum()}/{len(d)} candidates where MACREL exceeds "
              f"AMPlify by more than 0.3")
    print()


if __name__ == "__main__":
    main()
