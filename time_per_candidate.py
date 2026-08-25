#!/usr/bin/env python
"""
Time the per-candidate work directly. Seconds, not a profiled run.

extract_results does exactly five things per candidate. Four are already
measured; roughly 0.22 s of the ~0.29 s is unaccounted for, and that term is
the largest single cost once hemolysis is sampled -- and it is CPU work, so a
GPU does not touch it. This times each piece on a real structure from a
previous run, with no evolution loop involved.

    python time_per_candidate.py
    python time_per_candidate.py results/v3/sorf-r1/structures   # explicit path
"""
import glob
import gzip
import os
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd                                     # noqa: E402
from score import get_nconts                            # noqa: E402
from psique import pypsique                             # noqa: E402

N = 20


def bench(label, fn, n=N):
    fn()                                                # warm
    t = time.perf_counter()
    for _ in range(n):
        fn()
    ms = (time.perf_counter() - t) / n * 1000
    print(f"  {label:34s} {ms:8.2f} ms")
    return ms


def main():
    roots = sys.argv[1:] or sorted(glob.glob("results/*/*/structures")) \
                          + sorted(glob.glob("pfes-results/results/*/*/structures"))
    files = []
    for r in roots:
        files = sorted(glob.glob(os.path.join(r, "*.pdb.gz")))
        if files:
            break
    if not files:
        sys.exit("no structures found -- pass a path to a run's structures/ directory")

    pdb_txt = gzip.open(files[-1], "rt").read()
    nres = sum(1 for l in pdb_txt.splitlines()
               if l.startswith("ATOM  ") and l[12:16].strip() == "CA")
    print(f"structure: {files[-1]}  ({nres} residues)\n")

    tmp = tempfile.mkdtemp(prefix="pfes_time_")
    path = os.path.join(tmp, "x.pdb")
    cols = {c: 0.5 for c in
            ("gndx seq_len prot_len_penalty max_alpha_penalty max_beta_penalty ptm "
             "mean_plddt num_conts iplddt num_inter_conts score").split()}
    cols.update({c: "x" * 27 for c in
                 "id sel_mode amp_prob amp_src hemo_prob sequence mutation prev_id ss".split()})

    total = 0.0
    total += bench("write .pdb",
                   lambda: open(path, "wb").write(pdb_txt.encode()))
    total += bench("get_nconts (Eq.5 contacts)",
                   lambda: get_nconts(pdb_txt, "A", 6.0, 0.5))
    total += bench("pypsique (secondary structure)",
                   lambda: pypsique(pdb_txt, "A"))
    total += bench("pd.DataFrame(1 row)",
                   lambda: pd.DataFrame(cols, index=[0]))
    total += bench("os.system('gzip ... &')",
                   lambda: os.system(f"gzip -f '{path}' >/dev/null 2>&1 &"))

    print(f"\n  {'TOTAL per candidate':34s} {total:8.2f} ms")
    print(f"  {'observed (v3, by subtraction)':34s} {290.0:8.2f} ms")
    gap = 290.0 - total
    print(f"  {'UNACCOUNTED':34s} {gap:8.2f} ms  ({100*gap/290:.0f}%)")
    print("\n  If TOTAL is close to 290 ms, the slow piece is named above and")
    print("  fixable here. If a large gap remains, the cost is in the folding")
    print("  call or the selection step, not in extract_results.")


if __name__ == "__main__":
    main()
