#!/usr/bin/env python
"""
Where does a generation actually go?

The per-generation budget is known except for one term. Measured on the
production node at pop 100: ESM3 folding 5.8 s, HemoPI2 57.7 s, pypsique
0.78 s, in-process MACREL 0.015 s -- and ~24 s unaccounted for, which is
~0.22 s per candidate. That term scales with population, is CPU-bound, and
once hemolysis is sampled it becomes the largest single cost. It also decides
whether a GPU is worth using: a GPU accelerates folding, which is 6.5% of a
generation, and nothing else.

This runs a short evolution under cProfile with hemolysis switched off, so the
remaining time is the unexplained term plus folding, and prints the functions
that own it.

    python profile_generation.py -ps 32 -ng 6       # use at least 6 generations

USE AT LEAST 6 GENERATIONS. At 2-3 the ESM3 weight load dominates and the
steady-state cost is invisible: ignore `uniform_`, `get_storage_from_record`,
`_imp.create_dynamic` and `_ssl._SSLSocket.read`, which are model loading and
the HuggingFace fetch, and happen once per run rather than once per generation.

What matters is what recurs per candidate. If the top steady-state entries are
torch operations, a GPU will help. If they are `_posixsubprocess.fork_exec`,
`select.poll`, `communicate` (psique or gzip), file writes or pandas, a GPU
will not help and the fix is in this repository, not in different hardware.
"""
import cProfile
import os
import pstats
import runpy
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    argv = sys.argv[1:]
    if not any(a in ("-ps", "--pop_size") for a in argv):
        argv += ["-ps", "8"]
    if not any(a in ("-ng", "--num_generations") for a in argv):
        argv += ["-ng", "3"]

    init = os.path.join(HERE, "init", "init_random.faa")
    if not os.path.exists(init):
        sys.exit(f"no init file at {init} -- run from the repository root")

    out = tempfile.mkdtemp(prefix="pfes_profile_")
    # --hemo-every larger than the run skips HemoPI2 entirely except the final
    # generation, so its fixed ~58 s does not swamp everything else.
    sys.argv = ["pfes.py", "--start", "file", "--start-file", init,
                "--norepeat", "--hemo-every", "9999", "-o", out] + argv
    print(f"profiling: {' '.join(sys.argv)}\n")

    stats_path = os.path.join(out, "profile.out")
    try:
        cProfile.runctx(
            "runpy.run_path(os.path.join(HERE, 'pfes.py'), run_name='__main__')",
            {"runpy": runpy, "os": os, "HERE": HERE}, {}, stats_path)
    except SystemExit:
        pass                      # pfes.py calls sys.exit on completion

    print("\n" + "=" * 78)
    print("TOP 25 BY SELF TIME  (tottime = time in the function itself)")
    print("=" * 78)
    pstats.Stats(stats_path).sort_stats("tottime").print_stats(25)
    print("=" * 78)
    print("TOP 15 BY CUMULATIVE TIME")
    print("=" * 78)
    pstats.Stats(stats_path).sort_stats("cumulative").print_stats(15)
    print(f"\nraw profile: {stats_path}")


if __name__ == "__main__":
    main()
