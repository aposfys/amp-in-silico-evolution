# Runbook — how to run a series, and what it costs

Measured on the production host, not estimated. Every figure here comes from the
v4 series of 30 August 2026: six runs, pop 100 × 600 generations, on an
RTX 5090 with the CPU work capped at 12 cores.

For what the series *found*, see [`RESULTS.md`](RESULTS.md). For what each term
of the objective does, [`OBJECTIVE.md`](OBJECTIVE.md).

## Cost

| arm | run | wall | s / generation | ms / candidate |
|---|---|---|---|---|
| `main` | random | 2072 s | 3.45 | 34.5 |
| `main` | fragments | 2056 s | 3.43 | 34.3 |
| `main` | orfs | 2048 s | 3.41 | 34.1 |
| `control-fold-only` | random | 2231 s | 3.72 | 37.2 |
| `control-fold-only` | fragments | 2280 s | 3.80 | 38.0 |
| `control-fold-only` | orfs | 2286 s | 3.81 | 38.1 |

**Six runs, 3 h 36 m, serial.** A single arm is under two hours.

Two things worth reading off that table. The spread within an arm is under 1 %,
so nothing was contending for the machine — if a run comes in noticeably slower,
something else is on the GPU or the cores. And the control arm costs only 10 %
more despite chains growing from 26 to 57–65 residues, so folding cost is far
less length-sensitive here than expected.

### Against the old numbers

The previous HPC series ran at **88.6 s/generation**. This is **26× faster**,
and the difference is not the GPU alone:

| | HPC (v3) | GPU (v4) |
|---|---|---|
| folding | CPU, 40 threads | RTX 5090 |
| MACREL | subprocess, model reload per generation | in-process, loaded once |
| HemoPI2 | every generation | `--hemo-every 25` |
| dedup | pandas scan, `pop² × gen²` | position index |

Every earlier cost estimate in this file has been deleted rather than corrected.
Three were badly wrong and one was wrong in a way worth remembering: a term of
**~0.22 s per candidate** was attributed to unexplained CPU work and used to
argue a GPU could not help. It does not exist. Profiling on this
host puts the non-folding per-candidate work at **6.4 ms**, and its top entries
are entirely torch and ESM3 operations — no subprocess, no poll, no pandas. The
0.22 s was contention on a shared 48-core node, measured by subtraction and
mistaken for a property of the code.

### Population

Cost is close to linear in population once HemoPI2 is sampled, so pop 400 should
land near 13–14 s/generation, or **~2.2 h a run** — a full six-run series in
about half a day. That is a projection from the per-candidate figure above, not
a measurement; run one arm before committing to six.

This reopens a question that was closed on cost grounds. Pop 100 was chosen
because a series was believed to take days; at 3.4 s/generation it does not, and
a larger population reduces the drift that a six-run design cannot absorb. The
counter-argument is that the corrected objective already converges tightly —
score spread across the three origins is 0.0019 ([`RESULTS.md`](RESULTS.md)) —
so extra population may buy little.

## Running a series

```bash
cd /path/to/amp-in-silico-evolution
conda activate pfes_amps
./preflight.sh                      # must pass; every launcher sources it
MAXTOK=4096 ./run_v4.sh main        # 3 runs, ~1h45m
```

Then the control arm, in a worktree so both objectives stay checked out:

```bash
git worktree add ../amps-ctrl control-fold-only
cd ../amps-ctrl
MAXTOK=8192 ./run_v4.sh ctrl        # 3 runs, ~2h05m
```

Use `tmux`; a series outlasts an SSH session. `PFES_CORES` overrides the
12-core cap, `PFES_CPUSET` the core list.

**The control branch must be merged from `main` before every series.** Its
`score.py` and guards have to be identical; the only intended difference is the
`score = np.prod([mean_plddt, ptm])` block. A stale control branch is how the
v2/v3 series ended up comparing arms scored by different tools.

### Monitoring

```bash
cat results/v4/main.master.log                              # START/END with durations
tail -n1 results/v4/main/random/progress.log | cut -f1      # current generation
tail -f results/v4/main/random.console.log                  # live, Ctrl-C to exit
```

Do not pipe `tail -f` into anything that waits for EOF — it never arrives.

### After a run, before believing anything

```bash
for r in random fragments orfs; do echo -n "$r: "; awk -F'\t' '$1=="gndx"{for(i=1;i<=NF;i++) if($i=="amp_src") c=i; next} c && $1 ~ /^gndx[0-9]/ {n[$c]++} END{for(k in n) printf "%s=%d ", k, n[k]; print ""}' results/v4/main/$r/progress.log; done
```

`macrel` only. Any `proxy` count means MACREL was undefined for those rows —
below 10 or above 100 residues — and the surrogate answered instead, changing
what `amp_prob` measures partway through. It matters most on the control arm,
which has no length penalty.

Row counts fall slightly short of `pop × generations`, concentrated in the first
few generations and one just after the strong-selection switch. That is
`--norepeat` dropping duplicates while the population is still near-clonal, and
is expected. In v4 the sORF arm started at 61 unique of 100 against 88 for
fragments and 96 for random — a real property of the origin, since translated
small ORFs carry repeated segments that single mutations collide on.

## Environment

`setup_gpu.sh` builds it from nothing and ends by running `preflight.sh`. It has
been run end to end on an RTX 5090 (sm_120, driver CUDA 13.2). Each of the
following stopped a launch at least once and is now guarded:

| failure | guard |
|---|---|
| `conda` is a shell function, invisible to a child script | locates `conda.sh` via `CONDA_EXE`/`CONDA_PREFIX` |
| Anaconda ToS gate on the `defaults` channels | builds with `-c conda-forge --override-channels` |
| `esm` imports `httpx` but does not declare it | pinned in `requirements.txt` |
| torchvision from PyPI mismatches torch's CUDA major | installed from the same index; `--no-deps` on repair |
| `numpy` 2 breaks HemoPI2's pickled models | asserted and repaired after every install step |
| `onnxruntime` > 1.25.1 makes MACREL return raw decision values | same |
| `HF_TOKEN` set to a documentation placeholder | rejected by length and content |
| HemoPI2 models 1, 2 and 4 return uncalibrated or constant scores | model 3, plus an orientation check |

That last one is the pattern to watch for. Three separate times a classifier has
returned something that was not a calibrated probability and looked entirely
plausible in a log. The `[0, 1]` range check in `score.py` caught all three;
clamping instead would have written `1.0000` for every candidate of every
generation. `preflight.sh` now probes both classifiers end to end — MACREL on
magainin 2, HemoPI2 on magainin 2 *and* melittin, because range alone passes a
model that returns a constant.

## Analysis

```bash
python visual_pfes.py -l results/v4/main/<run>/progress.log \
       -s results/v4/main/<run>/structures -o results/v4/main/<run>/analysis --notraj
python analysis/compare_runs.py results/v4/*/*/ -o results/v4/comparison
python analysis/score_posthoc.py results/v4/*/*/ -o results/v4/comparison --lineage
```

`score_posthoc.py --lineage` is now the slow step, and deliberately so: HemoPI2
model 3 is an ESM2 language model reloaded per subprocess call, where model 1
was a random forest. Inside a run it barely registers, because `--hemo-every 25`
fires on 24 of 600 generations. Post hoc it dominates. `PFES_HEMO_MODEL=1` is
much faster and writes uncalibrated values, so use it only when the hemolysis
column does not matter.

AMPlify needs its own environment (Python 3.6, old TensorFlow) and is the
independent audit that turns a MACREL score of 0.99 from a restatement of the
objective into evidence. It has not yet been run against the v4 series.

## What is archived

One compressed archive per run is kept by the author: `progress.log`, the
console log with the exact invocation, and the final generation's structures.
These are **not published in this repository** and are available on request —
apostolosfysekidis1@gmail.com. The 60,000 intermediate structures per run are
not kept at all — 237 MB each, 1.4 GB for a series, re-derivable from the
sequences in `progress.log`, and needed for no figure or table.

Runs are not individually reproducible: neither this fork nor upstream PFES
seeds any RNG. The starting population is fixed by file and recorded by
checksum in every run banner; everything else is a draw. The evidence this
design admits is replication across independent runs.
