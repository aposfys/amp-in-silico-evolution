# Runbook — the corrected production series

Everything here is measured on this project's own runs, not estimated. The v2
series is an accidental controlled experiment: it executed the identical
pipeline with MACREL and HemoPI2 **not installed**, falling back to in-process
proxies, which isolates what those two subprocesses cost.

## Where the time actually goes

Per generation, v3 `sorf-r1`, pop 100, 600 generations, 14.77 h wall:

| Component | s/gen | share | scales with |
|---|---|---|---|
| ESM3 folding | 5.8 | **6.5 %** | population — *the only thing a GPU touches* |
| MACREL + HemoPI2 subprocesses | 60.6 | 68.4 % | per-invocation model reload, **not** population |
| ├─ MACREL (measured) | 0.35 | 0.4 % | subprocess; now in-process, 71.5× |
| └─ HemoPI2 (by subtraction) | ~60.2 | ~68 % | **the actual bottleneck — untouched** |
| PSIQUE + contacts + bookkeeping | 22.2 | 25.0 % | population (PSIQUE is one subprocess **per structure**) |

Derivation: three clean v2 runs averaged **28.0 s/generation** with proxies
(`fragments-r1` 28.7, `random-r2` 27.7, `fragments-r2` 27.6) against **88.6 s**
for v3 with the real classifiers, at the same 5.7–5.8 s folding cost.

### A GPU cannot fix this

```
GPU  5× faster at folding  →  whole run 1.06× faster
GPU 50× faster at folding  →  whole run 1.07×
GPU  ∞  faster at folding  →  whole run 1.07×
```

Amdahl's law on a 6.5 % term. **Buy engineering, not GPU hours.** At pop 400
folding falls to roughly 1.5 % of the generation, so the ceiling gets *worse*
as the population grows.

> v2 `random-r1` took **80.7 h** against 4.6 h for its siblings, purely from
> contention on 48 cores. With PSIQUE and both classifiers CPU-bound,
> concurrency is the first thing that will hurt you, not the last.

## What was fixed, and what each fix buys

| Fix | Commit | Effect |
|---|---|---|
| MACREL model loaded once per process, not per generation | `macrel in-process: fix against the real API…` | 71.5× on MACREL, but only **0.35 s/gen** — MACREL was never the bottleneck |
| `--norepeat` dedup indexed by sequence | `dedup: index the ancestral memory…` | 7.4 ms → 0.022 µs per test at 400 k rows (**330,000×**); removes a `pop² × gen²` wall |
| `preflight.sh` on every launcher | `Add the preflight guard…` | refuses to start if a classifier is missing — the v2 failure |

**The in-process MACREL path validates itself.** It calls library internals
rather than the documented CLI, and this project has twice shipped a classifier
that silently returned wrong numbers, so: it is checked against the subprocess
on magainin 2 before use, any disagreement beyond 1e-3 disables it permanently,
a mid-run exception reverts to the subprocess, `PFES_NO_INPROC=1` forces the
slow path, and `preflight.sh` prints which path the run will use. It cannot
quietly produce numbers the subprocess would not.

### Measured on the production node (HPC7)

The estimates above were replaced by direct measurement. One was badly wrong.

| | measured | my earlier estimate |
|---|---|---|
| `pypsique` | **7.8 ms / structure** | ~215 ms — **off by 27×** |
| HemoPI2 | **57.7 s / generation** | ~60 s ✓ |
| MACREL subprocess | **973 ms / generation** | 350 ms (laptop) |
| MACREL in-process | **15.2 ms / generation** | — (**64×**, saves ~1 s/gen) |

`preflight.sh` on the production node confirms `macrel path: in-process (fast)`
after the magainin-2 probe, so the fast path is live and validated there.

### Cost is mostly FIXED, which argues for a large population

Two points from the same node — pop 8 at ~62 s/generation, pop 100 at 88.6 s —
decompose as **~60 s fixed plus ~0.29 s per candidate**. HemoPI2 accounts for
~97 % of the fixed part. Treat this as indicative rather than settled: the two
runs used different `--threads`, and it needs confirming with a proper sweep.

| pop | s / generation | ms per candidate |
|---|---|---|
| 8 | 62.0 | 7,750 |
| 100 | 88.6 | 886 |
| 400 | 175 (projected) | 438 |

**A larger population is genuinely cheaper per candidate**, because the fixed
cost amortises: pop 400 evaluates 4× the candidates of pop 100 for roughly 2×
the wall time. That is the real argument for raising it, not the GPU.

### The unexplained 0.22 s per candidate

Of the ~0.29 s per candidate, folding is 0.058 s, `pypsique` 0.008 s and
in-process MACREL 0.0001 s. **~0.22 s is unaccounted for** — at pop 400 that is
89 s of every generation, 25 h per 1000-generation run, and after HemoPI2 is
sampled it becomes the largest single term. Profile it before committing six
runs:

```bash
python -c "
import cProfile,pstats,sys; sys.path.insert(0,'.')
sys.argv='pfes.py --start file --start-file init/init_random.faa -ps 8 -ng 2 --norepeat --hemo-every 99 -o /tmp/prof'.split()
cProfile.run(open('pfes.py').read(),'/tmp/prof.out')
pstats.Stats('/tmp/prof.out').sort_stats('cumulative').print_stats(25)"
```

**There is no psique problem.** At 7.8 ms it is 0.78 s of a generation at
pop 100 — under 1 %. The ~215 ms figure came from attributing everything
unexplained to psique, which was never sound reasoning.

**HemoPI2 is 65 % of a generation** — 57.7 s of 88.6 s, or 16 h of a
1000-generation run, for a column that never enters the fitness.

`--hemo-every N` samples it instead. One caveat that matters: `hemo_prob` is
cached with the sequence, so a molecule first evaluated in an unmeasured
generation **keeps NaN for as long as it survives, including into the final
population** — and screening candidates by hemolysis is the point. Recover it
post-hoc: `score_posthoc.py` scores HemoPI2 for the audited candidates, and
with `--lineage` along the whole ancestral line, so the trajectory survives a
sparsely-sampled run.

## The series

Production arm only (`main`), three origins, two replicates — **6 runs**.
The control arm (`fitness-esm3`) follows once these look right.

```bash
./preflight.sh || exit 1        # must pass; every launcher sources it

for origin in random fragments orfs; do
  for rep in r1 r2; do
    python pfes.py --start file --start-file init/init_${origin}.faa \
      -ps 400 -ng 1000 -sm weak -b 20 --strong_selection_after_n_gen 800 \
      --norepeat --max-tokens-per-batch 4096 \
      -o results/v4/${origin}-${rep} 2>&1 | tee results/v4/${origin}-${rep}.console.log
  done
done
```

`-b 20` is deliberate — the value derived in thesis §2.4.8.3, **not** the 32
that the corrected score scale would imply. See the open question below.

### Cost

| Scenario | s/gen | per run | 6 runs serial |
|---|---|---|---|
| before the fixes | 192–232 | 53–64 h | **13.3–16.1 days** |
| with MACREL in-process | ~122 | ~34 h | **~8.5 days** |

The second row is **wrong** and is kept only to show what the assumption cost.
It assumed moving MACREL in-process roughly halves the classifier budget.
Measured with macrel installed: the MACREL subprocess is **350 ms**, so the
in-process path saves 0.35 s of an 88.6 s generation — 0.4 %. By subtraction
**HemoPI2 is ~60 s per generation, essentially the entire classifier budget**,
and it is untouched. One v3 run logged eleven 600-second HemoPI2 timeouts,
which is consistent.

**HemoPI2 is the fix.** Two routes, and the choice is a design decision rather
than an optimisation: cache its model in-process the way MACREL now is, or stop
scoring it every generation. Hemolysis is an *attribute* that never enters the
fitness, so the second is available on the same argument that puts AMPlify
post-hoc — except §2.4.7.3 defends per-generation logging for trajectory
claims, so it costs something real. `PFES_SKIP_HEMO=1` already exists and turns
an 88.6 s generation into roughly 28 s.

For reference the v3 series was 8 runs, pop 100 × 600 generations, **78.6 h
(3.27 days) wall** with runs overlapping in streams.

## Two things to decide before this runs

**β against the corrected score scale.** Restoring Eq. 5 lowers scores by
roughly 0.6×, and Boltzmann weights depend on β × score *spread*, so the
selection strength matching the old β = 20 is nearer β = 32. §2.4.8.3 derives
20 and the series above keeps it. Deliberate, but re-derive it or accept
measurably weaker pressure than v3 had.

**The contact term is now near-inert at this chain length.** Correctly
implemented it contributes a log₁₀ spread of 0.0072, below the helix penalty's
0.0093, against 0.0374 before — it was the third-strongest selection pressure
in the objective and is now the weakest. A single α-helix has no tertiary
contacts, so 1.0 is the *right* answer for this class; Eq. 5 was designed for
250-residue globular folds. Keeping it is defensible (it costs nothing and it
is the published definition); dropping it from the AMP objective and saying why
is arguably more honest. Decide deliberately rather than by inheritance.

## Reading the output

`progress.log` carries **`amp_src`** per row — `macrel` or `proxy`. MACREL is
defined for 10–100 residues and the surrogate is substituted silently outside
it. Check this column before interpreting `amp_prob` from any run whose
population approached 100 aa.

```bash
python visual_pfes.py -l results/v4/<run>/progress.log -s results/v4/<run>/structures \
       -o results/v4/<run>/analysis --notraj
python analysis/score_posthoc.py results/v4/*/ -o results/v4/comparison --lineage
```

`--lineage` scores the ancestral line with both classifiers and plots agreement
against generation. The line is short — 28 to 49 members on the v3 runs — so it
is one extra batched call per run. A MACREL–AMPlify gap that is present from
generation zero means the two models disagree about this sequence family; a gap
that **opens** over the trajectory is evidence of specification gaming, and the
script flags any run where it widens by more than 0.15.
