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
| MACREL model loaded once per process, not per generation | `macrel: load the model once…` | the 68 % term; self-validating, see below |
| `--norepeat` dedup indexed by sequence | `dedup: index the ancestral memory…` | 7.4 ms → 0.022 µs per test at 400 k rows (**330,000×**); removes a `pop² × gen²` wall |
| `preflight.sh` on every launcher | `Add the preflight guard…` | refuses to start if a classifier is missing — the v2 failure |

**The in-process MACREL path validates itself.** It calls library internals
rather than the documented CLI, and this project has twice shipped a classifier
that silently returned wrong numbers, so: it is checked against the subprocess
on magainin 2 before use, any disagreement beyond 1e-3 disables it permanently,
a mid-run exception reverts to the subprocess, `PFES_NO_INPROC=1` forces the
slow path, and `preflight.sh` prints which path the run will use. It cannot
quietly produce numbers the subprocess would not.

Still un-fixed and worth doing next: **PSIQUE spawns one subprocess per
structure.** At pop 400 that is 400 spawns per generation and it is the largest
remaining term after the classifiers.

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
