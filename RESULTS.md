# Results — series v4, the objective ablation

Six runs, three starting populations × two objectives, from identical seeds.
Completed 30 August 2026 on an RTX 5090. Design and reasoning in the
[README](README.md#the-experiment); every constant in [`OBJECTIVE.md`](OBJECTIVE.md).

## The runs were comparable

This is stated first because the previous series was not, and the failure was
invisible in its output ([`VOID-RUNS.md`](VOID-RUNS.md)). Both arms passed an
identical preflight immediately before launching:

```
hemopi2:     model 3, magainin 2 0.229 < melittin 0.764
probe:       magainin-2 -> amp_prob 0.9505
macrel path: in-process (fast)
init md5:    08f234d4 / 81a0990d / e254933d
```

Same classifiers, same versions, same starting sequences by checksum. The
control branch was **merged** from `main` rather than cherry-picked, so the only
difference between the arms is the `score = np.prod([...])` block itself.

| | s / run | s / generation |
|---|---|---|
| `main` | 2072, 2056, 2048 | 3.41 – 3.45 |
| `control-fold-only` | 2231, 2280, 2286 | 3.72 – 3.81 |

## The result

Final generation, mean over the surviving population of 100.

| arm | origin | length | `amp_prob` | score | net charge | `amp_src` |
|---|---|---|---|---|---|---|
| **main** | fragments | 26.0 | **0.990** | 0.3366 | +3.26 | macrel, no proxy |
| **main** | orfs | 26.0 | **0.980** | 0.3347 | +2.89 | macrel, no proxy |
| **main** | random | 26.0 | **0.987** | 0.3363 | +4.04 | macrel, no proxy |
| ctrl | fragments | 65.0 | **0.091** | 0.8330 | −4.42 | macrel, no proxy |
| ctrl | orfs | 65.0 | **0.114** | 0.8349 | −2.90 | macrel, no proxy |
| ctrl | random | 56.6 | **0.028** | 0.8278 | +3.68 | macrel, no proxy |

**The two `score` columns are not comparable and must never be plotted
together.** They are different functions: `main` multiplies seven terms, the
control two. A control score of 0.83 against 0.34 says nothing except that a
product of two numbers near 0.9 is larger than a product of seven.

## What it establishes

### The objective is necessary, and the effect is not small

`amp_prob` drives selection on `main` and is a held-out measurement on the
control. It reads **0.980–0.990 against 0.028–0.114** — an order of magnitude,
with no overlap between the arms and none between any pair of runs.

The direction matters as much as the size. The control does not merely fail to
find antimicrobial peptides; it moves **away** from them. Optimising ESM3 fold
confidence alone drives a population that starts at MACREL ≈ 0.2 down to 0.03,
while the same optimiser from the same seeds under the composite objective drives
it to 0.99. Whatever `pLDDT × pTM` rewards is close to the opposite of what
MACREL recognises.

This is the controlled contrast §4 of the thesis lists as missing, where the
structural terms were "justified by the audit of §2.4.7.1 and by exploratory
observation rather than by a controlled contrast".

### The failure mode, demonstrated rather than argued

`pLDDT × pTM` is a reasonable objective — close to what
[Sahakyan et al. (2025)](https://doi.org/10.1073/pnas.2509015122) published, and
what they published it *for* was evolving globular folds from random sequence.
Optimised alone here it produces the wrong molecule class:

- **Chains more than double**, 26 → 57–65 residues, because fold confidence
  rises with length and nothing bounds it.
- **Two of three origins go anionic**, to −2.90 and −4.42, against +2.89 to +4.04
  under the full objective. Antimicrobial peptides are cationic; this is the
  defining property, and the control abandons it.

Same optimiser, same operators, same starting sequences, same 600 generations.
Only the objective differs, and it decides whether the run returns a peptide or
a protein. That is the clearest available statement of the methodological
argument this thesis makes.

### The length penalty binds, and the correction did not break it

`P_L` holds the population at exactly 26.0 residues in all three `main` runs,
against 56.6–65.0 with the term absent. The natural comparison: anuran AMPs
average 24 residues and APD6's synthetic entries 19, so the constraint lands the
search inside the class it was retuned for ([`OBJECTIVE.md`](OBJECTIVE.md)).

### Net charge is lower than the previous series, and closer to nature

The `main` arm reaches **+2.89 to +4.04**, against +6.05 to +7.29 reported for
the old condition A. That difference is consistent with the confound documented
in [`VOID-RUNS.md`](VOID-RUNS.md): the four arms driving those figures ran the
biophysical `calculate_samp` surrogate, which is an explicit function of net
charge, so optimising it drove charge up. MACREL's random forest saturates
without needing +7. The new figures sit much nearer the +2.5 mean of the anuran
set and well inside the +2 to +9 reported for the class
([Zhang et al. 2021](https://doi.org/10.1186/s40779-021-00343-2)).

### The origin effect does not reproduce

Under the corrected objective and a single scorer, the three origins converge:

| | score spread across origins | `amp_prob` spread |
|---|---|---|
| `main` | **0.0019** (0.3347 – 0.3366) | 0.010 |
| previous condition A | 0.05 – 0.07, reported as an origin effect | — |

A spread of 0.0019 is smaller than the within-replicate spread of the previous
series (0.002–0.022). Three origins that were said to reach different optima
reach the same one, to three decimal places, once they are scored by the same
classifier. This is what the confound predicted, and it is the reason §3.2.3 and
§4.3.1 of the thesis need revising.

The control arm is more dispersed — length 56.6 to 65.0, charge −4.42 to +3.68,
`amp_prob` 0.028 to 0.114 — which is consistent with the mechanism §4.3.1
proposed: with no length constraint the population has room to move, and where
it lands is less determined. But see the limits below before reading anything
into it.

## Limits

**One replicate per cell.** Six runs is enough for the ablation, whose effects
are an order of magnitude larger than any plausible run-to-run noise. It is
**not** enough for any statement about differences *between origins within an
arm*: those differences are 0.0019 in score on `main`, and there is no estimate
of the within-cell variance to compare them against. The control arm's larger
dispersion is likewise not yet distinguishable from noise. Replicating the
control is the next run, and it is cheap: an arm is under two hours.

**Nothing was measured in a laboratory.** Every quantity here is a model output.
A MACREL probability of 0.99 may reflect proximity to the classifier's training
distribution as readily as biological activity, and these sequences are novel by
construction.

**The held-out classifiers are not yet reported.** AMPlify needs its own Python
3.6 environment and was not installed for these runs; HemoPI2 is present and
sampled every 25 generations. The independent audit in
[`analysis/README.md`](analysis/README.md) is what turns MACREL 0.99 from a
restatement of the objective into evidence, and it is still outstanding.

**HemoPI2 model 3 puts poly-alanine at 0.638**, above its own 0.55 threshold —
a false positive on an inert homopolymer. Predicted hemolysis on non-natural
evolved peptides carries more uncertainty than the published AUROC of 0.921
suggests, and these peptides are non-natural by construction.

## Reproducing

```bash
./preflight.sh                       # must pass before anything is believed
MAXTOK=4096 ./run_v4.sh main         # 3 runs, ~1h45m total

git worktree add ../amps-ctrl control-fold-only
cd ../amps-ctrl && MAXTOK=8192 ./run_v4.sh ctrl
```

Runs are not reproducible individually — neither this fork nor upstream PFES
seeds any RNG, so a trajectory cannot be regenerated and must be preserved from
its `progress.log`. What is fixed is the starting population, by file and by
checksum. The evidence this design admits is replication across independent
runs, not reproduction of one.
