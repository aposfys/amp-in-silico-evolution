# Figure and report scripts

Every script that produces a graph in this project, and where it came from.

The objective these scripts describe is the one in [`OBJECTIVE.md`](../OBJECTIVE.md),
which traces every constant to its source and measures the pressure each term
applies. Two of its findings change how output should be read, so they are
repeated here rather than left to be discovered: **pTM is the strongest term in
the fitness for most of a run** and is evaluated below the length floor of its
own derivation, and **β = 20 is strong selection** by the source paper's own
criterion, not the weak selection the flag names.

> `make_word.py` and `make_ppt.py` were removed. They were v1-era report
> generators hard-coded to `results_macrel_300/` and `results_macrel_pfes_300/`,
> both since deleted, and they documented the fitness as containing
> `(1 − hemo_proxy)`. Hemolysis has not been in the score since; it is an
> attribute, and putting a safety term in the objective is the thing this design
> deliberately does not do. Recover them from git history if the layout is ever
> useful: `git show fba3369:analysis/make_word.py`.

## Inherited from upstream PFES

Both live at the repository root, at their original paths, so the upstream
invocation still works unchanged.

### `visual_pfes.py` — the main analysis

Reads a run's `progress.log` and `structures/`, reconstructs the winning
lineage, and writes the plot set. Extended in this fork; the upstream functions
are untouched.

```bash
python visual_pfes.py -l results/<run>/progress.log \
                      -s results/<run>/structures \
                      -o results/<run>/analysis --notraj
```

| Output | Origin |
|---|---|
| `Summary.png` | upstream |
| `Secondary_structures.png` | upstream |
| `plots/<column>.png` (one per logged metric) | upstream |
| `Evolution.png` | this fork |
| `Score_components.png` | this fork |
| `Score_distribution.png` | this fork |
| `Fitness_landscape.png` | this fork |
| `AA_composition.png` | this fork |
| `lineage.tsv`, `bestlog.tsv` | this fork |

`plots/amp_prob.png` and `plots/hemo_prob.png` exist only because this fork adds
those columns to the log.

Drop `--notraj` to also write the superposed backbone trajectory (needs
MDAnalysis). That trajectory is the input to the movie script below.

### `pymol_vstraj.py` — trajectory movie

Renders every *n*th structure with PyMOL and concatenates the frames into a
video, which is how the upstream paper produced its supplementary movies.
Requires `pymol` and `moviepy`, neither of which is in `requirements.txt`.

```bash
python pymol_vstraj.py results/<run>/structures
```

## Written for this thesis

These are not part of the simulator and nothing in the run path imports them.

### `compare_runs.py` — across runs

`visual_pfes.py` describes one run and knows nothing about any other, so every
cross-run statement otherwise has to be assembled by hand. This takes any number
of run directories, groups them by objective, and plots each metric as a mean
with an across-repeat band.

```bash
python analysis/compare_runs.py results/v4/*/*/ -o results/v4/comparison
```

| Output | Shows |
|---|---|
| `cmp_trajectories` | chain length and fitness per generation, band across repeats |
| `cmp_inputs` | pLDDT, pTM, contacts, MACREL probability |
| `cmp_helix` | longest α-helix, against the `-hl0` threshold |
| `cmp_hemolysis` | HemoPI2 over generations (skipped if `PFES_SKIP_HEMO` was set) |
| `cmp_variance` | endpoint metrics as individual runs plus mean ± SD |
| `per_run.tsv`, `summary.tsv` | per-run rows and per-arm mean and SD |

Works with a single run per arm, drawing no band and labelling the legend `n=1`.

Arms are identified **from the data**: both objectives are reconstructed from
the logged columns and whichever reproduces `score` wins, so a renamed or
mislabelled directory is reported rather than trusted.

Two counting conventions, because both are easy to get wrong:

- `lineage_nodes` counts unique members of the winning ancestral line.
  `visual_pfes.extract_lineage` writes its starting node twice, so
  `lineage.tsv` has one more row than the line has members.
- `winner_*` describes the **globally** best individual, which is not always the
  best of the final generation, because a run can peak and drift down. The
  lineage is traced from the final generation's best, matching `lineage.tsv`;
  the global best lies on that line but is not necessarily its endpoint.

### `score_posthoc.py` — the independent audit

MACREL drives selection, so a high MACREL score on an evolved peptide restates
the objective rather than corroborating it. AMPlify is architecturally unrelated
and never sees the search, so its verdict is independent evidence.

```bash
export AMPLIFY_CMD="conda run -n amplify AMPlify"   # own env: py3.6, old TF
python analysis/score_posthoc.py results/v4/*/*/ -o results/v4/comparison
```

| Output | Shows |
|---|---|
| `posthoc_classifiers` | MACREL against AMPlify per candidate. Points below the diagonal are peptides the objective rewarded that an independent model does not recognise. |
| `posthoc_charge_length` | length against net charge with natural AMPs as background |
| `posthoc_scores.tsv` | length, charge, hydrophobic fraction, helix, MACREL, AMPlify, HemoPI2 |

Takes the global best plus the top `-n` of the final generation from each run.
Net charge counts K and R as +1, D and E as −1, and H as +0.1 for its protonated
fraction at pH 7.4, reproducing the +5.2 reported for the structured winner.

Natural references come from `data/dbaasp.faa` if `amp_db.py --fetch` has been
run, otherwise from a small curated set, and the figure states which. Without
AMPlify installed it still runs, reports MACREL and HemoPI2, and says which
figure it could not draw.

### The rest

| Script | Produces |
|---|---|
| `make_loop_diagram.py` | `loop_diagram.png`, the evolutionary-loop schematic (matplotlib) |
| `make_init_sets.py` | the three starting populations; see [`init/README.md`](../init/README.md) |
| `build_sorf_source.py` | `source_sorfs.faa`, the sampling frame for the sORF arm |

Run them from the repository root so relative paths resolve:

```bash
python analysis/make_loop_diagram.py
```

Root-level `*.png` output is gitignored, being regenerable from the logs.

## What the analysis assumes, and where those choices come from

Conventions that appear in more than one script, stated once so a number is not
computed two ways in two places.

**Net charge** counts K and R as +1, D and E as −1, and H as +0.1 for its
approximate protonated fraction at pH 7.4 (`score_posthoc.py:58`). That
reproduces the +5.2 reported for the structured winner. The evolved populations
reach +6 to +7, inside but at the top of the +2 to +9 range reported for the
class ([Zhang et al. 2021](https://doi.org/10.1186/s40779-021-00343-2)) and well
above the +2.5 mean of the anuran set.

**Hydrophobic fraction** is `AVLIMFWYC` in `score_posthoc.py:66` — including
tyrosine. [`init/README.md`](../init/README.md) reports starting compositions on
`AVLIMFWC`, *excluding* it, so its own table and this one are not directly
comparable; the matched recomputation is given there. Check which definition a
figure used before comparing a starting value with a final one.

**AMPlify is the auditor, never the objective.** MACREL drives selection, so a
high MACREL score on an evolved peptide restates the objective rather than
corroborating it; AMPlify shares neither MACREL's input representation, its
inductive bias, nor its training set
([Li et al. 2022](https://doi.org/10.1186/s12864-022-08310-4)), so its errors are
usefully independent. The asymmetry is permanent, not provisional: AMPlify
saturates at exactly 1.0000 for a large fraction of AMP-like candidates, and a
term pinned at its ceiling supplies no gradient for selection to climb.
Generating under classifier guidance and screening with filters the generator
never saw follows
[Das et al. (Nat Biomed Eng 2021)](https://www.nature.com/articles/s41551-021-00689-x);
the reason it is necessary is that optimising a learned predictor drives the
search to where the predictor is unreliable
([Brookes, Park & Listgarten, ICML 2019](https://proceedings.mlr.press/v97/brookes19a.html)).

**Hemolysis is read, never optimised.** HemoPI2
([Chaudhary et al. 2016](https://doi.org/10.1038/srep22843)) is logged and never
enters the fitness, which is what lets it report on where a search went rather
than on how well the search satisfied it. The quantity it speaks to is the
therapeutic index, the ratio of minimum haemolytic to minimum inhibitory
concentration ([Cardoso et al. 2021](https://doi.org/10.1007/s12551-021-00784-y)).

**`--hemo-every N` leaves gaps that must be recovered.** `hemo_prob` is cached
with the sequence, so a molecule first evaluated in an unmeasured generation
keeps NaN for as long as it survives, including into the final population.
`score_posthoc.py --lineage` scores the ancestral line post hoc and closes them.

**`amp_src` before `amp_prob`.** MACREL is defined for 10–100 residues
([Santos-Júnior et al. 2020](https://peerj.com/articles/10555/)); outside that
window the biophysical surrogate is substituted per sequence and the column
changes identity. Filter on `amp_src == 'macrel'` — this matters most on the
`fitness-esm3` arm, which has no length penalty and whose chains grow past 100.

## Note on reproducibility

Plots regenerate exactly from a `progress.log`, because the log is the record.
The runs themselves do not: neither upstream PFES nor this fork seeds any random
number generator, so a run cannot be repeated. Keep `progress.log` and
`structures/` for anything you intend to cite.
