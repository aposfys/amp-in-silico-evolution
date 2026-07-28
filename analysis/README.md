# Figure and report scripts

Every script that produces a graph in this project, and where it came from.

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
python analysis/compare_runs.py results/final/*/ -o results/comparison
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
  best of the final generation. The structured arm peaks at generation 241 and
  drifts down. The lineage is traced from the final generation's best, matching
  `lineage.tsv`; the global best lies on that line but is not its endpoint.

### `score_posthoc.py` — the independent audit

MACREL drives selection, so a high MACREL score on an evolved peptide restates
the objective rather than corroborating it. AMPlify is architecturally unrelated
and never sees the search, so its verdict is independent evidence.

```bash
export AMPLIFY_CMD="conda run -n amplify AMPlify"   # own env: py3.6, old TF
python analysis/score_posthoc.py results/final/*/ -o results/comparison
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
| `make_ppt.py` | arm-comparison figures and a `.pptx` deck (matplotlib, python-pptx) |
| `make_word.py` | a `.docx` analysis report from one or more run logs (python-docx) |

Run them from the repository root so relative paths to `results/` resolve:

```bash
python analysis/make_loop_diagram.py
```

Their output (`*.docx`, `*.pptx`, root-level `*.png`) is gitignored, being
regenerable from the logs.

## Note on reproducibility

Plots regenerate exactly from a `progress.log`, because the log is the record.
The runs themselves do not: neither upstream PFES nor this fork seeds any random
number generator, so a run cannot be repeated. Keep `progress.log` and
`structures/` for anything you intend to cite.
