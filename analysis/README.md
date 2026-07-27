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
