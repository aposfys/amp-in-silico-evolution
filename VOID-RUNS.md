# Every run before v4 is void

There is no results directory in this repository, and this file is why. All
previous production output has been removed. The record of what was run, and of
what was wrong with it, is kept here because the failures are the reusable part.

The deleted archives remain in git history at commit `0942520` and can be
recovered with `git show 0942520:production/runs/<name>.tar.gz > <name>.tar.gz`,
so every measurement below is still checkable.

## What was deleted

| | | |
|---|---|---|
| `pfes-results/` | 3.1 GB | raw HPC output — v1, v2, v3, and the `final` set |
| `production/runs/` | 10 MB | twelve archives, one per run of the v2+v3 design study |
| `production/winners.tsv` | — | the twelve winners |
| `results_macrel_pfes_300/` | 13 MB | the v1 300-generation run and its figures |
| `runs/` | 1.1 MB | older per-origin summaries |

The starting populations survive at [`init/`](init/) and
[`init_varlen/`](init_varlen/), and the launch scripts at the repository root.
`production/init`, `production/init_varlen` and `production/scripts` were
byte-identical duplicates of those, so nothing was lost with them.

## Why each series is void

### v1 — a different fitness function

`(num_CA_contacts_within_8Å + L) / L` as the contact term, and the hemolysis
proxy multiplied into the score. Neither is the objective this project runs. The
`ANALYSIS.md` that described it has been superseded by
[`OBJECTIVE.md`](OBJECTIVE.md).

### v2 — the classifiers were not installed

Four runs, about 60 hours, completed with neither MACREL nor HemoPI2 present.
`score.py` fell back per sequence to the biophysical `calculate_samp` and
`calculate_hemo_proxy` surrogates and warned only on stderr, while each run's
header kept printing `[MACREL AMP + PFES]`. The cause was `conda activate
pfes_amps` resolving to whichever environment of that name came first on the
HPC — there were two, and only one had macrel.

Each of the four console logs carried **2,401 fallback warnings**:

```
Warning: macrel not installed — conda install -c bioconda macrel
           falling back to biophysical hemo proxy
```

`preflight.sh` exists because of this, and every launcher sources it. It probes
magainin-2 end to end and refuses to start if the surrogate answers.

### v3 — sound, but on a fitness function that has since been corrected

Eight runs with both classifiers live. Three defects were found afterwards by
auditing the code against Sahakyan et al. 2025 and its SI Appendix
(`0d06bb8`, `9c0ba93`):

- `get_nconts` counted the α-helical contacts Eq. 5 excludes. Measured on this
  project's own output at 8 Å, **95–98 % of counted contacts were |i − j| = 4** —
  the term was scoring helical turns, not tertiary compactness.
- `cbiplddt` used the wrong interface definition, a 100× error on ipLDDT.
- `uniprotrates` was scaled 20× wrong.

The empirical output stands as a record of what those runs did. It is not
comparable to anything run after 2026-08-23, because the objective changed.

### The condition-A origin comparison was also confounded

Worth separating from the above, because it is a design failure rather than a
code one, and it is the reason the v4 series exists.

`production/runs/` held twelve archives presented as three origins × two penalty
conditions × two replicates. Four of the six **condition-A** cells were v2
archives — both `random` replicates and both `fragments` replicates — while the
`sORF` cell was v3. `results/v3/` never contained a condition-A run for random
or fragments; they were never re-run after the classifiers were installed, and
the v2 archives were used to fill the cells.

| cell | series | gen-0 `amp_prob` max | > 0.5 of 100 | median `hemo_prob` |
|---|---|---|---|---|
| `random-r1` | **v2** | 0.980 | **32** | 0.879 |
| `random-r2` | **v2** | 0.928 | **27** | 0.868 |
| `fragments-r1` | **v2** | 0.800 | **24** | 0.870 |
| `fragments-r2` | **v2** | 0.853 | **25** | 0.903 |
| `sorf-r1` | v3 | 0.594 | 4 | 0.445 |
| `sorf-r2` | v3 | 0.634 | 3 | 0.415 |

The starting populations are MACREL-screened to **0/100 above 0.5, maximum
0.495**, so no screened set can yield 24–32 of 100 above 0.5 after one mutation
round. That is the second, independent confirmation of the fallback.

**Origin was therefore confounded with scorer**, which accounts for every
separation the comparison reported: the score gap (two different fitness
functions), hemolysis 0.36–0.49 against 0.74–0.76 (`calculate_hemo_proxy`
against HemoPI2 — two different predictors, not one measurement), and net charge
+6.1–7.3 against +3.2–3.9 (`calculate_samp` is an explicit function of charge;
MACREL's random forest saturates without needing +7).

Condition B was internally consistent — all six `long-*` runs were v3 — and it
showed **no** origin effect. So the one clean condition found nothing, and the
condition that found something could not separate origin from scorer. Thesis
§3.2.3 and §4.3.1 are written against the earlier reading and need revising.

## The pattern

Three failures, one shape: a component adequate for the job it was built for
became load-bearing in a job it was not built for, and nothing in the pipeline
announced the substitution. A surrogate stood in for a classifier. A contact
definition written for 250-residue globular folds was applied to 25-mers. An
archive from one series filled a cell in another.

The same shape is still present in the current objective and is documented
rather than fixed: pTM is the strongest term in the fitness and is evaluated 14
residues below the floor of its own derivation, and β = 20 is strong selection
by the source paper's criterion while the flag that sets it is named `weak`. See
[`OBJECTIVE.md`](OBJECTIVE.md).

Where a learned score, an inherited metric or a published constant becomes part
of an objective rather than part of a report, it needs what any other dependency
needs: a version pin, a test against a case whose answer is known by
construction, and a record in version control of the value that actually ran.

## What replaces them

The v4 series — three origins × two objectives, six runs, described under
**The experiment** in the [README](README.md#the-experiment) and launched by
[`run_v4.sh`](run_v4.sh). It asserts its own branch and refuses to start on a
modified `pfes.py` or `score.py`, and it greps each console log for a surrogate
fallback at the end of every run rather than at the end of the analysis.
