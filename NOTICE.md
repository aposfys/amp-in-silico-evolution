# Third-party components

This repository is MIT licensed ([`LICENSE`](LICENSE)), copyright 2026
Apostolos Fysekidis. That covers the work here and nothing else.

`LICENSE` is deliberately the bare MIT text with no additions, so that GitHub
and automated licence scanners classify it correctly; everything qualifying it
lives in this file. The pipeline bundles one binary and depends on several
models and tools, each under its own terms, and a permissive licence on this
repository does not extend to any of them.

## Required, but not distributed here

| Component | Licence | Source |
|---|---|---|
| `psique` | MIT, © 2023 Francisco Adasme | https://github.com/fadasme/psique |

`psique` assigns secondary structure and is called once per predicted structure
by `pfes.py`. The compiled binary is **not redistributed in this repository** —
it is a Linux x86-64 executable, and shipping third-party binaries in a source
repository serves nobody. Build or download it from the source above and place
it at `bin/psique`. `setup_gpu.sh` expects it there and `preflight.sh` checks
for it.

> Adasme-Carreño F., Caballero J. & Ireta J. (2021). PSIQUE: Protein Secondary
> Structure Identification on the Basis of Quaternions and Electronic Structure
> Calculations. *J. Chem. Inf. Model.* **61**, 1789–1800.
> [doi:10.1021/acs.jcim.0c01343](https://doi.org/10.1021/acs.jcim.0c01343)

`bin/dssp.sh` is a wrapper written for this project; it calls `dssp`, which is
not bundled and must be installed separately if used.

**Removed:** `bin/iqtree2` and `bin/muscle5` were inherited from upstream, were
not referenced anywhere in this repository, and are distributed under
GPL-family licences that the Unlicense this repository previously carried had
no authority to override. They were deleted rather than re-licensed, and are
not recoverable from this repository's history. If a phylogenetics or alignment
step is ever added, obtain them from upstream with their own licences.

## Models the pipeline fetches at runtime

None of these are redistributed here. Each is downloaded by the user, under
whatever terms its provider sets at the time.

| Component | Role | Terms |
|---|---|---|
| **ESM3** `esm3-sm-open-v1` | folds every candidate | the `esm` package is MIT; the **weights are gated on HuggingFace** and require accepting EvolutionaryScale's terms before download. Read what the gate presents — that acceptance, not this repository, governs your use of the weights. |
| **MACREL** | the activity term in the fitness | installed from bioconda; see its own repository |
| **HemoPI2** | hemolysis attribute, never in the fitness | installed from PyPI; see its own distribution |
| **AMPlify** | post-hoc independent audit | needs a separate environment (Python 3.6, old TensorFlow); see its own repository |

The ESM3 gate is the one that matters practically. The pipeline cannot run
without accepting it, and its terms may be narrower than this repository's MIT
licence — if you intend any use beyond research, read the gate text rather than
inferring permission from this file.

## Upstream: PFES, and what the MIT licence here does and does not claim

This project is a fork of [PFES](https://github.com/sahakyanhk/pfes) by Harutyun
Sahakyan, Sofya Babajanyan, Yuri Wolf and Eugene Koonin, released by its authors
into the public domain under the Unlicense — still its licence at the time of
writing, with no additional conditions stated in its README.

> Sahakyan H., Babajanyan S. G., Wolf Y. I. & Koonin E. V. (2025). In silico
> evolution of globular protein folds from random sequences. *PNAS* **122**,
> e2509015122. [doi:10.1073/pnas.2509015122](https://doi.org/10.1073/pnas.2509015122)

The Unlicense places no conditions on derivative works, so this fork may carry
any licence, and carries MIT. **Three things follow, and the copyright line in
`LICENSE` should be read against all of them.**

**The upstream code is public domain and stays public domain.** Nothing here
withdraws it. Anyone may take PFES from its own repository under the Unlicense,
and the MIT licence on this repository has no bearing on that.

**The MIT copyright claim is over this project's contributions**, not over the
upstream work it builds on. Six files are forked from PFES and still carry
substantial upstream code — `pfes.py`, `score.py`, `evolution.py`,
`visual_pfes.py`, `psique.py`, `bin/dssp.sh`. Across the
core simulator files the commit history is 65 commits by Harutyun Sahakyan
against 62 by Apostolos Fysekidis. Those files are derivative works: the
modifications are this project's, the material they modify is nobody's.

**Everything else in the repository is original**, including the objective and
its parameterisation, `preflight.sh` and the guards, `analysis/` in its
entirety, the starting-population construction, and the documentation. Which
terms of the objective are inherited and which are not is recorded per-term in
[`OBJECTIVE.md`](OBJECTIVE.md); which files are forked and which are not is the
list above.

If you want the unencumbered version, take PFES from upstream. If you want this
fork's objective, audit, starting populations or analysis, MIT applies and
attribution is the only condition.

## Data

Starting populations in [`init/`](init/) derive from UniProt/Swiss-Prot
(CC BY 4.0) and from translated small ORFs; provenance for every sequence is in
`init/init_fragments.tsv` and `init/source_sorfs.tsv`. Natural AMP reference
sequences used in figures come from DBAASP via `amp_db.py --fetch`, which is not
vendored.

## A note on scope

Licence terms are stated here as they were read at the time of writing, with a
link to each source. They change. Anyone redistributing this work, or using it
beyond research, should verify each one against its upstream rather than relying
on this table.
