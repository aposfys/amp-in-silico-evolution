# Third-party components

This repository is MIT licensed ([`LICENSE`](LICENSE)). That covers the work
here and nothing else. The pipeline bundles one binary and depends on several
models and tools, each under its own terms, and a permissive licence on this
repository does not extend to any of them.

## Bundled in this repository

| Component | Licence | Notice |
|---|---|---|
| `bin/psique` | MIT, © 2023 Francisco Adasme | [`bin/LICENSE.psique`](bin/LICENSE.psique) |

`psique` assigns secondary structure and is called once per predicted structure
by `pfes.py`. MIT requires the copyright and permission notice to accompany
every copy, so it is kept beside the binary rather than only referenced here.

> Adasme-Carreño F., Caballero J. & Ireta J. (2021). PSIQUE: Protein Secondary
> Structure Identification on the Basis of Quaternions and Electronic Structure
> Calculations. *J. Chem. Inf. Model.* **61**, 1789–1800.
> [doi:10.1021/acs.jcim.0c01343](https://doi.org/10.1021/acs.jcim.0c01343)

`bin/dssp.sh` is a wrapper written for this project; it calls `dssp`, which is
not bundled and must be installed separately if used.

**Removed:** `bin/iqtree2` and `bin/muscle5` were inherited from upstream, were
not referenced anywhere in this repository, and are distributed under
GPL-family licences that the Unlicense this repository previously carried had
no authority to override. They were deleted rather than re-licensed. Recover
them from git history if a phylogenetics or alignment step is ever added —
`git show 1d4f2a8:bin/iqtree2` — but obtain them from upstream with their own
licences instead.

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

## Upstream

This project is a fork of [PFES](https://github.com/sahakyanhk/pfes) by Harutyun
Sahakyan, Sofya Babajanyan, Yuri Wolf and Eugene Koonin, released into the
public domain under the Unlicense. The Unlicense imposes no conditions on
derivative works, so this fork is MIT licensed; the upstream code remains public
domain regardless.

> Sahakyan H., Babajanyan S. G., Wolf Y. I. & Koonin E. V. (2025). In silico
> evolution of globular protein folds from random sequences. *PNAS* **122**,
> e2509015122. [doi:10.1073/pnas.2509015122](https://doi.org/10.1073/pnas.2509015122)

The objective, its parameterisation, the audit against the published equations,
the starting-population construction and the analysis in this repository are
original work. What is inherited and what is not is recorded per-term in
[`OBJECTIVE.md`](OBJECTIVE.md).

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
