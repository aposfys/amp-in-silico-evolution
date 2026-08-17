# PFES-AMPs

In silico evolution of antimicrobial peptides. **ESM3** folds each candidate, **MACREL**
scores it for antimicrobial activity, **HemoPI2** predicts hemolysis, and the **PFES**
structural penalties constrain length and secondary structure. A population of one
hundred peptides is mutated and selected for six hundred generations.

Built for the MSc thesis *Investigating novel antimicrobial peptides using machine
learning, structure prediction and in silico evolution* (Apostolos Fysekidis,
Bioinformatics, National and Kapodistrian University of Athens), on top of
[PFES](https://github.com/sahakyanhk/pfes) by Sahakyan et al.

## Where to look

| Branch | Objective | |
|---|---|---|
| **`fitness-pfes-macrel`** | pLDDT × pTM × length × helix × β × **MACREL** × contacts | the full objective — **the production branch** |
| `fitness-fold-macrel` | pLDDT × pTM × **MACREL** | activity and fold only, no structural penalties |
| `fitness-pfes` | pLDDT × pTM × length × helix × β × contacts | the upstream objective, no activity term |
| `upstream` | — | PFES as published, **unmodified**; the reference the parameterisation audit is against (was `alpha`) |
| `main` | — | this page |

[`production/`](production/) holds the twelve-run design study the thesis reports:
three starting populations × two penalty conditions × two replicates, one compressed
archive per run, plus the starting populations, the launch scripts and the winners.
Its [README](production/README.md) describes the design and the result.

## The objective

```
score = pLDDT × pTM × length_penalty × helix_penalty × beta_penalty × AMP_probability × contact_density
```

| Term | Source | Range | Note |
|---|---|---|---|
| pLDDT | ESM3 | 0–1 | per-residue fold confidence |
| pTM | ESM3 | 0–1 | overall-fold confidence |
| AMP_probability | MACREL | 0–1 | probability the peptide is antimicrobial |
| length_penalty | PFES | 0–1 | logistic, → 0 above the `-pl0` target |
| helix_penalty | PFES (PSIQUE) | 0–1 | → 0 for α-helices longer than `-hl0` |
| beta_penalty | PFES (PSIQUE) | 0–1 | → 0 for β-strands longer than `-bl0` |
| contact_density | PFES | ≥ 1 | rewards compact folds |
| (1 − hemolysis) | HemoPI2 | 0–1 | **not in the score** unless `--hemo-in-score` |

> **Hemolysis is an attribute by default.** HemoPI2 runs every generation and writes
> `hemo_prob` for every peptide, but it does not affect selection. Screen candidates on
> it afterwards, or add `--hemo-in-score` to fold `× (1 − hemolysis)` into the fitness.

## Install

```bash
git clone https://github.com/aposfys/PFES-AMPs.git && cd PFES-AMPs
git checkout fitness-pfes-macrel

conda create -n pfes_amps python=3.11 -y && conda activate pfes_amps
pip install torch --index-url https://download.pytorch.org/whl/cpu   # drop --index-url for GPU
pip install -r requirements.txt
pip install "numpy<2" "pillow<12"          # HemoPI2/scikit-learn need numpy < 2

conda install -c bioconda -c conda-forge macrel -y

# psique (secondary structure) is bundled — bin/psique + psique.py, nothing to install
# ESM3 needs a token; accept the licence at
#   huggingface.co/EvolutionaryScale/esm3-sm-open-v1
export HF_TOKEN=your_token_here
```

Verify:

```bash
python -c "from esm.models.esm3 import ESM3; print('ESM3 ok')"
python -c "import score; print(score.macrel_score_batch(['GIGKFLHSAKKFGKAFVGEIMNS']))"
```

## Run

**1. Seed database and hemolysis labels** (once):

```bash
export PFES_AMP_DB=$PWD/data/dbaasp.faa
python amp_db.py --fetch          # ~10,700 validated AMPs from DRAMP
python amp_db.py --annotate-hemo  # a HemoPI2 label per database peptide
```

**2. A starting population.** The production series does not use `amp_db.py` for this —
its three sets are built by the scripts in `production/` and shipped in
`production/init/` and `production/init_varlen/`. For a quick run:

```bash
python amp_db.py --make-init --pop 100 --frac-existing 0.5 --seed 42 -o init_mix.faa
```

**3. Smoke test** before committing to a full run:

```bash
python pfes.py --start file --start-file init_mix.faa -ps 8 -ng 3 -sm weak -b 20 -o /tmp/smoke
```

No "falling back" warnings, and a `hemo_prob` column that varies, means MACREL and
HemoPI2 both engaged.

**4. The full run**, as the production series was launched:

```bash
python pfes.py \
  --start file --start-file production/init/init_fragments.faa \
  -ps 100 -ng 600 -sm weak -b 20 \
  -pl0 30 -hl0 30 -bl0 12 \
  --norepeat --max-tokens-per-batch 4096 \
  -o results/fragments-r1
```

Run it inside `tmux`, so a disconnect cannot kill it.

**5. Analyse:**

```bash
python visual_pfes.py -l results/fragments-r1/progress.log \
  -s results/fragments-r1/structures -o results/fragments-r1/analysis --notraj
```

## Options that matter

| Flag | Default | Meaning |
|---|---|---|
| `--start {random,randoms,existing,mix,file}` | random | starting population; `file` loads `--start-file` |
| `--start-file PATH` | — | a fixed init FASTA — how every controlled comparison is set up |
| `-ng` / `-ps` | 100 / 10 | generations / population size |
| `-sm` / `-b` | weak / 1 | selection: `weak` (Gibbs, β via `-b`), `weak2`, `strong` |
| `-pl0` / `-hl0` / `-bl0` | 30 / 20 / 12 | length / helix / β-strand penalty thresholds |
| `--hemo-in-score` | off | put hemolysis back into the fitness |
| `--norepeat` | off | never re-evaluate a sequence already seen |
| `--max-tokens-per-batch` | 512 | fold batch size; ~4096 minimises classifier reloads |

## Output

`<outdir>/progress.log` is tab-separated with a `#`-prefixed banner recording the exact
invocation. One row per surviving individual per generation, nineteen columns:

`gndx` `id` `seq_len` · `prot_len_penalty` `max_alpha_penalty` `max_beta_penalty` ·
`ptm` `mean_plddt` `num_conts` `iplddt` `num_inter_conts` · `sel_mode` ·
`amp_prob` (MACREL) `hemo_prob` (HemoPI2) · `score` `sequence` `mutation` `prev_id` `ss`

Note that `gndx` is a **string** (`gndx0`, `gndx1`, …), not an integer. Structures are
written as gzipped PDB to `<outdir>/structures/` — sixty thousand of them for a
100 × 600 run, about 244 MB.

## Reproducibility

Runs are **not** reproducible, matching upstream PFES, which seeds no random number
generator anywhere. Mutation, survivor sampling and random-sequence generation all draw
from OS entropy, so the same command produces a different trajectory every time. Repeat
runs are independent samples of one process, which is what makes replication across runs
meaningful — but an individual result cannot be regenerated and has to be preserved from
its log.

The starting population is the one thing that *is* fixed, and it is fixed by a file
rather than a seed: build it once, then pass the same `--start-file` to every arm.

## Known constraints

**MACREL is defined for 10–100 residues.** Outside that window `macrel_score_batch`
falls back to the `calculate_samp` biophysical proxy, silently and per sequence, so the
objective can change identity mid-run if chains grow past 100 residues. This is not
hypothetical: the fold-only arm reached 93 residues by generation 300. Either cap length
or log which scorer produced each value.

**A length target far above where the population settles is inert.** The penalty is a
logistic in `-pl0` with steepness hard-coded at `pfes.py:237`; set the target high enough
and it saturates at 0.999 and stops selecting on anything. Condition B of the production
series is exactly this case, deliberately — see `production/README.md`.

**Deduplication cost grows quadratically.** `--norepeat` scans the whole
`ancestral_memory` for every candidate, and that table gains `-ps` rows per generation,
so runtime is not linear in `-ng`.

## Hardware

Upstream PFES was developed and tested on Rocky Linux 8.7 with NVIDIA V100 and A100
GPUs. **The twelve production runs here ran on CPU** — every one logs `ready [cpu]` at
startup — so the fork is exercised on a path upstream did not report.

This does not change what is computed. Device selection at `pfes.py:930` only chooses
where ESM3 runs; MACREL, HemoPI2 and PSIQUE are CPU tools on either platform. The one
place the code branches on device is `use_threads = device.type != 'cpu'`, which
overlaps the scoring thread of one batch with the folding of the next: a pipelining
optimisation on GPU, serial on CPU, identical arithmetic either way.

What it does change is throughput, and that shaped the protocol — six hundred
generations rather than the several thousand of the original, with selection sharpened
to compensate. Anyone reproducing on a GPU should expect the same distribution of
outcomes but not the same numbers, since ESM3 matmuls run in reduced precision under
TF32 on Ampere while CPU stays at full fp32. That is moot in practice: nothing here is
bit-reproducible on any hardware, because no random number generator is seeded.

A third path exists and is a fork addition: Apple MPS is checked *before* CUDA, and a
monkeypatch at the top of `pfes.py` disables ESM3's fp32 autocast there because MPS does
not implement it. That path is untested at scale and was not used for any reported run.

## Licence

Public domain (Unlicense), as upstream. See [LICENSE](LICENSE).
