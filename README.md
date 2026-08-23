# PFES-AMPs — `fitness-pfes-macrel`

**The production branch.** In silico evolution of antimicrobial peptides:
**ESM3** folds each candidate, **MACREL** scores it for antimicrobial activity,
the **PFES** structural terms constrain length and secondary structure,
**HemoPI2** records hemolysis as an attribute, and **AMPlify** audits the
survivors afterwards.

```
score = pLDDT × pTM × length_penalty × helix_penalty × beta_penalty
        × AMP_probability × contact_density
```

| Term | Source | Range | Note |
|---|---|---|---|
| pLDDT | ESM3 | 0–1 | mean per-residue fold confidence |
| pTM | ESM3 | 0–1 | overall-fold confidence |
| length_penalty | PFES | 0–1 | →0 above `-pl0` (default 30) |
| helix_penalty | PFES / PSIQUE | 0–1 | →0 for α-helices > `-hl0` |
| beta_penalty | PFES / PSIQUE | 0–1 | →0 for β-strands > `-bl0` |
| AMP_probability | MACREL | 0–1 | probability the peptide is antimicrobial |
| contact_density | PFES Eq. 5 | ≥ 1 | tertiary compactness — see below |
| (1 − hemolysis) | HemoPI2 | 0–1 | **OFF by default**; `--hemo-in-score` to enable |

Paired against [`fitness-esm3`](../../tree/fitness-esm3), which is pLDDT × pTM
alone, this isolates what the structural terms and MACREL together contribute.
Run both from the same `--start-file` and the difference is attributable to the
objective.

## Contact density measures compactness, not helicity

`contact_density = (n_contacts + l) / l`, with contacts defined exactly as in
Eq. 5 of Sahakyan et al. 2025: **Cβ atoms, within 6 Å, more than 5 residues
apart in sequence**, both above the pLDDT floor.

The `|i − j| > 5` rule is the whole point — it excludes the i→i+4 α-helical
register, so a long helix does not score as a compact fold. Upstream PFES loops
`range(i + 4, …)`, which admits |i − j| = 4 and 5 and therefore counts exactly
those helical contacts. At upstream's 6 Å the geometry mostly hides it; at a
wider cutoff it dominates. Measured on this project's own output at 8 Å,
**95–98 % of counted contacts were |i − j| = 4** — the term was counting helical
turns. Fixed in `0d06bb8`.

A single α-helix has no tertiary contacts, so the corrected term sits near
**1.0** for 25–30 aa peptides and rises only for genuinely globular ones
(measured: 1.01 for a 26-mer helix, 1.08 for a 38-mer). That is the correct
answer, not a regression.

> **Open question.** The correction lowers scores by roughly 0.6×. Boltzmann
> selection depends on `β(sᵢ − s_max)`, so selective pressure scales with the
> spread of the score. β is held at **20** per the derivation in thesis
> §2.4.8.3, but that derivation predates this correction — worth re-checking
> against the new score distribution before the next production series.

ESM3 emits **backbone-only** structures (N, CA, C, O), so Cβ is reconstructed
with the standard virtual-Cβ formula (Yang et al. 2020); verified at
CA–CB 1.529 Å, sd 0.0006.

## Hemolysis is an attribute, not a selection pressure

HemoPI2 runs every generation and writes `hemo_prob` per candidate, but it does
**not** drive selection. Screen and rank candidates by `hemo_prob` afterwards,
or pass `--hemo-in-score` to fold `× (1 − hemo_prob)` into the fitness.
`PFES_SKIP_HEMO=1` drops the per-generation call entirely.

## `amp_src`: which scorer answered

MACREL is defined for **10–100 residues**. Outside that window
`macrel_score_batch` substitutes the biophysical `calculate_samp` surrogate per
sequence, which changes what the objective measures. `progress.log` therefore
carries an **`amp_src`** column — `macrel` or `proxy` — per row. Check it
before interpreting `amp_prob` from any run whose population approached 100 aa.

## Before you launch: `preflight.sh`

Every launcher sources `preflight.sh`, which aborts if MACREL, HemoPI2 or ESM3
are not actually callable, and probes magainin-2 end-to-end to confirm MACREL
answered rather than the surrogate.

This exists because the v2 production series — four runs, ~60 h — completed with
neither classifier installed. `score.py` fell back per sequence, warned only on
stderr, and each run's header still printed "MACREL AMP + PFES". Every
`amp_prob` in that series was the surrogate. **`conda activate pfes_amps`
resolves to whichever env of that name comes first, and only one has macrel.**

## Install

```bash
conda create -n pfes_amps python=3.11 -y && conda activate pfes_amps
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install -r requirements.txt
conda install -c bioconda -c conda-forge macrel -y
export HF_TOKEN=...      # huggingface.co/EvolutionaryScale/esm3-sm-open-v1
./preflight.sh           # must pass before any run
```

`requirements.txt` pins `onnxruntime<=1.25.1`: 1.26 changed the shape of ONNX
`output_probability`, so MACREL returned raw decision values instead of
calibrated probabilities — magainin-2 came back at **−0.050** and was classified
NOT an AMP, sending every candidate to the surrogate. `psique` is bundled
(`bin/psique`); no install needed.

## Run

```bash
python pfes.py --start file --start-file init/init_fragments.faa \
  -ps 100 -ng 600 -sm weak -b 20 --strong_selection_after_n_gen 480 \
  --norepeat --max-tokens-per-batch 4096 \
  -o results/frag-r1 2>&1 | tee results/frag-r1.console.log
```

Starting populations are built by `analysis/make_init_sets.py` (`--random`,
`--uniprot`, `--orfs`); each arm is AMP-screened, topped up to `--pop`, and
writes a provenance TSV. Reuse the same file across branches so every arm starts
identically.

## Analyse

```bash
python visual_pfes.py -l results/frag-r1/progress.log -s results/frag-r1/structures \
       -o results/frag-r1/analysis --notraj
python analysis/score_posthoc.py results/frag-r1/ -o results/comparison
```

`score_posthoc.py` re-scores survivors with **AMPlify** (Li et al. 2022), an
attentive BiLSTM unrelated to MACREL's random forest that never saw the search,
so its agreement is independent evidence in a way MACREL's own score is not.
It needs its own environment: `export AMPLIFY_CMD="conda run -n amplify AMPlify"`.

## Reproducibility

Runs are **not** reproducible, matching upstream PFES, which seeds no RNG.
Mutation, survivor sampling and random-sequence generation draw from OS entropy.
Repeat runs are independent samples of the same process — which is what makes
replication meaningful — but an individual trajectory cannot be regenerated and
must be preserved from its log and `structures/`. The starting population is the
one thing that is fixed, and it is fixed by a file, not a seed.
