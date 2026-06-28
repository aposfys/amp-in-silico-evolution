# PFES — fitness-macrel-structured

Evolves antimicrobial peptides using **MACREL** (ML AMP classifier), **HemoPI2** (ML hemolysis predictor), **ESM3** (fold quality), and a full set of **structural penalties** (length, secondary structure, contact density). This is the recommended branch for generating drug candidates — it constrains sequences to therapeutically realistic lengths (~20–30 aa).

> Renamed from `fitness-macrel-pfes` (both branches run on the PFES framework; the suffix now reflects that this one applies the full PFES **structural** fitness terms). The length-unconstrained variant is `fitness-macrel-foldonly`.

## Fitness formula

```
score = pLDDT × pTM × length_penalty × helix_penalty × beta_penalty
        × AMP_probability × (1 − hemolytic_probability)
        × (num_contacts + seq_len) / seq_len
```

| Term | Source | Notes |
|------|--------|-------|
| pLDDT | ESM3 per-residue confidence | [0, 1] |
| pTM | ESM3 predicted TM-score | [0, 1] |
| length_penalty | Sigmoid penalty for sequences > target length | [0, 1] |
| helix_penalty | Penalty for excessive α-helix content | [0, 1] |
| beta_penalty | Penalty for excessive β-strand content | [0, 1] |
| AMP_probability | MACREL ML classifier | [0, 1] |
| 1 − hemolytic_probability | HemoPI2 ML predictor (falls back to biophysical proxy) | [0, 1] |
| (contacts + len) / len | Contact density — rewards compact folds | > 1.0 for helices |

The contact density term is the only term that can exceed 1.0. A compact 26 aa helix with ~28 contacts gives (28+26)/26 = **2.08**, which compensates for the structural penalties and keeps short peptides competitive.

---

## Installation

```bash
# 1. Clone and switch to this branch
git clone https://github.com/aposfys/PFES-AMPs.git && cd PFES-AMPs
git checkout fitness-macrel-structured

# 2. Create environment
conda create -n pfes_amps python=3.11
conda activate pfes_amps

# 3. Python dependencies (includes peptides + hemopi2)
pip install -r requirements.txt

# 4. MACREL (Bioconda)
conda install -c bioconda -c conda-forge macrel

# 5. PSIQUE (secondary structure)
pip install git+https://github.com/sahakyanhk/psique

# 6. ESM3 access token
#    Accept the licence at: https://huggingface.co/EvolutionaryScale/esm3-sm-open-v1
export HF_TOKEN=your_token_here
```

Verify:
```bash
macrel --version
python -c "from score import macrel_score_batch; print('OK')"
```

---

## Quick test (3 generations)

```bash
python pfes.py \
  -em single_chain -sm weak -b 20 \
  -ps 4 -ng 3 \
  --random_seq_len 20 \
  --max-tokens-per-batch 256 \
  -o test_macrel_pfes
```

Each generation prints: `score  pLDDT  pTM  AMP  hemo  len  mutation  sequence`

Scores will be low early on (0.001–0.05 range) — the product of 8 terms is inherently small. They rise as the optimizer finds sequences that satisfy all constraints simultaneously.

Verify:
```bash
wc -l test_macrel_pfes/progress.log   # > 15 lines
head -2 test_macrel_pfes/progress.log # check column headers
```

---

## Full production run

```bash
# Mac MPS / GPU (recommended)
python pfes.py \
  -em single_chain -sm weak -b 20 \
  -ps 100 -ng 500 \
  --random_seq_len 24 \
  --norepeat \
  -o results_macrel_pfes

# CPU only
python pfes.py \
  -em single_chain -sm weak -b 20 \
  -ps 20 -ng 200 \
  --random_seq_len 24 \
  --max-tokens-per-batch 256 \
  -o results_macrel_pfes_cpu
```

**Structural penalty flags:**

| Flag | Default | Controls |
|------|---------|---------|
| `-pl0 N` | 30 | Max sequence length before penalty (aa) |
| `-hl0 N` | 20 | Max helix length before penalty (residues) |
| `-bl0 N` | 12 | Max beta strand length before penalty (residues) |

**Selection flags:**

| Flag | Meaning |
|------|---------|
| `-b 20` | Boltzmann β — use 20–50 for meaningful selection |
| `-sm weak` | Boltzmann-weighted sampling |
| `-sm strong` | Deterministic top-N selection |
| `--norepeat` | Reject duplicate sequences |

---

## Output

All results in `<outdir>/progress.log` (tab-separated):

```
gndx              generation index
seq_len           length in amino acids
prot_len_penalty  length penalty [0,1]
max_alpha_penalty helix fraction penalty [0,1]
max_beta_penalty  beta fraction penalty [0,1]
ptm               ESM3 predicted TM-score [0,1]
mean_plddt        ESM3 per-residue confidence [0,1]
num_conts         Cα contacts within 8 Å
amp_prob          AMP probability from MACREL [0,1]
hemo_prob         hemolytic proxy [0,1]
score             total fitness score
sequence          amino acid sequence
mutation          mutation that produced this sequence
ss                secondary structure string (PSIQUE)
```

Structures saved as PDB files under `<outdir>/structures/`.

---

## Analyse results

```bash
python visual_pfes.py \
  -l results_macrel_pfes/progress.log \
  -s results_macrel_pfes/structures/ \
  -o results_macrel_pfes/analysis/
```

Produces `Summary.png`, `Evolution.png`, `Score_components.png`, `Secondary_structures.png`, `Fitness_landscape.png`, `AA_composition.png`, `lineage.tsv`, `pfestraj.pdb`.
