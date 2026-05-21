# PFES: Protein Fold Evolution Simulation
## Branch: `fitness-macrel-pfes` — MACREL ML AMP Score + Full Structural Penalties + ESM3

This branch evolves antimicrobial peptides using **MACREL** (ML AMP classifier) combined with the full set of PFES structural penalties and ESM3 fold quality. It is the most selective scoring strategy: both AMP fitness and structural quality drive evolution simultaneously.

---

### Fitness score

**single_chain mode:**
```
score = pLDDT × pTM × len_penalty × helix_penalty × beta_penalty
        × AMP_probability × (1 - hemolytic_probability) × contact_term
```

**inter_chain mode:**
```
score = pLDDT × pTM × iPLDDT × len_penalty × helix_penalty × beta_penalty
        × AMP_probability × (1 - hemolytic_probability) × contact_term × inter_contact_term
```

- **pLDDT / pTM**: ESM3 per-residue confidence and predicted TM-score [0, 1]
- **AMP_probability**: MACREL prediction that the sequence is antimicrobial [0, 1]
- **hemolytic_probability**: MACREL prediction of hemolytic activity — penalised as `(1 - hemo)` [0, 1]
- **len_penalty**: discourages sequences longer than `--prot_len_penalty` (default 30 AA)
- **helix_penalty**: discourages helices longer than `--helix_len_penalty` (default 20 residues)
- **beta_penalty**: discourages beta strands longer than `--beta_len_penalty` (default 12 residues)
- **contact_term**: rewards compact folds with more intra-chain contacts
- **inter_contact_term** / **iPLDDT**: inter-chain contacts and interface pLDDT (inter_chain mode only)

MACREL is valid for 10–100 AA sequences. Sequences outside this range fall back to biophysical `s_amp`.

---

### Installation

```bash
# 1. Clone and enter the repo
git clone https://github.com/aposfys/PFES-AMPs.git
cd PFES-AMPs
git checkout fitness-macrel-pfes

# 2. Create environment (conda recommended)
conda create -n pfes-macrel python=3.10
conda activate pfes-macrel

# 3. Install Python dependencies
pip install -r requirements.txt

# 4. Install MACREL (requires Bioconda channel)
conda install -c bioconda -c conda-forge macrel

# 5. Install PSIQUE (secondary structure assignment)
pip install git+https://github.com/sahakyanhk/psique

# 6. Set your HuggingFace token (required for ESM3)
#    Accept the model license at https://huggingface.co/EvolutionaryScale/esm3-sm-open-v1
export HF_TOKEN=your_token_here
```

Verify MACREL works:
```bash
macrel --version
```

---

### Quick test (CPU, 5 generations)

```bash
python pfes.py \
  -em single_chain -sm weak -b 20 \
  -ps 4 -ng 5 \
  --random_seq_len 20 \
  --max-tokens-per-batch 256 \
  -o test_macrel_pfes
```

Expected output: 5 generations each printing `score`, `pLDDT`, `pTM`, `AMP`, `hemo`, `len`, `mutation`, `sequence`.

Verify the run succeeded:
```bash
wc -l test_macrel_pfes/progress.log         # should be > 20 lines
grep "s_amp" test_macrel_pfes/progress.log | head -1  # AMP probabilities should appear
```

---

### Full AMP evolution run

```bash
# GPU (recommended)
python pfes.py \
  -em single_chain -sm weak -b 20 \
  -ps 50 -ng 200 \
  --random_seq_len 24 \
  -pl0 30 -hl0 20 -bl0 12 \
  --norepeat \
  -o amp_macrel_pfes_run

# CPU (reduce population and generations)
python pfes.py \
  -em single_chain -sm weak -b 20 \
  -ps 8 -ng 50 \
  --random_seq_len 24 \
  --max-tokens-per-batch 256 \
  -o amp_macrel_pfes_cpu
```

**Note on scoring pressure:** The product of 8–10 terms all in [0,1] produces scores in the 0.001–0.01 range. Use `-b 20` to `-b 100` for strong Boltzmann selection, or `-sm strong` for deterministic top-N. Default `-b 1` is nearly neutral.

---

### Analyse results

```bash
python visual_pfes.py \
  -l amp_macrel_pfes_run/progress.log \
  -s amp_macrel_pfes_run/structures/ \
  -o amp_macrel_pfes_run/analysis/
```

Produces:
- `analysis/Summary.png` — pLDDT, pTM, score, length, AMP probability trajectories
- `analysis/Secondary_structures.png` — secondary structure along the lineage
- `analysis/lineage.tsv` — best evolutionary path
- `analysis/pfestraj.pdb` — backbone trajectory (open in PyMOL)

---

### Branch comparison

| Branch | AMP scoring | Structural penalties | External tool |
|---|---|---|---|
| `fitness-samp` | Biophysical s_amp | Full (length, SS, contacts) | None |
| `fitness-macrel` | MACREL ML classifier | None (ESM3 confidence only) | MACREL |
| **fitness-macrel-pfes** (this) | MACREL ML classifier | Full (length, SS, contacts) | MACREL |

---

### Hardware

- NVIDIA V100 / A100 (GPU, recommended — `--max-tokens-per-batch 2048`)
- Apple Silicon M-series via MPS
- CPU (functional but slow — use `-ps 4 -ng 20 --max-tokens-per-batch 256`)

Reference: [In silico evolution of globular protein folds from random sequences](https://www.pnas.org/doi/10.1073/pnas.2509015122)

Extended data: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14061036.svg)](https://doi.org/10.5281/zenodo.14061036)
