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

## Parameters (full reference)

All options are command-line flags — `python pfes.py --help` prints them.

### Run size & control
| Flag | Default | Meaning |
|------|---------|---------|
| `-ng`, `--num_generations` | 100 | number of generations |
| `-ps`, `--pop_size` | 10 | sequences per generation (population size) |
| `--seed` | none | RNG seed. Same seed **and** same `-iseq` ⇒ identical start across branches |
| `--norepeat` | off | never re-evaluate / re-select a sequence already seen (diversity + saves folds) |
| `-o`, `--outpath` | `output` | output directory |
| `-l`, `--log` | `progress.log` | log file name |
| `--nobackup` | off | overwrite existing output instead of backing up |

### Starting population — **`--start`** (the easy way to choose random / existing / mix)
This is the high-level switch for *what the population starts from*:

| `--start` | What you get | extra flags |
|-----------|--------------|-------------|
| `random` | one random sequence copied across the population | `--random_seq_len N` (24) |
| `randoms` | a different random sequence in every slot | `--random_seq_len N` |
| `existing` | known AMPs from the database (diverse) | — |
| `mix` | **existing AMPs + random sequences** at the existing seeds' mean length | `--mix-frac F` (0.5 = 50/50) |
| `file` | a **fixed** population loaded from a FASTA (identical start for every branch) | `--start-file PATH` |
| `seq` | a literal sequence | give it via `-iseq <SEQUENCE>` |

Example: `--start mix --mix-frac 0.5` → 50% existing AMPs + 50% random (your 50/50 plan).

<details><summary>Underlying <code>-iseq</code> values (advanced; <code>--start</code> just sets these)</summary>

`random` · `randoms` · `<SEQUENCE>` · `db`/`db:random` · `db:diverse` · `db:low`/`db:high` · `db:<name>` (e.g. `db:magainin_2`) · `db:mix[:frac]` · `file:<path>`
</details>

### Mutation & selection
| Flag | Default | Meaning |
|------|---------|---------|
| `-ed`, `--evoldict` | `flatrates` | mutation-rate set: `flatrates`, `codonrates`, `flatoptim`, `uniprotrates` |
| `-sm`, `--selection_mode` | `weak` | **`strong`** = deterministic top-N (fittest survive); `weak` = fitness-proportional sampling; `weak2` |
| `-b`, `--beta` | 1 | selection sharpness for the `weak` modes (20–50 typical). **Ignored when `-sm strong`.** |

### Structural penalties
*On the **structured** branches these drive the score; on the **foldonly** branches they are logged but not part of the score.*
| Flag | Default | Meaning |
|------|---------|---------|
| `-pl0`, `--prot_len_penalty` | 30 | sequence length (aa) above which the length penalty kicks in |
| `-hl0`, `--helix_len_penalty` | 20 | α-helix length threshold |
| `-bl0`, `--beta_len_penalty` | 12 | β-strand length threshold |

### ESM3 / performance / mode
| Flag | Default | Meaning |
|------|---------|---------|
| `--num-recycles` | 1 | ESM3 denoising steps per fold (higher = better fold, slower) |
| `--max-tokens-per-batch` | 512 | fold batch size; lower if you hit OOM (256 for CPU) |
| `-em`, `--evolution_mode` | `single_chain` | `single_chain`, `inter_chain`, `multimer` |

---

## 1. Set up the AMP seed database (once)

```bash
pip install hemopi2
python amp_db.py --fetch            # DRAMP ~10.7k AMPs -> data/dbaasp.faa
python amp_db.py --annotate-hemo    # adds the hemo_risk attribute (uses HemoPI2)
```

## 2. Build ONE shared starting population

```bash
python amp_db.py --make-init --pop 100 --frac-existing 0.5 --seed 42 -o init_pop.faa
# 50 existing AMPs + 50 random at the existing mean length; reuse this file everywhere
```

## 3. Smoke test first (~5 min) — confirm the real tools engage

```bash
python pfes.py --start file --start-file init_pop.faa -ps 8 -ng 3 -sm strong --seed 42 -o /tmp/smoke
```
Check: **no "falling back" warnings**, `hemo_prob` varies in the log, the AMP classifier + HemoPI2 actually ran.

## 4. Full run — strong evolution

```bash
python pfes.py \
  -em single_chain -sm strong \
  --start file --start-file init_pop.faa \
  -ps 100 -ng 300 \
  --seed 42 --norepeat \
  -o results/macrel-structured
```
`-sm strong` keeps the deterministic top-`pop_size` sequences each generation (greediest selection). For CPU add `--max-tokens-per-batch 256` and use a smaller `-ps`/`-ng`.

> Don't want a fixed file? Use `--start mix --mix-frac 0.5` to generate 50 existing + 50 random on the fly — but for a *fair 4-branch comparison* use the shared `--start file` so every branch starts identically.

## 5. Compare the four branches (2×2)

Run the **identical** command on each branch (same `init_pop.faa` + `--seed`), changing only `-o`:
```bash
for B in macrel-foldonly macrel-structured amplify-foldonly amplify-structured; do
  git checkout fitness-$B
  python pfes.py --start file --start-file init_pop.faa -ps 100 -ng 300 -sm strong --seed 42 --norepeat -o results/$B
done
# amplify branches also need:  export AMPLIFY_CMD="conda run -n amplify AMPlify"
```
Because the start is identical, differences are attributable to the two factors: classifier (MACREL vs AMPlify) and structural penalties (off vs on).

---

## Output (`<outdir>/progress.log`, tab-separated)

```
gndx               generation index            amp_prob   AMP probability (MACREL/AMPlify)
id                 sequence id                 hemo_prob  hemolysis probability (HemoPI2)
seq_len            length (aa)                 score      total fitness score
prot_len_penalty   length penalty [0,1]        sequence   amino-acid sequence
max_alpha_penalty  helix penalty [0,1]         mutation   mutation that produced it
max_beta_penalty   beta penalty  [0,1]         prev_id    parent id
ptm / mean_plddt   ESM3 confidence [0,1]       ss         secondary structure (PSIQUE)
num_conts          Cα contacts within 8 Å
```
Structures are saved as gzipped PDB under `<outdir>/structures/`.

## Analyse / graphs

```bash
python visual_pfes.py \
  -l results/macrel-structured/progress.log \
  -s results/macrel-structured/structures/ \
  -o results/macrel-structured/analysis/
```
Produces `Summary.png`, `Evolution.png`, `Score_components.png`, `Fitness_landscape.png`, `AA_composition.png`, `Secondary_structures.png`, per-column plots, and `lineage.tsv`.
