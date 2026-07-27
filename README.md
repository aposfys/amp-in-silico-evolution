# PFES — fitness-macrel-structured

Evolves antimicrobial peptides using **MACREL** (AMP classifier), **HemoPI2** (hemolysis, logged as an attribute), **ESM3** (fold), and the full **PFES structural penalties** (length, secondary structure, contact density).

> **2×2 matrix** (classifier × structural penalties): `fitness-macrel-foldonly` · `fitness-macrel-structured` · `fitness-amplify-foldonly` · `fitness-amplify-structured`. All share ESM3 folding + HemoPI2 + DRAMP seeding.

## Fitness score
```
score = pLDDT × pTM × length_penalty × helix_penalty × beta_penalty × AMP_probability × contact_density
```
| Term | Source | Range | Note |
|---|---|---|---|
| pLDDT | ESM3 | 0–1 | per-residue fold confidence |
| pTM | ESM3 | 0–1 | overall-fold confidence |
| AMP_probability | MACREL | 0–1 | probability the peptide is antimicrobial |
| length_penalty | PFES | 0–1 | →0 above ~30 aa (caps length) |
| helix_penalty | PFES (PSIQUE) | 0–1 | →0 for α-helices > ~20 residues |
| beta_penalty | PFES (PSIQUE) | 0–1 | →0 for β-strands > ~12 residues |
| contact_density | PFES | > 1 | rewards compact folds |
| (1 − hemolysis) | HemoPI2 | 0–1 | **OFF by default** — only in the score with `--hemo-in-score` |

> **Hemolysis is an ATTRIBUTE by default:** HemoPI2 runs every generation and writes `hemo_prob` for each peptide, but it does **not** affect selection. Screen/rank candidates by `hemo_prob` afterwards, **or** add `--hemo-in-score` to make `× (1 − hemolysis)` part of the fitness.

## Install
```bash
# 1. clone + branch
git clone https://github.com/aposfys/PFES-AMPs.git && cd PFES-AMPs
git checkout fitness-macrel-structured

# 2. environment
conda create -n pfes_amps python=3.11 -y && conda activate pfes_amps
pip install torch --index-url https://download.pytorch.org/whl/cpu    # CPU build (drop --index-url for GPU)
pip install -r requirements.txt
pip install "numpy<2" "pillow<12"        # HemoPI2/scikit-learn need numpy < 2

# 3. MACREL (Bioconda)  — the AMP classifier
conda install -c bioconda -c conda-forge macrel -y

# 4. psique (secondary structure) is BUNDLED (bin/psique + psique.py) — no install needed

# 5. ESM3 access token (accept licence at huggingface.co/EvolutionaryScale/esm3-sm-open-v1)
export HF_TOKEN=your_token_here
```

Verify:
```bash
python -c "from esm.models.esm3 import ESM3; print('ESM3 ok')"
python -c "import score; print(score.macrel_score_batch(['GIGKFLHSAKKFGKAFVGEIMNS']))"
```

---

## How to run — step by step

### 1. One-time setup: seed database + hemolysis labels
```bash
export PFES_AMP_DB=$PWD/data/dbaasp.faa
python amp_db.py --fetch          # ~10,700 validated AMPs from DRAMP -> data/dbaasp.faa
python amp_db.py --annotate-hemo  # HemoPI2 hemo_risk label per DB peptide (needs hemopi2)
```

### 2. Build the starting population — pick an arm
```bash
# MIX (50 existing AMPs + 50 random at the existing mean length)
python amp_db.py --make-init --pop 100 --frac-existing 0.5 --seed 42 -o init_mix.faa
# EXISTING ONLY (100 known AMPs)
python amp_db.py --make-init --pop 100 --frac-existing 1.0 --seed 42 -o init_existing.faa
# RANDOM ONLY (100 random sequences)
python amp_db.py --make-init --pop 100 --frac-existing 0.0 --rand-len 25 --seed 42 -o init_random.faa
```
Using your OWN AMPs instead of DRAMP: `python amp_db.py --import your_amps.faa` first, then `--make-init --frac-existing 1.0`.
The same init file is reused across branches so every arm starts identically (a fair comparison).

### 3. Smoke test first (a few minutes — confirms the tools engage)
```bash
python pfes.py --start file --start-file init_mix.faa -ps 8 -ng 3 -sm weak -b 20 -o /tmp/smoke
```
Check: no "falling back" warnings; the `hemo` column varies; MACREL + HemoPI2 actually ran.

### 4. Full run
```bash
nice -5 python pfes.py \
  --start file --start-file init_mix.faa \
  -ps 100 -ng 100 -sm weak -b 20 \
  --norepeat --max-tokens-per-batch 4096 \
  -o results/macrel-structured 2>&1 | tee results/macrel-structured.log
```
- **Hemolysis:** attribute-only by default (add `--hemo-in-score` to put it in the fitness).
- **Selection:** `-sm weak -b 20` (fitness-proportional, keeps diversity) or `-sm strong` (greedy top-N).
- **CPU:** with `--max-tokens-per-batch 4096` a generation is ~1–2 min on ~48 cores → ~2–3 h for 100 generations. **Run inside `tmux`** so a disconnect can't kill it (detach: Ctrl-b then d).

### 5. Analyse (graphs)
```bash
python visual_pfes.py -l results/macrel-structured/progress.log -s results/macrel-structured/structures -o results/macrel-structured/analysis --notraj
```
→ `Summary.png`, `Evolution.png`, `Score_components.png`, `Fitness_landscape.png`, `AA_composition.png`, plus per-column plots (`plots/score.png`, `amp_prob.png`, `hemo_prob.png`, `seq_len.png`).

### 6. Compare arms / branches
Run the same command with a different `--start-file` (existing / random / mix) or a different branch (change `git checkout` + `-o`). Identical start ⇒ differences are attributable to the design choice.

---

## Key options
| Flag | Default | Meaning |
|---|---|---|
| `--start {random,randoms,existing,mix,file}` | random | starting population; `file` loads `--start-file` |
| `--start-file PATH` / `--mix-frac F` | — / 0.5 | fixed init FASTA / existing-vs-random ratio for mix |
| `--hemo-in-score` | off | put HemoPI2 hemolysis back into the fitness score |
| `-ng` / `-ps` | 100 / 10 | generations / population size |
| `-sm` / `-b` | weak / 1 | selection: `weak` (β via `-b`), `weak2`, `strong` |
| `--norepeat` | off | never re-evaluate a sequence already seen |
| `-pl0 / -hl0 / -bl0` | 30 / 20 / 12 | length / helix / β penalty thresholds
| `--max-tokens-per-batch` | 512 | fold batch size (use ~4096 to minimise classifier reloads) |

## Reproducibility

Runs are **not** reproducible, matching upstream PFES, which seeds no random number
generator anywhere. Mutation, survivor sampling and random-sequence generation all
draw from OS entropy, so the same command produces a different trajectory each time.
Repeat runs are independent samples of the same process, which is what makes
replication across runs meaningful, but an individual result cannot be regenerated
and must be preserved from its log and `structures/` output.

The starting population is the one thing that *is* fixed, and it is fixed by a file
rather than by a seed: build it once with `amp_db.py --make-init` and pass the same
`--start-file` to every arm.

## Known constraints

**MACREL is defined for 10–100 residues.** Outside that window `macrel_score_batch`
falls back to the `calculate_samp` biophysical proxy, silently, per sequence. The
objective therefore changes identity mid-run if chains grow past 100 aa. This is not
hypothetical: the `foldonly` arm reached **93 residues by generation 300**. Fold-only
runs should not be extended much beyond 300 generations without either capping length
or logging which scorer produced each value.

**Deduplication cost grows quadratically.** `--norepeat` scans the full
`ancestral_memory` for every candidate, and that table gains `pop_size` rows per
generation. Doubling the generation count roughly quadruples the dedup work, so
runtime does not scale linearly in `-ng`.

**Measured CPU throughput** (12 threads, `--max-tokens-per-batch 4096`, pop 100):

| Arm | s / candidate | h / 300 generations |
|---|---|---|
| structured | 1.70 | 14.2 |
| fold-only | 2.14 | 17.8 |

Fold-only is the more expensive arm because folding cost rises with chain length.

## Output (`<outdir>/progress.log`, tab-separated)
`gndx` generation · `seq_len` · `prot_len_penalty` `max_alpha_penalty` `max_beta_penalty` · `ptm` `mean_plddt` (ESM3) · `num_conts` · `amp_prob` (MACREL) · **`hemo_prob` (HemoPI2 — the attribute)** · `score` · `sequence` · `mutation` · `ss`. Structures: gzipped PDB in `<outdir>/structures/`.
