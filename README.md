# PFES-AMPs

In silico evolution of antimicrobial peptides. **ESM3** folds each candidate,
**MACREL** scores it for antimicrobial activity, and the **PFES** structural
terms constrain length and secondary structure. **HemoPI2** records hemolysis
and **AMPlify** re-scores survivors, both as attributes that never touch
selection. A population of one hundred peptides is mutated and selected for six
hundred generations.

Built for the MSc thesis *Investigating novel antimicrobial peptides using
machine learning, structure prediction and in silico evolution* (Apostolos
Fysekidis, Bioinformatics, National and Kapodistrian University of Athens), on
top of [PFES](https://github.com/sahakyanhk/pfes) by Sahakyan et al.

## The objective

`main` **is** the production arm. Its fitness is:

```
score = pLDDT × pTM × length_penalty × helix_penalty × beta_penalty
        × AMP_probability × contact_density
```

| Term | Source | Range | In the score? |
|---|---|---|---|
| pLDDT | ESM3 | 0–1 | **yes** |
| pTM | ESM3 | 0–1 | **yes** |
| length_penalty | PFES Eq. 6 | 0–1 | **yes** — →0 above `-pl0` (default 30) |
| helix_penalty | PFES Eq. 7 / PSIQUE | 0–1 | **yes** |
| beta_penalty | PFES Eq. 8 / PSIQUE | 0–1 | **yes** |
| AMP_probability | MACREL | 0–1 | **yes** — the only classifier driving selection |
| contact_density | PFES Eq. 5 | ≥ 1 | **yes** |
| hemo_prob | HemoPI2 | 0–1 | **no** — logged; `--hemo-in-score` to enable |
| AMPlify probability | AMPlify | 0–1 | **no** — post-hoc only, never in the loop |

Hemolysis is measured but never selected on. It is measured at all because it
is half of what separates an antimicrobial peptide from a detergent: an ideal
candidate has *"high therapeutic activity and minimum hemolytic activity"*
([Chaudhary 2016][c16]), the therapeutic index being the ratio of minimum
haemolytic to minimum inhibitory concentration ([Cardoso 2021][c21]). It is not
selected on because optimising a predictor makes its output a restatement of the
objective rather than evidence about the design — the same reason AMPlify audits
post-hoc rather than driving the search.

## Where to look

Two arms, run from identical starting populations so the difference between them
is attributable to the objective alone.

| Branch | Objective | |
|---|---|---|
| **`main`** | pLDDT × pTM × length × helix × β × **MACREL** × contacts | the full objective — **the production arm** |
| **`fitness-esm3`** | pLDDT × pTM | **the control** — fold confidence alone, no structural terms, no classifier |
| `pfes-original` | — | PFES exactly as published, **unmodified**; the reference this parameterisation is audited against (was `upstream`, originally `alpha`) |

Both arms compute and log every quantity either one selects on, so their
`progress.log` files share a schema and compare column by column — including
`amp_prob` on the control, which MACREL scores but never drives.

[`production/`](production/) holds the design study the thesis reports: the
starting populations, the launch scripts, and one compressed archive per run.

## Contact density measures compactness, not helicity

Contacts follow Eq. 5 of Sahakyan et al. exactly: **Cβ atoms, within 6 Å, more
than 5 residues apart in sequence**, both above the pLDDT floor.

The `|i − j| > 5` rule is the whole point — it excludes the i→i+4 α-helical
register, so a long helix cannot score as a compact fold. Upstream PFES loops
`range(i + 4, …)`, admitting |i − j| = 4 and 5 and therefore counting exactly
those helical contacts; its 6 Å cutoff mostly hides this, but a wider cutoff
un-hides it. Measured on this project's own output at 8 Å, **95–98 % of counted
contacts were |i − j| = 4** — the term was counting helical turns. Corrected in
`0d06bb8`.

A single α-helix has no tertiary contacts, so the corrected term sits near
**1.0** for 25–30 aa peptides and rises only for genuinely globular ones
(measured: 1.01 for a 26-mer helix, 1.08 for a 38-mer). That is the right
answer, not a regression — and the reason is biological, not arithmetic.
Contact density measures *tertiary packing*, the defining property of the
globular domains PFES was built to evolve. Membrane-active AMPs never occupy
that state: they are disordered in aqueous solution and fold to an amphipathic
helix only on contact with a lipid bilayer ([Cardoso 2021][c21],
[Zhang 2021][z21]). A compactness term asks this class for a property it should
not have.

**No amphipathicity term replaces it, deliberately.** Hydrophobic moment and net
charge are what separate active from inactive peptides at scale
([Wang 2017][w17]) — but MACREL already carries the hydrophobic moment among its
features, so an explicit term would weight one descriptor twice and hand the
optimiser a single scalar to game. It is also unnecessary: measured with
MACREL's own implementation the winners average **μH 0.826** against **0.780**
for magainin 2, melittin, LL-37, pexiganan and cecropin B, at net charges of +4
to +8 — inside the +2 to +9 reported for the class ([Zhang 2021][z21]).

> **Open:** the correction lowers scores by roughly 0.6×, and Boltzmann
> selection depends on `β(sᵢ − s_max)`, so pressure scales with score spread. β
> is held at **20** per thesis §2.4.8.3, but that derivation predates the
> correction and is worth re-checking against the new distribution.

ESM3 emits **backbone-only** structures (N, CA, C, O), so Cβ is reconstructed
with the standard virtual-Cβ formula (Yang et al. 2020) — verified at CA–CB
1.529 Å, sd 0.0006.

## `amp_src`: which scorer answered

MACREL is defined for **10–100 residues**. Outside that window the biophysical
`calculate_samp` surrogate is substituted per sequence, changing what the
objective measures. `progress.log` therefore carries an **`amp_src`** column —
`macrel` or `proxy` — per row. Check it before interpreting `amp_prob` from any
run whose population approached 100 aa.

## Before you launch: `preflight.sh`

Every launcher sources `preflight.sh`, which aborts if MACREL, HemoPI2 or ESM3
are not callable and probes magainin-2 end-to-end to confirm MACREL answered
rather than the surrogate.

This exists because the v2 series — four runs, ~60 h — completed with neither
classifier installed. `score.py` fell back per sequence, warned only on stderr,
and each run's header still printed "MACREL AMP + PFES". Every `amp_prob` in
that series was the surrogate. **`conda activate pfes_amps` resolves to
whichever env of that name comes first, and only one has macrel.**

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
calibrated probabilities — magainin-2 came back at **−0.050**, classified NOT an
AMP, sending every candidate to the surrogate. `psique` is bundled
(`bin/psique`); no install needed.

## Run

```bash
python pfes.py --start file --start-file init/init_fragments.faa \
  -ps 100 -ng 600 -sm weak -b 20 --strong_selection_after_n_gen 480 \
  --norepeat --max-tokens-per-batch 4096 \
  -o results/frag-r1 2>&1 | tee results/frag-r1.console.log
```

Starting populations come from `analysis/make_init_sets.py` (`--random`,
`--uniprot`, `--orfs`); each arm is AMP-screened, topped up to `--pop`, and
writes a provenance TSV. Reuse the same file across both branches so every arm
starts identically.

## Analyse

```bash
python visual_pfes.py -l results/frag-r1/progress.log -s results/frag-r1/structures \
       -o results/frag-r1/analysis --notraj
python analysis/score_posthoc.py results/frag-r1/ -o results/comparison
```

`score_posthoc.py` re-scores survivors with **AMPlify** (Li et al. 2022), an
attentive BiLSTM unrelated to MACREL's random forest that never saw the search,
so its agreement is independent evidence in a way MACREL's own score is not. It
needs its own environment:
`export AMPLIFY_CMD="conda run -n amplify AMPlify"`.

## Reproducibility

Runs are **not** reproducible, matching upstream PFES, which seeds no RNG.
Mutation, survivor sampling and random-sequence generation draw from OS entropy.
Repeat runs are independent samples of the same process — which is what makes
replication meaningful — but an individual trajectory cannot be regenerated and
must be preserved from its log and `structures/`. The starting population is the
one thing that is fixed, and it is fixed by a file, not a seed.

## References

[c21]: https://doi.org/10.1007/s12551-021-00784-y
[z21]: https://doi.org/10.1186/s40779-021-00343-2
[w17]: https://doi.org/10.3390/molecules22112037
[c16]: https://doi.org/10.1038/srep22843

- **Sahakyan H., Babajanyan S. G., Wolf Y. I. & Koonin E. V. (2025).** In silico
  evolution of globular protein folds from random sequences. *PNAS* **122**,
  e2509015122. — the objective, Eqs. 4–8, and the contact definition this fork
  restores.
- **Hayes T. et al. (2025).** Simulating 500 million years of evolution with a
  language model. *Science* **387**, 850–858. — ESM3; Table S9 is the
  structure-prediction benchmark quoted above.
- **Santos-Júnior C. D., Pan S., Zhao X.-M. & Coelho L. P. (2020).** Macrel:
  antimicrobial peptide screening in genomes and metagenomes. *PeerJ* **8**,
  e10555. — the activity classifier, its 22 features, and its 10–100 residue
  domain.
- **Li C. et al. (2022).** AMPlify: attentive deep learning model for discovery
  of novel antimicrobial peptides effective against WHO priority pathogens.
  *BMC Genomics* **23**, 77. — the architecturally unrelated post-hoc auditor.
- **[Cardoso P. et al. (2021)][c21].** Molecular engineering of antimicrobial
  peptides. *Biophys Rev* **13**, 35–69. — AMPs disordered in water, helical
  only on the membrane; the therapeutic index.
- **[Zhang Q.-Y. et al. (2021)][z21].** Antimicrobial peptides: mechanism of
  action, activity and clinical potential. *Mil Med Res* **8**, 48. — the same
  conformational point, and the <40 aa / +2 to +9 range behind `-pl0`.
- **[Wang C.-K., Shih L.-Y. & Chang K. Y. (2017)][w17].** Large-scale analysis of
  antimicrobial activities in relation to amphipathicity and charge.
  *Molecules* **22**, 2037. — why amphipathicity is the property that matters.
- **[Chaudhary K. et al. (2016)][c16].** A web server and mobile app for
  computing hemolytic potency of peptides. *Sci Rep* **6**, 22843. — HemoPI;
  why a hemolysis readout is kept at all.
- **Yang J. et al. (2020).** Improved protein structure prediction using
  predicted interresidue orientations. *PNAS* **117**, 1496–1503. — the
  virtual-Cβ construction used where ESM3 emits backbone only.
