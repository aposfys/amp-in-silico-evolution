# ANALYSIS — fitness-macrel-pfes

MACREL + structural penalties: fitness = pLDDT × pTM × len × SS × AMP × (1−hemo) × contacts/len. Recommended branch for drug candidate generation.

---

## 1. Design rationale

### Why structural penalties?

Without length constraints, directed evolution exploits the `d` (full duplication) mutation operator to grow sequences to 60–80 aa. Longer sequences accumulate higher absolute AMP probability but are impractical as drug candidates: expensive to synthesise, potentially immunogenic, outside the canonical AMP therapeutic range.

Three penalties constrain the optimizer:

**Length penalty:** Sigmoid-based decay for sequences exceeding the target window (default ≤30 aa). Prevents runaway growth via duplication.

**Helix / beta penalties:** Penalise excessive secondary structure fractions. In practice both stay near 1.0 — the cationic amphipathic helix the optimizer converges on naturally satisfies these constraints.

### The contact density term: why > 1.0?

```
contact_term = (num_Cα_contacts_within_8Å + seq_len) / seq_len
```

This is the only term in the fitness that can exceed 1.0. For a compact 26 aa α-helix with ~28 contacts: (28+26)/26 = **2.077**. Without this multiplier, the best pfes sequence would score ~0.35 — too low to survive Boltzmann selection (β=20). The contact density booster bridges the gap, rewarding compact well-packed structures and making short peptides competitive.

### MACREL hemolytic failure and the proxy

MACREL's hemolytic classifier outputs exactly 0.000 for all evolved cationic sequences (see `fitness-macrel/ANALYSIS.md` for root cause). The biophysical proxy used here:

```
hemo = sigmoid(hydro_fraction × 10 − min(charge_pH7.4, 8.0) × 0.5 − 2.0)
```

The charge is capped at +8. Without the cap, sequences with charge +12 (typical for pfes sequences) get hemo ≈ 0.047, which is realistic. With an uncapped formula and extreme charge (+19, as in the MACREL-only run), hemo collapses to 0.001 — an artificial free pass. The cap forces a genuine balance between charge and hydrophobicity.

---

## 2. Production run results

### v1 run (original proxy, no charge cap)

**Parameters:** pop_size=100, generations=500, β=20, start_len=24 aa, Mac MPS.

| Metric | Value |
|--------|-------|
| Best score | **0.7240** |
| At generation | gndx483 |
| Mean population score at gen 499 | ~0.69 |

Convergence is slower than the unconstrained branch (~400 generations vs ~150) because the optimizer must satisfy eight constraints simultaneously. Local optima are common — improving AMP probability often worsens a structural penalty.

### v1 best sequence

```
RRVYKKFFAPRVRLKRLAKAIKLVRK
```

| Property | Value |
|----------|-------|
| Length | 26 aa |
| pLDDT | 0.980 |
| pTM | 0.550 |
| AMP probability | 0.993 |
| Hemolytic proxy | 0.047 (v1) / 0.269 (v2 proxy) |
| Net charge | +12 (K=6, R=6) |
| Hydrophobic fraction | 50% (F, V, A, I, L) |
| Contact term | (28+26)/26 = 2.077 |

**Score breakdown:** 0.980 × 0.550 × 0.70 × 1.0 × 1.0 × 0.993 × 0.953 × 2.077 = **0.724**

**Why pTM = 0.55 (lower than the unconstrained run)?** pTM measures global structural similarity. A 26 aa sequence naturally scores lower on pTM than a 70 aa sequence — there is simply less structure to align. For short peptides (< 30 aa), pTM of 0.5–0.6 is normal and acceptable. The low pTM is a consequence of the length constraint, not poor fold quality.

**Drug candidate assessment:**
- ✅ 26 aa — practical synthesis cost (~€100–200/mg)
- ✅ Charge +12 — strong electrostatic interaction with bacterial membranes
- ✅ 50% hydrophobic — appropriate for amphipathic helix
- ⚠️ K+R fraction 46% — more cationic than most natural AMPs; may benefit from iterative charge reduction
- Comparison: magainin-2 (23 aa, charge +4), LL-37 (37 aa, charge +6), pexiganan (22 aa, charge +9)

### v2 run (improved proxy, charge cap at +8)

A v2 re-run was launched after fixing the hemolytic proxy. With the charge cap, the optimizer can no longer suppress hemo by extreme K/R accumulation alone — it must genuinely balance hydrophobicity and charge. Expected outcome: sequences with slightly lower K/R content and more defined amphipathic structure.

---

## 3. Sequence length dynamics

The population stabilises at **~26 aa** within 30 generations and stays there throughout the run. The small drift from 24 aa (start) to 26 aa reflects a subtle selection pressure: at 26 aa, one extra helical turn adds ~2 contacts, boosting the contact density term.

This is in contrast to the unconstrained `fitness-macrel` branch where sequences grow to ~71 aa by generation 499.

---

## 4. Amino acid composition (final population)

**Enriched:**
- R: ~22% (from 5% baseline)
- K: ~20%
- F: ~8% (hydrophobic, amphipathic face)
- V: ~7% (compact hydrophobic)

**Depleted:**
- D, E, G, P: near zero

Total K+R ≈ 42% — the optimizer converges on the most cationic solution compatible with the structural penalties. This exceeds most natural AMPs (magainin-2: 17% K+R, defensins: ~22%). The high cationicity reflects the optimizer using charge as the primary mechanism to suppress the hemolytic proxy term.

---

## 5. Contact density drives short-peptide competitiveness

The mathematical mechanism that enables this branch to work:

| Sequence | Length | Contacts | Contact term | Score without | Score with |
|----------|--------|----------|-------------|---------------|------------|
| Best pfes (v1) | 26 aa | 28 | 2.077 | 0.349 | 0.724 |
| Hypothetical 24 aa | 24 | 22 | 1.917 | 0.349 | 0.669 |
| Hypothetical 30 aa | 30 | 35 | 2.167 | 0.349 | 0.756 |

Without the contact density multiplier, all short peptides would score below 0.40 and be eliminated from the population within a few generations. The term is essential for this fitness function to work.

---

## 6. Structural penalties in practice

**Length penalty:** Active constraint throughout — stabilises at 0.65–0.70. The optimizer pushes against the length ceiling because longer sequences gain more contacts. The penalty prevents this.

**Helix penalty / beta penalty:** Near 1.0 throughout. The evolved cationic amphipathic helix naturally satisfies both — no excessive helix, zero beta. These constraints are non-binding for this sequence class.

---

## 7. Known limitations

1. **Hemolytic proxy is imprecise.** The biophysical proxy cannot distinguish magainin-2 (non-hemolytic) from mastoparan (hemolytic) using composition descriptors alone. The proxy provides selection pressure but is not a reliable hemolysis predictor. Experimental hemolysis assays are required before clinical evaluation.

2. **MACREL AMP training distribution.** High AMP probabilities (>0.99) may reflect similarity to MACREL training data rather than guaranteed in vitro activity. The evolved sequences are novel — experimental MIC measurement is the only true validation.

3. **ESM3 is a predictor, not a solver.** pLDDT=0.98 does not mean the sequence folds in a membrane environment as ESM3 predicts. AMPs operate at the lipid-water interface — a context not well represented in ESM3's training.

4. **Population convergence.** By generation 500 the population is near-monomorphic. For candidate selection, use the top 10–20 sequences, not just the single best.

---

## 8. Recommended next steps

**Computational:**
- Use the v2 run results (charge-capped proxy) for candidate selection
- Run `visual_pfes.py` on the v2 results to compare with v1
- Benchmark top candidates against AMP databases (APD3, DRAMP) for novelty check

**Experimental (top candidates):**
1. Synthesise 5–10 top sequences from the v2 run (~€150/mg each)
2. MIC against E. coli, S. aureus, P. aeruginosa
3. Hemolysis assay (0.5% human RBC, 100 μg/mL peptide)
4. CD spectroscopy in TFE/SDS micelles (confirm helical structure)
5. If positive: serum stability, cytotoxicity panel (HEK293, HepG2)
