# The objective, term by term

Every term in the fitness, what it is for, where its constants come from, and
how much selection pressure it actually applies. Four of the seven terms turn
out to be inert at this chain length; the activity classifier does its work in
the first ten generations and then saturates; and for the remaining 590 the
strongest voice in the objective is a metric used outside its published domain
of validity.

Measured on `sorf-r1` and `sorf-r2`, the only two runs of the previous series
that had MACREL and HemoPI2 genuinely installed. Those runs have since been
deleted along with the rest ([`VOID-RUNS.md`](VOID-RUNS.md)); they are
recoverable from git history at `0942520`, and the numbers below are what they
measured. Spread is log₁₀(max/min) across the 100 members of a generation — a
pure number comparable across terms, because the fitness is a product, so
`log(score) = Σ log(terms)` and these are additive contributions to the quantity
selection ranks on.

**The ordering changes over a run, and the table below is the endpoint.** The
term that dominates depends on when you look:

| gen | pTM | `P_L` | **MACREL** | contacts | pLDDT | mean `amp_prob` |
|---|---|---|---|---|---|---|
| 0 | 0.938 | 0.295 | **1.473** | 0.378 | 0.362 | 0.256 |
| 10 | 0.459 | 0.246 | **0.296** | 0.156 | 0.098 | 0.607 |
| 50 | 0.158 | 0.055 | **0.137** | 0.066 | 0.041 | 0.892 |
| 200 | 0.220 | 0.114 | **0.108** | 0.063 | 0.037 | 0.934 |
| 450 | 0.115 | 0.042 | **0.101** | 0.079 | 0.023 | 0.944 |
| 599 | 0.069 | 0.042 | **0.036** | 0.050 | 0.004 | 0.971 |

MACREL dominates generation 0 overwhelmingly — 1.47 against pTM's 0.94 — and is
what drags a random population into AMP space, lifting mean `amp_prob` from
0.256 to 0.607 within ten generations and to 0.918 by generation 100. Then it
saturates. **From generation 10 onward pTM is the largest single term for the
remaining 590 generations**, comparable to MACREL through the mid-run and ahead
of it at the end.

This does not mean the peptides are not antimicrobial; they are, at 0.92–0.97
from generation 100. It means that what selects *among* AMP-like candidates for
98 % of the run is a superposition score whose tolerance at this chain length is
0.96 Å — below the coordinate error of any single-sequence predictor. Two
candidates both at `amp_prob` 0.95 are separated because one has pTM 0.58 and
the other 0.55, and at 26 residues that is not a fold-quality difference.

| Term | Constants | Published? | log₁₀ spread | Status |
|---|---|---|---|---|
| **pTM** | — | Eq. 4 | **0.069 / 0.050** | **largest pressure — outside its domain** |
| **length penalty** `P_L` | `L0=30, C=0.12` | **changed** from `250, 0.2` | 0.042 / 0.044 | binding, and justified below |
| **MACREL** `amp_prob` | — | added by this fork | 0.036 / 0.022 | dominant at gen 0 (1.47), saturated by gen 100 |
| helix penalty `P_α` | `L0=30, C=0.5` | Eq. 7, unchanged | 0.009 / 0.022 | inert — `P_L` binds first |
| contact density | Cβ, 6 Å, \|i−j\|>5 | Eq. 5, restored | 0.007 | near-inert by construction |
| pLDDT | — | Eq. 4 | 0.004 | saturated at 0.990 |
| beta penalty `P_β` | `L0=12, C=0.5` | Eq. 8, unchanged | **0.000** | completely inert |

---

## pTM is the largest term and it is out of range

TM-score normalises by a length-dependent distance scale
(Zhang & Skolnick 2004):

```
d0 = 1.24 · (L − 15)^(1/3) − 1.8
```

| L (residues) | 19 | 25 | **26** | 30 | 40 | 100 | 250 |
|---|---|---|---|---|---|---|---|
| d0 (Å) | 0.17 | 0.87 | **0.96** | 1.26 | 1.83 | 3.65 | 5.85 |

At the 26 residues this search converges to, **d0 is 0.96 Å** — sub-Ångström. A
Cα displaced by one Ångström already contributes under half weight. PFES was
built for 250-residue proteins, where d0 is 5.85 Å; the same score at 26
residues is a six-fold tighter question.

Two published boundaries are relevant. The original TM-score derivation
**contained no structure shorter than 40 residues**, and `d0` goes **negative
below 19** — the metric is undefined there, and this objective operates 14
residues below the first boundary and 7 above the second
([Dunbrack 2025](https://www.biorxiv.org/content/10.1101/2025.02.10.637595v2)).

This is why the winners sit at pTM 0.55–0.60, a figure the now-deleted
`ANALYSIS.md` explained away as "normal and acceptable for short peptides". It
is neither normal nor a fold-quality statement: it is what a global
superposition score returns when its tolerance is smaller than the coordinate
error of any predictor. And because the term still *varies* (sd 0.010–0.015),
that variation is the single strongest thing selection sees for most of a run.

**Kept for the final series, and reported as a limitation.** Removing a term
from the published objective is a larger claim than this thesis can support
from one measurement, and `control-fold-only` — which is pLDDT × pTM alone — already
measures what these two terms select for on their own. What changes is the
claim: the objective is not "fold quality × activity", it is "a
length-mismatched superposition score × activity", and the ranking above is the
evidence.

## β = 20 is not weak selection

`-sm weak -b 20` is used throughout, and §2.4.8.3 derives 20. The source paper
is explicit about what that number means:

> With β = [0, ∞], where β = 0 means no selection … β = 1 resembles stochastic
> selection as in Eq. 2, and **β > 5 starts behaving as strong deterministic
> selection as in Eq. 1** (SI Appendix, Fig. S10).
> — Sahakyan et al. 2025, *Methods*

β = 20 is four times past the threshold at which the authors say the Gibbs
selection collapses onto deterministic truncation. The runs are therefore under
near-deterministic selection from generation 1, and
`--strong_selection_after_n_gen 480` switches from strong selection to strong
selection.

**This matters specifically for the origin experiment.** A search that cannot
drift cannot leave the basin it started in, so near-deterministic selection is
*itself* a mechanism that would produce an apparent origin effect. Any result
of the form "where the search begins decides where it ends" has to be stated
against the selection strength that produced it, and at β = 20 the honest
statement is that the search was never given the opportunity to forget its
starting point.

**Kept at 20 for the final series** — it is the value the thesis derives and the
runs must match the thesis — but described as strong selection throughout, and
§2.4.8.3 re-checked against Fig. S10.

## The length penalty is the one constant this fork changed, and it is justified

Published Eq. 6 is `C = 0.2, L0 = 250`. This fork runs `C = 0.12, L0 = 30`,
derived in §2.4.7.1 by the rule *`L0` is the midpoint of the published length
range of the target class, and `C` is set so the 90→10 % transition spans that
range*.

The AMP length literature supports both:

- Anuran AMPs: **99.9 % below 50 residues, mean length 24**, mean net charge
  +2.5 ([Amphibian AMP analysis, 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7459754/)).
- APD6: synthetic AMPs average **19.09 aa**; the database's own inclusion
  boundary is **under 100 residues**
  ([APD6](https://academic.oup.com/nar/article/54/D1/D363/8250474),
  [APD3](https://academic.oup.com/nar/article/44/D1/D1087/2503090)).
- MACREL is defined over **10–100 residues**
  ([Santos-Júnior et al. 2020](https://peerj.com/articles/10555/)), which is
  where the long arm's `pl0 = 100` cap comes from.
- The class is described at **<40 residues, net charge +2 to +9**
  ([Zhang et al. 2021](https://doi.org/10.1186/s40779-021-00343-2)).

The populations converge to **26–27 residues**, against a natural mean of 24 and
a synthetic mean of 19. The length term is doing what it was retuned to do.

One caveat the same literature supplies: the runs reach net charge +6 to +7,
against +2.5 for the anuran set and +2 to +9 for the class. Inside the range,
but at its top — consistent with the search maximising the most reliable
descriptor MACREL carries.

## The two secondary-structure penalties are inherited unchanged, and both are inert

`P_α` (`C = 0.5, L0 = 30`) and `P_β` (`C = 0.5, L0 = 12`) are the published
constants, untouched.

| helix length | 15 | 20 | 22 | 25 | 28 | 30 | 35 |
|---|---|---|---|---|---|---|---|
| `P_α` | 0.999 | 0.993 | 0.982 | 0.924 | 0.731 | 0.500 | 0.076 |

Measured: the longest helical run in the final populations is **18.9 and 22.6
residues** in chains of 26 and 27, so `P_α` returns 0.98–0.99 and applies no
pressure. It cannot: the length penalty caps the chain at 26–27 before a helix
long enough to be penalised can exist. **`P_α` is redundant given `P_L` at this
chain length.**

That redundancy ends on the long arm. At `pl0 = 100` chains reach 35–40
residues, a 30-residue helix takes `P_α = 0.5` and a 35-residue helix 0.076 —
which would forbid exactly the architecture the class is built on, LL-37 being
37 residues with a helix spanning about 30
([Cardoso et al. 2021](https://doi.org/10.1007/s12551-021-00784-y)). Raising
`hl0` to 40 there was correct and is now quantitatively supported.

`P_β` returns **exactly 1.000, sd 0.000** — no member of either final
population has a strand longer than 12. Natural β-sheet AMPs have strands of
5–7 residues and depend on disulfide bonds that neither the fitness nor ESM3
models reliably, so that subclass is unreachable regardless, and the term
serves only to forbid extended sheets. It is kept because it is free and
because it is the published definition, not because it does anything.

## Contact density measures compactness, and that is why it is near-inert

Eq. 5 as published: **Cβ atoms, within 6 Å, more than 5 residues apart**. The
`|i − j| > 5` rule excludes the i→i+4 α-helical register; the previous
implementation admitted it, and 95–98 % of counted contacts were helical turns
(`0d06bb8`). Corrected, the term contributes 0.007 against 0.037 before.

That is the right answer, and the reason is biological. Contact density
measures tertiary packing — the defining property of the globular domains PFES
was built for. Membrane-active AMPs are disordered in aqueous solution and fold
to an amphipathic helix only on contact with a bilayer
([Cardoso 2021](https://doi.org/10.1007/s12551-021-00784-y),
[Zhang 2021](https://doi.org/10.1186/s40779-021-00343-2)), so a compactness
term asks this class for a property it should not have.

ESM3 emits backbone-only structures, so Cβ is reconstructed with the standard
virtual-Cβ formula ([Yang et al. 2020](https://doi.org/10.1073/pnas.1914677117)),
verified at CA–CB 1.529 Å, sd 0.0006.

## pLDDT is saturated

Spread 0.004 at a mean of 0.990. It stopped discriminating early and now
functions as a floor rather than a gradient. Worth stating rather than hiding,
because a high pLDDT on this class is not straightforwardly good news: AMPs are
disordered in water, so a confident single-conformer prediction is a statement
about the model's prior, not about the molecule in solution.

## What never enters the score, and why

**Hemolysis** (HemoPI2, [Chaudhary et al. 2016](https://doi.org/10.1038/srep22843))
is logged every sampled generation and never selected on. The therapeutic index
is the ratio of minimum haemolytic to minimum inhibitory concentration
([Cardoso 2021](https://doi.org/10.1007/s12551-021-00784-y)), so a hemolysis
readout is half of what separates an antimicrobial from a detergent — but an
optimised safety term returns a number about how well the search satisfied it,
whereas a held-out one reports on where the search went.

**AMPlify** ([Li et al. 2022](https://doi.org/10.1186/s12864-022-08310-4)) audits
post hoc. It saturates at exactly 1.0000 for a third to two thirds of AMP-like
candidates, so as an objective it supplies no gradient, and it needs Python 3.6
with an old TensorFlow. The asymmetry with MACREL is permanent by construction,
not a stage of work left undone.

**No amphipathicity term.** Hydrophobic moment and net charge are what separate
active from inactive peptides at scale
([Wang et al. 2017](https://doi.org/10.3390/molecules22112037)), but MACREL
already carries the hydrophobic moment among its 22 features, so an explicit
term would weight one descriptor twice and hand the optimiser a single scalar to
game. Measured, the winners reach μH 0.826 against 0.780 for magainin 2,
melittin, LL-37, pexiganan and cecropin B — natural values, without the term.

## Why the objective is a product of models that fail differently

Optimising a learned predictor drives the search to where the predictor is
unreliable; the canonical statement is
[Brookes, Park & Listgarten (ICML 2019)](https://proceedings.mlr.press/v97/brookes19a.html),
whose example failure mode is sequences that will not fold. Multiplying an
activity term by fold-quality and geometry terms forces a candidate to satisfy
models that fail in uncorrelated ways. The auditor split — generate under
classifier guidance, screen with filters the generator never saw — follows
[Das et al. (Nat Biomed Eng 2021)](https://www.nature.com/articles/s41551-021-00689-x).
Classifier guidance is separately known to narrow generated peptide diversity
([Brief Bioinform 2025](https://academic.oup.com/bib/article/26/5/bbaf500/8301249)),
which is the second reason AMPlify and HemoPI2 audit rather than drive.

## Summary for the final series

Nothing is changed for the production runs. What changes is what may be claimed:

1. **pTM is the strongest selection pressure**, at roughly twice the activity
   term, and it is evaluated 14 residues below the floor of its own derivation.
2. **β = 20 is strong selection by the source paper's criterion**, not weak, so
   no run in this project has had an exploratory phase.
3. **Four of seven terms are inert** — `P_β` exactly, `P_α` by redundancy with
   `P_L`, contact density by construction, pLDDT by saturation.
4. **The objective is effectively `pTM × P_L × MACREL`**, in that order of
   strength, with the activity term third.

Points 1 and 2 are the two open questions the next series should be designed to
close, and neither needs new compute to state.
