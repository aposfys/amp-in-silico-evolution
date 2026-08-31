# Starting populations — the origin comparison

Three sets of 100 sequences differing in where the sequences came from. The
question they answer is whether antimicrobial activity emerges as readily from
sequence with no evolutionary history, from novel protein-coding sequence, or
from sequence that already encodes something else.

> **The sequence files themselves are not published in this repository.** What
> follows documents exactly how the three sets were built, screened and
> verified, and `analysis/make_init_sets.py` is the code that builds them. The
> delivered sets used for the reported runs are available from the author —
> apostolosfysekidis1@gmail.com. They are held back so that anyone reproducing
> the published results makes contact first.

| File | Origin | Status |
|---|---|---|
| `init_random.faa` | completely random, uniform over the 20 residues | on request |
| `init_orfs.faa` | small ORFs from metazoan transcriptomes | on request |
| `init_fragments.faa` | fragments of conserved metazoan proteins that are not AMPs | on request |

Each is fed to a run as a fixed starting population:

```bash
python pfes.py --start file --start-file init/init_fragments.faa -ps 100 ...
```

## The AMP screen — applied

These sets are meant to start *without* antimicrobial activity, so anything
MACREL already calls an AMP is dropped. Both directories were built with
`--max-amp-prob 0.5`, and the delivered files were verified against
`macrel_score_batch`, the same function the run itself calls:

| File | above 0.5 | median | max |
|---|---|---|---|
| `init/init_random.faa` | **0 / 100** | 0.188 | 0.495 |
| `init/init_fragments.faa` | **0 / 100** | 0.139 | 0.495 |
| `init/init_orfs.faa` | **0 / 100** | 0.202 | 0.495 |
| `init_varlen/init_random.faa` | **0 / 100** | — | — |
| `init_varlen/init_fragments.faa` | **0 / 100** | — | — |

> **Verified 2026-08-30.** The sORF arm was re-scored with
> `score.macrel_score_batch`, the same function the run calls, and comes back
> **0 / 100 above 0.5, median 0.202, max 0.495** — the same screen signature as
> the other two. The three delivered sets are therefore all verified against
> the function that scores them at run time, and the row above is measured
> rather than assumed.
>
> Median probability does differ across the arms — 0.139 fragments, 0.188
> random, **0.202 sORF** — so small ORFs do start marginally nearer the class
> than either comparator. That is a property of small ORFs, not of the
> assembly, and after screening no set contains a sequence MACREL would call an
> AMP.
>
> ```bash
> python -c "
> import score
> s=[l.strip() for l in open('init/init_orfs.faa') if not l.startswith('>')]
> d=score.macrel_score_batch(s); v=[t[0] for t in d.values()]
> print(f'above 0.5: {sum(x>0.5 for x in v)}/{len(v)}  max {max(v):.3f}')"
> ```
>
> Two files this document refers to do not exist: `init_orfs.tsv` (the
> per-fragment provenance record) and `init_orfs.draw1_prescreen.faa`. The
> `--orfs` provenance TSV was added in `cd673f6`, after these sets were built.
> **Do not regenerate the set to produce them** — screening changes which draws
> are accepted and therefore the whole RNG stream, so `make_init_sets.py` would
> return a different set of 100 and invalidate every run seeded from this file.
> The TSV has to be reconstructed from the existing FASTA headers against
> `source_sorfs.tsv`.

The maximum of exactly 0.495 in both `init/` files is the screen's signature:
MACREL is a random forest, so its probabilities fall on a discrete grid, and a
`< 0.5` cut leaves 0.495 as the largest admissible value.

Screening changes which draws are accepted, which shifts the whole RNG stream
downstream of the first rejection. A screened set is therefore **not a subset**
of the unscreened draw made under the same seed — it is a different set of 100.
The pre-screen draws are retained as `*.draw1_prescreen.faa` for comparison and
are not used by any run. They hold 2 (random) and 4 (fragments) sequences above
the threshold, maxima 0.564 and 0.604.

> **Historical note.** An earlier revision of this file stated that the screen
> had *not* been applied and that the sets were built with `--no-screen`. That
> was true when it was written and was not updated when the sets were
> regenerated with MACREL available. The production series of §3.2 ran on the
> screened files listed above.

## Two length variants

The sets exist in two versions, cut by the same code from the same sampling
frame and differing only in the length of the window.

| Directory | Cut | Median length | Median `P_len` at `pl0 = 30` | at `pl0 = 100` |
|---|---|---|---|---|
| `init/` | fixed 25 aa | 25 | 0.65 | 1.00 |
| `init_varlen/` | uniform random 10–100 aa | 55 (frag), 60 (rand) | 0.05 / 0.03 | 0.99 |

The structured fitness carries `P_len = 1 − sigmoid(L, pl0, c)` (`pfes.py:237`),
which is `1 / (1 + e^{c(L − pl0)})`, a sigmoid equal to one half **at** `pl0`:

| L | 15 | 25 | 30 | 40 | 50 | 70 | 100 |
|---|---|---|---|---|---|---|---|
| `pl0 = 30` | 0.86 | 0.65 | 0.50 | 0.231 | 0.083 | 0.008 | 0.0002 |
| `pl0 = 100` | 1.000 | 1.000 | 1.000 | 0.999 | 0.998 | 0.973 | 0.500 |

**The steepness is `c = 0.12`, not the upstream `0.2`.** A logistic falls from
0.9 to 0.1 over `2·ln9/c` residues — 22 residues at `c = 0.2`, 37 at `c = 0.12`
— *regardless of* `pl0`. Against the published 250-residue target of upstream
PFES, 22 residues is a narrow wall at the threshold; against a 30-residue target
it spans the whole class and clips its upper half. The value here follows the
rule derived in §2.4.7.1 of the thesis: `pl0` is the midpoint of the published
length range of the target class and `c` is set so the 90→10% transition spans
that range, which for antimicrobial peptides at 12–50 residues gives 31 and
0.116, rounded to 30 and 0.12. The β-strand steepness is likewise 0.5 rather
than 0.6, restoring the published value.

Both constants are hard-coded, not exposed on the command line, so a clone that
does not carry the working-tree state will silently run a different objective.
They are frozen under the tag `v2-production`.

**Which directory to use is determined by `pl0`, and the two must be matched.**

At `pl0 = 30`, `init/` is the only usable version. In `init_varlen/` the median
fragment is 55 aa, so its median `P_len` is 0.05 — a thirteen-fold handicap
against the 0.65 a 25-residue sequence carries, before any other term is
considered. Since the fitness is a product, and
`--selection_mode strong` is a straight truncation to the top `pop_size`, half
the set is gone in the first round. Two things then break. The run stops
measuring where sequence came from and starts measuring which fragments happened
to be short. And because the long ones go at once, the arm collapses to a few
dozen genuinely competitive individuals — an unequal starting *sample size* that
no later conditioning on length repairs.

At `pl0 = 100` the relation inverts. A 55-residue fragment now carries
`P_len = 0.995`, so the length term is inert across the whole starting
distribution and the run is free to explore to about 75 residues before the
penalty begins to bite. This is what `init_varlen/` was built for, and the
threshold coincides with the definitional boundary of the peptide class: APD3
admits sequences up to 100 aa and MACREL is defined over 10–100.

Note that this uses `pl0` with **cap semantics rather than target semantics**,
and so departs deliberately from the rule of §2.4.7.1. That rule derives `pl0`
as the midpoint of a target class; here there is no target, only a boundary past
which the activity predictor is undefined. Holding `c = 0.12` while raising
`pl0` to 100 places the 90→10% transition at 82–118 residues, which is a wall at
the definitional edge rather than a preference for any particular length —
which is the intent. Re-deriving under the rule for an extended class of 12–100
would instead give `pl0 = 56, c = 0.05`, a length term that is again actively
selective, merely centred higher. The two are different experiments and the
choice should be stated wherever the long arm is reported.

Raising `pl0` **requires raising `hl0` with it.** The helix penalty uses a
steeper logistic (`c = 0.5`), so `hl0 = 30` assigns a factor of 0.076 to a
35-residue helix. A long peptide under that constraint cannot form the single
extended amphipathic helix that characterises the class — LL-37 is 37 residues
with a helix spanning about 30 — and the optimiser will instead build short
disconnected elements. The matched setting is:

```
short arm:  -pl0 30   -hl0 30  -bl0 12
long arm:   -pl0 100  -hl0 40  -bl0 12
```

`bl0` stays at 12 in both. Natural β-sheet AMPs have strands of 5–7 residues
(protegrin-1, tachyplesin) and are stabilised by disulfide bonds that neither
the fitness function nor ESM3 models reliably, so that subclass is out of reach
regardless and the threshold serves only to forbid extended sheets.

## How each was built

### `init_random.faa`

```bash
python analysis/make_init_sets.py --random --pop 100 --len 25 --seed 42 --max-amp-prob 0.5 -o init/
```

Uniform sampling over the 20 standard residues. Mean pairwise identity 5.0%,
which is what independent random draws give.

### `init_fragments.faa`

```bash
python analysis/make_init_sets.py --uniprot --pop 100 --len 25 --seed 42 --max-amp-prob 0.5 -o init/
```

**Conserved, defined as a UniRef50 cluster old enough to predate the split
between the animal phyla.** UniRef50 groups sequences at 50% identity, so the
size of a cluster counts how often a protein recurs across genomes rather than
how often it has been studied, and the cluster's *common taxon* is the last
common ancestor of every member it holds. The sampling frame is the 12,841
clusters (of 74,622 metazoan clusters) that hold at least 100 UniProtKB entries
and whose common ancestor is Bilateria or older. A cluster whose ancestor is
Mammalia or Euteleostomi is rejected however large it is, because it is
conserved within one lineage rather than across many species.

**Sampling is family-first.** Clusters are shuffled under the seed and walked in
order; one Swiss-Prot entry is taken from each. Sampling entries directly does
**not** work — it returns whatever UniProt has curated most, which is human, and
lets one well-studied family contribute a dozen fragments. Within a cluster the
species is drawn before the entry, for the same reason one level down. No two
fragments come from the same family, let alone the same protein.

Excluded by keyword: Antimicrobial (KW-0929), Immunity (KW-0391), Innate
immunity (KW-0399), Adaptive immunity (KW-1064), Defensin (KW-0211), Antibiotic
(KW-0044), Antiviral protein (KW-0930), Antiviral defense (KW-0051), Fungicide
(KW-0295), Bacteriolytic enzyme (KW-0081), Toxin (KW-0800), Cytolysis (KW-0204),
Hemolysis (KW-0354), and annotated fragments.

Excluded by GO term, which matches descendants and so covers the subtree in each
case: defense response (GO:0006952), immune system process (GO:0002376), killing
of cells of another organism (GO:0031640), antimicrobial humoral response
(GO:0019730).

Excluded by name, because keywords and GO terms miss things: a regex over the
named AMP families and over the parent proteins whose *fragments* are known
AMPs — histones, haemoglobin, haemocyanin, lactoferrin, chromogranin, thrombin,
lysozyme, the ribosomal proteins — plus interferon, interleukin, chemokine,
complement, Toll and lectin. GAPDH is a worked example of why the GO filter
earns its place: it passes every keyword, and is removed by GO:0006952, which it
carries because its own fragments are antimicrobial in fish.

One random window per protein. Headers keep the provenance:
`frag_000_P10768_Homo_sapiens_164-188` gives accession, organism and the residue
range cut. `init_fragments.tsv` carries the same 100 rows with the source protein
name, phylum, UniRef50 cluster, cluster size, common taxon, source length and
window, and is the authoritative record of what the set contains.

What comes out is housekeeping machinery — S-formylglutathione hydrolase,
GTP-binding nuclear protein Ran, V-type ATPase subunits, proteasome subunits.

**Composition of `init/init_fragments.faa` as delivered** (from
`init_fragments.tsv`, verified against the FASTA):

| | |
|---|---|
| Distinct UniRef50 families | 100 |
| Distinct source proteins | 100 |
| Distinct source species | 29 |
| Phyla | Chordata 84, Arthropoda 12, Nematoda 3, Platyhelminthes 1 |
| Common taxon | Bilateria 45, cellular organisms 18, Eukaryota 11, Eumetazoa 10, Opisthokonta 8, Metazoa 8 |
| Cluster size | median 387, maximum 6,554 |

The chordate share reflects Swiss-Prot's curation, not the sampling frame; the
*families* are bilaterian-wide by construction, whichever species the sequence
was read from.

The window is random within the protein, so a fragment is a piece of a conserved
protein rather than necessarily a piece of a conserved *region* — some windows
land in fast-evolving loops. Restricting the cut to annotated Pfam domains is the
change that would tighten this, if the distinction ever matters to a result.

### `init_varlen/`

Same code, same seed, `--len 10-100`:

```bash
python analysis/make_init_sets.py --random --uniprot --pop 100 --len 10-100 --seed 42 --max-amp-prob 0.5 -o init_varlen/
```

Delivered: `init_random.faa` 100 sequences, 11–100 aa, mean 60, screen dropped
1/100. `init_fragments.faa` 100 sequences, 10–99 aa, mean 55, screen dropped
3/100; 352 clusters were queried to fill 100, the other 252 having no reviewed
metazoan entry that survived the filters. `init_varlen/init_fragments.tsv`
carries the per-fragment provenance.

### `init_orfs.faa`

Built and screened: 100 fragments, one window per source ORF, drawn from the
2,999 small ORFs in `source_sorfs.faa` (provenance in `source_sorfs.tsv`, built
by `analysis/build_sorf_source.py`). Per-fragment provenance is in
`init_orfs.tsv`.

```bash
python analysis/make_init_sets.py --orfs init/source_sorfs.faa --pop 100 --len 25     --seed 42 --max-amp-prob 0.5 -o init/
python analysis/make_init_sets.py --orfs init/source_sorfs.faa --pop 100 --len 10-100 --seed 42 --max-amp-prob 0.5 -o init_varlen/
```

The screen drops any window MACREL already calls antimicrobial and then tops up
from ORFs not already used, so the set is 100 fragments from 100 distinct source
ORFs. Without that top-up a screened ORF set comes out smaller than the random
and fragments arms and the three origins stop being size-matched.

Same cutting logic, one window per ORF, and whichever `--len` the other two sets
in that directory were cut to. An arm cut to a different distribution from the
arms it is compared against is measuring length, not origin.

## Composition

| Set | Length | Mean charge | Charge range | Hydrophobic | Pairwise identity |
|---|---|---|---|---|---|
| random | 25 | +0.02 | −4 to +5 | 41% | 4.9% |
| fragments | 25 | −0.19 | −7 to +5 | 37% | 5.8% |
| sORF | 25 | +0.49 | −8 to +7 | 35% | — |

Hydrophobic is the fraction of AVLIMFWC; charge counts K and R as +1 and D and E
as −1. (§2.4.8 of the thesis reports +0.05 and −0.19 from a pH-7 calculation that
also gives histidine a partial charge; the two agree to within that difference.)

### Matched against the analysis definitions

The table above uses this file's own conventions, which are **not** the ones
`analysis/score_posthoc.py` applies to the final populations: it counts
histidine at +0.1 and includes tyrosine among the hydrophobic residues
(`AVLIMFWYC`). Recomputed under those definitions, so that a starting value and
a final value mean the same thing:

| arm | n | unique | net charge | sd | hydrophobic | sd |
|---|---|---|---|---|---|---|
| random | 100 | 100 | +0.13 | 2.03 | 0.464 | 0.093 |
| fragments | 100 | 100 | −0.11 | 2.38 | 0.398 | 0.093 |
| sORF | 100 | 100 | +0.57 | 2.90 | 0.381 | 0.120 |

**The spread between arms is smaller than the spread within any one of them** —
0.68 in charge against standard deviations of 2.0–2.9, and 0.083 in hydrophobic
fraction against 0.09–0.12. The three origins are therefore indistinguishable in
the bulk composition the objective rewards, and differ only in provenance, which
is what the comparison requires. All three are 100 sequences, all distinct, all
exactly 25 residues.

One consequence worth stating, because it constrains how a result may be read:
these starting differences are far too small to carry a compositional difference
through to the final population by inheritance. Any endpoint separation of the
size reported in §3.2 — several charge units — has to be produced by the search,
not transported from the seed.

The two sets are matched in length and close in mean charge and hydrophobic
content, so neither starts with the amphipathic composition an AMP needs, which
is the property the runs are meant to evolve. They differ where real sequence
differs from random: the fragments carry natural amino-acid usage — L, A, E, K,
S, V most frequent — against a flat distribution in the random set, and a charge
range skewed negative where the random set's is symmetric, because real protein
regions are locally charged in a way independent draws are not.

Pairwise identity near 5% in both confirms no family is over-represented, and for
the fragments the construction guarantees it anyway: 100 distinct UniRef50
families and 100 distinct proteins.
