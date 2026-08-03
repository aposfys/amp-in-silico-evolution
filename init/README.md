# Starting populations — the origin comparison

Three sets of 100 sequences differing in where the sequences came from. The
question they answer is whether antimicrobial activity emerges as readily from
sequence with no evolutionary history, from novel protein-coding sequence, or
from sequence that already encodes something else.

| File | Origin | Status |
|---|---|---|
| `init_random.faa` | completely random, uniform over the 20 residues | ready |
| `init_orfs.faa` | small ORFs from metazoan transcriptomes | **awaiting the source FASTA** |
| `init_fragments.faa` | fragments of conserved metazoan proteins that are not AMPs | ready |

Each is fed to a run as a fixed starting population:

```bash
python pfes.py --start file --start-file init/init_fragments.faa -pop 100 ...
```

## Two length variants

The sets exist in two versions, cut by the same code from the same sampling
frame and differing only in the length of the window.

| Directory | Cut | Median `P_len` at generation 0 |
|---|---|---|
| `init/` | fixed 25 aa | 0.731 |
| `init_varlen/` | uniform random 10–100 aa | 0.015 |

`init/` is the version to run a comparison on. `init_varlen/` implements random
cutting in the fullest sense, at random length as well as random position, over
the range MACREL is defined for. The reason it is kept separate is below, and it
is not the general point that length is a confound.

The structured fitness carries `P_len = 1 − sigmoid(L, pl0, 0.2)` with
`pl0 = 30` (`pfes.py:242`), which is `1 / (1 + e^{0.2(L − 30)})`:

| L | 15 | 25 | 30 | 40 | 50 | 70 | 100 |
|---|---|---|---|---|---|---|---|
| `P_len` | 0.95 | 0.73 | 0.50 | 0.119 | 0.018 | 0.00034 | 0.0000008 |

In `init_varlen/` the median fragment is 51 aa, so its median `P_len` is 0.015
and 49 of the 100 fragments start below 0.01. Since the fitness is a product,
half the set scores essentially zero at generation zero, and under
`--selection_mode strong`, a straight truncation to the top `pop_size`, it is
gone in the first round. Two things then break. The run stops measuring where
sequence came from and starts measuring which fragments happened to be short.
And because the long ones go at once, the arm collapses to a few dozen genuinely
competitive individuals, an unequal starting *sample size* that no later
conditioning on length repairs. Its random set is drawn over the same 10–100 so
the two arms at least match each other, but both are crippled in the same way,
which fixes the comparison and not the experiment.

If the requirement is that cutting be random and the run still has to be
informative, the middle range is 15–30, where the penalty is present but
functional, `P_len` from 0.95 to 0.50:

```bash
python analysis/make_init_sets.py --random --uniprot --pop 100 --len 15-30 --seed 42 --no-screen -o init_1530/
```

## How each was built

### `init_random.faa`

```bash
python analysis/make_init_sets.py --random --pop 100 --len 25 --seed 42 -o init/
```

Uniform sampling over the 20 standard residues. Mean pairwise identity 5.0%,
which is what independent random draws give.

### `init_fragments.faa`

```bash
python analysis/make_init_sets.py --uniprot --pop 100 --len 25 --seed 42 --no-screen -o init/
```

**Conserved, defined as a UniRef50 cluster old enough to predate the split
between the animal phyla.** UniRef50 groups sequences at 50% identity, so the
size of a cluster counts how often a protein recurs across genomes rather than
how often it has been studied, and the cluster's *common taxon* is the last
common ancestor of every member it holds. The sampling frame is the 12,841
clusters that have a metazoan member, hold at least 100 UniProtKB entries, and
whose common ancestor is Bilateria or older — Bilateria, Eumetazoa, Metazoa,
Opisthokonta, Eukaryota or cellular organisms. A cluster whose ancestor is
Mammalia or Euteleostomi is rejected however large it is, because it is
conserved within one lineage rather than across many species. The 100 clusters
actually used have a median size of 360 members and reach 3,016.

**Sampling is family-first.** Clusters are shuffled under the seed and walked in
order; one Swiss-Prot entry is taken from each. Sampling entries directly does
**not** work — it returns whatever UniProt has curated most, which is human, and
lets one well-studied family contribute a dozen fragments. Within a cluster the
species is drawn before the entry, for the same reason one level down. 267
clusters were queried to fill 100, the other 167 having no reviewed metazoan
entry that survived the filters. No two fragments come from the same family, let
alone the same protein.

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
`frag_001_O00159_Homo_sapiens_803-827` gives accession, organism and the residue
range cut. `init_fragments.tsv` carries the same 100 rows with the source protein
name, phylum, UniRef50 cluster, cluster size, common taxon, source length and
window.

What comes out is housekeeping machinery — V-type ATPase subunit d2, COP9
signalosome subunit 3, dynein axonemal heavy chain 1, GMP reductase 1,
multifunctional protein CAD, 26S proteasome subunit Rpn3. 100 species-diverse families across six phyla:
Chordata 84, Arthropoda 12, and one each of Mollusca, Echinodermata, Cnidaria
and Nematoda, from 24 source species. The chordate share reflects Swiss-Prot's
curation, not the sampling frame; the *families* are bilaterian-wide by
construction, whichever species the sequence was read from.

The window is random within the protein, so a fragment is a piece of a conserved
protein rather than necessarily a piece of a conserved *region* — some windows
land in fast-evolving loops. Restricting the cut to annotated Pfam domains is the
change that would tighten this, if the distinction ever matters to a result.

### `init_orfs.faa`

Not yet built. Once the source FASTA arrives:

```bash
python analysis/make_init_sets.py --orfs <source.faa> --pop 100 --len 25 --seed 42 -o init/
```

Same cutting logic, one window per ORF, and whichever `--len` the other two
sets in that directory were cut to — 25 for `init/`, `10-100` for
`init_varlen/`. An arm cut to a different distribution from the arms it is
compared against is measuring length, not origin.

## Outstanding: the AMP screen

These sets are meant to start *without* antimicrobial activity, so anything
MACREL already calls an AMP should be dropped. That screen needs MACREL on PATH
and has **not** been applied to the files as committed, because it was not
available where they were generated.

`init_fragments.faa` was therefore built with `--no-screen` deliberately. Without
MACREL, `score.py` falls back to the biophysical surrogate, and that surrogate is
a different and much harsher standard — it dropped 6 of the first 7 candidates in
a trial run. Screening one set on the surrogate while `init_random.faa` sits
unscreened would confound the comparison at the point it is meant to be cleanest.

Regenerate on a machine with MACREL before running:

```bash
python analysis/make_init_sets.py --random --uniprot --max-amp-prob 0.5 -o init/
```

Or screen the existing files and remove anything above the threshold. Either way
the same screen must be applied to all three sets, or the comparison is between
sets filtered to different standards.

## Composition

| Set | Length | Mean charge | Charge range | Hydrophobic | Pairwise identity |
|---|---|---|---|---|---|
| random | 25 | +0.22 | −4.0 to +6.1 | 41% | 4.9% |
| fragments | 25 | −0.24 | −5.0 to +5.1 | 38% | 5.8% |

Hydrophobic is the fraction of AVLIMFWC and charge is net charge at pH 7; both
rows are computed the same way. The two sets are matched in length and close in
mean charge and hydrophobic content, so neither starts with the amphipathic
composition an AMP needs, which is the property the runs are meant to evolve.
They differ where real sequence differs from random: the fragments carry natural
amino-acid usage, L, A, E, K, S, V most frequent, against a flat distribution in
the random set, and a charge range of −5.0 to +5.1 that is skewed
negative where the random set's is not, because real protein regions are locally
charged in a way independent draws are not.

Pairwise identity near 5% in both confirms no family is over-represented, and for
the fragments the construction guarantees it anyway: 100 distinct UniRef50
families and 100 distinct proteins.
