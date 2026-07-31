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

## Length

| Set | Length |
|---|---|
| `init_random.faa` | fixed 25 aa |
| `init_fragments.faa` | uniform random, 10–100 aa (median 51, observed 11–99) |

The fragments are cut at a random length as well as a random position, over the
range MACREL is defined for. The window length is drawn *before* the source
protein is chosen, so the distribution is the one asked for rather than one bent
towards long windows by which proteins happened to be long enough to supply them.

Note what this costs, since it is the reason the earlier version of these sets
was cut to one length. Length is a confound: the structured objective penalises
anything past `-pl0` (30 aa), so roughly two thirds of these fragments start with
a length penalty the 25 aa random set never pays, and a difference between the
two sets will be part origin and part length. Reading the comparison therefore
means conditioning on length rather than comparing the two sets wholesale, or
rebuilding the random set over the same range:

```bash
python analysis/make_init_sets.py --random --pop 100 --len 10-100 --seed 42 -o init/
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
python analysis/make_init_sets.py --uniprot --pop 100 --len 10-100 --seed 42 --no-screen -o init/
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
actually used have a median size of 362 members and reach 3,016.

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

One random window per protein, at a random length. Headers keep the provenance:
`frag_001_O00159_Homo_sapiens_803-827` gives accession, organism and the residue
range cut. `init_fragments.tsv` carries the same 100 rows with the source protein
name, phylum, UniRef50 cluster, cluster size, common taxon, source length and
fragment length.

What comes out is housekeeping machinery — V-type ATPase subunit d2, COP9
signalosome subunit 3, dynein axonemal heavy chain 1, GMP reductase 1,
multifunctional protein CAD. 100 species-diverse families across six phyla:
Chordata 84, Arthropoda 12, and one each of Mollusca, Echinodermata, Cnidaria
and Nematoda, from 25 source species. The chordate share reflects Swiss-Prot's
curation, not the sampling frame; the *families* are bilaterian-wide by
construction, whichever species the sequence was read from.

The window is random within the protein, so a fragment is a piece of a conserved
protein rather than necessarily a piece of a conserved *region* — some windows
land in fast-evolving loops. Restricting the cut to annotated Pfam domains is the
change that would tighten this, if the distinction ever matters to a result.

### `init_orfs.faa`

Not yet built. Once the source FASTA arrives:

```bash
python analysis/make_init_sets.py --orfs <source.faa> --pop 100 --len 10-100 --seed 42 -o init/
```

Same cutting logic, one window per ORF, and `--len` should be whatever range the
fragment set was cut to, so the two sets that both come from real sequence are at
least comparable to each other.

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

| Set | Length | Mean charge | Charge range | Hydrophobic |
|---|---|---|---|---|
| random | 25 | +0.22 | −4.0 to +5.2 | 41% |
| fragments | 11–99, median 51 | +0.61 | −7.6 to +14.0 | 40% |

Hydrophobic is the fraction of AVLIMFWC; both rows are computed the same way.
Charge is net charge at pH 7, and is not normalised by length, so the fragments'
wider range is partly the longer sequences and partly the fact that real protein
regions are locally charged in a way independent draws are not. Hydrophobic
content is near identical, so neither set starts with the amphipathic composition
an AMP needs — which is the property the runs are meant to evolve.

The fragments carry natural amino-acid usage, L, A, E, K, S, V most frequent,
against a flat distribution in the random set. Pairwise identity is not quoted
for the fragments because the sequences differ in length; the guarantee that
replaces it is stronger, namely that the 100 fragments come from 100 distinct
UniRef50 families and 100 distinct proteins.
