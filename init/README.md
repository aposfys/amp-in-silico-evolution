# Starting populations — the origin comparison

Three sets of 100 sequences, all cut to 25 residues, differing only in where the
sequences came from. The question they answer is whether antimicrobial activity
emerges as readily from sequence with no evolutionary history, from novel
protein-coding sequence, or from sequence that already encodes something else.

| File | Origin | Status |
|---|---|---|
| `init_random.faa` | completely random, uniform over the 20 residues | ready |
| `init_orfs.faa` | small ORFs from metazoan transcriptomes | **awaiting the source FASTA** |
| `init_fragments.faa` | fragments of conserved metazoan proteins that are not AMPs | ready |

Each is fed to a run as a fixed starting population:

```bash
python pfes.py --start file --start-file init/init_fragments.faa -pop 100 ...
```

## Why one length

All three are cut to 25 aa. Length is otherwise a confound: the structured
objective penalises anything past `-pl0` (30 aa), and MACREL is defined only over
10–100 residues. A set that happened to be longer would score worse for reasons
having nothing to do with its origin, which is precisely the comparison being
made.

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

One random window per protein. Headers keep the provenance:
`frag_001_O00159_Homo_sapiens_803-827` gives accession, organism and the residue
range cut. `init_fragments.tsv` carries the same 100 rows with the source
protein name, phylum, UniRef50 cluster, cluster size and common taxon.

What comes out is housekeeping machinery — V-type ATPase subunit d2, COP9
signalosome subunit 3, dynein axonemal heavy chain 1, GMP reductase 1,
multifunctional protein CAD. 100 species-diverse families across six phyla:
Chordata 84, Arthropoda 12, and one each of Mollusca, Echinodermata, Cnidaria
and Nematoda, from 21 source species. The chordate share reflects Swiss-Prot's
curation, not the sampling frame; the *families* are bilaterian-wide by
construction, whichever species the sequence was read from.

### `init_orfs.faa`

Not yet built. Once the source FASTA arrives:

```bash
python analysis/make_init_sets.py --orfs <source.faa> --pop 100 --len 25 --seed 42 -o init/
```

Same cutting logic, one window per ORF.

## Outstanding: the AMP screen

These sets are meant to start *without* antimicrobial activity, so anything
MACREL already calls an AMP should be dropped. That screen needs MACREL on PATH
and has **not** been applied to the files as committed, because it was not
available where they were generated.

Regenerate on a machine with MACREL before running:

```bash
python analysis/make_init_sets.py --random --uniprot --max-amp-prob 0.5 -o init/
```

Or screen the existing files and remove anything above the threshold. Either way
the same screen must be applied to all three sets, or the comparison is between
sets filtered to different standards.

## Composition

| Set | Mean charge | Range | Hydrophobic | Pairwise identity |
|---|---|---|---|---|
| random | +0.22 | −4.0 to +5.2 | 41% | 5.0% |
| fragments | −0.08 | −8.9 to +8.0 | 39% | 5.8% |

Hydrophobic is the fraction of AVLIMFWC; both columns are computed the same way.
The two sets are matched in length and close in mean charge and hydrophobic
content, so neither starts with the composition an AMP needs. They differ where
real sequence differs from random: the fragments carry natural amino-acid usage
(L, A, E, K, S, V most frequent, against a flat distribution in the random set)
and a far wider charge range, −8.9 to +8.0, because real protein regions are
locally charged in a way independent draws are not. Pairwise identity near 5% in
both confirms that no family is over-represented.
