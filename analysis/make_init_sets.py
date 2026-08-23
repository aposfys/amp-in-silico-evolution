#!/usr/bin/env python3
"""
Build the starting populations for the origin comparison.

The question is whether antimicrobial activity evolves equally readily from
sequences of three different origins:

  random     completely random sequences, no evolutionary history
  orfs       small ORFs from metazoan transcriptomes, novel protein-coding
  fragments  pieces of conserved metazoan proteins that already do something else

All three are cut to the same length, because length is otherwise a confound:
the structured objective penalises anything past `-pl0` (30 aa by default) and
MACREL is only defined over 10-100 residues, so a set that happens to be longer
would look worse for reasons that have nothing to do with its origin.

The fragment set is drawn family-first. UniRef50 clusters whose members share a
common ancestor at Bilateria or older are the sampling frame, one protein is
taken from each sampled cluster and one random window is cut from that protein.
Sampling entries directly instead would return whatever Swiss-Prot has curated
most, which is human, and would let one well-studied family contribute a dozen
fragments.

    python analysis/make_init_sets.py --random --uniprot -o init/
    python analysis/make_init_sets.py --orfs nikos_orfs.faa -o init/

Candidates already predicted to be antimicrobial are dropped, since the point is
to watch activity *emerge*. That screen needs MACREL on PATH; without it the
script says so and writes the set unscreened rather than failing.

Author: Apostolos Fysekidis
"""

import argparse
import gzip
import os
import random
import re
import sys
import time
import urllib.parse
import urllib.request

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))

AA = "ACDEFGHIKLMNPQRSTVWY"

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIREF_URL = "https://rest.uniprot.org/uniref/search"
METAZOA = "33208"

# --- what counts as conserved ------------------------------------------------
# A UniRef50 cluster is a set of sequences sharing at least 50% identity, so a
# large cluster is a protein that recurs across genomes rather than one that
# happens to be well studied. Two properties are asked of it.
#
#   size          how many UniProtKB entries fall in the cluster, which is a
#                 count over all of UniProt and so is not limited to the few
#                 species Swiss-Prot curates deeply.
#   common taxon  the last common ancestor of every member. If that ancestor is
#                 Bilateria or older then the protein is present in vertebrates,
#                 insects, nematodes and molluscs alike, which is the operational
#                 meaning of "we find it in many species". A cluster whose common
#                 ancestor is Mammalia or Euteleostomi is excluded however large
#                 it is, because it is conserved only within one lineage.
CONSERVED_TAXA = {
    "Bilateria", "Eumetazoa", "Metazoa", "Holozoa", "Filozoa", "Choanozoa",
    "Opisthokonta", "Obazoa", "Amorphea", "Eukaryota", "cellular organisms",
}

# --- what disqualifies a protein --------------------------------------------
# Keywords and GO terms are asked of UniProt, which is cheaper than filtering
# 110,000 entries locally. GO terms match their descendants, so the four below
# cover the whole defence and immunity subtrees.
EXCLUDE_KEYWORDS = [
    "KW-0929",   # Antimicrobial
    "KW-0391",   # Immunity
    "KW-0399",   # Innate immunity
    "KW-1064",   # Adaptive immunity
    "KW-0211",   # Defensin
    "KW-0044",   # Antibiotic
    "KW-0930",   # Antiviral protein
    "KW-0051",   # Antiviral defense
    "KW-0295",   # Fungicide
    "KW-0081",   # Bacteriolytic enzyme
    "KW-0800",   # Toxin
    "KW-0204",   # Cytolysis
    "KW-0354",   # Hemolysis
]
EXCLUDE_GO = [
    "0006952",   # defense response
    "0002376",   # immune system process
    "0031640",   # killing of cells of another organism
    "0019730",   # antimicrobial humoral response
]

# Keywords and GO terms miss things, so names are screened too. The list is the
# families that are known sources of antimicrobial peptides even when the parent
# protein is annotated as something else entirely: histones (buforin, parasin),
# haemoglobin and haemocyanin (their proteolytic fragments), lactoferrin
# (lactoferricin), chromogranin, thrombin and the ribosomal proteins.
NAME_EXCLUDE = re.compile(
    r"antimicrob|antibact|antifung|antivir|antibiot|bacteriocin|bactericid|"
    r"microbicid|defensin|cathelicidin|cecropin|magainin|bombinin|dermcidin|"
    r"buforin|penaeidin|crustin|drosomycin|attacin|diptericin|moricin|"
    r"apidaecin|abaecin|hepcidin|histatin|thionin|granulysin|perforin|"
    r"NK-lysin|saposin|amoebapore|lysozyme|lactoferr|transferrin|histone|"
    r"h(a?e)moglobin|h(a?e)mocyanin|chromogranin|secretogranin|thrombin|"
    r"casein|interferon|interleukin|chemokine|cytokine|complement|"
    r"immunoglobulin|histocompatibility|toll-like|toll-interacting|"
    r"tumor necrosis factor|ribosomal protein|ubiquitin-40S|ubiquitin-60S|"
    r"immune|defen[cs]e|lectin|peptidoglycan recognition|"
    r"lipopolysaccharide-binding|bactericidal",
    re.I,
)

PHYLUM_RE = re.compile(r"([A-Za-z ]+) \(phylum\)")


# --------------------------------------------------------------------------- #
def parse_fasta(text):
    out = []
    for block in text.split(">")[1:]:
        head, *rest = block.split("\n")
        seq = "".join(rest).strip().upper()
        if seq:
            out.append((head.split()[0], seq))
    return out


def write_fasta(path, records):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")
    if not records:
        print(f"  wrote {path} (empty)")
        return
    lens = [len(s) for _, s in records]
    span = (f"{min(lens)} aa" if min(lens) == max(lens)
            else f"{min(lens)}-{max(lens)} aa, mean {sum(lens) / len(lens):.0f}")
    print(f"  wrote {path}  ({len(records)} sequences, {span})")


def clean(seq):
    """Standard 20 residues only; ambiguity codes break the folder and MACREL."""
    return "".join(c for c in seq if c in AA)


def parse_length(spec):
    """`25` for one fixed length, `10-100` for a uniform range."""
    if "-" in str(spec):
        lo, hi = str(spec).split("-", 1)
        lo, hi = int(lo), int(hi)
        if lo < 1 or hi < lo:
            raise ValueError(f"bad length range {spec}")
        return lo, hi
    return int(spec), int(spec)


def cut(seq, length, rng):
    """A random window of `length`, or None if the sequence is too short.

    `length` is (lo, hi); the window length is drawn first and uniformly, so the
    set's length distribution is the one asked for rather than one bent by which
    proteins happened to be long enough.
    """
    lo, hi = length
    seq = clean(seq)
    n = rng.randint(lo, hi)
    if len(seq) < n:
        return None
    start = rng.randint(0, len(seq) - n)
    return seq[start:start + n]


# --------------------------------------------------------------------------- #
def screen_amp(records, max_prob):
    """
    Drop anything MACREL already calls antimicrobial.

    Two properties this screen must have, and previously did not:

    Fail CLOSED. An unscored sequence used to default to 0.0, which is <=
    max_prob, so anything MACREL did not return silently PASSED the screen --
    the exact opposite of what a screen is for. The starting population must
    contain no known AMPs, so a sequence whose AMP probability is unknown is
    dropped, not kept, and the count is reported.

    Say which scorer answered. macrel_score_batch substitutes the biophysical
    calculate_samp surrogate whenever MACREL is unavailable or the sequence
    falls outside its 10-100 residue window. A screen run entirely on the
    surrogate is not a MACREL screen, and printing "dropped N with MACREL >
    0.5" in that case is simply false. The v2 production series ran ~60 h that
    way. Now the source is counted and a surrogate-only screen warns loudly.
    """
    try:
        import score
    except Exception as e:
        sys.stderr.write(f"  AMP screen skipped, cannot import score.py ({e})\n")
        return records, None
    seqs = [s for _, s in records]
    try:
        scored = score.macrel_score_batch_src(seqs)
    except Exception as e:
        sys.stderr.write(f"  AMP screen skipped, MACREL failed ({e})\n")
        return records, None
    if not scored:
        sys.stderr.write("  AMP screen skipped, MACREL returned nothing\n")
        return records, None

    kept, unscored = [], 0
    for n, seq in records:
        hit = scored.get(seq)
        if hit is None:
            unscored += 1          # fail closed: unknown is not "not an AMP"
            continue
        if hit[0] <= max_prob:
            kept.append((n, seq))
    dropped = len(records) - len(kept)

    by_macrel = sum(1 for v in scored.values() if v[2] == 'macrel')
    src = (f"{by_macrel}/{len(scored)} by MACREL"
           if by_macrel else "biophysical surrogate ONLY")
    print(f"  AMP screen: dropped {dropped}/{len(records)} with prob > {max_prob}"
          f"  [{src}]")
    if unscored:
        sys.stderr.write(f"  AMP screen: {unscored} sequence(s) had no score and "
                         f"were dropped rather than assumed non-AMP\n")
    if not by_macrel:
        sys.stderr.write("  *** AMP screen ran on the biophysical surrogate, not "
                         "MACREL. The init set is NOT MACREL-screened. ***\n")
    return kept, dropped


# --------------------------------------------------------------------------- #
def build_random(n, length, rng):
    out = []
    for i in range(n):
        k = rng.randint(*length)
        out.append((f"rand{k}aa_{i}", "".join(rng.choice(AA) for _ in range(k))))
    return out


def get(url, timeout=120, tries=5):
    """One GET with retries, since UniProt occasionally answers 502."""
    for attempt in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"Accept-Encoding": "identity"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode(), r.headers.get("Link", "")
        except Exception as e:
            if attempt == tries - 1:
                raise
            sys.stderr.write(f"  retry {attempt + 1} ({type(e).__name__}: {e})\n")
            time.sleep(5 * (attempt + 1))


def tsv_rows(text):
    lines = text.rstrip("\n").split("\n")
    head = lines[0].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[1:] if l]


# --------------------------------------------------------------------------- #
def fetch_clusters(min_size, cache, timeout=180):
    """Every UniRef50 cluster with a metazoan member and at least `min_size`
    entries, cached, because the download is 150 pages and the sampling below
    has to be re-runnable without it."""
    if os.path.exists(cache):
        with gzip.open(cache, "rt") as fh:
            rows = tsv_rows(fh.read())
        print(f"  {len(rows)} UniRef50 clusters from cache {cache}")
        return rows

    query = f"(taxonomy_id:{METAZOA}) AND (identity:0.5) AND (count:[{min_size} TO *])"
    url = (f"{UNIREF_URL}?query={urllib.parse.quote(query)}"
           "&fields=id,name,common_taxon,count&format=tsv&size=500")
    print(f"  querying UniRef50 for clusters of {min_size}+ members "
          "(one page of 500 per request)...")
    rows, text_all = [], []
    while url:
        text, link = get(url, timeout=timeout)
        text_all.append(text if not rows else "\n".join(text.split("\n")[1:]))
        page = tsv_rows(text)
        rows += page
        print(f"\r    {len(rows)} clusters", end="", flush=True)
        m = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = m.group(1) if m else None
    print()

    os.makedirs(os.path.dirname(os.path.abspath(cache)) or ".", exist_ok=True)
    with gzip.open(cache, "wt") as fh:
        fh.write("\n".join(t.rstrip("\n") for t in text_all) + "\n")
    print(f"  cached {len(rows)} clusters to {cache}")
    return rows


def cluster_members(cluster_id, timeout=60):
    """Reviewed metazoan entries of one cluster that survive every annotation
    filter, whole sequences only."""
    excl = " ".join(f"NOT (keyword:{k})" for k in EXCLUDE_KEYWORDS)
    excl += " " + " ".join(f"NOT (go:{g})" for g in EXCLUDE_GO)
    query = (f"(uniref_cluster_50:{cluster_id}) AND (reviewed:true) "
             f"AND (taxonomy_id:{METAZOA}) AND (fragment:false) {excl}")
    url = (f"{UNIPROT_URL}?query={urllib.parse.quote(query)}"
           "&fields=accession,protein_name,organism_name,organism_id,lineage,"
           "sequence&format=tsv&size=50")
    text, _ = get(url, timeout=timeout)
    return tsv_rows(text)


def sample_conserved(n, length, rng, min_size, cache, tag="frag",
                     exclude=(), offset=0):
    """One fragment per conserved family.

    Clusters are shuffled and walked in order, so the sample is uniform over
    families rather than over entries. Sampling over entries would count a
    family once for every species Swiss-Prot happens to have curated, which is
    what makes a naive metazoan query return page after page of human protein.
    Within a cluster the species is drawn first and the entry second, for the
    same reason at the next level down.
    """
    clusters = fetch_clusters(min_size, cache)
    deep = [c for c in clusters if c["Common taxon"] in CONSERVED_TAXA
            and c["Cluster ID"] not in exclude]
    print(f"  {len(deep)} of {len(clusters)} clusters have a common ancestor "
          f"at Bilateria or older")
    rng.shuffle(deep)

    out, meta, seen_species, tried, empty = [], [], {}, 0, 0
    for c in deep:
        if len(out) >= n:
            break
        tried += 1
        frag_len = rng.randint(*length)          # drawn before the entry, so the
        members = cluster_members(c["Cluster ID"])   # distribution stays uniform
        members = [m for m in members
                   if not NAME_EXCLUDE.search(m.get("Protein names", ""))
                   and len(clean(m.get("Sequence", ""))) >= frag_len]
        if not members:
            empty += 1
            continue

        by_species = {}
        for m in members:
            by_species.setdefault(m["Organism (ID)"], []).append(m)
        species = rng.choice(sorted(by_species))
        entry = rng.choice(by_species[species])

        seq = clean(entry["Sequence"])
        start = rng.randint(0, len(seq) - frag_len)
        frag = seq[start:start + frag_len]

        organism = entry["Organism"].split(" (")[0].replace(" ", "_")
        phylum = PHYLUM_RE.search(entry.get("Taxonomic lineage", ""))
        name = (f"{tag}_{len(out) + offset:03d}_{entry['Entry']}_{organism}_"
                f"{start + 1}-{start + frag_len}")
        out.append((name, frag))
        meta.append({
            "fragment": name,
            "accession": entry["Entry"],
            "protein": entry.get("Protein names", "").split(" (")[0],
            "organism": entry["Organism"].split(" (")[0],
            "phylum": phylum.group(1).strip() if phylum else "",
            "uniref50": c["Cluster ID"],
            "cluster_size": c["Size"],
            "common_taxon": c["Common taxon"],
            "source_length": len(seq),
            "fragment_length": frag_len,
            "window": f"{start + 1}-{start + frag_len}",
            "sequence": frag,
        })
        seen_species[entry["Organism"].split(" (")[0]] = \
            seen_species.get(entry["Organism"].split(" (")[0], 0) + 1
        print(f"\r    {len(out)}/{n} fragments ({tried} clusters queried)",
              end="", flush=True)

    print(f"\n  {tried} clusters queried, {empty} had no eligible reviewed "
          f"entry, {len(seen_species)} source species")
    return out, meta


def build_fragments(source, n, length, rng, tag, exclude=None, offset=0):
    """
    One random window per protein, so fragments are maximally independent.

    Returns (records, meta). `meta` mirrors what sample_conserved produces for
    the fragments arm, so the ORF arm can write the same kind of provenance
    TSV: without it there is no record of which source ORF each fragment came
    from, and a set that cannot be traced cannot be audited.

    `exclude` is a set of source names already used, and `offset` continues the
    fragment numbering, so a screened set can be topped up without reusing a
    source ORF or colliding with existing fragment names.
    """
    exclude = exclude or set()
    pool = [(nm, s) for nm, s in source
            if len(clean(s)) >= length[0] and nm not in exclude]
    if len(pool) < n:
        sys.stderr.write(
            f"  only {len(pool)} sequences are at least {length[0]} aa; "
            f"windows will be reused across the shortfall\n")
    rng.shuffle(pool)
    out, meta, i = [], [], 0
    while len(out) < n and pool:
        nm, s = pool[i % len(pool)]
        frag = cut(s, length, rng)
        if frag:
            short = nm.split("|")[1] if "|" in nm else nm
            name = f"{tag}_{offset + len(out)}_{short}"
            out.append((name, frag))
            meta.append({"fragment": name,
                         "source": nm,
                         "source_len": len(clean(s)),
                         "fragment_len": len(frag)})
        i += 1
        if i > 20 * n:
            break
    return out, meta


def write_meta(path, meta, kept):
    """Provenance for every fragment that survived the screen."""
    kept = {n for n, _ in kept}
    rows = [m for m in meta if m["fragment"] in kept]
    if not rows:
        return
    cols = list(rows[0])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for m in rows:
            fh.write("\t".join(str(m[c]) for c in cols) + "\n")
    print(f"  wrote {path}  ({len(rows)} rows)")


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--random", action="store_true", help="build the random set")
    ap.add_argument("--uniprot", action="store_true",
                    help="build the conserved-protein fragment set from UniProt")
    ap.add_argument("--orfs", metavar="FASTA",
                    help="build the small-ORF set from a supplied FASTA")
    ap.add_argument("--pop", type=int, default=100)
    ap.add_argument("--len", default="25", dest="length",
                    help="fragment length: one number for a fixed length, or "
                         "MIN-MAX for a uniform random length per sequence "
                         "(e.g. 10-100, the range MACREL is defined over)")
    ap.add_argument("--min-cluster-size", type=int, default=100,
                    help="a UniRef50 cluster must hold at least this many "
                         "entries to count as conserved (default 100)")
    ap.add_argument("--cluster-cache", metavar="TSV.GZ",
                    help="where the UniRef50 cluster list is cached "
                         "(default <outdir>/.cache/)")
    ap.add_argument("--max-amp-prob", type=float, default=0.5,
                    help="drop candidates MACREL scores above this (default 0.5)")
    ap.add_argument("--no-screen", action="store_true",
                    help="skip the MACREL screen entirely")
    ap.add_argument("--seed", type=int, default=42,
                    help="seed for selection and cut positions; the sets are "
                         "reproducible even though the runs are not")
    ap.add_argument("-o", "--outdir", default="init")
    args = ap.parse_args()

    if not (args.random or args.uniprot or args.orfs):
        ap.error("choose at least one of --random, --uniprot, --orfs")
    try:
        args.length = parse_length(args.length)
    except ValueError as e:
        ap.error(str(e))

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    if args.random:
        print("\n[1] random sequences")
        recs = build_random(args.pop, args.length, rng)
        if not args.no_screen:
            recs, _ = screen_amp(recs, args.max_amp_prob)
            while len(recs) < args.pop:          # top up after screening
                extra = build_random(args.pop - len(recs), args.length, rng)
                extra, _ = screen_amp(extra, args.max_amp_prob)
                if not extra:
                    break
                recs += extra
        write_fasta(os.path.join(args.outdir, "init_random.faa"), recs[:args.pop])

    if args.uniprot:
        print("\n[3] fragments of conserved non-AMP metazoan proteins")
        cache = args.cluster_cache or os.path.join(
            args.outdir, ".cache", "uniref50_metazoa_clusters.tsv.gz")
        recs, meta = sample_conserved(args.pop, args.length, rng,
                                      args.min_cluster_size, cache)
        if not args.no_screen:
            kept, dropped = screen_amp(recs, args.max_amp_prob)
            # Top up from clusters not already used, so a screened set is still
            # `pop` fragments from `pop` distinct families.
            while dropped and len(kept) < args.pop:
                used = {m["uniref50"] for m in meta}
                extra, extra_meta = sample_conserved(
                    args.pop - len(kept), args.length, rng,
                    args.min_cluster_size, cache,
                    exclude=used, offset=len(kept))
                if not extra:
                    break
                extra, dropped = screen_amp(extra, args.max_amp_prob)
                kept += extra
                meta += extra_meta
            recs = kept
        write_fasta(os.path.join(args.outdir, "init_fragments.faa"), recs[:args.pop])
        write_meta(os.path.join(args.outdir, "init_fragments.tsv"),
                   meta, recs[:args.pop])

    if args.orfs:
        print("\n[2] small ORFs from metazoan transcriptomes")
        source = parse_fasta(open(args.orfs).read())
        print(f"  read {len(source)} ORFs from {args.orfs}")
        recs, meta = build_fragments(source, args.pop, args.length, rng, "orf")
        if not args.no_screen:
            kept, dropped = screen_amp(recs, args.max_amp_prob)
            # Top up from ORFs not already used, so a screened set is still
            # `pop` fragments from `pop` distinct source ORFs. Without this the
            # ORF arm silently ends up smaller than the random and fragments
            # arms whenever the screen drops anything, and the three origins
            # stop being comparable.
            while dropped and len(kept) < args.pop:
                used = {m["source"] for m in meta}
                extra, extra_meta = build_fragments(
                    source, args.pop - len(kept), args.length, rng, "orf",
                    exclude=used, offset=len(kept))
                if not extra:
                    break
                extra, dropped = screen_amp(extra, args.max_amp_prob)
                kept += extra
                meta += extra_meta
            recs = kept
        if len(recs) < args.pop:
            sys.stderr.write(
                f"  WARNING: ORF set is {len(recs)}/{args.pop}; it is NOT "
                f"size-matched to the random and fragments arms\n")
        write_fasta(os.path.join(args.outdir, "init_orfs.faa"), recs[:args.pop])
        write_meta(os.path.join(args.outdir, "init_orfs.tsv"),
                   meta, recs[:args.pop])

    print()


if __name__ == "__main__":
    main()
