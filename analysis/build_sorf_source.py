#!/usr/bin/env python3
'''
Build the sORF arm of the origin comparison from SmProt 2.0.

Two modes, because the two length regimes need different treatment.

    --emit pool   a sampling pool for `make_init_sets.py --orfs`, which cuts one
                  window per ORF. Used for the fixed-width arm.

    --emit init   a finished 100-sequence starting population of *whole* ORFs,
                  MACREL-screened here. Used for the variable-length arm.

Why the variable-length arm is not cut. `cut()` draws the window width first and
rejects sources shorter than it, so that the realised length distribution is the
one asked for rather than one bent by which sources happened to be long. That
holds when sources are proteins of several hundred residues. sORFs are 10-100
residues themselves, so acceptance becomes a function of source length, longer
ORFs are over-sampled and the distribution bends toward short windows — the
exact artefact the function exists to prevent. And a window cut out of a
56-residue sORF is a fragment of a small protein, which is the *other* arm's
treatment. For a sORF the biological unit is the whole ORF. The fixed-width arm
is cut regardless because all three sets must be exactly the same length there,
and a 25-residue window of a 56-residue ORF is a genuine window.

Three corrections are applied to the raw SmProt sets in both modes.

**Stop codons.** SmProt writes the terminating `*` into the sequence. Left in,
it is a residue MACREL and ESM3 have never seen.

**Nested ORFs.** SmProt lists one entry per start codon, so a locus appears
several times as progressively truncated N-termini sharing one C-terminus:

    SPRODRE1  MFQHFLTLQQVSRRQISSTVRRHMANKV...VPQKKA
    SPRODRE2            ISSTVRRHMANKV...VPQKKA

Around 60% of in-range human and mouse entries are redundant this way. Sampled
unfiltered, one locus would contribute several near-identical sequences; the
fragments arm guarantees 100 distinct UniRef50 families and 100 distinct
proteins, and an arm without comparable independence is a smaller experiment
wearing the same n. Entries are grouped on their C-terminal k-mer, which is the
invariant under a start-codon change, and only the longest of each group kept.

**Taxonomic composition.** SmProt's per-species depth follows how much ribosome
profiling each organism has received; human and mouse dominate. The fragments
arm's composition follows Swiss-Prot curation. Left alone the arms differ in
taxonomy as well as in origin and the comparison measures both. Sequences are
therefore drawn to the fragments arm's own proportions over the six species
SmProt covers, which account for 60 of its 100 sequences.

SmProt's high-confidence set is recorded as an annotation, not applied as a
filter: it is 95% human, and using it would trade the taxonomic match for a
translation-evidence tier.

Usage:
    python build_sorf_source.py --emit pool -o init/source_sorfs.faa
    python build_sorf_source.py --emit init --min-len 10 \
           --max-amp-prob 0.5 -o init_varlen/init_orfs.faa
'''
import argparse
import gzip
import os
import random
import sys
from collections import defaultdict

STANDARD = set('ACDEFGHIKLMNPQRSTVWY')

# SmProt file stem -> (binomial, count in init_fragments.tsv)
# Mus 20, Homo 13, Rattus 9, Drosophila 9, Danio 7, Caenorhabditis 2, of the
# 100 fragments; these six are the species SmProt and the fragments arm share.
SPECIES = {
    'mouse':     ('Mus_musculus',            20),
    'human':     ('Homo_sapiens',            13),
    'rat':       ('Rattus_norvegicus',        9),
    'fruitfly':  ('Drosophila_melanogaster',  9),
    'zebrafish': ('Danio_rerio',              7),
    'Celegans':  ('Caenorhabditis_elegans',   2),
}


def read_fasta_gz(path):
    name, buf = None, []
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name:
                    yield name, ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name:
        yield name, ''.join(buf)


def clean(seq):
    seq = seq.strip().upper().rstrip('*')
    if not seq or set(seq) - STANDARD:
        return None
    return seq


def collapse_nested(records, k):
    '''Keep the longest entry of each group sharing a C-terminal k-mer.

    Grouping on the C-terminus is what makes this linear rather than all-pairs:
    SmProt's redundancy is always a shared stop reached from a different start,
    so the tail is the invariant. Sequences shorter than k key on themselves.
    '''
    best = {}
    for sid, seq in records:
        key = seq[-k:] if len(seq) >= k else seq
        if key not in best or len(seq) > len(best[key][1]):
            best[key] = (sid, seq)
    return list(best.values())


def load_species(smprot_dir, min_len, max_len, ctail, quiet=False):
    kept = {}
    for stem, (binomial, _) in SPECIES.items():
        path = os.path.join(smprot_dir, f'SmProt2_{stem}_Ribo.fa.gz')
        if not os.path.isfile(path):
            sys.exit(f'missing: {path}')
        raw, bad, oor = [], 0, 0
        for sid, seq in read_fasta_gz(path):
            c = clean(seq)
            if c is None:
                bad += 1
            elif min_len <= len(c) <= max_len:
                raw.append((sid, c))
            else:
                oor += 1
        collapsed = collapse_nested(raw, ctail)
        kept[stem] = collapsed
        if not quiet:
            print(f'{binomial:26s} {len(collapsed):6d} kept  '
                  f'({len(raw)} in range, {len(raw) - len(collapsed)} nested, '
                  f'{bad} non-standard, {oor} out of length range)')
    return kept


def allocate(kept, total):
    '''Quotas at the fragments arm's proportions; a species short of its quota
    gives up the shortfall to the others rather than shrinking the draw.'''
    w = {s: SPECIES[s][1] for s in SPECIES}
    tw = sum(w.values())
    quota, short = {}, 0
    for s in SPECIES:
        want, have = round(total * w[s] / tw), len(kept[s])
        quota[s] = min(want, have)
        if have < want:
            short += want - have
            print(f'  note: {SPECIES[s][0]} has {have}, quota was {want}')
    if short:
        donors = [s for s in SPECIES if len(kept[s]) > quota[s]]
        dw = sum(w[s] for s in donors) or 1
        for s in donors:
            quota[s] += min(round(short * w[s] / dw), len(kept[s]) - quota[s])
    return quota


def draw(kept, quota, rng, seen):
    rows = []
    for s in SPECIES:
        binomial = SPECIES[s][0]
        for sid, seq in rng.sample(kept[s], quota[s]):
            if seq in seen:          # identical peptides across species
                continue
            seen.add(seq)
            rows.append((sid, binomial, seq))
    rng.shuffle(rows)
    return rows


class LengthMatcher:
    '''Draw ORFs whose lengths reproduce a reference set's distribution.

    The variable-length fragments are random windows of uniform width over
    10-100, so their length distribution is uniform by construction. sORFs are
    not: taken as they come they are markedly shorter, median 34 against 55.
    Since length is the strongest single predictor of antimicrobial activity and
    MACREL is itself length-sensitive, an arm with a different length profile
    would let length stand in for origin, which is the one thing this comparison
    must not do.

    So the reference lengths are read off the arm being matched, shuffled,
    dealt out to the species in their quota proportions, and each one filled by
    the nearest unused ORF of that species. Nearest rather than exact because
    the rarer long lengths are not present in every species, and a hard
    requirement would silently distort the species proportions instead.
    '''

    def __init__(self, kept, rng):
        self.rng = rng
        self.by_species = {}
        for s, recs in kept.items():
            order = sorted(recs, key=lambda r: len(r[1]))
            self.by_species[s] = order
        self.used = {s: set() for s in kept}

    def take(self, species, target):
        pool = self.by_species[species]
        used = self.used[species]
        best, best_d = None, None
        for i, (sid, seq) in enumerate(pool):
            if i in used:
                continue
            d = abs(len(seq) - target)
            if best_d is None or d < best_d:
                best, best_d = i, d
            elif len(seq) - target > best_d:
                break            # sorted by length: nothing further is closer
        if best is None:
            return None
        used.add(best)
        return pool[best]


def match_lengths(kept, quota, targets, rng, seen):
    '''Fill each species' quota against the shuffled reference lengths.'''
    matcher = LengthMatcher(kept, rng)
    rng.shuffle(targets)
    rows, pos = [], 0
    for s in SPECIES:
        binomial = SPECIES[s][0]
        for _ in range(quota[s]):
            if pos >= len(targets):
                break
            got = matcher.take(s, targets[pos])
            pos += 1
            if got is None or got[1] in seen:
                continue
            seen.add(got[1])
            rows.append((got[0], binomial, got[1]))
    rng.shuffle(rows)
    return rows, matcher


def read_lengths(path):
    lens, cur = [], 0
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            if cur:
                lens.append(cur)
            cur = 0
        else:
            cur += len(line)
    if cur:
        lens.append(cur)
    return lens


def screen(rows, max_prob, pfes_dir):
    '''Drop anything MACREL already calls antimicrobial, using the same
    function the run itself calls so the threshold means the same thing.

    Returns (kept, rejected) rather than a count, because a length-matched draw
    has to refill the exact slot it lost.
    '''
    sys.path.insert(0, pfes_dir)
    from score import macrel_score_batch
    sc = macrel_score_batch([r[2] for r in rows])
    keep, drop = [], []
    for r in rows:
        (drop if sc.get(r[2], (0.0, 0.0))[0] > max_prob else keep).append(r)
    return keep, drop


def write_out(rows, path, tag, hc):
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'w') as fh:
        for i, (sid, binomial, seq) in enumerate(rows):
            fh.write(f'>{tag}_{i:03d}_{sid}_{binomial}\n{seq}\n')
    meta = os.path.splitext(path)[0] + '.tsv'
    with open(meta, 'w') as fh:
        fh.write('name\tsmprot_id\torganism\tlength\thigh_confidence\tsequence\n')
        for i, (sid, binomial, seq) in enumerate(rows):
            fh.write(f'{tag}_{i:03d}_{sid}_{binomial}\t{sid}\t{binomial}\t'
                     f'{len(seq)}\t{"yes" if sid in hc else "no"}\t{seq}\n')
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--emit', choices=('pool', 'init'), default='pool',
                    help='pool: source for make_init_sets.py --orfs. '
                         'init: finished screened population of whole ORFs.')
    ap.add_argument('--smprot-dir', default='smprot_raw')
    ap.add_argument('--high-confidence', default=None,
                    help='SmProt2_highConfidenceSet.txt.gz, annotated not filtered')
    ap.add_argument('--n', type=int, default=None,
                    help='sequences to emit (default 3000 for pool, 100 for init)')
    ap.add_argument('--min-len', type=int, default=25,
                    help='shortest ORF kept; for a pool this must be at least '
                         'the window width, for an init set it is the low end '
                         'of the length distribution being matched')
    ap.add_argument('--max-len', type=int, default=100,
                    help='MACREL is undefined past 100')
    ap.add_argument('--ctail', type=int, default=20)
    ap.add_argument('--match-lengths', metavar='FASTA', default=None,
                    help='--emit init only; reproduce this set\'s length '
                         'distribution, so the arms differ in origin and not '
                         'in length')
    ap.add_argument('--max-amp-prob', type=float, default=0.5,
                    help='--emit init only; the screen make_init_sets applies')
    ap.add_argument('--no-screen', action='store_true')
    ap.add_argument('--pfes-dir', default='/data/apostolos/pfes')
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('-o', '--out', required=True)
    args = ap.parse_args()

    n = args.n or (3000 if args.emit == 'pool' else 100)
    rng = random.Random(args.seed)

    hc = set()
    if args.high_confidence:
        op = gzip.open if args.high_confidence.endswith('.gz') else open
        with op(args.high_confidence, 'rt') as fh:
            next(fh, None)
            for line in fh:
                p = line.split('\t')
                if len(p) > 1:
                    hc.add(p[1].strip())
        print(f'high-confidence annotation: {len(hc)} ids loaded')

    kept = load_species(args.smprot_dir, args.min_len, args.max_len, args.ctail)

    seen = set()

    if args.emit == 'pool':
        rows = draw(kept, allocate(kept, n), rng, seen)

    elif not args.match_lengths:
        target = n if args.no_screen else int(n * 1.3) + 10
        rows = draw(kept, allocate(kept, target), rng, seen)
        if not args.no_screen:
            rows, drop = screen(rows, args.max_amp_prob, args.pfes_dir)
            print(f'  AMP screen: dropped {len(drop)}/{len(drop) + len(rows)}')
        if len(rows) < n:
            sys.exit(f'only {len(rows)} available, needed {n}')
        rows = rows[:n]

    else:
        # Deal the reference lengths out to the species quotas, then refill any
        # slot the screen empties with the next-nearest ORF of the same species
        # and the same target length. Refilling the slot rather than topping up
        # from anywhere is what keeps both distributions intact at once.
        targets = read_lengths(args.match_lengths)
        print(f'matching lengths of {args.match_lengths} '
              f'({len(targets)} sequences)')
        if len(targets) < n:
            targets = (targets * (n // len(targets) + 1))[:n]
        targets = rng.sample(targets, n)

        quota = allocate(kept, n)
        matcher = LengthMatcher(kept, rng)
        slots, pos = [], 0
        for s in SPECIES:
            for _ in range(quota[s]):
                if pos < len(targets):
                    slots.append((s, targets[pos]))
                    pos += 1

        rows, dropped_total = [], 0
        for attempt in range(8):
            cand = []
            for s, t in slots:
                got = matcher.take(s, t)
                if got and got[1] not in seen:
                    seen.add(got[1])
                    cand.append((got[0], SPECIES[s][0], got[1], s, t))
            if not cand:
                break
            if args.no_screen:
                rows += cand
                break
            keep, drop = screen([(c[0], c[1], c[2]) for c in cand],
                                args.max_amp_prob, args.pfes_dir)
            kept_ids = {r[0] for r in keep}
            rows += [c for c in cand if c[0] in kept_ids]
            dropped_total += len(drop)
            slots = [(c[3], c[4]) for c in cand if c[0] not in kept_ids]
            if not slots:
                break
        if not args.no_screen:
            print(f'  AMP screen: dropped {dropped_total} across '
                  f'{attempt + 1} round(s), refilled at the same lengths')
        if len(rows) < n:
            sys.exit(f'only {len(rows)} filled, needed {n}')
        rng.shuffle(rows)
        rows = [(a, b, c) for a, b, c, _, _ in rows][:n]

        ach = sorted(len(r[2]) for r in rows)
        tgt = sorted(targets)
        print(f'  target lengths  min {tgt[0]}, median {tgt[len(tgt)//2]}, '
              f'max {tgt[-1]}')
        print(f'  achieved        min {ach[0]}, median {ach[len(ach)//2]}, '
              f'max {ach[-1]}')

    tag = 'sorf' if args.emit == 'pool' else 'orf'
    meta = write_out(rows, args.out, tag, hc)

    lens = sorted(len(r[2]) for r in rows)
    print(f'\nwrote {args.out}  ({len(rows)} sequences, '
          f'{lens[0]}-{lens[-1]} aa, median {lens[len(lens)//2]})')
    print(f'wrote {meta}')
    if hc:
        print(f'high-confidence: {sum(1 for r in rows if r[0] in hc)} of {len(rows)}')
    comp = defaultdict(int)
    for r in rows:
        comp[r[1]] += 1
    for k, v in sorted(comp.items(), key=lambda x: -x[1]):
        print(f'  {k:26s} {v:5d}  ({100 * v / len(rows):.1f}%)')


if __name__ == '__main__':
    main()
