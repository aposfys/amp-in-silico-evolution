"""
amp_db.py — known antimicrobial peptide (AMP) seed database.

Lets PFES start evolution from a *real, validated* AMP instead of a random
sequence or a manually typed one. Use it via the --initial_seq / -iseq flag:

    -iseq db                     random known AMP (one seed, whole population)
    -iseq db:random              same as above
    -iseq db:diverse             population seeded with many different AMPs
    -iseq db:low                 random AMP annotated low hemolysis-risk (good for
                                 lead optimisation: minimise hemolysis further)
    -iseq db:<name>              a specific AMP, e.g. -iseq db:magainin_2

    ("dbaasp" is accepted as an alias for "db" for backward compatibility.)

Where the sequences come from
-----------------------------
There are two layers:

  1. CURATED (below): 17 canonical, hand-verified AMPs with reliable hemolytic
     annotations. Always available, fully offline. This is the fallback.

  2. A local cache FASTA (default: <repo>/data/dbaasp.faa) holding the *full*
     DBAASP set. When that file exists, seeding draws from it (thousands of
     peptides). Build it once, e.g.:

         python amp_db.py --fetch                  # DRAMP (default; ~11k peptides,
                                                   #   one static download, no API gate)
         python amp_db.py --fetch --source dbaasp  # DBAASP REST API (often gated)
         python amp_db.py --import export.fasta     # any FASTA/TSV/CSV you downloaded

     DRAMP is the default because its bulk set is a single reliable file download,
     unlike DBAASP's rate-limited/User-Agent-gated API. Use --import for a DBAASP
     website export (FASTA or TSV) if you specifically want DBAASP records.

Citations: DRAMP — Shi et al., Nucleic Acids Res. 2022 (http://dramp.cpu-bioinfor.org);
DBAASP — Pirtskhalava et al., Nucleic Acids Res. 2021 (https://dbaasp.org).
"""

import os
import re
import sys
import random as _random


# ---------------------------------------------------------------------------
# Curated, high-confidence AMP seeds (standard 20 aa only).
# hemolytic: 'low' | 'moderate' | 'high' | 'unknown'  (literature consensus)
# note: source organism / structural class. 'disulfide' = needs S-S bonds for
#       full activity, so it is a weaker seed for linear-helix optimisation.
# ---------------------------------------------------------------------------
CURATED = [
    {"name": "magainin_2",      "seq": "GIGKFLHSAKKFGKAFVGEIMNS",                 "hemolytic": "low",      "note": "Xenopus laevis; helical"},
    {"name": "pexiganan",       "seq": "GIGKFLKKAKKFGKAFVKILKK",                  "hemolytic": "low",      "note": "MSI-78, synthetic magainin analog; helical"},
    {"name": "melittin",        "seq": "GIGAVLKVLTTGLPALISWIKRKRQQ",             "hemolytic": "high",     "note": "Apis mellifera venom; helical (toxic reference)"},
    {"name": "ll_37",           "seq": "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES",  "hemolytic": "moderate", "note": "human cathelicidin; helical"},
    {"name": "cecropin_a",      "seq": "KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK",  "hemolytic": "low",      "note": "Hyalophora cecropia; helical"},
    {"name": "indolicidin",     "seq": "ILPWKWPWWPWRR",                           "hemolytic": "moderate", "note": "bovine; extended/Trp-rich"},
    {"name": "aurein_1_2",      "seq": "GLFDIIKKIAESF",                           "hemolytic": "low",      "note": "Litoria aurea; helical"},
    {"name": "buforin_2",       "seq": "TRSSRAGLQFPVGRVHRLLRK",                   "hemolytic": "low",      "note": "toad; helical, DNA-binding"},
    {"name": "temporin_a",      "seq": "FLPLIGRVLSGIL",                           "hemolytic": "moderate", "note": "Rana temporaria; helical"},
    {"name": "mastoparan",      "seq": "INLKALAALAKKIL",                          "hemolytic": "high",     "note": "wasp venom; helical (toxic reference)"},
    {"name": "piscidin_1",      "seq": "FFHHIFRGIVHVGKTIHRLVTG",                  "hemolytic": "high",     "note": "fish; helical"},
    {"name": "polybia_mp1",     "seq": "IDWKKLLDAAKQIL",                          "hemolytic": "low",      "note": "wasp; helical"},
    {"name": "camel",           "seq": "KWKLFKKIGAVLKVL",                         "hemolytic": "moderate", "note": "cecropin-A/melittin hybrid; helical"},
    {"name": "histatin_5",      "seq": "DSHAKRHHGYKRKFHEKHHSHRGY",               "hemolytic": "low",      "note": "human saliva; His-rich"},
    {"name": "lactoferricin_b", "seq": "FKCRRWQWRMKKLGAPSITCVRRAF",              "hemolytic": "low",      "note": "bovine lactoferrin; loop"},
    {"name": "protegrin_1",     "seq": "RGGRLCYCRRRFCVCVGR",                      "hemolytic": "high",     "note": "porcine; beta-hairpin, disulfide"},
    {"name": "tachyplesin_1",   "seq": "KWCFRVCYRGICYRRCR",                       "hemolytic": "high",     "note": "horseshoe crab; beta-hairpin, disulfide"},
]

# Default local cache: <repo>/data/dbaasp.faa (override with $PFES_AMP_DB).
CACHE_PATH = os.environ.get(
    "PFES_AMP_DB",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "dbaasp.faa"),
)

# DRAMP general-AMP set: a single static FASTA (~11k peptides), no API gate.
DRAMP_GENERAL_FASTA = ("http://dramp.cpu-bioinfor.org/downloads/"
                       "download_data/DRAMP3.0_new/general_amps.fasta")

_AA = set("ACDEFGHIKLMNPQRSTVWY")
_pool_cache = None  # memoised loaded pool


def _norm(s):
    return s.strip().lower().replace("-", "_").replace(" ", "_")


def _clean_seq(seq):
    """Uppercase, keep only standard 20 aa; return '' if anything else remains."""
    s = re.sub(r"\s+", "", seq).upper()
    return s if s and all(c in _AA for c in s) else ""


# ---------------------------------------------------------------------------
# Pool loading (cache file if present, else curated)
# ---------------------------------------------------------------------------
def parse_fasta(path):
    """Parse a FASTA file into records. Header tokens like 'hemolytic=low' are
    captured; the first whitespace-delimited token is the name."""
    records, name, hemo, note, buf = [], None, "unknown", "", []

    def _flush():
        if name and buf:
            seq = _clean_seq("".join(buf))
            if seq and 5 <= len(seq) <= 100:
                records.append({"name": name, "seq": seq,
                                "hemolytic": hemo, "note": note})

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                _flush()
                header = line[1:].strip()
                parts = header.split()
                name = _norm(parts[0]) if parts else f"amp_{len(records)}"
                hemo, note_tokens = "unknown", []
                for tok in parts[1:]:
                    m = re.match(r"(?:hemo_risk|hemolytic)=(\w+)", tok, re.I)
                    if m:
                        hemo = m.group(1).lower()
                    else:
                        note_tokens.append(tok)
                note = " ".join(note_tokens)
                buf = []
            else:
                buf.append(line)
    _flush()
    return records


def load_pool(force_curated=False):
    """Return the active AMP pool: the cache FASTA if it exists and is usable,
    otherwise the curated 17. Result is memoised."""
    global _pool_cache
    if force_curated:
        return list(CURATED)
    if _pool_cache is not None:
        return _pool_cache
    if os.path.isfile(CACHE_PATH) and os.path.getsize(CACHE_PATH) > 0:
        try:
            recs = parse_fasta(CACHE_PATH)
            if recs:
                # merge curated hemolytic labels onto matching cache entries
                cur = {_norm(d["name"]): d for d in CURATED}
                for r in recs:
                    if r["hemolytic"] == "unknown" and _norm(r["name"]) in cur:
                        r["hemolytic"] = cur[_norm(r["name"])]["hemolytic"]
                _pool_cache = recs
                sys.stderr.write(f"  amp_db: loaded {len(recs)} AMPs from cache "
                                 f"{CACHE_PATH}\n")
                return _pool_cache
        except Exception as e:
            sys.stderr.write(f"  amp_db: cache unreadable ({e}); using curated\n")
    _pool_cache = list(CURATED)
    return _pool_cache


# ---------------------------------------------------------------------------
# Accessors
# ---------------------------------------------------------------------------
def list_amps(pool=None):
    return [d["name"] for d in (pool or load_pool())]


def get_amp(key, pool=None):
    pool = pool or load_pool()
    k = _norm(key)
    for d in pool:
        if _norm(d["name"]) == k:
            return d
    # fall back to curated by name even if a big cache is active
    for d in CURATED:
        if _norm(d["name"]) == k:
            return d
    raise KeyError(f"AMP '{key}' not in database "
                   f"({len(pool)} loaded). Try one of: "
                   f"{', '.join(list_amps(CURATED))}")


def random_amp(hemolytic=None, min_len=None, max_len=None,
               exclude_disulfide=False, rng=None, pool=None):
    rng = rng or _random
    pool = pool or load_pool()

    def _filter(src):
        out = []
        for d in src:
            if hemolytic and d.get("hemolytic") != hemolytic:
                continue
            if min_len and len(d["seq"]) < min_len:
                continue
            if max_len and len(d["seq"]) > max_len:
                continue
            if exclude_disulfide and "disulfide" in d.get("note", ""):
                continue
            out.append(d)
        return out

    cands = _filter(pool)
    if not cands and hemolytic:
        # cache entries are mostly hemolytic='unknown'; honour the filter using
        # the curated set, which has reliable labels
        cands = _filter(CURATED)
        if cands:
            sys.stderr.write(f"  amp_db: no '{hemolytic}' label in cache — "
                             f"drawing from curated set\n")
    if not cands:
        raise ValueError("no AMP matches the requested filters")
    return rng.choice(cands)


def seeds_for_population(spec, pop_size, rng=None, pool=None):
    """Resolve a DBAASP seed spec into exactly `pop_size` (name, sequence) tuples.

    spec (the part after 'dbaasp:'):
      '', 'random'    -> one random AMP, replicated across the population
      'diverse'       -> pop_size different AMPs (with replacement if needed)
      'low'/'safe'    -> one random low-hemolytic AMP, replicated
      'high'          -> one random high-hemolytic AMP, replicated
      <name>          -> the named AMP, replicated
    """
    rng = rng or _random
    pool = pool or load_pool()
    spec = (spec or "random").strip().lower()

    if spec in ("", "random"):
        d = random_amp(rng=rng, pool=pool)
        return [(d["name"], d["seq"])] * pop_size

    if spec == "diverse":
        chosen, bag = [], list(pool)
        rng.shuffle(bag)
        while len(chosen) < pop_size:
            if not bag:
                bag = list(pool)
                rng.shuffle(bag)
            d = bag.pop()
            chosen.append((d["name"], d["seq"]))
        return chosen[:pop_size]

    if spec in ("low", "safe", "low_hemo"):
        d = random_amp(hemolytic="low", rng=rng, pool=pool)
        return [(d["name"], d["seq"])] * pop_size

    if spec in ("high", "high_hemo"):
        d = random_amp(hemolytic="high", rng=rng, pool=pool)
        return [(d["name"], d["seq"])] * pop_size

    d = get_amp(spec, pool=pool)
    return [(d["name"], d["seq"])] * pop_size


def mixed_population(pop_size, frac_existing=0.5, randomseq=None,
                     default_len=24, rng=None, pool=None):
    """Build a mixed starting population: `frac_existing` of pop_size drawn from
    the known-AMP database (diverse), the remainder random sequences whose length
    equals the existing seeds' MEAN length.

    `randomseq` is a callable(length) -> sequence (e.g. Evolver.randomseq); if
    None, a uniform 20-aa generator is used. Returns a list of (name, sequence)
    of length pop_size; random entries are named 'rand<len>aa'.
    """
    rng = rng or _random
    pool = pool or load_pool()
    n_exist = round(pop_size * frac_existing)
    n_rand = pop_size - n_exist
    existing = seeds_for_population('diverse', n_exist, rng=rng, pool=pool) if n_exist > 0 else []
    if existing:
        mean_len = max(1, round(sum(len(s) for _, s in existing) / len(existing)))
    else:
        mean_len = default_len
    if randomseq is None:
        randomseq = lambda n: "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n))
    randoms = [(f"rand{mean_len}aa", randomseq(mean_len)) for _ in range(n_rand)]
    return existing + randoms


# ---------------------------------------------------------------------------
# Cache building: import a downloaded file, or fetch from the DBAASP API
# ---------------------------------------------------------------------------
def save_fasta(records, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        for d in records:
            fh.write(f">{d['name']} hemo_risk={d.get('hemolytic','unknown')} "
                     f"{d.get('note','')}\n{d['seq']}\n")
    return path


# Hemolysis-risk tiers derived from a predicted hemolytic probability in [0, 1].
HEMO_LOW_MAX = 0.30    # p < 0.30  -> low risk
HEMO_HIGH_MIN = 0.60   # p > 0.60  -> high risk; in between -> moderate


def risk_tier(prob):
    if prob is None:
        return "unknown"
    if prob < HEMO_LOW_MAX:
        return "low"
    if prob > HEMO_HIGH_MIN:
        return "high"
    return "moderate"


def annotate_hemo(records, batch=2000):
    """Compute the hemolysis-risk attribute for each record by scoring its
    sequence with HemoPI2 (via score.hemopi2_score_batch, which falls back to the
    biophysical proxy if HemoPI2 is not installed). Sets record['hemolytic'] to
    'low'/'moderate'/'high' and record['hemo_prob'] to the raw probability.
    Returns the records (and prints a tier distribution)."""
    from score import hemopi2_score_batch  # lazy import (pulls numpy/pandas)
    seqs = [r["seq"] for r in records]
    scores = {}
    for i in range(0, len(seqs), batch):
        chunk = seqs[i:i + batch]
        scores.update(hemopi2_score_batch(chunk))
        sys.stderr.write(f"  amp_db: hemolysis-scored {min(i + batch, len(seqs))}"
                         f"/{len(seqs)}\n")
    dist = {"low": 0, "moderate": 0, "high": 0, "unknown": 0}
    for r in records:
        p = scores.get(r["seq"])
        r["hemo_prob"] = p
        r["hemolytic"] = risk_tier(p)
        dist[r["hemolytic"]] += 1
    sys.stderr.write(f"  amp_db: hemo_risk distribution {dist}\n")
    return records


def import_file(path):
    """Build records from a user-downloaded DBAASP export (FASTA or TSV/CSV).

    FASTA: parsed directly. TSV/CSV: any column named like 'sequence' is used;
    a 'name'/'id' column and a 'hemolytic' column are used if present.
    """
    if path.lower().endswith((".fa", ".faa", ".fasta")):
        return parse_fasta(path)

    import csv
    recs = []
    with open(path, newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        reader = csv.DictReader(fh, delimiter=delim)
        cols = {c.lower(): c for c in (reader.fieldnames or [])}
        seq_c = next((cols[c] for c in cols if "seq" in c), None)
        name_c = next((cols[c] for c in cols if c in ("name", "id", "peptide_id")), None)
        hemo_c = next((cols[c] for c in cols if "hemolytic" in c or "hemo" in c), None)
        if not seq_c:
            raise ValueError(f"no sequence column in {path} (cols: {reader.fieldnames})")
        for i, row in enumerate(reader):
            seq = _clean_seq(row.get(seq_c, ""))
            if not seq or not (5 <= len(seq) <= 100):
                continue
            name = _norm(str(row.get(name_c) or f"dbaasp_{i}")) if name_c else f"dbaasp_{i}"
            hemo = "unknown"
            if hemo_c and row.get(hemo_c):
                val = str(row[hemo_c]).strip()
                if re.search(r"non|^no$|negative|false|0", val, re.I):
                    hemo = "low"
                elif re.search(r"hemo|toxic|positive|true|yes|1", val, re.I):
                    hemo = "high"
            recs.append({"name": name, "seq": seq, "hemolytic": hemo, "note": "DBAASP import"})
    return recs


def fetch_dramp(url=DRAMP_GENERAL_FASTA, timeout=120):
    """Download the DRAMP general-AMP set (one static FASTA, ~11k peptides) and
    return parsed records. This is the recommended bulk source: a single reliable
    file download with no gated API. DRAMP: Shi et al., Nucleic Acids Res. 2022.
    """
    import tempfile
    from urllib.request import urlopen, Request
    headers = {"User-Agent": "Mozilla/5.0",
               "Referer": "http://dramp.cpu-bioinfor.org/downloads/"}
    try:
        with urlopen(Request(url, headers=headers), timeout=timeout) as resp:
            data = resp.read().decode("utf-8", "replace")
    except Exception as e:
        sys.stderr.write(f"  amp_db: DRAMP fetch failed ({type(e).__name__}: {e})\n")
        return []
    if not data.lstrip().startswith(">"):
        sys.stderr.write("  amp_db: DRAMP response was not FASTA (download URL may "
                         "have changed) — try --import with a manual download\n")
        return []
    tmp = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)
    try:
        tmp.write(data)
        tmp.close()
        recs = parse_fasta(tmp.name)
    finally:
        os.unlink(tmp.name)
    for r in recs:
        r["note"] = "DRAMP general"
    return recs


def fetch_dbaasp(id_start=1, id_end=20000, delay=0.05, timeout=30, max_fail=200):
    """Best-effort live fetch from the DBAASP REST API by iterating peptide cards.

    Returns a list of records. The DBAASP API can be rate-limited or require a
    browser User-Agent; if it fails wholesale, prefer downloading an export from
    https://dbaasp.org and using import_file()/--import instead.
    """
    import json
    import time
    from urllib.request import urlopen, Request

    ua = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
          "AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.0 Safari/605.1.15")
    base = "https://dbaasp.org/api/v1"
    out, fails = [], 0

    def _harvest(obj, pid):
        seq = None
        if isinstance(obj, dict):
            seq = obj.get("sequence") or obj.get("seq")
            if not (isinstance(seq, str) and _clean_seq(seq)):
                for v in obj.values():
                    _harvest(v, pid)
                return
        elif isinstance(obj, list):
            for v in obj:
                _harvest(v, pid)
            return
        cs = _clean_seq(seq) if isinstance(seq, str) else ""
        if cs and 5 <= len(cs) <= 100:
            out.append({"name": f"dbaasp_{pid}", "seq": cs,
                        "hemolytic": "unknown", "note": "DBAASP API"})

    for pid in range(id_start, id_end + 1):
        url = f"{base}?query=peptide_card&peptide_id={pid}&format=json"
        try:
            req = Request(url, headers={"User-Agent": ua, "Accept": "application/json"})
            with urlopen(req, timeout=timeout) as resp:
                data = json.loads(resp.read().decode("utf-8", "replace"))
            before = len(out)
            _harvest(data, pid)
            if len(out) == before:
                fails += 1
            else:
                fails = 0
        except Exception:
            fails += 1
        if fails >= max_fail:
            sys.stderr.write(f"  amp_db: stopping after {max_fail} consecutive "
                             f"misses near id {pid}\n")
            break
        if delay:
            time.sleep(delay)
        if pid % 500 == 0:
            sys.stderr.write(f"  amp_db: fetched {len(out)} sequences "
                             f"(scanned to id {pid})\n")
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="AMP seed database utilities")
    ap.add_argument("--fetch", action="store_true",
                    help="download a bulk AMP set into the cache (see --source)")
    ap.add_argument("--source", choices=["dramp", "dbaasp"], default="dramp",
                    help="bulk source for --fetch (default: dramp, a reliable static download)")
    ap.add_argument("--import", dest="import_path", metavar="FILE",
                    help="build the cache from a downloaded export (FASTA/TSV/CSV)")
    ap.add_argument("--annotate-hemo", action="store_true",
                    help="(re)compute the hemo_risk attribute for the cache and save it")
    ap.add_argument("--no-annotate", action="store_true",
                    help="skip hemolysis annotation when fetching/importing")
    ap.add_argument("--out", default=CACHE_PATH, help=f"cache path (default {CACHE_PATH})")
    ap.add_argument("--id-end", type=int, default=20000, help="max peptide id to scan (--fetch --source dbaasp)")
    ap.add_argument("--list", action="store_true", help="list the active pool and exit")
    args = ap.parse_args(argv)

    def _maybe_annotate(recs):
        if not args.no_annotate:
            sys.stderr.write(f"  amp_db: scoring hemolysis risk for {len(recs)} AMPs"
                             " (HemoPI2 if installed, else biophysical proxy)…\n")
            annotate_hemo(recs)
        return recs

    if args.annotate_hemo:
        src = args.out if os.path.isfile(args.out) else CACHE_PATH
        if not os.path.isfile(src):
            print(f"no cache at {src} — run --fetch or --import first")
            return
        recs = annotate_hemo(parse_fasta(src))
        save_fasta(recs, args.out)
        print(f"annotated {len(recs)} AMPs with hemo_risk -> {args.out}")
        return
    if args.import_path:
        recs = _maybe_annotate(import_file(args.import_path))
        save_fasta(recs, args.out)
        print(f"imported {len(recs)} AMPs -> {args.out}")
        return
    if args.fetch:
        recs = fetch_dramp() if args.source == "dramp" else fetch_dbaasp(id_end=args.id_end)
        if not recs:
            print(f"fetch from {args.source} returned nothing — "
                  "try --import with a downloaded export")
            return
        _maybe_annotate(recs)
        save_fasta(recs, args.out)
        print(f"fetched {len(recs)} AMPs from {args.source} -> {args.out}")
        return

    pool = load_pool()
    print(f"{len(pool)} AMPs in active pool "
          f"({'cache: ' + CACHE_PATH if os.path.isfile(CACHE_PATH) else 'curated'}):")
    for d in (pool if args.list else pool[:17]):
        print(f"  {d['name']:18s} {len(d['seq']):3d} aa  "
              f"hemo={d.get('hemolytic','?'):8s} {d['seq']}")


if __name__ == "__main__":
    _main()
