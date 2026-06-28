"""
amp_db.py — known antimicrobial peptide (AMP) seed database.

Lets PFES start evolution from a *real, validated* AMP instead of a random
sequence or a manually typed one. Use it via the --initial_seq / -iseq flag:

    -iseq dbaasp                 random known AMP (one seed, whole population)
    -iseq dbaasp:random          same as above
    -iseq dbaasp:diverse         population seeded with many different AMPs
    -iseq dbaasp:low             random AMP annotated low-hemolytic (good for
                                 lead optimisation: minimise hemolysis further)
    -iseq dbaasp:<name>          a specific AMP, e.g. -iseq dbaasp:magainin_2

Where the sequences come from
-----------------------------
There are two layers:

  1. CURATED (below): 17 canonical, hand-verified AMPs with reliable hemolytic
     annotations. Always available, fully offline. This is the fallback.

  2. A local cache FASTA (default: <repo>/data/dbaasp.faa) holding the *full*
     DBAASP set. When that file exists, seeding draws from it (thousands of
     peptides). Build it once with either:

         python amp_db.py --fetch                  # DBAASP REST API
         python amp_db.py --import dbaasp_export.fasta   # a file you downloaded

     (Use --import if the API is rate-limited/blocked on your network: export
     the peptide list from https://dbaasp.org as FASTA or TSV, then import it.)

DBAASP citation: Pirtskhalava et al., Nucleic Acids Res. 2021 (DBAASP v3),
https://dbaasp.org
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
                    m = re.match(r"hemolytic=(\w+)", tok, re.I)
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


# ---------------------------------------------------------------------------
# Cache building: import a downloaded file, or fetch from the DBAASP API
# ---------------------------------------------------------------------------
def save_fasta(records, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        for d in records:
            fh.write(f">{d['name']} hemolytic={d.get('hemolytic','unknown')} "
                     f"{d.get('note','')}\n{d['seq']}\n")
    return path


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
                    help="fetch the full DBAASP set from the API into the cache")
    ap.add_argument("--import", dest="import_path", metavar="FILE",
                    help="build the cache from a downloaded DBAASP export (FASTA/TSV/CSV)")
    ap.add_argument("--out", default=CACHE_PATH, help=f"cache path (default {CACHE_PATH})")
    ap.add_argument("--id-end", type=int, default=20000, help="max peptide id to scan (--fetch)")
    ap.add_argument("--list", action="store_true", help="list the active pool and exit")
    args = ap.parse_args(argv)

    if args.import_path:
        recs = import_file(args.import_path)
        save_fasta(recs, args.out)
        print(f"imported {len(recs)} AMPs -> {args.out}")
        return
    if args.fetch:
        recs = fetch_dbaasp(id_end=args.id_end)
        if not recs:
            print("fetch returned nothing — try --import with a downloaded export")
            return
        save_fasta(recs, args.out)
        print(f"fetched {len(recs)} AMPs -> {args.out}")
        return

    pool = load_pool()
    print(f"{len(pool)} AMPs in active pool "
          f"({'cache: ' + CACHE_PATH if os.path.isfile(CACHE_PATH) else 'curated'}):")
    for d in (pool if args.list else pool[:17]):
        print(f"  {d['name']:18s} {len(d['seq']):3d} aa  "
              f"hemo={d.get('hemolytic','?'):8s} {d['seq']}")


if __name__ == "__main__":
    _main()
