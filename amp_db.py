"""
amp_db.py — known antimicrobial peptide (AMP) seed database.

Lets PFES start evolution from a *real, validated* AMP instead of a random
sequence or a manually typed one. Use it via the --initial_seq / -iseq flag:

    -iseq dbaasp                 random curated AMP (one seed, whole population)
    -iseq dbaasp:random          same as above
    -iseq dbaasp:diverse         population seeded with many different AMPs
    -iseq dbaasp:low             random AMP annotated low-hemolytic (good for
                                 lead optimisation: minimise hemolysis further)
    -iseq dbaasp:<name>          a specific AMP, e.g. -iseq dbaasp:magainin_2

Sequences are canonical AMPs curated from DBAASP v3 (https://dbaasp.org) and the
primary literature. The curated set below always works offline. fetch_dbaasp()
can additionally pull a larger live set from the DBAASP REST API.

DBAASP citation: Pirtskhalava et al., Nucleic Acids Res. 2021 (DBAASP v3).
"""

import random as _random
import sys


# Curated, high-confidence AMP seeds (standard 20 aa only).
# hemolytic: 'low' | 'moderate' | 'high' | 'unknown'  (literature consensus)
# note: source organism / structural class. 'disulfide' = needs S-S bonds for
#       full activity, so it is a weaker seed for linear-helix optimisation.
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


def _norm(s: str) -> str:
    return s.strip().lower().replace("-", "_").replace(" ", "_")


def list_amps():
    """Return the list of curated AMP names."""
    return [d["name"] for d in CURATED]


def get_amp(key: str, pool=None):
    """Look up one AMP by name (case-insensitive, '-'/' ' == '_'). Returns the
    record dict. Raises KeyError if not found."""
    pool = pool or CURATED
    k = _norm(key)
    for d in pool:
        if _norm(d["name"]) == k:
            return d
    raise KeyError(
        f"AMP '{key}' not in database. Available: {', '.join(list_amps())}"
    )


def random_amp(hemolytic=None, min_len=None, max_len=None,
               exclude_disulfide=False, rng=None, pool=None):
    """Pick one random AMP record, optionally filtered.

    hemolytic: keep only records with this annotation ('low'/'moderate'/'high').
    min_len/max_len: length bounds.
    exclude_disulfide: drop disulfide-dependent scaffolds.
    """
    rng = rng or _random
    pool = pool or CURATED
    cands = []
    for d in pool:
        if hemolytic and d.get("hemolytic") != hemolytic:
            continue
        if min_len and len(d["seq"]) < min_len:
            continue
        if max_len and len(d["seq"]) > max_len:
            continue
        if exclude_disulfide and "disulfide" in d.get("note", ""):
            continue
        cands.append(d)
    if not cands:
        raise ValueError("no AMP matches the requested filters")
    return rng.choice(cands)


def seeds_for_population(spec: str, pop_size: int, rng=None, pool=None):
    """Resolve a DBAASP seed spec into exactly `pop_size` (name, sequence) tuples.

    spec (the part after 'dbaasp:'):
      '', 'random'    -> one random AMP, replicated across the population
      'diverse'       -> pop_size different AMPs (with replacement if needed)
      'low'/'safe'    -> one random low-hemolytic AMP, replicated
      'high'          -> one random high-hemolytic AMP, replicated
      <name>          -> the named AMP, replicated
    """
    rng = rng or _random
    pool = pool or CURATED
    spec = (spec or "random").strip().lower()

    if spec in ("", "random"):
        d = random_amp(rng=rng, pool=pool)
        return [(d["name"], d["seq"])] * pop_size

    if spec == "diverse":
        chosen = []
        bag = list(pool)
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

    # otherwise: a specific AMP name
    d = get_amp(spec, pool=pool)
    return [(d["name"], d["seq"])] * pop_size


def to_fasta(path):
    """Write the curated set to a FASTA file (headers carry the metadata)."""
    with open(path, "w") as fh:
        for d in CURATED:
            fh.write(f">{d['name']} hemolytic={d['hemolytic']} {d['note']}\n{d['seq']}\n")
    return path


def fetch_dbaasp(max_records=500, timeout=30):
    """Best-effort live fetch of AMP sequences from the DBAASP REST API.

    Returns a list of records [{name, seq, hemolytic, note}] on success, or the
    curated set on any failure (network/schema). The DBAASP API
    (https://dbaasp.org/api/v1?query=...&format=json) is queried defensively:
    any object exposing an amino-acid 'sequence' field is harvested.

    Note: the DBAASP API schema can change; this parser is intentionally lenient
    and the curated set is always a safe fallback.
    """
    import json
    from urllib.request import urlopen, Request
    from urllib.parse import urlencode

    base = "https://dbaasp.org/api/v1"
    params = {"query": "search", "format": "json", "limit": str(max_records)}
    url = f"{base}?{urlencode(params)}"
    try:
        req = Request(url, headers={"Accept": "application/json"})
        with urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8", "replace"))
    except Exception as e:
        sys.stderr.write(
            f"  Warning: DBAASP fetch failed ({type(e).__name__}: {e}) — "
            "using curated set\n"
        )
        return list(CURATED)

    out = []

    def _harvest(obj):
        if isinstance(obj, dict):
            seq = obj.get("sequence") or obj.get("seq")
            if isinstance(seq, str) and seq.isalpha() and 5 <= len(seq) <= 100:
                name = str(obj.get("name") or obj.get("id") or f"dbaasp_{len(out)}")
                out.append({"name": _norm(name), "seq": seq.upper(),
                            "hemolytic": "unknown", "note": "DBAASP API"})
            for v in obj.values():
                _harvest(v)
        elif isinstance(obj, list):
            for v in obj:
                _harvest(v)

    _harvest(data)
    return out or list(CURATED)


if __name__ == "__main__":
    print(f"{len(CURATED)} curated AMPs:")
    for d in CURATED:
        print(f"  {d['name']:18s} {len(d['seq']):3d} aa  hemo={d['hemolytic']:8s} {d['seq']}")
