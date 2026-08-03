# Starting populations — the variable-length cut

The same two sets as `init/`, from the same code, the same sampling frame and the
same seed, cut at a **random length as well as a random position**: uniform over
10–100 aa, the range MACREL is defined for.

```bash
python analysis/make_init_sets.py --uniprot --pop 100 --len 10-100 --seed 42 --no-screen -o init_varlen/
python analysis/make_init_sets.py --random  --pop 100 --len 10-100 --seed 42 --no-screen -o init_varlen/
```

| File | Cut | Observed |
|---|---|---|
| `init_fragments.faa` | uniform 10–100 | 11–99 aa, median 51 |
| `init_random.faa` | uniform 10–100 | 11–100 aa, median 60 |

The window length is drawn *before* the source protein is chosen, so the
distribution is uniform rather than bent towards long windows by which proteins
happened to be long enough to supply them. Both arms are cut over the same range,
so they are comparable to each other but not to anything in `init/`.

Everything else is identical to `init/` and documented there: conservation
defined by UniRef50 common taxon at Bilateria or older, family-first sampling,
one fragment per family, and the keyword, GO and name exclusions. 100 distinct
UniRef50 families, 100 distinct proteins, 25 source species, six phyla.
`init_fragments.tsv` carries the provenance including `fragment_length`.

## What this costs, and why `init/` exists

The structured fitness carries `P_len = 1 / (1 + e^{0.2(L − 30)})` with
`pl0 = 30` (`pfes.py:242`). Over this range that penalty is not a tilt, it is a
cliff:

| L | 15 | 25 | 30 | 40 | 50 | 70 | 100 |
|---|---|---|---|---|---|---|---|
| `P_len` | 0.95 | 0.73 | 0.50 | 0.119 | 0.018 | 0.00034 | 0.0000008 |

Median `P_len` is 0.015 for the fragments and 0.002 for the random set; 49 and 59
of the 100 respectively start below 0.01. The fitness is a product, so half of
each population scores essentially zero at generation zero and is eliminated in
the first round under `--selection_mode strong`. What survives is not the
sequences whose *origin* suits the objective, it is the sequences that happened
to be short, and each arm is left with a few dozen genuinely competitive
individuals rather than 100.

Use `init/` (fixed 25 aa, `P_len` 0.731 for every member) for any run whose
result is meant to be read as a statement about origin. Use this directory when
the requirement is specifically that the cut be random. `--len 15-30` is the
compromise if both matter: random cutting, `P_len` still between 0.95 and 0.50.

## Outstanding: the AMP screen

Built with `--no-screen`, as `init/` is, and for the same reason. See
`init/README.md`.
