# PFES-AMPs — `fitness-esm3`

**The control arm.** Selection is driven by nothing but the two numbers ESM3
emits. No structural penalties, no secondary structure, no activity classifier.

```
score = mean_pLDDT × pTM
```

Paired against [`main`](../../tree/main), the production arm, which is
this plus the full PFES structural objective and MACREL, it isolates what all of
that machinery actually contributes. Run it with the same `--start-file` as the
production branch and the difference is attributable to the objective alone.

## What is logged but does not select

Everything else is still computed and written to `progress.log`, so this branch's
log has the **same schema** as the production branch's and the two are directly
comparable:

| Attribute | Source | Note |
|---|---|---|
| `prot_len_penalty`, `max_alpha_penalty`, `max_beta_penalty` | PFES / PSIQUE | computed, not in score |
| `num_conts` | Eq. 5 contacts (Cβ, 6 Å, \|i−j\| > 5) | computed, not in score |
| `amp_prob` | MACREL | **computed, not in score** |
| `hemo_prob` | HemoPI2 | computed, not in score |

Keeping `amp_prob` is the point of the arm. The question it answers — *what does
fold confidence alone select for, and are those peptides antimicrobial?* — can
only be read off if MACREL scores peptides it never drove. Set
`PFES_SKIP_HEMO=1` to drop the per-generation HemoPI2 call if you don't want the
hemolysis attribute at that price.

## Known behaviour: chains grow

There is no length penalty here, and fold confidence rises with chain length, so
the population drifts well past the 30 aa AMP window. Two consequences:

- **`amp_prob` becomes unreliable above 100 aa.** MACREL is defined for 10–100
  residues; outside that window `macrel_score_batch` substitutes the biophysical
  `calculate_samp` proxy per sequence. Selection is unaffected (MACREL is not in
  the score), but the *attribute you are analysing* silently changes identity.
  The run warns on stderr when this happens — read those warnings before
  interpreting `amp_prob` from a long run.
- **Folding cost rises with length**, so this arm is slower per candidate than
  the production branch despite the simpler score.

Cap length with `-pl0` if you want the comparison held at a fixed size, but note
that reintroduces a term this arm exists to exclude.

## Run

```bash
python pfes.py \
  --start file --start-file init/init_random.faa \
  -ps 100 -ng 600 -sm weak -b 32 \
  --strong_selection_after_n_gen 480 \
  --norepeat --max-tokens-per-batch 4096 \
  -o results/esm3-random-r1 2>&1 | tee results/esm3-random-r1.console.log
```

**Use `-b 32`, not `-b 20`.** Boltzmann selection depends on `β(sᵢ − s_max)`, so
the selective pressure scales with the spread of the score. This arm's scores
run lower than the production branch's (no contact term, no penalties), and β
must rise to keep the pressure comparable. Matching β across arms without
matching score scale makes the comparison meaningless.

## Everything else

Install, starting populations, analysis and the post-hoc AMPlify audit are
identical to the production branch — see
[`main`](../../tree/main).
