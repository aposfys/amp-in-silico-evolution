# amp-in-silico-evolution

**In silico evolution of antimicrobial peptides.** ESM3 folds each candidate,
MACREL scores it for antimicrobial activity, and the PFES structural terms
constrain length and secondary structure. HemoPI2 records hemolysis and AMPlify
re-scores survivors, both as attributes that never touch selection, so they
report on the search rather than restate it. A population of one hundred
peptides is mutated and selected for six hundred generations.

The reported work is a controlled ablation: **three starting populations × two
objectives**, testing whether the AMP character of the winners is caused by the
objective or is simply what ESM3 fold confidence rewards on its own.

Built on [PFES](https://github.com/sahakyanhk/pfes) by Sahakyan et al.

> **You are on `control-fold-only`.** This branch is `main` with one change:
> selection is driven by ESM3 fold confidence alone (mean pLDDT × pTM), with no
> structural terms and no MACREL. Every other quantity is still computed and
> logged as an attribute, so this branch's `progress.log` has the same schema as
> `main`'s and the two are directly comparable. It is the control arm of the
> ablation. For the experiment itself, see `main`.

## Before you run this

This code accompanies an MSc thesis. **If you intend to run it, please contact me
first** — apostolosfysekidis1@gmail.com. I would like to know who is using it.

The run outputs and the exact input populations behind the reported results are
**not published here**. They are available from me on request. Without them you
can read and adapt the method, but you will not reproduce the figures in the
thesis.

This repository is MIT licensed, so the licence does not oblige you to make
contact. The above is a request, not a condition.

## Where to look

| | |
|---|---|
| [`OBJECTIVE.md`](OBJECTIVE.md) | The fitness function, term by term, with every constant traced to its source and the selection pressure each term actually applies |
| [`RESULTS.md`](RESULTS.md) | Series v4, the objective ablation: six runs, what they found |
| [`RUNBOOK.md`](RUNBOOK.md) | How to run a series and what it costs, measured on the production host |
| [`VOID-RUNS.md`](VOID-RUNS.md) | Every series before v4, what was wrong with each, and why they were discarded |
| [`init/README.md`](init/README.md) | How the three starting populations were built and screened |
| [`analysis/README.md`](analysis/README.md) | The figure and report scripts |
| [`NOTICE.md`](NOTICE.md) | What is forked, what is original, and every component the pipeline requires |

## The experiment

Six runs: three starting populations (random sequence, fragments of conserved
metazoan proteins, small ORFs) crossed with two objectives (the full fitness,
and fold confidence alone). Population 100, 600 generations, strong selection
from generation 480. Identical seeds across arms.

The full objective holds chains at 26 residues with `amp_prob` near 0.99 and pTM
near 0.56. Fold confidence alone lets chains grow to 57–65 residues and reaches
pTM near 0.84 while `amp_prob` collapses below 0.12. The origin of the starting
population changes almost nothing. Numbers and interpretation in
[`RESULTS.md`](RESULTS.md).

Runs are not individually reproducible: neither this fork nor upstream PFES
seeds any RNG. The starting population is fixed by file and recorded by checksum
in every run banner. The evidence this design admits is replication across
independent runs.

## Install and run

`psique` is **not shipped here** — build it from
[its repository](https://github.com/fadasme/psique) and place it at
`bin/psique`. `requirements.txt` pins `onnxruntime<=1.25.1`; 1.26 changed the
shape of ONNX `output_probability` and MACREL silently returns raw decision
values instead of calibrated probabilities.

```bash
./setup_gpu.sh          # builds the environment on a CUDA node
./preflight.sh          # must pass before any run
./run_v4.sh             # the reported series
```

`preflight.sh` aborts before a run starts if the scoring tools are not answering
correctly. Do not skip it. Full procedure and measured costs in
[`RUNBOOK.md`](RUNBOOK.md).

## Licence

MIT, © 2026 Apostolos Fysekidis, over this project's contributions only. This is
a fork of [PFES](https://github.com/sahakyanhk/pfes), which is public domain
under the Unlicense. ESM3 weights are gated and governed by the terms you accept
on HuggingFace. See [`NOTICE.md`](NOTICE.md) for the full breakdown.
