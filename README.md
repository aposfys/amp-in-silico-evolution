# amp-in-silico-evolution
In silico evolution of antimicrobial peptides.

ESM3 folds each candidate, MACREL scores it for antimicrobial activity, and the PFES
structural terms constrain length and secondary structure. A population of 100 peptides
is mutated and selected for 600 generations. Built on
[PFES](https://github.com/sahakyanhk/pfes) (Sahakyan et al., *PNAS* 2025,
[e2509015122](https://doi.org/10.1073/pnas.2509015122)).

> **You are on `control-fold-only`.** This is `main` with one change: selection is
> driven by ESM3 fold confidence alone (mean pLDDT × pTM), with no structural terms
> and no MACREL. Every other quantity is still computed and logged as an attribute,
> so this branch's `progress.log` has the same schema as `main`'s and the two are
> directly comparable. It is the control arm of the ablation.

### Before you run this
This code accompanies an MSc thesis. **If you intend to run it, please contact me first**
— apostolosfysekidis1@gmail.com. The run outputs and the input populations are not
published here and are available on request. MIT licensed, so this is a request, not a
condition.

### Running
```
./setup_gpu.sh          # build the environment on a CUDA node
./preflight.sh          # must pass before any run
./run_v4.sh             # the reported series
```
`psique` is not shipped here; build it from [upstream](https://github.com/fadasme/psique)
and place it at `bin/psique`. Full procedure and measured costs in [RUNBOOK.md](RUNBOOK.md).

### Documentation
- [OBJECTIVE.md](OBJECTIVE.md) — the fitness function, term by term, every constant traced to its source
- [RESULTS.md](RESULTS.md) — series v4, the objective ablation
- [RUNBOOK.md](RUNBOOK.md) — how to run a series and what it costs
- [VOID-RUNS.md](VOID-RUNS.md) — every series before v4 and what was wrong with each
- [init/README.md](init/README.md) — how the starting populations were built and screened
- [analysis/README.md](analysis/README.md) — the figure and report scripts
- [NOTICE.md](NOTICE.md) — what is forked, what is original, what the pipeline requires

### Branches
- `main` — the full objective
- `control-fold-only` — ESM3 fold confidence alone, the control arm
- `pfes-original` — unmodified upstream PFES, the fork point

### Licence
MIT, © 2026 Apostolos Fysekidis, over this project's contributions only. See
[NOTICE.md](NOTICE.md).
