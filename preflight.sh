#!/usr/bin/env bash
# Abort a run before it starts if the scoring tools are not actually callable.
#
# Why this exists: the v2 production series (four runs, ~60 h) completed with
# neither MACREL nor HemoPI2 installed. score.py falls back to the biophysical
# calculate_samp / calculate_hemo_proxy surrogates per sequence and writes a
# warning to stderr, so the runs looked healthy, their headers still printed
# "MACREL AMP + PFES", and every amp_prob in progress.log was the proxy. The
# cause was `conda activate pfes_amps` resolving to whichever env of that name
# came first on the HPC -- there are two, and only one has macrel.
#
# A fitness function that silently swaps identity is worse than one that
# crashes. Source this from every launcher, after activating the environment.
#
#   source "$(dirname "$0")/preflight.sh"
#
# Set PFES_SKIP_HEMO=1 to make the HemoPI2 check advisory instead of fatal.

_pf_fail=0

echo "--- preflight ---"
echo "  python:     $(command -v python || echo MISSING)"
echo "  conda env:  ${CONDA_PREFIX:-<none>}"

if command -v macrel >/dev/null 2>&1; then
    echo "  macrel:     $(command -v macrel)  ($(macrel --version 2>&1 | head -1))"
else
    echo "  macrel:     *** NOT FOUND ***  conda install -c bioconda macrel"
    _pf_fail=1
fi

if command -v hemopi2_classification >/dev/null 2>&1; then
    echo "  hemopi2:    $(command -v hemopi2_classification)"
elif [ "${PFES_SKIP_HEMO:-0}" = "1" ]; then
    echo "  hemopi2:    not found, but PFES_SKIP_HEMO=1 — hemo_prob will be 0.0"
else
    echo "  hemopi2:    *** NOT FOUND ***  pip install hemopi2"
    echo "              (or set PFES_SKIP_HEMO=1 to run without the attribute)"
    _pf_fail=1
fi

if _pf_esm=$(python -c "from esm.models.esm3 import ESM3" 2>&1); then
    echo "  esm3:       ok"
else
    # Print the actual error. "check the env and HF_TOKEN" is not a diagnosis,
    # and the import can fail for reasons that have nothing to do with either --
    # a missing transitive dependency, a numpy or torch version conflict, or
    # `esm` resolving to a different PyPI package of the same name.
    echo "  esm3:       *** IMPORT FAILED ***"
    printf '%s\n' "$_pf_esm" | tail -5 | sed 's/^/              | /'
    echo "              full traceback:  python -c 'from esm.models.esm3 import ESM3'"
    _pf_fail=1
fi

# HF_TOKEN is needed to FETCH the weights, not to import the package. Report it
# separately so a token problem is never mistaken for an install problem.
if [ -z "${HF_TOKEN:-}" ]; then
    if [ -d "${HF_HOME:-$HOME/.cache/huggingface}/hub" ] && \
       find "${HF_HOME:-$HOME/.cache/huggingface}/hub" -maxdepth 1 -name '*esm3*' 2>/dev/null | grep -q .; then
        echo "  HF_TOKEN:   not set, but esm3 weights are already cached — ok"
    else
        echo "  HF_TOKEN:   *** NOT SET *** and no cached esm3 weights"
        echo "              accept the licence at huggingface.co/EvolutionaryScale/esm3-sm-open-v1"
        echo "              then: export HF_TOKEN=hf_...   (add it to ~/.bashrc)"
        _pf_fail=1
    fi
else
    echo "  HF_TOKEN:   set"
fi

# End-to-end: does macrel actually return a probability for a known AMP?
# `command -v macrel` succeeding is not the same as macrel working.
if [ "$_pf_fail" = "0" ]; then
    _pf_probe=$(python - <<'PY' 2>/dev/null
import score
d = score.macrel_score_batch(['GIGKFLHSAKKFGKAFVGEIMNS'])   # magainin-2
print(round(list(d.values())[0][0], 4))
PY
)
    if [ -z "$_pf_probe" ]; then
        echo "  probe:      *** macrel_score_batch failed ***"
        _pf_fail=1
    else
        echo "  probe:      magainin-2 -> amp_prob $_pf_probe"
        # calculate_samp('GIGKFLHSAKKFGKAFVGEIMNS') == 0.8148. If the probe
        # returns exactly that, MACREL did not run and the proxy answered.
        if [ "$_pf_probe" = "0.8148" ]; then
            echo "              *** that is the biophysical proxy value, not MACREL ***"
            _pf_fail=1
        fi
    fi
fi

# Which MACREL path will this run use? The in-process path loads the model once
# instead of once per generation -- worth ~60 s of an 88.6 s generation at
# pop 100 -- but it calls library internals rather than the documented CLI, so
# it is used only after reproducing the subprocess on a known peptide. Report
# the answer here rather than leaving it to be inferred from the runtime.
if [ "$_pf_fail" = "0" ]; then
    _pf_path=$(python - <<'PY' 2>/dev/null
import score
print("in-process (fast)" if score.macrel_inproc_agrees() else "subprocess (per-generation model reload)")
PY
)
    echo "  macrel path: ${_pf_path:-unknown}"
fi

if [ "$_pf_fail" != "0" ]; then
    echo "--- preflight FAILED — refusing to start ---"
    exit 1
fi
echo "--- preflight ok ---"
