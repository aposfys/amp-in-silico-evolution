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
    # The executable existing is not the same as it working, and this one has
    # failed in three distinct ways: a torch/torchvision CUDA mismatch aborting
    # the transformers import, numpy past its pin, and models emitting values
    # outside [0,1]. In every case score.py falls back per sequence to
    # calculate_hemo_proxy and warns only on stderr, so a whole run logs a
    # biophysical surrogate in hemo_prob while the tool reports as present.
    #
    # Probe end to end, and check ORDER as well as range: model 4 returns
    # exactly 0.0 for every sequence, which is inside [0,1], perfectly useless,
    # and calls melittin non-hemolytic. See score.hemopi2_agrees.
    _pf_hemo=$(python - <<'PY' 2>/dev/null
import score
ok, why = score.hemopi2_agrees()
print(f"ok|model {score.HEMOPI2_MODEL}, {why}" if ok else f"bad|{why}")
PY
)
    case "$_pf_hemo" in
        ok\|*)   echo "  hemopi2:    ${_pf_hemo#ok|}" ;;
        "")      echo "  hemopi2:    *** CALL FAILED *** $(command -v hemopi2_classification)"
                 _pf_fail=1 ;;
        *)       echo "  hemopi2:    *** INSTALLED BUT NOT WORKING ***"
                 echo "              ${_pf_hemo#bad|}"
                 # If the failure is not one this understands, show what score.py
                 # itself wrote when it fell back, rather than guessing.
                 python -c "
import score
score.hemopi2_score_batch(['GIGKFLHSAKKFGKAFVGEIMNS'])" 2>&1 >/dev/null \
                     | grep -v '^\s*$' | tail -4 | sed 's/^ *//; s/^/              | /'
                 echo "              try another model:  PFES_HEMO_MODEL=1 ./preflight.sh"
                 _pf_fail=1 ;;
    esac
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
elif [ "${#HF_TOKEN}" -lt 20 ] || case "$HF_TOKEN" in *...*) true ;; *) false ;; esac; then
    # A placeholder pasted from documentation sets the variable without setting
    # a token, and every check that only tests for non-empty then passes. Real
    # tokens are hf_ followed by roughly 34 characters.
    echo "  HF_TOKEN:   *** SET TO A PLACEHOLDER *** ('$HF_TOKEN')"
    echo "              that is the example string, not a token. Replace it, and"
    echo "              check ~/.bashrc if you appended it there."
    _pf_fail=1
else
    echo "  HF_TOKEN:   set (${#HF_TOKEN} chars)"
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
