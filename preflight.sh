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
    # The executable existing is not the same as it working. On rucker it was
    # installed and on PATH, and every invocation raised at import -- transformers
    # pulls in torchvision, and a torchvision built for a different CUDA major
    # version than torch aborts the import. score.py then falls back per sequence
    # to calculate_hemo_proxy and warns only on stderr, so an entire run logs a
    # biophysical surrogate in the hemo_prob column while preflight reports the
    # tool as present. Probe it end to end, the way macrel is probed.
    _pf_hemo=$(python - <<'PY' 2>/dev/null
import score, sys
d = score.hemopi2_score_batch(['GIGKFLHSAKKFGKAFVGEIMNS'])   # magainin-2
v = list(d.values())[0]
# calculate_hemo_proxy is deterministic; if HemoPI2 did not answer, the value
# is exactly what the surrogate returns for this sequence.
print('proxy' if abs(v - score.calculate_hemo_proxy('GIGKFLHSAKKFGKAFVGEIMNS')) < 1e-9
      else f'{v:.4f}')
PY
)
    case "$_pf_hemo" in
        "")      echo "  hemopi2:    *** CALL FAILED *** $(command -v hemopi2_classification)"
                 _pf_fail=1 ;;
        proxy)   echo "  hemopi2:    *** INSTALLED BUT NOT WORKING ***"
                 echo "              it returned the biophysical surrogate, not HemoPI2."
                 # Name the cause rather than guessing at it. Both known failures
                 # are import-time and both are version skew in the torch stack.
                 python - <<'PY' 2>/dev/null | sed 's/^/              /'
from packaging.version import Version
try:
    import numpy
    if Version(numpy.__version__) >= Version("2"):
        print(f"cause: numpy {numpy.__version__} against the numpy<2 pin —")
        print("       HemoPI2's pickled scikit-learn models cannot load.")
        print("       fix:   pip install 'numpy<2'")
except Exception:
    pass
try:
    import torch, torchvision
    if torch.version.cuda.split(".")[0] != torchvision.version.cuda.split(".")[0]:
        print(f"cause: torch CUDA {torch.version.cuda} against torchvision "
              f"{torchvision.version.cuda} — the transformers import aborts.")
        print("       fix:   pip install --force-reinstall --no-deps torchvision \\")
        print("                  --index-url https://download.pytorch.org/whl/cu128")
except Exception:
    pass
PY
                 echo "              full error:  hemopi2_classification --help"
                 _pf_fail=1 ;;
        *)       echo "  hemopi2:    $(command -v hemopi2_classification)  (magainin-2 -> $_pf_hemo)" ;;
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
