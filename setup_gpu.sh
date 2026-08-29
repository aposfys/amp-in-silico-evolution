#!/usr/bin/env bash
# Build the amp-in-silico-evolution environment on a CUDA machine, from nothing.
#
#   ./setup_gpu.sh [env_name]        # default: pfes_amps
#
# Ends by running preflight.sh, which refuses to pass unless MACREL, HemoPI2
# and ESM3 are all genuinely callable. If preflight fails, the environment is
# not ready and no run should be started -- the v2 production series lost ~60 h
# to exactly that, completing with neither classifier installed while its own
# log header claimed otherwise.
#
# Verified on an RTX 5090 (sm_120, driver CUDA 13.2) with torch 2.11.0+cu128:
# kernel launch ok, macrel 1.6.1 from bioconda, onnxruntime pinned at 1.25.1,
# preflight green on the magainin-2 probe. Everything is idempotent, so
# re-running is safe.
set -u
ENV_NAME="${1:-pfes_amps}"

say() { printf '\n\033[1m── %s\033[0m\n' "$*"; }
die() { printf '\n*** %s\n' "$*" >&2; exit 1; }

say "prerequisites"

# `conda` is normally a SHELL FUNCTION installed by `conda init`, and shell
# functions are not inherited by a child process, so `command -v conda` fails
# inside this script even when conda works perfectly in the shell that launched
# it. Locate conda.sh and source it instead of testing for the command.
# Override with CONDA_ROOT=/path/to/conda if the search below misses.
find_conda_sh() {
    local p
    if [ -n "${CONDA_ROOT:-}" ] && [ -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]; then
        echo "$CONDA_ROOT/etc/profile.d/conda.sh"; return 0
    fi
    # CONDA_EXE is exported by conda init and survives into a child process.
    if [ -n "${CONDA_EXE:-}" ] && [ -x "${CONDA_EXE}" ]; then
        p="$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"
        [ -f "$p" ] && { echo "$p"; return 0; }
    fi
    # An active environment gives the root two levels up from <root>/envs/<name>.
    if [ -n "${CONDA_PREFIX:-}" ]; then
        for p in "$CONDA_PREFIX/etc/profile.d/conda.sh" \
                 "$CONDA_PREFIX/../../etc/profile.d/conda.sh"; do
            [ -f "$p" ] && { echo "$p"; return 0; }
        done
    fi
    for p in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3" \
             "$HOME/mambaforge" "$HOME/micromamba" /opt/conda \
             /usr/local/miniconda3 /usr/local/anaconda3 /opt/miniforge3; do
        [ -f "$p/etc/profile.d/conda.sh" ] && { echo "$p/etc/profile.d/conda.sh"; return 0; }
    done
    return 1
}

CONDA_SH="$(find_conda_sh)" || {
    echo "  could not locate conda.sh. Tried CONDA_ROOT, CONDA_EXE, CONDA_PREFIX," >&2
    echo "  and the usual install prefixes. In the shell where conda works, run:" >&2
    echo "      echo \"\$CONDA_EXE\"" >&2
    echo "  then re-run this script as:" >&2
    echo "      CONDA_ROOT=/path/to/miniconda3 ./setup_gpu.sh" >&2
    die "conda not found"
}
echo "  conda: $CONDA_SH"
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader | sed 's/^/  /'
else
    echo "  WARNING: nvidia-smi not found. Continuing, but this will build a CPU"
    echo "           environment and there is no point running it on this host."
fi

# Pick the torch wheel from the GPU's COMPUTE CAPABILITY first, then the
# driver. Compute capability is the binding constraint: a wheel without kernels
# for the architecture cannot drive the card no matter how new the driver is.
# Blackwell (sm_120, e.g. RTX 5090) needs cu128 or later -- cu124 has no
# Blackwell kernels and fails at launch, which the driver version alone does not
# reveal, since a Blackwell card ships with a CUDA 13 driver that would
# otherwise select cu124.
CU_TAG=cu121
if command -v nvidia-smi >/dev/null 2>&1; then
    CC=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d ' ')
    DRV=$(nvidia-smi 2>/dev/null | grep -o 'CUDA Version: [0-9.]*' | head -1 | grep -o '[0-9.]*')
    CC_MAJOR=${CC%%.*}
    if [ -n "$CC_MAJOR" ] && [ "$CC_MAJOR" -ge 12 ] 2>/dev/null; then
        CU_TAG=cu128                       # Blackwell and later
    else
        case "$DRV" in
            12.[6-9]*|13.*) CU_TAG=cu124 ;;
            12.[1-5]*)      CU_TAG=cu121 ;;
            11.*)           CU_TAG=cu118 ;;
        esac
    fi
    echo "  compute capability ${CC:-unknown}, driver CUDA ${DRV:-unknown} -> torch wheel $CU_TAG"
fi

say "conda environment: $ENV_NAME"
# shellcheck disable=SC1090
source "$CONDA_SH" || die "could not source $CONDA_SH"
# conda-forge only, --override-channels, and never the Anaconda `defaults`
# channels. Two reasons, and the first one will stop the script dead otherwise:
#
#   Recent conda refuses to solve against repo.anaconda.com/pkgs/main and
#   /pkgs/r until their Terms of Service are accepted, and that acceptance is a
#   licensing decision for whoever runs this, not something a setup script
#   should make on their behalf. Overriding the channels sidesteps it entirely.
#
#   Bioconda requires conda-forge, listed first, and mixing `defaults` into a
#   bioconda solve is a known source of broken environments. macrel comes from
#   bioconda, so the whole environment is built this way for consistency.
conda env list | grep -qE "^${ENV_NAME}[[:space:]]" \
    && echo "  exists, reusing" \
    || conda create -n "$ENV_NAME" python=3.11 -y \
           -c conda-forge --override-channels
conda activate "$ENV_NAME" || die "could not activate $ENV_NAME"
echo "  python: $(command -v python)"

# A run must never inherit packages from whichever environment happened to be
# active when this was launched. The v2 series lost ~60 h to `conda activate
# pfes_amps` resolving to the wrong environment of that name.
case "${CONDA_PREFIX:-}" in
    */"$ENV_NAME") : ;;
    *) die "activated '$CONDA_PREFIX', expected an env named '$ENV_NAME'" ;;
esac

say "torch (CUDA build)"
# Not --quiet: this is a ~3 GB download and the longest step in the script.
# Silence here reads as a hang.
pip install torch --index-url "https://download.pytorch.org/whl/${CU_TAG}" \
    || die "torch install failed -- try a different CU_TAG, see https://pytorch.org/get-started/locally/"
python - <<'PY'
import torch, sys
print(f"  torch {torch.__version__}  compiled for CUDA {torch.version.cuda}")
if not torch.cuda.is_available():
    print("  *** torch cannot see a GPU. The CPU wheel is installed, or the driver")
    print("      is too old for this wheel. Fix this before continuing. ***")
    sys.exit(1)
p = torch.cuda.get_device_properties(0)
print(f"  sees: {p.name}  {p.total_memory/1e9:.1f} GB  sm_{p.major}{p.minor}")
# is_available() can be True while the wheel lacks kernels for this
# architecture, which only shows up at launch. Force an actual kernel.
try:
    (torch.randn(64, 64, device="cuda") @ torch.randn(64, 64, device="cuda")).sum().item()
    print("  kernel launch: ok")
except Exception as e:
    print(f"  *** kernel launch FAILED: {e}")
    print(f"  *** the wheel has no kernels for sm_{p.major}{p.minor}. Install a newer CUDA build.")
    sys.exit(1)
PY
[ $? -eq 0 ] || die "torch cannot see the GPU"

say "pipeline dependencies"
# requirements.txt pins onnxruntime<=1.25.1: 1.26 changed the shape of ONNX
# output_probability, so MACREL returns raw decision values and magainin 2
# comes back at -0.050, classified NOT an AMP.
pip install --quiet -r requirements.txt || die "requirements.txt install failed"
pip install --quiet hemopi2 || echo "  WARNING: hemopi2 failed; set PFES_SKIP_HEMO=1 or fix before running"

say "MACREL"
# conda-forge before bioconda, per bioconda's own channel-order requirement.
if conda install -y -c conda-forge -c bioconda --override-channels macrel >/dev/null 2>&1; then
    echo "  $(macrel --version 2>&1 | head -1)  [bioconda]"
elif pip install --quiet macrel; then
    # Fallback: bioconda does not always carry a build for the newest Python.
    # macrel is on PyPI, and preflight probes it end to end either way.
    echo "  $(macrel --version 2>&1 | head -1)  [pip fallback]"
else
    die "macrel install failed from both bioconda and pip"
fi

say "onnxruntime pin"
# Installing macrel after requirements.txt can pull onnxruntime past the pin.
# 1.26 changed the shape of ONNX output_probability, so MACREL returns raw
# decision values instead of calibrated probabilities and magainin-2 comes back
# at -0.050, classified NOT an AMP -- sending every candidate to the surrogate,
# silently. Repair it here rather than reporting it.
if ! python -c "
import sys, onnxruntime as o
from packaging.version import Version
sys.exit(0 if Version(o.__version__) <= Version('1.25.1') else 1)" 2>/dev/null; then
    echo "  above the pin after the macrel install — reinstalling"
    pip install --quiet 'onnxruntime<=1.25.1' || die "could not pin onnxruntime"
fi
python - <<'PY'
import onnxruntime as o
from packaging.version import Version
bad = Version(o.__version__) > Version("1.25.1")
print(f"  onnxruntime {o.__version__}" + ("  *** STILL ABOVE THE PIN — MACREL WILL BE WRONG ***" if bad else "  (ok)"))
PY

say "ESM3 access"
if [ -z "${HF_TOKEN:-}" ]; then
    echo "  HF_TOKEN not set. ESM3 weights cannot be fetched."
    echo "  Accept the licence at huggingface.co/EvolutionaryScale/esm3-sm-open-v1,"
    echo "  then:  export HF_TOKEN=hf_..."
else
    echo "  HF_TOKEN set"
fi

say "preflight"
./preflight.sh || die "preflight failed — do not start a run"

say "next"
echo "  ./gpu_inventory.sh                     # confirms --max-tokens-per-batch"
echo "  python profile_generation.py -ps 32 -ng 6"
echo "  python pfes.py --start file --start-file init/init_random.faa \\"
echo "      -ps 8 -ng 3 -sm weak -b 20 --norepeat -o /tmp/smoke   # expect: ready [cuda]"
