#!/usr/bin/env bash
# Build the PFES-AMPs environment on a CUDA machine, from nothing.
#
#   ./setup_gpu.sh [env_name]        # default: pfes_amps
#
# Ends by running preflight.sh, which refuses to pass unless MACREL, HemoPI2
# and ESM3 are all genuinely callable. If preflight fails, the environment is
# not ready and no run should be started -- the v2 production series lost ~60 h
# to exactly that, completing with neither classifier installed while its own
# log header claimed otherwise.
#
# NOT TESTED ON CUDA. It was written and syntax-checked on a machine with no
# NVIDIA GPU, so treat the torch install step as the likely failure point and
# read what it prints. Everything is idempotent: re-running is safe.
set -u
ENV_NAME="${1:-pfes_amps}"

say() { printf '\n\033[1m── %s\033[0m\n' "$*"; }
die() { printf '\n*** %s\n' "$*" >&2; exit 1; }

say "prerequisites"
command -v conda >/dev/null 2>&1 || die "conda not on PATH"
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader | sed 's/^/  /'
else
    echo "  WARNING: nvidia-smi not found. Continuing, but this will build a CPU"
    echo "           environment and there is no point running it on this host."
fi

# Pick the torch wheel from the driver's CUDA version. The driver is backward
# compatible, so a cu121 wheel runs on a 12.4 driver, but not the reverse.
CU_TAG=cu121
if command -v nvidia-smi >/dev/null 2>&1; then
    DRV=$(nvidia-smi 2>/dev/null | grep -o 'CUDA Version: [0-9.]*' | head -1 | grep -o '[0-9.]*')
    case "$DRV" in
        12.[6-9]*|13.*) CU_TAG=cu124 ;;
        12.[1-5]*)      CU_TAG=cu121 ;;
        11.*)           CU_TAG=cu118 ;;
    esac
    echo "  driver CUDA ${DRV:-unknown} -> torch wheel $CU_TAG"
fi

say "conda environment: $ENV_NAME"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda env list | grep -qE "^${ENV_NAME}\s" \
    && echo "  exists, reusing" \
    || conda create -n "$ENV_NAME" python=3.11 -y
conda activate "$ENV_NAME" || die "could not activate $ENV_NAME"
echo "  python: $(command -v python)"

say "torch (CUDA build)"
pip install --quiet torch --index-url "https://download.pytorch.org/whl/${CU_TAG}" \
    || die "torch install failed -- try a different CU_TAG, see https://pytorch.org/get-started/locally/"
python - <<'PY'
import torch, sys
print(f"  torch {torch.__version__}  compiled for CUDA {torch.version.cuda}")
if not torch.cuda.is_available():
    print("  *** torch cannot see a GPU. The CPU wheel is installed, or the driver")
    print("      is too old for this wheel. Fix this before continuing. ***")
    sys.exit(1)
p = torch.cuda.get_device_properties(0)
print(f"  sees: {p.name}  {p.total_memory/1e9:.1f} GB")
PY
[ $? -eq 0 ] || die "torch cannot see the GPU"

say "pipeline dependencies"
# requirements.txt pins onnxruntime<=1.25.1: 1.26 changed the shape of ONNX
# output_probability, so MACREL returns raw decision values and magainin 2
# comes back at -0.050, classified NOT an AMP.
pip install --quiet -r requirements.txt || die "requirements.txt install failed"
pip install --quiet hemopi2 || echo "  WARNING: hemopi2 failed; set PFES_SKIP_HEMO=1 or fix before running"

say "MACREL"
conda install -y -c bioconda -c conda-forge macrel >/dev/null 2>&1 \
    && echo "  $(macrel --version 2>&1 | head -1)" \
    || die "macrel install failed"

say "onnxruntime pin"
python - <<'PY'
import onnxruntime as o
from packaging.version import Version
bad = Version(o.__version__) > Version("1.25.1")
print(f"  onnxruntime {o.__version__}" + ("  *** ABOVE THE PIN — MACREL WILL BE WRONG ***" if bad else "  (ok)"))
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
