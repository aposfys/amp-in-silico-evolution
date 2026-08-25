#!/usr/bin/env bash
# Read-only inventory of a GPU node. Changes nothing, safe to run anywhere.
# Answers the three questions that decide how to configure a run:
#   what GPU and how much VRAM, is torch actually using it, and what is missing.
#
#   ./gpu_inventory.sh

export KMP_DUPLICATE_LIB_OK=TRUE      # keep OpenMP warnings out of the report

echo "════ host ════"
echo "  $(uname -srm)   $(hostname)"
echo "  cores: $(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo '?')   RAM: $(free -g 2>/dev/null | awk '/^Mem:/{print $2" GB"}' || echo '?')"

echo
echo "════ GPU ════"
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=index,name,memory.total,memory.used,driver_version,compute_cap \
               --format=csv,noheader 2>/dev/null | sed 's/^/  /'
    echo "  CUDA (driver): $(nvidia-smi 2>/dev/null | grep -o 'CUDA Version: [0-9.]*' | head -1)"
else
    echo "  nvidia-smi NOT FOUND — no NVIDIA GPU, or drivers not installed"
fi

echo
echo "════ torch ════"
python - <<'PY' 2>/dev/null | sed 's/^/  /'
try:
    import torch
    print(f"torch {torch.__version__}")
    print(f"compiled for CUDA: {torch.version.cuda}")
    print(f"cuda.is_available(): {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        for i in range(torch.cuda.device_count()):
            p = torch.cuda.get_device_properties(i)
            print(f"  [{i}] {p.name}  {p.total_memory/1e9:.1f} GB  sm_{p.major}{p.minor}")
    else:
        print("*** torch cannot see a GPU — the CPU wheel is probably installed ***")
        print("    fix: pip install torch --index-url https://download.pytorch.org/whl/cu121")
except ImportError:
    print("torch NOT INSTALLED")
PY

echo
echo "════ pipeline tools ════"
for t in macrel hemopi2_classification; do
    printf "  %-24s " "$t"
    command -v $t >/dev/null 2>&1 && echo "$(command -v $t)" || echo "MISSING"
done
printf "  %-24s " "psique (bundled)"
if [ -x bin/psique ]; then
    if ./bin/psique --version >/dev/null 2>&1; then echo "runs"; else echo "present, but will not exec on this arch"; fi
else echo "MISSING"; fi
printf "  %-24s " "esm3"
if python -c "from esm.models.esm3 import ESM3" >/dev/null 2>&1; then echo ok; else echo "IMPORT FAILED"; fi
printf "  %-24s " "onnxruntime"
_ort=$(python -c "import onnxruntime as o; print(o.__version__)" 2>/dev/null)
if [ -n "$_ort" ]; then
    echo "$_ort $(python - <<'PY' 2>/dev/null
from packaging.version import Version
import onnxruntime as o
print("  <-- ABOVE 1.25.1: MACREL will return raw decision values, not probabilities" if Version(o.__version__) > Version("1.25.1") else "(within the <=1.25.1 pin)")
PY
)"
else echo "MISSING"; fi
printf "  %-24s " "HF_TOKEN"
[ -n "$HF_TOKEN" ] && echo set || echo "NOT SET — ESM3 weights cannot be fetched"

echo
echo "════ suggested --max-tokens-per-batch ════"
python - <<'PY' 2>/dev/null | sed 's/^/  /'
try:
    import torch
    if torch.cuda.is_available():
        gb = torch.cuda.get_device_properties(0).total_memory/1e9
        # ESM3-open 1.4B is ~5.6 GB in fp32; the rest is activations, which
        # scale with tokens in the batch. Conservative: ~1024 tokens per free GB.
        rec = int(max(2048, (gb - 7) * 1024))
        print(f"{gb:.0f} GB VRAM -> try --max-tokens-per-batch {rec}")
        print(f"(current default is 512; the CPU runs used 4096)")
        print("Raise it until you OOM, then halve. A whole pop-400 generation of")
        print("30-mers is ~12,000 tokens, so one batch per generation is the goal.")
    else:
        print("no GPU visible — cannot advise; 4096 is the tested CPU value")
except ImportError:
    print("torch not installed")
PY
