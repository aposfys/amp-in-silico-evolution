#!/usr/bin/env bash
# amp-in-silico-evolution, series v4 — the objective ablation. Six runs, serial, on a GPU host.
#
#   ./run_v4.sh main     # three runs on `main`          (structure x MACREL)
#   ./run_v4.sh ctrl     # three runs on `control-fold-only`  (pLDDT x pTM)
#
# Three origins x two objectives from identical starting populations. What the
# comparison tests, and why one replicate per cell is enough for it, is in the
# top-level README under "The experiment".
#
# The three origins run STRICTLY SEQUENTIALLY, one process at a time -- random,
# then fragments, then orfs -- each starting only after the previous exits.
# They share one GPU and one CPU budget, so overlapping them would divide the
# core cap three ways and contend for VRAM.
#
# Run the `main` arm first. It is the cheap half (~4-5 h a run against 10-15 h
# for the control, whose chains grow), and if the machine is taken back it is
# the half that still yields the primary contrast.
#
# CPU cap defaults to 12 cores of the 32 on rucker; PFES_CORES overrides it.
#
# The two arms live on different branches, so check the control out beside this
# one rather than switching back and forth in place:
#
#   git worktree add ../amps-ctrl control-fold-only
#   ./run_v4.sh main                       # from here
#   cd ../amps-ctrl && ./run_v4.sh ctrl
#
# Both arms must be seeded from the SAME init files. The control worktree's
# init/ is a separate checkout of the same committed bytes, which is fine; what
# is not fine is regenerating either set (see init/README.md).

set -u

ARM="${1:-}"
case "$ARM" in
    main) WANT_BRANCH=main         ; TAG="main" ;;
    ctrl) WANT_BRANCH=control-fold-only ; TAG="ctrl" ;;
    *)    echo "usage: $0 main|ctrl" >&2; exit 2 ;;
esac

R="$(cd "$(dirname "$0")" && pwd)"
cd "$R" || exit 1

# --------------------------------------------------------------------------- #
# Guard: is this checkout actually the objective the arm name claims?
#
# The v2 series ran four production runs with a header printing
# "[MACREL AMP + PFES]" while macrel was not installed. An arm that cannot
# verify its own objective is how that happens. Refuse rather than mislabel.
BRANCH="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
if [ "$BRANCH" != "$WANT_BRANCH" ]; then
    echo "ABORT: arm '$ARM' expects branch '$WANT_BRANCH', this checkout is on '$BRANCH'" >&2
    exit 1
fi
if ! git diff --quiet HEAD -- pfes.py score.py; then
    echo "ABORT: pfes.py or score.py is modified in the working tree." >&2
    echo "       The objective constants are hard-coded, so an uncommitted edit" >&2
    echo "       changes what runs without changing anything a log records." >&2
    exit 1
fi
echo "arm $ARM  branch $BRANCH  commit $(git rev-parse --short HEAD)"

# Guard: the right environment, not merely a working one.
#
# preflight.sh catches an environment with no classifiers. It cannot catch one
# that belongs to a different project and happens to import esm -- and this
# pipeline has now been launched twice from `deeppeptide`, whose onnxruntime is
# not pinned and whose torch is shared with the other half of the thesis.
# Override with PFES_ENV if the environment is named something else.
WANT_ENV="${PFES_ENV:-pfes_amps}"
if [ -z "${CONDA_PREFIX:-}" ]; then
    die_env="no conda environment is active"
elif [ "$(basename "$CONDA_PREFIX")" != "$WANT_ENV" ]; then
    die_env="active environment is '$(basename "$CONDA_PREFIX")', expected '$WANT_ENV'"
fi
if [ -n "${die_env:-}" ]; then
    echo "ABORT: $die_env" >&2
    echo "       conda activate $WANT_ENV      (build it with ./setup_gpu.sh $WANT_ENV)" >&2
    echo "       or set PFES_ENV=<name> if you meant a different one." >&2
    exit 1
fi
echo "env $(basename "$CONDA_PREFIX")  python $(command -v python)"

# --------------------------------------------------------------------------- #
# Shared host. Rucker has 32 cores and one GPU and is not exclusively ours.
# Folding is on the GPU; PSIQUE, MACREL and HemoPI2 stay on the CPU, and thread
# oversubscription there is what convoyed the v2 series (26 s/generation on a
# free machine against 415 s under contention).
#
# CORES is a hard cap, not a hint. Thread-count environment variables only bind
# libraries that read them, and this pipeline spawns processes that do not:
# psique is a subprocess per structure, HemoPI2 is a separate program, and
# onnxruntime picks its own pool. taskset bounds the whole process tree,
# children included, so the number below is the number.
CORES="${PFES_CORES:-12}"
CPUSET="${PFES_CPUSET:-0-$((CORES-1))}"

export KMP_BLOCKTIME=0
export OMP_WAIT_POLICY=PASSIVE     # sleep rather than spin at the OMP barrier
export OMP_NUM_THREADS="$CORES"
export MKL_NUM_THREADS="$CORES"
export OPENBLAS_NUM_THREADS="$CORES"
export NUMEXPR_NUM_THREADS="$CORES"
export TOKENIZERS_PARALLELISM=false

TASKSET=""
if command -v taskset >/dev/null 2>&1; then
    TASKSET="taskset -c $CPUSET"
    echo "cpu cap: $CORES cores, pinned to $CPUSET"
else
    echo "cpu cap: $CORES cores by thread count only — taskset not found, so"
    echo "         subprocesses are NOT bounded. Install util-linux or accept it."
fi

source "$R/preflight.sh" || exit 1   # macrel/hemopi2/esm3 must really be callable

OUT="$R/results/v4/$TAG"
LOG="$R/results/v4/$TAG.master.log"

NGEN=600
STRONG=$(( NGEN * 80 / 100 ))
POP=100
HEMO_EVERY=25            # HemoPI2 is 65 % of a generation for a column that
                         # never enters the fitness; score_posthoc.py --lineage
                         # recovers the gaps. See RUNBOOK.md.
MAXTOK="${MAXTOK:-4096}" # override from ./gpu_inventory.sh, sized on FREE VRAM

COMMON="-ps $POP -ng $NGEN -sm weak -b 20 \
        --strong_selection_after_n_gen $STRONG --norepeat \
        --hemo-every $HEMO_EVERY --max-tokens-per-batch $MAXTOK \
        -pl0 30 -hl0 30 -bl0 12"

# --------------------------------------------------------------------------- #
missing=0
for f in init/init_random.faa init/init_fragments.faa init/init_orfs.faa; do
    if [ ! -s "$R/$f" ]; then
        echo "MISSING: $R/$f" >&2; missing=1
    else
        printf '  ok  %-28s %s sequences  md5 %s\n' "$f" \
            "$(grep -c '^>' "$R/$f")" "$(md5sum "$R/$f" 2>/dev/null | cut -c1-8 || md5 -q "$R/$f" | cut -c1-8)"
    fi
done
[ "$missing" -eq 0 ] || { echo "aborting: init sets are not present" >&2; exit 1; }

mkdir -p "$OUT"
{
  echo "### v4/$TAG start $(date '+%F %T')"
  echo "### branch $BRANCH  commit $(git rev-parse HEAD)"
  echo "### pop $POP, $NGEN generations, strong from $STRONG, hemo every $HEMO_EVERY"
  echo "### max-tokens-per-batch $MAXTOK"
  echo "### cpu cap $CORES cores${TASKSET:+, pinned to $CPUSET}"
  echo "### load at start: $(uptime | sed 's/.*load average/load average/')"
  nvidia-smi --query-gpu=name,memory.used,memory.total,utilization.gpu \
             --format=csv,noheader 2>/dev/null | sed 's/^/### gpu: /'
} >> "$LOG"

# Strictly sequential. Each call blocks until its run exits, so only one
# pfes.py is ever resident: they share one GPU, and the CPU-side work is capped
# at $CORES for the series as a whole rather than per run. Overlapping them
# would divide that cap by three and contend for VRAM -- the v2 series lost a
# run to 80.7 h against 4.6 h for its siblings, purely from concurrency.
N_RUNS=3
run_i=0

run () {           # $1 = init fasta   $2 = run name
    local init="$R/$1" name="$2"
    run_i=$((run_i + 1))
    local t0=$(date +%s)
    printf '\n[%d/%d] %s  started %s\n' "$run_i" "$N_RUNS" "$name" "$(date '+%F %T')"
    echo "=== $(date '+%F %T') START $name ($run_i/$N_RUNS)" >> "$LOG"

    $TASKSET python "$R/pfes.py" --start file --start-file "$init" \
        $COMMON -o "$OUT/$name" \
        > "$OUT/$name.console.log" 2>&1
    local rc=$?

    local el=$(( $(date +%s) - t0 ))
    printf '[%d/%d] %s  %s after %dh%02dm  (exit %d)\n' \
        "$run_i" "$N_RUNS" "$name" \
        "$([ $rc -eq 0 ] && echo finished || echo FAILED)" \
        $((el/3600)) $(((el%3600)/60)) "$rc"
    echo "=== $(date '+%F %T') END   $name (exit $rc, ${el}s)" >> "$LOG"

    # A run that fell back to the surrogate is not a result. Say so at the end
    # of the run rather than at the end of the analysis, six runs later.
    if grep -qi 'macrel not installed\|falling back' "$OUT/$name.console.log"; then
        echo "*** $name FELL BACK TO A SURROGATE — discard it ***" | tee -a "$LOG"
    fi
    # Carry on rather than abandoning the series: the remaining origins are
    # independent, and one failure should not cost a night of the other two.
    [ $rc -eq 0 ] || echo "    continuing to the next origin; $name needs re-running" >&2
}

run init/init_random.faa     random
run init/init_fragments.faa  fragments
run init/init_orfs.faa       orfs

echo "### v4/$TAG done $(date '+%F %T')" >> "$LOG"

cat <<EOF

done. before reading anything:

  grep -c macrel $OUT/*/progress.log        # amp_src must be macrel, not proxy
  python analysis/score_posthoc.py $OUT/*/ -o $OUT/comparison --lineage
EOF
