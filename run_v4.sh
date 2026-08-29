#!/usr/bin/env bash
# PFES-AMPs, series v4 — the objective ablation. Six runs, serial, on a GPU host.
#
#   ./run_v4.sh main     # three runs on `main`          (structure x MACREL)
#   ./run_v4.sh ctrl     # three runs on `fitness-esm3`  (pLDDT x pTM)
#
# Three origins x two objectives from identical starting populations. What the
# comparison tests, and why one replicate per cell is enough for it, is in the
# top-level README under "The experiment".
#
# Run the `main` arm first. It is the cheap half (~4-5 h a run against 10-15 h
# for the control, whose chains grow), and if the machine is taken back it is
# the half that still yields the primary contrast.
#
# The two arms live on different branches, so check the control out beside this
# one rather than switching back and forth in place:
#
#   git worktree add ../PFES-AMPs-ctrl fitness-esm3
#   ./run_v4.sh main                       # from here
#   cd ../PFES-AMPs-ctrl && ./run_v4.sh ctrl
#
# Both arms must be seeded from the SAME init files. The control worktree's
# init/ is a separate checkout of the same committed bytes, which is fine; what
# is not fine is regenerating either set (see init/README.md).

set -u

ARM="${1:-}"
case "$ARM" in
    main) WANT_BRANCH=main         ; TAG="main" ;;
    ctrl) WANT_BRANCH=fitness-esm3 ; TAG="ctrl" ;;
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

# --------------------------------------------------------------------------- #
# Shared host. Rucker has 32 cores and one GPU, and is not exclusively ours.
# Folding is on the GPU; PSIQUE, MACREL and HemoPI2 stay on the CPU, and thread
# oversubscription there is what convoyed the v2 series (26 s/generation on a
# free machine against 415 s under contention). Leave half the cores.
export KMP_BLOCKTIME=0
export OMP_WAIT_POLICY=PASSIVE
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export TOKENIZERS_PARALLELISM=false

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
  echo "### load at start: $(uptime | sed 's/.*load average/load average/')"
  nvidia-smi --query-gpu=name,memory.used,memory.total,utilization.gpu \
             --format=csv,noheader 2>/dev/null | sed 's/^/### gpu: /'
} >> "$LOG"

run () {           # $1 = init fasta   $2 = run name
    local init="$R/$1" name="$2"
    echo "=== $(date '+%F %T') START $name" >> "$LOG"
    python "$R/pfes.py" --start file --start-file "$init" \
        $COMMON -o "$OUT/$name" \
        > "$OUT/$name.console.log" 2>&1
    local rc=$?
    echo "=== $(date '+%F %T') END   $name (exit $rc)" >> "$LOG"
    # A run that fell back to the surrogate is not a result. Say so at the end
    # of the run rather than at the end of the analysis, six runs later.
    if grep -qi 'macrel not installed\|falling back' "$OUT/$name.console.log"; then
        echo "*** $name FELL BACK TO A SURROGATE — discard it ***" | tee -a "$LOG"
    fi
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
