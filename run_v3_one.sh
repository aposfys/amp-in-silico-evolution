#!/usr/bin/env bash
# One run, one stream.  Usage: run_v3_one.sh <init> <name> <geometry...>
set -u
export KMP_BLOCKTIME=0
export OMP_WAIT_POLICY=PASSIVE
export OMP_NUM_THREADS=12
export MKL_NUM_THREADS=12
export TOKENIZERS_PARALLELISM=false
source "$(dirname "$0")/preflight.sh"   # abort if macrel/hemopi2/esm3 are not really there
R=/data/apostolos/pfes
OUT=$R/results/v3
LOG=$R/results/v3.master.log
INIT="$1"; NAME="$2"; shift 2
[ -s "$R/$INIT" ] || { echo "MISSING $R/$INIT" >&2; exit 1; }
if [ -e "$OUT/$NAME" ]; then echo "REFUSING: $OUT/$NAME exists" >&2; exit 1; fi
echo "=== $(date '+%F %T') START $NAME" >> "$LOG"
python "$R/pfes.py" --start file --start-file "$R/$INIT" \
    -ps 100 -ng 600 -sm weak -b 20 --strong_selection_after_n_gen 480 \
    --norepeat --max-tokens-per-batch 4096 --threads 12 "$@" \
    -o "$OUT/$NAME" > "$OUT/$NAME.console.log" 2>&1
echo "=== $(date '+%F %T') END   $NAME (exit $?)" >> "$LOG"
