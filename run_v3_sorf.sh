#!/usr/bin/env bash
# The sORF stream: the short replicate first, then both long ones.
# Chained in one script rather than gated on a pgrep, because a pgrep pattern
# naming the run matches the waiting shell's own command line and never clears.
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
COMMON="-ps 100 -ng 600 -sm weak -b 20 --strong_selection_after_n_gen 480 \
        --norepeat --max-tokens-per-batch 4096 --threads 12"

go () {   # $1 init   $2 name   $3 geometry
  echo "=== $(date '+%F %T') START $2" >> "$LOG"
  python "$R/pfes.py" --start file --start-file "$R/$1" $COMMON $3 \
      -o "$OUT/$2" > "$OUT/$2.console.log" 2>&1
  echo "=== $(date '+%F %T') END   $2 (exit $?)" >> "$LOG"
}

go init/init_orfs.faa         sorf-r2      "-pl0 30  -hl0 30 -bl0 12"
go init_varlen/init_orfs.faa  long-sorf-r1 "-pl0 100 -hl0 40 -bl0 12"
go init_varlen/init_orfs.faa  long-sorf-r2 "-pl0 100 -hl0 40 -bl0 12"
