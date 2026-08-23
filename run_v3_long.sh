  #!/usr/bin/env bash
  # One origin of the long arm, both replicates, serially within the stream.
  # Three of these run concurrently: the per-generation scoring loop is
  # single-threaded for roughly 60% of each generation, so concurrent streams
  # fill each other's gaps rather than competing.
  # Usage: run_v3_long.sh {random|fragments|orfs}
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
  NGEN=600
  STRONG=$(( NGEN * 80 / 100 ))
  LONG="-pl0 100 -hl0 40 -bl0 12"
NGEN=600
STRONG=$(( NGEN * 80 / 100 ))
LONG="-pl0 100 -hl0 40 -bl0 12"
COMMON="-ps 100 -ng $NGEN -sm weak -b 20 \
        --strong_selection_after_n_gen $STRONG --norepeat \
        --max-tokens-per-batch 4096 --threads 12"

case "${1:-}" in
  random)    INIT=init_varlen/init_random.faa;    TAG=long-random ;;
  fragments) INIT=init_varlen/init_fragments.faa; TAG=long-fragments ;;
  orfs)      INIT=init_varlen/init_orfs.faa;      TAG=long-sorf ;;
  *) echo "usage: $0 {random|fragments|orfs}" >&2; exit 1 ;;
esac
[ -s "$R/$INIT" ] || { echo "MISSING $R/$INIT" >&2; exit 1; }

for rep in r1 r2; do
  name="$TAG-$rep"
  echo "=== $(date '+%F %T') START $name" >> "$LOG"
  python "$R/pfes.py" --start file --start-file "$R/$INIT" \
      $COMMON $LONG -o "$OUT/$name" > "$OUT/$name.console.log" 2>&1
  echo "=== $(date '+%F %T') END   $name (exit $?)" >> "$LOG"
done
