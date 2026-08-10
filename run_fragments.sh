#!/usr/bin/env bash
set -u
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate pfes_amps
R=/data/apostolos/pfes; N=/data/apostolos/pfes-noamp; OUT=$R/results/final
INIT=$R/init/init_fragments.faa
mkdir -p "$OUT"

run () {
  echo "=== $(date '+%F %T')  START  $2"
  cd "$1" || return 1
  nice -n 10 python pfes.py --start file --start-file "$INIT" \
    -ps 100 -ng 300 -sm weak -b 20 --strong_selection_after_n_gen 240 \
    --norepeat --max-tokens-per-batch 4096 --threads 12 \
    -o "$OUT/$2" > "$OUT/$2.console.log" 2>&1
  echo "=== $(date '+%F %T')  DONE   $2  (exit $?)"
}

echo "###  4 runs, fragments, ~70 ώρες.  Έναρξη $(date '+%F %T')"
run "$R" fragments-pfesmacrel-r1
run "$R" fragments-pfesmacrel-r2
run "$N" fragments-pfes-r1
run "$N" fragments-pfes-r2
echo "###  ΤΕΛΟΣ  $(date '+%F %T')"
