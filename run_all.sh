#!/usr/bin/env bash
set -u
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate pfes_amps
R=/data/apostolos/pfes; N=/data/apostolos/pfes-noamp; OUT=$R/results/final
mkdir -p "$OUT"

run () {   # $1 φάκελος  $2 tag
  echo "=== $(date '+%F %T')  START  $2"
  cd "$1" || return 1
  nice -n 10 python pfes.py --start file --start-file "$R/init/init_random.faa" \
    -ps 100 -ng 300 -sm weak -b 20 --strong_selection_after_n_gen 240 \
    --norepeat --max-tokens-per-batch 4096 --threads 12 \
    -o "$OUT/$2" > "$OUT/$2.console.log" 2>&1
  echo "=== $(date '+%F %T')  DONE   $2  (exit $?)"
}

echo "###  4 runs, ~56 ώρες.  Έναρξη $(date '+%F %T')"
run "$R" random-pfesmacrel-r1
run "$R" random-pfesmacrel-r2
run "$N" random-pfes-r1
run "$N" random-pfes-r2
echo "###  ΤΕΛΟΣ  $(date '+%F %T')"
