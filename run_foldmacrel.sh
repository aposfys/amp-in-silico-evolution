#!/usr/bin/env bash
set -u
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate pfes_amps
source "$(dirname "$0")/preflight.sh"   # abort if macrel/hemopi2/esm3 are not really there
R=/data/apostolos/pfes; F=/data/apostolos/pfes-foldmacrel; OUT=$R/results/final
mkdir -p "$OUT"

run () {
  echo "=== $(date '+%F %T')  START  $2"
  cd "$F" || return 1
  nice -n 10 python pfes.py --start file --start-file "$R/init/$1" \
    -ps 100 -ng 300 -sm weak -b 20 --strong_selection_after_n_gen 240 \
    --norepeat --max-tokens-per-batch 4096 --threads 12 \
    -o "$OUT/$2" > "$OUT/$2.console.log" 2>&1
  echo "=== $(date '+%F %T')  DONE   $2  (exit $?)"
}

echo "###  MACREL + ESM3 arm.  Έναρξη $(date '+%F %T')"
run init_random.faa     random-foldmacrel-r1
run init_random.faa     random-foldmacrel-r2
run init_fragments.faa  fragments-foldmacrel-r1
run init_fragments.faa  fragments-foldmacrel-r2
echo "###  ΤΕΛΟΣ  $(date '+%F %T')"
