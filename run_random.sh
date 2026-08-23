#!/usr/bin/env bash
# Runs 1-4: random set, both arms, 2 replicates each. Strictly serial.
set -u

# Το tmux ξεκινά καθαρό shell και δεν κληρονομεί το conda environment.
# Χωρίς αυτό το `python` είναι του base και λείπουν esm, macrel, hemopi2.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate pfes_amps
source "$(dirname "$0")/preflight.sh"   # abort if macrel/hemopi2/esm3 are not really there
python -c "import esm, score" 2>/dev/null || { echo "!!! λάθος environment"; exit 1; }
echo "environment: $CONDA_DEFAULT_ENV"


R=/data/apostolos/pfes
F=/data/apostolos/pfes-foldonly
INIT="$R/init/init_random.faa"
OUT="$R/results/final"

mkdir -p "$OUT"

run () {
  local dir=$1 tag=$2
  echo "=== $(date '+%F %T')  START  $tag"
  cd "$dir" || { echo "!!! δεν υπάρχει $dir"; return 1; }
  nice -n 10 python pfes.py \
      --start file --start-file "$INIT" \
      -ps 100 -ng 300 -sm weak -b 20 --strong_selection_after_n_gen 240 \
      --norepeat --max-tokens-per-batch 4096 --threads 12 \
      -o "$OUT/$tag" > "$OUT/$tag.console.log" 2>&1
  local rc=$?
  echo "=== $(date '+%F %T')  DONE   $tag  (exit $rc)"
  return $rc
}

echo "###  4 runs, σειριακά, ~64 ώρες.  Έναρξη $(date '+%F %T')"
run "$R" random-structured-r1
run "$R" random-structured-r2
run "$F" random-foldonly-r1
run "$F" random-foldonly-r2
echo "###  ΟΛΑ ΤΕΛΕΙΩΣΑΝ  $(date '+%F %T')"
