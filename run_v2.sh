#!/usr/bin/env bash
# Production series v2: one objective, three origins, two replicates each.
#
#   objective   pLDDT . pTM . P_AMP . P_len . P_helix . P_beta . contact_density
#
# The experimental variable is the ORIGIN of the starting population: sequence
# with no evolutionary history, sequence that already encodes something else,
# and (when available) sequence that is protein-coding but evolutionarily young.
# All three sets are cut to the same length, so length cannot masquerade as
# origin, and all are screened with MACREL below 0.5 so activity has to emerge.
#
# Horizon is 600 generations, double the earlier series: that series reached 95%
# of its final best only at generation 216 of 300, so the shorter horizon was
# truncating lineages that were still improving.  The original PFES protocol
# runs 4,000 stochastic + 1,000 strong; the 80/20 split is kept, the total is
# not, because this runs on CPU.
#
# Output goes to results/v2 so the earlier series -- run with -hl0 20 and the
# steeper length penalty -- survives intact as a sensitivity comparison.
set -u
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate pfes_amps
source "$(dirname "$0")/preflight.sh"   # abort if macrel/hemopi2/esm3 are not really there

R=/data/apostolos/pfes
OUT=$R/results/v2
mkdir -p "$OUT"

NGEN=600
STRONG=$(( NGEN * 80 / 100 ))       # strong selection for the final 20%

# Runs execute one at a time, each with the whole machine.  Total wall time is
# about the same as running two at half width, but the first result arrives in
# half the time, which is when a misconfiguration would show up.
CORES=$(nproc 2>/dev/null || echo 24)
THREADS=$CORES;  [ "$THREADS" -lt 4 ] && THREADS=4

# Thresholds are passed explicitly rather than left to defaults, so the command
# line is a complete record of the parameterisation:
#   -pl0 30   midpoint of the published AMP length range, 12-50 aa
#             (Mookherjee 2020; APD3 Wang 2016; DBAASP v3)
#   -hl0 30   the value specified in Sahakyan et al., PNAS 2025 -- the reference
#             implementation's default of 20 contradicts its own paper
#   -bl0 12   likewise from the paper
# Steepness constants live in pfes.py: 0.12 (length, set so the 90->10%
# transition spans the 12-50 aa class range), 0.5 (helix), 0.5 (beta, restored
# from the reference implementation's undocumented 0.6).
COMMON="-ps 100 -ng $NGEN -sm weak -b 20 --strong_selection_after_n_gen $STRONG
        --norepeat --max-tokens-per-batch 4096 --threads $THREADS
        -pl0 30 -hl0 30 -bl0 12"

run () {          # $1 = init fasta   $2 = run name
  local t0=$SECONDS
  echo "=== $(date '+%F %T')  START  $2"
  cd "$R" || { echo "!!! cannot cd $R"; return 1; }
  nice -n 10 python pfes.py --start file --start-file "$R/init/$1" \
       $COMMON -o "$OUT/$2" > "$OUT/$2.console.log" 2>&1
  echo "=== $(date '+%F %T')  DONE   $2  exit=$?  $(( (SECONDS-t0)/60 )) min"
}

# Origins alternate rather than running both replicates of one set back to back,
# so that an interruption leaves one replicate of each origin rather than one
# origin complete and the other untouched.
echo "### v2 series start $(date '+%F %T')"
echo "### $NGEN generations, strong selection from $STRONG, $THREADS threads, serial"

run init_random.faa     random-r1
run init_fragments.faa  fragments-r1
run init_random.faa     random-r2
run init_fragments.faa  fragments-r2

echo "### v2 series done  $(date '+%F %T')"
