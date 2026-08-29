#!/usr/bin/env bash
# amp-in-silico-evolution, series v3. Eight runs, serial, on a free node.
#
# Completes the design the v2 series opened. v2 ran random and fragments at the
# 25-residue cut; v3 adds the sORF arm there, and opens a second length regime
# in which all three origins are run again with the length term acting as a cap
# at the definitional boundary of the class rather than as a target.
#
#   1-2  sORF       25 aa fixed     directly comparable to the four v2 runs
#   3-4  random     10-100 aa       new regime
#   5-6  fragments  10-100 aa
#   7-8  sORF       10-100 aa
#
# Differences from run_v2.sh, both deliberate:
#
#   No `nice -n 10`. Under the contention that dominated the v2 series (load
#   ~120 against 48 cores) that niceness cost the first run a factor of
#   seventeen, 80 hours against 4h40m for the three that ran on a quiet node.
#   Equal footing, not a claim on the machine.
#
#   -hl0 40 in the long arm. The helix penalty uses c=0.5, so hl0=30 assigns a
#   factor of 0.076 to a 35-residue helix. Raising the chain limit without
#   raising this one does not remove a constraint, it moves it: a long peptide
#   would be forced into short disconnected elements, which is the opposite of
#   the single extended amphipathic helix the class is built on. LL-37 is 37
#   residues with a helix spanning about 30, so 40 admits what exists in nature
#   and closes above it.
#
# Check the node is quiet before starting:  uptime

set -u

# Thread oversubscription. torch intra-op, torch inter-op, MKL and OpenMP nest,
# so --threads 48 produced 191 threads on 48 cores. Under any contention at all
# the OpenMP spin barrier then convoys and throughput collapses non-linearly:
# 26 s per generation on a free machine against 415 s at a run queue of 52.
# Sleep rather than spin, pin every layer to one count, leave 8 cores of slack.
export KMP_BLOCKTIME=0
export OMP_WAIT_POLICY=PASSIVE
export OMP_NUM_THREADS=40
export MKL_NUM_THREADS=40
export TOKENIZERS_PARALLELISM=false
source "$(dirname "$0")/preflight.sh"   # abort if macrel/hemopi2/esm3 are not really there

R=/data/apostolos/pfes
OUT=$R/results/v3
LOG=$R/results/v3.master.log
NGEN=600
STRONG=$(( NGEN * 80 / 100 ))
THREADS=40

SHORT="-pl0 30  -hl0 30 -bl0 12"
LONG="-pl0 100 -hl0 40 -bl0 12"
COMMON="-ps 100 -ng $NGEN -sm weak -b 20 \
        --strong_selection_after_n_gen $STRONG --norepeat \
        --max-tokens-per-batch 4096 --threads $THREADS"

# --------------------------------------------------------------------------- #
# Refuse to start on a missing input rather than fail four hours in.
missing=0
for f in init/init_orfs.faa \
         init_varlen/init_random.faa \
         init_varlen/init_fragments.faa \
         init_varlen/init_orfs.faa; do
    if [ ! -s "$R/$f" ]; then
        echo "MISSING: $R/$f" >&2
        missing=1
    else
        printf '  ok  %-34s %s sequences\n' "$f" "$(grep -c '^>' "$R/$f")"
    fi
done
[ "$missing" -eq 0 ] || { echo "aborting: build the init sets first" >&2; exit 1; }

mkdir -p "$OUT"
{
  echo "### v3 series start $(date '+%F %T')"
  echo "### $NGEN generations, strong selection from $STRONG, $THREADS threads, serial"
  echo "### short arm: $SHORT"
  echo "### long arm:  $LONG"
  echo "### load at start: $(uptime | sed 's/.*load average/load average/')"
} >> "$LOG"

run () {           # $1 = init fasta (path under $R)   $2 = run name   $3 = arm
    local init="$R/$1" name="$2" arm="$3"
    echo "=== $(date '+%F %T') START $name" >> "$LOG"
    python "$R/pfes.py" --start file --start-file "$init" \
        $COMMON $arm \
        -o "$OUT/$name" \
        > "$OUT/$name.console.log" 2>&1
    echo "=== $(date '+%F %T') END   $name (exit $?)" >> "$LOG"
}

# Short arm first: it is the one that completes an existing comparison, so if
# the node is taken back it is the half that still yields a result.
run init/init_orfs.faa            sorf-r1        "$SHORT"
run init/init_orfs.faa            sorf-r2        "$SHORT"

run init_varlen/init_random.faa    long-random-r1     "$LONG"
run init_varlen/init_fragments.faa long-fragments-r1  "$LONG"
run init_varlen/init_orfs.faa      long-sorf-r1       "$LONG"
run init_varlen/init_random.faa    long-random-r2     "$LONG"
run init_varlen/init_fragments.faa long-fragments-r2  "$LONG"
run init_varlen/init_orfs.faa      long-sorf-r2       "$LONG"

echo "### v3 series done $(date '+%F %T')" >> "$LOG"
