#!/usr/bin/env bash
# Waits for sorf-r2 to close, then replaces the serial sorf chain with a
# dedicated stream for long-sorf-r1.  Counts a NEW END line rather than any
# END line, so a leftover from the killed first attempt cannot trigger it.
set -u
R=/data/apostolos/pfes
LOG=$R/results/v3.master.log
KEY='END   sorf-r2'
BASE=$(grep -c "$KEY" "$LOG" || true)
until [ "$(grep -c "$KEY" "$LOG" || true)" -gt "$BASE" ]; do sleep 60; done
tmux kill-session -t sorf 2>/dev/null
sleep 10
rm -rf "$R/results/v3/long-sorf-r1"
echo "### $(date '+%F %T') sorf chain retired, long-sorf-r1 moved to own stream" >> "$LOG"
exec bash "$R/run_v3_one.sh" init_varlen/init_orfs.faa long-sorf-r1 -pl0 100 -hl0 40 -bl0 12
