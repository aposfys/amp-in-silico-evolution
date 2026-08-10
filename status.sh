#!/usr/bin/env bash
# Κατάσταση της σειράς PFES v2. Read-only.
cd /data/apostolos/pfes || exit 1

echo "=== Διεργασίες ==="
ps -eo pid,user,etime,pcpu,args | grep -E "[p]fes\.py|[r]un_v2\.sh" || echo "  (καμία — η σειρά δεν τρέχει)"

echo
echo "=== Σειρά (master log) ==="
tail -6 v2.master.log

echo
echo "=== Πρόοδος ανά run ==="
for d in results/v2/*/; do
  f="$d/progress.log"; [ -f "$f" ] || continue
  name=$(basename "$d")
  gen=$(tail -1 "$f" | awk '{g=$1; gsub(/[^0-9]/,"",g); print g}')
  start=$(grep "START $name\$" v2.master.log | tail -1 | awk '{print $2" "$3}')
  end=$(stat -c %Y "$f")
  if [ -n "$start" ] && [ -n "$gen" ] && [ "$gen" -gt 0 ]; then
    s=$(date -d "$start" +%s)
    rate=$(awk -v a=$s -v b=$end -v g=$gen 'BEGIN{printf "%.1f", (b-a)/60/g}')
    left=$(awk -v r=$rate -v g=$gen 'BEGIN{printf "%.1f", r*(600-g)/60/24}')
    printf "  %-14s gen %4s/600  %5s min/gen  ~%s μέρες ακόμη  (τελ. %s)\n" \
      "$name" "$gen" "$rate" "$left" "$(date -d @$end '+%m-%d %H:%M')"
  else
    printf "  %-14s gen %4s  (χωρίς χρόνο έναρξης)\n" "$name" "${gen:-?}"
  fi
done

echo
echo "=== Τάση fitness (τρέχον run) ==="
cur=$(ls -td results/v2/*/ | head -1)
echo "  $(basename $cur)   [γενιά  μέσο_score  μέσο_amp_prob  μέσο_plddt  μέσο_len]"
awk '$15 ~ /^-?[0-9.]+$/ {g=$1; gsub(/[^0-9]/,"",g);
     sc[g]+=$15; ap[g]+=$13; pl[g]+=$8; ln[g]+=$3; n[g]++}
     END {for (k in n) printf "%d %.4f %.4f %.2f %.1f\n", k, sc[k]/n[k], ap[k]/n[k], pl[k]/n[k], ln[k]/n[k]}' \
     "$cur/progress.log" | sort -n > /tmp/pfes_gen_means
(head -3 /tmp/pfes_gen_means; echo "  ..."; tail -5 /tmp/pfes_gen_means) \
  | awk '{ if (NF==1) print "     ..."; else printf "     %-6s %-10s %-10s %-8s %s\n", $1,$2,$3,$4,$5 }'

echo
echo "=== Φόρτος μηχανήματος ($(nproc) πυρήνες) ==="
uptime
