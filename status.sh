#!/usr/bin/env bash
# Κατάσταση μιας σειράς PFES.  Χρήση: ./status.sh [v3]
S=${1:-v3}
cd /data/apostolos/pfes || exit 1
NG=$(grep -oE 'NGEN=[0-9]+' run_${S}.sh 2>/dev/null | head -1 | tr -dc '0-9'); NG=${NG:-600}

echo "=== Διεργασία ==="
ps -eo pid,ni,pcpu,etime,args | grep -E "[p]fes\.py.*results/$S" || echo "  (καμία)"

echo; echo "=== Σειρά ==="
tail -6 "results/${S}.master.log" 2>/dev/null

echo; echo "=== Πρόοδος ==="
for d in results/$S/*/; do
  f="$d/progress.log"; [ -f "$f" ] || continue
  n=$(basename "$d")
  g=$(tail -1 "$f" | awk '{x=$1; gsub(/[^0-9]/,"",x); print x}')
  st=$(grep "START $n\$" "results/${S}.master.log" | tail -1 | awk '{print $2" "$3}')
  end=$(stat -c %Y "$f")
  if [ -n "$st" ] && [ -n "$g" ] && [ "$g" -gt 0 ]; then
    s=$(date -d "$st" +%s)
    r=$(awk -v a=$s -v b=$end -v g=$g 'BEGIN{printf "%.2f",(b-a)/60/g}')
    l=$(awk -v r=$r -v g=$g -v N=$NG 'BEGIN{printf "%.1f", r*(N-g)/60}')
    printf "  %-20s gen %4s/%s  %5s λ/γενιά  ~%s ώρες ακόμη\n" "$n" "$g" "$NG" "$r" "$l"
  fi
done

echo; echo "=== Τάση (τρέχον run) ==="
cur=$(ls -td results/$S/*/ 2>/dev/null | head -1)
[ -n "$cur" ] && awk '$15 ~ /^[0-9.]+$/ {g=$1; gsub(/[^0-9]/,"",g)
  s[g]+=$15; a[g]+=$13; p[g]+=$8; L[g]+=$3; c[g]++}
  END {mx=0; for(k in c) if(k+0>mx) mx=k+0
    for(i=0;i<=mx;i+=50) if(i in c)
      printf "   gen %4d  score %.4f  amp %.4f  plddt %.3f  len %.1f\n",i,s[i]/c[i],a[i]/c[i],p[i]/c[i],L[i]/c[i]
    printf "   gen %4d  score %.4f  amp %.4f  plddt %.3f  len %.1f  <- τώρα\n",mx,s[mx]/c[mx],a[mx]/c[mx],p[mx]/c[mx],L[mx]/c[mx]}' \
  "$cur/progress.log"

echo; echo "=== Μηχάνημα ($(nproc) πυρήνες) ==="; uptime
