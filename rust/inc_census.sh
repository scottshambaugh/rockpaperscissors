#!/usr/bin/env bash
# Full inclusive(N) census runner, N-parametric (N even, 8..16): the SAME code
# path serves the n=8 rehearsal gate and the n>=10 production runs. L_inc(N) =
# sum of all even nullity strata L2..L_{N-2} (the high-strata loop below is what
# makes it complete at N>=10; see README "Inclusive census").
#
#   bash inc_census.sh N            (env: WORK, WORKERS, SHARDS, RUST)
#
# Phases (each resumable; all gates fatal):
#   1. strata: one pass over the (N-2)-vertex grandparent stream per shard,
#      tee'd into inc10 (stratum 2), f3x d=3/d=5 unique parents, and
#      sigma_sweep (the (2,2,1^{N-4}) statistics + L_inc(N-2) GATE).
#      Then inc4 over the concatenated parent streams (strata 4 and 6).
#   2. sigma: k-tuple sweeps for (2,2,2,1^{N-6}) and, at N=10, (2,2,2,2,1^2);
#      sharded sigma_fix brutes for the remaining heavy types; small types and
#      coprime types are computed inside sigma_gate.py.
#   3. assembly: sigma_gate.py N with the L_inc table; integrality gate; at
#      N=8 the result must equal 1,198,013 exactly.
#
# Known values (gates + coprime-rule inputs). L_inc(9) from the labsum run.
set -euo pipefail
N="${1:?usage: inc_census.sh N}"
RUST="${RUST:-$(cd "$(dirname "$0")" && pwd)}"
WORK="${WORK:-$HOME/inc_census_$N}"
WORKERS="${WORKERS:-3}"
SHARDS="${SHARDS:-64}"
LINK="/tmp/incc_bshim.o -L/usr/local/lib -lnauty"
NAUTY_INC="$(dirname "$(find /usr/include /usr/local/include -name nauty.h 2>/dev/null | head -1)")"
Q=$((N-2))
# LINC_Q = L_inc(N-2), the sweep gate. It is the FULL inclusive-labeled count of
# the (N-2)-vertex games, i.e. it must already include every nullity stratum for
# N-2 (see the high-strata loop below). L_inc(10) here includes L8 = 8,884,050.
case "$N" in
  8)  LINC_Q=130950                 # L_inc(6)
      EXPECT=1198013 ;;
  10) LINC_Q=46778967018            # L_inc(8): complete (no nullity-8 at n=8)
      EXPECT="" ;;
  12) LINC_Q=1297146370885586550    # L_inc(10) = L2+L4+L6+L8 (incl. nullity-8)
      EXPECT="" ;;
  *)  echo "N must be even in 8..=16 (dimension caps are width-16)"; exit 1 ;;
esac
LINC_TABLE="2=0,3=2,4=42,5=978,6=130950,7=49473198,8=46778967018,9=235837146265362,10=1297146370885586550"
mkdir -p "$WORK/strata" "$WORK/sigma"

gcc -O2 -c "$RUST/balanced_shim.c" -I"$NAUTY_INC" -o /tmp/incc_bshim.o
rustc -O -C target-cpu=native -C overflow-checks=on "$RUST/inc10.rs"        -o /tmp/c_inc10  -C link-args="$LINK"
rustc -O -C target-cpu=native -C overflow-checks=on "$RUST/inc_stratum.rs"         -o /tmp/c_incstrat   -C link-args="$LINK"
rustc -O -C target-cpu=native -C overflow-checks=on "$RUST/f3x.rs"          -o /tmp/c_f3x    -C link-args="$LINK"
rustc -O -C target-cpu=native -C overflow-checks=on "$RUST/sigma_sweep.rs"  -o /tmp/c_ssweep -C link-args="$LINK"
rustc -O -C target-cpu=native -C overflow-checks=on "$RUST/sigma_ktuple.rs" -o /tmp/c_ktup   -C link-args="$LINK"
rustc -O -C target-cpu=native -C overflow-checks=on "$RUST/sigma_fix.rs"    -o /tmp/c_sfix

# ---------------- phase 1: strata ----------------
cd "$WORK/strata"
export N Q SHARDS
seq 0 $((SHARDS-1)) | xargs -P"$WORKERS" -I@ bash -c '
  s=@
  [ -f shard_$s.done ] && exit 0
  fd="fifo_$s"; rm -rf "$fd"; mkdir "$fd"
  mkfifo "$fd/a" "$fd/b" "$fd/d"
  /tmp/c_inc10 "$N" < "$fd/a" 2>/dev/null | grep -o "L_nullity2_labeled=[0-9]*" > s2_$s.txt & p1=$!
  # fused f3x: d=3 stream to stdout, d=5 stream to the path arg -- one kernel
  # computation per orientation instead of a second full stream pass
  /tmp/c_f3x "$Q" 3 u p5_$s.d6 < "$fd/b" 2>/dev/null > p3_$s.d6 & p2=$!
  /tmp/c_ssweep "$Q" < "$fd/d" 2>/dev/null | grep -E "^n=$Q " > sw_$s.txt & p4=$!
  nauty-geng "$Q" $s/'"$SHARDS"' 2>/dev/null | nauty-directg -o 2>/dev/null \
    | tee "$fd/a" "$fd/b" > "$fd/d"
  wait $p1 $p2 $p4
  rm -rf "$fd"
  if [ ! -s s2_$s.txt ] || [ ! -s sw_$s.txt ]; then
    echo "shard $s: consumer output missing -- NOT marking done" >&2
    exit 1
  fi
  # strata 4/6 per shard: f3x-unique partitions parent classes across shards,
  # so partial sums add exactly (like L2) and the parent files can be deleted
  # immediately -- the full-stream parent set would not fit on disk
  wc -l < p3_$s.d6 > pc3_$s.txt
  wc -l < p5_$s.d6 > pc5_$s.txt
  /tmp/c_incstrat "$N"   < p3_$s.d6 | grep -o "L_nullity4_labeled=[0-9]*" > s4_$s.txt
  /tmp/c_incstrat "$N" 6 < p5_$s.d6 | grep -o "L_nullity6_labeled=[0-9]*" > s6_$s.txt
  if [ ! -s s4_$s.txt ] || [ ! -s s6_$s.txt ]; then
    echo "shard $s: stratum-4/6 output missing -- NOT marking done" >&2
    exit 1
  fi
  rm -f p3_$s.d6 p5_$s.d6
  touch shard_$s.done
  echo "shard $s done, parents $(cat pc3_$s.txt)/$(cat pc5_$s.txt) ($(date +%H:%M:%S))"
'
# EXACT integer sum (awk uses f64 doubles -> silently rounds values > 2^53;
# L2 ~1.3e18 was corrupted by ~112). Python bignum has no precision ceiling.
isum() { python3 -c "import glob,sys; print(sum(int(open(f).read().split('=')[-1]) for f in glob.glob(sys.argv[1]) if open(f).read().strip()))" "$1"; }
isum 's2_*.txt' > "$WORK/L2.txt"
echo "L2 = $(cat "$WORK/L2.txt")"

# sigma sweep totals: L_inc(Q) GATE + the (2,2,1^{N-4}) formula
LINC_Q="$LINC_Q" Q="$Q" python3 - <<'PY' > "$WORK/fix_pairsweep.txt"
import glob, os, re
acc = {"L_inc": 0, "rawJ1": 0, "rawJ0": 0, "rawK0b": 0, "rawH1": 0}
for f in glob.glob("sw_*.txt"):
    line = open(f).read()
    for k in acc:
        m = re.search(rf"{k}=(\d+)", line)
        if m:
            acc[k] += int(m.group(1))
want = int(os.environ["LINC_Q"])
assert acc["L_inc"] == want, f"sweep L_inc GATE FAILED: {acc['L_inc']} != {want}"
q = int(os.environ["Q"])
raw = 2 * acc["rawJ1"] + acc["rawJ0"] + 4 * acc["rawH1"] + 2 * acc["rawK0b"]
assert raw % (q * (q - 1)) == 0, "not integral"
print(raw // (q * (q - 1)))
PY
echo "Fix_(2,2,1^$((N-4))) = $(cat "$WORK/fix_pairsweep.txt")  [sweep L_inc($Q) gate passed]"

pc3sum=$(python3 -c "import glob; print(sum(int(open(f).read()) for f in glob.glob('pc3_*.txt') if open(f).read().strip()))")
pc5sum=$(python3 -c "import glob; print(sum(int(open(f).read()) for f in glob.glob('pc5_*.txt') if open(f).read().strip()))")
echo "parents: d=3 $pc3sum, d=5 $pc5sum"
isum 's4_*.txt' > "$WORK/L4.txt"
isum 's6_*.txt' > "$WORK/L6.txt"
echo "L4 = $(cat "$WORK/L4.txt")  L6 = $(cat "$WORK/L6.txt")"

# ---- high-nullity strata L8, L10, ... (D = 8, 10, ..., N-2) ----
# A rank-(N-D) inclusive game (e.g. the rank-2 rock-paper-scissors blow-ups) has
# nullity D = N - rank. At N=8 the top such family is nullity-6 (in L6); at N>=10
# it reaches nullity 8, 10, ... which L2+L4+L6 alone MISS. These families are
# rare (low rank), so one extra non-fused f3x pass per stratum is cheap. Each
# child nullity D comes from a parent of nullity D-1 (odd), i.e. f3x d=D-1, and
# is counted by inc_stratum N D. f3x-unique partitions parent classes across shards,
# so per-shard sums add exactly (same as L4/L6). The zero matrix (nullity N) is
# never inclusive, so the top stratum is D = N-2.
HI_SUM=""
for ((D=8; D<=N-2; D+=2)); do
  dd=$((D-1))
  export D dd
  seq 0 $((SHARDS-1)) | xargs -P"$WORKERS" -I@ bash -c '
    s=@
    [ -f hi_${D}_$s.done ] && exit 0
    nauty-geng "$Q" $s/'"$SHARDS"' 2>/dev/null | nauty-directg -o 2>/dev/null \
      | /tmp/c_f3x "$Q" "$dd" u 2>/dev/null > pD_${D}_$s.d6
    /tmp/c_incstrat "$N" "$D" < pD_${D}_$s.d6 2>/dev/null \
      | grep -o "L_nullity${D}_labeled=[0-9]*" > sD_${D}_$s.txt
    if [ ! -s sD_${D}_$s.txt ]; then echo "shard $s stratum $D: no output" >&2; exit 1; fi
    rm -f pD_${D}_$s.d6
    touch hi_${D}_$s.done
  '
  isum "sD_${D}_*.txt" > "$WORK/L$D.txt"
  echo "L$D = $(cat "$WORK/L$D.txt")"
  HI_SUM="$HI_SUM + $(cat "$WORK/L$D.txt")"
done

L10=$(python3 -c "print($(cat "$WORK/L2.txt")+$(cat "$WORK/L4.txt")+$(cat "$WORK/L6.txt")$HI_SUM)")
echo "$L10" > "$WORK/Linc_N.txt"
echo "L_inc($N) = L2+L4+L6${HI_SUM:+ + higher strata} = $L10"

# ---------------- phase 2: heavy sigma types ----------------
cd "$WORK/sigma"
# (2,2,2,1^{N-6}) via ktup k=3 over the (N-3)-vertex stream
if [ ! -f fix_k3.txt ]; then
  nauty-geng $((N-3)) 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/c_ktup $((N-3)) 3 2>/dev/null \
    | grep -o "= [0-9]*$" | grep -o "[0-9]*" > fix_k3.txt
fi
echo "(2,2,2,1^$((N-6))) = $(cat fix_k3.txt) [ktup]"
FIXES=(--fix "2,2,2$(printf ',1%.0s' $(seq 1 $((N-6))))=$(cat fix_k3.txt)")
FIXES+=(--fix "2,2$(printf ',1%.0s' $(seq 1 $((N-4))))=$(cat "$WORK/fix_pairsweep.txt")")
if [ "$N" -ge 10 ]; then
  if [ ! -f fix_k4.txt ]; then
    nauty-geng $((N-4)) 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/c_ktup $((N-4)) 4 2>/dev/null \
      | grep -o "= [0-9]*$" | grep -o "[0-9]*" > fix_k4.txt
  fi
  echo "(2,2,2,2,1^$((N-8))) = $(cat fix_k4.txt) [ktup]"
  FIXES+=(--fix "2,2,2,2$(printf ',1%.0s' $(seq 1 $((N-8))))=$(cat fix_k4.txt)")
  # remaining heavy brutes at N=10, sharded sigma_fix
  for spec in "4,2,1,1,1,1:729" "3,2,2,1,1,1:729" "3,3,1,1,1,1:729" "2,2,2,2,2:729"; do
    TYPE="${spec%:*}"; OF="${spec#*:}"; safe="${TYPE//,/_}"
    if [ ! -f fix_$safe.txt ]; then
      echo "=== brute $TYPE ($OF shards) ==="
      seq 0 $((OF-1)) | xargs -P"$WORKERS" -I@ sh -c "/tmp/c_sfix $TYPE @ $OF 2>/dev/null | grep -o 'fix=[0-9]*'" \
        | python3 -c "import sys; print(sum(int(l.split('=')[1]) for l in sys.stdin if l.strip()))" > fix_$safe.txt
    fi
    echo "$TYPE = $(cat fix_$safe.txt) [brute]"
    FIXES+=(--fix "$TYPE=$(cat fix_$safe.txt)")
  done
fi

# ---------------- phase 3: assembly ----------------
python3 "$RUST/sigma_gate.py" "$N" \
  --linc "$LINC_TABLE,$N=$(cat "$WORK/Linc_N.txt")" \
  "${FIXES[@]}" | tee "$WORK/assembly.txt"
RESULT=$(grep -o "inclusive($N) = [0-9]*" "$WORK/assembly.txt" | grep -o "[0-9]*$" || true)
if [ -n "$EXPECT" ]; then
  if [ "$RESULT" = "$EXPECT" ]; then
    echo "=== REHEARSAL GATE PASSED: inclusive($N) = $RESULT ==="
  else
    echo "=== REHEARSAL GATE FAILED: got '$RESULT' want $EXPECT ==="; exit 1
  fi
else
  echo "=== inclusive($N) = $RESULT ==="
fi
