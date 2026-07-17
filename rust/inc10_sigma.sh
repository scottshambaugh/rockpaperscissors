#!/usr/bin/env bash
# Phase 2 of the inclusive(10) production run: the 7 heavy sigma-correction
# types. Small brutes and coprime types are computed by sigma_gate.py itself;
# this script only produces the --fix values for the heavy ones.
#
# Each heavy type is sharded `sfix TYPE s OF` (OF a power of 3); the sum over
# shards is the Fix value (shard-sum identity validated at n=8). (2,2,1^6) is
# instead the sigma_sweep marked-pair pass over the grandparent stream.
#
#   RUST=/mnt/.../rust WORK=~/inc10 WORKERS=3 bash inc10_sigma.sh
set -euo pipefail
RUST="${RUST:-$(cd "$(dirname "$0")" && pwd)}"
WORK="${WORK:-$HOME/inc10}"
WORKERS="${WORKERS:-3}"
SHARDS="${SHARDS:-32}"
LINK="/tmp/inc10_bshim.o -L/usr/local/lib -lnauty"
NAUTY_INC="$(dirname "$(find /usr/include /usr/local/include -name nauty.h 2>/dev/null | head -1)")"
mkdir -p "$WORK/sigma"
cd "$WORK/sigma"
gcc -O2 -c "$RUST/balanced_shim.c" -I"$NAUTY_INC" -o /tmp/inc10_bshim.o
rustc -O -C target-cpu=native "$RUST/sigma_fix.rs"    -o /tmp/p_sfix
rustc -O -C target-cpu=native "$RUST/sigma_ktuple.rs" -o /tmp/p_ktup -C link-args="$LINK"
rustc -O -C target-cpu=native "$RUST/sigma_sweep.rs"  -o /tmp/p_ssweep -C link-args="$LINK"

# (2^k,1^f) types via the k-tuple class sweep (anchored: (2,2,2)=66,
# (2,2,2,1,1)=32310, (2,2,2,2)=20298); minutes instead of core-hours
if [ ! -f fix_2_2_2_1_1_1_1.txt ]; then
  nauty-geng 7 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/p_ktup 7 3 2>/dev/null \
    | grep -o "Fix_(2^3,1^4) = [0-9]*" | grep -o "[0-9]*$" > fix_2_2_2_1_1_1_1.txt
fi
echo "2,2,2,1,1,1,1 = $(cat fix_2_2_2_1_1_1_1.txt) [ktup]"
if [ ! -f fix_2_2_2_2_1_1.txt ]; then
  nauty-geng 6 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/p_ktup 6 4 2>/dev/null \
    | grep -o "Fix_(2^4,1^2) = [0-9]*" | grep -o "[0-9]*$" > fix_2_2_2_2_1_1.txt
fi
echo "2,2,2,2,1,1 = $(cat fix_2_2_2_2_1_1.txt) [ktup]"

# remaining heavy brute types: TYPE:OF (OF = 3^k shard count)
for spec in "4,2,1,1,1,1:729" "3,2,2,1,1,1:729" "3,3,1,1,1,1:729" \
            "2,2,2,2,2:19683"; do
  TYPE="${spec%:*}"; OF="${spec#*:}"
  safe="${TYPE//,/_}"
  [ -f fix_$safe.txt ] && { echo "$TYPE = $(cat fix_$safe.txt) [cached]"; continue; }
  echo "=== $TYPE  ($OF shards) ==="
  seq 0 $((OF-1)) | xargs -P"$WORKERS" -I@ sh -c "/tmp/p_sfix $TYPE @ $OF 2>/dev/null | grep -o 'fix=[0-9]*'" \
    | awk -F= '{s+=$2} END{print s}' > fix_$safe.txt
  echo "$TYPE = $(cat fix_$safe.txt)"
done

# (2,2,1^6): normally produced by the tee'd strata run (inc10_fix_2_2_1x6.txt);
# recompute here only if missing
if [ -f "$WORK/inc10_fix_2_2_1x6.txt" ] && [ ! -f fix_2_2_1_1_1_1_1_1.txt ]; then
  cp "$WORK/inc10_fix_2_2_1x6.txt" fix_2_2_1_1_1_1_1_1.txt
fi
if [ ! -f fix_2_2_1_1_1_1_1_1.txt ]; then
  echo "=== 2,2,1^6 via sigma_sweep ==="
  seq 0 $((SHARDS-1)) | xargs -P"$WORKERS" -I@ sh -c '
    s=@
    nauty-geng 8 $s/'"$SHARDS"' 2>/dev/null | nauty-directg -o 2>/dev/null \
      | /tmp/p_ssweep 8 2>/dev/null | grep -E "^n=8 " > sweep_$s.txt
    echo "sweep shard $s done ($(date +%H:%M:%S))"
  '
  # sum rawJ1,rawJ0,rawK0b,rawH1 and L_inc across shards, then apply the formula
  python3 - <<'PY' > fix_2_2_1_1_1_1_1_1.txt
import glob, re
acc = {"L_inc":0,"rawJ1":0,"rawJ0":0,"rawK0b":0,"rawH1":0}
for f in glob.glob("sweep_*.txt"):
    line = open(f).read()
    for k in acc:
        m = re.search(rf"{k}=(\d+)", line)
        if m: acc[k] += int(m.group(1))
# The S_10 type (2,2,1^6) has m = n-2 = 8 supervertices; the sweep runs over
# 8-vertex classes and pairs = m(m-1) = 56. (Sweep-at-8 == S_10 type, exactly
# as sweep-at-6 == the S_8 type (2,2,1^4) in the n=6 anchor.) Sum raw first,
# then divide once (per-shard raws need not be divisible by 56 individually).
pairs = 8*7
fix = (2*acc["rawJ1"] + acc["rawJ0"] + 4*acc["rawH1"] + 2*acc["rawK0b"]) // pairs
assert (2*acc["rawJ1"] + acc["rawJ0"] + 4*acc["rawH1"] + 2*acc["rawK0b"]) % pairs == 0
assert acc["L_inc"] == 46778967018, f"sweep L_inc(8) cross-check failed: {acc['L_inc']}"
print(fix)
PY
  echo "2,2,1^6 = $(cat fix_2_2_1_1_1_1_1_1.txt)  (also cross-check sweep L_inc(8)=46778967018)"
fi

echo ""
echo "Assemble with:"
echo "  python3 $RUST/sigma_gate.py 10 \\"
echo "    --linc 2=0,3=2,4=42,5=978,6=130950,7=49473198,8=46778967018,9=235837146265362,10=\$L10 \\"
for f in fix_*.txt; do
  t="${f#fix_}"; t="${t%.txt}"; t="${t//_/,}"
  echo "    --fix $t=$(cat "$f") \\"
done
