#!/usr/bin/env bash
# Phase 1 of the inclusive(10) production run, one-pass edition: each shard of
# the 8-vertex grandparent stream is generated ONCE and tee'd into all four
# consumers simultaneously:
#   * inc10 10        -> stratum-2 labeled sum (the heavyweight)
#   * f3x 8 3 u       -> 9-vertex F3' parents, each class emitted exactly once
#                        (argmax acceptance; no labelg, no sort; shard outputs
#                        concatenate because a child class determines its
#                        grandparent class, which lives in exactly one shard)
#   * f3x 8 5 u       -> F5' parents (stratum 6)
#   * sigma_sweep 8   -> the (2,2,1^6) sigma statistics + L_inc(8) cross-check
# Then stratum 4/6 run over the concatenated parent files.
#
# Resumable per shard (.done markers). Gates: the sweep's L_inc(8) must sum to
# 46,778,967,018 exactly across shards.
#
#   RUST=/mnt/.../rust WORK=~/inc10 WORKERS=3 SHARDS=32 bash inc10_strata.sh
set -euo pipefail
RUST="${RUST:-$(cd "$(dirname "$0")" && pwd)}"
WORK="${WORK:-$HOME/inc10}"
WORKERS="${WORKERS:-3}"
SHARDS="${SHARDS:-64}"
LINK="/tmp/inc10_bshim.o -L/usr/local/lib -lnauty"
NAUTY_INC="$(dirname "$(find /usr/include /usr/local/include -name nauty.h 2>/dev/null | head -1)")"
mkdir -p "$WORK/strata"
cd "$WORK/strata"

gcc -O2 -c "$RUST/balanced_shim.c" -I"$NAUTY_INC" -o /tmp/inc10_bshim.o
rustc -O -C target-cpu=native "$RUST/inc10.rs"       -o /tmp/p_inc10  -C link-args="$LINK"
rustc -O -C target-cpu=native "$RUST/inc4.rs"        -o /tmp/p_inc4   -C link-args="$LINK"
rustc -O -C target-cpu=native "$RUST/f3x.rs"         -o /tmp/p_f3x    -C link-args="$LINK"
rustc -O -C target-cpu=native "$RUST/sigma_sweep.rs" -o /tmp/p_ssweep -C link-args="$LINK"
export SHARDS

seq 0 $((SHARDS-1)) | xargs -P"$WORKERS" -I@ bash -c '
  s=@
  [ -f shard_$s.done ] && exit 0
  fd="fifo_$s"; rm -rf "$fd"; mkdir "$fd"
  mkfifo "$fd/a" "$fd/b" "$fd/c" "$fd/d"
  /tmp/p_inc10 10 < "$fd/a" 2>/dev/null | grep -o "L_nullity2_labeled=[0-9]*" > s2_$s.txt & p1=$!
  /tmp/p_f3x 8 3 u < "$fd/b" 2>/dev/null > p3_$s.d6 & p2=$!
  /tmp/p_f3x 8 5 u < "$fd/c" 2>/dev/null > p5_$s.d6 & p3=$!
  /tmp/p_ssweep 8 < "$fd/d" 2>/dev/null | grep -E "^n=8 " > sw_$s.txt & p4=$!
  nauty-geng 8 $s/'"$SHARDS"' 2>/dev/null | nauty-directg -o 2>/dev/null \
    | tee "$fd/a" "$fd/b" "$fd/c" > "$fd/d"
  wait $p1 $p2 $p3 $p4
  rm -rf "$fd"
  # inc10 and ssweep print a summary line even on an empty stream, so an
  # empty output file means a consumer died mid-shard
  if [ ! -s s2_$s.txt ] || [ ! -s sw_$s.txt ]; then
    echo "shard $s: consumer output missing -- NOT marking done" >&2
    exit 1
  fi
  touch shard_$s.done
  echo "shard $s done ($(date +%H:%M:%S))"
'

echo "=== all shards done ==="
awk -F= '{s+=$2} END{print s}' s2_*.txt > "$WORK/inc10_L2.txt"
echo "L2 = $(cat "$WORK/inc10_L2.txt")"

# sigma sweep totals: sum raw accumulators, verify L_inc(8), apply the formula
python3 - <<'PY' | tee "$WORK/inc10_fix_2_2_1x6.txt"
import glob, re
acc = {"L_inc": 0, "rawJ1": 0, "rawJ0": 0, "rawK0b": 0, "rawH1": 0}
for f in glob.glob("sw_*.txt"):
    line = open(f).read()
    for k in acc:
        m = re.search(rf"{k}=(\d+)", line)
        if m:
            acc[k] += int(m.group(1))
assert acc["L_inc"] == 46778967018, f"sweep L_inc(8) GATE FAILED: {acc['L_inc']}"
raw = 2 * acc["rawJ1"] + acc["rawJ0"] + 4 * acc["rawH1"] + 2 * acc["rawK0b"]
assert raw % 56 == 0, "not integral"
print(raw // 56)
PY
echo "Fix_(2,2,1^6) = $(cat "$WORK/inc10_fix_2_2_1x6.txt")  [sweep L_inc(8) gate passed]"

# strata 4 and 6 over the concatenated (already-unique) parent streams
cat p3_*.d6 > parents_9_3.d6
cat p5_*.d6 > parents_9_5.d6
echo "parents: d=3 $(wc -l < parents_9_3.d6), d=5 $(wc -l < parents_9_5.d6)"
/tmp/p_inc4 10   < parents_9_3.d6 | grep -o "L_nullity4_labeled=[0-9]*" | cut -d= -f2 > "$WORK/inc10_L4.txt"
/tmp/p_inc4 10 6 < parents_9_5.d6 | grep -o "L_nullity6_labeled=[0-9]*" | cut -d= -f2 > "$WORK/inc10_L6.txt"
echo "L4 = $(cat "$WORK/inc10_L4.txt")"
echo "L6 = $(cat "$WORK/inc10_L6.txt")"
echo "L_inc(10) = $(( $(cat "$WORK/inc10_L2.txt") + $(cat "$WORK/inc10_L4.txt") + $(cat "$WORK/inc10_L6.txt") ))"
