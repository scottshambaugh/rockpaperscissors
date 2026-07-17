#!/usr/bin/env bash
# Phase 1 of the inclusive(10) production run: the three nullity strata over
# the 8-vertex grandparent stream. Emits L2, L4, L6 to inc10_L{2,4,6}.txt.
# Sharded by geng res/mod; WORKERS parallel, leaving a core spare on a 4-core
# box. Resumable: each shard writes a .done marker.
#
#   RUST=/mnt/.../rust WORK=~/inc10 WORKERS=3 SHARDS=32 bash inc10_strata.sh
set -euo pipefail
RUST="${RUST:-$(cd "$(dirname "$0")" && pwd)}"
WORK="${WORK:-$HOME/inc10}"
WORKERS="${WORKERS:-3}"
SHARDS="${SHARDS:-32}"
LINK="/tmp/inc10_bshim.o -L/usr/local/lib -lnauty"
NAUTY_INC="$(dirname "$(find /usr/include /usr/local/include -name nauty.h 2>/dev/null | head -1)")"
mkdir -p "$WORK/strata"
cd "$WORK/strata"

gcc -O2 -c "$RUST/balanced_shim.c" -I"$NAUTY_INC" -o /tmp/inc10_bshim.o
rustc -O -C target-cpu=native "$RUST/inc10.rs" -o /tmp/p_inc10 -C link-args="$LINK"
rustc -O -C target-cpu=native "$RUST/inc4.rs"  -o /tmp/p_inc4  -C link-args="$LINK"
rustc -O -C target-cpu=native "$RUST/f3x.rs"   -o /tmp/p_f3x
export SHARDS

# --- Stratum 2: sum L_nullity2_labeled directly from the grandparent stream ---
seq 0 $((SHARDS-1)) | xargs -P"$WORKERS" -I@ sh -c '
  s=@
  [ -f s2_$s.done ] && exit 0
  nauty-geng 8 $s/'"$SHARDS"' 2>/dev/null | nauty-directg -o 2>/dev/null \
    | /tmp/p_inc10 10 2>/dev/null | grep -o "L_nullity2_labeled=[0-9]*" > s2_$s.txt
  touch s2_$s.done
  echo "stratum2 shard $s done ($(date +%H:%M:%S))"
'
awk -F= '{s+=$2} END{print s}' s2_*.txt > "$WORK/inc10_L2.txt"
echo "L2 = $(cat "$WORK/inc10_L2.txt")"

# --- Stratum 4 & 6 parents: f3x per shard, then global labelg | sort -u ---
for D in 3 5; do
  seq 0 $((SHARDS-1)) | xargs -P"$WORKERS" -I@ sh -c '
    s=@; D='"$D"'
    [ -f p${D}_$s.done ] && exit 0
    nauty-geng 8 $s/'"$SHARDS"' 2>/dev/null | nauty-directg -o 2>/dev/null \
      | /tmp/p_f3x 8 $D 2>/dev/null > p${D}_$s.d6
    touch p${D}_$s.done
    echo "f3x d=$D shard $s: $(wc -c < p${D}_$s.d6) bytes ($(date +%H:%M:%S))"
  '
  LC_ALL=C sort -S 2G -T "$WORK/strata" -m -u \
    <(cat p${D}_*.d6 | nauty-labelg 2>/dev/null | LC_ALL=C sort -S 2G -T "$WORK/strata" -u) \
    > parents_9_$D.d6
  echo "parents d=$D: $(wc -l < parents_9_$D.d6) classes"
done

# --- Stratum 4 and 6 counts ---
/tmp/p_inc4 10   < parents_9_3.d6 | grep -o "L_nullity4_labeled=[0-9]*" | cut -d= -f2 > "$WORK/inc10_L4.txt"
/tmp/p_inc4 10 6 < parents_9_5.d6 | grep -o "L_nullity6_labeled=[0-9]*" | cut -d= -f2 > "$WORK/inc10_L6.txt"
echo "L4 = $(cat "$WORK/inc10_L4.txt")"
echo "L6 = $(cat "$WORK/inc10_L6.txt")"
echo "L_inc(10) = L2+L4+L6 = $(( $(cat "$WORK/inc10_L2.txt") + $(cat "$WORK/inc10_L4.txt") + $(cat "$WORK/inc10_L6.txt") ))"
