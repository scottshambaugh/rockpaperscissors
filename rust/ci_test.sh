#!/usr/bin/env bash
# Build the nauty-FFI census tools and check them against known small-n counts.
# Self-contained: locates nauty.h, compiles the C shim + both Rust binaries, runs
# them on n where the (total, twin-free, prime) answers are known, and fails on any
# mismatch. Build artifacts go to /tmp/ci_* so a local run never touches binaries an
# overnight job may be using. Needs: gcc, rustc, libnauty2-dev.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# nauty.h from the distro package (/usr/include/...) or a source build (/usr/local/include)
NAUTY_H="$(find /usr/include /usr/local/include -name nauty.h 2>/dev/null | head -1 || true)"
[ -n "$NAUTY_H" ] || { echo "nauty.h not found — install libnauty2-dev or build nauty"; exit 1; }
NAUTY_INC="$(dirname "$NAUTY_H")"
echo "nauty headers: $NAUTY_INC"

LINK="/tmp/ci_bshim.o -L/usr/local/lib -lnauty"  # -L covers a source-built libnauty.a
gcc -O2 -c rust/balanced_shim.c -I"$NAUTY_INC" -o /tmp/ci_bshim.o
rustc -O rust/balanced.rs -o /tmp/ci_balanced -C link-args="$LINK"
rustc -O rust/regular.rs  -o /tmp/ci_regular  -C link-args="$LINK"

fail=0
expect () {  # $1 label, $2 expected substring, $3.. command
  local label="$1" expected="$2"; shift 2
  local out; out="$("$@" 2>/dev/null)"
  if echo "$out" | grep -qF "$expected"; then
    echo "ok   $label"
  else
    echo "FAIL $label: got [$out] want [$expected]"; fail=1
  fi
}

# balanced (A308239 totals; twin-free / prime from the brute-force range)
expect "balanced n=5" "total=4 twin_free=3 prime=3"          /tmp/ci_balanced 5 3 1 0
expect "balanced n=6" "total=16 twin_free=13 prime=13"       /tmp/ci_balanced 6 3 1 0
expect "balanced n=7" "total=175 twin_free=152 prime=152"    /tmp/ci_balanced 7 5 1 0
expect "balanced n=8" "total=5274 twin_free=4921 prime=4917" /tmp/ci_balanced 8 5 1 0
# regular (strict tier, summed over degree strata)
expect "regular n=6"  "total=5 twin_free=4 prime=4"          /tmp/ci_regular  6 4 1 0
expect "regular n=7"  "total=13 twin_free=12 prime=12"       /tmp/ci_regular  7 4 1 0
expect "regular n=8"  "total=82 twin_free=76 prime=76"       /tmp/ci_regular  8 5 1 0
expect "regular n=9"  "total=2016 twin_free=1973 prime=1972" /tmp/ci_regular  9 5 1 0
# sharding must partition: 3 shards of balanced n=7 sum to the whole
s0=$(/tmp/ci_balanced 7 5 3 0 2>/dev/null); s1=$(/tmp/ci_balanced 7 5 3 1 2>/dev/null); s2=$(/tmp/ci_balanced 7 5 3 2 2>/dev/null)
tot=$(printf '%s\n%s\n%s\n' "$s0" "$s1" "$s2" | awk '{for(i=1;i<=NF;i++)if($i~/^total=/){split($i,a,"=");t+=a[2]}}END{print t}')
[ "$tot" = "175" ] && echo "ok   balanced n=7 shard-sum=175" || { echo "FAIL shard-sum: got $tot want 175"; fail=1; }

[ $fail -eq 0 ] && echo "ALL RUST CHECKS PASSED" || { echo "RUST CHECKS FAILED"; exit 1; }
