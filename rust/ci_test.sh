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
[ -n "$NAUTY_H" ] || { echo "nauty.h not found — sudo apt install libnauty-dev (see rust/README.md)"; exit 1; }
NAUTY_INC="$(dirname "$NAUTY_H")"
echo "nauty headers: $NAUTY_INC"

LINK="/tmp/ci_bshim.o -L/usr/local/lib -lnauty"  # -L covers a source-built libnauty.a
gcc -O2 -c rust/balanced_shim.c -I"$NAUTY_INC" -o /tmp/ci_bshim.o
rustc -O rust/balanced.rs -o /tmp/ci_balanced -C link-args="$LINK"
rustc -O rust/regular.rs  -o /tmp/ci_regular  -C link-args="$LINK"
rustc -O rust/wbal.rs     -o /tmp/ci_wbal     -C link-args="$LINK"
rustc -O rust/cm_filter.rs -o /tmp/ci_cmf
rustc -O rust/factorcrit.rs -o /tmp/ci_fc
rustc -O rust/cm_extend.rs -o /tmp/ci_cmx
rustc -O rust/prime_filter.rs -o /tmp/ci_prime
rustc -O rust/inc_fast.rs -o /tmp/ci_incf
rustc -O rust/inc_extend.rs -o /tmp/ci_incx
rustc -O rust/inc10.rs -o /tmp/ci_inc10 -C link-args="$LINK"
rustc -O rust/inc4.rs -o /tmp/ci_inc4 -C link-args="$LINK"
rustc -O rust/inc_strata.rs -o /tmp/ci_incs -C link-args="$LINK"
rustc -O -C overflow-checks=on rust/burnside_regular.rs -o /tmp/ci_burn

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

# weighted-balanced core enumerator (the balanced-bracket complement method):
# identity weights reproduce tf-balanced(5); the blow-up families reproduce the
# non-twin-free counts (bal(5): 4 total - 3 tf = 1, from the (2,1,1,1) family)
expect "wbal identity n=5"  "cores=3" /tmp/ci_wbal 1,1,1,1,1
expect "wbal blowup (2,1,1,1)" "cores=1" /tmp/ci_wbal 2,1,1,1
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
# completely mixed (Kaplansky): Pfaffian-cofactor filter over nauty-directg.
# n=5/n=7 anchor to the Python census; candidates anchor to A001174. The
# factor-critical prefilter (fc, used by the big n=9 run) must lose no CM game.
if command -v nauty-geng >/dev/null && command -v nauty-directg >/dev/null; then
  cmf ()  { nauty-geng "$1" 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/ci_cmf "$1"; }
  cmff () { nauty-geng "$1" 2>/dev/null | /tmp/ci_fc "$1" 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/ci_cmf "$1"; }
  expect "compl-mixed n=5"          "candidates=582 completely_mixed=7"        cmf  5
  expect "compl-mixed n=7"          "candidates=2142288 completely_mixed=7268" cmf  7
  expect "compl-mixed n=5 +fc"      "completely_mixed=7"                       cmff 5
  expect "compl-mixed n=7 +fc"      "completely_mixed=7268"                    cmff 7  # fc lossless
  # extension method (cm_extend): build n-vertex CM games from (n-1)-vertex
  # parents, canonicalize + dedup. An INDEPENDENT algorithm from the cm_filter
  # scan above -- agreement on 7 / 7268 cross-validates both. Needs labelg.
  if command -v nauty-labelg >/dev/null; then
    cmx () { nauty-geng "$1" 2>/dev/null | nauty-directg -o 2>/dev/null \
             | /tmp/ci_cmx "$2" 2>/dev/null | nauty-labelg 2>/dev/null | sort -u | wc -l; }
    expect "cm-extend n=5 (from 4)"  "7"    cmx 4 5
    expect "cm-extend n=7 (from 6)"  "7268" cmx 6 7
    # prime subcount of the CM games (rust/prime_filter.rs, reusing balanced.rs's
    # is_prime); must match the Python census [prime] brackets 6 and 7240.
    cmprime () { nauty-geng "$1" 2>/dev/null | nauty-directg -o 2>/dev/null \
                 | /tmp/ci_cmx "$2" 2>/dev/null | nauty-labelg 2>/dev/null | sort -u | /tmp/ci_prime "$2"; }
    expect "cm-prime n=5"  "total=7 prime=6"       cmprime 4 5
    expect "cm-prime n=7"  "total=7268 prime=7240" cmprime 6 7
    # inclusive: the fast Pfaffian-classified filter and the extension method are
    # independent algorithms; both must give the census 15 / 10525, and the
    # extension's CM sub-stream must reproduce the cm_extend counts exactly.
    incf () { nauty-geng "$1" 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/ci_incf "$1"; }
    incx () { nauty-geng "$1" 2>/dev/null | nauty-directg -o 2>/dev/null \
              | /tmp/ci_incx "$2" 2>/dev/null | nauty-labelg 2>/dev/null | sort -u | wc -l; }
    expect "inclusive-fast n=5" "inclusive=15 (cm/nullity1=7"       incf 5
    expect "inclusive-fast n=7" "inclusive=10525 (cm/nullity1=7268" incf 7
    expect "inc-extend n=5 (from 4)" "15"    incx 4 5
    expect "inc-extend n=7 (from 6)" "10525" incx 6 7
    # inclusive(10) stratum engines: labeled nullity-2 sums via the fused
    # grandparent-stream engine (two-sided endpoint rule) and labeled nullity-4
    # sums via the Motzkin-fused counter; anchors are the n=6/n=8 ground truth
    inc10a () { nauty-geng "$1" 2>/dev/null | nauty-directg -o 2>/dev/null | /tmp/ci_inc10 "$2"; }
    inc4a ()  { nauty-geng "$1" 2>/dev/null | nauty-directg -o 2>/dev/null \
                | /tmp/ci_incs "$1" f3-emit 3 2>/dev/null | /tmp/ci_inc4 "$2"; }
    expect "inc10 nullity-2 n=6" "L_nullity2_labeled=126900"      inc10a 4 6
    expect "inc10 nullity-2 n=8" "L_nullity2_labeled=45897886776" inc10a 6 8
    expect "inc4 nullity-4 n=6"  "L_nullity4_labeled=4050"        inc4a 5 6
    # split-stream mode: cm to stdout (7268), nullity>=3 to file (3257)
    hi7="$(mktemp)"
    cmside="$(nauty-geng 6 2>/dev/null | nauty-directg -o 2>/dev/null \
              | /tmp/ci_incx 7 "$hi7" 2>/dev/null | nauty-labelg 2>/dev/null | sort -u | wc -l)"
    hiside="$(nauty-labelg 2>/dev/null < "$hi7" | sort -u | wc -l)"
    if [ "$cmside" = "7268" ] && [ "$hiside" = "3257" ]; then
      echo "ok   inc-extend n=7 split-stream cm=7268 hi=3257"
    else
      echo "FAIL inc-extend split-stream: cm=$cmside hi=$hiside want 7268/3257"; fail=1
    fi
  else
    echo "skip cm-extend checks (nauty-labelg not on PATH)"
  fi
else
  echo "skip completely-mixed checks (nauty-geng/nauty-directg not on PATH)"
fi

# burnside counting method (no enumeration): totals must reproduce the known
# regular column through n=10, per-stratum (A096368: 15 regular tournaments at n=9)
expect "burnside regular n=10" "regular(n) totals: [0, 0, 1, 1, 2, 5, 13, 82, 2016, 154831]" /tmp/ci_burn 10
expect "burnside stratum d=4"  "d=4: connected iso by n: [0, 0, 0, 0, 0, 0, 0, 0, 15, 3987]" /tmp/ci_burn 10

# sharding must partition: 3 shards of balanced n=7 sum to the whole
s0=$(/tmp/ci_balanced 7 5 3 0 2>/dev/null); s1=$(/tmp/ci_balanced 7 5 3 1 2>/dev/null); s2=$(/tmp/ci_balanced 7 5 3 2 2>/dev/null)
tot=$(printf '%s\n%s\n%s\n' "$s0" "$s1" "$s2" | awk '{for(i=1;i<=NF;i++)if($i~/^total=/){split($i,a,"=");t+=a[2]}}END{print t}')
[ "$tot" = "175" ] && echo "ok   balanced n=7 shard-sum=175" || { echo "FAIL shard-sum: got $tot want 175"; fail=1; }

[ $fail -eq 0 ] && echo "ALL RUST CHECKS PASSED" || { echo "RUST CHECKS FAILED"; exit 1; }
