#!/bin/bash
# Two-paradox (S2 / Erdos-Schutte) tournament count via nauty + the Rust S2 filter.
#
#   nauty/two_paradox.sh 11 [mod]
#
# nauty-gentourng generates one tournament per isomorphism class; rust/s2_filter
# applies the S2 predicate (every pair of vertices has a common dominator). The
# unfiltered total printed per shard must sum to A000568(n) -- a built-in
# completeness checksum. `mod` (default = nproc) splits the work across cores via
# gentourng's res/mod; the S2 counts are summed. NB: do NOT use gentourng -c
# (strongly-connected); 5/226 of the n=9 S2 tournaments are not strongly
# connected, so -c undercounts.
set -euo pipefail
N=${1:?usage: two_paradox.sh n [mod]}
MOD=${2:-$(nproc)}
FILTER=${S2FILTER:-/tmp/s2filter}
[ -x "$FILTER" ] || { rustc -O "$(dirname "$0")/../rust/s2_filter.rs" -o "$FILTER"; }

tmp=$(mktemp -d)
pids=()
for r in $(seq 0 $((MOD-1))); do
  ( nauty-gentourng "$N" "$r/$MOD" 2>/dev/null | "$FILTER" "$N" > "$tmp/shard_$r" ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done
awk '{for(i=1;i<=NF;i++){if($i~/^total=/){split($i,a,"=");t+=a[2]}
       if($i~/^S2/){split($i,b,"=");s+=b[2]}}}
     END{printf "n=%s: tournaments=%d  S2(two-paradox)=%d\n", N, t, s}' N="$N" "$tmp"/shard_*
rm -rf "$tmp"
