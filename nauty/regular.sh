#!/bin/bash
# regular(n) via nauty -- isomorphism classes of connected oriented graphs in
# which every vertex has out-degree = in-degree = d, summed over d.
#
#   nauty/regular.sh 10
#
# For each d, the decisive edges form a connected 2d-regular graph (geng), and a
# regular orientation is an Eulerian orientation with out=in=d -- exactly what
# `watercluster2 i$d o$d S` produces (max in- and out-degree <= d on a 2d-regular
# graph forces in=out=d). watercluster2 emits one representative per isomorphism
# class, so the line count is the class count. Validated against the Python
# `search_regular`: 2, 5, 13, 82, 2016 at n = 5..9.
set -euo pipefail
N=${1:?usage: regular.sh n}
dmax=$(( (N - 1) / 2 ))
total=0
brk=""
for d in $(seq 1 "$dmax"); do
  c=$(nauty-geng "$N" -d$((2*d)) -D$((2*d)) -c 2>/dev/null \
        | nauty-watercluster2 "i$d" "o$d" S Z 2>/dev/null | wc -l)
  total=$((total + c))
  brk="$brk d$d:$c"
done
echo "regular($N) = $total  [$brk ]"
