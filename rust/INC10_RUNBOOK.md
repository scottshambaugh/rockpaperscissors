# inclusive(10) production runbook

Everything below is anchor-gated: every engine reproduces its n=8 (or smaller)
ground truth exactly before its n=10 output is trusted, and the assembled n=8
pipeline must give exactly inclusive(8) = 1,198,013 before launch.

## The identity

    inclusive(10) = ( L2 + L4 + L6 + SIGMA ) / 10!

where L_D = labeled count of nullity-D inclusive 10-vertex games and
SIGMA = sum over the 41 nontrivial cycle types of S_10 of #sigma x Fix_sigma.
The identity term L_inc(10) = L2 + L4 + L6.

## Prerequisites (ALL COMPUTED AND CONFIRMED)

| quantity | value | source |
|---|---|---|
| L_inc(2) | 0 | no paradoxical 2-vertex game |
| L_inc(3..6) | 2 / 42 / 978 / 130,950 | sigma_fix brute = strata sums (3 routes agree at n=6) |
| L_inc(7) | **49,473,198** | n=7 census (incx): CM 36,273,536 + hi 13,199,662; = blind gate prediction |
| L_inc(8) | 46,778,967,018 | strata: 45,897,886,776 + 880,869,360 + 210,882 |
| L_inc(9) | **235,837,146,265,362** | labsum: CM 211,720,417,352,944 + hi 24,116,728,912,418 |
| n=8 sigma gate | corrections = 1,524,917,142 | 20 brute + Fix_(2,1^6)=L_inc(7) (coprime rule) |
| **n=8 assembly** | **(46,778,967,018 + 1,524,917,142)/8! = 1,198,013 EXACT** | **GATE PASSED — validates the entire pipeline** |

## n=10 sigma type strategy (41 nontrivial types)

- **16 coprime types**: Fix = 3^(within) * L_inc(m), all L_inc known. Instant.
  Notable: (2,1^8) = L_inc(9); (3,1^7) = 3*L_inc(8); (4,1^6) = (3,2,1^5) =
  3*L_inc(7); (9,1) = (7,3) = 3^4*L_inc(2) = 0.
- **18 small brutes** (<= 15 bundle slots, <= 3^15 leaves): sigma_fix inline.
- **7 heavy types** (shard sigma_fix `TYPE s OF`, OF a power of 3, or sweep):

  | type | bundle slots | route |
  |---|---|---|
  | 4,2,1^4 | 16 | shard 3^6 (minutes) |
  | 3,2,2,1^3 | 16 | shard 3^6 |
  | 3,3,1^4 | 17 | shard 3^6 |
  | 2,2,2,2,2 | 20 | shard 3^9 (~1-2h sharded) |
  | 2,2,2,2,1,1 | 21 | shard 3^9 |
  | 2,2,2,1^4 | 24 | shard 3^12 (~hours; pruning cuts 3^24 hard) |
  | **2,2,1^6** | **29** | **sigma_sweep (3^29 infeasible to brute; sweep anchored at n=6)** |

  Pruning caveat: row-boundary + sign prune cut explored leaves far below
  3^slots at n=8 (e.g. (2,2,1^4): 3^16 raw -> 11M explored -> 355k balanced).
  Calibrate each heavy shard on a single s before launching the fan-out.

## Phase 1 -- strata (the 8-vertex grandparent stream, ~575M classes)

The stream `nauty-geng 8 | nauty-directg -o` is regenerated per consumer
(generation is cheap relative to the consumers). Shard with geng res/mod for
parallelism; shards are unbalanced, so use >= 32 shards under xargs -P4.

1. **Stratum 2** (~12 core-h, ~3h wall): per shard s:
   `nauty-geng 8 s/32 | nauty-directg -o | inc10 10` and sum the
   L_nullity2_labeled values. (Engine: fused grandparent-adjugate, two-sided
   endpoint rule, mid-DFS degree-domination prune; anchors 126,900 /
   45,897,886,776 exact.)
2. **Stratum 4 parents**: `nauty-geng 8 s/32 | nauty-directg -o | f3x 8 3`
   concatenated, then `nauty-labelg | sort -S 2G -T <scratch> -u
   > parents_9_3.d6`. (f3x anchor: reproduces the inc_strata f3-emit streams
   byte-identically at n=7.) Volume estimate ~100-300M raw records / a few GB.
3. **Stratum 4**: split parents_9_3.d6 into 4 chunks, `inc4 10 4` each, sum.
   (Calibrate first on 100k parents; n=8 rate was ~200us/parent.)
4. **Stratum 6 parents**: same with `f3x 8 5` -> parents_9_5.d6 (tiny).
5. **Stratum 6**: `inc4 10 6 < parents_9_5.d6` (minutes).

## Phase 2 -- sigma corrections

1. **Sweep types**: one pass `nauty-geng 8 s/32 | nauty-directg -o |
   sigma_sweep 8` per shard; sum the raw accumulators, divide by 56 at the
   end. Gives Fix_(2,2,1^6) = 2*J1 + J0 + 4*H1 + 2*K0b AND the L_inc(8)
   cross-check (must equal 46,778,967,018 exactly).
   (Sweep anchor: reproduces brute Fix_(2,2,1^4) = 339,138 at n=6.)
2. **Heavy brutes** (sigma_fix, shard arg `s OF` with OF a power of 3):
   - `2,2,2,1,1,1,1` (24 slots): 27-way shard, sum. The big one (~10 core-h).
   - `2,2,2,2,1,1` (21 slots): 9-way shard.
   - `2,2,2,2,2` (20 slots): 9-way shard.
   Shard-sum identity is validated (9-way = unsharded at n=8).
3. **Everything else**: `sigma_gate.py 10 --linc ... --fix ...` computes
   coprime types via Fix = 3^k * L_inc(m) (validated rule, m >= 2) and small
   brutes inline, then assembles.

## Phase 3 -- assembly

    sigma_gate.py 10 \
      --linc 2=0,3=2,4=42,5=978,6=130950,7=49473198,8=46778967018,9=235837146265362,10=<L2+L4+L6> \
      --fix 2,2,1,1,1,1,1,1=<sweep> --fix 2,2,2,1,1,1,1=<shardsum> \
      --fix 2,2,2,2,1,1=<shardsum> --fix 2,2,2,2,2=<shardsum> \
      --fix 4,2,1,1,1,1=<shardsum> --fix 3,3,1,1,1,1=<shardsum> \
      --fix 3,2,2,1,1,1=<shardsum>

The driver computes the 16 coprime types and 18 small brutes itself; only the
7 heavy types need --fix. L_inc(10) = L2+L4+L6 is the identity term (passed as
--linc 10=). Final gate: total % 10! == 0 AND matches the divisibility check.

Final integrality check (total % 10! == 0) is a real gate: any missed or
double-counted term breaks divisibility with probability ~1.

## Budget

| item | core-h |
|---|---|
| stratum 2 | ~12 |
| stratum 4 (f3x + sort + inc4) | ~4-7 |
| stratum 6 | ~0.3 |
| sigma sweep (stream-bound) | ~3 |
| sigma heavy brutes | ~12 |
| L_inc(9) labsum | ~1 |
| **total** | **~32-35 core-h, ~8-9h wall on 4 cores** |

## Post-run (separate, optional)

The (twin-free) [prime] brackets for inclusive(10) via the bracket-complement
method: non-tf inclusive = blow-ups of inclusive cores (the weight condition
folds into w, so plain inclusive cores suffice -- same simplification as the
sigma quotients); prime gap via module substitution. Not part of this run.
