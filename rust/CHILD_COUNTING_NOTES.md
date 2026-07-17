# Child-counting investigation: reducing the Σz enumeration in inclusive(10)

Fork investigation into whether the nullity-2 child enumeration (currently
Σz ≈ 3.5×10¹¹ weighted leaves at n=10, z̄ ≈ 3.7 constructions per child) can be
structurally reduced. **Outcome: a validated replacement identity (the
"two-sided endpoint-local rule") that leaves the enumeration count unchanged
but collapses the per-leaf cost and makes ~half the leaves mid-DFS prunable —
estimated 1.6–2.2× on the production engine, which combined with SMT reaches
the ~2h wall target. Prototype `rust/cc_proto.rs`, both anchors exact.**

## Baseline

`inc10.rs` counts labeled nullity-2 inclusive n-games via
L₂ = n·Σ_{(P,r)} 1/z(child): every valid (parent, extension) pair walks to a
leaf, computes the child's second equilibrium endpoint (λ*-argmax over supp
fractions — the dominant leaf cost), and contributes 1/z, z = |Z₁|+|Z₂|.
Anchors: n=6 → 126,900 (805 leaves), n=8 → 45,897,886,776 (4,306,838 leaves).
Engine is branch-bound at ~68ns/event with a near-optimal tree
(1.83 nodes/leaf).

## Avenues considered

### 1+2. Canonical-construction rules → the two-sided endpoint-local rule (VALIDATED)

Key observation: a construction (P, r) enters its child through endpoint
E₁ = (v, 0) — the parent kernel extended by 0 — and the deleted-vertex
competition within that endpoint's zero set is **parent-known**:
Z₁ = Z(v) ∪ {new}. The child's other endpoint E₂ never contains the new vertex
(its new-coordinate is 1).

**Rule:** accept (P, r) iff the new vertex attains the maximum of a
deterministic iso-invariant vertex function over Z₁ (computed in the child);
accepted constructions weigh 1/(2·T₁), T₁ = number of tied maxima.

**Identity proof:** each child C has endpoints p, q with disjoint zero sets.
Constructions of C via u ∈ Z_p all have Z₁ = Z_p; exactly the T₁ sig-argmax
members are accepted, contributing T₁·(1/(2T₁)) = 1/2. The Z_q side
contributes the other 1/2. Total 1 per labeled child. Ties need no
canonicalization: fractional weights absorb them exactly, and any
iso-invariant vertex function yields the same argmax set in every labeled
construction of the same child. (This is cleaner than construction
*deduplication*, which does need canon — here we only weight.)

**Why it wins:**
- **Endpoint-2 disappears.** No supp tracking rows, no fraction argmax — the
  single largest leaf cost in the production engine. z is never computed.
- **CM parents (Z(v) = ∅, 28.9% of all leaves at n=8) become trivial:**
  Z₁ = {new}, always accepted at weight 1/2; leaf = paradox check + accumulate.
- **The other 71% are mid-DFS prunable:** acceptance requires the new vertex to
  (weakly) degree-dominate every Z(v) vertex in the child, and both sides'
  child-degrees are linear in r — sound suffix bounds can kill non-canonical
  subtrees early (not implemented in the prototype; the acceptance rate is
  53.6%, so up to ~46% of leaves plus their subtrees are prunable).

**Prototype validation** (`cc_proto.rs`, derived from the validated
`inc_count.rs`, signature = (od, id, out/in neighbour-degree digests)):

    n=6: parents=127  walked=836      leaves(accepted)=499        L = 126900            EXACT
    n=8: parents=179477 walked=4365161 cm_walked=1260526
         accepted=2307428 (53.6% ≈ 2/z̄ = 54% predicted)           L = 45897886776       EXACT

Same totals from a completely different weighting — a strong dual validation
of both identities. Enumeration volume unchanged (walked ≈ old leaves + the
paradox-rejected), as predicted: the win is cheap leaves + prunability, not
fewer constructions.

### Endpoint-choice-only weighting (REJECTED, proof)

Choosing a canonical endpoint (e.g. larger zero set) and weighting 1/z_side
still enumerates every (P, r): whether a leaf's side is canonical depends on
|Z₂| (an exact-tie count) known only at the leaf, so nothing prunes. Verified
analytically; no reduction.

### 3. Different parent families (REJECTED)

- Min-z / one-sided-endpoint parents: same leaf-level z₂ dependence, no
  early pruning, and loses the parent-known competition set.
- Two-vertex extensions from (n−2)-grandparents directly to children:
  multiplicity becomes ≈ z·y ≈ 26 per child — strictly worse.

### 4. Algebraic per-parent closed forms (REJECTED for the typical case)

Σ_r 1/z per parent = N/(z₁+1) + tie corrections, and N (lattice points of
{r ∈ {−1,0,1}^p : v·r = 0, Ar < 0, paradox}) is computable by a
transfer-matrix DP over coordinates with state = partial row sums. But the
DP's state count (≈ range of partial sums × levels ≈ thousands) exceeds the
DFS's ~180 events for a typical parent (64 leaves), and the tie corrections
need one extra equality-row DP per supp pair (~28 more). Only a heavy-tail
hybrid (DP for parents with ≥10³–10⁴ leaves) could pay; tail mass unmeasured.
The paradox conditions do force coordinates (no-win vertex ⇒ r_i = +1,
no-loss ⇒ r_i = −1) — worth fixing before the DFS in the production engine
regardless (currently walked 3-ways and leaf-rejected).

## Recommendation

Integrate the two-sided rule into `inc10.rs`:
1. Delete the endpoint-2 machinery (supp_v, supp_rows, fraction argmax, z).
2. Leaf: paradox masks → if Z(v) = ∅: accept, weight (n−1)!/|Aut|·LCM/2.
   Else compute child vertex signatures once (shared degree arrays), argmax
   over Z(v) ∪ {new}, weight LCM/(2T₁) on acceptance. LCM = lcm(1..2n)
   already covers 2T₁ ≤ 2n.
3. Add the mid-DFS degree-domination prune for Z(v) ≠ ∅ parents (partial
   (od,id) of new and each Z-vertex are linear in r; suffix bounds sound).
   [DONE: prune on cap = nod + remaining vs each Z-vertex's minimum child
   out-degree, with an od-tie in-degree refinement and a cap ≤ zmaxpod+1 gate
   that skips the loop where nothing can fire. Kills 29% of surviving leaves
   at n=8 (4.31M → 3.07M) and 59% of walked leaves on the dense n=10 band;
   anchors and accepted counts unchanged.]
4. Fix paradox-forced coordinates before the DFS.
5. Re-run the full anchor ladder (126,900 / 45,897,886,776) — mandatory.

Expected: 1.6–2.2× over the current engine (leaf cost collapse + pruning),
taking stratum-2 from ~5h to ~2.3–3h wall on 4 cores; with WSL processors=8
(SMT on the same silicon), ~1.8–2.4h — at the 2-hour target. The identity
change is validated; the remaining risk is only implementation, gated as
always by the anchors.
