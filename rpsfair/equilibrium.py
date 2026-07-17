"""Symmetric zero-sum equilibrium computation.

The payoff matrix A = M is antisymmetric (A = -Aᵀ), so the game value is 0 and
pᵀAp = 0 for every p. The equilibrium set is the polytope

    O = { p in the simplex : A p <= 0 }

and by complementary slackness any *fully-mixed* (strictly interior) equilibrium
satisfies A p = 0 exactly -- it lives in ker(A). Since rank(A) is always even,
dim ker(A) has the same parity as n: a unique interior equilibrium (dim ker = 1)
is only possible at odd n. At even n a fully-mixed equilibrium is never isolated;
it sits in a family of dimension dim ker(A) - 1 >= 1 (e.g. the n=4 cop game, whose
equilibrium set is a 1-D segment). So "the equilibrium odds" are only well-defined
once we pick a canonical representative. We use the leximin (max-min) point of O
(`maxmin_equilibrium`): maximise the smallest play rate, then the next, ... . It is
unique, lands on clean rationals for these games (cop game -> 0.2/0.2/0.2/0.4), and
coincides with the usual answer whenever O is a single point. `max_entropy_equilibrium`
offers a structure-weighted alternative.
"""

from itertools import combinations

import numpy as np
from scipy.optimize import linprog, minimize


def equilibrium(M):
    """Return the canonical leximin symmetric Nash equilibrium.

    This is the backwards-compatible entry point for callers that want one
    equilibrium. It delegates to `maxmin_equilibrium`, which enforces the actual
    Nash polytope constraints p >= 0, sum(p) = 1, and M p <= 0. In particular it
    remains valid when no normalized kernel vector exists and the equilibrium is
    on the simplex boundary.
    """
    return maxmin_equilibrium(M)


def kernel_dim(M, tol=1e-8):
    """dim ker(A) = n - rank(A). Equals the parity of n (rank is always even)."""
    s = np.linalg.svd(np.asarray(M, dtype=float), compute_uv=False)
    return int((s < tol * max(s[0], 1.0)).sum())


def has_fully_mixed(M, tol=1e-8):
    """Return (is_fully_mixed, witness_distribution_or_None).

    Tests whether a *fully-mixed* equilibrium EXISTS, i.e. whether ker(A) contains
    a strictly positive vector (equivalently O has an interior point). This is the
    correct existence test -- it does not depend on any single returned solution.

    The witness is a particular interior point (non-canonical when dim ker(A) >= 2);
    do NOT read odds/metrics off it -- use `maxmin_equilibrium` instead.
    """
    n = len(M)
    _, S, Vt = np.linalg.svd(M.astype(float))
    null = Vt[tol * max(S[0], 1.0) > S]
    if len(null) == 0:
        return False, None
    if len(null) == 1:
        v = null[0]
        if (v > tol).all():
            return True, v / v.sum()
        if (v < -tol).all():
            return True, (-v) / (-v).sum()
        return False, None
    # dim ker(A) >= 2: the on-kernel system is rank-deficient, so test positivity
    # with an LP (a single linear solve would silently pick one face). We maximise
    # the minimum coordinate; > 0 means a strictly interior point exists.
    k = len(null)
    A_ub = np.zeros((n, k + 1))
    A_ub[:, :k] = -null.T
    A_ub[:, k] = 1
    A_eq = np.zeros((1, k + 1))
    A_eq[0, :k] = null.sum(1)
    c = np.zeros(k + 1)
    c[k] = -1
    r = linprog(
        c,
        A_ub=A_ub,
        b_ub=np.zeros(n),
        A_eq=A_eq,
        b_eq=[1.0],
        bounds=[(None, None)] * (k + 1),
        method="highs",
    )
    if r.success and -r.fun > tol:
        x = null.T @ r.x[:k]
        return True, x / x.sum()
    return False, None


def is_completely_mixed(M, tol=1e-8):
    """Whether EVERY Nash equilibrium is fully mixed (a *completely mixed game*).

    Kaplansky (1945) named these: a matrix game is *completely mixed* when every
    optimal strategy (of both players) is completely mixed, and proved this forces
    the optimal strategy to be UNIQUE. For skew-symmetric games it also forces odd
    n (the parity theorem: rank is even, so at even n a fully-mixed equilibrium
    always sits in a continuum and the continuum's boundary points are equilibria
    that drop a strategy). Kaplansky (1995) characterised the skew-symmetric case
    via principal Pfaffians; operationally the test is just:

        dim ker(M) = 1  and  the kernel vector is strictly one-signed.

    That this implies "every equilibrium is fully mixed" (not merely "a fully-mixed
    equilibrium exists") is a two-line argument: let p > 0 span ker(M) and let q be
    any symmetric equilibrium (M q <= 0). Then pᵀ(M q) = -(M p)ᵀq = 0, and a
    positive combination of non-positive terms vanishes only if all do, so
    M q = 0 and q ∈ ker(M) = span(p), i.e. q = p. Conversely nullity >= 2 gives a
    fully-mixed *family* whose closure hits the simplex boundary in equilibria of
    partial support. So completely mixed <=> inclusive with a unique equilibrium,
    the strict top of the equilibrium-fairness ladder: every strategy is not just
    *playable* (inclusive) but *required* -- played in every equilibrium.
    """
    _, S, Vt = np.linalg.svd(np.asarray(M, dtype=float))
    null = Vt[tol * max(S[0], 1.0) > S]
    if len(null) != 1:
        return False
    v = null[0]
    return bool((v > tol).all() or (v < -tol).all())


def required_strategies(M, tol=1e-7):
    """Indices of strategies played in EVERY equilibrium (positive at every vertex of O).

    The complement of the "droppable" moves: strategy i is required iff
    min_{p in O} p_i > 0, and the minimum of a linear function over the polytope O
    is attained at a vertex, so we just scan the extreme equilibria. Contrast with
    the *essential* strategies (played in at least one equilibrium, the union of
    vertex supports): inclusive means every strategy is essential, completely mixed
    means every strategy is required. E.g. the 5-move elemental game (twins of the
    cop game) is inclusive, but Fire and Water are capped and only Grass+Clay's
    *combined* weight is bounded below -- its required set is a proper subset.
    """
    V = equilibrium_vertices(M, tol)
    if not V:
        return []
    P = np.array(V)
    return [i for i in range(len(M)) if (P[:, i] > tol).all()]


def maxmin_equilibrium(M, tol=1e-7):
    """Canonical equilibrium: the leximin point of O = {p simplex : A p <= 0}.

    Repeatedly maximises the smallest play rate, then the next-smallest, ... .
    This is the "max-min" criterion done canonically: a single max-min LP can
    leave a flat optimal face (re-introducing an arbitrary choice), so we fix
    each coordinate that is forced to the current minimum and recurse on the
    rest. The result is unique and lands on clean rationals for these {-1,0,+1}
    games -- e.g. the cop game gives exactly (0.2, 0.2, 0.2, 0.4).

    A single max-min LP can leave a flat optimal face, so at each stage we (a)
    compute the current minimum level t*, (b) pin every free coordinate that cannot
    exceed t* while the others are held at the floor t*, then (c) recurse on the rest
    with t* as a hard floor. Pinning against the *numeric* floor t* (not a free LP
    variable) is what makes it balance every independent group at its own scale --
    two unrelated free directions each get evened out rather than one absorbing all
    the slack.

    For a fully-mixed game the optimum is interior (A p = 0); when no fully-mixed
    equilibrium exists it returns the leximin boundary equilibrium. Returns None
    only if O is empty (never happens for a connected paradoxical game).
    """
    M = np.asarray(M, dtype=float)
    n = len(M)
    fixed = [None] * n  # coords pinned to their leximin value

    def pins():
        rows, vals = [np.ones(n)], [1.0]  # sum p = 1
        for i in range(n):
            if fixed[i] is not None:
                row = np.zeros(n)
                row[i] = 1.0
                rows.append(row)
                vals.append(fixed[i])
        return np.array(rows), np.array(vals)

    def max_min_level(free):
        # maximise t over { p in O, sum p = 1, pins, p_j >= t for every free j }
        c = np.zeros(n + 1)
        c[n] = -1.0
        A_ub = [np.hstack([M, np.zeros((n, 1))])]
        b_ub = [np.zeros(n)]
        for j in free:  # t - p_j <= 0
            row = np.zeros(n + 1)
            row[j] = -1.0
            row[n] = 1.0
            A_ub.append(row[None])
            b_ub.append([0.0])
        a_eq, b_eq = pins()
        a_eq = np.hstack([a_eq, np.zeros((len(a_eq), 1))])  # no t in the equalities
        return linprog(
            c,
            A_ub=np.vstack(A_ub),
            b_ub=np.concatenate(b_ub),
            A_eq=a_eq,
            b_eq=b_eq,
            bounds=[(0.0, 1.0)] * n + [(None, None)],
            method="highs",
        )

    def max_coord(j, free, floor):
        # maximise p_j over { p in O, sum p = 1, pins, every OTHER free coord >= floor }
        c = np.zeros(n)
        c[j] = -1.0
        a_eq, b_eq = pins()
        lo = max(0.0, floor - 1e-9)
        bounds = [(0.0, 1.0)] * n
        for k in free:
            if k != j:
                bounds[k] = (lo, 1.0)
        return linprog(
            c, A_ub=M, b_ub=np.zeros(n), A_eq=a_eq, b_eq=b_eq, bounds=bounds, method="highs"
        )

    free = list(range(n))
    while free:
        res = max_min_level(free)
        if not res.success:
            return None
        t = res.x[n]
        # pin every free coord that cannot strictly exceed t* with the others at the floor
        newly = []
        for j in list(free):
            r = max_coord(j, free, t)
            if (not r.success) or (-r.fun) <= t + tol:
                fixed[j] = t
                newly.append(j)
        if not newly:  # numerical safety: pin the achieved minimum coordinate
            j = min(free, key=lambda k: res.x[k])
            fixed[j] = res.x[j]
            newly = [j]
        free = [k for k in free if fixed[k] is None]
    p = np.clip(np.array(fixed, dtype=float), 0.0, None)
    return p / p.sum()


def max_entropy_equilibrium(M, tol=1e-8):
    """Alternative canonical equilibrium: the max-entropy point of O.

    Maximises -Σ pᵢ log pᵢ subject to 1ᵀp = 1 and A p <= 0. The objective is
    strictly concave and O is convex, so the maximiser is unique. It is the
    most-uniform equilibrium, reduces to the usual point when O is a single point,
    and (unlike a symmetric split) accounts for the game's structure -- e.g. the
    cop game gives (0.216, 0.216, 0.177, 0.392), not (0.2, 0.2, 0.2, 0.4).

    Returns None when O has no point at all (never happens here: every connected
    paradoxical game has at least a boundary equilibrium). For a fully-mixed game
    the optimum is interior (A p = 0), so we optimise within ker(A); otherwise we
    optimise over the full polytope O.
    """
    M = np.asarray(M, dtype=float)
    n = len(M)
    fm, witness = has_fully_mixed(M, tol)

    def negent(p):
        pp = np.clip(p, 1e-300, 1.0)
        return float((pp * np.log(pp)).sum())

    if fm:
        # interior optimum lives in ker(A); parametrise p = N c (N = kernel basis)
        _, S, Vt = np.linalg.svd(M)
        N = Vt[tol * max(S[0], 1.0) > S].T  # n x k
        c0 = np.linalg.lstsq(N, witness, rcond=None)[0]

        def neg_c(c):
            return negent(N @ c)

        def grad_c(c):
            p = np.clip(N @ c, 1e-300, 1.0)
            return N.T @ (np.log(p) + 1.0)

        cons = [
            {"type": "eq", "fun": lambda c: (N @ c).sum() - 1, "jac": lambda c: N.sum(0)},
            {"type": "ineq", "fun": lambda c: N @ c, "jac": lambda c: N},
        ]
        res = minimize(
            neg_c,
            c0,
            jac=grad_c,
            method="SLSQP",
            constraints=cons,
            options={"ftol": 1e-14, "maxiter": 1000},
        )
        p = np.clip(N @ res.x, 0.0, None)
        return p / p.sum()

    # boundary case: optimise over O = {p>=0, 1ᵀp=1, A p <= 0} directly
    def grad(p):
        pp = np.clip(p, 1e-300, 1.0)
        return np.log(pp) + 1.0

    cons = [
        {"type": "eq", "fun": lambda p: p.sum() - 1, "jac": lambda p: np.ones(n)},
        {"type": "ineq", "fun": lambda p: -M @ p, "jac": lambda p: -M},
    ]
    res = minimize(
        negent,
        np.ones(n) / n,
        jac=grad,
        method="SLSQP",
        bounds=[(0.0, 1.0)] * n,
        constraints=cons,
        options={"ftol": 1e-14, "maxiter": 1000},
    )
    if not res.success:
        return None
    p = np.clip(res.x, 0.0, None)
    return p / p.sum()


def equilibrium_vertices(M, tol=1e-7):
    """Extreme equilibria (vertices) of the symmetric-Nash polytope

        O = { p in the simplex : M p <= 0 }.

    O is the set of all symmetric Nash equilibria (p is a best response to itself
    iff (M p)_i <= pᵀM p = 0 for every move i). It is a convex polytope, so the
    full solution set is the convex hull of these vertices: 1 vertex = a unique
    equilibrium, 2 = a segment (the cop game), more = a polygon/polytope.

    When a fully-mixed equilibrium exists, positivity proves every equilibrium
    lies in ker(M), so we use the fast kernel-coordinate enumeration. Otherwise
    boundary equilibria need only satisfy M p <= 0. We then enumerate vertices of
    the full polytope by selecting n-1 active inequalities from p >= 0 and M p <= 0,
    adjoining sum(p)=1, solving, and retaining feasible full-rank intersections.
    The general path is O(C(2n,n-1)); the fairness censuses stay on the much smaller
    kernel path.
    """
    M = np.asarray(M, dtype=float)
    n = len(M)
    fm, _ = has_fully_mixed(M, tol)
    if not fm:
        # Inequalities C p <= 0: the first n rows encode p >= 0 as -I p <= 0,
        # followed by the Nash inequalities M p <= 0. A vertex has n-1 linearly
        # independent active inequalities in addition to sum(p)=1.
        C = np.vstack([-np.eye(n), M])
        verts, seen = [], set()
        for combo in combinations(range(2 * n), n - 1):
            sysA = np.vstack([np.ones(n), C[list(combo)]])
            if np.linalg.matrix_rank(sysA, tol=1e-10) < n:
                continue
            rhs = np.zeros(n)
            rhs[0] = 1.0
            p = np.linalg.solve(sysA, rhs)
            if (C @ p <= tol).all():
                p = np.clip(p, 0.0, None)
                p /= p.sum()
                key = tuple(np.round(p, 9))
                if key not in seen:
                    seen.add(key)
                    verts.append(p)
        return verts

    _, s, vt = np.linalg.svd(M)
    thresh = tol * max(s[0], 1.0)
    N = vt[s < thresh].T  # n x d kernel basis
    d = N.shape[1]
    if d == 0:
        return []  # defensive: fm implies a nontrivial kernel
    ones_n = np.ones(n) @ N  # 1ᵀN, length d
    verts, seen = [], set()
    for combo in combinations(range(n), d - 1):
        sysA = np.vstack([N[list(combo)], ones_n]) if combo else ones_n[None]
        rhs = np.zeros(d)
        rhs[-1] = 1.0
        if abs(np.linalg.det(sysA)) < 1e-12:
            continue  # chosen facets + normalization not independent
        p = N @ np.linalg.solve(sysA, rhs)
        if (p >= -tol).all():
            p = np.clip(p, 0.0, None)
            p = p / p.sum()
            key = tuple(np.round(p, 6))
            if key not in seen:
                seen.add(key)
                verts.append(p)
    return verts


def num_equilibria(M, tol=1e-7):
    """Number of extreme equilibria (vertices of O). 1 = unique Nash equilibrium."""
    return len(equilibrium_vertices(M, tol))


def equilibrium_dim(M, tol=1e-7):
    """Dimension of the equilibrium polytope O (0 = a single point / unique).

    Computed from the affine hull of the vertices, so it reflects the WHOLE Nash
    set -- including any boundary equilibria -- not just the fully-mixed kernel
    family (whose dimension is kernel_dim - 1, generally a lower bound on this).
    """
    V = equilibrium_vertices(M, tol)
    if len(V) <= 1:
        return 0
    P = np.array(V)
    return int(np.linalg.matrix_rank(P - P[0], tol=1e-6))


def equilibrium_info(M, tol=1e-8):
    """Audit summary of the equilibrium set O for an antisymmetric game.

    Returns a dict:
      nullity     -- dim ker(A) = n - rank(A)
      fully_mixed -- whether a strictly-interior equilibrium exists
      family_dim  -- dimension of the fully-mixed equilibrium family
                     (= nullity - 1 when fully_mixed, else None)
      unique      -- whether the full equilibrium polytope O is a single point
      xs          -- canonical (leximin max-min) equilibrium, or None
    """
    nd = kernel_dim(M, tol)
    fm, _ = has_fully_mixed(M, tol)
    fam = nd - 1 if fm else None
    return {
        "nullity": nd,
        "fully_mixed": fm,
        "family_dim": fam,
        "unique": num_equilibria(M, tol) == 1,
        "xs": maxmin_equilibrium(M, tol),
    }
