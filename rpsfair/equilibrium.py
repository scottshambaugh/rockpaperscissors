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

import numpy as np
from scipy.optimize import linprog, minimize


def equilibrium(M):
    """Optimal mixed strategy via constrained least squares: min ||Mx|| s.t. sum(x)=1.

    NOTE: returns one particular solution; when the equilibrium set is not a single
    point (even n with a fully-mixed family) this is non-canonical. Prefer
    `max_entropy_equilibrium` for any reported odds or metrics.
    """
    n = len(M)
    A = np.vstack([M.astype(float), np.ones(n)])
    b = np.zeros(n + 1)
    b[n] = 1
    x, *_ = np.linalg.lstsq(A, b, rcond=None)
    return x


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


def maxmin_equilibrium(M, tol=1e-7):
    """Canonical equilibrium: the leximin point of O = {p simplex : A p <= 0}.

    Repeatedly maximises the smallest play rate, then the next-smallest, ... .
    This is the "max-min" criterion done canonically: a single max-min LP can
    leave a flat optimal face (re-introducing an arbitrary choice), so we fix
    each coordinate that is forced to the current minimum and recurse on the
    rest. The result is unique and lands on clean rationals for these {-1,0,+1}
    games -- e.g. the cop game gives exactly (0.2, 0.2, 0.2, 0.4).

    For a fully-mixed game the optimum is interior (A p = 0); when no fully-mixed
    equilibrium exists it returns the leximin boundary equilibrium. Returns None
    only if O is empty (never happens for a connected paradoxical game).
    """
    M = np.asarray(M, dtype=float)
    n = len(M)
    fixed = [None] * n  # coords pinned to their leximin value

    def solve(free, minimise_for=None):
        # maximise t (or, if minimise_for=j, maximise p_j) over
        #   p_i = fixed[i] (pinned), p_j >= t (free), p >= 0, sum p = 1, A p <= 0
        # variables: [p_0..p_{n-1}, t]
        c = np.zeros(n + 1)
        c[n if minimise_for is None else minimise_for] = -1.0
        A_ub = [np.hstack([M, np.zeros((n, 1))])]  # A p <= 0
        b_ub = [np.zeros(n)]
        for j in free:  # t - p_j <= 0
            row = np.zeros(n + 1)
            row[j] = -1.0
            row[n] = 1.0
            A_ub.append(row[None])
            b_ub.append([0.0])
        eq = [np.append(np.ones(n), 0.0)]  # sum p = 1
        beq = [1.0]
        for i in range(n):  # pin already-fixed coords
            if fixed[i] is not None:
                row = np.zeros(n + 1)
                row[i] = 1.0
                eq.append(row)
                beq.append(fixed[i])
        bounds = [(0.0, 1.0)] * n + [(None, None)]
        return linprog(
            c,
            A_ub=np.vstack(A_ub),
            b_ub=np.concatenate(b_ub),
            A_eq=np.vstack(eq),
            b_eq=beq,
            bounds=bounds,
            method="highs",
        )

    free = list(range(n))
    while free:
        r = solve(free)
        if not r.success:
            return None
        t = r.x[n]
        # a free coord is "stuck" at t iff it cannot exceed t under the current pins
        newly = []
        for j in list(free):
            rj = solve([k for k in free if k != j], minimise_for=j)
            if (not rj.success) or rj.x[j] <= t + tol:
                fixed[j] = t
                newly.append(j)
        if not newly:  # numerical safety: pin the achieved minimum
            j = min(free, key=lambda k: r.x[k])
            fixed[j] = r.x[j]
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


def equilibrium_info(M, tol=1e-8):
    """Audit summary of the equilibrium set O for an antisymmetric game.

    Returns a dict:
      nullity     -- dim ker(A) = n - rank(A)
      fully_mixed -- whether a strictly-interior equilibrium exists
      family_dim  -- dimension of the fully-mixed equilibrium family
                     (= nullity - 1 when fully_mixed, else None)
      unique      -- whether O's fully-mixed set is a single point
      xs          -- canonical (leximin max-min) equilibrium, or None
    """
    nd = kernel_dim(M, tol)
    fm, _ = has_fully_mixed(M, tol)
    fam = nd - 1 if fm else None
    return {
        "nullity": nd,
        "fully_mixed": fm,
        "family_dim": fam,
        "unique": (fam == 0) if fm else None,
        "xs": maxmin_equilibrium(M),
    }
