"""Hand-built toy games with known properties, shared across test modules.

Each matrix is an n-by-n skew-symmetric {-1,0,+1} game built from a list of
(winner, loser) decisive edges; unlisted pairs tie.
"""

import numpy as np


def beat(M, w, l):
    M[w, l] = 1
    M[l, w] = -1


def game(n, edges):
    M = np.zeros((n, n), dtype=np.int8)
    for w, l in edges:
        beat(M, w, l)
    return M


# Rock(0) > Scissors-as-1 ... we just use the cyclic 0>1>2>0 (isomorphic to RPS).
RPS = game(3, [(0, 1), (1, 2), (2, 0)])

# Rock-Paper-Scissors-Lizard-Spock: each move beats the next two (mod 5).
RPSLS = game(5, [(i, (i + 1) % 5) for i in range(5)] + [(i, (i + 2) % 5) for i in range(5)])


def ring(n):
    """Each move beats the next (mod n), ties the rest -- the circulant 'ring'."""
    return game(n, [(i, (i + 1) % n) for i in range(n)])


# The n=4 "cop" game. Cop=0, K-9=1, Perp=2, Witness=3.
# Witness>Cop, Cop>Perp, Cop>K-9, K-9>Perp, Perp>Witness; K-9 ties Witness.
COP = game(4, [(3, 0), (0, 2), (0, 1), (1, 2), (2, 3)])

# Brick(0)+Boulder(1)+Paper(2)+Scissors(3): Brick & Boulder are tie-twins.
# Paper beats both rocks, loses to Scissors; Scissors beats Paper, loses to rocks.
BRICK = game(4, [(2, 0), (2, 1), (3, 2), (0, 3), (1, 3)])

# Elemental n=5: Water=0, Fire=1, Clay=2, Sand=3, Grass=4. Clay & Grass are twins.
ELEM = game(5, [(0, 1), (0, 3), (1, 2), (1, 4), (2, 0), (3, 1), (4, 0)])

# Paley tournament Q_7: i beats j iff (j-i) mod 7 is a quadratic residue {1,2,4}.
# The unique smallest two-paradox (P2) tournament.
PALEY7 = game(7, [(i, (i + d) % 7) for i in range(7) for d in (1, 2, 4)])

# RPS plus a strictly dominated 4th move (loses to all) -> NOT paradoxical,
# no fully-mixed equilibrium; unique boundary equilibrium drops the 4th move.
DOMINATED = game(4, [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)])

# Paradoxical and connected but without a fully-mixed equilibrium. Its unique
# equilibrium is the boundary point (1/3, 1/3, 1/3, 0), disproving the tempting
# but false extension of the inclusive kernel/parity theorem to all fair games.
BOUNDARY = game(4, [(0, 1), (1, 2), (2, 0), (1, 3), (2, 3), (3, 0)])

# A prime, rigid (|Aut|=1) *regular* n=9 game (profile (3,2,3)) with 11 extreme
# equilibria -- the first counterexample to `prime => n_eq <= n`. Found by the
# regular n=9 enumeration; hard-coded here so the regression test is self-contained.
PRIME_NEQ11 = np.array(
    [
        [0, 1, 1, 1, 0, 0, -1, -1, -1],
        [-1, 0, 0, -1, 1, 0, 1, 1, -1],
        [-1, 0, 0, 1, 1, 1, -1, 0, -1],
        [-1, 1, -1, 0, -1, 1, 0, 1, 0],
        [0, -1, -1, 1, 0, 1, -1, 0, 1],
        [0, 0, -1, -1, -1, 0, 1, 1, 1],
        [1, -1, 1, 0, 1, -1, 0, -1, 0],
        [1, -1, 0, -1, 0, -1, 1, 0, 1],
        [1, 1, 1, 0, -1, -1, 0, -1, 0],
    ],
    dtype=np.int8,
)
