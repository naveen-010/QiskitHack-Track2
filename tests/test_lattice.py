"""Unit tests for the classical lattice code: lll, hnf, extract."""

import math
import random
from fractions import Fraction

import pytest


def det(M):
    """Exact determinant by fraction Gaussian elimination."""
    A = [[Fraction(x) for x in r] for r in M]
    n, d = len(A), Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c]), None)
        if p is None:
            return 0
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return d


def in_lattice(v, basis):
    """Is v an integer combination of the rows of the (square, full-rank) basis?"""
    n = len(basis)
    # solve c . basis = v over Q by Cramer, check integrality
    D = det(basis)
    for i in range(n):
        M = [list(r) for r in basis]
        M[i] = list(v)
        if (det(M) / D).denominator != 1:
            return False
    return True


def same_lattice(B1, B2):
    return all(in_lattice(r, B2) for r in B1) and all(in_lattice(r, B1) for r in B2)


def gram_schmidt(B):
    bs = []
    for r in B:
        v = [Fraction(x) for x in r]
        for u in bs:
            mu = sum(a * b for a, b in zip(r, u)) / sum(a * a for a in u)
            v = [a - mu * b for a, b in zip(v, u)]
        bs.append(v)
    return bs


def random_basis(rng, n, size):
    while True:
        B = [[rng.randrange(-size, size + 1) for _ in range(n)] for _ in range(n)]
        if det(B):
            return B


@pytest.mark.parametrize("seed", range(8))
def test_lll_same_lattice_and_reduced(solver, seed):
    rng = random.Random(seed)
    B = random_basis(rng, 5, 10**6)
    R = solver.lll(B)
    assert abs(det(R)) == abs(det(B))                       # unimodular change of basis
    assert same_lattice(R, B)
    bs = gram_schmidt(R)
    norms = [sum(x * x for x in v) for v in bs]
    for k in range(1, len(R)):                              # Lovasz condition (delta = 3/4)
        mu = sum(a * b for a, b in zip(R[k], bs[k - 1])) / norms[k - 1]
        assert norms[k] >= (Fraction(3, 4) - mu * mu) * norms[k - 1]
    for i in range(len(R)):                                 # size reduction |mu_ij| <= 1/2
        for j in range(i):
            mu = sum(a * b for a, b in zip(R[i], bs[j])) / norms[j]
            assert abs(mu) <= Fraction(1, 2)


@pytest.mark.parametrize("seed", range(8))
def test_hnf_triangular_same_lattice(solver, seed):
    rng = random.Random(100 + seed)
    B = random_basis(rng, 4, 50)
    H = solver.hnf(B)
    assert H is not None
    for i in range(4):
        assert H[i][i] > 0
        for j in range(i):
            assert H[i][j] == 0                             # upper triangular
    assert math.prod(H[i][i] for i in range(4)) == abs(det(B))
    assert same_lattice(H, B)


def test_hnf_singular_returns_none(solver):
    assert solver.hnf([[1, 2, 3, 4], [2, 4, 6, 8], [0, 0, 1, 0], [0, 0, 0, 1]]) is None


def relation_lattice_basis(moduli, elems):
    """Brute-force-free exact basis of {x : sum x_i e_i = 0} for elements of Z_d1 x ... (test oracle).

    Built as the kernel of Z^4 -> G via HNF of [e_i | I] stacked with the moduli relations."""
    m = len(moduli)
    rows = [list(e) + [1 if j == i else 0 for j in range(4)] for i, e in enumerate(elems)]
    rows += [[d if j == k else 0 for j in range(m)] + [0] * 4 for k, d in enumerate(moduli)]
    # integer row reduction on the first m columns; rows whose first m entries vanish span the kernel
    M = [r[:] for r in rows]
    piv_rows = []
    for c in range(m):
        while True:
            nz = [r for r in M if r[c] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda r: abs(r[c]))
            p = nz[0]
            for r in nz[1:]:
                q = r[c] // p[c]
                for k in range(len(r)):
                    r[k] -= q * p[k]
        piv = [r for r in M if r[c] != 0]
        if piv:
            M.remove(piv[0])
            piv_rows.append(piv[0])
    return [r[m:] for r in M if any(r[m:])]


def group_case(moduli, g, A, B, C):
    from silentgroup.group import AbelianGroup, compute_answer
    G = AbelianGroup(tuple(moduli))
    els = {"g": G.reduce(g), "A": G.reduce(A), "B": G.reduce(B), "C": G.reduce(C)}
    return [els["A"], els["B"], els["C"], els["g"]], compute_answer(G, els), G, els


@pytest.mark.parametrize("case", [
    ([97], [5], [10], [15], [96]),                 # cyclic, logs 2, 3, 19*... all in <g>
    ([12, 6], [1, 1], [2, 2], [3, 3], [0, 1]),     # C = (0,1) not in <(1,1)> -> None
    ([8, 4], [1, 2], [3, 2], [0, 0], [7, 2]),      # A in <g>, B = e, C in <g>
    ([2**20, 3**5], [3, 1], [6, 2], [9, 3], [3 * 1000, 1000]),
])
def test_extract_hand_built(solver, case):
    moduli, g, A, B, C = case
    elems, answer, G, els = group_case(moduli, g, A, B, C)
    gens = relation_lattice_basis(moduli, elems)
    H = solver.hnf(gens)                            # hnf accepts >= 4 generating rows
    assert H is not None
    n, logs = solver.extract(H)
    assert n == answer["order"]
    for k, lg in zip("ABC", logs):
        true = G.discrete_log(els["g"], els[k])
        assert (lg is None) == (true is None)
        if true is not None:
            assert lg % n == true % n
