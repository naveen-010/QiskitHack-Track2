"""Silent Group solver: iterative (semiclassical) phase estimation + relation-lattice recovery.

Pipeline
  1. One-qubit iterative QPE (Griffiths-Niu / Kitaev) measures, in ONE shot and for ONE random
     character chi, b bits of phi_e = chi(e) for e in (A, B, C, g).
  2. t shots give t points of the dual lattice; LLL finds the relation lattice
        L = {x in Z^4 : A^x0 B^x1 C^x2 g^x3 = e}.
  3. Hermite normal form of L -> ord(g), discrete logs of A, B, C (or forgery), secret = abc mod n.
  4. Cheap Hadamard tests verify every claimed relation before we answer.
"""

import math
from fractions import Fraction

from qiskit import QuantumCircuit

from silentgroup import oracle_gate

ELS = ["A", "B", "C", "g"]          # lattice coordinates; g last so HNF gives ord(g)
MAX_LOG2_GROUP = 60


# ----------------------------------------------------------------------------- phase estimation
def build_pe_circuit(b, trunc=12):
    """Iterative QPE with a single recycled qubit; bit x_m of each phase goes to its own clbit."""
    qc = QuantumCircuit(1, 4 * b)
    for ei, e in enumerate(ELS):
        base = ei * b
        for m in range(b, 0, -1):                       # LSB first: power 2^(m-1)
            qc.h(0)
            qc.append(oracle_gate(e, 1 << (m - 1)), [0])
            for l in range(m + 1, min(b, m + trunc) + 1):  # feed back already-measured bits
                with qc.if_test((qc.clbits[base + (b - l)], 1)):
                    qc.p(-2 * math.pi / 2 ** (l - m + 1), 0)
            qc.h(0)
            qc.measure(0, base + (b - m))
            qc.reset(0)
    return qc


def decode(key, b):
    bits = key[::-1]                                     # bits[i] = classical bit i
    out = []
    for ei in range(4):
        base = ei * b
        out.append(sum(int(bits[base + (b - m)]) << (b - m) for m in range(1, b + 1)))
    return out                                           # ~ 2^b * phi_e


# ----------------------------------------------------------------------------- integer lattices
def lll(basis, delta=Fraction(3, 4)):
    b = [list(map(int, r)) for r in basis]
    n = len(b)
    dot = lambda u, v: sum(x * y for x, y in zip(u, v))
    bs, mu, B = [], [[Fraction(0)] * n for _ in range(n)], []
    for i in range(n):
        v = [Fraction(x) for x in b[i]]
        for j in range(i):
            mu[i][j] = dot(b[i], bs[j]) / B[j]
            v = [x - mu[i][j] * y for x, y in zip(v, bs[j])]
        bs.append(v)
        B.append(dot(v, v))
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                b[k] = [x - q * y for x, y in zip(b[k], b[j])]
                for i in range(j):
                    mu[k][i] -= q * mu[j][i]
                mu[k][j] -= q
        if B[k] >= (delta - mu[k][k - 1] ** 2) * B[k - 1]:
            k += 1
        else:
            b[k], b[k - 1] = b[k - 1], b[k]
            m_ = mu[k][k - 1]
            Bn = B[k] + m_ * m_ * B[k - 1]
            mu[k][k - 1] = m_ * B[k - 1] / Bn
            B[k] = B[k - 1] * B[k] / Bn
            B[k - 1] = Bn
            for j in range(k - 1):
                mu[k][j], mu[k - 1][j] = mu[k - 1][j], mu[k][j]
            for i in range(k + 1, n):
                t = mu[i][k]
                mu[i][k] = mu[i][k - 1] - m_ * t
                mu[i][k - 1] = t + mu[k][k - 1] * mu[i][k]
            k = max(k - 1, 1)
    return b


def hnf(rows):
    """Row-style upper-triangular HNF of a full-rank square integer matrix (None if singular)."""
    M = [list(r) for r in rows]
    n = len(M[0])
    out = []
    for c in range(n):
        while True:
            nz = [r for r in M if r[c] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda r: abs(r[c]))
            p = nz[0]
            for r in nz[1:]:
                q = r[c] // p[c]
                for k in range(n):
                    r[k] -= q * p[k]
        piv = [r for r in M if r[c] != 0]
        if not piv:
            return None
        p = piv[0]
        M.remove(p)
        out.append([-x for x in p] if p[c] < 0 else p)
    return out


def relation_lattice(samples, b, max_fail=2):
    """Short integer x with x . X_i ~ 0 (mod 2^b) for all samples -> 4 generators of the lattice."""
    t = len(samples)
    N = 4 + t
    rows = []
    for c in range(4):
        rows.append([1 if j == c else 0 for j in range(4)] + [s[c] for s in samples])
    for i in range(t):
        rows.append([0] * 4 + [(1 << b) if j == i else 0 for j in range(t)])
    red = lll(rows)
    tol = 1 << (b - 12)
    good = []
    for r in red:
        tail = [((v + (1 << (b - 1))) % (1 << b)) - (1 << (b - 1)) for v in r[4:]]
        if max(abs(v) for v in tail) <= tol and any(r[:4]):
            good.append(r[:4])
    if len(good) < 4:
        return None
    H = hnf(good[:4]) or hnf(good)
    if H is None:
        return None
    det = 1
    for i in range(4):
        det *= H[i][i]
    if det > (1 << MAX_LOG2_GROUP):
        return None
    return H, good[:4]


def consensus(samples, b):
    """Allow a few corrupted shots: try all samples, then leave-one-out, then leave-two-out."""
    t = len(samples)
    H = relation_lattice(samples, b)
    if H is not None:
        return H, samples
    for i in range(t):
        sub = samples[:i] + samples[i + 1:]
        H = relation_lattice(sub, b)
        if H is not None:
            return H, sub
    for i in range(t):
        for j in range(i + 1, t):
            sub = [s for k, s in enumerate(samples) if k not in (i, j)]
            if len(sub) >= 4:
                H = relation_lattice(sub, b)
                if H is not None:
                    return H, sub
    return None, samples


def extract(H):
    """(order of g, [log_g A, log_g B, log_g C] with None = not in <g>)."""
    n = H[3][3]
    logs = []
    for target in range(3):
        c = [0, 0, 0]
        ok = True
        for i in range(3):
            want = (1 if i == target else 0) - sum(c[j] * H[j][i] for j in range(i))
            if want % H[i][i]:
                ok = False
                break
            c[i] = want // H[i][i]
        if not ok:
            logs.append(None)
        else:
            a = sum(c[i] * H[i][3] for i in range(3))
            logs.append((-a) % n)
    return n, logs


# ----------------------------------------------------------------------------- verification
def relation_fails(oracle, x, shots):
    """Hadamard test of the claimed relation prod e^x_e = identity. Returns #outcome-1 shots."""
    qc = QuantumCircuit(1, 1)
    qc.h(0)
    for e, k in zip(ELS, x):
        if k:
            qc.append(oracle_gate(e, int(k)), [0])
    qc.h(0)
    qc.measure(0, 0)
    return oracle.run(qc, shots=shots).get("1", 0)


def small_primes(n, limit=50):
    return [p for p in range(2, limit + 1) if all(p % q for q in range(2, int(p ** 0.5) + 1)) and n % p == 0]


def verify(oracle, n, short, shots=16):
    """Every short generator of the lattice must be a true relation; then order is not too big for tiny primes."""
    for x in short:
        if relation_fails(oracle, x, shots) >= 3:
            return False
    for p in small_primes(n, 7):
        if relation_fails(oracle, [0, 0, 0, n // p], 14) == 0:
            return False
    return True


# ----------------------------------------------------------------------------- main
def solve(oracle, b=36, min_samples=4):
    qc = build_pe_circuit(b)
    per_shot = oracle.cost(qc, 1)
    samples = []
    while True:
        left = oracle.budget_left()
        if left is not None and left < per_shot + 4000:
            break
        key = next(iter(oracle.run(qc, shots=1)))
        samples.append(decode(key, b))
        if len(samples) < min_samples:
            continue
        res, used = consensus(samples, b)
        if res is None:
            continue
        H, short = res
        n, logs = extract(H)
        if verify(oracle, n, short):
            return {"order": n, "secret": "none" if None in logs else math.prod(logs) % n}
    res, _ = consensus(samples, b)
    if res is None:
        return {"order": 1, "secret": "none"}
    n, logs = extract(res[0])
    return {"order": n, "secret": "none" if None in logs else math.prod(logs) % n}
