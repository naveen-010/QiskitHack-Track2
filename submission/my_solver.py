"""Silent Group solver: iterative (semiclassical) phase estimation + relation-lattice recovery.

Pipeline
  1. One-qubit iterative QPE (Griffiths-Niu / Kitaev) measures, in ONE shot and for ONE random
     character chi, b bits of phi_e = chi(e) for e in (A, B, C, g).
  2. t shots give t points of the dual lattice; LLL finds short integer vectors x with
     x . phi ~ 0 (mod 1) for every shot. They generate a lattice L' that should equal the relation
     lattice   L = {x in Z^4 : A^x0 B^x1 C^x2 g^x3 = e}.
  3. Two-sided check with cheap Hadamard tests:
       L' <= L  every LLL basis vector of L' is a true relation         (verify_generators)
       L <= L'  no element of prime order p in Z^4/L' is trivial in G   (complete)
     A noise-flipped QPE stage can make L' a proper SUBLATTICE of L (small groups); the second
     check finds the missing relations and adds them, so the answer is never built on it.
  4. Hermite normal form of L -> ord(g), discrete logs of A, B, C (or forgery), secret = abc mod n.

Guards: wall-clock (TIME_GUARD s), query budget (never asks the oracle for more than is left),
bounded LLL work, a precision ladder b = 36 -> 44 -> 52. "none" is only returned when the
lattice has passed BOTH checks; on any exhaustion a best-effort integer answer is returned.
"""

import math
import random
import time
from fractions import Fraction

from qiskit import QuantumCircuit

from silentgroup import oracle_gate

ELS = ["A", "B", "C", "g"]          # lattice coordinates; g last so HNF gives ord(g)
MAX_LOG2_GROUP = 60                 # |H| <= |G| <= 2^60, so det(L) <= 2^60
LADDER = (36, 44, 52)               # bits of phase per element; escalate when b is not enough
MIN_SAMPLES = 4                     # 4 unknown coordinates need >= 4 dual points
SAMPLE_CAP = {36: 7, 44: 8, 52: 8}  # samples per precision level before escalating (measured:
                                    # a 2-D part with det ~2^57 needs b=44; waiting longer at 36 wastes budget)
WINDOW = 10                         # consensus only looks at the 10 most recent samples
MAX_REJECTS = 3                     # rejected lattices per level before escalating
NONCYCLIC_MIN = 6                   # non-cyclic quotient: need >= 6 samples and a stable lattice
TIME_GUARD = 400.0                  # seconds; the judge's limit is 600
MAX_LLL_CALLS = 400                 # bound on classical work (each call < 0.1 s)
VERIFY_RESERVE = 8000               # queries kept back so a found lattice can still be verified
MAX_COMPLETE_STEPS = 64             # each step divides det(L') by a prime >= 2, det <= 2^60
LINE_CAP = 40                       # p-torsion lines tested exhaustively up to this many


class _Exhausted(Exception):
    """Raised internally when time or budget would be exceeded; solve() then answers best-effort."""


class _Ctx:
    """Per-solve bookkeeping: oracle, deadline, LLL-call counter and the cache of tested vectors."""

    def __init__(self, oracle):
        self.oracle = oracle
        self.deadline = time.monotonic() + TIME_GUARD
        self.lll_calls = 0
        self.skew_tried = 0                              # skew attempts per solve (each ~15k queries)
        self.relations = set()                           # vectors already verified as relations
        self.non_relations = []                          # vectors already shown NOT to be relations

    def left(self):
        left = self.oracle.budget_left()
        return math.inf if left is None else left   # None means an unlimited budget

    def check_time(self):
        if time.monotonic() > self.deadline:
            raise _Exhausted("time")

    def run(self, qc, shots, per_shot):
        """oracle.run, but never asks for more queries than remain (BudgetExceeded is an error)."""
        self.check_time()
        if per_shot * shots > self.left():
            raise _Exhausted("budget")
        return self.oracle.run(qc, shots=shots)


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
                with qc.if_test((qc.clbits[base + (b - l)], 1)):   # (only the 12 nearest: error < 2^-12)
                    qc.p(-2 * math.pi / 2 ** (l - m + 1), 0)
            qc.h(0)
            qc.measure(0, base + (b - m))
            qc.reset(0)
    return qc


def build_single_circuit(b, element, trunc=12):
    """Same iterative QPE, but for one element alone (used to find a very long relation e^n = 1)."""
    qc = QuantumCircuit(1, b)
    for m in range(b, 0, -1):
        qc.h(0)
        qc.append(oracle_gate(element, 1 << (m - 1)), [0])
        for l in range(m + 1, min(b, m + trunc) + 1):
            with qc.if_test((qc.clbits[b - l], 1)):
                qc.p(-2 * math.pi / 2 ** (l - m + 1), 0)
        qc.h(0)
        qc.measure(0, b - m)
        qc.reset(0)
    return qc


def decode(key, b, nel=4):
    bits = key[::-1]                                     # bits[i] = classical bit i
    out = []
    for ei in range(nel):
        base = ei * b
        out.append(sum(int(bits[base + (b - m)]) << (b - m) for m in range(1, b + 1)))
    return out                                           # ~ 2^b * phi_e


# ----------------------------------------------------------------------------- integer lattices
def lll(basis, delta=Fraction(3, 4)):
    """Exact-rational LLL (Lenstra-Lenstra-Lovasz) on the rows of `basis`."""
    b = [list(map(int, r)) for r in basis]
    n = len(b)
    dot = lambda u, v: sum(x * y for x, y in zip(u, v))
    bs, mu, B = [], [[Fraction(0)] * n for _ in range(n)], []
    for i in range(n):                                   # Gram-Schmidt
        v = [Fraction(x) for x in b[i]]
        for j in range(i):
            mu[i][j] = dot(b[i], bs[j]) / B[j]
            v = [x - mu[i][j] * y for x, y in zip(v, bs[j])]
        bs.append(v)
        B.append(dot(v, v))
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):                   # size-reduce b_k
            q = round(mu[k][j])
            if q:
                b[k] = [x - q * y for x, y in zip(b[k], b[j])]
                for i in range(j):
                    mu[k][i] -= q * mu[j][i]
                mu[k][j] -= q
        if B[k] >= (delta - mu[k][k - 1] ** 2) * B[k - 1]:   # Lovasz condition holds
            k += 1
        else:                                            # swap b_k, b_{k-1}; update GS data
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
    """Row-style upper-triangular HNF of an integer matrix of full column rank (None if not).

    Accepts any number (>= 4) of generating rows; returns the 4 basis rows."""
    M = [list(r) for r in rows]
    n = len(M[0])
    out = []
    for c in range(n):
        while True:                                      # Euclid on column c across rows
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
        out.append([-x for x in p] if p[c] < 0 else p)  # positive pivots
    for i in range(n):                                   # canonical form: 0 <= H[i][j] < H[j][j]
        for j in range(i + 1, n):                        # (so equal lattices give equal matrices)
            q = out[i][j] // out[j][j]
            if q:
                out[i] = [x - q * y for x, y in zip(out[i], out[j])]
    return out


def det_of(H):
    return math.prod(H[i][i] for i in range(4))         # HNF is triangular: det = product of pivots


def short_relations(ctx, samples, b, ndim=4, slack=12):
    """LLL on [I | X ; 0 | 2^b I]; returns the rows (head part) whose tails are small, in LLL order."""
    ctx.check_time()
    if ctx.lll_calls >= MAX_LLL_CALLS:
        raise _Exhausted("lll")
    ctx.lll_calls += 1
    t = len(samples)
    rows = []
    for c in range(ndim):
        rows.append([1 if j == c else 0 for j in range(ndim)] + [s[c] for s in samples])
    for i in range(t):
        rows.append([0] * ndim + [(1 << b) if j == i else 0 for j in range(t)])
    red = lll(rows)
    tol = 1 << (b - slack)                               # true relations: |x . err| << 2^(b-12)
    good = []
    for r in red:
        tail = [((v + (1 << (b - 1))) % (1 << b)) - (1 << (b - 1)) for v in r[ndim:]]   # centred mod 2^b
        if max(abs(v) for v in tail) <= tol and any(r[:ndim]):
            good.append(r[:ndim])
    return good


def independent_prefix(rows, k=4):
    """Greedily keep rows (in LLL order) that are linearly independent of those kept, up to k."""
    kept, echelon = [], []                               # echelon: reduced copies over Q
    for r in rows:
        v = [Fraction(x) for x in r]
        for e in echelon:
            piv = next(i for i, x in enumerate(e) if x)
            if v[piv]:
                f = v[piv] / e[piv]
                v = [a - f * b_ for a, b_ in zip(v, e)]
        if any(v):
            kept.append(r)
            echelon.append(v)
            if len(kept) == k:
                break
    return kept


def relation_lattice(ctx, samples, b):
    """Short relations -> (H, short basis) of the lattice they generate, or None."""
    good = independent_prefix(short_relations(ctx, samples, b))
    if len(good) < 4:
        return None
    # The FIRST 4 independent good rows only: later "good" rows are often false (head ~ tail ~ 2^21
    # also passes the tail test; measured: hnf(all good) -> many rejections), and in skewed lattices
    # several short true relations can be linearly dependent (measured on a=2,b=3,c=5).
    H = hnf(good)
    if H is None or det_of(H) > (1 << MAX_LOG2_GROUP):   # rank < 4 or index too big: not L
        return None
    return H, lll(H)                                     # short basis = cheap Hadamard tests


def candidates(ctx, samples, b, deep):
    """Lattices from all samples, then leave-one-out, then (if deep) leave-two-out.

    A corrupted shot usually makes the all-samples lattice fail (big groups) or shrink to a
    sublattice (small groups, repaired by `complete`). Duplicates are removed."""
    win = samples[-WINDOW:]                              # stall guard: at most 10 samples
    t = len(win)
    subsets = [win]
    subsets += [win[:i] + win[i + 1:] for i in range(t)] if t - 1 >= MIN_SAMPLES else []
    if deep and t - 2 >= MIN_SAMPLES:
        subsets += [[s for k, s in enumerate(win) if k not in (i, j)] for i in range(t) for j in range(i + 1, t)]
    seen, out = set(), []
    for sub in subsets:
        res = relation_lattice(ctx, sub, b)
        if res is None:
            continue
        key = tuple(map(tuple, res[0]))
        if key not in seen:
            seen.add(key)
            out.append(res)
        if len(out) == 1 and sub is win:
            return out                                   # all samples agree: nothing to vote on
    return out


def extract(H):
    """(order of g, [log_g A, log_g B, log_g C] with None = not in <g>)."""
    n = H[3][3]                                          # last pivot = smallest n with g^n = e
    logs = []
    for target in range(3):
        # want c . H[:3] = e_target in the first 3 coordinates (triangular back-substitution)
        c = [0, 0, 0]
        ok = True
        for i in range(3):
            want = (1 if i == target else 0) - sum(c[j] * H[j][i] for j in range(i))
            if want % H[i][i]:
                ok = False                               # no relation X^1 g^y = e: X not in <g>
                break
            c[i] = want // H[i][i]
        if not ok:
            logs.append(None)
        else:
            a = sum(c[i] * H[i][3] for i in range(3))   # X * g^a = e  =>  log_g X = -a
            logs.append((-a) % n)
    return n, logs


# ----------------------------------------------------------------------------- number theory
def is_prime(n):
    """Deterministic Miller-Rabin for n < 3.3e24 (bases = first 12 primes)."""
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def pollard_rho(n, ctx):
    """A non-trivial factor of composite n (Brent's variant); deadline-capped."""
    if n % 2 == 0:
        return 2
    rng = random.Random(n)                               # deterministic per n (reproducible)
    while True:
        ctx.check_time()
        y, c, m = rng.randrange(1, n), rng.randrange(1, n), 128
        g, r, q = 1, 1, 1
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += m
            r *= 2
        if g == n:                                       # overshoot: step one at a time
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
        if g != n:
            return g


def prime_factors(n, ctx):
    """Set of primes dividing n (n <= 2^60 here: trial division to 1000, then Pollard rho)."""
    out = set()
    for p in range(2, 1000):
        while n % p == 0:
            out.add(p)
            n //= p
    stack = [n] if n > 1 else []
    while stack:
        m = stack.pop()
        if is_prime(m):
            out.add(m)
        else:
            f = pollard_rho(m, ctx)
            stack += [f, m // f]
    return sorted(out)


def left_kernel_mod_p(M, p):
    """Basis of {c in F_p^4 : c . M = 0 (mod p)} for a 4x4 integer matrix M."""
    A = [[M[i][j] % p for i in range(4)] for j in range(4)]   # A = M^T, solve A c = 0
    pivots, row = [], 0
    for col in range(4):
        r = next((i for i in range(row, 4) if A[i][col]), None)
        if r is None:
            continue
        A[row], A[r] = A[r], A[row]
        inv = pow(A[row][col], -1, p)
        A[row] = [x * inv % p for x in A[row]]
        for i in range(4):
            if i != row and A[i][col]:
                f = A[i][col]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[row])]
        pivots.append(col)
        row += 1
    free = [c for c in range(4) if c not in pivots]
    basis = []
    for f in free:                                       # one kernel vector per free column
        c = [0] * 4
        c[f] = 1
        for i, pc in enumerate(pivots):
            c[pc] = (-A[i][f]) % p
        basis.append(c)
    return basis


def p_rank(basis, p):
    return len(left_kernel_mod_p(basis, p))              # = number of invariant factors divisible by p


def reduce_mod_lattice(x, basis):
    """Babai rounding: x - round(x B^-1) B. Same coset of L', much smaller entries (cheaper test)."""
    for _ in range(2):
        # solve y . basis = x over Q by Cramer's rule-free elimination (4x4, exact)
        M = [[Fraction(basis[i][j]) for i in range(4)] + [Fraction(x[j])] for j in range(4)]
        for c in range(4):
            r = next(i for i in range(c, 4) if M[i][c] != 0)
            M[c], M[r] = M[r], M[c]
            for i in range(4):
                if i != c and M[i][c] != 0:
                    f = M[i][c] / M[c][c]
                    M[i] = [a - f * b_ for a, b_ in zip(M[i], M[c])]
        y = [M[i][4] / M[i][i] for i in range(4)]
        x = [x[j] - sum(round(y[i]) * basis[i][j] for i in range(4)) for j in range(4)]
    return x


def torsion_candidates(basis, p):
    """Vectors x with p*x in L', x not in L': one per line of (Z^4/L')[p] (capped at LINE_CAP)."""
    ker = left_kernel_mod_p(basis, p)
    r = len(ker)
    if r == 0:
        return []
    lines = (p ** r - 1) // (p - 1)
    rng = random.Random(p * 1000003 + r)                 # deterministic choice when sampling
    combos = []
    if lines <= LINE_CAP:
        # every line has a unique representative whose last non-zero coefficient is 1:
        # choose that position k, then the k earlier coefficients freely (p^k choices)
        for k in range(r):
            for idx in range(p ** k):
                combos.append([(idx // p ** j) % p for j in range(k)] + [1] + [0] * (r - k - 1))
    else:
        combos = [[rng.randrange(p) for _ in range(r)] for _ in range(LINE_CAP)]
        combos = [co for co in combos if any(co)]
    out = []
    for co in combos:
        c = [sum(co[k] * ker[k][i] for k in range(r)) % p for i in range(4)]
        v = [sum(c[i] * basis[i][j] for i in range(4)) for j in range(4)]
        out.append(reduce_mod_lattice([x // p for x in v], basis))   # exact: c.basis = 0 mod p
    return out


# ----------------------------------------------------------------------------- Hadamard tests
def hadamard_circuit(x):
    """H - controlled (A^x0 B^x1 C^x2 g^x3) - H - measure. P(1) = 0 if a relation, 1/2 otherwise."""
    qc = QuantumCircuit(1, 1)
    qc.h(0)
    for e, k in zip(ELS, x):
        if k:
            qc.append(oracle_gate(e, int(k)), [0])
    qc.h(0)
    qc.measure(0, 0)
    return qc


NOISE = 0.002                       # competition noise (turns per sqrt(query)), from the rules
TEST_ERR = 1e-4                     # max false-accept AND false-reject probability per relation test


def test_plan(per_shot):
    """Smallest (n, k) for the rule 'relation iff fewer than k ones in n shots' meeting TEST_ERR.

    A true relation still gives P(1) = E[sin^2(pi*N(0, s))] = (1 - exp(-2 pi^2 s^2))/2 with
    s = NOISE*sqrt(cost) (all gates' noise adds up); a non-relation gives P(1) = 1/2 exactly."""
    p_true = (1 - math.exp(-2 * math.pi ** 2 * NOISE ** 2 * per_shot)) / 2
    for n in (16, 20, 24, 28, 32, 40, 48, 64, 96, 128):
        for k in range(1, n + 1):
            false_reject = 1 - sum(math.comb(n, j) * p_true ** j * (1 - p_true) ** (n - j) for j in range(k))
            false_accept = sum(math.comb(n, j) for j in range(k)) / 2 ** n
            if false_accept > TEST_ERR:
                break                                    # larger k only accepts more: try a larger n
            if false_reject <= TEST_ERR:
                return n, k
    return 128, 40                                       # unreachable for costs below ~10^4


def relation_test(ctx, x):
    """Sequential version of 'run n shots, x is a relation iff fewer than k ones'.

    Stops as soon as the fixed-size outcome is decided (k ones seen, or k impossible), so the error
    rates are exactly those of the fixed test; a non-relation is usually rejected after ~2k shots."""
    qc = hadamard_circuit(x)
    per_shot = ctx.oracle.cost(qc, 1)
    n_max, k_rej = test_plan(per_shot)
    ones = shots = 0
    batch = k_rej                                        # cannot reject with fewer than k shots
    while True:
        if ones >= k_rej:
            return False                                 # already rejected
        if ones + (n_max - shots) < k_rej:
            return True                                  # cannot reach the threshold any more
        step = min(batch, n_max - shots)
        ones += ctx.run(qc, step, per_shot).get("1", 0)
        shots += step
        batch = 4


def in_lattice(x, H):
    """Is x an integer combination of the HNF rows? (triangular back-substitution)"""
    x = list(x)
    for i in range(4):
        if x[i] % H[i][i]:
            return False
        q = x[i] // H[i][i]
        x = [a - q * h for a, h in zip(x, H[i])]
    return True


def verify_generators(ctx, short, H=None):
    """L' <= L: every basis vector must be a relation (relation_test: both error rates <= 1e-4;
    typically 'reject if >= 3 ones in 24'). Results are cached per vector: a known
    relation is not re-tested, and a lattice containing a known non-relation is rejected for free."""
    if H is not None and any(in_lattice(x, H) for x in ctx.non_relations):
        return False
    for x in sorted(short, key=lambda v: tuple(v) in ctx.relations):   # untested vectors first
        key = tuple(x)
        if key in ctx.relations:
            continue
        if not relation_test(ctx, x):
            ctx.non_relations.append(key)
            return False
        ctx.relations.add(key)
    return True


def downward_repair(ctx, short):
    """A rejected lattice L' may still be close to L: if exactly one basis vector y is not a relation
    but k*y is (small k), then L'' = <other basis vectors, k*y> is a sublattice of L, which `complete`
    can enlarge to L. Measured case: 15 shots all had chi(h) = 0 for one order-2 element h."""
    bad = [x for x in short if tuple(x) in ctx.non_relations]
    if len(bad) != 1:
        return None
    others = [x for x in short if x is not bad[0]]
    if not verify_generators(ctx, others):
        return None
    for k in (2, 3, 4, 5, 6, 7, 8):
        y = [k * v for v in bad[0]]
        if verify_generators(ctx, [y]):
            return hnf(others + [y])
    return None


def is_missing_relation(ctx, x):
    """Is the candidate x (x not in L') a relation? Same noise-aware test as verify_generators."""
    return relation_test(ctx, x)


def complete(ctx, H):
    """L <= L': add any relation x with p*x in L' (p prime | det) until none is left.

    If L' is a proper sublattice of L then L/L' is a non-trivial finite group, so it contains an
    element x of prime order p, and p divides [L:L'] which divides det(L'). Testing one vector per
    line of the p-torsion of Z^4/L' therefore detects every missing relation."""
    basis = lll(H)
    cleared = set()                                      # primes p with no missing relation of order p
    last = None
    for _ in range(MAX_COMPLETE_STEPS):
        found = None
        primes = [q for q in prime_factors(det_of(H), ctx) if q not in cleared]
        primes.sort(key=lambda q: q != last)             # after a repair try the same prime first
        for p in primes:
            for x in torsion_candidates(basis, p):
                if is_missing_relation(ctx, x):
                    found = x
                    break
            if found:
                last = p
                break
            # [L:L'] is coprime to p now, hence also for every larger L'' between L' and L
            cleared.add(p)
        if found is None:
            return H                                     # both checks passed: L' = L
        H = hnf(basis + [found])                         # enlarge L' by the missing relation
        basis = lll(H)
    return None                                          # did not converge: do not trust


def is_noncyclic(ctx, H):
    """Z^4/L' is non-cyclic iff some prime p | det has p-rank >= 2."""
    basis = lll(H)
    return any(p_rank(basis, p) >= 2 for p in prime_factors(det_of(H), ctx))


# ----------------------------------------------------------------------------- skewed lattices
SKEW_B = 80                         # single-element precision: n <= 2^60 needs b > 60 + 12 (tolerance)
SKEW_SAMPLES = (4, 8)               # min / max single-element samples


def order_of(ctx, ei):
    """ord(ELS[ei]) (or a multiple, fixed later by `complete`) from single-element 80-bit QPE samples.

    1-D version of the lattice step: the shortest h with h * phi_i ~ 0 (mod 1) for every shot."""
    qc = build_single_circuit(SKEW_B, ELS[ei])
    per_shot = ctx.oracle.cost(qc, 1)
    ys, tested = [], set()
    unit = [1 if j == ei else 0 for j in range(4)]
    while len(ys) < SKEW_SAMPLES[1]:
        if ctx.left() < per_shot + VERIFY_RESERVE:
            raise _Exhausted("budget")
        ys.append(decode(next(iter(ctx.run(qc, 1, per_shot))), SKEW_B, nel=1))
        if len(ys) < SKEW_SAMPLES[0]:                    # fewer than 4 shots: spurious short vectors
            continue
        # ~10 % of 80-bit shots have a harmful stage flip, so also try every leave-one-out subset
        t = len(ys)
        subsets = [ys] + ([ys[:i] + ys[i + 1:] for i in range(t)] if t > SKEW_SAMPLES[0] else [])
        if t >= SKEW_SAMPLES[0] + 2:                     # and two bad shots: leave-two-out
            subsets += [[y for k, y in enumerate(ys) if k not in (i, j)] for i in range(t) for j in range(i + 1, t)]
        heads = set()
        for sub in subsets:
            # true tail = n * err <= 2^63, so 2^(80-16) still accepts it; a false h passes one shot's
            # test with P ~ 2^-16 (less if the order is a power of 2: phases then have only log2 n bits)
            good = short_relations(ctx, sub, SKEW_B, ndim=1, slack=16)
            heads |= {abs(v[0]) for v in good if 0 < abs(v[0]) <= 1 << MAX_LOG2_GROUP}
        for h in sorted(heads - tested)[:3]:             # smallest untested candidates, each tested
            tested.add(h)
            if verify_generators(ctx, [[h * u for u in unit]]):
                return h                                 # e^h = 1 (a multiple of the order is fixed by complete)
    return None


def skew_lattice(ctx, M):
    """3 short verified relations M + ord(e) * unit_e: a full-rank sublattice of L of finite index.

    e is the first of g, A, B, C whose unit vector is independent of M (g if possible: e.g. tiny logs;
    A when the long direction is a forged A of large order modulo <B, C, g>)."""
    for ei in (3, 0, 1, 2):
        unit = [1 if j == ei else 0 for j in range(4)]
        if len(independent_prefix([list(r) for r in M] + [unit])) == 4:
            break
    else:
        return None
    n = order_of(ctx, ei)
    if n is None:
        return None
    H = hnf([list(r) for r in M] + [[n * u for u in unit]])
    # det may exceed 2^60 here: M need not be saturated (measured index 262145); `complete` closes
    # the finite index [L : L'], whose primes divide det. 2^90 keeps Pollard rho fast.
    if H is None or det_of(H) > (1 << 90):
        return None
    return H


def small_span(ctx, samples, b):
    """Do the 3 shortest relations span a tiny volume, so that the 4th basis vector of a lattice with
    det <= 2^60 could be longer than b bits resolve (~2^(b-14))? Gram det = volume^2."""
    M = independent_prefix(short_relations(ctx, samples[-WINDOW:], b), 3)
    if len(M) < 3:
        return False
    gram = [[sum(x * y for x, y in zip(u, v)) for v in M] for u in M]
    g = (gram[0][0] * (gram[1][1] * gram[2][2] - gram[1][2] * gram[2][1])
         - gram[0][1] * (gram[1][0] * gram[2][2] - gram[1][2] * gram[2][0])
         + gram[0][2] * (gram[1][0] * gram[2][1] - gram[1][1] * gram[2][0]))
    return g <= 1 << (2 * (MAX_LOG2_GROUP - (b - 14)))


def try_skew(ctx, samples, b, rank3):
    """Skewed lattice: 3 short relations, the 4th basis vector far too long for any b here (e.g. tiny
    logs a, b, c with ord(g) ~ 2^55, or a forged A of large order modulo <B, C, g>). After 2 hints
    (rank-3 data, or a rejected lattice whose 3 shortest relations span a tiny volume) measure the order
    of one element directly; the SAME two-sided tests then decide.
    Returns None (not tried), "fail", or an answer."""
    if rank3 < 2 or ctx.skew_tried >= 2 or len(samples) < MIN_SAMPLES + 1:
        return None
    M = independent_prefix(short_relations(ctx, samples[-WINDOW:], b), 3)
    if len(M) < 3 or not verify_generators(ctx, M):
        return "fail"
    H = skew_lattice(ctx, M)
    if H is None or not verify_generators(ctx, lll(H)):
        return "fail"
    Hc = complete(ctx, H)
    return "fail" if Hc is None else answer(Hc, trusted=True)


# ----------------------------------------------------------------------------- main
def answer(H, trusted):
    n, logs = extract(H)
    if None in logs:
        # "none" only on positive, verified evidence; otherwise a best-effort integer
        return {"order": n, "secret": "none" if trusted else 0}
    return {"order": n, "secret": math.prod(logs) % n}


def solve(oracle):
    ctx = _Ctx(oracle)
    best = None                                          # (H, trusted) best-effort fallback
    try:
        for b in LADDER:
            qc = build_pe_circuit(b)
            per_shot = oracle.cost(qc, 1)
            samples, history, rejected, n_rej = [], {}, set(), 0
            repaired = set()                             # lattices already given a downward repair
            rank3 = 0
            while len(samples) < SAMPLE_CAP[b] and n_rej < MAX_REJECTS:
                if ctx.left() < per_shot + VERIFY_RESERVE:
                    raise _Exhausted("budget")
                key = next(iter(ctx.run(qc, 1, per_shot)))
                samples.append(decode(key, b))
                if len(samples) < MIN_SAMPLES:
                    continue
                cands = candidates(ctx, samples, b, deep=False)
                if not cands and len(samples) >= MIN_SAMPLES + 2:
                    cands = candidates(ctx, samples, b, deep=True)   # tolerate 2 corrupted shots
                if not cands:
                    # skew hint: exactly 3 independent short relations AND they span a tiny volume
                    # (without the volume test a corrupted shot also gives 3 rows: measured on public_2)
                    if len(independent_prefix(short_relations(ctx, samples[-WINDOW:], b))) == 3 and                             small_span(ctx, samples, b):
                        rank3 += 1
                    res = try_skew(ctx, samples, b, rank3)
                    if res is not None:
                        ctx.skew_tried += 1
                        rank3 = 0                        # a retry needs 2 fresh hints
                        if res != "fail":
                            return res
                    continue
                t = len(samples)
                keys = [tuple(map(tuple, H)) for H, _ in cands]
                history[t] = set(keys)
                # vote: lattices also found at the two previous sample counts go first
                votes = {k: (k in history.get(t - 1, ())) + (k in history.get(t - 2, ())) for k in keys}
                order = sorted(range(len(cands)), key=lambda i: -votes[keys[i]])
                for i in order:
                    (H, short), hkey = cands[i], keys[i]
                    if hkey in rejected:
                        # a rejected lattice that the data keeps producing (3 counts in a row): the
                        # samples are fooled about one low-order element -> try the downward repair once
                        if votes[hkey] >= 2 and hkey not in repaired:
                            repaired.add(hkey)
                            Hd = downward_repair(ctx, short)
                            if Hd is not None:
                                Hc = complete(ctx, Hd)   # L'' <= L verified; complete gives L
                                if Hc is not None:
                                    return answer(Hc, trusted=True)
                        continue
                    if best is None:
                        best = (H, False)
                    # rank-2 / non-cyclic quotients are where 4 samples are least reliable: demand
                    # >= 6 samples and the same lattice at the two previous counts (budget permitting)
                    if is_noncyclic(ctx, H) and t < SAMPLE_CAP[b] and (t < NONCYCLIC_MIN or votes[hkey] < 2):
                        break                            # sample more before spending on tests
                    if n_rej >= MAX_REJECTS:
                        break                            # escalate precision instead
                    if not verify_generators(ctx, short, H):
                        rejected.add(hkey)
                        n_rej += 1
                        if small_span(ctx, samples, b):  # rejected 4th vector + tiny 3-volume: skew hint
                            rank3 += 1
                        res = try_skew(ctx, samples, b, rank3)
                        if res is not None:
                            ctx.skew_tried += 1
                            rank3 = 0                    # a retry needs 2 fresh hints
                            if res != "fail":
                                return res
                        continue
                    best = (H, False)                    # L' <= L holds
                    Hc = complete(ctx, H)
                    if Hc is None:
                        rejected.add(hkey)
                        n_rej += 1
                        continue
                    return answer(Hc, trusted=True)
    except _Exhausted:
        pass
    if best is None:
        return {"order": 1, "secret": 0}                 # nothing recovered: still never "none"
    return answer(best[0], trusted=False)
