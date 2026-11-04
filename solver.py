"""Track 2 -- The Silent Group.

    solve(oracle) -> {"order": ord(g), "secret": abc mod ord(g), or "none"}

How the device actually behaves
-------------------------------
The sealed register never entangles with our qubits.  For one shot the device
draws a single character chi of G uniformly at random and holds it fixed for the
whole circuit, and ``oracle_gate(h, k)`` is then exactly the phase gate

    diag(1, exp(2*pi*i * k * chi(h)))

on its control qubit, with Gaussian jitter of 0.002*sqrt(cost) turns added to
the exponent.  So the problem is phase estimation of the angles theta_h =
chi(h), and one qubit is enough.

Two primitives follow.

1.  Zero test.  H . (oracle gates) . H . measure, over many shots.  Character
    orthogonality gives  E_chi[chi(h^k)] = 1 if h^k = e, else 0.  So every shot
    returns 0 when h^k = e, while an h^k != e gives P(0) = 1/2 per shot
    independently.  An all-zeros column certifies h^k = e.  Chaining gates on
    one control tests any word g^k1 A^k2 B^k3 C^k4 = e.  The test is exact and
    essentially noise-immune (jitter enters only as cos(2*pi*eps)), so it is the
    final arbiter: nothing is returned that it has not confirmed.

2.  Iterative phase estimation.  One qubit, reset between rounds, with
    semiclassical Griffiths-Niu feedback through ``if_test``.  It has to live
    inside a single circuit, because chi is redrawn on every shot -- splitting
    the rounds across runs would estimate bits of different angles.  Feedback
    rotations are not oracle gates, so the budget does not charge them.

    Round j applies h^(2^(bits-1-j)) and measures bit j, least significant
    first, so the costliest and noisiest round fixes the least important bit.
    The flip side is that feedback makes a bit error cascade upward -- one wrong
    bit throws the next correction by a quarter turn -- so a sweep is
    all-or-nothing, and must be validated rather than trusted.

Why that solves the problem
---------------------------
g has order n, so chi(g) is always a multiple of 1/n, with chi(g) = s/n for s
uniform on Z_n.  Estimating theta_g and taking the best rational approximation
with denominator <= 2^60 recovers q = n/gcd(s,n), a divisor of n, and the lcm
over a few sweeps reaches n.

Sweeping at 2*60+24 bits rather than the 2*60+3 Shor needs makes each sweep
self-certifying: a clean sweep's convergent explains theta to ~2^-bits, while a
cascaded sweep is an essentially random real whose best denominator-<=2^60
rational sits ~2^-(bits-24) away.  Demanding agreement to 2^-(bits-8)
therefore admits clean sweeps -- they clear it by ~256x -- while letting a
corrupted one through only with probability ~2^-16.

Note carefully that g^L = e proves only n | L, never L = n -- an lcm polluted by
one bad sweep passes that test while being a multiple of the true order.  So the
order is finished by dividing primes out of L under the zero test, which ends at
the unique minimal L with g^L = e.  That is what makes the answer exact rather
than merely likely.

If A = g^a then chi(A) = a*chi(g) for *every* character.  So estimating theta_g
and theta_A in the same shot (hence under one chi) gives s and t with
a*s = t (mod n); a handful of shots pins a down, and the zero test on A*g^-a
confirms it.  A forged A satisfies no single a, so no candidate ever confirms --
this is strictly stronger than testing A^n = e, which the forged broadcast in
public_2 passes.
"""

import math
import random
from fractions import Fraction

from qiskit import QuantumCircuit

from silentgroup import oracle_gate

MAX_ORDER = 1 << 60  # |G| <= 2^60, so ord(g) <= 2^60
FEEDBACK_WINDOW = 12  # truncated correction costs < 2^-13 turns of phase error
# Precision ladder.  A rung at b bits certifies orders up to 2^((b-24)/2), so
# the cheap rungs finish small groups outright and cost less than one full
# sweep between them; cascade risk grows as b^2, so they also fail far less.
ORDER_LADDER = (64, 88, 112, 144, 144, 144, 144, 144)
ORDER_SLACK = 8  # accept a convergent explaining theta to 2^-(bits-8):
# a clean sweep clears that by ~256x, a cascaded one slips through at only 2^-16
MULT_SEARCH = 48  # backstop: look for n = L*m with m this small
LOG_SLACK = 12  # dlog precision past log2(n), giving the same kind of filter
LOG_RESIDUAL = 0.01  # a clean dlog estimate lands this close to a multiple of 1/n
LOG_SHOTS = 8
QUICK_SHOTS = 8  # stage 1 of the zero test: rejects a non-identity cheaply
FULL_SHOTS = 12  # stage 2 tops it up to 20
ZERO_ONES = 2  # accept "= e" only at <=2 ones in 20: false accept 2e-4, reject 2e-5
RESERVE = 16000  # held back so a full dlog shot plus confirmations always fits
MIN_RESERVE = 3000  # minimising the order outranks the dlogs: a wrong order scores 0

DEBUG = False


def _log(*a):
    if DEBUG:
        print("   [solver]", *a)


# --------------------------------------------------------------------------- #
# integer helpers
# --------------------------------------------------------------------------- #
def _congruences(eqs):
    """Solve a*x = b (mod m) simultaneously.  Returns (x0, period) or None."""
    x0, period = 0, 1
    for a, b, m in eqs:
        A, B = a * period % m, (b - a * x0) % m
        d = math.gcd(A, m)
        if B % d:
            return None
        md = m // d
        t = (B // d) * pow(A // d, -1, md) % md if md > 1 else 0
        x0, period = (x0 + period * t) % (period * md), period * md
    return x0, period


def _is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d, r = d // 2, r + 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _factor(n):
    """Distinct prime factors of n (n <= 2^60, so this is instant)."""
    out, rng = set(), random.Random(0)
    stack = [n]
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if _is_prime(m):
            out.add(m)
            continue
        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
            if m % p == 0:
                out.add(p)
                while m % p == 0:
                    m //= p
                break
        else:  # Pollard rho
            while True:
                c, x, y, d = rng.randrange(1, m), 2, 2, 1
                while d == 1:
                    x = (x * x + c) % m
                    y = (y * y + c) % m
                    y = (y * y + c) % m
                    d = math.gcd(abs(x - y), m)
                if d != m:
                    stack += [d, m // d]
                    break
            continue
        stack.append(m)
    return sorted(out)


def _affordable(oracle, cost, reserve=0):
    left = oracle.budget_left()
    return left is None or cost + reserve <= left


# --------------------------------------------------------------------------- #
# primitive 1: the exact zero test
# --------------------------------------------------------------------------- #
def _is_identity(oracle, word, reserve=0):
    """True if the product of element**power over ``word`` is the identity.

    Two stages: a short batch rejects a non-identity word at ~6 shots, and only
    a word that survives it pays for the long batch that accepts.
    """
    word = [(h, k) for h, k in word if k]
    if not word:
        return True
    qc = QuantumCircuit(1, 1)
    qc.h(0)
    for h, k in word:
        qc.append(oracle_gate(h, k), [0])
    qc.h(0)
    qc.measure(0, 0)
    per = oracle.cost(qc, 1)
    if not _affordable(oracle, per * QUICK_SHOTS, reserve):
        return False
    ones = oracle.run(qc, shots=QUICK_SHOTS).get("1", 0)
    if ones > ZERO_ONES:
        return False
    if not _affordable(oracle, per * FULL_SHOTS, reserve):
        return False
    ones += oracle.run(qc, shots=FULL_SHOTS).get("1", 0)
    return ones <= ZERO_ONES


# --------------------------------------------------------------------------- #
# primitive 2: single-qubit iterative phase estimation
# --------------------------------------------------------------------------- #
def _ipe(oracle, names, bits, reserve=RESERVE):
    """Estimate chi(h) for each name, all under one character.

    Returns one integer m per name with chi(h) ~ m / 2**bits, or None if the run
    does not fit in the remaining budget.
    """
    qc = QuantumCircuit(1, bits * len(names))
    for slot, name in enumerate(names):
        off = slot * bits
        for j in range(bits):
            qc.reset(0)
            qc.h(0)
            qc.append(oracle_gate(name, 1 << (bits - 1 - j)), [0])
            for i in range(max(0, j - FEEDBACK_WINDOW), j):
                with qc.if_test((qc.clbits[off + i], 1)):
                    qc.p(-2.0 * math.pi * 2.0 ** (i - j - 1), 0)
            qc.h(0)
            qc.measure(0, off + j)
    if not _affordable(oracle, oracle.cost(qc, 1), reserve):
        return None
    out = next(iter(oracle.run(qc, shots=1)))[::-1]
    return [
        sum(1 << j for j in range(bits) if out[slot * bits + j] == "1")
        for slot in range(len(names))
    ]


# --------------------------------------------------------------------------- #
# part 1: the order of g
# --------------------------------------------------------------------------- #
def _minimise(oracle, L, reserve=MIN_RESERVE):
    """Given g^L = e, divide primes out of L to reach the true ord(g).

    Writing M = L/n, the test g^(L/p^f) = e holds exactly when p^f divides M, so
    it is monotone in f and the largest valid f is found by bisection in
    O(log f) tests rather than f of them.  That matters: a polluted lcm can sit
    2^18 above n, and stepping down one prime at a time both costs far more
    than the budget allows and stops silently in the middle, which is how a
    non-minimal order used to escape.
    """
    for p in _factor(L):
        high = 0
        while L % p**(high + 1) == 0:
            high += 1
        lo, hi = 0, high  # f = 0 holds by assumption
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if _is_identity(oracle, [("g", L // p**mid)], reserve):
                lo = mid
            else:
                hi = mid - 1
        if lo:
            L //= p**lo
            _log(f"minimise: dropped {p}^{lo} -> {L}")
    return L


def _order(oracle):
    """ord(g), proved by the zero test.  Falls back to the best divisor found."""
    acc = 1
    for sweep, bits in enumerate(ORDER_LADDER):
        m = _ipe(oracle, ["g"], bits)
        if m is None:
            break
        theta = Fraction(m[0], 1 << bits)
        approx = theta.limit_denominator(min(MAX_ORDER, 1 << ((bits - 24) // 2)))
        # Self-certification: only a sweep with no cascaded bit explains theta
        # to very nearly the full precision it was measured at.  A cascaded
        # sweep is near-random, and its best small-denominator rational sits
        # ~2^-(bits-24) away, so this admits it with probability only ~2^-12.
        if abs(theta - approx) > Fraction(1, 1 << (bits - ORDER_SLACK)):
            _log(f"sweep {sweep} ({bits}b): rejected (cascaded)")
            continue
        cand = math.lcm(acc, approx.denominator)
        _log(f"sweep {sweep} ({bits}b): q={approx.denominator} lcm={cand}")
        if cand == 1 or cand > MAX_ORDER:
            continue
        if cand != acc and _is_identity(oracle, [("g", cand)], RESERVE):
            return _minimise(oracle, cand)
        acc = cand  # a proper divisor of n: fold in another sweep

    # Backstop: acc divides n but no sweep completed it.  n = acc*m for some m
    # whose prime factors we missed, and that m is small in practice.
    for mult in range(2, MULT_SEARCH + 1):
        cand = acc * mult
        if cand > MAX_ORDER:
            break
        if _is_identity(oracle, [("g", cand)], RESERVE // 2):
            return _minimise(oracle, cand)
    return max(acc, 1)


# --------------------------------------------------------------------------- #
# part 2: the three discrete logs, and forgery
# --------------------------------------------------------------------------- #
def _split(m, bits, n):
    """Read chi = m/2**bits as t/n; return (t, distance to the nearest t/n)."""
    t, rem = divmod(m * n, 1 << bits)
    frac = rem / (1 << bits)
    if frac > 0.5:
        t, frac = t + 1, 1.0 - frac
    return t % n, frac


def _candidates(pairs, n):
    """Exponents consistent with the (s, t) pairs, most-supported first.

    A single character with gcd(s, n) = 1 pins a down on its own, but when it
    does not, several characters together still do: each contributes
    a = t*s^-1 (mod n/gcd(s, n)), and the intersection closes in on one value.
    The solution set comes back as x0 modulo a period; once that period reaches
    n there is at most one exponent below n in it, namely x0 itself.
    """
    out = []
    subsets = [pairs]
    if len(pairs) > 2:  # one corrupted pair should not veto the rest
        subsets += [pairs[:i] + pairs[i + 1:] for i in range(len(pairs))]
    for sub in subsets:
        joint = _congruences([(s, t, n) for s, t in sub])
        if joint is not None and joint[1] >= n and joint[0] < n and joint[0] not in out:
            out.append(joint[0])
    for s, t in reversed(pairs):  # newest character first
        if math.gcd(s, n) == 1:
            a = t * pow(s, -1, n) % n
            if a not in out:
                out.append(a)
    return out


def _dlogs(oracle, n):
    """Return ({"A": a, "B": b, "C": c}, forged).  Short-circuits on a forgery."""
    names = ["A", "B", "C"]
    bits = n.bit_length() + LOG_SLACK
    pairs = {k: [] for k in names}
    tried = {k: set() for k in names}
    found = {k: None for k in names}

    for shot in range(LOG_SHOTS):
        todo = [k for k in names if found[k] is None]
        if not todo:
            break
        m = _ipe(oracle, ["g"] + todo, bits, reserve=RESERVE // 2)
        if m is None:
            break
        s, res = _split(m[0], bits, n)
        if res > LOG_RESIDUAL:
            _log(f"dlog shot {shot}: g cascaded, discarding")
            continue
        for name, mi in zip(todo, m[1:]):
            t, rt = _split(mi, bits, n)
            if rt <= LOG_RESIDUAL:
                pairs[name].append((s, t))
                continue
            # Either this sub-sweep cascaded, or chi(X) is genuinely not a
            # multiple of 1/n -- and the latter happens only when X^n != e.
            # The exact test tells the two apart; guessing conflates a 9%
            # cascade rate with a forgery.
            if not _is_identity(oracle, [(name, n)], RESERVE // 4):
                _log(f"dlog shot {shot}: {name}^n != e, so {name} is not in <g>")
                return found, True
            _log(f"dlog shot {shot}: {name} cascaded, discarding")

        for name in todo:
            budget_tests = 2  # cap confirmations per shot
            for a in _candidates(pairs[name], n):
                if a in tried[name] or budget_tests == 0:
                    continue
                tried[name].add(a)
                budget_tests -= 1
                if _is_identity(oracle, [(name, 1), ("g", -a)], RESERVE // 4):
                    found[name] = a
                    break
            _log(f"dlog shot {shot}: {name} found={found[name]} "
                 f"pairs={len(pairs[name])} tried={len(tried[name])}")
            # A genuine power of g confirms on the first candidate a clean
            # character offers.  Two independent characters that each yield a
            # candidate the zero test refuses means X is not in <g>.
            if found[name] is None and len(pairs[name]) >= 3 and len(tried[name]) >= 2:
                return found, True
    return found, any(v is None for v in found.values())


# --------------------------------------------------------------------------- #
# submission entry point
# --------------------------------------------------------------------------- #
def solve(oracle) -> dict:
    n = _order(oracle)
    _log("order =", n, "| queries used", oracle.queries_used)
    if n <= 1:  # g = e, so the only honest broadcasts are e as well
        forged = any(not _is_identity(oracle, [(x, 1)]) for x in ("A", "B", "C"))
        return {"order": 1, "secret": "none" if forged else 0}

    logs, forged = _dlogs(oracle, n)
    _log("logs =", logs, "forged =", forged, "| queries used", oracle.queries_used)
    if forged or any(v is None for v in logs.values()):
        return {"order": n, "secret": "none"}
    return {"order": n, "secret": logs["A"] * logs["B"] * logs["C"] % n}
