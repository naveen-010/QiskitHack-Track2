"""Heavy validation, weighted toward the shapes that actually break things.

The public meetings cover two group shapes; the hidden ones "vary in shape".
The hard cases are those with small phi(n)/n -- a character's s shares a factor
with n often, so no single character inverts and the exponent stays ambiguous.
A pure power of two is the worst: s is odd only half the time.
"""

import math
import sys
import time
import traceback

from silentgroup import MockGroup, Oracle, load_public
from silentgroup.check import score

import solver

BUDGET = 80_000

CASES = [(n, load_public(n)) for n in ("public_1", "public_2")] + [
    # pure powers of two: phi(n)/n = 1/2, the ambiguous-coset stress test
    ("pow2-2^40", MockGroup([2**40, 2**18], seed=23)),
    ("pow2-2^50", MockGroup([2**50], seed=41)),
    ("pow2-2^59", MockGroup([2**59], seed=43)),
    ("pow2-2^30", MockGroup([2**30, 2**7], seed=47)),
    # low phi(n)/n from many small primes
    ("lowphi-3factor", MockGroup([2**20, 2**20, 60], seed=11)),
    ("lowphi-smooth", MockGroup([2**18 * 3**5, 2**6 * 5**3], seed=53)),
    # ordinary shapes
    ("prime-2^30", MockGroup([10**9 + 7], seed=5)),
    ("prime-2^56", MockGroup([72057594037927931], seed=17)),
    ("cyclic-2^59", MockGroup([576460752303423487], seed=3)),
    ("ord-g-small", MockGroup([2**30, 2**20], g=(0, 1), seed=4)),
    ("tiny", MockGroup([7, 5], seed=2)),
    # forgeries, including over a hard shape
    ("forged-B", MockGroup([12345678, 60], forged="B", seed=5)),
    ("forged-A-pow2", MockGroup([2**31, 2**10], forged="A", seed=9)),
    ("forged-C-lowphi", MockGroup([2**20, 2**20, 60], forged="C", seed=11)),
    ("forged-AB", MockGroup([2**25, 2**12], forged="AB", seed=31)),
]


def wilson_upper(fails, n, z=1.96):
    if n == 0:
        return 1.0
    p, d = fails / n, 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return min(1.0, centre + half)


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    only = sys.argv[2] if len(sys.argv) > 2 else None
    tot = runs = bad_order = bad_secret = 0
    tot = 0.0
    slowest = 0.0
    for name, inst in CASES:
        if only and only not in name:
            continue
        pts, qs, fails = [], [], []
        for _ in range(reps):
            orc = Oracle(inst, budget=BUDGET)
            t0 = time.time()
            try:
                got = solver.solve(orc)
            except Exception:
                got = None
                fails.append(traceback.format_exc().strip().splitlines()[-1])
            el = time.time() - t0
            slowest = max(slowest, el)
            p = score(got, inst.answer, BUDGET, orc.queries_used)
            if el > 600:
                p, _ = 0.0, fails.append("OVER TIME LIMIT")
            if p == 0.0:
                bad_order += 1
                fails.append(f"ORDER got={got} want={inst.answer} q={orc.queries_used}")
            elif p < 1:
                bad_secret += 1
                fails.append(f"SECRET got={got} want={inst.answer} q={orc.queries_used}")
            pts.append(p)
            qs.append(orc.queries_used)
            tot += p
            runs += 1
        ok = sum(1 for x in pts if x >= 1)
        flag = "" if ok == reps else "   <-- FAILURES"
        print(f"{name:18s} n=2^{math.log2(inst.answer['order']):4.1f} "
              f"pts={sum(pts)/len(pts):.3f} q={sum(qs)//len(qs):6d} "
              f"max_q={max(qs):6d} ok={ok}/{reps}{flag}", flush=True)
        for f in dict.fromkeys(fails):
            print(f"      ! {f}", flush=True)
    bad = bad_order + bad_secret
    print(f"\nmean {tot/runs:.3f} over {runs} runs | wrong order {bad_order} | "
          f"wrong secret {bad_secret} | failure rate {100*bad/runs:.2f}% "
          f"(95% CI <= {100*wilson_upper(bad, runs):.2f}%) | slowest {slowest:.1f}s")


if __name__ == "__main__":
    main()
