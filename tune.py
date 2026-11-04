"""Compare solver configurations on score and, more importantly, failure rate.

Order correctness is worth 1.0 vs 0.0 and the secret 1.0 vs 0.3, while leftover
budget is worth at most 0.5 -- so a config that wastes queries but always lands
beats a frugal one that sometimes misses.
"""

import itertools
import math
import sys

from silentgroup import MockGroup, Oracle, load_public
from silentgroup.check import score

import solver

BUDGET = 80_000
CASES = [(n, load_public(n)) for n in ("public_1", "public_2")] + [
    ("cyclic-prime-2^30", MockGroup([10**9 + 7], seed=5)),
    ("cyclic-prime-2^56", MockGroup([72057594037927931], seed=17)),
    ("cyclic-2^20x6", MockGroup([2**20, 6], seed=1)),
    ("big-cyclic-2^59", MockGroup([576460752303423487], seed=3)),
    ("3-factor-lowphi", MockGroup([2**20, 2**20, 60], seed=11)),
    ("smooth-2^40", MockGroup([2**40, 2**18], seed=23)),
    ("ord-g-small", MockGroup([2**30, 2**20], g=(0, 1), seed=4)),
    ("tiny", MockGroup([7, 5], seed=2)),
    ("forged-B", MockGroup([12345678, 60], forged="B", seed=5)),
    ("forged-A", MockGroup([2**31, 2**10], forged="A", seed=9)),
    ("forged-C-lowphi", MockGroup([2**20, 2**20, 60], forged="C", seed=11)),
    ("forged-AB", MockGroup([2**25, 2**12], forged="AB", seed=31)),
]

LADDERS = {
    "64-88-112": (64, 88, 112) + (144,) * 5,
    "72only": (72,) + (144,) * 5,
    "none": (144,) * 6,
    "88only": (88,) + (144,) * 5,
}


def run(ladder, reserve, reps):
    solver.ORDER_LADDER = LADDERS[ladder]
    solver.RESERVE = reserve
    tot, bad_order, bad_secret, runs, qs = 0.0, 0, 0, 0, []
    for _, inst in CASES:
        for _ in range(reps):
            orc = Oracle(inst, budget=BUDGET)
            try:
                got = solver.solve(orc)
            except Exception:
                got = None
            p = score(got, inst.answer, BUDGET, orc.queries_used)
            if p == 0.0:
                bad_order += 1
            elif p < 1:
                bad_secret += 1
            tot += p
            runs += 1
            qs.append(orc.queries_used)
    return tot / runs, bad_order, bad_secret, runs, sum(qs) / len(qs)


if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    print(f"{'ladder':12s} {'reserve':>8s} {'mean':>7s} {'badord':>7s} "
          f"{'badsec':>7s} {'runs':>6s} {'avg q':>8s}")
    for ladder, reserve in itertools.product(LADDERS, (16000, 24000, 30000)):
        mean, bo, bs, runs, aq = run(ladder, reserve, reps)
        print(f"{ladder:12s} {reserve:8d} {mean:7.3f} {bo:7d} {bs:7d} "
              f"{runs:6d} {aq:8.0f}", flush=True)
