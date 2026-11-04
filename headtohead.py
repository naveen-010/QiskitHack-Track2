"""Head-to-head on the two surviving configs, at enough reps to mean something.

The sweep in tune.py showed that dropping the cheap rungs is clearly worse, but
could not separate the remaining configs: at 56 runs, 0 vs 1 vs 2 failures is
noise.  This runs far more reps on fewer configs and reports a Wilson upper
bound on the failure rate, which is the number that actually matters -- a wrong
order scores 0 for that meeting.
"""

import math
import sys

from tune import CASES, LADDERS, BUDGET
from silentgroup import Oracle
from silentgroup.check import score

import solver

CONFIGS = [
    ("64-88-112", 16000),
    ("64-88-112", 30000),
    ("72only", 16000),
]


def wilson_upper(fails, n, z=1.96):
    """Upper end of a 95% CI for the failure rate."""
    if n == 0:
        return 1.0
    p = fails / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return min(1.0, centre + half)


def run(ladder, reserve, reps):
    solver.ORDER_LADDER = LADDERS[ladder]
    solver.RESERVE = reserve
    tot, bad_order, bad_secret, runs, q = 0.0, 0, 0, 0, 0
    worst = []
    for name, inst in CASES:
        for _ in range(reps):
            orc = Oracle(inst, budget=BUDGET)
            try:
                got = solver.solve(orc)
            except Exception as exc:
                got = None
                worst.append(f"{name}: raised {type(exc).__name__}")
            p = score(got, inst.answer, BUDGET, orc.queries_used)
            if p == 0.0:
                bad_order += 1
                worst.append(f"{name}: order wrong ({got})")
            elif p < 1:
                bad_secret += 1
                worst.append(f"{name}: secret wrong ({got})")
            tot += p
            runs += 1
            q += orc.queries_used
    return tot / runs, bad_order, bad_secret, runs, q / runs, worst


if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    print(f"{'ladder':12s} {'reserve':>8s} {'mean':>7s} {'bad_ord':>8s} "
          f"{'bad_sec':>8s} {'fail%<=':>8s} {'runs':>6s} {'avg q':>8s}", flush=True)
    for ladder, reserve in CONFIGS:
        mean, bo, bs, n, aq, worst = run(ladder, reserve, reps)
        print(f"{ladder:12s} {reserve:8d} {mean:7.3f} {bo:8d} {bs:8d} "
              f"{100*wilson_upper(bo + bs, n):7.1f}% {n:6d} {aq:8.0f}", flush=True)
        for w in worst[:4]:
            print(f"      ! {w}", flush=True)
