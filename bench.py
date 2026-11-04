"""Score solver.py the way the judge will, over many instances and repeats."""
import sys, time, traceback, math
from silentgroup import MockGroup, Oracle, load_public
from silentgroup.check import score
import solver

BUDGET = 80_000
cases = [(n, load_public(n)) for n in ("public_1", "public_2")] + [
    ("cyclic-prime-2^30",  MockGroup([10**9 + 7], seed=5)),
    ("cyclic-prime-2^56",  MockGroup([72057594037927931], seed=17)),
    ("cyclic-2^20x6",      MockGroup([2**20, 6], seed=1)),
    ("big-cyclic-2^59",    MockGroup([576460752303423487], seed=3)),
    ("3-factor-lowphi",    MockGroup([2**20, 2**20, 60], seed=11)),
    ("smooth-2^40",        MockGroup([2**40, 2**18], seed=23)),
    ("ord-g-small",        MockGroup([2**30, 2**20], g=(0, 1), seed=4)),
    ("tiny",               MockGroup([7, 5], seed=2)),
    ("forged-B",           MockGroup([12345678, 60], forged="B", seed=5)),
    ("forged-A",           MockGroup([2**31, 2**10], forged="A", seed=9)),
    ("forged-C-lowphi",    MockGroup([2**20, 2**20, 60], forged="C", seed=11)),
    ("forged-AB",          MockGroup([2**25, 2**12], forged="AB", seed=31)),
]
reps = int(sys.argv[1]) if len(sys.argv) > 1 else 3
tot, runs, allfail = 0.0, 0, []
for name, inst in cases:
    pts, qs, fails = [], [], []
    for _ in range(reps):
        orc = Oracle(inst, budget=BUDGET)
        t0 = time.time()
        try:
            got = solver.solve(orc)
        except Exception:
            got = None; fails.append(traceback.format_exc().strip().splitlines()[-1])
        el = time.time() - t0
        p = score(got, inst.answer, BUDGET, orc.queries_used)
        if el > 600: p, _ = 0.0, fails.append("TIMEOUT")
        if p < 1: fails.append(f"p={p} got={got} want={inst.answer}")
        pts.append(p); qs.append(orc.queries_used); tot += p; runs += 1
    ok = sum(1 for x in pts if x >= 1)
    print(f"{name:20s} n=2^{math.log2(inst.answer['order']):4.1f} "
          f"pts={sum(pts)/reps:.3f} q={sum(qs)//reps:6d} ok={ok}/{reps}")
    for f in fails[:2]:
        print(f"      ! {f}")
    if ok < reps: allfail.append(name)
print(f"\nmean score {tot/runs:.3f} over {runs} runs; imperfect: {allfail or 'none'}")
