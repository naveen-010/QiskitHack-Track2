"""Mock-group regression harness for the Silent Group solver.

    python scripts/regression.py                          # my_solver.py, 20 seeds
    python scripts/regression.py --seeds 40 --tag combined_final_stress
    python scripts/regression.py --only tiny_1000003 --seeds 10

Uses answers ONLY here, outside the submission file. Every |G| <= 2^60.
Writes results/regress_<tag>.csv (one row per run) and prints a summary table with
Wilson 95 % intervals for "fully correct".
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import multiprocessing as mp
import statistics
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from silentgroup import MockGroup, Oracle  # noqa: E402

# name -> (moduli, MockGroup keyword arguments)
SHAPES = {
    "cyclic_2^58":        ([(1 << 58) - 27], {}),
    "p1e9_x_p1e9":        ([10**9 + 7, 10**9 + 9], {}),
    "2^20x6x3^8":         ([2**20, 6, 3**8], {}),
    "cyclic_2^55":        ([2**55], {}),
    "2^29_sq":            ([2**29, 2**29], {}),
    "tiny_1000003":       ([1000003], {}),
    "2^30x3^18":          ([2**30, 3**18], {}),
    "tiny_7":             ([7], {}),
    "tiny_64":            ([64], {}),
    "tiny_1001":          ([1001], {}),
    "tiny_65537":         ([65537], {}),
    "forged1_C":          ([2**25 * 3**3, 2**20 * 5**3], {"forged": "C"}),
    "forged2_AB":         ([10**8, 3 * 10**8], {"forged": "AB"}),
    "forged3_ABC":        ([2**28, 3 * 2**28], {"forged": "ABC"}),
    "forged_small_12x6":  ([12, 6], {"forged": "C"}),
    "public1_shape":      ([6145295, 53834531687], {}),
    "public2_shape":      ([10628367750, 12297285], {}),
    "public2_shape_forgA": ([10628367750, 12297285], {"forged": "A"}),
    "skew_2^55_a2b3c5":   ([2**55], {"a": 2, "b": 3, "c": 5}),
    "skew_2^58_a1b1c1":   ([(1 << 58) - 27], {"a": 1, "b": 1, "c": 1}),
    "skew_p1e9sq_small":  ([10**9 + 7, 10**9 + 9], {"a": 3, "b": 7, "c": 11}),
}
GROUPS = {
    "core": ["cyclic_2^58", "p1e9_x_p1e9", "2^20x6x3^8", "cyclic_2^55", "public1_shape", "public2_shape"],
    "tiny": ["tiny_7", "tiny_64", "tiny_1001", "tiny_65537", "tiny_1000003", "forged_small_12x6"],
}


def load(path):
    spec = importlib.util.spec_from_file_location("solver_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_SOLVER = None


def _init(path):
    global _SOLVER
    _SOLVER = load(path)


def run_one(task):
    """One (shape, seed) solve. Instance seed = seed, oracle seed = 100 + seed (reproducible)."""
    name, seed, noise, budget = task
    moduli, kw = SHAPES[name]
    inst = MockGroup(moduli, seed=seed, **kw)
    assert math.prod(moduli) <= 1 << 60
    orc = Oracle(inst, budget=budget, seed=100 + seed, noise=noise)
    t0 = time.perf_counter()
    err, res = "", None
    try:
        res = _SOLVER.solve(orc)
    except Exception:
        err = traceback.format_exc().strip().splitlines()[-1][:200]
    dt = time.perf_counter() - t0
    ans = inst.answer
    ok_order = isinstance(res, dict) and res.get("order") == ans["order"]
    ok_secret = ok_order and res.get("secret") == ans["secret"]
    got_none = isinstance(res, dict) and res.get("secret") == "none"
    return {
        "shape": name, "seed": seed, "noise": noise,
        "order_ok": int(ok_order), "secret_ok": int(ok_secret),
        "false_forgery": int(got_none and ans["secret"] != "none"),
        "missed_forgery": int(isinstance(res, dict) and ans["secret"] == "none" and not got_none),
        "queries": orc.queries_used, "seconds": round(dt, 3), "error": err,
        "got": repr(res), "want": repr(ans),
    }


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, mid - half), min(1.0, mid + half))


def summarise(rows, out=sys.stdout):
    by = {}
    for r in rows:
        by.setdefault(r["shape"], []).append(r)
    hdr = (f"{'shape':22s} {'n':>3s} {'ok':>3s} {'wilson95':>13s} {'ordOK':>5s} {'falseF':>6s} "
           f"{'missF':>5s} {'exc':>3s} {'medQ':>7s} {'maxQ':>7s} {'maxT':>6s}")
    print(hdr, file=out)
    allpass = True
    for name, rs in by.items():
        n = len(rs)
        k = sum(r["secret_ok"] for r in rs)
        lo, hi = wilson(k, n)
        q = [r["queries"] for r in rs]
        t = [r["seconds"] for r in rs]
        exc = sum(bool(r["error"]) for r in rs)
        ff = sum(r["false_forgery"] for r in rs)
        passed = k / n >= 0.95 and lo >= 0.90 and exc == 0 and ff == 0 and max(q) <= 80_000 and max(t) < 60
        allpass &= passed
        print(f"{name:22s} {n:3d} {k:3d} [{lo:.3f},{hi:.3f}] {sum(r['order_ok'] for r in rs):5d} "
              f"{ff:6d} {sum(r['missed_forgery'] for r in rs):5d} {exc:3d} "
              f"{statistics.median(q):7.0f} {max(q):7d} {max(t):6.1f} {'PASS' if passed else 'FAIL'}", file=out)
    return allpass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solver", default=str(ROOT / "my_solver.py"))
    ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--only", nargs="*", help="shape names or group names (core, tiny)")
    ap.add_argument("--noise", type=float, default=0.002)
    ap.add_argument("--budget", type=int, default=80_000)
    ap.add_argument("--tag", default="current")
    ap.add_argument("--workers", type=int, default=max(1, mp.cpu_count() - 2))
    ap.add_argument("--deadline", type=float, default=3600, help="global wall-clock cap (s)")
    args = ap.parse_args()

    names = []
    for x in args.only or SHAPES:
        names += GROUPS.get(x, [x])
    tasks = [(nm, s, args.noise, args.budget) for nm in names for s in range(args.seed0, args.seed0 + args.seeds)]
    rows = []
    t0 = time.time()
    with mp.Pool(args.workers, initializer=_init, initargs=(args.solver,)) as pool:
        it = pool.imap_unordered(run_one, tasks)
        for _ in tasks:
            try:
                rows.append(it.next(timeout=max(1, args.deadline - (time.time() - t0))))
            except mp.TimeoutError:
                print("GLOBAL DEADLINE hit; unfinished tasks recorded as timeouts", file=sys.stderr)
                done = {(r["shape"], r["seed"]) for r in rows}
                for nm, s, nz, _b in tasks:
                    if (nm, s) not in done:
                        rows.append({"shape": nm, "seed": s, "noise": nz, "order_ok": 0, "secret_ok": 0,
                                     "false_forgery": 0, "missed_forgery": 0, "queries": -1,
                                     "seconds": float("inf"), "error": "harness timeout", "got": "", "want": ""})
                pool.terminate()
                break
    rows.sort(key=lambda r: (names.index(r["shape"]), r["seed"]))
    (ROOT / "results").mkdir(exist_ok=True)
    path = ROOT / "results" / f"regress_{args.tag}.csv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(ROOT / "results" / f"regress_{args.tag}.txt", "w") as f:
        f.write(f"solver={args.solver} seeds={args.seed0}..{args.seed0 + args.seeds - 1} noise={args.noise}\n")
        ok = summarise(rows, f)
        for r in rows:
            if not r["secret_ok"]:
                f.write(f"FAIL {r['shape']} seed={r['seed']} got={r['got']} want={r['want']} "
                        f"q={r['queries']} err={r['error']}\n")
    ok = summarise(rows)
    for r in rows:
        if not r["secret_ok"]:
            print(f"FAIL {r['shape']} seed={r['seed']} got={r['got']} want={r['want']} q={r['queries']} {r['error']}")
    print(f"ALL PASS: {ok}   ({time.time() - t0:.0f}s, wrote {path})")


if __name__ == "__main__":
    main()
