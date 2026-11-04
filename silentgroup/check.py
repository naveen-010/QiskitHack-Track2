"""Rehearse a submission on the public instances, scored exactly as the judge scores hidden ones.

    python -m silentgroup.check my_solver.py
    python -m silentgroup.check my_solver.py --repeat 5      # several runs: noise makes results vary
"""

from __future__ import annotations

import argparse
import importlib.util
import time
import traceback

from .group import load_public
from .oracle import BUDGET, TIME_LIMIT, Oracle


def score(result, answer, budget: int, used: int) -> float:
    """Per instance: 0 if wrong; 0.3 for the order alone; 1 + up to 0.5 for budget left."""
    if not isinstance(result, dict) or result.get("order") != answer["order"]:
        return 0.0
    if result.get("secret") != answer["secret"]:
        return 0.3
    return 1 + 0.5 * (budget - used) / budget


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("solver", help="your submission: a .py file defining solve(oracle)")
    ap.add_argument("--instances", nargs="+", default=["public_1", "public_2"])
    ap.add_argument("--repeat", type=int, default=1, help="runs per instance")
    ap.add_argument("--budget", type=int, default=BUDGET)
    args = ap.parse_args()

    spec = importlib.util.spec_from_file_location("submission", args.solver)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    for name in args.instances:
        instance = load_public(name)
        for _ in range(args.repeat):
            oracle = Oracle(instance, budget=args.budget)
            t0 = time.monotonic()
            try:
                result, error = module.solve(oracle), None
            except Exception:
                result, error = None, traceback.format_exc().strip().splitlines()[-1]
            elapsed = time.monotonic() - t0
            pts = score(result, instance.answer, args.budget, oracle.queries_used)
            if elapsed > TIME_LIMIT:
                pts, error = 0.0, f"over the {TIME_LIMIT}s time limit"
            status = error or ("correct" if pts >= 1 else f"wrong: got {result}, expected {instance.answer}")
            print(f"{name}: {pts:.3f} pts  queries={oracle.queries_used}  {elapsed:.1f}s  {status}")


if __name__ == "__main__":
    main()
