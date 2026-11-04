"""Seeded benchmark shared by the combined solver and any supplied reference file.

python engineering/benchmark_combined.py --solver submission/my_solver.py --tag combined
Each solver uses the same instances, oracle seeds, noise, budget and scoring rule.
The selected suite is not an estimate of the hidden judge's instance distribution.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import multiprocessing as mp
from pathlib import Path
import statistics
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from silentgroup import MockGroup, Oracle, load_public
from silentgroup.check import score

CASES = {
    "prime_56": ([72057594037927931], {}),
    "pow2_59": ([2**59], {}),
    "prime_product": ([10**9 + 7, 10**9 + 9], {}),
    "smooth": ([2**18 * 3**5, 2**6 * 5**3], {}),
    "g_small": ([2**30, 2**20], {"g": (0, 1)}),
    "tiny_7": ([7], {}),
    "tiny_1000003": ([1000003], {}),
    "forged_B": ([12345678, 60], {"forged": "B"}),
    "forged_AB": ([10**8, 3 * 10**8], {"forged": "AB"}),
    "forged_ABC": ([2**28, 3 * 2**28], {"forged": "ABC"}),
    "skew_small_logs": ([2**55], {"g": (1,), "a": 2, "b": 3, "c": 5}),
    "skew_equal_logs": ([2**58 - 27], {"g": (1,), "a": 1, "b": 1, "c": 1}),
    "identity_g_honest": ([1009], {"g": (0,)}),
    "identity_g_forged": ([1009], {"g": (0,), "forged": "A"}),
    "zero_A": ([2**59], {"g": (1,), "a": 0}),
    "exact_2pow60": ([2**60], {"g": (1,)}),
    "prime_rank4_forgery": ([1009] * 4, {"g": (1, 0, 0, 0), "forged": "ABC"}),
    "two_long_directions": ([2**30] * 2, {"g": (1, 0), "a": 1, "b": 1, "c": 1, "forged": "A"}),
}


def init(path):
    global SOLVER
    spec = importlib.util.spec_from_file_location("bench_submission", path)
    SOLVER = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(SOLVER)


def run(task):
    name, seed, offset = task
    inst = load_public(name) if name.startswith("public_") else MockGroup(CASES[name][0], seed=seed, **CASES[name][1])
    assert inst.group.order <= 2**60
    oracle = Oracle(inst, budget=80000, seed=offset + seed, noise=0.002)
    start = time.monotonic()
    result, error = None, ""
    try:
        result = SOLVER.solve(oracle)
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
    elapsed = time.monotonic() - start
    return dict(case=name, seed=seed, oracle_seed=offset + seed, ok=int(result == inst.answer),
                order_ok=int(isinstance(result, dict) and result.get("order") == inst.answer["order"]),
                queries=oracle.queries_used, seconds=round(elapsed, 3),
                score=score(result, inst.answer, 80000, oracle.queries_used) if elapsed <= 600 else 0,
                error=error, result=json.dumps(result), answer=json.dumps(inst.answer))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", default=str(ROOT / "submission/my_solver.py"))
    parser.add_argument("--tag", default="combined")
    parser.add_argument("--seeds", type=int, default=8)
    parser.add_argument("--seed0", type=int, default=0)
    parser.add_argument("--public-repeats", type=int, default=20)
    parser.add_argument("--offset", type=int, default=10000)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--only", nargs="+")
    args = parser.parse_args()
    names = args.only or ["public_1", "public_2", *CASES]
    tasks = [(name, seed, args.offset) for name in names
             for seed in range(args.seed0, args.seed0 + (args.public_repeats if name.startswith("public_") else args.seeds))]
    rows, start = [], time.monotonic()
    target = ROOT / "results" / f"benchmark_{args.tag}"
    target.parent.mkdir(exist_ok=True)
    with mp.Pool(args.workers, initializer=init, initargs=(str(Path(args.solver).resolve()),)) as pool:
        with target.with_suffix(".csv").open("w", newline="") as f:
            writer = None
            for row in pool.imap_unordered(run, tasks):
                rows.append(row)
                if writer is None:
                    writer = csv.DictWriter(f, fieldnames=list(row))
                    writer.writeheader()
                writer.writerow(row)
                f.flush()
                if not row["ok"]:
                    print("FAIL", json.dumps(row), flush=True)
                if len(rows) % 20 == 0:
                    print(f"{len(rows)}/{len(tasks)} runs, {time.monotonic() - start:.0f}s", flush=True)
    summaries = []
    for name in names:
        rs = [r for r in rows if r["case"] == name]
        item = dict(case=name, runs=len(rs), correct=sum(r["ok"] for r in rs),
                    correct_order=sum(r["order_ok"] for r in rs),
                    mean_queries=round(statistics.mean(r["queries"] for r in rs)),
                    max_queries=max(r["queries"] for r in rs),
                    mean_score=round(statistics.mean(r["score"] for r in rs), 6),
                    max_seconds=max(r["seconds"] for r in rs))
        summaries.append(item)
        print(json.dumps(item), flush=True)
    report = dict(solver=str(Path(args.solver).resolve()),
                  solver_sha256=hashlib.sha256(Path(args.solver).read_bytes()).hexdigest(),
                  offset=args.offset, seed0=args.seed0,
                  runs=len(rows), correct=sum(r["ok"] for r in rows),
                  mean_queries=statistics.mean(r["queries"] for r in rows),
                  mean_score=statistics.mean(r["score"] for r in rows), cases=summaries)
    target.with_suffix(".json").write_text(json.dumps(report, indent=2))
    print("TOTAL", json.dumps({k: v for k, v in report.items() if k != "cases"}), flush=True)


if __name__ == "__main__":
    main()
