# Combined submission checklist

Current branch: `combined-solver` (based on Pulkit, with Naveen's order/log ideas).

Submit **`submission/my_solver.py`**. The root **`my_solver.py`** is byte-identical.
Both are standalone files defining `solve(oracle) -> dict`.

SHA256 (both files): `5d3635eb8439b62a07ef19d70eb1922976ee451e50f1b549570dc94d15c81e85`.
The two files use LF endings, enforced by `.gitattributes`.
`results/combined_sha256.txt` records their checksums.

## Validation on the final file

| Check | Result |
|---|---:|
| Unit tests and reproduced failure regressions | 118 passed |
| Same seeded comparison used for both original branches | 184/184 correct |
| Independent-seed comparison | 328/328 correct |
| Original 21-shape, 40-seed regression | 840/840 correct |
| Official unseeded public checker, 20 repeats per instance | 40/40 correct |
| Total solver runs | 1,392/1,392 correct |
| Maximum observed queries | 71,094 / 80,000 |
| Maximum observed wall time | 8.991 s / 600 s |

Environment: Python 3.14.6, Qiskit 2.6.1, NumPy 2.5.4; Windows.

The solver uses one recycled qubit. Runtime imports are only `math`, `random`,
`time`, `fractions`, `qiskit.QuantumCircuit`, and `silentgroup.oracle_gate`.
Device interaction uses only `run`, `cost`, and `budget_left`; no instance access,
file access, networking, nonstandard dependencies, or output printing.

Finite-shot checks and phase reconstruction remain probabilistic. The selected
suite is not the unknown hidden-instance distribution. On total exhaustion the
required interface receives a best-effort integer guess, not a claimed forgery.

## Reproduce

```powershell
.\.venv\Scripts\python.exe -m silentgroup.check submission/my_solver.py --repeat 20
.\.venv\Scripts\python.exe -m pytest tests -q
.\.venv\Scripts\python.exe engineering/benchmark_combined.py --tag combined_final
.\.venv\Scripts\python.exe engineering/benchmark_combined.py --tag combined_final_fresh --seed0 8 --seeds 16 --public-repeats 20 --offset 20000
.\.venv\Scripts\python.exe engineering/regress.py --seeds 40 --tag combined_final_stress --workers 4
```

See `COMBINED_SOLVER.md` for the derivations, optimizations, and limitations.
See `results/combined_validation.json` for machine-readable measurements and
original branch provenance. `REPORT_NOTES.md`, `results/final_checker.txt`, and
`results/frozen_sha256.txt` are historical Pulkit-version artifacts, not validation
of this combined implementation.
