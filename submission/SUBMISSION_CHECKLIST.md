# Submission checklist — `submission/my_solver.py`

* Date: 2026-11-04
* Frozen file sha256: `2b0b1b4da9104cf93f822b84e79ba38970fb3b7a7ad015d6deba45036917ba8f`
  (`results/frozen_sha256.txt`; re-check with `sha256sum submission/my_solver.py`)
* Repo: `ProAltro/qff-hackathon` @ `89265613ec267045e5f36be745fb2fd4d029dd10` (unmodified `silentgroup/`)
* Environment (`results/environment.txt`): Windows 11, Python 3.14.6, qiskit 2.6.1, numpy 2.5.4
  (+ pytest 9.2.0, matplotlib 3.11.2, scipy 1.18.1 for tests/plots only — not used by the solver)

## Compliance audit (frozen file)
| Rule | Status | How checked |
|---|---|---|
| `solve(oracle) -> {"order": int, "secret": int or "none"}` | ✅ | `answer()` builds the dict; tests assert types |
| Imports only stdlib / qiskit / silentgroup | ✅ | `grep -n import` → `math, random, time, fractions, qiskit.QuantumCircuit, silentgroup.oracle_gate` (numpy not needed) |
| Oracle used only via `run`, `cost`, `budget_left` | ✅ | `grep -n "oracle\." submission/my_solver.py` |
| No instance files / answers / `_instance` / MockGroup internals | ✅ | `grep -n "_instance\|instances\|answer\[" ` → none; the only file I/O-free code |
| No prints, no file access | ✅ | `grep -n "print(\|open("` → none |
| ≤ 16 qubits, no loops in circuits, parameters bound | ✅ | every circuit is `QuantumCircuit(1, …)`; Python unrolls; all angles are floats |
| Budget ≤ 80,000 | ✅ | `_Ctx.run` refuses runs that exceed `budget_left()`; max observed 76,413 (checker), 70,055 (regression) |
| Time < 600 s | ✅ | 400 s monotonic guard; max observed 3.5 s (checker), 10.2 s (regression) |
| `"none"` only with verified lattice | ✅ | `answer(H, trusted)`; trusted only after `verify_generators` + `complete`; 0 false forgeries in 840 regression runs |

## Commands run and outputs
```
python -m venv .venv && source .venv/Scripts/activate
pip install -e . && pip install numpy qiskit pytest matplotlib scipy
python -m silentgroup.check my_solver.py --repeat 5                 # baseline -> results/baseline.txt (10/10)
python -m pytest tests -q                                           # 96 passed
python engineering/regress.py --seeds 40 --tag final_40seeds        # results/regress_final_40seeds.txt
python -m silentgroup.check submission/my_solver.py --repeat 20     # results/final_checker.txt
```

### Final checker (`results/final_checker.txt`, frozen file, 20 repeats each)
| instance | correct | mean queries | max queries | mean pts | max time |
|---|---|---|---|---|---|
| public_1 | 20/20 | 17,025 | 19,516 | 1.393 | 3.5 s |
| public_2 | 20/20 | 40,161 | 76,413 | 1.249 | 3.5 s |

### Final regression (`results/regress_final_40seeds.txt`, 40 seeds × 21 shapes, noise 0.002)
840/840 fully correct; every shape 40/40, Wilson 95 % lower bound 0.912 ≥ 0.90; 0 false forgeries,
0 missed forgeries, 0 exceptions; max 70,055 queries; max 10.2 s.

## Re-check command
```
python -m silentgroup.check submission/my_solver.py --repeat 20
python engineering/regress.py --seeds 40 --tag recheck
```
