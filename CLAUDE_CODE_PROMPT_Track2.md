# PROMPT FOR CLAUDE CODE (paste everything below the line)

Before pasting: put `QFF_Track2_FINAL_Plan.md` and `my_solver.py` inside your cloned `qff-hackathon` folder, open the terminal there, and start `claude`. Suggested model: Opus at medium effort. Remember the deadline is today; the plan says to submit the baseline first.

---

You are my senior quantum-software engineer. We are in a clone of `ProAltro/qff-hackathon` (Qiskit Fall Fest 2026, BITS Pilani, Track 2 "The Silent Group"). Your job: take the existing reference solver to a hardened, frozen, submission-ready single file, with tests, measurements, open-track experiments and report material, working end to end with minimal back-and-forth.

## 0. Read first (do this before any code)
1. Read `QFF_Track2_FINAL_Plan.md` completely. It is the spec for design decisions, measured results, known bugs and priorities. Where it says **[measured]** or **[source]** you may trust it; where it says **[hypothesis]** you must test it.
2. Read `my_solver.py` completely, then `silentgroup/simulator.py`, `oracle.py`, `program.py`, `group.py`, `check.py`, `examples/solver_template.py`. Re-verify the facts the plan relies on (one random character per shot, fixed across the circuit; both `if_test` branches charged; cost = `max(1, bit_length(|k|))`; noise σ = 0.002·√cost turns). If anything contradicts the plan, stop and tell me.
3. Summarise in ≤ 15 lines what the solver does and what you plan to change. Then proceed without waiting for approval unless a rule below says to ask.

## 1. Hard rules (never violate)
* The submitted file `submission/my_solver.py` must define `solve(oracle) -> {"order": int, "secret": int | "none"}`, import only the Python standard library, numpy, qiskit and `silentgroup` (as the template does), and use the hidden meeting **only** through the `oracle` object (`run`, `cost`, `budget_left`, `elements`).
* Inside `solve`, never read `silentgroup/instances/`, public answers, `MockGroup` internals or `oracle._instance`. Answers may be used only in test/harness code outside the submission file.
* ≤ 16 qubits per circuit; no loops inside circuits (unroll in Python); all gate parameters bound.
* Budget ≤ 80,000 queries and wall-clock < 600 s per solve. Add a time guard (~400 s) and a budget guard.
* Output `"none"` only on positive evidence that some log does not exist **and** the lattice was verified. Never output `"none"` merely because of a timeout or failed recovery.
* Do not weaken verification to save queries. Reliability beats query savings (the whole solve is ≈14k of 80k; the maximum gain from squeezing further is ≈0.04 points).
* Do not tune only to the two public instances. Use many seeds and shapes.
* If a permission prompt or sandbox blocks a command, do not work around it; tell me what you were trying to do and why, and continue with other tasks.
* Do not commit or push unless I ask. Keep all work in the working tree.

## 2. Project layout to create
```
submission/my_solver.py        # final frozen single file (self-contained)
engineering/regress.py         # mock-group regression harness
engineering/experiments.py     # open-track experiments + plots
tests/test_oracle_semantics.py # F0: one shot = one character; measure/reset/if_test; cost formula
tests/test_small_exact.py      # F1: noise=0 exact answers on tiny groups
tests/test_lattice.py          # unit tests for lll / hnf / extract
results/                       # CSVs, PNG plots, logs
REPORT_NOTES.md                # material for the 6-page report and the 8-min talk
```
Keep `my_solver.py` at the repo root as the untouched baseline until the final file passes everything.

## 3. Work packages, in this order

**WP0 – Baseline (≤ 10 min).** Create/activate a venv, `pip install -e .`, make sure qiskit and numpy are installed. Run `python -m silentgroup.check my_solver.py --repeat 5`. Record the output in `results/baseline.txt`. Expected: all correct, ≈14k queries, ≈1.4 pts. Report any deviation.

**WP1 – Tests of facts (F0, F1).**
* `test_oracle_semantics.py`: (a) with a group where phases differ per character, a single circuit that measures one qubit twice with the same oracle gate gives perfectly correlated outcomes in the noiseless case while two separate `run(shots=1)` calls do not; (b) `measure → reset → if_test` feedback works; (c) `oracle.cost` equals `shots × Σ max(1, bit_length(|k|))` and charges both branches of an `if_test`; (d) identity relation gives P(1)≈0, non-identity ≈0.5.
* `test_small_exact.py`: with `noise=0`, the solver returns the exact answer on ℤ₈, ℤ₁₅, ℤ₈×ℤ₃, ℤ₄×ℤ₈, ℤ₇×ℤ₁₁ and on forged variants (use `MockGroup(..., forged="C")` etc., needing non-cyclic moduli).
* `test_lattice.py`: LLL output is unimodularly equivalent to input (same determinant) and reduced; HNF is upper triangular with positive pivots and generates the same lattice; `extract` returns the right order/logs on hand-built lattices, including `None` for an element outside ⟨g⟩.

**WP1b – Regression harness `engineering/regress.py`.** For a list of group shapes × N seeds (default 20), build `MockGroup`, run `solve` with `Oracle(inst, budget=80000, seed=...)`, and log: correct order, correct secret, false forgery, missed forgery, queries, seconds, exception. Print a table with Wilson 95 % intervals. Keep every |G| ≤ 2^60 (an earlier invalid mock with |G|≈2^61 stalled the solver). Shapes must include: cyclic ≈2^58; ℤ_{10^9+7}×ℤ_{10^9+9}; ℤ_{2^20}×ℤ_6×ℤ_{3^8}; ℤ_{2^55}; ℤ_{2^29}²; ℤ_{1000003}; ℤ_{2^30}×ℤ_{3^18}; several tiny groups (orders 7, 64, 1001, 65537); forged variants with 1, 2 and 3 forgeries using non-cyclic moduli with |G| ≤ 2^58; the shapes of `public_1` and `public_2`; and skew cases where a, b, c are small (pass `a=`, `b=`, `c=` to `MockGroup`).

**WP2 – P0 fixes in the solver (work in `submission/my_solver.py`, copied from the baseline).**
1. *Stall guard.* Cap the leave-one/two-out consensus search to the 10 most recent samples, bound total LLL calls, add a wall-clock guard (~400 s) and a budget guard. On exhaustion return a best-effort answer (never an invented `"none"`).
2. *Tiny-group bug.* Known failures: `MockGroup([1000003], seed=1)` (wrong order 2048006144) and `seed=4` (right order, secret `"none"`), with `Oracle(inst, seed=100+seed)`. Reproduce, print per-sample-count HNF diagonals, determinant, number of "good" short vectors and tails, then find the root cause. Hypotheses to test, not assume: (a) more than four near-relations make `good[:4]` pick a dependent or wrong set, so select a maximal independent subset or the 4 shortest with full rank and then check that remaining good vectors are in the lattice they span; (b) tolerance `2^(b-12)` badly tuned, so scale the tolerance with the observed norms or `log2` of the determinant; (c) verification threshold too weak for small groups. Fix the cause, not the symptoms.
3. *Rank-2 p-groups* (ℤ_{2^29}² had 1/10 wrong). Require more samples (≥ 6) when the HNF quotient is non-cyclic, and keep sampling until the lattice is unchanged after two further samples if budget allows (budget guard permitting).

After each fix, rerun `engineering/regress.py` and the checker. Record before/after in `results/`.

**WP3 – P1 robustness.**
1. *Precision ladder.* If verification fails, or the minimality test fails, or no lattice is found after the sample cap, retry with b = 44, then 52, discarding old samples of lower precision. Count the extra queries in the budget guard.
2. *Skewed lattices.* Add regression cases with very small logs (e.g. a=2, b=3, c=5 in a 2^55 group). If `n` can come out as a proper multiple of the true order, add a certification step: factor `n` with Pollard-rho (pure Python, deadline-capped) and test `g^(n/p) ≠ e` for every prime `p | n` with the 14-shot Hadamard test, within a query cap. Make sure certification is skipped when it is not needed (no extra queries in the common case).
3. *Sequential verification.* Replace fixed 16-shot tests with a sequential test that stops early when it sees ≥ 3 ones (reject) or completes 16 shots with ≤ 2 ones (accept); keep the same error rates. Only keep this change if the regression shows no loss of reliability.

**WP4 – Acceptance criteria for the frozen file.**
* Official checker: `--repeat 20` on both public instances, all correct.
* Regression: ≥ 95 % correct (lower Wilson bound ≥ 90 %) on every shape listed above with 30 seeds each; zero exceptions; max wall time < 60 s; median queries ≲ 16k on ordinary shapes.
* No run exceeds 80,000 queries. No `"none"` produced on valid instances (false forgery rate = 0 in the regression).
* If any criterion fails, document why and what you tried; do not hide failures, and do not lower the criteria.

**WP5 – Open-track experiments (`engineering/experiments.py`, results in `results/`, PNGs).** Title: *Noise–precision–query laws for relation-lattice recovery in a noisy black-box abelian group.*
1. *Cost law.* For log₂|G| ∈ {20, 30, 40, 50, 58}, find the smallest b with ≥ 95 % success at 4 samples; compare with the heuristic b ≈ 16 + log₂|H|/t; plot predicted vs measured queries.
2. *Noise sweep.* `Oracle(..., noise=ν)` for ν ∈ {0.002, 0.004, 0.008, 0.02, 0.05}; plot success rate and queries; test the hypothesis that corrupted-shot rate ∝ ν².
3. *Mitigation.* Repeat the cheap high-weight stages (m ≤ 19) twice with a vote; report success per query vs the baseline.
4. *Group structure.* Add a function (outside the submission path, or inside only if it costs no queries) that outputs the full structure ℤ_{d1}×… from the Smith normal form of Λ; test it on mock groups.
Use ≥ 20 seeds per point, Wilson intervals, and say clearly which statements are theorem, simulator fact, or empirical observation. Report negative results.

**WP6 – Freeze and audit.** Merge everything into `submission/my_solver.py`; verify imports; no prints; no file access; run the checker on the frozen file with `--repeat 20`; run the full regression once more. Then produce `submission/SUBMISSION_CHECKLIST.md` with the exact commands run, outputs, date, environment (`python --version`, `pip freeze`).

**WP7 – Report and talk material (`REPORT_NOTES.md`).** Fill in: the 6-page report outline from the plan with concrete numbers and figures; the 8-minute talk script with timings; a "what failed" section (b=28/32, stability-rule overspend, Shor-style 120-bit cost, anything new); a function-by-function explanation table; 15 viva questions with short answers. Write it so a student who did not write the code can explain and modify it without AI.

## 4. Way of working
* Use a task list; mark items complete only when their acceptance test passes.
* Make small, reviewable edits. After each work package, give me a 3–5 line summary: what changed, test results, anything unexpected.
* Run long experiments in the background and keep the main loop moving. Keep each solver run under the wall-clock limit during tests.
* Prefer simple, explainable code over cleverness. Every non-obvious line gets a one-line comment, since I must defend it in a viva without AI.
* If you must choose between a faster and a more reliable option, choose reliable and tell me the cost.
* Stop and ask me only if: a rule in §1 would be violated, the plan contradicts the code, an acceptance criterion cannot be met after two honest attempts, or a command is blocked.

## 5. Final response
When finished, reply with: (1) the pass/fail table for WP4, (2) the final query/time statistics, (3) the list of known remaining weaknesses, (4) the exact command to re-check the frozen file, (5) the paths of all deliverables. Be honest about anything unproven.
