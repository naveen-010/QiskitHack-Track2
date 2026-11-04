# The Silent Group — FINAL Plan (merged)
Qiskit Fall Fest 2026 · BITS Pilani · Track 2 · 4 Nov 2026
Merges: your "Research Engineering Plan" (governance, gates, compliance, report) + my plan and prototype (working method, measurements).
Tags: **[measured]** = run on the organisers' simulator earlier today; **[hypothesis]** = untested; **[source]** = read directly in the repo code.

---

## 0. Verdict on your plan, and what changed

**Keep as is:** the stage-gate discipline, compliance boundary (no peeking at instances/answers inside `solve`), the research-integrity stance, the evaluation rubric, trial logging with confidence intervals, the 6-page report outline, the 8-minute presentation timing and the viva checklist. These are the parts that protect your presentation and viva marks.

**Resolved by facts you flagged as open:**

| Your open question | Answer |
|---|---|
| Does one shot keep one character across gates? (F0) | **Yes [source]**: `simulator.py` draws one uniform character per shot and holds it for the whole circuit; every oracle gate becomes a phase `exp(2πi·(k·χ(h)+noise))` |
| Do mid-circuit `measure`/`reset`/`if_test` work, and are both branches charged? | **Yes [source/measured]**: iterative QPE with 144 measured bits in one circuit ran fine; `oracle.cost` matched my formula |
| Is 16-bit QPE enough? Need ~120 bits for continued fractions? | Not with this method: lattice reduction needs only **b≈36 bits per element [measured]** (§3) |
| "Package source not downloaded" | I read it all; the facts in this document come from the real code |

**Corrections to your plan:**
1. **Scoring weights conflict.** The website's Track 2 tab says auto **35 %**, open track **20 %**, presentation **45 %**; the GitHub README (your plan) says **60/20/20**. Ask the organisers; until then treat the presentation as worth a lot.
2. The "20k order / 40k logs / 12k verify / 8k reserve" split is unnecessary: the whole solve costs **≈14k** because order, logs and forgery come from one lattice.
3. The Step 5 candidates P1–P3 reduce to one design: **single-qubit iterative QPE on all four elements in one shot (P2+P3 combined)**; fixed-width QPE (P1) is not needed.
4. Your references R2, R4, R6–R9 I did **not** open; I confirmed R1 (Kitaev), R3 (Griffiths–Niu), R5 (Kimmel–Low–Yoder). The Bayesian-phase-estimation papers (R6–R8) are not needed for the chosen method, since the per-shot random character makes them a poor fit; cite them only as background after you have checked them.

---

## 1. Final technical decision

**Relation-lattice method.** Each shot gives one random point of the dual of H = ⟨A,B,C,g⟩. Measure the phases of all four elements to b bits inside one shot, collect t≈4 shots, and use LLL to recover the lattice

  Λ = { x ∈ ℤ⁴ : A^x0 · B^x1 · C^x2 · g^x3 = e }.

Λ has index |H| ≤ 2^60 in ℤ⁴, so generically it has a basis of norm ≈ 2^15, which is why 36 bits suffice. From the Hermite normal form (coordinates ordered A,B,C,g):
* `ord(g)` = last pivot,
* `log_g A, log_g B, log_g C` by back-substitution (or "not in ⟨g⟩" = forged),
* `secret = abc mod ord(g)` or `"none"`.

Cheap Hadamard tests then confirm each short relation (testing a relation is cheap; searching for one is what costs queries).

---

## 2. Hard constraints (merged)

| Item | Rule |
|---|---|
| Qubits | ≤ 16 per circuit (my solver uses **1**) |
| Gate cost | `max(1, bit_length(|k|))` per gate per shot; both `if_test` branches charged |
| Noise | Gaussian, σ = 0.002·√cost turns per oracle call |
| Budget / time | 80,000 queries / 600 s per hidden meeting |
| Imports | stdlib, numpy, qiskit (and `silentgroup` as in the template) |
| Interface | one file, `def solve(oracle) -> {"order": int, "secret": int | "none"}` |
| Auto score per meeting | 0 wrong order · 0.3 order only · 1 + 0.5·(budget left) both right |
| Forbidden | reading instance files/answers or judge internals in `solve` (development/testing only) |

---

## 3. Algorithm

```
repeat until verified (or budget low):
  SAMPLE   one circuit, shots=1: iterative QPE of φ_A, φ_B, φ_C, φ_g to b=36 bits each, 1 recycled qubit
  LATTICE  LLL on [I4 | X_i ; 0 | 2^b·I_t]; keep 4 short vectors whose tails are small; det ≤ 2^60
           (leave-one/two-out consensus tolerates 1–2 corrupted shots)
  HNF      -> order, logs or None
  VERIFY   Hadamard test (16 shots, reject if ≥3 ones) on each short relation;
           minimality test g^(n/p) != e for primes p ≤ 7 dividing n
return {"order": n, "secret": "none" if any log is None else (la*lb*lc) % n}
```

Per-circuit cost = 4·b(b+1)/2 = **2,664 queries** at b=36. Feedback is truncated to the 12 nearest bits (error < 2^-12 turn, below the noise floor). Reference code: `my_solver.py` (~240 lines).

Why one circuit per sample: the character changes between shots, so every bit that must refer to the same character has to be measured inside the same circuit, using feedback (Griffiths–Niu).

---

## 4. Evidence so far [measured, competition noise, organisers' simulator]

| Test | Result |
|---|---|
| Official `silentgroup.check`, public_1 / public_2, 6 repeats each | **12/12 correct**, ≈14.0k / 14.7k queries (one 31k run), 0.6–2 s, ≈1.41 pts each |
| Precision sweep, b=28/32/36/40 | order correct on public_1: 0/12, 9/12, **12/12**, 12/12 |
| Mock cyclic ≈2^58 · ℤ_{10^9+7}×ℤ_{10^9+9} · ℤ_{2^20}×ℤ_6×ℤ_{3^8} · ℤ_{2^55} | 10/10 each, median 13–15k queries |
| Mock ℤ_{2^29}² (rank-2 p-group) | 9/10 |
| Mock ℤ_{1000003} (tiny) | **8/10** |

The first version (b=40, "wait for a stable lattice") was also 10/10 but used 26k–79k queries: a reliability-for-cost lesson worth showing in the presentation.

**What this does not show:** performance on hidden groups of other shapes, noise levels above 0.002, or two bugs below.

---

## 5. Stage gates (your F0–F4, with status)

| Gate | Status | Remaining work |
|---|---|---|
| **F0** simulator truth table | Done by source reading; confirm with a 20-line test file (`tests/test_oracle_semantics.py`: one shot, same character, measure/reset/if_test) | write the test so you can demo it |
| **F1** small exact groups (ℤ₈, ℤ₁₅, ℤ₈×ℤ₃, ℤ₄×ℤ₈, ℤ₇×ℤ₁₁) | Not done | run solver with `noise=0` on these; assert exact answers |
| **F2** information scaling vs order bits | Partly (precision sweep, shapes) | sweep log₂|G| ∈ {20,30,40,50,58}; record b needed, samples, queries |
| **F3** logs and forgery on held-out groups | Partly (public_2, forged mocks) | 20 seeds × 8 shapes incl. 1, 2 and 3 forgeries; log false/missed forgery |
| **F4** full budget/time/official checker | Done for public instances | repeat ≥ 20 times after any code change |

---

## 6. Hardening backlog (do in this order)

1. **P0 — stall guard.** If no lattice is ever found, the leave-one/two-out search grows ~t² LLL calls and can exceed 600 s. Cap it to the 10 most recent samples, add a ~400 s wall-clock guard, and return the best-effort answer. *Not implemented in the file I gave you.*
2. **P0 — tiny groups (2/10 wrong).** Seed 1: wrong order 2048006144; seed 4: right order but secret `"none"`. Cause unknown. [hypotheses] (a) many near-relations in tiny groups make `good[:4]` pick a wrong or dependent set; (b) tolerance `2^(b−12)` badly tuned; (c) verification threshold lets a false relation pass. Debug by printing HNF diagonals and determinant per sample count. This debug run was blocked by the session's safety classifier earlier, so it is un-run.
3. **P1 — rank-2 / non-cyclic p-groups (1/10 wrong).** Require ≥6 samples when the HNF quotient is non-cyclic.
4. **P1 — skewed lattices** (very small logs make a ~2^56-norm basis vector that b=36 cannot resolve; order would come out as a multiple of the truth). Precision ladder 36→44→52 when minimality or verification fails; factor n (Pollard rho is instant below 2^60) and test `g^(n/p) ≠ e` for all primes within a cost cap.
5. **P2 — never answer `"none"` on timeout** (your rule is right): only output forgery when the lattice positively says the log does not exist and verification passed.

Keep the three internal states from your Step 8 (`VALID / FORGED / UNRESOLVED`); in this method they correspond to "log found and verified", "HNF says no log and lattice verified", "lattice not yet verified".

---

## 7. Budget math (replaces the 20/40/12/8 split)

| Item | Queries |
|---|---|
| 4 samples × 2,664 | 10,656 |
| verify 4 short relations (16 shots × ~60) | ≈ 3,800 |
| small-prime minimality | 0–3,000 |
| **Total [measured]** | **≈ 14k → 82 % left → ≈ 1.41 pts** |

Do not optimise queries further; the maximum extra is ≈ 0.04 pts. Reliability is worth far more.

---

## 8. Critical path tonight (deadline today; cut-off hour not stated — check the form)

- [ ] **P0** Clone, install, run `python -m silentgroup.check my_solver.py --repeat 5`; confirm it matches §4
- [ ] **P0** Submit that file via the form **now** as a safety net
- [ ] **P0** Read and re-derive the solver until each of you can explain it (§10 table); delete or rewrite anything you cannot explain
- [ ] **P0** Stall guard (§6.1), re-run checker ×20
- [ ] **P0** Tiny-group bug (§6.2) with a mock regression: 8 shapes × 20 seeds, log failures
- [ ] **P1** Rank-2 fix and precision ladder (§6.3–6.4); re-run regression
- [ ] **P1** Freeze one file; audit imports; remove prints; re-run checker on the frozen file; **re-submit**
- [ ] **P1** Open-track experiment and plots (§9)
- [ ] **P2** 6-page report and 8-minute deck (§10)

Stop conditions (yours, kept): a distribution is falsified by exact small tests; the simulator lacks a required operation; held-out failures show overfitting. In that case revert to the last frozen version that passed the checker.

---

## 9. Open track (20 %)

Merge of your title with my experiment: **"Noise–precision–query laws for relation-lattice recovery in a noisy black-box abelian group."**
1. *Cost law:* per sample 2b(b+1); samples ≈ 4 (cyclic dual), more for higher rank; required precision **b ≈ 16 + log₂|H| / t** (my heuristic count: ~2^60 candidate vectors each passing t tolerance tests with probability ~2^{16−b}). Compare predicted vs measured over log₂|G| ∈ {20,30,40,50,58}.
2. *Noise sweep:* `Oracle(..., noise=ν)` for ν ∈ {0.002, 0.004, 0.008, 0.02, 0.05}; plot success and queries. **[hypothesis]** corrupted-shot rate ∝ ν² (a stage flips with probability ≈ (πσ)², and stages with weight ≳ 2^-19 break a norm-2^15 relation), predicted collapse near ν ≈ 0.02–0.05.
3. *Mitigation:* repeat the cheap high-weight stages (m ≤ 19) twice and vote; report whether it helps per query.
4. *Bonus:* output ℤ_{d1}×ℤ_{d2}×… via Smith normal form of Λ (same lattice, no extra queries).
Report negative results (b=28/32 failures, stability-rule overspend); separate theorem / simulator fact / empirical observation, as your plan says.

---

## 10. Report, presentation and viva

**6-page report (your outline, filled in):** p1 problem and the per-shot-character insight; p2 character math and why the lattice has short vectors; p3 HNF → order/logs/forgery and the verification tests; p4 circuit, cost formula, noise, compliance; p5 measurements with confidence intervals (§4, §5); p6 open-track result, limitations, references.

**8-minute talk (your timings):** 0–1 problem · 1–2:15 character model · 2:15–3:45 circuit + lattice · 3:45–5 membership and noise · 5–6:30 results and demo (`check` output) · 6:30–8 open track, **what failed**, limits.

**Function map (be able to modify any of these):**

| Function | Job |
|---|---|
| `build_pe_circuit` | 1-qubit iterative QPE, four elements, feedback phases |
| `decode` | classical bits → integers ≈ 2^b·φ |
| `lll` | exact-rational LLL |
| `relation_lattice`, `consensus` | embed samples, reduce, filter small tails, HNF, drop corrupted shots |
| `hnf`, `extract` | order, logs, forgery by back-substitution |
| `relation_fails`, `verify` | Hadamard tests of relations |

**Viva questions to rehearse:** why a fresh character per shot; why `U_h|χ̂⟩ = e^{2πiχ(h)}|χ̂⟩`; why one qubit suffices; what the index of Λ is; why 2^60 gives a determinant check; what happens when A ∉ ⟨g⟩; why failing to find a log is not proof of forgery; why testing a relation is cheaper than finding it; false-relation pass probability of 16 shots (≈0.2 %); what breaks first as noise grows; modular equations with non-coprime coefficients.

---

## 11. Compliance checklist (from your plan)
- [ ] Only `oracle.run/cost/budget_left/elements`; no instance files, no public answers, no judge introspection
- [ ] ≤ 16 qubits; no loops in circuits; no unbound parameters
- [ ] No debug prints, no hard-coded public answers
- [ ] Returns exact types: `int` order, `int` or `"none"` secret
- [ ] Time < 600 s with a guard; budget never exceeded
- [ ] Frozen file re-checked with `--repeat 5` or more

---

## 12. References
**Checked by me:** Kitaev quant-ph/9511026 · Griffiths–Niu PRL 76, 3228 (arXiv quant-ph/9511007) · Mosca–Ekert quant-ph/9903071 (one control qubit) · Cheung–Mosca cs/0101004 · Seifert ePrint 2000/036 · Ekerå–Håstad ePrint 2017/077 · Ekerå J. Math. Cryptol. 15(1) (ePrint 2018/797) · Ekerå arXiv:1905.09084 · Kimmel–Low–Yoder arXiv:1502.02677 (note their erratum) · arXiv:2609.36480 (abstract only) · Bermejo-Vega–Lin–Van den Nest arXiv:1409.4800.
**From your plan, not checked by me:** Cleve et al. quant-ph/9708016 · IET Inf. Security 2023 survey · Smith et al. PRA 109, 042412 · de Neeve et al. arXiv:2405.08930 · Hurant et al. arXiv:2407.18339 · IBM dynamic-circuits docs.
The precision rule b ≳ 16 + 60/t and the ν² corruption scaling are **my derivations**, not paper claims.

---

## 13. Honest status
The prototype passes the public instances and most mock shapes, but it has known failures (tiny groups, one p-group case, stall risk), has never seen a hidden instance, and the guard code is not yet written. AI help is allowed for the build; anything you cannot explain or modify in the presentation counts as not submitted.
