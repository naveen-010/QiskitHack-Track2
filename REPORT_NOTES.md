# REPORT_NOTES — The Silent Group (QFF 2026, Track 2)

Material for the 6-page report and the 8-minute talk. Written so that a team member who did not write
the code can explain and change it. Labels: **[theorem]** = mathematically proven; **[simulator fact]**
= read in `silentgroup/` source and confirmed by `tests/test_oracle_semantics.py`; **[empirical]** =
measured with the harnesses in `engineering/` (numbers, seeds and CSVs in `results/`).

Frozen file: `submission/my_solver.py` (≈770 lines, imports: `math, random, time, fractions`,
`qiskit.QuantumCircuit`, `silentgroup.oracle_gate`). Baseline kept untouched: `my_solver.py`.

---

## 1. The 6-page report outline (filled in)

### Page 1 — Problem and the key insight
* Task: the hidden group G ≅ ℤ_{d1}×…, |G| ≤ 2^60, is only reachable through controlled `oracle_gate(h, k)`
  (multiply by h^k). Output ord(g) and abc mod ord(g) where A=g^a, B=g^b, C=g^c, or `"none"` if one of
  A, B, C ∉ ⟨g⟩. Score per meeting: 0 if order wrong, 0.3 order only, 1 + 0.5·(budget left) if both right.
* **[simulator fact]** The group register starts in |e⟩ = |G|^{-1/2} Σ_χ |χ̂⟩ and is never measured, so
  each shot samples ONE uniformly random character χ and every oracle gate becomes a phase
  e^{2πi(k·χ(h) + noise)} on its control qubit (`simulator.py:79-82`). The character is fixed for the
  whole circuit and fresh for the next shot (test `test_one_character_per_shot`: two Hadamard tests in
  one circuit agree 400/400 in ℤ₂; across shots they disagree ~50 %).
* Consequence: everything that must refer to the same character has to be measured **inside one
  circuit**. We measure the phases of all four elements in one shot.

### Page 2 — Character maths: why a lattice, why it has short vectors
* χ(h) ∈ ℚ/ℤ. For x ∈ ℤ⁴, A^{x0}B^{x1}C^{x2}g^{x3} = e ⇔ x·(χ(A),χ(B),χ(C),χ(g)) ≡ 0 (mod 1) for every χ
  **[theorem: characters separate points of a finite abelian group]**.
* The relation lattice Λ = {x : A^{x0}B^{x1}C^{x2}g^{x3} = e} has ℤ⁴/Λ ≅ H = ⟨A,B,C,g⟩, so
  det Λ = |H| ≤ |G| ≤ 2^60 **[theorem]**. A "random" lattice of determinant 2^60 in dimension 4 has a
  basis of norm ≈ 2^15 (Minkowski/Gaussian heuristic) — this is why b ≈ 36 phase bits suffice instead of
  the ≈120 bits continued fractions would need.
* Each shot gives a point φ of the dual lattice Λ* (to b bits). With t shots, LLL on
  `[I₄ | X ; 0 | 2^b I_t]` (X_i ≈ 2^b φ_i) finds short x whose "tails" x·X_i mod 2^b are small.

### Page 3 — From the lattice to the answer, and the two-sided check
* HNF of Λ (upper triangular, coordinates A,B,C,g): **ord(g) = last pivot**; back-substitution gives
  log_g A etc., or proves "no relation A·g^y = e exists" ⇒ A ∉ ⟨g⟩ ⇒ forged **[theorem]**.
* Testing a relation is cheap, finding one is expensive: a Hadamard test on h = A^{x0}…g^{x3} gives
  P(1) = (1 − Re χ(h))/2, averaged over χ: **0 if h = e, exactly 1/2 otherwise [theorem]**.
* **Check 1, L′ ⊆ Λ** (`verify_generators`): every LLL basis vector of the recovered L′ is tested.
* **Check 2, Λ ⊆ L′** (`complete`) — new in this version. If L′ ⊊ Λ then Λ/L′ is a non-trivial finite
  group, so it contains x of prime order p, p | [Λ:L′] | det L′, and x is (up to L′) on one of the
  (p^r − 1)/(p − 1) lines of the p-torsion of ℤ⁴/L′, computed as the left kernel of the basis mod p
  **[theorem]**. We test one vector per line; a relation found is added and the step repeats.
  A prime that is cleared stays cleared (if [Λ:L′] is coprime to p, so is [Λ:L″] for L′ ⊆ L″ ⊆ Λ).
* Why Check 2 matters **[empirical]**: one noise-flipped QPE stage (phase error exactly ±2^j) makes L′ a
  proper SUBLATTICE of true relations. In a group of order 2^60 that sublattice has det > 2^60 and is
  discarded; in a group of order 10^6 it is accepted. Baseline: 80–90 % on tiny groups, false `"none"`.
* `"none"` is output only after both checks pass ("VALID / FORGED / UNRESOLVED": unresolved ⇒ best-effort
  integer secret, never `"none"`).

### Page 4 — Circuit, cost, noise, compliance
* One qubit, recycled with `reset`. Iterative QPE (Kitaev; Griffiths–Niu semiclassical QFT): for
  m = b…1 (LSB first) apply H, controlled h^{2^{m−1}}, phase corrections −2π/2^{l−m+1} conditioned on the
  12 nearest already-measured bits, H, measure, reset. Truncating feedback at 12 bits leaves an error
  < 2^{-12} turn per stage — below the noise.
* **[simulator fact]** cost of `oracle_gate(h,k)` = max(1, bit_length|k|) per shot; both if/else branches
  are charged (`program.py:166-176`, test `test_cost_formula_both_branches`). Our feedback ifs contain no
  oracle gates, so they are free. One 4-element sample at b bits costs 4·b(b+1)/2 = 2,664 at b=36.
* **[simulator fact]** noise: each oracle call's phase gets N(0, (0.002·√cost)²) turns.
* Compliance: only `oracle.run / cost / budget_left` are used (`oracle.elements` is not needed: names are
  fixed A,B,C,g); no file access, no prints, no instance/answer access; ≤ 16 qubits (we use 1); no
  loops in circuits (Python unrolls them); all parameters bound. Guards: 400 s wall clock, budget check
  before every `run`, LLL-call cap 400, completion-step cap 64.

### Page 5 — Measurements (with 95 % Wilson intervals)
See §3 below (tables copied from `results/`).

### Page 6 — Open track, limitations, references
See §4 (open-track results), §5 (what failed), §7 (limitations). References: plan §12.

---

## 2. Algorithm in one screen (for slides)

```
for b in 36, 44, 52:                         # precision ladder
  repeat (≤ 7/8/8 samples):
    SAMPLE  one shot of the 1-qubit QPE circuit -> X_A, X_B, X_C, X_g  (b bits each)
    LATTICE LLL on [I4 | X; 0 | 2^b I]  (all samples, then leave-1-out, then leave-2-out)
            keep first 4 independent rows with small tails; HNF; det <= 2^60
    VOTE    prefer lattices also found at the two previous sample counts
            non-cyclic quotient: need >= 6 samples and 2 votes
    CHECK1  Hadamard-test each short basis vector (curtailed, error <= 1e-4 both ways)
    CHECK2  complete(): p-torsion lines for every prime p | det, add missing relations
    return order = last pivot, secret = abc mod n  or "none"
  skew mode (3 short relations only): 80-bit QPE of ONE element e -> ord(e); L' = M + ord(e) e_e; CHECK1+2
  downward repair (data keeps producing one rejected lattice): replace the bad vector y by k*y; CHECK2
best effort on budget/time exhaustion (never "none")
```

---

## 3. Measurements

### 3.1 Official checker (judge scoring, unseeded OS randomness)
| file | instance | runs | correct | mean queries | mean pts |
|---|---|---|---|---|---|
| baseline `my_solver.py` | public_1 | 5 | 5 | 13,997 | 1.413 |
| baseline | public_2 | 5 | 5 | 21,358 (14.7k–37.1k) | 1.366 |
| **frozen** `submission/my_solver.py` | public_1 | 20 | 20 | 17,025 (max 19.5k) | **1.393** |
| **frozen** | public_2 | 20 | 20 | 40,161 (33k–76k) | **1.249** |

[empirical] The frozen file costs more queries (≈ −0.02 pts on public_1, ≈ −0.12 on public_2) in exchange for
the reliability below. Max wall time 3.5 s (limit 600 s).

### 3.2 Mock regression, baseline (20 seeds) — `results/regress_baseline.txt`
```
shape                    n  ok      wilson95 ordOK falseF missF exc    medQ    maxQ   maxT
cyclic_2^58             20  20 [0.839,1.000]    20      0     0   0   14120   28360   11.2 FAIL
p1e9_x_p1e9             20  20 [0.839,1.000]    20      0     0   0   14224   20248    1.6 FAIL
2^20x6x3^8              20  20 [0.839,1.000]    20      0     0   0   13252   16576    1.4 FAIL
cyclic_2^55             20  20 [0.839,1.000]    20      0     0   0   14699   21126    1.7 FAIL
2^29_sq                 20  19 [0.764,0.991]    20      1     0   0   12663   31238    2.3 FAIL
tiny_1000003            20  17 [0.640,0.948]    19      3     0   0   11704   36704    2.7 FAIL
2^30x3^18               20  20 [0.839,1.000]    20      0     0   0   15685   24366    2.3 FAIL
tiny_7                  20  16 [0.584,0.919]    19      3     0   0   10886   75294    5.6 FAIL
tiny_64                 20  18 [0.699,0.972]    20      2     0   0   11037   16374    1.6 FAIL
tiny_1001               20  16 [0.584,0.919]    20      4     0   0   11312   12128    1.3 FAIL
tiny_65537              20  18 [0.699,0.972]    20      2     0   0   11504   12176    1.3 FAIL
forged1_C               20  17 [0.640,0.948]    17      0     2   0   18835   75160 1386.1 FAIL
forged2_AB              20  20 [0.839,1.000]    20      0     0   0   15046   36174    2.8 FAIL
forged3_ABC             20  20 [0.839,1.000]    20      0     0   0   14932   34680    3.0 FAIL
forged_small_12x6       20  20 [0.839,1.000]    20      0     0   0   10988   19188    1.7 FAIL
public1_shape           20  20 [0.839,1.000]    20      0     0   0   14888   21024    2.1 FAIL
public2_shape           20  17 [0.640,0.948]    20      3     0   0   14435   20112    1.8 FAIL
public2_shape_forgA     20  16 [0.584,0.919]    16      0     3   1   20159   35788    inf FAIL
skew_2^55_a2b3c5        20   0 [0.000,0.161]     0      5     0   9   12568   75312    inf FAIL
skew_2^58_a1b1c1        20   0 [0.000,0.161]     0      0     0  18      -1   27070    inf FAIL
skew_p1e9sq_small       20   0 [0.000,0.161]     0      0     0  20      -1      -1    inf FAIL
```
(The PASS/FAIL column applies the WP4 rule; with 20 seeds even 20/20 has Wilson lower bound 0.839.)

### 3.3 Mock regression, frozen file (40 seeds) — `results/regress_final_40seeds.txt`
```
shape                    n  ok      wilson95 ordOK falseF missF exc    medQ    maxQ   maxT
cyclic_2^58             40  40 [0.912,1.000]    40      0     0   0   15962   21248    3.8 PASS
p1e9_x_p1e9             40  40 [0.912,1.000]    40      0     0   0   16540   21496    2.6 PASS
2^20x6x3^8              40  40 [0.912,1.000]    40      0     0   0   13522   31515    2.1 PASS
cyclic_2^55             40  40 [0.912,1.000]    40      0     0   0   15754   30659    3.6 PASS
2^29_sq                 40  40 [0.912,1.000]    40      0     0   0   13043   39917    2.4 PASS
tiny_1000003            40  40 [0.912,1.000]    40      0     0   0   12009   21633    1.4 PASS
2^30x3^18               40  40 [0.912,1.000]    40      0     0   0   16308   52550    8.0 PASS
tiny_7                  40  40 [0.912,1.000]    40      0     0   0   10898   15386    1.4 PASS
tiny_64                 40  40 [0.912,1.000]    40      0     0   0   11065   34908    3.2 PASS
tiny_1001               40  40 [0.912,1.000]    40      0     0   0   11330   18236    1.7 PASS
tiny_65537              40  40 [0.912,1.000]    40      0     0   0   11734   18855    1.4 PASS
forged1_C               40  40 [0.912,1.000]    40      0     0   0   23724   63681    6.5 PASS
forged2_AB              40  40 [0.912,1.000]    40      0     0   0   23758   55308   10.2 PASS
forged3_ABC             40  40 [0.912,1.000]    40      0     0   0   22519   56573    9.7 PASS
forged_small_12x6       40  40 [0.912,1.000]    40      0     0   0   16335   22962    2.0 PASS
public1_shape           40  40 [0.912,1.000]    40      0     0   0   17060   22474    3.7 PASS
public2_shape           40  40 [0.912,1.000]    40      0     0   0   14416   39087    2.0 PASS
public2_shape_forgA     40  40 [0.912,1.000]    40      0     0   0   38069   67810    7.4 PASS
skew_2^55_a2b3c5        40  40 [0.912,1.000]    40      0     0   0   36567   64410    6.1 PASS
skew_2^58_a1b1c1        40  40 [0.912,1.000]    40      0     0   0   33520   54291    4.2 PASS
skew_p1e9sq_small       40  40 [0.912,1.000]    40      0     0   0   37206   70055    7.0 PASS
```
[empirical] 840/840; Wilson 95 % lower bound 0.912 on every shape; 0 false and 0 missed forgeries;
0 exceptions; max 70,055 queries, max 10.2 s. Note: the WP4 rule "lower bound ≥ 0.90 with 30 seeds" is
unreachable at n = 30 (30/30 gives 0.886), so 40 seeds were used.

### 3.4 Unit and semantics tests — `python -m pytest tests -q` → 96 passed
* F0 (`test_oracle_semantics.py`): one character per shot; measure→reset→if_test feedback; cost formula with
  both branches charged; identity vs non-identity P(1); σ ∝ √cost.
* F1 (`test_small_exact.py`): noise = 0, exact answers on ℤ8, ℤ15, ℤ8×ℤ3, ℤ4×ℤ8, ℤ7×ℤ11, ℤ7, ℤ64, ℤ1001,
  ℤ65537 and forged ℤ12×ℤ6 (C), ℤ4×ℤ8 (AB), ℤ6×ℤ6×ℤ10 (ABC), ℤ2×ℤ4×ℤ8 (B), ℤ7×ℤ7 (A); 5 seeds each.
* `test_lattice.py`: LLL output has the same |det| and lattice, satisfies Lovász and size reduction;
  HNF is triangular with positive pivots and the same lattice; `extract` on hand-built lattices incl. None.


---

## 4. Open track: noise–precision–query laws

Title: **Noise–precision–query laws for relation-lattice recovery in a noisy black-box abelian group.**
Code: `engineering/experiments.py` (`e1`…`e4`, `plots`); data `results/e*.csv`; figures `results/e*.png`.
Seeds ≥ 20 per point, Wilson 95 % intervals.

### E1 — Cost law (precision needed at t = 4 samples) — `e1_cost_law.png`, `e1_cost_law.csv`
Cyclic mock groups of order ≈ 2^k (random g, A, B, C), 30 seeds per (k, b), success = the 4-sample LLL
lattice equals the true Λ exactly (no consensus, no tests).

| log₂|H| | b₅₀ measured | heuristic 16 + log₂|H|/4 | plateau (top b) | b for ≥ 95 % |
|---|---|---|---|---|---|
| 20.0 | 17.3 | 21.0 | 0.96 [0.91, 0.98] | 20 |
| 29.9 | 20.1 | 23.5 | 0.96 [0.91, 0.98] | 26 |
| 39.9 | 23.1 | 26.0 | 0.86 [0.79, 0.91] | 32 |
| 50.0 | 28.0 | 28.5 | 0.80 [0.72, 0.86] | not reached |
| 57.8 | 31.6 | 30.4 | 0.78 [0.70, 0.85] | not reached |

* **[empirical]** Fit: b₅₀ ≈ 9.0 + 0.38·log₂|H|. The heuristic has the right magnitude near |H| ≈ 2^50–2^58
  but the wrong slope (0.25 vs 0.38): it is too pessimistic for small groups and slightly optimistic at 2^58.
  Partially **negative** result for the heuristic.
* **[empirical]** Above the onset, success does not reach 1: it plateaus at 0.78–0.96, decreasing with |H|.
  The plateau is set by corrupted shots (E2), not by precision — raising b from 36 to 42 does not help.
  This is why the solver needs consensus (leave-one/two-out) and why "≥ 95 % at 4 samples" is unreachable
  for |H| ≥ 2^50. Queries for the 4 samples at b₅₀: 2.5k (2^20) → 8.3k (2^58); predicted 3.7k → 7.7k.
* **[theorem]** per-sample cost 4·b(b+1)/2 (sum of bit lengths of 2^{m−1}, m = 1…b, four elements).

### E2 — Noise sweep — `e2_noise.png`, `e2_corruption.png`
Full frozen solver, 20 seeds per point, budget 80k:

| ν | cyclic 2^58 ok | public1-shape ok | median queries | false "none" |
|---|---|---|---|---|
| 0.002 (competition) | 20/20 | 20/20 | 16.1k / 17.3k | 0 |
| 0.004 | 20/20 | 20/20 | 18.5k / 17.3k | 0 |
| 0.008 | 5/20 [0.11, 0.47] | 9/20 [0.26, 0.66] | 66.9k (budget guard) | 0 |
| 0.02 | 0/20 | 0/20 | 66.9k | 0 |
| 0.05 | 0/20 | 0/20 | 66.9k | 0 |

Corrupted-shot rate at b = 36 (|G| ≈ 2^40, 300 shots per point; a shot is corrupted if a TRUE short relation
has a tail > 2^{b−12}): 0.033 [0.018, 0.060], 0.14, 0.46, 0.98, 1.00.
* **[empirical]** log–log slope on the unsaturated points = **1.90**, consistent with the ν² hypothesis
  (a stage flips with probability ≈ (πσ)², σ = ν√m). Hypothesis **supported** below saturation.
* **[empirical]** The plan predicted collapse near ν ≈ 0.02–0.05; it actually collapses at **ν ≈ 0.008**
  (4× competition noise), where ≈ 46 % of shots are corrupted. Prediction **falsified** (too optimistic).
* **[empirical]** Failure mode under high noise is graceful: no false forgeries and no exceptions in 200 runs;
  the solver stops at the budget guard (66.9k) and returns a best-effort answer.

### E3 — Mitigation: 3-fold majority vote on the high-weight stages m ≤ 19 — `e3_mitigation.png`
Stages m ≤ 19 are measured 3 times; feedback uses the majority (register `if_test` on values 3, 5, 6, 7 —
free, no oracle gates). Cost 4,184 vs 2,664 queries per shot (+57 %). 300 shots / 40 seeds per point.

| ν | corrupted shots plain → voted | 4-sample recovery plain → voted | recoveries per 10k queries |
|---|---|---|---|
| 0.002 | 3.3 % → 0.33 % | 0.83 → 0.85 | 0.77 → **0.51** |
| 0.008 | 46 % → 13 % | 0.00 → 0.10 | 0 → 0.06 |
| 0.02 | 98 % → 68 % | 0 → 0 | 0 → 0 |

* **[empirical] Negative result at competition noise:** voting removes 90 % of corrupted shots but does not
  improve 4-sample recovery (the remaining failures are precision/“false short vector” failures at |G| = 2^58),
  so it is **worse per query**. Not adopted in the solver.
* **[empirical]** At 4× noise it is the only variant that recovers at all (10 %) — a mitigation for a noisier
  oracle, not for this one.

### E4 — Group structure from the Smith normal form (no extra queries) — `e4_structure.csv`
`group_structure(H)` (in `experiments.py`, outside the submission) computes invariant factors
s_k = D_k / D_{k−1} from determinantal divisors D_k (gcd of k×k minors) **[theorem]** of the verified lattice
captured from the solver. 10 seeds × 6 shapes: structure equal to the truth **60/60**; the largest invariant
equals lcm(ord A, ord B, ord C, ord g) computed independently from the group **60/60**. Examples:
cyclic ℤ_288230376151711717; forged2: ℤ_100000000 × ℤ_300000000; forged3: ℤ_134217728 × ℤ_805306368.


---

## 5. What failed (negative results — show these)

| Attempt | What happened | Lesson |
|---|---|---|
| b = 28 / 32 (plan §4) | order correct on public_1 in 0/12 and 9/12 runs | precision must exceed ≈16 + log₂|H|/4; see E1 |
| "wait for a stable lattice" at b=40 (first prototype) | 10/10 correct but 26k–79k queries | stability rules overspend; vote only where needed (non-cyclic) |
| Shor-style continued fractions | needs ≈120 phase bits per element (cost ≈ 4·120·121/2 ≈ 29k per shot) | the lattice uses ≈36 bits because Λ has short vectors |
| baseline verification (L′ ⊆ Λ only) | tiny groups 80–90 %, false `"none"` (14 runs in 120), one 1386 s stall | verification must be two-sided: added `complete` |
| `hnf(all good LLL rows)` | false rows with head ≈ tail ≈ 2^21 pass the tail test → public_1 cost 44k | use only the first 4 linearly independent good rows |
| "same lattice as 2 samples ago" (literal) | HNF was not canonical → never equal; after canonical form, a corrupted 7th shot still broke it | canonical HNF + a vote over candidate sets |
| first p-torsion enumeration | looped `range(p)` for p ≈ 2^50 → hang | enumerate line representatives directly |
| fixed "≥2 ones in 20 ⇒ not a relation" in completion | a true relation of cost 220 has P(1) ≈ 0.9 % → falsely rejected, wrong order 2^58 | noise-aware (n, k) from the known σ, both error rates ≤ 1e-4 |
| order finding from 4 samples at 80 bits | if ord(g) = 2^54 phases carry only 54 bits → spurious short h | tighter tolerance, leave-one/two-out, Hadamard-test each candidate |
| skew mode using only g | forged A with order 2^25.9 modulo ⟨B,C,g⟩: M + n·e_g is rank-deficient | use the first element whose unit vector is independent of M |
| mitigation by 3-fold voting (E3) | corrupted shots 3.3 % → 0.33 % but recovery per 10k queries 0.77 → 0.51 at ν = 0.002 | not adopted; helps only at 4× noise |
| heuristic b ≈ 16 + log₂|H|/4 (E1) | slope 0.38 measured, not 0.25 | use the measured law |
| predicted noise collapse at ν ≈ 0.02–0.05 (E2) | collapse at ν ≈ 0.008 | |
| "wait for 2 more samples" for skew mode (t ≥ 6) | skew shapes fell to 36/40 (budget) | reverted |

---

## 6. Function-by-function map (`submission/my_solver.py`)

| Function | Job | Change it if… |
|---|---|---|
| `_Ctx` | per-solve state: deadline (400 s), LLL counter, cache of tested vectors; `run()` refuses any run that would exceed the budget | you change the guards |
| `build_pe_circuit(b)` | 1-qubit iterative QPE of A,B,C,g, b bits each, 12-bit truncated feedback | you change precision/feedback |
| `build_single_circuit(b, e)` | same for one element (skew mode, b=80) | |
| `decode(key, b)` | bitstring → integers X_e ≈ 2^b·χ(e) | you change the clbit layout |
| `lll` | exact-rational LLL, δ = 3/4 | never needed |
| `hnf` | canonical row HNF (0 ≤ H[i][j] < H[j][j]), accepts ≥ 4 rows | |
| `short_relations` | builds `[I | X; 0 | 2^b I]`, LLL, keeps rows with tails ≤ 2^(b−12) | tolerance |
| `independent_prefix` | first k linearly independent rows in LLL order | |
| `relation_lattice` | 4 independent good rows → (HNF, short basis); rejects det > 2^60 | |
| `candidates` | all samples, leave-one-out, leave-two-out (window of 10) | consensus depth |
| `extract` | HNF → ord(g), logs or None | coordinate order |
| `is_prime`, `pollard_rho`, `prime_factors` | factor det ≤ 2^90 (deterministic MR, Brent rho) | |
| `left_kernel_mod_p`, `p_rank`, `torsion_candidates`, `reduce_mod_lattice` | p-torsion of ℤ⁴/L′, one short vector per line (≤ 40 lines, else 40 random) | LINE_CAP |
| `hadamard_circuit` | H – oracle gates – H – measure | |
| `test_plan`, `relation_test` | smallest (n, k) with false accept/reject ≤ 1e-4 for this cost; curtailed (stops when decided) | TEST_ERR |
| `in_lattice`, `verify_generators` | Check 1 with cache; free rejection if a known non-relation lies in L′ | |
| `downward_repair` | rejected lattice that keeps coming back: one bad y → k·y | |
| `complete` | Check 2 (missing relations), cleared primes kept | |
| `is_noncyclic` | some p with p-rank ≥ 2 | |
| `order_of`, `skew_lattice`, `small_span`, `try_skew` | skew mode | SKEW_B |
| `answer` | dict; `"none"` only if trusted | |
| `solve` | ladder / sampling loop / vote / checks / fallbacks | |

---

## 7. Known limitations (honest)

* **Not a proof of correctness on hidden instances.** All reliability numbers are on `MockGroup` shapes we chose (21 shapes × 40 seeds) and the two public instances. Hidden instances may have shapes we did not imagine.
* **Cost on non-cyclic / forged groups.** public_2: ≈34–39k queries typically (1.28 pts) and ≈60–75k in roughly 10–20 % of runs (skew mode triggered by a mildly skewed true lattice — B forged with order 585585 — still correct). Baseline averaged ≈21k on public_2 but had 15 % false forgeries on public_2-shaped mocks.
* **Completion is exhaustive only up to 40 lines per prime.** For p-rank ≥ 2 and p ≥ 41, 40 random p-torsion elements are tested instead (not a proof). Corruption-induced sublattices observed so far had 2-power index.
* **Test errors are ≤ 1e-4 per test assuming the competition noise σ = 0.002·√cost.** The constant `NOISE = 0.002` is hard-coded (the oracle's noise attribute is not in the allowed interface). At higher noise the false-reject rate of true relations rises (see E2).
* **Skew mode needs a 3-volume test and 2 hints;** cases with two long directions (rank-2 short part) rely on the precision ladder (b ≤ 52), not on skew mode.
* **Best-effort answers on exhaustion** get the order right in all observed cases, but the secret is 0 (0.3 pts) rather than a guess.
* **Seeded-RNG artifact:** forged2_AB seed 13 needed the downward repair because 15 consecutive shots had χ(h)=0 for one order-2 element (P≈2^-15 under fair sampling).


---

## 8. The 8-minute talk script

| Time | Slide | What to say (and show) |
|---|---|---|
| 0:00–1:00 | Problem | Diffie–Hellman in a group we never see. Output ord(g) and abc mod ord(g), or "none" if a broadcast is forged. Score: correctness first, 0.5 bonus for budget left. |
| 1:00–2:15 | Character model | Show `simulator.py` lines 79–82: one random character per shot, fixed for the circuit. Demo `tests/test_oracle_semantics.py::test_one_character_per_shot` (400/400 agree inside a shot, ~50 % across shots). "So all four phases must come from ONE shot." |
| 2:15–3:45 | Circuit + lattice | 1 qubit, iterative QPE, 4 × 36 bits, 2,664 queries. Each shot = a point of the dual lattice; LLL on [I | X; 0 | 2^b I] gives the relation lattice Λ, det = |H| ≤ 2^60 → short vectors ≈ 2^15 → 36 bits suffice (vs ≈120 for continued fractions). HNF → ord(g) = last pivot, logs by back-substitution. |
| 3:45–5:00 | Membership + noise | Hadamard test: P(1)=0 for a relation, exactly ½ otherwise. Two-sided check: generators are relations (L′ ⊆ Λ) AND no p-torsion element of ℤ⁴/L′ is a relation (Λ ⊆ L′). Show the bug it fixed: one flipped QPE stage (error −2^25) gave a sublattice of index 2^11 → baseline answered order 2048·1000003. |
| 5:00–6:30 | Results + demo | Run `python -m silentgroup.check submission/my_solver.py --repeat 3` live (≈10 s). Table: 840/840 on 21 mock shapes × 40 seeds, 0 false forgeries; baseline 300/420 with stalls. Cost: ≈16.6k on public_1 (1.40 pts), ≈35k typical on public_2 (non-cyclic, forged). |
| 6:30–8:00 | Open track, what failed, limits | E1 law b50 ≈ 9 + 0.38·log₂|H| and the corrupted-shot plateau; E2 noise sweep and ν² test; E3 voting is a negative result (see §4). What failed: b=28/32, stability rule overspend, hnf(all good rows), fixed test thresholds. Limits: §7. |


---

## 9. Viva questions (15) with short answers

1. **Why does each shot see a fresh character?** The register starts in the uniform superposition of character states and is never measured, so different χ never interfere; the simulator samples χ per shot (`simulator.py:79`).
2. **Why is U_h|χ̂⟩ = e^{2πiχ(h)}|χ̂⟩?** Character states are the Fourier basis of the regular representation; multiplication by h is diagonal in it with eigenvalue χ(h).
3. **Why does one qubit suffice?** Iterative QPE measures one bit at a time and feeds it back classically (Griffiths–Niu); `reset` recycles the qubit; ≤16 qubits is never an issue.
4. **What is the index of Λ in ℤ⁴?** |H| where H = ⟨A,B,C,g⟩, because ℤ⁴/Λ ≅ H. Hence det Λ ≤ |G| ≤ 2^60.
5. **Why does 2^60 give a determinant check?** A lattice computed from a corrupted shot (or with a false vector) usually has det ≠ |H|; det > 2^60 is impossible for the true Λ. It is one-sided: a sublattice with det ≤ 2^60 passes — that is the tiny-group bug, fixed by `complete`.
6. **What happens when A ∉ ⟨g⟩?** No vector (1, 0, 0, y) is in Λ; in the HNF, back-substitution for e_A hits a non-divisible pivot, so `extract` returns None for log A and we answer "none" — only after both checks passed.
7. **Why is failing to find a log not proof of forgery?** Our lattice might be a proper sublattice (missing relations) — exactly what the baseline got wrong on tiny groups. Proof requires Λ ⊆ L′, which `complete` checks.
8. **Why is testing a relation cheaper than finding one?** A test is a single Hadamard circuit costing Σ bit_length(x_i) ≈ 60 queries per shot; ~20 shots decide with error ≤ 1e-4. Finding needs b-bit phases of 4 elements in one shot (2,664 queries) and several shots.
9. **False-relation pass probability?** A non-relation has P(1) = ½ exactly; we accept a relation only with < k ones in n shots, chosen so P(Bin(n, ½) < k) ≤ 1e-4 (e.g. k=3, n=24 → 1.8e-5).
10. **What breaks first as noise grows?** Single stage flips in the high-weight stages (errors 2^j with large j), which corrupt whole shots; then the tests' P(1) for true relations grows (≈π²σ²·cost). See E2.
11. **Modular equations with non-coprime coefficients?** We never divide modulo n: everything is integer HNF/back-substitution; divisibility of pivots decides existence of logs.
12. **What is p-torsion of ℤ⁴/L′ and why test lines?** Elements x with p·x ∈ L′. If L′ ⊊ Λ, Λ/L′ has an element of prime order p, which is a non-zero p-torsion element; up to scalars (and L′) there are (p^r−1)/(p−1) of them (r = p-rank).
13. **Why can a cleared prime stay cleared?** If no p-torsion line is a relation, [Λ:L′] is coprime to p; for L′ ⊆ L″ ⊆ Λ, [Λ:L″] divides [Λ:L′].
14. **When does skew mode run and why is it still safe?** When only 3 short relations exist (4th basis vector longer than b bits can resolve, e.g. tiny logs). We measure one element's order with 80-bit QPE, form M + ord(e)·e_e (a sublattice of Λ) and let the SAME two checks finish; nothing is trusted unchecked.
15. **What guarantees no timeout or budget overrun?** `_Ctx.run` refuses any run costing more than the remaining budget, a 400 s monotonic deadline is checked before every run and in Pollard rho, LLL calls are capped at 400, completion steps at 64; on exhaustion we return a best-effort integer, never "none".

