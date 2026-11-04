# WP2 / WP3 change log (before → after), 4 Nov 2026

All regressions: `engineering/regress.py`, MockGroup seed = s, Oracle seed = 100 + s, noise 0.002,
budget 80k. "ok" = order and secret both right. Files: `results/regress_<tag>.txt/.csv`.

## Baseline (`my_solver.py`, tag `baseline`, 20 seeds)
| shape | ok/20 | failure kind |
|---|---|---|
| core shapes (cyclic 2^58, p·q, 2^20×6×3^8, 2^55, 2^30×3^18, public1) | 20 | – |
| 2^29² | 19 | 1 false "none" |
| tiny 1000003 / 7 / 64 / 1001 / 65537 | 17 / 16 / 18 / 16 / 18 | false "none" (most), wrong order (2) |
| forged1_C | 17 | order 1 returned; one 1386 s stall |
| public2_shape | 17 | 3 false "none" |
| public2_shape_forgA | 16 | missed forgery, 1 harness timeout (stall) |
| skew (a,b,c tiny) ×3 shapes | 0 | stalls / wrong order |

## Root cause of the tiny-group bug (WP2.2)
Measured with `engineering/debug_tiny.py` and `engineering/debug_phase_err.py` (brute-force the
character of every shot in Z_1000003):
* seed 1: shot 0 has phase error exactly −2^25 on g (one noise-flipped QPE stage). All short LLL
  vectors are TRUE relations, but they span a sublattice of index 2^11 (HNF diag [1,1,1,2048·1000003]).
* seed 4: shot 1 has error +2^21 on C → sublattice of index 2^15 in the C coordinate → C looks
  outside <g> → baseline answers "none" (false forgery) because verification passed.
* Hypotheses (a) "good[:4] dependent/wrong" and (b) "tolerance badly tuned" were **falsified** for this
  failure (all good vectors were true relations). (c) was right in a stronger form: verification only
  checked L' ⊆ L, never L ⊆ L'. In big groups the det ≤ 2^60 bound happens to catch sublattices; in
  small groups it cannot.
* Baseline also never fed a verification failure back: `consensus` kept returning the same lattice,
  and the budget-exhaustion fallback returned it unverified.

**Fix:** `complete()` — for every prime p | det(L'), compute the p-torsion of Z^4/L' (left kernel of the
basis mod p) and Hadamard-test one vector per line. A relation found is added to L' and the step
repeats. Theorem: L' ⊊ L ⇒ L/L' has an element of prime order p, p | det(L'), and it is one of the
tested lines (exhaustive while the number of lines ≤ 40).

## Other changes, each motivated by a measured failure
| change | evidence |
|---|---|
| canonical HNF (reduced above pivots) | stability rule never fired: same lattice had different keys |
| `torsion_candidates` enumerates line representatives directly | first version looped `range(p)` for p ~ 2^50 (hang) |
| keep only first 4 *independent* good LLL rows | `hnf(all good)` admitted false rows (head≈tail≈2^21) → 44k queries on public_1; skew data gave 4 good but rank-3 rows |
| non-cyclic quotient: ≥ 6 samples and the lattice seen at the 2 previous counts (vote) | WP2.3; literal "same as t−2" missed a stable lattice after a corrupted 7th shot |
| precision ladder 36 → 44 → 52, caps 7/8/8 samples | forged A with order 2^25.9 mod <B,C,g> needs b ≥ 44 |
| skew mode: 3 verified short relations + ord(e)·e_e from 80-bit single-element QPE, then the same two checks | skew shapes 0/20 → 20/20 |
| cache of tested vectors; free rejection of lattices containing a known non-relation | public2_forgA spent 2k per re-tested candidate |
| completion: a cleared prime stays cleared (theorem: [L:L'] coprime to p is inherited) | completion cost grew to 13k per step |
| noise-aware relation test (n, k) with false-accept and false-reject ≤ 1e-4 | skew seed 3: true relation of cost 220 rejected (P(1)≈0.9 %) → wrong order 2^58 |
| curtailed (sequential) tests | identical decisions to the fixed-n test by construction (stop only when the outcome is decided) |
| time guard 400 s, LLL-call cap 400, budget check before every run, VERIFY_RESERVE 8000 | baseline stalls of 1386 s and harness timeouts |

## After (tag `wp3e`, 20 seeds; plus `forged_extra`, 100 extra seeds on the 3 hardest forged shapes)
* 419/420 on the 21 shapes, 0 false forgeries, 0 exceptions, max 10 s per solve.
* The one miss (forged2_AB seed 13): all 15 samples had χ(h)=0 for one order-2 element h (P≈2^-15
  for fair characters; the seeded stream for this seed starts with a run of even parities); best-effort
  answer, never "none". 100 fresh seeds of the same shape: 100/100.
* forged extra seeds: forged2_AB 100/100, forged1_C 100/100, public2_shape_forgA 99/100 (budget
  exhausted, best-effort with correct order).
* Cost: median queries 11–17k on ordinary/tiny shapes (baseline 11–15k), 22–24k on forged non-cyclic
  shapes, 32–38k on skew and forged-A shapes.
