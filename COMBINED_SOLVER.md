# Combined Track 2 solver

The standalone submission is `submission/my_solver.py`. The root `my_solver.py` is
an identical copy, so the default README command also runs the current solver.
Both use one recycled qubit and only the standard library, Qiskit, and
`silentgroup.oracle_gate`. Neither reads files or inspects the hidden instance.
The only oracle methods used are `run`, `cost`, and `budget_left`.

Base: Pulkit's `submission/my_solver.py` at `a067e0ad3eec13904e92e784e5ff61580b4fdaa7`.
Additional algorithmic ideas: Naveen's `solver.py` at
`82e17e1be64df3d24956e60bc756f09f3f4da428`. Both original branches remain unchanged.
The combined implementation is on the local branch `combined-solver`.

## Algorithm

1. Check the identity case before phase estimation. Confirming `g=e` gives order 1;
   then test whether all broadcasts are also the identity.
2. Try a cheap 48-bit, single-element rational order measurement. A denominator is
   only a proposal: two independent device tests must establish that it annihilates
   `g`. This path targets orders up to roughly `2^20` without paying for a lattice.
3. Minimize any verified order multiple. Test `n/p` first, because most proposed
   orders are already minimal. If a prime can be removed, use Naveen's prime-power
   bisection, accepting a reduction only after two independent confirmations.
   Unlike the original main branch, minimization can reach 1.
4. For a known small order, use same-shot phase pairs and modular congruences to
   recover the logs. Verify each positive log against the device. Omit logs already
   solved from subsequent circuits. Treat insufficient samples/budget as unresolved.
5. Otherwise use Pulkit's joint four-element relation lattice, with phase precision
   36, 44, then 52 bits. Use exact-rational LLL, canonical integer triangular bases,
   sample consensus, independently verified generators, and outlier subsets.
6. A verified lattice can be incomplete. Minimize its proposed order of `g` before
   answering. If all three positive logs are available, their verified relations
   plus a minimal order suffice: full lattice saturation is unnecessary. Add a fresh
   eight-shot audit of each extracted log; suspicious audits receive a full test.
7. Detect forgeries with a twice-rejected `X^n=e` check, or with the non-cyclic witness
   below. These often avoid certifying the complete subgroup structure.
8. If necessary, complete the relation lattice, with a 20,000-query allowance per
   attempt and an 8,000-query retry reserve. Keep verified repairs and cleared
   primes between attempts. Exhaustively check small torsion spaces; a capped
   space returns unresolved, never a false completeness certificate. Reuse partial
   progress immediately for direct log recovery before buying more measurements.
9. Very skewed lattices need a longer direction. Find three tiny, verified relations
   across outlier subsets, measure a complementary element at 80 bits, minimize its
   order, and insert the resulting relation. After unsuccessful low-precision joint
   recovery, also try 88-bit rational order recovery. Combine denominators and close
   cheap missing factors with independently tested multiples. A bounded 144-bit
   rational fallback remains available if no order survived.

Each shot must measure all phases it compares inside a single circuit. The sampled
character changes between shots, not between successive oracle calls in a shot.
The implementation uses this operational equivalence; it makes no claim that the
physical group register never becomes entangled with a control qubit.

## Why the new forgery witness works

Write the candidate relation lattice as L', and the actual one as L. Once all
generators of L' are verified, L' is contained in L, conditional on the statistical
checks. Thus Q = Z^4/L' maps onto the actual subgroup H = <A,B,C,g>.

If all broadcasts are honest, H = <g> is cyclic. A cyclic group cannot contain two
independent elements of the same prime order p.

Compute a two-dimensional slice of Q[p] using the left kernel of the basis modulo
p. It has p+1 one-dimensional lines. If a representative of every line maps to a
nonidentity element of H, the slice maps injectively into H. Consequently H contains
Z_p x Z_p, so at least one broadcast is forged.

The solver checks **every** line in this slice and independently repeats each
nonidentity decision. It skips primes whose p+1 lines exceed the cap. It does not
sample some lines and pretend to have checked them all. It does not need to prove
that L'=L before accepting this witness. The order of g is established separately.

This is particularly useful when a forged broadcast satisfies X^ord(g)=e: the
usual annihilation check alone cannot establish membership in <g>.

## Statistical confidence and limitations

For a relation word costing W queries per shot, at the competition noise
sigma=0.002 the outcome-1 probability is

    p_identity = (1 - exp(-2*pi*pi*sigma*sigma*W))/2.

For a nonrelation, character orthogonality gives probability 1/2. The calibrated
test selects a shot count and threshold targeting both error types at most 1e-4
per full test, then stops early once that test's decision is fixed. Repeated
confirmations protect order reductions, added completion relations, and negative
membership witnesses. Final positive-log audits are additional checks of relations
already implied by the tested basis; they are not standalone 1e-4 tests.

Finite sampling is not zero-error certification. Errors can accumulate across
tests. Phase residual filters, LLL reconstruction, outlier consensus, and the direct
log fallback remain probabilistic. Incompatible direct-log data is allowed to
establish a forgery only after at least five accepted phase pairs remain
incompatible even when any two outliers are removed.

The interface has no abstention sentinel. If every recovery route is exhausted,
the solver preserves its best order information and returns integer secret 0 as
an explicit best-effort guess. It does not label exhaustion as a forgery.

The calibrated thresholds target the competition's stated noise. Historical
noise-sweep figures in `results/` do not measure this new solver. Likewise,
historical full-structure results belong to the previous version. The adapted
`engineering/experiments.py e4` separately attempts full saturation after a fast
forgery decision, and honestly reports unresolved structure if its budget or
exhaustive-search limit prevents certification. Finding the key or a forgery does
not by itself recover the entire subgroup structure.

## Measured results

The final files have SHA256
`5d3635eb8439b62a07ef19d70eb1922976ee451e50f1b549570dc94d15c81e85`.

Same 184-case comparison, including failed runs in query/score averages:

| Solver | Fully correct | Mean queries | Mean score / 1.5 |
|---|---:|---:|---:|
| Naveen original | 177/184 | 37,215 | 1.230510 |
| Pulkit original final | 183/184 | 23,167 | 1.351402 |
| Combined final | **184/184** | **20,236** | **1.373527** |

Additional final-version validation:

- Independent-seed suite: **328/328**, mean 19,548 queries, mean score 1.377824.
- Original 21-shape, 40-seed regression: **840/840**, mean 19,489 queries.
- Official unseeded checker: **40/40**.
- Unit tests: **118 passed**, including the reproduced verification counterexamples.
- Across these **1,392 solver runs**, the maximum cost was **71,094 queries** and
  the maximum measured elapsed time was **8.991 seconds** under concurrent testing.

The combined solver is not cheaper on every individual shape. It pays for stronger
order checks and positive-log audits on ordinary large honest instances, while its
small-order handling and focused forgery witnesses improve the selected suite as
a whole. Raw rows, per-case summaries, and exact provenance are in `results/`.

## Reproduction

Use the existing environment or install the package dependencies from
`pyproject.toml`. No extra runtime dependency is needed by the submission.

```powershell
# Official checker, OS-randomized measurements:
.\.venv\Scripts\python.exe -m silentgroup.check submission/my_solver.py --repeat 20

# Unit tests and regression counterexamples:
.\.venv\Scripts\python.exe -m pytest tests -q

# Same 184-case suite used to compare the original branches:
.\.venv\Scripts\python.exe engineering/benchmark_combined.py --tag combined_final

# Independent seeds, 328 runs:
.\.venv\Scripts\python.exe engineering/benchmark_combined.py --tag combined_final_fresh --seed0 8 --seeds 16 --public-repeats 20 --offset 20000

# Original 21-shape, 40-seed regression suite:
.\.venv\Scripts\python.exe engineering/regress.py --seeds 40 --tag combined_final_stress --workers 4
```

The benchmark accepts `--solver path/to/other_solver.py` for a reference file. The
original comparison rows are preserved in `results/benchmark_originals.csv`.
Its historical label `large_prime_rank2` denotes the same Z_1009^4 construction
now correctly labelled `prime_rank4_forgery`; use the explicit harness definition.
Instance answers are used only by test and benchmark code, never by the solver.

Final measured results and artifact hashes are recorded in
`results/combined_validation.json` and `submission/SUBMISSION_CHECKLIST.md`.
The selected case distribution is not the unknown judge distribution, and passing
these experiments is not a guarantee of perfect hidden-instance performance.
