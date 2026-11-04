# Algorithm and correctness

The standalone implementation is [`my_solver.py`](../my_solver.py). It uses one recycled qubit and interacts with the device only through its public oracle interface.

## Algorithm

1. Check the identity case before phase estimation. Confirming `g=e` gives order 1;
   then test whether all broadcasts are also the identity.
2. Try a cheap 48-bit, single-element rational order measurement. A denominator is
   only a proposal: two independent device tests must establish that it annihilates
   `g`. This path targets orders up to roughly `2^20` without paying for a lattice.
3. Minimize any verified order multiple. Test `n/p` first, because most proposed
   orders are already minimal. If a prime can be removed, use Naveen's prime-power
   bisection, accepting a reduction only after two independent confirmations.
   Minimization includes order 1.
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

The calibrated thresholds target the competition's stated noise. Finding a key or
proving a forgery does not necessarily recover the entire subgroup structure.
Full-structure and elevated-noise experiments from earlier prototypes are not
presented as results of this submission.

## Open extension: one-qubit implementation and resource analysis

Every circuit uses one control qubit, recycled by measurement and reset. All
same-character measurements remain within one circuit. This meets the device's
16-qubit limit with a factor of sixteen fewer allocated control qubits.

At least one control qubit is necessary in this oracle model: an oracle gate acts
on a control qubit, so a zero-qubit circuit cannot query the device. Public costs
and element names alone cannot distinguish arbitrary hidden group instances.
Thus one is the minimum possible number of allocated control qubits for a general
solver using this interface. This statement concerns the provided implicit-register
model; it does not count or eliminate the physical group register.

For b phase bits on d named elements, the exact oracle cost per shot is

    d * sum(m for m in 1..b) = d*b*(b+1)/2.

The four-element 36-bit circuit therefore costs 2,664 queries per shot, compared
with 10,440 for a single-element 144-bit circuit. The 12-stage feedback window uses
at most `d * (12*b - 78)` conditional rotations when `b >= 12`. Feedback rotations
do not add oracle queries. The circuit stores `d*b` classical bits.

For relation verification, the per-shot cost is

    sum(max(1, bit_length(abs(x_i))) for nonzero x_i).

The calibrated sequential test stops once its fixed-sample decision is determined.
The forgery witness checks p+1 lines in a two-dimensional slice instead of
`(p**r - 1)/(p - 1)` lines in the entire rank-r torsion space. This reduction does
not relax the slice's exhaustive injectivity check.

Total cost also depends on rejected candidates, repairs, and sampled characters;
the solver checks the remaining budget before every run. Actual per-instance
query counts are available in the [recorded results](../results/README.md).
