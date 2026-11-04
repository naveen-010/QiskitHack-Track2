# The Silent Group — Qiskit Fall Fest 2026, Track 2

Submission: **`solver.py`** (one file, defines `solve(oracle)`).

## The observation the whole solution rests on

Reading `silentgroup/simulator.py`, the sealed register never entangles with your
qubits. Per **shot** the device draws one character `χ` of `G` uniformly at random
and holds it fixed for the entire circuit, and `oracle_gate(h, k)` is then exactly
a phase gate on its control qubit:

```
diag(1, exp(2πi · k · χ(h)))
```

So this is not a hidden-subgroup problem needing a big register — it is phase
estimation of the angles `θ_h = χ(h)`, and **one qubit** suffices (the 16-qubit cap
is never approached). Two primitives follow.

### 1. An exact zero test

`H → oracle gates → H → measure`, over many shots. Character orthogonality gives
`E_χ[χ(h^k)] = 1` if `h^k = e` and `0` otherwise. So when `h^k = e` every shot
returns 0; otherwise each shot is an independent fair coin. Chaining gates on one
control tests any word `g^k₁ A^k₂ B^k₃ C^k₄ = e`.

This is exact and essentially noise-immune — the phase jitter enters only as
`cos 2πε`. It is the **final arbiter**: the solver returns nothing it has not
confirmed with this test. Accepting at ≤2 ones in 20 shots gives a false-accept
rate of 2·10⁻⁴ and a false-reject rate of 2·10⁻⁵.

### 2. Single-qubit iterative phase estimation

One qubit, `reset` between rounds, semiclassical Griffiths–Niu feedback via
`if_test`. It must live inside **one circuit**, because `χ` is redrawn every shot —
splitting rounds across runs would estimate bits of different angles. Feedback
rotations are not oracle gates, so the budget never charges them.

Round `j` applies `h^(2^(bits-1-j))` and measures bit `j`, least significant first,
so the costliest and noisiest round fixes the least important bit. The flip side is
that feedback makes an error **cascade upward** — one wrong bit throws the next
correction by a quarter turn — so a sweep is all-or-nothing and gets validated
rather than trusted.

## Order of g

`g` has order `n`, so `χ(g)` is always a multiple of `1/n`, with `χ(g) = s/n` for
`s` uniform on `Z_n`. The best rational approximation of `θ_g` recovers
`q = n/gcd(s,n)`, a divisor of `n`; the lcm over a few characters reaches `n`.

Two details matter:

- **Self-certifying sweeps.** Sweeping at `2·60+24` bits rather than the `2·60+3`
  Shor requires leaves 24 bits of headroom. A clean sweep's convergent explains `θ`
  to ~`2^-bits`; a cascaded sweep is near-random and its best denominator-≤`2^60`
  rational sits ~`2^-(bits-24)` away. Demanding agreement to `2^-(bits-8)` admits
  clean sweeps (they clear it by ~256x) and lets a cascaded one through only with
  probability ~`2^-16`.
- **`g^L = e` proves only `n | L`, never `L = n`.** An lcm polluted by one bad sweep
  passes that test while being a *multiple* of the true order — this was a real bug,
  producing answers of exactly `n·2^k`. The order is therefore finished by dividing
  primes out of `L` under the zero test, which terminates at the unique minimal `L`
  with `g^L = e`. That is what makes the answer exact rather than merely likely.

A precision ladder (64 → 88 → 112 → 144 bits) runs first because cascade risk grows
as `bits²`: the cheap rungs finish small groups outright, cost less than one full
sweep combined, and cut the all-sweeps-cascaded failure mode from ~0.5% to ~10⁻⁴.

## Discrete logs and forgery

If `A = g^a` then `χ(A) = a·χ(g)` for **every** character. Estimating `θ_g` and
`θ_A` in the same shot — hence under one `χ` — gives `s` and `t` with
`a·s ≡ t (mod n)`. One character with `gcd(s,n) = 1` pins `a` down alone; when it
does not, several characters intersect to one value, which matters for groups with
small `φ(n)/n`. Every candidate is confirmed by the zero test on `A·g^{-a}`.

Forgery is detected two ways, and **not** by the obvious test:

- `X^n ≠ e` ⟹ `X ∉ ⟨g⟩`, decided exactly by the zero test.
- Otherwise, two independent characters each yielding a candidate the zero test
  refuses means no single exponent fits — `X ∉ ⟨g⟩`.

The naive `X^n = e` check alone is **not sufficient**: the forged broadcast in
`public_2` satisfies it. That group is non-cyclic (moduli `[10628367750, 12297285]`,
exponent `2^36` ≪ order `2^57`), so `B` sits outside `⟨g⟩` while still being killed
by `n`. Only the cross-character consistency check catches it.

A large residual is also ambiguous — it means *either* a cascaded sub-sweep (~9% per
element per shot) *or* a genuinely non-`1/n` character. Guessing conflates the two
and cries forgery on honest instances; the `X^n` test separates them.

## Reproducing

```sh
git clone --depth 1 https://github.com/ProAltro/qff-hackathon.git upstream
python3 -m venv .venv && .venv/bin/pip install -e ./upstream
.venv/bin/python -m silentgroup.check solver.py --repeat 5   # official rehearsal
.venv/bin/python bench.py 10                                 # 13 group shapes
```

`bench.py` scores with the judge's own `silentgroup.check.score` over the two public
meetings plus eleven `MockGroup` shapes — cyclic and non-cyclic, prime and smooth
order, `ord(g)` far below `|G|`, and four forgery layouts. Set `solver.DEBUG = True`
for a per-sweep trace.

## Qubit and gate budget

| | |
|---|---|
| Qubits per circuit | **1** (cap is 16) |
| Classical bits | up to `4 × (log₂n + 12)`, ~280 |
| Oracle queries | ~5k–50k of 80,000, depending on group shape |
| Gates on our qubits | `h`, `p`, `reset`, `measure`, `if_test` only |

Using a single qubit also answers one of the open-track prompts ("solve with as few
qubits as possible") outright.
