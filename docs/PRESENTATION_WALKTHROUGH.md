# The Silent Group Development Walkthrough

Preparation guide for the Track 2 class presentation

Prepared for Naveen and Pulkit on 5 November 2026. This document explains the submitted solution, the reasoning behind its two starting approaches, the failures reproduced during review, the development timeline, and the purpose of the tests. It is a presentation guide, not the slide deck itself.

## Executive summary

We had to recover a three-party Diffie-Hellman secret inside a finite abelian group whose elements and register we could not inspect. Our only useful access was controlled multiplication through an oracle. The challenge combined algebra, noisy quantum measurements, a strict query budget, and detection of dishonest public broadcasts.

Both original implementations recognized the central opportunity: within one circuit shot, the observable statistics are equivalent to using one randomly sampled group character. Controlled multiplication therefore supplies a phase. This lets us estimate algebraic relationships using a single recycled control qubit.

Naveen's implementation found the order of `g` first, then solved modular equations for the three discrete logs. Pulkit's implementation estimated four phases together and recovered their relation lattice. Naveen's approach was cheaper for small orders and easier to verify claim by claim. Pulkit's final implementation was usually more economical on large instances and provided the stronger base for joint reconstruction.

Neither original was fully reliable. We reproduced an order-one bug and false forgery decisions in Naveen's solver. We reproduced budget exhaustion in Pulkit's solver and constructed a counterexample to its completeness check. The first combined versions also failed: stricter verification and repair limits initially consumed too much budget or stopped useful work too early. We used those failures to redesign the control flow rather than simply concatenate the two files.

The submitted solver keeps joint lattice recovery, adds cheap direct order recovery and minimization, preserves verified partial results, and distinguishes an unresolved attempt from a forgery. It also adds a useful forgery witness: demonstrating two independent elements of the same prime order proves that the generated subgroup is not cyclic, so the broadcasts cannot all be powers of `g`.

| Same comparison suite | Correct answers | Mean queries | Mean score out of 1.5 |
|---|---:|---:|---:|
| Naveen original | 177 of 184 | 37,215 | 1.230510 |
| Pulkit original final | 183 of 184 | 23,167 | 1.351402 |
| Final submission | 184 of 184 | 20,236 | 1.373527 |

The frozen submission then passed 328 independent-seed runs, 840 regression runs, and 40 official public-checker runs. Together with the 184-run comparison, that is **1,392 correct final-version solver runs**. The unit suite has **118 passing cases**. These measurements establish tested performance, not a zero-error guarantee for hidden instances.

For a short pitch, the central story is: **we turned phase measurements into candidate algebraic relations, then made the solver economical and careful about what those relations actually prove**.

## The exact Track 2 task

The official title is **The Silent Group**. Its stated goal is: “Break a three-party Diffie–Hellman exchange in a group you can never look at.” [S1]

Let `G` be a finite abelian group with `|G| <= 2^60`. Public broadcasts are `g`, `A=g^a`, `B=g^b`, and `C=g^c`. For each meeting, output:

- `order`: the smallest positive `n` with `g^n=e`.
- `secret`: `a*b*c mod n`, representing the key `K=g^(abc)`.
- `secret="none"` if any broadcast is outside `<g>`. The order is still required.

The oracle supplies controlled multiplication by a named element raised to an integer power. The hidden register cannot be read or manipulated directly. The rules allow at most 16 allocated qubits, 80,000 queries, and 600 seconds per meeting. A power `k` costs `max(1, bit_length(abs(k)))` per shot, including both conditional branches. Each call adds Gaussian phase noise with standard deviation `0.002*sqrt(cost)` turns. Adaptivity is allowed; only public oracle interaction is permitted. [S1]

Per meeting, a wrong order earns zero; a correct order with a wrong secret earns 0.3; a fully correct answer earns `1 + 0.5*(budget left / 80000)`. The website weights automated scoring at 35 percent, the open extension at 20 percent, and presentation at 45 percent. [S1]

This is not a request to print the hidden group element `K`, factor the group by inspecting simulator internals, or solve only cyclic groups. The required numerical exponent identifies `K` relative to `g`. The local checker and tests may read public/mock answers to evaluate us; the submitted solver may not use them.

## Why these two approaches were natural

### Start with the oracle rather than a named quantum algorithm

For an abelian group, multiplication operators share a character eigenbasis. Write the phase of character `chi` on element `h` as `phi_h`, so controlled multiplication contributes `exp(2*pi*i*k*phi_h)` to the control's one-state amplitude. The organizer's simulator implements an equivalent measurement model: one random character is fixed throughout a shot, and a new character is sampled for the next shot. [S2]

This immediately gives two useful tools. Iterative phase estimation measures `phi_h`, and a Hadamard test checks whether a product of element powers is the identity. It also explains a crucial circuit rule: phases that must refer to the same character must be measured within the same circuit. Measuring `g` in one run and `A` in an unrelated run would not generally give a valid equation between them.

The original branches already existed when our review began. The derivations below explain why their methods fit the oracle and match their code; they are not an invented account of the authors' earlier brainstorming.

### Naveen used order finding followed by modular equations

If `n=ord(g)`, a character takes the value `chi(g)=s/n`, with `s` uniformly distributed modulo `n`. A sufficiently precise rational approximation reveals a denominator

```text
q = n / gcd(s,n).
```

Thus one good measurement may give only a divisor of `n`. Combining denominators with an LCM can recover the missing factors. A proposed `N` is tested by checking `g^N=e`. That establishes `n | N`, so unnecessary prime factors must then be removed to obtain the minimal order.

Once `n` is known, same-shot phases for `g` and `A` give numerators `s` and `t` satisfying `a*s=t mod n`. If `s` is invertible, one equation yields `a`. Otherwise, several congruences or a short list of remaining candidates can identify it. Test `A*g^(-a)=e` before trusting the log. Repeat for `B` and `C`, then multiply the logs modulo `n`.

This approach is attractive because its intermediate claims are direct: an order multiple, a minimal order, and three verified discrete logs. The original solver also stopped measuring a broadcast once its log had been found.

### Pulkit recovered relations jointly

Instead of treating the broadcasts separately, define the integer relation lattice

```text
Lambda = {x in Z^4 : A^x0 * B^x1 * C^x2 * g^x3 = e}.
```

Every character turns a genuine relation into a phase equation: `x dot phi = 0 mod 1`. Several noisy four-phase samples therefore constrain short integer vectors. Pulkit used LLL reduction to find candidate relations and an integer triangular basis, called HNF in the code, to extract the order and logs.

The quotient `Z^4/Lambda` is isomorphic to the subgroup generated by the four public elements. Its size is at most `2^60`. For an ordinary well-balanced four-dimensional lattice, useful relations can be much shorter than the full group order. That makes relatively short phase measurements productive. This is a geometric heuristic about typical instances, not a guarantee that every basis direction is short.

The cost advantage is concrete. With `b` bits for `d` elements, a phase-estimation shot costs `d*b*(b+1)/2`. Four elements at 36 bits cost **2,664 queries**. A 144-bit sweep of `g` alone costs **10,440**. Four of the shorter joint shots can sometimes recover the whole answer for approximately the cost of one long order-only sweep.

### What both solutions already shared

Both used the same character interpretation, one-qubit iterative phase estimation, classical feedback, resets, same-shot phase comparisons, classical reconstruction, and Hadamard-style relation tests. Neither had a qubit-count advantage over the other. Their important difference was how they organized the information and verification budget.

It is misleading to say that the physical group register never entangles with the control. The valid statement is that the unobserved character branches give the same accessible measurement statistics as the sampled-character simulation. Resetting our control does not sample a new character within the shot.

## What was wrong with the starting implementations

### First distinguish Pulkit baseline from Pulkit final

The `pulkit` branch contained a short root `my_solver.py` baseline and a substantially hardened `submission/my_solver.py`. We compared Naveen against the latter. Treating the root baseline as Pulkit's final solver would give a misleading comparison.

The baseline's archived regression results already showed incomplete-lattice errors on tiny groups and failures on skewed cases with small secret exponents. Its earlier record had 0 of 20 successes on each of three skew constructions. These were historical prototype measurements, not new measurements of the final combined code. An old summary line marked `FAIL` could also mean its confidence-bound threshold was unmet: even 20 of 20 successes has a Wilson lower bound of only about 0.839.

### Naveen flaw one was the order one boundary

In `Z_1009`, set `g=e` and make all broadcasts honest. For instance seed 0 and oracle seed 10000, the original solver returned `order=2, secret=0`; the correct order is 1. A corrupted phase sample proposed a denominator `2^32`. The minimizer removed powers of two but had a hard floor of 2, so it could never finish the reduction to 1.

The final code checks `g=e` before expensive measurements and permits order minimization to reach 1. We test both honest and forged broadcasts with identity `g`, and separately test minimization of a large verified multiple down to 1.

### Naveen flaw two was treating an unfinished log search as forgery

For an honest instance built with modulus `2^59`, instance seed 1 and oracle seed 10001, the original solver spent **68,029 queries** recovering the correct order. It then could not afford its first full discrete-log measurement with the reserved confirmation budget. It returned `"none"`, despite having no evidence that a broadcast was forged.

Other fresh comparison failures appeared with small exponents, a zero exponent, and the exact `2^60` order boundary. The failure pattern was primarily unfinished recovery being reported as a trap. The final control flow keeps unresolved work separate from a forgery conclusion and uses the actual next circuit's cost to decide what is affordable.

### Pulkit flaw one was an incomplete completeness check

Verifying every generator of a candidate lattice `Lprime` establishes `Lprime` is contained in the true relation lattice. It does not establish equality. Missing relations can inflate the order or make an honest broadcast appear to have no discrete log.

We supplied the verifier with a valid but incomplete lattice in `Z_1009`, with `A=g`, `B=g^2`, and `C=g^3`. Its generators were `A^1009=e`, `B*g^-2=e`, `C*g^-3=e`, and `g^1009=e`. Every generator was genuine, but the relation `A*g^-1=e` was missing.

There were **1,010** relevant torsion lines. The original routine sampled only **40**, missed the essential direction, and still called the lattice complete. Its trusted output was `order=1009, secret="none"`, while the correct secret is **6**.

This was a constructed verifier-path counterexample. We did not claim that the ordinary solver produced that exact candidate during the benchmark. It nevertheless disproved the claim that the verification path always certified completeness. The final routine returns unresolved when its exhaustive search limit is exceeded.

### Pulkit flaw two was an expensive repair consuming the meeting budget

On `public_2` with oracle seed 10014, the original final solver entered a skew repair followed by lattice completion. It consumed **79,949 queries** and returned the correct order with integer secret 0 instead of the correct forgery answer. Its earlier 840 of 840 regression record did not cover every possible measurement sequence.

We bounded each completion attempt, saved verified repairs, cached completed prime checks, and used partial progress immediately. Most importantly, we added a way to prove a forgery without insisting on full subgroup reconstruction.

### Both descriptions overstated exactness

An identity check made from finitely many random shots is not logically infallible. Under noise, even a genuine relation can produce ones. A nonrelation can also produce a run of zeros by chance. Naveen's fixed threshold had false acceptance probability about `2.01e-4`; Pulkit calibrated tests to the word's cost and the stated noise.

The final solver retains calibrated statistical tests and repeats consequential decisions. Algebraic implications are exact conditional on correct evidence; measurement and reconstruction success remain probabilistic. A per-test error target is not an overall solver failure bound.

## The chronological development record

The dates verified by Git are 4 November 2026 for the reviewed Pulkit commit, the combined implementation commit, and the repository cleanup. The individual experiments below are ordered by the retained records and conversation sequence; we do not assign invented wall-clock times to them. The guide itself was prepared on 5 November.

### Stage one was review without changing code

We read the official task, inspected the oracle and simulator, and identified the two different architectures. We verified that the local simulator matched the branch source and snapshotted both solvers. We ran the original 96-test suite and built a common seeded comparison rather than comparing unrelated success anecdotes.

The comparison used 20 runs on each public meeting plus eight runs on each of eighteen additional constructions: **184 runs per solver**. Both used the same instance seeds, oracle-seed schedule, noise, budget, and scoring function. Different circuits consume different random draws, so identical seeds do not mean identical characters at every later measurement.

Naveen finished 177 of 184 correctly and Pulkit 183 of 184. This justified choosing Pulkit's final implementation as the base, while the small-order results justified preserving Naveen's direct method. We then reproduced the order-one, exhausted-log, exhausted-repair, and incomplete-certificate failures separately.

The choice was supported by differences within the same suite, not just its overall average:

| Construction | Naveen correct and mean queries | Pulkit correct and mean queries | What we learned |
|---|---:|---:|---|
| Order-seven group | 8 of 8; 4,576 | 8 of 8; 10,895 | Direct recovery avoids lattice overhead at tiny orders |
| Order-1000003 group | 8 of 8; 7,388 | 8 of 8; 11,996 | A cheap small-order probe is worth retaining |
| Official honest meeting | 20 of 20; 39,510 | 20 of 20; 18,074 | Joint recovery is a strong base for ordinary large cases |
| Binary modulus `2^59` | 7 of 8; 52,515 | 8 of 8; 16,312 | A long order-first sweep can leave too little budget for logs |
| Two long directions | 8 of 8; 31,708 | 8 of 8; 53,731 | Direct recovery remains useful when the relation geometry is unfavorable |
| Official forged meeting | 20 of 20; 45,689 | 19 of 20; 40,078 | Lower query cost alone cannot compensate for an unsupported final answer |

Query averages include failed runs and are rounded to whole queries. These observations motivated a shared controller with complementary routes; they do not establish which method wins on every unseen instance.

### Stage two combined the methods and initially regressed

We created `combined-solver`. The first combined candidate added direct order/log recovery, order-one handling, strict capped-search behavior, bounded completion, and direct positive-log checks. Its inherited 96 tests passed, but the larger comparison exposed a serious regression: **163 of 184** correct, including failures on all twenty `public_2` runs.

This was an important lesson: fixing a proof gap and adding safeguards can make the practical solver worse if the budget policy is too restrictive. A 12,000-query completion allowance stopped legitimate rank-two verification. Simply making checks stricter was not enough.

### Stage three adjusted verification and allocation

We raised the local completion allowance to 20,000 queries, cached cleared primes and known logs, and reduced repeated work. We tested `n/p` first before starting prime-power bisection: a minimal order can usually reject that single division immediately. We shortened the small-order probe from 64 to 48 bits, keeping independent device confirmation instead of trusting the rational residual.

For logs already implied by a verified lattice, we replaced repeated full tests with fresh eight-shot audits, escalating suspicious outcomes to full tests. Those audits supplement the tested basis; they are not standalone high-confidence certificates. This reduced ordinary large-instance overhead.

We also tried directing measurements to a single unresolved broadcast. Public forgery failures decreased, but some repair paths still discarded useful partial progress or spent too long rediscovering it.

### Stage four repaired geometry and reused partial results

A corrupted phase can distort the apparent geometry of the three short relations, not just the final fourth vector. We extended outlier-subset handling to the skew detector itself and required genuinely tiny verified spans before paying for a long-direction measurement. We allowed more early candidate rejections before escalating precision.

When a bounded completion attempt made progress and then stopped, we started using its repaired lattice immediately for log recovery. Previously, the solver could buy an entire new precision level and rediscover the same relations only after too much of the budget was gone.

### Stage five added a focused mathematical forgery witness

If every broadcast is a power of `g`, the generated subgroup must be cyclic. A cyclic group cannot contain two independent elements of the same prime order. We therefore looked for an injective two-dimensional prime-torsion slice in the image of a verified candidate quotient.

Testing every one of its `p+1` lines, with repeated nonidentity decisions, can establish that injectivity. This proves a forgery without certifying the entire lattice or checking every prime component. It addressed the main practical weakness of the earlier strict-completion versions.

We also added medium-precision direct order recovery after unsuccessful low-precision lattice work, so a difficult rank-two relation lattice did not force us to spend the whole meeting budget at higher joint precision.

### Retained intermediate measurements

| Candidate | Runs and scope | Correct | Main lesson |
|---|---|---:|---|
| V1 | Full 184-case suite | 163 of 184 | Repair cap too tight; all public forgery runs failed |
| V2 | 40 public runs | 37 of 40 | Larger allowance helped but did not solve partial-progress handling |
| V3 | Full 184-case suite | 178 of 184 | Lower verification cost helped; skew and forgery gaps remained |
| V4 | Public and two-long-direction subset, 48 runs | 46 of 48 | Targeting the remaining broadcast helped |
| V5 | Full 184-case suite | 181 of 184 | Robust geometry fixed skew equality; forgery and rank-two gaps remained |
| V6 | Public forgery and hard-shape subset, 36 runs | 35 of 36 | Reusing repaired partial results fixed the public forgery subset |
| V7 | Full 184-case suite | 184 of 184 | Focused non-cyclic witness and medium order fallback closed known gaps |

The rows are development checkpoints, not a controlled ablation study. Several changes accumulated between them, and the subset rows have different case distributions. Their mean scores should not be ranked against one another as though they were the same experiment.

Raw run records survived for these checkpoints, but a separate immutable source snapshot for every V1 through V7 candidate did not. The ledger preserves observed inputs, outputs, costs and failures; it does not make every intermediate binary independently reproducible. The original prototype commit identities and the frozen final source are separately identified.

### Stage six challenged the apparent success with fresh seeds

Passing the original 184 runs was not the stopping condition. An independent-seed 328-run suite exposed one more failure: an honest phase denominator was a proper divisor of the true order in a difficult two-direction forgery construction. The result was **327 of 328**, even though the earlier suite had been perfect.

The final fix combined useful rational denominators and tried a bounded set of cheap missing-factor multiples. Every proposed order multiple still required independent identity confirmation and minimization. This came directly from the arithmetic meaning of a sampled denominator; it was not a special case for the failing seed.

The final version passed the original 184 runs, all 328 fresh-seed runs, all 840 regression runs, and 40 official checker runs. The solver's frozen SHA256 is `5d3635eb8439b62a07ef19d70eb1922976ee451e50f1b549570dc94d15c81e85`.

### Stage seven prepared the submitted repository

Commit `ed6da6f` records the combined implementation; commit `76f9726` records the clean submission tree. We reduced the repo to one entry point, the unmodified organizer simulator, tests, reproducible scripts, documentation, and final evidence. We removed duplicate solver files, prompts, scratch code, old plans, and obsolete experimental outputs. Their development history remains in Git.

We checked installation, local documentation links, solver checksum, and unchanged simulator sources. We ran all 118 tests in the working tree and in a fresh clone, six public-checker smoke runs, four benchmark smoke runs, and forty regression smoke runs. The repo remained private, as requested. The current submission is the root `my_solver.py` on `combined-solver`.

## What we took from each solution

| Component | Source | Why it belongs in the final solver |
|---|---|---|
| Character interpretation and recycled-qubit phase estimation | Both | Establishes the legal measurement model and same-shot constraints |
| Hadamard relation testing | Both, with Pulkit calibration | Converts reconstructed claims into independent device checks |
| Joint four-element relation lattice | Pulkit | Makes short measurements useful on ordinary large instances |
| LLL, triangular extraction, and outlier consensus | Pulkit | Reconstructs order and membership from correlated measurements |
| Rational order proposals and denominator combination | Naveen | Provides economical small-order and fallback routes |
| Prime-power minimization | Naveen, corrected and optimized | Removes extra order factors and handles order 1 |
| Modular log equations and omitting solved broadcasts | Naveen | Focuses queries on the uncertainty that remains |
| Budget and time bookkeeping | Pulkit, expanded | Bounds attempts while preserving retry capacity |
| Reusing repaired partial lattices and verified log caches | Combined redesign | Avoids buying the same information twice |
| Two-dimensional non-cyclic forgery witness | Added during combination | Proves some traps without completing the entire lattice |

The final file is not two complete solvers run one after another with two separate budgets. It is one controller sharing the meeting's 80,000-query budget, verified relations, known logs, candidate orders, and stopping decisions.

## The submitted algorithm step by step

1. Test the identity boundary. If `g=e`, establish order 1 and check whether the broadcasts are also identities.
2. Try a cheap small-order phase measurement. Regard its rational denominator as a proposal, not a certificate.
3. Confirm a proposed annihilating exponent and minimize it. Test `n/p` first; use bisection only when a prime can actually be removed. Repeat accepted reductions independently.
4. If a useful order is known, recover logs from same-shot modular equations. Check positive log claims and stop measuring logs already solved.
5. Otherwise collect joint four-phase samples at 36, 44, then 52 bits. Recover candidate short relations with LLL and integer reduction; tolerate outlier samples.
6. Verify a candidate basis. Preserve only verified relation information for later fallbacks. Minimize its proposed order before interpreting the answer.
7. If the basis yields three positive logs, audit those claims and return their product modulo the verified order. Full lattice saturation is unnecessary for this positive conclusion.
8. If membership is unresolved, test cheap trap witnesses, including `X^n != e` and the non-cyclic slice witness. A capped search cannot clear an entire torsion space.
9. When needed, perform bounded completion or measure a complementary long direction of a skewed lattice. Save progress and reserve room for another attempt.
10. On total exhaustion, return the best available order information and an integer guess, rather than reporting a forgery solely because a search stopped. The interface has no abstention value; this fallback is not a verified success.

### Why the non cyclic witness is sound

Let `Lprime` be the verified candidate lattice and `Lambda` the true relation lattice. Conditional on the tests, `Lprime` is contained in `Lambda`, so `Q=Z^4/Lprime` maps onto the actual generated subgroup `H`.

Choose two independent order-`p` directions in `Q`. Their two-dimensional space has `p+1` one-dimensional lines. Every nonzero kernel would contain one such line. Therefore, if the image of a representative of every line is nonidentity, the map is injective on that slice. Its image contains `Z_p x Z_p`, so `H` is not cyclic. All-honest broadcasts would generate the cyclic subgroup `<g>`, giving the contradiction that establishes a forgery.

This checks a complete small slice, not a random sample of a large space. The order of `g` is verified separately. This distinction is what makes the witness useful without repeating the original completeness mistake.

## Why each kind of test was necessary

The complete unit catalog is [TEST_CATALOG.csv](evidence/TEST_CATALOG.csv). It has one row per collected pytest case, including the parameterized seeds. The [run ledger](evidence/TEST_RUN_LEDGER.csv) records retained solver-run results across versions. The following explanations give the reasoning behind the individual families and constructions.

### Five oracle semantics tests

| Test function | Logic and why we needed it |
|---|---|
| `test_one_character_per_shot` | In `Z_2`, two Hadamard blocks inside one noiseless shot must agree, while separate shots can disagree. This tests the same-character assumption on which both algorithms depend. |
| `test_measure_reset_if_feedback` | Measure, reset, and conditionally flip using the measured classical bit. The next outcome must reproduce the stored bit. This validates the recycled-qubit feedback primitive. |
| `test_cost_formula_both_branches` | Use positive, negative, zero, and large powers plus both branches of an `if_test`. Compare analytical and reported costs, and the budget decrement. This prevents undercharging adaptive circuits. |
| `test_identity_vs_nonidentity` | A known identity word should produce almost no ones at normal noise; a nonidentity averages one half. This checks the discrimination used by every relation test. |
| `test_noise_sigma_scales_with_sqrt_cost` | At an intentionally elevated diagnostic noise of 0.02, compare observed ones to `(1-exp(-2*pi^2*sigma^2*cost))/2`. This tests the noise formula, not full-solver robustness at elevated noise. |

### Twenty one integer algebra tests

Eight seeded LLL cases use independent random five-dimensional integer bases. We verify unchanged lattice membership and determinant, the Lovasz condition, and size reduction. Short vectors are not useful if the reduction changes which integer relations they generate.

Eight seeded HNF cases verify upper-triangular shape, positive pivots, determinant preservation, and equal lattices. One singular-input case must return `None` rather than fabricate a full-rank basis.

Four hand-built extraction cases compare against exact group arithmetic: a cyclic group of order 97; a forged broadcast in `Z_12 x Z_6`; a zero broadcast in `Z_8 x Z_4`; and a mixed-modulus large example. These isolate extraction logic from phase noise. They check both the final order and each log's membership/value, so a correct final product cannot hide compensating log errors.

### Seventy exact small group tests

Fourteen constructions run with instance seeds 0 through 4 and zero oracle noise. Nine honest choices use moduli `[8]`, `[15]`, `[8,3]`, `[4,8]`, `[7,11]`, `[7]`, `[64]`, `[1001]`, and `[65537]`. They cover powers of two, odd orders, products, and different subgroup shapes.

Five forged choices are `[12,6]` with `C` forged; `[4,8]` with `AB`; `[6,6,10]` with `ABC`; `[2,4,8]` with `B`; and `[7,7]` with `A`. They test that membership is not confused with a single element-order condition. Every case checks the exact answer, allowed output types, and the query budget.

Zero noise removes device phase jitter, helping separate algebraic or decoding defects from noise handling. Phase measurements and sampled characters still involve randomness; passing these cases is not a theorem about all noiseless outcomes.

### Twenty two combined solver regressions

| Test family | Cases | Purpose and expected conclusion |
|---|---:|---|
| Frozen artifact hash | 1 | The reviewed root solver must match the validation checksum. This replaced the old duplicate-file equality check after cleanup. |
| Identity `g` with honest or forged broadcasts | 2 | Recover order 1 correctly at the original failing noisy seed, with fewer than 1,000 queries. |
| Minimization down to 1 | 1 | Reduce a verified `2^32` multiple for identity `g` all the way to 1. |
| Rational proper-divisor completion | 1 | Inject a clean phase denominator half the true `2^30` order and recover the full order using verified multiples, within 8,000 queries. |
| Capped torsion certificate counterexample | 1 | The `Z_1009` incomplete lattice must remain unresolved or recover the genuine answer; it must not become a false trusted forgery. |
| Repair allowance with finite or unlimited oracle budget | 2 | A long sequence of valid repairs must stop within its own allowance, save progress, and retain retry queries. Unlimited external budget must not disable the internal cap. |
| Actual non-cyclic witness at primes 2, 3, and 7 | 3 | Genuine independent prime-order directions must establish non-cyclicity. |
| Incomplete lattice over an actually cyclic image, same primes | 3 | A candidate quotient that looks non-cyclic must not fool the witness when a slice direction is truly trivial in the group. |
| Insufficient log budget | 1 | Return unresolved without running a circuit; lack of budget must not assert forgery. |
| Previous failures and large torsion | 6 | Replay large-order, small-exponent, zero-exponent, maximum-order, and high-rank forgery cases at specific failure seeds. |
| Original public repair exhaustion seed | 1 | `public_2` at oracle seed 10014 must now return the right answer with fewer than 75,000 queries. |

The six replay cases are a random `Z_(2^59)` instance at seed 1; fixed-generator `Z_(2^55)` with exponents 2, 3, 5 at seed 1; fixed-generator `Z_(2^59)` with `a=0` at seed 1; fixed-generator `Z_(2^60)` at seeds 1 and 6; and `Z_1009^4` with three forged broadcasts at seed 3. The full catalog records each distinct test ID.

### The twenty constructions in the comparison suite

Each public meeting ran twenty times; each of the eighteen additional constructions ran eight times in the 184-run suite and sixteen times in the 328-run independent suite.

| Construction | Logic behind choosing it |
|---|---|
| `public_1` | Honest full-size rehearsal with the official checker semantics |
| `public_2` | Official forgery that cannot be identified just by checking `X^n=e` |
| `prime_56` | Modulus 72057594037927931 challenges high order with few useful small factors |
| `pow2_59` | Modulus `2^59` stresses noninvertible sampled numerators and missing powers of two |
| `prime_product` | Coprime moduli 1000000007 and 1000000009 test very large order without assuming a single coordinate |
| `smooth` | Moduli `2^18*3^5` and `2^6*5^3` stress shared factors and multiple congruences |
| `g_small` | In `Z_(2^30) x Z_(2^20)`, fix `g=(0,1)` to separate its order from the ambient group size |
| `tiny_7` | Detect unnecessary overhead and tiny-order correctness problems |
| `tiny_1000003` | Test a small prime order where a direct approach should be economical |
| `forged_B` | Moduli 12345678 and 60 with a forged `B` test a trap in a composite setting |
| `forged_AB` | Moduli 100000000 and 300000000 with two forgeries test multiple invalid broadcasts |
| `forged_ABC` | Moduli `2^28` and `3*2^28` with all three broadcasts forged test shared-prime structure |
| `skew_small_logs` | Fix a generator of order `2^55` and exponents 2, 3, 5 to create three very short relations and one long direction |
| `skew_equal_logs` | Fix a generator with modulus `2^58-27` and all exponents 1 to create extreme relation dependence |
| `identity_g_honest` | `g=0` in additive `Z_1009` tests the order-one boundary |
| `identity_g_forged` | Same identity generator with forged `A` tests membership when `<g>` contains only the identity |
| `zero_A` | Fixed generator of order `2^59` with `a=0` makes the honest secret zero |
| `exact_2pow60` | A fixed generator of `Z_(2^60)` reaches the allowed order boundary exactly |
| `prime_rank4_forgery` | `Z_1009^4`, `g=(1,0,0,0)`, and three forged broadcasts challenge large torsion spaces |
| `two_long_directions` | `Z_(2^30) x Z_(2^30)` with other honest logs 1 and forged `A` challenges rank-two geometry and proper-divisor order recovery |

Unless a generator is explicitly fixed, a modulus such as `2^59` does not mean every sampled `g` has that exact order. Also, a product of two cyclic groups with coprime moduli is itself cyclic; two coordinates alone do not prove a non-cyclic group. The equal-modulus constructions provide actual non-cyclic ambient groups.

### The twenty one regression shapes

The 840-run suite uses forty seeds per shape and a different oracle-seed schedule. These shapes overlap conceptually with the comparison suite but are independent execution records.

| Shape | Reason for the test |
|---|---|
| `cyclic_2^58` | Large modulus `2^58-27`; ordinary high-order relation recovery |
| `p1e9_x_p1e9` | Large coprime prime factors and potentially near-limit element order |
| `2^20x6x3^8` | Mixed prime powers and several coordinates |
| `cyclic_2^55` | Large binary modulus with noninvertible character samples |
| `2^29_sq` | Equal large binary moduli give a genuinely non-cyclic ambient group |
| `tiny_1000003` | Small prime order and low-cost fast path |
| `2^30x3^18` | Large coprime prime powers and long combined order |
| `tiny_7` | Very small order and random identity-generator possibilities |
| `tiny_64` | Small binary order and repeated factors |
| `tiny_1001` | Multiple odd factors and congruence intersections |
| `tiny_65537` | Small prime just beyond a common 16-bit boundary |
| `forged1_C` | One forged broadcast in a mixed composite product |
| `forged2_AB` | Two forged broadcasts with common factors in the ambient moduli |
| `forged3_ABC` | All broadcasts forged in a large shared-prime group |
| `forged_small_12x6` | Tiny forgery whose relation lattice must still be complete enough for the claim |
| `public1_shape` | New random elements with the first public group's moduli |
| `public2_shape` | Honest random elements with the second public group's moduli |
| `public2_shape_forgA` | A new forgery in that second shape, not merely replaying the public answer |
| `skew_2^55_a2b3c5` | Small exponents creating a long missing lattice direction |
| `skew_2^58_a1b1c1` | Equal exponents creating several dependent short relations |
| `skew_p1e9sq_small` | Small exponents 3, 7, 11 at large odd product order |

Forty seeds also make the harness's confidence criterion meaningful: even forty successes give a two-sided 95 percent Wilson lower bound of only about 0.912. Thirty successes out of thirty would give about 0.886, below its 0.90 threshold. This is why a small perfect sample was not sufficient for that gate.

### Why we used all four final validation layers

The 184-run suite made the before-and-after comparison fair. The 328-run suite challenged tuning to the first seed schedule and found the last missing-factor bug. The 840-run suite broadened shape coverage and repeated each shape enough times for the chosen confidence gate. The official checker used OS-randomized measurements and the organizer's own scoring path, preventing our custom harness from being the only evidence.

All four layers check exact answers, not merely absence of exceptions. Correct order with a wrong secret is recorded as partial score, not a successful solve. Time and budget are also checked. Four worker processes accelerate independent runs; they do not give any solver extra queries for a meeting.

### One off diagnostics and packaging checks

The four original counterexample diagnostics separated root causes: order-one minimization, log-budget exhaustion, repair-budget exhaustion, and the supplied incomplete lattice. A further numerical diagnostic compared relation-test error rates at costs 1, 20, 60, 120, and 240, confirming why a fixed threshold cannot have one universal false-rejection rate.

We also smoke-tested a legacy full-structure extension. Two selected examples recovered the expected structure, while another correctly remained unresolved after the new exhaustive-certification restriction. We did not turn that unresolved result into a claimed structure theorem, and removed the legacy extension outputs from the final submission. The chosen presentation extension is one-qubit implementation with resource analysis.

Repository checks served a different purpose from algorithm tests. The solver checksum prevented accidental code changes during cleanup. Comparing simulator sources ensured we were not validating against a modified oracle. Local-link checks caught documentation packaging issues. Editable installation validated the README setup. An initial diagnostic install without build isolation failed because the existing environment lacked a build backend; the ordinary documented isolated installation succeeded. This was an environment/setup distinction, not an algorithm failure.

Finally, tests in a fresh Git clone ensured the submission did not depend on deleted scratch files or ignored local artifacts. Additional smoke runs established that the renamed benchmark and regression scripts still loaded the single canonical solver correctly.

## How to interpret the final numbers

| Final version evidence | Correct solver runs | Mean queries | Mean score |
|---|---:|---:|---:|
| Original comparison schedule | 184 of 184 | 20,235.71 | 1.373527 |
| Independent schedule | 328 of 328 | 19,548.13 | 1.377824 |
| Forty seeds across twenty one shapes | 840 of 840 | 19,489.38 | 1.378191 |
| Official public checker | 40 of 40 | 25,394.55 | 1.341284 |

The maximum final-validation query count was **71,094**, and the longest measured solver execution was **8.991 seconds**. The unit-test runtime is a separate measurement and is not the per-meeting solver runtime.

The total 1,392 refers to the frozen final solver's four validation groups only. It does not include original-solvers' comparisons, unsuccessful intermediate versions, repeated unit invocations, or later packaging smoke tests. The larger companion ledger includes retained development and packaging runs too, clearly labelled by stage.

The final approach is not cheaper on every shape. Additional order checks and positive-log audits add cost to some large honest instances. Small-order handling and focused forgery witnesses improve the selected overall suite. We should say that the final version performed best on our recorded comparison, not that it is universally optimal.

## Suggested eight minute presentation

| Slide | Time | Main point and evidence |
|---|---:|---|
| The task | 45 seconds | Hidden abelian group, the order and exponent outputs, and traps |
| The oracle insight | 60 seconds | One fixed character per shot turns multiplication into phase information |
| Our two starting methods | 60 seconds | Sequential order/log recovery versus joint relation recovery |
| What broke | 75 seconds | Order 1, unresolved versus forged, and incomplete lattice certification |
| The final algorithm | 90 seconds | Joint reconstruction with direct recovery, verified minimization, and saved progress |
| The forgery witness | 60 seconds | Two independent prime-order directions cannot fit inside a cyclic subgroup |
| Tests and results | 60 seconds | Explain the four validation layers before showing 1,392 of 1,392 |
| Contribution and limits | 30 seconds | One control qubit, explicit query accounting, and honest statistical limits |

These times total eight minutes. Keep the 118-case catalog and detailed timeline as backup material; the talk should explain representative tests rather than read every seed aloud.

### A short opening you can adapt

“Our challenge was to recover a shared secret from a group we were never allowed to inspect. We could only request controlled multiplications, and every shot spent a limited noisy-query budget. We used the phases produced by those multiplications to recover algebraic relations with one control qubit. The hardest part was deciding which conclusions were actually justified when measurements were noisy or the search was incomplete.”

### Questions to prepare for

- **Why not compare phases from separate runs?** They usually correspond to different sampled characters, so the modular equation would not be valid.
- **Why is `g^N=e` not enough for the order?** It proves only that the actual order divides `N`; minimality requires ruling out removable prime factors.
- **Why is `X^n=e` not enough for membership?** In a non-cyclic group, an element outside `<g>` may still have order dividing `n`.
- **Why can a verified lattice still be wrong for a negative claim?** Its generators may be genuine while an essential relation is missing.
- **How did one qubit estimate many bits?** Measurement, reset, classical storage, and conditional feedback reuse the control within one shot.
- **Did we inspect answers to solve hidden meetings?** No. Only the separate test harness reads public/mock answers; the solver uses public oracle methods.
- **Does 1,392 successes prove zero failure probability?** No. The tests exercise selected cases; finite shots, reconstruction, and unseen instances remain probabilistic.
- **Did we recover the full structure of every subgroup?** No. A key or forgery can be established with less information; we do not claim universal full-structure recovery.
- **What is the open extension?** A one-control-qubit implementation, its minimum-qubit argument in this oracle model, and explicit query-cost analysis.

## Evidence and source guide

- **S1** [Official Track 2 problem statement](https://qff-website-pi.vercel.app/#hackathon). The website is used for the task and current scoring weights; an earlier participant README had outdated weights.
- **S2** [Organizer simulator](../silentgroup/simulator.py), [oracle](../silentgroup/oracle.py), and [official checker](../silentgroup/check.py).
- **S3** [Submitted solver](../my_solver.py) and [algorithm notes](ALGORITHM.md).
- **S4** [Final validation manifest](../results/combined_validation.json) and [results index](../results/README.md).
- **S5** [Original comparison records](../results/benchmark_originals.csv), with prototype commit identities in the manifest.
- **S6** [Unit case catalog](evidence/TEST_CATALOG.csv), [retained run ledger](evidence/TEST_RUN_LEDGER.csv), and [development checkpoints](evidence/DEVELOPMENT_TIMELINE.csv).
- **S7** [Reproduction commands](REPRODUCIBILITY.md), [comparison harness](../scripts/benchmark.py), and [regression harness](../scripts/regression.py).

The retained ledger documents individual solver executions for which raw records remain available, including failed intermediate candidates. It is not a count of every shell command or every repeated pytest invocation. Those repetitions are described in the timeline and test-family explanations. Original pre-review development is distinguished from our review and combination work throughout this guide.

The ledger contains **3,832 retained runs**: 368 original comparisons, 2,028 intermediate-version runs, 1,392 frozen final-validation runs, and 44 packaging smoke runs. The additional six public packaging checks are recorded in the validation manifest, but their per-run raw output was not retained as a file, so they are not added to this ledger. See the [evidence reading guide](evidence/README.md) for field definitions and provenance limits.
