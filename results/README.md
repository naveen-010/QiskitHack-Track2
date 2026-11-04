# Validation results

All final-solver records below refer to the unchanged `my_solver.py` checksum in
[`combined_sha256.txt`](combined_sha256.txt).

| Validation | Correct | Evidence |
|---|---:|---|
| Seeded comparison suite | 184/184 | [Summary](benchmark_combined_final.json), [raw runs](benchmark_combined_final.csv) |
| Independent-seed suite | 328/328 | [Summary](benchmark_combined_final_fresh.json), [raw runs](benchmark_combined_final_fresh.csv) |
| 21 shapes, 40 seeds each | 840/840 | [Summary](regress_combined_final_stress.txt), [raw runs](regress_combined_final_stress.csv) |
| Official public checker, 20 repeats each | 40/40 | [Checker output](combined_checker.txt) |
| Unit tests | 118 passed | [Test output](combined_tests.txt) |

Total: **1,392/1,392 solver runs correct**, with a maximum observed cost of
**71,094 queries** and maximum measured elapsed time of **8.991 seconds**.
These observations are not guarantees about unseen instances.

[`combined_validation.json`](combined_validation.json) collects the measurements,
source checksum, dependencies, and provenance in one machine-readable record.

## Development comparison

On the same 184-case suite, averages include unsuccessful runs:

| Version | Correct | Mean queries | Mean score / 1.5 |
|---|---:|---:|---:|
| Naveen prototype | 177/184 | 37,215 | 1.230510 |
| Pulkit prototype | 183/184 | 23,167 | 1.351402 |
| Final submission | 184/184 | 20,236 | 1.373527 |

The prototype rows are preserved in [`benchmark_originals.csv`](benchmark_originals.csv).
Its older label `large_prime_rank2` refers to the same group `Z_1009^4` now named
`prime_rank4_forgery` in the benchmark harness; it does not assert rank two.

Reproduction commands and the exact environment are in
[`docs/REPRODUCIBILITY.md`](../docs/REPRODUCIBILITY.md). The archived unit-test log
includes the old duplicate-copy packaging check; the current equivalent check
verifies the single solver against its frozen checksum instead.
