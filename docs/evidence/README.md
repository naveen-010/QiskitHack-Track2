# Walkthrough evidence

These files support the [development walkthrough](../PRESENTATION_WALKTHROUGH.md).
They preserve the reasoning behind individual tests and the recorded development
results without adding experimental solver variants to the submission entry point.

| File | Contents |
|---|---|
| [TEST_CATALOG.csv](TEST_CATALOG.csv) | All 118 current pytest cases, their real parameter values, purpose, expected assertion, and source file |
| [TEST_RUN_LEDGER.csv](TEST_RUN_LEDGER.csv) | 3,832 retained solver executions across original, intermediate, frozen final, and packaging versions |
| [DEVELOPMENT_TIMELINE.csv](DEVELOPMENT_TIMELINE.csv) | 28 ordered checkpoints, changes, failures, reasons, and evidence references |

## Reading the counts

| Ledger group | Runs |
|---|---:|
| Two original solvers, 184 runs each | 368 |
| Intermediate comparisons and stress tests | 2,028 |
| Frozen final validation | 1,392 |
| Packaging benchmark and regression smoke tests | 44 |
| Total retained ledger rows | 3,832 |

Only rows with `counts_in_final_1392=1` belong to the headline final validation.
Repeated instances across versions are separate executions, not independent new
problem constructions. The 118 unit cases and repeated pytest invocations are
separate from this ledger. Six additional public packaging checks are confirmed
by the validation manifest but excluded because their raw per-run file was not
retained.

The catalog was regenerated with pytest collection on 5 November 2026. Collection
enumerates the cases; it does not claim that a new test execution occurred. The
passing-suite evidence is in the existing [validation manifest](../../results/combined_validation.json)
and [test log](../../results/combined_tests.txt).

## Field meanings and provenance

The ledger keeps exact stored correctness flags, query counts, results, expected
answers, seeds and timings where available. `logic` explains why that construction
was selected. `source_record` names the original CSV or checker output and its
row/line number; the first data row of a CSV is line 2. Row order is grouped by
checkpoint, not the concurrent workers' completion timestamps.

Sources under `results/` remain in the repository. Sources labelled `archived
development/` or `archived packaging/` were local experiment CSVs consolidated
into this ledger. Their normalized results are retained here; those labels are
provenance descriptions, not links to files in this tree. A separate immutable
source-code snapshot for each intermediate candidate was not retained. Candidate
names and the cumulative development narrative therefore do not constitute a
controlled ablation or a promise that every intermediate executable can be rerun.

Original comparison CSVs omit the oracle seed column. The ledger fills it using
the recorded schedule `10000 + seed`. Regression oracle seeds are `100 + seed`.
The independent comparison schedule uses `20000 + seed`. Public comparison cases
use a fixed instance, so their seed chooses oracle randomness rather than new
public group elements. Official checker runs have no recorded seed.

The old label `large_prime_rank2` is normalized to `prime_rank4_forgery`; the
original label remains in its own column. The construction is `Z_1009^4`.
Regression result dictionaries are normalized from Python representations to
JSON; their values are unchanged.

Comparison scores are retained directly. Regression CSVs did not store scores,
so their scores are derived from the official rule using correctness and queries.
Official-checker scores are likewise calculated from exact query counts because
the printed scores have only three decimal places. `score_basis` distinguishes
these sources. Official output does not contain result dictionaries, so the
ledger leaves those fields empty rather than reconstructing unrecorded output.
Reported timings retain the precision in the source files.

Some unit cases deliberately inject a denominator or an incomplete candidate
lattice. These isolate a proof or control-flow property; they do not claim the
top-level solver naturally produced that exact injected state in a benchmark.
Timeline rows marked as session records summarize observations from development;
they are distinguished from surviving raw logs and Git timestamps.

CSV files use UTF-8 with a byte-order mark so names and explanations open cleanly
in common spreadsheet applications. No spreadsheet application is required to
read them.
