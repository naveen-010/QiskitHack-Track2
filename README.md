# The Silent Group

**Qiskit Fall Fest 2026 | Track 2 | Final submission**

Recover the order of a hidden group element and break a three-party Diffie-Hellman
exchange using **one recycled control qubit**. Detect forged broadcasts without
reading, measuring, or inspecting the hidden group register.

Developed by **Naveen and Pulkit**. [Official problem statement](https://qff-website-pi.vercel.app/#hackathon).

## Submitted solution

**[`my_solver.py`](my_solver.py) is the complete standalone solution.** It defines
`solve(oracle)` and returns the required `order` and `secret` fields. A forged
exchange returns `"none"` as the secret. Only this file is needed by the judge;
the rest of the repository supports explanation and local reproduction.

The method combines joint relation-lattice recovery with rational order finding,
prime-power minimization, and adaptive discrete-log equations. Independent relation
tests validate candidates. A non-cyclic subgroup witness can establish a forgery
without completing the entire relation lattice.

## Quick start

From a checkout of this repository, create a virtual environment:

```sh
python -m venv .venv
```

Activate it with `.venv\Scripts\Activate.ps1` on Windows PowerShell, or
`source .venv/bin/activate` on Linux/macOS. Then install and run:

```sh
python -m pip install -e ".[test]"
python -m silentgroup.check my_solver.py --repeat 5
python -m pytest tests -q
```

The recorded experiments used **Python 3.14.6**. For the exact tested package
versions, install `requirements-tested.txt` instead. The checker prints score,
queries, execution time, and correctness for each public instance.

## Validation

| Check | Result |
|---|---:|
| Seeded comparison suite | 184/184 correct |
| Independent-seed suite | 328/328 correct |
| Regression across 21 group shapes | 840/840 correct |
| Official public checker | 40/40 correct |
| Unit tests | 118 passed |

Across **1,392 solver runs**, the maximum observed consumption was **71,094 of
80,000 queries**, and the longest measured run was **8.991 seconds** against a
600-second limit. The implementation uses **one qubit** against the 16-qubit cap.
Finite-shot verification remains probabilistic; these results do not guarantee
performance on the unknown hidden-instance distribution.

The final solver has SHA256
`5d3635eb8439b62a07ef19d70eb1922976ee451e50f1b549570dc94d15c81e85`.

## Guide for reviewers

- [Development walkthrough and presentation guide](docs/PRESENTATION_WALKTHROUGH.md),
  with the problem, original approaches, failure timeline, final design, and the
  logic behind every test. A [formatted reading copy](docs/PRESENTATION_WALKTHROUGH.html)
  and detailed evidence catalogs are included.
- [Algorithm and correctness arguments](docs/ALGORITHM.md), including the forgery
  witness and the one-qubit resource analysis.
- [Reproduction instructions](docs/REPRODUCIBILITY.md), exact environment, and all
  benchmark commands.
- [Recorded results](results/README.md), with raw runs and machine-readable summaries.
- [Attribution](docs/ATTRIBUTION.md), including the organizer-provided simulator.

## Repository layout

```text
my_solver.py               Standalone submission: solve(oracle)
pyproject.toml             Local simulator and testing dependencies
requirements-tested.txt    Exact package versions used in recorded validation
docs/                      Algorithm, walkthrough, test reasoning, and reproduction
silentgroup/               Unmodified organizer simulator and public instances
scripts/                   Reproducible benchmark and regression harnesses
tests/                     Unit tests and failure-mode regressions
results/                   Final validation evidence and comparison records
```

Prompts, scratch scripts, obsolete experiments, and duplicate solver files are
excluded from this submission tree. Development history remains in Git.
