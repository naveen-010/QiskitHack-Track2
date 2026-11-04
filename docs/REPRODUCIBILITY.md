# Reproduce the submission

Run commands from the repository root after setting up the environment described
in the [README](../README.md). The standalone entry point is `my_solver.py`.

## Environment

The recorded runs used Windows, Python 3.14.6, Qiskit 2.6.1, NumPy 2.5.4, and
pytest 9.2.0. Install these package versions with:

```sh
python -m pip install -r requirements-tested.txt
```

For another supported Python version, `python -m pip install -e ".[test]"` lets
pip resolve compatible versions. The recorded measurements refer to the tested
environment above.

## Verify the artifact

The solver has SHA256
`5d3635eb8439b62a07ef19d70eb1922976ee451e50f1b549570dc94d15c81e85`.
Its code is unchanged from the 1,392-run validation. Git keeps this file in LF
format to preserve its checksum across operating systems.

```sh
python -c "import hashlib; from pathlib import Path; print(hashlib.sha256(Path('my_solver.py').read_bytes()).hexdigest())"
```

## Tests and official checker

```sh
python -m pytest tests -q
python -m silentgroup.check my_solver.py --repeat 20
```

The checker uses fresh measurement randomness, so exact query counts will vary.
The unit suite includes small exact groups, oracle semantics, integer-lattice
invariants, identity order, incomplete certificates, and repair-budget regressions.

## Seeded benchmarks

```sh
# The 184-run comparison suite:
python scripts/benchmark.py --tag combined_final

# The independent 328-run suite:
python scripts/benchmark.py --tag combined_final_fresh --seed0 8 --seeds 16 --public-repeats 20 --offset 20000

# The 840-run regression suite (21 shapes, 40 seeds):
python scripts/regression.py --seeds 40 --tag combined_final_stress --workers 4
```

These commands recreate their CSV/JSON/text outputs under `results/`. The two
benchmark suites use noise 0.002 and an 80,000-query budget. Their oracle seeds are
`offset + instance_seed`; the regression harness uses `100 + instance_seed`.
Both harnesses accept `--solver path/to/file.py` to evaluate another implementation.

The benchmark's worker count defaults to four. Reduce `--workers` on a small
machine. Recorded elapsed times include concurrent execution.

## Interpreting the evidence

See the [results index](../results/README.md) for raw records and summary metrics.
The selected cases cover multiple group shapes and failure modes, but are not a
statistical model of the hidden judge distribution. Finite-shot verification and
phase reconstruction remain probabilistic.
