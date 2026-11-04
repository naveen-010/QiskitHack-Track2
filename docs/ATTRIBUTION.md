# Attribution and source provenance

The solution was developed by **Naveen and Pulkit** for Qiskit Fall Fest Track 2,
The Silent Group.

- Pulkit's prototype supplies joint relation-lattice recovery, exact-rational LLL,
  integer triangular bases, noise-calibrated relation tests, and outlier handling.
- Naveen's prototype supplies direct rational order recovery, prime-power order
  minimization, modular discrete-log equations, and adaptive measurements.
- The final implementation combines these methods with bounded repair attempts,
  verified non-cyclic forgery witnesses, identity handling, and recovery from
  proper-divisor order samples.

The organizer-provided [`silentgroup/`](../silentgroup/) simulator and public
instances are retained so judges can reproduce the solution locally. Their source
is the [official participant package](https://github.com/ProAltro/qff-hackathon),
commit `89265613ec267045e5f36be745fb2fd4d029dd10`. The simulator is not part of our
submitted solver and its source files have not been modified.

Original prototype commits in the development repository:

- Naveen: `82e17e1be64df3d24956e60bc756f09f3f4da428`, `solver.py`.
- Pulkit: `a067e0ad3eec13904e92e784e5ff61580b4fdaa7`, `submission/my_solver.py`.

The final solver's validation and provenance are recorded in
[`results/combined_validation.json`](../results/combined_validation.json).
Only `my_solver.py` is passed to the oracle interface. Test and benchmark programs
use public/mock answers exclusively to evaluate results.
