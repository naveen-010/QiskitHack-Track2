"""Shot-by-shot simulation of validated programs.

The group register starts in |e> = |G|^{-1/2} sum_chi |chi^>, and every U_h is
diagonal in the character basis with eigenvalue exp(2 pi i chi(h)). The register is
never measured, so the branches for different chi never interfere. Each shot
therefore samples chi uniformly and turns every oracle gate into a phase gate on
its control qubit. This is exact for any circuit, at any group size.
"""

from __future__ import annotations

import cmath
import math

import numpy as np

from .program import oracle_cost


class _Shot:
    def __init__(self, nq, nc, phase, rng):
        self.nq = nq
        self.psi = np.zeros((2,) * nq, dtype=complex)
        self.psi[(0,) * nq] = 1.0
        self.bits = [0] * nc
        self.phase = phase
        self.rng = rng
        self._slices = [
            tuple(tuple(v if i == q else slice(None) for i in range(nq)) for v in (0, 1)) for q in range(nq)
        ]

    def apply(self, mat, qs):
        k = len(qs)
        targets = list(reversed(qs))  # the matrix's most significant bit is the last qubit
        out = np.tensordot(mat.reshape((2,) * (2 * k)), self.psi, axes=(list(range(k, 2 * k)), targets))
        self.psi = np.moveaxis(out, list(range(k)), targets)

    def measure(self, q) -> int:
        s0, s1 = self._slices[q]
        p1 = min(max(float(np.vdot(self.psi[s1], self.psi[s1]).real), 0.0), 1.0)
        bit = int(self.rng.random() < p1)
        if bit:
            self.psi[s0] = 0
            self.psi /= math.sqrt(p1)
        else:
            self.psi[s1] = 0
            self.psi /= math.sqrt(1 - p1)
        return bit

    def run(self, ops):
        for op in ops:
            kind = op[0]
            if kind == "u":
                self.apply(op[1], op[2])
            elif kind == "oracle":
                _, element, power, q = op
                self.psi[self._slices[q][1]] *= cmath.exp(2j * math.pi * self.phase(element, power))
            elif kind == "measure":
                self.bits[op[2]] = self.measure(op[1])
            elif kind == "reset":
                if self.measure(op[1]):
                    s0, s1 = self._slices[op[1]]
                    self.psi[s0], self.psi[s1] = self.psi[s1].copy(), 0
            elif kind == "if":
                _, cbits, value, true_ops, false_ops = op
                v = sum(self.bits[c] << i for i, c in enumerate(cbits))
                self.run(true_ops if v == value else false_ops)


def simulate(instance, nq, nc, ops, shots, rng, noise: float = 0.0) -> dict[str, int]:
    """Run ``shots`` independent shots; return Qiskit-style counts (clbit 0 rightmost).

    ``noise``: each oracle call's kicked-back phase is off by an independent Gaussian
    with standard deviation noise * sqrt(oracle_cost(power)) turns.
    """
    group, elements = instance.group, instance.elements
    D = group.exponent
    counts: dict[str, int] = {}
    for _ in range(shots):
        chi = group.random_character(rng)
        nums = {name: group.character_numerator(chi, h) for name, h in elements.items()}

        def phase(element, power, nums=nums):
            # Reduce exactly mod 1 before going to float: power can be ~2^140.
            exact = (power * nums[element] % D) / D
            return exact + rng.gauss(0.0, noise * math.sqrt(oracle_cost(power))) if noise else exact

        shot = _Shot(nq, nc, phase, rng)
        shot.run(ops)
        key = "".join(str(b) for b in reversed(shot.bits))
        counts[key] = counts.get(key, 0) + 1
    return counts
