"""The oracle handle teams program against, locally and in the judge."""

from __future__ import annotations

import random

from qiskit.circuit import QuantumCircuit

from .group import Instance
from .program import InvalidCircuit, compile_circuit, count_queries, validate
from .simulator import simulate

MAX_QUBITS = 16
NOISE = 0.002  # turns of phase jitter per sqrt(query); the competition value
BUDGET = 80_000  # queries per hidden instance
TIME_LIMIT = 600  # seconds of wall-clock per hidden instance


class BudgetExceeded(RuntimeError):
    pass


class Oracle:
    """Runs circuits against one instance and charges a query budget.

    One ``oracle_gate(h, k)`` costs bit_length(|k|) queries (min 1) per shot, counting
    both branches of every if. Its kicked-back phase carries Gaussian noise with standard
    deviation ``noise * sqrt(cost)`` turns. ``budget=None`` means unlimited; pass
    ``noise=0`` to debug against a perfect oracle.
    """

    def __init__(
        self,
        instance: Instance,
        budget: int | None = None,
        seed=None,
        noise: float = NOISE,
        max_qubits: int = MAX_QUBITS,
    ):
        self._instance = instance
        self._budget = budget
        # Unseeded runs (the judge) use OS randomness: measurement outcomes leak the generator's
        # outputs, and a recoverable Mersenne Twister would reveal which character each shot drew.
        self._rng = random.Random(seed) if seed is not None else random.SystemRandom()
        self.noise = noise
        self.max_qubits = max_qubits
        self.queries_used = 0

    @property
    def elements(self) -> list[str]:
        return sorted(self._instance.elements)

    def budget_left(self) -> int | None:
        return None if self._budget is None else self._budget - self.queries_used

    def cost(self, circuit: QuantumCircuit, shots: int = 1) -> int:
        _, _, ops = validate(compile_circuit(circuit), self.elements, self.max_qubits)
        return count_queries(ops) * shots

    def run(self, circuit: QuantumCircuit, shots: int = 1) -> dict[str, int]:
        """Simulate ``circuit``; returns counts keyed by bitstring, clbit 0 rightmost."""
        return self.run_program(compile_circuit(circuit), shots)

    def run_program(self, program: dict, shots: int = 1) -> dict[str, int]:
        if not isinstance(shots, int) or shots < 1:
            raise InvalidCircuit("shots must be a positive int")
        nq, nc, ops = validate(program, self.elements, self.max_qubits)
        cost = count_queries(ops) * shots
        left = self.budget_left()
        if left is not None and cost > left:
            raise BudgetExceeded(f"run costs {cost} queries but only {left} remain")
        self.queries_used += cost
        return simulate(self._instance, nq, nc, ops, shots, self._rng, self.noise)
