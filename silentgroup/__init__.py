"""silentgroup: Diffie-Hellman in a group you never see.

    from silentgroup import Oracle, oracle_gate, MockGroup, load_public

    oracle = Oracle(load_public("public_1"))
    qc = QuantumCircuit(1, 1)
    qc.h(0); qc.append(oracle_gate("g", 1), [0]); qc.h(0); qc.measure(0, 0)
    oracle.run(qc, shots=100)
"""

from .gate import OracleGate, oracle_gate
from .group import Instance, MockGroup, load_public
from .oracle import BUDGET, MAX_QUBITS, NOISE, TIME_LIMIT, BudgetExceeded, Oracle
from .program import InvalidCircuit

__all__ = [
    "BUDGET",
    "MAX_QUBITS",
    "NOISE",
    "TIME_LIMIT",
    "BudgetExceeded",
    "Instance",
    "InvalidCircuit",
    "MockGroup",
    "Oracle",
    "OracleGate",
    "load_public",
    "oracle_gate",
]
