"""Quickstart: phase kickback from the sealed group register.

Each shot samples one character chi of the hidden group, and oracle_gate(h, k) kicks
the phase exp(2 pi i k chi(h)) back onto its control qubit. Here we read chi(g) to
8 bits with textbook phase estimation on a small mock group.
"""

from fractions import Fraction

from qiskit import QuantumCircuit
from qiskit.circuit.library import QFTGate

from silentgroup import MockGroup, Oracle, load_public, oracle_gate

t = 8
qc = QuantumCircuit(t, t)
qc.h(range(t))
for j in range(t):
    qc.append(oracle_gate("g", 2**j), [j])
qc.append(QFTGate(t).inverse(), range(t))
qc.measure(range(t), range(t))

mock = MockGroup([256], g=[3], seed=0)  # Z_256, g = 3 has order 256
oracle = Oracle(mock)
counts = oracle.run(qc, shots=10)
print("mock answer:", mock.answer)
for key, c in counts.items():
    print(f"  chi(g) = {Fraction(int(key, 2), 2**t)}  x{c}")
print("queries used:", oracle.queries_used)

# The public instances are 2^40 - 2^60 in size: 8 bits won't cut it, and you only have 16 qubits.
public = Oracle(load_public("public_1"))
print("public_1 elements:", public.elements, "- answer is published in silentgroup/instances/public_1.json")
