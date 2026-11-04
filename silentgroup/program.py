"""Lower a Qiskit circuit to a small JSON-safe program, and validate programs.

The judge receives programs, never pickles, so a submitted solver cannot run code
in the process that holds the hidden group. In-process runs use the same path.

Program format::

    {"num_qubits": n, "num_clbits": m, "ops": [op, ...]}
    op = ["u", [[re, im], ...], [q, ...]]      unitary, row-major, Qiskit qubit order
       | ["oracle", element, power, q]
       | ["measure", q, c]
       | ["reset", q]
       | ["if", [c, ...], value, [op, ...], [op, ...]]   value read little-endian from the bits
"""

from __future__ import annotations

import numpy as np
from qiskit.circuit import ClassicalRegister, Clbit, IfElseOp, QuantumCircuit
from qiskit.quantum_info import Operator

from .gate import OracleGate

MAX_GATE_QUBITS = 4
MAX_DEPTH = 8


class InvalidCircuit(ValueError):
    pass


def compile_circuit(circuit: QuantumCircuit) -> dict:
    qmap = {q: i for i, q in enumerate(circuit.qubits)}
    cmap = {c: i for i, c in enumerate(circuit.clbits)}
    return {
        "num_qubits": circuit.num_qubits,
        "num_clbits": circuit.num_clbits,
        "ops": _lower(circuit, qmap, cmap),
    }


def _lower(circuit, qmap, cmap):
    ops = []
    for inst in circuit.data:
        op = inst.operation
        qs = [qmap[q] for q in inst.qubits]
        cs = [cmap[c] for c in inst.clbits]
        if isinstance(op, OracleGate):
            ops.append(["oracle", op.element, op.power, qs[0]])
        elif op.name == "measure":
            ops.append(["measure", qs[0], cs[0]])
        elif op.name == "reset":
            ops.append(["reset", qs[0]])
        elif op.name in ("barrier", "delay"):
            continue
        elif isinstance(op, IfElseOp):
            target, value = op.condition
            if isinstance(target, Clbit):
                bits = [cmap[target]]
            elif isinstance(target, ClassicalRegister):
                bits = [cmap[b] for b in target]
            else:
                raise InvalidCircuit("if_test conditions must be (clbit, value) or (register, value)")
            branches = []
            for body in (op.blocks + (None,))[:2]:
                if body is None:
                    branches.append([])
                    continue
                inner_q = dict(zip(body.qubits, qs))
                inner_c = dict(zip(body.clbits, cs))
                branches.append(_lower(body, inner_q, inner_c))
            ops.append(["if", bits, int(value), *branches])
        elif op.name == "oracle":
            raise InvalidCircuit("use silentgroup.oracle_gate to build oracle gates")
        elif (mat := _matrix(op, len(qs))) is not None:
            ops.append(["u", [[float(z.real), float(z.imag)] for z in mat.ravel()], qs])
        elif getattr(op, "definition", None) is not None:
            # Library and user-defined composite gates (QFTGate, qc.to_gate(), ...).
            body = op.definition
            ops.extend(_lower(body, dict(zip(body.qubits, qs)), dict(zip(body.clbits, cs))))
        else:
            raise InvalidCircuit(f"unsupported operation {op.name!r}")
    return ops


def _matrix(op, num_qubits):
    if getattr(op, "num_clbits", 0) or num_qubits > MAX_GATE_QUBITS:
        return None
    try:
        return np.asarray(op.to_matrix() if hasattr(op, "to_matrix") else Operator(op).data)
    except Exception:
        return None


def validate(program: dict, elements, max_qubits: int):
    """Check a program and return (num_qubits, num_clbits, ops) with numpy matrices.

    Raises InvalidCircuit on anything malformed. Never trusts the input.
    """
    try:
        nq, nc, ops = int(program["num_qubits"]), int(program["num_clbits"]), program["ops"]
    except (KeyError, TypeError, ValueError) as exc:
        raise InvalidCircuit("malformed program") from exc
    if not 1 <= nq <= max_qubits:
        raise InvalidCircuit(f"circuit uses {nq} qubits; the cap is {max_qubits}")
    if not 0 <= nc <= 100_000:
        raise InvalidCircuit("too many classical bits")
    elements = set(elements)

    def qubit(q):
        if not isinstance(q, int) or not 0 <= q < nq:
            raise InvalidCircuit(f"bad qubit index {q!r}")
        return q

    def clbit(c):
        if not isinstance(c, int) or not 0 <= c < nc:
            raise InvalidCircuit(f"bad clbit index {c!r}")
        return c

    def walk(ops, depth):
        if depth > MAX_DEPTH or not isinstance(ops, list):
            raise InvalidCircuit("control flow nested too deeply or malformed")
        out = []
        for op in ops:
            kind = op[0] if isinstance(op, list) and op else None
            if kind == "u" and len(op) == 3:
                qs = [qubit(q) for q in op[2]]
                if not qs or len(set(qs)) != len(qs) or len(qs) > MAX_GATE_QUBITS:
                    raise InvalidCircuit("bad gate qubits")
                dim = 2 ** len(qs)
                mat = np.array([complex(re, im) for re, im in op[1]], dtype=complex)
                if mat.size != dim * dim:
                    raise InvalidCircuit("gate matrix has the wrong size")
                mat = mat.reshape(dim, dim)
                if not np.allclose(mat @ mat.conj().T, np.eye(dim), atol=1e-8):
                    raise InvalidCircuit("gate matrix is not unitary")
                out.append(("u", mat, tuple(qs)))
            elif kind == "oracle" and len(op) == 4:
                if op[1] not in elements:
                    raise InvalidCircuit(f"unknown element {op[1]!r}; this instance has {sorted(elements)}")
                if not isinstance(op[2], int):
                    raise InvalidCircuit("oracle power must be an int")
                out.append(("oracle", op[1], op[2], qubit(op[3])))
            elif kind == "measure" and len(op) == 3:
                out.append(("measure", qubit(op[1]), clbit(op[2])))
            elif kind == "reset" and len(op) == 2:
                out.append(("reset", qubit(op[1])))
            elif kind == "if" and len(op) == 5:
                bits = tuple(clbit(c) for c in op[1])
                out.append(("if", bits, int(op[2]), walk(op[3], depth + 1), walk(op[4], depth + 1)))
            else:
                raise InvalidCircuit(f"malformed op {str(op)[:60]}")
        return tuple(out)

    try:
        return nq, nc, walk(ops, 0)
    except InvalidCircuit:
        raise
    except (TypeError, ValueError, IndexError) as exc:
        raise InvalidCircuit(f"malformed program: {exc}") from exc


def oracle_cost(power: int) -> int:
    """Queries charged for one oracle gate: h^k takes bit_length(|k|) group operations."""
    return max(1, abs(power).bit_length())


def count_queries(ops) -> int:
    """Query cost of a validated program for one shot, counting both branches of every if."""
    total = 0
    for op in ops:
        if op[0] == "oracle":
            total += oracle_cost(op[2])
        elif op[0] == "if":
            total += count_queries(op[3]) + count_queries(op[4])
    return total
