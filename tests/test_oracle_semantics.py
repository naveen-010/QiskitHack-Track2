"""F0: simulator truth table. Each test checks one fact the solver's design depends on."""

import math

from qiskit import QuantumCircuit

from silentgroup import MockGroup, Oracle, oracle_gate


def hadamard_block(qc, power, clbit, element="g"):
    """H - controlled g^power - H - measure - reset on qubit 0."""
    qc.h(0)
    qc.append(oracle_gate(element, power), [0])
    qc.h(0)
    qc.measure(0, clbit)
    qc.reset(0)


def test_one_character_per_shot():
    """(a) In Z_2 the phase chi(g) is 0 or 1/2, so a Hadamard test is deterministic for a fixed
    character. Two tests in ONE circuit always agree; two separate shots agree only ~50%."""
    inst = MockGroup([2], g=[1], a=1, b=1, c=1, seed=0)
    orc = Oracle(inst, seed=1, noise=0)
    qc = QuantumCircuit(1, 2)
    hadamard_block(qc, 1, 0)
    hadamard_block(qc, 1, 1)
    counts = orc.run(qc, shots=400)
    assert set(counts) <= {"00", "11"}, counts             # perfectly correlated inside a shot
    assert min(counts.get("00", 0), counts.get("11", 0)) > 120  # and both characters occur

    single = QuantumCircuit(1, 1)
    hadamard_block(single, 1, 0)
    pairs = [(next(iter(orc.run(single))), next(iter(orc.run(single)))) for _ in range(400)]
    disagree = sum(x != y for x, y in pairs)
    assert 120 < disagree < 280                             # independent characters across shots


def test_measure_reset_if_feedback():
    """(b) measure -> reset -> if_test: copy the first outcome into a second clbit via feedback."""
    inst = MockGroup([2], g=[1], a=1, b=1, c=1, seed=0)
    orc = Oracle(inst, seed=2, noise=0)
    qc = QuantumCircuit(1, 2)
    hadamard_block(qc, 1, 0)                                # random bit c0, qubit reset to |0>
    with qc.if_test((qc.clbits[0], 1)):
        qc.x(0)                                             # flip only if c0 == 1
    qc.measure(0, 1)
    counts = orc.run(qc, shots=300)
    assert set(counts) <= {"00", "11"}, counts
    assert len(counts) == 2


def test_cost_formula_both_branches():
    """(c) cost = shots * sum max(1, bit_length(|k|)), with BOTH if/else branches charged."""
    inst = MockGroup([1009], seed=0)
    orc = Oracle(inst, seed=3, noise=0)
    powers = [1, 5, -1024, 0, 2**40 + 1]
    qc = QuantumCircuit(1, 1)
    for k in powers:
        qc.append(oracle_gate("g", k), [0])
    qc.measure(0, 0)
    with qc.if_test((qc.clbits[0], 1)) as else_:
        qc.append(oracle_gate("A", 7), [0])                 # cost 3
    with else_:
        qc.append(oracle_gate("B", 255), [0])               # cost 8
    expected = sum(max(1, abs(k).bit_length()) for k in powers) + 3 + 8
    assert orc.cost(qc, 1) == expected
    assert orc.cost(qc, 7) == 7 * expected
    before = orc.queries_used
    orc.run(qc, shots=7)
    assert orc.queries_used - before == 7 * expected


def test_identity_vs_nonidentity():
    """(d) g^ord(g) = e gives P(1) ~ 0 (noise only); a random non-identity element gives ~0.5."""
    inst = MockGroup([1000003], seed=5)
    n = inst.answer["order"]
    orc = Oracle(inst, seed=4)                              # competition noise 0.002
    qc_id = QuantumCircuit(1, 1)
    qc_id.h(0)
    qc_id.append(oracle_gate("g", n), [0])
    qc_id.h(0)
    qc_id.measure(0, 0)
    ones_id = orc.run(qc_id, shots=400).get("1", 0)
    # P(1) = sin^2(pi * noise) with sigma = 0.002*sqrt(20) ~ 0.009 turns -> ~1e-3
    assert ones_id <= 4

    qc_g = QuantumCircuit(1, 1)
    qc_g.h(0)
    qc_g.append(oracle_gate("g", 1), [0])
    qc_g.h(0)
    qc_g.measure(0, 0)
    ones_g = orc.run(qc_g, shots=400).get("1", 0)
    assert abs(ones_g / 400 - 0.5) < 0.1


def test_noise_sigma_scales_with_sqrt_cost():
    """sigma = noise * sqrt(cost): an identity relation costing ~60 has P(1) ~ (pi*sigma)^2."""
    inst = MockGroup([1000003], seed=5)
    n = inst.answer["order"]
    orc = Oracle(inst, seed=9, noise=0.02)                  # large noise so the effect is visible
    qc = QuantumCircuit(1, 1)
    qc.h(0)
    qc.append(oracle_gate("g", n * (1 << 40)), [0])         # identity, cost = 61
    qc.h(0)
    qc.measure(0, 0)
    shots = 4000
    p = orc.run(qc, shots=shots).get("1", 0) / shots
    sigma = 0.02 * math.sqrt(61)
    pred = (1 - math.exp(-2 * (math.pi * sigma) ** 2)) / 2  # E[sin^2(pi*N(0,sigma))]
    assert abs(p - pred) < 4 * math.sqrt(pred * (1 - pred) / shots) + 0.005
