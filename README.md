# The Silent Group

*Break a three-party Diffie–Hellman exchange in a group you can never look at.*

## The story

Three members of the Illuminati, known only as **Galileo**, **Bernini** and **Milton**, are planning a secret meeting. The time and place are never written down. Instead, the three of them agree on a shared secret K, and K is the meeting.

They don't trust paper, and they don't trust ordinary arithmetic. Every key exchange the Order runs goes through the *Silent Group*. This is a sealed quantum device that holds one element of a finite abelian group in a quantum register. You hand it a control qubit and an instruction such as "multiply by h to the power k". If the control is |1⟩, the device multiplies its register by h^k. If the control is |0⟩, nothing happens. That's all it does. It has no display and no readout, and it can't be opened. Outside the Order, all anyone knows is that the group is finite and abelian, with at most 2^60 elements.

To set a meeting, the three of them agree on a public element g. Galileo, Bernini and Milton each pick secret exponents a, b and c, and broadcast A = g^a, B = g^b and C = g^c. Then they pass messages around the circle over private channels: Galileo sends B^a to Milton, Bernini sends C^b to Galileo, and Milton sends A^c to Bernini. Each raises what they receive to their own secret, and all three arrive at the same key, K = g^(abc).

You've intercepted the public broadcasts g, A, B and C for three upcoming meetings. The private messages went over channels you couldn't tap. Tonight, for a few hours, you also have the device itself. Find each meeting's K.

You have four things to keep in mind:

- **The device is old.** It raises elements to powers by repeated squaring, so h^k takes about log₂ k internal operations. Every operation drains its power cell.
- **The device is noisy.** Each multiplication leaves a small random error in the phase it kicks back. The bigger the power, the bigger the error.
- **The Order suspects a mole.** Before some meetings, one or more members broadcast a forgery: an element that isn't a power of g at all. No real key comes out of that exchange, and the meeting is a trap. Spot it.
- **The power cell is per meeting.** Every shot of every circuit you run against a meeting's broadcasts drains that meeting's budget.

## Your task

Each meeting gives you four element names: `g`, `A`, `B` and `C`. You can't print group elements, so every answer is a number:

| Field | Meaning |
| --- | --- |
| `order` | ord(g): the smallest n > 0 with g^n = e |
| `secret` | abc mod ord(g), so that K = g^secret |
| `secret = "none"` | if any of A, B or C is not a power of g (a forgery) |

Always report `order`, even for a forgery.

## What you get

This repo is a Python package, `silentgroup`, that simulates the device on your laptop:

| Name | What it does |
| --- | --- |
| `oracle_gate(element, power)` | A one-qubit Qiskit gate. Its qubit is the control, and it multiplies the device's register by `element ** power`. `power` can be any int, including negative or very large values. The register itself is implicit: you never allocate, touch or measure it. It starts in the identity e. |
| `Oracle(instance)` | `.run(circuit, shots)` returns counts. Keys are bitstrings with classical bit 0 rightmost. Also: `.cost(circuit, shots)`, `.budget_left()`, `.queries_used`, `.elements` |
| `load_public("public_1")` | One of two practice meetings at full size. Their answers are in `silentgroup/instances/`. |
| `MockGroup(moduli, ...)` | A group Z_d1 × Z_d2 × … that you choose, for debugging. Its answer is known, and `forged="C"` (or `"AB"`, …) makes those broadcasts forgeries. |

```python
from qiskit import QuantumCircuit
from silentgroup import MockGroup, Oracle, load_public, oracle_gate

qc = QuantumCircuit(1, 1)
qc.h(0)
qc.append(oracle_gate("g", 5), [0])   # controlled multiplication by g^5
qc.h(0)
qc.measure(0, 0)

Oracle(MockGroup([2**20, 6], seed=1), noise=0).run(qc, shots=100)   # perfect device, for debugging
Oracle(load_public("public_1")).run(qc, shots=100)                   # realistic device
```

The source is all here, and you're welcome to read it.

## Rules

1. **Qubits.** At most 16 qubits per circuit. The device's register doesn't count.
2. **Gates.** Anything Qiskit can express on your own qubits: standard and library gates, `measure`, `reset`, and `if_test` conditioned on a classical bit or register. `for_loop`, `while_loop` and `switch_case` aren't supported, so unroll loops in Python. Gates must be fully bound, with no free parameters.
3. **Cost.** Each `oracle_gate(h, k)` costs max(1, bit_length(|k|)) queries per shot. Oracle gates in both branches of an `if_test` are charged. A run costs the circuit's cost per shot times the number of shots.
4. **Noise.** Each oracle call adds an independent Gaussian phase error to the phase it kicks back. Its standard deviation is 0.002 × √(that call's cost) turns, where 1 turn = 2π radians.
5. **Budget.** Each hidden meeting has 80,000 queries. Once they run out, the device stops responding.
6. **Time.** Each hidden meeting gets 10 minutes of wall-clock time.
7. **Adaptivity.** Your solver can look at the results of one run before building the next.
8. **Fair play.** Your solver may interact with the device only through the `oracle` object. An attempt to inspect, import or attack the judge scores 0 for that meeting.

## Scoring

| Component | Weight | How it's scored |
| --- | --- | --- |
| Auto score | 60% | Summed over the 3 hidden meetings, scored per meeting as below |
| Open track | 20% | One extension of your choice, judged on depth |
| Presentation | 20% | 8 minutes: your algorithm, why it works, and what failed |

The auto score for one meeting:

- **0** if `order` is wrong.
- **0.3** if `order` is right but `secret` is wrong.
- **1 + 0.5 × (fraction of budget left)** if both are right.

Ideas for the open track:

- Find the full structure of the subgroup generated by g, A, B and C, as a product Z_n1 × Z_n2 × ….
- Derive a query-cost bound for your method and compare it with what you measure.
- Characterise how your solver degrades as the noise grows.
- Solve with as few qubits as possible.

## Submitting

Submit one Python file that defines:

```python
def solve(oracle) -> dict:
    return {"order": ..., "secret": ...}
```

`examples/solver_template.py` is a starting point. After the deadline the organizers run your `solve` once on each hidden meeting. Your solver runs in a sandbox, and its `oracle` has the same interface as `Oracle` here. It can use the Python standard library, `numpy` and `qiskit`. If you need another package, ask the organizers before the deadline.

To rehearse, score your solver on the public meetings exactly as the judge will score the hidden ones:

```bash
python -m silentgroup.check my_solver.py --repeat 5
```

The noise is random, so run it more than once. A solver that's right 4 times out of 5 can still fail on the night.

## Getting started

```bash
git clone https://github.com/ProAltro/qff-hackathon.git && cd qff-hackathon
python -m venv .venv && source .venv/bin/activate
pip install -e .
python examples/quickstart.py
python -m silentgroup.check examples/solver_template.py
```

## Notes on the physics

The register starts in the identity state |e⟩, and the device never measures it. G is abelian, so all the multiplication operators U_h share one eigenbasis. There is one eigenvector |χ̂⟩ for each character χ of G, with

U_h |χ̂⟩ = e^(2πi χ(h)) |χ̂⟩,

and |e⟩ is the uniform superposition of all these eigenvectors. Work out what that means for what your measurements can tell you. The simulator in `silentgroup/simulator.py` does exactly this calculation.

## FAQ

**Can I see the group?** Only for `MockGroup` and the two public meetings. The hidden groups exist only on the judge's machine.

**Is the hidden group cyclic? How big is it?** Assume nothing beyond "finite, abelian, at most 2^60 elements". Hidden groups vary in shape.

**Why is my answer different each run?** Every shot is random, and so is the noise. Design for that.

**Do failed runs cost queries?** A run that's refused for breaking the rules (too many qubits, an unknown element, an unsupported gate) costs nothing. A run that executes costs its full price, whatever you do with the result.
