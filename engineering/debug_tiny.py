"""Diagnostics for the tiny-group failure: replay solve() step by step and print lattice data."""
import math, sys, importlib.util
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from silentgroup import MockGroup, Oracle
spec = importlib.util.spec_from_file_location("s", sys.argv[1] if len(sys.argv) > 1 else ROOT / "my_solver.py")
s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
mod = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [1000003]
inst = MockGroup(mod, seed=seed)
print("answer", inst.answer, "elements", inst.elements)
orc = Oracle(inst, budget=80000, seed=100 + seed)
b = 36
qc = s.build_pe_circuit(b)
G = inst.group; D = G.exponent
samples = []
for it in range(12):
    key = next(iter(orc.run(qc, shots=1)))
    samples.append(s.decode(key, b))
    if len(samples) < 4: continue
    t = len(samples)
    rows = [[1 if j == c else 0 for j in range(4)] + [x[c] for x in samples] for c in range(4)]
    rows += [[0]*4 + [(1 << b) if j == i else 0 for j in range(t)] for i in range(t)]
    red = s.lll(rows)
    tol = 1 << (b - 12)
    print(f"--- t={t} used={orc.queries_used}")
    for r in red:
        tail = [((v + (1 << (b-1))) % (1 << b)) - (1 << (b-1)) for v in r[4:]]
        x = r[:4]
        true_rel = G.reduce([0]*len(mod)) == G.reduce([sum(xi*e[k] for xi, e in zip(x, [inst.elements[n] for n in s.ELS])) for k in range(len(mod))])
        print(f"   x={x}  maxtail=2^{math.log2(max(1,max(abs(v) for v in tail))):.1f} good={max(abs(v) for v in tail)<=tol and any(x)} TRUE={true_rel}")
    res = s.consensus(samples, b)
    if res[0] is not None:
        H, short = res[0]
        n, logs = s.extract(H)
        print("   HNF diag", [H[i][i] for i in range(4)], "det", math.prod(H[i][i] for i in range(4)), "n", n, "logs", logs, "used", len(res[1]))
