"""For a cyclic mock group, recover each shot's character by brute force and print per-element phase errors."""
import sys, importlib.util, math
import numpy as np
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from silentgroup import MockGroup, Oracle
spec = importlib.util.spec_from_file_location("s", sys.argv[1]); s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
seed = int(sys.argv[2]); N = int(sys.argv[3]) if len(sys.argv) > 3 else 1000003
nshots = int(sys.argv[4]) if len(sys.argv) > 4 else 12
inst = MockGroup([N], seed=seed); orc = Oracle(inst, budget=None, seed=100 + seed)
b = 36; qc = s.build_pe_circuit(b)
el = np.array([inst.elements[e][0] for e in s.ELS], dtype=object)
cs = np.arange(N, dtype=np.int64)
for i in range(nshots):
    X = s.decode(next(iter(orc.run(qc, shots=1))), b)
    # error for each candidate character c, element e: X_e - 2^b * (c*e mod N)/N, wrapped
    tot = np.zeros(N)
    errs = []
    for k in range(4):
        ph = (cs * int(el[k]) % N) / N
        d = (X[k] / 2**b - ph + 0.5) % 1 - 0.5
        tot += np.abs(d)
        errs.append(d)
    c = int(np.argmin(tot))
    e = [round(errs[k][c] * 2**b) for k in range(4)]
    print(f"shot {i}: c={c} err(units of 2^-36) = " + " ".join(f"{s.ELS[k]}:{e[k]:+d}(2^{math.log2(abs(e[k])+1):.1f})" for k in range(4)))
