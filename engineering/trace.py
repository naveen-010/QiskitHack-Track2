"""Wrap an Oracle to log every run (kind, shots, cost). Development only."""
import sys, importlib.util
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from silentgroup import Oracle, MockGroup, load_public

class TracingOracle(Oracle):
    def __init__(self, *a, **k):
        super().__init__(*a, **k); self.log = []
    def run(self, qc, shots=1):
        before = self.queries_used
        out = super().run(qc, shots)
        kind = "sample" if qc.num_clbits > 1 else "htest"
        self.log.append((kind, shots, self.queries_used - before, out.get("1", 0) if kind == "htest" else None))
        return out

def summary(log):
    s = sum(c for k, _, c, _ in log if k == "sample"); h = sum(c for k, _, c, _ in log if k == "htest")
    return f"samples={sum(1 for k,*_ in log if k=='sample')} ({s}q) htests={sum(1 for k,*_ in log if k=='htest')} ({h}q)"

if __name__ == "__main__":
    spec = importlib.util.spec_from_file_location("s", sys.argv[1]); s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
    name = sys.argv[2]; seeds = range(int(sys.argv[3]))
    for seed in seeds:
        inst = load_public(name) if name.startswith("public") else MockGroup([int(x) for x in name.split(",")], seed=seed)
        o = TracingOracle(inst, budget=80000, seed=1000 + seed)
        r = s.solve(o)
        print(seed, "OK" if r == inst.answer else f"BAD {r} {inst.answer}", o.queries_used, summary(o.log))
        if len(sys.argv) > 4: print("   ", [(k, sh, c, one) for k, sh, c, one in o.log if k == "htest"])
