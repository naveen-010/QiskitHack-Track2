"""How close does a CLEAN sweep's convergent sit to theta? That licenses ORDER_SLACK."""
import math
from fractions import Fraction
from silentgroup import MockGroup, Oracle
import solver

for label, moduli in (("2^59", [2**59]), ("prime-2^56", [72057594037927931]), ("2^35ish", [37199287125])):
    inst = MockGroup(moduli, seed=7)
    n = inst.group.element_order(inst.elements["g"])
    for B in (136, 144):
        clean, resid, cascaded = 0, [], 0
        for _ in range(60):
            orc = Oracle(inst, budget=None)
            m = solver._ipe(orc, ["g"], B)
            theta = Fraction(m[0], 1 << B)
            approx = theta.limit_denominator(min(solver.MAX_ORDER, 1 << ((B - 24) // 2)))
            if n % approx.denominator == 0 and approx.denominator > 1:
                clean += 1
                d = abs(theta - approx)
                resid.append(math.log2(float(d)) + B if d else -99)
            else:
                cascaded += 1
        worst = max(resid) if resid else None
        print(f"{label:11s} B={B}: clean {clean}/60  cascaded {cascaded}  "
              f"worst clean residual = 2^-(B-{-worst:.1f})" if worst is not None
              else f"{label} B={B}: no clean sweeps")
