"""F1: with a perfect oracle (noise=0) the solver must return the exact answer on tiny groups."""

import pytest

from silentgroup import MockGroup, Oracle

CASES = [
    ([8], ""), ([15], ""), ([8, 3], ""), ([4, 8], ""), ([7, 11], ""),
    ([7], ""), ([64], ""), ([1001], ""), ([65537], ""),
    # forged variants need non-cyclic moduli so that some element lies outside <g>
    ([12, 6], "C"), ([4, 8], "AB"), ([6, 6, 10], "ABC"), ([2, 4, 8], "B"), ([7, 7], "A"),
]


@pytest.mark.parametrize("moduli,forged", CASES, ids=[f"{m}{f}" for m, f in CASES])
@pytest.mark.parametrize("seed", range(5))
def test_exact(solver, moduli, forged, seed):
    inst = MockGroup(moduli, forged=forged, seed=seed)
    orc = Oracle(inst, budget=80_000, seed=100 + seed, noise=0)
    res = solver.solve(orc)
    assert res == inst.answer
    assert isinstance(res["order"], int)
    assert res["secret"] == "none" or isinstance(res["secret"], int)
    assert orc.queries_used <= 80_000
