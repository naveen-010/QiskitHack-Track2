"""Regressions for the comparison failures and the combined solver's evidence handling."""
import pytest
from pathlib import Path
import hashlib
import json

from silentgroup import MockGroup, Oracle, load_public


def test_submission_matches_frozen_validation_hash():
    root = Path(__file__).resolve().parents[1]
    validation = json.loads((root / "results/combined_validation.json").read_text())
    assert hashlib.sha256((root / "my_solver.py").read_bytes()).hexdigest() == validation["source_sha256"]


@pytest.mark.parametrize("forged", ["", "A"])
def test_identity_order_with_noisy_phase_seed(solver, forged):
    inst = MockGroup([1009], g=(0,), forged=forged, seed=0)
    oracle = Oracle(inst, budget=80000, seed=10000)
    assert solver.solve(oracle) == inst.answer
    assert oracle.queries_used < 1000


def test_minimization_can_reach_one(solver):
    oracle = Oracle(MockGroup([1009], g=(0,), seed=0), budget=80000, seed=3, noise=0)
    ctx = solver._Ctx(oracle)
    assert solver.minimize_order(ctx, 2**32) == 1


def test_rational_proper_divisor_completion(solver, monkeypatch):
    inst = MockGroup([2**30], g=(1,), seed=0)
    oracle = Oracle(inst, budget=80000, seed=0, noise=0)
    # Emulate a clean character whose phase denominator is only half the order.
    monkeypatch.setattr(solver, "decode", lambda *args: [2**88 // 2**29])
    assert solver.cheap_order(solver._Ctx(oracle), bits=88, bound=2**32) == 2**30
    assert oracle.queries_used < 8000


def test_capped_torsion_search_is_not_a_certificate(solver):
    p = 1009
    inst = MockGroup([p], g=(1,), a=1, b=2, c=3, seed=0)
    oracle = Oracle(inst, budget=80000, seed=12345, noise=0)
    ctx = solver._Ctx(oracle)
    H = solver.hnf([[p, 0, 0, 0], [0, 1, 0, -2], [0, 0, 1, -3], [0, 0, 0, p]])
    assert solver.verify_generators(ctx, solver.lll(H))
    assert solver.torsion_candidates(solver.lll(H), p) is None
    assert solver.complete(ctx, H) is None
    result = solver.finish_lattice(ctx, H)
    assert result is None or result == inst.answer
    assert inst.answer == {"order": p, "secret": 6}


@pytest.mark.parametrize("budget", [80000, None])
def test_completion_preserves_budget_for_another_attempt(solver, budget):
    inst = MockGroup([2], g=(1,), a=0, b=0, c=0, seed=0)
    oracle = Oracle(inst, budget=budget, seed=1, noise=0)
    ctx = solver._Ctx(oracle)
    H = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 2**30]]
    solver.complete(ctx, H)
    assert oracle.queries_used <= solver.COMPLETE_BUDGET
    if budget is not None:
        assert oracle.budget_left() >= 80000 - solver.COMPLETE_BUDGET
    assert ctx.best is not None
    assert solver.det_of(ctx.best) < solver.det_of(H)


@pytest.mark.parametrize("p", [2, 3, 7])
def test_noncyclic_witness_survives_in_actual_group(solver, p):
    inst = MockGroup([p, p], g=(1, 0), a=1, b=1, c=1, seed=0)
    # Set the public forgery explicitly in this test, never inside the solver.
    inst.elements["A"] = (0, 1)
    ctx = solver._Ctx(Oracle(inst, budget=80000, seed=12, noise=0))
    H = [[p, 0, 0, 0], [0, 1, 0, -1], [0, 0, 1, -1], [0, 0, 0, p]]
    assert solver.verify_generators(ctx, solver.lll(H))
    assert solver.noncyclic_witness(ctx, H)


@pytest.mark.parametrize("p", [2, 3, 7])
def test_incomplete_lattice_cannot_fake_noncyclic_witness(solver, p):
    inst = MockGroup([p], g=(1,), a=1, b=1, c=1, seed=0)
    ctx = solver._Ctx(Oracle(inst, budget=80000, seed=12, noise=0))
    H = [[p, 0, 0, 0], [0, 1, 0, -1], [0, 0, 1, -1], [0, 0, 0, p]]
    assert solver.verify_generators(ctx, solver.lll(H))
    assert not solver.noncyclic_witness(ctx, H)


def test_insufficient_log_budget_is_unresolved(solver):
    inst = MockGroup([2**59], g=(1,), seed=0)
    oracle = Oracle(inst, budget=100, seed=0)
    assert solver.direct_logs(solver._Ctx(oracle), inst.answer["order"]) is None
    assert oracle.queries_used == 0


@pytest.mark.parametrize("moduli,kwargs,seed", [
    ([2**59], {}, 1),
    ([2**55], {"g": (1,), "a": 2, "b": 3, "c": 5}, 1),
    ([2**59], {"g": (1,), "a": 0}, 1),
    ([2**60], {"g": (1,)}, 1),
    ([2**60], {"g": (1,)}, 6),
    ([1009] * 4, {"g": (1, 0, 0, 0), "forged": "ABC"}, 3),
])
def test_previous_failures_and_large_torsion(solver, moduli, kwargs, seed):
    inst = MockGroup(moduli, seed=seed, **kwargs)
    oracle = Oracle(inst, budget=80000, seed=10000 + seed)
    assert solver.solve(oracle) == inst.answer
    assert oracle.queries_used <= 80000


def test_public2_previously_exhausted_seed(solver):
    inst = load_public("public_2")
    oracle = Oracle(inst, budget=80000, seed=10014)
    assert solver.solve(oracle) == inst.answer
    assert oracle.queries_used < 75000
