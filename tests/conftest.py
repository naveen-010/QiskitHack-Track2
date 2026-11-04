"""Shared test helpers: load the submission file as a module (it is not a package)."""

import importlib.util
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
# SOLVER=path/to/file.py lets the same tests run against the baseline or the frozen file.
SOLVER_PATH = Path(os.environ.get("SOLVER", ROOT / "submission" / "my_solver.py"))


def load_solver(path=SOLVER_PATH):
    spec = importlib.util.spec_from_file_location("solver_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="session")
def solver():
    return load_solver()


sys.path.insert(0, str(ROOT))
