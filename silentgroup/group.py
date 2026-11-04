"""Finite abelian groups Z_d1 x ... x Z_dm, their characters, and problem instances.

Holding one of these objects reveals everything about the group. Teams only ever
hold mock and public instances; hidden instances live in the judge process.
"""

from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass
from pathlib import Path

from ._arith import solve_linear_congruences

Element = tuple[int, ...]


@dataclass(frozen=True)
class AbelianGroup:
    moduli: tuple[int, ...]

    def __post_init__(self):
        if not self.moduli or any(d < 1 for d in self.moduli):
            raise ValueError("moduli must be positive integers")

    @property
    def order(self) -> int:
        return math.prod(self.moduli)

    @property
    def exponent(self) -> int:
        return math.lcm(*self.moduli)

    def reduce(self, h) -> Element:
        return tuple(x % d for x, d in zip(h, self.moduli, strict=True))

    def add(self, x: Element, y: Element) -> Element:
        return self.reduce(a + b for a, b in zip(x, y, strict=True))

    def scale(self, k: int, x: Element) -> Element:
        return self.reduce(k * a for a in x)

    def element_order(self, h: Element) -> int:
        return math.lcm(*(d // math.gcd(x, d) for x, d in zip(h, self.moduli)))

    def discrete_log(self, g: Element, h: Element) -> int | None:
        """Smallest x >= 0 with x*g = h, or None if h is not in <g>."""
        sol = solve_linear_congruences(zip(g, h, self.moduli))
        return None if sol is None else sol[0]

    def random_element(self, rng: random.Random) -> Element:
        return tuple(rng.randrange(d) for d in self.moduli)

    def random_character(self, rng: random.Random) -> Element:
        # Characters of Z_d1 x ... x Z_dm are indexed by c in the same group:
        # chi_c(h) = sum_i c_i h_i / d_i  (mod 1).
        return tuple(rng.randrange(d) for d in self.moduli)

    def character_numerator(self, chi: Element, h: Element) -> int:
        """chi(h) as an exact fraction numerator / exponent, numerator in [0, exponent)."""
        D = self.exponent
        return sum(c * x * (D // d) for c, x, d in zip(chi, h, self.moduli)) % D


@dataclass
class Instance:
    name: str
    group: AbelianGroup
    elements: dict[str, Element]
    answer: dict | None = None

    def to_json(self) -> dict:
        return {
            "name": self.name,
            "moduli": list(self.group.moduli),
            "elements": {k: list(v) for k, v in self.elements.items()},
            "answer": self.answer,
        }

    @classmethod
    def from_json(cls, data: dict) -> Instance:
        group = AbelianGroup(tuple(data["moduli"]))
        elements = {k: group.reduce(v) for k, v in data["elements"].items()}
        return cls(data["name"], group, elements, data.get("answer"))


PUBLIC_KEYS = ("A", "B", "C")  # Galileo's, Bernini's and Milton's broadcasts


def compute_answer(group: AbelianGroup, elements: dict[str, Element]) -> dict:
    """ord(g), and abc mod ord(g) (the exponent of K = g^(abc)), or "none" if a broadcast is forged."""
    g = elements["g"]
    n = group.element_order(g)
    logs = [group.discrete_log(g, elements[k]) for k in PUBLIC_KEYS]
    secret = "none" if None in logs else math.prod(logs) % n
    return {"order": n, "secret": secret}


def MockGroup(moduli, *, g=None, a=None, b=None, c=None, forged: str = "", seed=None) -> Instance:
    """A debugging instance in Z_d1 x ... x Z_dm whose answer you know.

    Unspecified pieces are random. ``forged`` names the broadcasts to push outside <g>,
    e.g. "C" or "AB" (this needs a non-cyclic choice of moduli, e.g. [12, 6]).
    """
    rng = random.Random(seed)
    group = AbelianGroup(tuple(int(d) for d in moduli))
    g = group.reduce(g) if g is not None else group.random_element(rng)
    n = group.element_order(g)
    exps = [x if x is not None else rng.randrange(1, max(n, 2)) for x in (a, b, c)]
    elements = {"g": g} | {k: group.scale(x, g) for k, x in zip(PUBLIC_KEYS, exps)}
    for k in forged:
        for _ in range(1000):
            t = group.random_element(rng)
            if group.discrete_log(g, group.add(elements[k], t)) is None:
                elements[k] = group.add(elements[k], t)
                break
        else:
            raise ValueError(f"could not place {k} outside <g>; use non-cyclic moduli")
    return Instance("mock", group, elements, compute_answer(group, elements))


_INSTANCE_DIR = Path(__file__).parent / "instances"


def load_public(name: str) -> Instance:
    """Load a public instance (``"public_1"`` or ``"public_2"``); its answer is included."""
    path = _INSTANCE_DIR / f"{name}.json"
    if not path.exists():
        names = sorted(p.stem for p in _INSTANCE_DIR.glob("*.json"))
        raise KeyError(f"unknown public instance {name!r}; available: {names}")
    return Instance.from_json(json.loads(path.read_text()))
