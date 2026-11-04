"""The oracle gate: controlled multiplication of the sealed group register by h^k."""

from qiskit.circuit import Gate


class OracleGate(Gate):
    """Controlled multiplication by ``element ** power`` on the sealed group register.

    The gate acts on a single control qubit. The group register it multiplies is
    implicit: it starts in |e>, only oracle gates touch it, and it is never measured.
    ``power`` may be any Python int, including negative and very large values.
    """

    def __init__(self, element: str, power: int = 1):
        if not isinstance(element, str):
            raise TypeError("element must be a name such as 'g', 'A' or 'B'")
        power = int(power)
        super().__init__("oracle", 1, [], label=f"{element}^k")
        self.element = element
        self.power = power


def oracle_gate(element: str, power: int = 1) -> OracleGate:
    """Return an :class:`OracleGate`; append it with ``qc.append(oracle_gate("g", k), [q])``."""
    return OracleGate(element, power)
