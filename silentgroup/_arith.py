"""Integer helpers used to build instances."""

import math


def solve_linear_congruences(equations):
    """Solve a*x = b (mod n) simultaneously for every (a, b, n).

    Returns (x0, m) such that the solutions are exactly x = x0 (mod m), or None.
    """
    x0, m = 0, 1
    for a, b, n in equations:
        # Substitute x = x0 + m*t and solve (a*m) t = b - a*x0 (mod n) for t.
        A, B = a * m % n, (b - a * x0) % n
        d = math.gcd(A, n)
        if B % d:
            return None
        nd = n // d
        t0 = (B // d) * pow(A // d, -1, nd) % nd if nd > 1 else 0
        x0, m = x0 + m * t0, m * nd
        x0 %= m
    return x0, m
