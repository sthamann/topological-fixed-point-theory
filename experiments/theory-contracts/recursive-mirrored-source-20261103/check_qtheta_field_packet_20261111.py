"""Exact tests of moment recovery, conditional on a positive scalar fibre."""
from fractions import Fraction as F


def moments(nodes):
    return tuple(sum(mass * a**k for a, mass in nodes) for k in range(3))


def exact_two_nodes(u, m0, m1, m2):
    assert m0 > 0 and m1**2 <= m0 * m2 <= u**2 * m0**2
    mean, second = m1 / m0, m2 / m0
    if mean == u:
        return ((u, m0),)
    x = (u * mean - second) / (u - mean)
    p = (mean - x) / (u - x)
    assert -u <= x <= mean and 0 <= p <= 1
    return ((x, (1 - p) * m0), (u, p * m0))


def main():
    u = F(13, 100)
    # Use arbitrary rational masses; these are algebra checks, not TFPT weights.
    m, r, eps, b = F(7, 5), u / 2, F(3, 7), u
    even = ((-r, m / 2), (r, m / 2))
    skew = ((-r, m * (1 - eps) / 2), (r, m * (1 + eps) / 2))
    for k in range(0, 14, 2):
        assert sum(w * a**k for a, w in even) == sum(w * a**k for a, w in skew)
    assert moments(even) == (m, 0, m*r*r)
    m0, m1, m2 = moments(skew)
    assert (m0, b*m1, m2/2) == (m, b*m*eps*r, m*r*r/2)

    for j in range(-12, 13):
        for k in range(-12, 13):
            for p in (F(0), F(1, 7), F(1, 2), F(1)):
                target = moments(((u*j/12, m*p), (u*k/12, m*(1-p))))
                assert moments(exact_two_nodes(u, *target)) == target

    m0, m1, m2 = moments(((u/2, m),))
    weight_at_minus_u = (m2/u**2 - m1/u)/2
    assert weight_at_minus_u == -m/8
    assert moments(exact_two_nodes(u, m0, m1, m2)) == (m0, m1, m2)
    print("PASS: parity ambiguity, exact positive reconstruction, fixed-node counterexample")


if __name__ == "__main__":
    main()
