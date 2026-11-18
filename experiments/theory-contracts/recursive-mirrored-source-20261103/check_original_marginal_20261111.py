"""Test a coordinate marginal on the unchanged QTheta diagnostic chart."""
from collections import defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path
import sys

import numpy as np


def add(left, right):
    out = left.copy()
    for key, (real, imag) in right.items():
        a, b = out.get(key, (0, 0))
        out[key] = a + real, b + imag
    return {key: value for key, value in out.items() if any(value)}


def matmul(left, right):
    rows = defaultdict(list)
    for (row, col), value in right.items():
        rows[row].append((col, value))
    out = {}
    for (row, mid), (a, b) in left.items():
        for col, (c, d) in rows[mid]:
            real, imag = out.get((row, col), (0, 0))
            out[row, col] = real + a*c - b*d, imag + a*d + b*c
    return {key: value for key, value in out.items() if any(value)}


def inner(left, right):
    return sum(a*c + b*d for key in left.keys() & right.keys()
               for a, b in (left[key],) for c, d in (right[key],))


def exact_lift(*matrices):
    """Lift all supplied binary matrix entries to the same dyadic lattice."""
    dyadic = []
    scale = 1
    for matrix in matrices:
        entries = {}
        for i, j in zip(*np.nonzero(matrix)):
            value = matrix[i, j]
            re = Fraction(float(value.real))
            im = Fraction(float(value.imag))
            scale = max(scale, re.denominator, im.denominator)
            entries[i, j] = (re, im)
        dyadic.append(entries)
    return scale, [{k: (int(re*scale), int(im*scale)) for k, (re, im) in m.items()}
                   for m in dyadic]


def main(repo):
    pkg = repo / ("experiments/theory-contracts/original-paper-source-reconciliation-20261102"
                  "/contribution_20261106Q")
    sys.path.insert(0, str(pkg / "runtime"))
    from assemble import Assembler

    native = Assembler(pkg / "domain/native-diagnostic-v1.json")
    profiles = {(x, y): native.assemble({"genes": [x+1, 2, y+1]})["h_N"]
                for x, y in product((-1, 0, 1), repeat=2)}
    h = profiles[0, 0]
    b = (profiles[1, 0] - profiles[-1, 0]) / 2
    c = (profiles[0, 1] - profiles[0, -1]) / 2
    for (x, y), value in profiles.items():
        np.testing.assert_allclose(value, h+x*b+y*c, atol=1e-14, rtol=0)
    for matrix in (h, b, c):
        np.testing.assert_array_equal(matrix, matrix.conj().T)

    scale, (h_exact, b_exact, c_exact) = exact_lift(h, b, c)
    assert inner(b_exact, c_exact) == inner(h_exact, c_exact) == 0
    c_squared_norm = inner(c_exact, c_exact)
    assert c_squared_norm > 0

    exact_second = {}
    exact_fourth = {}
    for x, y in profiles:
        bx = {k: (x*re, x*im) for k, (re, im) in b_exact.items()}
        cy = {k: (y*re, y*im) for k, (re, im) in c_exact.items()}
        matrix = add(add(h_exact, bx), cy)
        square = matmul(matrix, matrix)
        exact_second[x, y] = inner(matrix, matrix)
        exact_fourth[x, y] = inner(square, square)

    for x, y in profiles:
        assert exact_second[x, y] - exact_second[x, 0] == y*y*c_squared_norm
    fourth_increment = {(x, y): value - exact_fourth[x, 0]
                        for (x, y), value in exact_fourth.items() if y}
    # Both hidden signs at each outer point lie strictly above both
    # hidden signs at the center, independent of any relative weights.
    c_min = Fraction(
        min(value for (x, _), value in fourth_increment.items() if x)
        - max(value for (x, _), value in fourth_increment.items() if not x),
        scale**4,
    )
    assert c_min > Fraction(24, 100)

    # Each H has norm < 12 by the maximum absolute row sum, bounded
    # using |re+i im| <= |re|+|im| on the exact dyadic entries.
    rows = defaultdict(int)
    for matrix in (h_exact, b_exact, c_exact):
        for (i, _), (re, im) in matrix.items():
            rows[i] += abs(re) + abs(im)
    bound = 1 + max(rows.values()) // scale
    d = h.shape[0]
    # The contrast has total absolute coefficient at most 4.
    # For exp(-t H^2), the remainder past t^2 is <= d (t ||H||^2)^3/6.
    # This gives a rigorous sign interval for the *action* contrast.
    time_upper = 3*c_min / (4*d*bound**6)
    print("size", d, "rowsum bound", bound)
    print("common second increment", float(Fraction(c_squared_norm, scale**2)))
    print("exact lower fourth separation", float(c_min))
    print("certified open interval", float(time_upper))

    eigs = {key: np.linalg.eigvalsh(value) for key, value in profiles.items()}
    fourth = {k: np.sum(e**4) for k, e in eigs.items()}
    def contrast(values, y):
        return (values[1, y]+values[-1, y])/2-values[0, y]
    for key, value in fourth_increment.items():
        x, y = key
        np.testing.assert_allclose(fourth[key]-fourth[x, 0],
                                   float(Fraction(value, scale**4)),
                                   atol=1e-9, rtol=0)
    for t in (0.001, 0.02, 1.0):
        action = {key: np.sum(np.expm1(-t*e**2)) for key, e in eigs.items()}
        quotients = {(x, y): np.expm1(-(action[x, y]-action[x, 0])) for x, y in profiles}
        gaps = [contrast(quotients, y) for y in (-1, 1)]
        print("time", t, "quotient contrast", gaps)
        assert all(g < 0 for g in gaps) if t < 1 else all(g > 0 for g in gaps)


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
