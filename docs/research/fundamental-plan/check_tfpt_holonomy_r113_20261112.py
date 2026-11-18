"""Exact E6/A2 Coxeter-factor census and R113 return discriminators."""
from collections import Counter
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import runpy
import sys

import numpy as np
import sympy as sp


repository = Path(sys.argv[1])
parent = repository / "docs/research/fundamental-plan/check_tfpt_coxeter_cartan_20261112.py"
with redirect_stdout(StringIO()):
    previous = runpy.run_path(str(parent))
roots = np.array(sorted(previous["roots"]), dtype=np.int64)
root_index = {tuple(root): i for i, root in enumerate(roots)}
identity = sp.eye(8)
root_identity = np.arange(len(roots))
coxeter, quarter, family, clock = (previous[key] for key in ("coxeter", "quarter", "family", "combined"))
axis_pairs = [(0, 1), (2, 3), (4, 5), (6, 7)]
basis_roots = [tuple(2 if k == i else 2 * sign if k == j else 0 for k in range(8))
               for i, j in axis_pairs for sign in (1, -1)]
basis_indices = [root_index[root] for root in basis_roots]
basis = roots[basis_indices].T
basis_inverse = sp.Matrix(basis).inv()
variable = sp.Symbol("x")


def permutation(matrix):
    multiplier = int(sp.ilcm(*[x.q for x in matrix]))
    integers = np.array(matrix * multiplier).astype(np.int64)
    images = roots @ integers.T
    assert not (images % multiplier).any()
    return np.array([root_index[tuple(v)] for v in images // multiplier])


def matrix_of(permutation_array):
    return sp.Matrix(roots[permutation_array[basis_indices]].T) * basis_inverse


def trace_powers(permutation_array):
    power = root_identity
    traces = []
    for _ in range(8):
        power = permutation_array[power]
        # The eight chosen orthogonal roots have squared norm eight.
        trace_numerator = np.sum(roots[power[basis_indices]].T * basis)
        assert trace_numerator % 8 == 0
        traces.append(int(trace_numerator // 8))
    return tuple(traces)


def orthogonal_projector(columns):
    return columns * (columns.T * columns).inv() * columns.T


def return_blocks(h, p):
    q = identity - p
    b, complement = q * h * p, q * h * q
    return p * h * p, b.T * b, b.T * complement * b, b.T * complement**2 * b


def coxeter_factor_census():
    reflections = [permutation(identity - sp.Matrix(r) * sp.Matrix(r).T / 4) for r in roots]
    candidates = {}
    for i, j in zip(*np.where(roots @ roots.T == -4)):
        perm = reflections[i][reflections[j]]
        candidates.setdefault(tuple(perm), (int(i), int(j)))
    assert len(candidates) == 2240
    m = permutation(coxeter)
    e6_traces = (1, 3, 4, -1, 1, 0, 1, -1)
    results = {}
    for b_key, pair in candidates.items():
        b = np.array(b_key)
        a = m[np.argsort(b)]
        if trace_powers(a) != e6_traces:
            continue
        a_matrix = matrix_of(a)
        fixed = sp.Matrix.hstack(*(a_matrix - identity).nullspace())
        moving_roots = [i for i, r in enumerate(roots) if fixed.T * sp.Matrix(r) == sp.zeros(2, 1)]
        assert len(moving_roots) == 72 and np.count_nonzero(a == root_identity) == 6
        # A cycle of length below 12 would miss the primitive E6 eigenspace.
        a_cubed = a[a[a]]
        assert all(i != a_cubed[i] for i in moving_roots)
        w = a[b[np.argsort(a)[np.argsort(b)]]]
        w_matrix = matrix_of(w)
        polynomial = sp.factor(w_matrix.charpoly(variable).as_expr())
        w_order = 4 if polynomial == (variable - 1)**4 * (variable**2 + 1)**2 else 3
        assert polynomial == ((variable - 1)**4 * (variable**2 + 1)**2 if w_order == 4
                              else (variable - 1)**6 * (variable**2 + variable + 1))
        assert w_matrix**w_order == identity
        results[b_key] = (w_order, pair, fixed, w_matrix)
    assert Counter(row[0] for row in results.values()) == {4: 5, 3: 15}

    remaining = set(results)
    orbit_lengths = []
    inverse = np.argsort(m)
    while remaining:
        start = remaining.pop()
        current = tuple(m[np.array(start)[inverse]])
        length = 1
        while current != start:
            assert results[current][0] == results[start][0]
            remaining.remove(current)
            length += 1
            current = tuple(m[np.array(current)[inverse]])
        orbit_lengths.append((results[start][0], length))
    assert sorted(orbit_lengths) == [(3, 15), (4, 5)]
    print("2240 A2 Weyl factors; 20 E6 partners. One Coxeter orbit of 15 order-three loops and one of 5 order-four loops.")
    return results


def physical_discriminators(results):
    isotropic = []
    for order, pair, fixed, w in results.values():
        p = orthogonal_projector(fixed)
        b_plane = orthogonal_projector(sp.Matrix(roots[list(pair)].T))
        h = 2 * identity - w - w.T
        a, g, k, l = return_blocks(h, p)
        if order == 4:
            assert h**2 == 2 * h and h.trace() == 8
            assert (a, g, k, l) == (p, p, p, p)
            assert return_blocks(h, b_plane) == (2 * b_plane, sp.zeros(8), sp.zeros(8), sp.zeros(8))
            assert p * b_plane * p == p / 3
            orbit = [clock**j * h * clock**(-j) for j in range(6)]
            assert len({tuple(m) for m in orbit}) == 6
            if sum(orbit, sp.zeros(8)) == 6 * identity:
                isotropic.append((pair, orbit))
        else:
            assert h**2 == 3 * h
            assert a.rank() == 1 and a.trace() == sp.Rational(4, 3)
            assert (g, k, l) == (5*a/3, 25*a/9, 125*a/27)
    assert len(isotropic) == 1
    pair, orbit = isotropic[0]
    assert clock**6 == -identity
    shift = sp.zeros(48)
    for j in range(6):
        row = (j + 1) % 6
        shift[8 * row:8 * row + 8, 8 * j:8 * j + 8] = clock
    registered = sp.diag(*orbit)
    assert registered * shift == shift * registered and shift**6 == -sp.eye(48)
    p1 = np.ones(8, dtype=int)
    p2 = np.array((-1, 1, 1, 1, 1, -1, 1, 1))
    gamma = p1 - p2
    # The D8 pair and its vacuum return are even under the old half parity.
    assert gamma[0] % 2 == 0 and p1[0] % 2 == p2[0] % 2 == 1
    pullback = sp.diag(*(clock**j for j in range(6)))
    assert pullback.T * registered * pullback == sp.diag(*(orbit[0],) * 6)
    print("All 5 quarter loops: PHP=G=K=L=P on Fix(E6); G=0 on the A2 path plane.")
    print("All 15 third-turn loops: rank-one return on Fix(E6), active energy 5/3.")
    print("Exactly one of the 5 has clock mean I in the pinned chart. Doubled A2 roots:", roots[list(pair)].tolist())
    print("Registered transfer is six copies in the comoving frame; the sixth clock power is -I.")
    print("The original D8 pair and neutral vacuum are even; an intertwiner preserving the identified half phase vanishes.")


def original_short_time_reader():
    a, b, c = (sp.Matrix(2, 2, sp.symbols(f"{letter}:4")) for letter in ("a", "b", "c"))
    a, c = a + a.T, c + c.T
    full = a.row_join(b.T).col_join(b.row_join(c))
    moment = [full**j for j in range(1, 5)]
    m1, m2, m3, m4 = [p[:2, :2] for p in moment]
    g = m2 - m1**2
    k = m3 - m1**3 - m1*g - g*m1
    l = m4 - m1**4 - m1**2*g - m1*g*m1 - g*m1**2 - g**2 - m1*k - k*m1
    assert sp.expand(g - b.T*b) == sp.zeros(2)
    assert sp.expand(k - b.T*c*b) == sp.zeros(2)
    assert sp.expand(l - b.T*c*c*b) == sp.zeros(2)
    print("Noncommuting fourth-order compressed-jet formulas for G,K,L checked exactly.")


if __name__ == "__main__":
    physical_discriminators(coxeter_factor_census())
    original_short_time_reader()
