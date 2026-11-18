"""Exact fixed-source clock obstruction and charge-compensated Ramond return."""
import ast
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import runpy

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
original = runpy.run_path(str(HERE / "check_ramond_full_return_20261110.py"))
native = runpy.run_path(str(HERE / "check_fable_slot_clock_20261110.py"))
M, q, V, Q, clean = (original[name] for name in ("M", "q", "V", "Q", "cl"))
unit = sp.eye(8)
c = sp.Matrix(native["c"])
R = M.T * c * M
r = M.T * (unit[:, 4] + unit[:, 5])
q_prime = R * q
gamma = q - q_prime
y = sp.Matrix(native["y6"]) / 6
center = sp.Matrix([1, -1, 0, -1, -1, -1, 0, -1])
alphas = (sp.Rational(1, 4), sp.Rational(3, 4))
assert gamma.dot(gamma) == 4 and gamma.dot(q_prime) == -2
assert M * gamma == 2 * unit[:, 7]
assert sorted(map(abs, gamma)) == [0] * 4 + [1] * 4
assert y.dot(M * gamma) == 0 and center.dot(gamma) % 3 == 0


def drift(alpha):
    return sp.Matrix([(i + alpha) / 8 - sp.Rational(1, 2) for i in range(8)])


def spectral_and_color():
    for alpha in alphas:
        d = drift(alpha)
        q_shifts = [d.dot(R**k * q - q) for k in range(12)]
        r_shifts = [d.dot(R**k * r - r) for k in range(12)]
        common = [k for k, (a, b) in enumerate(zip(q_shifts, r_shifts))
                  if a % 1 == b % 1 == 0]
        assert common == [0] and q_shifts[1] == sp.Rational(1, 2)
        assert d.dot(gamma) == -sp.Rational(1, 2)
        print("quarter", alpha, "q and r K shifts", q_shifts, r_shifts, flush=True)
        for moved_alpha in alphas:
            moved = drift(moved_alpha)
            simultaneous = [k for k in range(12)
                            if all((moved.dot(R**k * p) - d.dot(p)) % 1 == 0
                                   for p in (q, r))]
            assert simultaneous == ([0] if moved_alpha == alpha else [])
    assert all(center.dot(R**k * q) % 3 == 0 for k in range(12))
    excluded = [k for k in range(12) if center.dot(R**k * r) % 3]
    assert excluded == [2, 4, 5, 8, 10, 11]
    assert M * R * r == unit[:, 1] - unit[:, 0]
    assert y.dot(M * r) == sp.Rational(1, 2) and y.dot(M * R * r) == 0
    print("r color-center excluded powers", excluded, flush=True)


def cover_support():
    # R072 U_alpha is a diagonal arc phase times this Fourier matrix.
    zeta = (1 + sp.I) / sp.sqrt(2)
    F = sp.Matrix(8, 8, lambda a, k: zeta**((a*k) % 8) / sp.sqrt(8))
    positive = [i for i in range(8) if gamma[i] == 1]
    negative = [i for i in range(8) if gamma[i] == -1]
    assert positive == [1, 5] and negative == [3, 7]
    pairs = list(combinations(range(8), 2))

    def minors(columns):
        return {p: sp.simplify(F.extract(p, columns).det()) for p in pairs}

    creation, annihilation = minors(positive), minors(negative)
    for values in (creation, annihilation):
        assert all(sp.simplify(value * sp.conjugate(value)) ==
                   (sp.Rational(1,16) if (a-b) % 2 else 0)
                   for (a,b), value in values.items())
    coefficient = creation[0,1] * sp.conjugate(annihilation[2,3])
    assert sp.simplify(coefficient * sp.conjugate(coefficient)) == sp.Rational(1,256)
    assert all(any(a in b and creation[b] != 0 and annihilation[d] != 0
                   and set(b).isdisjoint(d) for b in pairs for d in pairs)
               for a in range(8))
    print("folded compensator: cover coefficient on disjoint arcs (0,1;2,3) has modulus 1/16",
          flush=True)


def cocycle_phase():
    path = ROOT / "experiments/theory-contracts/seam-joint-ward-state-20261014/operator_transport_check.py"
    nodes = [n for n in ast.parse(path.read_text()).body if
             (isinstance(n, ast.FunctionDef) and n.name in {"eps_car", "split_even_lattice", "eps_e8"})
             or (isinstance(n, ast.Assign) and any(getattr(t, "id", "") in {"e", "basis"} for t in n.targets))]
    env = {"S": sp, "q": sp.ones(8, 1), "s": sp.ones(8, 1) / 2}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), env)
    eps, basis = env["eps_e8"], env["basis"]
    coefficients = basis.inv() * q

    def phase(frame):
        point, result = sp.zeros(8, 1), 1
        for b, k in zip((frame[:, j] for j in range(8)), coefficients, strict=True):
            assert k.q == 1
            inverse_sign = eps(b, -b) if k < 0 else 1
            direction = -b if k < 0 else b
            for _ in range(abs(k)):
                result *= inverse_sign * eps(point, direction)
                point += direction
        return result

    # Generator convention E_b -> E_(Rb), reused from the original J lift.
    return phase(basis) * phase(R * basis) * eps(gamma, q_prime)


def gram_comparison():
    q_axes = sp.Matrix.hstack(V, q)
    axes = sp.Matrix.hstack(q_axes, R * q_axes, gamma)
    index_gamma = axes.cols - 1
    pair_metric = clean(axes.H * axes)
    metrics = [[Q.from_sympy(x) for x in row] for row in pair_metric.tolist()]
    grams = original["inner"].__globals__["inner_q"]
    grams[:] = [metrics, metrics]
    original["contractions"].cache_clear()
    dot = [Q.from_sympy(x) for x in clean(axes.T * gamma)]
    phase = Q.from_sympy(cocycle_phase())
    print("chosen ordered-lattice phase", Q.to_sympy(phase), flush=True)

    def dressed(poly):
        # On input momentum q_prime the z^(-2) coefficient preserves oscillator grade.
        # This is the original lattice vertex annihilation/creation expansion.
        result = defaultdict(lambda: Q.zero)
        for word, val in poly.items():
            word = tuple((n, i + q_axes.cols) for n, i in word)
            for count in range(len(word) + 1):
                for removed in combinations(range(len(word)), count):
                    value = phase * val
                    for j in removed:
                        value *= -dot[word[j][1]]
                    if not value:
                        continue
                    kept = tuple(w for j, w in enumerate(word) if j not in removed)
                    order = sum(word[j][0] for j in removed)
                    for added, coeff in original["exp_y"](order):
                        key = tuple(sorted(kept + tuple((n, index_gamma) for n, _ in added)))
                        result[key] += value * coeff
        return {k: v for k, v in result.items() if v}

    polynomials, _ = original["m"].code()
    # Use the five unchanged preparations and their full q vertex outputs.
    inner = original["inner"]

    def pairing(left, right):
        return clean(sp.Matrix(len(left), len(right), lambda i, j: inner(left[i], right[j])))

    def rank(matrix):
        return matrix.to_DM().convert_to(Q).rank()

    for alpha in alphas:
        d = [Q.from_sympy(x) for x in clean(V.T * drift(alpha))]
        states = [original["vertex_zero"](original["series_response"](p, d, 4))
                  for p in polynomials]
        images = [dressed(v) for v in states]
        initial = pairing(states, states)
        dressed_cross = pairing(states, images)
        dressed_norm = pairing(images, images)
        shifted = [{tuple((n, i + q_axes.cols) for n, i in word): a
                    for word, a in state.items()} for state in states]
        assert pairing(shifted, shifted) == initial
        assert all(sum(n for n, _ in word) == 3 for image in images for word in image)
        assert rank(initial) == 5
        assert dressed_cross * initial[0, 0] != initial * dressed_cross[0, 0]
        assert dressed_norm * initial[0, 0] != initial * dressed_norm[0, 0]
        print("quarter", alpha, "rank undressed Gram / dressed cross / dressed norm",
              rank(initial), rank(dressed_cross), rank(dressed_norm), flush=True)
        print("packet 0 and 1 dressed cross", dressed_cross[:2, :2].tolist(), flush=True)
        print("packet 0 and 1 original Gram", initial[:2, :2].tolist(), flush=True)
        print("individual dressed/initial squared norm ratios",
              [sp.simplify(dressed_norm[i, i] / initial[i, i]) for i in range(5)], flush=True)


if __name__ == "__main__":
    spectral_and_color()
    cover_support()
    gram_comparison()
