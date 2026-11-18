"""Exact checks for the Coxeter lift and the two family actions."""
from collections import Counter
from contextlib import redirect_stdout
from itertools import combinations_with_replacement
from pathlib import Path
import io
import json
import runpy
import sys

import sympy as sp


repository = Path(sys.argv[1])
plan = repository / "docs/research/fundamental-plan"
with redirect_stdout(io.StringIO()):
    original = runpy.run_path(str(plan / "check_fable_origin_rule_20261112.py"))
roots = original["roots"]
dot = original["dot"]
glue = original["glue_reader"]
certificate = json.loads((repository / "experiments/theory-contracts/compiler-integral-triality-20260918/certificate.json").read_text())
clocks = {k: sp.Matrix(v["vector"]).applyfunc(sp.Rational) for k, v in certificate["clock_matrices"].items()}
identity = sp.eye(8)
coxeter, quarter, family, combined = (clocks[k] for k in ("C", "J", "sigma", "c"))
z = sp.Symbol("z")


def image(matrix, root):
    return tuple(matrix * sp.Matrix(root))


def root_orbits(matrix):
    remaining = set(roots)
    lengths = Counter()
    while remaining:
        start = remaining.pop()
        current = image(matrix, start)
        length = 1
        while current != start:
            remaining.remove(current)
            current = image(matrix, current)
            length += 1
        lengths[length] += 1
    return lengths


def phase_table(root_set, delta24, sigma24, spin_offset=0):
    # Roots use doubled coordinates; the torus vectors use denominator 12.
    counts = Counter(((dot(delta24, p) + spin_offset * (p[0] % 2)) % 24,
                      dot(sigma24, p) % 24) for p in root_set)
    counts[0, 0] += 8
    return tuple(tuple(counts[d, s] for s in (0, 8, 16)) for d in (0, 6, 12, 18))


def seifert_factor(size):
    return sp.Matrix(
        size, size, lambda i, j: int(i == j) - int(j == i + 1))


def lift_check():
    seifert = sp.kronecker_product(seifert_factor(2), seifert_factor(4))
    monodromy = -seifert.inv() * seifert.T
    metric = seifert + seifert.T
    assert metric.det() == 1 and monodromy.T * metric * monodromy == metric
    assert coxeter.T * coxeter == identity
    for operator in (monodromy, coxeter):
        assert operator.charpoly(z).as_expr() == sp.cyclotomic_poly(30, z)
        assert operator**15 == -identity and (identity - operator).det() == 1
        assert sum((operator**j for j in range(15)), sp.zeros(8)) == 2 * (identity - operator).inv()
        gold = operator**7 - operator**3 - operator**2
        assert gold**2 + gold == identity
    assert root_orbits(coxeter) == {30: 8}
    assert quarter**2 == -identity and family**3 == identity
    assert combined == quarter * family == family * quarter
    assert combined**4 == family and combined**9 == quarter
    assert {image(quarter, p) for p in roots} == roots
    assert {image(family, p) for p in roots} == roots
    print("Original Seifert M and native C: Phi30, M^15=-I, det(1-M)=1; 8 root orbits of length 30.")
    print("Torusal norm of 15 steps is 2*(1-M)^(-1); parity-matched lifts form one orbit of size 2^8.")


def marked_half_period_check():
    p1 = (1,) * 8
    p2 = (-1, 1, 1, 1, 1, -1, 1, 1)
    middle = tuple(x - y for x, y in zip(p1, p2))
    spin_phase = lambda p: (-1) ** (p[0] % 2)
    assert (spin_phase(p2), spin_phase(middle), spin_phase(p1)) == (-1, 1, -1)
    assert all((2 * dot(glue, p) / 2 - sp.Rational(p[0] % 2, 2)) % 1 == 0 for p in roots)
    surviving = []
    for before, after in ((0, 0), (1, 0), (0, 1), (1, 1)):
        charge = tuple((-1)**after * ((-1)**before * a + b) - c
                       for a, b, c in zip(p2, middle, p1))
        surviving.append(not any(charge))
    assert surviving == [True, False, False, False]
    print("At fixed original end marks, old-Cartan momentum inversion kills every inserted word.")
    print("The new-Cartan lift gives original relative factors +1,-1,-1,+1 on the four words.")


def family_check():
    assert root_orbits(family) == {1: 12, 3: 76}
    assert len((family - identity).nullspace()) == 4
    # In the fixed-S vector model D cycles each of four Gaussian groups.
    # Across the four groups sigma has eigenvalues 1,1,omega,omega^-1
    # at each D-eigenvalue. The two torus rows below select one line
    # from each conjugate pair in that 16-dimensional vector action.
    d = (0, 0, 6, 6, 3, 3, 3, 3)
    sigma = (0, 4, 0, 4, 0, 0, 4, -4)
    target = ((22, 19, 19), (24, 20, 20), (22, 19, 19), (24, 20, 20))
    assert phase_table(roots, d, sigma) == target
    for vector_lift_sign in (0, 6):
        shifted = tuple(a + vector_lift_sign for a in d)
        for spin_shift in (0, 12):
            assert phase_table(roots, shifted, sigma, spin_shift) == target
    assert sum(row[0] for row in target) == 92
    # Classify order-three factors in the original Spin(10) x SU(4)
    # cover. The common central kernel has order 4.
    matches = []
    pure_family = None
    for exponents in combinations_with_replacement(range(3), 4):
        if sum(exponents) % 3:
            continue
        a, b, c, d = exponents
        d -= sum(exponents)
        wedge = [a + b, a + c, a + d]
        for nonzero in range(6):
            for final_sign in ((1, -1) if nonzero == 5 else (1,)):
                carrier = [int(i < nonzero) for i in range(5)]
                carrier[-1] *= final_sign
                torus = carrier + wedge
                if sum(torus) % 2:
                    torus[0] += 3
                old_d = tuple(int(a * 12) for a in glue)
                table = phase_table(roots, old_d, [4 * a for a in torus])
                if table == target:
                    matches.append((nonzero, exponents))
                if nonzero == 0 and exponents == (0, 0, 1, 2):
                    pure_family = table
    assert matches == [(3, (1, 1, 2, 2))]
    assert sum(row[0] for row in pure_family) == 134
    first, second = (1, 2, 2, 1), (2, 1, 2, 1)
    assert Counter(first) == Counter(second) == {1: 2, 2: 2}
    assert tuple((a + b) % 3 for a, b in zip(first, second)) == (0, 0, 1, 2)
    # The original exp(4 pi i Y) on the full spinor lattice is the
    # carrier class with three nontrivial defining vector planes.
    y6 = (-2, -2, -2, 3, 3, 0, 0, 0)
    carrier = (4, 1, 1, 0, 0, 0, 0, 0)
    assert all((dot(y6, p) - dot(carrier, p)) % 6 == 0 for p in roots)
    p1, p2 = (1,) * 8, (-1, 1, 1, 1, 1, -1, 1, 1)
    gamma = tuple(a - b for a, b in zip(p1, p2))
    assert tuple((dot(y6, p) % 6) // 2 for p in (p1, p2, gamma)) == (0, 2, 1)
    print("The D5 companions can be the original opposite color-center elements; a separated middle insertion has phase omega^(-1).")
    print("Native sigma lift joint (D,sigma) multiplicities:", target)
    print("Its unique original factor class rotates three D5 planes and has SU(4) multiplicities 2+2.")
    print("Pure meridian 3-cycle fixes 134 directions, native lift 92; two transported native-class factors give the pure 3-cycle.")


def transported_family_orbit_check():
    phases = [k for k in range(30) if coxeter**k * quarter == quarter * coxeter**k]
    orientations = [k for k in range(15)
                    if (coxeter**k * family * coxeter**(-k)) * quarter
                    == quarter * (coxeter**k * family * coxeter**(-k))]
    assert phases == [0, 15] and orientations == [0]
    print("In the original pinned Coxeter orbit only its initial sigma orientation returns to the fixed quarter.")


def time_moment_check():
    m = tuple(map(sp.Rational, ("1/2", "5/8", "7/8", "21/16", "131/64", "209/64")))
    hankel = sp.Matrix(3, 3, lambda i, j: m[i+j])
    coefficients = -hankel.inv() * sp.Matrix(m[3:6])
    recurrence = z**3 + sum(coefficients[i] * z**i for i in range(3))
    assert hankel.det() == sp.Rational(1, 4096)
    assert sp.factor(recurrence) == (2*z - 3)*(4*z**2 - 10*z + 5)/8
    energies = (sp.Rational(3, 2), (5-sp.sqrt(5))/4, (5+sp.sqrt(5))/4)
    weights = (sp.Rational(1, 4), (5+sp.sqrt(5))/40, (5-sp.sqrt(5))/40)
    assert all(sp.simplify(sum(w * e**i for e, w in zip(energies, weights))) == v
               for i, v in enumerate(m))
    print("Reported return moments imply three energies:", energies, "with weights", weights)


def transported_hamiltonian_check():
    r = sp.eye(8)
    r.row_swap(0, 1)  # original E8 root reflection in e1-e2
    polys = []
    for m in (coxeter, r * coxeter * r):
        source_h = identity - (m + m.inv()) / 2
        h0 = sum((combined**j * source_h * combined**(-j) for j in range(12)), sp.zeros(8)) / 12
        assert h0.T == h0 and h0 * quarter == quarter * h0 and h0 * family == family * h0
        assert all(h0[:i, :i].det() > 0 for i in range(1, 9))
        polys.append(sp.factor(h0.charpoly(z).as_expr()))
        for n in (0, 1, 14, 15):
            un = m**n
            hn = un * h0 * un.inv()
            assert m * hn * m.inv() == m**(n + 1) * h0 * m**(-n - 1)
    assert polys[0] != polys[1]
    print("The same native clock and parity allow two distinct positive averaged spectra:", *polys)


if __name__ == "__main__":
    lift_check()
    marked_half_period_check()
    family_check()
    transported_family_orbit_check()
    time_moment_check()
    transported_hamiltonian_check()
