"""Exact controls on actual original insertions. No physical-source promotion."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "verification"))
SOURCE = ROOT / "verification/v252_full_finite_triple.py"
spec = importlib.util.spec_from_file_location("original_v252", SOURCE)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


def exact(array):
    """Reject noninteger entries instead of silently rationalizing floats."""
    def entry(z):
        z = complex(z)
        if not (z.real.is_integer() and z.imag.is_integer()):
            raise ValueError("Constructor control left Gaussian integers")
        return sp.Integer(int(z.real)) + sp.I * sp.Integer(int(z.imag))
    return sp.MutableSparseMatrix(sp.Matrix(np.asarray(array).tolist()).applyfunc(entry))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def equal(a, b, message):
    require(a == b, message)


def comm(a, b):
    return a * b - b * a


def zero(a, message):
    equal(a, sp.zeros(*a.shape), message)


def run():
    # Controls of the real constructor, not measured or derived mass values.
    ys = [np.diag(x).astype(complex) for x in
          ([1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6])]
    mr = np.array([[2, 1, 0], [1, 3, 1], [0, 1, 4]], complex)
    phi = exact(v.full_D(*ys, mr))
    gf = exact(v.grading())
    equal(phi.shape, (96, 96), "actual finite dimension")
    equal(phi, phi.H, "actual D_F adjoint")
    zero(phi * gf + gf * phi, "actual internal oddness")
    equal(gf.trace(), 0, "actual grading balance")
    rank = phi.rank()
    equal(rank, 96, "invertible-input constructor rank")

    z2, z3 = np.zeros((2, 2), complex), np.zeros((3, 3), complex)
    qbasis = [np.eye(2, dtype=complex),
              np.array([[0, 1j], [1j, 0]]),
              np.array([[0, 1], [-1, 0]], complex),
              np.diag([1j, -1j])]
    bf = exact(v.rep(1, z2, z3))
    for q in qbasis:
        af = exact(v.rep(0, q, z3))
        zero(af * bf, "one-form derivative term vanishes")
        ins = af * comm(phi, bf)
        expected = sp.MutableSparseMatrix(96, 96, {})
        # Independently populate the claimed particle LR quaternion/Y block.
        for gen in range(3):
            for colour, up, down, yu, yd in [(0, "nu", "e", ys[3], ys[2])] + [
                    (c, "u", "d", ys[0], ys[1]) for c in range(3)]:
                for row, name in enumerate((up, down)):
                    for col, target in enumerate((up, down)):
                        for gj in range(3):
                            ii = v.IDX[(name, "L", colour, gen)]
                            jj = v.IDX[(target, "R", colour, gj)]
                            y = (yu, yd)[col]
                            expected[ii, jj] = exact([[q[row, col] * y[gen, gj]]])[0, 0]
        equal(ins, expected, "actual quaternion one-form image")
        require(ins != sp.zeros(96), "nonzero actual one-form")

    gamma = sp.diag(1, -1)
    jn = sp.Matrix([[0, -1], [1, 0]])
    s = jn * gamma
    pr = (sp.eye(2) + gamma) / 2
    equal(jn * jn, -sp.eye(2), "normal Clifford square")
    equal(s, sp.Matrix([[0, 1], [1, 0]]), "product factor sign")
    equal(comm(s, pr), jn, "bare-right mixing sign")
    p = sp.Symbol("p", real=True)
    n = sp.Symbol("n", integer=True)
    # Use sparse actual 192-dimensional product, keeping actual gamma_F.
    bp = sp.kronecker_product(gamma, sp.eye(96)) * p - sp.kronecker_product(s, phi)
    pbare = sp.kronecker_product(pr, sp.eye(96))
    equal(comm(bp, pbare), -sp.kronecker_product(jn, phi), "actual bare-right mixing")
    equal(bp.subs(p, p-n)*pbare-pbare*bp,
          -n*pbare-sp.kronecker_product(jn, phi), "all-mode same-time Ward defect")
    total = sp.kronecker_product(gamma, gf)
    zero(comm(bp, total), "actual TOTAL grading preserved")
    rp = sp.kronecker_product(s, gf)
    # Reflection reverses partial_r: test its symbol and both remaining terms.
    jprod = sp.kronecker_product(jn, sp.eye(96))
    tprod = sp.kronecker_product(s, sp.eye(96))
    pot = sp.kronecker_product(gamma, phi)
    equal(rp * jprod * rp, -jprod, "product reflection normal symbol")
    equal(rp * tprod * rp, tprod, "product reflection tangent symbol")
    equal(rp * pot * rp, pot, "product reflection finite potential")

    dims = []
    for masses, expected_dim in [([1, 2], 2), ([1, -1], 4),
                                 ([1, 0], 3), ([1, -1, 0], 6)]:
        ph = sp.diag(*masses)
        size = 2 * len(masses)
        gg = sp.kronecker_product(gamma, sp.eye(len(masses)))
        ss = sp.kronecker_product(s, ph)
        basis = [sp.Matrix(size, size, lambda r, c: int(r == i and c == j))
                 for i in range(size) for j in range(size)]
        columns = [sp.Matrix(list(comm(gg, x)) + list(comm(ss, x))) for x in basis]
        constraints = sp.Matrix.hstack(*columns)
        null = constraints.nullspace()
        equal(len(null), expected_dim, "complete constant linear commutant")
        for vector in null:
            x = sp.Matrix(size, size, list(vector))
            zero(comm(gg, x), "commutant basis grading")
            zero(comm(ss, x), "commutant basis finite coupling")
        dims.append({"masses": masses, "commutant_complex_dimension": len(null)})

    # A genuinely nondiagonal reducing projection, including opposite masses.
    ph = sp.diag(1, -1)
    q = sp.Matrix([[9, 12], [12, 16]]) / 25
    el = ph * q * ph
    e = sp.diag(q, el)
    b = sp.kronecker_product(gamma, sp.eye(2)) * p - sp.kronecker_product(s, ph)
    equal(e * e, e, "nondiagonal orthoprojection")
    equal(e.H, e, "nondiagonal selfadjoint projection")
    zero(comm(b, e), "nondiagonal projection time reduction")
    nprime = sp.kronecker_product(jn, ph)
    equal(nprime * b * nprime.H, -b, "massive complement pairing")
    zero(comm(e, nprime), "massive projection respects pairing")
    equal(q.rank() - el.rank(), 0, "massive chiral balance")
    ph0 = sp.diag(1, 0)
    e0 = sp.kronecker_product(pr, sp.diag(0, 1))
    b0 = sp.kronecker_product(gamma, sp.eye(2)) * p - sp.kronecker_product(s, ph0)
    zero(comm(b0, e0), "unpaired actual kernel channel")
    equal(e0.rank(), 1, "unit kernel level witness")

    # The pre-existing normal E2 example has index1 for every k, but its
    # whole scalar current has a single right Ward only at k=1.
    kernel_examples = []
    for k in range(1, 5):
        g0 = sp.diag(*([1]*k + [-1]*(k-1)))
        equal(g0.trace(), 1, "old E2 normal kernel net level")
        single_right_ward = (-n*g0 == -n*sp.eye(2*k-1))
        equal(single_right_ward, k == 1, "net index does not force time Ward")
        kernel_examples.append({"positive": k, "negative": k-1,
                                "net_level": 1, "single_right_time_ward": single_right_ward})

    # Singular finite odd examples retain zero TOTAL graded kernel index.
    gtest = sp.diag(1, 1, -1, -1)
    for cross in (sp.diag(1, 0), sp.zeros(2)):
        odd = sp.BlockMatrix([[sp.zeros(2), cross], [cross.H, sp.zeros(2)]]).as_explicit()
        kernel = sp.Matrix.hstack(*odd.nullspace())
        equal((kernel.H*gtest*kernel).trace(), 0, "balanced finite TOTAL kernel index")

    c = -2 * sp.I
    va = sp.Matrix([[1, 1], [1, -1]])
    k = va * sp.diag(c, -c) * va.T
    equal(k, sp.Matrix([[0, 2*c], [2*c, 0]]), "mixed VA cocycles")
    right = sp.Matrix([[sp.Rational(1, 2), sp.Rational(1, 2)]])
    equal((right * k * right.T)[0, 0], c, "right current cocycle recovered")
    return {
        "verdict": "EXACT_CONTROLS_PASS_WITH_ANALYTIC_CLASSIFICATION_IN_PROOF",
        "original_source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "finite_constructor_control": {"dimension": 96, "rank": rank,
             "positive_internal_grades": 48, "negative_internal_grades": 48,
             "quaternion_insertion_witnesses": 4, "input_status": "integer controls, not physical Yukawas"},
        "constant_commutants": dims,
        "existing_normal_kernel_controls": kernel_examples,
        "mixed_va_cocycle_at_original_profiles": str(2*c),
        "right_cocycle_at_original_profiles": str(c),
        "actual_product_grading_and_reflection_symbols": "PASS",
        "state": "retained pure two-face CAR graph; right face thermal",
        "actual_physical_P1_source_image_derived": False,
        "full_TFPT_solution": False,
    }


if __name__ == "__main__":
    print(json.dumps(run(), ensure_ascii=False, indent=2))
