"""Replay the charge-conjugation pairing on the original full96 spatial assembler."""

from itertools import product
from fractions import Fraction
from math import comb
import runpy
import sys

import numpy as np
import sympy as sp

from reconstruction import require
from r113_paired_response import native_jet
from source import NATIVE


ASSEMBLER = NATIVE.parent.parent/"runtime/assemble.py"


def load_assembler(root):
    path = root/ASSEMBLER
    original_path = sys.path[:]
    try:
        sys.path.insert(0, str(path.parent))
        return runpy.run_path(str(path))["Assembler"](root/NATIVE)
    finally:
        sys.path[:] = original_path


def fixed_form_spectral_certificate(couplings, momenta):
    """A finite-field certificate for exact nonresonance of reversed modes."""
    prime = 65537

    def residue(value):
        rational = Fraction(value)
        return rational.numerator*pow(rational.denominator, -1, prime) % prime

    def decode(name, conjugate=False):
        sign = -1 if conjugate else 1
        return sp.Matrix([[(residue(real)+sign*256*residue(imag)) % prime for real, imag in row]
                          for row in couplings[name]])

    species = ("Yu", "Yd", "Ye", "Ynu")
    determinants = [int(decode(key).det()) % prime for key in species]
    x = sp.Symbol("x")

    def charpoly(matrix):
        return sp.Poly((matrix.applyfunc(lambda value: int(value) % prime)).charpoly(x).as_expr(),
                       x, modulus=prime)

    total = sp.Poly(1, x, modulus=prime)
    for key in species[:-1]:
        total *= charpoly(decode(key, conjugate=True).T * decode(key))
    y, ybar = decode("Ynu"), decode("Ynu", conjugate=True)
    m, mbar = decode("MR"), decode("MR", conjugate=True)
    zero = sp.zeros(3)
    b = sp.BlockMatrix([[y, zero], [mbar.T, y.T]]).as_explicit()
    bbar = sp.BlockMatrix([[ybar, zero], [m.T, ybar.T]]).as_explicit()
    total *= charpoly(bbar.T*b)
    differences = sorted({abs(p*p-q*q) for p, q in zip(momenta, reversed(momenta))})
    gaps = {str(gap): total.gcd(total.compose(sp.Poly(x+residue(gap), x, modulus=prime))).degree()
            for gap in differences}
    return {"prime": prime, "dirac_determinants": determinants, "reversed_mode_squared_gcd_degrees": gaps}


def controls(root):
    assembler = load_assembler(root)
    fixture = assembler.domain
    choices = fixture["field_proposal"]["allowed_multipliers"]
    slots = fixture["field_proposal"]["slots"]
    count = len(fixture["modes"])
    gamma, swap = assembler.gamma, assembler.K
    internal = gamma@swap
    pairing = np.kron(np.eye(count), internal)
    dimension = len(pairing)
    checks = []

    def check(name, condition):
        require(condition, name)
        checks.append(name)

    def assemble(genes):
        return assembler.assemble({"schema": "tfpt.field-configuration.v1",
                                   "domain_sha256": fixture["_sha256"], "genes": list(genes)})["h_N"]

    def reflected(genes):
        return tuple(choices.index(str(-int(choices[g]))) if slot["kind"] == "sin" else g
                     for g, slot in zip(genes, slots))

    check("v252 grading and real structure give a skew unitary quaternionic form",
          np.array_equal(internal.T, -internal) and np.array_equal(internal@internal, -np.eye(96)))
    check("the declared domain explicitly uses a spatial and not OS coordinate",
          fixture["coordinate_role"] == "spatial_tangential_circle_NOT_OS_time")
    check("the signed permutation is a full rank skew form on the Fourier cut",
          np.count_nonzero(pairing) == dimension and np.array_equal(pairing.T@pairing, np.eye(dimension)))

    matrices = {g: assemble(g) for g in product(range(len(choices)), repeat=len(slots))}
    failures = []
    for genes, h in matrices.items():
        other = matrices[reflected(genes)]
        check(f"all entries of the reflected infinitesimal pairing cancel for {genes}",
              np.array_equal(h@pairing, -pairing@other.T))
        defect = float(np.max(np.abs(h@pairing+pairing@h.T)))
        if defect:
            failures.append((genes, defect))
        check(f"same profile preserves the pairing precisely on reflection-even genes {genes}",
              (defect == 0) == (genes == reflected(genes)))

    # Connected internal gauge generators computed from the actual real
    # bimodule lift rather than treating the left representation alone as gauge.
    zero2, zero3 = np.zeros((2, 2), complex), np.zeros((3, 3), complex)
    pauli = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))
    variations = [(1j, zero2, zero3)] + [(0, 1j*axis, zero3) for axis in pauli]
    for i in range(3):
        for j in range(i+1, 3):
            e = zero3.copy()
            e[i, j] = 1
            variations.extend(((0, zero2, e-e.T), (0, zero2, 1j*(e+e.T))))
    variations.extend((0, zero2, 1j*np.diag(axis)) for axis in ((1, -1, 0), (0, 1, -1)))
    for i, legs in enumerate(variations):
        left = assembler.original.rep(*legs)
        generator = left+swap@left.conjugate()@swap
        check(f"connected v252 gauge generator {i} preserves the internal pairing",
              np.array_equal(generator@internal, -internal@generator.T))

    cos_slot = next(slot for slot in slots if slot["kind"] == "cos" and slot["harmonic"] == 1)
    sin_slot = next(slot for slot in slots if slot["kind"] == "sin" and slot["harmonic"] == 1)
    check("actual cosine and sine field directions are the first two quaternion axes",
          (cos_slot["component"], sin_slot["component"]) == (0, 1)
          and tuple(map(Fraction, fixture["field_proposal"]["base_constant"])) == (0, 0, 1, 0))
    left = assembler.original.rep(1, np.diag([1j, -1j]), np.eye(3))
    align = left@swap@left.conjugate()@swap
    check("actual weak gauge rotation aligns those quaternion axes",
          np.array_equal(align@assembler.V[0]@align.conj().T, assembler.V[1]))
    certificate = fixed_form_spectral_certificate(fixture["couplings"],
        [mode+Fraction(fixture["quarter"]) for mode in fixture["modes"]])
    check("exact finite field certificate excludes reversed Fourier block resonances",
          all(certificate["dirac_determinants"]) and certificate["reversed_mode_squared_gcd_degrees"]
          == {"1/2": 0, "3/2": 0})
    grading = np.kron(np.eye(count), gamma)
    check("the matched wedge is odd under doubled internal grading", np.array_equal(grading@pairing@grading.T, -pairing))

    def transfer(h):
        eigenvalues, vectors = np.linalg.eigh(h)
        return (vectors*np.exp(-eigenvalues/10))@vectors.conj().T

    early, late = failures[0][0], failures[-1][0]
    r0, r1 = (transfer(matrices[g]) for g in (early, late))
    s0, s1 = (transfer(matrices[reflected(g)]) for g in (early, late))
    same = float(np.max(np.abs((r1@r0)@pairing@(s1@s0).T-pairing)))
    reverse = float(np.max(np.abs((r1@r0)@pairing@(s0@s1).T-pairing)))
    check("noncommuting reflected two-window replay depends on order", same < 1e-11 and reverse > 1e-8)

    # Vacuum and one normalized invariant wedge vector form the two axes.
    a = sp.diag(1, -1, 1)
    b = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    projector = sp.diag(1, 1, 0)
    t = sp.symbols("t", positive=True)
    t0, t1, t2, t3, t4 = sp.symbols("r:5", positive=True)
    segments = [sp.diag(t, t, value) for value in (t0, t1, t2, t3, t4)]
    r0, r1, r2, r3, r4 = segments
    word = r4*a*r3*b*r2*a*r1*b*r0
    bare = r4*r3*r2*r1*r0
    check("on the reducing vacuum and pair cap all scalar history weights cancel in the word ratio",
          a*a == b*b == sp.eye(3) and projector*word*projector == -projector*bare*projector)
    complement = sp.eye(3)-projector
    check("a cap on the complement has the opposite sign", complement*word*complement == complement*bare*complement)
    cap = sp.Matrix([[1, 1, 1]])
    check("even a mixed outer cap leaves only one selected overlap to compute",
          (cap*word*cap.T)[0] == sp.expand((cap*bare*cap.T-2*cap*projector*bare*cap.T)[0]))

    late_assembled = {"h_N":matrices[late]}
    yukawa_response = native_jet(assembler, late_assembled, 0.1)
    return {"status": "EXACT_SCOPED_CONTROLS_PASSED_F1_OPEN", "checks": checks,
            "physical_solution": False, "native_domain_digest": fixture["_sha256"],
            "one_particle_dimension": dimension, "profile_count": len(matrices),
            "same_branch_obstructions": len(failures), "strongest_one_leg_defect": max(d for _, d in failures),
            "matched_order_max_residual": same, "reversed_order_max_residual": reverse,
            "paired_wedge_norm_squared": dimension, "common_cap_word_ratio": -1,
            "paired_quartic_half_wedge_norm_squared": comb(dimension, 2),
            "general_cap_word_ratio": "1 - 2 * selected_paired_cap_response / bare_response",
            "fixed_same_branch_form_certificate": certificate,
            "paired_yukawa_response":yukawa_response,
            "scope": "actual original spatial matrices; paired exterior comparison, no OS time or physical cap binding"}
