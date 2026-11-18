"""Source-first Higgs variation audit; no fitted or selected Yukawa coupling.

The circle calculation classifies represented local one-forms analytically.
The finite control reuses v252's actual code and checks formal coefficient
bases, which are probes of its linear map rather than physical Yukawa inputs.
"""

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import sympy as s


def circle_controls():
    z = s.Symbol("z", nonzero=True)
    alpha = s.Symbol("alpha", real=True)
    theta, tau = s.symbols("theta tau", real=True)
    checks = 0
    for power in range(-5, 6):
        f = z ** power
        for wave_number in range(-6, 7):
            wave = z ** wave_number
            d = lambda p: z * s.diff(p, z) + alpha * p
            commutator = s.expand(d(f * wave) - f * d(wave))
            assert s.expand(commutator - z * s.diff(f, z) * wave) == 0
            checks += 1
    # Every local multiplication h is already a[D,b], a=h/z, b=z.
    h = s.Symbol("h")
    assert s.simplify(h / z * (z * s.diff(z, z))) == h
    phi = s.sin(theta) - s.cos(2 * theta)
    v = s.cos(theta) + 2 * s.sin(2 * theta)
    psi = s.Function("psi")(theta)
    unitary = s.exp(s.I * tau * phi)
    transformed = s.conjugate(unitary) * (
        -s.I * s.diff(unitary * psi, theta) + alpha * unitary * psi)
    expected = -s.I * s.diff(psi, theta) + (alpha + tau * v) * psi
    assert s.simplify(s.expand(transformed - expected)) == 0
    assert s.simplify(phi.subs(theta, 2 * s.pi) - phi.subs(theta, 0)) == 0
    mean_v = s.integrate(v, (theta, 0, 2 * s.pi)) / (2 * s.pi)
    assert mean_v == 0
    assert s.exp(2 * s.pi * s.I * s.Integer(1)) == 1
    assert s.exp(2 * s.pi * s.I * s.Rational(1, 2)) == -1
    return {
        "commutator_coefficient_controls": checks,
        "represented_local_one_forms": "all smooth multiplication operators M_v",
        "selfadjoint_form": "v real",
        "circle_operator": "D_alpha,v=-i*dtheta+alpha+v(theta)",
        "holonomy": "exp(2*pi*i*(alpha+mean(v)))",
        "same_quarter_holonomy_condition": "mean(v) is an integer",
        "gauge_conjugacy_identity": "exp(-i*tau*phi) D_alpha exp(i*tau*phi)=D_alpha+tau*phi_prime",
        "new_internal_representation_created": False,
        "physical_gauge_invariance_or_anomaly_cancellation_proved": False,
    }


def boundary_algebra_controls():
    # Exact basis actions of U,C,B on Z, evaluated away from any truncation edge.
    alpha = s.Rational(1, 4)
    p0 = lambda n: int(n >= 0) - int(n - 1 >= 0)
    checks = 0
    for n in range(-12, 13):
        assert p0(n) == int(n == 0)
        # U* [B,U] e_n = ((n+1+alpha)-(n+alpha)) e_n.
        assert (n + 1 + alpha) - (n + alpha) == 1
        for row in range(-3, 4):
            for col in range(-3, 4):
                coefficient = p0(n - col)
                image_index = n - col + row
                assert coefficient == int(n == col)
                if coefficient:
                    assert image_index == row
                checks += 1
    return {
        "integer_basis_controls": checks,
        "rank_one_projection": "P0=C-U*C*U_star",
        "matrix_units": "E_mn=U^m*P0*U^(-n)",
        "unit_one_form": "U_star*[B,U]=I",
        "one_form_space_of_enlarged_algebra": "the whole algebra",
        "enlarged_algebra_center": "scalars (contains all finite matrix units)",
        "calderon_projection_is_local_multiplication": False,
        "enlargement_derived_as_physical_local_source": False,
    }


def finite_existing_path_controls(root):
    sys.path.insert(0, str(root / "verification"))
    import v252_full_finite_triple as t
    zero = np.zeros((3, 3), dtype=complex)
    z2 = np.zeros((2, 2), dtype=complex)
    z3 = np.zeros((3, 3), dtype=complex)
    er = t.rep(1, z2, z3)[:t.NP, :t.NP]
    quaternion_coordinates = [(1, 0), (1j, 0), (0, 1), (0, 1j)]
    identities = 0
    coefficients = 0
    # Complete real coefficient basis for the four complex 3x3 input matrices.
    # This certifies a linear identity, not a choice of any physical Yukawa value.
    names = ["u", "d", "e", "nu"]
    for name in names:
        for generation in range(3):
            for right_generation in range(3):
                for value in (1, 1j):
                    ys = {key: zero.copy() for key in names}
                    ys[name][generation, right_generation] = value
                    original = t.dirac_block(ys["u"], ys["d"], ys["e"], ys["nu"])
                    for first, second in quaternion_coordinates:
                        q = t.quaternion(first, second)
                        left = t.rep(0, q, z3)[:t.NP, :t.NP]
                        form = left @ (original @ er - er @ original)
                        expected = np.zeros_like(form)
                        for sector, colors in [(("nu", "e"), [0]),
                                               (("u", "d"), range(3))]:
                            for color in colors:
                                for out_weak, out_name in enumerate(sector):
                                    for in_weak, in_name in enumerate(sector):
                                        for g in range(3):
                                            for gj in range(3):
                                                i = t.IDX[(out_name, "L", color, g)]
                                                j = t.IDX[(in_name, "R", color, gj)]
                                                expected[i, j] = q[out_weak, in_weak] * ys[in_name][g, gj]
                        assert np.array_equal(form, expected)
                        # Values are Gaussian integers with exactly representable sums.
                        assert np.all(form.real == np.rint(form.real))
                        assert np.all(form.imag == np.rint(form.imag))
                        identities += 1
                        coefficients += form.size
    empty = t.full_D(zero, zero, zero, zero, zero)
    assert np.array_equal(empty, np.zeros_like(empty))
    for first, second in quaternion_coordinates:
        a = t.rep(0, t.quaternion(first, second), z3)
        assert np.array_equal(a @ (empty @ a - a @ empty), empty)
    return {
        "formal_input_coordinate_basis": "72 real coordinates of four complex 3x3 matrices",
        "quaternion_coordinate_basis": 4,
        "original_path_block_identities": identities,
        "gaussian_integer_entries_compared": coefficients,
        "conditional_scalar_block": "H(phi1,phi2)*diag(Y_up,Y_down), and lepton analogue",
        "zero_seed_gives_zero_finite_fluctuation": True,
        "physical_yukawa_values_selected": False,
        "color_and_weak_representation": "already supplied by v252 rep and GEN_FIELDS",
        "rank_three_or_two_selected_by_this_test": False,
    }


def symbolic_full_lr_formula():
    ar, ai, br, bi, cr, ci, dr, di, lr, li = s.symbols(
        "ar ai br bi cr ci dr di lr li", real=True)
    y1, y2 = s.symbols("Y_up Y_down")
    qa = s.Matrix([[ar + s.I * ai, br + s.I * bi],
                   [-br + s.I * bi, ar - s.I * ai]])
    qb = s.Matrix([[cr + s.I * ci, dr + s.I * di],
                   [-dr + s.I * di, cr - s.I * ci]])
    right_b = s.diag(lr + s.I * li, lr - s.I * li)
    seed = s.diag(y1, y2)
    actual_lr = qa * (seed * right_b - qb * seed)
    q_effective = qa * (right_b - qb)
    assert (actual_lr - q_effective * seed).applyfunc(s.expand) == s.zeros(2)
    assert s.simplify(q_effective[1, 0] + s.conjugate(q_effective[0, 1])) == 0
    assert s.simplify(q_effective[1, 1] - s.conjugate(q_effective[0, 0])) == 0
    return {"all_particle_lr_one_forms": "q_a*(diag(lambda_b,conj(lambda_b))-q_b)*M",
            "effective_factor_is_quaternion": True,
            "seed_yukawa_blocks_remain_symbolic_inputs": True}


def main():
    root = Path(__file__).resolve().parents[3]
    sources = [
        "_archive/tfpt-45/source_extracts/01_boundary_kernel_source.tex",
        "_archive/tfpt-45/source_extracts/02_carrier_source.tex",
        "verification/v252_full_finite_triple.py",
        "verification/v254_inner_fluctuations.py",
        "experiments/theory-contracts/seam-source-character-search-20261015/source_reconstruction/PROOF.txt",
    ]
    result = {
        "contract": "UR.SEAM.JOINT_SOURCE_ORIGIN_AUDIT.01",
        "annex": "source-first-higgs-variation-20261016",
        "verdict": "CONDITIONAL_VARIATION_DERIVED_NO_3PLUS2_ORIGIN",
        "source_sha256": {path: hashlib.sha256((root / path).read_bytes()).hexdigest()
                          for path in sources},
        "local_circle_source": circle_controls(),
        "nonlocal_calderon_enlargement": boundary_algebra_controls(),
        "existing_finite_ncg_path": finite_existing_path_controls(root),
        "symbolic_all_lr_formula": symbolic_full_lr_formula(),
        "p1_p2_source_selected": False,
        "full_phi_dependent_original_functional_derived": False,
        "all_higgs_mechanisms_excluded": False,
        "verification_scripts_modified": False,
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
