"""Distinguish native mixed response, coordinate drift and common jumps."""

from fractions import Fraction as F
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import expm
import sympy as sp


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
old = runpy.run_path(str(base / "check_tfpt_conditioned_r113_weights_20261110.py"))
source, native = old["source"], old["native"]
gamma, swap, identity = old["gamma"], old["swap"], old["identity"]
higgs, majorana = old["clifford_higgs"], old["clifford_majorana"]
colors = old["blocks"]
omega = np.exp(2j * np.pi / 3)
z2, z3 = np.zeros((2, 2)), np.zeros((3, 3))


def real_generator(p):
    return -1j * (p + swap @ p.conj() @ swap)


weak_axes = [old["sigma_x"], np.array([[0, -1j], [1j, 0]]), old["sigma_z"]]
weak = [real_generator(native.rep(0, .5j * a, z3)) for a in weak_axes]
y = real_generator(native.rep(3j, z2, -1j * np.eye(3))) / 6
q = weak[2] + y


def comm(a, b):
    return a @ b - b @ a


def color_jet(r0, r_i, r_j, r_ij):
    terms = []
    for power in range(3):
        log_value = 0.
        sign = 1.
        li, lj, lij = 0., 0., 0.
        for idx, charge in colors:
            sub = np.ix_(idx, idx)
            phase = omega**(charge * power)
            matrix = np.eye(len(idx)) + phase * r0[sub]
            s, v = np.linalg.slogdet(matrix)
            sign *= s
            log_value += v
            ui, uj, uij = [np.linalg.solve(matrix, phase * t[sub]) for t in (r_i, r_j, r_ij)]
            li += np.trace(ui)
            lj += np.trace(uj)
            lij += np.trace(uij - ui @ uj)
        # Constant exp(log_value) is small enough in this finite native chart.
        terms.append(sign * np.exp(log_value) * np.array([1, li, lj, lij + li * lj]))
    z, zi, zj, zij = np.mean(terms, axis=0)
    connected = zij / z - zi * zj / z**2
    assert abs(connected.imag) < 1e-11
    return float(connected.real)


def scalar_log(r0):
    vals = []
    for power in range(3):
        sign, magnitude = 1., 0.
        for idx, charge in colors:
            d, l = np.linalg.slogdet(np.eye(len(idx)) + omega**(charge * power) * r0[np.ix_(idx, idx)])
            sign *= d
            magnitude += l
        vals.append(sign * np.exp(magnitude))
    mean = np.mean(vals)
    assert abs(mean.imag) < 1e-12 * mean.real and mean.real > 0
    return np.log(mean.real)


def native_matrix(radius_log, active):
    p = .25 * np.exp(-radius_log)
    h0 = p * gamma + (majorana if active else 0)
    h_left = h0 + (.35 * higgs if active else 0)
    h_right = h0 + (.55 * higgs if active else 0)
    left = expm(-.2 * h_left)
    middle = expm(-.4 * h_right)
    assert max(np.max(abs(comm(q, h))) for h in (h_left, h_right)) < 1e-13
    r0 = left @ middle @ left

    def susceptibility(a, b):
        assert np.max(abs(comm(a, b))) < 1e-13
        ri = 1j * left @ comm(middle, a) @ left
        rj = 1j * left @ comm(middle, b) @ left
        rij = -left @ comm(comm(middle, a), b) @ left
        return color_jet(r0, ri, rj, rij)

    gauge = [weak[2], y]
    matrix = np.array([[susceptibility(a, b) for b in gauge] for a in gauge])
    transverse = [susceptibility(a, y) for a in weak[:2]]
    assert np.max(abs(matrix @ np.ones(2))) < 1e-11
    assert max(map(abs, transverse)) < 1e-11
    if active:
        assert abs(matrix[0,1]) > 1e-3
    else:
        assert np.max(abs(matrix)) < 1e-11

    def vary(a, b):
        link = expm(1j * (a * weak[2] + b * y))
        return scalar_log(left @ link.conj().T @ middle @ link @ left)

    eps = 1e-3
    fd = (vary(eps, eps) + vary(-eps, -eps) - vary(eps, -eps) - vary(-eps, eps)) / (4 * eps**2)
    assert abs(fd - matrix[0,1]) < 2e-6
    assert abs(vary(.2, .2) - vary(0., 0.)) < 1e-12
    return dict(hessian_theta3_phi=matrix.tolist(), transverse_cross=transverse,
                finite_difference_cross=float(fd), fixed_EM_generator_exact=True)


def old_trace_susceptibility():
    # Marker in the original full finite triple, not a heat-kernel a6
    # coefficient by itself. T3 Y vanishes after the full H=0 trace.
    assert abs(np.trace(weak[2] @ y)) < 1e-13
    full_phi = source["inner_form"](np.array([0., 1.]))
    phi_with_majorana = full_phi + source["majorana"]
    direct = np.trace(weak[2] @ y @ phi_with_majorana @ phi_with_majorana).real
    u, d, e, nu = source["qref"].yukawa
    # v252 full_D positional order is verified below independently by the trace.
    norms = {name: float(np.linalg.norm(a)**2) for name, a in zip(("u", "d", "e", "nu"), (u,d,e,nu))}
    predicted = (norms["u"] - norms["d"] + norms["e"] - norms["nu"]) / 2
    assert abs(direct - predicted) < 1e-12
    return dict(trace_T3_Y_at_zero=0, trace_T3_Y_Phi_squared=float(direct),
                predicted_from_four_Yukawa_norms=predicted, yukawa_frobenius_squared=norms)


def infinitesimal_tests():
    # θ convention is exp(i θ T3). On D^j_mn and exp(i q φ),
    # the old Y generator at fixed U becomes ∂φ − X3 at fixed X.
    a, b, m, charge, k, rho = sp.symbols("a b m q k rho", real=True)
    block = sp.Matrix([[a+b, -b],[-b,b]])
    assert sp.factor(block.det()) == a*b
    expected = b*(charge-m)**2
    assert sp.expand(expected - (b*charge**2+b*m**2-2*b*m*charge)) == 0
    mixed_diffusion = k*rho**2*m*charge
    assert mixed_diffusion.subs(m,0) == 0
    assert mixed_diffusion.subs(m,-m) == -mixed_diffusion
    assert mixed_diffusion.subs(charge,-charge) == -mixed_diffusion
    positive_condition = sp.factor(sp.Matrix([[a,k*rho**2/2],[k*rho**2/2,b]]).det())
    assert sp.expand(positive_condition-(a*b-k*k*rho**4/4)) == 0
    jump_minus_jump = F(1,36) - F(-1,36)
    assert jump_minus_jump == F(1,18)

    # The short-time cross increment removes arbitrary scalar potential.
    x, angle = sp.symbols("x angle", real=True)
    v = sp.Function("v")(x, angle)
    def generator(f):
        return a*sp.diff(f,x,2) + b*sp.diff(f,angle,2) + k*rho**2*sp.diff(f,x,angle) + v*f
    covariance_probe = (generator(x*angle)-x*generator(angle)-angle*generator(x)
                        +x*angle*generator(sp.Integer(1)))/2
    assert sp.simplify(covariance_probe - k*rho**2/2) == 0

    # Bit proposal: the marginal flipping rate is a in each direction;
    # replace b of each single flip rate with one common flip rate b.
    flip_a, shared = sp.symbols("flip_a shared", nonnegative=True)
    eigen_one = -2 * (flip_a-shared) - 2 * shared
    eigen_joint = -4 * (flip_a-shared)
    assert sp.expand(eigen_one) == -2*flip_a
    assert eigen_joint.subs(shared,flip_a/2) == -2*flip_a

    metric, eps, time = sp.symbols("metric epsilon time", real=True)
    difference = -time*eps*metric**4
    assert all(sp.diff(difference, metric, i).subs(metric,0)==0 for i in range(4))
    return dict(coordinate_shifted_diffusion_matrix=str(block), coordinate_mixed_term="-2 b X3 d_phi",
                coordinate_weak_diagonal_addition="+b X3^2",
                intrinsic_mixed_diffusion_determinant=str(positive_condition),
                Plan232_jump_generator_difference_on_j1_m0_q2=str(jump_minus_jump),
                local_kappa_difference_on_same_mode=0,
                short_time_cross_increment_without_scalar_potential=str(sp.simplify(covariance_probe)),
                binary_marginal_generator=str(sp.expand(eigen_one)), binary_joint_generator=str(eigen_joint),
                absolute_reference_scale_example="c(m)=epsilon*m^4 preserves all center metric jets of orders 0..3, but changes second jets away from center")


print(json.dumps({"native_fixed_H_relative_link": {
    str(g): {"original": native_matrix(g, True), "S_zero_H_reference": native_matrix(g, False)}
    for g in (-.001, 0., .001)},
    "Higgs_dependent_original_trace": old_trace_susceptibility(),
    "exact_generator_tests": infinitesimal_tests(),
    "scope": "fixed native R032 96-mode, actual color trace. No full P_adm or stochastic generator inferred from the deterministic response."}, indent=2))
