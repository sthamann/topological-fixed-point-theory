"""Independent scoped checks of the Yang–Mills trace and one Majorana threshold."""

import json
from fractions import Fraction

import mpmath as mp
import sympy as sp

from tfpt_constants import Mbar, N_fam
from v3_em_alpha import make_F
from v245_spectral_matter_content import GEN
from v246_unification_data_crosscheck import M_Z, ainv_1loop


def check_spin_trace():
    pauli = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]),
             sp.diag(1, -1)]
    gamma = [sp.kronecker_product(pauli[0], matrix) for matrix in pauli]
    gamma.append(sp.kronecker_product(pauli[1], sp.eye(2)))
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    omega = sp.symbols('o01 o02 o03 o12 o13 o23')
    ef = sum((gamma[i] * gamma[j] * entry for (i, j), entry in zip(pairs, omega)),
             start=sp.zeros(4))
    omega_sq = 2 * sum(entry**2 for entry in omega)
    assert sp.expand(sp.trace(ef**2) + 2 * omega_sq) == 0
    coefficient = (-sp.Integer(2) / 2 + sp.Integer(4) / 12) / (4 * sp.pi)**2
    assert sp.simplify(coefficient + 1 / (24 * sp.pi**2)) == 0


def check_gauge_trace():
    trace_su3 = sum(Fraction(weak, 2) for _, color, weak, _ in GEN if color == 3)
    trace_su2 = sum(Fraction(color, 2) for _, color, weak, _ in GEN if weak == 2)
    trace_y = sum(color * weak * charge**2 for _, color, weak, charge in GEN)
    full = tuple(2 * N_fam * trace for trace in (trace_su3, trace_su2, trace_y))
    assert full == (12, 12, 20)
    assert full[2] * Fraction(3, 5) == full[1]
    f0, g_squared = sp.symbols('f0 g_squared', positive=True)
    assert sp.solve(f0 * g_squared * full[0] / (24 * sp.pi**2) - sp.Rational(1, 4),
                    f0 * g_squared) == [sp.pi**2 / 2]


def check_profile():
    x, y = sp.symbols('x y', positive=True)
    f = lambda x: (sp.exp(-x) + sp.exp(-2 * x)) / 2
    gap = (sp.exp(-x) - sp.exp(-2 * x)) * (sp.exp(-y) - sp.exp(-2 * y)) / 4
    assert sp.simplify(f(x + y) - f(x) * f(y) - gap) == 0
    rho = sp.diag(sp.Rational(3, 4), sp.Rational(1, 4))
    modular_operator = sp.kronecker_product(rho, rho.inv())
    assert modular_operator.eigenvals() == {sp.Integer(1): 2, sp.Integer(3): 1,
                                            sp.Rational(1, 3): 1}
    return mp.mpf(str(sp.N(gap.subs({x: 1, y: 1}), mp.mp.dps)))


def check_mass(mass, planck, cutoff):
    ratio = mass / cutoff
    spectral_heat = 2 * cutoff**2 * mp.expm1(-ratio**2) / (48 * mp.pi**2)
    spectral_poisson = 4 * cutoff**2 * mp.expm1(mp.log1p(ratio) - ratio) / (48 * mp.pi**2)
    inverse_cutoff_sq = cutoff**-2
    relative_integral = (mp.expm1(-ratio**2) / inverse_cutoff_sq
                         - mass**2 * mp.e1(ratio**2))
    cutoff_loop = relative_integral / (192 * mp.pi**2)
    einstein = planck**2 / 2
    return {"cutoff_GeV": cutoff, "M_over_cutoff": ratio,
            "heat_C_over_EH_f0_1": spectral_heat / einstein,
            "poisson_C_over_EH_f0_1": spectral_poisson / einstein,
            "proper_time_loop_C_over_EH": cutoff_loop / einstein}


def main():
    mp.mp.dps = 55
    check_spin_trace()
    check_gauge_trace()
    gap = check_profile()
    c3 = 1 / (8 * mp.pi)
    alpha = mp.findroot(make_F(c3, 41), mp.mpf('0.0073'))
    alpha_as_uv_g_squared = 32 * mp.pi * alpha / 3
    alpha_as_uv_f0 = mp.pi**2 / (2 * alpha_as_uv_g_squared)
    planck = Mbar
    scalaron = c3**mp.mpf('3.5') * planck
    mass = 3 * scalaron
    cutoffs = {'M_scal': scalaron, 'Mbar': planck}
    gauge_scales = {'M_Z': M_Z, **cutoffs}
    gauge = {}
    for name, scale in gauge_scales.items():
        ai_1 = ainv_1loop(float(scale))
        gauge[name] = {'alpha_inv_1loop_GUT': ai_1.tolist(),
                       'f0_needed_1loop': (float(mp.pi) / 8 * ai_1).tolist()}
    # Compare the closed expressions to convergent quadrature at a non-small mass.
    x = mp.mpf('0.37')
    direct_poisson = mp.quad(lambda s: mp.exp(-1 / (4 * s) - s * x**2)
                            / (2 * mp.sqrt(mp.pi) * s**mp.mpf('2.5')), [0, 1, mp.inf])
    assert mp.almosteq(direct_poisson, 2 * (1 + x) * mp.exp(-x))
    direct_loop = mp.quad(lambda t: mp.expm1(-t * x**2) / t**2, [1, mp.inf])
    assert mp.almosteq(direct_loop, mp.expm1(-x**2) - x**2 * mp.e1(x**2))
    high_cutoff_limit = -2 * mass**2 / (24 * mp.pi**2 * planck**2)
    loop_msbar_at_m = -mass**2 / (96 * mp.pi**2 * planck**2)
    print(json.dumps({
        'status': 'conditional ordinary Spin(4) product / one free Majorana',
        'gauge_traces_SU3_SU2_U1Y': [12, 12, 20],
        'profile_gap_at_x_eq_y_eq_1': gap,
        'thomson_alpha_inverse_v3': 1 / alpha,
        'if_thomson_alpha_falsely_identified_with_same_scale_g_squared': alpha_as_uv_g_squared,
        'if_thomson_alpha_falsely_identified_with_same_scale_f0': alpha_as_uv_f0,
        'f0_gauge_RG': gauge,
        'Mbar_GeV': planck,
        'M_scal_GeV': scalaron,
        'v986_conditional_M_R_GeV': mass,
        'M_R_squared_over_Mbar_squared': (mass / planck)**2,
        'massive_minus_massless_thresholds': {
            name: check_mass(mass, planck, cutoff) for name, cutoff in cutoffs.items()},
        'spectral_high_cutoff_C_over_EH_f0_1': high_cutoff_limit,
        'spectral_high_cutoff_C_over_EH_f0_from_thomson_alpha': high_cutoff_limit * alpha_as_uv_f0,
        'one_majorana_loop_MSbar_at_mu_eq_M_C_over_EH': loop_msbar_at_m,
    }, default=lambda value: mp.nstr(value, 17), indent=2))


if __name__ == '__main__':
    main()
