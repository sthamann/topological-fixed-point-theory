"""Exact joint source integration for the supplied charged collar.

This is a conditional complex-fermion quantization, not a P1/P2 selection proof.
The primitive kinetic operator, source coupling and Berezin orientation are
explicit inputs. No determinant of the archival D_rel comparison part is used.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as s


def clean(poly):
    return {m: s.simplify(c) for m, c in poly.items() if s.simplify(c) != 0}


def add(a, b):
    result = dict(a)
    for mask, value in b.items():
        result[mask] = result.get(mask, 0) + value
    return clean(result)


def scale(a, value):
    return clean({m: c * value for m, c in a.items()})


def multiply(a, b):
    result = {}
    for left, ca in a.items():
        for right, cb in b.items():
            if left & right:
                continue
            inversions = sum((right & ((1 << i) - 1)).bit_count()
                             for i in range(left.bit_length()) if left >> i & 1)
            mask = left | right
            result[mask] = result.get(mask, 0) + (-1) ** inversions * ca * cb
    return clean(result)


def variable(index):
    return {1 << index: s.Integer(1)}


def exponential(poly, dimension):
    if poly.get(0, 0) != 0:
        raise ValueError("use an explicit scalar weight for the constant action")
    result, term = {0: s.Integer(1)}, {0: s.Integer(1)}
    for k in range(1, dimension + 1):
        term = scale(multiply(term, poly), s.Rational(1, k))
        if not term:
            break
        result = add(result, term)
    return result


def bilinear(matrix, bar_start=0, field_start=None):
    n = matrix.rows
    field_start = n if field_start is None else field_start
    out = {}
    for i in range(n):
        for j in range(n):
            out = add(out, scale(multiply(variable(bar_start + i),
                                          variable(field_start + j)), matrix[i, j]))
    return out


def berezin_source(kinetic):
    """Integrate exp(-barpsi K psi + bareta psi + barpsi eta) exactly.

    Internal order is all barpsi then all psi. The orientation is explicitly
    chosen so integral exp(-sum barpsi_j psi_j)=1. External order is bareta,eta.
    """
    kinetic = s.Matrix(kinetic)
    n = kinetic.rows
    if kinetic.cols != n or not 1 <= n <= 3:
        raise ValueError("bounded exact checker accepts square size 1..3")
    action = scale(bilinear(kinetic), -1)
    for j in range(n):
        action = add(action, multiply(variable(2*n + j), variable(n + j)))
        action = add(action, multiply(variable(j), variable(3*n + j)))
    integrand = exponential(action, 4*n)
    internal = (1 << (2*n)) - 1
    orientation = (-1) ** (n * (n+1) // 2)
    return clean({mask >> (2*n): orientation * value
                  for mask, value in integrand.items() if mask & internal == internal})


def determinant_source(kinetic):
    kinetic = s.Matrix(kinetic)
    if kinetic.det() == 0:
        raise ValueError("inverse chart requires a separate zero-mode treatment")
    return scale(exponential(bilinear(kinetic.inv()), 2*kinetic.rows), kinetic.det())


def schur_source(kinetic, retained):
    """Return both the effective boundary action and the interior scalar."""
    kinetic = s.Matrix(kinetic)
    retained = list(retained)
    if len(set(retained)) != len(retained) or any(i not in range(kinetic.rows) for i in retained):
        raise ValueError("retained indices must be unique valid field indices")
    inner = [i for i in range(kinetic.rows) if i not in retained]
    if not inner:
        return kinetic.extract(retained, retained), s.Integer(1)
    ii = kinetic.extract(inner, inner)
    if ii.det() == 0:
        raise ValueError("interior chart contains zero modes")
    eff = (kinetic.extract(retained, retained)
           - kinetic.extract(retained, inner)*ii.inv()*kinetic.extract(inner, retained))
    return eff.applyfunc(s.simplify), s.factor(ii.det())


def profile_obstruction(profile, x):
    """Necessary commuting-sector test for S_f=-logdet K+constant.

    K must first be the actual fermionic kinetic operator. It is not identified
    with D_rel by this diagnostic. A generic heat profile fails this identity.
    """
    return s.simplify(x*s.diff(profile, x) + 1)


def sea_energy_zeta(alpha):
    if not 0 < alpha < 1:
        raise ValueError("this chart excludes zero-mode domains")
    a = 1-alpha
    return s.simplify((a*a - a + s.Rational(1, 6))/2)


def circle_source(alpha):
    """Zeta determinant and Green response of the GIVEN angular Dirac source.

    This angular spectral amplitude is not a positive thermal density matrix.
    Negative eigenvalues have spectral argument +pi. There are no zero modes
    for 0<alpha<1. The diagnostic conjugate domain is not archival D_ref.
    """
    alpha = s.sympify(alpha)
    if not 0 < alpha < 1:
        raise ValueError("the zero-mode domains need another determinant chart")
    h = s.simplify(s.expand_complex(s.exp(2*s.pi*s.I*alpha)))
    eta = 1-2*alpha
    modulus = 2*s.sin(s.pi*alpha)
    determinant = s.simplify(s.expand_complex(modulus*s.exp(-s.I*s.pi*eta/2)))
    before = s.simplify(s.I/(h-1))
    after = s.simplify(h*before)
    return {'holonomy': h, 'eta': eta, 'modulus': modulus,
            'determinant': determinant, 'green_before': before, 'green_after': after,
            'pv_inverse_trace': s.pi*s.cot(s.pi*alpha),
            'zeta_logdet_derivative': s.pi*s.cot(s.pi*alpha)+s.pi*s.I}


def spectral_state_ward(alpha, beta, modes=1):
    """Joint spectral/thermal Ward identities in the same specified zeta scheme.

    The source coordinate alpha is held in a single zero-free chart. The
    response is a spectral charge for this supplied Hamiltonian, not a physical
    event or charge-selection derivation from P1/P2. The full sea counterterm
    is used even when the thermal excitation product is truncated.
    """
    alpha, beta = s.sympify(alpha), s.sympify(beta)
    if not 0 < alpha < 1 or beta.is_positive is not True:
        raise ValueError('need 0<alpha<1 and positive beta')
    if not isinstance(modes, int) or modes < 1:
        raise ValueError('need at least one particle/hole pair')
    x = s.Symbol('ward_alpha', real=True)
    energy = -s.zeta(-1, 1-x)
    eta = 1-2*x
    sea_derivative = s.diff(energy, x)
    logz = -beta*energy
    occupation_charge = sea_derivative
    for n in range(modes):
        plus, minus = n+x, n+1-x
        logz += s.log(1+s.exp(-beta*plus))+s.log(1+s.exp(-beta*minus))
        occupation_charge += 1/(1+s.exp(beta*plus))-1/(1+s.exp(beta*minus))
    source_charge = -s.diff(logz, x)/beta
    at = lambda expr: s.simplify(expr.subs(x, alpha))
    return {
        'energy_eta_residual': at(energy-eta**2/8+s.Rational(1,24)),
        'sea_derivative_eta_residual': at(sea_derivative+eta/2),
        'thermal_ward_residual': at(source_charge-occupation_charge),
        'sea_spectral_charge': at(sea_derivative),
        'charge_from_log_partition': at(source_charge),
        'charge_from_occupations': at(occupation_charge),
    }


def control():
    # The old source is reused, including its signed graph and particle/hole map.
    from charged_pair_control import verify_charged_pair
    pair = verify_charged_pair()
    r = s.diag(s.Rational(1, 4), s.Rational(1, 64))
    # Closing TWO collar amplitudes realizes the earlier beta=2*ell state.
    kinetic = s.eye(2) + r*r
    integral = berezin_source(kinetic)
    expected = determinant_source(kinetic)
    assert integral == expected
    eff, interior_scalar = schur_source(kinetic, [0])
    assert s.factor(interior_scalar*eff.det()) == kinetic.det()
    assert eff.inv() == kinetic.inv().extract([0], [0])
    n = s.eye(2)-kinetic.inv()
    assert n == s.diag(s.Rational(1, 17), s.Rational(1, 4097))
    x = s.Symbol('x', positive=True)
    assert profile_obstruction(-s.log(x), x) == 0
    assert profile_obstruction(s.exp(-x*x), x) != 0
    alpha = s.Rational(1, 4)
    vacuum_energy = sea_energy_zeta(alpha)
    # This is a SPECIFIED zeta/Abel prescription, not a P1/P2 consequence.
    beta = 8*s.log(4)
    vacuum_factor = s.exp(-beta*vacuum_energy).simplify()
    assert vacuum_factor == 2**s.Rational(1, 6)
    angular, conjugate = circle_source(alpha), circle_source(1-alpha)
    assert angular['determinant'] == 1-s.I
    assert conjugate['determinant'] == 1+s.I
    ward = spectral_state_ward(alpha, beta)
    assert ward['energy_eta_residual'] == ward['sea_derivative_eta_residual'] == ward['thermal_ward_residual'] == 0
    return {
        'schema': 'tfpt.joint-quadratic-source.v1',
        'verdict': 'EXACT_JOINT_INTEGRATION_UNDER_EXPLICIT_FERMION_MEASURE',
        'original_pair_verified': pair['ok'],
        'source_inputs': {
            'charged_quarter': 'alpha=1/4, same n=0,-1 chiral component',
            'collar_length': 'ell=4 log4', 'closed_length': 'beta=2 ell',
            'kinetic': 'I+R_ell^2 under ordinary fermionic trace closure',
            'quadratic_Berezin_measure_assumed': True,
            'orientation': 'integral exp(-sum barpsi_j psi_j)=1',
            'physical_K_identified_from_P1_P2': False,
            'archival_D_rel_is_not_automatically_full_kinetic_K': True,
        },
        'joint_outputs': {
            'determinant_scalar': str(kinetic.det()),
            'coherent_source_covariance_diag': [str(v) for v in kinetic.inv().diagonal()],
            'normal_ordered_occupation_diag': [str(v) for v in n.diagonal()],
            'integrated_mode_determinant_must_be_retained': str(interior_scalar),
            'effective_retained_mode_kinetic': str(eff[0, 0]),
            'partition_scope': 'retained particle mode and integrated hole mode; not a geometric collar cut',
            'geometric_boundary_pairing_newly_derived': False,
            'dropping_integrated_factor_changes_source_normalization': True,
            'Berezin_and_determinant_sources_identical_exactly': True,
            'all_field_source_coefficients': {str(k): str(v) for k, v in integral.items()},
        },
        'primitive_profile_test': {
            'required_identity': 'x f_prime(x)+1=0 for S_f=-logdet K',
            'solutions': 'f(x)=-log(x)+constant on the tested positive interval',
            'log_profile_passes': True,
            'heat_profile_residual': str(profile_obstruction(s.exp(-x*x), x)),
            'original_profile_or_fermionic_measure_selected_by_this_checker': False,
        },
        'conditional_absolute_zeta_control': {
            'prescription': 'Hurwitz-zeta / exponential-Abel subtraction for the given charged sea',
            'sea_energy': str(vacuum_energy),
            'closed_vacuum_factor': str(vacuum_factor),
            'two_mode_closed_source_zero_insertions': str(s.factor(vacuum_factor*kinetic.det())),
            'q_beta': 'exp(-beta)=1/65536; not q_ell=1/256',
            'full_zero_source_product': 'exp(-beta E_sea_zeta)*product_n(1+q_beta^(n+alpha))*(1+q_beta^(n+1-alpha))',
            'full_functional_for_finite_mode_sources': 'Z[0]*exp(bareta_S (I+R_beta)^(-1)_SS eta_S)',
            'higher_modes_tail_bounds': 'existing PROOF section8, with length beta',
            'field_quarter_remains': 'i on particles, -i on holes',
            'absolute_loop_central_phase': 'exp(2pi i E_sea_zeta)=exp(-i pi/48)',
            'central_phase_is_not_the_field_quarter': True,
            'eta_phase_cancellation_or_P1_identification_proved': False,
            'P1_P2_selection_of_prescription_proved': False,
        },
        'conditional_angular_dirac_source': {
            'coordinate': 'angular theta in [0,2pi]; not normal collar time or a thermal state',
            'operator': '-i partial_theta, f(2pi)=i f(0)',
            'alpha_variation_comparison_frame': 'periodic frame: -i partial_theta+alpha; Green values above are in the twisted frame',
            'complex_fermion_zeta_quantization_assumed': True,
            'negative_spectral_argument': '+pi',
            'zeta_D_squared_at_zero': '0',
            'eta_at_zero': str(angular['eta']),
            'determinant_modulus': str(s.simplify(angular['modulus'])),
            'determinant': str(angular['determinant']),
            'green_before_insertion': str(angular['green_before']),
            'green_after_insertion': str(angular['green_after']),
            'green_jump': str(s.simplify(angular['green_after']-angular['green_before'])),
            'field_holonomy': str(angular['holonomy']),
            'diagnostic_conjugate_domain_reference': 'alpha_ref=3/4; not identified with archival D_ref',
            'diagnostic_relative_determinant': str(s.simplify(angular['determinant']/conjugate['determinant'])),
            'paired_determinant_norm_squared': str(s.simplify(angular['determinant']*s.conjugate(angular['determinant']))),
            'pv_inverse_trace': str(angular['pv_inverse_trace']),
            'zeta_logdet_derivative': str(angular['zeta_logdet_derivative']),
            'eta_phase_contact_term': 'i*pi; must not be dropped in the alpha Ward variation',
            'source_functional': 'det_zeta(D_alpha)*exp(bareta D_alpha^(-1) eta) on smooth finite Grassmann sources',
            'log_abs_profile_matches_original_eta_form_only_if_kinetic_and_reference_are_identified': True,
            'P1_P2_quantization_or_reference_identification_proved': False,
        },
        'shared_spectral_thermal_ward': {
            'same_prescription_required': 'same Hurwitz-zeta scheme for eta and negative sea; 0<alpha<1',
            'energy_eta_identity': 'E_sea=eta_D(0)^2/8-1/24',
            'sea_response_identity': 'partial_alpha E_sea=-eta_D(0)/2',
            'thermal_identity': '-partial_alpha logZ/beta=-eta/2+sum_n(f_beta(n+alpha)-f_beta(n+1-alpha))',
            'lowest_pair_values': {k: str(v) for k, v in ward.items()},
            'scope': 'full prescribed sea counterterm plus truncated excitation product; spectral response, not P1/P2 event selection',
            'angular_contact_term_is_a_different_response': True,
        },
        'absolute_original_TFPT_source_recovered': False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    result = control()
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k: result[k] for k in ('verdict', 'joint_outputs',
                                            'absolute_original_TFPT_source_recovered')}, sort_keys=True))


if __name__ == '__main__':
    main()
