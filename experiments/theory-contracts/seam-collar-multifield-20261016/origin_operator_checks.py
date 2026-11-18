"""Decisive original-operator and fixed-field tests, not a new TFPT source.

The compact symbol control keeps the quarter domain and reflection parity.
The CAR control retains the given source spectrum/thermal Ward identity while
testing whether these data fix the field dictionary. Its Fourier-mode unitary
is not asserted to be a local geometric P1 realization.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sympy as s
from source_sewing import creator


def original_perturbation_check():
    a = s.Symbol('a', real=True)
    k = s.Symbol('k', real=True)
    m = s.Symbol('m', positive=True)
    sx, sy = s.Matrix([[0,1],[1,0]]), s.Matrix([[0,-s.I],[s.I,0]])
    anchor = k*sx+m*sy
    perturbation = a*m*s.eye(2)
    total = anchor+perturbation
    # Reflection r -> -r with fiber lift sigma_x. The theta twist is fixed.
    even = sx*anchor.subs(m,-m)*sx-anchor
    odd = sx*perturbation.subs(m,-m)*sx+perturbation
    ratio = m/s.sqrt(m*m+s.Rational(1,16))
    slow = m-s.sqrt(m*m+s.Rational(1,16))
    return {
        'quarter_domain': 'psi(theta+2pi,r)=i psi(theta,r); periodic r',
        'geometry': 'compact doubled product torus; diagnostic of the claimed general perturbation lemma',
        'anchor_symbol_squared': str((anchor*anchor)[0,0]),
        'anchor_reflection_even': even == s.zeros(2),
        'perturbation_reflection_odd': odd == s.zeros(2),
        'relative_bound_ratio_limit_per_unit_amplitude': str(s.limit(ratio,m,s.oo)),
        'total_symbol_determinant': str(s.factor(total.det())),
        'at_amplitude_one_slow_eigenvalue_limit': str(s.limit(slow,m,s.oo)),
        'fixed_metric_normal_Clifford_residual': str(s.expand((sy+a*s.eye(2))**2-s.eye(2))),
        'correct_relative_bound': '|a|, not zero; |a|<1 suffices for Kato-Rellich on the anchor domain',
        'domain_failure_at_a1': 'unit modes with eigenvalues m-sqrt(m^2+1/16) approach zero; sum_m v_m/m is L2 but not H1, so the H1 realization is not closed',
        'actual_original_source_excluded': False,
        'same_metric_Dirac_type_for_this_family_requires': 'a=0',
        'repair': 'bounded symmetric order-zero perturbation; or proved relative bound<1; or independently verified elliptic total symbol and compatible selfadjoint boundary domain',
    }


def fixed_field_ward_check():
    # Actual first TWO excitation pairs of the supplied single chiral source.
    # These four field modes are not the P2 five-slot internal carrier.
    alpha = s.Symbol('alpha', real=True)
    dim = 16
    charges = (1,1,-1,-1)
    energies = (alpha,1+alpha,1-alpha,2-alpha)
    h = s.diag(*(sum(energies[j] for j in range(4) if state&(1<<j)) for state in range(dim)))
    q = s.diag(*(sum(charges[j] for j in range(4) if state&(1<<j)) for state in range(dim)))
    u = s.eye(dim)
    first, second = 5,10  # p0 h0 and p1 h1, neutral energies 1 and 3.
    u[first,first]=u[second,second]=s.Rational(3,5)
    u[first,second]=-s.Rational(4,5)
    u[second,first]=s.Rational(4,5)
    hu = s.simplify(u*h*u.adjoint())
    eye = s.eye(dim)
    low = [0,1,2,4,8]
    quarter = s.diag(*(s.I**int(q[j,j]) for j in range(dim)))
    shift = s.diff(-s.zeta(-1,1-alpha),alpha)
    ward = s.simplify(hu.diff(alpha)-q)
    # At t=log2, use the exact pair block only, retaining the common sea factor
    # separately. There is no fitted many-body matrix in this check.
    transfer = s.diag(*(2**(-h[j,j]) for j in range(dim)))
    tu = u*transfer*u.adjoint()
    pair = creator(4,0)*creator(4,2)
    vacuum = eye[:,0]
    original = (vacuum.adjoint()*pair.adjoint()*transfer*pair*vacuum)[0]
    changed = (vacuum.adjoint()*pair.adjoint()*tu*pair*vacuum)[0]
    c = creator(4,0)
    field_residual = hu*c-c*hu-alpha*c
    unitary_ok = u.adjoint()*u == eye
    low_fixed = all(u[:,j] == eye[:,j] for j in low)
    charge_ok = u*q == q*u
    similarity_ok = hu == u*h*u.adjoint()
    return {
        'scope': 'exact finite inverse-data diagnostic; not a selected geometric-local source or an inequivalent theory after transporting all fields',
        'actual_four_excitation_modes': [str(e) for e in energies],
        'unitary_exact': unitary_ok,
        'vacuum_and_entire_one_particle_sector_fixed': low_fixed,
        'same_full_charge_operator': charge_ok,
        'same_charged_quarter_operator': u*quarter == quarter*u,
        'quarter_field_phases': ['i','i','-i','-i'],
        'same_joint_energy_charge_spectrum_for_all_alpha': similarity_ok and charge_ok,
        'same_partition_for_all_beta_and_alpha': unitary_ok and similarity_ok,
        'same_sea_eta_identity_and_absolute_loop_phase': 'the identical scalar E_sea(alpha) and exp(2pi i E_sea) multiply both',
        'thermal_alpha_Ward_operator_residual_zero': ward == s.zeros(dim),
        'same_twist_response_with_full_sea': 'partial_alpha H_absolute = Q + partial_alpha E_sea in both models',
        'sea_response_at_quarter': str(shift.subs(alpha,s.Rational(1,4))),
        'same_all_vacuum_two_field_time_functions': unitary_ok and similarity_ok and low_fixed,
        'first_neutral_pair_return_at_t_log2': {'original':str(original),'fixed_fields_conjugated_time':str(changed)},
        'fixed_field_linear_equation_fails_after_conjugation': field_residual != s.zeros(dim),
        'geometric_locality_of_the_many_body_conjugation_proved': False,
        'P1_P2_source_selection_proved': False,
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    result={'schema':'tfpt.original-source-operator-checks.v1',
            'perturbation':original_perturbation_check(),
            'fixed_field_ward':fixed_field_ward_check(),
            'complete_P1_P2_solution':False}
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
