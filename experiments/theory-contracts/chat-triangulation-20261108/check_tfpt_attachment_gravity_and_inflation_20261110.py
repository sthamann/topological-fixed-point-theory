"""Scoped attachment checks. No canonical TFPT files are changed."""
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

repo = Path(sys.argv[1])
sys.path.insert(0, str(repo / 'verification'))
from v86_nstar_reheating import nstar_chain


def tt_coefficients():
    q0, q3, a, b, rank = sp.symbols('q0 q3 A B N', real=True)
    q = sp.Matrix([q0, 0, 0, q3])
    h = sp.Matrix([[0, 0, 0, 0], [0, a, b, 0], [0, b, -a, 0], [0, 0, 0, 0]])
    assert h*q == sp.zeros(4, 1)
    assert sp.trace(h) == 0
    q2 = q.dot(q)
    rh = (q2*h - q*(h*q).T - (h*q)*q.T + sp.trace(h)*q*q.T)/2
    scalar_r = sp.trace(rh)
    assert scalar_r == 0
    # On closed/interior perturbations, drop integrated Euler and divergence.
    # Weyl² = 2 Ric² -(2/3) R² modulo Euler in four dimensions.
    a4 = -18*rank*(2*sp.trace(rh.T*rh)-sp.Rational(2, 3)*scalar_r**2)/(5760*sp.pi**2)
    assert sp.simplify(a4 + rank*q2**2*sp.trace(h.T*h)/(640*sp.pi**2)) == 0

    # R046 Dirac heat-kernel sign convention E=-R/4; N is total trace rank.
    R, positive_cutoff, rank_ref = sp.symbols('R omega2 N_ref', real=True)
    a2 = positive_cutoff*rank*(-R/4+R/6)/(4*sp.pi)**2
    a2_ref = a2.subs(rank, rank_ref)
    assert (a2-a2_ref).subs(rank_ref, rank) == 0
    return {'linear_scalar_curvature': str(scalar_r), 'ordinary_a4_TT': str(sp.factor(a4)),
            'weighted_ordinary_relative_a2': str(sp.factor(a2-a2_ref))}


def eh_sign():
    A, B, C, d, w2, chi, n, vr = sp.symbols('A B C d w2 chi n v_ref', positive=True)
    xi = sp.Rational(1, 6)
    mu_sq = 2*w2*chi**2/d-4*C/A
    quartic = 8*sp.pi**2*B/(d*A**2)
    v_sq = mu_sq/(2*quartic)
    absolute = w2*chi**2*n/(48*sp.pi**2)-xi*v_sq
    decomposed = w2*chi**2*(n-A**2/B)/(48*sp.pi**2) + d*A*C/(24*sp.pi**2*B)
    assert sp.simplify(absolute-decomposed) == 0
    rel_tree = -xi*(v_sq-vr**2)
    assert sp.simplify(rel_tree.subs(vr, sp.sqrt(v_sq))) == 0
    assert sp.simplify(rel_tree.subs(vr, 0)+xi*v_sq) == 0
    return {
        'ordinary_unsubtracted_positive_at_finite_cutoff':
            'w2*chi**2*(n-A**2/B)/(48*pi**2) + d*A*C/(24*pi**2*B)',
        'inequality': 'A**2 <= n*B and C >= 0; strict first bound for unequal DY**2 eigenvalues',
        'equal_rank_equal_geometry_relative_tree': '-xi*(v**2-v_ref**2)'
    }


def pull(measurement, ns):
    center, sigma = map(mp.mpf, measurement)
    return float((center-ns)/sigma)


def evaluate_tilt():
    datasets = {
        'ACT DR6 incl lowE optical-depth prior': ('0.9666', '0.0077'),
        'P-ACT': ('0.9709', '0.0038'),
        'P-ACT-LB DESI DR1': ('0.9743', '0.0034'),
        'P-ACT-LB2 DESI DR2': ('0.9752', '0.0030'),
        'CMB-SPA': ('0.9679', '0.0033'),
    }
    low, high = (1-2/mp.mpf(n) for n in (50, 60))
    # v86's own N_eff=4 and, under the same H=Gamma pivot prescription,
    # the conventional 4 minimally coupled species with denominator 192pi.
    orig_n, orig_t, orig_gamma = nstar_chain(4, '0.05')
    corr_n, corr_t, corr_gamma = nstar_chain(1, '0.05')
    assert mp.almosteq(orig_gamma/corr_gamma, 4)
    assert mp.almosteq(orig_t/corr_t, 2)
    xi = sp.Rational(1, 6)
    assert 4*(1-6*xi)**2 == 0
    assert abs(pull(datasets['P-ACT-LB2 DESI DR2'], high) - 2.8444444444444446) < 1e-12
    return {
        'band': [float(low), float(high)],
        'pull_sigma_by_dataset_lower_upper_sharp_v86': {
            key: [pull(value, n_s) for n_s in (low, high, 1-2/orig_n)]
            for key, value in datasets.items()},
        'minimal_comparison': [
            {'width_denominator': 48, 'N_star': float(orig_n), 'T_GeV': float(orig_t),
             'n_s_LO': float(1-2/orig_n), 'gamma_GeV': float(orig_gamma)},
            {'width_denominator': 192, 'N_star': float(corr_n), 'T_GeV': float(corr_t),
             'n_s_LO': float(1-2/corr_n), 'gamma_GeV': float(corr_gamma)}],
        'R048_xi': str(xi),
        'Higgs_width_factor_in_units_M3_over_192pi_Mbar2': 0
    }


mp.mp.dps = 35
print(json.dumps({'tt': tt_coefficients(), 'eh_sign': eh_sign(), 'inflation': evaluate_tilt()}, indent=2))
