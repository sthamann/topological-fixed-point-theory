"""Decisive symbolic checks of actual proposed source-origin bridges.

These checks reject particular proposed inferences. They do not generate a new
P1 measure, replace its field space or certify complete TFPT.
"""
import json
from pathlib import Path
import sympy as sp


def need(test, message):
    if not test:
        raise RuntimeError(message)


def noncommuting_spectral_hessian():
    t = sp.symbols('t', real=True)
    # D(t)=diag(1,3)+t X, actual eigenvalues 2+-sqrt(1+t^2).
    actual_trace = 2*sp.exp(-2)*sp.cosh(sp.sqrt(1+t*t))
    actual = sp.expand(sp.diff(actual_trace, t, 2).subs(t, 0).rewrite(sp.exp))
    d, v = sp.diag(1, 3), sp.Matrix([[0, 1], [1, 0]])
    naive = sp.trace(sp.diag(sp.exp(-1), sp.exp(-3))*v*v)
    proposed_anticommutator = 2*naive
    need(sp.simplify(actual-(sp.exp(-1)-sp.exp(-3))) == 0, 'Exact matrix trace derivative failed')
    need(sp.simplify(actual-naive) == -2*sp.exp(-3), 'Noncommuting residual failed')
    return {'verdict': 'EXCLUDED_PROPOSED_GENERIC_HESSIAN_FORMULA',
            'actual_second_derivative': str(actual),
            'trace_f_second_times_V_squared': str(naive),
            'proposed_anticommutator_formula': str(proposed_anticommutator),
            'nonzero_residual': str(sp.simplify(actual-naive)),
            'scope': 'Exact e^(-u) trace variation in a two-dimensional noncommuting direction; not a full TFPT model or a chosen original insertion.'}


def modular_generator_roles():
    # The independently established ACTUAL quarter source has minimum positive
    # excitation energy min(alpha,1-alpha), for alpha=1/4 or3/4.
    cases=[]
    for alpha in (sp.Rational(1, 4), sp.Rational(3, 4)):
        gap = min(alpha, 1-alpha)
        scaled_gap = gap/2
        need(0 < scaled_gap < gap, 'Positive gap contraction did not contradict minimum')
        cases.append({'alpha': str(alpha), 'actual_positive_gap': str(gap),
                      'required_gap_after_unitary_scaling': str(scaled_gap)})
    return {'verdict': 'EXCLUDED_GAPPED_ROTATION_AS_NONTRIVIAL_BORCHERS_TRANSLATION',
            'cases': cases,
            'analytic_argument': 'For the full all-parameter dilation law, Stone uniqueness gives V H V*=H/2. Unitary conjugacy preserves spectrum and its least positive point. Scaling H halves that point, a contradiction. The spectral-set implication is separately kernel-checked in Lean.',
            'scope': 'The proposed use of the SAME actual gapped source H. Other gapless translation generators and correctly supplied continuum standard pairs remain allowed.'}


def kms_profile_typing():
    u = sp.symbols('u', nonnegative=True)
    ratios=[]
    for scale in (1, 2):
        cutoff=sp.exp(-scale*u)
        ratios.append(sp.integrate(cutoff, (u,0,sp.oo))/cutoff.subs(u,0))
    need(ratios == [1, sp.Rational(1,2)], 'Cutoff moments failed')
    return {'verdict': 'MISSING_MODULAR_TO_SPECTRAL_OPERATOR_IDENTIFICATION',
            'cutoff_moment_ratios': list(map(str,ratios)),
            'argument': 'KMS beta1 normalizes the given modular flow. It does not, by itself, equate its generator/weight with a separately declared spectral-action cutoff. Both positive Laplace profiles above exist mathematically. Their admissibility as complete P1 models is not asserted.',
            'source_check': 'v259 sets exp(-u) before integrating and records the decisive identification via check(..., True). Correct moment integrals are retained as conditional results.',
            'scope': 'Only the KMS-normalization-alone inference; no refutation of a correctly proved original modular/spectral dictionary.'}


def main():
    result={'spectral_hessian':noncommuting_spectral_hessian(),
            'modular_roles':modular_generator_roles(),
            'cutoff_role':kms_profile_typing(),
            'full_tfpt_solution':False,
            'next_target':'Derive the actual odd field Hessian from original S_adm and its local insertion map; then derive the full source identity, including allowed interactions, rather than input target moments.'}
    (Path(__file__).parent/'route-check-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v['verdict'] for k,v in result.items() if isinstance(v,dict)},indent=2))


if __name__=='__main__':
    main()
