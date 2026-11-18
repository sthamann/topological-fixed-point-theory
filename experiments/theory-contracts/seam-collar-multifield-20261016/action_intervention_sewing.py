"""Exact checks of source propagation vs response-kernel sewing.

No intervention law is installed. Infinite shift products are checked by
their coefficients, with no finite/cyclic replacement of the circle shift.
"""
import json
import sys
from pathlib import Path
import sympy as s
from boundary_operation_response import occupation_chart, ratio


def check():
    rows=[]
    for n in range(-3,4):
        a=occupation_chart(n)
        raw_plus=a*occupation_chart(n+1)
        raw_minus=a*occupation_chart(n-1)
        sewn_plus=s.cancel(raw_plus/occupation_chart(n+1))
        sewn_minus=s.cancel(raw_minus/occupation_chart(n-1))
        assert sewn_plus==a and sewn_minus==a
        assert raw_plus!=raw_minus
        r=ratio(n)
        a_double=1/(1+r**4)
        assert r*r==ratio(n)**2
        assert a_double != a*a
        rows.append({'n':n,'ordinary_Lstar_L':str(raw_plus),
                     'ordinary_L_Lstar':str(raw_minus),
                     'inverse_Gram_sewn_both_orders':str(a),
                     'A_double_length':str(a_double),
                     'A_squared':str(a*a)})
    # Identity kernel fixes the sewing matrix in an exact finite regulator.
    a=s.diag(*(occupation_chart(n) for n in (-1,0,1)))
    variables=s.symbols('w:9');w=s.Matrix(3,3,variables)
    solution=s.solve(list(a*w*a-a),variables)
    solved=w.subs(solution)
    assert solved==a.inv()
    t,lam=s.symbols('t lam',real=True)
    scalar=1/(1+s.exp(-2*t*lam))
    derivative=s.simplify(s.diff(s.log(scalar),t,2).subs(t,0))
    assert derivative==-lam**2
    curvatures=[derivative.subs(lam,x) for x in (s.Rational(1,4),s.Rational(5,4))]
    assert curvatures==[-s.Rational(1,16),-s.Rational(25,16)]
    # A domain-safe all-bounded identity also records the form sewing:
    # F^-1 K_g = M_g F; no raw unbounded A^-1 is required in the proof.
    return {'verdict':'LF_NOT_DERIVED_AS_ACTUAL_INTERVENTION',
            'ordinary_response_noncommutativity_retained':True,
            'inverse_Gram_sewing_conditional_on_kernel_frame_interpretation':True,
            'inverse_Gram_sewing_selected_as_physical_law':False,
            'same_side_scalar_insertions_commute_under_that_sewing':True,
            'identity_sewing_matrix_regulator':str(solved),
            'log_A_second_derivative_at_zero':str(derivative),
            'quarter_source_two_mode_curvatures':[str(x) for x in curvatures],
            'scalar_normalization_of_A_to_source_semigroup_possible':False,
            'source_response_mode_checks':rows,
            'physical_source_perturbation_rule_derived':False,
            'measurement_projection_or_coupling_assumed':False}


if __name__=='__main__':
    result=check()
    target=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('action_intervention_sewing_results.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
