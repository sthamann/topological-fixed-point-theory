"""Exact reflected modulus score/stiffness gate for finite source histories.

For fixed measure and domain Z(x)=sum_h mu_h exp(-S_h(x)),
d2 logZ/dx2 = Var(S'_h) - E(S''_h).
The original TFPT moving-domain, selector, eta and reference derivatives
remain extra physical terms, NOT silently removed.
"""
import sympy as sp

def controls():
    x,a,c=sp.symbols("x a c",real=True)
    Splus=c*x*x/2+a*x
    Sminus=c*x*x/2-a*x
    energies=(Splus,Sminus)
    Z=sum(sp.exp(-S) for S in energies)
    p=[sp.exp(-S)/Z for S in energies]
    avg=lambda v:sp.simplify(sum(p[i]*v[i] for i in range(2)))
    scores=[sp.diff(S,x) for S in energies]
    mean=avg(scores)
    variance=sp.simplify(avg([k*k for k in scores])-mean*mean)
    stiffness=avg([sp.diff(S,x,2) for S in energies])
    log_curvature=sp.simplify(sp.diff(sp.log(Z),x,2))
    if sp.simplify(log_curvature-variance+stiffness)!=0:
        raise ArithmeticError("native finite source-score identity violated")
    if sp.simplify(sp.diff(Z,x).subs(x,0))!=0:
        raise ArithmeticError("mirror-paired histories should be stationary")
    if sp.simplify(log_curvature.subs(x,0)-(a*a-c))!=0:
        raise ArithmeticError("source curvature a^2-c changed")
    witnesses={}
    for label,values in (
        ("square_maximum",{a:1,c:2}),
        ("square_minimum",{a:2,c:2}),
        ("square_flat_to_second_order",{a:1,c:1})):
        witnesses[label]=str(log_curvature.subs(values).subs(x,0))
    if witnesses!={"square_maximum":"-1","square_minimum":"2",
                  "square_flat_to_second_order":"0"}:
        raise ArithmeticError("positive mirror source did not yield all curvatures")
    return {
      "status":"EXACT_WEIGHTED_MODULUS_STRESS_VS_FLUCTUATION",
      "symmetric_coordinate":"x=log(t) with t=tan(delta/2), square x=0",
      "two_mirror_histories":"S_plus=c*x^2/2+a*x, S_minus=c*x^2/2-a*x",
      "exact_second_derivative":"d2log(Z)/dx2=Var(Sprime)-E(Ssecond)",
      "necessary_square_maximum":"E(Ssecond)>Var(Sprime)",
      "square_curvature":"a^2-c",
      "original_modulus_stress_fluctuations_computed":False,
      "original_moving_domain_projector_and_eta_variations_included":False,
      "mirror_only_selects_square":False,
      "positive_source_counterexamples":witnesses,
      "original_R113_four_values_produced":False,
      "F1_solved":False
    }

if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
