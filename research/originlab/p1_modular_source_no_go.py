"""Exact E8_1 modular and OS-positive P1 modulus nonselection controls.

The chiral E8 vacuum character is E4/eta^8, not the original full TFPT
Q0. Synthetic positive source weights are countermodels, not native R113
outputs. A local physical source would need additional modulus curvature.
"""
from itertools import combinations
import sympy as sp
import mpmath as mp

def require(ok,message):
    if not ok:
        raise ArithmeticError(message)

def chi_E8(tau,terms=55):
    q=mp.exp(2j*mp.pi*tau)
    e4=1+240*sum(k**3*q**k/(1-q**k) for k in range(1,terms+1))
    e6=1-504*sum(k**5*q**k/(1-q**k) for k in range(1,terms+1))
    eta=mp.exp(mp.pi*1j*tau/12)*mp.fprod(1-q**k for k in range(1,terms+1))
    return e4/eta**8,e4,e6

def modular_saddle():
    mp.mp.dps=60
    value,e4,e6=chi_E8(1j)
    curvature=4*mp.pi**2*e4/3
    require(abs(value-12)<mp.mpf('1e-45') and
            abs(e6)<mp.mpf('1e-45') and
            curvature.real>0 and abs(curvature.imag)<mp.mpf('1e-45'),
            "E8 chiral square j/derivative check failed")
    x=mp.mpf('0.01')
    horizontal=abs(chi_E8(1j+x)[0])
    vertical=abs(chi_E8(1j*(1+x))[0])
    require(horizontal<12<vertical,
            "square tau=i is not the predicted chiral saddle")
    # Ramanujan derivative:
    # d log(chi)/dτ=-(2 pi i /3) E6/E4;
    # E6(i)=0, dE6(i)/dτ=-pi i E4(i)^2.
    # Thus d2 log chi(i)/dτ2 = -(2 pi^2/3)E4(i).
    return {
      "status":"EXACT_MODULAR_DERIVATIVE_NUMERICAL_60_DIGIT",
      "E8_vacuum_character":"E4(tau)/eta(tau)^8=j(tau)^(1/3)",
      "chi_E8_at_tau_i":"12",
      "E6_at_tau_i":"0",
      "log_chi_holomorphic_second_derivative_at_i":"-(2*pi^2/3)*E4(i)",
      "log_abs_chi_squared_Hessian_at_i":"diag(-4*pi^2*E4(i)/3,+4*pi^2*E4(i)/3)",
      "extra_source_curvature_required_on_imaginary_axis_if_factorized":
          "d_y^2 log(W_rest)(i)<-4*pi^2*E4(i)/3",
      "conditional_rest_curvature_threshold_numerical":mp.nstr(curvature.real,22),
      "abs_chi_i_plus_0p01":mp.nstr(horizontal,22),
      "abs_chi_1p01i":mp.nstr(vertical,22),
      "chiral_E8_square_selected_by_maximizing_character":False,
      "actual_TFPT_factorization_into_E8_character_times_rest_proved":False,
    }

def os_positive_source_models():
    t=sp.symbols("t",positive=True)
    shape=4*t**2/(1+t**2)**2
    require(sp.factor(shape.subs(t,1/t)-shape)==0 and
            sp.diff(shape,t,2).subs(t,1)==-2,
            "original antipodal mirror shape changed")
    rows={}
    for name,sgn in (("selects_square",1),("rejects_square",-1)):
        w=1+sgn*shape/2
        curvature=sp.factor(sp.diff(sp.log(w),t,2).subs(t,1))
        # H=diag(0,1), v=(1,1); K(τ)=v*e^(-τH)v=1+e^-τ.
        # On discrete times k ln(2), K_ij=w(t)(1+2^(-i-j)).
        Q=sp.Matrix(3,3,lambda i,j:w.subs(t,sp.Rational(1,2))*
                    (1+sp.Rational(1,2)**(i+j)))
        require(Q.rank()==2 and
                all(Q.extract(ix,ix).det()>=0
                    for n in (1,2,3) for ix in combinations(range(3),n)),
                "declared common OS-positive semigroup lost reflection positivity")
        require(sp.simplify(w.subs(t,1/t)-w)==0 and
                sp.diff(w,t).subs(t,1)==0,
                "declared source breaks antipodal mirror symmetry")
        rows[name]={
          "square_log_weight_curvature":str(curvature),
          "w_half_over_w_one":str(sp.factor(w.subs(t,sp.Rational(1,2))/w.subs(t,1))),
          "reflection_positive_time_Gram_rank":2,
          "same_time_semigroup":"diag(1,exp(-tau))",
          "positive_source_weight":True,
        }
    require(rows["selects_square"]["square_log_weight_curvature"]=="-2/3" and
            rows["rejects_square"]["square_log_weight_curvature"]=="2",
            "two symmetric positive sources no longer disagree on the square")
    return {
      "status":"EXACT_MIRROR_OS_POSITIVE_SOURCE_NONIDENTIFIABILITY",
      "geometric_shape":"4*t*t/(1+t*t)^2=sin(delta)^2",
      "positive_symmetry_compatible_weights":["1+shape/2","1-shape/2"],
      "two_sources_same_time_generator_and_symmetric_marks":True,
      "source_models":rows,
      "symmetry_reflection_positivity_and_sewing_select_square":False,
      "actual_original_weighted_Q0_source_used":False,
    }

def controls():
    return {
      "E8_chiral_modular":modular_saddle(),
      "positive_modulus_countermodels":os_positive_source_models(),
      "original_blocked_fiber_packet_produced":False,
      "original_R113_Z0_ZU_ZY_ZUY_produced":False,
      "F1_solved":False
    }

if __name__=="__main__":
    import json
    print(json.dumps(controls(),sort_keys=True,indent=2))
