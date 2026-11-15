"""Exact bridge from the original v512 nonsquare j witness to E8_1 character.

For D_t={+1,-1,+u,-u}, u=(1+i*t)/(1-i*t), the anharmonic
elliptic j-invariant is j(t). Its real positive cube root is the
absolute E8_1 character on the associated elliptic double cover.
The native physical Q0 and R113 are not derived.
"""
import sympy as sp

def controls():
    t=sp.symbols("t",positive=True)
    j=sp.factor(256*(1+t*t+t**4)**3/(t**4*(1+t*t)**2))
    dlog=sp.factor(sp.diff(j,t)/j)
    target=2*(t*t-1)*(t*t+2)*(2*t*t+1)/(t*(1+t*t)*(1+t*t+t**4))
    if sp.simplify(dlog-target)!=0:
        raise ArithmeticError("elliptic j source derivative changed")
    points={str(x):str(sp.factor(j.subs(t,x))) for x in
            (sp.Rational(1,3),sp.Rational(1,2),sp.Rational(2,3),sp.S.One)}
    if points!={"1/3":"48228544/2025","1/2":"148176/25",
               "2/3":"37642192/13689","1":"1728"}:
        raise ArithmeticError("original v512 exact rational j witness changed")
    return {
       "status":"EXACT_V512_J_TO_CHIRAL_E8_VACUUM_BRIDGE",
       "source_four_marks":"(1,u,-1,-u), u=(1+i*t)/(1-i*t)",
       "elliptic_j":"256*(1+t^2+t^4)^3/(t^4*(1+t^2)^2)",
       "j_log_derivative":"2*(t^2-1)*(t^2+2)*(2*t^2+1)/(t*(1+t^2)*(1+t^2+t^4))",
       "j_strictly_decreases_for_t_in_0_to_1":True,
       "square_j":"1728",
       "square_abs_chi_E8":"12",
       "nonsquare_v512_t_half_j":"148176/25",
       "exact_heldout_j":points,
       "chiral_E8_magnitude_minimum_not_maximum_at_square":True,
       "original_weighted_source_factors_as_chiral_E8_proved":False,
       "original_R113_Z0_ZU_ZY_ZUY_computed":False,
       "F1_solved":False
    }

if __name__=="__main__":
    import json
    print(json.dumps(controls(),sort_keys=True,indent=2))
