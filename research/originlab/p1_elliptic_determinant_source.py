"""Unit-area rectangular elliptic cover spectral determinant P1 source gate.

The original 4 mark cross ratio x=t^2/(1+t^2) uniformizes a
rectangular double-cover torus tau=i*y via y=K(1-x)/K(x).
The known zeta Laplacian determinant y*abs(eta(i*y))^4 is uniquely
maximal at the square among all fixed-area rectangular tori.
Any ORIGINAL physical source determinant exponent and factorization
with the chiral E8_1 vacuum character remain UNPROVED.
"""
import mpmath as mp
import sympy as sp

def eta_iy(y,n=75):
    q=mp.exp(-2*mp.pi*y)
    return mp.exp(-mp.pi*y/12)*mp.fprod(1-q**k for k in range(1,n+1))

def y_of_t(t):
    x=t*t/(1+t*t)
    return mp.ellipk(1-x)/mp.ellipk(x)

def torus_det(y):
    return y*eta_iy(y)**4

def controls():
    mp.mp.dps=60
    t=sp.symbols("t",positive=True)
    x=t*t/(1+t*t)
    jx=256*(1-x+x*x)**3/(x*x*(1-x)**2)
    jt=256*(1+t*t+t**4)**3/(t**4*(1+t*t)**2)
    if sp.factor(jx-jt)!=0:
        raise ArithmeticError("elliptic double-cover j != native v512 j")
    sample={}
    for tv in (sp.Rational(1,3),sp.Rational(1,2),sp.Rational(2,3),sp.S.One):
        p=mp.mpf(int(sp.numer(tv)))/int(sp.denom(tv))
        y=y_of_t(p)
        sample[str(tv)]={"elliptic_lambda":mp.nstr(p*p/(1+p*p),22),
                         "rectangular_modulus_y":mp.nstr(y,22),
                         "unit_area_det":mp.nstr(torus_det(y),22)}
    value=torus_det(mp.mpf(1))
    if abs(value-mp.mpf("0.3483009824214192147961"))>mp.mpf("1e-21"):
        raise ArithmeticError("unit-area square determinant changed")
    for y in (mp.mpf("0.8"),mp.mpf("1.25"),mp.mpf("2")):
        if abs(torus_det(y)-torus_det(1/y))>mp.mpf("1e-47"):
            raise ArithmeticError("rectangular modular S determinant fails")
    if not all(torus_det(mp.mpf(y))<value
               for y in ("0.8","0.9","1.1","1.2","2","3")):
        raise ArithmeticError("square determinant extremum missing")
    # The exact logarithmic second derivative is STRICTLY negative for
    # all y>0: -1/y^2-16*pi^2 sum(n^2*q_n/(1-q_n)^2).
    q=mp.exp(-2*mp.pi)
    scalar_curvature=-1-16*mp.pi**2*sum(
       k*k*q**k/(1-q**k)**2 for k in range(1,76))
    E4=1+240*sum(k**3*q**k/(1-q**k) for k in range(1,76))
    chiral_curvature=4*mp.pi**2*E4/3
    threshold=chiral_curvature/(-scalar_curvature)
    if not (mp.mpf("14.75")<threshold<mp.mpf("14.76") and
            chiral_curvature+14*scalar_curvature>0 and
            chiral_curvature+15*scalar_curvature<0):
        raise ArithmeticError("E8 times determinant Hessian source gate failed")
    return {
      "status":"CONDITIONAL_RECTANGULAR_DOUBLE_COVER_DETERMINANT_SELECTOR",
      "elliptic_cover_parameter":"x=t^2/(1+t^2), y=K(1-x)/K(x)",
      "native_v512_j_equals_elliptic_j_of_x":True,
      "unit_area_torus_det":"D(y)=y*abs(eta(i*y))^4",
      "square_det":mp.nstr(value,25),
      "det_log_second_y_derivative_at_square":mp.nstr(scalar_curvature,25),
      "chiral_E8_abs_character_squared_log_second_y_derivative_at_square":mp.nstr(chiral_curvature,25),
      "conditional_chiral_E8_times_D_power_min_integer_for_local_square_max":15,
      "conditional_real_power_threshold":mp.nstr(threshold,22),
      "positive_det_power_maximizes_square_in_rectangular_class":True,
      "free_boson_negative_half_det_power_maximizes_square":False,
      "source_physical_elliptic_Dirac_Laplacian_identified":False,
      "source_net_ghost_fermion_boson_exponent_identified":False,
      "source_E8_character_factorization_proved":False,
      "native_Q0_and_R113_produced":False,
      "numeric_samples":sample,
      "F1_solved":False
    }

if __name__=="__main__":
    import json
    print(json.dumps(controls(),sort_keys=True,indent=2))
