"""Original-source necessary F1 tests: universal antipodal seam, determinant and Hardy Cauchy data.

Every result is an exact scoped calculation, not a native weighted R113
producer or a full physical no-go theorem. Works with SymPy only.
"""
import sympy as sp


def require(condition, text):
    if not condition:
        raise ArithmeticError(text)


def four_mark_family():
    """Original v512 one-parameter antipodal family with no square assumption."""
    z,u=sp.symbols("z u",nonzero=True)
    polynomial=sp.expand((z*z-1)*(z*z-u*u))
    r2=sp.diag(1,-1,1,-1,1)
    s=sp.Matrix(5,5,lambda i,j:u**(j-2) if i+j==4 else 0)
    inclusion=sp.Matrix(5,3,lambda i,j:int(i==j+1))
    evaluation=sp.Matrix([[1/u**2,0,0,0,0],[0,0,0,0,1]])
    family_r2=sp.diag(-1,1,-1)
    family_s=sp.Matrix(3,3,lambda i,j:-u**(j-1) if i+j==2 else 0)
    endpoint_s=sp.Matrix([[0,1],[1,0]])
    require(sp.simplify(polynomial.subs(z,-z)-polynomial)==0 and
            sp.simplify(polynomial.subs(z,u/z)-u*u*polynomial/z**4)==0,
            "antipodal divisor lacks the expected V4 actions")
    require(s*s==sp.eye(5) and s*r2==r2*s and
            s*inclusion==-inclusion*family_s and
            r2*inclusion==inclusion*family_r2 and
            evaluation*s==endpoint_s*evaluation and
            evaluation*inclusion==sp.zeros(2,3) and
            (evaluation.rank(),inclusion.rank())==(2,3),
            "equivariant 3+2 exact sequence requires the square")
    require(sp.simplify(polynomial.subs(z,sp.I*z)-polynomial-
                        2*z*z*(1+u*u))==0,
            "the genuine quarter symmetry condition changed")
    witness=(3+4*sp.I)/5
    require(sp.simplify(witness*witness+1)!=0 and
            sp.simplify(1/witness**2)!=0,
            "non-square witness unexpectedly has quarter symmetry")
    # For ALL u (not merely u=i), eta=dlog(z) has unit winding
    # and reflection orientation -1 under s:z->u/z.
    require(sp.residue(1/z,z,0)==1 and
            sp.simplify((-u/z**2)/(u/z)+1/z)==0,
            "logarithmic orientation changes with the modulus")
    return {
        "status":"EXACT_ANTIPODAL_FAMILY_COUNTERMODEL",
        "parameter":"u nonzero, u^2!=1; v512 restricted to unit circle",
        "one_sided_involution":"z -> -z",
        "pair_exchange":"z -> u/z",
        "three_plus_two_equivariant_for_all_noncolliding_moduli":True,
        "unit_dlog_winding_for_all_moduli":True,
        "D4_quarter_rotation_iff":"u^2=-1 (for z->i z)",
        "non_square_witness":"u=(3+4i)/5",
        "square_modulus_selected":False,
        "physical_R113_weight_computed":False
    }


def gaussian_arc_determinant():
    """Finite Gaussian over mean-zero vertex fluctuations on the v512 arc graph."""
    a,b=sp.symbols("a b",positive=True)
    d=sp.symbols("delta",positive=True)
    L=sp.Matrix([[a+b,-a,0,-b],[-a,a+b,-b,0],
                 [0,-b,a+b,-a],[-b,0,-a,a+b]])
    projector=sp.ones(4,4)/4
    detp=sp.factor((L+projector).det())
    require(detp==8*a*b*(a+b), "four-cycle Laplacian pseudodeterminant changed")
    geometric=sp.factor(detp.subs({a:1/d,b:1/(sp.pi-d)}))
    require(sp.simplify(geometric-8*sp.pi/(d*d*(sp.pi-d)**2))==0,
            "geometric determinant changed")
    # Gaussian partition function on the 3 nonzero modes:
    # Z ~ det'(L)^(-1/2) = delta (pi-delta) / sqrt(8 pi).
    zshape=d*(sp.pi-d)
    require(sp.diff(zshape,d).subs(d,sp.pi/2)==0 and
            sp.diff(zshape,d,2)==-2 and
            sp.simplify(sp.pi**2/4-zshape)==(d-sp.pi/2)**2,
            "Gaussian square extremum not unique")
    cw=2*sp.atan(sp.Rational(1,2))
    ratio=sp.simplify(4*cw*(sp.pi-cw)/sp.pi**2)
    require(0<float(ratio)<1,
            "non-square covariance witness does not disfavor the saddle")
    return {
        "status":"EXACT_DECLARED_GAUSSIAN_ARC_MODEL",
        "primitive_arc_generator":"L = d1 diag(1/delta,1/(pi-delta),1/delta,1/(pi-delta)) d1*",
        "nonzero_eigenvalues":["2/delta","2/(pi-delta)","2*pi/(delta*(pi-delta))"],
        "pseudodeterminant":"8*pi/(delta^2*(pi-delta)^2)",
        "relative_gaussian_weight":"delta*(pi-delta) / (pi*pi/4)",
        "unique_maximum_in_0_delta_pi":"delta=pi/2 (square)",
        "physical_original_boson_action_identified":False,
        "fermion_ghost_eta_factors_included":False,
        "sign_control":"reversing the net determinant exponent reverses square selection",
        "original_R113_computed":False,
    }


def calderon_hardy_obstruction():
    """Exact half-cylinder Dirac Fourier modes, before any physical identification."""
    sigma3=sp.diag(1,-1)
    gamma=sp.Matrix([[0,-1],[1,0]])
    require(gamma*gamma==-sp.eye(2) and gamma*sigma3==-sigma3*gamma,
            "Dirac product normal Clifford relations failed")
    # D=gamma(derivative_n+B), B_k=k sigma3; an L2 collar
    # solution e^{-n B_k}v exists only for positive B_k eigenvalue.
    bplus=sigma3
    bminus=-sigma3
    pplus=(sp.eye(2)+sigma3)/2
    pminus=(sp.eye(2)-sigma3)/2
    require(bplus*pplus==pplus and bminus*pminus==pminus and
            pplus*pminus==sp.zeros(2) and pplus+pminus==sp.eye(2),
            "positive-frequency Calderon projections changed")
    rplus=2*pplus-sp.eye(2)
    rminus=2*pminus-sp.eye(2)
    require(rplus!=rminus and rplus==-rminus,
            "a single local deck matrix unexpectedly equals both Calderon polarizations")
    return {
      "status":"EXACT_PRODUCT_DIRAC_CALDERON_COUNTERMODEL",
      "collar_operator":"gamma_n(partial_n + (-i partial_theta)*sigma3)",
      "P_C_at_k_plus_1":[[1,0],[0,0]],
      "P_C_at_k_minus_1":[[0,0],[0,1]],
      "required_2P_minus_I_at_k_plus_1":[[1,0],[0,-1]],
      "required_2P_minus_I_at_k_minus_1":[[-1,0],[0,1]],
      "no_frequency_independent_boundary_deck_matrix":True,
      "scope":"symbol countermodel to treating the Calderon polarization as an automatic geometric reflection",
      "full_TFPT_original_boundary_conditional_relation_proved":False,
    }


def controls():
    return {
      "four_mark_family":four_mark_family(),
      "gaussian_arc_determinant":gaussian_arc_determinant(),
      "calderon_hardy":calderon_hardy_obstruction(),
      "physical_source_identity_proved":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),sort_keys=True,indent=2))
