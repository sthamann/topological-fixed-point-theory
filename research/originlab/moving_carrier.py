"""Replay the Hardy carrier and test whether its compression closes in time."""
import sympy as sp

from reconstruction import require


def replay():
    r, x, q = sp.symbols("r x q", positive=True)
    # r=|a|^8, x=q^4 and q=exp(-kappa*tau). Each nonconstant
    # normalized section has squared coefficient (1-r)*r^m at k+4*m.
    moment = lambda step: (1-r)/(1-r*x**step)
    defect = sp.factor(moment(2)-moment(1)**2)
    positive_defect = r*(1-r)*(1-x)**2/((1-r*x**2)*(1-r*x)**2)
    require(sp.cancel(defect-positive_defect)==0, "invalid compression identity")
    reference = {r:sp.S.Half,x:sp.S.Half}
    require(positive_defect.subs(reference)==sp.Rational(8,63), "reference defect changed")
    require(positive_defect.subs(r,0)==0 and positive_defect.subs(x,1)==0,
            "the defect must vanish at zero radius and at zero time")

    # On normalized moving sections the weight is q^k sqrt((1-r)/(1-r*x^2)).
    # At the moved radius r becomes r*x^2; composition therefore telescopes.
    normalizer = (1-r)/(1-r*x**2)
    composed_normalizer = normalizer*normalizer.subs(r,r*x**2)
    require(sp.cancel(composed_normalizer-(1-r)/(1-r*x**4))==0,
            "moving-carrier transports must compose")

    # Finite controls for the all-band statement: on each k mod 4,
    # successive images have distinct geometric ratios a^4*x^l.
    # Every initial finite Fourier window sees a Vandermonde matrix.
    a4 = sp.symbols("a4", nonzero=True)
    window = sp.Matrix(4,4,lambda m,l:(a4*x**l)**m)
    determinant = sp.prod((a4*x**j-a4*x**i) for i in range(4) for j in range(i+1,4))
    require(sp.expand(window.det()-determinant)==0, "incorrect Fourier control")

    # For coefficient t^m, Im<psi,d psi> depends only on the geometric
    # mean index r/(1-r), hence it is identical across the four classes.
    mean_index = sp.cancel(r*sp.diff(1/(1-r),r)*(1-r))
    require(sp.cancel(mean_index-r/(1-r))==0, "invalid Berry coefficient")
    return {
        "status":"EXACT_DECLARED_HARDY_CARRIER_CONTROL",
        "fixed_carrier_defect_rank":4,
        "defect_scalar_per_sector":str(positive_defect),
        "defect_at_reference":["0",*[str(sp.Rational(8,63)*2**(-sp.Rational(k,2))) for k in range(1,5)]],
        "moving_weight_squared":[str(q**(2*k)*normalizer) for k in range(1,5)],
        "moving_carrier_composition_exact":True,
        "cyclic_completion":"Full Hardy space for 0<|a|<1 and positive step. Orthogonality to all iterates "
                            "forces a power series in each residue class to vanish at a^4*x^l, which accumulate at zero.",
        "finite_stationary_transfer_reproduces_all_moments":False,
        "all_moment_obstruction":"In each nonconstant sector the spectral measure gives positive mass to "
                                 "every distinct q^(k+4*m). A nonzero polynomial cannot vanish on that support; "
                                 "the Hankel Gram matrix has full rank at every order.",
        "berry_sector_structure":"0 on the constant line and a common scalar on the four other classes",
        "physical_source_identity_proved":False,
        "first_original_test":"Use both fixed and transported projections at two ordered times on the "
                              "same charged R113 preparation and its adjoint, with the same absolute reference.",
    }
