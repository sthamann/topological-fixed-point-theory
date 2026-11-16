"""Weak-centre projection of the original Yukawa pair response."""

import numpy as np
import sympy as sp

from reconstruction import require


def reversible_sewing_control():
    """Same local sea marks and transfer spectrum, different two-window contractions."""
    # For a charged R032 pair with energy sqrt(p^2 + sigma r^2), use
    # sigma r^2 = 0, 1, 9/4. This factors out the common 1/sigma in U.
    u = sp.Matrix([0, 1, sp.Rational(9,4)])
    p = sp.Rational(1,4)
    y = u.applyfunc(lambda value: value/sp.sqrt(p*p+value))
    constant, curvature, edge = (sp.Matrix(values) for values in
                                  ((1,1,1),(1,-2,1),(1,0,-1)))
    q = curvature*curvature.T/6
    r = edge*edge.T/2
    first, second = q+2*r, 2*q+r
    require(all(g*constant == sp.zeros(3,1) and g.eigenvals() == {0:1,1:1,2:1}
                and all(g[i,j]<0 for i in range(3) for j in range(3) if i!=j)
                for g in (first,second)),
            "the two original-score controls must be irreducible reversible generators")
    a = (u.T*q*y)[0]/3
    b = (u.T*r*y)[0]/3
    z = sp.simplify(-a/b)
    require(a.is_negative and b.is_positive and z.is_positive and (1-z).is_positive,
            "the shared local marks must separate the two sewings")
    root = sp.simplify(a*z+b*z*z)
    changed = sp.simplify(a*z*z+b*z)
    require(root == 0 and changed.is_positive,
            "same stationary marks and spectrum did not yield different contractions")
    equal_time = a+b
    require(equal_time.is_positive, "the sample must share a strictly positive equal-time covariance")
    # If a positive physical generator is supplied, 0 <= I-exp(-tH) <= tH
    # yields a short-time lower bound using only its equal-time forms.
    safe_times = []
    for generator in (first,second):
        energy_u = (u.T*generator*u)[0]/3
        energy_y = (y.T*generator*y)[0]/3
        require(energy_u.is_positive and energy_y.is_positive, "invalid Dirichlet forms")
        safe_times.append(float(equal_time/sp.sqrt(energy_u*energy_y)))
    require(max(safe_times) < float(-sp.log(z)), "the crossing must be outside the certified interval")
    return {
        "status":"ISOSPECTRAL_REVERSIBLE_SEWING_NOT_IDENTIFIED_BY_LOCAL_MARKS",
        "local_mark_u_times_sigma":[str(x) for x in u],
        "local_mark_y":[str(x) for x in y],
        "same_transfer_eigenvalues_exponents":[0,1,2],
        "same_normalized_stationary_law":["1/3"]*3,
        "exp_minus_time":str(z),
        "connected_zero":"0",
        "connected_other_times_sigma":str(changed),
        "connected_other_times_sigma_numeric":float(changed),
        "same_equal_time_connected_times_sigma":float(equal_time),
        "strict_positive_times_less_than":safe_times,
        "scope":"one original charged-pair sea-score shape; two auxiliary bosonic sewings, "
                "neither identified with the actual original endpoint kernel",
    }


def weak_center_response(assembler, h, yukawa, tangent, width, pair):
    """Compare the joint gauge action and a fermion-only projection at fixed Higgs."""
    left = assembler.original.rep(1, -np.eye(2), np.eye(3))
    real_structure = assembler.K
    center = left @ real_structure @ left.conj() @ real_structure
    signs = np.diag(center).real
    require(np.array_equal(center, np.diag(signs)) and np.all(signs**2 == 1),
            "expected the original involutive weak-centre action")

    def conjugate(matrix, signature):
        return (signature[:, None]*matrix)*signature[None, :]

    require(np.array_equal(conjugate(assembler.DM, signs), assembler.DM)
            and all(np.array_equal(conjugate(v, signs), -v) for v in assembler.V)
            and np.array_equal(conjugate(pair, signs), pair),
            "wrong weak-centre action on Majorana, Higgs directions or pair")

    size = h.shape[0]
    full_signs = np.tile(signs, len(assembler.domain["modes"]))
    baseline = h-yukawa
    require(np.array_equal(conjugate(baseline, full_signs), baseline)
            and np.array_equal(conjugate(yukawa, full_signs), -yukawa),
            "the whole native matrix does not have the expected Higgs parity")
    reflected = conjugate(tangent, full_signs)
    even, odd = (tangent+reflected)/2, (tangent-reflected)/2

    # The pair is invariant, so inner products between its marked legs
    # reduce to the normalized one-particle Hilbert–Schmidt pairing.
    overlap = float(-np.vdot(yukawa, tangent).real/size)
    lower_bound = float(width*np.vdot(yukawa, yukawa).real/size)
    require(lower_bound > 0 and overlap >= lower_bound-1e-12,
            "the canonical Yukawa cap must detect the tangent")
    require(abs(np.vdot(yukawa, even)) < 1e-11,
            "the fermion-only weak-centre average must erase this pairing")
    return {
        "status":"WEAK_CENTRE_RETAINS_JOINT_YUKAWA_CAP",
        "internal_center_multiplicities":[int(np.sum(signs == -1)), int(np.sum(signs == 1))],
        "joint_projector":"(I + field_sign_flip ⊗ fermion_weak_center)/2",
        "joint_cap_pairing":overlap,
        "joint_cap_pairing_lower_bound":lower_bound,
        "fermion_only_cap_pairing":"0",
        "even_tangent_norm_squared":float(np.vdot(even,even).real/size),
        "odd_tangent_norm_squared":float(np.vdot(odd,odd).real/size),
        "reversible_sewing_control":reversible_sewing_control(),
        "whole_original_selector_evaluated":False,
    }
