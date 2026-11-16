"""Yukawa response of the doubled reflected exterior cap."""

from itertools import combinations
from math import comb

import numpy as np
import sympy as sp

from reconstruction import require
from r113_weak_center import weak_center_response


def exterior_diagonal_response(size, degree, first_trace, second_trace, first_square_trace):
    """First and second derivatives of e_degree(M)/binomial(size, degree) at M=I."""
    singles = sp.Rational(degree, size)
    pairs = sp.Rational(degree*(degree-1), size*(size-1)) if size > 1 else sp.S.Zero
    return (singles*first_trace,
            singles*second_trace + pairs*(first_trace**2-first_square_trace))


def exact_controls():
    # A traceless change on one reflected leg. Expand determinants of the
    # actual principal minors independently of the coefficient formula.
    source = sp.symbols("y")
    direction = sp.diag(1, -1, 0, 0)
    marked = sp.diag(sp.exp(source), sp.exp(-source), 1, 1)
    measured = {}
    for degree in range(5):
        minors = [marked.extract(axes, axes).det()
                  for axes in combinations(range(4), degree)]
        overlap = sum(minors)/comb(4, degree)
        first, second = (sp.diff(overlap, source, order).subs(source, 0) for order in (1, 2))
        t2 = sp.trace(direction**2)
        require((first, second) == exterior_diagonal_response(4, degree, 0, t2, t2),
                "paired exterior derivatives disagree with principal minors")
        measured[str(degree)] = [str(first), str(second)]
    require(measured["1"] == ["0","1/2"] and measured["2"] == ["0","2/3"],
            "first response cancels but the second responds")
    # Dropping a partner changes the cap. The coordinate readout can see
    # the exact trace-free variation that the balanced cap averages out.
    require(sp.diff(marked[0,0], source).subs(source, 0) == 1,
            "coordinate cap should detect a linear mark")
    return measured


def native_jet(assembler, assembled, width):
    """Use the original tangential matrix with Majorana held under a Yukawa variation."""
    h = assembled["h_N"]
    modes = assembler.domain["modes"]
    bare = np.kron(np.diag([float(sp.Rational(n)+sp.Rational(assembler.domain["quarter"]))
                            for n in modes]), assembler.gamma)
    bare = bare + np.kron(np.eye(len(modes)), -1j*assembler.gamma@assembler.DM)
    v = h-bare
    size = h.shape[0]
    require(np.trace(h) == 0 and np.trace(v) == 0 and np.array_equal(v, v.conj().T),
            "original Yukawa variation must be Hermitian and traceless")
    eigenvalues, basis = np.linalg.eigh(h)
    yukawa_in_basis = basis.conj().T@v@basis
    # B = (d R) R^-1 is similar to the Hermitian matrix
    # R^-1/2 (d R) R^-1/2. In its eigenbasis the divided-difference
    # multiplier is -width*sinh(x)/x, including x=0.
    x = width*(eigenvalues[:,None]-eigenvalues[None,:])/2
    factor = np.ones_like(x)
    np.divide(np.sinh(x), x, out=factor, where=x!=0)
    hermitian = -width*factor*yukawa_in_basis
    norm_sq = float(np.sum(np.abs(hermitian)**2))
    require(abs(np.trace(hermitian)) < 1e-11 and norm_sq > 0,
            "native paired derivative replay")
    # On the degree-one cross pair the tangent leaves the paired line:
    # its component along that line is tr(B)/size, while its norm squared
    # is ||B||_F^2/size. Here B is similar, not equal, to hermitian.
    tangent_in_basis = np.exp(-x)*hermitian
    leakage_sq = float(np.sum(np.abs(tangent_in_basis)**2))/size
    tangent = basis@tangent_in_basis@basis.conj().T
    center = weak_center_response(assembler,h,v,tangent,width,assembler.gamma@assembler.K)
    gap = float(np.min(np.abs(eigenvalues)))
    # The sewn sea in the two-window control uses exp(-width Tr|h|/2).
    # Keep its derivative explicit and separate from the frozen-sea Y mark.
    varied_sea = None
    if gap > 1e-10:
        sea_score = float(-width/2*np.dot(np.sign(eigenvalues),np.diag(yukawa_in_basis).real))
        require(sea_score <= 1e-11, "weak-centre symmetry makes the sea action even and convex")
        varied_sea = {"derivative_log_sea":sea_score,"sample_gap":gap,
                      "scope":"differentiate only the late Yukawa in the sewn sea; comparison branch"}
    return {
        "width_role":"tangential transfer comparison width, not physical OS time",
        "width":width, "size":size, "original_yukawa_trace":"0",
        "first_y_response_frozen_sea":"0",
        "first_uy_response_frozen_sea_on_pure_paired_cap":"0",
        "second_y_response_paired_cap_before_sea":{
            str(k): k*(size-k)/(size*(size-1))*norm_sq for k in (1,2)
        },
        "first_y_tangent_norm_squared_outside_degree_one_pair_line":leakage_sq,
        "weak_center":center,
        "varied_sea_comparison":varied_sea,
        "scope":"matched spatial reflected legs, both external caps in the paired span, "
                "one-leg Yukawa variation and unchanged history scalar",
        "explicit_principal_minor_controls":exact_controls(),
        "offered_as_physical_source":False,
    }
