"""Test central sign and phase symmetries of the actual R032 Yukawa vertex."""
import contextlib
import io
import itertools
import json
import runpy
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
import sympy as sp

repo = Path(sys.argv[1]).resolve()
directory = repo / "experiments/theory-contracts/chat-triangulation-20261108"
path = directory / "check_original_whole_two_strips_20261109.py"
argv = sys.argv
try:
    sys.argv = [str(path), str(repo)]
    with contextlib.redirect_stdout(io.StringIO()):
        existing = runpy.run_path(str(path))
finally:
    sys.argv = argv

native, gamma = (existing[k] for k in ("native", "grading"))
source = existing["source"]
higgs = source["inner_form"]
majorana = existing["majorana"]
original_yukawas = source["qref"].yukawa
history,durations = (existing[k] for k in ("history", "durations"))


def exact_complex(x):
    return sp.Rational(float(x.real))+sp.I*sp.Rational(float(x.imag))


def family_commutant(first, second):
    """Read the joint Gram commutant, exact on the stored binary inputs."""
    eye = np.eye(3)
    grams = [y@y.conj().T for y in (first,second)]
    normalized = [g/np.linalg.norm(g,2) for g in grams]
    map_numeric = np.vstack([np.kron(g.T,eye)-np.kron(eye,g) for g in normalized])
    singular = np.linalg.svd(map_numeric,compute_uv=False)
    matrices = [sp.Matrix([[exact_complex(x) for x in row] for row in y])
                for y in (first,second)]
    grams_exact = [y*y.conjugate().T for y in matrices]
    basis = [sp.zeros(3) for _ in range(9)]
    for i,mat in enumerate(basis):
        mat[i%3,i//3]=1
    # Only tracefree directions can be additional commuting symmetries.
    # A nonzero 8x8 minor already proves exact rank 8, avoiding a large
    # rational symbolic SVD or rref.
    import scipy.linalg
    _,_,column_permutation = scipy.linalg.qr(map_numeric,pivoting=True)
    columns = list(map(int,column_permutation[:8]))
    _,_,row_permutation = scipy.linalg.qr(map_numeric[:,columns].T,pivoting=True)
    rows = list(map(int,row_permutation[:8]))
    exact_columns = [[
        sp.expand((b*g-g*b)[k%3,k//3])
        for g in grams_exact for k in range(9)] for b in basis]
    certificate = sp.Matrix([[exact_columns[c][r] for c in columns] for r in rows])
    # The entries are Gaussian rationals; use Gaussian rational domain
    # elimination rather than default symbolic complex zero heuristics.
    exact_nonzero = certificate.det(method="domain-ge") != 0
    assert exact_nonzero
    assert singular[-2] > .01 and singular[-1] < 1e-13
    assert all(y.det(method="domain-ge") != 0 for y in matrices)
    return {"exact_commutant_dimension":1,
            "normalized_commutator_gap":float(singular[-2]),
            "least_yukawa_singular":min(float(np.linalg.svd(y,compute_uv=False)[-1])
                                      for y in (first,second))}


up,down,e,nu = original_yukawas
families = {"quark":family_commutant(up,down), "lepton":family_commutant(e,nu)}


def full_action(q,l,zeta):
    phases = {("u","L"):q,("d","L"):q,("u","R"):zeta*q,("d","R"):zeta.conjugate()*q,
              ("e","L"):l,("nu","L"):l,("e","R"):zeta.conjugate()*l,("nu","R"):zeta*l}
    particle = np.array([phases[(name,ch)] for name,ch,_,_ in native.P])
    return np.diag(np.r_[particle,particle.conj()])


def size(x):
    return float(np.linalg.norm(x,2))


def local(h,p):
    return p*gamma-1j*gamma@(higgs(h)+majorana)


def whole(fields,p):
    result = existing["identity"].copy()
    for h,t in zip(fields,durations):
        vals,vecs=eigh(local(h,p))
        result=((vecs*np.exp(-t*vals))@vecs.conj().T)@result
    return result


# All order-two point symmetries commuting with the native gauge and
# grading, under the exact scalar Gram-commutant check above.
sign_cases = list(itertools.product((-1,1),repeat=3))
i_le = native.IDX[("e","L",0,0)]
i_re = native.IDX[("e","R",0,0)]
i_rn = native.IDX[("nu","R",0,0)]
order4 = full_action(-1j,-1j,1j)
assert np.allclose(np.diag(order4)[[i_le,i_re,i_rn]],[-1j,-1,1])
assert np.allclose(np.diag(order4@order4)[[i_le,i_re,i_rn]],[-1,1,1])
assert size(order4@majorana@order4.conj().T-majorana) < 1e-14
rows = {}
for m in (-.001,.001):
    p=.25*np.exp(-m)
    strip = whole(history,p)
    worst_sign_local = 0.
    worst_sign_whole = 0.
    for q,l,h in sign_cases:
        symmetry=full_action(q,l,complex(h))
        assert symmetry[i_re,i_re]==symmetry[i_rn,i_rn]
        transformed = [h*np.array(x) for x in history]
        worst_sign_local=max(worst_sign_local,max(
            size(symmetry@local(x,p)@symmetry.conj().T-local(y,p))
            for x,y in zip(history,transformed)))
        worst_sign_whole=max(worst_sign_whole,
                             size(symmetry@strip@symmetry.conj().T-whole(transformed,p)))
    quarter = whole([1j*np.array(x) for x in history],p)
    quarter_error=size(order4@strip@order4.conj().T-quarter)
    assert max(worst_sign_local,worst_sign_whole,quarter_error)<1e-10
    rows[str(m)]={"sign_local_covariance_error":worst_sign_local,
                  "sign_ordered_covariance_error":worst_sign_whole,
                  "order_four_ordered_covariance_error":quarter_error}

# Even dropping the Gram commutation constraints does not rescue a
# central opposite sign eR=-1, nuR=+1. Optimize the common left unitary.
ge=e@e.conj().T
gn=nu@nu.conj().T
val,vec=eigh(gn-ge)
optimal=(vec*np.sign(val))@vec.conj().T
least_squared=2*np.trace(ge+gn).real-2*np.abs(val).sum()
direct_squared=(np.linalg.norm(optimal@e+e,"fro")**2
                + np.linalg.norm(optimal@nu-nu,"fro")**2)
assert abs(direct_squared-least_squared)<1e-11
assert least_squared>0
assert size(optimal@optimal-np.eye(3))<1e-13
print(json.dumps({
    "scope":"original finite R032 diagnostic Yukawa data, gauge-commuting linear family and scalar Higgs transformations only; not QTheta identification",
    "family_checks":families,
    "central_order_two_sign_count":len(sign_cases),
    "right_e_and_right_nu_agree_for_all_signs":True,
    "order_four_example":{"H":"i","left_e":"-i","right_e":"-1","right_nu":"+1",
                          "left_e_square":"-1"},
    "two_metrics":rows,
    "least_squared_frobenius_defect_for_opposite_right_signs":float(least_squared),
    "least_joint_defect_relative_to_yukawa_pair_frobenius_norm":
        float(np.sqrt(least_squared/np.trace(ge+gn).real)),
    "full_original_QTheta_and_gap_evaluated":False,
},indent=2))
