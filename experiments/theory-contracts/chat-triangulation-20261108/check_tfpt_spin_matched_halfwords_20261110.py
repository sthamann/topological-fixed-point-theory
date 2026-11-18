"""Original Ramond graph lift, neutral weak control, and native marked Higgs Cauchy trace."""
from contextlib import redirect_stdout
from itertools import combinations
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import log_expit

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    old = runpy.run_path(str(base / "check_tfpt_ramond_covariance_and_cover_20261110.py"))
    native = runpy.run_path(str(base / "check_tfpt_conditioned_r113_weights_20261110.py"))

sp, Q, conj = old["sp"], old["Q"], old["m"].conj
M, q = old["M"], old["q"]
d = M.T * (sp.eye(8)[:,3]-sp.eye(8)[:,4])
assert d.dot(d)==2 and d.dot(q)==0 and q.dot(q)==2
assert all(sp.exp(2*sp.pi*sp.I*x)==-1 for x in q)
assert all(sp.exp(2*sp.pi*sp.I*x)==1 for x in d)
color_adjoint = {}
g_color = (-1,-1,1)
for i,j in combinations(range(3),2):
    pulled = M.T*(sp.eye(8)[:,i]-sp.eye(8)[:,j])
    signs = {sp.exp(2*sp.pi*sp.I*x) for x in pulled}
    assert signs == {g_color[i]*g_color[j]}
    assert pulled.dot(q)==0
    color_adjoint[f"{i+1}-{j+1}"] = int(next(iter(signs)))
assert list(color_adjoint.values())==[1,-1,-1]
raise_channel = list(d).index(1)
lower_channel = list(d).index(-1)
bit = old["bit"]
levels = old["levels"]


def hop(states, dest, src):
    out = {}
    for bits, coeff in states.items():
        for n in levels:
            a, b = bit(dest,n), bit(src,n)
            if bits>>a&1 or not bits>>b&1:
                continue
            sign = (bits&((1<<b)-1)).bit_count()
            rem = bits^(1<<b)
            sign += (rem&((1<<a)-1)).bit_count()
            new = rem^(1<<a)
            out[new] = out.get(new,Q.zero)+((-1)**sign)*coeff
    return {k:v for k,v in out.items() if v}


def inner(a,b):
    return Q.to_sympy(sum((conj(v)*b.get(k,Q.zero) for k,v in a.items()),Q.zero))


def normalized_gram(states):
    gram = sp.Matrix(5,5,lambda i,j:inner(states[i],states[j]))
    scale = sp.diag(*[1/sp.sqrt(x) for x in old["norms"]])
    return (scale*gram*scale).applyfunc(sp.simplify)


def encoded(mat):
    return [[str(x) for x in row] for row in mat.tolist()]


ramond = {}
for alpha in (sp.Rational(1,4),sp.Rational(3,4)):
    delta = sp.Matrix([(j+alpha)/8-sp.Rational(1,2) for j in range(8)])
    assert delta.dot(q)==0
    drift = [Q.from_sympy(x) for x in old["cl"](old["V"].T*delta)]
    prepared = [old["old"]["series_response"](f,drift,4) for f in old["polynomials"]]
    packets = [old["as_car"](old["vac_r"],old["old"]["vertex_zero"](f)) for f in prepared]
    g = normalized_gram(packets)
    assert g==sp.Matrix(old["original"][str(alpha)]["full_return_Gram"])
    up = normalized_gram([hop(x,raise_channel,lower_channel) for x in packets])
    down = normalized_gram([hop(x,lower_channel,raise_channel) for x in packets])
    assert up==down  # Cartan charge d.q vanishes on all five packets.
    rank = up.rank()
    assert rank>0
    ramond[str(alpha)] = {
        "exact_ramond_graph_gram_rank":g.rank(),
        "ramond_energy_both_faces_old_H":32,
        "bare_same_spin_collar_with_one_ramond_and_even_higgs_exact_zero":True,
        "old_weak_root_shift_in_H_units":str(8*delta.dot(d)),
        "balanced_old_weak_control_gram_rank":rank,
        "balanced_old_weak_control_gram":encoded(up),
        "balanced_weak_norm_over_q_norm":[str(sp.simplify(up[i,i]/g[i,i])) for i in range(5)],
        "balanced_weak_norm_over_q_norm_float":[float((up[i,i]/g[i,i]).evalf()) for i in range(5)],
    }

gamma, mr, make_phi = native["gamma"], native["majorana"], native["inner"]
weak_lift = native["weak_lift"]
p, ell = .25, 2*np.pi
background = .55*np.array([0.,1.])
mark = np.array([0.,1.])
g = np.array([[1,1j],[1j,1]],dtype=complex)/np.sqrt(2)
g96 = weak_lift(g)
nodes, weights = leggauss(96)


def higgs_cut(field, direction, first, last):
    h = p*gamma-1j*gamma@(mr+make_phi(field))
    v = -1j*gamma@make_phi(direction)
    np.testing.assert_allclose(h,h.conj().T,atol=2e-14)
    np.testing.assert_allclose(v,v.conj().T,atol=2e-14)
    e, u = np.linalg.eigh(h)
    s_log = .5*log_expit(2*ell*e)
    c_log = .5*log_expit(-2*ell*e)
    s,c = np.exp(s_log),np.exp(c_log)
    w = np.vstack((u*s,u*c))
    perpendicular = np.vstack((-u*c,u*s))
    a,b = ell*first,ell*last
    r = (a+b)/2+(b-a)*nodes/2
    left = np.exp(s_log[:,None]-e[:,None]*(ell-r)[None,:])
    right = np.exp(s_log[:,None]-e[:,None]*r[None,:])
    kernel = -((u.conj().T@v@u)*((left*(b-a)*weights/2)@right.T))
    tangent = perpendicular@kernel@w.conj().T
    return h, v, tangent, kernel, w


h, v, mark_left, kl, w = higgs_cut(background,mark,1/8,1/4)
_, _, mark_right, kr, _ = higgs_cut(background,mark,3/4,7/8)
np.testing.assert_allclose(kr,kl.conj().T,atol=2e-12)
h_rot, v_rot, mark_rot, _, _ = higgs_cut(g@background,g@mark,1/8,1/4)
both_faces = np.kron(np.eye(2),g96)
np.testing.assert_allclose(h_rot,g96@h@g96.conj().T,atol=2e-13)
np.testing.assert_allclose(v_rot,g96@v@g96.conj().T,atol=2e-13)
cov_defect = np.linalg.norm(mark_rot-both_faces@mark_left@both_faces.conj().T)
assert cov_defect<1e-10
_, _, whole, _, _ = higgs_cut(background,mark,0.,1.)
pieces = sum((higgs_cut(background,mark,a,b)[2] for a,b in ((0.,1/8),(1/8,1/4),(1/4,3/4),(3/4,7/8),(7/8,1.))),np.zeros_like(whole))
assert np.linalg.norm(whole-pieces)<2e-11


def projection(phi):
    e, u = np.linalg.eigh(p*gamma-1j*gamma@(mr+make_phi(phi)))
    w = np.vstack((u*np.exp(.5*log_expit(2*ell*e)),u*np.exp(.5*log_expit(-2*ell*e))))
    return w@w.conj().T


dx = 1e-5
numerical = (projection(background+dx*mark)-projection(background-dx*mark))/(2*dx)
derivative_defect = np.linalg.norm(numerical-whole-whole.conj().T)
assert derivative_defect<2e-8
assert np.linalg.norm(w.conj().T@mark_left@w)<2e-12
print(json.dumps({
    "original_q":list(map(str,q)), "original_weak_root_pullback":list(map(str,d)),
    "old_spin_character_on_actual_color_roots":color_adjoint,
    "inner_color_representative_diagonal":list(g_color),
    "old_ramond_graph_and_balanced_weak_controls":ramond,
    "original_native_conditional_higgs_poisson_mark": {
        "quarter_p":p,"radial_background":.55,"collar_half_width":ell,
        "higgs_window_fraction":[1/8,1/4],
        "reflected_window_fraction":[3/4,7/8],
        "native_tangent_hilbert_schmidt_norm":float(np.linalg.norm(mark_left)),
        "native_tangent_reflected_hilbert_schmidt_norm":float(np.linalg.norm(mark_right)),
        "conditional_tangent_gram_trace":float(np.linalg.norm(mark_left)**2),
        "constant_whole_projector_derivative_fd_error":float(derivative_defect),
        "native_weak_covariance_error":float(cov_defect),
        "scope":"Two-face Cauchy trace of original marked coefficients, fixed background. No physical Padm assignment, determinant integration or embedding of old q."},
    "physical_full_Padm_mixed_response":None,
},indent=2))
