"""Exact original Fuchs endpoints, finite direct pole controls and flat TT test."""
import importlib.util
import json
import re
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
from tfpt_constants import phi0
from v117_monodromy_weyl_a3 import M0_EXACT, U_EXACT
from v252_full_finite_triple import IDX, NH, Kswap, full_D, grading, quaternion, rep
from v420_fpole_qed_running import ME, MMU, MTAU

mp.mp.dps = 42
aps = repo / "experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review"
spec = importlib.util.spec_from_file_location("hardy_lift", aps / "check_hardy_bmu_weak_clifford_lift.py")
lift = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lift)


def unit(row, col):
    matrix = sp.zeros(3)
    matrix[row, col] = 1
    return matrix


a, b, y, delta = sp.symbols("a b y delta", real=True)
i = sp.I
g = sp.Matrix([[0, 0, -1+i], [1, i, 0], [i, 1, 0]]) / sp.sqrt(2)
sign = sp.diag(1, -1, 1)
p0 = (sp.eye(3)+U_EXACT**2)/2
f = a*unit(0,2)+b*unit(2,1)
n = -i*(a*unit(1,2)+b*unit(2,0))
fhat, nhat = (sp.simplify(g*letter*g.H) for letter in (f,n))
assert sp.simplify(g*g.H) == sp.eye(3)
assert sp.simplify(fhat-a*M0_EXACT*p0+b*p0*M0_EXACT) == sp.zeros(3)
assert sp.simplify(nhat+sign*fhat*sign) == sp.zeros(3)
assert p0 == unit(0,0)
d = y*sp.eye(3)-delta*M0_EXACT
gram_derivative = sp.simplify(sp.diff(d.H*d,delta))
assert gram_derivative == 2*delta*sp.eye(3)-y*(M0_EXACT+M0_EXACT.H)
assert sp.simplify(gram_derivative*M0_EXACT-M0_EXACT*gram_derivative) == sp.zeros(3)
fcomm = sp.simplify(M0_EXACT*fhat-fhat*M0_EXACT)
ncomm = sp.simplify(M0_EXACT*nhat-nhat*M0_EXACT)
assert fcomm.subs({a:1,b:0}) != sp.zeros(3)
assert ncomm.subs({a:1,b:0}) != sp.zeros(3)

period = (aps / "HARDY-BMU-FULL-MU01-READER.txt").read_text()
rel_period = mp.mpf(re.search(r"I21/I02=(-[0-9]+[.][0-9]+)",period)[1])
phase = mp.exp(-5j*mp.pi/12)
assert max(abs(complex(z)) for z in fcomm.subs({a:1,b:float(rel_period)})) > 1
assert max(abs(complex(z)) for z in ncomm.subs({a:1,b:float(rel_period)})) > 1
w = lift.whole96_dictionary()
creator = w.T @ np.kron(np.eye(3, dtype=int), lift.creators()[4]) @ w
gamma, swap = grading(), Kswap()
even = (np.eye(NH)+gamma)/2
left = [IDX[("e","L",0,k)] for k in range(3)]
right = [IDX[("e","R",0,k)] for k in range(3)]
charged_entries = {}
letters = {
    "single_B02": phase*unit(0,2),
    "whole_F": phase*f.subs({a:1,b:float(rel_period)}),
    "whole_N": phase*n.subs({a:1,b:float(rel_period)}),
}
for label, letter in letters.items():
    mf = np.array(letter,complex)
    coefficient = w.T @ np.kron(mf, np.eye(32)) @ w @ even
    directed = creator@coefficient-coefficient@creator
    adj = directed + directed.conj().T
    ko = adj + swap@adj.conj()@swap
    assert np.allclose(ko,ko.conj().T) and np.allclose(ko,swap@ko.conj()@swap)
    lr = ko[np.ix_(left,right)]
    assert np.allclose(lr, -(mf+mf.T))
    assert np.allclose(ko[np.ix_(right,left)],lr.conj().T)
    assert np.allclose(np.diag(lr),0)
    charged_entries[label] = {"L0_R2": str(lr[0,2]), "L1_R2": str(lr[1,2]),
                              "reverse_is_adjoint": True, "diagonal_zero_in_G_chart": True}

# Conditional direct insertion: G labels in order correspond to e, mu, tau
# of the fixed diagonal source. This is not a physical charged endpoint map.
u = mp.mpf(53)/54*phi0
m = [mp.mpf(16)/7*u**5, mp.mpf(4)/3*u**3, mp.mpf(7)/6*u**2]
e, mu, tau = (x/m[-1] for x in m)
poles = [ME/MTAU,MMU/MTAU,mp.mpf(1)]
target_q = sum(poles)/sum(mp.sqrt(x) for x in poles)**2
anchor_e, anchor_tau = mu*ME/MMU, mu*MTAU/MMU
required_square = (anchor_e**2+anchor_tau**2-e**2-tau**2)/2
det_lower_bound = abs(required_square-e*tau)
assert det_lower_bound > 9*anchor_e*anchor_tau

# Real-basis audit of the actual finite SM one-form caller. The original
# representation is family blind; the charged LR output can only be a
# complex scalar times the given Y_e, for any spectator matrices.
q0m = np.zeros((2,2),complex)
z3 = np.zeros((3,3),complex)
finite_basis = [rep(1,q0m,z3),rep(1j,q0m,z3)]
finite_basis += [rep(0,quaternion(c,d),z3) for c,d in ((1,0),(1j,0),(0,1),(0,1j))]
for rr in range(3):
    for cc in range(3):
        color_matrix = np.array(unit(rr,cc),complex)
        finite_basis += [rep(0,q0m,color_matrix),rep(0,q0m,1j*color_matrix)]
source_ye = np.diag([float(x) for x in (e,mu,tau)])
spectator = np.array([[1,2j,3],[4,5,6j],[7,8j,9]],complex)
source_d = full_D(spectator,spectator.T,source_ye,spectator.conj(),np.eye(3))
one_form_cases = 0
for bmat in finite_basis:
    bracket = source_d@bmat-bmat@source_d
    for amat in finite_basis:
        block = amat[left,:] @ bracket[:,right]
        assert np.array_equal(block,np.diag(np.diag(block)))
        assert np.allclose(block,block[2,2]*source_ye,rtol=1e-13,atol=1e-14)
        one_form_cases += 1

# Exhaust the common real strength of F and N separately, retaining original
# relative periods and phase. I2 is the sum of squared absolute 2x2 minors.
target_i1 = sum(x**2 for x in poles)
target_i2 = sum(poles[j]**2*poles[k]**2 for j,k in ((0,1),(0,2),(1,2)))
target_i3 = mp.fprod(poles)**2
target_shape2 = target_i2 / target_i1**2
target_shape3 = target_i3 / target_i1**3
moment_tests = {}
for label, (v1,v2,cos2,letter_phase) in {
    "whole_F": (mp.mpf(1),rel_period,-mp.sqrt(3)/2,phase),
    "whole_N": (rel_period,mp.mpf(1),mp.sqrt(3)/2,-1j*phase),
}.items():
    def moments(x):
        y1,y2 = v1*v1*x,v2*v2*x
        i1 = e*e+mu*mu+tau*tau+2*(y1+y2)
        i2 = ((e*mu)**2+(e*tau)**2+(mu*tau)**2 + (y1+y2)**2
              +2*y1*(mu*mu-e*tau*cos2)+2*y2*(e*e-mu*tau*cos2))
        j = mu*y1+e*y2
        i3 = (e*mu*tau)**2+j*j-2*e*mu*tau*j*cos2
        return i1,i2,i3

    for x in (mp.mpf("0.001"),mp.mpf("0.2")):
        z1,z2 = mp.sqrt(x)*letter_phase*v1, mp.sqrt(x)*letter_phase*v2
        masses = np.linalg.svd(np.array([[e,0,z1],[0,mu,z2],[z1,z2,tau]],complex),compute_uv=False)
        sq = masses*masses
        check = [sum(sq),sq[0]*sq[1]+sq[0]*sq[2]+sq[1]*sq[2],np.prod(sq)]
        assert np.allclose(moments(x),check,rtol=1e-11,atol=1e-14)

    def residual2(x):
        i1,i2,_ = moments(x)
        return i2-target_shape2*i1*i1

    constant = residual2(0)
    ca = (residual2(1)+residual2(-1))/2-constant
    cb = (residual2(1)-residual2(-1))/2
    roots = mp.polyroots([ca,cb,constant], maxsteps=200)
    positive = []
    for x in roots:
        if not isinstance(x,mp.mpf) or x < 0:
            continue
        i1,i2,i3 = moments(x)
        assert abs(i2/i1**2-target_shape2) < mp.mpf("1e-36")
        mismatch = i3/i1**3/target_shape3-1
        assert abs(mismatch) > mp.mpf("1e-4")
        zp = mp.sqrt(x)*letter_phase
        masses = sorted(mp.svd_c(mp.matrix([[e,0,v1*zp],[0,mu,v2*zp],
                                           [v1*zp,v2*zp,tau]]),compute_uv=False))
        t1,t2,t3 = masses
        positive.append({"x":x, "normalized_det_relative_residual":mismatch,
                         "mu_over_e":t2/t1,"tau_over_mu":t3/t2,
                         "Q":sum(masses)/sum(mp.sqrt(z) for z in masses)**2})
    moment_tests[label] = {"I2_match_roots":roots,"positive_root_readouts":positive}

# Flat Euclidean TT plane wave with the full ordinary R046 spinor coefficient.
q0,q3,A,B,N = sp.symbols("q0 q3 A B N",real=True)
q = [q0,0,0,q3]
h = sp.Matrix([[0,0,0,0],[0,A,B,0],[0,B,-A,0],[0,0,0,0]])
assert h*sp.Matrix(q) == sp.zeros(4,1) and sp.trace(h) == 0
q2 = q0*q0+q3*q3
riem = {}
for aa in range(4):
    for bb in range(4):
        for cc in range(4):
            for dd in range(4):
                riem[aa,bb,cc,dd] = (q[aa]*q[cc]*h[bb,dd]+q[bb]*q[dd]*h[aa,cc]
                                     -q[aa]*q[dd]*h[bb,cc]-q[bb]*q[cc]*h[aa,dd])/2
ric = sp.Matrix(4,4,lambda bb,dd: sum(riem[aa,bb,aa,dd] for aa in range(4)))
assert sp.simplify(ric-q2*h/2) == sp.zeros(4)
curvature = sp.trace(ric)
curvature_sq = sp.expand(sum(x*x for x in riem.values()))
ric_sq = sp.expand(sp.trace(ric*ric))
assert curvature == 0
assert sp.simplify(curvature_sq-4*ric_sq) == 0
a4_tt = sp.factor(N*(5*curvature**2-8*ric_sq-7*curvature_sq)/(5760*sp.pi**2))
assert sp.simplify(a4_tt+N*(A*A+B*B)*q2**2/(320*sp.pi**2)) == 0

output = {
    "same_original_transport_marks": {
        "F": "a M0 P0 - b P0 M0", "N": "-diag(1,-1,1) F diag(1,-1,1)",
        "global_real_delta_gram_tangent_commutes_with_M0": True,
        "F_commutator": fcomm, "N_commutator": ncomm,
    },
    "actual_finite_two_charged_ends":charged_entries,
    "v252_finite_one_form_real_basis_pairs":one_form_cases,
    "fixed_G_order_single_B02_test_any_phase": {
        "required_off_diagonal_norm_squared":required_square,
        "minimum_possible_absolute_det":det_lower_bound,
        "required_absolute_det":anchor_e*anchor_tau,
        "lower_bound_over_required":det_lower_bound/(anchor_e*anchor_tau),
    },
    "fixed_G_order_whole_letters_individual_strength_test":moment_tests,
    "reference_pole_Q":target_q,
    "flat_ordinary_R046_TT": {"scalar_curvature_linear":curvature,
                              "a4_TT":a4_tt,"higgs_R048_exchange_TT":0}
}
print(json.dumps(output,indent=2,default=lambda x: mp.nstr(x,24) if isinstance(x,mp.mpf) else str(x)))
