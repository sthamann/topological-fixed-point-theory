"""Use original v252, old affine Ramond packets, and the same finite time controls."""
from contextlib import redirect_stdout
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np

repo = Path(sys.argv[1]).resolve()
base = repo/"experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    x = runpy.run_path(str(base/"check_tfpt_conditioned_r113_weights_20261110.py"))
    words = runpy.run_path(str(base/"check_tfpt_ramond_null_and_native_products_20261110.py"))
    cut = runpy.run_path(str(base/"check_tfpt_spin_matched_halfwords_20261110.py"))

native, swap = x["native"],x["swap"]
Y = x["inner"](np.array([0.,1.]))
MR = x["majorana"]
I2,I3 = np.eye(2),np.eye(3)
Z2,Z3 = np.zeros((2,2)),np.zeros((3,3))
quaternions = [native.quaternion(1,0),native.quaternion(1j,0),
               native.quaternion(0,1),native.quaternion(0,1j)]
colors = [s*np.outer(I3[:,i],I3[:,j]) for s in (1,1j) for i in range(3) for j in range(3)]
sm_basis = ([native.rep(z,Z2,Z3) for z in (1,1j)]+
            [native.rep(0,q,Z3) for q in quaternions]+[native.rep(0,Z2,m) for m in colors])
# v252 rep_PS fixes a scalar on antileptons even at zero arguments.
# Isolate that scalar and take differences to span the actual real matrices.
anti_lepton = native.rep_PS(Z2,Z2,Z3)
ps_basis = ([anti_lepton]+[native.rep_PS(q,Z2,Z3)-anti_lepton for q in quaternions]+
            [native.rep_PS(Z2,q,Z3)-anti_lepton for q in quaternions]+
            [native.rep_PS(Z2,Z2,m)-anti_lepton for m in colors])
assert len(sm_basis)==24 and len(ps_basis)==27
right_leptons = [i+k*native.NP for i,row in enumerate(native.P) if row[0] in ("nu","e") and row[1]=="R" for k in range(2)]
mask = np.zeros((96,96),dtype=bool)
mask[np.ix_(right_leptons,right_leptons)] = True
particle_antiparticle = np.zeros_like(mask)
particle_antiparticle[:native.NP,native.NP:] = True
particle_antiparticle[native.NP:,:native.NP] = True
mask &= particle_antiparticle


def hat(a):
    return swap@a.conj()@swap


def comm(a,b):
    return a@b-b@a


def check_basis(basis):
    hats = list(map(hat,basis))
    assert all(not comm(a,b).any() for a in basis for b in hats)
    readouts = []
    for D in (Y,MR):
        blocks = [comm(comm(D,a),b) for a in basis for b in hats]
        assert all(not v[~mask].any() for v in blocks)
        nonzero = [v[mask] for v in blocks if v.any()]
        rank = int(np.linalg.matrix_rank(np.array(nonzero))) if nonzero else 0
        readouts.append({"number_of_nonzero_basis_pairs":len(nonzero),
                         "complex_defect_span_rank":rank,
                         "largest_HS_norm":float(max(map(np.linalg.norm,nonzero),default=0))})
    return readouts


sm,ps = check_basis(sm_basis),check_basis(ps_basis)
assert not any(t["number_of_nonzero_basis_pairs"] for t in sm)
assert ps[0]["number_of_nonzero_basis_pairs"]==0 and ps[1]["number_of_nonzero_basis_pairs"]>0
D = Y+MR
unitary_right = (I2+1j*np.array([[0,1],[1,0]]))/np.sqrt(2)
u = native.rep_PS(I2,unitary_right,I3)
u_hat = hat(u)
A1 = u@comm(D,u.conj().T)
A2 = u_hat@comm(A1,u_hat.conj().T)
transformed = u@u_hat@D@u_hat.conj().T@u.conj().T
assert np.linalg.norm(A2[~mask])<1e-14
assert np.linalg.norm(D+A1+hat(A1)+A2-transformed)<1e-13
assert np.linalg.norm(A2)>.1

sp,Q = words["sp"],words["Q"]
old,lift = words["old"],words["lift"]
m = old["m"]
q, M = lift["q"],lift["M"]
d0 = sp.Matrix([0,1,1,-1,0,0,1,0])
d1 = sp.Matrix([1,0,0,0,1,-1,0,1])
S = sp.eye(8)-q*q.T
mid = d0-q
assert S*q==-q and S*d0==d1 and d0-d1==2*q
assert [q.dot(d) for d in (d0,mid,d1)]==[2,0,-2]
assert [d.dot(d)/2 for d in (d0,mid,d1)]==[2,1,2]
assert [sp.exp(2*sp.pi*sp.I*d[0]) for d in (d0,mid,d1)]==[1,-1,1]
# The q ladder on the three normalized old states, from weight +2
# and the original adjoint. Its pi/2 rotation has endpoint overlap 1/2.
lowering = sp.Matrix([[0,0,0],[sp.sqrt(2),0,0],[0,sp.sqrt(2),0]])
rotation_generator = (lowering-lowering.T)/2
assert rotation_generator**3 == -rotation_generator
theta = sp.Symbol("theta", real=True)
rotation = sp.eye(3)+sp.sin(theta)*rotation_generator+(1-sp.cos(theta))*rotation_generator**2
assert sp.simplify(rotation[2,0]-sp.sin(theta/2)**2)==0
for root in lift["E"][:,0]-lift["E"][:,1], lift["E"][:,1]-lift["E"][:,2], lift["weak"]:
    assert root.dot(M*q)==0
for d in (d0,d1):
    p = M*d
    assert p[0]==p[1]==p[2] and p.dot(lift["weak"])==-1 and p.dot(lift["Y"])==sp.Rational(1,2)
    E = lift["E"]
    assert max(abs(p.dot(E[:,i]-E[:,j])) for i,j in ((0,1),(0,2),(1,2),(3,4)))==1
# In the old lattice representation these weights are affine highest
# under the gauge currents: every semisimple gauge root has scalar product
# 0 or +/-1 with the two momenta. Positive current modes annihilate.
V = old["V"]
reflection = [Q.from_sympy(z) for z in list(V.T*q)]
def transform_oscillators(poly):
    result = {}
    for monomial,value in poly.items():
        term = {():value}
        for n,a in monomial:
            choices = [((n,a),Q.one),((n,4),-reflection[a])] if a!=4 else [((n,4),-Q.one)]
            out = {}
            for w,coeff in term.items():
                for mark,s in choices:
                    key = tuple(sorted(w+(mark,)))
                    out[key] = out.get(key,Q.zero)+coeff*s
            term = out
        for w,z in term.items():
            result[w] = result.get(w,Q.zero)+z
    return {w:z for w,z in result.items() if z}


family = {}
for alpha,states in words["history"]["quarter_states"].items():
    delta = sp.Matrix([(j+sp.Rational(alpha))/8-sp.Rational(1,2) for j in range(8)])
    assert S*delta==delta
    drift = [Q.from_sympy(z) for z in old["cl"](V.T*delta)]
    polys = [old["old"]["vertex_zero"](old["old"]["series_response"](f,drift,4)) for f in old["polynomials"]]
    vac = old["sea"]([z<0 for z in q])
    reflected = [old["as_car"](vac,transform_oscillators(p)) for p in polys]
    assert words["gmat"](reflected)==words["gmat"](states)
    quartet = words["four_field_zero_mode"]
    direct0 = words["gmat"]([quartet(v,d0) for v in states])
    direct1 = words["gmat"]([quartet(v,d1) for v in states])
    switched0 = words["gmat"]([quartet(v,d0) for v in reflected])
    switched1 = words["gmat"]([quartet(v,d1) for v in reflected])
    assert direct0==switched1==sp.zeros(5) and direct1==switched0
    assert words["exact_rank"](direct1)==5
    projected_at_half_angle = sp.simplify(rotation[2,0].subs(theta,sp.pi/2)**2)*direct1
    assert 4*projected_at_half_angle==direct1
    family[alpha] = {"fixed_q_ranks":[0,5],"Weyl_transported_q_ranks":[5,0],
                     "entire_nonzero_5x5_Gram_preserved_exactly":True,
                     "after_old_Ramond_projection_half_angle_Gram_ratio":"1/4",
                     "P0_norm_over_original_q_norm":str(sp.simplify(direct1[0,0]/words["gmat"](states)[0,0]))}

t = sp.Rational(1,8)
probabilities = [sp.exp(-2*t*E) for E in (32,33)]
conditional_weight = sp.simplify(probabilities[0]/sum(probabilities))
assert conditional_weight != sp.Rational(1,2)
gamma = x["gamma"]
plus = np.flatnonzero(np.diag(gamma)>0)
minus = np.flatnonzero(np.diag(gamma)<0)
native_minimal = {}
for radius in (0.,.55):
    h = .25*gamma-1j*gamma@(MR+radius*Y)
    b = h[np.ix_(minus,plus)]
    rank = int(np.linalg.matrix_rank(b))
    assert rank == (3 if radius==0. else 48)
    native_minimal[str(radius)] = {"short_time_rank":rank,"minimum_total_state_dimension":len(plus)+rank}
native_h = .25*gamma-1j*gamma@(MR+.55*Y)
native_B = native_h[np.ix_(minus,plus)]
Y_B = (-1j*gamma@Y)[np.ix_(minus,plus)]
quadratic_yukawa = Y_B.conj().T@Y_B
unique_yukawa_squares = np.unique(np.round(np.linalg.eigvalsh(quadratic_yukawa),9))
assert len(unique_yukawa_squares)==12 and np.linalg.matrix_rank(quadratic_yukawa)==48
z = 2j
resolver = np.linalg.inv(z*np.eye(96)-native_h)[np.ix_(plus,plus)]
effective = np.linalg.inv((z-.25)*np.eye(48)-native_B.conj().T@native_B/(z+.25))
assert np.linalg.norm(resolver-effective)<2e-13
energy = np.linalg.eigvalsh(cut["h"])
K = cut["kl"]
shift = np.exp(-.05*energy)
asymmetry = np.linalg.norm((shift[:,None]-shift[None,:])*K)
np.testing.assert_allclose(asymmetry,.0230894685024,rtol=2e-9)
print(json.dumps({
    "full_algebra_first_order_basis":{
        "SM_basis_size":24,"SM_Y_and_MR":sm,"given_PS_basis_size":27,"given_PS_Y_and_MR":ps,
        "PS_pure_gauge_A2_HS_norm":float(np.linalg.norm(A2)),
        "PS_pure_gauge_covariance_error":float(np.linalg.norm(D+A1+hat(A1)+A2-transformed)),
        "A2_support_only_R_lepton_particle_antiparticle":True},
    "old_family_SU2_triplet":{
        "target_d0_mid_d1":[list(map(str,M*r)) for r in (d0,mid,d1)],
        "old_spin_grades":[1,-1,1],"exact_H_q_weights":[2,0,-2],
        "conformal_weight_of_triplet_all":2,
        "after_old_Ramond_projection_fixed_q_rotation_Gram":"sin(theta/2)^4 times G_d1",
        "comparison":family},
    "original_old_process_control":{
        "t":str(t),"energies_32_and_33_filter_traces":list(map(float,probabilities)),
        "conditional_mixture_P32_weight":float(conditional_weight),
        "mixture_of_separately_normalized_P32_weight":.5,
        "deterministic_trace_preserving_Euclidean_map_with_common_scalar":False},
    "native_time_and_shift_comparison":{
        "minimal_total_dimension_on_actual_native_clifford_compression":native_minimal,
        "Yukawa_quadratic_term_distinct_values":list(map(float,unique_yukawa_squares)),
        "Yukawa_quadratic_term_rank":48,
        "same_native_Schur_resolvent_error":float(np.linalg.norm(resolver-effective)),
        "same_generator_Higgs_mark_shift_HS_norm":float(asymmetry)},
    "physical_original_Padm_values":None,
},indent=2))
