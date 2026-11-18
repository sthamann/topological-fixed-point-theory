"""Test the old E8 current and the original eightfold native-CAR Higgs lift."""
from contextlib import redirect_stdout
from itertools import combinations, product
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np
import sympy as sp

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    data = runpy.run_path(str(base / "check_tfpt_conditioned_r113_weights_20261110.py"))
    old = runpy.run_path(str(base / "check_tfpt_ramond_covariance_and_cover_20261110.py"))

native = data["native"]
gamma, phi = data["gamma"], data["inner"]
tau = [(-1j*gamma@phi(np.array(x))) for x in ((1,0),(1j,0),(0,1),(0,1j))]
np.testing.assert_allclose(tau, np.asarray(tau).conj().transpose(0,2,1), atol=1e-14)
A = tau[2]+1j*tau[3]
B = A@A.conj().T-A.conj().T@A
assert not set(np.nonzero(A)[0]) & set(np.nonzero(A)[1])
# In all formulas the full 96 coefficients are matrices, not a choice of 96
# independent fermion species. The irreducible-copy test below uses particles.
mass_blocks = {}
singular = []
for field in ("u","d","e","nu"):
    ids = [native.IDX[(field,"L",0,k)] for k in range(3)]
    js = [native.IDX[(field,"R",0,k)] for k in range(3)]
    ym = phi(np.array([0,1]))[np.ix_(ids,js)]
    sm = np.linalg.svd(ym,compute_uv=False)
    singular.extend(sm)
    mass_blocks[field] = ym


def exact_complex_matrix(a):
    return sp.Matrix([[sp.Rational(float(x.real))+sp.I*sp.Rational(float(x.imag))
                       for x in row] for row in a])


exact_squares = {}
for name,matrix in mass_blocks.items():
    y = exact_complex_matrix(matrix)
    exact_squares[name] = (4*y*y.H).applyfunc(sp.expand)
z = sp.Symbol("z")
polynomials = [sp.Poly(mat.charpoly(z).as_expr(), z, domain=sp.QQ)
               for mat in exact_squares.values()]
assert all(p.gcd(p.diff()).degree()==0 and p.eval(0)!=0 for p in polynomials)
assert all(p.gcd(qp).degree()==0 for p,qp in combinations(polynomials,2))

A_spectrum = np.linalg.svd(A,compute_uv=False)
active = A_spectrum[A_spectrum>1e-12]
strengths = 2*np.array(singular)
assignment_error = max(min(abs(x-strengths)) for x in active)
assert assignment_error<2e-13
assert all(min(abs(x-active))<2e-13 for x in strengths)
pair_separation = min(np.diff(np.sort(strengths**2)))
assert len(strengths)==12 and pair_separation>1e-6

# Construct the twelve commuting Hermitian directions without unstable
# raw Vandermonde powers. Spectral interpolation of ad_B is a Lie polynomial
# in A,A*, since every monomial is an iterated commutator.
a_left, sigma, a_right = np.linalg.svd(A)
tori, eig_residuals = [], []
for strength in strengths:
    indices = np.flatnonzero(abs(sigma-strength)<1e-10)
    component = (a_left[:,indices]*sigma[indices])@a_right[indices,:]
    ad = B@component-component@B
    eig_residuals.append(np.linalg.norm(ad-2*strength**2*component))
    t = component@component.conj().T-component.conj().T@component
    tori.append(t/np.linalg.norm(t))
gram = np.array([[np.vdot(x,y) for y in tori] for x in tori])
torus_gram_error = np.linalg.norm(gram-np.eye(12))
torus_commutator_max = max(np.linalg.norm(x@y-y@x) for x,y in combinations(tori,2))
assert max(eig_residuals)<2e-12 and torus_gram_error<2e-12 and torus_commutator_max<2e-12

# Gauge-compatible particle-preserving reducing subspaces for the full
# native Higgs family. Distinct left Yukawa eigenvalues plus a connected
# graph of the other species in that basis force the generation projection
# to be scalar. Nonzero right blocks are then forced by invertibility.
def family_test(a,b):
    ya,yb = (mass_blocks[x] for x in (a,b))
    first = ya@ya.conj().T
    second = yb@yb.conj().T
    evals, vectors = np.linalg.eigh(first)
    other = vectors.conj().T@second@vectors
    offdiag = [abs(other[i,j]) for i,j in combinations(range(3),2)]
    min_gap = float(min(np.diff(evals)))
    min_link = float(min(offdiag))
    assert min_gap>1e-7 and min_link>1e-7
    x,y = exact_squares[a],exact_squares[b]
    eye3 = sp.eye(3)
    equations = (sp.kronecker_product(x,eye3)-sp.kronecker_product(eye3,x.T)).col_join(
                sp.kronecker_product(y,eye3)-sp.kronecker_product(eye3,y.T))
    exact_rank = equations.to_DM().convert_to(sp.QQ_I).rank()
    assert exact_rank==8
    return {"left_distinct_gap":min_gap,"other_min_generation_link":min_link,
            "exact_dimension_of_common_generation_commutant":9-exact_rank,
            "spectral_residual":float(np.linalg.norm(first-(vectors*evals)@vectors.conj().T))}

families = {"quarks":family_test("u","d"),"leptons":family_test("nu","e")}

fold = repo / "experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json"
M = sp.Matrix(json.loads(fold.read_text())["explicit_Weyl_matrix"])
E = sp.eye(8)
q_target = E[:,7]-E[:,6]
q = M.T*q_target
weak = E[:,3]-E[:,4]
Y = sp.Matrix([-2,-2,-2,3,3,0,0,0])/6
integer_roots = [a*E[:,i]+b*E[:,j] for i,j in combinations(range(8),2)
                 for a,b in product((-1,1),repeat=2)]
half_roots = [sp.Matrix(s)/2 for s in product((-1,1),repeat=8) if s.count(-1)%2==0]
roots = {tuple(r) for r in integer_roots+half_roots}
assert len(roots)==240 and {tuple(M*r) for r in integer_roots+half_roots}==roots


def grade(x):
    p = M.T*x
    assert len({sp.exp(2*sp.pi*sp.I*z) for z in p})==1
    return int(sp.exp(2*sp.pi*sp.I*p[0]))


doublets = []
for r0 in sorted(roots):
    r = sp.Matrix(r0)
    if not (r[0]==r[1]==r[2] and Y.dot(r)==sp.Rational(1,2) and weak.dot(r)==-1):
        continue
    upper = r+weak
    assert tuple(upper) in roots and grade(upper)==grade(r)
    assert all(tuple(r+sgn*(E[:,i]-E[:,j])) not in roots
               for i,j in combinations(range(3),2) for sgn in (-1,1))
    steps = [sgn for sgn in (-1,1) if tuple(r+sgn*q_target) in roots]
    assert all(grade(r+sgn*q_target)==-grade(r) for sgn in steps)
    doublets.append({"lower_root_target":list(map(str,r)),
                     "upper_root_target":list(map(str,upper)),
                     "lower_root_old_charge":list(map(str,M.T*r)),
                     "Xi":grade(r),"q_root_steps":steps})
assert len(doublets)==10
summary = {str(sign): {"count":sum(x["Xi"]==sign for x in doublets),
                      "nontrivial_q_steps":sum(x["Xi"]==sign and bool(x["q_root_steps"]) for x in doublets)}
           for sign in (1,-1)}

# Only the four integer-current doublets are used here. These are complete
# old current controls with both q ends and a paired weak-neutral component,
# not replacements for the native Higgs force excluded above.
Q, levels, bit = old["Q"], old["levels"], old["bit"]


def zero_mode(state, charge, ramond=True):
    (i,si),(j,sj) = [(k,int(x)) for k,x in enumerate(charge) if x]
    out = {}
    for n in levels:
        other = -si*sj*n - int(not ramond and si==sj)
        if other not in levels:
            continue
        pairs = ((bit(j,other),sj),(bit(i,n),si))
        for bits, coefficient in state.items():
            for k, create in pairs:
                if bool(bits>>k&1) == (create>0):
                    break
                if (bits&((1<<k)-1)).bit_count()%2:
                    coefficient = -coefficient
                bits ^= 1<<k
            else:
                out[bits] = out.get(bits,Q.zero)+coefficient
    return {k:v for k,v in out.items() if v}


def inner(left,right):
    return Q.to_sympy(sum((old["m"].conj(x)*right.get(k,Q.zero) for k,x in left.items()),Q.zero))


controls = {}
for alpha in (sp.Rational(1,4),sp.Rational(3,4)):
    delta = sp.Matrix([(j+alpha)/8-sp.Rational(1,2) for j in range(8)])
    drift = [Q.from_sympy(x) for x in old["cl"](old["V"].T*delta)]
    inputs = [old["old"]["series_response"](f,drift,4) for f in old["polynomials"]]
    ns = [old["as_car"](old["vac_ns"],old["old"]["input_word_map"](f)) for f in inputs]
    packets = [old["as_car"](old["vac_r"],old["old"]["vertex_zero"](f)) for f in inputs]
    control = []
    for cand in (x for x in doublets if x["Xi"]==1):
        charge = sp.Matrix(cand["lower_root_old_charge"])
        images = [zero_mode(x,charge) for x in packets]
        opposite = [zero_mode(x,-charge) for x in packets]
        gram = sp.Matrix(5,5,lambda i,j:inner(images[i],images[j]))
        reverse = sp.Matrix(5,5,lambda i,j:inner(opposite[i],opposite[j]))
        original = sp.Matrix(5,5,lambda i,j:inner(packets[i],packets[j]))
        assert gram-reverse == -q.dot(charge)*original
        shift = 8*delta.dot(charge)
        ratios = [sp.simplify(gram[i,i]/original[i,i]) for i in range(5)]
        ns_images = [zero_mode(x,charge,ramond=False) for x in ns]
        ns_ratios = [sp.simplify(inner(y,y)/inner(x,x)) for x,y in zip(ns,ns_images,strict=True)]
        time = sp.Rational(1,8)
        control.append({"old_charge":list(map(str,charge)),
                        "balanced_gram_rank":gram.rank(),
                        "original_H_intermediate_energy":str(32+shift),
                        "exact_norm_ratios":list(map(str,ratios)),
                        "exact_ns_control_norm_ratios":list(map(str,ns_ratios)),
                        "closed_q_minus_ns_at_old_time_1_over_8": [float((x-y)*sp.exp(-time*shift)) for x,y in zip(ratios,ns_ratios,strict=True)],
                        "ratios_at_old_time_1_over_8":[float(x*sp.exp(-time*shift)) for x in ratios]})
    controls[str(alpha)] = control

print(json.dumps({
    "old_neutral_q":list(map(str,q)),
    "native_Higgs_force_commutators": {
        "pair_singular_values_squared_by_species": {
            f:[float(4*y*y) for y in np.linalg.svd(m,compute_uv=False)]
            for f,m in mass_blocks.items()},
        "smallest_distinct_gap":float(pair_separation),
        "exact_binary_input_yukawa_characteristic_polynomials_squarefree_and_pairwise_coprime":True,
        "ad_eigenvector_max_residual":float(max(eig_residuals)),
        "torus_Gram_residual":float(torus_gram_error),
        "torus_commutator_max":float(torus_commutator_max),
        "native_independent_commuting_directions":len(tori),
        "old_E8_compact_rank":8,
    },
    "fixed_particle_fold_and_spin_lift": {
        "left_reducing_subspace_test":families,
        "minimal_particle_reducing_rank_preserving_weak_color_Y_and_all_Higgs":12,
        "full_particle_copy_rank":48,
        "old_q_weight":1,
        "minimum_weight_for_same_q_on_nonzero_reducing_particle_copies_per_residue":12,
        "full_copy_q_weight":48,
        "scope":"Same diagonal half-charge q on a reducing particle subspace of every original residue. No identification with physical R113 selector.",
    },
    "old_E8_color_singlet_weight_one_Higgs_doublet_inventory": {
        "counts_by_Xi":summary,"lower_and_upper_root_pairs":doublets},
    "actual_old_balanced_q_current_controls":controls,
    "physical_balanced_R113_Ramond_Higgs_value":None,
},indent=2))
