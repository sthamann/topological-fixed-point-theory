"""Exact family curvature and independent finite-loop tests on original v252 Cauchy data."""
from contextlib import redirect_stdout
from itertools import combinations
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np
from scipy.linalg import expm_frechet, polar
from scipy.special import expit
import sympy as sp

repo = Path(sys.argv[1]).resolve()
base = repo/"experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    original = runpy.run_path(str(base/"check_tfpt_conditioned_r113_weights_20261110.py"))

rep, gamma, mass, inner, gauge = (original[k] for k in ("native","gamma","majorana","inner","weak_lift"))
fields = ((1,0),(1j,0),(0,1),(0,1j))
ell = 2*np.pi
prime = 65537
imaginary = 256
assert imaginary**2%prime==prime-1


def exact(M):
    return sp.Matrix([[sp.Rational(float(v.real))+sp.I*sp.Rational(float(v.imag))
                       for v in row] for row in M])


def modular(M):
    result = []
    for row in M.applyfunc(sp.expand).tolist():
        values = []
        for x in row:
            a,b = x.as_real_imag()
            assert a.is_Rational and b.is_Rational
            values.append((int(a.p)*pow(int(a.q),-1,prime)
                           +imaginary*int(b.p)*pow(int(b.q),-1,prime))%prime)
        result.append(values)
    return np.array(result,dtype=np.int64)


def lie_rank(x,y):
    I = sp.eye(3)
    seeds = [modular(sp.I*(a-sp.trace(a)*I/3)) for a in (x,y)]
    basis, echelon = [],{}


    def add(matrix):
        row = matrix.ravel().copy()
        for pivot,old in sorted(echelon.items()):
            row = (row-row[pivot]*old)%prime
        nz = np.flatnonzero(row)
        if not len(nz):
            return
        k = nz[0]
        echelon[k] = (row*pow(int(row[k]),-1,prime))%prime
        basis.append(matrix)


    for seed in seeds:
        add(seed)
    rounds = [len(basis)]
    while True:
        for A,B in combinations(tuple(basis),2):
            add((A@B-B@A)%prime)
        if rounds[-1]==len(basis):
            return rounds
        rounds.append(len(basis))


def family_system(up,down):
    modes = [rep.IDX[(f,s,0,g)] for s in ("L","R") for f in (up,down) for g in range(3)]
    ix = np.ix_(modes,modes)
    if up=="u":
        assert np.linalg.norm(mass[ix])==0
    grad = gamma[ix]
    H = [(-1j*gamma@inner(np.array(h)))[ix] for h in fields]
    assert all(np.linalg.norm(grad@a+a@grad)<1e-14 for a in H)
    yf = inner(np.array([0,1]))[ix]
    y1,y2 = exact(yf[:3,6:9]),exact(yf[3:6,9:12])
    A,B = y1@y1.H,y2@y2.H
    fa = H[0]@H[1]-H[1]@H[0]
    fb = H[2]@H[3]-H[3]@H[2]
    np.testing.assert_allclose(fa[:3,:3],-2j*np.array(B.evalf(),complex),atol=1e-13)
    np.testing.assert_allclose(fb[:3,:3],2j*np.array(A.evalf(),complex),atol=1e-13)
    ranks=lie_rank(A,B)
    assert ranks[-1]==8
    commutator = (A@B-B@A).applyfunc(sp.expand)
    cubic = sp.expand(sp.trace((commutator**2).applyfunc(sp.expand)@commutator))
    assert sp.re(cubic)==0 and sp.im(cubic)!=0
    return modes,grad,H,fa,fb,ranks,float(sp.im(cubic))


modes,G,H,FA,FB,ranks,quark_cubic = family_system("u","d")
# The leptonic generation matrices use the original Yukawa coefficients;
# the clean Cauchy curvature formula below is confined to quarks,
# where Majorana vanishes exactly.
_,_,_,_,_,lepton_ranks,lepton_cubic = family_system("nu","e")


def cauchy_frame(h):
    E,U = np.linalg.eigh(h)
    return np.vstack(((U*np.sqrt(expit(2*ell*E)))@U.conj().T,
                      (U*np.sqrt(expit(-2*ell*E)))@U.conj().T))


def make_vertices(i,j,width,subdivisions):
    coords = np.zeros((5,4))
    coords[1:3,i] = width
    coords[2:4,j] = width
    return np.vstack([np.linspace(a,b,subdivisions,endpoint=False)
                      for a,b in zip(coords[:-1],coords[1:],strict=True)]
                     +[coords[-1:]])


def wilson(p,i,j,width,subdivisions,change_frames=False,raw_projection=False):
    vertices = make_vertices(i,j,width,subdivisions)
    frames = [cauchy_frame(p*G+sum(c*T for c,T in zip(v,H,strict=True))) for v in vertices]
    total = np.eye(len(G),dtype=complex)
    for k,(A,B) in enumerate(zip(frames[:-1],frames[1:],strict=True)):
        if change_frames:
            x,y = [0.2*np.sin(2*np.pi*m/(len(frames)-1)) for m in (k,k+1)]
            g1,g2 = [gauge(np.array([[np.cos(a),np.sin(a)],[-np.sin(a),np.cos(a)]]))[np.ix_(modes,modes)]
                     for a in (x,y)]
            pair1,pair2,link = [np.kron(np.eye(2),q) for q in (g1,g2,g2@g1.conj().T)]
            M = (pair2@B).conj().T@link@(pair1@A)
        else:
            M = B.conj().T@A
        total = (M if raw_projection else polar(M)[0])@total
    return total


quarter = {}
for p in (.25,.75):
    c = (np.tanh(ell*p)/(2*p))**2
    target_A,target_B = c*FA,c*FB
    target_com = target_A@target_B-target_B@target_A
    assert np.linalg.norm(target_com[:3,:3])>1
    loops={}
    for width in (.02,.01,.005):
        U,V = [wilson(p,i,j,width,16) for i,j in ((0,1),(2,3))]
        identity = np.eye(len(G))
        for loop,target in ((U,target_A),(V,target_B)):
            np.testing.assert_allclose(loop.conj().T@loop,identity,atol=1e-12)
            assert np.linalg.norm((loop-identity)/width**2+target)/np.linalg.norm(target)<.06
        commutator = U@V@U.conj().T@V.conj().T
        compare = np.linalg.norm((commutator-identity)/width**4-target_com)/np.linalg.norm(target_com)
        assert compare<.10 and abs(np.linalg.det(commutator)-1)<1e-11
        loops[str(width)] = {
            "group_commutator_distance_from_identity":float(np.linalg.norm(commutator-identity)),
            "determinant_residual":float(abs(np.linalg.det(commutator)-1)),
            "relative_curvature_commutator_error":float(compare),
        }
    original_loop = wilson(p,2,3,.01,16)
    same_loop_gauge = wilson(p,2,3,.01,16,True)
    assert np.linalg.norm(original_loop-same_loop_gauge)<1e-11
    np.testing.assert_allclose(wilson(p,2,3,.01,32),original_loop,atol=3e-9)
    D = .2
    Ua,Ub = [expm_frechet(-D*p*G,-D*X,compute_expm=False) for X in H[2:]]
    actual = Ub@Ua-Ua@Ub
    time_c=(np.sinh(D*p)/p)**2
    np.testing.assert_allclose(actual,-time_c*FB,atol=2e-13)
    quarter[str(p)] = {
        "geometric_curvature_scalar":float(c),
        "open_ordered_mixed_time_scalar_at_window_width_point_2":float(time_c),
        "window_width_that_matches_geometric_scalar_at_this_one_p":float(np.arcsinh(.5*np.tanh(ell*p))/p),
        "actual_up_family_curvature_commutator_HS":float(np.linalg.norm(target_com[:3,:3])),
        "passive_weak_gauge_transport_residual":float(np.linalg.norm(original_loop-same_loop_gauge)),
        "geometric_loop_tests":loops,
    }
assert quarter["0.25"]["window_width_that_matches_geometric_scalar_at_this_one_p"]>quarter["0.75"]["window_width_that_matches_geometric_scalar_at_this_one_p"]

# A chain of repeated projections onto the instantaneous Cauchy graph
# converges to the geometric transport, even with no physical generator.
# This is a different operation from the original whole strip.
projection_chain={}
for p in (.25,.75):
    errors=[]
    for subdivisions in (8,16,32):
        U=wilson(p,2,3,.03,subdivisions)
        raw=wilson(p,2,3,.03,subdivisions,raw_projection=True)
        errors.append(float(np.linalg.norm(U-raw)))
    assert all(1.9 < a/b < 2.1 for a,b in zip(errors[:-1],errors[1:],strict=True))
    projection_chain[str(p)]={"subdivisions":[8,16,32],
                              "difference_from_geometric_unitary":errors}


# Extend the original quark graph to its actual Fourier circle and
# compare one uniform lower-Higgs mark with a first Fourier harmonic.
# This is a differential of the same graph projection, not a new H.
def angle(E):
    return np.arctan(np.exp(-ell*np.clip(E,-100,100)))


def tangent_factor(E,F):
    if E==F:
        return -ell/(2*np.cosh(ell*E)) if abs(ell*E)<700 else 0.
    return np.sin(angle(E)-angle(F))/(E-F)


def first_fourier_curvature(pm,pn):
    t = tangent_factor(pm,-pn)
    a,b = tangent_factor(pm,-pm),tangent_factor(pn,-pn)
    return .25*t*(a+b)*FB[:3,:3]


def finite_projector_check(p):
    # Real harmonic cos x couples the two adjacent momenta with 1/2.
    A = np.kron(np.eye(2),H[2])
    B = np.kron(np.array([[0.,.5],[.5,0.]]),H[3])
    h = np.kron(np.diag([p,p+1]),G)
    W = cauchy_frame(h)
    def proj(a):
        frame=cauchy_frame(a)
        return frame@frame.conj().T
    width=1e-5
    dA=(proj(h+width*A)-proj(h-width*A))/(2*width)
    dB=(proj(h+width*B)-proj(h-width*B))/(2*width)
    F=W.conj().T@(dA@dB-dB@dA)@W
    cross=F[12:15,:3]
    direct=first_fourier_curvature(p+1,p)
    np.testing.assert_allclose(cross,direct,atol=1e-8,rtol=5e-8)
    return float(np.linalg.norm(cross-direct))


def relative_metric_response(p):
    x=ell*abs(p)
    return 4*(1-2*x/np.sinh(2*x)) if x<300 else 4.


def traced_geometric_intensity(alpha,cutoff):
    p=np.arange(-cutoff,cutoff)+alpha
    c=(np.tanh(ell*p)/(2*p))**2
    squared=c*c
    D=np.array([relative_metric_response(k) for k in p])
    return float(squared.sum()),float(squared@D)


spatial={}
for alpha in (.25,.75):
    first=first_fourier_curvature(alpha+1,alpha)
    second=first_fourier_curvature(alpha,alpha-1)
    cartan_commutator = second@second.conj().T-first.conj().T@first
    assert np.linalg.norm(cartan_commutator)>1e-2
    for n in (100,200):
        pm=alpha+n
        F=first_fourier_curvature(pm,pm-1)
        asym=np.linalg.norm(pm**2*F-FB[:3,:3]/8)/np.linalg.norm(FB[:3,:3]/8)
        assert asym<.02
    energy_intensity,metric=traced_geometric_intensity(alpha,10000)
    bound=1/(24*(10000-1)**3)
    spatial[str(alpha)]={
        "first_harmonic_curvature_norm":float(np.linalg.norm(first)),
        "first_harmonic_commutator_noncentral_at_anchor_HS":float(np.linalg.norm(cartan_commutator)),
        "independent_projector_derivative_error":finite_projector_check(alpha),
        "full_tower_geometric_intensity_scalar":energy_intensity,
        "log_radius_derivative_of_intensity_scalar":metric,
        "intensity_tail_upper_bound":bound,
    }
np.testing.assert_allclose(
    [spatial["0.25"][k] for k in ("full_tower_geometric_intensity_scalar","log_radius_derivative_of_intensity_scalar")],
    [spatial["0.75"][k] for k in ("full_tower_geometric_intensity_scalar","log_radius_derivative_of_intensity_scalar")],
    atol=1e-12,rtol=0)


def full_geometric_lie_rank():
    charges = [exact(T) for T in H]
    grading = exact(G)
    generators = [modular(A@B-B@A) for A,B in combinations(charges,2)]
    generators += [modular(grading@A-A@grading) for A in charges]
    basis, echelon = [],{}


    def add(matrix):
        row = matrix.ravel().copy()
        for pivot,old in sorted(echelon.items()):
            row = (row-row[pivot]*old)%prime
        nz = np.flatnonzero(row)
        if len(nz):
            k = nz[0]
            echelon[k] = (row*pow(int(row[k]),-1,prime))%prime
            basis.append(matrix)


    for generator in generators:
        add(generator)
    initial = len(basis)
    i = 0
    while i<len(basis):
        a = basis[i]
        for b in generators:
            add((a@b-b@a)%prime)
        i+=1
    return initial,len(basis)


def mixed_geometry_check(p):
    # The metric variable changes p->p*exp(-m) with the original ell fixed.
    # This is exactly the radius convention used in the earlier C58 test.
    W = cauchy_frame(p*G)
    def projection(h):
        V=cauchy_frame(h)
        return V@V.conj().T
    eps=1e-5
    dM=(projection(p*np.exp(-eps)*G)-projection(p*np.exp(eps)*G))/(2*eps)
    dH=(projection(p*G+eps*H[2])-projection(p*G-eps*H[2]))/(2*eps)
    F=W.conj().T@(dM@dH-dH@dM)@W
    coef = -ell*np.tanh(ell*p)/(4*np.cosh(ell*p))
    expected = coef*(G@H[2]-H[2]@G)
    np.testing.assert_allclose(F,expected,atol=2e-9,rtol=1e-8)
    assert np.linalg.norm(G@expected+expected@G)<1e-13
    # The lower complex direction has the original weak charge. It
    # cannot itself be the SM-neutral old family root.
    charged=coef*(G@(H[2]-1j*H[3])-(H[2]-1j*H[3])@G)
    T3 = np.diag([.5]*3+[-.5]*3+[0.]*6)
    assert np.linalg.norm(T3@charged-charged@T3+charged/2)<1e-13
    return {"odd_metric_Higgs_curvature_scalar":float(coef),
            "odd_metric_Higgs_curvature_norm":float(np.linalg.norm(expected)),
            "relative_weak_charge_of_lower_complex_component":"-1/2",
            "independent_projector_finite_difference_residual":float(np.linalg.norm(F-expected))}


initial_rank, full_rank = full_geometric_lie_rank()
assert full_rank==143
# Two neutral geometric generators on the original quark weak doublet.
# The sum cancels weak charge on the left; the difference is purely
# left weak Cartan. Two brackets with the difference restore neutrality.
S = FA+FB
D = FB-FA
comm1 = D@S-S@D
R = D@comm1-comm1@D
np.testing.assert_allclose(S[:3,:3],S[3:6,3:6],atol=1e-13)
np.testing.assert_allclose(R[:3,:3],R[3:6,3:6],atol=1e-11)
np.testing.assert_allclose(R[6:,:],0.,atol=1e-11)
weak_raise = np.zeros((12,12))
weak_raise[:3,3:6]=np.eye(3)
Y_quark = np.diag([1/6]*6+[2/3]*3+[-1/3]*3)
assert max(np.linalg.norm(a@b-b@a) for a in (S,R) for b in (weak_raise,Y_quark))<1e-11
# Compute the modular rank from the dyadic inputs, with no rounded
# intermediate matrices in this certificate.
SA,SB = [exact(H[i])@exact(H[j])-exact(H[j])@exact(H[i]) for i,j in ((0,1),(2,3))]
SE,DE = (SA+SB).applyfunc(sp.expand),(SB-SA).applyfunc(sp.expand)
E1=(DE@SE-SE@DE).applyfunc(sp.expand)
RE=(DE@E1-E1@DE).applyfunc(sp.expand)
neutral_ranks=lie_rank(-sp.I*SE[:3,:3],-sp.I*RE[:3,:3])
assert neutral_ranks[-1]==8
massive_quark= .25*G+.55*H[2]
neutral_massive_commutator = float(np.linalg.norm(massive_quark@R-R@massive_quark))
assert neutral_massive_commutator>1

mixed={str(p):mixed_geometry_check(p) for p in (.25,.75)}

print(json.dumps({
    "quark_traceless_family_lie_closure_dimensions":ranks,
    "lepton_Yukawa_square_traceless_lie_closure_dimensions":lepton_ranks,
    "exact_modular_certificate_prime":prime,
    "diagnostic_input_imaginary_trace_of_cubed_Yukawa_square_commutator":{
        "quark":quark_cubic,"lepton":lepton_cubic,
        "both_nonzero_exactly":True,
        "scope":"Provided diagnostic Yukawas; no prediction of physical CP phase.",
    },
    "quark_Cauchy_curvature":"tanh(ell*p)^2/(4*p^2) [T_a,T_b]",
    "quarters":quarter,
    "no_proper_flatly_parallel_generation_projector_on_this_up_family_fiber":True,
    "original_local_Fourier_harmonic_test":spatial,
    "first_harmonic_curvature_asymptotic":"p^2 F_(p,p-1) -> [T_Re,T_Im]/8 on the quark up family",
    "finite_Fourier_curvature_is_trace_class_on_full_original_circle":True,
    "curvature_pair_fails_affine_Cartan_by_nonzero_one_particle_bracket":True,
    "original_mixed_metric_Higgs_curvature":mixed,
    "mixed_Higgs_and_metric_geometric_lie_closure":{
        "initial_curvature_span_mod_prime":initial_rank,
        "full_Lie_span_mod_prime":full_rank,
        "maximum_possible_traceless_rank":143,
    },
"gauge_neutral_quark_family_geometry": {
        "sum_of_curvatures_and_double_bracket_commute_with_weak_and_hypercharge":True,
        "left_traceless_lie_dimensions_mod_prime":neutral_ranks,
        "actual_double_bracket_commutator_with_fixed_nonzero_Higgs_HS":neutral_massive_commutator,
    },
    "repeated_instantaneous_projection_without_a_physical_generator":projection_chain,
    "physical_Padm_and_weighted_source":None,
},indent=2))
