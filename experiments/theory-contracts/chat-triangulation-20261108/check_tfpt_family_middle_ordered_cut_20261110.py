"""Original old weight-two triplet and ordered native Higgs diagnostics."""
from contextlib import redirect_stdout
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    old = runpy.run_path(str(base / "check_ramond_full_return_20261110.py"))
    native = runpy.run_path(str(base / "check_tfpt_conditioned_r113_weights_20261110.py"))

sp, Q = old["sp"], old["Q"]
q, M, V = old["q"], old["M"], old["V"]
d0 = sp.Matrix([0,1,1,-1,0,0,1,0])
d1 = sp.Matrix([1,0,0,0,1,-1,0,1])
r = d0-q
assert d1 == r-q and q.dot(r)==0 and q.dot(q)==r.dot(r)==2
axes = sp.Matrix.hstack(V,q,r)
metric = [[Q.from_sympy(x) for x in row] for row in (axes.H@axes).applyfunc(sp.simplify).tolist()]
generator = {tuple(q):((4,1),),tuple(-q):((4,-1),),tuple(r):((5,1),),
             tuple(d0):((4,1),(5,1)),tuple(d1):((4,-1),(5,1))}
polys, norms = old["m"].code()


def combine(terms):
    result = {}
    for word, value in terms:
        result[word] = result.get(word,Q.zero)+value
    return {word:value for word,value in result.items() if value}


def grade(word):
    return sum(n for n,_ in word)


def exponential(charge, order):
    if order==0:
        yield (),Q.one
        return
    options = [(n,k,c) for n in range(1,order+1) for k,c in generator[tuple(charge)]]
    for counts in product(*(range(order//n+1) for n,_,_ in options)):
        if sum(count*n for count,(n,_,_) in zip(counts,options,strict=True))!=order:
            continue
        word, value = [],Q.one
        for count,(n,k,c) in zip(counts,options,strict=True):
            word.extend([(n,k)]*count)
            value *= Q.from_sympy(sp.Rational(c**count,n**count*factorial(count)))
        yield tuple(word),value


def vertex(poly, charge, momentum, energy_step=0):
    # Coefficient of z^(energy_step-|charge|^2/2) in the lattice vertex.
    # The omitted cocycle is constant for each fixed input and output momentum.
    dot = [Q.from_sympy(x) for x in (axes.T@charge).applyfunc(sp.simplify)]
    shift = int(energy_step-charge.dot(charge)/2-charge.dot(momentum))
    terms = []
    for word,value in poly.items():
        for n in range(len(word)+1):
            for deleted in combinations(range(len(word)),n):
                removed = sum(word[i][0] for i in deleted)
                if removed+shift<0:
                    continue
                kept = tuple(word[i] for i in range(len(word)) if i not in deleted)
                c = value
                for i in deleted:
                    c *= -dot[word[i][1]]
                if c:
                    terms.extend((tuple(sorted(kept+new)),c*a)
                                 for new,a in exponential(charge,removed+shift))
    return combine(terms)


def middle(poly):
    # (q(-1)|r>)_0 = (:J_q V_r:)_0. q.r=0 makes the two modes commute.
    terms = [(word,2*a) for word,a in vertex(poly,r,q).items()]
    for n in range(1,4):
        cut = combine((w[:i]+w[i+1:],n*metric[4][k]*a)
                      for w,a in poly.items() for i,(degree,k) in enumerate(w)
                      if degree==n)
        # q is real. metric[4][k] is the bilinear contraction q dot axis_k.
        terms.extend(vertex(cut,r,q,n).items())
        terms.extend((tuple(sorted(w+((n,4),))),a)
                     for w,a in vertex(poly,r,q,-n).items())
    return combine(terms)


def inner(left,right):
    total = Q.zero
    for a,x in left.items():
        for b,y in right.items():
            if len(a)!=len(b):
                continue
            for perm in permutations(b):
                if any(n!=m for (n,_),(m,_) in zip(a,perm,strict=True)):
                    continue
                value = old["m"].conj(x)*y
                for (n,i),(_,j) in zip(a,perm,strict=True):
                    value *= n*metric[i][j]
                total += value
    return Q.to_sympy(total)


def pairing(left,right):
    return sp.Matrix(len(left),len(right),lambda i,j:inner(left[i],right[j]))


def gram(packets):
    return pairing(packets,packets)


def rank(mat):
    return mat.to_DM().convert_to(Q).rank()


def encode(mat):
    return [[str(v) for v in row] for row in mat.tolist()]


unit = sp.eye(8)
weak = unit[:,3]-unit[:,4]
hypercharge = sp.Matrix([-2,-2,-2,3,3,0,0,0])/6
t, delta_t = sp.Rational(1,8),sp.Rational(1,20)
families = {}
for alpha in (sp.Rational(1,4),sp.Rational(3,4)):
    drift = sp.Matrix([j+alpha-4 for j in range(8)])
    inputs = [old["series_response"](f,[Q.from_sympy(z) for z in old["cl"](V.T*drift/8)],4)
              for f in polys]
    states = [old["vertex_zero"](f) for f in inputs]
    assert [vertex(old["input_word_map"](f),q,sp.zeros(8,1)) for f in inputs]==states
    assert all(not vertex(v,d0,q) for v in states)
    raw = gram(states)
    minus_images = [vertex(v,d1,q) for v in states]
    g_minus = gram(minus_images)
    mids = [middle(v) for v in states]
    # Independent check using the actual family lowering commutator:
    # the first term V_-q V_d0 vanishes on these states.
    lowers = [vertex(v,-q,q) for v in states]
    comm = [vertex(u,d0,sp.zeros(8,1)) for u in lowers]
    assert all(a==b for a,b in zip(mids,comm,strict=True)) or all(
        a=={word:-c for word,c in b.items()} for a,b in zip(mids,comm,strict=True))
    g_middle = gram(mids)/2   # unit norm of q(-1)|r> is sqrt(2)
    # A neutral family current inserted at the cut makes the two spin
    # paths meet. The two orderings have the same cocycle phase since
    # q dot d1 = -2, so the relative sign of the commutator is fixed.
    assert q.dot(d1)==-2
    q_after = [vertex(w,q,r) for w in minus_images]
    q_before = [vertex(vertex(v,q,q),d1,2*q) for v in states]
    differences = [combine([*a.items(),*((w,-v) for w,v in b.items())])
                   for a,b in zip(q_after,q_before,strict=True)]
    assert differences==mids or differences==[
        {w:-a for w,a in v.items()} for v in mids]
    sign = 1 if differences==mids else -1
    after_pair = sign*pairing(mids,q_after)
    before_pair = sign*pairing(mids,q_before)
    assert after_pair-before_pair==2*g_middle
    assert sp.simplify(g_minus[0,0]/raw[0,0])==sp.Rational(19017,4277)
    assert rank(g_minus)==5
    assert all(grade(word)==2 for result in mids for word in result)
    final = [q+d for d in (d0,r,d1)]
    assert [int(sp.exp(2*sp.pi*sp.I*p[0])) for p in final]==[-1,1,-1]
    assert [p.dot(p)/2 for p in final]==[5,2,1]
    energy = 32+drift.dot(r)
    assert energy==30+2*alpha
    assert all(drift.dot(d)==drift.dot(r) for d in (d0,d1))
    folded_drift = M@drift
    color = folded_drift[:3]
    weak_coefficient = folded_drift.dot(weak)/2
    gauge = sp.Matrix([x-sum(color)/3 for x in color]+
                      [weak_coefficient,-weak_coefficient,0,0,0])
    assert gauge.dot(M@q)==0
    for d in (d0,r,d1):
        charge = M@d
        assert charge[0]==charge[1]==charge[2]==0
        assert charge.dot(weak)==-1 and charge.dot(hypercharge)==sp.Rational(1,2)
        assert gauge.dot(charge)==2
    # At pi/2 the normalized spin-one operator has amplitudes 1/sqrt(2)
    # and 1/2 into the two nonzero, mutually orthogonal charge spaces.
    full_gram = g_middle/2+g_minus/4
    ramond_gram = g_minus/4

    def returns(matrix):
        c = sp.simplify(matrix[0,0]/raw[0,0])
        values = [c*sp.exp(-energy*t), c*sp.exp(-energy*2*t),
                  c*sp.exp(-energy*t-32*delta_t), c*sp.exp(-energy*(t+delta_t))]
        assert sp.simplify(c*values[1]-values[0]**2)==0
        assert values[2]<values[3]
        return {"norm_ratio":str(c),"t_2t_shift_R_shift_H":list(map(float,values)),
                "shift_R_minus_H":float(values[2]-values[3]),
                "shift_R_over_H_exact":str(sp.simplify(values[2]/values[3])),
                "two_time_defect_exact":0}
    # Non-proportionality matters for determining a missing spin contribution
    # from the full five-packet response.
    ratio = sp.simplify(g_middle[0,0]/g_minus[0,0])
    assert g_middle != ratio*g_minus
    assert rank(g_middle)>0
    families[str(alpha)] = {
        "relative_energy":str(energy),
        "initial_reference_energy":32,
        "retained_weak_transport_rate":2,
        "middle_NS_rank":rank(g_middle),
        "middle_Gram":encode(g_middle),
        "outer_R_rank":rank(g_minus),
        "outer_Gram":encode(g_minus),
        "middle_is_proportional_to_outer":False,
        "neutral_family_current_between_the_two_spin_paths": {
            "after_pairing_with_unnormalized_middle":encode(after_pair),
            "before_pairing_with_unnormalized_middle":encode(before_pair),
            "after_rank":rank(after_pair),
            "before_rank":rank(before_pair),
            "difference_equals_twice_normalized_middle_Gram":True,
            "phase_convention":"relative sign fixed by family commutator",
        },
        "pi_over_2_both_old_spin_sectors":returns(full_gram),
        "pi_over_2_old_Ramond_filter_only":returns(ramond_gram),
    }

gamma, mass, inner_form = (native[k] for k in ("gamma","majorana","inner"))
plus, minus = (np.flatnonzero(sgn*np.diag(gamma)>0) for sgn in (1,-1))
marks = [(-1j*gamma@inner_form(np.array(h)))[np.ix_(minus,plus)]
         for h in ((0,1),(0,1j),(1,0),(1j,0))]
assert all(np.linalg.matrix_rank(x)==48 for x in marks)
anchor = marks[0]
quadratic = anchor.conj().T@anchor
other = marks[1]
cross = anchor.conj().T@other
assert np.linalg.norm(cross+cross.conj().T)<1e-12
assert np.linalg.norm(other.conj().T@other-quadratic)<1e-12
ordered_difference = np.linalg.norm(cross-cross.conj().T)
assert ordered_difference>1
# Reconstruct all coefficients in one common frame using the complete ordered Gram.
# Only unitary orientation relative to independently specified end vectors remains.
e,v = np.linalg.eigh(quadratic)
inv_sqrt = (v*e**(-.5))@v.conj().T
recovered = [inv_sqrt@anchor.conj().T@b for b in marks]
ordered_error = max(np.linalg.norm(a.conj().T@b-u.conj().T@v)
                    for a,u in zip(marks,recovered,strict=True)
                    for b,v in zip(marks,recovered,strict=True))
assert ordered_error<1e-11
majorana = (-1j*gamma@mass)[np.ix_(minus,plus)]
polar_recovered = inv_sqrt@anchor.conj().T@majorana
majorana_error = max(np.linalg.norm(a.conj().T@majorana-u.conj().T@polar_recovered)
                    for a,u in zip(marks+[majorana],recovered+[polar_recovered],strict=True))
assert majorana_error<1e-11 and np.linalg.matrix_rank(majorana)==3

# A concrete residual relative orientation survives the given grading,
# real structure, SM representation and Majorana. It is a passive change
# of native frame. Keeping independent external ends fixed distinguishes it.
rep = native["native"]
quark = np.array([f in ("u","d") for f,*_ in rep.P]*2)
frame = np.diag(np.exp(-1j*np.pi/4*np.diag(gamma)*quark))
swap = native["swap"]
assert np.linalg.norm(frame-swap@frame.conj()@swap)<1e-13
assert np.linalg.norm(frame@mass@frame.conj().T-mass)<1e-13
zeros2,zeros3 = np.zeros((2,2)),np.zeros((3,3))
quaternions = [rep.quaternion(1,0),rep.quaternion(1j,0),
               rep.quaternion(0,1),rep.quaternion(0,1j)]
color_matrices = [k*np.outer(np.eye(3)[:,i],np.eye(3)[:,j])
                  for k in (1,1j) for i in range(3) for j in range(3)]
basis = ([rep.rep(z,zeros2,zeros3) for z in (1,1j)]+
         [rep.rep(0,q,zeros3) for q in quaternions]+
         [rep.rep(0,zeros2,c) for c in color_matrices])
assert all(np.linalg.norm(frame@a-a@frame)<1e-13 for a in basis)
rotated = [frame[np.ix_(minus,minus)]@a@frame[np.ix_(plus,plus)].conj().T for a in marks]
all_marks,all_rotated = marks+[majorana],rotated+[majorana]
frame_gram_error = max(np.linalg.norm(a.conj().T@b-u.conj().T@v)
                       for a,u in zip(all_marks,all_rotated,strict=True)
                       for b,v in zip(all_marks,all_rotated,strict=True))
i,j = np.unravel_index(np.argmax(np.abs(anchor*quark[minus,None])),anchor.shape)
assert quark[minus[i]] and quark[plus[j]] and abs(anchor[i,j])>1
assert abs(rotated[0][i,j]-1j*anchor[i,j])<1e-13 and frame_gram_error<1e-11
print(json.dumps({
    "old_exact_family_middle_and_time":families,
    "same_native_96_ordered_quadratic": {
        "pure_Yukawa_lower_complex_plane_gram_equals_radius_square_times_common_Gram": True,
        "lower_components_ordered_Gram_difference_HS_norm":float(ordered_difference),
        "native_anchoring_Yukawa_block_rank":int(np.linalg.matrix_rank(anchor)),
        "four_marks_ordered_reconstruction_error":float(ordered_error),
        "including_existing_Majorana_ordered_Gram_error":float(majorana_error),
        "native_Majorana_block_rank":int(np.linalg.matrix_rank(majorana)),
        "gauge_and_reality_compatible_quark_rephasing_Gram_error":float(frame_gram_error),
        "quark_amplitude_ratio_with_fixed_basis_ends":"+i",
        "common_output_unitary_undetermined_by_Grams":True,
    },
    "physical_full_Padm_values":None,
},indent=2))
