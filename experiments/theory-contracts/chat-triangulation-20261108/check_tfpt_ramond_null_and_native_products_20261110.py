"""Exact native Higgs products and ordered tests of the four old Ramond nulls."""
from contextlib import redirect_stdout
from itertools import product
from pathlib import Path
import io
import json
import sys
import runpy

import numpy as np

repo = Path(sys.argv[1]).resolve()
base = repo/"experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    lift = runpy.run_path(str(base/"check_tfpt_ramond_native_higgs_lifts_20261110.py"))
    history = runpy.run_path(str(base/"check_tfpt_ramond_marked_histories_20261110.py"))

sp,Q = lift["sp"],lift["Q"]
old = lift["old"]
zero_mode, inner = lift["zero_mode"],lift["inner"]
exact_word = history["exact_word"]
native = lift["native"]
n = native.NP
A = lift["A"][:n,:n]
B = A@A.conj().T-A.conj().T@A
u, ss, vh = np.linalg.svd(A)
products = []
for field, y in lift["mass_blocks"].items():
    for index,strength in enumerate(2*np.linalg.svd(y,compute_uv=False)):
        ix = np.flatnonzero(abs(ss-strength)<1e-10)
        c = (u[:,ix]*ss[ix])@vh[ix,:]
        h = (c@c.conj().T-c.conj().T@c)/strength**2
        residual = np.linalg.norm(h@c-c@h-2*c)
        assert residual<1e-12
        products.append({"species":field,"generation_singular_index":index,
                         "singular_square":float(strength**2),"number_of_channels":len(ix),
                         "local_CAR_su2_bracket_residual":float(residual)})
        assert np.linalg.norm(B@c-c@B-2*strength**2*c)<2e-12
assert len(products)==12


def mark(state,index,keys,phase):
    if index==4:
        return exact_word(state,phase)
    return {bits:a for bits,a in state.items() if all(bits>>k&1 for k in keys[index])}


def mix(words, coeffs):
    out = {}
    for word,c in zip(words,coeffs,strict=True):
        if not c:
            continue
        for bits,a in word.items():
            out[bits] = out.get(bits,Q.zero)+c*a
    return {k:v for k,v in out.items() if v}


def gmat(vectors):
    return sp.Matrix(len(vectors),len(vectors),lambda i,j:inner(vectors[i],vectors[j]))


def exact_rank(matrix):
    return matrix.to_DM().convert_to(Q).rank()


q = lift["q"]
q_coeff = [Q.from_sympy(x) for x in q]


def cartan(state, mode, ramond=True):
    if mode==0:
        offset = sp.Rational(11,2) if ramond else 5
        width = len(old["levels"])
        mask = (1<<width)-1
        out = {}
        for bits,a in state.items():
            occupations = [((bits>>(j*width))&mask).bit_count() for j in range(8)]
            charge = sum(x*(number-offset) for x,number in zip(q,occupations,strict=True))
            if charge:
                out[bits] = a*Q.from_sympy(charge)
        return out
    out = {}
    for bits,a in state.items():
        for j,v in enumerate(q_coeff):
            for new,sign in old["step"](bits,-mode,j):
                out[new] = out.get(new,Q.zero)+sign*v*a
    return {k:v for k,v in out.items() if v}


def difference_root_mode(state,charge,mode):
    i,j = list(charge).index(1),list(charge).index(-1)
    output = {}
    for level in old["levels"]:
        if level-mode not in old["levels"]:
            continue
        pair = ((old["bit"](i,level-mode),True),(old["bit"](j,level),False))
        for k,v in exact_word(state,pair).items():
            output[k] = output.get(k,Q.zero)+v
    return {k:v for k,v in output.items() if v}


def local_normal_product(state,charge,ramond=True):
    # Weight-two zero mode of :J_(q.H) J_r:. With q.r=0 the two
    # currents have no mutual singular OPE. Positive modes act first.
    terms = []
    for mode in range(-4,5):
        if mode<0:
            terms.append(cartan(difference_root_mode(state,charge,-mode),mode,ramond))
        else:
            terms.append(difference_root_mode(cartan(state,mode,ramond),charge,-mode))
    return mix(terms,[Q.one]*len(terms))


def four_field_zero_mode(state,charge):
    channels = [(j,int(x)) for j,x in enumerate(charge) if x]
    levels = old["levels"]
    output = {}
    for bits,a in state.items():
        available = [{n for n in levels if bool(bits>>old["bit"](j,n)&1)!=(sign>0)}
                     for j,sign in channels]
        for start in product(*(sorted(x) for x in available[:3])):
            last = -channels[3][1]*sum(n*sign for n,(_,sign) in zip(start,channels[:3],strict=True))
            if last not in available[3]:
                continue
            final, coefficient = bits, a
            for n,(j,_) in reversed(list(zip((*start,last),channels,strict=True))):
                index = old["bit"](j,n)
                if (final&((1<<index)-1)).bit_count()%2:
                    coefficient = -coefficient
                final ^= 1<<index
            output[final] = output.get(final,Q.zero)+coefficient
    return {k:v for k,v in output.items() if v}


reports = {}
for alpha, states in history["quarter_states"].items():
    known = history["joint"][alpha]
    keys = known["number_word_modes"]
    (i,j),(k,l) = known["phase_sensitive_word_modes"]
    phase = ((i,True),(j,True),(l,False),(k,False))
    G = sp.Matrix(known["all_five_same_words_exact_gram"])
    marks = [mark(v,i,keys,phase) for v in states for i in range(5)]
    assert gmat(marks)==G
    # Normalize each exact null vector by its last nonzero coordinate.
    rows = G.to_DM().convert_to(Q).nullspace().to_Matrix().tolist()
    nulls = [[Q.from_sympy(sp.expand(x)) for x in row] for row in rows]
    nulls = [[x/next(y for y in reversed(row) if y) for x in row] for row in nulls]
    assert len(nulls)==4
    assert all(not mix(marks,row) for row in nulls)

    charges = [sp.Matrix(x["lower_root_old_charge"]) for x in lift["doublets"] if x["Xi"]==1]
    charges.append(lift["E"][:,3]-lift["E"][:,7])
    activated = []
    joint = sp.zeros(4)
    delta = sp.Matrix([(j+sp.Rational(alpha))/8-sp.Rational(1,2) for j in range(8)])
    for charge in charges:
        after = [zero_mode(word,charge) for word in marks]
        assert all(not mix(after,row) for row in nulls)
        early = [mark(zero_mode(v,charge),i,keys,phase) for v in states for i in range(5)]
        revived = [mix(early,row) for row in nulls]
        revived_gram = gmat(revived)
        if all(x==0 for x in revived_gram):
            rank = 0
        else:
            rank = exact_rank(revived_gram)
        joint += revived_gram
        activated.append({"charge_old":list(map(str,charge)),
                          "late_same_operation_all_null_exact":True,
                          "early_null_subspace_rank":rank,
                          "null_norms_squared_exact":[str(x) for x in revived_gram.diagonal()],
                          "old_H_intermediate_energy":str(32+8*delta.dot(charge))})
    extra = {}
    root = next(charge for charge in charges[:-1] if charge.dot(q)==0)
    for key,operation in (("original_neutral_family_Cartan_minus_one",lambda state:cartan(state,-1)),
                          ("original_local_normal_product",lambda state:local_normal_product(state,root))):
        assert all(not mix([operation(word) for word in marks],row) for row in nulls)
        early = [mark(operation(v),i,keys,phase) for v in states for i in range(5)]
        revived_gram = gmat([mix(early,row) for row in nulls])
        extra[key] = {"early_rank":exact_rank(revived_gram),
                      "norms_squared_exact":[str(x) for x in revived_gram.diagonal()]}
    composite = [local_normal_product(v,root) for v in states]
    plain = [zero_mode(v,root) for v in states]
    product_gram = gmat(composite)
    raw = gmat(states)
    drift = [Q.from_sympy(x) for x in old["cl"](old["V"].T*delta)]
    inputs = [old["old"]["series_response"](f,drift,4) for f in old["polynomials"]]
    ns = [old["as_car"](old["vac_ns"],old["old"]["input_word_map"](f)) for f in inputs]
    ns_composite = [local_normal_product(v,root,False) for v in ns]
    for mode in (3,4):
        assert all(not difference_root_mode(v,root,mode) for v in states)
    assert all(not cartan(v,4) for v in states)
    assert all(not difference_root_mode(v,root,4) for v in ns)
    shift = 8*delta.dot(root)
    q_ratios = [sp.simplify(product_gram[i,i]/raw[i,i]) for i in range(5)]
    ns_ratios = [sp.simplify(inner(y,y)/inner(x,x)) for x,y in zip(ns,ns_composite,strict=True)]
    extra["original_local_normal_product"].update({
        "charge":list(map(str,root)),"conformal_weight":2,
        "same_old_time_H_intermediate_energy":str(32+shift),
        "closed_q_Gram_rank":exact_rank(product_gram),
        "rank_root_and_composite_together":exact_rank(gmat(plain+composite)),
        "exact_closed_q_norm_ratios":list(map(str,q_ratios)),
        "exact_ns_norm_ratios":list(map(str,ns_ratios)),
        "closed_q_minus_ns_at_old_time_1_over_8":[float((x-y)*sp.exp(-shift/8)) for x,y in zip(q_ratios,ns_ratios,strict=True)]})
    charges4 = [sp.Matrix([0,1,1,-1,0,0,1,0]),sp.Matrix([1,0,0,0,1,-1,0,1])]
    same_frequency = []
    for charge in charges4:
        drift4 = 8*delta.dot(charge)
        assert drift4 == -2+2*sp.Rational(alpha)
        g4 = gmat([four_field_zero_mode(v,charge) for v in states])
        same_frequency.append({"charge":list(map(str,charge)),
                               "unshifted_ground_energy":str((q+charge).dot(q+charge)/2),
                               "old_H_shift":str(drift4),"rank":exact_rank(g4),
                               "exact_P0_norm_over_q_norm":str(sp.simplify(g4[0,0]/raw[0,0]))})
    assert [x["rank"] for x in same_frequency]==[0,5]
    assert same_frequency[1]["exact_P0_norm_over_q_norm"]=="19017/4277"
    reports[alpha] = {
        "exact_old_rank":exact_rank(G),
        "null_vectors_sparse_packet_word_order": [
            {str(i):str(Q.to_sympy(v)) for i,v in enumerate(row) if v} for row in nulls],
        "old_current_mark_before_vs_after_same_halfword":activated,
        "combined_early_orthogonal_charge_control_rank_on_four_nulls":exact_rank(joint),
        "other_actual_ordered_controls":extra,
        "independent_same_frequency_four_CAR_product_check":same_frequency,
    }

# Existing original clock comparison under the declared equal-period-fraction
# convention. This is a multiplier of a marked cross, not a value for it.
clock = history["fields"]
E = history["old_energy_in_pin_x"]
step = .05
M2 = [{"radius":r["radius"],"d":step,
       "lower_bound_abs_relative_shift_difference":
           float(np.exp(-step*r["normal_zero_abs_K_max"])-np.exp(-step*E)),
       "upper_bound_abs_relative_shift_difference":
           float(np.exp(step*r["normal_zero_abs_K_max"])-np.exp(-step*E))}
      for r in clock]

print(json.dumps({
    "twelve_native_particle_CAR_product_components":products,
    "exact_old_ramond_null_tests":reports,
    "stationary_old_semigroup_two_time_defect_exact_zero":{
        "old_energy":32,"rank_of_old_word_span":21,
        "scope":"Common autonomous self-adjoint time with same cut Gram and old diagonal at t and 2t. Excludes leakage to orthogonal channels."},
    "marked_seam_relative_shift_multipliers_if_same_old_clocks":M2,
    "physical_full_Padm_Ramond_Higgs_Ramond":None,
},indent=2))
