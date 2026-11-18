"""Exact determinant reduction of two native DET16 letters with a gauge twist."""
import contextlib
import io
import json
import runpy
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

repo = Path(sys.argv[1]).resolve()
own = Path(__file__).resolve().with_name("check_original_det16_selector_discriminants_20261109.py")
argv = sys.argv
try:
    sys.argv = [str(own), str(repo)]
    with contextlib.redirect_stdout(io.StringIO()):
        original = runpy.run_path(str(own))
finally:
    sys.argv = argv

native, n, dimension = (original[k] for k in ("native", "n", "dimension"))
projector, masks = (original[k] for k in ("projector", "positive"))
e, nu = (original[k] for k in ("e_index", "nu_index"))
s = np.diag(original["carrier_sign"]).astype(complex)
one = np.eye(n, dtype=complex)
fock_one = sparse.eye(dimension, dtype=complex, format="csr")
states = np.arange(dimension)
occupation = np.array([int(x).bit_count() for x in states])
contrast = sparse.diags(occupation/8-1, dtype=complex, format="csr")
mode = {i: native.mode_operator(n,i,False) for i in range(n)}
fock_det_sign = 2*projector-fock_one


def adjugate(x):
    # Use cofactors because r=-I is an actual center twist. I+r is singular.
    size = len(x)
    return np.array([[
        (-1)**(i+j)*np.linalg.det(np.delete(np.delete(x,j,axis=0),i,axis=1))
        for j in range(size)] for i in range(size)])


def coefficients(time):
    a = np.exp(-2*time)
    return a, -np.expm1(-2*time)


def empty(r, marker, beta):
    a,b = coefficients(beta)
    if marker is None:
        return 2-empty(r,one,beta)
    return a*np.linalg.det(one+r@marker)+b


def charged(r, marker, beta, time):
    a,b = coefficients(beta-time)
    c,d = coefficients(time)
    inverse = r.conj().T
    if marker is None:
        return a*inverse-charged(r,one,beta,time)
    return (a*c*(marker@adjugate(one+r@marker))
            + a*d/2*inverse + b*c/2*marker)


def neutral(r, marker, beta, time):
    a,b = coefficients(beta-time)
    c,d = coefficients(time)
    if marker is None:
        return 2*a-neutral(r,one,beta,time)
    polynomial = np.array([1+0j])
    for root in np.linalg.eigvals(r@marker):
        polynomial = np.convolve(polynomial,[1,root])
    return a*c*np.dot(polynomial,(np.arange(n+1)/8-1)**2)+a*d+b*c


def exterior_monomial(r):
    """Independent direct action on all 65536 basis states for these twists."""
    columns = np.arange(n)
    positions = np.argmax(np.abs(r), axis=0)
    phases = r[positions, columns]
    assert len(set(positions)) == n
    assert np.max(np.abs(r-np.eye(n)[:,positions]*phases)) < 1e-13
    rows = np.zeros(dimension,dtype=np.int64)
    values = np.ones(dimension,dtype=complex)
    occupied = [(states>>i&1) == 1 for i in range(n)]
    for i in range(n):
        rows[occupied[i]] |= 1<<positions[i]
        values[occupied[i]] *= phases[i]
        for j in range(i):
            if positions[j] > positions[i]:
                values[occupied[i]&occupied[j]] *= -1
    return sparse.csr_matrix((values,(rows,states)),shape=(dimension,dimension))


def wedge_even(carrier):
    coordinates = [[j for j in range(5) if (int(bit)>>(4-j))&1] for bit in masks]
    return np.array([[np.linalg.det(carrier[np.ix_(rows,cols)])
                      if len(rows)==len(cols) else 0
                      for cols in coordinates] for rows in coordinates])


def heat(time):
    a,b = coefficients(time)
    return a*fock_one+b*projector


def trace(x):
    return x.diagonal().sum()


# An actual weak SU(2) representative, acting on the native even wedge.
# It does not commute with the determinant-minus direction chosen in s.
carrier = np.eye(5,dtype=complex)
carrier[3:,3:] = [[0,1],[-1,0]]
weak = wedge_even(carrier)
assert abs(np.linalg.det(weak)-1) < 1e-12
assert np.linalg.norm(weak@s-s@weak) > 1
markers = {"no_mark":one,"gauge_sign":s,"fermion_parity":-one,"rank_one_DET":None}
twists = {str(q):q*one for q in (1,1j,-1,-1j)}
twists["noncommuting_weak_rotation"] = weak
beta,time = 1.,.8
max_errors = {"empty":0.,"charged":0.,"neutral":0.}
fock_markers = {key:(fock_det_sign if val is None else exterior_monomial(val))
                for key,val in markers.items()}
raw = {}
for key,r in twists.items():
    gr = exterior_monomial(r)
    assert np.max(np.abs((gr@projector-projector).data),initial=0.) < 1e-13
    assert np.max(np.abs((gr@projector-projector@gr).data),initial=0.) < 1e-13
    raw[key] = {}
    for label,m in markers.items():
        mid = heat(time/2)@fock_markers[label]@heat(time/2)
        left = gr@heat(beta-time)
        direct_empty = trace(left@mid)
        direct_neutral = trace(left@contrast@mid@contrast)
        model_empty = empty(r,m,beta)
        model_neutral = neutral(r,m,beta,time)
        max_errors["empty"] = max(max_errors["empty"],abs(direct_empty-model_empty))
        max_errors["neutral"] = max(max_errors["neutral"],abs(direct_neutral-model_neutral))
        model_charged = charged(r,m,beta,time)
        positions = [(e,e),(nu,nu),(e,nu),(nu,e)]
        if key=="noncommuting_weak_rotation":
            positions += list(zip(*np.nonzero(np.abs(weak)>1e-12)))
        for i,j in set(positions):
            direct_charged = trace(left@mode[i]@mid@mode[j].getH())
            max_errors["charged"] = max(max_errors["charged"],
                                        abs(direct_charged-model_charged[i,j]))
        raw[key][label] = np.array([model_empty,model_charged[e,e],
                                    model_charged[nu,nu],model_neutral])
assert max(max_errors.values()) < 2e-9

center = {key:sum(raw[t][key] for t in list(twists)[:4])/4 for key in markers}
baseline = center["no_mark"]
assert np.min(baseline.real) > 0
ratios = {key:center[key]/baseline for key in markers}
assert all(np.max(np.abs(x.imag)) < 1e-12 for x in ratios.values())
assert np.max(np.abs(ratios["fermion_parity"]-[1,-1,-1,1])) < 1e-12

# The original DET4 content shadow is exp(2pi i Y). Its weak determinant
# is +1. Compare the two possible determinant-sign extensions within the
# original native S(U(3)xU(2)) carrier, with full group averages in DET16.
weights = original["cartan_weights"]
z = np.diag(np.exp(1j*weights@(np.pi*np.array([-2/3,-2/3,-2/3,1,1]))))
central = np.diag(np.exp(1j*weights@(np.pi*np.array([-1,-1,1,.5,.5]))))
e_left = int(np.flatnonzero(masks==29)[0])
assert np.max(np.abs(central@central-np.linalg.matrix_power(z,3))) < 1e-12
assert abs(z[e,e]-1) < 1e-12
assert abs(s[e,e]+1) < 1e-12


def group_read(mark):
    # Distinct z^k and z^k mark, k=0,...,5; det(weak) distinguishes the cosets.
    subgroup = [np.linalg.matrix_power(z,k)@x for k in range(6) for x in (one,mark)]
    read_modes = [e,nu,e_left]
    def block(x):
        diagonal = sum(np.diag(charged(r,x,beta,time)) for r in subgroup)/12
        vacuum = sum(empty(r,x,beta) for r in subgroup)/12
        return np.array([vacuum,*diagonal[read_modes]])
    plain = block(one)
    marked = block(mark)/plain
    double = block(mark@mark)/plain
    chars = np.array([1,*np.diag(mark)[read_modes]])
    assert np.max(np.abs(marked-chars)) < 1e-12
    assert np.max(np.abs(double-chars**2)) < 1e-12
    return {"native_projected_empty":float(plain[0].real),
            "relative_mark":[[float(x.real),float(x.imag)] for x in marked],
            "relative_double_mark":[[float(x.real),float(x.imag)] for x in double]}


group_results = {"split_C6_times_C2":group_read(s),
                 "central_C12":group_read(central)}

def total_residue(diagonal, order):
    one_body_charge = np.rint(np.angle(diagonal)*order/(2*np.pi)).astype(int)
    result = np.zeros(dimension,dtype=int)
    for i,q in enumerate(one_body_charge):
        result += q*(states>>i&1)
    return result%order


common = total_residue(np.diag(z),6)==0
split_mask = common & (total_residue(np.diag(s),2)==0)
central_mask = common & (total_residue(np.diag(central),4)==0)
rank_split = np.count_nonzero(split_mask)
rank_central = np.count_nonzero(central_mask)
assert rank_split == rank_central
assert np.any(split_mask != central_mask)
one_particle_states = 1 << np.arange(n)
one_particle_split = split_mask[one_particle_states]
one_particle_central = central_mask[one_particle_states]
assert np.array_equal(one_particle_split,one_particle_central)
assert np.array_equal(split_mask,split_mask[::-1])
assert np.array_equal(central_mask,central_mask[::-1])
degree_second_difference = int(((occupation[central_mask]-8)**2).sum()
                               - ((occupation[split_mask]-8)**2).sum())
assert degree_second_difference == 624
first_pair = next((i,j) for i in range(n) for j in range(i+1,n)
                  if np.count_nonzero(central_mask & (states>>i&1==1) & (states>>j&1==1))
                  != np.count_nonzero(split_mask & (states>>i&1==1) & (states>>j&1==1)))
i,j = first_pair
occupied_pair = (states>>i&1==1) & (states>>j&1==1)
pair_difference = int(np.count_nonzero(central_mask & occupied_pair)
                      - np.count_nonzero(split_mask & occupied_pair))
a,b = coefficients(beta)
assert abs(group_results["central_C12"]["native_projected_empty"]-(a*rank_central+b)) < 1e-10
assert abs(group_results["split_C6_times_C2"]["native_projected_empty"]-(a*rank_split+b)) < 1e-10
joint = int(np.count_nonzero(split_mask & central_mask))
first_difference = int(np.flatnonzero(split_mask!=central_mask)[0])
group_results["exact_census"] = {
    "rank_of_each":int(rank_central),"intersection_rank":joint,
    "native_one_particle_invariant_wedge_masks":masks[one_particle_split].tolist(),
    "identical_full_ordered_unmarked_two_point_matrices_for_all_beta_and_t":True,
    "central_minus_split_sum_of_N_minus_8_squared":degree_second_difference,
    "central_minus_split_normalized_neutral_total_number_word":
        float(a*degree_second_difference/64/(a*rank_split+b)),
    "first_unmarked_density_pair_difference":{
        "carrier_wedge_masks":[int(masks[i]),int(masks[j])],
        "central_minus_split_count":pair_difference},
    "split_occupation_degree_histogram":np.bincount(occupation[split_mask],minlength=n+1).tolist(),
    "central_occupation_degree_histogram":np.bincount(occupation[central_mask],minlength=n+1).tolist(),
    "one_different_occupation":{
        "mask":first_difference,
        "occupied_internal_wedge_masks":[int(masks[i]) for i in range(n) if (first_difference>>i)&1],
        "in_split_selection":bool(split_mask[first_difference]),
        "in_central_selection":bool(central_mask[first_difference]),
    },
}

print(json.dumps({
    "scope": "native original DET16 mirror factor; natural supplied carrier twists only, not the missing v1013 gauge embedding or actual QTheta",
    "fock_dimension": dimension,
    "time_parameters":{"beta":beta,"between_field_ends":time},
    "maximum_absolute_errors_vs_direct_fock":max_errors,
    "native_spin10_center_projected_raw_unmarked":baseline.real.tolist(),
    "readout_order":["empty","e_conjugate","nu_conjugate","total_number_contrast"],
    "native_spin10_center_projected_marked_over_unmarked":
        {key:val.real.tolist() for key,val in ratios.items()},
    "tested_noncommuting_carrier_pair_commutator_norm":float(np.linalg.norm(weak@s-s@weak,2)),
    "native_carrier_subgroup_comparisons":{
        "readout_order":["empty","e_conjugate","nu_conjugate","left_e"],
        "results":group_results},
    "actual_v1013_Gauss_embedding_or_QTheta_supplied":False,
},indent=2))
