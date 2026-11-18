"""Compare supplied TFPT hypotheses on the unchanged Weyl/DET16 carriers."""
import itertools
import json
import runpy
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

repo = Path(sys.argv[1]).resolve()
sys.path[:0] = [str(repo), str(repo / "verification")]
import v1002_spin10_mirror_projector as native

n = native.SPIN10_MODES
gens, pairs = native.positive_weyl_generators()
cartans = [gens[pairs.index((2*j,2*j+1))] for j in range(5)]
twice_weights = np.array([2*np.diag(x).real for x in cartans], dtype=int).T
bits = ((twice_weights+1)//2).astype(int)
masks = bits @ (2**np.arange(4,-1,-1))
assert np.array_equal(np.sort(masks), [k for k in range(32) if k.bit_count()%2==0])
# Units of pi/6 on Lambda^even C^5. Sum of carrier charges is zero, so
# subtracting the common half weight in the Spin lift makes no difference.
charge = bits @ np.array([-2,-2,-2,3,3])
det_g = (twice_weights @ np.array([-1,-1,1,1,0]))//2
assert np.array_equal((twice_weights @ np.array([-1,-1,1,1,0]))%2, np.zeros(n))
k = np.arange(1<<n)
occ = ((k[:,None] >> np.arange(n)) & 1).astype(int)
total = occ @ charge
common = (total%6==0)
split = common & ((occ @ det_g)%2==0)
central = common & (total%12==0)
assert (sum(split),sum(central),sum(split & central)) == (5464,5464,3072)
assert set(charge) == {-4,-3,0,1,2,6}
assert np.array_equal((4*charge)%12, 4*(bits[:,:3].sum(axis=1)%3))
assert np.array_equal((6*charge)%12, 6*(bits[:,3:].sum(axis=1)%2))
charges = [{"mask":int(m),"Y":str(Fraction(int(q),6)),
            "M_exponent_mod12":int(q%12),"M9_exponent_mod12":int(9*q%12),
            "split_sign":int((-1)**int(g%2))}
           for m,q,g in zip(masks,charge,det_g)]

def mode_class(q_mod6, second, residue):
    return np.flatnonzero((charge%6 == q_mod6) & (second == residue)).tolist()


def parity_slice(baseline, changes, selector):
    # Changes replace disjoint bundles of original occupations. Every
    # occupation is affine in exactly one Boolean variable or constant.
    used = [i for bundle in changes for i in bundle]
    assert len(used)==len(set(used))
    start = sum(1<<i for i in baseline)
    toggles = [sum(1<<i for i in bundle) for bundle in changes]
    checked = 0
    for active in itertools.product((0,1), repeat=len(changes)):
        state = start
        for flag, toggle in zip(active,toggles):
            if flag:
                state ^= toggle
        assert common[state]
        assert selector[state] == (sum(active)%2==0)
        checked += 1
    return {"independent_parity_variables":len(changes),"checked_vertices":checked,
            "start_wedge_masks":masks[baseline].tolist(),
            "toggle_wedge_mask_groups":[masks[bundle].tolist() for bundle in changes]}


gs = det_g%2
cs = (3*charge)%4
e = int(np.flatnonzero(masks==3)[0])
split_start, split_switches = [], []
for q in (1,2,3):
    even = mode_class(q,gs,0)
    odd = mode_class(q,gs,1)
    assert len(even)==len(odd)
    split_start.extend(even)
    split_switches.extend([list(pair) for pair in zip(even,odd)])
split_switches.append([e])
central_start = mode_class(2,cs,0) + mode_class(3,cs,3) + [e]
central_switches = [list(pair) for pair in zip(mode_class(2,cs,0),mode_class(2,cs,2))]
odd_three, quarks = mode_class(3,cs,3), mode_class(1,cs,3)
assert len(odd_three)==2 and len(quarks)==6
central_switches.extend([[i]+quarks[3*j:3*j+3] for j,i in enumerate(odd_three)])
central_switches.append([e])
split_slice = parity_slice(split_start,split_switches,split)
central_slice = parity_slice(central_start,central_switches,central)
assert (split_slice["independent_parity_variables"],central_slice["independent_parity_variables"])==(8,6)

color_count = occ[:,charge%3==1].sum(axis=1)-occ[:,charge%3==2].sum(axis=1)
weak_parity = occ[:,charge%2==1]
split_parity = occ[:,gs==1]
assert weak_parity.shape[1]==split_parity.shape[1]==8
assert (color_count.min(),color_count.max())==(-6,6)
assert (total.min(),total.max())==(-18,18)
color_penalty = color_count**2*(color_count**2-9)**2*(36-color_count**2)
split_penalty = (color_penalty + (1-np.prod(1-2*weak_parity,axis=1))
                              + (1-np.prod(1-2*split_parity,axis=1)))
central_penalty = total**2*(total**2-144)**2
assert split_penalty.min()==central_penalty.min()==0
assert np.array_equal(split_penalty==0,split)
assert np.array_equal(central_penalty==0,central)

clock_path = repo / "experiments/theory-contracts/joint-native-algebra-selection-20261016/marked_lift_checker.py"
old_clock = runpy.run_path(str(clock_path))
clean = old_clock["clean"]
audit = old_clock["audit"]
marker = old_clock["build_marker_data"]()["data"]
_, supports = old_clock["_code_basis"]()
sigma, _, _ = old_clock["_source_operators"]()
sigma5 = audit.quartic_compression(sigma, supports)
rotation = old_clock["_parse_matrix"](marker["basis_alignment"]["rotation_matrix"])
frame = clean(rotation * old_clock["_process_frame"]().T).T
tau = clean(frame.T * audit.permutation_matrix((1,0,2,3,4,5)) * frame)
i5 = sp.eye(5)
y5 = old_clock["_parse_matrix"](marker["native_polar_transport"]["native_charge_matrix"])
p2 = clean((6*y5 + 2*i5)/5)
p3 = i5-p2
f = clean((i5+sigma5+sigma5**2)/3)
r5 = clean(p3*(i5-2*f)+sp.I*p2*tau)
c5 = clean(-p3+sp.I*p2)
assert clean(r5*r5)==clean(c5*c5)==clean(i5-2*p2)
assert sp.simplify(r5.det()) == -1 and sp.simplify(c5.det()) == 1
weak_r_det = sp.simplify(clean(p3+r5*p2).det())
weak_c_det = sp.simplify(clean(p3+c5*p2).det())
assert (weak_r_det,weak_c_det)==(1,-1)

def exterior_multiplicities(carrier):
    raw = clean(carrier).eigenvals()
    spectrum = [phase for phase in (1,sp.I,-1,-sp.I)
                for _ in range(raw.get(phase,0))]
    assert len(spectrum)==5
    products = [sp.prod(spectrum[i] for i in s) for d in (0,2,4)
                for s in itertools.combinations(range(5),d)]
    return [products.count(phase) for phase in (1,sp.I,-1,-sp.I)]

r_spec, c_spec = exterior_multiplicities(r5), exterior_multiplicities(c5)
assert (r_spec,c_spec)==([4,4,4,4],[4,0,4,8])
print(json.dumps({
    "scope":"actual v1002 Weyl weights and two existing DET16 diagonal comparison selectors",
    "one_particle_rows":charges,
    "original_selection_ranks":[int(sum(split)),int(sum(central)),int(sum(split & central))],
    "split_selector_parity_slice":split_slice,
    "central_selector_parity_slice":central_slice,
    "minimal_exact_ground_parent_field_degrees":{"split":16,"central":12},
    "explicit_nonnegative_diagonal_penalty_checked_all_configurations":int(len(k)),
    "exact_original_quartic_logical_clock": {
        "R5_and_M9_same_carrier_square":True,
        "det_R5":int(r5.det()),
        "det_M9":int(c5.det()),
        "relative_e_conjugate_over_nu_conjugate_for_exterior_lifts":[int(weak_r_det),int(weak_c_det)],
        "relative_even_exterior_phase_multiplicities_R5_then_M9":[r_spec,c_spec]
    }
},indent=2))
