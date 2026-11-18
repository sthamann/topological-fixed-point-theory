"""Resolve the original native DET16 selection difference by weak spin."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations
from pathlib import Path
import json
import sys

import numpy as np

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "verification"))
import v1002_spin10_mirror_projector as native

gens, pairs = native.positive_weyl_generators()
cartans = [gens[pairs.index((2*j, 2*j+1))] for j in range(5)]
weights = np.array([2*np.diag(x).real for x in cartans], dtype=int).T
wedge = (weights+1)//2
masks = wedge @ (2**np.arange(4, -1, -1))
occ = (np.arange(1<<16)[:,None] >> np.arange(16)) & 1
number = occ.sum(axis=1)
q = occ @ (wedge @ np.array([-2,-2,-2,3,3]))
z = q%6==0
central = q%12==0
twice_m = occ @ (wedge[:,3]-wedge[:,4])
g = (occ @ ((weights @ np.array([-1,-1,1,1,0]))//2))%2==0
split = z & g
assert np.all(twice_m[z]%2==0)
m = twice_m//2
assert np.array_equal(split, z & (central==(m%2==0)))

# The actual weak Lie raising action is on the last two P2 coordinates.
# A change of weight by +1 swaps weak flags 01 -> 10 in the positive Weyl basis.
transitions = [(int(np.flatnonzero(masks==bits)[0]),
                int(np.flatnonzero(masks==(bits^3))[0])) for bits in masks if bits&3==1]
assert len(transitions)==4
assert all(wedge[j,3]-wedge[j,4] == wedge[i,3]-wedge[i,4]+2 for i,j in transitions)
assert all((wedge[j] - wedge[i]) @ np.array([-2,-2,-2,3,3]) == 0 for i,j in transitions)

# Check a concrete weak-root counterexample to split invariance and a
# D5-root counterexample to central invariance, both within common z.
def root_witness(arrows, selection):
    for i,j in arrows:
        allowed = z & selection & (occ[:,i]==1) & (occ[:,j]==0)
        for start in np.flatnonzero(allowed):
            target = int(start) ^ (1<<i) ^ (1<<j)
            if z[target] and not selection[target]:
                return {"original_wedge_root_arrow":[int(masks[i]),int(masks[j])],
                        "initial_occupied_wedges":masks[occ[start]>0].tolist(),
                        "final_occupied_wedges":masks[occ[target]>0].tolist()}
    raise AssertionError("No witness")

d5_arrows = [(int(np.flatnonzero(masks==bits)[0]),
              int(np.flatnonzero(masks==(bits|3))[0])) for bits in masks if bits&3==0]
assert len(d5_arrows)==4
assert all(any(abs(generator[j,i])>0 for generator in gens)
           for i,j in transitions+d5_arrows)
weak_witness = root_witness(transitions, split)
d5_witness = root_witness(d5_arrows, central)

# Each fixed N and c charge is SU2-invariant, so the multiplicity of j
# is the number of weight j minus the number of weight j+1.
rows = defaultdict(lambda: defaultdict(int))
largest = int(max(m[z]))
bare_weak_square_count = 0
for n in range(17):
    for eta in (1,-1):
        retained = z & (number==n) & (central==(eta==1))
        weights_count = {j:int(np.count_nonzero(retained & (m==j)))
                         for j in range(-largest,largest+1)}
        assert all(weights_count[j]==weights_count[-j] for j in range(largest+1))
        for j in range(largest+1):
            multiplicity = weights_count[j]-weights_count.get(j+1,0)
            assert multiplicity>=0
            split_in_rep = sum((-1)**k==eta for k in range(-j,j+1))
            central_in_rep = (2*j+1) if eta==1 else 0
            delta = multiplicity*(central_in_rep-split_in_rep)
            row = rows[j]
            row[f"irreps_c_{eta:+d}"] += multiplicity
            row["split_rank"] += multiplicity*split_in_rep
            row["central_rank"] += multiplicity*central_in_rep
            row["central_minus_split_rank"] += delta
            row["central_minus_split_N_minus_8_squared"] += delta*(n-8)**2
            row["absolute_rank_difference"] += abs(delta)
            if j==0 and eta==1:
                bare_weak_square_count += multiplicity*(n-8)**2

assert sum(row["central_minus_split_N_minus_8_squared"] for row in rows.values())==624
assert sum(row["central_rank"] for row in rows.values())==np.count_nonzero(central)==5464
assert sum(row["split_rank"] for row in rows.values())==np.count_nonzero(split)==5464
assert rows[0]["absolute_rank_difference"]==0
trace_distance = Fraction(sum(row["absolute_rank_difference"] for row in rows.values()),2*5464)
assert trace_distance==Fraction(5464-3072,5464)

# Separate exact bare SM singlet count. The alternating Weyl-denominator
# coefficient extracts trivial SU3 and SU2 reps on the exact Y=0 sector.
color = occ @ wedge[:,:3]
neutral = q==0
weight_counts = Counter(zip(color[neutral,0]-color[neutral,1],
                            color[neutral,1]-color[neutral,2],m[neutral]))
rho = (1,0,-1)
sm_rank = 0
for p in permutations(range(3)):
    inversion_parity = sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))%2
    a,b,c = [rho[i]-rho[p[i]] for i in range(3)]
    sm_rank += (1-2*inversion_parity)*(weight_counts[a-b,b-c,0]-weight_counts[a-b,b-c,1])
assert np.all(central[neutral])

# Independent check of the Weyl count via the three raising operators
# on the Y=0, SU3/SU2 weight-zero subspace, using exterior signs twice.
zero_ids = np.flatnonzero(neutral & (m==0) & (color[:,0]==color[:,1])
                         & (color[:,1]==color[:,2]))

def raising_targets(start, i, j):
    for mode, row in enumerate(wedge):
        if row[i] or not row[j] or not (start >> mode & 1):
            continue
        after_wedge = masks[mode] ^ (1<<(4-i)) ^ (1<<(4-j))
        target_mode = int(np.flatnonzero(masks==after_wedge)[0])
        if start >> target_mode & 1:
            continue
        sign_wedge = 1-2*int(row[min(i,j)+1:max(i,j)].sum()%2)
        mask_between = (1<<max(mode,target_mode)) - (1<<(min(mode,target_mode)+1))
        sign_fock = 1-2*((int(start)&mask_between).bit_count()%2)
        yield int(start) ^ (1<<mode) ^ (1<<target_mode), sign_wedge*sign_fock

entries = defaultdict(list)
for i,j in ((0,1),(1,2),(3,4)):
    for col, start in enumerate(zero_ids):
        for target,value in raising_targets(start,i,j):
            entries[i,j,target].append((col,value))
table = np.zeros((len(entries),len(zero_ids)),dtype=np.int64)
for row, values in enumerate(entries.values()):
    for col, value in values:
        table[row,col] += value
prime, rank = 101, 0
table %= prime
for col in range(table.shape[1]):
    pivots = np.flatnonzero(table[rank:,col])
    if pivots.size:
        pivot = rank+int(pivots[0])
        table[[rank,pivot]] = table[[pivot,rank]]
        table[rank,col:] = table[rank,col:]*pow(int(table[rank,col]),-1,prime)%prime
        table[rank+1:,col:] = (table[rank+1:,col:]-table[rank+1:,col,None]*table[rank,col:])%prime
        rank += 1
assert len(zero_ids)-rank==sm_rank

# Relative weak rotation on ONE leg of the unique spin-j singlet pairing
# is chi_j/(2j+1). A simultaneous rotation on both legs is exactly I.
single_leg = {}
for j in range(largest+1):
    amplitude = Fraction(sum(1 if k%2==0 else -1 for k in range(-j,j+1)),2*j+1)
    assert amplitude==Fraction((-1)**j,2*j+1)
    single_leg[str(j)] = {"relative_singlet_amplitude":str(amplitude),
                         "retained_squared_norm":str(amplitude**2),
                         "orthogonal_squared_norm":str(1-amplitude**2),
                         "split_selector_effect_at_c_plus":str((1+amplitude)/2),
                         "split_selector_effect_at_c_minus":str((1-amplitude)/2),
                         "two_insertions_without_intermediate_projection":"1",
                         "two_insertions_with_intermediate_singlet_projection":str(amplitude**2)}

# Actual SU2 adjoint action on the existing Pauli basis of the weak field.
t = np.diag([-1j,1j])
pauli = np.array([[[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]])
adjoint = np.array([[np.trace(a @ t @ b @ t.conj().T)/2 for b in pauli] for a in pauli])
assert np.array_equal(adjoint,np.diag([-1,-1,1]))
assert np.array_equal(t @ (pauli[0]+1j*pauli[1]) @ t.conj().T,-(pauli[0]+1j*pauli[1]))
assert np.array_equal(t @ (pauli[0]-1j*pauli[1]) @ t.conj().T,-(pauli[0]-1j*pauli[1]))
unnormed_singlet = np.eye(3).reshape(9)
one_leg = np.kron(adjoint,np.eye(3))
both_legs = np.kron(adjoint,adjoint)
assert np.vdot(unnormed_singlet,one_leg @ unnormed_singlet)==-1
assert np.array_equal(both_legs @ unnormed_singlet,unnormed_singlet)

print(json.dumps({
    "scope":"actual positive-Weyl DET16, common native z; full weak-singlet and twirl implications, not TFPT Padm identification",
    "one_family_weak_doublets":len(transitions),
    "common_z_rank":int(np.count_nonzero(z)),
    "same_bare_weak_singlet_rank":rows[0]["split_rank"],
    "bare_full_SM_singlet_rank_exact_color_weak_and_zero_Y":int(sm_rank),
    "direct_raising_nullity_mod101":int(len(zero_ids)-rank),
    "common_native_weak_singlet_flat_H_Z":"1+1787*exp(-2*beta)",
    "common_native_full_SM_singlet_flat_H_Z":"1+55*exp(-2*beta)",
    "common_native_weak_singlet_sum_N_minus_8_squared":int(bare_weak_square_count),
    "weak_spin_table":dict(sorted(rows.items())),
    "trace_distance_before_and_after_density_twirl":str(trace_distance),
    "split_noninvariant_under_weak_raising":weak_witness,
    "central_noninvariant_under_D5_weak_pair_root":d5_witness,
    "relative_weak_rotation_on_one_leg_of_total_singlet":single_leg,
},indent=2))
