"""B7 native physical cube: operator identification and normalized time control."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import expm_multiply

path = Path(sys.argv[1])
spec = importlib.util.spec_from_file_location("b7", path)
native = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = native
spec.loader.exec_module(native)

coupling = 0.8  # Frozen B7 test point.
tree_a, tree_b = native.TREE_A, native.TREE_B
gauge_a, gauge_b = native.TreeGauge(tree_a), native.TreeGauge(tree_b)
H_a = native.coupled_hamiltonian(coupling, tree_a)
H_b = native.coupled_hamiltonian(coupling, tree_b)
dim = H_a.shape[0]
assert dim == 4096
cycles = 1 << len(gauge_a.chord_edges)
rows, signs = [], []
for chord in range(cycles):
    full = gauge_a.representative_bits(chord)
    for n in native.EVEN_OCCUPANCIES:
        target, sign = gauge_b.canonicalize(full,n)
        rows.append(native.matrix_index(target,n))
        signs.append(sign)
assert len(set(rows)) == dim
transport = sparse.csr_matrix((signs,(rows,np.arange(dim))),shape=(dim,dim))


def max_entry(A):
    return float(np.max(np.abs(A.data),initial=0))


h_residue = max_entry(H_b@transport-transport@H_a)
assert h_residue < 1e-12
numbers = np.tile([n.bit_count() for n in native.EVEN_OCCUPANCIES],cycles)
face_numbers = np.tile([(n&native.PHYSICAL_MASK).bit_count()
                        for n in native.EVEN_OCCUPANCIES],cycles)
entries = H_a.tocoo()
commutator_total = entries.data*(numbers[entries.col]-numbers[entries.row])
commutator_face = entries.data*(face_numbers[entries.col]-face_numbers[entries.row])
total = float(np.max(abs(commutator_total)))
face = float(np.max(abs(commutator_face)))
assert abs(total-2*native.DET_STRENGTH)<1e-12
assert abs(face-native.HOPPING)<1e-12
vacuum = native.matrix_index(0,0)
full_mirror = native.matrix_index(0,native.MIRROR_MASK)
direct_mirror = abs(H_a[full_mirror,vacuum])
assert abs(direct_mirror-native.DET_STRENGTH/2)<1e-12
occupancies = len(native.EVEN_OCCUPANCIES)
matter = H_a[:occupancies,:occupancies].copy()
matter.setdiag(0)
matter.eliminate_zeros()
component_count, components = connected_components(abs(matter),directed=False)
assert component_count==2
component_sizes = sorted(int(np.sum(components==j)) for j in range(component_count))
assert component_sizes==[56,72]
for j in range(component_count):
    residues = set(numbers[:occupancies][components==j]%4)
    assert len(residues)==1


def native_marks(tree):
    gauge = native.TreeGauge(tree)
    xy = next(p for p in native.PLAQUETTES if p.axes==(0,1) and p.base[2]==0)
    w = [(-1)**((gauge.representative_bits(chord)&xy.mask).bit_count())
         for chord in range(cycles)]
    W = sparse.diags(np.repeat(w,len(native.EVEN_OCCUPANCIES)),format="csr",dtype=float)
    full = native.coupled_hamiltonian(coupling,tree)
    without_hop = native.coupled_hamiltonian(coupling,tree,hopping=0.0)
    B = (without_hop-full)/native.HOPPING
    assert max_entry(B-B.getH())<1e-12
    return W,B


def ground(H):
    values, states = native.lowest(H,count=2,vectors=True)
    assert values[1]-values[0]>1e-8
    state = states[:,0]
    resid = np.linalg.norm(H@state-values[0]*state)
    assert resid<1e-8
    return float(values[0]),float(values[1]-values[0]),state


def normalized_window(H,E0,state,A,B,order=12):
    # Keep the same fractions as the previous two-window test, now in B7 time.
    early,late=(0.05,0.10),(0.20,0.60)
    A0,B0=A@state,B@state
    A0-=np.vdot(state,A0)*state
    B0-=np.vdot(state,B0)*state
    M=-(H-E0*sparse.eye(H.shape[0],format="csr"))
    nodes,weights=np.polynomial.legendre.leggauss(order)
    # e^{(s-t)M} splits into two commuting factors about the first late time.
    def vector(window,shift,direction,v):
        a,b=window
        times=(a+b)/2+(b-a)*nodes/2
        return sum(w*expm_multiply(direction*(t-shift)*M,v)
                   for t,w in zip(times,weights))/2
    early_vector=vector(early,late[0],-1,A0)
    late_vector=vector(late,late[0],1,B0)
    return np.vdot(late_vector,early_vector)


W_a,B_a=native_marks(tree_a)
W_b,B_b=native_marks(tree_b)
mark_residue=max(max_entry(W_b@transport-transport@W_a),max_entry(B_b@transport-transport@B_a))
assert mark_residue<1e-12
E_a,gap_a,state_a=ground(H_a)
E_b,gap_b,state_b=ground(H_b)
overlap=abs(np.vdot(state_b,transport@state_a))
assert abs(overlap-1)<1e-10
assert max(abs(E_a-E_b),abs(gap_a-gap_b))<1e-8
result_a=normalized_window(H_a,E_a,state_a,W_a,B_a)
result_b=normalized_window(H_b,E_b,state_b,W_b,B_b)
coarse=normalized_window(H_a,E_a,state_a,W_a,B_a,order=6)
assert abs(result_a-result_b)<1e-9
assert abs(result_a-coarse)<1e-10
print(json.dumps(dict(
    cube_dim=dim,coupling=coupling,native_ground_energy=E_a,native_gap=gap_a,
    exact_tree_hamiltonian_residual=h_residue,tree_native_marks_residual=mark_residue,
    ground_overlap=overlap,
    total_number_commutator_max=total,physical_face_number_commutator_max=face,
    mirror_vacuum_full_transition=float(direct_mirror),
    occupation_components_sizes=component_sizes,
    native_wilson_hopping_window_response_tree_A=[result_a.real,result_a.imag],
    native_wilson_hopping_window_response_tree_B=[result_b.real,result_b.imag],
    quadrature_difference=abs(coarse-result_a),physical_higgs_yukawa_response=None
),indent=2))
