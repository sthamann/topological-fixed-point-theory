"""Original quaternion/weak action versus the unchanged B7 face map."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse


def load(name,path):
    spec = importlib.util.spec_from_file_location(name,Path(path))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


b7 = load("b7",sys.argv[1])
r032 = load("r032",sys.argv[2])
v252 = r032.native


def q(h1,h2):
    return np.array([[h2.conjugate(),h1],[-h1.conjugate(),h2]],dtype=complex)


def residual(matrix):
    return float(np.max(np.abs(matrix.data if sparse.issparse(matrix) else matrix),initial=0))


def original_field_algebra():
    charts = [(1+0j,0j),(1j,0j),(0j,1+0j),(0j,1j)]
    squares = [q(*h) for h in charts]
    assert all(np.array_equal(a.conj().T@b+b.conj().T@a,
                              2*int(i==j)*np.eye(2))
               for i,a in enumerate(squares) for j,b in enumerate(squares))
    ref = q(0j,1+0j)
    rel1,rel2 = ref.conj().T@squares[0],ref.conj().T@squares[3]
    comm = np.linalg.norm(rel1@rel2-rel2@rel1)
    assert abs(comm-2*np.sqrt(2))<1e-14

    p = v252.rep(1,np.zeros((2,2),complex),np.zeros((3,3),complex))
    k = v252.Kswap()
    dy,dm,gamma = r032.higgs_matrix,r032.mass_matrix,r032.grading
    weak = np.diag([-1 if tag[1]=="L" else 1 for tag in v252.P]*2)
    assert residual(weak@dm@weak-dm)<1e-14
    A = np.linalg.norm(dy)**2
    rows = []
    for h1,h2 in charts+[(0.2+0.3j,0.4-0.1j),(0j,0j)]:
        a = v252.rep(0,q(h1,h2),np.zeros((3,3),complex))
        one = a@(dy@p-p@dy)
        one += one.conj().T
        one += k@one.conj()@k
        phi = dm+one
        h = 0.25*gamma-1j*gamma@phi
        expected = abs(h1)**2+abs(h2)**2
        U = (np.linalg.norm(phi)**2-np.linalg.norm(dm)**2)/A
        err = residual((h-weak@h@weak)/2+1j*gamma@one)
        assert abs(U-expected)<1e-13 and err<1e-13
        rows.append([float(expected),float(U),err])
    return dict(quaternion_anticommutator_exact=True,
                relative_quaternion_commutator_norm=float(comm),
                expected_reconstructed_norm_and_weak_odd_residual=rows)


def native_face_field():
    vertical = [l for l,edge in enumerate(b7.EDGES) if edge.axis==2]
    assert len(vertical)==4
    left = {site:k for k,site in enumerate(b7.PHYSICAL_FACE)}
    right = {site:k for k,site in enumerate(b7.MIRROR_FACE)}

    def block(bits):
        V = np.zeros((4,4))
        for l in vertical:
            edge = b7.EDGES[l]
            assert b7.edge_phase(l,0)==1
            V[left[edge.source],right[edge.target]] = -b7.HOPPING*(-1)**((bits>>l)&1)
        return V

    signed = [block(sum(1<<vertical[j] for j in range(4) if mask>>j&1))
              for mask in range(16)]
    U = [np.trace(V.T@V)/4 for V in signed]
    assert max(U)-min(U)==0 and abs(U[0]-b7.HOPPING**2)<1e-15
    relative = [np.linalg.solve(signed[0],V) for V in signed]
    commutators = [residual(A@B-B@A) for A in relative for B in relative]
    assert max(commutators)==0
    return dict(covariant_vertical_link_configurations=len(signed),
                native_norm_trace=float(U[0]),native_norm_trace_spread=float(max(U)-min(U)),
                max_relative_vertical_commutator=max(commutators))


def gauge_and_windows():
    ga,gb = (b7.TreeGauge(T) for T in (b7.TREE_A,b7.TREE_B))
    occ = b7.EVEN_OCCUPANCIES
    count = len(occ)
    dims = (1<<len(ga.chord_edges))*count
    signs = [np.prod([1-2*((n>>i)&1) for i in b7.PHYSICAL_FACE]) for n in occ]
    W = sparse.diags(np.tile(signs,dims//count),dtype=float,format="csr")
    Hs, Ys, states = [], [], []
    readouts = []
    for tree in (b7.TREE_A,b7.TREE_B):
        H = b7.coupled_hamiltonian(0.8,tree)
        Y = (H-W@H@W)/2
        Y.eliminate_zeros()
        off = Y.tocoo()
        vertical_masks = {1<<edge.source | 1<<edge.target for edge in b7.EDGES if edge.axis==2}
        assert all((occ[i%count]^occ[j%count]) in vertical_masks for i,j in zip(off.row,off.col))
        assert np.all(off.row//count==off.col//count)
        assert max(abs(off.data))==b7.HOPPING
        first_form_counts = []
        second_contact_has_matter = False
        for k in range(5):
            Wk = sparse.diags(np.repeat([(-1)**((q>>k)&1) for q in range(dims//count)],count),
                              dtype=float,format="csr")
            jet = H@Wk-Wk@H
            jet.eliminate_zeros()
            entries = jet.tocoo()
            assert np.all(entries.row//count!=entries.col//count)
            assert np.all(entries.row%count==entries.col%count)
            assert Y.conj().multiply(jet).nnz==0
            first_form_counts.append(jet.nnz)
            second = (H@jet-jet@H).tocoo()
            same_link = second.row//count==second.col//count
            different_matter = second.row%count!=second.col%count
            assert not np.any(same_link&different_matter&(np.abs(second.data)>1e-12))
            second_contact_has_matter |= bool(np.any(
                different_matter&(np.abs(second.data)>1e-12)))
        assert all(n>0 for n in first_form_counts) and second_contact_has_matter
        Hs.append(H)
        Ys.append(Y)
        val, vectors = b7.lowest(H,count=2,vectors=True)
        assert val[1]-val[0]>1e-8
        state = vectors[:,0]
        states.append(state)
        mean_y = np.vdot(state,Y@state)
        U = b7.HOPPING**2*state
        centered_u = U-b7.HOPPING**2*state
        assert residual(centered_u)==0
        readouts.append(dict(ground_energy=float(val[0]),ground_gap=float(val[1]-val[0]),
                             late_vertical_mark_mean=float(mean_y.real),
                             first_link_forms_nonzero_entries=first_form_counts,
                             second_time_contact_changes_matter=second_contact_has_matter,
                             second_time_contact_changes_matter_within_link_block=False,
                             centered_radial_candidate_norm=float(np.linalg.norm(centered_u)),
                             connected_candidate_window_response=0.0))
    rows,phases = [],[]
    for chord in range(dims//count):
        for n in occ:
            q_b,phase = gb.canonicalize(ga.representative_bits(chord),n)
            rows.append(b7.matrix_index(q_b,n))
            phases.append(phase)
    assert len(set(rows))==dims
    S = sparse.csr_matrix((phases,(rows,np.arange(dims))),shape=(dims,dims))
    transport = [residual(Hs[1]@S-S@Hs[0]),residual(Ys[1]@S-S@Ys[0]),residual(W@S-S@W)]
    assert max(transport)==0
    overlap = abs(np.vdot(states[1],S@states[0]))
    assert abs(overlap-1)<1e-12
    assert abs(readouts[0]["late_vertical_mark_mean"]-readouts[1]["late_vertical_mark_mean"])<1e-10
    return dict(transport_residues_H_Y_W=transport,
                independent_ground_state_overlap=float(overlap),
                early_window=[0.05,0.10],late_window=[0.20,0.60],
                gauges=readouts,physical_higgs_yukawa_response=None)


print(json.dumps(dict(original_field=original_field_algebra(),
                      native_face_field=native_face_field(),
                      same_transition=gauge_and_windows()),indent=2))
