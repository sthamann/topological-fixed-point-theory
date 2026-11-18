"""Classify conserved B7 operators preserving native link coordinates."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csgraph

source_path = Path(sys.argv[1])
spec = importlib.util.spec_from_file_location("b7",source_path)
native = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = native
spec.loader.exec_module(native)


def certificate(tree):
    gauge = native.TreeGauge(tree)
    H = native.coupled_hamiltonian(0.8,tree)
    occupations = np.array(native.EVEN_OCCUPANCIES)
    matter_dim = len(occupations)
    link_dim = 1 << len(gauge.chord_edges)

    # Each electric move connects a distinct link block. Its gauge restoration
    # acts on occupations by the original parity sign.
    moves, chars = [], []
    for edge in range(len(native.EDGES)):
        cases = [gauge.canonicalize(1<<edge,int(n)) for n in occupations]
        assert len({q for q,_ in cases})==1
        moves.append(cases[0][0])
        chars.append(np.array([sign for _,sign in cases]))
    assert len(set(moves))==12 and 0 not in moves
    assert set(moves[j] for j in gauge.chord_edges)=={1<<k for k in range(5)}
    assert all(np.all(chars[j]==1) for j in gauge.chord_edges)

    # Confirm directly that the full original off-diagonal link blocks contain
    # exactly these electric moves, with no hidden matter transition.
    blocks = H.tocoo()
    q_in,n_in = divmod(blocks.col,matter_dim)
    q_out,n_out = divmod(blocks.row,matter_dim)
    mask = q_in!=q_out
    assert int(mask.sum())==link_dim*matter_dim*len(moves)
    assert np.array_equal(n_in[mask],n_out[mask])
    edge_for_move = np.full(link_dim,-1)
    edge_for_move[moves] = np.arange(len(moves))
    expected = -0.5*0.8**2*np.array(chars)[edge_for_move[q_in[mask]^q_out[mask]],n_in[mask]]
    assert np.all(edge_for_move[q_in[mask]^q_out[mask]]>=0)
    electric_residue = float(np.max(abs(blocks.data[mask]-expected)))
    assert electric_residue<1e-14

    # Chord flips imply a single common matter matrix. Tree-restoration
    # signs distinguish every one of the 128 occupation vectors.
    joint_character = np.array(chars)[list(tree)].T
    n_signs = np.unique(joint_character,axis=0).shape[0]
    assert n_signs==matter_dim==128

    first = H[:matter_dim,:matter_dim].copy()
    first.setdiag(0)
    first.eliminate_zeros()
    count, labels = csgraph.connected_components(abs(first),directed=False)
    number = np.array([int(n).bit_count() for n in occupations])
    assert count==2
    residues = []
    for i in range(count):
        unique = np.unique(number[labels==i]%4)
        assert len(unique)==1
        residues.append([int(unique[0]),int(np.sum(labels==i))])
    assert sorted(residues)==[[0,72],[2,56]]

    # Number of a subsystem with m active coordinate modes, with total
    # parity fixed. Include either parity to cover particle-hole charts.
    incompatible = {}
    for parity in (0,1):
        sector = [n for n in range(256) if n.bit_count()%2==parity]
        for m in range(1,9):
            values = [(n&((1<<m)-1)).bit_count() for n in sector]
            ranks = np.unique(values,return_counts=True)[1].tolist()
            assert sorted(ranks)!=[56,72]
            incompatible[f"parity_{parity}_modes_{m}"]=ranks

    values, vectors = native.lowest(H,count=2,vectors=True)
    assert values[1]-values[0]>1e-8
    state = vectors[:,0]
    grade = np.tile(np.where(number%4==0,1.0,-1.0),link_dim)
    commutator = float(np.max(abs(blocks.data*(grade[blocks.row]-grade[blocks.col]))))
    assert commutator==0
    grade_mean = np.vdot(state,grade*state)
    centered_residual = float(np.linalg.norm((grade-grade_mean)*state))
    assert abs(abs(grade_mean)-1)<1e-10 and centered_residual<1e-9
    return dict(dimension=H.shape[0],electric_moves=moves,
                electric_transport_residual=electric_residue,
                distinct_joint_parity_characters=n_signs,
                matter_commutant_component_sizes_by_number_mod4=sorted(residues),
                other_number_spectra_ranks=incompatible,
                ground_energy=float(values[0]),ground_gap=float(values[1]-values[0]),
                ground_N_mod4=0 if grade_mean.real>0 else 2,
                N_mod4_commutator_residual=commutator,
                centered_residual_for_N_mod4=centered_residual,
                connected_ground_window_response_for_any_conserved_coordinate_mark=0,
                physical_higgs_yukawa_response=None)


a = certificate(native.TREE_A)
b = certificate(native.TREE_B)
assert a["ground_N_mod4"]==b["ground_N_mod4"]
assert abs(a["ground_energy"]-b["ground_energy"])<1e-9

# For any two centered ground-state vectors the native gap gives a
# dimensionless bound for the normalized window average, not their norms.
d = a["ground_gap"]
window_decay = (-np.expm1(-d*0.05)/(d*0.05)
                * -np.expm1(-d*0.40)/(d*0.40)*np.exp(-d*0.10))
assert 0 < window_decay < 1
print(json.dumps(dict(tree_A=a,tree_B=b,native_window_decay_bound=window_decay),indent=2))
