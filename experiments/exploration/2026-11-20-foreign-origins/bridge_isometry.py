#!/usr/bin/env python3
"""Conditionally identify the five facets of a 4-simplex with five qubit legs.
Verify that every geometrically allowed 2|3 cut yields an exact isometry.
This identification is an extra hypothesis, not a TFPT derivation.
"""
from itertools import combinations
from sympy import Matrix, eye, simplify
from five_qubit_qec import GENS, exact_stabilizer_state
import json


def reshape_for_cut(psi, positions):
    positions=tuple(positions)
    outpos=tuple(i for i in range(5) if i not in positions)
    M=Matrix.zeros(8,4)
    for idx in range(32):
        bits=tuple((idx>>(4-i))&1 for i in range(5))
        inp=sum(bits[p]<<(1-j) for j,p in enumerate(positions))
        out=sum(bits[p]<<(2-j) for j,p in enumerate(outpos))
        M[out,inp]=psi[idx]
    return M


def test():
    psi=exact_stabilizer_state(GENS)
    verified=0
    projection_ranks=set()
    for two in combinations(range(5),2):
        M=reshape_for_cut(psi,two)
        V=2*M
        assert V.conjugate().T*V == eye(4)
        P=V*V.conjugate().T
        assert P*P==P
        projection_ranks.add(P.rank())
        assert P!=eye(8)
        verified+=1
    assert verified==10 and projection_ranks=={4}
    return {'Pachner_facet_2_to_3_splits':10,
            'all_ten_AMEmaps_isometric':True,
            'input_Hilbert_dimension':4,'output_Hilbert_dimension':8,
            'physical_output_projector_rank':4,
            'full_output_unitarity':False,
            'extra_assumption':'facets of a 4-simplex = five physical qubits carrying this stabilizer state',
            'critical_missing_condition':'Pachner pentagon state sum with coherent amplitudes and positive physical time'}


if __name__=='__main__':print(json.dumps(test(),indent=2,ensure_ascii=False))
