"""Read the original W coupling and test total family transport on Q=2 and Q=4."""
from fractions import Fraction
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.sparse import csr_matrix

repo=Path(sys.argv[1]).resolve()
common=runpy.run_path(str(repo/'experiments/theory-contracts/universalraum-native-operations-ground-response-20260915/common.py'))
W=common['load_tensor']()
Ws=csr_matrix(W)
PAIRS=common['PAIRS']
pair_ix=common['PAIR_INDEX']
data={}
for i in range(4):
    for j in range(4):
        if i==j:
            continue
        X=np.zeros((64,64),dtype=np.int64)
        X[4*np.arange(16)+i,4*np.arange(16)+j]=1
        Nf=common['exterior_square_derivation'](X)
        assert Nf.shape==(2016,2016)
        y=common['boson_generator'](W,X)
        assert np.all(y.imag==0) and np.array_equal(y.real,np.rint(y.real))
        Nb=np.rint(y.real).astype(np.int64)
        assert np.array_equal((Ws@Nf).toarray(),Nb@W)
        assert np.array_equal((Ws@Nf.T).toarray(),Nb.T@W)
        assert (Nb@Nb==0).all()
        changed=np.flatnonzero(np.count_nonzero(Nb,axis=0))
        assert len(changed)==20
        assert np.all(np.count_nonzero(Nb[:,changed],axis=0)==1)
        assert np.all(np.max(abs(Nb[:,changed]),axis=0)==1)
        fermion_triplet_count=0
        for a,b in PAIRS:
            if a%4!=j or b%4!=j:
                continue
            k=pair_ix[(a,b)]
            n0=csr_matrix(([1],([k],[0])),shape=(2016,1))
            n1,n2=Nf@n0,Nf@(Nf@n0)
            assert [int(u.conjugate().multiply(u).sum().real) for u in (n0,n1,n2)]==[1,2,4]
            assert all((Ws@u).nnz==0 for u in (n0,n1,n2))
            assert (Nf@n2).nnz==0
            fermion_triplet_count+=1
        assert fermion_triplet_count==120
        # Q=4 normalized two-boson states aa, ab, bb in a doublet
        # with Nb e_a = +/- e_b. Boson lowering gives the stated
        # norms and couples via the exact two fermion pair W row.
        for A in changed:
            B=int(np.flatnonzero(Nb[:,A])[0])
            assert not np.any(Nb[:,B])
            assert W[A]@W[A]==W[B]@W[B]==8 and W[A]@W[B]==0
            assert A!=B
        data[f'{i+1}<-{j+1}']={
            'fermion_triplets_Q2':fermion_triplet_count,
            'every_tested_Q2_triplet_annihilated_by_original_W':True,
            'one_boson_doublets':len(changed),
            'example_original_boson_rows':[int(changed[0]),int(np.flatnonzero(Nb[:,changed[0]])[0])],
            'boson_triplet_Q4_0_1_2_power_norm2':[1,2,4],
            'interaction_action_squared_norm_on_each_Q4_normalized_triplet_state_in_units_g_squared':16,
            'total_charge4_current_commutes_with_original_H_by_W_intertwiner':True,
        }

# Both unnormalized readouts use one and the same Q=2 block. The
# 60 singular channels have coupling sqrt(8)*g. The other 1956
# fermion pair directions have zero energy. Use the declared |g|/Delta
# range; the numerical point below is an illustrative value, not a fit.
Delta=1.
g=1/20
beta=1.
split=np.sqrt(Delta**2+32*g*g)
energies=((Delta-split)/2,(Delta+split)/2)
coupled=60*sum(np.exp(-beta*e) for e in energies)
empty=1956+coupled
reference=2016+60*np.exp(-beta*Delta)
p1_counts=json.loads((base:=repo/'experiments/theory-contracts/chat-triangulation-20261108/tfpt-p1-closed-return-results-20261110.json').read_text())['closed_original']['old_P1_same_reference']
ident={}
for a,v in p1_counts.items():
    w=[Fraction(x) for x in v['Cartan_count_k0_through_k3_weight_over_input']]
    contrast=(1,-3,3,-1)
    eps=min(x/(2*abs(c)) for x,c in zip(w,contrast,strict=True))
    plus=[x+eps*c for x,c in zip(w,contrast,strict=True)]
    minus=[x-eps*c for x,c in zip(w,contrast,strict=True)]
    assert min(plus+minus)>0
    mom=lambda u,k: sum(x*i**k for i,x in enumerate(u))
    assert all(mom(plus,k)==mom(minus,k)==mom(w,k) for k in range(3))
    parity=lambda u:sum(x*(-1)**i for i,x in enumerate(u))
    assert parity(plus)-parity(minus)==16*eps
    ident[a]={'epsilon':str(eps),'original_unmarked_and_first_two_moments_preserved':True,
              'parity_difference_exact':str(16*eps),'parity_difference':float(16*eps)}

print(json.dumps({'original_W':data,
                  'absolute_conditional_Q2':{'Delta':Delta,'g':g,'beta':beta,
                      'empty_Z':empty,'coupled_marker_Z':coupled,'fixed_g0_reference':float(reference),
                      'empty_relative':float(empty/reference),'coupled_marker_relative':float(coupled/reference),
                      'coupled_normalized':float(coupled/empty)},
                  'P1_cartan_identifiability':ident,
                  'physical_original_TFPT_values':None},indent=2))
