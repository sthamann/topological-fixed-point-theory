"""Selector algebra and a diagnostic of field-embedding identifiability.

No choice below is asserted to be the primitive TFPT determinant lift.
The two Klein placements test whether the printed selector/CAR relations
alone decide the mixed determinant-marked field leg.
"""
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import sympy as s

root=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(root/'verification'))
spec=importlib.util.spec_from_file_location('joint_original_v1002',root/'verification/v1002_spin10_mirror_projector.py')
v=importlib.util.module_from_spec(spec);sys.modules[spec.name]=v;spec.loader.exec_module(v)
weights=v.spin10_chiral_weights_doubled();n=len(weights);dim=1<<n
checks={}
def ck(name,val):
 if not val:raise RuntimeError(name)
 checks[name]=True
parity=np.array([1 if int(i).bit_count()%2==0 else -1 for i in range(dim)],dtype=np.int8)
nu=weights.index((1,1,1,1,1))
ck('actual_carrier16',n==16)
ck('actual_singlet_field_color_neutral',(weights[nu][0]+weights[nu][1]-2*weights[nu][2])==0)
for i in range(n):
 a=v.mode_operator(n,i,False).tocoo()
 ck(f'actual_car_parity_anticommutes_mode_{i}',np.all(parity[a.row]+parity[a.col]==0))

# Same fixed determinant involution Q=diag(+1,-1) on a spectral doublet.
# The spectator doublet is an algebraic control, not a selected physical
# theta Hilbert space, state, source weight or TFPT completion.
I=sp.eye(dim,format='csr');X=sp.csr_matrix([[0,1],[1,0]])
Z=sp.csr_matrix([[1,0],[0,-1]])
Q=sp.kron(Z,I,format='csr')
a=v.mode_operator(n,nu,False)
c0=sp.kron(sp.eye(2),a,format='csr')
c1=sp.kron(X,a,format='csr')
def zero(A):
 A.eliminate_zeros()
 return A.nnz==0
ck('fixed_determinant_Q_involution',zero(Q@Q-sp.eye(2*dim)))
ck('field_embedding0_CAR',zero(c0@c0.getH()+c0.getH()@c0-sp.eye(2*dim)))
ck('field_embedding1_CAR',zero(c1@c1.getH()+c1.getH()@c1-sp.eye(2*dim)))
ck('same_even_number_word',zero(c0.getH()@c0-c1.getH()@c1))
ck('embedding0_Q_even',zero(Q@c0@Q-c0))
ck('embedding1_Q_odd',zero(Q@c1@Q+c1))
psi=sp.csr_matrix(([1],([0],[0])),shape=(2*dim,1))
g0=(psi.getH()@c0@Q@c0.getH()@psi)[0,0]
g1=(psi.getH()@c1@Q@c1.getH()@psi)[0,0]
ck('mixed_original_coordinate_leg_separates_embeddings',g0==1 and g1==-1)
ck('ordinary16_determinant_blind_to_central_sign',(-1)**n==1)

# Universal simultaneous decomposition needs only the three commuting
# original projections, not commutation with the complete Z3 generator.
binary=[(p,c,t) for p in [0,1] for c in [0,1] for t in [-1,1]]
bits=[(a,b,c) for a in [0,1] for b in [0,1] for c in [0,1]]
F=np.array([[(2*p-1)**a*(2*c-1)**b*t**d for a,b,d in bits]
            for p,c,t in binary],dtype=np.int64)
ck('universal_eight_component_Fourier_inverse',np.array_equal(F@F.T,8*np.eye(8,dtype=int)))
avg=F.sum(axis=1)//8
ck('universal_joint_projector_selects_original_sector',
   avg.tolist()==[int(p==1 and c==1 and t==1) for p,c,t in binary])

# Universal joint functional calculus of the three original commuting
# selectors: Pprim, determinant involution, color center. This is its
# character algebra, not a candidate physical realization.
z=(-1+s.sqrt(3)*s.I)/2
labels=[(p,t,r) for p in [0,1] for t in [-1,1] for r in range(3)]
E=s.zeros(12)
for a,(pa,ta,ra) in enumerate(labels):
 for b,(pb,tb,rb) in enumerate(labels):
  color=s.simplify(sum(z**(k*(rb-ra)) for k in range(3))/3)
  E[a,b]=s.simplify((pb if pa else 1-pb)*(1+ta*tb)*color/2)
ck('joint_spectral_partition_is_exact',E==s.eye(12))

out={'result_id':'CRAIG20261103D','passed':len(checks),'checks':checks,
 'actual_original_mode_count':n,'neutral_coordinate_index':nu,
 'mixed_determinant_marked_coordinate_legs':[int(g0.real),int(g1.real)],
 'universal_joint_selector_labels':[list(x) for x in binary],
 'joint_selector_labels':[list(x) for x in labels],
 'scope':'Joint selector functional calculus; exact CAR/color/even-word diagnostic only. No two complete primitive theories claimed; actual joint bosonic algebra, time, original state, real structure and determinant-line lift can decide the embedding and are not assumed equal.',
 'full_original_source_integral':'not_evaluated',
 'full_primitive_TFPT_kernel':'not_derived'}
Path(__file__).with_name('joint-selector-results.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))
