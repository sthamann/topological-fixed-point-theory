#!/usr/bin/env python3
"""Independent exact audit of physical 5-qubit code and 3-share operator span."""
from itertools import combinations, permutations
import json
import sympy as s

I=s.eye(2); X=s.Matrix([[0,1],[1,0]]); Z=s.diag(1,-1)
paulis={'I':I,'X':X,'Z':Z}
rows=['XZZXI','IXZZX','XIXZZ','ZXIXZ']
Gs=[s.kronecker_product(*(paulis[c] for c in row)) for row in rows]
v0=s.zeros(32,1);v0[0]=1
for G in Gs:v0=(v0+G*v0)/2
v0*=4
v1=s.kronecker_product(*([X]*5))*v0
assert (v0.T*v0)[0]==(v1.T*v1)[0]==1 and (v0.T*v1)[0]==0
P=v0*v0.T+v1*v1.T

def permute(index,perm):
 return sum(((index>>i)&1)<<perm[i] for i in range(5))
physical=[]
for perm in permutations(range(5)):
 ids=[permute(i,perm) for i in range(32)]
 if P.extract(ids,ids)==P:physical.append(perm)
assert len(physical)==10

def bits(index,sites):
 return sum(((index>>site)&1)<<k for k,site in enumerate(sites))
def reduced(M,sites):
 env=tuple(i for i in range(5) if i not in sites)
 out=s.zeros(2**len(sites))
 for i in range(32):
  for j in range(32):
   if M[i,j] and bits(i,env)==bits(j,env):out[bits(i,sites),bits(j,sites)]+=M[i,j]
 return out
ranks=[]
for triple in combinations(range(5),3):
 Es=[reduced(a*b.T,triple) for a in (v0,v1) for b in (v0,v1)]
 M=s.Matrix.hstack(*(E.reshape(64,1) for E in Es))
 rank=M.rank();assert rank==4
 herm=[Es[0]+Es[3],Es[1]+Es[2],s.I*(Es[2]-Es[1]),Es[0]-Es[3]]
 trless=s.Matrix.hstack(*(E.reshape(64,1) for E in herm[1:]))
 assert trless.rank()==3
 assert herm[0]==s.eye(8)/4
 eig=Es[0].eigenvals();assert eig=={s.Rational(1,4):4,s.Integer(0):4}
 ranks.append({'triple':triple,'operator_image_rank':rank,'traceless_image_rank':trless.rank(),
               'pure_logical_reduced_spectrum':{str(k):v for k,v in eig.items()}})
print(json.dumps({'phase_sensitive_physical_code_permutations':len(physical),
 'all_3_share_reductions':ranks,
 'full_3_qubit_operator_space_dimension':64,'full_3_qubit_traceless_dimension':63,
 'actual_encoded_state_operator_span':4,'actual_encoded_traceless_span':3,
 'scope':'Fixed encoded one-qubit operator space. Arbitrary physical errors, syndrome records and extra dynamics enlarge the reachable space and require new assumptions.'},indent=2))
