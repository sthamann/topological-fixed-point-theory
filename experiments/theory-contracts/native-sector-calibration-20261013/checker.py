from pathlib import Path
import importlib.util,sys,json,hashlib
import numpy as np
import sympy as sp
REPO=Path(__file__).resolve().parents[3]
source=REPO/'verification/v1033_charged_disorder.py'
sys.path.insert(0,str(REPO/'verification'))
spec=importlib.util.spec_from_file_location('v1033_calibration_source',source)
mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
n=6;d=2**n;z=1j
cs=[]
for j in range(n):
 c=np.zeros((d,d),complex)
 for k in range(d):
  if k&(1<<j):c[k^(1<<j),k]=(-1)**((k&((1<<j)-1)).bit_count())
 cs.append(c)
num=np.array([k.bit_count() for k in range(d)])
S=np.roll(np.eye(4,dtype=complex),1,axis=0)
Z=np.diag(z**np.arange(4));q=np.array([1,1,0,0,0,0])
Gamma=np.diag([z**sum(q[j] for j in range(n) if k&(1<<j)) for k in range(d)])
T=np.zeros((4*d,4*d),complex);H=T.copy()
for r in range(4):
 h=mod.qwz_cylinder(3,1,1.0,r)
 hr=sum((h[a,b]*(cs[a].conj().T@cs[b]) for a in range(n) for b in range(n)),np.zeros((d,d),complex))
 H[r*d:(r+1)*d,r*d:(r+1)*d]=hr
 T[r*d:(r+1)*d,r*d:(r+1)*d]=np.linalg.matrix_power(Gamma,r)
C=np.kron(Z,np.diag(z**num));D=np.kron(S,Gamma)
Fs=[T@np.kron(S,c)@T.conj().T for c in cs]
Ps=[np.diag(np.array([int((r+num[k])%4==s) for r in range(4) for k in range(d)])) for s in range(4)]
cs4=np.array([0,.5,1.75,-1.5]);delta=sum((v*P for v,P in zip(cs4,Ps)),np.zeros_like(H))
comm=lambda A,B:A@B-B@A
norm=lambda A:float(np.max(np.abs(A)))
K=D.conj().T@comm(delta,D)
b=np.array([float(np.trace(P@K).real/d) for P in Ps]);expected=np.roll(cs4,-1)-cs4
B=sp.zeros(4)
for s in range(4):B[s,s]=-1;B[s,(s+1)%4]=1
B2=sp.zeros(4)
for s in range(4):B2[s,s]=-1;B2[s,(s+2)%4]=1
recovered=np.r_[0,np.cumsum(b[:3])];recovered-=recovered.mean()
# A constructive response protocol: use the existing source ground and only
# its translates by native D. Preparation and complex response access are
# stated conditions, not a P1-derived laboratory implementation.
eigvals,eigvecs=np.linalg.eigh(H)
psi=eigvecs[:,0];s0=int(np.argmax([np.linalg.norm(P@psi) for P in Ps]))
ground_sector_leakage=float(np.linalg.norm(psi-Ps[s0]@psi))
response_by_sector={}
for k in range(4):
 v=np.linalg.matrix_power(D,k)@psi
 s=(s0+k)%4
 response_by_sector[s]=float(np.vdot(v,K@v).real)
protocol_error=max(abs(response_by_sector[s]-expected[s]) for s in range(4))
checks={
 'source_H_preserves_C':norm(comm(H,C))==0,
 'source_D_shifts_C':norm(C@D-z*D@C)==0,
 'sector_offsets_preserve_C':norm(comm(delta,C))==0,
 'sector_offsets_preserve_all_completed_CAR_generators':max(norm(comm(delta,F)) for F in Fs)==0,
 'native_D_detects_nonconstant_offsets':norm(comm(delta,D))>0,
 'source_operator_recovers_exact_adjacent_differences':bool(np.array_equal(b,expected)),
 'four_cycle_consistency_sum_zero':float(b.sum())==0,
 'exact_offset_reconstruction_mod_global_scalar':bool(np.array_equal(recovered,cs4-cs4.mean())),
 'native_D_rank_three_kernel_dimension_one':B.rank()==3 and len(B.nullspace())==1,
 'keeping_D_squared_only_leaves_two_constants':B2.rank()==2 and len(B2.nullspace())==2,
 'native_D_translated_ground_protocol_recovers_all_offsets':protocol_error<1e-12,
 'ground_is_in_reported_C_sector':ground_sector_leakage<1e-12,
 'ground_source_matches_prior_E_minus4_unique':abs(eigvals[0]+4)<1e-12 and eigvals[1]-eigvals[0]>1e-8,
 'global_scalar_stays_invisible_to_full_source':norm(comm(np.eye(4*d),D))==0 and all(norm(comm(np.eye(4*d),F))==0 for F in Fs)
}
checks={key:bool(value) for key,value in checks.items()}
out={'research_id':'UR.SOURCE.NATIVE_SECTOR_CALIBRATION.01','scope':'Existing v1033 nx=3 ny=1 mass=+1; previously declared conditional C4 tensor Fock lift and completed CAR dictionary; no new dynamics or physical-source selection.', 'source_path':str(source.relative_to(REPO)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'dimension':4*d,'offsets':cs4.tolist(),'measured_adjacent_differences':b.tolist(),'centered_reconstruction':recovered.tolist(),'native_cycle_rank':int(B.rank()),'native_cycle_nullspace':[[str(v) for v in x] for x in B.nullspace()],'D_squared_rank':int(B2.rank()),'max_CAR_commutator_error':max(norm(comm(delta,F)) for F in Fs),'native_D_commutator_max':norm(comm(delta,D)),'checks':checks,'passed':sum(checks.values()),'failed':sum(not x for x in checks.values()),'physical_goals_closed':[],'native_response_protocol':{'initial_state':'unique ground of original H; same fixed preparations compared, not reselected shifted ground','initial_C_sector':s0,'ground_energy':float(eigvals[0]),'ground_sector_leakage':ground_sector_leakage,'four_preparations':'D^k psi for k=0,1,2,3','required_readout':'baseline-subtracted complex first derivative [G_prime(0;H_c)-G_prime(0;H)]/i in each fixed translated preparation; known native H response is required','difference_by_sector':response_by_sector,'max_error':protocol_error}}
print(json.dumps(out,indent=2));assert all(checks.values())
