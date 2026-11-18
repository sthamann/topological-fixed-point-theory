"""Direct use of unchanged original native source; no new physical model."""
import sys,importlib.util,json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply
root=Path(__file__).resolve().parents[3]
if not (root/'verification').is_dir():root=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
p=root/'experiments/tfpt-discovery/tfpt3p1d_coupled_wall_probe.py'
spec=importlib.util.spec_from_file_location('native_coupled_source',p)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
g=0.8 # unchanged first declared source coupling, not a derived primitive value
H=m.coupled_hamiltonian(g); gauge=m.TreeGauge(m.TREE_A)
f=len(m.EVEN_OCCUPANCIES);nb=H.shape[0]//f;dim=H.shape[0]
def maxabs(A):return float(np.max(np.abs(A.data))) if A.nnz else 0.0
Xs=[]
for ell in range(len(m.EDGES)):
 rows=[];cols=[];data=[]
 for q in range(nb):
  bits=gauge.representative_bits(q)
  for occ in m.EVEN_OCCUPANCIES:
   qq,sgn=gauge.canonicalize(bits^(1<<ell),occ)
   cols.append(m.matrix_index(q,occ));rows.append(m.matrix_index(qq,occ));data.append(sgn)
 Xs.append(sparse.csr_matrix((data,(rows,cols)),shape=H.shape,dtype=complex))
I=sparse.eye(dim,format='csr');c=g*g/2
HE=c*(len(Xs)*I-sum(Xs));R=(H-HE).tocsr()
coo=R.tocoo();off=coo.data[(coo.row//f)!=(coo.col//f)]
checks={
 'dimension':dim,'boson_coordinate_values':nb,'retained_inner_dimension':f,
 'H_nnz':H.nnz,'hermiticity':maxabs(H-H.getH()),
 'electric_involution':max(maxabs(X@X-I) for X in Xs),
 'electric_commutation':max(maxabs(X@Y-Y@X) for i,X in enumerate(Xs) for Y in Xs[i+1:]),
 'remaining_cross_boson_block':float(np.max(np.abs(off))) if len(off) else 0.,
 'operator_reconstruction':maxabs(H-R-HE)}
if any(checks[k]>1e-12 for k in ['hermiticity','electric_involution','electric_commutation','remaining_cross_boson_block','operator_reconstruction']):raise RuntimeError(checks)
evals=[];evecs=[]
for q in range(nb):
 w,V=eigh(R[q*f:(q+1)*f,q*f:(q+1)*f].toarray());evals.append(w);evecs.append(V)
def exp_R(s,v):
 v=v.reshape(nb,f);u=np.empty_like(v)
 for q,(w,V) in enumerate(zip(evals,evecs)):
  u[q]=V@(np.exp(-s*w)*(V.conj().T@v[q]))
 return u.ravel()
def exp_E(s,v):
 u=v.copy();ch=np.cosh(c*s);sh=np.sinh(c*s)
 for X in Xs:u=ch*u+sh*(X@u)
 return np.exp(-c*s*len(Xs))*u
def slice_evolution(t,n,v):
 s=t/n;u=v.copy()
 for _ in range(n):u=exp_R(s/2,exp_E(s,exp_R(s/2,u)))
 return u
# Deterministic original physical field-basis leg. This is a representation
# check on a whole vector, not selection of the physical OS vacuum.
v=np.zeros(dim,dtype=complex);v[0]=1
t=0.2;exact=expm_multiply(-t*H,v)
errors=[]
for n in [4,8,16]:
 u=slice_evolution(t,n,v)
 errors.append({'slices':n,'whole_vector_error':float(np.linalg.norm(u-exact))})
if not all(errors[i+1]['whole_vector_error']<errors[i]['whole_vector_error']/3.5 for i in range(2)):raise RuntimeError(errors)
# Nontrivial actual boson field and even matter field marks. Source comparison
# uses the same endpoints, time gaps and absolute amplitudes on both paths.
markA=np.repeat(np.array([1. if (q&1)==0 else -1. for q in range(nb)]),f)
markB=np.tile(np.array([float((occ>>0)&1) for occ in m.EVEN_OCCUPANCIES]),nb)
exact_word=expm_multiply(-t*H,markA*expm_multiply(-t*H,markB*v))
word_errors=[]
# use A then B on nonzero propagated vector, to avoid a trivially zero source
exact_word=expm_multiply(-t*H,markB*expm_multiply(-t*H,markA*expm_multiply(-t*H,v)))
for n in [4,8,16]:
 u=slice_evolution(t,n,markB*slice_evolution(t,n,markA*slice_evolution(t,n,v)))
 word_errors.append({'slices_per_gap':n,'marked_word_vector_error':float(np.linalg.norm(u-exact_word))})
if not all(word_errors[i+1]['marked_word_vector_error']<word_errors[i]['marked_word_vector_error']/3.5 for i in range(2)):raise RuntimeError(word_errors)
b0=H[6*f:7*f,0:f].toarray();diag=np.diag(b0)
checks['electric_block_q0_to_q6_diagonal_values']=sorted(set(np.round(diag.real,14)))
checks['empty_occupation_scalar_replacement_error']=float(np.linalg.norm(b0-diag[0]*np.eye(f),2))
result={'scope':'Exact original finite native coupled cube source; g0.8 is existing model input. No primitive P1/Sadm identity or physical spectrum inferred. Direct matrix integration with all inner marks, no scalar boson replacement.','checks':checks,'time_diagnostic':t,'trotter_errors':errors,'marked_word_errors':word_errors,'primitive_source_identification_proved':False}
Path(__file__).with_name('native-matrix-source-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
