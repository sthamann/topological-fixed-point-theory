"""Exact linear family certificate for the actual v252 compatibility checks.

This audits the search objective; it does not claim to impose all TFPT physics.
"""
from pathlib import Path
import sys,json,time
import numpy as np
import sympy as s
R=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(R/'verification'))
import v252_full_finite_triple as f
start=time.time();z=np.zeros((3,3),complex);I=np.eye(96)
g=f.grading();swap=f.Kswap();perm=np.argmax(swap,axis=1)
checks={};basis=[];coordinate_labels=[]
for sector in range(5):
 for i in range(3):
  for j in range(3):
   if sector==4 and j<i:continue
   for phase in (1,1j):
    args=[z.copy() for _ in range(5)];args[sector][i,j]=phase
    if sector==4:args[sector][j,i]=phase
    basis.append(f.full_D(*args));coordinate_labels.append([['Yu','Yd','Ye','Ynu','MR'][sector],i,j,'real' if phase==1 else 'imag'])
D=np.array(basis);assert len(basis)==84
checks['entire84_family_hermitian']=np.array_equal(D,D.conj().transpose(0,2,1))
checks['entire84_family_gamma_odd']=not np.any(D*(np.diag(g)[None,:,None]+np.diag(g)[None,None,:]))
checks['entire84_family_J_real']=np.array_equal(D.conj()[:,perm][:,:,perm],D)
def coords(mat):
 v=[]
 for name in ['u','d','e','nu']:
  for i in range(3):
   for j in range(3):
    x=mat[f.IDX[name,'L',0,i],f.IDX[name,'R',0,j]];v.extend([x.real,x.imag])
 for i in range(3):
  for j in range(i,3):
   x=mat[f.IDX['nu','R',0,i],48+f.IDX['nu','R',0,j]];v.extend([x.real,x.imag])
 return np.array(v)
checks['84_independent_original_coordinates']=np.array_equal(np.array([coords(b) for b in basis]),np.eye(84))
# Algebra acts identically on generations. Every ordered family block can
# therefore be certified on the actual32-slot module, with no discarded modes.
slots=[list(range(16*a,16*(a+1)))+list(range(48+16*a,48+16*(a+1))) for a in range(3)]
alg=[]
for phase in (1,1j):alg.append(f.rep(phase,np.zeros((2,2)),z))
for a,b in [(1,0),(1j,0),(0,1),(0,1j)]:alg.append(f.rep(0,f.quaternion(a,b),z))
for i in range(3):
 for j in range(3):
  for phase in (1,1j):
   m=z.copy();m[i,j]=phase;alg.append(f.rep(0,np.zeros((2,2)),m))
checks['actual_algebra24_family_diagonal']=all(not np.any(a[np.ix_(slots[i],slots[j])]) for a in alg for i in range(3) for j in range(3) if i!=j)
checks['actual_algebra24_family_blind']=all(np.array_equal(a[np.ix_(slots[0],slots[0])],a[np.ix_(slots[j],slots[j])]) for a in alg for j in (1,2))
left=np.array([a[np.ix_(slots[0],slots[0])] for a in alg]);right=np.array([(swap@a.conj()@swap)[np.ix_(slots[0],slots[0])] for a in alg])
checks['SM_order_zero_all24_by24']=not np.any(left[:,None]@right[None,:]-right[None,:]@left[:,None])
unique={}
for d in basis:
 for si in slots:
  for sj in slots:
   b=d[np.ix_(si,sj)]
   if np.any(b):unique[b.tobytes()]=b
violations=[]
for n,b in enumerate(unique.values()):
 comm=b@left-left@b
 residual=comm[:,None]@right[None,:]-right[None,:]@comm[:,None]
 if np.any(residual):violations.append(n)
checks['SM_order_one_entire84_linear_family']=not violations
print('All84 input directions checked, ordered family block types',len(unique),'elapsed',round(time.time()-start,2),flush=True)
# All calculations above involve Gaussian integers and sums/products well
# within2^53. Certify integrality rather than treating this as tolerance math.
checks['original_basis_gaussian_integer_exact']=all(np.array_equal(v.real,np.rint(v.real)) and np.array_equal(v.imag,np.rint(v.imag)) for v in [D,left,right])
# Reuse the actual six-multiplet coordinate-frame construction.
positive={}
for i in range(96):
 anti=i>=48;name,ch,col,family=f.P[i%48]
 if g[i,i]!=1:continue
 label=('Q' if name in ('u','d') else 'L') if not anti else {'u':'uc','d':'dc','e':'ec','nu':'nc'}[name]
 positive.setdefault(label,{}).setdefault((name,col),{})[family]=i
H=[];Hlabels=[]
for i in range(3):
 h=z.copy();h[i,i]=1;H.append(h);Hlabels.append(f'diag{i}')
for i in range(3):
 for j in range(i+1,3):
  h=z.copy();h[i,j]=h[j,i]=1;H.append(h);Hlabels.append(f'sym{i}{j}')
  h=z.copy();h[i,j]=1j;h[j,i]=-1j;H.append(h);Hlabels.append(f'imag{i}{j}')
xgens=[];glabels=[]
for label,components in positive.items():
 for h,hl in zip(H,Hlabels):
  xp=np.zeros((96,96),complex)
  for comp,idx in components.items():
   for i in range(3):
    for j in range(3):xp[idx[i],idx[j]]=1j*h[i,j]
  x=xp+xp.conj()[perm][:,perm];xgens.append(x);glabels.append(label+':'+hl)
Y=[np.array([[1+1j,2,0],[0,3-1j,1],[1j,0,4]],complex)+j*np.eye(3) for j in range(4)]
MR=np.array([[2+1j,1,0],[1,3,1j],[0,1j,4-1j]],complex);point=f.full_D(*Y,MR)
deltas=[x@point-point@x for x in xgens];Mnp=np.array([coords(d) for d in deltas]).T
checks['orbit_entries_integer_exact']=np.array_equal(Mnp,np.rint(Mnp))
checks['orbit_preserves_actual84_family']=all(np.array_equal(sum((c*b for c,b in zip(coords(d),basis)),np.zeros((96,96),complex)),d) for d in deltas)
M=s.Matrix(np.rint(Mnp).astype(np.int64).tolist());reduced,pivcols=M.rref(simplify=False);rank=len(pivcols)
free=[i for i in range(54) if i not in pivcols];null=[]
for k in free:
 v=s.zeros(54,1);v[k]=1
 for row,col in enumerate(pivcols):v[col]=-reduced[row,k]
 null.append(v)
checks['generic_orbit_rank53']=rank==53
checks['generic_stabilizer_one_dimension']=len(null)==1 and M*null[0]==s.zeros(84,1)
# Exhibit the exact nonzero minor, avoiding a numerical rank assertion.
_,rows=M[:,list(pivcols)].T.rref(simplify=False)
minor=M.extract(list(rows),list(pivcols));det=minor.det(method='domain-ge')
checks['exact53_by53_minor_nonzero']=det!=0 and minor.shape==(53,53)
checks['null_is_original_baryon_rotation']=len(null)==1 and all(((sum((float(c)*x for c,x in zip(null[0],xgens)),np.zeros((96,96),complex)))[i,i]==(0 if f.P[i%48][0] in ('nu','e') else (1j if i<48 else -1j))*complex((sum((float(c)*x for c,x in zip(null[0],xgens)),np.zeros((96,96),complex)))[f.IDX['u','L',0,0],f.IDX['u','L',0,0]]).imag) for i in range(96))
result=dict(contract='FUNDAMENTAL-SEARCH-OBJECTIVE-AUDIT-20261030',all_passed=all(checks.values()),checks={k:bool(v) for k,v in checks.items()},real_parameter_dimension=84,frame_dimension=54,generic_orbit_rank=rank,generic_stabilizer_dimension=len(null),generic_quotient_dimension=84-rank,ordered_family_block_types=len(unique),exact_minor=dict(rows=list(rows),columns=list(pivcols),determinant=str(det)),stabilizer_support=[(lab,str(c)) for lab,c in zip(glabels,null[0]) if c],elapsed_seconds=time.time()-start,search_verdict='NON_SELECTIVE_COMPATIBILITY_CLASS',final_solution_found=False,binary_subfamily_size=str(2**84),binary_subfamily_enumerated=False,coverage_method='Exact basis plus linearity, all real coefficient values',independent_source_selection_equation_supplied=False,scope='All real linear combinations in the ORIGINAL v252 finite-triple family, not all TFPT source/action/continuum constraints. Gaussian-integer finite algebra plus exact rational rank. No final physical solution selected.')
out=R/'experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review'
(out/'SEARCH-OBJECTIVE-AUDIT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
(out/'SEARCH-OBJECTIVE-ORBIT-MATRIX.json').write_text(json.dumps(dict(coordinates=coordinate_labels,generators=glabels,matrix=Mnp.astype(int).tolist()),indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2));assert result['all_passed']
