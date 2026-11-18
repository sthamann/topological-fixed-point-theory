"""Decisive native tensor test: quartic recovers the actual E7 adjoint channel.
Success: exact sparse projector and equality with all original E7 generators.
No physical state, action or contraction is selected by this calculation.
"""
import ast,itertools,json,math
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import numpy as np
import sympy as sp
from scipy.sparse import csr_matrix,eye
OUT=Path(__import__('os').environ.get('TFPT_E_SERIES_DATA_DIR','work/quartic-descent'));REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
qpath=OUT/'quartic-descent-results.json'
if not qpath.exists():qpath=REPO/'experiments/theory-contracts/e-series-compiler-closure-20261023/quartic-descent-results-20261106.json'
data=json.loads(qpath.read_text())
bootstrap=REPO/'experiments/theory-contracts/e-series-compiler-closure-20261023/e7_full_source_response_20261106.py'
t=ast.parse(bootstrap.read_text());nodes=[]
for n in t.body:
    if isinstance(n,ast.For):break
    if any(isinstance(z,ast.Name) and z.id in ['parser','OUT'] for z in ast.walk(n)):continue
    nodes.append(n)
env={};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(bootstrap),'exec'),env)
ch,ck,ip=(env[k] for k in ['ch','ck','ip'])
alpha=tuple(data['normal_root_doubled']);ap=ch.ridx[alpha]
G=[ch.ridx[tuple(r)] for r in data['native56_roots_doubled']];idx={r:i for i,r in enumerate(G)}
omega=np.zeros((56,56),dtype=np.int64)
for i,j in itertools.product(range(56),repeat=2):omega[i,j]=ch.bracket(G[i],G[j]).get(ap,0)
ck(np.array_equal(omega@omega,-np.eye(56,dtype=int)),'native omega is signed symplectic permutation')
oi=-omega
pairs=list(itertools.combinations_with_replacement(range(56),2));pi={p:i for i,p in enumerate(pairs)}
rows=[];cols=[];values=[];hist=Counter()
for term in data['quartic_monomials']:
    mon=term['indices'];perms=set(itertools.permutations(mon));v=F(term['coefficient'],len(perms))
    ck((2*v).denominator==1,'twice normalized quartic tensor entries integral');hist[str(v)]+=1
    for i,j,k,l in perms:
        if i>j:continue
        m=int(np.flatnonzero(oi[k])[0]);n=int(np.flatnonzero(oi[l])[0])
        rows.append(pi[(i,j)]);cols.append(pi[tuple(sorted((m,n)))]);values.append(int(2*v)*int(oi[k,m])*int(oi[l,n]))
K=csr_matrix((values,(rows,cols)),shape=(1596,1596),dtype=np.int64);K.sum_duplicates();K.eliminate_zeros()
K2=K@K
i,j=next((i,j) for i,j in zip(*K.nonzero()) if i!=j)
a=F(int(K2[i,j]),int(K[i,j]));b=F(int(K2[0,0]))-a*int(K[0,0])
ck(a.denominator==b.denominator==1,'quadratic tensor algebra has integral coefficients')
I=eye(1596,dtype=np.int64,format='csr')
ck((K2-int(a)*K-int(b)*I).nnz==0,'exact whole1596 tensor minimal polynomial')
z=sp.Symbol('z');eigs=sp.solve(z*z-int(a)*z-int(b),z)
trace=int(K.diagonal().sum());mult=sp.solve([sp.Symbol('m')+sp.Symbol('n')-1596,sp.Symbol('m')*eigs[0]+sp.Symbol('n')*eigs[1]-trace],['m','n'])
lam=next(e for e in eigs if F(trace-1596*next(f for f in eigs if f!=e),e-next(f for f in eigs if f!=e))==133)
other=next(e for e in eigs if e!=lam);den=int(lam-other);PN=K-int(other)*I
ck((PN@PN-den*PN).nnz==0 and F(int(PN.diagonal().sum()),den)==133,'exact rank133 projector from omega and quartic')

# Original E7 action, converted to a symmetric covariant tensor by omega.
E7=[r for r in range(240) if ip(alpha,ch.roots[r])==0]
cartan=sp.Matrix([[ch.roots[r][j]-F(alpha[j],2) for j in range(8)] for r in G]);cartan_basis=cartan.columnspace()
gens=[]
for r in E7:
    A=np.zeros((56,56),dtype=np.int64)
    for j,g in enumerate(G):
        for h,c in ch.bracket(r,g).items():A[idx[h],j]=int(c)
    gens.append(A)
for v in cartan_basis:gens.append(np.diag([int(x) for x in v]))
ck(len(gens)==133,'all126 original root actions and seven Cartans')
for A in gens:
    ck(np.array_equal(A.T@omega+omega@A,np.zeros((56,56),dtype=int)),'actual E7 action symplectic')
    B=omega@A;vec=np.array([B[i,j] for i,j in pairs],dtype=np.int64)
    ck(np.array_equal(K@vec,int(lam)*vec),'quartic projector includes actual E7 generator')

# Symplectic cubic gradient: no black-hole or physical-clock identification.
qterms=[(tuple(t['indices']),int(t['coefficient'])) for t in data['quartic_monomials']]
def qval(x):return sum(c*math.prod(x[i] for i in m) for m,c in qterms)
def sharp(x):
    grad=[F(0)]*56
    for m,c in qterms:
        for i,n in Counter(m).items():
            mm=list(m);mm.remove(i);grad[i]+=c*n*math.prod(x[j] for j in mm)
    return [sum(F(int(oi[i,j]))*grad[j] for j in range(56)) for i in range(56)]
witnesses=[]
for offset in range(4):
    x=[F((i*i+offset*i+2*offset)%5-2) for i in range(56)];qx=qval(x)
    if not qx:continue
    sx=sharp(x);ss=sharp(sx);rat={ss[i]/(qx*qx*x[i]) for i in range(56) if x[i]}
    ck(len(rat)==1 and all(not ss[i] for i in range(56) if not x[i]),'Freudenthal gradient composition whole56 point witness')
    witnesses.append({'q':str(qx),'composition_factor':str(next(iter(rat)))})
result={'claim':'MATH.ESERIES.QUARTIC.ADJ.PROJECTOR.20261106','tensor_definition':'Q(x)=T_ijkl x_i x_j x_k x_l; K_ij^mn=2 T_ijkl omegaInv_km omegaInv_ln on Sym2 covariant tensors; factor2 clears half-integral T entries','tensor_entry_histogram_by_monomial':dict(hist),'K_sparse_entries':K.nnz,'K_minimal_polynomial':f'K^2=({a})K+({b})I','K_eigenvalues':list(map(str,eigs)),'K_trace':trace,'adjoint_eigenvalue':str(lam),'other_eigenvalue':str(other),'projector':f'P133=(K-({other})I)/({den})','projector_rank':133,'actual_E7_generator_images_checked':len(gens),'Sym2_complement_dimension':1596-133,'sharp_definition':'omegaInv * gradient Q, without rescaling','sharp_composition_point_witnesses':witnesses,'scope':'Exact sparse full1596 projector identity and native E7 generator identification. Cubic-gradient composition here checked at four exact56 points, general FTS identity supplied by original-source construction and cited literature. No new physical measure, EM insertion, dynamical rule or clock identified.','checks':dict(env['checks'])}
(OUT/'quartic-projector-results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
