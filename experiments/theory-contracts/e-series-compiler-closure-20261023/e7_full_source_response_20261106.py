"""Same original E8 affine source: E7+A1 Ward response and full glued trace."""
import ast
import argparse
import hashlib
import itertools
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import mpmath as mp
import numpy as np
import sympy as sp

REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
parser=argparse.ArgumentParser()
parser.add_argument('--out',required=True,type=Path)
OUT=parser.parse_args().out
OUT.mkdir(exist_ok=True,parents=True)
PINS={}
checks=Counter()
def ck(x,name):
    if not bool(x): raise RuntimeError(name)
    checks[name]+=1
class Always(ast.NodeTransformer):
    def visit_Assert(self,n):
        return ast.copy_location(ast.Expr(ast.Call(ast.Name('ck',ast.Load()),[n.test,ast.Constant('original assertion')],[])),n)
def original(path,name,env):
    raw=(REPO/path).read_bytes()
    PINS[path]=hashlib.sha256(raw).hexdigest()
    tree=ast.parse(raw)
    node=next(n for n in tree.body if isinstance(n,(ast.ClassDef,ast.FunctionDef)) and n.name==name)
    env['ck']=ck
    exec(compile(ast.fix_missing_locations(Always().visit(ast.Module(body=[node],type_ignores=[]))),str(REPO/path),'exec'),env)
    return env[name]
def ip(a,b):return sum(x*y for x,y in zip(a,b))
def add(a,b):return tuple(x+y for x,y in zip(a,b))
raw=original('verification/v1_e8_glue.py','e8_roots',{'np':np,'itertools':itertools})()
roots=[tuple(int(2*x) for x in r) for r in raw]
triality_path='experiments/theory-contracts/compiler-integral-triality-20260918/certificate.json'
triality_raw=(REPO/triality_path).read_bytes()
PINS[triality_path]=hashlib.sha256(triality_raw).hexdigest()
triality=json.loads(triality_raw)
R=sp.Matrix(triality['integral_triality']['vector_basis']);Ri=R.inv()
coords={r:tuple(F(x) for x in Ri*sp.Matrix(r)/2) for r in roots}
class Adapter:
    n=dim=8
    def __init__(self):
        self.roots=roots
        self.simple=[tuple(int(x) for x in 2*R[:,i]) for i in range(8)]
    def alpha_coords(self,r):return coords[r]
src='verification/v498_celestial_wp5b_singular_vector.py'
ch=original(src,'Chevalley',{'F':F,'ip':ip,'vadd':add})(Adapter())
af=original(src,'Affine',{'F':F})(ch,1)
vac={():F(1)}
ck(ch.sgn==ch.kappa_root==-1,'original compact convention')
def star(i):return (ch.opp[i],-1) if i<240 else (i,1)
for i in range(240):
    for n in [1,2,3]:
        j,c=star(i)
        value=c*af.act(j,n,af.act(i,-n,vac)).get((),0)
        ck(value==n,'inherited level-one root response')

phys=[sp.Matrix(r)/2 for r in ch.roots]
rho=sp.Matrix([F(7-2*j,2) for j in range(8)])
q=sp.ones(8,1)
Y=sp.Matrix([F(-1,3)]*3+[F(1,2)]*2+[0]*3)
yc=Ri*Y
Ystate={}
for i in range(8):
    if yc[i]:Ystate[((-1,240+i),)]=F(yc[i])
Ygram=F(0)
for i in range(8):
    if yc[i]:Ygram+=F(yc[i])*af.act(240+i,1,Ystate).get((),F(0))
ck(Ygram==F(5,6),'original hypercharge current Gram')
normals=[]
for i,j in itertools.combinations([5,6,7],2):
    for sign in [-1,1]:
        a=sp.zeros(8,1);a[i]=1;a[j]=sign
        normals.append(('D5_'+str(i)+'_'+str(j)+'_'+str(sign),a))
normals.append(('Q_neutral',q/2))
decomp=[]
for name,a in normals:
    e7=[i for i,r in enumerate(phys) if ip(a,r)==0]
    a1=[i for i,r in enumerate(phys) if r==a or r==-a]
    mixed=[i for i in range(240) if i not in set(e7+a1)]
    ck((len(e7),len(a1),len(mixed))==(126,2,112),'full E7 A1 glue branching')
    P1=a*a.T/2;P7=sp.eye(8)-P1
    ck(sum((r*r.T for r in [phys[i] for i in e7]),sp.zeros(8))==36*P7,'E7 adjoint trace tensor')
    ck(sum((r*r.T for r in [phys[i] for i in a1]),sp.zeros(8))==4*P1,'A1 adjoint trace tensor')
    ck(sum((r*r.T for r in [phys[i] for i in mixed]),sp.zeros(8))==24*P7+56*P1,'mixed112 trace tensor')
    ck(sum((r*r.T for r in phys),sp.zeros(8))==60*sp.eye(8),'whole E8 invariant response')
    norms=Counter((F(ip(P7*phys[i],P7*phys[i])),F(ip(P1*phys[i],P1*phys[i]))) for i in mixed)
    ck(norms=={(F(3,2),F(1,2)):112},'glue weights 3/4 plus 1/4')
    ck(P7*Y==Y and P1*Y==sp.zeros(8,1),'actual hypercharge retained in E7')
    ytr=[sum(F(ip(Y,phys[i]))**2 for i in part) for part in [e7,a1,mixed]]
    ck(ytr==[F(30),F(0),F(20)],'full hypercharge adjoint trace30 plus20')
    decomp.append({'cut':name,'root_counts':[126,2,112],'E7_current_level':1,'A1_current_level':1,'mixed_projection_squared_norms':['3/2','1/2'],'hypercharge_current_Gram':'5/6','hypercharge_adjoint_trace_E7_A1_glue':list(map(str,ytr))})
print('Original affine levels and seven full branchings checked.',flush=True)

def addstate(target,state,c=F(1)):
    for key,v in state.items():target[key]=target.get(key,F(0))+c*v
    return target
def sugawara(a=None,which='all'):
    P=sp.eye(8) if a is None else (a*a.T/2 if which=='a1' else sp.eye(8)-a*a.T/2)
    inds=list(range(240)) if a is None else [i for i,r in enumerate(phys) if (ip(a,r)!=0 if which=='a1' else ip(a,r)==0)]
    if which=='a1':inds=[i for i,r in enumerate(phys) if r==a or r==-a]
    C=Ri*P*Ri.T
    ans={}
    for i in inds:addstate(ans,af.act(i,-1,af.act(ch.opp[i],-1,vac)),F(-1))
    for i,j in itertools.product(range(8),repeat=2):
        if C[i,j]:addstate(ans,af.act(240+i,-1,af.act(240+j,-1,vac)),F(C[i,j]))
    h={'all':30,'e7':18,'a1':2}[which]
    return {k:v/F(2*(1+h)) for k,v in ans.items() if v}
@lru_cache(None)
def wordgram(left,right):
    state={right:F(1)};c=F(1)
    for m,i in left:
        j,ci=star(i);c*=ci;state=af.act(j,-m,state)
    ck(set(state)<={()},'grade-two PBW vacuum contraction')
    return c*state.get((),F(0))
def inner(left,right):return sum(c*d*wordgram(k,l) for k,c in left.items() for l,d in right.items())
T8=sugawara()
stress=[]
# One native D5-compatible cut and the source-charge-selected cut are checked
# at the actual operator level. All seven invariant trace tensors above agree.
for name,a in [normals[0],normals[-1]]:
    T7=sugawara(a,'e7');T1=sugawara(a,'a1')
    D=dict(T8);addstate(D,T7,F(-1));addstate(D,T1,F(-1));D={k:v for k,v in D.items() if v}
    values=[inner(T7,T7),inner(T1,T1),inner(T7,T1),inner(T8,T8),inner(D,D)]
    ck(values==[F(7,2),F(1,2),F(0),F(4),F(0)],'actual Sugawara Gram and null coset')
    stress.append({'cut':name,'Gram_E7_A1':['7/2','1/2','0'],'Gram_total':'4','norm_T8_minus_T7_minus_T1':'0'})
    print('Actual PBW stress response:',name,list(map(str,values)),flush=True)

mp.mp.dps=85
M=18
def ms(x):return mp.mpf(str(x))
def one(t,w,half=False,signed=False,scale=4):
    return mp.fsum(((-1)**n if signed else 1)*mp.exp(-scale*t*(mp.mpf(n)+(mp.mpf('.5') if half else 0))**2+t*w*(mp.mpf(n)+(mp.mpf('.5') if half else 0))) for n in range(-M,M+1))
def theta8(t,w):
    return (mp.fprod(one(t,x) for x in w)+mp.fprod(one(t,x,signed=True) for x in w)+mp.fprod(one(t,x,True) for x in w)+mp.fprod(one(t,x,True,True) for x in w))/2
def theta_coordinate_cut(t,w,i,j,sign):
    # Turn alpha=e_i+e_j into a difference normal by flipping the j variable.
    # An odd sign flip changes the half-integer parity; retain it explicitly.
    ww=list(w);parity=1
    if sign==1:ww[j]=-ww[j];parity=-1
    rest=[k for k in range(8) if k not in [i,j]]
    S=mp.fprod(one(t,ww[k]) for k in rest);A=mp.fprod(one(t,ww[k],signed=True) for k in rest)
    H=mp.fprod(one(t,ww[k],True) for k in rest);B=parity*mp.fprod(one(t,ww[k],True,True) for k in rest)
    Ui=one(t,ww[i]+ww[j],scale=8);Uh=one(t,ww[i]+ww[j],True,scale=8)
    E0=(S+A)/2*Ui+(H+B)/2*Uh
    E56=(S-A)/2*Uh+(H-B)/2*Ui
    A0=one(t,ww[i]-ww[j],scale=8);A2=one(t,ww[i]-ww[j],True,scale=8)
    return E0,A0,E56,A2
def coefficient(t,rho,sum_value,half=False):
    # Coefficient extraction of the actual eight-coordinate lattice product.
    # Keys track twice the sum of coordinates, so all arithmetic labels are exact.
    acc={0:mp.mpf(1)}
    for w in rho:
        nxt={}
        for n in range(-M,M+1):
            v=mp.mpf(n)+(mp.mpf('.5') if half else 0);key=2*n+(1 if half else 0)
            val=mp.exp(-4*t*v*v+t*w*v)
            for k,c in acc.items():nxt[k+key]=nxt.get(k+key,mp.mpf(0))+c*val
        acc=nxt
    return acc[2*sum_value]
def theta_neutral(t,w):
    rr=[ms(x) for x in rho]
    E0=coefficient(t,rr,0)+coefficient(t,rr,0,True)
    E56=mp.exp(2*t)*(coefficient(t,rr,2)+coefficient(t,rr,2,True))
    wa=mp.fsum(w)/2
    A0=one(t,wa,scale=8);A2=one(t,wa,True,scale=8)
    return E0,A0,E56,A2
def trace(t,w):
    osc=mp.fprod(1-mp.exp(-8*t*n) for n in range(1,401))
    return theta8(t,w)/osc**8
samples=[]
for a in [mp.mpf('.25'),mp.mpf('.75')]:
    w=[ms(x)-(2*a-1)/2 for x in rho]
    for t in [mp.mpf('.25'),mp.mpf('.5'),mp.mpf(1),mp.mpf(2)]:
        th=theta8(t,w);z=trace(t,w)
        row={'a':str(a),'t':str(t),'Z_full':mp.nstr(z,35),'cuts':[]}
        for name,normal in normals:
            if name=='Q_neutral':E0,A0,E56,A2=theta_neutral(t,w)
            else:
                _,si,sj,sgn=name.split('_')
                E0,A0,E56,A2=theta_coordinate_cut(t,w,int(si),int(sj),int(sgn))
            val=E0*A0+E56*A2
            err=abs(val/th-1)
            ck(err<mp.mpf('1e-65'),'full glued original heat trace identity')
            row['cuts'].append({'name':name,'relative_error':mp.nstr(err,5),'glue_fraction':mp.nstr(E56*A2/th,25),'E7_A1_product_fraction':mp.nstr(E0*A0/th,25)})
        samples.append(row)
        print('Full original heat trace:',str(a),str(t),mp.nstr(z,18),'neutral glue',row['cuts'][-1]['glue_fraction'],flush=True)

# Exact source-time response, retaining all positive affine modes by its
# geometric series. This is a vacuum current kernel, separate from heat traces.
time_response=[]
for a in [F(1,4),F(3,4)]:
    w=rho-(2*a-1)*q/2
    row={'a':str(a),'cuts':[]}
    for name,normal in normals:
        e7=[r for r in phys if ip(normal,r)==0]
        a1=[r for r in phys if r==normal or r==-normal]
        mixed=[r for r in phys if r not in e7 and r not in a1]
        hist=[Counter(F(8-ip(w,r)) for r in part) for part in [e7,a1,mixed]]
        row['cuts'].append({'name':name,'grade_one_root_energy_histograms':[{str(k):v for k,v in sorted(h.items())} for h in hist]})
    time_response.append(row)

result={'claim':'MATH.SEAM.E7.FULL.RESPONSE.20261106','scope':'The given unitary level-one E8 lattice VOA, its existing vacuum/current normalization and the original H8 generator. No identification with the fully integrated primitive physical S_adm or a new measured coupling.','checks':dict(checks),'branchings':decomp,'PBW_stress':stress,'heat_trace_samples':samples,'source_time_current_response':time_response,'exact_formulas':{'module':'V8=(V7^0 tensor V1^0) direct_sum (V7^56 tensor V1^2)','trace':'Z8=Z7^0 Z1^0+Z7^56 Z1^2','H8':'8 L0-rho_A7.a0+(2a-1)Q/2','vacuum_current_kernel':'K8(t)=exp(-8t)/(1-exp(-8t))^2 * [8+sum_roots exp(t*(rho.r-(2a-1)Q(r)/2))]','local_stress_Gram':'[[7/2,0],[0,1/2]], total4','current_levels':'k_E7=k_A1=k_E8=1','inherited_seam_normalization':'Full coordinate rewrite changes existing c3 by factor1; absolute P1 origin not newly derived.'}}
result['source_sha256']=PINS
(OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS',sum(checks.values()),'exact/numerical assertions; no new physical constant derived.',flush=True)
