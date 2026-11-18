"""Exact ordered return of original mixed E8 source response on native code.
Sparse Gaussian-rational Heisenberg Fock arithmetic; no intermediate projection.
"""
from pathlib import Path
import sys,json,math,itertools as it
from collections import Counter
from functools import lru_cache
import sympy as s
from sympy.polys.domains import QQ_I as Q
R=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4');sys.path.insert(0,str(R))
BASE=R/'experiments/theory-contracts/joint-native-algebra-selection-20261016'
from tfpt_explorer.cartan_source import _lattice_data,_TARGET_SIMPLE
from tfpt_explorer.current_block_geometry import _actual_current_source
from tfpt_explorer.neutral_source_response import _native_neutral_branch
from tfpt_explorer.process import _code_basis

def need(ok,msg):
    if not ok:raise ArithmeticError(msg)
def clean(A):return A.applyfunc(s.simplify)
def zero(A):return clean(A)==s.zeros(*A.shape)
def mat(A):return [[str(v) for v in row] for row in A.tolist()]
def conj(c):return Q(c.x,-c.y)
def add(*ps):
    out={}
    for a,p in ps:
        for m,c in p.items():out[m]=out.get(m,Q.zero)+a*c
    return {m:c for m,c in out.items() if c}
def mul(a,b):
    out={}
    for m,c in a.items():
        for n,d in b.items():
            w=tuple(sorted(m+n));out[w]=out.get(w,Q.zero)+c*d
    return {m:c for m,c in out.items() if c}
def deriv(f,b):
    out={}
    for m,c in f.items():
        for j,k in enumerate(m):
            need(k<4,'Derivative acts only on original positive code variables')
            n=m[:j]+m[j+1:];out[n]=out.get(n,Q.zero)+c*b[k]
    return {m:c for m,c in out.items() if c}
@lru_cache(None)
def weight(m):
    return math.prod(math.factorial(n) for n in Counter(m).values())*math.prod(1+k//8 for k in m)
def inner(a,b):return sum((conj(c)*b.get(m,Q.zero)*weight(m) for m,c in a.items()),Q.zero)
def linear(b,mode):
    out={}
    for j,c in enumerate(b):
        if c:out[((mode-1)*8+j,)]=conj(c);out[((mode-1)*8+4+j,)]=c
    return out
def code():
    patterns=[{(i,i,i,i):Q.one for i in range(4)},
              {(0,0,1,1):Q.one,(2,2,3,3):Q.one},
              {(0,0,2,2):Q.one,(1,1,3,3):Q.one},
              {(0,0,3,3):Q.one,(1,1,2,2):Q.one},
              {(0,1,2,3):Q.one}]
    norms=[96,8,8,8,1];_,supports=_code_basis()
    for j,sp in enumerate(supports):
        actual=Counter(tuple(sorted(w)) for w in sp)
        need(set(actual)==set(patterns[j]),'Actual original code monomials')
        cs=set(actual.values());need(len(cs)==1,'Constant support multiplicity')
        c=next(iter(cs));need(s.Rational(c*c,24*len(sp))==s.Rational(1,norms[j]),'Original tensor/Fock normalization')
        need(inner(patterns[j],patterns[j])==Q(norms[j]),'Fock norm of actual unnormalized code polynomial')
    return patterns,norms
def response(f,b,c):
    degree=max(map(len,f))
    need(degree in (3,4) and all(len(m)==degree for m in f),'Actual homogeneous quartic/trimer source input')
    B=[{},linear(b,1),linear(b,2),linear(b,3)]
    V=[{},linear(c,1),linear(c,2),linear(c,3)]
    E=[{():Q.one},B[1],add((Q(1)/2,mul(B[1],B[1])),(Q(1)/2,B[2])),
       add((Q(1)/6,mul(mul(B[1],B[1]),B[1])),(Q(1)/2,mul(B[1],B[2])),(Q(1)/3,B[3]))]
    F=[f]
    for k in range(1,degree+1):F.append(add((-Q.one/k,deriv(F[-1],b))))
    G=[deriv(f,c)]
    for k in range(1,degree):G.append(add((-Q.one/k,deriv(G[-1],b))))
    p=add(*[(Q.one,mul(V[m],mul(F[k],E[k-m-1]))) for m in range(1,degree) for k in range(m+1,degree+1)],
          *[(Q.one,mul(G[k],E[k])) for k in range(degree)])
    need(all(sum(1+k//8 for k in m)==degree-1 for m in p),'Actual charged oscillator grade and total conformal grade')
    # Independent previously derived partition-(3) coefficient, for all8
    # orthonormal positive/negative one-particle directions.
    part3=add((Q.one,mul(V[degree-1],F[degree])),(Q.one/(degree-1),mul(G[degree-1],B[degree-1])))
    need({m:a for m,a in p.items() if len(m)==1}==part3,'Exact first-response PBW partition3 agrees')
    return p

def main():
    L=_lattice_data();S=L['basis'];T=s.Matrix.hstack(*map(s.Matrix,_TARGET_SIMPLE));O=clean(T*S.inv()/s.sqrt(2))
    need(zero(O.T*O-s.eye(8)),'Original displayed orthogonal isometry')
    J=s.zeros(8);U=s.zeros(8,4)
    for j in range(4):J[2*j,2*j+1]=-1;J[2*j+1,2*j]=1;U[2*j,j]=1/s.sqrt(2);U[2*j+1,j]=-s.I/s.sqrt(2)
    need(zero(O*J-J*O),'Original map preserves native J')
    source=_actual_current_source();roots=source['chevalley'].roots
    need({tuple(T*s.Matrix(c)) for c in L['root_coordinates']}==set(roots),'All240 original source roots match, no fitted map')
    AR=_native_neutral_branch()['R_matrix'];AN=clean(O.T*AR*O)
    polys,norms=code();D=s.diag(*[s.sqrt(n) for n in norms]);rows=[];sum_raw=s.zeros(5);rank_counts=Counter();scalar_count=0;commP=0
    old=json.loads((BASE/'polar_fibre_results.json').read_text());e0=s.Matrix(old['fixed_basis']['e0']);e1=s.Matrix(old['fixed_basis']['e1']);a=s.Matrix(old['fixed_basis']['a']);cp=clean(-s.sqrt(s.Rational(3,5))*e0+s.sqrt(s.Rational(2,5))*e1);Fix=clean(a*a.H+e0*e0.H+e1*e1.H);P=clean(Fix-cp*cp.H);tau=s.eye(5)-2*a*a.H
    raw=s.Matrix([[2,2,-1,-1,-1,-1],[0,0,-1,-1,1,1],[0,0,-1,1,-1,1],[0,0,1,-1,-1,1],[1,-1,0,0,0,0]])
    N=clean(s.diag(*[s.sqrt(v) for v in (4,12,12,12,24)])*raw/s.sqrt(48));perm=s.zeros(6)
    for j,k in enumerate((0,1,3,4,2,5)):perm[k,j]=1
    sigma=clean(N*perm*N.H);Praw=clean(D.inv()*P*D);sigraw=clean(D.inv()*sigma*D);tauraw=clean(D.inv()*tau*D)
    for ix,r in enumerate(roots):
        beta=clean(O.T*s.Matrix(r)/2);v=clean(-AN*beta+((beta.T*AN*beta)[0]/2)*beta)
        b=[Q.from_sympy(s.simplify(t)) for t in U.T*beta];c=[Q.from_sympy(s.simplify(t)) for t in U.T*v]
        W=[response(f,b,c) for f in polys]
        G=s.Matrix(5,5,lambda i,j:Q.to_sympy(inner(W[i],W[j])));need(G==G.H,'Actual positive-adjoint return Gram')
        Kraw=s.diag(*[s.Rational(1,n) for n in norms])*G
        rk=G.rank();rank_counts[rk]+=1;scalar=zero(Kraw-s.trace(Kraw)*s.eye(5)/5);scalar_count+=scalar
        compatible=zero(Kraw*Praw-Praw*Kraw);commP+=compatible
        rows.append({'root_doubled':list(r),'response_norm_squared':str((v.T*v)[0]),'rank':rk,'scalar':scalar,'commutes_original_P2':compatible,'raw_gram':mat(G)})
        sum_raw+=G
        if ix%60==59:print('completed',ix+1,flush=True)
    K=clean(D.inv()*sum_raw*D.inv());x=s.Symbol('x');char=s.factor(K.charpoly(x).as_expr());print('SUM',mat(K),flush=True);print('CHAR',char,flush=True)
    out={'status':'PASS_EXACT_ORIGINAL_MIXED_ORDERED_RETURN','original_isometry':mat(O),'original_R_native_chart':mat(AN),'root_count':len(rows),'nonzero_response_count':sum(v['rank']>0 for v in rows),'individual_rank_counts':dict(rank_counts),'individual_scalar_count':scalar_count,'individual_commutes_original_P2':commP,'response_basis':'Derivatives indexed by original unit-norm root-current basis; W_beta are NOT independently renormalized. Sum is an original-R-weighted response Gram, not a physical event-channel choice.','ordered_return_sum':mat(K),'sum_characteristic_polynomial':str(char),'sum_commutes_P2':zero(K*P-P*K),'sum_commutes_sigma':zero(K*sigma-sigma*K),'sum_commutes_tau':zero(K*tau-tau*K),'sum_traceless_norm_squared':str(s.simplify(s.trace((K-s.trace(K)*s.eye(5)/5)**2))),'scope':'Actual source mode-insertion return matrix elements on code, without an intermediate projection. No physical coupling, detector/reset rule, selected history or polar/time selection derived.','individual_return_grams':rows}
    Path(__file__).with_name('mixed_return_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='individual_return_grams'},indent=2))
if __name__=='__main__':main()
