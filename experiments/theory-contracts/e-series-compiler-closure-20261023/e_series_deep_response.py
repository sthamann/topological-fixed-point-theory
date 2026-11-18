"""Full Cartan polynomial discriminators and actual source contraction controls."""
import ast,itertools,json,math
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
import sympy as sp
REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4');OUT=Path(__import__('os').environ.get('TFPT_E_SERIES_DATA_DIR','work/quartic-descent'))
boot=REPO/'experiments/theory-contracts/e-series-compiler-closure-20261023/e7_full_source_response_20261106.py'
t=ast.parse(boot.read_text());nodes=[]
for n in t.body:
    if isinstance(n,ast.For):break
    if any(isinstance(z,ast.Name) and z.id in ['parser','OUT'] for z in ast.walk(n)):continue
    nodes.append(n)
env={};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(boot),'exec'),env)
ch,af,ck,ip=(env[k] for k in ['ch','af','ck','ip']);vac={():F(1)}
alpha=(0,0,0,0,0,2,-2,0);G=[i for i,r in enumerate(ch.roots) if ip(alpha,r)==4]
native48=[i for i,r in enumerate(ch.roots) if all(abs(x)==1 for x in r) and sum(x<0 for x in r[:5])%2==1 and r[5:]!=(-1,-1,-1)]
x=sp.symbols('x0:5');s2=sum(z*z for z in x);s4=sum(z**4 for z in x)
def moment(space,n):return sp.expand(sum(sum(sp.Rational(ch.roots[r][i],2)*x[i] for i in range(5))**n for r in space))
expected={(56,2):12*s2,(48,2):12*s2,(56,4):6*s2**2,(48,4):9*s2**2-6*s4,(56,5):sp.Integer(0),(48,5):-180*sp.prod(x)}
for (dim,n),v in expected.items():ck(sp.expand(moment(G if dim==56 else native48,n)-v)==0,'full D5 Cartan polynomial trace identity')
# Independent original E8 adjoint weight identities, all8 Cartan variables.
y=sp.symbols('y0:8');r2=sum(z*z for z in y)
hy=list(y);hy[6]=hy[5];norm7=sum(z*z for z in hy)
for n,v in [(2,12*norm7),(4,6*norm7**2)]:
    actual=sp.Poly(sum(sum(sp.Rational(ch.roots[r][i],2)*hy[i] for i in range(8))**n for r in G),*y)
    ck(actual==sp.Poly(v,*y),'full seven-dimensional E7 Cartan trace identity')
for n,expected8 in [(2,60*r2),(4,36*r2**2),(6,30*r2**3)]:
    actual=sp.Poly(sum(sum(sp.Rational(r[i],2)*y[i] for i in range(8))**n for r in ch.roots),*y)
    ck(actual==sp.Poly(expected8,*y),'full E8 Cartan trace identity')

def linear(terms,n,state):
    result={}
    for i,c in terms:
        if c:
            for word,v in af.act(i,n,state).items():result[word]=result.get(word,F(0))+c*v
    return {w:c for w,c in result.items() if c}
Y=sp.Matrix([F(-1,3)]*3+[F(1,2)]*2+[0]*3);yc=env['Ri']*Y
yt=[(240+i,F(yc[i])) for i in range(8)]
def C(state):
    a=linear(yt,1,state);b=linear(yt,-1,state)
    for w,c in b.items():a[w]=a.get(w,F(0))+c
    return {w:c for w,c in a.items() if c}
state=vac;mom={}
for n in range(1,5):state=C(state);mom[n]=state.get((),F(0))
ck(mom[2]==F(5,6) and mom[4]==F(25,12),'actual Y affine vacuum moments')
connected=mom[4]-3*mom[2]**2;ck(connected==0,'actual Y fourth connected oscillator cumulant zero')
ap=ch.ridx[alpha];am=ch.opp[ap];r=G[0];b=ch.bracket(am,r)
bstar={ch.opp[r]:F(-1)}
def word(B):
    state=linear(list(B.items()),-1,vac);state=linear([(r,F(1))],0,state)
    state=linear(list(B.items()),0,state);return linear([(r,F(1))],1,state).get((),F(0))
lower=word(b);compact=word(bstar)
ck(lower==0 and compact==2,'A1 lowering and actual compact adjoint contractions differ')
result={'claim':'MATH.ESERIES.HIGHER.RESPONSE.DISCRIMINATORS.20261106','full_D5_Cartan_trace_identities':{f'Tr{dim}_H{n}':str(sp.factor(v)) for (dim,n),v in expected.items()},'fourth_trace_difference_56_minus48':'3*(2*sum_i x_i^4-(sum_i x_i^2)^2)','fifth_moment_chirality':'Tr48 H^5=-180*x0*x1*x2*x3*x4, whereas Tr56 H^5=0','E7_56_quartic_trace_factorization':'Tr56 H^4=(Tr56 H^2)^2/24, checked as full polynomial on all seven native E7 Cartan directions','E8_adjoint_trace_identities':{'Tr248_H2':'60*sum_i y_i^2','Tr248_H4':'36*(sum_i y_i^2)^2=(Tr248_H2)^2/100','Tr248_H6':'30*(sum_i y_i^2)^3=(Tr248_H2)^3/7200'},'actual_Y_oscillator_vacuum':{'moment2':str(mom[2]),'moment4':str(mom[4]),'connected4':str(connected)},'single_root_contraction_control':{'root_doubled':ch.roots[r],'A1_lowering_R4':str(lower),'compact_adjoint_R4':str(compact)},'scope':'Internal trace polynomials, original affine oscillator and ordered source words. Holomorphic quartic invariant on56 is a different type from a fourth power trace of a Lie algebra element and from a connected physical EM response. This does not exclude interacting or full integrated EM responses.','checks':dict(env['checks'])}
(OUT/'deep-response-results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
