"""Neutral ordered four-current response realizing the native E7 quartic."""
import ast,json,itertools,argparse
from pathlib import Path
from collections import Counter,defaultdict
from fractions import Fraction as F
REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True)
OUT=parser.parse_args().out;OUT.mkdir(parents=True,exist_ok=True)
data=json.loads((OUT/'quartic-descent-results.json').read_text())
bootstrap=REPO/'experiments/theory-contracts/e-series-compiler-closure-20261023/e7_full_source_response_20261106.py'
tree=ast.parse(bootstrap.read_text());nodes=[]
for n in tree.body:
    if isinstance(n,ast.For):break
    if any(isinstance(z,ast.Name) and z.id in ['parser','OUT'] for z in ast.walk(n)):continue
    nodes.append(n)
env={};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(bootstrap),'exec'),env)
ch,af,ck=(env[k] for k in ['ch','af','ck']);vac={():F(1)}
alpha=tuple(data['normal_root_doubled']);ap=ch.ridx[alpha];am=ch.opp[ap]
G=[ch.ridx[tuple(r)] for r in data['native56_roots_doubled']]
minus=[]
for i in G:
    b=ch.bracket(am,i);ck(len(b)==1,'actual A1 lowering intertwiner')
    j,c=next(iter(b.items()));ck(abs(c)==1,'same-source lowering normalization');minus.append((j,c))

# Y(x)=[Xplus(x),Xminus(x)] and exact coefficient-wise kappa(Y,Y).
Y=defaultdict(Counter)
for i,j in itertools.product(range(56),repeat=2):
    r,c=minus[j]
    for k,v in ch.bracket(G[i],r).items():Y[k][tuple(sorted((i,j)))]+=c*v
Y={k:{m:c for m,c in p.items() if c} for k,p in Y.items() if any(p.values())}
response=Counter()
for i,pi in Y.items():
    for j,pj in Y.items():
        g=ch.kappa(i,j)
        if not g:continue
        for mi,ci in pi.items():
            for mj,cj in pj.items():response[tuple(sorted(mi+mj))]+=g*ci*cj
response={m:c for m,c in response.items() if c}
Q={tuple(m['indices']):F(m['coefficient']) for m in data['quartic_monomials']}
ck(response=={m:-c for m,c in Q.items()},'all1036 neutral four-current coefficients equal minus Q')
def linact(terms,n,state):
    out={}
    for i,c in terms:
        if not c:continue
        for k,v in af.act(i,n,state).items():out[k]=out.get(k,F(0))+c*v
    return {k:v for k,v in out.items() if v}
tests=[]
points=[]
chosen=next(m for m,c in Q.items() if len(set(m))==4 and abs(c)==24)
points.append([int(i in chosen) for i in range(56)])
points.append([(i%3)-1 for i in range(56)])
for x in points:
    plus=[(G[i],F(x[i])) for i in range(56)]
    neg=[(j,c*x[i]) for i,(j,c) in enumerate(minus)]
    state=linact(neg,-1,vac)
    state=linact(plus,0,state)
    state=linact(neg,0,state)
    state=linact(plus,1,state)
    value=sum(c*F(__import__('math').prod(x[i] for i in m)) for m,c in response.items())
    ck(state==({():value} if value else {}),'actual original affine four-current vacuum evaluation')
    tests.append({'x_pattern':'four-support' if x is points[0] else 'ternary-full56','affine_vacuum_response':str(value)})
result={'claim':'MATH.ESERIES.NEUTRAL.QUARTIC.RESPONSE.20261106','definition':'Xminus=[e_-alpha,Xplus]; R4(x)=<Omega|J_Xplus(1) J_Xminus(0) J_Xplus(0) J_Xminus(-1)|Omega>','exact_identity':'R4(x)=kappa([Xplus,Xminus],[Xplus,Xminus])=-Q(x)','matched_quartic_monomials':len(Q),'A1_grade_sequence_from_right':[-1,0,-1,0],'net_A1_charge':0,'actual_affine_witnesses':tests,'E6_derivative_identity':'d/dt R4(t s+x)|t=0=-24 P(x), with P=d(x,x,x)/6; normalized symmetric quartic polarization is minus d(suppressed singlet,x,x,x) in the stated root conventions.','scope':'Exact neutral ordered source correlator, using original affine vacuum and current algebra. Not a fully integrated primitive physical response or a determination of the electromagnetic delta_top coefficient48.','checks':dict(env['checks'])}
(OUT/'neutral-quartic-response-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
