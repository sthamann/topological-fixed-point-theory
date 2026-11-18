"""Actual SM-preserving two-reflection transport between the two E7 windows."""
import ast,json,argparse,itertools
from collections import Counter,defaultdict
from pathlib import Path
from fractions import Fraction as F
REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
parser=argparse.ArgumentParser();parser.add_argument('--out',required=True,type=Path)
OUT=parser.parse_args().out;OUT.mkdir(parents=True,exist_ok=True)
bootstrap=REPO/'experiments/theory-contracts/e-series-compiler-closure-20261023/e7_full_source_response_20261106.py'
t=ast.parse(bootstrap.read_text());nodes=[]
for n in t.body:
    if isinstance(n,ast.For):break
    if any(isinstance(z,ast.Name) and z.id in ['parser','OUT'] for z in ast.walk(n)):continue
    nodes.append(n)
env={};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(bootstrap),'exec'),env)
ch,ck,sp=(env[k] for k in ['ch','ck','sp'])
I=sp.eye(8);alpha=sp.Matrix([0]*5+[1,-1,0]);aQ=sp.ones(8,1)/2
gamma=sp.Matrix([0]*6+[-1,-1]);beta=sp.Matrix([0]*5+[1,0,1])-aQ
M=(I-beta*beta.T)*(I-gamma*gamma.T)
ck(M*alpha==aQ and M.T*M==I,'literal two-reflection E7 normal transport')
Y=sp.Matrix([F(-1,3)]*3+[F(1,2)]*2+[0]*3)
sm=[]
for i in range(3):
    for j in range(3):
        if i!=j:r=sp.zeros(8,1);r[i]=1;r[j]=-1;sm.append(r)
r=sp.zeros(8,1);r[3]=1;r[4]=-1;sm.extend([r,-r])
ck(M*Y==Y and all(M*r==r for r in sm),'fixed original SM charges and nonabelian roots')
comm=[sp.Matrix(r)/2 for r in ch.roots if (sp.Matrix(r)/2).dot(Y)==0 and all((sp.Matrix(r)/2).dot(z)==0 for z in sm)]
ck(len(comm)==20 and sp.Matrix.hstack(*comm).rank()==4,'actual A4 root commutant of the marked SM')
def ridx(r):return ch.ridx[tuple(int(2*x) for x in r)]
def expad(root,v):
    out=dict(v);term=dict(v)
    for n in range(1,8):
        nxt={}
        for j,c in term.items():
            for k,cf in ch.bracket(root,j).items():nxt[k]=nxt.get(k,F(0))+c*cf/F(n)
        term={k:c for k,c in nxt.items() if c}
        if not term:return {k:c for k,c in out.items() if c}
        for k,c in term.items():out[k]=out.get(k,F(0))+c
    raise RuntimeError('nilpotent source exponential did not terminate')
def tits(r,v):
    i=ridx(r)
    return expad(i,expad(ch.opp[i],expad(i,v)))
def lift(v):return tits(beta,tits(gamma,v))
ck(all(lift({ridx(r):F(1)})=={ridx(r):F(1)} for r in sm),'actual original Tits lift fixes all8 SM root currents exactly')
yc=env['Ri']*Y;Yv={240+i:F(yc[i]) for i in range(8) if yc[i]}
ck(lift(Yv)==Yv,'actual source lift fixes hypercharge current exactly')
normal=lift({ridx(alpha):F(1)})
ck(set(normal)=={ridx(aQ)} and abs(normal[ridx(aQ)])==1,'normal root fields transported with actual cocycle sign')
for i,r in enumerate(ch.roots):
    image=M*sp.Matrix(r)/2
    v=lift({i:F(1)})
    ck(set(v)=={ridx(image)} and abs(next(iter(v.values())))==1,'all240 root field images match source Weyl map')
old=json.loads((OUT/'quartic-descent-results.json').read_text())
Gnew=[i for i,r in enumerate(ch.roots) if sum(r)==4]
new_index={i:k for k,i in enumerate(Gnew)}
new_alpha=ridx(aQ);V={ch.opp[new_alpha]:{():1}}
for degree in range(1,5):
    nxt=defaultdict(Counter)
    for a,root in enumerate(Gnew):
        for b,poly in V.items():
            for target,c in ch.bracket(root,b).items():
                for mon,v in poly.items():nxt[target][tuple(sorted(mon+(a,)))]+=int(c)*v
    V={b:{m:c for m,c in poly.items() if c} for b,poly in nxt.items() if any(poly.values())}
mapped=[];signs=[]
for r in old['native56_roots_doubled']:
    value=lift({ch.ridx[tuple(r)]:F(1)});j,c=next(iter(value.items()))
    mapped.append(new_index[j]);signs.append(int(c))
predicted={tuple(sorted(mapped[i] for i in mon['indices'])):int(mon['coefficient'])*__import__('math').prod(signs[i] for i in mon['indices']) for mon in old['quartic_monomials']}
ck(V=={new_alpha:predicted},'all1036 quartic coefficients transported by the same SM-preserving source lift')
rho=sp.Matrix([F(7-2*i,2) for i in range(8)]);q=sp.ones(8,1)
changes=[]
for a in [F(1,4),F(3,4)]:
    w=rho-F(2*a-1,2)*q
    ck(M*w!=w,'SM-preserving transport is not a fixed-H8 symmetry')
    changes.append({'a':str(a),'w_fixed':False,'M_w':list(map(str,M*w))})
result={'claim':'MATH.ESERIES.SM.PRESERVING.E7.TRANSPORT.20261106','native_E7_normal':list(map(str,alpha)),'Q_neutral_E7_normal':list(map(str,aQ)),'reflection_roots':[list(map(str,gamma)),list(map(str,beta))],'source_Tits_normal_sign':str(normal[ridx(aQ)]),'all240_root_fields_transported':True,'all8_SM_root_currents_fixed':True,'hypercharge_current_fixed':True,'transformed_time_markings':changes,'scope':'Exact two-reflection/Tits transport inside the original E8 source. It fixes the SM current algebra, but transports the source clock marking; not a symmetry of the same fixed H8 or identification of56 with48.','checks':dict(env['checks'])}
result['SM_commutant_root_count_rank_type']=[20,4,'A4']
result['quartic_coefficients_transported']=1036
(OUT/'sm-preserving-transport-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
