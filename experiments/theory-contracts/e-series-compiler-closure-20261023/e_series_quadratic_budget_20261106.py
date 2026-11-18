"""Exact native E7 56 versus existing three-family48 response and orbit flag."""
import ast,itertools,json,argparse
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import sympy as sp
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
ch,ck,ip=(env[k] for k in ['ch','ck','ip'])
R=[sp.Matrix(r)/2 for r in ch.roots]
alpha=sp.Matrix([0]*5+[1,-1,0])
beta=sp.Matrix([0]*5+[1,0,-1])
lam=sp.Matrix([0]*5+[1,0,1])
G=[r for r in R if r.dot(alpha)==1]
carrier=[r for r in R if all(abs(x)==sp.Rational(1,2) for x in r) and sum(bool(x<0) for x in r[:5])%2==1]
Q=[r for r in carrier if tuple(2*x for x in r[5:])!=(-1,-1,-1)]
ck(len(G)==56 and len(carrier)==64 and len(Q)==48,'same native56 64 48')
M56=sum((r[:5,0]*r[:5,0].T for r in G),sp.zeros(5))
M48=sum((r[:5,0]*r[:5,0].T for r in Q),sp.zeros(5))
ck(M56==M48==12*sp.eye(5),'all D5 Cartan quadratic responses match')
Y=sp.Matrix([F(-1,3)]*3+[F(1,2)]*2+[0]*3)
traces={str(n):{'E7_56':str(sum(r.dot(Y)**n for r in G)),'native_48':str(sum(r.dot(Y)**n for r in Q))} for n in [1,2,3,4,5,6]}
ck(traces['2']=={'E7_56':'10','native_48':'10'},'same actual hypercharge quadratic trace')
ck(traces['4']=={'E7_56':'25/6','native_48':'95/18'},'nonlinear hypercharge discriminator')
ck(traces['5']=={'E7_56':'0','native_48':'5/3'},'fifth hypercharge moment detects chirality')
half56=[r for r in G if all(abs(x)==sp.Rational(1,2) for x in r)]
vec56=[r for r in G if any(r[:5,0]) and r not in half56]
sing56=[r for r in G if not any(r[:5,0])]
ck((len(half56),len(vec56),len(sing56))==(32,20,4),'actual56 D5 content32+20+4')
chiral56=Counter(sum(bool(x<0) for x in r[:5])%2 for r in half56)
chiral48=Counter(sum(bool(x<0) for x in r[:5])%2 for r in Q)
ck(chiral56=={0:16,1:16} and chiral48=={1:48},'quadratic equality is not chirality equality')
for rr in [r for r in R if not any(r[5:])]:
    ri=ch.ridx[tuple(int(2*x) for x in rr)];rm=ch.opp[ri]
    values=[]
    for space in [G,Q]:
        members={ch.ridx[tuple(int(2*x) for x in r)] for r in space}
        val=F(0)
        for j in members:
            b=ch.bracket(rm,j)
            ck(all(k in members for k in b),'actual D5 root action preserves both modules')
            for k,c in b.items():val+=c*ch.bracket(ri,k).get(j,0)
        values.append(val)
    ck(values==[F(-12),F(-12)],'all40 actual D5 root-pair charge traces match')

E7=[r for r in R if r.dot(alpha)==0]
E6=[r for r in E7 if r.dot(beta)==0]
D5=[r for r in E6 if r.dot(lam)==0]
ck((len(E7),len(E6),len(D5))==(126,72,40),'actual marked orbit stabilizer chain')
N=sp.Matrix.hstack(alpha,beta)
w7=beta-alpha*(alpha.dot(beta)/2)
w6=lam-N*(N.T*N).inv()*N.T*lam
w5=sp.Matrix([sp.Rational(-1,2)]*5+[0]*3)
def orbit(v,roots):
    # Scale by12 so all weights and actual root reflections are integral.
    rr=[tuple(int(2*x) for x in r) for r in roots]
    first=tuple(int(12*x) for x in v);seen={first};todo=[first]
    while todo:
        x=todo.pop()
        for r in rr:
            dot=sum(a*b for a,b in zip(x,r));ck(dot%4==0,'exact reflected weight integrality')
            y=tuple(a-dot//4*b for a,b in zip(x,r))
            if y not in seen:seen.add(y);todo.append(y)
    return seen
orbit_sets=[orbit(v,rs) for v,rs in [(alpha,R),(w7,E7),(w6,E6),(w5,D5)]]
orbits=list(map(len,orbit_sets))
ck(orbits==[240,56,27,16],'actual orbit ladder240 56 27 16')
ck(tuple(int(12*x) for x in lam-alpha*(alpha.dot(lam)/2)) in orbit_sets[1],'next marked weight belongs to the same56 orbit')
ck(tuple(int(12*x) for x in w5) in {v[:5]+(0,0,0) for v in orbit_sets[2]},'actual27 contains the selected same-chirality D5 spinor projection')
ck(Y.dot(w5)==0,'existing hypercharge lies in the final standard A4 Cartan')
norms=[v.dot(v) for v in [alpha,w7,w6,w5]]
ck(norms==list(map(sp.Rational,['2','3/2','4/3','5/4'])),'native projection norm ladder')
ck(len([r for r in D5 if r.dot(w5)==0])==20,'actual final A4 stabilizer roots')
result={'claim':'MATH.ESERIES.QUADRATIC.BUDGET.20261106','quadratic_response_matrix56_and48':'12 I5','E7_56_D5_branching':'16+16bar+2*10+4*1','hypercharge_traces':traces,'chirality_counts56':dict(chiral56),'chirality_counts48':dict(chiral48),'actual_Weyl_orbit_sizes':orbits,'projected_weight_norms':list(map(str,norms)),'stabilizer_root_counts':[126,72,40,20],'scopes':['Internal source current characters and actual root reflections, not physical fermion species.','Equal quadratic responses do not identify chiral representations; fourth moments distinguish them.','Known orbit hierarchy and norm formulas, newly compared with native TFPT48; no new physical constant.'],'checks':dict(env['checks'])}
result['response_definition']='Internal character chi_R(s)=Tr_R exp(s Y0); these are character moments, not the affine vacuum current OPE level or a spacetime fermion loop without an additional physical identification.'
result['actual_D5_root_pair_trace_both_modules']='-12 for every opposite root pair in the original sign convention'
result['orbit_flag_scope']='Nested stabilizer flag with a new marked weight from the preceding orbit at each stage, not repeated projection of the same single vector. Final D5 spinor is the odd chirality present in this27.'
Gneutral=[r for r in R if sum(r)==2]
pair_charges=[Y[i]+Y[j] for i,j in itertools.combinations(range(8),2)]
ck(Counter(r.dot(Y) for r in Gneutral)==Counter(pair_charges+[-v for v in pair_charges]),'source-selected neutral E7 56 is literal A7 pair charges28 plus28bar')
int28=[r for r in Gneutral if all(x.q==1 for x in r)]
half28=[r for r in Gneutral if all(x.q==2 for x in r)]
for block in [int28,half28]:
    members={ch.ridx[tuple(int(2*x) for x in r)] for r in block}
    ck(len(members)==28,'actual A7 pair module size')
    for ii,jj in itertools.permutations(range(8),2):
        r=sp.zeros(8,1);r[ii]=2;r[jj]=-2;ri=ch.ridx[tuple(r)]
        ck(all(k in members for j in members for k in ch.bracket(ri,j)),'all actual A7 root steps retain28 or dual28')
result['Q_neutral_56_hypercharge_moments']={str(n):str(sum(r.dot(Y)**n for r in Gneutral)) for n in [2,4,5,6]}
result['eight_channel_index_identity']='sum_(i<j)(y_i+y_j)^2=(N-2)sum y_i^2+(sum y_i)^2; for traceless N8 each28 has index6, while complete56 has coefficient12.'
result['native48_original_source_Q_histogram']={str(k):v for k,v in Counter(sum(r) for r in Q).items()}
result['source_sha256']=env['PINS']
result['all_three_lepton_factor_trace_recoding']={'c_e':'128/56=16/7','c_mu':'64/48=4/3','c_tau':'56/48=7/6','scope':'Existing equal-current-Gram projector traces. This is a representation of existing coefficients, not a derivation of the mass readout rule.'}
(OUT/'quadratic-budget-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
