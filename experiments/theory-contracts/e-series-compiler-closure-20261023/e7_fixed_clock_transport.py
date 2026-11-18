"""Can native E7 cuts be transported while both SM and the original H8 stay fixed?"""
import ast,itertools,json
from pathlib import Path
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
ch,ck,ip=(env[k] for k in ['ch','ck','ip']);R=[sp.Matrix(r)/2 for r in ch.roots]
Y=sp.Matrix([F(-1,3)]*3+[F(1,2)]*2+[0]*3);rho=sp.Matrix([F(7-2*j,2) for j in range(8)]);q=sp.ones(8,1)
sm=[]
for i,j in [(0,1),(1,2),(3,4)]:
    v=sp.zeros(8,1);v[i]=1;v[j]=-1;sm.append(v)
comm=[r for r in R if all(r.dot(v)==0 for v in sm+[Y])]
alpha=sp.Matrix([0]*5+[1,-1,0]);target=sp.ones(8,1)/2
result=[]
for a in [F(1,4),F(3,4)]:
    w=rho-(2*a-1)*q/2;zero=[r for r in comm if r.dot(w)==0]
    ck(all(r.dot(w)==0 for r in zero),'centralizer roots preserve fixed original H8')
    start=tuple(alpha);goal=tuple(target);todo=[start];words={start:[]}
    while todo:
        x=sp.Matrix(todo.pop())
        for r in zero:
            y=tuple(x-x.dot(r)*r)
            if y not in words:words[y]=words[tuple(x)]+[tuple(r)];todo.append(y)
    result.append({'quarter':str(a),'SM_commutant_roots':len(comm),'fixed_H_commutant_roots':len(zero),'fixed_H_commutant_rank':sp.Matrix.hstack(*zero).rank() if zero else 0,'alpha_H_charge':str(alpha.dot(w)),'target_H_charge':str(target.dot(w)),'normal_orbit_size':len(words),'target_in_orbit':goal in words,'word_roots':[[str(z) for z in r] for r in words.get(goal,[])]})
print(json.dumps(result,indent=2));(OUT/'fixed-clock-transport-results.json').write_text(json.dumps(result,indent=2)+'\n')
