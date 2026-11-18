"""Canonical symmetric four-current probe of the original E7 quartic.
This is a declared charged current-state probe, not a physical selection law.
"""
import ast,itertools,json,argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from pathlib import Path

REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
parser=argparse.ArgumentParser();parser.add_argument('--out',required=True,type=Path)
OUT=parser.parse_args().out;OUT.mkdir(parents=True,exist_ok=True)
bootstrap=REPO/'experiments/theory-contracts/e-series-compiler-closure-20261023/e7_full_source_response_20261106.py'
tree=ast.parse(bootstrap.read_text());nodes=[]
for n in tree.body:
    if isinstance(n,ast.For):break
    if any(isinstance(t,ast.Name) and t.id in ['parser','OUT'] for t in ast.walk(n)):continue
    nodes.append(n)
env={};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(bootstrap),'exec'),env)
ch,af,ip,ck=(env[k] for k in ['ch','af','ip','ck'])
R=ch.roots;vac={():F(1)}
rho=tuple(F(7-2*j,2) for j in range(8))
out=[]
for direction in [1,-1]:
    alpha=(direction,)*8;ap=ch.ridx[alpha];am=ch.opp[ap]
    G1=[i for i,r in enumerate(R) if ip(alpha,r)==4]
    V={am:{():1}}
    for degree in range(1,5):
        nxt=defaultdict(Counter)
        for a,root in enumerate(G1):
            for basis,poly in V.items():
                for target,coeff in ch.bracket(root,basis).items():
                    for m,c in poly.items():nxt[target][tuple(sorted(m+(a,)))]+=int(coeff)*c
        V={b:{m:c for m,c in p.items() if c} for b,p in nxt.items() if any(p.values())}
    ck(set(V)=={ap},'source quartic output')
    Q=V[ap];state={}
    for mon,c in Q.items():
        perms=set(itertools.permutations(mon))
        for word in perms:
            v=vac
            for i in reversed(word):v=af.act(G1[i],-1,v)
            for key,value in v.items():state[key]=state.get(key,F(0))+F(c,len(perms))*value
    state={k:c for k,c in state.items() if c}
    unit=af.act(ap,-3,af.act(ap,-1,vac))
    def pair_unit(st):return af.act(am,1,af.act(am,3,st)).get((),F(0))
    unit_norm=pair_unit(unit);amplitude=pair_unit(state)
    ck(unit_norm==1,'same-source unit lattice head state')
    ck(amplitude!=0,'quartic sym current probe survives integrable quotient')
    energies={}
    for a in [F(1,4),F(3,4)]:
        values=set()
        for key in state:
            L=-sum(m for m,_ in key)
            p=tuple(F(x,2) for x in af.weight(key))
            ck(L==4 and p==tuple(alpha),'whole probe has grade4 and charge2alpha')
            values.add(8*L-ip(rho,p)+F(2*a-1,2)*sum(p))
        ck(len(values)==1,'actual H8 eigenstate, not average energy')
        energies[str(a)]=str(next(iter(values)))
    out.append({'root_orientation':direction,'quartic_monomials':len(Q),'sym_current_PBW_terms':len(state),'unit_reference_norm':str(unit_norm),'quartic_probe_unit_amplitude':str(amplitude),'source_H8_eigenvalues':energies,'physical_identification':'none; specified symmetric current probe, inherited source time only'})
    print(out[-1],flush=True)
sp=env['sp'];basis=env['R'];I=sp.eye(8);C=I
for i in range(8):
    root=basis[:,i];C=(I-root*root.T)*C
v=alpha0=sp.ones(8,1)/2;orbit=[]
for k in range(1,61):
    orbit.append(tuple(v));v=C*v
    if v==alpha0:break
ck(len(orbit)==30,'actual source-basis Coxeter root orbit has30 points')
ck(C*alpha0 not in [alpha0,-alpha0],'Coxeter action does not preserve the two charged probes')
counter={'root_orbit_size':len(orbit),'preserves_two_charged_probe_directions':False,'H8_energies_along_orbit_at_quarter1_4':list(map(str,sorted({32-2*ip(rho,r)-F(1,2)*sum(r) for r in orbit})))}
(OUT/'quartic-clock-results.json').write_text(json.dumps({'claim':'MATH.ESERIES.QUARTIC.CLOCK.PROBE.20261106','probe_definition':'Q mapped by fully symmetrized product of four inherited mode-minus-one currents; this is a declared charged current-state probe, not a self-adjoint observable or an original physical insertion selection theorem.','results':out,'Coxeter_action_countertest':counter,'checks':dict(env['checks']),'source_sha256':env['PINS']},indent=2)+'\n')
