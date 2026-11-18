"""Actual E8 bracket -> E7 quartic -> the pre-existing source E6 cubic.
Question: does polarized quartic reproduce the existing signed45 cubic?
Success: exact root/cocycle dictionary and all coefficients proportional.
Stop: failed dictionary or failed proportionality; no alpha-cubic inference.
"""
import ast
import argparse
import hashlib
import json
import itertools
from pathlib import Path
from collections import Counter,defaultdict
from fractions import Fraction as F
import sympy as sp

REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
parser=argparse.ArgumentParser();parser.add_argument('--out',required=True,type=Path)
OUT=parser.parse_args().out;OUT.mkdir(exist_ok=True,parents=True)
bootstrap=REPO/'experiments/theory-contracts/e-series-compiler-closure-20261023/e7_full_source_response_20261106.py'
tree=ast.parse(bootstrap.read_text());nodes=[]
for n in tree.body:
    if isinstance(n,ast.For):break
    if any(isinstance(t,ast.Name) and t.id in ['parser','OUT'] for t in ast.walk(n)):continue
    nodes.append(n)
env={}
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(bootstrap),'exec'),env)
ch=env['ch'];ck=env['ck'];ip=env['ip']
R=ch.roots
alpha=(0,0,0,0,0,2,-2,0)
ap=ch.ridx[alpha];am=ch.opp[ap]
G1=[i for i,r in enumerate(R) if ip(alpha,r)==4]
E7=[i for i,r in enumerate(R) if ip(alpha,r)==0]
ck(len(G1)==56 and len(E7)==126,'actual56 and126 root grading')
index={r:i for i,r in enumerate(G1)}
omega=sp.zeros(56)
for a,b in itertools.product(range(56),repeat=2):
    bracket=ch.bracket(G1[a],G1[b])
    ck(set(bracket)<={ap},'grade2 bracket in original root line')
    omega[a,b]=bracket.get(ap,0)
ck(omega.T==-omega and omega.rank()==56,'native symplectic56 form nondegenerate')

# Polynomial coefficients of (ad x)^4 e_-alpha, retaining original basis.
V={am:{():1}}
stages=[]
for degree in range(1,5):
    nxt=defaultdict(Counter)
    for a,root in enumerate(G1):
        for basis,poly in V.items():
            bracket=ch.bracket(root,basis)
            for target,coeff in bracket.items():
                ck(coeff.denominator==1,'integral original bracket')
                coeff=int(coeff)
                for mon,v in poly.items():
                    key=tuple(sorted(mon+(a,)))
                    nxt[target][key]+=coeff*v
    V={b:{m:c for m,c in p.items() if c} for b,p in nxt.items() if any(p.values())}
    stages.append({'degree':degree,'basis_directions':len(V),'polynomial_entries':sum(len(p) for p in V.values())})
    print('stage',stages[-1],flush=True)
ck(set(V)=={ap},'quartic ends in alpha root only')
Q=V[ap]
ck(bool(Q),'quartic nonzero')

actions={}
for s in E7:
    steps={}
    for i,r in enumerate(G1):
        b=ch.bracket(s,r)
        ck(all(t in index for t in b),'E7 preserves native56')
        if b:steps[i]=[(index[t],int(v)) for t,v in b.items()]
    actions[s]=steps
for s,steps in actions.items():
    deriv=Counter()
    for mon,c in Q.items():
        for old,mult in Counter(mon).items():
            for new,v in steps.get(old,[]):
                mm=list(mon);mm.remove(old);mm.append(new)
                deriv[tuple(sorted(mm))]+=c*mult*v
    ck(not any(deriv.values()),'all126 E7 Ward equations for quartic')
print('quartic',len(Q),'coefficients',dict(Counter(Q.values())),'all126 invariance PASS',flush=True)

# E6 plane inside the same E7: family A2 has roots e5-e6 and e6-e7.
u=(0,0,0,0,0,1,1,-2)
charge={i:F(ip(u,R[r]),2) for i,r in enumerate(G1)}
hist=Counter(charge.values());ck(hist=={F(-3):1,F(-1):27,F(1):27,F(3):1},'56 as1+27+27bar+1 by actual E6 charge')
sing={k:next(i for i in range(56) if charge[i]==k) for k in [F(-3),F(3)]}
polarized={}
for sc,si in sing.items():
    p={}
    for mon,c in Q.items():
        if mon.count(si)==1 and all(charge[j]==-sc/3 for j in mon if j!=si):
            mm=list(mon);mm.remove(si);p[tuple(mm)]=c
    ck(len(p)==45,'quartic polarization yields45 cubic monomials')
    polarized[str(sc)]=p
print('E6 charge pieces',dict(hist),'cubic polarizations',[(k,len(v)) for k,v in polarized.items()],flush=True)

# Existing signed cubic certificate, without a newly solved invariant fit.
old=json.loads((REPO/'experiments/theory-contracts/compiler-vacuum-current-cubic-20260918/certificate.json').read_text())['cubic']
base=[tuple(r) for r in old['base_roots_doubled']]
parent=json.loads((REPO/'experiments/theory-contracts/compiler-e6-seam-readout-20260918/certificate.json').read_text())['e6']
blocks=[[tuple(r) for r in b] for b in parent['three_27_blocks']]
base_family=tuple(sum(r[k] for r in base) for k in [5,6,7])
matches=[]
for sc,poly in polarized.items():
    target_set={R[G1[i]] for i in range(56) if charge[i]==-F(sc)/3}
    for orientation in [1,-1]:
        for block in blocks:
            oriented={tuple(orientation*x for x in r) for r in block}
            if target_set!=oriented:continue
            target_family=[F(sum(r[k] for r in oriented),27) for k in [5,6,7]]
            base_family=[F(sum(r[k] for r in base),27) for k in [5,6,7]]
            # Transport the original base by its actual A2 root operation.
            delta=tuple([0]*5+[int(target_family[k]*orientation-base_family[k]) for k in range(3)])
            mapped=[];signs=[]
            for r in base:
                if not any(delta):image=r;sign=1
                else:
                    image=tuple(r[k]+delta[k] for k in range(8))
                    b=ch.bracket(ch.ridx[delta],ch.ridx[r]);ck(b=={ch.ridx[image]:b.get(ch.ridx[image])} and abs(b.get(ch.ridx[image],0))==1,'literal source A2 transport to E7 27')
                    sign=int(b[ch.ridx[image]])
                if orientation==-1:
                    # Compact conjugation contributes - on each root vector.
                    image=tuple(-x for x in image);sign=-sign
                mapped.append(index[ch.ridx[image]]);signs.append(sign)
            expected={tuple(sorted(mapped[i] for i in m['indices'])):int(m['coefficient'])*sp.prod(signs[i] for i in m['indices']) for m in old['monomials']}
            ck(set(poly)==set(expected),'all45 original cubic supports match polarization')
            ratios={F(poly[m],int(expected[m])) for m in poly}
            ck(len(ratios)==1,'all45 original signed cubic coefficients match one source scale')
            matches.append({'singlet_charge':sc,'orientation':orientation,'quartic_to_existing_cubic_scale':str(next(iter(ratios))),'matched_monomials':45,'A2_root_doubled':delta})
ck(bool(matches),'source cubic is actually identified, not inferred from dimension')
print('OLD CUBIC MATCH',matches,flush=True)

result={'claim':'MATH.ESERIES.NATIVE.QUARTIC.DESCENT.20261106','scope':'Original source E8 bracket and internal E7/E6 representations. Known Freudenthal construction realized on native source; no physical coupling or alpha-polynomial identification.','normal_root_doubled':alpha,'native56_roots_doubled':[R[i] for i in G1],'symplectic_rank':56,'stages':stages,'quartic_monomials':[{'indices':m,'coefficient':c} for m,c in sorted(Q.items())],'quartic_coefficient_histogram':dict(Counter(Q.values())),'E7_root_Ward_tests':126,'E6_charge_histogram':{str(k):v for k,v in hist.items()},'source_cubic_matches':matches,'checks':dict(env['checks'])}
result['source_sha256']={**env['PINS'],str(bootstrap.relative_to(REPO)):hashlib.sha256(bootstrap.read_bytes()).hexdigest(),'experiments/theory-contracts/compiler-vacuum-current-cubic-20260918/certificate.json':hashlib.sha256((REPO/'experiments/theory-contracts/compiler-vacuum-current-cubic-20260918/certificate.json').read_bytes()).hexdigest()}
(OUT/'quartic-descent-results.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS',sum(env['checks'].values()),flush=True)
