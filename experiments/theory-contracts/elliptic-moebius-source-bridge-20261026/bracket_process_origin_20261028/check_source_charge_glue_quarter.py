"""Exact relation of the given charged fold holonomy and internal E8 glue clock."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
FOLD=REPO/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json'
WARD=REPO/'experiments/theory-contracts/seam-joint-ward-state-20261014/PROOF.txt'
PINS={FOLD:'253f370a9c84b70c28169ed0aad83d1866245804227e44aa0189e823382f5639',
      WARD:'7c0cd92dfc5eca88260baae0af9edad66f958a27ce0d15d7965414b3971f8da0'}
checks=[]
def need(value,label):
    if not value:raise RuntimeError(label)
    checks.append(label)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def mv(A,x):return tuple(dot(row,x) for row in A)
def ine8(x):
    integral=all(z.denominator==1 for z in x)
    half=all(z.denominator==2 for z in x)
    return (integral or half) and sum(x)%2==0
def counter(values):return {str(k):v for k,v in sorted(Counter(values).items())}

for path,pin in PINS.items():need(hashlib.sha256(path.read_bytes()).hexdigest()==pin,'source pin '+path.name)
fold=json.loads(FOLD.read_text())
M=tuple(tuple(F(z) for z in row) for row in fold['explicit_Weyl_matrix'])
MT=tuple(zip(*M))
need(all(dot(M[i],M[j])==int(i==j) for i in range(8) for j in range(8)),'given M is orthogonal')
roots=[]
for i,j in combinations(range(8),2):
    for a,b in product((-1,1),repeat=2):
        root=[F(0)]*8;root[i]=F(a);root[j]=F(b);roots.append(tuple(root))
for signs in product((-1,1),repeat=8):
    if signs.count(-1)%2==0:roots.append(tuple(F(z,2) for z in signs))
rootset=set(roots)
need(len(rootset)==240 and all(ine8(r) and dot(r,r)==2 for r in roots),'standard E8 roots')
need({mv(M,r) for r in roots}==rootset,'given M preserves the entire E8 root lattice')
q=mv(M,(F(1),)*8)
ward=(F(0),F(0),F(1),F(0),F(0),F(1),F(0),F(0))
need(q==tuple(2*z for z in ward),'same-source charge q maps to twice the Ward root')
need(ward in rootset and dot(q,ward)==4,'Ward field has original source charge four')
v=(F(0),)*5+(F(1,2),)*3
c0=mv(M,tuple(F(j,4) for j in range(8)))
need(c0==tuple(map(F,(-1,0,2,-1,0,F(3,2),F(-1,2),F(-1,2)))),'explicit geometric-quarter cocharacter')
integer_error=tuple(a+b for a,b in zip(c0,v))
need(integer_error==tuple(map(F,(-1,0,2,-1,0,2,0,0))) and ine8(integer_error),
     'c0 plus internal glue cocharacter is an E8 lattice vector')
need(all((dot(c0,r)+dot(v,r))%1==0 for r in roots),'C0 equals inverse internal glue clock on all roots')
need(sum((dot(c0,r)-dot(v,r))%1!=0 for r in roots)==128,'C0 differs from internal clock on all128 halfroots')
need(ine8(tuple(4*z for z in v)) and not ine8(tuple(2*z for z in v)),'internal clock has exact order four')
x=tuple(map(F,(-2,-2,-2,3,3,0,0,0)))
res3={r for r in roots if dot(x,r)%5==3}
creation={r for r in res3 if dot(v,r)%1==F(1,2)}
need(len(res3)==50 and len(creation)==30,'full U5 sector50 has joint glue-grade-two subspace30')
need({r for r in res3 if dot(c0,r)%1==F(1,2)}==creation,'geometric-quarter joint field component is exactly creation30')
need(counter(dot(c0,r)%1 for r in res3)=={'1/4':20,'1/2':30},'source geometric quarter separates30 and20')
cases=[]
for alpha in (F(1,4),F(3,4)):
    mu=mv(M,tuple((F(j)+alpha)/8 for j in range(8)))
    need(ine8(tuple(16*z for z in mu)) and not ine8(tuple(8*z for z in mu)),
         'given charged fold holonomy has order16 '+str(alpha))
    difference=tuple(2*mu[i]-alpha*q[i]/4-c0[i] for i in range(8))
    need(all(z==0 for z in difference),'exact charged/charge/geometric factorization '+str(alpha))
    # No state, time, background or holonomy is changed. This is a relation
    # among commuting characters already defined on the stipulated module.
    need(all((2*dot(mu,r)-alpha*dot(q,r)/4+dot(v,r))%1==0 for r in roots),
         'g_alpha squared equals charge phase times inverse glue on every root '+str(alpha))
    raw_ward=(2*dot(mu,ward))%1
    raw={r for r in res3 if (2*dot(mu,r))%1==raw_ward}
    correct={r for r in res3 if (2*dot(mu,r)-alpha*dot(q,r)/4)%1==F(1,2)}
    need(correct==creation,'joint actual holonomy and charge recovers creation component '+str(alpha))
    need(raw!=creation,'discarding source charge phase gives a different component '+str(alpha))
    cases.append({'alpha':str(alpha),'raw_g2_ward_phase_exponent':str(raw_ward),
                  'raw_g2_ward_eigenspace_inside_U5_residue3':len(raw),
                  'raw_selected_creation_roots':len(raw&creation),
                  'retained_charge_joint_selection':len(correct),
                  'source_charge_histogram_on_creation30':counter(dot(q,r) for r in creation)})
need([x['raw_g2_ward_eigenspace_inside_U5_residue3'] for x in cases]==[1,11],
     'raw quarter comparison has exact one and eleven dimensional components')
result={'verdict':'GIVEN_CHARGED_FOLD_CLOCK_AND_SOURCE_CHARGE_RECOVER_INVERSE_INTERNAL_GLUE',
        'formula':'g_alpha^2 = exp(2 pi i alpha Q_source/4) C_internal^(-1)',
        'geometric_cocharacter_target':list(map(str,c0)),
        'lattice_error_c0_plus_v':list(map(str,integer_error)),
        'proof_for_all_lattice_states':'c0+v belongs to E8=E8*, not only a finite root test',
        'creation_field_component':'U5 residue3 and (g_alpha^2 exp(-2 pi i alpha Q_source/4)) eigenvalue -1',
        'full_U5_residue3_dimension':50,'creation_joint_dimension':30,'cases':cases,
        'checks':len(checks),'check_names':checks,
        'sources':{str(p.relative_to(REPO)):pin for p,pin in PINS.items()},
        'scope':'Characters of the given marked E8 Fock/vertex representation. Full charged g_alpha is retained. No new physical clock, modified state, dropped phase, affine winding operator, NS/R transition, local P1 field map, executed measurement or physical OS time is inferred.'}
(HERE/'source-charge-glue-quarter-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
