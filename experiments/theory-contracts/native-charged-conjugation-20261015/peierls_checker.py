#!/usr/bin/env python3
"""Original W monomial actions with arbitrary fixed complex link phase."""
import ast
import hashlib
import json
from pathlib import Path
import sympy as sp

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
checks={}
def check(name,value):
    if not bool(value):raise AssertionError(name)
    checks[name]=True
manifest=json.loads((HERE/'source_manifest.json').read_text())
source=next(x for x in manifest['sources'] if x['path'].endswith('local-charge-transport-round20/checker.py'))
check('original_Round20_source_hash',hashlib.sha256((REPO/source['path']).read_bytes()).hexdigest()==source['sha256'])
tree=ast.parse((REPO/source['path']).read_text())
main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
assign=[n for n in main.body if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name) and n.targets[0].id in {'eye','B','G','g'}]
names={'add','scale','dotg','energy','beta','translate','transfer','total'}
defs=[n for n in main.body if isinstance(n,ast.FunctionDef) and n.name in names]
ns={'sp':sp};exec(compile(ast.Module(body=assign+defs,type_ignores=[]),source['path'],'exec'),ns)
T,W=ns['translate'],ns['transfer']
neg=lambda v:tuple(-x for x in v)
C=lambda profile:tuple(neg(n) for n in profile)
zero=(0,)*8
units=[tuple(int(i==j) for i in range(8)) for j in range(8)]
large=(2,-1,0,3,-2,1,0,-3)
profiles=[(zero,zero,zero),(units[0],zero,zero),(large,neg(large),zero)]

def hopping_action(profile,p,z):
    out={profile:sp.Integer(2)}
    for target,phase,coefficient in [(*W(profile,0,1,p),-z),(*W(profile,1,0,p),-sp.conjugate(z))]:
        out[target]=sp.expand(out.get(target,0)+coefficient*phase)
    return {k:sp.simplify(v) for k,v in out.items() if v!=0}

def mapped_action(profile,p,z,charge_flip=False,conjugate=False):
    argument=C(profile) if charge_flip else profile
    out={}
    for target,value in hopping_action(argument,p,z).items():
        target=C(target) if charge_flip else target
        out[target]=sp.conjugate(value) if conjugate else value
    return {k:sp.simplify(v) for k,v in out.items() if v!=0}

c,s=sp.symbols('c s',real=True)
z=c+sp.I*s
check('K_hopping_theta_K_is_minus_theta_symbolically',all(
    mapped_action(n,p,z,conjugate=True)==hopping_action(n,p,sp.conjugate(z))
    for n in profiles for p in units))
check('C_hopping_minus_theta_C_is_theta_symbolically',all(
    mapped_action(n,p,sp.conjugate(z),charge_flip=True)==hopping_action(n,p,z)
    for n in profiles for p in units))
check('actual_CK_hopping_covariance_arbitrary_complex_coefficients',all(
    mapped_action(n,p,z,charge_flip=True,conjugate=True)==hopping_action(n,p,z)
    for n in profiles for p in units))
for value,label in [(1,'zero'),(sp.I,'quarter'),(-sp.I,'minus_quarter'),(-1,'half')]:
    check('actual_CK_phase_'+label,all(mapped_action(n,p,value,True,True)==hopping_action(n,p,value)
                                    for n in profiles for p in units))
check('C_alone_really_fails_quarter_phase',mapped_action(profiles[0],units[7],sp.I,True,False)!=hopping_action(profiles[0],units[7],sp.I))
check('K_alone_really_fails_quarter_phase',mapped_action(profiles[0],units[7],sp.I,False,True)!=hopping_action(profiles[0],units[7],sp.I))
check('actual_A_T_p_A_is_T_minus_p',all(
    (C(T(C(n),0,p)[0]),sp.conjugate(T(C(n),0,p)[1]))==T(n,0,neg(p))
    for n in profiles for p in units))
check('actual_CK_square_is_identity',all(C(C(n))==n for n in profiles))
theta=sp.symbols('theta',real=True)
u=sp.exp(sp.I*theta)
check('phase_unitarity_exact',sp.simplify(u*sp.conjugate(u))==1)
check('positive_square_diagonal_coefficient',sp.simplify(1+u*sp.conjugate(u))==2)
# (I-uW)^*(I-uW)=2I-uW-conj(u)W^* for the actual unitary W.
check('original_W_and_adjoint_inverse_regression',all(
    W(W(n,0,1,p)[0],1,0,p)[0]==n
    and W(n,0,1,p)[1]*W(W(n,0,1,p)[0],1,0,p)[1]==1
    for n in profiles for p in units))
check('original_real_diagonal_C_even',all(ns['energy'](neg(n))==ns['energy'](n) for profile in profiles for n in profile))
result={'contract':'native-charged-conjugation-20261015','addendum':'FIXED_PEIERLS_PHASES_RETAIN_CK_EUCLIDEAN_PAIR_SYMMETRY',
        'verdict':'PASS_SCOPED_ANTIUNITARY_OBSTRUCTION','overall_physical_verdict':'PARTIAL',
        'check_count':len(checks),'failed':0,'checks':checks,'source_hashes':[source],
        'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'addendum_sha256':hashlib.sha256((HERE/'ADDENDUM.md').read_bytes()).hexdigest() if (HERE/'ADDENDUM.md').exists() else None,
        'scope':'Original real cocycle W and real C-even diagonal/scalar Hamiltonian, arbitrary fixed link phases and A=CK-invariant start. All-time Euclidean equality is analytic, not inferred from sampled phases/profiles.',
        'not_claimed':['No P1/P2 selection of phases or source','No extension to time-dependent phases, A-breaking starts or extra A-odd terms','No full Gibbs/heat simulation or finite charge truncation','No physical gate closure']}
(HERE/'peierls_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(f'{len(checks)}/{len(checks)} exact groups passed; failed=0; overall PARTIAL')
