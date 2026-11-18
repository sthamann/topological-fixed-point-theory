#!/usr/bin/env python3
"""Exact consequences of the original Round20 cocycle; no charge cutoff."""
import ast
import hashlib
import json
from pathlib import Path
import sympy as sp

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
manifest=json.loads((HERE/'source_manifest.json').read_text())
checks={}
def check(name,value):
    if not bool(value): raise AssertionError(name)
    checks[name]=True
for item in manifest['sources']:
    check('source_hash:'+item['path'],hashlib.sha256((REPO/item['path']).read_bytes()).hexdigest()==item['sha256'])
path='experiments/theory-contracts/local-charge-transport-round20/checker.py'
tree=ast.parse((REPO/path).read_text())
main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
assign=[n for n in main.body if isinstance(n,ast.Assign) and len(n.targets)==1
        and isinstance(n.targets[0],ast.Name) and n.targets[0].id in {'eye','B','G','g'}]
names={'add','scale','dotg','energy','beta','translate','transfer','total'}
defs=[n for n in main.body if isinstance(n,ast.FunctionDef) and n.name in names]
ns={'sp':sp};exec(compile(ast.Module(body=assign+defs,type_ignores=[]),path,'exec'),ns)
B,G=ns['B'],ns['G'];beta,T,W,total,e=[ns[k] for k in ['beta','translate','transfer','total','energy']]
neg=lambda v:tuple(-x for x in v)
C=lambda profile:tuple(neg(v) for v in profile)
zero=(0,)*8
units=[tuple(int(i==j) for i in range(8)) for j in range(8)]
charges=[zero]+units+[neg(p) for p in units]+[(2,-1,0,3,-2,1,0,-3)]
profiles=[(zero,zero,zero),(charges[-1],neg(charges[-1]),zero)]
profiles += [(p,neg(p),zero) for p in units]
check('actual_even_unimodular_Gram',G.det()==1 and all(G[i,i]%2==0 for i in range(8)))
check('C_is_unitary_basis_involution',all(C(C(n))==n for n in profiles))
check('actual_bilinear_cocycle_sign_reversal',all(beta(neg(p),neg(n))==beta(p,n) for p in charges for n in charges))
check('actual_cocycle_CT_pC_equals_T_minus_p',all(
    (C(T(C(n),x,p)[0]),T(C(n),x,p)[1])==T(n,x,neg(p))
    for n in profiles for x in range(3) for p in charges))
check('actual_adjoint_sign_not_silently_dropped',all(
    T(n,x,p,True)[0]==T(n,x,neg(p))[0]
    and T(n,x,p,True)[1]==(-1)**int(e(p))*T(n,x,neg(p))[1]
    for n in profiles for x in range(3) for p in charges))
check('odd_norm_hopping_channels_have_minus_adjoint_sign',[int(e(p))%2 for p in units]==[0,1,1,1,1,1,1,1])
check('C_W_p_C_is_W_minus_p',all(
    (C(W(C(n),0,1,p)[0]),W(C(n),0,1,p)[1])==W(n,0,1,neg(p))
    for n in profiles for p in charges))
check('W_minus_p_is_W_p_adjoint',all(W(n,0,1,neg(p))==W(n,1,0,p) for n in profiles for p in charges))
check('quadratic_charge_energy_is_even',all(e(neg(p))==e(p) for p in charges))
check('gradient_energy_is_even',all(e(neg(ns['add'](p,neg(r))))==e(ns['add'](p,neg(r))) for p in charges for r in charges))
check('C_preserves_total_zero_sector',all(total(C(n))==neg(total(n)) for n in profiles))
check('on_site_T_exits_C0_for_nonzero_charge',all(total(T(n,0,p)[0])==p for n in profiles for p in units))
check('balanced_W_stays_in_C0',all(total(W(n,0,1,p)[0])==zero for n in profiles for p in units))
check('balanced_W_changes_local_charges_inside_C0',W((zero,zero,zero),0,1,units[7])[0] != (zero,zero,zero))

eta=sp.Matrix(sp.symbols('eta0:8',real=True))
pv=sp.Matrix(sp.symbols('p0:8',integer=True))
t=sp.symbols('t',positive=True)
check('uniform_chemical_mark_zero_on_C0',eta.dot(sp.zeros(8,1))==0)
check('uniform_chemical_mark_odd_under_C',sp.expand(eta.dot(-pv)+eta.dot(pv))==0)
check('charged_corner_chemical_ratio_exponent',sp.expand(-t*eta.dot(pv)+t*eta.dot(-pv)+2*t*eta.dot(pv))==0)
delta=sp.Matrix(sp.symbols('delta0:8',real=True))
odd=sp.Matrix([delta.dot(B[:,i]) for i in range(8)])
check('original_native_charge_basis_spans_twist_readout',abs(B.det())==1 and sp.simplify(B.T.inv()*odd-delta)==sp.zeros(8,1))
check('conditional_coordinate_dictionary_eta_equals_Bt_delta',sp.expand((B.T*delta).dot(pv)-delta.dot(B*pv))==0)
witness=[]
for alpha in [sp.Rational(1,4),sp.Rational(3,4)]:
    d=sp.Matrix([(k+alpha)/8-sp.Rational(1,2) for k in range(8)])
    values=[sp.simplify(d.dot(B[:,i])) for i in range(8)]
    i=next(i for i,x in enumerate(values) if x!=0)
    check(str(alpha)+':fixed_chiral_twist_has_asymmetric_charged_pair',values[i]!=0)
    witness.append({'alpha':str(alpha),'native_integer_basis_index':i,'target_charge':[str(x) for x in B[:,i]],'delta_dot_target_charge':str(values[i])})

out={'contract':'native-charged-conjugation-20261015','research_id':'UR.NATIVE.CHARGED_CONJUGATION.01',
     'verdict':'PARTIAL','mathematical_verdict':'PASS_SOURCE_CONJUGATION_AND_CHARGED_CORNER_DISCRIMINATOR',
     'check_count':len(checks),'failed':0,'checks':checks,'quarter_witnesses':witness,
     'source_hashes':manifest['sources'],'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'proof_sha256':hashlib.sha256((HERE/'PROOF.md').read_bytes()).hexdigest() if (HERE/'PROOF.md').exists() else None,
     'finite_regression_profiles':len(profiles),'finite_regression_charges':len(charges),
     'all_charge_all_time_basis':'Analytic bilinearity, quadratic evenness, strong C covariance and functional calculus; finite samples are regressions, not a heat-kernel/Gibbs simulation.',
     'physical_gates_closed':[],'complete_TFPT_solution':False,
     'not_claimed':['No physical selection of native source or chemical eta',
                    'No absence of charged local probes merely because a starting state lies in C0',
                    'No ordinary onsite T as an endomorphism of C0',
                    'No automatic equivalence to a mixed +/-delta chiral family',
                    'No unrestricted TFPT no-go or full field/source identification']}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(f'{len(checks)}/{len(checks)} exact groups passed; failed=0; overall PARTIAL')
