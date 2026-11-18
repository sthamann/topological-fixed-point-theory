"""One actual E8 root: integrated current and charge-weighted pair response.

The geometric seam class is kept separate from the internal root loop.
"""
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import sys
import sympy as s
sys.path.insert(0,'/Users/stefanhamann/Projekte/tfpt-theoryv4')
from tfpt_explorer.current_block_geometry import _actual_current_source,_kappa,_bracket
from tfpt_explorer.neutral_source_response import _native_neutral_branch
from tfpt_explorer.source_charge_response import hypercharge_vector
from tfpt_explorer.flavor_path_transport import exact_clock_dictionary

src=_actual_current_source();alg=src['chevalley'];n=_native_neutral_branch()
root=tuple([-1]*8);lam=s.Matrix(root)/2
assert root in alg.ridx and (lam.T*lam)[0]==2
# Genuine SM gauge invariance uses brackets, not charge counting alone.
Y=hypercharge_vector()
simple=[(2,-2,0,0,0,0,0,0),(0,2,-2,0,0,0,0,0),(0,0,0,2,-2,0,0,0)]
gauge_roots=simple+[tuple(-x for x in r) for r in simple]
singlets=[]
for r in alg.roots:
    beta=s.Matrix(r)/2
    if (Y.T*beta)[0] or any(_bracket(src,{alg.ridx[a]:Fraction(1)},{alg.ridx[r]:Fraction(1)}) for a in gauge_roots):
        continue
    singlets.append(r)
assert len(singlets)==20
family_roots=[tuple([-1]*5+list(tail)) for tail in n['family']]
opposite_family=[tuple(-x for x in r) for r in family_roots]
family_A3=[r for r in alg.roots if not any(r[:5])]
assert set(singlets)==set(family_roots+opposite_family+family_A3)
M=exact_clock_dictionary()['M'];U=s.diag(1,s.I,-s.I)
M4=s.diag(1,M);U4=s.diag(1,U)
fixed=(M4-s.eye(4)).col_join(U4-s.eye(4)).nullspace()
assert fixed==[s.Matrix([1,0,0,0])]
# The A3 root-only subspace consists of zero-diagonal sl4 matrices.
# Any simultaneous invariant there must also be in the full commutant.
variables=s.symbols('z:12')
Z=s.zeros(4);k=0
for a in range(4):
    for b in range(4):
        if a!=b:Z[a,b]=variables[k];k+=1
fixed_root_solution=s.linsolve(list(Z*M4-M4*Z)+list(Z*U4-U4*Z),variables)
assert next(iter(fixed_root_solution))==(0,)*12
H=alg.hvec_index(root)
assert _kappa(src,H,H)==2
# The ACTUAL Cartan zero-mode action has integer lambda.beta charges.
charges={}
for r in alg.roots:
    beta=s.Matrix(r)/2;q=(lam.T*beta)[0]
    assert q.q==1
    index=alg.ridx[r]
    assert _bracket(src,H,{index:Fraction(1)})==({index:Fraction(int(q))} if q else {})
    charges[r]=q
full=Counter(charges.values())
assert full==Counter({-2:1,-1:56,0:126,1:56,2:1})
assert sum(charges.values())==0
moment=s.zeros(8)
for r in n['spinor']:
    beta=s.Matrix(r)/2
    moment+=charges[r]*beta*beta.T
assert moment==-n['R_matrix']
assert s.trace(moment)==0 and s.trace(moment**2)/2==48
nu=[tuple([-1]*5+list(tail)) for tail in n['family']]
assert [charges[r] for r in nu]==[2,1,1,1]
selected=Counter(charges[r] for r in n['retained'])
assert selected==Counter({-1:15,0:30,1:3})
# J=H/sqrt2: affine mode norm n. Primitive X_n=J_n/n has 1/n.
m=s.symbols('m',integer=True,positive=True)
assert s.simplify(2*m/2/m**2)==1/m
result={
 'native_root_lambda':list(map(str,lam)),
 'SM_gauge_root_centralizer_dimension':20,
 'SM_root_centralizer_decomposition':'12 A3 roots + 4 C + 4 conjugate C',
 'common_clock_fixed_charged_root_lines':['C0','C0 adjoint'],
 'uniqueness_scope':'the declared native C4 chirality and U4/M4 lift; not the whole Cartan current space',
 'actual_Cartan_kappa':str(_kappa(src,H,H)),
 'all_root_zero_mode_action_checked':len(alg.roots),
 'log_kernel_mode':'J=H_lambda/sqrt2, <J_n J_-n>=n, X_n=J_n/n, <X_n X_-n>=1/n',
 'radial_log_kernel':'-log(1-r^2*exp(i*(theta-phi))), r<1',
 'root_pair_weighting':'sum_64 (lambda.beta) beta beta^T = -R',
 'R_norm_squared':'48',
 'integer_charge_distributions':{
   '240':{str(q):c for q,c in sorted(full.items())},
   '64':{str(q):c for q,c in sorted(Counter(charges[r] for r in n['spinor']).items())},
   '48':{str(q):c for q,c in sorted(selected.items())}},
 'finite_determinant_windings':{'adjoint':0,'64':0,'48':str(sum(charges[r] for r in n['retained'])),'four_neutrino_currents':5},
 'finite_winding_scope':'These finite representation determinant windings do not determine an independently specified boundary Fredholm spectral flow.',
 'first_proposed_identification_rejected':'The internal root loop, interpreted as the same pointwise seam clutching operator on the listed source fibers, does not have unit determinant winding.',
 'source_global_guard':'A1/E7 charge glue and zero modes cannot be discarded when interpreting the mean-zero oscillator kernel.',
 'not_derived':['raw P1 Calderon-to-current map','common physical coupling/time normalization','physical 3+1D execution','prime law','AI superiority'],
 'verdict':'EXACT_NATIVE_TWO_OBSERVABLES; INTERNAL_LOOP_IS_NOT_P1_UNIT_CLUTCHING_ON_THE_TESTED_FIBERS; FUNDAMENTAL_SOLUTION_OPEN'
}
Path(__file__).with_name('native_root_joint_operation_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2))
