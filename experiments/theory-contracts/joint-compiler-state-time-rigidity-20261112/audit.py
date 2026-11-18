from pathlib import Path
import argparse
import importlib.util
import hashlib
import json
import sys
import sympy as s

parser = argparse.ArgumentParser()
parser.add_argument('--repo', required=True, type=Path)
parser.add_argument('--output', type=Path, default=Path(__file__).with_name('results.json'))
args = parser.parse_args()
source = args.repo / 'verification/v1_e8_glue.py'
sys.path.insert(0, str(source.parent))
spec = importlib.util.spec_from_file_location('native_e8', source)
native = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)
roots = [s.Matrix([s.Rational(str(x)) for x in r]) for r in native.e8_roots()]
simple = s.Matrix([[s.Rational(str(x)) for x in r] for r in native.e8_simple_roots()])
assert len(roots) == 240 and all(r.dot(r) == 2 for r in roots)
assert sum((r*r.T for r in roots), s.zeros(8)) == 60*s.eye(8)
coordinates = [simple.T.inv()*r for r in roots]
highest = max((r for r in coordinates if all(x >= 0 for x in r)), key=lambda r: sum(r))
assert sorted(highest) == [2, 2, 3, 3, 4, 4, 5, 6]
assert sum(highest)+1 == 30
assert s.Rational(248,31) == 8

n,m = s.symbols('n m', integer=True, nonnegative=True)
energy = n + 2*m*n*(n-1)
assert s.expand(energy.subs(n,n+1)-energy) == 1+4*m*n
assert all((2*j*k*(k-1)) % 4 == 0 for j in range(3) for k in range(17))
# The all-integer quarter-turn identity follows analytically from n(n-1) even.
table=[]
for j in range(3):
    energies=[int(energy.subs({m:j,n:k})) for k in range(4)]
    assert energies[0:2] == [0,1] and energies[2] == 2+4*j
    table.append({'m':j,'energies_grade_0_1_2_3':energies,
                  'first_normalized_return':'exp(-t)',
                  'second_normalized_return':f'exp(-{energies[2]}*t)'})
# One Cartan oscillator is a literal subalgebra of the current module.
# Four states suffice for the two specified words; no boundary is crossed.
creation=s.zeros(4)
for k in range(3): creation[k+1,k]=s.sqrt(k+1)
vac=s.Matrix([1,0,0,0]); annihilation=creation.T
assert ((creation*vac).T*(creation*vac))[0] == 1
assert ((creation**2*vac).T*(creation**2*vac))[0] == 2
assert (annihilation*creation-creation*annihilation)*vac == vac
assert (annihilation*creation-creation*annihilation)*creation*vac == creation*vac
dims=[4,5,20]
weights=[[s.Rational(1,3)]*3,
         [s.Rational(1,2),s.Rational(1,3),s.Rational(1,6)],
         [s.Rational(1,6),s.Rational(1,3),s.Rational(1,2)]]
purities=[s.factor(sum(w*w/d for w,d in zip(ws,dims))) for ws in weights]
assert len(set(purities))==3 and all(sum(ws)==1 for ws in weights)
for ws in weights:
    assert sum(d*(w/d) for w,d in zip(ws,dims))==1
    # H_d=-log(w_d/d), hence exp(-H_d)=w_d/d and Z=1 exactly.
    assert all(s.simplify(s.exp(s.log(w/d))-w/d)==0 for w,d in zip(ws,dims))
q=s.symbols('q', integer=True)
quarter_controls=[]
for alpha in [s.Rational(1,4),s.Rational(3,4)]:
    F=q*(q-1)/2+alpha*q
    h=(q+alpha-s.Rational(1,2))**2/2
    assert s.expand(h-s.Rational(1,32)-F)==0
    assert s.expand(F.subs(q,q+1)-F)==q+alpha
    assert s.expand(F.subs(q,q+2)-F)==2*q+1+2*alpha
    assert F.subs(q,0)==0
    assert F.subs(q,1)==alpha and F.subs(q,-1)==1-alpha
    assert all(F.subs(q,j)>0 for j in range(-5,6) if j)
    quarter_controls.append({'alpha':str(alpha),'h0':'1/32',
        'charge_energies_minus3_to3':[str(F.subs(q,j)) for j in range(-3,4)],
        'winding_increment':str(q+alpha),'even_winding_increment':str(2*q+1+2*alpha)})
# These finite diagnostics do not replace the all-charge odd-square proof.
assert 2*s.Rational(1,32)==s.Rational(1,16)
assert 4*s.Rational(1,32)+s.Rational(1,2)==s.Rational(5,8)
results={'native_v1_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'root_count':len(roots),'root_second_moment':'60 I8',
         'highest_root_comarks':[int(x) for x in highest],
         'dual_coxeter':30,'sugawara_c':8,'dynamics_controls':table,
         'cartan_word_norms_squared':[1,2],
         'native_block_state_purities':[str(x) for x in purities],
         'actual_quarter_source_controls':quarter_controls,
         'quarter_all_charge_stress_cyclicity':'analytic odd-square/Kac proof, not finite scan',
         'quarter_three_intertwiner_uniqueness':'analytic dense-core proof, not physical P1 evaluation',
         'all_history_theorem':'analytic proof; not inferred from finite checks',
         'full_P1_countermodels':False,'original_transfer_identity_computed':False,
         'F1_complete':False}
args.output.write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
