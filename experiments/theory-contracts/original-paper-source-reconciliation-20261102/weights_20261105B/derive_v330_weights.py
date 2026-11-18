"""Exact original v330 spectral/history weights; no primitive S_adm claim."""
from pathlib import Path
import importlib.util
import json
import sys

import numpy as np
import sympy as s

ROOT = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0, str(ROOT/'verification'))
spec = importlib.util.spec_from_file_location('original_v330', ROOT/'verification/v330_qgamb_admissible_measure.py')
v330 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v330)
T_actual, omega_actual, r_actual, s_actual = v330.symmetric_transfer()

r, q = s.Rational(64,729), s.Rational(1,729)
u0 = s.ones(3,1)/s.sqrt(3)
u1 = s.Matrix([1,-1,0])/s.sqrt(2)
u2 = s.Matrix([1,1,-2])/s.sqrt(6)
P = [u*u.T for u in [u0,u1,u2]]
T = P[0] + r*P[1] + q*P[2]
a,b,c = s.symbols('a b c', real=True)
D = s.diag(a,b,c)
weights = [s.factor((u0.T*D*pj*D*u0)[0]) for pj in P]
expected = [(a+b+c)**2/9,(a-b)**2/6,(a+b-2*c)**2/18]

checks = {}
checks['projectors_sum_identity'] = sum(P,s.zeros(3)) == s.eye(3)
checks['projectors_are_orthogonal'] = all(P[i]*P[j] == (P[i] if i==j else s.zeros(3)) for i in range(3) for j in range(3))
checks['spectral_weights_exact_all_diagonal_insertions'] = all(s.expand(x-y)==0 for x,y in zip(weights,expected))
checks['actual_v330_constructor_matches_exact_matrix'] = bool(np.max(np.abs(T_actual-np.array(T,dtype=float))) < 2e-16)
checks['actual_v330_perron_vector_matches'] = bool(np.max(np.abs(omega_actual-np.array(u0,dtype=float).ravel())) < 2e-16)
checks['actual_v330_rates_match'] = abs(float(r)-r_actual)<1e-16 and abs(float(q)-s_actual)<1e-16
checks['strictly_positive_transition_entries'] = all(x>0 for x in T)
checks['rows_sum_one'] = T*s.ones(3,1)==s.ones(3,1)
checks['stationary_label_weights_one_third'] = u0.applyfunc(lambda x:x*x)==s.ones(3,1)/3
checks['all_moments_symbolic_0_to_7'] = all(s.simplify((u0.T*D*T**n*D*u0)[0] - (weights[0]+weights[1]*r**n+weights[2]*q**n))==0 for n in range(8))
checks['empty_mean_subtraction_keeps_two_connected_weights'] = s.simplify(weights[0]-(s.trace(D)/3)**2)==0

# The original v330 correlator is independently called on one deterministic
# admissible diagonal test; exact general coefficients above cover all a,b,c.
Dn=np.diag([1.,2.,-3.]); substitution={a:1,b:2,c:-3}
checks['original_correlator_matches_closed_connected_formula'] = all(abs(v330.correlator(T_actual,omega_actual,Dn,n)-float((weights[1]*r**n+weights[2]*q**n).subs(substitution)))<3e-15 for n in range(8))

# Identity shifts in the insertion change only the zero-rate coefficient.
d=s.symbols('d',real=True)
checks['identity_shift_leaves_connected_weights'] = all(s.expand(w.subs({a:a+d,b:b+d,c:c+d},simultaneous=True)-w)==0 for w in weights[1:])

# Label insertions are not transition weights.  V330 has the complete diagonal
# coordinate projectors.  Their finite ordered history law is pi_i0 product T.
E=[s.diag(*[int(i==j) for i in range(3)]) for j in range(3)]
checks['coordinate_diagonals_separate_original_three_basis_states'] = all(sum(E,s.zeros(3))==s.eye(3) for _ in [0]) and all(s.trace(e)==1 and e*e==e for e in E)
checks['whole_two_step_label_law'] = all(s.simplify((u0.T*E[i]*T*E[j]*T*E[k]*u0)[0] - T[i,j]*T[j,k]/3)==0 for i in range(3) for j in range(3) for k in range(3))

assert all(checks.values()), {k:v for k,v in checks.items() if not v}
out={
 'source':'verification/v330_qgamb_admissible_measure.py:symmetric_transfer,correlator',
 'scope':'Entire existing three-state symmetric seam source. Not the primitive S_adm, its field map or 3+1D dynamics.',
 'checks':checks,'passed':sum(checks.values()),'total':len(checks),
 'T_exact':[[str(x) for x in row] for row in T.tolist()],
 'weights_spectral_A_B_C':[str(x) for x in weights],
 'stationary_coordinate_weights':['1/3']*3,
 'whole_history_law':'omega(E_i0 T^n1 E_i1 ... T^nk E_ik)=(1/3) product_j (T^nj)[i_(j-1),i_j]',
 'not_derived':['No identity with primitive S_adm','No independently chosen new state or measure','No arbitrary diagonal label weight interpreted as an interaction coefficient'],
}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
