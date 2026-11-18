"""Small independent exact checks; not an original physical source derivation."""
import json
from pathlib import Path
import sympy as s

out = {}
rho = (1 + s.I*s.sqrt(3))/2
omega = s.simplify(rho**2)
t = s.symbols('t')
P = s.Matrix(6, 6, lambda i,j: int((j-i) % 6 == 1))
D = s.diag(*[(-1)**j for j in range(6)])
assert D*(-P)*D == P
assert P.charpoly(t).as_expr() == t**6-1
assert (-P).charpoly(t).as_expr() == t**6-1
proj = sum((s.simplify(rho**(-j))*P**j for j in range(6)), s.zeros(6))/6
assert (proj*proj-proj).applyfunc(s.simplify) == s.zeros(6)
assert (P*proj-rho*proj).applyfunc(s.simplify) == s.zeros(6)
assert s.simplify(s.trace(proj)) == 1
assert s.simplify(s.trace(proj*P)-rho) == 0
sigma, sheet = P**4, P**3
triality_projector = (s.eye(6)+omega*sigma+omega**2*sigma**2)/3
joint_q = (s.eye(6)-sheet)*triality_projector/2
joint_l = (s.eye(6)+sheet)*triality_projector/2
assert (joint_q-proj).applyfunc(s.simplify) == s.zeros(6)
assert s.simplify(s.trace(joint_l)) == 1
assert (P*joint_l-rho**4*joint_l).applyfunc(s.simplify) == s.zeros(6)
out['six_cycle'] = {'charpoly':str(t**6-1), 'rho_projector_rank':1,
                    'rho_trace':str(rho), 'alternating_car_intertwiner':True,
                    'full_trace':str(s.trace(P)),
                    'old_triality_and_sheet_joint_rank': [1,1],
                    'old_triality_and_sheet_eigenvalues': ['rho','rho^4'],
                    'all_intertwiners': 'D times the six-dimensional commutant of P; D is unique only in the diagonal orbit-line comparison'}

r = P.copy()
r[5,0] = -1
assert r**6 == -s.eye(6)
assert r.charpoly(t).as_expr() == t**6+1
J, sig = r**9, r**4
assert J**2 == -s.eye(6) and sig**3 == s.eye(6)
assert J*sig == r
q = [0,2,4]
sigq = sig.extract(q,q)
assert s.trace(sigq) == 0 and sigq.charpoly(t).as_expr() == t**3-1
assert s.simplify(2*omega**2+omega) != 0
out['spin_clock'] = {'charpoly':str(t**6+1),'sixth_power': '-I',
                     'geometric_sigma_q_trace':0,
                     'native_sigma_q_trace':str(s.simplify(2*omega**2+omega)),
                     'native_literal_clock_match':False}

Yu=s.diag(1,2,3)
Yd=s.Matrix([[1,1,0],[0,1,1],[1,0,rho]])
Hu,Hd=Yu*Yu.H,Yd*Yd.H
comm=Hu*Hd-Hd*Hu
cp=s.simplify(s.im(s.trace(comm**3)))
assert cp == -360*s.sqrt(3)
out['relative_phase_counterexample']={'Im_trace_commutator_cubed':str(cp),
                                      'scope':'algebraic counterexample to universal CP inference; no TFPT Yukawa prediction'}

Rh=s.diag(rho,omega)
Q=-s.I*s.Matrix([[0,1],[1,0]])
Ph=s.diag(1,0)
assert Q*Ph*Q.inv() == s.eye(2)-Ph
assert s.simplify(s.trace(Ph*Rh*Ph)-rho)==0
assert s.simplify(Rh.det()) == -1
w=s.symbols('w')
assert s.cancel((1-w**6)/(1+w**3)) == 1-w**3
assert s.expand(((1+w**3)+(1-w**3))/2) == 1
out['hodge']={'quotient_projection_trace':str(rho),'whole_hodge_determinant':-1,
               'Q_exchanges_hodge_lines':True,'spin_square_sum_identity':True}
out['physical_source_completed']=False
Path(__file__).with_name('ROOT-MIXED-SPIN-U6-CHECK.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(out,ensure_ascii=False))
