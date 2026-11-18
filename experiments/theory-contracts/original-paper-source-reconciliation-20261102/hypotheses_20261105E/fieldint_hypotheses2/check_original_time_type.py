"""Exact kinematic test of the literal unchanged quarter-circle -> line candidate.
No new source, state, parameter search, or physical field assignment.
"""
from pathlib import Path
import sympy as s
import json
n=s.symbols('n',integer=True);x=s.symbols('x',real=True)
checks=[]
def need(name,ok):
 checks.append({'name':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
T=8*s.pi
phases=[]
for alpha in [s.Rational(1,4),s.Rational(3,4)]:
 need(f'original_Balpha_period_all_n_{alpha}',s.simplify(s.exp(s.I*T*(n+alpha)))==1)
 phase=s.simplify(s.exp(s.I*2*s.pi*(n+alpha)))
 need(f'original_Balpha_2pi_phase_{alpha}',phase==s.exp(2*s.pi*s.I*alpha).expand(complex=True))
 phases.append({'alpha':str(alpha),'2pi_one_particle_phase':str(phase),'8pi_one_particle_phase':'1'})
 for gamma in [-1,1]:
  need(f'free_original_Balpha_gamma_period_all_n_a{alpha}_g{gamma}',s.simplify(s.exp(s.I*T*gamma*(n+alpha)))==1)
need('faithful_line_translation_at_same_original_period_is_not_identity',s.simplify((x+T)-x)!=0)
# PSL2 conjugation preserves determinant; real nonzero rescaling preserves sign.
J=s.Matrix([[0,-1],[1,0]])
P=s.Matrix([[0,1],[0,0]])
D=s.diag(s.Rational(1,2),-s.Rational(1,2))
need('elliptic_rotation_positive_det',J.det()==1)
need('parabolic_translation_zero_det',P.det()==0)
need('hyperbolic_dilation_negative_det',D.det()==-s.Rational(1,4))
need('rotation_generator_square_negative',J*J==-s.eye(2))
need('translation_generator_nilpotent',P*P==s.zeros(2))
need('dilation_generator_square_positive',D*D==s.eye(2)/4)
# Cayley/stereographic chart transforms circle rotation to (1+x²)/2 dx,
# never to dx or x dx by this fixed chart alone.
theta=s.symbols('theta',real=True)
need('circle_rotation_stereographic_vector_is_elliptic',s.simplify(s.diff(s.tan(theta/2),theta)-(1+s.tan(theta/2)**2)/2)==0)
out={'scope':'Literal original free periodic quarter-circle flow only. Original general P1 tau_t and nonzero-Phi coupled circle Hamiltonian are not assumed periodic or identified with this free flow. No global holography exclusion.',
 'checked':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks,
 'period':'8*pi','phases':phases,
 'psl2_types':{'rotation':'elliptic, det positive','line_translation':'parabolic, det zero','modular_halfline_dilation':'hyperbolic, det negative'},
 'candidate_verdict':'Literal unchanged globally periodic circle flow cannot be faithfully intertwined with line translation, nor conjugated to dilation by a global time-preserving change of variables.',
 'general_P1_time_classification':'Original operational seed supplies tau_t as data; its equality to a PSL2 generator is not established by inspected source.',
 'primitive_source_evaluated':False}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['checked','total','period','phases','candidate_verdict','scope']},indent=2))
