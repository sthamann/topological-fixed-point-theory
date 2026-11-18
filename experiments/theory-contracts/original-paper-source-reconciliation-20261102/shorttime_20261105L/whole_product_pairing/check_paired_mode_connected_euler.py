"""Exact connected global Euclidean Euler limit of the actual paired polynomial.
Assumptions: original four-real-H isotropic IR exp(-r*u/2), held rest with
dominated degree-zero ray limit and nonzero leading angular integral.
No numerical p/Y, no independent state, no fixed-time Yham identification.
"""
from pathlib import Path
import json,sympy as s
p2,Sigma=s.symbols('p_squared Sigma',real=True)
d=4;m=192;n=96
moment=lambda k:2**k*s.rf(s.Rational(d,2),k)
checks=[]
def need(name,b):
 checks.append({'name':name,'passed':bool(b)})
 if not b:raise AssertionError(name)
ratio=s.factor(moment(n-1)/moment(n))
need('whole_paired_Hnorm_moment_ratio',ratio==s.Rational(1,194))
Au_A=s.factor(moment(n+1)/moment(n))
need('whole_paired_leading_U_mark_ratio',Au_A==196)
B_A=s.factor(4*p2*Sigma*ratio)
need('actual_subleading_integrated_ratio',s.expand(B_A-2*p2*Sigma/97)==0)
Bu_A=s.factor(4*p2*Sigma)
need('actual_subleading_U_mark_ratio',s.expand(Bu_A-194*B_A)==0)
connected=s.factor(2*(Bu_A-Au_A*B_A))
need('connected_leading_Euler_scalar_cancels_before_limit',s.expand(connected+4*B_A)==0)
need('fully_evaluated_whole_paired_connected_Euler',s.expand(connected+8*p2*Sigma/97)==0)
need('zero_momentum_returns_exact_mass_homogeneity_zero',connected.subs(p2,0)==0)
out={'passed':sum(c['passed']for c in checks),'total':len(checks),'checks':checks,
 'source_polynomial':'Genuine384 skew opposite-mode Pf = C_Y²*u^96+4*p²*C_Y²*Sigma*u^95+O(u^94)',
 'actual_four_real_H_IR':'exp(-r*u/2), r down to0; same original coordinates and all held weights kept',
 'held_rest_premise':'dominated degree-zero ray limit w_infinity with leading angular integral nonzero; rest independent of composite Euler source g',
 'moment_ratio':str(ratio),'B_over_A':str(B_A),'connected_global_Euclidean_Y_limit':str(connected),
 'scope':'Global affine fermionic Euclidean action Euler mark of the actual free Spin4 constant-H opposite-momentum source; not fixed-time R113 Yham, whole mode tower, full history, coupled gravity or uncut primitive source.',
 'not_R048_positive_quartic_ray_class':True,'native_physical_parameters_chosen':False}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
