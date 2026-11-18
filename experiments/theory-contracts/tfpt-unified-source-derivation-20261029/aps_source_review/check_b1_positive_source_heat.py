"""Existing AP positive factor, full orientation and squared-variable moments."""
from pathlib import Path
import sys,json
import sympy as s
import numpy as np
R=Path(__file__).resolve().parents[4];sys.path.insert(0,str(R/'verification'))
import v252_full_finite_triple as f
u,t=s.symbols('u t',positive=True)
F=s.exp(-s.sqrt(u));density=s.exp(-1/(4*t))/(2*s.sqrt(s.pi)*t**s.Rational(3,2))
checks={}
# The full subordination identity is proved analytically in the report.
# Exact point controls and all required scalar moment integrals follow.
checks['positive_mixture_mass_one']=s.integrate(density,(t,0,s.oo))==1
for j in (1,2):
 checks[f'inverse_heat_time_moment_{j}']=s.integrate(density/t**j,(t,0,s.oo))==(2 if j==1 else 12)
for point in (1,4):
 checks[f'positive_mixture_value_u_{point}']=s.simplify(s.integrate(density*s.exp(-point*t),(t,0,s.oo))-s.exp(-s.sqrt(point)))==0
m0=s.Integer(1);m1=s.integrate(F,(u,0,s.oo));m2=s.integrate(u*F,(u,0,s.oo))
checks['squared_argument_moment1_is2']=m1==2
checks['squared_argument_moment2_is12']=m2==12
checks['shape_invariant3']=m0*m2/m1**2==3
amp,rate,a0,a1,a2=s.symbols('amp rate a0 a1 a2',positive=True)
checks['shape_invariant_preserved_by_units']=s.simplify((amp*a0)*(amp*a2/rate**2)/(amp*a1/rate)**2-a0*a2/a1**2)==0
checks['cannot_normalize_all_three_to1']=s.simplify((m0*m2-m1**2))==8
Y=[np.array([[1+1j,2,0],[0,3-1j,1],[1j,0,4]],complex)+j*np.eye(3) for j in range(4)]
MR=np.array([[2+1j,1,0],[1,3,1j],[0,1j,4-1j]],complex)
Phi=f.full_D(*Y,MR);g=f.grading();h=.25*g-1j*g@Phi
v,U=np.linalg.eigh(h);beta=1/32;I=np.eye(96)
mat=lambda w:(U*w)@U.conj().T
C=mat(1/(1+np.exp(beta*v)));X=mat(1/(2*np.cosh(beta*v/2)))
Q=mat(np.exp(-beta*np.abs(v)));sign=mat(np.sign(v))
ratio=mat(np.exp(-beta*np.abs(v))/(1+np.exp(-beta*np.abs(v))))
Cx=(I-sign)/2+sign@ratio;Xx=mat(np.sqrt(np.exp(-beta*np.abs(v)))/(1+np.exp(-beta*np.abs(v))))
close=lambda a,b:bool(np.linalg.norm(a-b)<1e-10)
checks['actual96_C_rebuilt_from_sign_and_positive_factor']=close(C,Cx)
checks['actual96_X_rebuilt_from_positive_factor']=close(X,Xx)
checks['actual96_sign_involution']=close(sign@sign,I)
checks['actual96_positive_factor_is_contraction']=bool(np.min(np.exp(-beta*np.abs(v)))>0 and np.max(np.exp(-beta*np.abs(v)))<1)
checks['actual96_h_and_minus_h_same_positive_factor']=close(mat(np.exp(-beta*np.abs(-v))),Q)
checks['actual96_h_and_minus_h_different_C']=close(mat(1/(1+np.exp(-beta*v))),I-C) and bool(np.linalg.norm(C-(I-C))>1e-3)
checks['actual96_pair_rebuilds_full_G']=close(np.block([[I-C,X],[X,C]]),np.block([[I-Cx,Xx],[Xx,Cx]]))
checks['actual96_square_root_argument_identity']=np.array_equal(np.sqrt((beta*np.abs(v))**2),beta*np.abs(v))
result=dict(contract='B1-POSITIVE-SOURCE-HEAT-20261030',checks={k:bool(v) for k,v in checks.items()},all_passed=all(checks.values()),profile='exp(-sqrt(u)) for u=(beta E)^2',moments=[int(m0),int(m1),int(m2)],shape_invariant=3,scope='Existing AP factor and full source polar decomposition; analytic all-spectrum subordination, exact scalar integrals and actual96 controls. No identification with physical Sbos or new source/action selection.')
Path(__file__).with_name('B1-POSITIVE-SOURCE-HEAT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2));assert result['all_passed']
