"""Exact action-level and contour tests for the existing scalaron/Higgs sector.
No quantum measure, real-rho contour or positive heat generator is selected.
"""
from pathlib import Path
import json
import sympy as s

rho,rho0,alpha,ZH,z=s.symbols('rho rho0 alpha ZH z',positive=True)
R,F,xi,phi,M=s.symbols('R F xi phi M',real=True)
checks=[]
def need(name,ok):
 checks.append({'name':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
aux=rho*R-(rho-F)**2/(4*alpha)
saddle=F+2*alpha*R
need('actual_R2_Legendre_stationarity',s.simplify(s.diff(aux,rho).subs(rho,saddle))==0)
need('actual_R2_classical_elimination',s.simplify(aux.subs(rho,saddle)-(F*R+alpha*R**2))==0)
need('Euclidean_vertical_contour_integrand',s.simplify((-aux).subs(rho,F+s.I*xi)-(-F*R-xi**2/(4*alpha)-s.I*xi*R))==0)
need('Euclidean_HS_square_completion',s.simplify(-xi**2/(4*alpha)-s.I*xi*R+(xi+2*s.I*alpha*R)**2/(4*alpha)+alpha*R**2)==0)
need('real_auxiliary_integral_has_wrong_Gaussian_sign',s.diff(-aux,rho,2)==1/(2*alpha))
# R_g = Omega²[R_E + 6 box_E lnOmega - 6(d lnOmega)²].
# The total derivative has constant rho0 coefficient after Omega²=rho/rho0.
weyl_rho_kinetic=s.simplify(-6*rho0*s.diff(s.log(rho/rho0)/2,rho)**2)
need('raw_plus_Euclidean_R_gives_negative_scalar_kinetic',weyl_rho_kinetic==-3*rho0/(2*rho**2))
Ghealthy=s.diag(3*rho0/rho**2,*([2*ZH*rho0/rho]*4))
Geuclid=s.diag(-3*rho0/rho**2,*([2*ZH*rho0/rho]*4))
J=s.diag(2*z,*([s.sqrt(6/ZH)]*4))
want=12*rho0/z**2*s.eye(5)
need('healthy_classical_H5_metric_pullback',s.simplify(J.T*Ghealthy.subs(rho,z*z)*J-want)==s.zeros(5))
need('raw_Euclidean_metric_stays_indefinite',s.simplify(J.T*Geuclid.subs(rho,z*z)*J-12*rho0/z**2*s.diag(-1,1,1,1,1))==s.zeros(5))
volume=4*s.sqrt(3)*ZH**2*rho0**s.Rational(5,2)/rho**3
need('target_metric_volume_only',s.simplify(volume**2-Ghealthy.det())==0)
need('target_volume_coordinate_change',s.simplify(volume.subs(rho,z*z)*J.det()-(12*rho0)**s.Rational(5,2)/z**5)==0)
V=rho0**2/(4*alpha)*(1-s.exp(-phi/s.sqrt(3*rho0)))**2
mass=s.simplify(s.diff(V,phi,2).subs(phi,0))
need('classical_scalaron_mass_with_actual_alpha',mass==rho0/(6*alpha))
need('Starobinsky_M_normalization_requires_alpha_rho0_over6M2',s.simplify(mass.subs(alpha,rho0/(6*M**2))-M**2)==0)
out={'scope':'Existing homogeneous residue/constant-chi0 action, symbolic classical Legendre/Weyl map. Healthy Lorentz sigma target differs from raw printed Euclidean plus-R target. No source contour or measure equivalence, no primitive physical kernel derived.',
 'checked':sum(c['passed']for c in checks),'total':len(checks),'checks':checks,
 'healthy_target_metric':'3 rho0 d(rho)^2/rho^2 + 2 ZH rho0 sum_a dH_a^2/rho',
 'raw_Euclidean_target_signature':'(-,+,+,+,+)',
 'H5_radius_squared':'12 rho0 = 6 Mpl^2',
 'H5_scalar_curvature':'-20/(12 rho0) = -5/(3 rho0)',
 'hypothetical_Laplace_Beltrami':'[z^2(partial_z^2+sum_a partial_ya^2)-3z partial_z]/(12rho0); NOT chosen as original quantum generator',
 'exact_Euclidean_auxiliary_contour':'rho=F+i xi, xi real, not independent real positive rho',
 'primitive_quantum_kernel_computed':False,'original_measure_replaced':False}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
