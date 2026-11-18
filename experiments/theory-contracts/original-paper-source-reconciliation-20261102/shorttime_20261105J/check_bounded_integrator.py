"""Pure algorithm controls for the R048 formula family, NOT TFPT native data."""
from pathlib import Path
import sys,json,math
import numpy as np
from scipy.special import roots_hermitenorm
import mpmath as mp
sys.path.insert(0,str(Path(__file__).resolve().parent/'shorttime20261105J_lean'))
from integrate_original_radial_tilt import integrate_original_integral
out=[]
# Numerical parameters here test the integrator; they are not selected physical inputs.
d={'original_quadratic_precision':[[100.]],'original_quadratic_mean':[.15],'original_cell_volumes':[1.],'original_lambda':.1,'original_vev':.5,'absolute_error_for_bounded_expectation':.02,'algorithm_max_evaluations':200000,'u_indices':[0]}
r=integrate_original_integral(d,execute=True)
mp.mp.dps=45;a=mp.sqrt(2)*mp.mpf('.5');beta=mp.mpf('.1');mean=mp.mpf('.15');sigma=mp.mpf('.1')
def ff(z):
 eta=mean+sigma*z
 W=beta/4*(eta+3*a)**2*((eta-a)**2+2*a*a)
 return mp.exp(-W-z*z/2)/mp.sqrt(2*mp.pi)
zz=mp.quad(ff,[-mp.inf,-8,0,8,mp.inf]);uu=mp.quad(lambda z:((mean+sigma*z+a)**2/2)*ff(z),[-mp.inf,-8,0,8,mp.inf])
for name,val in [('Z',zz),('u[0]',uu)]:
 q=complex(*r['unnormalized_original_insertions'][name]['unnormalized_bounded_estimate']);err=r['analytic_absolute_errors'][name];e=abs(q-complex(val));assert e<=err
 out.append({'case':'one-dimensional original-formula algorithm control','insertion':name,'evaluations':r['evaluations'],'absolute_error_vs_independent_45digit_quad':float(e),'analytic_quadrature_bound_before_rounding':err,'passed':True})
# Non-diagonal whole covariance: no independent mode simplification.
d2={**d,'original_quadratic_precision':[[100.,15.],[15.,80.]],'original_quadratic_mean':[.15,-.1],'original_cell_volumes':[1.,.8],'absolute_error_for_bounded_expectation':.2,'u_indices':[],'algorithm_max_evaluations':250000}
r2=integrate_original_integral(d2,execute=True)
cov=np.linalg.inv(np.array(d2['original_quadratic_precision']));L=np.linalg.cholesky(cov);m=np.array(d2['original_quadratic_mean']);aa=math.sqrt(2)*.5;bb=.1*np.array([1.,.8])
refs=[]
for order in [48,80]:
 nodes,weights=roots_hermitenorm(order);z=np.stack(np.meshgrid(nodes,nodes,indexing='ij'),axis=-1).reshape(-1,2);eta=m+z@L.T
 W=np.sum(bb/4*(eta+3*aa)**2*((eta-aa)**2+2*aa*aa),axis=1);w=np.outer(weights,weights).ravel()/(2*math.pi);refs.append(float(np.dot(w,np.exp(-W))))
err=abs(r2['bounded_expectation_midpoint']-refs[-1]);assert err<=r2['analytic_absolute_errors']['Z'];assert abs(refs[0]-refs[1])<1e-12
out.append({'case':'correlated two-dimensional algorithm control','evaluations':r2['evaluations'],'absolute_error_vs_independent_Hermite80':err,'Hermite48_80_difference':abs(refs[0]-refs[1]),'analytic_quadrature_bound_before_rounding':r2['analytic_absolute_errors']['Z'],'passed':True})
# Whole zero denominator must remain unnormalized and cannot cause division.
d0={**d,'u_indices':[],'original_polynomial_insertions':{'zero':[]},'original_denominator_polynomial':'zero','normalize_insertions':['zero'],'absolute_error_for_bounded_expectation':.2}
r0=integrate_original_integral(d0,execute=True);assert not r0['denominator_disc_excludes_zero'] and not r0['normalized_original_insertions'];out.append({'case':'exact zero denominator preserved, no division','passed':True})
from integrate_original_radial_tilt import integrate_until_original_tolerance
ad={**d,'algorithm_max_evaluations':1000000,'normalized_absolute_tolerance':.005,'normalized_relative_tolerance':0.,'algorithm_max_refinements':6}
ar=integrate_until_original_tolerance(ad)
aq=ar['normalized_original_insertions']['u[0]'];actual=abs(complex(*aq['normalized_estimate'])-complex(uu/zz));assert actual<=aq['analytic_ratio_error']<=.005
assert len(ar['refinement_history'])>1 and ar['rounding_enclosed'] is False
out.append({'case':'adaptive neutral moment tolerance, unchanged algorithm-control inputs','rungs':len(ar['refinement_history']),'final_evaluations':ar['evaluations'],'absolute_error_vs_independent_quad':actual,'analytic_ratio_error_before_rounding':aq['analytic_ratio_error'],'passed':True})
bad={**d,'normalize_insertions':['UY_not_supplied']}
try:
 integrate_original_integral(bad,execute=False)
except ValueError as exc:
 assert 'no supplied original polynomial' in str(exc)
 out.append({'case':'requested but missing UY insertion rejected before integration','passed':True})
else:raise AssertionError('Missing source insertion must not report tolerance success')
res={'scope':'Algorithm controls only, NOT supplied native TFPT field data and NOT physical predictions','rounding_enclosed':False,'checks':out,'passed':len(out),'total':len(out),'native_primitive_source_value_evaluated':False}
Path('work/shorttime20261105J_gravity/integrator-math-controls.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res))
