"""All-parameter DtN -> full heat trace on the actual quarter circle/full96.
Only unchanged R055 diagnostic Yukawas. Proper time is not physical source time.
"""
from pathlib import Path
import sys,json
import numpy as np
from scipy.integrate import quad_vec
ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(ROOT/'verification'))
import v252_full_finite_triple as v
rng=np.random.default_rng(252)
Ys=[rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for _ in range(4)]
mr=rng.normal(size=(3,3));MR=mr+mr.T
Phi=v.full_D(*Ys,MR)  # Actual Hstar full block, no dropped Majorana.
lam=np.linalg.eigvalsh(Phi@Phi)
L=2*np.pi;a=L/4;tau=1.
checks=[];rows=[];heatrows=[]
def need(name,ok):
 checks.append({'name':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
for alpha in [.25,.75]:
 U=np.zeros((4,4),complex)
 for i in range(3):U[i,i+1]=1
 U[3,0]=np.exp(2j*np.pi*alpha)
 need('actual_quarter_holonomy_'+str(alpha),np.max(abs(np.linalg.matrix_power(U,4)-np.exp(2j*np.pi*alpha)*np.eye(4)))<1e-14)
 for z in [0.,1.,2.+.7j]:
  kappas=np.sqrt(lam+z+0j)
  trace_dtn=0j;log_det_error=0.
  for k in kappas:
   sh=np.sinh(k*a);ch=np.cosh(k*a)
   R=k*(2*ch/sh*np.eye(4)-(U+U.conj().T)/sh)
   Rk=2*(ch/sh-k*a/(sh*sh))*np.eye(4) - (1/sh-k*a*ch/(sh*sh))*(U+U.conj().T)
   trace_dtn+=np.trace(np.linalg.solve(R,Rk/(2*k)))
   detR=np.linalg.det(R)
   gluing=(2*np.sinh(k*a)/k)**4*detR/16
   expected=2*np.cosh(k*L)  # cos(2pi alpha)=0 exactly for both original quarters.
   log_det_error=max(log_det_error,float(abs(gluing/expected-1)))
  circle=(L/(2*kappas)*np.tanh(L*kappas)).sum()
  cut=(4*(a/(2*kappas)/np.tanh(a*kappas)-1/(2*kappas*kappas))).sum()
  error=abs(trace_dtn-(circle-cut))
  need('native96_quarter_BFK_all_modes_a'+str(alpha)+'_z'+str(z),log_det_error<1e-9)
  need('native96_quarter_resolvent_gluing_a'+str(alpha)+'_z'+str(z),error<1e-8)
  rows.append({'alpha':alpha,'z':[float(np.real(z)),float(np.imag(z))],
   'full96_logdet_relative_error':log_det_error,'resolvent_absolute_error':float(error),
   'dtn_resolvent_trace':[float(trace_dtn.real),float(trace_dtn.imag)]})
 # Independent original complete circle spectrum and four Dirichlet pieces.
 N=16;n=np.arange(-N,N+1);k=np.arange(1,N+1)
 internal=float(np.exp(-tau*lam).sum())
 heat_circle=internal*float(np.exp(-tau*(n+alpha)**2).sum())
 heat_cut=4*internal*float(np.exp(-tau*(np.pi*k/a)**2).sum())
 circle_tail=96*2*np.exp(-tau*(N+1-alpha)**2)/(1-np.exp(-tau*(2*(N+1-alpha)+1)))
 # Whole DtN contour reconstruction, independent of truncating circle modes.
 left=float((lam[0]+min(alpha,1-alpha)**2)/2)
 right=float(lam[-1]+(N+1)**2+.123);height=2.
 corners=[left-1j*height,right-1j*height,right+1j*height,left+1j*height]
 def boundary_resolvent(w):
  kk=np.sqrt(lam-w+0j)
  return -(2/(kk*kk)+L/(2*kk)*np.tanh(L*kk)-2*a/kk/np.tanh(a*kk)).sum()
 contour=0j;estimates=[]
 for j in range(4):
  start=corners[j];step=corners[(j+1)%4]-start
  val,err=quad_vec(lambda u:np.exp(-tau*(start+u*step))*boundary_resolvent(start+u*step)*step,
                   0.,1.,epsabs=1e-9,epsrel=1e-10)
  contour+=val;estimates.append(float(err))
 boundary_heat=contour/(2j*np.pi)
 heat_error=float(abs(boundary_heat-(heat_circle-heat_cut)))
 need('actual_whole96_circle_heat_from_entire_DtN_'+str(alpha),heat_error<1e-7)
 reference_heat=96*float(np.exp(-tau*(n+alpha)**2).sum())
 # Full pointwise kernel across two distinct original quarter arcs, not a trace.
 x=.2;y=np.pi+.4
 def point_boundary_green(w):
  kk=np.sqrt(lam-w+0j);sh=np.sinh(kk*a);ch=np.cosh(kk*a)
  R=kk[:,None,None]*(2*(ch/sh)[:,None,None]*np.eye(4)-(U+U.conj().T)[None,:,:]/sh[:,None,None])
  px=np.zeros((96,4),complex);py=np.zeros((96,4),complex)
  px[:,0]=np.sinh(kk*(a-x))/sh;px[:,1]=np.sinh(kk*x)/sh
  yy=y-2*a
  py[:,2]=np.sinh(kk*(a-yy))/sh;py[:,3]=np.sinh(kk*yy)/sh
  return np.einsum('ni,ni->n',px,np.linalg.solve(R,py[...,None])[...,0])
 # Minus sign changes (K-w)^-1 to the Cauchy convention (w-K)^-1.
 point_contour=np.zeros(96,complex);point_estimates=[]
 for j in range(4):
  start=corners[j];step=corners[(j+1)%4]-start
  value,error=quad_vec(lambda u:-np.exp(-tau*(start+u*step))*point_boundary_green(start+u*step)*step,
                       0.,1.,epsabs=1e-10,epsrel=1e-10)
  point_contour+=value;point_estimates.append(float(error))
 point_boundary=point_contour/(2j*np.pi)
 scalar_circle=np.sum(np.exp(-tau*(n+alpha)**2+1j*(n+alpha)*(x-y)))/L
 point_native=np.exp(-tau*lam)*scalar_circle
 point_error=float(np.linalg.norm(point_boundary-point_native))
 need('whole96_cross_quarter_pointwise_heat_kernel_'+str(alpha),point_error<1e-8)
 heatrows.append({'alpha':alpha,'proper_time_diagnostic':tau,'physical_source_time_identified':False,
  'full96_circle_heat':heat_circle,'four_Dirichlet_piece_heat':heat_cut,
  'boundary_reconstructed_heat_difference':[float(boundary_heat.real),float(boundary_heat.imag)],
  'heat_absolute_error':heat_error,'circle_mode_tail_bound':float(circle_tail),
  'quadrature_estimates':estimates,'original_S_reference_heat':reference_heat,
  'relative_S_heat':heat_circle-reference_heat,
  'pointwise_cross_quarter_kernel':{'x':x,'y':y,'native_trace':[float(point_native.sum().real),float(point_native.sum().imag)],
   'boundary_trace':[float(point_boundary.sum().real),float(point_boundary.sum().imag)],
   'whole96_Frobenius_error':point_error,'quadrature_estimates':point_estimates}})
out={'scope':'Actual fixed circle L2pi, original quarters1/4,3/4, full96 Hstar v252 with unchanged R055 seed252. No spatial/proper/source-time identification, new physical Yukawas, chosen profile/measure or primitive H integral. Four geometric cut points are the existing quarter marks.',
 'checked':sum(c['passed']for c in checks),'total':len(checks),'checks':checks,
 'all_z_controls':rows,'whole_circle_heat_controls':heatrows,
 'eta_or_pfaffian_phase_replaced':False,'primitive_common_integral_evaluated':False}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k]for k in ['checked','total','whole_circle_heat_controls']},indent=2))
