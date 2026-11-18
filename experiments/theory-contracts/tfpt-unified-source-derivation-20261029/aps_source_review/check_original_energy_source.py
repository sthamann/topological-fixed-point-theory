"""R042/R043: direct controls on the existing quarter/AP source.
No source selection, geometry fit, new vertex or verification-suite change.
"""
from pathlib import Path
import json
import mpmath as mp
mp.mp.dps=70
I=mp.eye(2);N=mp.matrix([[0,1j],[1j,0]])
B=mp.matrix([[mp.mpf(1)/4,-mp.j/3],[mp.j/3,-mp.mpf(1)/4]])
E=mp.mpf(5)/12;ell=mp.mpf(1);gap=mp.mpf(1)/4
checks={};err={}
def norm(a):return mp.sqrt(sum(abs(x)**2 for x in a))
def need(k,ok):
 checks[k]=bool(ok)
 if not ok:raise AssertionError(k)
def near(k,a,b,tol='1e-60'):
 e=norm(a-b) if isinstance(a,mp.matrix) else abs(a-b)
 err[k]=float(e);need(k,e<mp.mpf(tol))
def block(a,b,c,d):
 o=mp.matrix(4)
 for i in range(2):
  for j in range(2):o[i,j]=a[i,j];o[i,j+2]=b[i,j];o[i+2,j]=c[i,j];o[i+2,j+2]=d[i,j]
 return o
def pair(z):
 A=B+z*N;T=mp.expm(-ell*A);a=(I+T*T)**-1
 C=block(a,T*a,T*a,T*T*a)
 aa=(I+T.H*T)**-1;Q=block(aa,aa*T.H,T*aa,T*aa*T.H)
 return A,T,C,Q
z=mp.mpf(1)/8
A,T,C,Q=pair(z)
kap=mp.sqrt(E*E-z*z);c=mp.cosh(ell*kap);s=mp.sinh(ell*kap)/kap;d=c*c+s*s*z*z
tt=mp.tanh(ell*kap)/kap
Cf=block((I+tt*A)/2,I/(2*c),I/(2*c),(I-tt*A)/2)
Qf=block(I/2+c*s*B/(2*d)+s*s*z*N*B/(2*d),(c*I+s*z*N)/(2*d),(c*I-s*z*N)/(2*d),I/2-c*s*B/(2*d)+s*s*z*N*B/(2*d))
near('raw_trace_closed_formula',C,Cf)
near('orthogonal_trace_closed_formula',Q,Qf)
near('raw_trace_idempotent',C*C,C)
near('orthogonal_trace_idempotent',Q*Q,Q)
near('orthogonal_trace_selfadjoint',Q,Q.H)
need('real_energy_raw_trace_is_not_a_state_projector',norm(C-C.H)>mp.mpf('0.01'))
_,T0,C0,Q0=pair(0);near('zero_energy_raw_equals_orthogonal',C0,Q0)
for zz in [1j*mp.mpf(2)/3,-1j*mp.mpf(2)/3]:
 az,tz,cz,qz=pair(zz)
 near('imaginary_energy_raw_is_orthogonal_'+str(zz.imag),cz,qz)
def Green(z,t):
 az=B+z*N
 return (mp.expm(-t*az) if t>=0 else -mp.expm(-(t+2*ell)*az))*(I+mp.expm(-2*ell*az))**-1
near('original_AP_Green_jump',Green(z,0)-(-mp.expm(-2*ell*A)*(I+mp.expm(-2*ell*A))**-1),I)
near('original_AP_Green_antiperiodicity',Green(z,-ell/3),-Green(z,2*ell-ell/3))
near('Poisson_both_actual_face_limits',block(Green(z,0),-Green(z,-ell),Green(z,ell),T*T*(I+T*T)**-1),C)
a=mp.tanh(ell*E)/E;b=1/mp.cosh(ell*E)
Cp=block(a*N/2,mp.zeros(2),mp.zeros(2),-a*N/2)
Qp=block(a*a*N*B/2,a*b*N/2,-a*b*N/2,a*a*N*B/2)
near('raw_energy_derivative',mp.diff(lambda zz:pair(zz)[2],0),Cp)
near('orthogonal_energy_derivative',mp.diff(lambda zz:pair(zz)[3],0),Qp)
near('orthogonal_derivative_from_raw_commutator',Qp,Cp*C0-C0*Cp)
near('exact_uniform_derivative_norm',Qp*Qp,(a*a/4)*mp.eye(4))
need('all_mode_quarter_bound_not_finite_cutoff',a/2<=min(ell/2,1/(2*gap)))
Tp=mp.diff(lambda zz:pair(zz)[1],0)
Pplus=(I+B/E)/2
moment=mp.exp(-ell*E)*sum((Pplus*N*Tp)[i,i] for i in range(2))
near('whole_Green_trace_gives_existing_flat_moment',moment,(1-mp.exp(-2*ell*E))/(2*E))
near('orthogonal_graph_reads_same_transfer',Q[:2,:2]**-1, (I+T.H*T))
near('orthogonal_graph_transfer_quotient',Q[2:4,:2]*Q[:2,:2]**-1,T)
# Keep the ACTUAL OLD endpoint projector; do not replace it by P_-(B).
P=mp.matrix([[0,0],[0,1]]);Qold=I-P
def fixed(z,r,s0):
 az=B+z*N;tr=mp.expm(-r*az);te=mp.expm(-ell*az)
 pinv=P/te[1,1]
 return (mp.expm(-(r-s0)*az) if r>s0 else mp.zeros(2))-tr*pinv*mp.expm(-(ell-s0)*az)
def incoming(z):
 az=B+z*N;te=mp.expm(-ell*az);pinv=P/te[1,1]
 R=Qold-pinv*te*Qold+pinv
 U=mp.matrix(4,2)
 U[:2,:]=R;U[2:,:]=te*R
 return U
U=incoming(z);Qu=U*(U.H*U)**-1*U.H
near('actual_fixed_APS_two_face_source_recovers_whole_graph',Qu,Q)
r=mp.mpf(2)/5
Gstart=mp.expm(-r*A)-mp.expm(-r*A)*(P/T[1,1])*T
Gend=-mp.expm(-r*A)*(P/T[1,1])
near('actual_fixed_APS_boundary_columns_match_Poisson_source',Gstart*Qold-Gend*P,mp.expm(-r*A)*U[:2,:])
near('actual_fixed_APS_left_domain',Qold*fixed(z,0,mp.mpf(1)/3),mp.zeros(2))
near('actual_fixed_APS_right_domain',P*fixed(z,ell,mp.mpf(1)/3),mp.zeros(2))
# Decisive counterexample to inferring full Tprime from a single cross kernel.
Pad=(I-B/E)/2
def adapted_cross(zz):
 az=B+zz*N;te=mp.expm(-ell*az)
 # M inverse on the one-dimensional negative eigenspace.
 mm=sum((Pad*te)[i,i] for i in range(2))
 return te-te*Pad/mm*te
near('single_APS_cross_derivative_vanishes',mp.diff(adapted_cross,0),mp.zeros(2))
need('whole_transfer_derivative_is_nonzero_despite_cross_zero',norm(Tp)>mp.mpf('0.1'))
# Original AP normal Fourier mode k=pi/(2ell), same quarter/Yukawa block.
k=mp.pi/(2*ell);K0=N*(1j*k*I+B);hh=mp.mpf(2)/3
V=N*mp.matrix([[0,-mp.j/3],[mp.j/3,0]])
near('actual_AP_normal_mode_K_selfadjoint',K0,K0.H)
zz=1j*mp.mpf(3)/7
def Ksigma(e,outer):
 f=mp.exp(-e*hh/2)
 return f*f*(K0-V)+ (V if outer else f*f*V)
for outer in [False,True]:
 KK=Ksigma(0,outer);G=(KK-zz*I)**-1
 dg=mp.diff(lambda e:(Ksigma(e,outer)-zz*I)**-1,0)
 W=hh*V if outer else mp.zeros(2)
 defect=hh*G+zz*hh*G*G-dg
 near('whole_scale_source_reads_actual_W_'+str(outer),defect,G*W*G)
 reconstructed=(KK-zz*I)*defect*(KK-zz*I)
 near('sandwich_readout_uniquely_recovers_W_'+str(outer),reconstructed,W)
# These checks validate actual-source formulas, not an original Brel identity.
result={'contract':'ORIGINAL-ENERGY-SOURCE-20261031','all_passed':all(checks.values()),'checks':checks,'errors':err,'scope':'Existing historical p=1/4,m=1/3 block and actual AP normal Fourier mode, unchanged old endpoint projector. Analytic proof treats all modes. No new Y/MR,physical scale,normal source or Brel map selected. No new Lean formalization or ledger promotion.'}
Path(__file__).with_name('ORIGINAL-ENERGY-SOURCE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_passed':result['all_passed'],'checks':len(checks),'largest_error':max(err.values())}))
