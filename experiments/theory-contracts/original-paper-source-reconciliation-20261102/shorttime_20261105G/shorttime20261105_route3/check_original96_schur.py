"""Whole-spectral-parameter Schur gluing on the unchanged original R055/v252 fixture.
This is an internal block cut, not a new geometric boundary or physical measure.
"""
from pathlib import Path
import sys, json, hashlib
import numpy as np
from scipy.integrate import quad_vec

ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(ROOT/'verification'))
import v252_full_finite_triple as v

rng=np.random.default_rng(252)
Ys=[rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for _ in range(4)]
mr=rng.normal(size=(3,3));MR=mr+mr.T
z3=np.zeros((3,3),complex);z2=np.zeros((2,2),complex)
DY=v.full_D(*Ys,z3);DM=v.full_D(z3,z3,z3,z3,MR)
algP=v.rep(1,z2,z3);Jswap=v.Kswap();gamma=v.grading()
# Original Hstar=(0,1), p=alpha=1/4, as used in the original R055 check.
q=np.eye(2,dtype=complex)
a=v.rep(0,q,z3);one=a@(DY@algP-algP@DY);one=one+one.conj().T
Phi=one+Jswap@one.conj()@Jswap+DM
p=.25;h=p*gamma-1j*gamma@Phi
K=h@h
A=K[:48,:48];B=K[:48,48:];C=K[48:,48:]
Ia=np.eye(48);Ik=np.eye(96)
evK=np.linalg.eigvalsh(K);evA=np.linalg.eigvalsh(A)
checks=[];resolvent_rows=[]
def need(name,ok):
 checks.append({'name':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
need('whole96_native_h_hermitian',np.max(abs(h-h.conj().T))<1e-12)
need('whole96_native_square',np.max(abs(K-(p*p*Ik+Phi@Phi)))<1e-11)
cross_norm=float(np.trace(B.conj().T@B).real)
moment2=float((np.trace(K@K)-np.trace(A@A)-np.trace(C@C)).real)
need('native96_first_cross_return_moment',abs(moment2-2*cross_norm)<1e-8)
# Source-specific invisible sector of this original particle/anti internal cut.
neutrinos=[i for i,(name,_,_,_) in enumerate(v.P) if name=='nu']
dark=[i for i in range(48) if i not in neutrinos]
need('actual_neutrino_subspace_reduces_A',np.max(abs(A[np.ix_(dark,neutrinos)]))<1e-11)
need('actual_cross_block_has_no_dark_rows',np.max(abs(B[dark,:]))<1e-11)
rankB=int(np.linalg.matrix_rank(B,tol=1e-9))
need('actual_original_cut_cyclic_rank_six',rankB==6)
need('actual_dark_dimension_42',len(dark)==42)
for z in [1.,2.+.7j,5.+2j]:
 R=np.linalg.inv(A+z*Ia)
 S=C+z*Ia-B.conj().T@R@B
 Sp=Ia+B.conj().T@R@R@B
 direct=np.trace(np.linalg.inv(K+z*Ik))-np.trace(R)
 schur=np.trace(np.linalg.solve(S,Sp))
 error=abs(direct-schur)
 need('entire_original_resolvent_'+str(z),error<1e-10)
 # Derivative includes the interior return factor, not just S inverse.
 static=np.trace(np.linalg.inv(S))
 resolvent_rows.append({'z':[float(np.real(z)),float(np.imag(z))],
  'direct':[float(direct.real),float(direct.imag)],
  'schur':[float(schur.real),float(schur.imag)],
  'error':float(error),'missing_return_if_Sprime_omitted':float(abs(direct-static))})
need('static_inverse_loses_nonzero_return',all(r['missing_return_if_Sprime_omitted']>1e-5 for r in resolvent_rows))
# Cauchy heat reconstruction from the whole Schur family only, independent of
# diagonalizing K to obtain the boundary resolvent. tau=1 is the preexisting
# R055 heat-component diagnostic, not source time or a new physical profile.
tau=1.;left=float(evK[0]/2);right=float(evK[-1]+1);height=2.
corners=[left-1j*height,right-1j*height,right+1j*height,left+1j*height]
def boundary_resolvent(w):
 R=np.linalg.inv(w*Ia-A)
 S=w*Ia-C-B.conj().T@R@B
 return np.trace(np.linalg.solve(S,Ia+B.conj().T@R@R@B))
contour=0j;quad_errors=[]
for j in range(4):
 start=corners[j];step=corners[(j+1)%4]-start
 value,error=quad_vec(lambda u:np.exp(-tau*(start+step*u))*boundary_resolvent(start+step*u)*step,
                      0.,1.,epsabs=1e-10,epsrel=1e-10)
 contour+=value;quad_errors.append(float(error))
heat_from_boundary=contour/(2j*np.pi)
heat_direct=float(np.exp(-tau*evK).sum()-np.exp(-tau*evA).sum())
need('native96_whole_heat_gluing_from_SpectralFamily',abs(heat_from_boundary-heat_direct)<1e-8)
# The real original S reference has no H-dependent Phi on even H; its squared
# whole p block is p²I. Retain its absolute spectral contribution separately.
reference_heat=96*np.exp(-tau*p*p)
out={'scope':'Unchanged whole96 R055/v252 seed252 diagnostic at original Hstar,p1/4. Internal particle/anti cut; no physical boundary identification, no new Yukawa, measure, spectral profile, source time or primitive field integral.',
 'checked':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks,
 'cross_BdaggerB_trace':cross_norm,'first_cross_return_coefficient_2TrBdaggerB':2*cross_norm,
 'moment2_direct':moment2,'actual_internal_cut_cyclic_dimension':rankB,'actual_internal_cut_dark_dimension':len(dark),
 'resolvent_controls':resolvent_rows,
 'heat_gluing':{'proper_time_diagnostic':tau,'source_time_identification':False,
  'boundary_reconstructed':[float(heat_from_boundary.real),float(heat_from_boundary.imag)],
  'native_direct':heat_direct,'absolute_error':float(abs(heat_from_boundary-heat_direct)),
  'quadrature_error_estimates':quad_errors,'native_full_heat':float(np.exp(-tau*evK).sum()),
  'original_S_reference_heat':float(reference_heat)},
 'provenance':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in [ROOT/'verification/v252_full_finite_triple.py',ROOT/'experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/check_exact_cutoff_higgs_measure.py']},
 'primitive_boson_body_evaluated':False}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['checked','total','cross_BdaggerB_trace','first_cross_return_coefficient_2TrBdaggerB','actual_internal_cut_cyclic_dimension','actual_internal_cut_dark_dimension','heat_gluing']},indent=2))
