#!/usr/bin/env python3
"""Regress the native Gaussian-transfer normalization, not a truncated charge model."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,itertools,math
import numpy as np
L=3;N=L**3
sites=list(itertools.product(range(L),repeat=3));index={s:i for i,s in enumerate(sites)}
lap=np.zeros((N,N))
for x,s in enumerate(sites):
 for axis in range(3):
  v=list(s);v[axis]=(v[axis]+1)%L;y=index[tuple(v)]
  lap[x,x]+=1;lap[y,y]+=1;lap[x,y]-=1;lap[y,x]-=1
# Original admissible inputs for diagnostic: a=m=1,u=1/3,nu=0; p=e8 has e(p)=1.
energies=np.zeros(N);energies[0]=energies[index[(1,0,0)]]=1
D=float(energies.sum()/N) # nu=0; retain the original charge energy in the block trace.
A=lap+np.eye(N)+(2/3)*np.diag(energies)
omega=np.sqrt(np.linalg.eigvalsh(A))
rows=[]
for T in (3,4,8,16):
 delta=1/T
 Lt=np.zeros((T,T))
 for j in range(T):
  k=(j+1)%T;Lt[j,j]+=1;Lt[k,k]+=1;Lt[j,k]-=1;Lt[k,j]-=1
 M=np.kron(Lt,np.eye(N))/delta+delta*np.kron(np.eye(T),A)
 sign,logdet=np.linalg.slogdet(M)
 logZ_integral=-D-N*T*math.log(delta)/2-logdet/2
 alpha=2*np.arcsinh(delta*omega/2)
 logZ_spectrum=-D-np.log(2*np.sinh(T*alpha/2)).sum()
 a=alpha/np.sinh(alpha);b=alpha*np.sinh(alpha)/delta**2
 assert sign>0 and abs(logZ_integral-logZ_spectrum)<2e-10
 assert np.all(a+1e-13>=1/(1+delta**2*omega**2/4)) and np.all(b+1e-12>=omega**2)
 rows.append({'T':T,'delta':delta,'logZ_gaussian_integral':float(logZ_integral),'logZ_reconstructed_spectrum':float(logZ_spectrum),'absolute_residual':float(abs(logZ_integral-logZ_spectrum)),'effective_frequencies_max_error':float(np.max(abs(alpha/delta-omega)))})
# Uniform resolvent bound diagnostics; the all-x theorem is analytic in PROOF.md.
x=np.logspace(-7,7,301)
res=[]
for shift in (.01,1.,100.):
 value=1/(-np.expm1(-x)+shift)-1/(x+shift)
 assert value.min()>-1e-12 and value.max()<1+1e-12
 res.append({'delta_times_z':shift,'min_ratio':float(value.min()),'max_ratio':float(value.max())})
result={'status':'PASS','scope':'Gaussian log and normalization in one original static neutral charge block; analytical operator proof handles moving charges and nonzero signed hopping. No full charge history simulation.','declared_diagnostic_inputs':{'L':3,'a':1,'m':1,'u':'1/3','beta':1,'nu':0,'D':D,'charge_profile':'p=e8 at (0,0,0), -p at (1,0,0), zero elsewhere'},'gaussian_checks':rows,'resolvent_bound_samples':res,'max_normalization_residual':max(r['absolute_residual'] for r in rows)}
Path(__file__).with_name('reconstruction-result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
