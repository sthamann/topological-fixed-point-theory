#!/usr/bin/env python3
"""Direct checks of specific equations in the user-supplied synthesis."""
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import sys

for variable in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[variable]='1'
import numpy as np
import mpmath as mp

ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(ROOT/'verification'))
spec=importlib.util.spec_from_file_location('original_v1024',ROOT/'verification/v1024_parent_internal_boundaries.py')
mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
captured=io.StringIO()
with contextlib.redirect_stdout(captured):
    data=mod.parent_internal_intertwiner()
Path(__file__).with_name('v1024-parent-replay.txt').write_text(captured.getvalue())
if not all(x.passed for x in mod.RESULTS): raise RuntimeError('original parent check failed')
H=data['h_intertwined']
energies,vecs=np.linalg.eigh(H)
probabilities=np.exp(-energies);probabilities/=probabilities.sum()
rho=(vecs*probabilities)@vecs.conj().T
ann=[mod.jordan_wigner_annihilator(6,k) for k in (0,1)]
numbers=[np.kron(np.eye(4),a.conj().T@a) for a in ann]
n0,n1=numbers;pair=n0@n1
av=lambda X:float(np.trace(X).real/256)
covH=lambda X:av(X@H)-av(X)*av(H)
slope=-covH(pair)+av(n0)*covH(n1)+av(n1)*covH(n0)
kappa=float(np.trace(rho@pair).real-np.trace(rho@n0).real*np.trace(rho@n1).real)
actual_block=np.array([[7/3,-2/15],[-2/15,83/150]])
block_eval=np.linalg.eigvalsh(actual_block)
delta=89/50;c=2/15;B0=83/150
pasted=np.array([B0+1/4+(delta-sign*np.sqrt(delta**2+4*c*c))/2 for sign in (1,-1)])
corrected=pasted-1/4
if not np.allclose(corrected,block_eval,atol=1e-13):raise RuntimeError('block correction failed')
mp.mp.dps=60
c3=1/(8*mp.pi);top=48*c3**4;phi0=1/(6*mp.pi)+top
def f(alpha):
    q=top*mp.exp(-2/alpha)
    phis=1/(6*mp.pi)+q*(1-q)**(-mp.mpf(5)/4)
    return alpha**3-2*c3**3*alpha**2-mp.mpf(4)/5*41*c3**6*mp.log(1/phis)
alpha=mp.findroot(f,(mp.mpf(1)/138,mp.mpf(1)/136))
em_spec=importlib.util.spec_from_file_location('original_v3',ROOT/'verification/v3_em_alpha.py')
em=importlib.util.module_from_spec(em_spec);em_spec.loader.exec_module(em)
original_f=em.make_F(c3,41)
original_alpha=mp.findroot(original_f,(mp.mpf(1)/138,mp.mpf(1)/136))
# Exact-scope counterexample to the missing MaxEnt compatibility premise.
tau=np.array([.75,.25]);sigma=np.array([.25,.75])
S=lambda p:float(-np.sum(p*np.log(p)))
D=float(np.sum(sigma*np.log(sigma/tau)))
result={
 'original_parent_checks_passed':len(mod.RESULTS),
 'actual_v1024_beta1_pair_cumulant':kappa,
 'actual_v1024_beta0_cumulant_derivative':slope,
 'beta0_derivative_matches_minus_1_over_300':abs(slope+1/300)<1e-13,
 'ground_pair_bright_block_actual_eigenvalues':block_eval.tolist(),
 'pasted_formula_k1_eigenvalues':pasted.tolist(),
 'pasted_formula_constant_error_for_zero_single_occupancy':.25,
 'correct_general_bright_formula':'k*B0 + s/4 + (delta +/- sqrt(delta^2+4*k*c^2))/2; k>=1, s=number of singly occupied families, k+s<=3',
 'vacuum_transition_amplitude_bound':str(mp.mpf(1600)/72889),
 'alpha_inverse_from_exact_pasted_function':str(1/alpha),
 'alpha_inverse_from_original_v3_function':str(1/original_alpha),
 'alpha_formula_correction':'q(alpha)=delta_top*exp(-2*alpha), not delta_top*exp(-2/alpha); the latter exponential belongs to a different cosmological readout.',
 'phi0':str(phi0),'lambda':str(mp.sqrt(phi0*(1-phi0))),
 'pasted_maxent_without_E_sigma_equal_tau':{'E':'identity','rho_tau':tau.tolist(),'sigma':sigma.tolist(),'S_tau_minus_S_sigma':S(tau)-S(sigma),'D_sigma_tau':D,'correction':'Require E(sigma)=tau=E(rho); otherwise the displayed identity is false.'},
 'scope':'Selected direct equation checks; not validation of the entire pasted document or original-physics derivation.'}
Path(__file__).with_name('pasted-formula-result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
