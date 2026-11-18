"""Untruncated original whole96 fourth heat vertex and weight combinatorics.

The old rng252 inputs and one diagnostic proper-time value check identities.
They do not select physical Yukawas, nu, a field measure or a new source.
"""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json
import numpy as np
from scipy.linalg import expm

ROOT = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0, str(ROOT/'verification'))
import v252_full_finite_triple as original

rng = np.random.default_rng(252)
Y = [rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for _ in range(4)]
mr0 = rng.normal(size=(3,3)); MR=mr0+mr0.T
z=np.zeros((3,3))
DY=original.full_D(*Y,z)
DM=original.full_D(z,z,z,z,MR)
A=DM@DM; L=DM@DY+DY@DM; B=DY@DY
n=96
left_sign=np.array([-1 if row[1]=='L' else 1 for row in original.P])
Uweak=np.diag(np.r_[left_sign,left_sign])
oddness_error=float(np.linalg.norm(Uweak@DY@Uweak+DY))
majorana_fixed_error=float(np.linalg.norm(Uweak@DM@Uweak-DM))
assert oddness_error==0 and majorana_fixed_error==0

def toeplitz_coefficient(heat_time, order):
    block=np.zeros(((order+1)*n,(order+1)*n),complex)
    for i in range(order+1):
        block[i*n:(i+1)*n,i*n:(i+1)*n]=A
        if i+1<=order:block[i*n:(i+1)*n,(i+1)*n:(i+2)*n]=L
        if i+2<=order:block[i*n:(i+1)*n,(i+2)*n:(i+3)*n]=B
    return np.trace(expm(-heat_time*block)[:n,order*n:(order+1)*n])

heat_time=.4
c0=np.trace(expm(-heat_time*A))
c2=toeplitz_coefficient(heat_time,2)
c4=toeplitz_coefficient(heat_time,4)
# The old exact full96 radial second derivative is reused only as a control
# on coefficient extraction; it is not counted as the new result.
Ach=np.trace(3*Y[0].conj().T@Y[0]+3*Y[1].conj().T@Y[1]+Y[2].conj().T@Y[2]).real
c2_closed=-4*heat_time*(Ach+np.trace(expm(-heat_time*MR@MR.conj().T)@Y[3].conj().T@Y[3]).real)

contours=[]
for radius in (.08,.12):
    N=64
    angles=2*np.pi*np.arange(N)/N
    traces=[]
    for theta in angles:
        r=radius*np.exp(1j*theta)
        traces.append(np.trace(expm(-heat_time*(DM+r*DY)@(DM+r*DY))))
    traces=np.array(traces)
    tc4=np.mean(traces*np.exp(-4j*angles))/radius**4
    weights=np.exp(-(traces-c0))
    wc4=np.mean(weights*np.exp(-4j*angles))/radius**4
    predicted=c2*c2/2-c4
    contours.append({'radius':radius,'nodes':N,'trace_fourvertex':float(tc4.real),
       'trace_relative_error':float(abs(tc4-c4)/max(1,abs(c4))),
       'real_profile_weight_fourvertex':float(wc4.real),
       'weight_relative_error':float(abs(wc4-predicted)/max(1,abs(predicted)))})

# Order-sensitive path subdivision preserves the complete matrices. Splitting
# into A/L/B commuting scalar means is not used anywhere above.
noncomm=float(np.linalg.norm(A@B-B@A))
sources=[
 'verification/v252_full_finite_triple.py',
 'experiments/theory-contracts/recursive-mirrored-source-20261103/ORIGINAL-PROFILE-DICTIONARY.txt',
 'experiments/theory-contracts/original-local-source-variable-background-20261102/CURRENT-ARGUMENT-MAP.txt',
 'documents/archive--paper-latex--tfpt-31/sections/007-subsection-occupancy-and-boundary-to-bulk-bridge.tex',
]
result={
 'scope':'Unchanged original finite96 matrices; new untruncated fourth proper-time coefficient and full real-profile weight fourth coefficient. No physical profile or field integration.',
 'diagnostic_heat_time':heat_time,'dimension':n,
 'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
 'noncommuting_A_B_norm':noncomm,
 'actual_weak_central_sign_DY_odd_error':oddness_error,
 'actual_weak_central_sign_DM_fixed_error':majorana_fixed_error,
 'odd_trace_coefficients_c1_c3':'exactly zero by original weak central sign conjugation',
 'c0':float(c0.real),'c2':float(c2.real),'c4':float(c4.real),
 'real_profile_weight_c4':float((c2*c2/2-c4).real),
 'old_second_coefficient_control_error':float(abs(c2-c2_closed)),
 'independent_complex_cauchy_controls':contours,
 'T_even_identical_whole_arguments_relative_c2_c4':[0,0],
 'S_even_fixed_reference_c4':'numerator_c4',
 'physical_nu_selected':False,'eta_and_later_terms_not_numerically_evaluated':True,
 'primitive_full_source_evaluated':False,
}
assert noncomm>1
assert abs(c2-c2_closed)<1e-9
assert all(q['trace_relative_error']<1e-8 and q['weight_relative_error']<1e-8 for q in contours)
assert max(abs(c0.imag),abs(c2.imag),abs(c4.imag))<1e-8
out=Path(__file__).with_name('original-loop-fourvertex-results.json')
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
