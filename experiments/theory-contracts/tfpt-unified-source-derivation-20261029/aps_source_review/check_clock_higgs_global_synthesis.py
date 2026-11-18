"""R083: conditional clock/Higgs unification tests on existing source data.

Shift and C5 carrier identifications are hypotheses, not physical TFPT maps.
"""
from pathlib import Path
import json, sys
import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy.special import expit

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'verification'))
import v252_full_finite_triple as original
from v311_gap_decoupled_measure import transfer

checks={}
def require(name, value):
    checks[name]=bool(value)
    if not value: raise RuntimeError(name)

# The actual R081 seeded input, not a new physical Yukawa fit.
rng=np.random.default_rng(252)
ys=[rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for _ in range(4)]
mr=rng.normal(size=(3,3))
D=original.full_D(*ys,mr+mr.T)
require('original_finite_triple_96',D.shape==(96,96))
s=float(np.linalg.svd(ys[0],compute_uv=False)[0]);ell=.2;rad=.3;m=rad*s;alpha=.25
def c(p, mass):
    e=np.sqrt(p*p+mass*mass)
    return (1+(p/e)*np.tanh(ell*e))/2
def dc(p,mass):
    e=np.sqrt(p*p+mass*mass);u=np.tanh(ell*e)
    return p*mass/(2*e**3)*(ell*e*(1-u*u)-u)
p=alpha;h=1e-5
K=2*ell*np.array([[p,m],[m,-p]])
e,u=np.linalg.eigh(K);numeric=((u*expit(e))@u.T)[0,0]
require('actual_singular_pair_occupation',abs(c(p,m)-numeric)<1e-14)
der=dc(p,m)
require('higgs_derivative_nonzero',abs(der)>1e-5)
require('higgs_derivative_finite_difference',abs(der-(c(p,m+h)-c(p,m-h))/(2*h))<1e-9)

# Conditional same scalar Hardy component and unilateral shift, not established
# equivalence to the original full geometric or E8 carrier.
z=.3; limits=[]
for nmax in [100,1000,10000]:
    ns=np.arange(-nmax,nmax+1,dtype=float)
    cs=c(ns+alpha,m); prev=c(ns+alpha-1,m)
    logdet=float(np.sum(np.log1p(-z*prev)-np.log1p(-z*cs)))
    boundary=float(np.log1p(-z*c(-nmax+alpha-1,m))-np.log1p(-z*c(nmax+alpha,m)))
    require('telescoping_identity_'+str(nmax),abs(logdet-boundary)<1e-11)
    trace_dist=float(np.sum(abs(cs-(ns+alpha>0))))
    limits.append({'cutoff':nmax,'relative_determinant':float(np.exp(logdet)),
                   'distance_to_target':float(abs(np.exp(logdet)-1/(1-z))),
                   'trace_distance_to_step':trace_dist})
require('shifted_determinant_converges_to_winding',limits[-1]['distance_to_target']<1e-7)
require('finite_cutoffs_improve',limits[-1]['distance_to_target']<limits[0]['distance_to_target'])
require('relative_trace_tail_bounded',limits[-1]['trace_distance_to_step']<2*limits[0]['trace_distance_to_step'])
# Compare topology against a genuinely varying matter response.
require('same_higgs_response_varies',abs(c(p,m)-c(p,0))>1e-5)

# AGL(1,5) on the *natural five-letter permutation carrier*. This is only a
# proposed common carrier; v311's 3x3 transfer is not asserted to live on it.
C=sp.zeros(5);G=sp.zeros(5)
for j in range(5): C[(j+1)%5,j]=1;G[(2*j)%5,j]=1
require('actual_clock_semidirect_relation',G*C*G.inv()==C**2)
xs=sp.symbols('x0:25');X=sp.Matrix(5,5,xs)
eq=list(X*C-C*X)+list(X*G-G*X)
M,_=sp.linear_eq_to_matrix(eq,xs);rank=M.rank(); nullity=25-rank
require('five_carrier_commutant_dimension_two',nullity==2)
require('all_ones_and_identity_in_commutant',C*sp.ones(5)==sp.ones(5)*C and G*sp.ones(5)==sp.ones(5)*G)
T,l2,l3=transfer();eigs=np.linalg.eigvalsh(T)
require('v311_three_distinct_levels',np.min(np.diff(eigs))>1e-6)

out={'scope':'conditional same scalar shift/Higgs carrier; natural five-letter clock representation; no identification with actual joint TFPT process',
     'singular_value':s,'radius':rad,'ell':ell,'alpha':alpha,'mass_control':m,
     'one_mode_occupation':float(c(p,m)),'one_mode_mass_derivative':float(der),
     'relative_determinant_target':1/(1-z),'shift_cutoffs':limits,
     'clock_commutant_dimension':nullity,'native_transfer_eigenvalues':eigs.tolist(),
     'checks':checks,'passed':len(checks),'failed':0,'new_lean':False,'physical_source_or_spacetime_derived':False}
Path(__file__).with_name('CLOCK-HIGGS-GLOBAL-SYNTHESIS.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
