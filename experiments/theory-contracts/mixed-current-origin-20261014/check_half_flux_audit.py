"""Independent controls of the pasted half-flux argument, not its missing package.
Analytic local nonnormality is proved in the accompanying review, not inferred
from this finite census. Original verification modules are only imported/read.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib, importlib, json, sys
import numpy as np
import sympy as S

REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(REPO/'verification'))
v746=importlib.import_module('v746_phys_gnet_local_functor')
v1033=importlib.import_module('v1033_charged_disorder')

def need(b,label):
    if not b:raise ValueError(label)

# Actual Z4 Weyl convention from v1033 symbolic_construction.
Z=S.diag(1,S.I,-1,-S.I)
X=S.zeros(4)
for j in range(4):X[(j+1)%4,j]=1
need(X**2*Z*(X**2).T==-Z,'register clock remains a half-flux witness')
need(X**4==S.eye(4),'actual order-four register')

# The v746 paired-window witness, using the CURRENT original shift/window.
clock_cases=[]
for N in [8,16,48]:
    H=np.linalg.matrix_power(v746.shift(N),N//2)
    sigma=np.diag(np.where(np.arange(N)<N//2,-1.,1.))
    need(np.array_equal(H@H,-np.eye(N)),'NS clock square')
    for j in range(N//2):
        pair=v746.window(N,j,1)
        need(pair==[j,j+N//2],'actual doubled window')
        K=np.zeros((N,N),complex)
        K[j,j+N//2]=1j;K[j+N//2,j]=-1j
        need(np.array_equal(H@K@H.T,K),'clock-invariant pair current')
        need(np.array_equal(sigma@K@sigma,-K),'half-flux flips invariant pair current')
    clock_cases.append({'N':N,'all_base_pair_windows':N//2,'passed':True})

# Current actual half-filled cylinder state convention, not a pure zero sea.
strip=[]
for sector in [0,2]:
    h=v1033.qwz_cylinder(8,8,1.,sector)
    _,C,zero,k=v1033.half_filled_ground_space(h)
    d=zero.shape[1]
    need((d,k)==((2,1) if sector==0 else (0,0)),'actual R/NS zero-mode census')
    strip.append({'sector':sector,'zero_modes':d,'zero_occupied':k,
                  'covariance_idempotence_defect':float(np.linalg.norm(C@C-C))})

rho=S.diag(0,S.Rational(1,2),S.Rational(1,2),0)
quasi=S.eye(4)/4
pair=S.diag(0,0,0,1)
need(S.trace(rho*pair)==0 and S.trace(quasi*pair)==S.Rational(1,4),'quartic state distinction')
need(S.trace(rho*rho)==S.Rational(1,2) and rho.rank()==2,'mixed R state purity')
for d in range(1,8):
    for k in range(d+1):
        for m in range(d+1):
            lhs=F(int(S.binomial(d-m,k-m)),int(S.binomial(d,k))) if k>=m else F(0)
            rhs=F(1)
            for i in range(m):rhs*=F(k-i,d-i)
            need(lhs==rhs,'fixed-filling full normally ordered moment')

# Exact analytic cutoff identity; divergence follows by delta ->0.
delta,eps=S.symbols('delta eps',positive=True)
# Integrate first in v, then u, using the elementary antiderivatives.
boundary=2*S.log(eps+delta)-S.log(2*eps)-S.log(2*delta)
claimed=S.log((eps+delta)**2/(4*eps*delta))
need(S.simplify(S.expand_log(boundary-claimed,force=True))==0,'cutoff logarithm')

# Favorable J=C counterexample: R is real-linear but not complex-linear.
I2=S.Matrix([[0,-1],[1,0]])
complex_i=S.diag(I2,I2)
real_R=S.diag(1,1,1,-1)
gamma=S.diag(1,1,-1,-1)
need(real_R*complex_i*real_R==gamma*complex_i,'sector-selective conjugation changes complex scalar')
need(real_R*complex_i!=complex_i*real_R,'not complex linear')
need(real_R*complex_i!=-complex_i*real_R,'not globally antilinear')

files=['verification/v746_phys_gnet_local_functor.py','verification/v1033_charged_disorder.py']
out={'date':'2026-10-22','verdict':'SHARP_HALF_FLUX_EXCLUSION_CONFIRMED_IN_STATED_NETS',
     'register_obstruction_exact':True,'quotient_clock_witness_cases':clock_cases,
     'original_strip_controls':strip,'R_fixed_filling_state_purity':'1/2',
     'R_zero_mode_pair_moment':'0','same_covariance_quasifree_pair_moment':'1/4',
     'local_cutoff_HS_integral':'log((eps+delta)^2/(4 eps delta)) -> infinity',
     'selective_antiunitary_candidate_complex_linear':False,
     'finite_or_analytic_control_scope':'No replay of unavailable 18+38 report package; analytic proof reviewed separately.',
     'all_P1_no_go':False,'source_hashes':{f:hashlib.sha256((REPO/f).read_bytes()).hexdigest() for f in files}}
Path(__file__).with_name('half-flux-audit-results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
