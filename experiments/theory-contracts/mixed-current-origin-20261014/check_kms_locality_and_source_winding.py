#!/usr/bin/env python3
"""Controls for the written source/KMS audit, not an E8 origin claim."""
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json
import math
import sys

import numpy as np
import sympy as s

HERE=Path(__file__).resolve().parent
REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')

def need(ok, why):
    if not bool(ok):
        raise RuntimeError(why)

# Exact two-mode subalgebra of the original CAR, not a finite-net substitute.
def car(n,j):
    a=s.zeros(2**n)
    for b in range(2**n):
        if (b>>j)&1:
            a[b^(1<<j),b]=(-1)**sum((b>>k)&1 for k in range(j))
    return a

a,b=car(2,0),car(2,1)
nf=a.H*a
exterior=a.H*b+b.H*a
comm=nf*exterior-exterior*nf
need(comm==a.H*b-b.H*a,'analytic even exterior commutator')
need(comm.rank()==2 and comm.H*comm==s.diag(0,1,1,0),
     'nonzero exact even-locality witness')

# Refined central character of the existing odd Phi, checked on source charges.
r=s.Matrix([s.Rational(1,2)]*4+[-s.Rational(1,2)]*4)
charges=[s.zeros(8,1)]
for i in range(8):
    for v in (-2,-1,1,2):
        p=s.zeros(8,1);p[i]=v;charges.append(p)
for i in range(7):
    p=s.zeros(8,1);p[i]=1;p[i+1]=1;charges.append(p)
for p in charges:
    parity=1 if int(sum(p))%2==0 else -1
    q=p if parity==1 else s.Matrix([p[(k+4)%8] for k in range(8)])+r
    need(sum(q)==sum(p),'existing Phi total charge')
    need(all(s.exp(2*s.pi*s.I*qk)==parity for qk in q),
         'target individual Cartan holonomy equals old parity')
    need(all(s.exp(2*s.pi*s.I*pk)==1 for pk in p),
         'old individual Cartan holonomy is identity')

# Exact truncated ORIGINAL edge spectrum, partition and parity trace.
# This is a finite control of analytic infinite products proved in the text.
x=s.symbols('x')
levels=[1,3,5,7,9,11,13,15]  # energies in units of1/4, eitherquarter.
zpoly=s.prod(1+x**e for e in levels)
spoly=s.prod(1-x**e for e in levels)
ztrace=0;strace=0
for bits in itertools.product((0,1),repeat=len(levels)):
    term=x**sum(v*e for v,e in zip(bits,levels))
    ztrace+=term;strace+=(-1)**sum(bits)*term
need(s.expand(ztrace-zpoly)==0 and s.expand(strace-spoly)==0,
     'full finite original-edge trace/product identity')
ze,zo=(zpoly+spoly)/2,(zpoly-spoly)/2
need(s.cancel((ze**2-zo**2)/(ze**2+zo**2)
              -2*zpoly*spoly/(zpoly**2+spoly**2))==0,
     'total-even two-edge central character formula')
need(s.cancel(1-2*zpoly*spoly/(zpoly**2+spoly**2)
              -(zpoly-spoly)**2/(zpoly**2+spoly**2))==0,
     'strictly positive two-edge thermal defect')

# Source charge-one winding parity, checked algebraically without making a
# truncated shift a supposedly unitary finite spectral-flow implementation.
u,v=s.symbols('u v')
P=s.diag(1,-1);W=s.Matrix([[0,u],[v,0]]);Pe=(s.eye(2)+P)/2
need(P*W+W*P==s.zeros(2),'odd source winding parity relation')
need(Pe*W*Pe==s.zeros(2),'left-even projection removes source W1 compression')
need(W.H*P*W*P==-W.H*W,'exact odd-winding/Cartan four-word sign')

thermal=[]
for alpha in (.25,.75):
    for beta in (.7,2*math.pi,8.):
        energies=np.array([e for m in range(96) for e in (m+alpha,m+1-alpha)])
        ratio=float(np.prod(np.tanh(beta*energies/2)))
        even_expect=2*ratio/(1+ratio*ratio)
        need(0<ratio<1 and 0<even_expect<1,'all-beta formulas finite controls')
        thermal.append({'alpha':alpha,'beta_source_H':beta,
            'single_edge_target_Cartan_holonomy':ratio,
            'total_even_two_edge_target_holonomy':even_expect,
            'old_Cartan_holonomy':1.})

# Actual original QWZ time moves an even number observable out of its region.
sys.path.insert(0,str(REPO/'verification'))
source=REPO/'verification/v1033_charged_disorder.py'
spec=importlib.util.spec_from_file_location('kms_audit_v1033',source)
mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod
spec.loader.exec_module(mod)
actual=[]
for sector in (1,3):
    nx,ny=8,8
    h=mod.qwz_cylinder(nx,ny,1.,sector)*nx/(2*math.pi)
    es,vs=np.linalg.eigh(h)
    t=.2
    U=(vs*np.exp(1j*t*es))@vs.conj().T
    j0=0                  # source site inside column0
    j1=2*ny              # exterior observable in adjacent column1
    f=U[:,j0]
    number_t=np.outer(f,f.conj())
    number_out=np.zeros_like(h);number_out[j1,j1]=1
    c=number_t@number_out-number_out@number_t
    cnorm=float(np.linalg.norm(c,2))
    need(cnorm>1e-4,'actual original-time fixed-region invariance fails')
    original=np.zeros_like(h);original[j0,j0]=1
    need(np.linalg.norm(original@number_out-number_out@original)<1e-14,
         'the original even local observables commute before time transport')
    # Single-mode reduced thermal state: modular flow fixes its number.
    beta=.7
    C=(vs*(1/(1+np.exp(beta*es))))@vs.conj().T
    occupation=float(C[j0,j0].real)
    need(0<occupation<1,'actual restricted finite thermal state is faithful')
    actual.append({'sector':sector,'Nx':nx,'Ny':ny,'physical_source_time':t,
        'outside_even_number_commutator_norm':cnorm,
        'local_thermal_occupation':occupation,
        'local_modular_one_mode_energy':math.log((1-occupation)/occupation)})

result={'date':'2026-10-23','checks_passed':True,
    'verdict':'KMS_GIVEN_EXTENSION_DOES_NOT_SUPPLY_MICROSCOPIC_LOCAL_SPINOR_MAP',
    'analytic_CAR_witness_rank':int(comm.rank()),
    'refined_charge_samples':len(charges),
    'finite_edge_fock_states':2**len(levels),
    'thermal_original_edge_controls':thermal,'actual_QWZ_time_controls':actual,
    'actual_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
    'global_rotations_equal_fixed_interval_modular_flow':False,
    'conditional_expectation_uniqueness_refuted':False,
    'source_W1_identified_with_primitive_P1_circle':False,
    'unitary_source_winding_Cartan_four_word_signs':{'old':1,'target':-1},
    'new_microscopic_spinor_sector':False,'all_P1_no_go':False}
(HERE/'kms-locality-and-source-winding-results.json').write_text(
    json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'checks_passed':True,'actual_quarters':2,
    'fixed_interval_invariance_failed':True,'thermal_dictionary_witness':True,
    'source_W1_survives_left_even_projection':False}))
