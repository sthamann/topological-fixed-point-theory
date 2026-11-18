#!/usr/bin/env python3
"""Actual root-frame duality, source metric, and overlap compatibility.
No new microscopic Hamiltonian or physical execution is assumed.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
from pathlib import Path
import sys
import numpy as np
import sympy as sp

REPO=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(REPO))
from tfpt_explorer.cartan_source import _roots,_rays
from tfpt_explorer.process import _code_basis
from check_collective_origin import need,sym_scaled


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=Path(__file__).with_name('dynamics_results.json'))
    args=parser.parse_args()
    rays=_rays(_roots())
    words=list(it.product(range(4),repeat=4))
    fr=np.zeros((256,60),dtype=np.int64)
    fi=np.zeros_like(fr)
    zs=[[(r[2*j],r[2*j+1]) for j in range(4)] for r in rays]
    for j,z in enumerate(zs):
        for i,w in enumerate(words):
            a,b=1,0
            for k in w:
                a,b=a*z[k][0]-b*z[k][1],a*z[k][1]+b*z[k][0]
            fr[i,j],fi[i,j]=a,b
    g=fr.T@fr+fi.T@fi
    gi=fr.T@fi-fi.T@fr
    need(not np.any(gi),'Actual fourth Gram is real')
    # Independent construction by fourth powers of one-copy inner products.
    direct=np.zeros_like(g)
    direct_i=np.zeros_like(g)
    for i,a in enumerate(zs):
        for j,b in enumerate(zs):
            re=sum(x*u+y*v for (x,y),(u,v) in zip(a,b))
            im=sum(x*v-y*u for (x,y),(u,v) in zip(a,b))
            x,y=1,0
            for _ in range(4):
                x,y=x*re-y*im,x*im+y*re
            direct[i,j],direct_i[i,j]=x,y
    need(np.array_equal(g,direct) and not np.any(direct_i),'Direct kernel matches tensor Gram')
    unit=np.eye(60,dtype=np.int64)
    need(not np.any(g@(g-384*unit)@(g-768*unit)),'Exact spectral polynomial')
    tr,tr2=int(np.trace(g)),int(np.trace(g@g))
    top=Q(tr2-384*tr,768*(768-384))
    middle=Q(tr-768*top,384)
    need(top==5 and middle==30,'Spectral multiplicities from exact traces')
    row_counts=[]
    for row in g:
        row_counts.append({str(Q(int(k),256)):v for k,v in sorted(Counter(map(int,row)).items())})
    need(all(c=={'-1/4':12,'0':15,'1/16':32,'1':1} for c in row_counts),'All60 interaction rows')
    triangle=next(t for t in it.combinations(range(60),3) if g[t[0],t[1]]*g[t[1],t[2]]*g[t[2],t[0]]<0)
    product=Q(int(g[triangle[0],triangle[1]])*int(g[triangle[1],triangle[2]])*int(g[triangle[2],triangle[0]]),256**3)
    need(product==Q(-1,64),'Gauge-invariant negative cycle')
    need(not np.any(g.sum(axis=1)),'Canonical constant label vector is null')
    ag=np.abs(g)
    need(np.all(ag.sum(axis=1)==1536) and np.all(ag@ag>0),'Magnitude-only Perron mode is unique at6')
    need(np.any(ag@(ag-384*unit)@(ag-768*unit)),'Magnitude-only negative control breaks spectrum')
    # K^+=(7K-2K^2)/9; K=g/256. F=(fr+i fi)/16.
    pin=7*256*g-2*g@g
    wr=fr@pin@fr.T+fi@pin@fi.T
    wi=fi@pin@fr.T-fr@pin@fi.T
    need(np.array_equal(wr,6291456*sym_scaled(4)) and not np.any(wi),'F K+ F* is entire Sym4')
    # This is the literal code overlap on five registers, not a chosen chain.
    _,supports=_code_basis()
    overlap=np.zeros((20,20),dtype=np.int64)
    target=np.zeros_like(overlap)
    for i,S in enumerate(supports):
        for j,T in enumerate(supports):
            for w in S:
                for v in T:
                    if w[1:]==v[:3]:
                        overlap[4*i+v[3],5*w[0]+j]+=1
        for a in range(4):
            target[4*i+a,5*a+i]=len(S)//4
    need(np.array_equal(overlap,target),'Adjacent code isometries have overlap1/4 times a permutation')
    # Actual normalized f0 and an original Gaussian root alpha=sqrt2 e0.
    need((2,0,0,0,0,0,0,0) in _roots(),'Witness is an original root')
    x,y,x2,x3=sp.symbols('x y x2 x3',real=True)
    t=(x+sp.I*y)/sp.sqrt(2)
    a1,a2,a3=sp.sqrt(2)*x,sp.sqrt(2)*x2,sp.sqrt(2)*x3
    response=sp.expand((-4*t**3+6*t**2*a1-4*t*(a2/2+a1**2/2)+a3/3+a1*a2/2+a1**3/6)/sp.sqrt(96))
    coeff=sp.simplify(response.coeff(x3))
    norm2=sp.simplify(3*coeff**2)
    need(norm2==sp.Rational(1,144),'Actual charged partition3 root-zero-mode leakage')
    files=['tfpt_explorer/process.py','tfpt_explorer/cartan_source.py',
        'experiments/theory-contracts/original-admissible-source-closure-20261025/QUANTUM-OCTIC-COMPLETION.txt']
    result={
        'verdict':'EXACT_FINITE_DUALITY_AND_COMPATIBILITY_BOUNDARIES; PHYSICAL_ORIGIN_OPEN',
        'arithmetic':'Gaussian/int64 tensor identities and exact rational/symbolic witness',
        'kernel':'K_ab=<z_a,z_b>^4=F*F',
        'kernel_spectrum':{'3':int(top),'3/2':int(middle),'0':60-int(top+middle)},
        'row_counts':row_counts[0],
        'negative_triangle':{'indices':list(triangle),'roots':[list(rays[i]) for i in triangle],'product':str(product)},
        'magnitude_only_kernel':{'perron_eigenvalue':6,'perron_multiplicity':1,'square_strictly_positive':True},
        'source_metric':{'pseudoinverse':'K+=(7K-2K^2)/9','whitened_frame':'F K+ F*=Pi_sym4','conformal_energy':'L0=4 on all35 one-mode quartics','generalized_problem':'4 K c=E K c; E=4 on rank35 support'},
        'overlap':{'adjacent_isometry_overlap':'(1/4) permutation20','adjacent_intersection':0,'projector_penalty_spectrum':{'3/4':20,'5/4':20,'2':984},'native_parent_pair_lower_bound':'3/8','all_identical_overlapping_blocks':'zero intersection for overlap sizes1,2,3 by symmetry transport; no independently rotated-code classification'},
        'root_zero_mode':{'input':'f0=sum_a(u_a.h(-1))^4/sqrt96','alpha':'sqrt2 e0 in original Gaussian chart','charged_partition3_coefficient':str(coeff),'charged_partition3_norm_squared':str(norm2),'outside_charge_zero_code':True},
        'analytic_consequence':'Any operator commuting with all horizontal E8 zero modes cannot have precisely this C5 as an isolated eigen/ground space in the full grade4 source; applies to the screened W8 zero mode if used as such a selector',
        'generator_nonuniqueness':'L_t=(3I-K)(I+tK),t>=0, shares top-five kernel; remaining eigenvalues3/2*(1+3t/2) and3. This is a finite polynomial-family witness, not full P1 countermodels.',
        'source_sha256':{f:hashlib.sha256((REPO/f).read_bytes()).hexdigest() for f in files},
        'not_claimed':['Actual physical60-site network','Physical identification of Gram with energy','All rotated-code or open-system exclusion','W8 cannot contribute to larger E8 multiplets or state-dependent dynamics','Complete TFPT origin or 3+2/clock/seed selection']
    }
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'kernel_spectrum':result['kernel_spectrum'],'negative_cycle':str(product),'metric_full35':True,'overlap_intersection':0,'root_witness_norm2':str(norm2),'magnitude_control_fires':True},sort_keys=True))

if __name__=='__main__':
    main()
