#!/usr/bin/env python3
"""Actual native source: relative-arc charge test and modular diagnostic.

No extra state, source weighting, coupling, spatial lattice or physical clock.
The finite ring is a regulator diagnostic of the retained one-circle sea.
"""
import hashlib
import argparse
import itertools as it
import json
from pathlib import Path
import sys

import numpy as np
import sympy as s

REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
SOURCE=REPO/'experiments/theory-contracts/joint-native-algebra-selection-20261016'
sys.path[:0]=[str(SOURCE),str(REPO)]
import check_native_quartic_multilocal_support as source
old=source.old
m=source.m

def need(ok,message):
    if not bool(ok):raise RuntimeError(message)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=Path(__file__).with_name('arc_current_completion_results.json'))
    args=parser.parse_args()
    need(hashlib.sha256(source.M_PATH.read_bytes()).hexdigest()==source.M_SHA,'Original M source')
    ld=m._lattice_data()
    T=s.Matrix.hstack(*map(s.Matrix,m._TARGET_SIMPLE))
    O=old.cl(T*ld['basis'].inv()/s.sqrt(2))
    M=s.Matrix(json.loads(source.M_PATH.read_text())['explicit_Weyl_matrix'])
    U=s.zeros(8,4)
    for a in range(4):U[2*a,a]=1/s.sqrt(2);U[2*a+1,a]=-s.I/s.sqrt(2)
    V=old.cl(M.T*O*U)
    need(old.cl(V.conjugate().T*V)==s.eye(4),'Original metric')
    need(old.cl(V.T*V)==s.zeros(4),'Original isotropy')
    polys,norms=m.code()
    subsets,C=source.principal_coefficients(polys,V)
    need(C.rank()==5,'Actual principal map')
    # For each original arc r, E_r=U_alpha^*|r><r|U_alpha.
    # On exterior^4 C8 the offdiagonal entries of dGamma(E_r) have
    # modulus 1/8 exactly when two four-subsets differ by one index.
    # The fourth-current highest-CAR coefficient is diag(C[:,i]).
    # Hence this Gram is the EXACT [N_r,T_i] principal coefficient Gram.
    edges=[(a,b) for a,K in enumerate(subsets) for b,L in enumerate(subsets)
           if len(set(K).symmetric_difference(L))==2]
    need(len(edges)==1120,'Johnson graph directed edge count')
    d=s.Matrix([list(C[a,:]-C[b,:]) for a,b in edges])
    G=old.cl(d.conjugate().T*d/64)
    rank=G.rank()
    principal_gram=old.cl(C.conjugate().T*C)
    total=s.Matrix([sum(C[:,i]) for i in range(5)])
    N=s.diag(*[1/s.sqrt(n) for n in norms])
    Gn=old.cl(N*G*N)
    # Direct agreement of one arc using the actual 8-Fourier phase,
    # without 70-by-70 floating tensor diagonalization. The phase for
    # r=0 is 1 and wedge replacement signs are +/-1.
    B=s.zeros(70,70)
    for a,K in enumerate(subsets):
        for j,k in enumerate(K):
            for ell in range(8):
                if ell in K:continue
                dest=list(K);dest[j]=ell
                inversions=sum(dest[x]>dest[y] for x in range(4) for y in range(x+1,4))
                b=subsets.index(tuple(sorted(dest)))
                B[b,a]=s.Rational((-1)**inversions,8)
        B[a,a]=s.Rational(1,2)
    direct=s.Matrix([[sum(s.conjugate(B[a,b]*(C[b,i]-C[a,i]))*
                          B[a,b]*(C[b,j]-C[a,j])
                          for a,b in edges) for j in range(5)] for i in range(5)])
    need(old.cl(direct)==G,'Direct actual arc zero exterior commutator')
    need(G==G.conjugate().T,'Hermitian charge witness')
    need(rank==5,'Every nonzero native quartic has relative-arc charge response')
    need(G.det()==s.Rational(498159,262144),'Exact reconstructed charge determinant')
    need(all(G[:k,:k].det()>0 for k in range(1,6)),'Positive definite exact witness')
    difference_classes=[]
    for step in range(1,5):
        chosen=[]
        for a,b in edges:
            removed=next(iter(set(subsets[a])-set(subsets[b])))
            added=next(iter(set(subsets[b])-set(subsets[a])))
            if (added-removed)%8 in {step,8-step}:chosen.append((a,b))
        ds=s.Matrix([list(C[a,:]-C[b,:]) for a,b in chosen])
        gs=old.cl(ds.conjugate().T*ds)
        need(gs.rank()==5,'Every nonconstant real relative phase detects every native quartic')
        difference_classes.append({'step':step,'directed_edges':len(chosen),'rank':gs.rank()})
    # Source-independent negative control: total density polynomial
    # e4(n_0,...,n_7) has constant principal vector and no relative charge.
    control=s.ones(70,1)
    dc=s.Matrix([control[a]-control[b] for a,b in edges])
    need(dc==s.zeros(1120,1),'Total-density invariant control')

    # Original local scalar one-particle potential is diagonal in the arc
    # chart. The positive native Cartan is A(v)=U_alpha diag(Vv)U_alpha^*.
    # Offdiagonal HS Gram is V*V - V*11^TV/8; independent of alpha.
    linear=old.cl(V.conjugate().T*V-V.conjugate().T*s.ones(8,8)*V/8)
    need(linear.rank()==4,'No positive native Cartan tangent is an arc-local scalar potential')

    diagnostics=[]
    for length in (64,128):
        # Eight disjoint equal arcs, each angle width pi/16, fixed as N grows.
        per=length//32
        ids=np.array([r*length//8+j for r in range(8) for j in range(per)])
        theta=2*np.pi*ids/length
        modes=np.arange(-length//2,0)
        for alpha in (0.25,0.75):
            F=np.exp(1j*np.outer(theta,modes+alpha))/np.sqrt(length)
            cov=F@F.conj().T
            lam,ev=np.linalg.eigh(cov)
            need(lam.min()>0 and lam.max()<1,'Finite faithful restriction without eigenvalue clipping')
            h=(ev*np.log((1-lam)/lam))@ev.conj().T
            mask=np.repeat(np.arange(8),per)
            off=np.where(mask[:,None]!=mask[None,:],h,0)
            L=np.diag(np.repeat([1.,0.,0.,0.,0.,0.,0.,0.],per))
            comm=h@L-L@h
            need(np.linalg.norm(comm)>1e-5,'Joint-source modular flow transfers arc charge')
            product=np.where(mask[:,None]==mask[None,:],cov,0)
            le,ve=np.linalg.eigh(product)
            hp=(ve*np.log((1-le)/le))@ve.conj().T
            need(np.linalg.norm(hp@L-L@hp)<1e-10,'Erasing cross-arc correlations removes this witness')
            # Unchanged actual fold matrix, not a fitted eigenvector chart.
            rr=np.arange(8)
            UF=np.exp(2j*np.pi*np.outer(rr,np.arange(8)+alpha)/8)/np.sqrt(8)
            native_cartan_errors=[]
            for a in range(4):
                coeff=np.array([complex(V[k,a]) for k in range(8)])
                A=np.kron((UF*coeff)@UF.conj().T,np.eye(per))
                native_cartan_errors.append(float(np.linalg.norm(h@A-A@h)))
            need(max(native_cartan_errors)<1e-9,'Joint modular matrix commutes with all original native Cartan matrices')
            diagnostics.append({'ring_size':length,'alpha':str(alpha),
                'retained_sites':len(ids),'occupied_source_modes':len(modes),'covariance_range':[float(lam.min()),float(lam.max())],
                'cross_arc_modular_hs_norm':float(np.linalg.norm(off)),
                'arc0_modular_commutator_hs_norm':float(np.linalg.norm(comm)),
                'product_covariance_commutator':float(np.linalg.norm(hp@L-L@hp)),
                'native_cartan_commutators':native_cartan_errors})
    result={'verdict':'EXACT_RELATIVE_ARC_CHARGE_OBSTRUCTION; FINITE_MODULAR_MIXING_OBSERVED; PHYSICAL_SOURCE_MAP_OPEN',
        'original_source_rank':5,'directed_subset_edges':len(edges),
        'charge_commutator_gram':source.rows(G),'charge_gram_det':str(G.det()),
        'normalized_charge_gram':source.rows(Gn),'normalized_charge_gram_det':str(Gn.det()),
        'charge_gram_rank_each_arc':[rank]*8,'coefficient_column_sums':[str(x) for x in total],
        'relative_phase_difference_classes':difference_classes,
        'positive_native_cartan_vs_local_scalar_gram':source.rows(linear),
        'positive_native_cartan_vs_local_scalar_rank':linear.rank(),
        'negative_control_constant_principal_vector':'zero relative-arc charge',
        'modular_diagnostic':diagnostics,
        'exact_scope':'Same M/O/quarter source, highest CAR coefficient, additive density/stress algebra on the eight proper disjoint preimage arcs.',
        'not_claimed':['P1 or full local CAR is exhausted by density/stress algebra','modular parameter is physical time','a modular matrix entry is a physical coupling','continuum convergence from two ring sizes','physical 3+2, alpha or 3+1D selected'],
        'source_sha256':{str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [source.M_PATH,SOURCE/'check_native_quartic_multilocal_support.py',SOURCE/'ORIGINAL-ORDERED-CURRENT-KERNEL-AND-SEAM-GRADING.txt']}}
    out=args.out
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'rank':rank,'det':str(G.det()),'normalized_det':str(Gn.det()),'linear_rank':linear.rank(),'modular_diagnostic':diagnostics},indent=2))

if __name__=='__main__':main()
