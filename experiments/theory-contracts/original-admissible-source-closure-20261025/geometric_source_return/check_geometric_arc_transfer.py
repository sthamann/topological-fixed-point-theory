"""Same geometric source: second-order exterior resolvent has arc transfer.

The nonlocal Schur term is a spectral response at z=0, not an autonomous
physical generator or the missing native quartic/P1 field identification.
"""
import json
from pathlib import Path
import argparse

import numpy as np
import sympy as sp

from source_geometric_machine import profile_values


def need(condition,message):
    if not condition:raise RuntimeError(message)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'runs/geometric_arc_transfer.json')
    args=parser.parse_args()
    k,a,z,hm=sp.symbols('k a z h_m',real=True)
    for b in [-1,1]:
        D=b*sp.Matrix([[0,k+sp.I*a],[k-sp.I*a,0]])
        B=-b*hm*sp.Matrix([k+sp.I*a/2,k-sp.I*a/2])/2
        actual=(sp.conjugate(B).T*(z*sp.eye(2)-D).inv()*B)[0]
        expected=hm**2*(z*(2*k*k+a*a/2)+b*(2*k**3+3*k*a*a/2))/(4*(z*z-k*k-a*a))
        need(sp.factor(actual-expected)==0,'Actual exterior source Schur block')
        static=-b*hm**2*k/2+b*hm**2*a*a*k/(8*(k*k+a*a))
        need(sp.factor(actual.subs(z,0)-static)==0,'Local plus nonlocal response split')
    # Same collar-preserving diagnostic geometry as the executed source.
    nodes,weights=np.polynomial.legendre.leggauss(1024)
    w=(nodes+1)/2;weights=weights/2
    h,_=profile_values(w,'collar_bump')
    m=np.arange(1,129);energies=np.pi*m
    coefficients=2*np.cos(energies[:,None]*w[None,:])@(weights*h)
    # Proper preimages of the base interval (1/4,3/4), with disjoint closures.
    length=1/16
    r=np.exp(-energies)
    positive_weights=(np.exp(-energies/2)*(1+r)/(1+r*r)
                      *np.sinh(energies*length/2)**2/(4*length))
    analytic=(1-1j)*np.sum(coefficients**2*positive_weights)
    need(abs(analytic)>1e-6,'Opposite native-arc transfer absent')
    trials=[]
    for nmax in [64,128,256,512,1024]:
        momenta=2*np.pi*(np.arange(-nmax,nmax+1)+.25)
        symbol=momenta/8*np.sum(coefficients[:,None]**2*energies[:,None]**2/
                  (momenta[None,:]**2+energies[:,None]**2),axis=0)
        value=np.sum(length*np.sinc(momenta*length/(2*np.pi))**2*
                     np.exp(-.5j*momenta)*symbol)
        trials.append(dict(longitudinal_bound=nmax,value=[float(value.real),float(value.imag)],
                           difference_from_analytic=float(abs(value-analytic))))
    need(trials[-1]['difference_from_analytic']<2e-11,'Independent Fourier/images source test')
    result=dict(actual_source='Original b=+1,c=+i quarter circle, same normal conformal probe',
                profile='Same collar_bump as source_geometric_machine',
                order='epsilon^2 coefficient of exterior Schur response at z=0',
                operator_symbol='R_h(k)=k/8 sum_m h_m^2 (pi*m)^2/(k^2+(pi*m)^2)',
                arcs={'I':[1/32,3/32],'J':[17/32,19/32],'smears':'normalized indicator functions on proper native preimage arcs; bounded nonlocal term only'},
                analytic_cross_arc_amplitude=[float(analytic.real),float(analytic.imag)],
                all_profile_fact='For q=i and opposite equal arcs each nonzero cosine coefficient contributes a positive multiple of (1-i); no cancellation for any nonconstant smooth even profile',
                compactness='The nonlocal remainder is compact: |R_h(k)| <= sum_m h_m^2 (pi*m)^2/(8 |k|). It cannot itself equal a nonzero folded local multiplication current.',
                trials=trials,
                limitation='An actual inter-arc transfer is necessary data; no native quartic composite/P1 insertion identification or profile selection is proved',
                physical_gates_closed=[])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
