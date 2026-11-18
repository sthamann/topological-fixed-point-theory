"""Put the native Higgs insertion between the two existing Pin domains."""

import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss


repo = Path(sys.argv[1]).resolve()
contract = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with contextlib.redirect_stdout(io.StringIO()):
    old = runpy.run_path(str(contract / "check_tfpt_conditioned_r113_weights_20261110.py"))
pin = runpy.run_path(str(repo / ("experiments/theory-contracts/"
    "original-admissible-source-closure-20261025/geometric_source_return/source_geometric_machine.py")))

x = np.array([[0.,1.],[1.,0.]])
y = np.array([[0.,-1j],[1j,0.]])
z_spin = np.diag([1.,-1.])
assert np.array_equal(z_spin@x,-x@z_spin) and np.array_equal(z_spin@y,-y@z_spin)
# Native R032 Z*Phi term in the actual Pin spin chart. Z exchanges
# (c,b) and (-c,-b), so it is a return via the other given domain.
m, v = old["majorana"], old["inner"](np.array([0.,1.]))
gamma_f = old["gamma"]
assert np.max(abs(gamma_f@m+m@gamma_f))<1e-14
assert np.max(abs(gamma_f@v+v@gamma_f))<1e-14

eig_m, u_m = np.linalg.eigh(m)
groups = [np.flatnonzero(abs(eig_m-value)<1e-9) for value in np.unique(eig_m.round(10))]
mass = np.array([eig_m[index].mean() for index in groups])
couplings = abs(u_m.conj().T@v@u_m)**2
weights = np.array([[couplings[np.ix_(a,b)].sum() for b in groups] for a in groups])
np.testing.assert_allclose(weights,weights.T,atol=1e-12)
np.testing.assert_allclose(weights.sum(),np.trace(v@v).real,atol=1e-12)

def rotated_weights(theta):
    link = old["weak_lift"](np.diag(np.exp(1j*theta*np.array([.5,-.5]))))
    np.testing.assert_allclose(link@m@link.conj().T,m,atol=1e-13)
    first = u_m.conj().T@v@u_m
    second = u_m.conj().T@link@v@link.conj().T@u_m
    entries = first*second.T
    return np.array([[entries[np.ix_(a,b)].sum() for b in groups] for a in groups])


def multiplication(modes, left, right):
    # Even extension to [-1,1]; both windows vanish near the physical boundary.
    nodes, quadrature = leggauss(128)
    coordinate = (right+left+(right-left)*nodes)/2
    profile = np.exp(1-1/(1-nodes**2))
    differences = np.arange(-2*modes[-1],2*modes[-1]+1)
    moments = np.cos(np.pi*differences[:,None]*coordinate)@(profile*quadrature*(right-left)/2)
    return moments[(modes[:,None]-modes[None,:])+2*modes[-1]]


def constant_return(lam, occupation, phi, z):
    squares = np.linalg.eigvalsh(phi)**2
    return np.sum(occupation[:,None]*(lam[:,None]+z)
                  /(lam[:,None]**2+squares-z*z))/len(squares)


def cross_jet(lam, u, f, g, z, weights=weights):
    # After the exact chirality transport of the opposite domain, the full
    # Dirac is diag(D,-D)+tau_1*Phi. Keep the original Majorana mass.
    numerator = lam[None,:,None,None]*z_spin + mass[:,None,None,None]*x + z*np.eye(2)
    res = numerator/(lam[None,:,None,None]**2+mass[:,None,None,None]**2-z*z)
    total = 0j
    for a in range(len(mass)):
        left = u.conj()[:,None]*res[a,:,0,:]
        right = u[:,None]*res[a,:,:,0]
        lf, lg = f.T@left, g.T@left
        rf, rg = f@right, g@right
        for b in range(len(mass)):
            middle = x@res[b]@x
            total += (weights[a,b]*np.einsum("ja,jab,jb->",lf,middle,rg)
                      + weights[a,b].conjugate()*np.einsum("ja,jab,jb->",lg,middle,rf))
    return total/len(m)


def independent_small_control():
    # A native massless color block checks the ordered resolvent insertions
    # against direct inversions of both domains at distinct mark amplitudes.
    modes,d,_,_ = pin["dirac_block"](8,.1,np.pi/2)
    f,g = (multiplication(modes,*w) for w in ((.125,.25),(.5,.625)))
    n=len(d)
    h=np.block([[d,np.zeros_like(d)],[np.zeros_like(d),-d]])
    vf,vg=(np.kron(x,mark) for mark in (f,g))
    z=2j
    p=np.eye(2*n)[:,n//2]
    r=np.linalg.inv(h-z*np.eye(2*n))
    expected=p.conj()@(r@vf@r@vg@r+r@vg@r@vf@r)@p
    delta=.002
    measured=sum(a*b*(p.conj()@np.linalg.solve(h+a*delta*vf+b*delta*vg-z*np.eye(2*n),p))
                 for a in (-1,1) for b in (-1,1))/(4*delta**2)
    assert abs(expected-measured)<2e-9
    return float(abs(expected-measured))


out={}
for epsilon in (0.,.1):
    results={}
    for cutoff in (32,64,128):
        modes,d,_,_ = pin["dirac_block"](cutoff,epsilon,np.pi/2)
        lam,u = np.linalg.eigh(d)
        source = u[cutoff,:].conj()
        f,g = (u.conj().T@multiplication(modes,*w)@u for w in ((.125,.25),(.5,.625)))
        zz=2j
        cross = cross_jet(lam,source,f,g,zz)
        # In S the pullback retains the odd field. Test the original even
        # windows as fields on the Pin cover; S is then H=M_R=0 and T copies A.
        absolute = constant_return(lam,abs(source)**2,m+.55*v,zz)
        ref_s = constant_return(lam,abs(source)**2,0*m,zz)
        results[str(cutoff)] = {"cross_resolvent_derivative":[cross.real,cross.imag],
                                "H_0_55_resolvent_A_minus_S":[(absolute-ref_s).real,(absolute-ref_s).imag],
                                "H_0_55_resolvent_A_minus_T":[0.,0.]}
    angles = {}
    for theta in (0.,np.pi,2*np.pi):
        rotated = cross_jet(lam,source,f,g,zz,rotated_weights(theta))
        np.testing.assert_allclose(rotated,np.cos(theta/2)*cross,atol=1e-12)
        angles[str(theta)] = [rotated.real,rotated.imag]
    results["relative_weak_Cartan_rotation"] = angles
    np.testing.assert_allclose(results["64"]["cross_resolvent_derivative"],
                               results["128"]["cross_resolvent_derivative"],atol=2e-9)
    out[str(epsilon)]=results

print(json.dumps({"two_domain_geometry":out,"independent_direct_inverse_error":independent_small_control(),
                  "scope":"Conditional direct sum of the original Pin domains (c,b),(-c,-b) with the actual native order-zero Higgs/Majorana block. Resolvent at z=2i, not a physical OS return or a choice of P_adm."},indent=2))
