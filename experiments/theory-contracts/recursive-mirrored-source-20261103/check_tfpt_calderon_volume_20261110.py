"""Check the real pullback metric of the actual full96 R032 graph."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

repo = Path(sys.argv[1]).resolve()
folder = repo/"experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    src = runpy.run_path(str(folder/"check_tfpt_conditioned_r113_weights_20261110.py"))

gamma = src["gamma"]
M = -1j*gamma@src["majorana"]
fields = [(1,0),(1j,0),(0,1),(0,1j)]
V = np.array([-1j*gamma@src["inner"](np.array(a)) for a in fields])
weak = src["weak_lift"]
blocks = [i for i,_ in src["blocks"]]
ell = 2*np.pi
ray = V[2]  # Actual existing H=(0,R) direction
assert gamma.shape == (96,96)
assert max(np.linalg.norm(x-x.conj().T) for x in [*V,M]) < 1e-13
assert max(np.linalg.norm(gamma@x+x@gamma) for x in [*V,M]) < 1e-13
for x in [gamma,*V,M]:
    y = x.copy()
    for i in blocks:
        y[np.ix_(i,i)] = 0
    assert np.linalg.norm(y) < 1e-13
assert np.linalg.matrix_rank(V.reshape(4,-1)) == 4
masses = np.linalg.eigvalsh(ray)
assert np.min(np.abs(masses)) > .05
assert np.sum(masses>0) == 48

for z in [np.array([1,0]),np.array([1j,1])/np.sqrt(2),
          np.array([.6j,.8]),np.array([.7,.5+.4j])/np.sqrt(.9)]:
    z = z/np.linalg.norm(z)
    a,b = z
    U = weak(np.array([[b.conjugate(),a],[-a.conjugate(),b]]))
    np.testing.assert_allclose(U@ray@U.conj().T, -1j*gamma@src["inner"](z), atol=2e-13)
    np.testing.assert_allclose(U@M@U.conj().T,M,atol=2e-13)
    np.testing.assert_allclose(U@gamma@U.conj().T,gamma,atol=2e-13)
phase = V[3]
assert np.linalg.norm(phase@ray+ray@phase) < 1e-12
np.testing.assert_allclose(phase@phase,ray@ray,atol=2e-13)
packed = [(gamma[np.ix_(i,i)],M[np.ix_(i,i)],V[:,i][:,:,i]) for i in blocks]


def dd(x,y,dy):
    dx = x[:,None]-x[None,:]
    diff = y[:,None]-y[None,:]
    diagonal = abs(dx)<1e-11
    divided = np.divide(diff,dx,out=np.zeros_like(dx),where=~diagonal)
    divided[diagonal] = ((dy[:,None]+dy[None,:])/2)[diagonal]
    return divided


def metric(p, R, sign=False, mr=True):
    out = np.zeros((4,4))
    for grad,mass,v in packed:
        e,U = np.linalg.eigh(p*grad+R*v[2]+(mass if mr else 0))
        X = (U.conj().T@v@U).reshape(4,-1)
        if sign:
            matrices = [dd(e,np.sign(e),np.zeros_like(e))]
        else:
            z = ell*e
            exp_z = np.exp(-np.abs(z))
            sech = 2*exp_z/(1+exp_z**2)
            tanh = np.tanh(z)
            matrices = [dd(e,tanh,ell*sech**2),dd(e,sech,-ell*sech*tanh)]
        for d in matrices:
            Y = X*d.reshape(-1)
            out += .5*np.real(Y.conj()@Y.T)
    return out


def circle(R,N,alpha=.25):
    return sum(metric(n+alpha,R) for n in range(-N,N+1))


def volume(G):
    sign, log = np.linalg.slogdet(G)
    assert sign>0
    return np.exp(log/2)

# First independent finite difference of the actual two-face projector.
def graph(p,H):
    h = p*gamma+M+np.einsum("a,aij->ij",H,V)
    e,U = np.linalg.eigh(h)
    z = ell*e
    u = np.exp(-np.abs(z))
    s = (U*(2*u/(1+u*u)))@U.conj().T
    t = (U*np.tanh(z))@U.conj().T
    I = np.eye(len(h))
    return np.block([[I+t,s],[s,I-t]])/2

R = .55
H = np.array([0.,0.,R,0.])
step = 1e-5
D = np.array([(graph(.25,H+step*v)-graph(.25,H-step*v)).ravel()/(2*step)
              for v in np.eye(4)])
err = np.max(abs(np.real(D.conj()@D.T)-metric(.25,R)))
assert err<1e-6

# For a fixed cutoff the radial derivative has an extra power of R.
fixed = []
for R in (16,32,64,128,256,512,1024):
    G = circle(R,2)
    fixed.append([R,volume(G),R**5*volume(G)])
assert abs(fixed[-1][2]/fixed[-2][2]-1)<.01

# The continuum-in-momentum constant follows from the original massive
# quarter circle. This is not a regulator chosen for the joint theory.
x,w = leggauss(260)
angle = x*np.pi/2
t = np.tan(angle)
weights = w*np.pi/2/np.cos(angle)**2
L = sum(wt*metric(p,1,sign=True,mr=False) for p,wt in zip(t,weights,strict=True))
m_radial = np.pi/2*np.sum(masses[masses>0])
assert abs(L[2,2]-m_radial)<1e-7
assert np.min(np.linalg.eigvalsh(L))>0
assert abs(L[3,3]-2*m_radial)<1e-7
tower = []
for R,N in [(4,100),(8,200),(16,400),(32,800)]:
    G = circle(R,N)
    tower.append([R,N,R**2*volume(G),float(R*G[2,2])])
# Finite ratio N/R=25 still omits a calculable large-momentum tail.
x,w = leggauss(640)
L25 = sum(25*wt*metric(25*p,1,sign=True,mr=False) for p,wt in zip(x,w,strict=True))
assert abs(tower[-1][2]/volume(L25)-1)<.02
assert abs(tower[-1][3]/L25[2,2]-1)<.02

# Positive powers and traces refer to the existing sea-corrected free
# reference. They only evaluate the conditional replacement factor.
# Radial total-mass coefficient of the original symmetric cutoff N.
# These quadratures are comparisons; the convergence proof is analytic.
def halfline_rule(size):
    x,w = leggauss(size)
    a = (x+1)*np.pi/4
    return np.tan(a), w*np.pi/4/np.cos(a)**2

xs,xweights = halfline_rule(44)
ts,tweights = leggauss(80)
mass_lim = 0.
quot_lim = 0.
for x,wx in zip(xs,xweights,strict=True):
    a = np.arctan(1/x)
    t = np.tan(a*ts)
    T = sum(wt*a/np.cos(a*z)**2*metric(p,1,sign=True,mr=False)
            for p,z,wt in zip(t,ts,tweights,strict=True))
    mass_lim += wx*x*volume(T)*2*np.pi**2
    inv = np.linalg.inv(T)
    quot_lim += wx/np.sqrt(x*inv[2,2])

mass_controls = []
for N in (3,7,15):
    total, quotient = 0.,0.
    for x,wx in zip(xs,xweights,strict=True):
        R=N*x
        G=circle(R,N)
        total += wx*x**3*N**2*volume(G)*2*np.pi**2
        quotient += wx*np.sqrt(N/np.linalg.inv(G)[2,2])
    mass_controls.append([N,total,quotient])
# The source diagnostic is [H, |H|²] applied to a test section that
# is constant near H. Scalar a cancels in the ratio printed below.
def radial_log_slope(R,N):
    eps=.001
    return (np.log(volume(circle(R*np.exp(eps),N)))-
            np.log(volume(circle(R*np.exp(-eps),N))))/(2*eps)

finite_slope = radial_log_slope(512,2)
large_tower_slope = radial_log_slope(8,800)
assert abs(finite_slope+5)<.002
assert abs(large_tower_slope+2)<.06
alpha=.25
t_ref=4*np.pi
ref_log=96*sum(np.log1p(np.exp(-t_ref*abs(n+alpha))) for n in range(-32,33))
print(json.dumps({
    "fixed_mode_cutoff":fixed,
    "tower_scaling_R_N_R2rho_Rgrr":tower,
    "full_limiting_metric_eigenvalues":np.linalg.eigvalsh(L).tolist(),
    "full_limiting_R2rho":volume(L),
    "limiting_R2rho_for_N_over_R_25":volume(L25),
    "exact_radial_limit_formula":float(m_radial),
    "mass_over_N_squared_rough_quadrature_limit":float(mass_lim),
    "quotient_mass_over_sqrt_N_rough_quadrature_limit":float(quot_lim),
    "mass_over_N_squared_and_quotient_mass_over_sqrt_N":mass_controls,
    "pure_S_heat_tail_ratio":float(np.exp(-ref_log)),
    "minimum_Yukawa_mass_for_original_ray":float(np.min(np.abs(masses))),
    "full_graph_derivative_independent_max_error":float(err),
    "r_times_log_density_derivative_at_fixed_N":float(finite_slope),
    "r_times_log_density_derivative_near_whole_tower":float(large_tower_slope),
    "commutator_contact_in_a_units_at_fixed_N":float(-8-2*finite_slope),
    "commutator_contact_in_a_units_near_whole_tower":float(-8-2*large_tower_slope),
    "primitive_F1_computed":False
},indent=2))
