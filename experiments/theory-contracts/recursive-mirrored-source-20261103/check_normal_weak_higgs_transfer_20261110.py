"""Check the charged normal Higgs/weak transfer on the existing degree-one split."""
import json
from pathlib import Path
import sys

import numpy as np
from mpmath import iv
from scipy.linalg import eigh_tridiagonal
import sympy as sp

sys.path.insert(0,str(Path(__file__).resolve().parents[3]/'verification'))
import v252_full_finite_triple as native

X=np.array([[0,1],[1,0]])
Y=np.array([[0,-1j],[1j,0]])
Z=np.diag([1,-1])
swap=native.Kswap()

def weak_lift(matrix):
    left=native.rep(0,matrix,np.zeros((3,3)))
    return left-swap@left.conj()@swap

weak_charge=weak_lift(Z)/2
charged=(weak_lift(X)+1j*weak_lift(Y))/2
assert np.count_nonzero(weak_charge@charged-charged@weak_charge-charged)==0
assert np.linalg.matrix_rank(charged)==24

center_left=native.rep(1,-np.eye(2),np.eye(3))
center=center_left@swap@center_left.conj()@swap
family=np.arange(1,10).reshape(3,3)+1j*np.eye(3)
majorana=3*np.eye(3)+np.ones((3,3))
reference=native.full_D(*(np.zeros((3,3)),)*4,majorana)
yukawa=native.full_D(family,2*family,3*family,4*family,majorana)-reference
assert np.array_equal(center@reference@center.conj().T,reference)
assert np.count_nonzero(yukawa)>0
assert np.array_equal(center@yukawa@center.conj().T,-yukawa)
assert np.linalg.det(-np.eye(2))==1

x=sp.symbols('x',real=True)
q=sp.Rational

def norm_moment(poly):
    p=sp.Poly(sp.expand(poly),x)
    return sum(coefficient/q(k+1) for (k,),coefficient in p.terms())

for sign in (1,-1):
    alpha,beta=(0,1) if sign==1 else (1,0)
    weight=x if sign==1 else 1-x
    polys=[]
    for n in range(6):
        f=sp.expand(sp.jacobi(n,alpha,beta,2*x-1))
        differential= -x*(1-x)*sp.diff(f,x,2)-(beta+1-3*x)*sp.diff(f,x)+f/2
        assert sp.expand(differential-(q(1,2)+n*(n+2))*f)==0
        polys.append(f/sp.sqrt(norm_moment(weight*f*f)))
    for n,f in enumerate(polys):
        for m,g in enumerate(polys):
            moment=norm_moment(weight*f*g)
            assert sp.simplify(moment-int(m==n))==0
            expected=0
            if m==n:
                expected=1+q(sign,(2*n+1)*(2*n+3))
            if abs(m-n)==1:
                j=min(m,n)
                expected=sp.sqrt((j+1)*(j+2))/(2*j+3)
            assert sp.simplify(norm_moment(weight*f*g*2*x)-expected)==0

assert norm_moment((2*x-1)**2)==q(1,3)
assert norm_moment(4*x*(1-x))==q(2,3)
assert norm_moment((2*x)**6)==q(64,7)
u0,u1=sp.symbols('u0 u1',nonnegative=True)
assert norm_moment((u0+2*x*u1)**2)==u0**2+2*u0*u1+q(4,3)*u1**2


def resolvent_bounds(sign, radius, strength, spectral_parameter, depth):
    """Positive tail bounds propagated through the exact Jacobi recursion."""
    iv.dps=48
    def diagonal(n):
        return (n*(n+2)+iv.mpf(1)/2)/(radius*radius)
    tail=iv.mpf([(1/(spectral_parameter+diagonal(depth)+2*strength)).a,
                (1/(spectral_parameter+diagonal(depth))).b])
    for n in range(depth-1,-1,-1):
        a_squared=iv.mpf((n+1)*(n+2))/(2*n+3)**2
        beta=iv.mpf(1)+iv.mpf(sign)/((2*n+1)*(2*n+3))
        denominator=spectral_parameter+diagonal(n)+strength*beta-strength**2*a_squared*tail
        assert denominator.a>0
        tail=1/denominator
    return tail


def direct(sign, size):
    n=np.arange(size,dtype=float)
    diagonal=n*(n+2)+.5+1+sign/((2*n+1)*(2*n+3))
    edge=np.sqrt((n[:-1]+1)*(n[:-1]+2))/(2*n[:-1]+3)
    energies,vectors=eigh_tridiagonal(diagonal,edge)
    weights=abs(vectors[0])**2
    covariance=np.sum(weights*np.exp(-np.sqrt(energies))/(2*np.sqrt(energies)))
    return float(np.sum(weights/(energies+1))),float(energies[0]),float(covariance)


results={}
for sign in (1,-1):
    bound=resolvent_bounds(sign,1,1,1,10)
    spectral,lowest,covariance=direct(sign,32)
    assert abs(float(bound.a)-spectral)<2e-15
    _,short_low,short_cov=direct(sign,16)
    assert abs(short_low-lowest)<3e-13 and abs(short_cov-covariance)<3e-14
    results[str(sign)]={'R_1_b_1_z_1_resolvent':spectral,
        'certified_interval_width':float(bound.delta),
        'lowest_normal_mass_squared_diagnostic':lowest,
        'projected_euclidean_covariance_at_t_1':covariance}


def log_complex_oscillator_determinant(sign, depth):
    n=np.arange(depth,dtype=float)
    diagonal=n*(n+2)+1.5+sign/((2*n+1)*(2*n+3))
    off_diagonal=np.sqrt((n[:-1]+1)*(n[:-1]+2))/(2*n[:-1]+3)
    eigenvalues=eigh_tridiagonal(diagonal,off_diagonal,eigvals_only=True)
    frequency=np.sqrt(eigenvalues)
    return np.sum(frequency+2*np.log1p(-np.exp(-frequency)))


depth=96
truncated=(log_complex_oscillator_determinant(1,depth)
           -log_complex_oscillator_determinant(-1,depth))
epsilon=1/(4*depth**2)
log_tail=1/(8*(depth-1)**2)+1/(6*(depth-1)**3)
log_upper=truncated/(1-2*epsilon)+log_tail
double=(log_complex_oscillator_determinant(1,2*depth)
        -log_complex_oscillator_determinant(-1,2*depth))
assert 0 < truncated < double < log_upper


def mixed_radial_weight_jet(time_difference):
    n=np.arange(1,2001,dtype=float)
    mass_squared=n*n-.5
    return np.sum(np.exp(-2*np.sqrt(mass_squared)*time_difference)/((4*n*n-1)*mass_squared))


zero_contact=-1-np.pi/(np.sqrt(2)*np.tan(np.pi/np.sqrt(2)))
assert abs(zero_contact-mixed_radial_weight_jet(0))<2e-11

print(json.dumps({'native_weak_positive_root_rank':int(np.linalg.matrix_rank(charged)),
    'weak_center_on_v252_yukawa':-1,
    'weak_center_on_determinant_line':1,
    'free_charged_normal_gap_on_unit_sphere':.5,
    'normal_harmonic_variance_for_same_charged_mode':str(q(1,3)),
    'charged_modes':results,
    'conditional_normal_azimuthal_ratio_at_R_1_b_1_time_1':{
        'relative_log_determinant_lower':float(truncated),
        'relative_log_determinant_upper':float(log_upper),
        'relative_weight_lower':float(np.exp(-log_upper)),
        'relative_weight_upper':float(np.exp(-truncated))},
    'relative_log_weight_mixed_b_jet_at_b_0':{
        'time_0':float(zero_contact),'time_0_point_5':float(mixed_radial_weight_jet(.5)),
        'time_1':float(mixed_radial_weight_jet(1.))}},indent=2))
