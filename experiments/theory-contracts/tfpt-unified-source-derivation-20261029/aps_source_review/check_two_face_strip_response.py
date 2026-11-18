"""Exact algebra and finite controls for the original R074 two-face graph.

The infinite proof is in TWO-FACE-STRIP-RESPONSE-20261102.txt. Finite sums
are controls, not an independent proof of the R066 theta identity.
"""
from pathlib import Path
import sys
import mpmath as mp
import sympy as s
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'experiments/theory-contracts/seam-collar-multifield-20261016'))
from boundary_amplitude import graph_projection

checks = []
def check(ok, label):
    if not ok:
        raise AssertionError(label)
    checks.append(label)

t, r = s.symbols('t r', positive=True)
G = graph_projection(s.Matrix([[1/t]]))
expected = s.Matrix([[t*t,t],[t,1]])/(1+t*t)
check(all(s.simplify(v) == 0 for v in G-expected), 'native original graph')
C, Cp = 1/(1+t*t), 1/(1+t*t/r**2)
X, Xp = t/(1+t*t), (t/r)/(1+t*t/r**2)
delta = Cp-C
check(s.factor((C-Cp)**2+X**2+Xp**2-(r*r+1)/(r*r-1)*delta) == 0,
      'pointwise diagonal telescoping identity')
check(s.factor(2*X*Xp-2*r/(r*r-1)*delta) == 0,
      'pointwise cross-face telescoping identity')
P = s.Matrix([[1,-1],[-1,1]])/2
check(P*P == P and P.T == P, 'normalized relative projector')
check(s.Matrix([[1,-1],[-1,1]])/(2*s.pi)-s.Matrix([[1,-1],[-1,1]])/(4*s.pi)
      == P/(2*s.pi), 'zero-mode rank-one coefficient')

mp.mp.dps = 70
def c(p,beta): return 1/(1+mp.exp(beta*p))
def x(p,beta): return 1/(2*mp.cosh(beta*p/2))
def response(alpha,beta,k,N=28):
    diag = cross = mp.mpf('0')
    for n in range(-N,N+1):
        a,b = c(n+alpha,beta),c(n+alpha-k,beta)
        u,v = x(n+alpha,beta),x(n+alpha-k,beta)
        diag += (a-b)**2+u**2+v**2
        cross -= 2*u*v
    return diag,cross
for alpha in (mp.mpf(1)/4, mp.mpf(3)/4):
    beta=4*mp.pi
    for k in (-3,-2,-1,1,2,3):
        d,o=response(alpha,beta,k)
        K=abs(k)
        check(abs(d-K*mp.coth(beta*K/2)) < mp.mpf('1e-55'), f'diag alpha={alpha} k={k}')
        check(abs(o+K/mp.sinh(beta*K/2)) < mp.mpf('1e-55'), f'cross alpha={alpha} k={k}')
        check(abs(2*(d+o)-2*K*mp.tanh(beta*K/4)) < mp.mpf('1e-55'), f'R074 common response alpha={alpha} k={k}')
    d,o=response(alpha,beta,0)
    check(abs(d-1/(4*mp.pi)) < mp.mpf('1e-55'), f'R066 quarter zero mode alpha={alpha}')
    check(abs(d+o) < mp.mpf('1e-65'), f'zero common mode alpha={alpha}')

# Independent direct harmonic energy for arbitrary complex boundary data.
L=2*mp.pi
for k in (1,2,3):
    a=mp.mpc('0.2','0.3'); b=mp.mpc('-0.1','0.4')
    f=lambda y: (a*mp.sinh(k*(L-y))+b*mp.sinh(k*y))/mp.sinh(k*L)
    df=lambda y: k*(-a*mp.cosh(k*(L-y))+b*mp.cosh(k*y))/mp.sinh(k*L)
    energy=mp.quad(lambda y: abs(df(y))**2+k*k*abs(f(y))**2,[0,L])
    rhs=k*mp.coth(k*L)*(abs(a)**2+abs(b)**2)-2*k/mp.sinh(k*L)*mp.re(a.conjugate()*b)
    check(abs(energy-rhs) < mp.mpf('1e-60'), f'harmonic strip energy k={k}')

# Independent diagnostic on the full coupled two-component graph. This is a
# fixed-background source with the same quarter, not an identified P1 map.
def coupled_response(m,k,N):
    beta=float(4*mp.pi); eye=np.eye(2)
    sigma_z=np.diag([1.,-1.]); sigma_x=np.array([[0.,1.],[1.,0.]])
    def blocks(p):
        E=np.sqrt(p*p+m*m)
        h=p*sigma_z+m*sigma_x
        C=(eye-h*np.tanh(beta*E/2)/E)/2
        X=eye/(2*np.cosh(min(beta*E/2,350)))
        return C,X
    diag=cross=0.; projection_error=0.
    for n in range(-N,N+1):
        C,X=blocks(n+.25); Cp,Xp=blocks(n+.25-k)
        diff=C-Cp
        diag+=np.trace(diff@diff+X@X+Xp@Xp)
        cross-=2*np.trace(X@Xp)
        projection_error=max(projection_error,np.linalg.norm(C@C-C+X@X))
    check(projection_error<1e-14,f'coupled block projection identity m={m} k={k}')
    # With j0=N-|k| and |p|,|p-k| >= j0 in the omitted tails,
    # ||C'(p)||_HS <= (m/p^2+4 beta exp(-beta |p|))/sqrt(2).
    # Summing a monotone integral bound on j^-4 and geometric tails gives:
    j0=N-abs(k)
    assert j0>0
    geom=lambda r: np.exp(-r*j0)/(1-np.exp(-r))
    bound=2*k*k*(m*m*(j0**-4+1/(3*j0**3))+16*beta*beta*geom(2*beta))+8*geom(beta)
    return float(diag),float(cross),float(bound)

for k in (0,1,2):
    massless=coupled_response(0.,k,80)
    massive=coupled_response(.3,k,80)
    # Even allowing the (individually rigorous) positive tail bounds, the
    # diagonal responses differ. No fit or extrapolation is involved.
    check(abs(massless[0]-massive[0])>massless[2]+massive[2]+1e-10,
          f'coupled mass changes whole response k={k}')
    print(f'coupled k={k} massless={massless[0]:.12f},{massless[1]:.12f} '
          f'massive={massive[0]:.12f},{massive[1]:.12f} '
          f'diag_tail_bounds={massless[2]:.3g},{massive[2]:.3g}')

print(f'PASS {len(checks)} checks; infinite free result uses the analytic telescoping proof and existing R066 identity; coupled example has explicit tail bound; no physical source identification')
