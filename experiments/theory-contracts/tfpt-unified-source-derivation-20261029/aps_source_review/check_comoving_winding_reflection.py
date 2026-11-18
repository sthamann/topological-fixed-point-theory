#!/usr/bin/env python3
"""Finite local-window checks of the exact polarized winding/reflection identities."""
import numpy as np

I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
P = np.diag([1., 0.])
M = I - P

def w(p, mass):
    t = np.arctan(mass / p) / 2
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])

def check(alpha, mass, nmax=12):
    ns = np.arange(-nmax, nmax)
    d = 2 * len(ns)
    U = np.zeros((d, d)); Ur = np.zeros_like(U)
    B = np.zeros_like(U); Br = np.zeros_like(U); R = np.zeros_like(U)
    for j, n in enumerate(ns):
        wp = w(n + alpha, mass)
        wref = w(n + 1-alpha, mass)
        B[2*j:2*j+2, 2*j:2*j+2] = wp @ (P if n < 0 else M) @ wp.T
        Br[2*j:2*j+2, 2*j:2*j+2] = wref @ (P if n < 0 else M) @ wref.T
        U[2*j:2*j+2, 2*j:2*j+2] = wp @ M @ wp.T
        Ur[2*j:2*j+2, 2*j:2*j+2] = wref @ P @ wref.T
        jr = int(-n-1+nmax)
        R[2*jr:2*jr+2,2*j:2*j+2]=X
        if j + 1 < len(ns):
            U[2*(j+1):2*(j+1)+2,2*j:2*j+2] = w(n+1+alpha,mass) @ P @ wp.T
        if j > 0:
            Ur[2*(j-1):2*(j-1)+2,2*j:2*j+2] = w(n-1+1-alpha,mass) @ M @ wref.T
    interior = np.r_[np.arange(2,d-2)]
    target=np.zeros_like(U)
    j0=nmax
    target[2*j0:2*j0+2,2*j0:2*j0+2]=w(alpha,mass) @ P @ w(alpha,mass).T
    test=(U@B@U.T-B-target)[np.ix_(interior,interior)]
    # The matrices are real, so R also represents the linear part of the chosen
    # antiunitary conjugation. Projection excludes finite-cutoff edge effects.
    assert np.max(np.abs(test)) < 2e-14
    assert abs(np.trace(target)-1) < 2e-14
    assert np.max(np.abs((R@B@R.T-Br)[np.ix_(interior,interior)])) < 2e-14
    assert np.max(np.abs((R@U@R.T-Ur)[np.ix_(interior,interior)])) < 2e-14
    # Winding cannot itself be an order-four quarter map: |1,+> goes to |5,+>.
    j=nmax+1
    v=np.zeros(d);v[2*j:2*j+2]=w(1+alpha,mass)[:,0]
    v4=np.linalg.matrix_power(U,4) @ v
    assert abs(np.vdot(v,v4)) < 2e-14
    assert abs(np.linalg.norm(v4)-1) < 2e-14
    print(f'alpha={alpha}, m={mass}: rank-one interior shift, reflection and non-quarter verified')

for a in (.25,.75):
    for m in (0., .3, 2.2746272169):
        check(a,m)
Y = np.array([[0.,-1j],[1j,0.]])
Z = np.diag([1.,-1.])
for p in (.25, -.75, 1.25):
    m, mdot, eps = .3, .7, 1e-6
    W = w(p,m)
    Wdot=(w(p,m+eps*mdot)-w(p,m-eps*mdot))/(2*eps)
    actual=W.T @ (p*Z+m*X) @ W - 1j*W.T @ Wdot
    expected=np.sign(p)*np.hypot(p,m)*Z-p*mdot/(2*(p*p+m*m))*Y
    assert np.max(np.abs(actual-expected)) < 1e-9
    assert abs(np.exp(2j*np.pi*np.hypot(p,m))-np.exp(2j*np.pi*abs(p))) > 1e-4
print('PASS: dynamical connection and conditional quarter-time mismatch, 3 mode samples')
print('PASS: 6 source-block cases; finite-window checks do not prove infinite-domain claims')
