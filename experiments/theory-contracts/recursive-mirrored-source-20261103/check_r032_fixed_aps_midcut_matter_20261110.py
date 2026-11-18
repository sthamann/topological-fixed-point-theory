"""Actual fixed APS exterior cap, with a matter mark inserted before gluing."""
import contextlib
import io
import json
from pathlib import Path
import runpy

import numpy as np
from scipy.linalg import eigh, solve

base = Path(__file__).resolve().parent
with contextlib.redirect_stdout(io.StringIO()):
    s = runpy.run_path(str(base / 'check_r113_mixed_metric_matter_20261110.py'))
I, grade, right = (s[k] for k in ('I','grade','right'))
n, r = len(I),len(right)
minus = np.diag(grade)<0
half = s['product'](s['factors'][:3])
E = tuple(half[:n,j*n:(j+1)*n][:,minus] for j in range(4))


def prod(x,y):
    a,b,c,d=x
    w,z,v,u=y
    return a@w,b@w+a@z,c@w+a@v,d@w+b@v+c@z+a@u


def inverse(x):
    A, m,h,mh=x
    G = np.linalg.inv(A)
    return G,-G@m@G,-G@h@G,G@m@G@h@G+G@h@G@m@G-G@mh@G


Et = tuple(x.conj().T for x in E)
A = prod(Et,E)
W = prod(prod(E,inverse(A)),Et)
C = tuple(x[np.ix_(right,right)] for x in W)
t = lambda M: np.array([np.trace(x).real for x in M])
tr = t(C)
paired = np.array([np.zeros(4),2*tr-np.array([r,0,0,0]),2*(tr-t(prod(C,C)))]).T


def cap(sigma,h,sign):
    R = I.astype(complex)
    p = sign*s['p0']*np.exp(-sigma)
    for H, width in zip(s['history'][:3],s['durations'][:3]):
        phi = s['inner_form']((1+h)*np.asarray(H))+s['M']
        d,U=eigh(p*grade-1j*grade@phi)
        R = ((U*np.exp(-width*d))@U.conj().T)@R
    V = R[:, sign*np.diag(grade)<0]
    overlap=V.conj().T@V
    G=V@solve(overlap,V.conj().T,assume_a='pos')
    return G, np.linalg.slogdet(overlap)[1]


def direct(sigma,h):
    P,l = cap(sigma,h,1)
    Q,l2 = cap(sigma,h,-1)
    assert np.max(abs(P+grade@Q@grade-I)) < 1e-12
    assert abs(l-l2)<1e-12
    a,b = (x[np.ix_(right,right)] for x in (P,Q))
    e = np.trace(a+b).real-r
    v = np.trace(a-a@a+b-b@b).real + e*e
    return np.array([e,np.trace(a-b).real,v])


compare=s['central_jet'](direct)
err = np.max(abs(paired-compare))
assert err < 4e-7
assert abs(paired[3,1]) > 1e-4 and abs(paired[3,2]) > 1e-4
print(json.dumps({'columns':['centered mirror-even occupation','oriented occupation',
                             'square of centered mirror-even occupation'],
                  'rows':['value','metric','Higgs','mixed'],
                  'fixed_aps_cap_jet':paired.tolist(),
                  'independent_finite_difference':compare.tolist(),
                  'max_difference':float(err)},indent=2))
