"""Read what the actual v330 transfer proves and test the v221/v369 Petz inference."""
from itertools import product
from pathlib import Path
import sys

import numpy as np
import sympy as sp

repo = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(repo / 'verification'))
from v330_qgamb_admissible_measure import symmetric_transfer

Q = sp.Rational
v = sp.Matrix([1,-1,0])
w = sp.Matrix([1,1,-2])
one = sp.ones(3,1)
I = sp.eye(3)
J = sp.ones(3)/3
T = J + Q(64,729)*v*v.T/2 + Q(1,729)*w*w.T/6
T_native, _, _, _ = symmetric_transfer()
assert np.max(abs(T_native-np.array(T,float))) < 1e-15
assert T*one == one and T==T.T and min(T)>0


def cylinder(P, xs):
    return Q(1,3)*sp.prod(P[i,j] for i,j in zip(xs,xs[1:]))


for P in (T, I, Q(1,1000)*T + Q(999,1000)*I):
    assert P*one == one and P==P.T and min(P)>=0
    for xs in product(range(3),repeat=3):
        p = cylinder(P,xs)
        assert sum(cylinder(P,xs+(j,)) for j in range(3)) == p
        assert sum(cylinder(P,(j,)+xs) for j in range(3)) == p

# The literal historical ungapped negative control is not a stochastic kernel.
bad = J + Q(999,1000)*v*v.T/2
assert bad[0,1] == -Q(997,6000)

# For the uniform faithful stationary reference, Petz on the actual classical
# algebra is the reverse transition, which equals T. This is not T^{-1}.
R = T.T
assert R*T*v == Q(4096,531441)*v and R*T*w == Q(1,531441)*w
assert T.inv()[0,1] > 0 and min(T.inv()) < 0
print('v330 native matrix verified; exact v330 cylinder consistency incl. the gapless endpoint.')
print('v330 literal negative-control edge:', bad[0,1])
print('v221/v330 native stochastic Petz round-trip deviation factors:', Q(4096,531441), Q(1,531441))
print('Min entry of formal inverse:', min(T.inv()))
print('Original v330 susceptibility of diag(1,-1,0):', Q(2,3)*(1+Q(64,729))/(1-Q(64,729)))
