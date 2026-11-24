#!/usr/bin/env python3
"""Exact conditional SU(5)_1 x SU(5)_1 conformal extension into E8_1.

Physics assumptions NOT derived: two chiral affine SU(5) level-one
current algebras, a local Z5 simple-current extension and its source.
Check h_k=k(5-k)/10, c=24/6=4, integral spin for k -> m k
mod 5, 248=24+24+4*50. This is a known conformal embedding.
"""
from fractions import Fraction
from math import comb
c = Fraction(24, 6)
assert c == 4 and 2*c == 8
h = lambda k: Fraction(k*(5-k), 10)
good = []
for m in (1,2,3,4):
    weights = []
    degrees = []
    for k in (1,2,3,4):
        j = (m*k) % 5
        weights.append(h(k)+h(j))
        degrees.append(comb(5,k)*comb(5,j))
    local = all(a.denominator == 1 for a in weights)
    print("m=",m,"h=",weights,"dims=",degrees,"bosonic=",local)
    if local:
        good.append(m)
        assert weights == [Fraction(1)] * 4
        assert degrees == [50] * 4
        assert 24+24+sum(degrees)==248
assert good == [2,3]
print("EXACT SU5_1^2 c=8, level-one current census 248, 2/4 even Z5 glues.")
print("Known conformal E8_1 extension, not a physical derivation from anyons.")
