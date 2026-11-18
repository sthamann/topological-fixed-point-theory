#!/usr/bin/env python3
"""Exact common discriminant/clock calculation; no physical-source assertion.

Input: the already reconstructed unitary lattice VOA
U=V_D5 tensor V_A3 and its marked diagonal C4 extension V_E8.
This checker uses integer arithmetic, including Z[zeta_24].
"""
from __future__ import annotations

import itertools as it
import json
from fractions import Fraction as F
from math import gcd, lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


A = tuple(it.product(range(4), repeat=2))
ZERO = (0, 0)


def add(x, y):
    return tuple((a+b) % 4 for a, b in zip(x, y))


def q8(x):
    a, b = x
    return (5*a*a+3*b*b) % 8


def bilinear4(x, y):
    return (5*x[0]*y[0]+3*x[1]*y[1]) % 4


def flux(x):
    return (x[0]-x[1]) % 4


def polynomial_reduce(coeff):
    """Reduction modulo Phi_24(X)=X^8-X^4+1."""
    v = list(coeff)+[0]*max(0, 8-len(coeff))
    for k in range(len(v)-1, 7, -1):
        c = v[k]
        v[k-4] += c
        v[k-8] -= c
    return tuple(v[:8])


RZERO = (0,)*8
RONE = (1,)+(0,)*7
POW = tuple(polynomial_reduce([0]*k+[1]) for k in range(24))


def radd(a, b):
    return tuple(x+y for x, y in zip(a, b))


def rmul(a, b):
    p = [0]*15
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    p[i+j] += x*y
    return polynomial_reduce(p)


def matmul(a, b):
    out = [[RZERO for _ in b[0]] for _ in a]
    for i, row in enumerate(a):
        for k, x in enumerate(row):
            if x != RZERO:
                for j, y in enumerate(b[k]):
                    if y != RZERO:
                        out[i][j] = radd(out[i][j], rmul(x, y))
    return out


def check():
    require(all((q8(add(x, y))-q8(x)-q8(y)) % 8 ==
                2*bilinear4(x, y) for x in A for y in A),
            "quadratic polarization")
    require([x for x in A if all(bilinear4(x,y)==0 for y in A)] == [ZERO],
            "discriminant nondegeneracy")

    isotropic = tuple(x for x in A if q8(x)==0)
    supports = []
    for bits in it.product((0,1), repeat=len(isotropic)-1):
        h = frozenset([ZERO]+[x for x,b in zip(isotropic[1:],bits) if b])
        if all(add(x,y) in h for x in h for y in h):
            supports.append(h)
    supports.sort(key=lambda h:(len(h), sorted(h)))
    hp = frozenset((k,k) for k in range(4))
    hm = frozenset((k,(-k)%4) for k in range(4))
    require(len(supports)==4 and {h for h in supports if len(h)==4} == {hp,hm},
            "complete local support census")
    require(hp & hm == frozenset((ZERO,(2,2))), "two closures share D8")
    require(not any(len(h)==4 and h & hp == {ZERO} for h in supports),
            "no transverse local holomorphic extension")

    # All sixteen U modules, arranged by the monodromy of the marked glue h.
    hd = (F(0),F(5,8),F(1,2),F(5,8))
    ha = (F(0),F(3,8),F(1,2),F(3,8))
    dd, da = (1,16,10,16), (1,4,6,4)
    fibers = []
    for f in range(4):
        labels = [x for x in A if flux(x)==f]
        weights = {x:hd[x[0]]+ha[x[1]] for x in labels}
        rho = min(weights.values())
        tops = [x for x in labels if weights[x]==rho]
        fibers.append({"monodromy_flux":f,"Li_kappa_power":(-f)%4,"rho":str(rho),
                       "top_dimension":sum(dd[a]*da[b] for a,b in tops),
                       "top_U_labels":[list(x) for x in tops],
                       "branches":[{"label":list(x),"weight":str(weights[x]),
                                    "top_dimension":dd[x[0]]*da[x[1]]} for x in labels]})
    require([f["rho"] for f in fibers]==["0","3/8","1/2","3/8"],
            "four twisted ground weights")
    require([f["top_dimension"] for f in fibers]==[1,4,16,4],
            "four twisted ground dimensions")

    # Explicit isometry to the type 4{2} cyclic orbifold quadratic module:
    # (i,j) -> (3i+j, 2i+j), inverse i=a-b, j=3b-2a.
    iso = {(i,j):((3*i+j)%4,(2*i+j)%4) for i,j in A}
    require(len(set(iso.values()))==16, "orbifold discriminant isomorphism")
    require(all(q8(iso[i,j]) == (i*i+2*i*j)%8 for i,j in A),
            "type 4{2} quadratic form, not just group order")
    require(all(flux(iso[i,j])==i for i,j in A), "twist projection respected")

    # Modular S numerator is a fourth root of unity; denominator is 4.
    sn = [[POW[(-6*bilinear4(x,y))%24] for y in A] for x in A]
    te = [(3*q8(x)-8)%24 for x in A]
    stn = [[rmul(sn[i][j],POW[te[j]]) for j in range(16)] for i in range(16)]
    s2 = matmul(sn,sn)
    st3 = matmul(matmul(stn,stn),stn)
    for i,x in enumerate(A):
        for j,y in enumerate(A):
            target = RONE if add(x,y)==ZERO else RZERO
            require(s2[i][j] == tuple(16*c for c in target), "S^2=C")
            require(st3[i][j] == tuple(64*c for c in target), "(ST)^3=C")
    order = lcm(*(24//gcd(24,e) for e in te))
    require(order==24, "faithful T has exact order24")
    require(lcm(*(8//gcd(8,q8(x)) for x in A))==8, "quadratic level8")
    anomaly = (16*F(3,8))%4
    require(anomaly==2, "cyclic orbifold type2")
    return {
        "status":"PASS_EXACT_COMMON_DISCRIMINANT_AND_CLOCK",
        "input":"marked U=V_D5 tensor V_A3; diagonal E8 extension and its dual clock",
        "quadratic_form":"q(a,b)=(5a^2+3b^2)/8 mod1 on C4^2",
        "isotropic_supports":[[list(x) for x in sorted(h)] for h in supports],
        "twist_projection":"a-b mod4",
        "orientation":"kappa=exp(+2pi*i*v_0); E8-a*v has monodromy_flux=-a mod4",
        "twisted_fibers":fibers,
        "orbifold_isometry":"(i,j)->(3i+j,2i+j) mod4; q=i^2/8+ij/4",
        "orbifold_type":"4{2}",
        "quadratic_level":8,
        "central_charge":8,
        "modular_T_exponents_over24":te,
        "modular_T_order":order,
        "modular_relations":"S^2=(ST)^3=charge conjugation, exact Z[zeta24]",
        "transverse_isotropic_complement":False,
        "coefficient_arithmetic":"integer, rational, cyclotomic integer",
        "not_claimed":["physical spacetime or microscopic source reconstruction",
                       "a 24-sheet base cover or a new physical time",
                       "a preferred E6 vertex","native4 equals twisted ground4"]
    }


if __name__ == "__main__":
    print(json.dumps(check(),sort_keys=True,indent=2))
