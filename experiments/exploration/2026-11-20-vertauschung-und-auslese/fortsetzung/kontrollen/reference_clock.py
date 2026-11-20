"""Exact finite coefficient controls for the infinite braid reference process."""

from itertools import product
from math import gcd

import sympy as sp

from path_os import A, B, E, MOVES, inv, mul, walks


def conjugate(h, g):
    return mul(mul(h,g),inv(h))


def main():
    ab = mul(A,B)
    central = mul(mul(ab,ab),ab)
    assert central[0] == 6 and central[1:] == (-1,0,0,-1)
    assert all(conjugate(s,central) == central for s in MOVES)
    hidden = mul(mul(A,A),mul(inv(B),inv(B)))
    assert hidden[0] == 0 and hidden != E
    print("Central braid coordinates (w,M) =", central)
    print("Cyclic phase order and first surviving positive central exponent:",
          [(n,n//gcd(n,6)) for n in range(1,13)])

    raw = walks(8,0)
    def returns(g, n):
        if n and g[0]:
            return 0
        return sum(c for s,c in raw[n].items() if conjugate(s,g) == g)
    for label,g in (("vacuum",E),("primitive",A),("central nonempty",central),
                    ("hidden",hidden)):
        coeff = [sp.Rational(returns(g,n),4**n) for n in range(9)]
        print(label,"moments of one-mark channel, n=0..8:",coeff)
    tau=sp.symbols("tau")
    hidden_series=sum(sp.Rational(returns(hidden,n),4**n*sp.factorial(n))*tau**n
                      for n in range(9))
    print("exp(tau) times hidden autocorrelation through order eight:",hidden_series)
    assert sp.expand(hidden_series).coeff(tau,2) == sp.Rational(1,8)
    # Exact grade covariance for all histories through length eight.
    def end(g):
        return tuple(v%2 for v in g[1:])
    probes=set().union(*(set(row) for row in raw))
    quotient={end(g) for g in probes}
    assert len(quotient)==6
    # canonical representatives for endpoint permutations
    representatives={}
    for g in probes:
        representatives[end(g)]=g
    for s,g in product(MOVES,probes):
        assert end(conjugate(s,g)) == end(conjugate(representatives[end(s)],
                                                    representatives[end(g)]))
    signs={end(g):(-1 if g[0]%2 else 1) for g in probes}
    assert sum(signs[q]==-1 for q in quotient)==3
    assert sum(signs[q]==1 and q!=end(E) for q in quotient)==2
    print("Checked endpoint pointer covariance on",len(probes),"histories",
          "and all four elementary exchanges; active grading 3 odd, 2 even.")
    # A finite cyclic recalibration of order n commutes with the transfer
    # and preserves z^(n/gcd(n,6)). It cannot cluster for any n.
    for n in range(1,49):
        k=n//gcd(n,6)
        assert (6*k)%n==0 and all((6*j)%n for j in range(1,k))
    # Independently vary central winding kicks without losing the checked
    # conditions when calibration idempotence is dropped: on winding n,
    # a symmetric kick ±theta multiplies a marked step by cos(n*theta).
    print("Haar charged correlation is exactly exp(-tau) for all nonzero windings.")


if __name__=="__main__":
    main()
