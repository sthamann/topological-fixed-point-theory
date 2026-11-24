"""Time-reflected braid histories. Run with Python 3 and SymPy."""
from collections import Counter
from itertools import product
from math import factorial

import sympy as sp


# The Artin representation paired with abelianization is faithful on B3.
# Coordinates: (w,a,b,c,d), with integer determinant-one matrix [[a,b],[c,d]].
E = (0, 1, 0, 0, 1)
A = (1, 1, 1, 0, 1)
B = (1, 1, 0, -1, 1)


def inv(g):
    w, a, b, c, d = g
    return -w, d, -b, -c, a


def mul(g, h):
    w, a, b, c, d = g
    v, p, q, r, s = h
    return w+v, a*p+b*r, a*q+b*s, c*p+d*r, c*q+d*s


MOVES = (A, B, inv(A), inv(B))


def walks(n, hold):
    # hold=4 yields one half identity, 4 generators of weight 1/8.
    # hold=0 is the embedded discrete chain of the Poisson process.
    layers = [Counter({E: 1})]
    for k in range(n):
        counts = Counter()
        for g, count in layers[k].items():
            if hold:
                counts[g] += count*hold
            for h in MOVES:
                counts[mul(g, h)] += count
        layers.append(counts)
    return layers


def brownian_clock():
    # Complex Gaussian moments of independent coordinates in units l=1.
    # The Vandermonde is translation-invariant, so the Gaussian center of
    # mass factors and does not change the norm ratios below.
    z = sp.symbols("z0:3")
    v = sp.prod(z[i]-z[j] for i in range(3) for j in range(i+1,3))
    def norm2(poly):
        return sum(abs(c)**2*sp.prod(factorial(k) for k in exponent)
                   for exponent,c in sp.Poly(poly,*z).terms())
    for m in (1,3,5):
        amplitude = v**m
        mass = norm2(amplitude)
        gradient = sum(norm2(sp.diff(amplitude,zi)) for zi in z)
        assert gradient == 3*m*mass
        # R_rel^2 is R^2 minus the center of mass Gaussian of mean one.
        radius = sum(norm2(zi*amplitude) for zi in z)/mass-1
        assert radius == 2+3*m
        clock = gradient/(m*m*mass)
        print("Repulsion m, pi^2 times winding rate(l=1), relative radius, product:",
              m, clock, radius, sp.factor(clock*radius))
    assert norm2(v) == 12 and sum(norm2(sp.diff(v,zi)) for zi in z) == 36
    # For free Brownian motion z1-z2~r with z3 separated, times a fixed s>0:
    # the density is bounded below and the local bracket diverges in expectation.
    r, eps = sp.symbols("r eps", positive=True)
    assert sp.limit(sp.integrate(1/r,(r,eps,1)),eps,0,dir="+") == sp.oo


def main():
    g = mul(mul(A, A), mul(inv(B), inv(B)))
    assert g == (0, 5, 2, 2, 1) and g != E
    raw = walks(8, 0)
    counts = [(k,raw[k][g],raw[k][E]) for k in range(9)]
    print("Nonlazy n, paths to g, paths to identity:", counts)
    assert next(k for k in range(9) if raw[k][g]) == 4
    print("At operational time tau, first hidden-loop coefficient in exp(tau)*p_tau(g):",
          sp.Rational(raw[4][g], 4**4*sp.factorial(4)))

    probes = [(1,E),(2,g),(3,inv(g)),(3,A)]
    lazy = walks(6, 4)
    def kernel(t,h):
        return sp.Rational(lazy[t][h], 8**t)
    gram = sp.Matrix([[kernel(t+s,mul(inv(h),j)) for s,j in probes]
                      for t,h in probes])
    # Calculate separately from the half-path amplitudes, without the
    # concatenation formula, including all braid lifts in both halves.
    halves = [{mul(h,x):sp.Rational(c,8**t) for x,c in lazy[t].items()}
              for t,h in probes]
    actual = sp.Matrix([[sum(value*v.get(x,0) for x,value in u.items())
                         for v in halves] for u in halves])
    assert gram == actual and gram == gram.T
    minors = [gram[:i,:i].det() for i in range(1,len(probes)+1)]
    assert all(v>0 for v in minors)
    print("Exact time-reflected principal minors:", minors)
    times = (2,4,6)
    print("Lazy discrete heat returns and hidden endpoints:",
          [(t,kernel(t,E),kernel(t,g)) for t in times])

    # Isotropic rates and inversion guarantee the winding projection
    # increments symmetrically at each event. This marker ignores generator names.
    h, tau = sp.symbols("h tau", real=True)
    mgf = sp.exp(tau*(sp.cosh(h)-1))
    even = [sp.expand(sp.diff(mgf,h,k).subs(h,0)) for k in (2,4,6)]
    print("Unconditioned winding moments orders 2,4,6:", even)
    assert even == [tau, 3*tau**2+tau, 15*tau**3+15*tau**2+tau]
    # Witness against survival under coarsening: identity and g have
    # precisely the same winding and same action modulo 2.
    assert g[0] == E[0] and tuple(v%2 for v in g[1:]) == tuple(v%2 for v in E[1:])
    # Compare every equally weighted degree-four word for the nonlazy readout.
    no_return, distinct = 0, set()
    for letters in product(MOVES,repeat=4):
        value = E
        for a in letters:
            value = mul(value,a)
        if value[0] == 0 and tuple(v%2 for v in value[1:]) == (1,0,0,1):
            distinct.add(value)
            no_return += value != E
    print("Length-four strings with trivial coarse observables but nonidentity:",
          no_return,"/",4**4,"; distinct coarse-null lift endpoints:",len(distinct))
    # Gauss--Bonnet numbers here multiply pi. Cone pullback has four defects of -2pi.
    curv_base = 2*2
    area_pullback = 2*curv_base
    total_intrinsic = area_pullback-4*2
    assert area_pullback == 8 and total_intrinsic == 0
    print("Pulled back curvature vs intrinsic conical curvature, in pi units:",
          area_pullback,total_intrinsic)
    mirror_g=mul(mul(inv(A),inv(A)),mul(B,B))
    invisible_orbit={g,inv(g),mirror_g,inv(mirror_g)}
    assert len(invisible_orbit) == 4
    assert all(v[0]==0 and tuple(k%2 for k in v[1:])==(1,0,0,1)
               for v in invisible_orbit)
    print("Added symmetric hidden-jump menu: new derivative at g = eta/",
          len(invisible_orbit),"; projected derivative identically zero.")
    brownian_clock()


if __name__ == "__main__":
    main()
