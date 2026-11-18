"""Close the actual compensated clock word in the neutral original CAR currents."""
from collections import defaultdict
from functools import cache
from itertools import combinations
from pathlib import Path
import runpy

import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
src = runpy.run_path(str(HERE / "check_ramond_clock_return_20261110.py"))
original = src["original"]
Q, V, R = src["Q"], src["V"], src["R"]
q, gamma, q_prime = src["q"], src["gamma"], src["q_prime"]
conj = original["m"].conj
zero_power = (0, 0, 0)


def shift(poly, terms):
    out = defaultdict(lambda: Q.zero)
    for power, a in poly.items():
        for inc, b in terms:
            out[tuple(x + y for x, y in zip(power, inc, strict=True))] += a * b
    return {p: a for p, a in out.items() if a}


def kernel(left_axes, right_axes, charges, exponent):
    """Neutral ordered kernel matrix element on the original zero-charge core."""
    left = (left_axes.H * sp.Matrix.hstack(*charges)).applyfunc(sp.simplify)
    right = (right_axes.T * sp.Matrix.hstack(*charges)).applyfunc(sp.simplify)
    metric = (left_axes.H * right_axes).applyfunc(sp.simplify)
    nl, nr = left.rows, right.rows
    ls = [[Q.from_sympy(left[i, j]) for j in range(len(charges))] for i in range(nl)]
    rs = [[-Q.from_sympy(right[i, j]) for j in range(len(charges))] for i in range(nr)]
    inner = [[Q.from_sympy(metric[i, j]) for j in range(nr)] for i in range(nl)]
    origin = zero_power[:len(charges)]

    def powers(n, i, values):
        return [(tuple(n * int(k == j) for k in range(len(charges))), a)
                for j, a in enumerate(values[i]) if a]

    @cache
    def wick(a, b):
        if not a:
            out = {origin: Q.one}
            for n, j in b:
                out = shift(out, powers(-n, j, rs))
            return out
        (n, i), rest = a[0], a[1:]
        out = defaultdict(lambda: Q.zero, shift(wick(rest, b), powers(n, i, ls)))
        for k, (m, j) in enumerate(b):
            if m == n and inner[i][j]:
                for power, value in wick(rest, b[:k] + b[k+1:]).items():
                    out[power] += n * inner[i][j] * value
        return {p: a for p, a in out.items() if a}

    @cache
    def read_word(a, b):
        return sum((exponent(p) * x for p, x in wick(a, b).items() if exponent(p)),
                   Q.zero)

    def read(a, b):
        return Q.to_sympy(sum((conj(x) * y * read_word(u, v)
                               for u, x in a.items() for v, y in b.items()), Q.zero))
    return read


def ground_packet_cross():
    charges = (q_prime, gamma, -q)  # Increasing radial coordinate order.
    assert sum(charges, sp.zeros(8, 1)) == sp.zeros(8, 1)
    assert [charges[i].dot(charges[j]) for i, j in ((0,1),(0,2),(1,2))] == [-2,0,-2]

    # Expand (z3-z2)^-2 (z2-z1)^-2. At equal outside grade the
    # required exponent (-1,-2,-1) picks unique terms in the two series.
    def triple_coefficient(p):
        p1, _, p3 = p
        return (-p1) * p3 if p1 <= -1 and p3 >= 1 else 0

    # The empty return word is (V_q)_0^* (V_q)_0.
    bare = kernel(V, V, (q, -q), lambda p: p[1] if p[1] >= 1 else 0)
    marked = kernel(V, R*V, charges, triple_coefficient)
    polynomials, _ = original["m"].code()
    phase = src["cocycle_phase"]()
    old = {
        sp.Rational(1,4): sp.Matrix([[-sp.Rational(77,16)+16*sp.I,9*sp.I/64],
                                     [-459*sp.I/512,sp.Rational(1539,65536)]]),
        sp.Rational(3,4): sp.Matrix([[-sp.Rational(77,16)+16*sp.I,sp.I/64],
                                     [-51*sp.I/512,sp.Rational(19,65536)]]),
    }
    for alpha in src["alphas"]:
        drift = src["drift"](alpha)
        d = [Q.from_sympy(x) for x in (V.T * drift).applyfunc(sp.simplify)]
        shifted = [original["series_response"](f, d, 4) for f in polynomials]
        ns = [original["input_word_map"](p) for p in shifted]
        G = sp.Matrix(5,5,lambda i,j: bare(ns[i],ns[j]))
        C = phase * sp.Matrix(5,5,lambda i,j: marked(ns[i],ns[j]))
        w = [original["vertex_zero"](p) for p in shifted]
        assert G == sp.Matrix(5,5,lambda i,j: original["inner"](w[i],w[j]))
        assert C[:2,:2] == old[alpha] and C.rank() == 5
        assert G.rank() == 5 and G == G.H
        # On the two intermediate fibers every 2pi charge character is -I.
        # On the neutral external sector and on all integer CAR sectors it is +I.
        for momentum in (q_prime, q):
            assert all(x % 1 == sp.Rational(1,2) for x in momentum)
        assert (C[:2,:2]+C[:2,:2].H)[0,1] != 0  # The two phase readings are distinct.
        print("quarter",alpha,"closed neutral Gram/cross rank",G.rank(),C.rank(),flush=True)
        print("original packet cross, two-state block",C[:2,:2].tolist(),flush=True)


def cut_support():
    """Actual two-packet source spans for outer and diagnostic inner forms."""
    a = sp.Matrix.hstack(V, q)
    axes = sp.Matrix.hstack(a, R*a, gamma)
    m = (axes.H*axes).applyfunc(sp.simplify)
    metric = [[Q.from_sympy(x) for x in row] for row in m.tolist()]
    dot = [Q.from_sympy(x) for x in (axes.T*gamma).applyfunc(sp.simplify)]
    store = original["inner"].__globals__["inner_q"]
    saved = list(store)
    store[:] = [metric, metric]
    original["contractions"].cache_clear()

    def shift_axes(state):
        return {tuple((n,j+5) for n,j in word):a for word,a in state.items()}

    def zero_mode(state, sign):
        out = defaultdict(lambda: Q.zero)
        for word, value in state.items():
            for size in range(len(word)+1):
                for omitted in combinations(range(len(word)), size):
                    amp = value
                    for j in omitted:
                        amp *= -sign*dot[word[j][1]]
                    if not amp:
                        continue
                    kept = tuple(x for j,x in enumerate(word) if j not in omitted)
                    degree = sum(word[j][0] for j in omitted)
                    for extra, coeff in original["exp_y"](degree):
                        factors = tuple(sorted(kept + tuple((n,10) for n,_ in extra)))
                        out[factors] += amp * coeff * sign**len(extra)
        return {k:v for k,v in out.items() if v}

    inner = original["inner"]
    phase = src["cocycle_phase"]()
    polynomials, _ = original["m"].code()
    try:
        for alpha in src["alphas"]:
            d = [Q.from_sympy(x) for x in (V.T * src["drift"](alpha)).applyfunc(sp.simplify)]
            states = [original["vertex_zero"](original["series_response"](f,d,4))
                      for f in polynomials[:2]]
            moved = [shift_axes(w) for w in states]
            dressed = [{k:Q.from_sympy(phase)*v for k,v in zero_mode(w,1).items()} for w in moved]
            returning = [{k:Q.from_sympy(phase)*v for k,v in zero_mode(w,-1).items()} for w in states]

            def gram(v):
                return sp.Matrix(4,4,lambda i,j:inner(v[i],v[j]))

            outside = gram(states + dressed)
            midpoint = gram(returning + moved)
            assert outside[:2,2:] == midpoint[:2,2:]  # True adjoint on the same cut.
            assert outside.rank() == midpoint.rank() == 4
            for form in (outside, midpoint):
                assert all(form[:k,:k].det() > 0 for k in range(1,5))
            # Removing just the common empty scalar cannot determine the physical cross form.
            print("quarter",alpha,"outer/inner four-vector span",outside.rank(),midpoint.rank(),flush=True)
    finally:
        store[:] = saved
        original["contractions"].cache_clear()


def thermal_center_read():
    """The old cylinder thermal comparison with its exact original center reference."""
    mp.mp.dps = 65
    beta = 4*mp.pi
    T = 8*beta
    x1, x2, x3 = 0, mp.mpf(2)/5, mp.mpf(7)/5  # Old H gaps 1/20, 1/8.
    nums = lambda v: [mp.mpf(int(a.p))/int(a.q) for a in v]
    X = [a*x1+b*x2-c*x3 for a,b,c in zip(nums(q_prime),nums(gamma),nums(q),strict=True)]
    center = list(map(int,src["center"]))

    def center_theta(delta, X, cutoff):
        out = [mp.mpf(1),mp.mpf(0),mp.mpf(0)]
        for w,d,x in zip(center,delta,X,strict=True):
            f = [mp.mpf(0)]*3
            for p in range(-cutoff,cutoff+1):
                f[(w*p)%3] += mp.exp(-T*(mp.mpf(p*p)/2+d*p)+p*x)
            out = [sum(out[j]*f[(k-j)%3] for j in range(3)) for k in range(3)]
        return out[0]

    def prime(t, cutoff):
        return 2*mp.sinh(t/2)*mp.fprod(
            ((1-mp.exp(-k*T+t))*(1-mp.exp(-k*T-t)))/(1-mp.exp(-k*T))**2
            for k in range(1,cutoff+1))

    def values(alpha, cutoff):
        d = nums(src["drift"](alpha))
        empty = center_theta(d,[0]*8,cutoff)
        posc = mp.fprod((1-mp.exp(-k*T))**(-8) for k in range(1,cutoff+1))
        charged = center_theta(d,X,cutoff)
        z0 = posc * empty
        cross = (int(src["cocycle_phase"]())*posc*charged*mp.exp(sum(a*b for a,b in zip(d,X))) /
                 (prime(x3-x2,cutoff)*prime(x2-x1,cutoff))**2)
        return z0, cross, cross/z0

    for alpha in src["alphas"]:
        a, b = values(alpha,4), values(alpha,7)
        assert max(abs(x-y) for x,y in zip(a,b,strict=True)) < mp.mpf("1e-58")
        print("quarter",alpha,"cylinder center empty, marked, ratio",*[mp.nstr(v,24) for v in b],flush=True)


if __name__ == "__main__":
    ground_packet_cross()
    cut_support()
    thermal_center_read()
