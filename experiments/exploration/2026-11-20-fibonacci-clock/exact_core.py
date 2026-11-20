#!/usr/bin/env python3
"""Fibonacci four-anyon braid clock over Q(zeta_10), no numeric tolerance.

This is an algebraic candidate for parts of TFPT, not a derivation of P1,
physical CAR fermions, 3+1D or gravity. Input: nonpointed rank-two
unitary Fibonacci fusion category and four tau excitations.

Matrices are obtained from Fibonacci F/R moves in the basis of four
tau-anyon fusion channels (two with vacuum charge, three with tau charge).
"""
from math import comb
from fractions import Fraction
from sympy import Matrix, Poly, Symbol, eye, zeros, expand, Rational

z = Symbol("z")
x = Symbol("x")
cycl = Poly(z**4 - z**3 + z**2 - z + 1, z, domain="QQ")

def red(y):
    return Poly(expand(y), z, domain="QQ").rem(cycl).as_expr()

def mat(A):
    return A.applyfunc(red)

def mult(*matrices):
    result = matrices[0]
    for another in matrices[1:]:
        result = mat(result * another)
    return result

def eq(A, B):
    return mat(A - B) == zeros(A.rows, A.cols)

def four_anyon_braids():
    # z=exp(i*pi/5), R1=-z, Rtau=z**3,
    # a = 2*cos(2*pi/5) = z**2-z**3, F^2=I.
    t, s = -z, z**3
    u, v, w, r = z**2-z, -z**3+z**2-z, -z**2, z**3-z**2
    return [
       [Matrix([[t, 0], [0, s]]),
        Matrix([[u, v], [w, r]]),
        Matrix([[t, 0], [0, s]])],
       [Matrix.diag(t, s, s),
        Matrix([[u, 0, v], [0, s, 0], [w, 0, r]]),
        Matrix([[s, 0, 0], [0, u, v], [0, w, r]])]
    ]

def main():
    braids = four_anyon_braids()
    cycles = []
    for dim, group in zip((2, 3), braids):
        assert eq(mult(group[0], group[1], group[0]),
                  mult(group[1], group[0], group[1]))
        assert eq(mult(group[1], group[2], group[1]),
                  mult(group[2], group[1], group[2]))
        assert eq(mult(group[0], group[2]), mult(group[2], group[0]))
        cycle = mult(*group)
        power = eye(dim)
        for _ in range(5):
            power = mult(power, cycle)
        J = -power
        assert eq(mult(J, J, J, J), eye(dim))
        if dim == 2:
            assert eq(mult(J, J), eye(2))
            assert red(J.trace()) == 0 and red(J.det()) == -1
        if dim == 3:
            assert eq(mult(J, J, J) + mult(J, J) + J + eye(3), zeros(3))
            assert red(J.trace()) == -1 and red(J.det()) == -1
        cycles.append((cycle, J))
        print("sector dimension:", dim, "; braid relations and J^4=I exact")

    assert red(cycles[0][1].det() * cycles[1][1].det()) == 1
    print("J = -(sigma1 sigma2 sigma3)^5 is SU(5),",
          "spectrum {+1,-1} (+) {-1,+i,-i} in the 2+3 sectors.")

    passes = []
    for k in range(1, 20):
        powers = []
        for sector in (0, 1):
            cycle = cycles[sector][0]
            p = eye(cycle.rows)
            for _ in range(k):
                p = mult(p, cycle)
            powers.append(p)
        for sign in (-1, 1):
            A, B = [sign * p for p in powers]
            if red(A.det() * B.det()) != 1:
                continue
            if not (eq(mult(A, A, A, A), eye(2)) and
                    eq(mult(B, B, B, B), eye(3))):
                continue
            if eq(mult(A, A), eye(2)) and eq(mult(B, B), eye(3)):
                continue
            if not eq(mult(B, B, B) + mult(B, B) + B + eye(3), zeros(3)):
                continue
            passes.append((k, sign))
    assert passes == [(5, -1), (15, -1)]
    print("Exact post-hoc baseline: 2/38 sign/power candidates", passes)

    def four_factor_multiplicities(n):
        a, b = 0, 1
        for _ in range(3):
            a, b = b, a + n*b
        return a, b
    assert four_factor_multiplicities(0) == (1, 0)
    assert four_factor_multiplicities(1) == (2, 3)
    print("rank-two fusion rules: pointed=1+0, Fibonacci=2+3.")

    def j(lam):
        return 256*(1-lam+lam**2)**3/(lam**2*(1-lam)**2)
    assert j(Rational(1, 2)) == 1728
    assert j(Rational(1, 3)) == Rational(21952, 9)
    print("P1 is not selected: same category, j=1728 or j=21952/9.")

    # Explicit EXTRA axiom: reinterpret all five fusion channels as odd
    # CAR one-particle modes with relative hypercharge normalization.
    weights = []
    for p in range(4):
        for q in range(3):
            if (p + q) % 2 == 0:
                weights += [Fraction(-p, 3) + Fraction(q, 2)] * (
                    comb(3, p)*comb(2, q))
    assert len(weights) == 16
    assert sum(weights) == 0
    assert sum(w**3 for w in weights) == 0
    print("conditional odd CAR even exterior: 16, TrY=TrY^3=0.")
    print("PASS EXACT ALGEBRA; PHYSICAL P1, TIME, 3+1, GRAVITY OPEN")

if __name__ == "__main__":
    main()
