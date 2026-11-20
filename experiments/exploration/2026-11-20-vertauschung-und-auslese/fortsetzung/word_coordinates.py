"""Exact positive-type OS kernel on the full braid group B3.

Signed words use +1=a, +2=b and negative integers for group inverses.
The infinite positivity proofs are analytic; the finite matrices are cross-checks.
"""
from functools import lru_cache
from itertools import product

from sympy import Matrix, Rational, S, diag, eye, simplify, sqrt

A = diag(1, -Rational(1, 2))
B = Matrix([[-1, sqrt(14)], [sqrt(14), 4]])/6
T = {1:A, 2:B}
r = Rational(1, 8)

# With X=aba and Y=ab, B3=<X,Y | X^2=Y^3>. This common power is
# central and has exponent six. Modding it out gives C2 * C3.
# Exponent and the reduced C2*C3 word therefore separate elements of B3.
letters = {
    1: (("Y", 2), ("X", 1)),  -1: (("X", 1), ("Y", 1)),
    2: (("X", 1), ("Y", 2)),  -2: (("Y", 1), ("X", 1)),
}


def reduced_quotient(word):
    stack = []
    for letter in word:
        for tag, power in letters[letter]:
            if stack and stack[-1][0] == tag:
                power = (stack.pop()[1] + power) % (2 if tag == "X" else 3)
            if power:
                stack.append((tag, power))
    return tuple(stack)


def degree(word):
    return sum(1 if letter > 0 else -1 for letter in word)


def group_inverse(word):
    return tuple(-g for g in reversed(word))


def reflection(word):
    return tuple(-g for g in word)


@lru_cache(maxsize=None)
def positives(n):
    # Finite, exact P-membership at the prescribed exponent. Useful for modest words.
    return {reduced_quotient(w):w for w in product((1,2), repeat=n)}


@lru_cache(maxsize=None)
def representation(positive_word):
    if not positive_word:
        return eye(2)
    return representation(positive_word[:-1])*T[positive_word[-1]]


def trace2(mat):
    return simplify(mat.trace()/2)


@lru_cache(maxsize=None)
def f(word):
    n = degree(word)
    if n < 0:
        return f(group_inverse(word))
    quotient = reduced_quotient(word)
    p = positives(n).get(quotient)
    if p is None:
        return S.Zero
    return r**n*trace2(representation(p))


def main():
    assert A == A.T and B == B.T
    assert A*B*A == B*A*B and A**2 == (A+eye(2))/2 and B**2 == (B+eye(2))/2
    assert A.det() == B.det() == -Rational(1,2) and (A*B)**3 == eye(2)/8
    for side in ((1,2,1,-2,-1,-2), (-1,1,-2,2), (1,2)*6):
        assert not reduced_quotient(side)
    assert f((1,-2)) == 0 and f((1,2,1,-2,-1,-2)) == 1
    assert f((1,2)*6) == r**12 / 64   # Quotient identity, but not group identity.

    # Infinite group Gram bound: ||f-delta_e||_1 <= 2 sum_(n>=1) (2r)^n.
    off_bound = 4*r/(1-2*r)
    assert off_bound == Rational(2,3)

    signed_words = ((),(1,),(2,),(-1,),(-2,),(1,1),(2,2),(1,2),(2,1),(1,-2),(2,-1))
    gram = Matrix([[f(group_inverse(u)+v) for v in signed_words] for u in signed_words])
    assert gram == gram.T and all(x==1 for x in gram.diagonal())
    assert all(sum(abs(gram[i,j]) for j in range(len(signed_words)) if i!=j) <= off_bound
               for i in range(len(signed_words)))

    observed = ((),(1,),(2,),(1,1),(1,2),(2,1),(2,2))
    reflected = Matrix([[f(group_inverse(reflection(u))+v) for v in observed] for u in observed])
    hs = Matrix([[r**(len(u)+len(v))*trace2(representation(u).T*representation(v))
                  for v in observed] for u in observed])
    assert hs == reflected and reflected.rank() == 4
    norm = simplify(r**4 * trace2((A**2 - B**2)**2))
    assert norm == Rational(7,65536)
    assert f((1,1)) == f((2,2)) == Rational(5,512)
    # A common continuation distinguishes the histories in the unbiased trace.
    probe = f((1,1,1))-f((2,2,1))
    assert probe == Rational(7,8192) and f((2,2,1))==0
    assert all(f(reflection(w)) == f(w) for w in signed_words)
    print("full Gram lower bound:",1-off_bound)
    print("rank of reflected Gram for words of length <=2:",reflected.rank())
    print("f(a^2) = f(b^2) =",f((1,1)))
    print("OS squared distance of a^2 and b^2:",norm)
    print("f(a^3), f(b^2 a), difference:",f((1,1,1)),f((2,2,1)),probe)
    print("Hecke central transfer (AB)^3:",(A*B)**3)
    print("reflection invariance, braid, full-group mixed-word samples and Gram bounds: OK")


if __name__=="__main__":
    main()
