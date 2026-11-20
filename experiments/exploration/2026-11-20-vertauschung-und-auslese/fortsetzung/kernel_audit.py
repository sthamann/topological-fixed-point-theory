"""The flag kernel on signed B3 words and an invisible perturbation of its vacuum response."""
from functools import lru_cache
from itertools import product
import sympy as sp

from flag_process import change, I, tr
from word_coordinates import reduced_quotient, degree, group_inverse, reflection, positives


def key(word):
    return degree(word), reduced_quotient(word)


def kernel_factory(r, epsilon):
    moves = {1: change(0), 2: change(1)}
    pure = (1, -2)*3
    orbit = {key(v) for v in (pure, group_inverse(pure), reflection(pure),
                              group_inverse(reflection(pure)))}
    assert key(()) not in orbit and all(n == 0 for n, _ in orbit)
    @lru_cache(maxsize=None)
    def transfer(word):
        return I if not word else transfer(word[:-1])*moves[word[-1]]

    @lru_cache(maxsize=None)
    def phi(word):
        k=key(word)
        if k in orbit:
            return epsilon
        n=k[0]
        if n<0:
            return phi(group_inverse(word))
        positive=positives(n).get(k[1])
        return sp.S.Zero if positive is None else r**n*tr(transfer(positive))
    return phi,len(orbit)


def run():
    r=sp.Rational(1,4)
    epsilon=sp.Rational(1,16)
    phi,orbit=kernel_factory(r,0)
    changed,_=kernel_factory(r,epsilon)
    pure=(1,-2)*3
    assert (phi(pure),changed(pure))==(0,epsilon)
    # Off-identity bound: all traces are <=1/2, and the length-one trace is 0.
    # Up to two signs of each positive word, at most 2^n spellings.
    lower=1-(2*r)**2/(1-2*r)
    assert lower==sp.Rational(1,2) and lower-orbit*epsilon>0
    words=[()]+[(i,) for i in (1,2,-1,-2)]+[(1,-2),(2,-1),pure]
    assert len(set(map(key,words)))==len(words)
    for f in (phi,changed):
        gram=sp.Matrix([[f(group_inverse(u)+v) for v in words] for u in words])
        assert gram.is_positive_definite
    short=[w for n in range(4) for w in product((1,2),repeat=n)]
    for u,v in product(short,repeat=2):
        history=group_inverse(reflection(u))+v
        assert phi(history)==changed(history)
    print("Ordinary full-group Gram lower bounds:",lower,lower-orbit*epsilon)
    print("Two kernels agree on the entire positive half; hidden vacuum response:",
          phi(pure),changed(pure))


if __name__=="__main__":
    run()
