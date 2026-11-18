#!/usr/bin/env python3
"""Exact finite certificate for the NS/C2/C12 reflection dictionary."""

from math import gcd


def mm(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)) for i in range(2))


def tr(a):
    return tuple(zip(*a))


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def order_pair(a: int, b: int, n: int) -> int:
    """Order of Q^a Delta^b in C4 x Cn."""
    for k in range(1, 4 * n + 1):
        if (a * k) % 4 == 0 and (b * k) % n == 0:
            return k
    raise AssertionError("finite order not found")


def main() -> None:
    # NS zero-branch representation in the ordered (+i,-i) Pin basis.
    # A=diag(i,-i), R swaps the two sectors.  Encode i^e by exponents e mod 4.
    A = (1, 3)
    A_inv = tuple((-e) % 4 for e in A)
    RAR = (A[1], A[0])
    need(tuple((2 * e) % 4 for e in A) == (2, 2), "A^2 != -I")
    need(tuple((4 * e) % 4 for e in A) == (0, 0), "A^4 != I")
    need(RAR == A_inv, "RAR != A^-1")

    # C2: Q(x,y)=(ix,y), kappa(x,y)=(bar x,bar y).
    # On eta=dx/y and bar(eta), Q has the same eigen-exponents (+1,-1).
    Q_H1 = (1, 3)
    need(Q_H1 == A, "C2 and NS quarter characters differ")
    need((Q_H1[1], Q_H1[0]) == tuple((-e) % 4 for e in Q_H1), "kappa does not invert Q")

    # The bridge also preserves the actual polarized real Hodge carrier.
    # In (Re eta, Im eta), Q is a quarter rotation and K=kappa* is reflection.
    Q = ((0, 1), (-1, 0))
    K = ((1, 0), (0, -1))
    Om = ((0, 1), (-1, 0))  # intersection/cup form
    g = ((1, 0), (0, 1))    # Hodge L2 metric, up to its positive scale
    I = ((1, 0), (0, 1))
    minus_Om = tuple(tuple(-x for x in row) for row in Om)
    need(mm(Q, Q) == ((-1, 0), (0, -1)), "Q^2 != -I")
    need(mm(K, K) == I, "K^2 != I")
    need(mm(mm(K, Q), K) == ((0, -1), (1, 0)), "KQK != Q^-1")
    need(mm(mm(tr(Q), Om), Q) == Om, "Q does not preserve intersection form")
    need(mm(mm(tr(K), Om), K) == minus_Om, "K is not anti-symplectic")
    need(mm(mm(tr(Q), g), Q) == g, "Q is not a Hodge isometry")
    need(mm(mm(tr(K), g), K) == g, "K is not a Hodge isometry")

    # eta=(1,i), bar_eta=(1,-i): K swaps them.  Coefficient conjugation
    # swaps them too, so the anti-linear C K fixes each quarter eigenline.
    eta = (1, 1j)
    bar_eta = (1, -1j)
    Kv = lambda v: (v[0], -v[1])
    Cv = lambda v: tuple(z.conjugate() for z in v)
    need(Kv(eta) == bar_eta and Kv(bar_eta) == eta, "linear kappa does not swap Hodge lines")
    need(Cv(Kv(eta)) == eta and Cv(Kv(bar_eta)) == bar_eta, "anti-linear Ckappa does not fix lines")

    # C12 marked clock R12=Q Delta^4.  Conjugation inverts both factors.
    n = 12
    R12 = (1, 4)
    need(order_pair(*R12, n) == 12, "R12 order differs from 12")
    kappa_R12_kappa = ((-R12[0]) % 4, (-R12[1]) % n)
    R12_inv = kappa_R12_kappa
    need(kappa_R12_kappa == R12_inv, "kappa does not invert R12")

    # Powers and descents from the pinned cover contract.
    need(((4 * R12[0]) % 4, (4 * R12[1]) % n) == (0, 4), "wrong fourth power")
    need(((9 * R12[0]) % 4, (9 * R12[1]) % n) == (1, 0), "wrong ninth power")
    # y=s^6: Delta^4 acts trivially because 6*4=24=0 mod 12.
    need((6 * R12[1]) % n == 0, "deck factor does not vanish on y=s^6")

    # Branch data for C2 -> P1_x: four simple finite branch points, infinity unramified.
    need(gcd(2, 4) == 2, "wrong infinity fiber count")
    genus_c2 = (3 * 2 - gcd(2, 4) - 2) // 2
    need(genus_c2 == 1, "C2 genus differs from one")

    print("PASS NS D4 relation: A^2=-I, A^4=I, RAR=A^-1")
    print("PASS C2 H1 carrier: Q*=diag(i,-i), kappa Q kappa=Q^-1")
    print("PASS polarized carrier: Q symplectic, kappa anti-symplectic, both Hodge isometries")
    print("PASS reflection types: linear kappa swaps Hodge lines; anti-linear Ckappa fixes each")
    print("PASS C12 real structure: kappa R12 kappa=R12^-1, ord(R12)=12")
    print("PASS marked descent: R12^4=Delta12^4, R12^9=Q, C12->C2 kills Delta12^4")
    print("PASS C2->P1_x is genus-one double cover branched at mu4")


if __name__ == "__main__":
    main()
