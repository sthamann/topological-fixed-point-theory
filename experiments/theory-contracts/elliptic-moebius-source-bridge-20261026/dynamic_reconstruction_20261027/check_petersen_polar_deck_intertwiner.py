#!/usr/bin/env python3
"""Exact narrow reproduction of the Petersen polar/deck dictionary."""

import sympy as sp


R = sp.Rational


def check(name, condition):
    if condition is not True and condition != True:
        raise RuntimeError("FAIL: " + name)
    print("PASS", name)


def block(a, b, c, d):
    return a.row_join(b).col_join(c.row_join(d))


def main():
    P = sp.Matrix([
        [0, 1, 0],
        [R(1, 3), 0, R(2, 3)],
        [0, R(1, 3), R(2, 3)],
    ])
    D = sp.diag(R(1, 10), R(3, 10), R(6, 10))
    sqrt_D = sp.diag(sp.sqrt(R(1, 10)), sp.sqrt(R(3, 10)),
                     sp.sqrt(R(6, 10)))

    f0 = sp.Matrix([1, 1, 1])
    q = sp.Matrix([6, -4, 1])
    r = sp.Matrix([3, 1, -1])
    B = sp.Matrix.hstack(f0, q / 3, sp.sqrt(5) * r / 3)
    check("marked source basis is H_pi orthonormal", B.T * D * B == sp.eye(3))

    e0 = sp.Matrix([1, 0, 0])
    ep = sp.Matrix([0, 1, 1]) / sp.sqrt(2)
    em = sp.Matrix([0, 1, -1]) / sp.sqrt(2)
    E = sp.Matrix.hstack(e0, ep, em)

    # U acts between Euclidean frames.  W = U sqrt(D) is the actual map
    # H_pi -> R^3 and satisfies W^T W = D.
    U = sp.simplify(E * (sqrt_D * B).T)
    W = sp.simplify(U * sqrt_D)
    check("U is orthogonal", sp.simplify(U * U.T) == sp.eye(3))
    check("W is an isometry H_pi to Euclidean R3", sp.simplify(W.T * W) == D)
    check("W sends q/3 to e+", sp.simplify(W * (q / 3) - ep) == sp.zeros(3, 1))

    abs_P = R(4, 5) * P**2 - R(1, 15) * P + R(4, 15) * sp.eye(3)
    J_P = -R(6, 5) * P**2 + R(8, 5) * P + R(3, 5) * sp.eye(3)
    M = sp.Matrix([
        [1, 0, 0],
        [0, R(1, 2), R(1, 6)],
        [0, R(1, 6), R(1, 2)],
    ])
    Delta = sp.Matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
    Z = sp.diag(1, -1, -1)
    K = Z * Delta

    check("polar decomposition", J_P**2 == sp.eye(3) and J_P * abs_P == P)
    check("positive transfer conjugacy", sp.simplify(W * abs_P * W.inv()) == M)
    check("polar sign becomes K=Z Delta", sp.simplify(W * J_P * W.inv()) == K)
    check("literal deck differs on the 2/3 line", K * ep == -ep and Delta * ep == ep)
    check("signed one-step dictionary", sp.simplify(W * P * W.inv()) == K * M)

    Y_shell = sp.diag(6, -4, 1)
    Y_target = sp.simplify(W * Y_shell * W.inv())
    check("charge spectrum 6,-4,1 is retained",
          sp.factor(Y_target.charpoly().as_expr()) ==
          (sp.Symbol("lambda") - 6) * (sp.Symbol("lambda") - 1) *
          (sp.Symbol("lambda") + 4))
    charge_commutator = Y_shell * P - P * Y_shell
    expected_commutator = sp.Matrix([
        [0, 10, 0],
        [-R(10, 3), 0, -R(10, 3)],
        [0, R(5, 3), 0],
    ])
    check("shell charge commutator has rank two",
          charge_commutator == expected_commutator and
          charge_commutator.rank() == 2)

    even_defect = sp.simplify(K * W * charge_commutator * W.inv())
    odd_defect = sp.simplify(W * charge_commutator * W.inv() * K)
    expected_even_defect = sp.Matrix([
        [0, -5 * sp.sqrt(2) / 2, -5 * sp.sqrt(2) / 2],
        [-5 * sp.sqrt(2) / 2, -5 * sp.sqrt(5) / 3, 0],
        [-5 * sp.sqrt(2) / 2, 0, 5 * sp.sqrt(5) / 3],
    ])
    check("moving charge defect alternates sign",
          even_defect == expected_even_defect and odd_defect == -even_defect)
    check("moving charge defect rank and Frobenius norm",
          even_defect.rank() == 2 and
          sp.trace(even_defect.T * even_defect) == R(700, 9))

    zero = sp.zeros(3)
    ident = sp.eye(3)
    P_hat = block(zero, P, P, zero)
    G = block(ident, zero, zero, J_P)
    W_hat = sp.diag(W, W) * G
    target = block(zero, M, M, zero)
    tau = block(zero, ident, ident, zero)
    tau_K = block(zero, K, K, zero)
    check("double-cover block conjugacy",
          sp.simplify(W_hat * P_hat * W_hat.inv()) == target)
    check("transported sheet exchange is K twisted",
          sp.simplify(W_hat * tau * W_hat.inv()) == tau_K)

    print("PASS_EXACT_PETERSEN_POLAR_DECK_REVIEW")


if __name__ == "__main__":
    main()
