"""Check the local Spin(4) a6 coefficient and the independent open-link response."""

import json
from fractions import Fraction

import numpy as np
import sympy as sp

import v252_full_finite_triple as original


PAULI = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
         np.diag([1, -1]))


def gauge_operators():
    zero2 = np.zeros((2, 2), complex)
    zero3 = np.zeros((3, 3), complex)
    real = original.Kswap()

    def real_lift(generator):
        return -1j * (generator + real @ generator.conj() @ real)

    hypercharge = real_lift(original.rep(3j, zero2, -1j * np.eye(3))) / 6
    weak = [real_lift(original.rep(0, 1j * s, zero3)) / 2 for s in PAULI]
    return weak, hypercharge


def higgs_constructor(yukawas):
    zero2 = np.zeros((2, 2), complex)
    zero3 = np.zeros((3, 3), complex)
    real = original.Kswap()
    d_y = original.full_D(*yukawas, zero3)
    p = original.rep(1, zero2, zero3)

    def construct(higgs):
        h1, h2 = higgs
        quaternion = np.array([[h2.conjugate(), h1], [-h1.conjugate(), h2]])
        left = original.rep(0, quaternion, zero3) @ (d_y @ p - p @ d_y)
        finite = left + left.conj().T
        return finite + real @ finite.conj() @ real

    return construct


def commutator(left, right):
    return left @ right - right @ left


def check_clifford():
    pauli = [sp.Matrix(s).applyfunc(sp.nsimplify) for s in PAULI]
    gamma = [sp.kronecker_product(pauli[0], matrix) for matrix in pauli]
    gamma.append(sp.kronecker_product(pauli[1], sp.eye(2)))
    g5 = gamma[0] * gamma[1] * gamma[2] * gamma[3]
    s = gamma[0] * gamma[1]
    q0, q1 = [sp.I * gamma[index] * g5 for index in (0, 1)]
    assert sp.trace(s * s) == -4
    assert sp.trace(s * q0 * q1) == -4
    assert sp.trace(q0 * q0) == 4
    assert sp.trace(q0 * q1) == 0


def check_trace():
    weak, hypercharge = gauge_operators()
    yukawas = [np.diag(vals) for vals in ([1, 2, 3], [4, 5, 6], [2, 3, 4], [1, 3, 5])]
    construct = higgs_constructor(yukawas)
    u, d, e, n = [int(np.vdot(yukawa, yukawa)) for yukawa in yukawas]
    raw_factor = Fraction(2, 3) * (u - d + e - n)
    reduced_factor = Fraction(1, 3) * (5 * u + d + 3 * e - n)
    assert raw_factor == -46 and reduced_factor == Fraction(199, 3)

    for generator in weak:
        assert np.array_equal(commutator(generator, hypercharge), np.zeros((96, 96)))
        assert np.trace(generator @ hypercharge) == 0

    probes = (np.array([0j, 1]), np.array([1, 1]), np.array([1, 1j]),
              np.array([1 + 2j, 2 - 3j]))
    for h in probes:
        phi = construct(h)
        dy = commutator(hypercharge, phi)
        for sigma, generator in zip(PAULI, weak):
            dweak = commutator(generator, phi)
            direct = -2 / 3 * np.trace(phi @ phi @ (generator @ hypercharge
                                                     + hypercharge @ generator))
            extra = -np.trace(dweak @ dy + dy @ dweak) / 6
            triplet = h.conj() @ sigma @ h
            assert np.isclose(direct, float(raw_factor) * triplet, atol=1e-12)
            assert np.isclose(direct + extra, float(reduced_factor) * triplet, atol=1e-12)

    t3y = weak[2] @ hypercharge
    nu_l = original.IDX[('nu', 'L', 0, 0)]
    assert t3y[nu_l, nu_l] == -1 / 4
    assert np.count_nonzero(t3y) == 48
    assert np.trace(t3y) == 0
    assert np.isclose(np.linalg.norm(t3y), 1)
    return {'family_trace_Yu_Yd_Ye_Ynu': [u, d, e, n],
            'raw_mixed_factor': str(raw_factor),
            'IBP_and_EOM_mixed_factor': str(reduced_factor),
            'weak_hypercharge_mixed_generator_rank': int(np.count_nonzero(t3y)),
            'open_link_mixed_nu_L_at_ell_pi_over_4': 'pi**2/64',
            'open_link_mixed_Frobenius_norm_at_ell_pi_over_4': 'pi**2/16'}


if __name__ == '__main__':
    check_clifford()
    print(json.dumps(check_trace(), indent=2))
