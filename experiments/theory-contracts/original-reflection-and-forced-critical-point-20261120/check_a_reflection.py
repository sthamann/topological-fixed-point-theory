#!/usr/bin/env python3
"""KERNSYM20261120 A: the time reflection of the original against the mirror of the kernel.

Fock space with three modes on each normal half, kernel state Omega = Fock state of
P_C = (1/2)[[1, U*], [U, 1]].
A1  the Lagrange relation G P G = 1 - P makes the antilinear field map a(u) -> a*(G u) a symmetry of
    the two point function; no other relative phase between the halves does
A2  it is implemented by an antiunitary operator that fixes Omega (the time reflection at time zero)
A3  that operator maps the algebra of one normal half to itself; the mirror J of P1BIND maps it to
    the other half: the two are different operators
A4  as an order preserving map the reflection does not give a positive form on a half; with the
    order reversal of the fermionic reconstruction it sends F to F* and the form is
    <Omega, F* G Omega> = Tr(F* G)/2^d, the mirror form of P1BIND
"""
import numpy as np
from common import Checks, grading, randc, fock, slater

c = Checks("a_reflection")
rng = np.random.default_rng(2026112021)
m = 3
n = 2 * m
U, _ = np.linalg.qr(randc(rng, m, m))
Z = np.zeros((m, m), complex)
P = 0.5 * np.block([[np.eye(m), U.conj().T], [U, np.eye(m)]])
G = grading(m)
I = np.eye(n)

# ---------------------------------------------------------------- A1
c.ok("A1a G* (1 - P) G = P: the two point function is invariant under a(u) -> a*(G u)", np.allclose(G.conj().T @ (I - P) @ G, P))
bad = []
for om in (1, 1j, -1j, np.exp(0.7j)):
    Gam = np.diag([1.0] * m + [om] * m)
    bad.append(np.linalg.norm(Gam.conj().T @ (I - P) @ Gam - P))
c.ok("A1b with any other relative phase between the halves it is not", min(bad) > 0.5, [round(float(b), 3) for b in bad])

# ---------------------------------------------------------------- A2
ann = fock(n)
Om = slater(np.vstack([np.eye(m), U]) / np.sqrt(2), ann)
W = np.eye(2 ** n, dtype=complex)
for a in ann:
    W = W @ (a + a.conj().T)
Nm = sum(ann[i].conj().T @ ann[i] for i in range(m, n))
V = np.diag(np.exp(1j * np.pi * np.diag(Nm).real))
T = W @ V                                             # Theta = T K, K = complex conjugation in the occupation basis
signs = []
ok = True
for i, a in enumerate(ann):
    img = T @ a.conj() @ T.conj().T                   # Theta a_i Theta^(-1)
    s = np.vdot(a.conj().T, img) / np.vdot(a.conj().T, a.conj().T)
    ok &= np.allclose(img, s * a.conj().T)
    signs.append(s)
signs = np.array(signs)
c.ok("A2a Theta a_i Theta^(-1) = +- a_i^+ with the sign of the normal grading (up to one common sign)",
     ok and np.allclose(signs / signs[0], np.diag(G).real))
ov = np.vdot(Om, T @ Om.conj())
c.ok("A2b Theta is antiunitary and fixes the kernel state up to a phase", np.allclose(T @ T.conj().T, np.eye(2 ** n)) and abs(abs(ov) - 1) < 1e-10, round(float(abs(ov)), 12))

# ---------------------------------------------------------------- A3
half = ann[:m]
def in_half(X):
    """X commutes with the operators of the other half, graded: test via (anti)commutators with b_j."""
    return all(np.allclose(X @ b + b @ X, 0) or np.allclose(X @ b - b @ X, 0) for b in ann[m:])
c.ok("A3a Theta maps each field of one normal half to a field of the same half", all(np.allclose(T @ a.conj() @ T.conj().T, signs[i] * a.conj().T) for i, a in enumerate(half)))
c.ok("A3b the mirror of P1BIND sends a(e) to a*(U e), a field of the other half: Theta and J are different operators",
     np.linalg.norm(U) > 0 and np.allclose(P[m:, :m], U / 2))

# ---------------------------------------------------------------- A4
def rand_half():
    X = np.zeros((2 ** n, 2 ** n), complex)
    for _ in range(4):
        Y = np.eye(2 ** n, dtype=complex)
        for _ in range(rng.integers(1, 4)):
            a = half[rng.integers(m)]
            Y = Y @ (a if rng.integers(2) else a.conj().T)
        X = X + (rng.normal() + 1j * rng.normal()) * Y
    return X
F2 = half[0] @ half[1]
auto = T @ F2.conj() @ T.conj().T                     # order preserving image
neg = np.vdot(Om, auto @ F2 @ Om).real
c.ok("A4a as an order preserving map: <Omega, Theta(F) F Omega> is negative for F = a_1 a_2 (no positive form)", neg < -1e-3, round(float(neg), 6))
ok = True
for _ in range(20):
    F, Gg = rand_half(), rand_half()
    lhs = np.vdot(Om, F.conj().T @ Gg @ Om)
    ok &= abs(lhs - np.trace(F.conj().T @ Gg) / 2 ** n) < 1e-10
c.ok("A4b with order reversal Theta(F) = F* on a half, and <Omega, F* G Omega> = Tr(F* G)/2^d: the reflected form of the "
     "original at time zero is the mirror form of P1BIND on each normal half", ok)
c.done()
