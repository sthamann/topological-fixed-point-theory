#!/usr/bin/env python3
"""R064 original 60-record channel: exact spectrum certificate and parent.

The original source builder is imported unchanged. Integer matrix identities
certify the Gram polynomial; a rank over F_101 proves full product span over C.
The time and process examples are numerical controls of explicit formulas.
This is not the distinct five-Kraus Hodge channel from R088.
"""
import importlib.util
from pathlib import Path

import numpy as np

source = Path(__file__).resolve().parents[1] / "source_channel.py"
spec = importlib.util.spec_from_file_location("original_source", source)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
rays = mod.source_rays()
U, _, _ = mod.reflection_actions(rays)
assert len(U) == 60 and all(np.array_equal(a @ a, np.eye(10)) for a in U)

# Every entry is a Gaussian integer of modulus at most one. The Gram matrix
# is real integral. Integer matrix multiplication below has a maximum entry
# far below the int64 limit, so the polynomial identity is exact arithmetic.
M = np.stack([a.reshape(-1) for a in U], axis=1)
Gcomplex = M.conj().T @ M
assert np.array_equal(Gcomplex.imag, np.zeros((60, 60)))
assert np.array_equal(Gcomplex.real, np.rint(Gcomplex.real))
G = Gcomplex.real.astype(np.int64)
I = np.eye(60, dtype=np.int64)
assert np.max(np.abs(G)) <= 10
assert 60**3 * 10 * 18 * 26 * 106 < np.iinfo(np.int64).max
P = ((G @ (G - 8 * I)) @ (G - 16 * I)) @ (G - 96 * I)
assert np.array_equal(P, np.zeros_like(P))
assert np.trace(G) == 600
assert np.trace(G @ G) == 14400
assert np.trace((G @ G) @ G) == 944640
# P has simple roots. The four exact power sums uniquely give multiplicities.
assert 0*5 + 8*45 + 16*9 + 96 == np.trace(G)
assert 8**2*45 + 16**2*9 + 96**2 == np.trace(G @ G)
assert 8**3*45 + 16**3*9 + 96**3 == np.trace((G @ G) @ G)
print("Gram eigenvalues and multiplicities: 0^5, 8^45, 16^9, 96^1")
print("normalized Choi spectrum: 0^45, (1/75)^45, (2/75)^9, (4/25)^1")

# Mod 101, 10 is a square root of -1. Full row rank in this finite field
# certifies full rank over Gaussian rationals and hence over C.
def rank_mod_101(a):
    a = np.asarray(a, dtype=np.int64) % 101
    row = 0
    for col in range(a.shape[1]):
        candidates = np.flatnonzero(a[row:, col])
        if not len(candidates):
            continue
        pivot = row + int(candidates[0])
        a[[row, pivot]] = a[[pivot, row]]
        a[row] = a[row] * pow(int(a[row, col]), -1, 101) % 101
        if row + 1 < len(a):
            a[row+1:] = (a[row+1:] - a[row+1:, col:col+1]*a[row]) % 101
        row += 1
        if row == len(a):
            break
    return row

products = np.stack([(a @ b).reshape(-1) for a in U for b in U], axis=1)
assert np.array_equal(products.real, np.rint(products.real))
assert np.array_equal(products.imag, np.rint(products.imag))
assert rank_mod_101(np.rint(products.real).astype(np.int64) +
                    10*np.rint(products.imag).astype(np.int64)) == 100
print("rank of Choi: 55; product-Kraus rank: 100 (mod 101 certificate)")

# The source mean makes the joint controlled involution's partial trace
# a mixture of Id and Phi. No 600-dimensional exponential is required.
assert np.array_equal(sum(U), 24*np.eye(10))
assert all(np.array_equal(a, a.conj().T) for a in U)
rho = np.zeros((10, 10), dtype=np.complex128)
rho[0, 0] = 1
def phi(r):
    return sum(a @ r @ a.conj().T for a in U) / 60

theta = 0.371
K = [np.exp(-1j*theta)*(np.cos(theta)*np.eye(10) +
                         1j*np.sin(theta)*a)/np.sqrt(60) for a in U]
evolved = sum(k @ rho @ k.conj().T for k in K)
assert np.allclose(sum(k.conj().T @ k for k in K), np.eye(10), atol=1e-14)
assert np.allclose(evolved, np.cos(theta)**2*rho + np.sin(theta)**2*phi(rho), atol=1e-14)
K_end = [np.exp(-.5j*np.pi)*(1j*a)/np.sqrt(60) for a in U]
assert all(np.allclose(k, a/np.sqrt(60), atol=1e-15) for k, a in zip(K_end, U))
reused = sum(a @ a @ rho @ a.conj().T @ a.conj().T for a in U)/60
fresh = phi(phi(rho))
assert np.allclose(reused, rho, atol=1e-14)
distance = np.linalg.norm(fresh - reused, ord="fro")
assert distance > 0.1
print(f"parent endpoint and reuse versus fresh records: Frobenius distance {distance:.12f}")
print("PASS: finite R064 mathematical channel; physical source and timescale not identified")
