#!/usr/bin/env python3
"""Conditional original-96 sector reduction and diagnostic number-phase complement."""
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import null_space

root = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(root / 'verification'))
import v252_full_finite_triple as v

rng = np.random.default_rng(252)
ys = [rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3)) for _ in range(4)]
mr = rng.normal(size=(3, 3)); majorana = mr + mr.T
z3 = np.zeros((3, 3), complex); z2 = np.zeros((2, 2), complex)
dy = v.full_D(*ys, z3); dm = v.full_D(z3, z3, z3, z3, majorana)
swap = v.Kswap(); projector = v.rep(1, z2, z3); gamma = v.grading()
evals, family = np.linalg.eigh(ys[0].conj().T @ ys[0])
y = np.sqrt(evals[-1]); right = np.zeros((96, 3), complex)
for color in range(3):
    for f in range(3):
        right[v.IDX[('u', 'R', color, f)], color] = family[f, -1]

p = .25; ell = 1.; z = 1j
for x in (.5, 1., 2.):
    h1, h2 = 0j, complex(x)
    q = np.array([[h2.conjugate(), h1], [-h1.conjugate(), h2]])
    a = v.rep(0, q, z3); b = a @ (dy @ projector - projector @ dy)
    b += b.conj().T
    phi = b + swap @ b.conj() @ swap + dm
    h = p*gamma - 1j*gamma @ phi
    left = (h + p*np.eye(96)) @ right / (x*y)
    sector = np.column_stack((right, left)); complement = null_space(sector.conj().T)
    ev, u = np.linalg.eigh(h)
    occ = 1/(1 + np.exp(2*ell*ev)); cov = (u*occ) @ u.conj().T
    c6 = sector.conj().T @ cov @ sector; c90 = complement.conj().T @ cov @ complement
    errors = [np.max(abs(sector.conj().T @ sector - np.eye(6))),
              np.max(abs(complement.conj().T @ h @ sector)),
              np.max(abs(complement.conj().T @ cov @ sector))]
    local = np.linalg.det(np.eye(6) - (1-z)*c6)
    factor = np.linalg.det(np.eye(90) - (1-z)*c90)
    full = np.prod(1-occ+z*occ)
    errors.append(abs(full - factor*local))
    if max(errors) > 2e-12:
        raise RuntimeError(f'original-96 reduction failed: {errors}')
    print(f'H={x}: errors={errors}; complement={factor}; local={local}; full={full}')
print('PASS: original 96 block and diagnostic complement; conditional source only')
