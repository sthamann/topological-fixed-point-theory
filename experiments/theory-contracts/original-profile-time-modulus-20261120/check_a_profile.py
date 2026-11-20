#!/usr/bin/env python3
"""KERNORIG20261120 A: the profile of the original.

The original reads its geometric branch from the "profile-free" family
    Z_rel(s) = Tr(|D_rel|^(-2s) - |D_ref|^(-2s)) + (i pi/2) Delta eta
(residue at s = 1, finite part at s = 0).  Checked here on finite matrices:
A1  half the s-derivative at s = 0 of that family is the primitive action with f = -log|x|
A2  so exp(-(1/2) dZ/ds(0) - phase) is the relative fermion determinant
A3  the value at s = 0 carries only the phase term
A4  the master functional of the original has + Re log det(D^2 + mu^2): at mu = 0 this is the
    profile +2 log|x|, the opposite sign and twice the weight of one fermion determinant
"""
import numpy as np
from common import Checks, randc

c = Checks("a_profile")
rng = np.random.default_rng(2026112011)

def herm(k):
    a = randc(rng, k, k)
    return (a + a.conj().T) / 2

def eta(D):
    w = np.linalg.eigvalsh(D)
    return int(np.sum(w > 0) - np.sum(w < 0))

def zrel(s, Drel, Dref):
    return np.sum(np.abs(np.linalg.eigvalsh(Drel)) ** (-2 * s)) - np.sum(np.abs(np.linalg.eigvalsh(Dref)) ** (-2 * s))

k = 8
a1 = a2 = a3 = a4 = True
for _ in range(20):
    Drel, Dref = herm(k), herm(k)
    h = 1e-5
    dz = (zrel(h, Drel, Dref) - zrel(-h, Drel, Dref)) / (2 * h)
    slog = -np.sum(np.log(np.abs(np.linalg.eigvalsh(Drel)))) + np.sum(np.log(np.abs(np.linalg.eigvalsh(Dref))))
    a1 &= abs(0.5 * dz - slog) < 1e-7
    phase = 0.5j * np.pi * (eta(Drel) - eta(Dref))
    ratio = np.linalg.det(Drel) / np.linalg.det(Dref)
    a2 &= abs(np.exp(-0.5 * dz - phase) - ratio) < 1e-6 * abs(ratio)
    a3 &= abs(zrel(0.0, Drel, Dref)) < 1e-12
    plus = np.sum(np.log(np.linalg.eigvalsh(Drel @ Drel))) - np.sum(np.log(np.linalg.eigvalsh(Dref @ Dref)))
    a4 &= abs(plus - (-2) * slog) < 1e-8
c.ok("A1 (1/2) d/ds Z_rel(s) at s = 0 equals tr f(D_rel) - tr f(D_ref) with f = -log|x|", a1)
c.ok("A2 with the phase term of the original this is minus the logarithm of det(D_rel)/det(D_ref)", a2)
c.ok("A3 at s = 0 the family itself carries only the phase term (equal ranks)", a3)
c.ok("A4 Re log det(D^2) of the master functional = -2 times that action: opposite sign, double weight", a4)
c.done()
