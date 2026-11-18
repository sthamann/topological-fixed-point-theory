from itertools import product
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss

base = Path(__file__).resolve().parent
masses = json.loads((base / "quarter-reference-result.json").read_text())["masses"]
time = json.loads((base / "r113-temporal-cross.json").read_text())["densities"]
mp.mp.dps = 90
p0 = mp.mpf(1) / 4
dt = mp.mpf(3) / 10
energies = [mp.sqrt(p0**2 + mp.mpf(m)**2) for m in masses]
overlaps = [(1 + p0 / e) / 2 for e in energies]
rates = sorted(2 * sum(e * b for e, b in zip(energies, bits)) for bits in product((0, 1), repeat=3))
assert len(set(rates)) == 8

def cap(t):
    return mp.fprod(a + (1 - a) * mp.exp(-2 * t * e) for a, e in zip(overlaps, energies))

moments = [cap(n * dt) for n in range(24)]
def hankel(n):
    return mp.matrix([[moments[i+j] for j in range(n)] for i in range(n)])
g8 = hankel(8)
eigenvalues = list(mp.eigsy(g8, eigvals_only=True))
assert eigenvalues[0] > 0
assert abs(mp.det(hankel(9))) < mp.mpf("1e-75")
coefficients = mp.lu_solve(g8, mp.matrix([-moments[i+8] for i in range(8)]))
predicted = moments[:16]
for n in range(16, 24):
    predicted.append(-sum(coefficients[j] * predicted[n-8+j] for j in range(8)))
held_out = [abs(moments[n] - predicted[n]) for n in range(16, 24)]
assert max(held_out) < mp.mpf("1e-70")
prior = json.loads((base / "whole-strategy-contacts.json").read_text())["reference_slater_cap"]
assert abs(float(mp.log(cap(dt))) - prior["native_first_mode_log_cap"]) < 2e-12
print("R032 one mode, three recorded massive pairs")
print("smallest rate separation", mp.nstr(min(np.diff(rates)), 12))
print("Hankel rank", len(rates), "first eigenvalue", mp.nstr(eigenvalues[0], 8))
print("held-out recurrence max error", mp.nstr(max(held_out), 8), "(high precision on diagnostic inputs)")
print("cap values", *(mp.nstr(v, 13) for v in moments[:3]))
print("two-cut return determinant", mp.nstr(moments[2] - moments[1]**2, 12))

T = time["time"]
B = time["squared_charged_weight"]
nodes, weights = leggauss(100)
width = T / 20
bump = weights * np.exp(-1/(1-nodes**2))
bump = bump * (T / bump.sum())
# Smooth reflection-even profiles with a strict gap, each integral exactly T
# in the quadrature. One has two reflected supports and the other has one.
x = np.concatenate([(T/5 + width*nodes), (4*T/5 + width*nodes)])
wx = np.concatenate([bump/2, bump/2])
z = T/2 + width*nodes
distance = np.abs(x[:, None] - z[None, :])
assert distance.min() > T / 6
outer = wx[:, None] * bump[None, :]

def score(p):
    # cosh(p(T-2d))/sinh(pT) in stable exponential form
    K = (np.exp(-2*p*distance)+np.exp(-2*p*(T-distance)))/(-np.expm1(-2*p*T))
    g = B * np.sum(outer * K)
    return g*np.tanh(T*p/2), -2*g/(1+np.exp(T*p))

for p in (1e-4, 1e-5, 1e-6):
    sewn, polar = score(p)
    print(f"p={p:g} sewn cross={sewn:.9f}  p*polar={p*polar:.9f}")
assert abs(score(1e-6)[0]/(B*T**2/2)-1)<1e-8
assert abs(1e-6*score(1e-6)[1]/(-B*T)-1)<1e-5
print("predicted sewn p->0", B*T**2/2, "polar coefficient before log radius (half jet, per r)", -B*T)
