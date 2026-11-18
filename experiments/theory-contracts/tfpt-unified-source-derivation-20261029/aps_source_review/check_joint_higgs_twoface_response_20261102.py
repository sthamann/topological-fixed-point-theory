"""Conditional two-face/Higgs identities: original constructor and exact controls.

The chosen weights and collar length are diagnostics, not a physical measure.
"""
from fractions import Fraction as F
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'verification'))
import v252_full_finite_triple as v

rng = np.random.default_rng(252)
ys = [rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3)) for _ in range(4)]
mr = rng.normal(size=(3, 3)); MR = mr + mr.T
z3 = np.zeros((3, 3), complex); z2 = np.zeros((2, 2), complex)
DY = v.full_D(*ys, z3); DM = v.full_D(z3, z3, z3, z3, MR)
K = v.Kswap(); P = v.rep(1, z2, z3); gamma = v.grading()
p, ell = .25, 1.
quarks = [v.IDX[(name, ch, col, gen)] for name in ('u', 'd')
          for ch in ('L', 'R') for col in range(3) for gen in range(3)]
left = [v.IDX[(name, 'L', 0, gen)] for name in ('u', 'd') for gen in range(3)]
sing = np.r_[np.linalg.svd(ys[0], compute_uv=False), np.linalg.svd(ys[1], compute_uv=False)]

def source(H):
    h1, h2 = H
    q = np.array([[np.conj(h2), h1], [-np.conj(h1), h2]])
    a = v.rep(0, q, z3)
    b = a @ (DY @ P - P @ DY); b = b + b.conj().T
    phi = b + K @ b.conj() @ K + DM
    h = p*gamma - 1j*gamma @ phi
    assert np.max(abs(h-h.conj().T)) < 1e-12
    ev, U = np.linalg.eigh(h)
    C = (U * (1/(1+np.exp(2*ell*ev)))) @ U.conj().T
    X = (U * (1/(2*np.cosh(ell*ev)))) @ U.conj().T
    G = np.block([[np.eye(96)-C, X], [X, C]])
    assert np.max(abs(G @ G-G)) < 1e-12
    u = sum(abs(z)**2 for z in H)
    g = float(np.trace(X[np.ix_(quarks, quarks)]).real)
    expected = float(6*np.sum(1/(2*np.cosh(ell*np.sqrt(p*p+u*sing*sing)))))
    assert abs(g-expected) < 2e-12
    return G, float(np.trace(C[np.ix_(left,left)]).real), g

samples = [source(H) for H in [(0j,0j),(0j,1+0j),(.6j,.8+0j),(0j,2+0j)]]
assert len(quarks) == 36 and len(left) == 6
assert abs(samples[1][2]-samples[2][2]) < 2e-12
assert samples[0][2] > samples[1][2] > samples[3][2]
assert samples[0][1] < samples[1][1] < samples[3][1]
G0,m0,g0 = samples[0]; G1,m1,g1 = samples[1]
Gbar = (G0+G1)/2
idx = [96+i for i in left]
defect = Gbar-Gbar@Gbar
lhs = float(np.trace(defect[np.ix_(idx,idx)]).real)
variance = (m0-m1)**2/4
assert lhs+1e-12 >= variance/6 and variance > 0
assert np.linalg.eigvalsh(defect).min() > -1e-12
f0,f1 = 0.,.5
mixed_covariance = (f0-f1)*(g0-g1)/4
assert mixed_covariance < 0

# Independent exact rational saturation of the rank-six inequality.  These
# are abstract 12-mode projectors, not an identification with the TFPT source.
a,b = F(1,5),F(4,5)
var = (6*a-6*b)**2/4
exact_lhs = 6*(F(1,2)-F(1,2)**2-F(2,5)**2)
assert exact_lhs == var/6
delta = 24*var
assert exact_lhs == delta/144
assert (24*a*a+24*b*b)/2 == 8*(3*a*a+3*b*b)/2
print({'result':'PASS','original96_g': [x[2] for x in samples],
       'original96_m': [x[1] for x in samples],
       'diagnostic_projected_defect':lhs,'diagnostic_variance_over_6':variance/6,
       'diagnostic_mixed_covariance':mixed_covariance,
       'exact_saturation':str(exact_lhs), 'physical_measure_selected':False})
