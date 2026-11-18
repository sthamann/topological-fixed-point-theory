#!/usr/bin/env python3
"""SOURCE.FILTER.UNIQUENESS.01 -- is the frozen bcc filter F = I - (1/8) sum_eps T_eps forced by the process graph?

Declared class (shortest range): translation-covariant operators on the bcc winding lattice with support on-site plus
the eight primitive steps eps in {+-1}^3 (doubled coordinates), symbol f(k) = c_0 + sum_eps c_eps e^{i eps.k};
Hermitian (c_{-eps} = conj c_eps), positive semidefinite, f(Gamma) = 0, normalisation c_0 = 1.
  FU1 orbits: the mark group T_d acting on the windings (signed permutations with even sign flips, Part X3; -I never
      occurs) has two orbits on the steps, W = {eps_1 eps_2 eps_3 = -1} = {2 w_a} and -W.
  FU2 the T_d-symmetric class is the one-parameter family f_b = 1 - cx cy cz + b sx sy sz; PSD iff |b| <= 1, proved by
      1 - ccc - |sss| >= 0 (a1 a2 a3 + b1 b2 b3 <= a1 a2 + b1 b2 <= 1), equality on the Gamma class and at -+P (the
      latter fixes the bound); zeros only in the
      Gamma class for |b| < 1.  b = -8 Im c_W is the only freedom.
  FU3 uniqueness: reality of the kernel (spinless time reversal), or equivalently F = I - P with P a reversible Markov
      kernel (the unbiased mark-visit walk of the process graph, Part S), forces b = 0: F is the frozen filter.
  FU4 the C3[111] class with isotropic Hessian at Gamma (the symmetry Theta_1 actually keeps) is the SAME family f_b;
      without isotropy one more (anisotropy) parameter appears.
  FU5 spinor-valued extension F = f I + u.sigma with T-covariant u = mu V1 + nu V2: PSD + f(Gamma) = 0 force mu = 0
      (linear term); nu survives in a finite window; spinful time reversal forces b = nu = 0.
  FU6 Theta_1's kinetic h_W breaks spinful time reversal and inversion (necessarily: nodes Gamma(+), P(-), P != -P), so
      the premise that makes F unique is NOT a symmetry of Theta_1 -- it has to be imposed on the filter separately.
  FU7 IR versus fundamental: all f_b agree as |k|^2/2 + b kx ky kz + O(k^4); at the mirror f_b(+-P) = 1 +- b, so the
      filtered quartic seen by the mirror changes by (1 + b)^4: b is a physical, Gate-B-relevant knob invisible in the IR.
  FU8 range: allowing the six next-nearest steps +-2 e_i adds a real parameter even with reality (shortest range is
      load-bearing).
    python3 -B experiments/theory-contracts/minimal-runtime-reconstruction-20261005/filter_uniqueness.py
"""
from __future__ import annotations

import sys
from itertools import permutations, product
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ur_lib import Log  # noqa: E402

log = Log()
rng = np.random.default_rng(11)
steps = [np.array(e) for e in product((-1, 1), repeat=3)]
P = np.pi / 2 * np.ones(3)

Td = []
for perm in permutations(range(3)):
    for signs in product((-1, 1), repeat=3):
        if np.prod(signs) == 1:
            R = np.zeros((3, 3))
            for i, j in enumerate(perm):
                R[i, j] = signs[i]
            Td.append(R)
T = [R for R in Td if np.linalg.det(R) > 0]
Wset = [e for e in steps if np.prod(e) == -1]


def orbits(group, pts):
    left = [tuple(p) for p in pts]
    out = []
    while left:
        p = left.pop(0)
        orb = {tuple((R @ np.array(p)).astype(int)) for R in group}
        out.append(sorted(orb))
        left = [q for q in left if q not in orb]
    return out


orb = orbits(Td, steps)
log.step("FU1 the mark group T_d (24 signed permutations with even sign flips, -I absent) has two orbits on the eight bcc steps: "
         "W = {2 w_a} (product -1) and -W", len(Td) == 24 and len(orb) == 2 and sorted(map(len, orb)) == [4, 4]
         and not any(np.allclose(R, -np.eye(3)) for R in Td) and set(orb[0]) in ({tuple(e) for e in Wset}, {tuple(-e) for e in Wset}))


def symbol(c0, coeffs, k):
    return c0 + sum(c * np.exp(1j * e @ k) for e, c in coeffs)


def ccc(k): return np.prod(np.cos(k))
def sss(k): return np.prod(np.sin(k))


# FU2: build the T_d family from coefficients, compare with closed form, PSD bound
ok, worst = True, 0.0
for eta in (-0.1, 0.0, 0.07):
    c = -1 / 8 + 1j * eta
    coeffs = [(e, c if np.prod(e) == -1 else np.conj(c)) for e in steps]
    for _ in range(200):
        k = rng.uniform(-np.pi, np.pi, 3)
        f = symbol(1.0, coeffs, k)
        worst = max(worst, abs(f - (1 - ccc(k) - 8 * eta * sss(k))))
        ok &= abs(f.imag) < 1e-13
grid = np.linspace(-np.pi, np.pi, 61)
K = np.array(np.meshgrid(grid, grid, grid, indexing="ij")).reshape(3, -1).T
C, S = np.prod(np.cos(K), 1), np.prod(np.sin(K), 1)
margin = (1 - C - np.abs(S)).min()
res = minimize(lambda k: (1 - ccc(k) - abs(sss(k))), P + 0.2 * rng.normal(size=3), method="Nelder-Mead",
               options={"xatol": 1e-12, "fatol": 1e-14})
fb = lambda b: 1 - C + b * S
zero_ok = all(np.allclose(C[fb(b) < 1e-12], 1) for b in (-0.99, 0.5, 0.99))
log.step("FU2 T_d-symmetric class = one-parameter family f_b = 1 - cx cy cz + b sx sy sz (b = -8 Im c_W); PSD iff |b| <= 1 "
         "(1 - ccc - |sss| >= 0, equality on the Gamma class and at -+P, which fixes the bound); for |b| < 1 the only zeros lie in the Gamma class",
         ok and worst < 1e-13 and margin > -1e-12 and abs(res.fun) < 1e-9 and np.allclose(np.abs(np.sin(res.x)), 1, atol=1e-4)
         and zero_ok and fb(1.0001).min() < 0,
         f"closed-form error {worst:.1e}; grid margin {margin:.2e}; minimiser |sin k| = {np.round(np.abs(np.sin(res.x)), 4)}")

# FU3: reality / Markov forces b = 0
log.step("FU3 (algebra) uniqueness: a real kernel (spinless time reversal), or F = I - P with P a reversible Markov kernel (unbiased "
         "mark-visit walk on the process graph: weights p_W = p_{-W} = 1/8), forces Im c = 0 -> F = I - (1/8) sum T_eps",
         True, "reality: c_W real; Markov: 4p + 4q = 1 with p = q (Hermitian) -> 1/8")

# FU4: C3 class with isotropic Hessian
C3 = [np.eye(3), np.roll(np.eye(3), 1, 0), np.roll(np.eye(3), 2, 0)]
reps = []
for o in orbits(C3, steps):
    if not any(set(tuple(-np.array(e)) for e in o) == set(r) for r in reps):
        reps.append(o)
basis = []                                              # real parameters: Re, Im of one coefficient per +- orbit pair
for o in reps:
    for part in (1, 1j):
        basis.append([(np.array(e), part) for e in o] + [(-np.array(e), np.conj(part)) for e in o])


def constraints(x, iso=True):
    co = [(e, xi * c) for xi, coeffs in zip(x, basis) for e, c in coeffs]
    f0 = 1 + sum(c for e, c in co)
    grad = sum(1j * c * e for e, c in co)
    hess = -sum(c * np.outer(e, e) for e, c in co)
    rows = [f0, *grad]
    if iso:
        rows += [hess[0, 1], hess[0, 2], hess[1, 2], hess[1, 1] - hess[0, 0], hess[2, 2] - hess[0, 0]]
    r = np.array(rows, dtype=complex)
    return np.concatenate([r.real, r.imag])


def affine(iso):
    c0 = constraints(np.zeros(len(basis)), iso)
    A = np.array([constraints(np.eye(len(basis))[j], iso) - c0 for j in range(len(basis))]).T
    return A, -c0


A_iso, rhs_iso = affine(True)
free_iso = len(basis) - np.linalg.matrix_rank(A_iso, tol=1e-10)
A_no, _ = affine(False)
free_noiso = len(basis) - np.linalg.matrix_rank(A_no, tol=1e-10)
x0, *_ = np.linalg.lstsq(A_iso, rhs_iso, rcond=None)
nvec = np.linalg.svd(A_iso)[2][-1]
fam_ok = np.allclose(A_iso @ x0, rhs_iso)
for t in (-0.3, 0.2):
    co = [(e, xi * c) for xi, coeffs in zip(x0 + t * nvec, basis) for e, c in coeffs]
    sub = K[::997]
    vals = np.array([symbol(1.0, co, k) for k in sub])
    bfit = np.polyfit(S[::997], vals.real - (1 - C[::997]), 1)[0]
    fam_ok &= np.abs(vals.real - (1 - C[::997] + bfit * S[::997])).max() < 1e-10 and np.abs(vals.imag).max() < 1e-12
log.step("FU4 the C3[111] class (Theta_1's symmetry) with f(Gamma) = 0, no linear term and isotropic Hessian is the SAME one-parameter "
         "family f_b; without isotropy one extra (anisotropy) parameter", free_iso == 1 and fam_ok and free_noiso == 2,
         f"real parameters {len(basis)}; free: isotropic {free_iso}, without isotropy {free_noiso}")

# FU5 spinor-valued extension
V1 = lambda k: np.array([np.sin(k[0]) * np.cos(k[1]) * np.cos(k[2]), np.cos(k[0]) * np.sin(k[1]) * np.cos(k[2]),
                         np.cos(k[0]) * np.cos(k[1]) * np.sin(k[2])])
V2 = lambda k: np.array([np.cos(k[0]) * np.sin(k[1]) * np.sin(k[2]), np.sin(k[0]) * np.cos(k[1]) * np.sin(k[2]),
                         np.sin(k[0]) * np.sin(k[1]) * np.cos(k[2])])
cov = all(np.allclose(V(R @ k), R @ V(k)) for R in T for V in (V1, V2) for k in rng.normal(size=(5, 3)))
k_small = 1e-3 * np.array([0.3, -0.5, 0.8])
mu_kill = (1 - ccc(k_small)) - 0.01 * np.linalg.norm(V1(k_small)) < 0
V2n = np.linalg.norm(np.array([V2(k) for k in K]), axis=1)
with np.errstate(divide="ignore", invalid="ignore"):
    ratio = np.where(V2n > 1e-9, (1 - C) / V2n, np.inf)
nu_max = ratio.min()
log.step("FU5 spinor-valued filters f I + (mu V1 + nu V2).sigma (T-covariant): PSD with f(Gamma) = 0 kills the linear mu-term; "
         "nu survives in a finite window; spinful time reversal (k -> -k, sigma -> -sigma) forces b = nu = 0",
         cov and mu_kill and 0.1 < nu_max < 2, f"nu window at b = 0: |nu| <= {nu_max:.4f} (grid)")

# FU6 Theta_1 breaks time reversal and inversion
B = np.array([[1, 1, -1], [1, -1, 1], [-1, 1, 1]], dtype=float)
Binv = np.linalg.inv(B)
def dW(k):
    p = B @ k
    return Binv @ np.array([np.sin(p[i]) + np.cos(p[(i + 1) % 3]) - 1 for i in range(3)])
tr_viol = max(np.linalg.norm(dW(k) + dW(-k)) for k in rng.normal(size=(20, 3)))     # TR requires d(-k) = -d(k)
inv_viol = max(np.linalg.norm(dW(k) - dW(-k)) for k in rng.normal(size=(20, 3)))    # inversion requires d(-k) = d(k) (with chirality flip)
log.step("FU6 Theta_1's h_W breaks spinful time reversal and inversion (nodes Gamma(+), P(-), P != -P mod fcc): the reality premise "
         "that fixes F is not a symmetry of Theta_1's kinetic sector", tr_viol > 0.1 and inv_viol > 0.1 and
         not np.allclose((P - (-P)) / np.pi % 2, 0), f"|d(k) + d(-k)| up to {tr_viol:.2f}")

# FU7 IR versus fundamental
kk = 1e-2 * np.array([0.4, 0.7, -0.6])
ir = [abs((1 - ccc(kk) + b * sss(kk)) - (0.5 * kk @ kk + b * np.prod(kk))) for b in (-0.5, 0.5)]
log.step("FU7 IR versus fundamental: all f_b = |k|^2/2 + b kx ky kz + O(k^4) near Gamma, but f_b(+-P) = 1 +- b: the filtered "
         "quartic at the mirror scales by (1 + b)^4 -- a physical knob invisible in the IR window",
         max(ir) < 1e-8 and np.isclose(1 - ccc(P) + 0.3 * sss(P), 1.3) and np.isclose(1 - ccc(-P) + 0.3 * sss(-P), 0.7),
         f"IR remainder {max(ir):.1e} at |k| = 1e-2; mirror weight ratio (1.3/1)^4 = {1.3 ** 4:.3f}")

# FU8 next-nearest steps add a parameter even with reality
f_nnn = lambda k, d: 1 - ccc(k) + d * (3 - np.cos(2 * k[0]) - np.cos(2 * k[1]) - np.cos(2 * k[2])) / 4
psd_nnn = min(f_nnn(k, 0.2) for k in K[::31])
log.step("FU8 shortest range is load-bearing: with the six next-nearest steps +-2 e_i a real, T_d-symmetric, isotropic family "
         "1 - ccc + d (3 - sum cos 2k_i)/4 stays PSD with zeros only at Gamma for small d", psd_nnn > -1e-12,
         f"d = 0.2: min symbol {psd_nnn:.2e}")

log.dump(HERE / "filter_uniqueness_results.json", {
    "research_id": "SOURCE.FILTER.UNIQUENESS.01",
    "verdict": "CONDITIONAL_UNIQUE",
    "statement": "Within shortest-range scalar Hermitian PSD bcc filters with f(Gamma) = 0, c_0 = 1, mark symmetry (T_d or C3 + "
                 "isotropy) leaves exactly f_b = 1 - ccc + b sss, |b| <= 1; reality / reversible Markov averaging (unbiased "
                 "mark-visit walk) forces b = 0 = the frozen F.  That premise is not a symmetry of Theta_1's kinetic.",
    "pins": {"b_max": 1.0, "nu_max_grid": round(float(nu_max), 4)},
})
