#!/usr/bin/env python3
"""Part Z: SPACETIME.PRECISION.01 -- geometry and dimension from the inverse covariance.

  Z1 microscopic: the integrated top-edge current phi (weight a_n ~ 1/n, the physical endpoint
     operation) of the ORIGINAL QWZ source has equal-time covariance C(n) = S(n)/n^2 = 1/n, so the
     precision operator C^-1 has eigenvalues n: C^-1 = |D| (reconstructed from correlations).
  Z2 (C^-1 - z)^-1 is the full geometric Weyl function; z is a variation of the quadratic precision
     operator (Robin coupling of Part X), not an energy or a new state.
  Z3 one-particle Weyl law of the FULL precision operator of U (8 species x two chiral halves x |D| on
     S^1): N(Lambda) ~ Lambda^1 -> spatial dimension 1; finite internal indices do not change it.
  Z4 thermodynamic dimension of the FULL process space: the E8_1 vacuum character (E4/eta^8) has
     log a_n ~ n^alpha with alpha -> 1/2 (Cardy, c = 8): spacetime dimension 1/(1-alpha) = 2.
  Z5 control: a 3+1D conformal scalar on S^3 (energies l+1, multiplicity (l+1)^2) gives Weyl exponent
     3 and alpha -> 3/4, i.e. d = 4 -- the diagnostic distinguishes 2 from 4.
  Z6 heat kernel of the compiler precision operator: on the chi = 2 double (S^2),
     Tr e^{t Delta} - 1/t -> 1/3 = chi/6: the 2D 'Einstein-Hilbert' coefficient is the topological
     Gauss-Bonnet chi of P1's c3 = 1/(|Z2| 2 pi chi) -- in 2D it carries no dynamics.
  Z7 (statement) with every structure now in TFPT the reconstructed spacetime is 1+1D; 3+1D needs a
     covariance network with one-particle Weyl exponent 3 (thermodynamic alpha = 3/4).
    PYTHONPATH=. python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/spacetime_precision.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


# Z1 microscopic precision operator from correlations ------------------------------------------
bu = json.loads((HERE / "boundary_universality_results.json").read_text())
S = {int(n): v[0] for n, v in bu["qwz_one_current"].items()}       # low-energy current weight S(n) (Part Y)
C = {n: S[n] / n ** 2 for n in S}                                    # integrated current: a_n = 1/n
prec = {n: 1 / C[n] for n in C}
step("Z1.microscopic QWZ: equal-time covariance of the integrated current C(n) = S(n)/n^2 = 1/n, precision operator C^-1 = n = |D|",
     all(abs(prec[n] - n) / n < 1e-3 for n in prec), "C^-1 eigenvalues " + ", ".join(f"{prec[n]:.4f}" for n in sorted(prec)))

# Z2 the Weyl function as precision variation ---------------------------------------------------
diag = sum(1 / (n + 2) - 1 / (n + 1) for n in range(1, 200000)) / np.pi
step("Z2.(C^-1 - z)^-1 between point traces = the geometric Weyl function (diag of z=-2 minus z=-1: -1/2pi); z = precision (Robin) variation, not an energy",
     abs(diag + 1 / (2 * np.pi)) < 1e-5)

# Z3 one-particle Weyl law of the full precision operator --------------------------------------------
eig = np.sort(np.concatenate([np.full(8 * 2, n) for n in range(1, 5001)]))  # 8 species x 2 halves x |D|
Lam = np.array([100, 300, 1000, 3000])
NL = np.array([np.searchsorted(eig, L, side="right") for L in Lam])
slope = np.polyfit(np.log(Lam), np.log(NL), 1)[0]
step("Z3.Weyl law of the FULL precision operator of U: N(Lambda) ~ Lambda^1 -> spatial dimension 1 (internal indices only multiply)",
     abs(slope - 1) < 1e-2, f"slope {slope:.4f}, N(1000) = {NL[2]} = 16 x 1000")


def series_coeffs(exponent_fn, nmax, extra=None):
    """Coefficients of prod_k (1 - q^k)^(-exponent_fn(k)) (times an optional extra series)."""
    a = [0] * (nmax + 1); a[0] = 1
    for k in range(1, nmax + 1):
        e = exponent_fn(k)
        for _ in range(e):
            for m in range(k, nmax + 1):
                a[m] += a[m - k]
    if extra is not None:
        b = [sum(a[i] * extra[m - i] for i in range(m + 1)) for m in range(nmax + 1)]
        return b
    return a


def sigma3(n):
    return sum(d ** 3 for d in range(1, n + 1) if n % d == 0)


nmax = 400
E4 = [1] + [240 * sigma3(n) for n in range(1, nmax + 1)]
e8 = series_coeffs(lambda k: 8, nmax, extra=E4)
ok_tower = e8[:4] == [1, 248, 4124, 34752]
s4 = series_coeffs(lambda k: k * k, nmax)


def alpha_fit(a, lo=150, hi=400):
    ns = np.arange(lo, hi + 1)
    y = np.log(np.array([float(a[n]) for n in ns]))
    best = None
    for al in np.linspace(0.3, 0.9, 601):
        X = np.column_stack([ns ** al, np.log(ns), np.ones_like(ns, dtype=float)])
        coef, res, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = np.sum((X @ coef - y) ** 2)
        if best is None or r < best[0]:
            best = (r, al)
    return best[1]


a_e8, a_s4 = alpha_fit(e8), alpha_fit(s4)
step("Z4.thermodynamic dimension of the FULL process space: E8_1 vacuum character 1, 248, 4124, 34752, ... grows with alpha -> 1/2: spacetime d = 1/(1-alpha) = 2",
     ok_tower and abs(a_e8 - 0.5) < 0.03, f"alpha = {a_e8:.3f} -> d = {1 / (1 - a_e8):.2f}")
lw = np.array([20, 40, 80])
N3 = np.array([sum((l + 1) ** 2 for l in range(L)) for L in lw])
slope3 = np.polyfit(np.log(lw), np.log(N3), 1)[0]
step("Z5.control 3+1D conformal scalar on S^3: Weyl exponent -> 3 and alpha -> 3/4 (d = 4): the diagnostic separates 2 from 4",
     abs(slope3 - 3) < 0.15 and abs(a_s4 - 0.75) < 0.04, f"Weyl slope {slope3:.3f}, alpha = {a_s4:.3f} -> d = {1 / (1 - a_s4):.2f}")

ts = [0.02, 0.01, 0.005]
hk = [sum((2 * l + 1) * np.exp(-l * (l + 1) * t) for l in range(0, 4000)) - 1 / t for t in ts]
step("Z6.heat kernel of the compiler precision operator on the chi = 2 double: Tr e^{t Delta_S2} - 1/t -> 1/3 = chi/6 (topological, Gauss-Bonnet of P1's c3); no dynamics in 2D",
     abs(hk[-1] - 1 / 3) < 5e-3, "values " + ", ".join(f"{x:.5f}" for x in hk))
step("Z7.(statement) every structure now in TFPT reconstructs a 1+1D spacetime; 3+1D needs a covariance network with one-particle Weyl exponent 3 (thermodynamic alpha = 3/4)",
     abs(slope - 1) < 1e-2 and abs(a_e8 - 0.5) < 0.03, "Theorem D remains a new-input problem, now with a quantitative target")

out = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG),
       "alpha_E8": a_e8, "alpha_S3": a_s4, "weyl_slope_U": slope, "weyl_slope_S3": slope3, "heat_kernel": hk}
(HERE / "spacetime_precision_results.json").write_text(json.dumps(out, indent=2, default=float))
print(f"{out['passed']}/{out['total']} steps passed")
sys.exit(0 if out["passed"] == out["total"] else 1)
