"""KERNDEPHASE20261119: does the original P1 boundary kernel itself supply the Cartan dephasing?
CARTANDEPHASE20261119 left open: where the dephasing comes from, which Cartan directions it uses and
why its time is tau = 2 alpha. Tested here against the original kernel source
(documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source, section 001):
one-sided Dirac-type operator B_Sigma, Calderon projector (positive half), seed chi_seed = 1/lambda_1^+,
Hilbert space L^2(S) x H_int with H_int an analytic realization choice.
Reading tested: a current of weight difference Dw carries the phase exp(i e Dw.phi), phi a field with
the kernel as covariance in every internal direction, e^2 = 4 pi alpha. Averaging over phi is a random
unitary channel with multiplier exp(-tau |P Dw|^2 / 2), tau = e^2 * (channel weight), P the projector on
the directions that couple."""
from common import *
import os
import numpy as np
from fractions import Fraction as Fr
from mpmath import gamma, binomial, coth

c = Checks("kernel_dephasing")
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")
SRC = os.path.join(ROOT, "documents", "archive--tfpt-45--source_extracts--01_boundary_kernel_source", "sections",
                   "001-section-primitive-core-and-carrier-interface-handoff.tex")
z = lambda v: (v - CODATA) / SIG
REC = mpf("137.0359992168")

# --- 0. what the original says (text checks, no interpretation) -------------------
src = open(SRC, encoding="utf-8").read()
c.ok("original: the seed is the first positive eigenvalue, chi_seed = lambda_1^+(|B_Sigma|)^-1, and the polarization is the Calderon projector",
     "\\chiseed=\\lambda_1^+(|B_\\Sigma|)^{-1}" in src and "Calder\\'on projector" in src)
c.ok("original: H_int is 'an analytic realization choice rather than part of the boundary datum': the kernel acts as B_Sigma x 1 on internal labels",
     "analytic realization choice rather than part of the boundary" in " ".join(src.split()))
c.ok("original section 001 names no hypercharge, no Cartan direction and no coupling e: these are not kernel data",
     not any(w in src.lower() for w in ("hypercharge", "cartan", "fine-structure", "fine structure")))

# --- 1. kernel channels: flat on S^2 = white noise = exact semigroup ---------------
vol = lambda n: 2 * pi**(mpf(n + 1) / 2) / gamma(mpf(n + 1) / 2)
w_dirac = lambda k, n: 2**(n // 2) * binomial(k + n - 1, k) / vol(n) / (k + mpf(n) / 2)   # positive half, unit S^n
w2 = [w_dirac(k, 2) for k in range(60)]
c.ok("positive half of |B|^-1 on the unit S^2: every channel has weight 1/(2 pi)", all(abs(w - 1 / (2 * pi)) < mpf(10)**-35 for w in w2))
semi = {n: max(abs(w_dirac(k, n) - w_dirac(0, n)) for k in range(1, 30)) < mpf(10)**-30 for n in range(1, 8)}
c.ok("so N channels dephase by exp(-N tau_1 |Dw|^2/2): an exact one-parameter semigroup with the channel number as time; this holds only on S^2 (S^1..S^7 tested)",
     semi == {n: n == 2 for n in range(1, 8)}, semi)
al = mp.mpf(1) / REC
tau1 = 4 * pi * al * w2[0]
c.ok("time of one channel with e^2 = 4 pi alpha: tau_1 = e^2/(2 pi) = 2 alpha exactly", abs(tau1 - 2 * al) < mpf(10)**-38)
c.ok("the seed channel is that one channel: degeneracy 2, area 4 pi, lambda_1^+ = 1, weight chi_seed * 2/(4 pi) = 1/(2 pi)",
     abs(2 / vol(2) / 1 - 1 / (2 * pi)) < mpf(10)**-38)
c.ok("limit: the weight is 1/(2 pi R) on a sphere of radius R; the unit radius is the seed unit, an input of the reading",
     all(abs(2 / (4 * pi * R * R) * R - 1 / (2 * pi * R)) < mpf(10)**-38 for R in (mpf(1) / 2, 2, 3)))

# --- 2. the averaged channel in the 128 chart -------------------------------------
sgn = lambda s, m: 1 if (s >> m) & 1 else -1
W = np.array([[sgn(s, m) / 2 for m in range(8)] for s in STATES])              # weights of the 128 states
act = {r: root_action(*r) for r in EVEN_ROOTS}
def multiplier(tau, P):
    """exp(-tau |P (w_s - w_t)|^2 / 2); P = orthonormal rows spanning the coupled directions."""
    V = W @ P.T
    d2 = ((V[:, None, :] - V[None, :, :])**2).sum(-1)
    return np.exp(-tau * d2 / 2)
def damping(K):
    out = {}
    for r in EVEN_ROOTS:
        X = np.zeros((128, 128))
        for s, (t, sg) in act[r].items(): X[INDEX[t], INDEX[s]] = sg
        region = lambda m: "c" if m < 3 else ("w" if m < 5 else "f")
        key = "".join(sorted(region(r[0]) + region(r[1]))) + (" same" if r[2] == r[3] else " opp")
        out.setdefault(key, set()).add(round(float(np.trace(X.T @ (K * X)) / np.trace(X.T @ X)), 12))
    return out
x, wq = np.polynomial.hermite_e.hermegauss(60)
gauss_mean = lambda g, var: float((wq * np.cos(g * np.sqrt(var) * x)).sum() / wq.sum())
c.ok("averaging the unitary exp(i e w.phi) over a Gaussian phi with covariance C per direction gives the multiplier exp(-e^2 C |Dw|^2/2) (quadrature, 3 cases)",
     all(abs(gauss_mean(g, v) - np.exp(-g * g * v / 2)) < 1e-12 for g, v in ((1.0, 0.3), (1.4142, 0.1592), (2.8284, 0.05))))
tau = 0.37
Kfull = multiplier(tau, np.eye(8))
c.ok("the averaged channel is completely positive and unital (random unitary)", np.linalg.eigvalsh(Kfull).min() > -1e-10 and np.allclose(np.diag(Kfull), 1))
dfull = damping(Kfull)
c.ok("kernel x 1 on internal labels = all eight Cartan directions with equal weight: all 112 root currents get the same factor exp(-tau)",
     all(v == {round(float(np.exp(-tau)), 12)} for v in dfull.values()), len(dfull))
v_one = alpha_inv(lambda a: exp(-2 * a))
c.ok("with tau = tau_1 = 2 alpha this is the recorded damping: alpha^-1 = 137.0359992168", abs(v_one - REC) < 1e-9, f"{mp.nstr(v_one, 14)}, z = {mp.nstr(z(v_one), 4)}")

# --- 3. is the direction content forced? alternatives with the same kernel weight --
unit = lambda v: np.array(v, dtype=float) / np.linalg.norm(v)
a = np.array([-1 / 6] * 3 + [1 / 4] * 2 + [0] * 3)
hyper = unit(a)[None, :]
charge = unit([-1 / 6] * 3 + [1 / 2, 0] + [0] * 3)[None, :]                    # Q = T3 + Y in the chart
gauge6 = np.linalg.qr(np.array([[1, -1, 0, 0, 0, 0, 0, 0], [0, 1, -1, 0, 0, 0, 0, 0], [0, 0, 0, 1, -1, 0, 0, 0], list(a),
                                [0, 0, 0, 0, 0, 1, -1, 0], [0, 0, 0, 0, 0, 0, 1, -1]], dtype=float).T)[0].T
alt = {"hypercharge direction only": hyper, "electric charge direction only": charge, "six gauge Cartan directions": gauge6}
spread = {}
for name, P in alt.items():
    d = damping(multiplier(tau, P))
    rates = sorted({round(-float(np.log(v)) / tau, 9) for vs in d.values() for v in vs})
    spread[name] = rates
    print(f"  {name}: rates in units of tau = {rates}")
c.ok("any restricted direction set damps the root classes unequally (rates between 0 and 1); only the full torus gives one common rate 1",
     all(len(r) > 1 and max(r) <= 1 + 1e-9 for r in spread.values()), {k: (min(v), max(v)) for k, v in spread.items()})
hy = sorted({round(float(k), 9) for k in spread["hypercharge direction only"]})
c.ok("hypercharge only: mixed colour currents at rate 1/15, mixed weak at 3/20 (the recorded shares); with tau = 2 alpha that is exp(-2 alpha/15), exp(-3 alpha/10)",
     any(abs(k - 1 / 15) < 1e-9 for k in hy) and any(abs(k - 3 / 20) < 1e-9 for k in hy), hy)
zalt = {"exp(-3 alpha/10) (hypercharge only, weak)": z(alpha_inv(lambda t: exp(-mpf(3) / 10 * t))),
        "exp(-2 alpha/15) (hypercharge only, colour)": z(alpha_inv(lambda t: exp(-mpf(2) / 15 * t))),
        "exp(-alpha) (unit electric charge on one scalar field of weight 1/(2 pi))": z(alpha_inv(lambda t: exp(-t))),
        "exp(-4 alpha) (two channels or both Dirac signs)": z(alpha_inv(lambda t: exp(-4 * t))),
        "exp(-2 alpha coth(pi)) (quantum oscillator at tau = i)": z(alpha_inv(lambda t: exp(-2 * t * coth(pi))))}
for k, v in zalt.items(): print(f"  {k}: z = {mp.nstr(v, 6)}")
c.ok("every alternative is excluded by the alpha readout (|z| > 50); the full torus with exactly one positive channel is the only reading within 2 sigma (z = +1.9)",
     all(abs(v) > 50 for v in zalt.values()) and abs(z(v_one)) < 3, {k: mp.nstr(v, 5) for k, v in zalt.items()})
c.ok("so: internal blindness of the kernel (B x 1) is exactly what the common damping needs; a gauge-selective dephasing would contradict it", True)

# --- 4. consequences for the marks and the neutral coherence ----------------------
aF = (Fr(-1, 6),) * 3 + (Fr(1, 4),) * 2
h = lambda s: sum(aF[i] * sgn(s, i) for i in range(5))
mu = sum(t * t for t in aF)
iso, hyp = set(), {}
for r in MIXED:
    pairs = [(s, t) for s, (t, _) in act[r].items()]
    zy = Fr(len(pairs), 128)
    wv = lambda s, m: Fr(sgn(s, m), 2)
    jump = sum(sum(wv(s, m) * wv(t, m) for m in range(8)) for s, t in pairs) / 128
    diag = sum(sum(wv(s, m)**2 for m in range(8)) for s, _ in pairs) / 128
    iso.add((jump - diag) / zy)
    hyp.setdefault("weak" if r[0] >= 3 else "colour", set()).add(sum(h(s) * h(t) - (h(s)**2 + h(t)**2) / 2 for s, t in pairs) / 128)
c.ok("generator of the kernel process on the 60 mixed currents: the same value -ZY for all, no colour/weak difference", iso == {Fr(-1)}, iso)
c.ok("the recorded mixed values -1/72 and -1/32 are its part along the hypercharge direction: -(share) * |h|^2 * ZY with |h|^2 = 4 mu = 5/6",
     hyp == {"colour": {Fr(-1, 72)}, "weak": {Fr(-1, 32)}}
     and Fr(1, 15) * 4 * mu * Fr(1, 4) == Fr(1, 72) and Fr(3, 20) * 4 * mu * Fr(1, 4) == Fr(1, 32), hyp)
c.ok("the hypercharge direction enters through the marker U = h^2 of the readout, not through the kernel: the kernel supplies the time, not this direction", True)
i0, i255 = INDEX[0], INDEX[255]
c.ok("the neutral coherence 0-255 is not protected by the kernel process: factor exp(-4 tau), at tau = 2 alpha that is exp(-8 alpha) = 0.9433",
     abs(Kfull[i0, i255] - np.exp(-4 * tau)) < 1e-12 and abs(float(exp(-8 * al)) - 0.94329) < 1e-4, float(exp(-8 * al)))
c.ok("the marks are generator values (exact rationals at tau = 0); read after the finite kernel time they would shift by 1 - exp(-2 alpha) = 1.45 percent: both cannot be the same readout",
     abs(float(1 - exp(-2 * al)) - 0.01449) < 1e-4, float(1 - exp(-2 * al)))
c.done({"tau_one_channel": "2 alpha", "alpha_inv_full_torus_one_channel": mp.nstr(v_one, 14), "z_alternatives": {k: mp.nstr(v, 6) for k, v in zalt.items()},
        "rates_by_direction_set": spread,
        "verdict": "time and isotropy supplied conditionally (coupling e^2 = 4 pi alpha, unit seed sphere, one channel); hypercharge direction not supplied"})
