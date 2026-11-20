"""DEBYEWALLER20261119 (idea 3): is exp(-2 alpha) the charged average of a root
insertion in the recorded self-dual E8 state at tau = i?
Declared before computing: phase exp(i e beta.p) with e^2 = 4 pi alpha
(Heaviside-Lorentz), beta an E8 root, two sheets (recorded cover degree d = 2).
Exact object: phi(z) = E_{4,1}(i,z)/E_4(i), z = e/(2 pi)."""
from common import *
from mpmath import jtheta, sqrt, taylor
c = Checks("debye_waller")
q = exp(-pi)                                           # nome at tau = i
def theta_e8(zz):                                      # root beta = e1+e2 in D8 coordinates
    t = lambda n, w: jtheta(n, pi * w, q)
    return (t(3, zz)**2 * t(3, 0)**6 + t(4, zz)**2 * t(4, 0)**6 + t(2, zz)**2 * t(2, 0)**6) / 2
phi = lambda zz: theta_e8(zz) / theta_e8(0)
E4i = theta_e8(0)
c.ok("lattice sum at tau=i equals E4(i) = 3 Gamma(1/4)^8/(2 pi)^6", abs(E4i - 3 * mp.gamma(mpf(1) / 4)**8 / (2 * pi)**6) < 1e-30, mp.nstr(E4i, 15))
tay = taylor(lambda w: log(phi(w)), 0, 8)
c.ok("second moment: log phi = -pi z^2 + ..., i.e. <p_i p_j> = delta_ij/(4 pi) = 2 c3 exactly", abs(tay[2] + pi) < 1e-25, mp.nstr(tay[2], 20))
c.ok("self-duality forces this: phi(z) e^{pi z^2} is invariant under z -> i z (no z^2, z^6 terms)", abs(tay[6]) < 1e-18 and abs(phi(mpf("0.3") * 1j) * exp(-pi * mpf("0.09")) - phi(mpf("0.3")) * exp(pi * mpf("0.09"))) < 1e-25)
a4 = tay[4]
print("  a4 =", mp.nstr(a4, 15), "  a8 =", mp.nstr(tay[8], 8))
base = lambda a: PHI_TREE + (DTOP * exp(-2 * a)) * (1 - DTOP * exp(-2 * a))**(-mpf(5) / 4)
def exact(a, sheets=2):
    x = DTOP * phi(sqrt(4 * pi * a) / (2 * pi))**sheets
    return PHI_TREE + x * (1 - x)**(-mpf(5) / 4)
a0 = 1 / alpha_inv(base)
c.ok("Gaussian part reproduces the recorded damping: phi(e/2pi)^2 = exp(-2 alpha)(1 + O(alpha^2))", abs(phi(sqrt(4 * pi * a0) / (2 * pi))**2 / exp(-2 * a0) - 1) < 1e-3, mp.nstr(phi(sqrt(4 * pi * a0) / (2 * pi))**2 / exp(-2 * a0) - 1, 6))
c.ok("one sheet alone gives exp(-alpha) to leading order: the factor 2 needs the cover degree, not the root norm alone", abs(log(phi(sqrt(4 * pi * a0) / (2 * pi))) / a0 + 1) < 1e-2)
zs = lambda f: sum(f(n) * exp(-pi * n * n) for n in range(-40, 41))
c.ok("not E8 selective: the plain Z lattice at tau=i has the same second moment 1/(4 pi)", abs(zs(lambda n: n * n) / zs(lambda n: 1) * 4 * pi - 1) < 1e-30)
v0, v1 = alpha_inv(base), alpha_inv(exact)
print(f"  recorded closure   {mp.nstr(v0, 15)}  z = {mp.nstr(z(v0), 5)}   (Rb: {mp.nstr(z(v0, RB, RB_SIG), 4)})")
print(f"  exact lattice mean {mp.nstr(v1, 15)}  z = {mp.nstr(z(v1), 5)}   (Rb: {mp.nstr(z(v1, RB, RB_SIG), 4)})")
shift = (v1 - v0) / SIG
c.ok("the non-Gaussian correction is parameter free and resolvable (|shift| > 0.5 sigma)", abs(shift) > 0.5, mp.nstr(shift, 5))
verdict = "SUPPORTED" if abs(z(v1)) < 2 else ("EXCLUDED at >5 sigma" if abs(z(v1)) > 5 else "TENSION 2-5 sigma")
print("  verdict for the strict lattice-average reading:", verdict)
c.ok("verdict recorded (strict reading tested against CODATA 2022)", True, verdict)
c.done({"a4": mp.nstr(a4, 15), "alpha_inv_recorded": mp.nstr(v0, 16), "alpha_inv_exact_lattice": mp.nstr(v1, 16),
        "z_recorded": mp.nstr(z(v0), 6), "z_exact_lattice": mp.nstr(z(v1), 6), "shift_sigma": mp.nstr(shift, 6), "verdict": verdict})
