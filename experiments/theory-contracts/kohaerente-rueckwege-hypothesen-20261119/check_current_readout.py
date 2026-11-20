"""EMCURRENT20261119 (hypothesis 2): is exp(-2 alpha) a protected Heisenberg current readout?
Declared before computing, as in DEBYEWALLER20261119: phase exp(i e beta.J), e^2 = 4 pi alpha,
beta an E8 root (beta^2 = 2). Lattice vertex algebra V_E8 = sum over lattice charges of
Heisenberg Fock modules; a smeared current J(f) = f_0 a_0 + sum_n f_n a_n lies in the
Heisenberg algebra. Three states are compared: the vacuum, the torus state at tau = i with
mean-free f (oscillators only), and the pure zero mode at tau = i (the recorded exclusion)."""
from common import *
import sympy as sp
from mpmath import jtheta, taylor, coth, sqrt

c = Checks("current_readout")
# --- oscillator part: exact Gaussianity (one mode, exact rational matrices) ---
n = 16
adag = sp.zeros(n)
for k in range(n - 1): adag[k + 1, k] = sp.sqrt(k + 1)
aop = adag.T
X = aop + adag
vac = sp.zeros(n, 1); vac[0] = 1
mom = [sp.simplify((vac.T * X**k * vac)[0]) for k in range(9)]
k2 = mom[2]; k4 = mom[4] - 3 * mom[2]**2
k6 = mom[6] - 15 * mom[4] * mom[2] + 30 * mom[2]**3
k8 = mom[8] - 28 * mom[6] * mom[2] - 35 * mom[4]**2 + 420 * mom[4] * mom[2]**2 - 630 * mom[2]**4
c.ok("vacuum, one oscillator mode: cumulants kappa_4 = kappa_6 = kappa_8 = 0 exactly", (k4, k6, k8) == (0, 0, 0), mom[:9])
c.ok("zero mode in the vacuum: the vacuum has lattice charge 0, so exp(i t a_0) acts as 1 (no lattice sum at all)", True,
     "vacuum vector lies in the charge-0 Fock module")
# --- thermal oscillator: still Gaussian, variance multiplied by coth -----------
q = sp.Rational(1, 7)                                   # any Boltzmann factor: exact check of Gaussianity
Z = sum(q**k for k in range(n))
th = lambda k: sum(q**m * (X**k)[m, m] for m in range(n - 6)) / sum(q**m for m in range(n - 6))
t2, t4 = sp.nsimplify(th(2)), sp.nsimplify(th(4))
c.ok("thermal oscillator stays Gaussian (kappa_4 -> 0 with the cutoff), variance (1 + q)/(1 - q) times the vacuum value",
     abs(float(t4 - 3 * t2**2)) < 1e-4 and abs(float(t2 - (1 + q) / (1 - q))) < 1e-5, (float(t2), float(t4 - 3 * t2**2)))

# --- zero mode at tau = i: the recorded non-Gaussian lattice mean -------------
qq = exp(-pi)
def theta_e8(zz):
    t = lambda m, w: jtheta(m, pi * w, qq)
    return (t(3, zz)**2 * t(3, 0)**6 + t(4, zz)**2 * t(4, 0)**6 + t(2, zz)**2 * t(2, 0)**6) / 2
phi = lambda zz: theta_e8(zz) / theta_e8(0)
tay = taylor(lambda w: log(phi(w)), 0, 4)
c.ok("pure zero mode at tau = i: log phi = -pi z^2 + a4 z^4 with a4 = 7.878 (recorded value reproduced)",
     abs(tay[2] + pi) < 1e-25 and abs(tay[4] - mpf("7.878")) < 1e-3, mp.nstr(tay[4], 10))

# --- consequences for the alpha root ------------------------------------------
gauss = lambda a: exp(-2 * a)
v_rec = alpha_inv(gauss)
zs = lambda v: (v - CODATA) / SIG
c.ok("recorded closure reproduced", abs(v_rec - mpf("137.0359992168")) < 1e-9, mp.nstr(v_rec, 14))
v_lat = alpha_inv(lambda a: phi(sqrt(4 * pi * a) / (2 * pi))**2)
c.ok("zero-mode reading at tau = i: excluded (recorded +141.7 sigma)", abs(zs(v_lat) - mpf("141.7")) < 0.5, mp.nstr(zs(v_lat), 6))
v_th = alpha_inv(lambda a: exp(-2 * a * coth(pi)))      # oscillator at the same tau = i: beta*omega = 2 pi
c.ok("oscillator-only reading at tau = i is Gaussian but its variance carries coth(pi): excluded as well",
     abs(zs(v_th)) > 5, f"alpha^-1 = {mp.nstr(v_th, 14)}, {mp.nstr(zs(v_th), 6)} sigma, coth(pi) - 1 = {mp.nstr(coth(pi) - 1, 5)}")
# normalisation: the recorded value needs <x^2> = <p^2> = 1/(4 pi) (h = 1); canonical [a, a^+] = 1 gives 1/2
v_hbar = alpha_inv(lambda a: exp(-4 * pi * a))
c.ok("with the canonical mode normalisation (hbar = 1) the same vacuum readout gives exp(-4 pi alpha): far off",
     abs(zs(v_hbar)) > 1000, f"alpha^-1 = {mp.nstr(v_hbar, 12)}, {mp.nstr(zs(v_hbar), 6)} sigma")
# the variance is a free number of the test function: kappa_2 = beta^2 * sum_n n |f_n|^2
need = 1 / (2 * pi)
c.ok("required smearing norm: sum_n n |f_n|^2 = 1/(2 pi) per unit beta^2/2; the current algebra does not supply this number", True,
     mp.nstr(need, 12))
c.done({"a4": mp.nstr(tay[4], 12), "alpha_inv_vacuum_gaussian": mp.nstr(v_rec, 15),
        "sigma_zero_mode_tau_i": mp.nstr(zs(v_lat), 6), "sigma_oscillator_tau_i": mp.nstr(zs(v_th), 6),
        "sigma_hbar_normalisation": mp.nstr(zs(v_hbar), 6), "required_norm": mp.nstr(need, 15)})
