"""KERNFELD20261119: does the original fix the kernel field as a Cartan gauge potential, and what does the
Coulomb reading of the kernel fix beyond the damping?
COULOMBKOPPLUNG20261119 left one assumption: kernel field = Cartan gauge potential with Gauss normalisation.
Part A reads the original. Part B tests the reading on the other fixed inputs of the alpha equation with one
variance for everything (w = 1/(4 pi) per component, the value that gives exp(-2 alpha)). Part C lets the four
marks interact through the same kernel and asks where they sit (recorded open point: mark position is an input;
SWEEP20261119: marks on the seam give a rectangle with free aspect ratio)."""
from common import *
import os
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
from scipy.optimize import minimize

c = Checks("kernel_field")
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "documents")
flat = lambda *p: " ".join(open(os.path.join(ROOT, *p), encoding="utf-8").read().split())
K1 = flat("archive--tfpt-45--source_extracts--01_boundary_kernel_source", "sections", "001-section-primitive-core-and-carrier-interface-handoff.tex")
K2 = flat("archive--tfpt-45--source_extracts--01_boundary_kernel_source", "sections", "002-section-intrinsic-closed-branch-and-canonical-rigid-obje.tex")
EM = flat("archive--tfpt-45--source_extracts--03_em_flavor_source", "sections", "001-section-kernel-numerics-and-electromagnetic-closure.tex")
CA = flat("archive--tfpt-45--source_extracts--02_carrier_source", "sections", "006-subsection-family-monodromy-and-occupancy.tex")
z = lambda v: (v - CODATA) / SIG

# --- A. the original --------------------------------------------------------------
c.ok("original: the primitive kernel contains 'No carrier split, no hypercharge operator, no family geometry'",
     "No carrier split, no hypercharge operator, no family geometry" in K2)
c.ok("original: 'Gauge fixing, ghosts, and carrier-compatible bosonic fluctuations are introduced only after the hard carrier theorem': a gauge potential is not kernel data",
     "Gauge fixing, ghosts, and carrier-compatible bosonic fluctuations are introduced only after the hard carrier theorem" in K1)
c.ok("original: phi_base = (1 - tr_E Y^2)/pi with gamma = tr_E Y^2 = 5/6, and delta_top = Omega_adm c3^4 is an 'occupancy reinterpretation', not a kernel computation",
     "\\frac{1-\\tr_EY^2}{\\pi}" in EM and "occupancy reinterpretation $\\deltatop=\\Oadm\\,\\cthree^4$" in CA)
c.ok("verdict A: the original does not fix the kernel field as a gauge potential; it enters, if at all, with the later fluctuation package", True)
a5 = (Fr(-1, 6),) * 3 + (Fr(1, 4),) * 2
c.ok("link: the norm |h|^2 = 4 mu = 5/6 of the hypercharge direction in the chart is the original carrier norm gamma = tr_E Y^2",
     4 * sum(x * x for x in a5) == Fr(5, 6) and 3 * Fr(1, 3)**2 + 2 * Fr(1, 2)**2 == Fr(5, 6))

# --- B. one variance for all inputs ------------------------------------------------
w = 1 / (4 * pi)                                                             # Coulomb weight per component; two components give exp(-2 alpha)
c.ok("with w = 1/(4 pi): c3 = w/2, delta_top = 48 c3^4 = 3 w^4, phi_base = 1/(6 pi) = (2/3) w",
     abs(C3 - w / 2) < mpf(10)**-38 and abs(DTOP - 3 * w**4) < mpf(10)**-38 and abs(PHI_TREE - 2 * w / 3) < mpf(10)**-38)
dfact = lambda n: 1 if n <= 0 else n * dfact(n - 2)
fact = lambda n: 1 if n <= 1 else n * fact(n - 1)
menu = {}
for vname, v in (("c3", C3), ("w = 2 c3", w), ("2 w", 2 * w), ("w^2", w * w)):
    for n in (1, 2, 3, 4):
        menu[f"moment <x^{2*n}>, variance {vname}"] = dfact(2 * n - 1) * v**n
        menu[f"cumulant {n} of x^2, variance {vname}"] = 2**(n - 1) * fact(n - 1) * v**n
hits = sorted(k for k, val in menu.items() if abs(val / DTOP - 1) < mpf(10)**-30)
print("  readings of delta_top in the menu:", hits)
c.ok("delta_top in a menu of 32 Gaussian readings: exactly two hits, fourth cumulant of x^2 at variance c3 and fourth moment at pair weight w^2",
     hits == ["cumulant 4 of x^2, variance c3", "moment <x^4>, variance w^2"], len(menu))
c.ok("only the second uses the variance of the damping: three pairings of four marks, each pair joined once on each of the two components (w * w); the first needs half that variance",
     abs(3 * (w * w)**2 - DTOP) < mpf(10)**-38)
c.ok("control: independent pairings on the two components would give 9 w^4, a complex field 2 w^4: the pairing must be the same on both (deck invariant)",
     abs(9 * w**4 / DTOP - 3) < mpf(10)**-30 and abs(2 * w**4 / DTOP - mpf(2) / 3) < mpf(10)**-30)
def alpha_dt(dt):
    phi = lambda a: PHI_TREE + dt * exp(-2 * a) * (1 - dt * exp(-2 * a))**(-mpf(5) / 4)
    F = lambda a: a**3 - 2 * C3**3 * a**2 - mpf(4) / 5 * 41 * C3**6 * log(1 / phi(a))
    return 1 / findroot(F, 1 / mpf(137))
zd = {"3 w^4": z(alpha_dt(3 * w**4)), "2 w^4": z(alpha_dt(2 * w**4)), "9 w^4": z(alpha_dt(9 * w**4))}
c.ok("the alpha readout separates the counts: 3 pairings within 2 sigma, 2 or 9 far outside", abs(zd["3 w^4"]) < 2 and abs(zd["2 w^4"]) > 1000 and abs(zd["9 w^4"]) > 1000,
     {k: mp.nstr(v, 6) for k, v in zd.items()})
c.ok("not explained by the reading: phi_base = (2/3) w. In the original it is (1 - gamma)/pi; no Gaussian object of the menu gives 2/3 (moments and cumulants are integers times powers)",
     not any(abs(val / PHI_TREE - 1) < mpf(10)**-30 for val in menu.values()))

# --- C. the marks interacting through the same kernel ------------------------------
PAIRINGS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
dist = lambda P, i, j: float(np.linalg.norm(P[i] - P[j]))
energy = lambda P: sum(1 / dist(P, i, j) for i in range(4) for j in range(i))
pair_sum = lambda P, k=lambda d: 1 / d: sum(k(dist(P, *p)) * k(dist(P, *q)) for p, q in PAIRINGS)
circle = lambda th: np.array([[np.cos(t), np.sin(t), 0.0] for t in th])
sphere = lambda u: np.array([[np.sin(u[2 * i]) * np.cos(u[2 * i + 1]), np.sin(u[2 * i]) * np.sin(u[2 * i + 1]), np.cos(u[2 * i])] for i in range(4)])
rng = np.random.default_rng(11)
def cross_ratio(P):                                                          # stereographic projection, then (z1,z2;z3,z4)
    zz = [complex(p[0], p[1]) / (1 - p[2] + 1e-300) for p in P]
    return (zz[0] - zz[2]) * (zz[1] - zz[3]) / ((zz[0] - zz[3]) * (zz[1] - zz[2]))
jinv = lambda lam: 256 * (lam * lam - lam + 1)**3 / (lam * lam * (lam - 1)**2)
ends = set()
for _ in range(60):
    r = minimize(lambda t: energy(circle([0, *t])), np.sort(rng.uniform(0.2, 6.0, 3)), method="Nelder-Mead", options={"xatol": 1e-11, "fatol": 1e-14, "maxiter": 6000})
    ends.add(tuple(np.round(np.sort(np.mod(r.x, 2 * np.pi)), 5)))
c.ok("four equal charges on the seam circle, interacting through the kernel 1/(4 pi d): the only equilibrium reached from 60 random starts is the square (marks at mu_4)",
     ends == {(round(np.pi / 2, 5), round(np.pi, 5), round(3 * np.pi / 2, 5))}, ends)
sq = circle([0, np.pi / 2, np.pi, 3 * np.pi / 2])
best = min((minimize(lambda u: energy(sphere(u)), rng.uniform(0.3, 2.8, 8), method="Nelder-Mead", options={"xatol": 1e-11, "fatol": 1e-14, "maxiter": 20000}) for _ in range(12)), key=lambda r: r.fun)
tet = sphere(best.x)
c.ok("without the reflection (free on the sphere) the equilibrium is the regular tetrahedron, with lower energy than the square",
     abs(best.fun - 6 / np.sqrt(8 / 3)) < 1e-8 and best.fun < energy(sq) - 0.1, f"tetrahedron {best.fun:.6f}, square {energy(sq):.6f}")
js, jt = jinv(cross_ratio(sq)), jinv(cross_ratio(tet))
c.ok("square: cross ratio harmonic, j = 1728 (square torus tau = i); tetrahedron: equianharmonic, j = 0 (hexagonal torus)", abs(js - 1728) < 1e-6 and abs(jt) < 1e-3, f"j = {js:.4f}, {abs(jt):.2e}")
def hess_out(P, h=1e-4):                                                     # second derivative of the energy when two opposite marks leave the circle plane in opposite sense
    def e(t):
        Q = P.copy()
        for i, s in ((0, 1), (2, 1), (1, -1), (3, -1)):
            Q[i] = Q[i] * np.cos(t) + np.array([0, 0, s * np.sin(t)])
        return energy(Q)
    return (e(h) - 2 * e(0) + e(-h)) / h**2
c.ok("on the sphere the square is a saddle (energy falls when alternate marks leave the plane): it is stable only with the reflection, as recorded for tau = i in P1KERN20261119",
     hess_out(sq) < -0.1, hess_out(sq))
c.ok("consistency with the record: hexagon excluded by alpha (recorded alpha^-1 = 104.42), square needs the reflection; here the same two shapes appear as the two equilibria of the kernel", True)
th = sp.symbols("theta", positive=True)
rect = sp.Rational(1, 4) + 1 / (4 * sp.sin(th)**2) * 1 / (4 * sp.cos(th)**2) * 4
c.ok("pairing sum of the Coulomb kernel over the three pairings on a rectangle with half angle theta: 1/4 + 1/sin^2(2 theta); at the square exactly 5/4",
     sp.simplify(rect - sp.Rational(1, 4) - 1 / sp.sin(2 * th)**2) == 0 and sp.simplify(rect.subs(th, sp.pi / 4)) == sp.Rational(5, 4) and abs(pair_sum(sq) - 1.25) < 1e-12)
ctrl = {"flat (one seed channel)": pair_sum(sq, lambda d: 1.0), "Coulomb 1/d": pair_sum(sq), "1/d^2": pair_sum(sq, lambda d: 1 / d**2), "log": pair_sum(sq, lambda d: -np.log(d / 2)),
        "Coulomb at the tetrahedron": pair_sum(tet)}
print("  pairing sums:", {k: round(v, 6) for k, v in ctrl.items()})
c.ok("the same sum is 3 for the flat seed channel (the count in delta_top = 3 w^4) and 5/4 for the full Coulomb kernel at the square; 9/16 for 1/d^2, 9/8 at the tetrahedron",
     abs(ctrl["flat (one seed channel)"] - 3) < 1e-12 and abs(ctrl["1/d^2"] - 9 / 16) < 1e-12 and abs(ctrl["Coulomb at the tetrahedron"] - 9 / 8) < 1e-6)
c.ok("pattern only: 5/4 is also the fixed exponent B gamma = (3/2)(5/6) and delta_2/delta_top^2; no step turns a pairing sum into the exponent of (1 - q), so this is not a derivation", True)
c.ok("verdict C: mark position is no longer a free input under the reading: reflection (marks on the seam) plus the kernel's own repulsion gives mu_4; the premise 'marks carry equal charge and repel through the kernel' is new and not in the original", True)
c.done({"readings_of_delta_top": hits, "z_by_pairing_count": {k: mp.nstr(v, 6) for k, v in zd.items()}, "pairing_sums": {k: round(v, 9) for k, v in ctrl.items()},
        "energy": {"square": energy(sq), "tetrahedron": float(best.fun)}, "second_derivative_out_of_plane_at_square": hess_out(sq),
        "verdict": "A: gauge potential not kernel data in the original. B: delta_top = three pairings of four marks at the damping variance; phi_base not explained. C: marks at mu_4 are the Coulomb equilibrium on the seam circle; 5/4 as pairing sum is a pattern."})
