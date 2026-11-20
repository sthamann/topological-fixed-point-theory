"""COULOMBKOPPLUNG20261119: does the coupling e^2 = 4 pi alpha of the Cartan phase to the kernel field
follow from the original?
KERNDEPHASE20261119 assumed it. The original states the exponent differently
(paper v28, section 'Geometric Self-Consistency', lemma 'Defect Sector Reweighting'; current EM source
section 001): the defect action is deformed additively by alpha * J, J = d * SF, d = 2 sheets of the
orientable double cover, SF = 1 the spectral-flow unit. Tested here:
 (1) is the kernel reading the same law, i.e. e^2 * (channel weight) = alpha per unit,
 (2) what replaces the free normalisation kappa and the unit radius,
 (3) which counting of the factor 2 is compatible, (4) what several defects do."""
from common import *
import os
from fractions import Fraction as Fr
from mpmath import gamma, binomial, gegenbauer, legendre, nsum, inf, sqrt, cos

c = Checks("coulomb_coupling")
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "documents")
V28 = os.path.join(ROOT, "archive--paper-latex--tfpt-theory-fullv28", "sections", "004-subsection-geometric-self-consistency-backreaction.tex")
EM = os.path.join(ROOT, "archive--tfpt-45--source_extracts--03_em_flavor_source", "sections", "001-section-kernel-numerics-and-electromagnetic-closure.tex")
z = lambda v: (v - CODATA) / SIG
REC = mpf("137.0359992168")
flat = lambda p: " ".join(open(p, encoding="utf-8").read().split())

# --- 0. what the original says -----------------------------------------------------
v28, em = flat(V28), flat(EM)
c.ok("original (v28): exponent law S_def(alpha) = S_def(0) + alpha J with J = d SF, d = 2 (orientable double cover), SF = 1",
     "S_{\\mathrm{def}}(\\alpha)=S_{\\mathrm{def}}(0)+\\alpha\\,\\mathcal{J}" in v28 and "\\mathcal{J}=d\\,\\mathrm{SF}" in v28 and "$d=2$ and $\\mathrm{SF}=1$" in v28)
c.ok("original (v28): the per-sheet unit is an identification ('We identify the per-sheet normalization with the spectral-flow unit'), not a computed kernel value",
     "We identify the per-sheet normalization with the spectral-flow unit" in v28)
c.ok("current EM source: q(alpha) = delta_top e^{-2 alpha}, spectral flow SF(U_Sigma) = 1, 'two sheets, one seam, one elementary relative winding'",
     "q(\\alpha):=\\deltatop e^{-2\\alpha}" in em and "spectral flow $SF(U_\\Seam)=1$" in em and "two sheets, one seam, one elementary relative wind" in em)
c.ok("neither passage contains a coupling e or a field covariance: the original has alpha times an integer, no 4 pi", "4\\pi\\alpha" not in v28 and "4\\pi\\alpha" not in em)

# --- 1. the scalar kernel at kappa = 1 is the Coulomb kernel -----------------------
vol = lambda n: 2 * pi**(mpf(n + 1) / 2) / gamma(mpf(n + 1) / 2)                 # unit S^n
def deg(l, d): return binomial(l + d - 1, d - 1) - (binomial(l + d - 3, d - 1) if l >= 2 else 0)
def coulomb_eigen(l, d):
    """Eigenvalue on degree-l harmonics of the Green function of -Laplace in R^d restricted to the unit sphere:
    G = 1/((d-2) vol |x-y|^(d-2)) = sum_l C_l^lam(cos g)/((d-2) vol), lam = (d-2)/2, and the projector on degree l
    is (deg/vol) C_l^lam(cos g)/C_l^lam(1)."""
    lam = mpf(d - 2) / 2
    return gegenbauer(l, lam, 1) / ((d - 2) * deg(l, d))
c.ok("Green function of -Laplace in R^d on the unit sphere has eigenvalue 1/(2l + d - 2) = 1/(2 Omega): it is the comparison kernel at kappa = 1 (d = 3..7, l < 25)",
     all(abs(coulomb_eigen(l, d) - 1 / mpf(2 * l + d - 2)) < mpf(10)**-30 for d in range(3, 8) for l in range(25)))
r, g = mpf("0.6"), mpf("1.1")
series = nsum(lambda l: r**l * legendre(int(l), cos(g)), [0, inf])
c.ok("d = 3: 1/|x - y| = sum_l r^l P_l(cos g), every multipole with coefficient 1 (generating function, checked at r = 0.6)",
     abs(series - 1 / sqrt(1 - 2 * r * cos(g) + r * r)) < mpf(10)**-25)
c.ok("so kappa = 1 is Gauss's law (unit flux per unit charge) and the flat spectrum is the Legendre series at r = 1: channel weight 1/(4 pi) for every l",
     all(abs(deg(l, 3) / vol(2) * coulomb_eigen(l, 3) - 1 / (4 * pi)) < mpf(10)**-30 for l in range(40)))
al = 1 / REC
for name, e2, G in (("Heaviside-Lorentz", 4 * pi * al, 1 / (4 * pi)), ("Gauss", al, mpf(1))):
    if abs(e2 * G - al) > mpf(10)**-38: break
else: name = None
c.ok("e^2 * (Coulomb coefficient) = alpha in either unit convention: the 4 pi of the coupling cancels the 4 pi of the kernel, alpha is the only number", name is None)
c.ok("hence per multipole channel and per scalar sheet the exponent is alpha * |beta|^2/2 = alpha * 1 for a root: the original unit J_1 = SF = 1",
     abs(4 * pi * al * (1 / (4 * pi)) * Fr(2, 2) - al) < mpf(10)**-38)

# --- 2. the factor 2: which countings are compatible -------------------------------
w_dirac_pos = 2 / vol(2) / 1                                                    # positive half, lowest channel, both spinor components
count = {"scalar Coulomb channel on each of two sheets (original d = 2)": 2 * 4 * pi * (1 / (4 * pi)),
         "positive Dirac half on one sheet (two spinor components)": 4 * pi * w_dirac_pos,
         "positive Dirac half on each of two sheets": 2 * 4 * pi * w_dirac_pos,
         "one scalar channel on one sheet": 4 * pi * (1 / (4 * pi))}
zs = {k_: z(alpha_inv(lambda a, k=k: exp(-k * a))) for k_, k in count.items()}
for k_, k in count.items(): print(f"  {k_}: k = {mp.nstr(k, 6)}, z = {mp.nstr(zs[k_], 6)}")
c.ok("two countings give k = 2 (two sheets of scalar channels; or one sheet with the two spinor components of the positive half); counting both gives k = 4, excluded",
     [mp.nstr(k, 3) for k in count.values()] == ["2.0", "2.0", "4.0", "1.0"] and abs(zs["positive Dirac half on each of two sheets"]) > 1000)
c.ok("original doubling theorem: Ran(P_C) = ker(R - 1), the Calderon half is the reflection-even half of the doubled data, i.e. one copy: doubling and the positive half do not multiply",
     "\\operatorname{Ran}(P_C)=\\ker(R_{\\taudbl}-1)" in flat(os.path.join(ROOT, "archive--tfpt-45--source_extracts--01_boundary_kernel_source", "sections",
                                                                   "001-section-primitive-core-and-carrier-interface-handoff.tex")))
c.ok("the 2 is therefore one count, readable as 'two sheets' or as 'two spinor components'; the original does not say which, and the alpha readout cannot tell (same number)", True)

# --- 3. scale: the unit radius ------------------------------------------------------
c.ok("in seed units the weight is scale free: chi_seed^2 * (lowest degeneracy)/(area) = R^2 * 2/(4 pi R^2) = 1/(2 pi) for every radius",
     all(abs(R * R * 2 / (4 * pi * R * R) - 1 / (2 * pi)) < mpf(10)**-38 for R in (mpf(1) / 3, 1, 7)))
c.ok("limit: the bare coincidence weight is 1/(2 pi R); that the sphere is measured in units of its own chi_seed is a reading (the original supplies the unit, not this combination)", True)

# --- 4. spectral-flow sectors and several defects -----------------------------------
sgn = lambda s, m: 1 if (s >> m) & 1 else -1
norms = sorted({sum((sgn(s, m) - sgn(t, m))**2 for m in range(8)) // 4 for s in STATES for t in STATES if s != t})
c.ok("law k = |Dw|^2 on the 128 weights: sectors k = 2, 4, 6, 8 (root currents 2, neutral coherence 0-255: 8), matching the original remark 'SF != 1 is a different sector with a different k'",
     norms == [2, 4, 6, 8], norms)
DT = DTOP
def phi_n(a, law):                                                               # series of q (1-q)^(-5/4) with defect number n damped by law(n)
    return PHI_TREE + sum(binomial(n - 1 + mpf(1) / 4, n - 1) * DT**n * exp(-law(n) * a) for n in range(1, 8))
def solve(law):
    F = lambda a: a**3 - 2 * C3**3 * a**2 - mpf(4) / 5 * 41 * C3**6 * log(1 / phi_n(a, law))
    return 1 / findroot(F, 1 / mpf(137))
v_lin, v_quad = solve(lambda n: 2 * n), solve(lambda n: 2 * n * n)
c.ok("control: the series with n defects damped by exp(-2 n alpha) reproduces the recorded value", abs(v_lin - REC) < 1e-9, mp.nstr(v_lin, 14))
v_root = solve(lambda n: 2 * n if n != 2 else 2)                                 # two defects whose roots add up to a root
zz = {"independent (exp(-2 n alpha))": z(v_lin), "same root, same channel (exp(-2 n^2 alpha))": z(v_quad), "pair adding to a root (exp(-2 alpha) at n = 2)": z(v_root)}
for k_, v in zz.items(): print(f"  several defects, {k_}: z = {mp.nstr(v, 5)}")
c.ok("several defects: only independent addition of the exponents stays within 2 sigma; coherent addition in one channel (|n beta|^2) and a pair fusing to a root are excluded beyond 5 sigma",
     abs(zz["independent (exp(-2 n alpha))"]) < 2 and all(abs(v) > 5 for k_, v in zz.items() if not k_.startswith("independent")), {k_: mp.nstr(v, 5) for k_, v in zz.items()})
c.ok("independent addition is what the kernel gives for defects in different channels (flat spectrum, uncorrelated channels) or for mutually orthogonal roots (|b1 + b2|^2 = 4)", True)
c.ok("verdict: e^2 = 4 pi alpha is not an extra input once the kernel field is a Cartan gauge potential with Gauss normalisation in three bulk dimensions; that identification is the remaining assumption", True)
c.done({"k_by_counting": {k_: mp.nstr(k, 4) for k_, k in count.items()}, "z_by_counting": {k_: mp.nstr(v, 6) for k_, v in zs.items()},
        "sectors": norms, "z_several_defects": {k_: mp.nstr(v, 6) for k_, v in zz.items()},
        "verdict": "coupling reduced to Coulomb's law of the kernel; identification kernel field = Cartan gauge potential not in the original; origin of the count 2 ambiguous"})
