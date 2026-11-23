"""Z16: Prüfung der drei externen Berichte vom 23. November (nach H6h8).

Die Berichte verweisen auf Pakete (PDF, ZIP), die nicht im Repo liegen. Geprüft wird, was sich aus dem
Text nachrechnen lässt:
 1. Einheitswindung als relative Chernzahl der Pfaffianlinie; explizite Masse M(z).
 2. Äquivariante Abbildung T vom Fünferträger auf die drei Familien.
 3. Verzweigte Doppelüberlagerung w^2 = s(z): vier Marken, zwei Blätter, Torus, mit Vierteldrehung tau = i.
 4. Zurückgezogene Bündel auf dem Nahttorus: gerade und ungerade Schnitte.
 5. Ginsparg-Wilson-Masse: verschwindet die Phasenbarriere aus Z9?
 6. Hopf-Modell: ein geladenes Ordnungsfeld auf S^3 erzwingt genau eine geschlossene Naht.
 7. Elternanomalie im Paarmodell; Sugawara-Rest; Gibbs-Variationsprinzip.
 8. Wilson-Parameter und die Untergrenze für a/l_P.
 9. Datenvergleiche (nur Arithmetik; die Messwerte stammen aus den Berichten).
"""
import json
from pathlib import Path

import numpy as np
import sympy as sp

OK, BEF = [], {}


def check(name, cond):
    OK.append((name, bool(cond)))
    return "ja" if cond else "NEIN"


print("1. Einheitswindung und Pfaffianlinie\n")
z = sp.symbols("z")
M = sp.Rational(1, 2) * sp.Matrix([[1 + z, sp.I * (z - 1)], [sp.I * (z - 1), -1 - z]])
detM = sp.simplify(M.det())
M2 = sp.simplify(M * M)
print(f"   det M(z) = {detM}, M(z)^2 = {M2.tolist()}")
check("det M = -z", sp.simplify(detM + z) == 0)
check("M^2 = z I", sp.simplify(M2 - z * sp.eye(2)) == sp.zeros(2))
for zz in (0.0, 0.3 + 0.4j, 1.0, 2.0j):
    Mn = np.array(M.subs(z, zz).evalf().tolist(), dtype=complex)
    sv = np.sort(np.linalg.svd(Mn, compute_uv=False))
    print(f"   z = {zz}: Singulärwerte {np.round(sv, 6)} (erwartet 1 und |z| = {abs(zz):.3f})")
sv_ok = all(np.allclose(np.sort(np.linalg.svd(np.array(M.subs(z, zz).evalf().tolist(), dtype=complex), compute_uv=False)),
                        np.sort([1.0, abs(zz)])) for zz in (0.0, 0.3 + 0.4j, 2.0j))
check("Singulärwerte {1, |z|}", sv_ok)
isy = np.array([[0, 1], [-1, 0]], complex)
Mn = np.array(M.subs(z, 0.6 + 0.2j).evalf().tolist(), dtype=complex)
Delta = np.kron(Mn, isy)
# Pfaffian of a 4x4 antisymmetric matrix
a = Delta
pf = a[0, 1] * a[2, 3] - a[0, 2] * a[1, 3] + a[0, 3] * a[1, 2]
print(f"   Pf(M x i sigma2) = {pf:.6f}, det M = {np.linalg.det(Mn):.6f}")
check("Pf = det M", abs(pf - np.linalg.det(Mn)) < 1e-12)
phis = np.linspace(0, 2 * np.pi, 400)
dets = [complex(detM.subs(z, np.exp(1j * p))) for p in phis]
wind = np.sum(np.diff(np.unwrap(np.angle(dets)))) / (2 * np.pi)
print(f"   Windung von det M auf dem Rand |z| = 1: {wind:+.3f}; im Inneren verliert M bei z = 0 genau einen Rang")
check("Windung 1 und ein Rangverlust", abs(wind - 1) < 1e-6 and np.linalg.matrix_rank(np.array(M.subs(z, 0).tolist(), dtype=complex)) == 1)
print("   Normierung bei 16 Paaren: primitiver Block 1, gesamte Pfaffianlinie 16, Determinante der vollen")
print("   Paarungsmatrix 32 (= Pf^2). Bestätigt: Innerhalb des Sektors Windung 1 ist der Rangverlust erzwungen;")
print("   die Wahl des Sektors bleibt offen.")

print("\n2. Äquivariante Abbildung T: H^0(O(D)) -> H^0(K(D)), D = mu4\n")
carrier = [z ** k / (z ** 4 - 1) for k in range(5)]                 # f = Q/(z^4 - 1), deg Q <= 4
fam = [z ** k / (z ** 4 - 1) for k in range(3)]                     # omega = P dz/(z^4 - 1), deg P <= 2


def P0(f):
    return sp.simplify(sum(f.subs(z, sp.I ** j * z) for j in range(4)) / 4)


def T(f):                                                            # coefficient of dz in (f - P0 f) dz/z
    return sp.simplify((f - P0(f)) / z)


def coords_fam(expr):
    num = sp.expand(sp.cancel(expr * (z ** 4 - 1)))
    poly = sp.Poly(num, z)
    assert poly.degree() <= 2, poly
    return [poly.coeff_monomial(z ** k) for k in range(3)]


Tmat = sp.Matrix([coords_fam(T(f)) for f in carrier]).T
rank = Tmat.rank()
kernel = Tmat.nullspace()
print(f"   Matrix von T (3 x 5): {Tmat.tolist()}, Rang {rank}, Kern {[list(v) for v in kernel]}")
check("T surjektiv, Kern 2-dimensional", rank == 3 and len(kernel) == 2)
Rc = sp.diag(*[sp.I ** k for k in range(5)])                         # f(z) -> f(iz) on carrier
Rf = sp.diag(*[sp.I * sp.I ** k for k in range(3)])                  # omega -> i P(iz) dz/(z^4-1)
equiv = sp.simplify(Tmat * Rc - Rf * Tmat) == sp.zeros(3, 5)
print(f"   T vertauscht mit der Vierteldrehung: {equiv}; Charaktere Träger {[sp.I ** k for k in range(5)]} -> Familien {list(Rf.diagonal())}")
check("T äquivariant", equiv)
print("   Die Familiencharaktere (i, -1, -i) sind der nichttriviale Charakteranteil des Fünferträgers.")

print("\n3. Verzweigte Doppelüberlagerung w^2 = s(z) mit s in H^0(O(4))\n")
g = sp.symbols("g")
r_needed = sp.solve(sp.Eq(2 * 1 - 2, 2 * (2 * 0 - 2) + g), g)
print(f"   Riemann-Hurwitz für g(E) = 1 über CP^1: Zahl der Verzweigungspunkte r = {r_needed}")
check("r = 4", r_needed == [4])
a_, b_, c_, d_, e_ = sp.symbols("a b c d e")
s = a_ * z ** 4 + b_ * z ** 3 + c_ * z ** 2 + d_ * z + e_
inv = sp.solve(sp.Poly(sp.expand(s.subs(z, sp.I * z) - s), z).coeffs(), [b_, c_, d_], dict=True)
print(f"   unter z -> iz invariante Quartiken: {inv} -> s = a z^4 + e")
check("invariante Quartik a z^4 + e", inv == [{b_: 0, c_: 0, d_: 0}])
roots = [1, sp.I, -1, -sp.I]
lam = sp.simplify((roots[0] - roots[2]) * (roots[1] - roots[3]) / ((roots[0] - roots[3]) * (roots[1] - roots[2])))
jval = sp.simplify(256 * (lam ** 2 - lam + 1) ** 3 / (lam ** 2 * (lam - 1) ** 2))
print(f"   Doppelverhältnis der Marken mu4 = {lam}, j = {jval}")
check("j = 1728", jval == 1728)
print("   omega = dz/w geht unter z -> iz (w fest) in i omega über: komplexe Multiplikation mit i, tau = i.")

print("\n4. Zurückgezogene Bündel auf dem Nahttorus E (pi: E -> CP^1 vom Grad 2)\n")
for name, dgr in (("Familienbündel O(2) = K(D)", 2), ("Trägerbündel O(4) = O(D)", 4)):
    even = dgr + 1                                                   # pullbacks f(z), deg f <= d
    odd = max(dgr - 2 + 1, 0)                                        # w * g(z), deg g <= d - 2
    print(f"   {name}: Grad auf E = {2 * dgr}, h^0 = {2 * dgr} = {even} gerade + {odd} ungerade")
    if dgr == 2:
        check("O(2): 3 gerade + 1 ungerade", (even, odd) == (3, 1))
    else:
        check("O(4): 5 gerade + 3 ungerade (Atlas 19.11.)", (even, odd) == (5, 3))
print("   Ein Elternoperator muss die leichten Fermionen an den geraden Teil koppeln, um drei und nicht vier")
print("   Familien zu bekommen. Auf dem Trägerbündel sind die 5 geraden der Träger, die 3 ungeraden w * H^0(O(2)).")

print("\n5. Ginsparg-Wilson-Masse gegen Wilson-Masse (freier 4D Dirac, Impulsraum)\n")
gam = [np.array(x, complex) for x in (
    [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]],
    [[0, 0, 0, -1j], [0, 0, -1j, 0], [0, 1j, 0, 0], [1j, 0, 0, 0]],
    [[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]],
    [[0, 0, -1j, 0], [0, 0, 0, 1j], [1j, 0, 0, 0], [0, -1j, 0, 0]])]
g5 = gam[0] @ gam[1] @ gam[2] @ gam[3]
I4 = np.eye(4)


def expm_g5(t):
    return np.cos(t) * I4 + 1j * np.sin(t) * g5


def D_ov(k, rho=1.0):
    DW = sum(1j * gam[m] * np.sin(k[m]) for m in range(4)) + (sum(1 - np.cos(k[m]) for m in range(4)) - rho) * I4
    X = DW
    sq = np.linalg.inv(np.real_if_close(np.linalg.cholesky(X.conj().T @ X)).conj().T)
    w, U = np.linalg.eigh(X.conj().T @ X)
    inv_sqrt = U @ np.diag(1 / np.sqrt(w)) @ U.conj().T
    return I4 + X @ inv_sqrt


rng = np.random.default_rng(1)
resid = 0.0
for _ in range(5):
    k = rng.uniform(-np.pi, np.pi, 4)
    Dk = D_ov(k)
    gw = g5 @ Dk + Dk @ g5 - Dk @ g5 @ Dk
    resid = max(resid, np.max(np.abs(gw)))
    g5hat = g5 @ (I4 - Dk)
    Cm = I4 - Dk / 2
    m, th = 0.3, 1.1
    lhs = Dk + m * expm_g5(th) @ Cm
    w_, U_ = np.linalg.eig(g5hat)
    rhs = expm_g5(th / 2) @ (Dk + m * Cm) @ (U_ @ np.diag(np.exp(1j * th / 2 * w_)) @ np.linalg.inv(U_))
    resid = max(resid, np.max(np.abs(lhs - rhs)))
print(f"   Ginsparg-Wilson-Relation und Drehidentität D_m(theta) = e^(i theta g5/2) D_m(0) e^(i theta g5hat/2): Rest {resid:.1e}")
check("Drehidentität exakt", resid < 1e-10)
n = 8
ks = 2 * np.pi * (np.arange(n) + 0.5) / n
grid = [(a1, a2, a3, a4) for a1 in ks for a2 in ks for a3 in ks for a4 in ks]
V_ov, V_w = [], []
for th in np.linspace(0, np.pi, 5):
    sov = sw = 0.0
    for k in grid:
        Dk = D_ov(np.array(k))
        sov += -np.log(abs(np.linalg.det(Dk + 0.3 * expm_g5(th) @ (I4 - Dk / 2))))
        DW = sum(1j * gam[m] * np.sin(k[m]) for m in range(4)) + sum(1 - np.cos(k[m]) for m in range(4)) * I4
        sw += -np.log(abs(np.linalg.det(DW + 0.3 * expm_g5(th))))
    V_ov.append(sov / len(grid))
    V_w.append(sw / len(grid))
span_ov, span_w = np.ptp(V_ov), np.ptp(V_w)
print(f"   Phasenpotential je Zelle: Overlap-Spannweite {span_ov:.1e}, Wilson-Spannweite {span_w:.4f}")
check("Overlap ohne Phasenbarriere", span_ov < 1e-10 and span_w > 0.01)
print("   Die Kritik hat recht: Die Barriere aus Z9 hängt am Wilson-Glied. Mit einem Ginsparg-Wilson-Operator")
print("   bleibt die Phasensymmetrie exakt. Für die Majoranapaarung von H6 ist das nicht gerechnet.")
BEF["Wand aus Z9 regulatorabhängig"] = True

print("\n6. Hopf-Modell: Ordnungsfeld mit Ladung 1 auf S^3\n")
x = sp.symbols("x0:4", real=True)
z0, z1 = x[0] + sp.I * x[1], x[2] + sp.I * x[3]
lap = lambda f: sum(sp.diff(f, xi, 2) for xi in x)
print(f"   z0, z1 harmonisch in R^4: {lap(z0) == 0 and lap(z1) == 0}; homogen vom Grad 1, also Eigenwert l(l+2) = 3 auf S^3")
check("lineare Moden harmonisch", lap(z0) == 0 and lap(z1) == 0)
print("   Funktionen mit Ladung 1 unter (z0, z1) -> e^(i a)(z0, z1) haben l >= 1; bei l = 1 sind es genau z0, z1.")
A, B = 0.6 + 0.3j, -0.2 + 0.7j
t = np.linspace(0, 2 * np.pi, 2000)
# zero set of A z0 + B z1 on S^3: (z0, z1) proportional to (B, -A): a great circle
v = np.array([B, -A]) / np.linalg.norm([B, -A])
w_perp = np.array([np.conj(A), np.conj(B)]) / np.linalg.norm([A, B])
eps = 0.05
loop = [(np.exp(1j * 0.3) * v * np.sqrt(1 - eps ** 2) + eps * np.exp(1j * tt) * w_perp) for tt in t]
vals = [A * p[0] + B * p[1] for p in loop]
wnd = np.sum(np.diff(np.unwrap(np.angle(vals)))) / (2 * np.pi)
print(f"   Nullmenge von A z0 + B z1: ein Großkreis; Windung um ihn: {wnd:+.3f}")
check("eine Naht mit Windung 1", abs(abs(wnd) - 1) < 1e-6)
print("   Verträglich mit Z11: H_2(S^3) = 0, also schneidet jede geschlossene Fläche den Kreis algebraisch null-")
print("   mal (gerade). Erzwungen wird die Naht durch die Ladung 1 unter der Hopf-U(1), nicht durch Topologie")
print("   des Raums allein. Offen: warum das Ordnungsfeld primitive Ladung trägt.")
BEF["Hopf-Modell erzwingt genau eine Naht bei Ladung 1"] = True

print("\n7. Elternanomalie, Sugawara, Gibbs\n")
for npairs in (1, 16):
    nu_pair = npairs * (1 - 1)
    print(f"   Paarmodell mit {npairs} Paaren (Ladungen +1 und -1): nu = n1 - n3 = {nu_pair} mod 16")
BEF["Paarmodell für jede Paarzahl Z16-anomaliefrei"] = True
BEF["Z8: 16 aus Z16 nur in der partnerfreien Fassung"] = True
print("   Korrektur zu Z8: Im Paarmodell erzwingt die Elternanomalie keine 16. Nur in der partnerfreien")
print("   Fassung (16 Weyl mit Ladung 1, Masse aus Phi mit Ladung 2) ist nu = 16 die Z16-Bedingung.")
for k_ in (1, 2, 3):
    csug = sp.Rational(248 * k_, k_ + 30)
    print(f"   {k_} Blöcke: c = {8 * k_}, diagonales Niveau {k_}, Sugawara {csug}, Rest {8 * k_ - csug}")
check("Sugawara-Rest 0, 1/2, 16/11", [8 * k_ - sp.Rational(248 * k_, k_ + 30) for k_ in (1, 2, 3)] == [0, sp.Rational(1, 2), sp.Rational(16, 11)])
coef = np.array([1, 248, 4124, 34752, 213126, 1057504, 4530744, 17333248, 60655377, 197230000, 603096260,
                 1749556736], float)                                  # chi_E8 = q^(-1/3) sum coef_n q^n
H = np.arange(len(coef)) - 1 / 3
pstar = coef * np.exp(-2 * np.pi * H)
pstar /= pstar.sum()
S_lvl = lambda d: -np.sum(d * np.log(d / coef))                      # entropy with level degeneracies
D_rel = lambda d, e: np.sum(d * np.log(d / e))
rng7 = np.random.default_rng(7)
gaps = []
for _ in range(5):
    u = rng7.normal(size=len(coef))
    e1 = np.ones_like(H)
    e2 = H - (pstar * H).sum()                                          # weighted Gram-Schmidt w.r.t. pstar
    for e_vec in (e1, e2):
        u = u - (pstar * u * e_vec).sum() / (pstar * e_vec * e_vec).sum() * e_vec
    rho = pstar * (1 + 0.5 * u / np.max(np.abs(u)))
    gaps.append((S_lvl(pstar) - S_lvl(rho), D_rel(rho, pstar), rho @ H))
print(f"   Mittlere Energie des KMS-Zustands bei beta = L: {pstar @ H:+.2e}")
for g_, d_, e_ in gaps[:3]:
    print(f"   Zustand mit <H> = {e_:+.1e}: S(rho*) - S(rho) = {g_:.6f}, D(rho || rho*) = {d_:.6f}")
check("Gibbs: S(rho*) - S(rho) = D >= 0 bei <H> = 0", all(abs(g_ - d_) < 1e-9 and g_ > 0 and abs(e_) < 1e-9 for g_, d_, e_ in gaps))
print("   Unter allen Zuständen mit mittlerer Energie null hat der KMS-Zustand bei beta = L eindeutig die größte")
print("   Entropie. Das macht K6 zu einem Variationsprinzip; eine exakte Nullenergieschale gibt es nicht (n - 1/3).")

print("\n8. Wilson-Parameter und Untergrenze für a/l_P\n")
src15 = Path(__file__).with_name("z15_a_lP_vollinhalt.py").read_text()
ns15 = {"__file__": str(Path(__file__).with_name("z15_a_lP_vollinhalt.py"))}
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    exec(src15[:src15.index("# --------------------------------------------------------------------------- Teil 1")], ns15)
kd = ns15["kappa_dirac"]
kap = {r: kd(0.0, r, nk=64, ring=120) for r in (1.0, 2.0, 4.0, 8.0)}
for r, kv in kap.items():
    print(f"   Wilson r = {r}: kappa je Diracfeld {kv:.5f}; nur 24 leichte Dirac: a/l_P = {2 * np.sqrt(24 * kv):.3f}")
check("r = 1 wie Z15 und Bericht", abs(kap[1.0] - 0.2209) < 0.003)
check("r = 4 wie Bericht", abs(kap[4.0] - 0.0434) < 0.003)
kap_s = 0.024262
light_bos = 12 * 2 * kap_s + 4 * kap_s
low = 2 * np.sqrt(light_bos)
print(f"   Bosonischer Sockel ohne Randterm (12 Eichbosonen x 2 Polarisationen + 4 Higgs): kappa = {light_bos:.3f},"
      f" a/l_P >= {low:.2f}")
BEF["Untergrenze 3,3 aus Z15 gilt nicht für großes Wilson-r"] = True
print("   Korrektur zu Z15: Die Untergrenze 3,3 gilt nur für die dort gerechneten Regulatoren. Für großes r sinkt")
print(f"   der Fermionanteil; ohne Randterm bleibt nur der bosonische Sockel, a/l_P >= {low:.2f}. Robust ist a/l_P > 1.")

print("\n9. Datenvergleich (Messwerte aus den Berichten, hier nicht geprüft)\n")
phi0 = 1 / (6 * np.pi) + 3 / (256 * np.pi ** 4)
s12, s13 = 1 / 3 - phi0 / 2, phi0 * np.exp(-5 / 6)
s2_2t13 = 4 * s13 * (1 - s13)
for name, pred, val, err in (("sin^2 th12 JUNO", s12, 0.3036, 0.0064), ("sin^2 th13 Bari", s13, 0.02218, 0.00055),
                             ("sin^2 th13 NuFIT 6.1", s13, 0.02249, 0.00057), ("sin^2 2th13 Daya Bay", s2_2t13, 0.0851, 0.0024)):
    print(f"   {name:22s}: TFPT {pred:.6f}, Messung {val} +- {err}, Abstand {(pred - val) / err:+.2f} sigma")
check("Arithmetik der Abstände", abs((s13 - 0.02249) / 0.00057 - 1.08) < 0.02 and abs((s2_2t13 - 0.0851) / 0.0024 - 2.17) < 0.02)
print("   delta = 240 Grad steht im TFPT-Register als zugewiesener Texturwert (REG.FREEZE.01), nicht hergeleitet.")
print("   Korrektur: Für theta13 sind Reaktorexperimente wie Daya Bay maßgeblich, nicht JUNO.")

print("\nZusammenfassung")
for n_, c in OK:
    print(f"   {'bestanden' if c else 'NICHT bestanden'}: {n_}")
Path(__file__).with_name("z16_ergebnis.json").write_text(json.dumps({
    "checks": {n_: c for n_, c in OK}, "befunde": BEF, "overlap_span": float(span_ov), "wilson_span": float(span_w),
    "kappa_wilson": {str(k): float(v) for k, v in kap.items()}, "bosonic_floor_a_lP": float(low)}, ensure_ascii=False, indent=1))
assert all(c for _, c in OK)
