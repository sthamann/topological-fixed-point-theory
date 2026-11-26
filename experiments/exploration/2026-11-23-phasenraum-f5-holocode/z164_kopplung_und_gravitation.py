"""Z164: Wege W1 (Kopplung) und W2 (Gravitationsvertex) im Rahmen P_TW (Twistorwirkung), erste Rechnungen.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Bestand: z110 (Einspuranteil = -k mal Parke-Taylor; g_eff = k^(-1/2) unter P_norm), z112 (native Klammern, Hodges;
*(I_pi - I_omega) = (1 + det x) <ij>), z123 (Kinetik nativ, c_V N = 1 aus Faktorisierung), z163 (Vertex = log Z,
Zähler aus dem Gewicht; Anomaliebedingung nur E8). Ledger: alpha-Gleichung (Paper 1), Wirkungsleiter, CELEST.SEAM.01.

Teil A  W1: Wie viele freie Zahlen hat der Rahmen? Skalierung von b; Zähler <12>^4 als Berezin-Integral (b als oberste
        Komponente eines Zusammenhangs)
Teil B  W1: P_norm (g = k^(-1/2) = 1 an der Nahtskala) gegen die gemessenen Kopplungen (Einschleifen-Lauf des
        Standardmodells bis zur reduzierten Planckmasse); gegen Zahlenspielerei markiert
Teil C  W2: Gewichte für Spin 2; Hodges-Determinante: Unabhängigkeit von der Wahl, Skalierung, Zähler <ij>^8
Teil D  W2: Nativer Unendlichkeitstwistor I_pi - I_omega: Nichtentartung, Krümmung, Kontraktion zum flachen Fall
Teil E  Bilanz und Urteil.
"""
import itertools
import json
import math
import os
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import numpy as np
import sympy as sp
from sympy.combinatorics import Permutation

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
PI = math.pi


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   -> {name}: {text}")


say("Z164: W1 (Kopplung) und W2 (Gravitationsvertex) im Rahmen der Twistorwirkung")
say("=" * 100)

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil A  W1: freie Zahlen des Rahmens")
say("   Wirkung S = A (b, dbar a + a a) + V Summe über Geraden (b K b K <12>^4); a ist ein Zusammenhang (keine freie")
say("   Normierung), b nicht. Ausintegrieren von b (quadratisch) gibt die Yang-Mills-Wirkung mit 1/g^2 ~ A^2/V.")
A_, V_, s_, F_, b_ = sp.symbols("A V s F b", positive=True)
S_ab = A_ * b_ * F_ - V_ * b_ ** 2
b_star = sp.solve(sp.diff(S_ab, b_), b_)[0]
S_eff = sp.simplify(S_ab.subs(b_, b_star))
check("A1 Ausintegrieren von b: S_eff = (A^2/(4 V)) F^2; unter b -> s b ist nur A^2/V invariant [E]",
      sp.simplify(S_eff - A_ ** 2 * F_ ** 2 / (4 * V_)) == 0
      and sp.simplify(((A_ * s_) ** 2 / (V_ * s_ ** 2)) - A_ ** 2 / V_) == 0)
# Berezin: chi^a(sigma) = theta^a_alpha sigma^alpha, a = 1..4; Koeffizient von prod theta in chi^4(s1) chi^4(s2)
s1 = sp.symbols("p1 p2")
s2 = sp.symbols("r1 r2")
br12 = s1[0] * s2[1] - s1[1] * s2[0]
# je Index a: chi(s1) chi(s2) = theta_1 theta_2 (s1^1 s2^2 - s1^2 s2^1); die vier Faktoren sind gerade und vertauschen
per_a = sp.expand(s1[0] * s2[1] - s1[1] * s2[0])       # Koeffizient von theta^a_1 theta^a_2
# Umordnung chi1(s1)..chi4(s1) chi1(s2)..chi4(s2) -> prod_a [chi_a(s1) chi_a(s2)]: Vorzeichen der Permutation
perm = [0, 4, 1, 5, 2, 6, 3, 7]
sign = Permutation(perm).signature()
coeff = sign * per_a ** 4
check("A2 Berezin-Integral: int d^8 theta chi^4(s1) chi^4(s2) = +<12>^4, Koeffizient eins [E]",
      sp.expand(coeff - br12 ** 4) == 0 and sign == 1)
say("   Lesart [L]: Schreibt man A = a + chi^4 b (b als oberste Komponente eines Zusammenhangs mit vier formalen")
say("   Gewicht-1-Variablen), dann ist der Vertex genau int d^4x d^8theta log det(dbar + A) mit W = k log det der Naht.")
say("   Der Koeffizient des Vertex ist dann das Niveau k, und die einzige freie Zahl ist der Koeffizient 1/h des")
say("   holomorphen Teils: g^2 ~ k h^2 (b kanonisch) bzw. g_eff^2 = 1/k in der Stromnormierung von z110 (P_norm).")
befund("A_freie_zahl", 1, "Der Rahmen hat genau eine freie reine Zahl (Kopplung des holomorphen Teils relativ zum "
       "Niveau der Naht). Verhältnisse von Kopplungen sind Niveauverhältnisse und fest; der Gesamtwert ist es "
       "nicht. Unter P_norm (z110) ist er g^2 = 1/k = 1 [E für die Skalierung, C unter P_norm].")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil B  W1: P_norm gegen die gemessenen Kopplungen")
say("   Eingaben von außen (PDG, gerundet): 1/alpha_em(M_Z) = 127,95; sin^2 theta_W = 0,2312; alpha_s(M_Z) = 0,1180;")
say("   M_Z = 91,19 GeV; reduzierte Planckmasse 2,435e18 GeV. Einschleifen-Lauf mit b = (41/10, -19/6, -7).")
AEM, SW2, AS, MZ, MBAR = 127.95, 0.2312, 0.1180, 91.19, 2.435e18
inv_mz = {"U(1) (GUT-Norm)": 0.6 * (1 - SW2) * AEM, "SU(2)": SW2 * AEM, "SU(3)": 1 / AS}
bcoef = {"U(1) (GUT-Norm)": 41 / 10, "SU(2)": -19 / 6, "SU(3)": -7.0}
L = math.log(MBAR / MZ)
inv_mbar = {k: inv_mz[k] - bcoef[k] / (2 * PI) * L for k in inv_mz}
for k in inv_mz:
    say(f"   {k:16s} 1/alpha(M_Z) = {inv_mz[k]:7.3f}   1/alpha(M_Planck,red) = {inv_mbar[k]:7.3f}   "
        f"1/g^2 = {inv_mbar[k] / (4 * PI):.3f}")
ZAHLEN.update(inv_alpha_mbar={k: round(v, 3) for k, v in inv_mbar.items()},
              inv_g2_mbar={k: round(v / (4 * PI), 3) for k, v in inv_mbar.items()})
check("B1 Der Koeffizient 41/10 der alpha-Gleichung ist der Einschleifen-Koeffizient von U(1) (Bestand, v434)",
      abs(bcoef["U(1) (GUT-Norm)"] - 4.1) < 1e-12)
check("B2 [N negativ] P_norm (1/g^2 = k = 1, also 1/alpha = 4 pi = 12,57 an der Nahtskala) verfehlt alle drei "
      "Kopplungen um den Faktor 2,7 bis 4,0", all(2.6 < v / (4 * PI) < 4.1 for v in inv_mbar.values()))
dev = {k: inv_mbar[k] / (16 * PI) - 1 for k in inv_mbar}
ZAHLEN["abw_von_16pi"] = {k: round(v, 4) for k, v in dev.items()}
say("   Abweichung von 1/alpha = 16 pi = 2/c3 (also g^2 = 1/4): "
    + ", ".join(f"{k}: {100 * v:+.1f} %" for k, v in dev.items()))
check("B3 [N, Hinweis] Die beiden nichtabelschen Kopplungen liegen bei M_Planck,red bei 1/g^2 = 3,9 bis 4,0 "
      "(16 pi auf 3,3 % und 0,7 %); U(1) nicht (2,73)", abs(dev["SU(2)"]) < 0.035 and abs(dev["SU(3)"]) < 0.01
      and abs(dev["U(1) (GUT-Norm)"]) > 0.3)
as_pred = 1 / (16 * PI + bcoef["SU(3)"] / (2 * PI) * L)
shift = {k: -bcoef[k] / (2 * PI) * math.log(math.sqrt(8 * PI)) for k in bcoef}
ZAHLEN.update(alpha_s_aus_16pi=round(as_pred, 4), verschiebung_bei_M_Planck={k: round(v, 2) for k, v in shift.items()})
say(f"   Gegenproben: alpha_s(M_Z) aus 16 pi bei M_Planck,red wäre {as_pred:.4f} (gemessen 0,1180, +4 %). Wahl der")
say(f"   Skala: bei M_Planck statt M_Planck,red verschieben sich die Werte um {shift['SU(2)']:+.2f} (SU(2)) und "
    f"{shift['SU(3)']:+.2f} (SU(3)).")
# Zweischleifen-Lauf (nur Eichanteil, ohne Yukawa-Terme und Schwellen), RK4 in ln mu
BIJ = np.array([[199 / 50, 27 / 10, 44 / 5], [9 / 10, 35 / 6, 12.0], [11 / 10, 9 / 2, -26.0]])
BI = np.array([41 / 10, -19 / 6, -7.0])
rhs = lambda y: -BI / (2 * PI) - BIJ @ (1 / y) / (8 * PI ** 2)
y2 = np.array([inv_mz[k] for k in ("U(1) (GUT-Norm)", "SU(2)", "SU(3)")])
NST = 4000
hstep = L / NST
for _ in range(NST):
    k1 = rhs(y2); k2 = rhs(y2 + hstep / 2 * k1); k3 = rhs(y2 + hstep / 2 * k2); k4 = rhs(y2 + hstep * k3)
    y2 = y2 + hstep / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
dev2 = y2 / (16 * PI) - 1
ZAHLEN.update(inv_alpha_mbar_2schleifen=[round(float(v), 3) for v in y2],
              abw_von_16pi_2schleifen=[round(float(v), 4) for v in dev2])
say("   Zweischleifen-Lauf (nur Eichanteil): 1/alpha(M_Planck,red) = "
    + ", ".join(f"{v:.3f}" for v in y2) + "; Abweichung von 16 pi: "
    + ", ".join(f"{100 * v:+.1f} %" for v in dev2))
check("B4 [N] Auch mit dem Zweischleifen-Eichanteil bleiben SU(2) und SU(3) innerhalb von 5 % bei 16 pi; U(1) nicht",
      abs(dev2[1]) < 0.05 and abs(dev2[2]) < 0.05 and abs(dev2[0]) > 0.25)
t_x = (inv_mz["SU(2)"] - inv_mz["SU(3)"]) / ((bcoef["SU(2)"] - bcoef["SU(3)"]) / (2 * PI))
inv_x = inv_mz["SU(3)"] - bcoef["SU(3)"] / (2 * PI) * t_x
sw2_pred = (16 * PI + bcoef["SU(2)"] / (2 * PI) * L) / AEM
ZAHLEN.update(schnitt_skala_GeV=float(f"{MZ * math.exp(t_x):.2e}"), schnitt_inv_alpha=round(inv_x, 2),
              sin2_aus_16pi=round(sw2_pred, 4))
say(f"   SU(2) und SU(3) schneiden sich (eine Schleife) bei {MZ * math.exp(t_x):.1e} GeV mit 1/alpha = {inv_x:.2f}; aus")
say(f"   16 pi bei M_Planck,red folgte sin^2 theta_W(M_Z) = {sw2_pred:.4f} (gemessen 0,2312, {100 * (sw2_pred / SW2 - 1):+.1f} %).")
check("B5 [N negativ] Als exakte Gleichung an einer gemeinsamen Skala ist 1/alpha_2 = 1/alpha_3 = 16 pi nicht "
      "haltbar: Die Kurven schneiden sich bei 1/alpha = 47, nicht 50,3; die Vorhersagen für alpha_s und "
      "sin^2 theta_W liegen 4 bis 6 % daneben", abs(inv_x - 16 * PI) > 3 and abs(sw2_pred / SW2 - 1) > 0.04)
say("   Zählung: eine ganze Zahl (4) und eine Skala gewählt; zwei von drei Kopplungen passen auf wenige Prozent; die")
say("   beiden sind korreliert (bekannte Beinahe-Vereinigung). Ohne Yukawa-Terme und Schwellen. Schwacher Hinweis.")
befund("B_pnorm", {"p_norm_widerlegt": True, "hinweis_g2": "1/4"},
       "P_norm in der Form g = 1 an der Nahtskala ist gegen die Daten widerlegt. Die Daten verlangen für die "
       "nichtabelschen Faktoren 1/g^2 nahe 4. Ein nativer Kandidat für die 4 ist die Quellladung |mu_4| = 4 aus "
       "CELEST.WP5E.M1.01 (nicht hergeleitet) [N; Hinweis L].")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil C  W2: Gewichte für Spin 2 und Hodges-Determinante")
w, hs = sp.symbols("w h")
check("C1 Serre-Paarung auf dem Twistorraum (K = O(-4)): O(w) mit O(-w-4); Helizität h hat Gewicht 2h - 2. "
      "Spin 1: 0 und -4; Spin 2: 2 und -6 [E]", [(2 * x - 2, -(2 * x - 2) - 4) for x in (1, 2)] == [(0, -4), (2, -6)]
      and all(-(2 * x - 2) - 4 == 2 * (-x) - 2 for x in (1, 2)))


def kin(n, seed):
    rng = np.random.default_rng(seed)
    lam = rng.normal(size=(n, 2)) + 1j * rng.normal(size=(n, 2))
    lt = rng.normal(size=(n, 2)) + 1j * rng.normal(size=(n, 2))
    M = np.array([[lam[0, 0], lam[1, 0]], [lam[0, 1], lam[1, 1]]])
    for d in range(2):                                  # Impulserhaltung: löse nach lt_1, lt_2
        rhs = -sum(lam[a] * lt[a, d] for a in range(2, n))
        lt[0, d], lt[1, d] = np.linalg.solve(M, rhs)
    return lam, lt


ang = lambda lam, a, b: lam[a, 0] * lam[b, 1] - lam[a, 1] * lam[b, 0]


def hodges(lam, lt, rem_r, rem_c, x, y):
    n = len(lam)
    Phi = np.zeros((n, n), complex)
    for a in range(n):
        for b in range(n):
            if a != b:
                Phi[a, b] = ang(lt, a, b) / ang(lam, a, b)
    for a in range(n):
        Phi[a, a] = -sum(ang(lt, a, b) * (lam[b] @ np.array([x[1], -x[0]])) * (lam[b] @ np.array([y[1], -y[0]]))
                         / (ang(lam, a, b) * (lam[a] @ np.array([x[1], -x[0]])) * (lam[a] @ np.array([y[1], -y[0]])))
                         for b in range(n) if b != a)
    rows = [i for i in range(n) if i not in rem_r]
    cols = [i for i in range(n) if i not in rem_c]
    sgn = (-1) ** (sum(rem_r) + sum(rem_c))
    c = lambda t: ang(lam, t[0], t[1]) * ang(lam, t[1], t[2]) * ang(lam, t[2], t[0])
    return sgn * np.linalg.det(Phi[np.ix_(rows, cols)]) / (c(rem_r) * c(rem_c))


worst_ind, worst_sc = 0.0, 0.0
for n in (5, 6):
    lam, lt = kin(n, 164 + n)
    x, y = np.array([0.3, 1.1 + 0.2j]), np.array([1.0, -0.7 + 0.5j])
    ref = hodges(lam, lt, (0, 1, 2), (0, 1, 2), x, y)
    for rr, cc in (((0, 1, 2), (2, 3, 4)), ((1, 2, 4), (0, 2, 3)), ((0, 3, 4), (0, 3, 4))):
        worst_ind = max(worst_ind, abs(hodges(lam, lt, rr, cc, x, y) / ref - 1))
    worst_ind = max(worst_ind, abs(hodges(lam, lt, (0, 1, 2), (0, 1, 2), y + 0.4, x - 0.2j) / ref - 1))
    t = 1.3 - 0.4j
    for a in range(n):
        l2, lt2 = lam.copy(), lt.copy()
        l2[a] *= t
        lt2[a] /= t
        worst_sc = max(worst_sc, abs(hodges(l2, lt2, (0, 1, 2), (0, 1, 2), x, y) / ref * t ** 4 - 1))
ZAHLEN.update(hodges_unabh=float(f"{worst_ind:.1e}"), hodges_skalierung=float(f"{worst_sc:.1e}"))
check("C2 Reduzierte Hodges-Determinante: unabhängig von den gestrichenen Zeilen/Spalten und den Referenzspinoren "
      "(n = 5, 6) [N exakt bis Rundung]", worst_ind < 1e-8)
check("C3 Sie skaliert je Bein mit t^-4 (Helizität +2). Für zwei Beine der Helizität -2 muss der Zähler t^+8 "
      "liefern: <ij>^8, ohne Supersymmetrie; allgemein <ij>^(4s) für Spin s (s = 1: 4, s = 2: 8)", worst_sc < 1e-8)
befund("C_spin2", {"zaehler": 8}, "Wie bei Spin 1 folgt der Zähler aus den Gewichten: Der Serre-Partner von O(2) ist "
       "O(-6), der Zähler ist <ij>^8. Eine eigene Prämisse (N = 8) ist dafür nicht nötig [E].")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil D  W2: nativer Unendlichkeitstwistor I_pi - I_omega")
xs = sp.symbols("x0:4", real=True)
pa, pb = sp.symbols("a1 a2"), sp.symbols("b1 b2")
eps = lambda u, v: u[0] * v[1] - u[1] * v[0]
sig = [sp.eye(2), sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
res_d = {}
for name, X in (("euklidisch (x0 + i x.sigma)", xs[0] * sig[0] + sp.I * sum((xs[i] * sig[i] for i in (1, 2, 3)), sp.zeros(2))),
                ("lorentzsch (t + x.sigma, hermitesch)", xs[0] * sig[0] + sum((xs[i] * sig[i] for i in (1, 2, 3)), sp.zeros(2)))):
    om = lambda p: sp.I * X * sp.Matrix(p)
    Ipi, Iom = eps(pa, pb), sp.expand(eps(om(pa), om(pb)))
    ratio = sp.simplify((Ipi - Iom) / Ipi)
    detx = sp.expand(X.det())
    eta = [1, 1, 1, 1] if name.startswith("eukl") else [-1, 1, 1, 1]      # Raumsignatur überwiegend plus
    Om = 1 + detx
    box = sum(e * sp.diff(Om, v, 2) for e, v in zip(eta, xs))
    grad2 = sum(e * sp.diff(Om, v) ** 2 for e, v in zip(eta, xs))
    sgn_metric = 1 if name.startswith("eukl") else -1    # det(dx) = -eta' im lorentzschen Schnitt
    # g = sgn * eta'/Omega^2; für konform flache Metrik eta/Omega^2 ist R = 6 Omega Box Omega - 12 |grad Omega|^2
    R = sp.simplify(sgn_metric * (6 * Om * box - 12 * grad2)) if sgn_metric == 1 else sp.simplify(6 * Om * box - 12 * grad2)
    res_d[name] = (sp.simplify(ratio - Om) == 0, detx, R)
    say(f"   {name}: det x = {detx}; (I_pi - I_omega)/I_pi auf L_x = 1 + det x; Skalarkrümmung R = {R}, Lambda = R/4 = {R / 4}")
check("D1 Auf L_x (omega = i x pi): I_omega = -det(x) I_pi, also (I_pi - I_omega) = (1 + det x) <12> wie in z112 [E]",
      all(v[0] for v in res_d.values()))
Rs = [v[2] for v in res_d.values()]
check("D2 Die Metrik det(dx)/(1 + det x)^2 hat konstante Krümmung: euklidischer Schnitt S^4 mit R = 48, "
      "lorentzscher Schnitt AdS_4 mit R = -48; |Lambda| = 12 in nativen Einheiten (Radius 1/2) [E]",
      Rs == [48, -48])
Imat = sp.Matrix([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])     # I_pi - I_omega auf (omega, pi)
e2 = sp.Symbol("epsilon2", positive=True)
Ieps = sp.Matrix([[0, -e2, 0, 0], [e2, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])
pf = lambda M: sp.simplify(M[0, 1] * M[2, 3] - M[0, 2] * M[1, 3] + M[0, 3] * M[1, 2])
check("D3 I_pi - I_omega ist nicht entartet (Pfaffsche -1): Einstein-Typ mit Lambda ungleich null. Die Familie "
      "I_pi - eps^2 I_omega hat Pfaffsche -eps^2 und wird genau bei eps = 0 entartet (flacher Fall, Rang 2) [E]",
      pf(Imat) == -1 and pf(Ieps) == -e2 and Ieps.subs(e2, 0).rank() == 2)
detx_s = sp.Symbol("D")
Om_e = 1 + e2 * detx_s
check("D4 Mit I_pi - eps^2 I_omega ist der Faktor 1 + eps^2 det x, also |Lambda| = 12 eps^2 [E]",
      sp.simplify(Om_e.subs(e2, 1) - (1 + detx_s)) == 0)
AINV = 137.0359992168407
lam_over_m2 = 3 / (4 * PI ** 2) * math.exp(-2 * AINV)          # Lambda/M^2 = rho_Lambda/M^4 (Ledger-Leiter)
eps2 = {a2g: lam_over_m2 * a2g / (8 * PI) / 12 for a2g in (0.2, 24.0)}
ZAHLEN.update(log10_Lambda_G8pi=round(math.log10(lam_over_m2), 2),
              log10_eps2={str(k): round(math.log10(v), 2) for k, v in eps2.items()})
say(f"   Mit der Ledger-Leiter ist Lambda/M^2 = 10^{ZAHLEN['log10_Lambda_G8pi']}. Für native Längen a^2/G = 0,2 bis 24 "
    f"(Bestand) wäre eps^2 = 10^{ZAHLEN['log10_eps2']['0.2']} bis 10^{ZAHLEN['log10_eps2']['24.0']}.")
befund("D_unendlich", {"R_eukl": 48, "R_lor": -48},
       "Der native Unendlichkeitstwistor ist vom Einstein-Typ (nicht entartet); das Abbruchkriterium 'nur konforme "
       "Gravitation' greift an dieser Stelle nicht. Nativ ist der Krümmungsradius aber die Nahtlänge (|Lambda| = 12). "
       "Flache Einstein-Gravitation ist die Kontraktion I_pi - eps^2 I_omega; die freie Zahl des Kandidaten ist "
       "genau eps^2, im Ledger die Lambda-Stufe. Dass der Vertex der konformen Theorie auf Einstein-Zuständen den "
       "Einstein-Vertex mit Koeffizient ~ Lambda G gibt, ist Literatur (Maldacena 2011, Adamo-Mason 2013) und hier "
       "nicht nachgerechnet [E für die Geometrie, L für die Reduktion].")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil E  Bilanz und Urteil")
bilanz = [
    ("E", "W1 Struktur", "Der Rahmen hat genau eine freie reine Zahl; Verhältnisse sind Niveauverhältnisse. Der "
     "Zähler <12>^4 ist ein Berezin-Integral mit Koeffizient eins."),
    ("N negativ", "W1 P_norm", "g = 1 an der Nahtskala verfehlt die gemessenen Kopplungen um den Faktor 2,7 bis 4,0."),
    ("N, Hinweis", "W1 Wert", "Die nichtabelschen Kopplungen liegen bei M_Planck,red bei 1/g^2 = 3,9 bis 4,0 "
     "(16 pi = 2/c3); eine ganze Zahl und eine Skala gewählt, U(1) passt nicht; als exakte Gleichung nicht haltbar "
     "(4 bis 6 % in alpha_s und sin^2 theta_W)."),
    ("O", "W1 alpha-Gleichung", "Für den kubischen Term und den Blatt-Term der alpha-Gleichung ist kein Gegenstück "
     "im Rahmen gerechnet; die Stationarität bleibt offen."),
    ("E", "W2 Gewichte", "Zähler <ij>^8 aus dem Serre-Partner O(-6); Hodges-Determinante skaliert richtig."),
    ("E", "W2 Unendlichkeitstwistor", "I_pi - I_omega ist nicht entartet, konstante Krümmung |Lambda| = 12 nativ; "
     "flach nur als Kontraktion mit eps^2."),
    ("L", "W2 Reduktion", "Einstein aus konformer Gravitation mit Koeffizient Lambda G: Literatur, nicht nachgerechnet."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Keiner der beiden Wege ist entschieden, aber beide sind schärfer gestellt. W1: Der Rahmen legt die Kopplung "
          "nicht fest. Er hat genau eine freie Zahl, und die naheliegende native Normierung (P_norm, g = 1 aus "
          "k = 1) ist gegen die Daten widerlegt. Die Daten verlangen für SU(2) und SU(3) an der reduzierten "
          "Planckmasse 1/g^2 nahe 4; das ist ein schwacher Hinweis mit einer gewählten ganzen Zahl, kein Treffer, und "
          "als exakte Gleichung ist er um 4 bis 6 % daneben. "
          "Dass die alpha-Gleichung die Stationarität desselben Funktionals ist, ist nicht gezeigt; das "
          "Erfolgskriterium von W1 ist nicht erfüllt, das Abbruchkriterium auch nicht. W2: Die Gewichte geben den "
          "Zähler <ij>^8 ohne Zusatzannahme, und der native Unendlichkeitstwistor ist vom Einstein-Typ. Damit ist "
          "die Form da. Die Zahl fehlt: Nativ ist die Krümmung von der Ordnung der Nahtlänge, und die flache "
          "Gravitation verlangt die Kontraktion eps^2, die genau die Lambda-Stufe der Leiter ist. Beide Wege enden "
          "bei derselben Art von Frage: Der Rahmen gibt Formen und Verhältnisse, die zwei reinen Zahlen (Kopplung "
          "und eps^2) kommen aus ihm allein nicht. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z164_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z164.txt").write_text("\n".join(OUT) + "\n")
