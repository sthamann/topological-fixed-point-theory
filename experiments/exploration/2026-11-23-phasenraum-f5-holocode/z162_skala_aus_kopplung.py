"""Z162: Vorprüfung der Neuausrichtung 'Skala aus Kopplung' (Diagnose nach den Wellen A bis C).

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Anlass: Die Wellen A bis C (z146 bis z161) haben zwölf Wege geprüft, eine der neun Prämissen aus der einzelnen Naht
herzuleiten; keiner hat eine Prämisse entfernt. Diese Rechnung prüft die Diagnose und die erste neue Strategie.

Diagnose (wird hier geprüft):
  D1  Die nativen Daten sind ein Gaußscher Zustand mit bilinearen Erzeugern. Daraus folgt Kinematik; jede Wechselwirkung
      ist eine Zusatzannahme.
  D2  P_viertel legt die Kompositionskette auf ihren selbstdualen Punkt. Dort gibt es nur eine Skala. Deshalb blieb
      G/R_vac^2 auf allen sieben Wegen eine freie Zahl der Ordnung eins in Gittereinheiten.
  D3  Die gelösten Zahlen (alpha, Wirkungsleiter e^(-A)) wurden nie als Eingabe der Dynamik benutzt.

Bestand: z120 (Gatter U = B^+ E B, U^2 = Tausch), z121 (Kette = Destri-de-Vega-Lichtkegelgitter, Theta = tau/2),
z155 (H_Q = (Q_A^2 + Q_B^2)/2 aus dem Gauß-Gesetz), z159 (geeichte Schrittkette = z121; nativer Zustand = e^(4 pi S^z)/Z),
z160 (elliptische Familie, Masse 8 e^(-pi K/lambda) [L]), Paper 1 (R_vac M = 2 pi e^(1/alpha), A_EW = 27,5259; Status I).

Teil A  Zuordnung Winkel -> Gitterparameter, mit und ohne Tausch (Spinkette, 6 Plätze)
Teil B  Fermionsprache: zwei Vierteltransfers = Tausch mal inverse Uhr; Durchgang = Tausch . Uhr^-1 . Sandwich
Teil C  Ein Magnon in der Durchgangskette: Dirac-Teilchen auf dem Lichtkegelgitter
Teil D  Gaußsches No-Go: bilineare Entwicklung erzeugt keine Kumulante, das Gatter schon
Teil E  Skalen: selbstdualer Punkt gegen schwache Kopplung; Zahlen gegen die Ledger-Leiter (alles gegen Zahlenspielerei
        markiert)
Teil F  Bilanz und Urteil.
"""
import json
import math
import os
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import numpy as np
import sympy as sp
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
TOL = 1e-10
PI = math.pi


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def de(x, n=4):
    return f"{x:.{n}f}".replace(".", ",")


def sci(x, n=3):
    return f"{x:.{n}e}".replace(".", ",")


# ---------------------------------------------------------------------------------------------------------------------
say("Z162: Vorprüfung 'Skala aus Kopplung'")
say("\nTeil A  Zuordnung Winkel -> Gitterparameter (Spinkette mit 6 Plätzen, Yang-R-Matrix R(w) = w + iP)")
L = 6
DIM = 2 ** L
ID = np.eye(DIM)


def perm(i, j):
    M = np.zeros((DIM, DIM))
    for s in range(DIM):
        t = s ^ ((1 << i) | (1 << j)) if ((s >> i) & 1) != ((s >> j) & 1) else s
        M[t, s] = 1
    return M


def gate(u, i, j):
    return (ID + 1j * u * perm(i, j)) / (1 + 1j * u)


def brick(u, swap):
    A, B = ID.astype(complex), ID.astype(complex)
    for i in range(0, L, 2):
        g = gate(u, i, (i + 1) % L)
        A = A @ (perm(i, (i + 1) % L) @ g if swap else g)
    for i in range(1, L, 2):
        g = gate(u, i, (i + 1) % L)
        B = B @ (perm(i, (i + 1) % L) @ g if swap else g)
    return A @ B, B @ A


def tmat(x, th):
    xi = [th * (-1) ** j for j in range(L)]
    bl = [[ID.astype(complex) if a == b else np.zeros((DIM, DIM), complex) for b in range(2)] for a in range(2)]
    for j in range(L):
        w = x - xi[j]
        Eab = [[np.zeros((DIM, DIM)) for _ in range(2)] for _ in range(2)]
        for s in range(DIM):
            b = (s >> j) & 1
            for a in range(2):
                Eab[a][b][(s & ~(1 << j)) | (a << j), s] = 1
        R = [[(w * ID if a == b else 0 * ID) + 1j * Eab[b][a] for b in range(2)] for a in range(2)]
        bl = [[sum(R[a][c] @ bl[c][b] for c in range(2)) for b in range(2)] for a in range(2)]
    return bl[0][0] + bl[1][1]


def prop_res(A, B):
    k = int(np.argmax(np.abs(B)))
    c = A.flat[k] / B.flat[k]
    return float(np.abs(A - c * B).max())


def ddv(th):
    return tmat(-th, th) @ np.linalg.inv(tmat(th, th))


PHIS = (0.4, 1.0, PI / 2, 2.2)
resA1, resA2, resA3 = [], [], []
for phi in PHIS:
    u = math.tan(phi / 2)
    resA1.append(prop_res(brick(u, False)[0], ddv(u / 2)))
    resA2.append(prop_res(brick(u, True)[1], np.linalg.inv(ddv(1 / (2 * u)))))
    if abs(phi - PI / 2) > 1e-9:
        resA3.append(min(prop_res(F, C) for F in brick(u, False)
                         for C in (ddv(1 / (2 * u)), np.linalg.inv(ddv(1 / (2 * u))))))
ZAHLEN["A_residuen"] = {"sandwich": f"{max(resA1):.1e}", "durchgang": f"{max(resA2):.1e}",
                        "falsche_zuordnung_min": de(min(resA3), 3)}
check("A1 [E] Sandwich-Kette (Gatter g(u) = (1 + iuP)/(1 + iu), u = tan(phi/2), ohne Tausch) = t(-Theta) t(Theta)^-1 "
      f"mit Theta = tan(phi/2)/2, für phi = 0,4; 1,0; pi/2; 2,2 (Residuum < 1e-10)", max(resA1) < TOL)
check("A2 [E] Durchgangs-Kette (Gatter Tausch . g(u)) = Inverses desselben Operators in der anderen Schichtfolge mit "
      "Theta = cot(phi/2)/2 (Residuum < 1e-10)", max(resA2) < TOL)
check("A3 [E] Korrektur zu z159 E1: Die Zuordnung Theta = cot(phi/2)/2 passt für die Sandwich-Kette nur bei phi = pi/2 "
      f"(kleinstes Residuum sonst {de(min(resA3), 3)}); z121 (Theta = tau/2) ist richtig", min(resA3) > 1e-2)
us = sp.symbols("u", real=True)
P4 = sp.Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
g4 = (sp.eye(4) + sp.I * us * P4) / (1 + sp.I * us)
dual = sp.simplify(P4 * g4 - g4.subs(us, -1 / us) * (sp.I * us * (1 - sp.I / us) / (1 + sp.I * us)))
check("A4 [E] Dualität: Tausch . g(u) = Phase . g(-1/u) (symbolisch); der Tausch bildet Theta auf 1/(4 Theta) ab",
      dual == sp.zeros(4, 4))
pg2 = sp.simplify((P4 * g4 * g4)[1, 2] * (1 + sp.I * us) ** 2)       # P g^2 proportional zu 1 <=> dieses Element = 0
sols = sorted(str(s) for s in sp.solve(sp.Eq(pg2, 0), us))
check("A5 [E] Tausch . g(u) ist proportional zu g(u)^-1 genau für u^2 = 1, also phi = +-pi/2: P_viertel ist der "
      f"selbstduale Punkt Theta = 1/2 (Lösungen u = {', '.join(sols)})", sols == ["-1", "1"])
befund("A", "Theta = tan(phi/2)/2 (Sandwich), cot(phi/2)/2 (Durchgang)",
       "Eine Familie, zwei Enden. Der Viertelwinkel ist der selbstduale Punkt. Großes Theta (kleine Skala) liegt beim "
       "Sandwich bei phi -> pi, beim Durchgang bei phi -> 0.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil B  Fermionsprache: zwei Nähte A, B mit je zwei Moden (Fockraum 16)")
FD = 16


def ann(j):
    M = np.zeros((FD, FD))
    for s in range(FD):
        if (s >> j) & 1:
            M[s ^ (1 << j), s] = (-1) ** bin(s & ((1 << j) - 1)).count("1")
    return M


a = [ann(j) for j in range(4)]                         # A1, A2, B1, B2
ad = [x.T for x in a]
num = [ad[j] @ a[j] for j in range(4)]
NA, NB = num[0] + num[1], num[2] + num[3]
NT = NA + NB
JX = 0.5 * sum(ad[k] @ a[k + 2] + ad[k + 2] @ a[k] for k in range(2))
BS = expm(-1j * PI / 2 * JX)
B2 = BS @ BS
RCL = expm(1j * PI / 2 * NT)                           # Uhr der Gesamtladung, R = i^N
W = B2 @ RCL
check("B1 [E] Zwei Vierteltransfers: B^2 a_k B^-2 = i b_k; W = B^2 R (R = i^N) erfüllt W a W^-1 = b, W b W^-1 = a, "
      "W^2 = 1. Also B^2 = Tausch . R^-1",
      max(np.abs(B2 @ a[k] @ B2.conj().T - 1j * a[k + 2]).max() for k in range(2)) < TOL
      and max(np.abs(W @ a[k] @ W.conj().T - a[k + 2]).max() for k in range(2)) < TOL
      and max(np.abs(W @ a[k + 2] @ W.conj().T - a[k]).max() for k in range(2)) < TOL
      and np.abs(W @ W - np.eye(FD)).max() < TOL)


def EQ(phi):
    return expm(-1j * phi * (NA @ NA + NB @ NB) / 2)


def usand(phi):
    return BS.conj().T @ EQ(phi) @ BS


def upass(phi):
    return BS @ EQ(phi) @ BS


check("B2 [E] Durchgang B E B = Tausch . R^-1 . (B^+ E B); die Uhr ist zentral und vertauscht mit dem Gatter",
      all(np.abs(upass(p) - W @ RCL.conj().T @ usand(p)).max() < TOL
          and np.abs(usand(p) @ RCL - RCL @ usand(p)).max() < TOL for p in PHIS))
NEUT = [0b1001, 0b0110]                                # A1 B2, A2 B1: je Naht ein Teilchen
okB3 = True
for p in PHIS:
    Us = usand(p)[np.ix_(NEUT, NEUT)]
    ev = sorted(np.round(np.mod(-np.angle(np.linalg.eigvals(Us)), 2 * PI), 9))
    okB3 &= np.abs(Us @ Us.conj().T - np.eye(2)).max() < TOL
    okB3 &= min(abs(abs(ev[1] - ev[0]) - q) for q in (p, 2 * PI - p)) < 1e-8
check("B3 [E] Im neutralen Sektor (je Naht ein Teilchen) ist das Sandwich-Gatter unitär und hat die relative Phase "
      "e^(-i phi) zwischen den zwei Tauschsektoren: g(u) mit u = tan(phi/2), wie in z159", okB3)
befund("B", "B^2 = Tausch . R^-1",
       "Die Durchgangs-Komposition unterscheidet sich vom Sandwich um den vollen Transfer; der ist nativ (zwei "
       "Viertelschritte) und trägt die inverse Uhr. Die Wahl B^+ E B gegen B E B ist die Prämisse P_sym aus z120.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Ein Magnon über dem Referenzzustand, Durchgangskette (Ring mit 24 Plätzen)")
LM = 24


def magnon_step(eps):
    u = math.tan(eps / 2)
    t, r = 1 / (1 + 1j * u), 1j * u / (1 + 1j * u)

    def layer(off):
        M = np.zeros((LM, LM), complex)
        for i in range(off, LM + off, 2):
            i0, i1 = i % LM, (i + 1) % LM
            M[i1, i0] = t
            M[i0, i0] = r
            M[i0, i1] = t
            M[i1, i1] = r
        return M
    return layer(1) @ layer(0)


devC = []
for eps in (0.05, 0.2, 1.1, PI / 2):
    F = magnon_step(eps)
    ev = np.sort(np.mod(-np.angle(np.linalg.eigvals(F)) + 1e-9, 2 * PI))
    K = 2 * PI * np.arange(LM // 2) / (LM // 2)
    kern = 2 * np.arccos(np.cos(eps / 2) * np.cos(K / 2))
    pr = np.sort(np.mod(np.concatenate([eps + kern, eps - kern]) + 1e-9, 2 * PI))
    devC.append(float(np.abs(ev - pr).max()))
check("C1 [E] Quasienergie eines Magnons: omega = eps +- 2 arccos(cos(eps/2) cos(K/2)), K = Impuls der Zweinahtzelle "
      "(eps = 0,05; 0,2; 1,1; pi/2; Abweichung < 1e-10)", max(devC) < TOL)
eps0 = 0.05
gap = 2 * math.acos(math.cos(eps0 / 2))
check("C2 [E] Für kleines eps ist das die Dirac-Dispersion omega - eps = +-sqrt(eps^2 + cos(eps/2) K^2 + ...): Masse "
      "m a = eps je Schritt, Lichtgeschwindigkeit eine Zelle je Schritt bis auf O(eps^2)",
      abs(gap - eps0) < 1e-12 and abs(2 * math.acos(math.cos(eps0 / 2) * math.cos(0.01 / 2))
                                      - math.sqrt(eps0 ** 2 + math.cos(eps0 / 2) * 0.01 ** 2)) < 1e-6)
n_mag = 1 / (1 + math.exp(4 * PI))
ZAHLEN["C_magnondichte_nativ"] = sci(n_mag)
check("C3 [E] Der native Zustand e^(4 pi S^z)/Z (z159) ist ein verdünntes Magnongas: Dichte 1/(1 + e^(4 pi)) = "
      f"{sci(n_mag)} je Naht", abs(n_mag - 3.487e-6) < 1e-8)
befund("C", "Dirac-Teilchen mit m a = eps",
       "Über dem Referenzzustand, bei dem der native Zustand liegt, ist die Anregung der Durchgangskette exakt ein "
       "relativistisches Teilchen auf dem Lichtkegelgitter. Die Sandwich-Kette ist dieselbe Formel mit eps = pi - phi: "
       "bei kleinem phi liegt die 'Masse' an der Gittergrenze, die Dispersion ist quadratisch (z159 B10).")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil D  Gaußsches No-Go (Fockraum 16, Kumulante kappa = <n_A1 n_B1> - <n_A1><n_B1> + |<a_1^+ b_1>|^2)")


def kappa(rho):
    ex = lambda X: np.trace(rho @ X)
    return float(np.real(ex(num[0] @ num[2]) - ex(num[0]) * ex(num[2]) + abs(ex(ad[0] @ a[2])) ** 2))


def gauss(levels):
    K = sum(levels[j] * num[j] for j in range(4))
    r = expm(-K)
    return r / np.trace(r)


rho_gleich = gauss([PI, 3 * PI, PI, 3 * PI])           # zwei gleiche native Nähte
rho_ungl = gauss([0.3, 1.1, 2.0, 0.7])                 # ungleiche Gaußsche Nähte (Probe)
rng = np.random.default_rng(162)
Hb = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
Hb = Hb + Hb.conj().T
Ubil = expm(-1j * sum(Hb[i, j] * ad[i] @ a[j] for i in range(4) for j in range(4)))
k0 = kappa(rho_ungl)
k_bil = kappa(Ubil @ rho_ungl @ Ubil.conj().T)
k_gate = [kappa(usand(p) @ rho_ungl @ usand(p).conj().T) for p in PHIS]
k_nat = [float(np.abs(usand(p) @ rho_gleich @ usand(p).conj().T - rho_gleich).max()) for p in PHIS]
ZAHLEN["D_kumulante_nach_gatter"] = [de(k, 5) for k in k_gate]
check("D1 [E] Gaußscher Zustand: kappa = 0; nach beliebiger bilinearer (quasifreier) Entwicklung weiter kappa = 0",
      abs(k0) < TOL and abs(k_bil) < TOL)
check("D2 [E] Das Gatter erzeugt aus ungleichen Gaußschen Nähten kappa != 0 (phi = 0,4; 1,0; pi/2; 2,2: "
      f"{', '.join(de(k, 5) for k in k_gate)}): Es ist kein quasifreier Schritt",
      min(abs(k) for k in k_gate) > 1e-4)
check("D3 [E] Auf dem Produkt gleicher nativer Nähte wirkt das Gatter für jedes phi trivial (z159): Der native "
      "Zustand kann den Winkel nicht messen", max(k_nat) < TOL)
befund("D", "Wechselwirkung liegt nicht in den nativen Daten",
       "Zustand, modulare Erzeuger, Uhr und so(3,2) sind bilinear oder Exponentiale von Bilinearen. Der erste Ort, an "
       "dem eine Kopplung eingeht, ist das Gatter; ihr Wert kann aus dem nativen Zustand nicht abgelesen werden.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil E  Skalen (Spinonmasse des Lichtkegelgitters m a = 4 e^(-pi Theta) für großes Theta [L, Destri-de Vega];"
    " exakt am Loch bei Rapidität 0: 4 arctan(e^(-pi Theta)) [L])")
# Kontrolle der Dichteformel, auf der die Massenformel ruht: rho + K*rho = Quelle, rho = 1/(2 cosh(pi x))
xs = np.linspace(-40, 40, 16001)
hx = xs[1] - xs[0]
rho = 1 / (2 * np.cosh(PI * xs))
x_test = np.array([0.0, 0.3, 1.0, 2.5])
lhs = [1 / (2 * math.cosh(PI * x)) + float(np.sum(rho / (PI * (1 + (x - xs) ** 2))) * hx) for x in x_test]
rhs = [(1 / (2 * PI)) / (x * x + 0.25) for x in x_test]
check("E1 [N] Bethe-Dichte des antiferromagnetischen Vakuums: rho(x) = 1/(2 cosh(pi x)) löst rho + K*rho = "
      "(1/2pi)/(x^2 + 1/4) mit K = (1/pi)/(1 + x^2) (vier Stützstellen, Abweichung < 1e-6)",
      max(abs(l - r) for l, r in zip(lhs, rhs)) < 1e-6)
ma_dual = 4 * math.atan(math.exp(-PI / 2))
ZAHLEN["E_ma_selbstdual"] = de(ma_dual, 4)
check(f"E2 [L] Am selbstdualen Punkt Theta = 1/2 ist m a = {de(ma_dual, 4)} (Kontinuumsformel 4 e^(-pi/2) = "
      f"{de(4 * math.exp(-PI / 2), 4)}): eine Skala, keine Trennung", 0.5 < ma_dual < 1.0)
ALPHA_INV = 137.03599921684                            # Ledger (gelöster Block, nur als Eingabe)
al = 1 / ALPHA_INV
eps_star = 2 * math.atan(PI * al / 2)
T_star = eps_star / (4 * PI * al)
ratio_pa = math.exp(-(PI / 2) / math.tan(PI * al / 2)) / math.exp(-ALPHA_INV)
ZAHLEN["E_eps_stern"] = de(eps_star, 8)
ZAHLEN["E_T_stern"] = de(T_star, 6)
ZAHLEN["E_verhaeltnis_bei_eps_pi_alpha"] = de(ratio_pa, 5)
check("E3 [N] Durchgangskette mit Theta = cot(eps/2)/2: pi Theta = 1/alpha verlangt eps = 2 arctan(pi alpha/2) = "
      f"{de(eps_star, 8)}; mit g^2 = 4 pi alpha ist eps/g^2 = {de(T_star, 6)} (ein Viertel bis auf 4e-5)",
      abs(T_star - 0.25) < 1e-4)
check(f"E4 [N] Mit exakt eps = pi alpha: m a = 4 e^(-1/alpha) mal {de(ratio_pa, 5)} (= e^(pi^2 alpha/12))",
      abs(ratio_pa - math.exp(PI ** 2 * al / 12)) < 1e-6)
tab = []
for nm, c in (("alpha", 1.0), ("2 alpha", 2.0), ("pi alpha", PI), ("4 alpha", 4.0), ("2 pi alpha", 2 * PI),
              ("4 pi alpha", 4 * PI)):
    A = (PI / 2) / math.tan(c * al / 2)
    tab.append((nm, A * al))
    say(f"      eps = {nm:<10}: Exponent pi Theta = {de(A, 3)} = {de(A * al, 4)}/alpha")
ZAHLEN["E_exponent_je_wahl"] = {nm: de(x, 4) for nm, x in tab}
check("E5 [N] Gegenprobe gegen Zahlenspielerei: Von sechs naheliegenden Wahlen eps = c alpha trifft genau eine "
      "(c = pi) den Exponenten 1/alpha; der Treffer ist die Wahl einer Zahl, kein Beweis",
      sum(abs(x - 1) < 1e-3 for _, x in tab) == 1)
d_log = (PI / 2) / math.sin(eps_star / 2) ** 2 * eps_star / 2
ZAHLEN["E_empfindlichkeit"] = de(d_log, 2)
check(f"E6 [N] Empfindlichkeit: d ln(m a)/d ln(eps) = {de(d_log, 2)}; ein Prozent in eps ändert die Skala um den "
      f"Faktor {de(math.exp(d_log / 100), 2)}. Prüfbar ist die Lesart nur, wenn eps hergeleitet wird",
      abs(d_log - ALPHA_INV) < 0.05)
a_M_1 = 8 * PI * ratio_pa                              # R_vac = 1/m
a_M_2 = 4 * ratio_pa                                   # R_vac = 2 pi/m (Compton-Wellenlänge)
a2G_1, a2G_2 = a_M_1 ** 2 * 8 * PI, a_M_2 ** 2 * 8 * PI
ZAHLEN["E_a2_durch_G"] = {"R_vac = 1/m": de(a2G_1, 0), "R_vac = 2 pi/m": de(a2G_2, 0)}
check("E7 [C: P_eps, P_see] Gegen das Ledger R_vac M = 2 pi e^(1/alpha): Gitterabstand a M = "
      f"{de(a_M_1, 2)} (R_vac = 1/m) oder {de(a_M_2, 3)} (R_vac = 2 pi/m), also a^2/G = {de(a2G_1, 0)} oder "
      f"{de(a2G_2, 0)}. Die sieben nativen Wege geben 0,2 bis 24: Der Vorfaktor ist nicht getroffen; aus der freien "
      "Zahl 10^123 würde eine Zahl zwischen 10 und 10^5", 300 < a2G_2 < 500 and 1.5e4 < a2G_1 < 1.7e4)
A_EW = 27.5259
N_EW = 2 * ALPHA_INV / A_EW
ZAHLEN["E_gruppenrang_fuer_EW"] = de(N_EW, 3)
check("E8 [N negativ] Leiter: Für SU(N)-Ketten ist der Exponent 2 pi Theta/N [L]. Bei festem Theta gäbe N = 2 den "
      f"Exponenten 1/alpha; A_EW = 27,5259 verlangte N = {de(N_EW, 3)}, keine ganze Zahl und kein natives N",
      abs(N_EW - round(N_EW)) > 0.03)
bits = (2 * ALPHA_INV + math.log(32 * PI ** 3)) / math.log(2)
ZAHLEN["E_zaehlweg_bits"] = de(bits, 1)
check(f"E9 [N] Zählweg zum Vergleich: R_vac^2/G = 32 pi^3 e^(2/alpha) ist die Dimension eines Raums von {de(bits, 1)} "
      "Qubits; als Anzahl von Nähten (P_N aus z158) wäre es eine freie Zahl, als Zustandszahl eine endliche Aufgabe",
      404 < bits < 407)
befund("E", "selbstdual: eine Skala; schwache Kopplung: Hierarchie",
       "P_viertel und eine Skalentrennung schließen sich in der Kompositionskette aus. Die Hierarchie des Ledgers hat "
       "die Form der Spinonmasse der Durchgangskette; der Exponent verlangt eps = pi alpha (P_eps), der Zustand das "
       "antiferromagnetische Vakuum (P_see), der Vorfaktor ist offen.")

# ---------------------------------------------------------------------------------------------------------------------
say("\nTeil F  Bilanz")
bilanz = [
    ("E", "Familie", "Sandwich und Durchgang sind eine Familie: Theta = tan(phi/2)/2 und cot(phi/2)/2. Sie "
     "unterscheiden sich um B^2 = Tausch . R^-1, zwei native Viertelschritte. z159 E1 ist korrigiert, z121 bestätigt."),
    ("E", "P_viertel", "Der Viertelwinkel ist genau der selbstduale Punkt (Tausch . Gatter = Gatter^-1). Das ist die "
     "dritte Fassung derselben Selbstdualität nach z155 (Villain) und z157 (komplexe Multiplikation)."),
    ("E/L", "Diagnose D2", "Am selbstdualen Punkt hat die Kette eine einzige Skala (m a = 0,82). G/R_vac^2 kann dort "
     "nativ nur von der Ordnung eins sein; das erklärt, warum sieben Wege die Zahl frei ließen."),
    ("E", "Diagnose D1", "Bilineare Entwicklung erzeugt keine Kumulante; das Gatter schon; der native Zustand sieht "
     "den Winkel nicht. Eine Kopplung folgt aus der einzelnen Naht nicht."),
    ("E", "Dirac", "In der Durchgangskette ist ein Magnon ein Dirac-Teilchen mit m a = eps; der native Zustand ist "
     "ein Magnongas der Dichte 3,5e-6 je Naht."),
    ("C", "Hierarchie", "Unter P_durch (Komposition B E B), P_eps (eps = pi alpha = g^2/4) und P_see "
     "(antiferromagnetisches Vakuum) ist die Vakuumskala m a = 4 e^(-1/alpha) mal 1,006 [L]. Der Vorfaktor gegen das "
     "Ledger verlangt a^2/G zwischen 4e2 und 1,6e4, nativ sind es 0,2 bis 24."),
    ("N negativ", "Leiter", "Die zweite Ledgerskala (A_EW) folgt aus keinem ganzzahligen Rang; die Leiter 1:5:10 ist "
     "damit nicht erklärt."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Die Diagnose trägt, die neue Strategie ist konsistent, aber kein Treffer. (1) Die Wellen A bis C konnten "
          "keine Prämisse entfernen, weil die nativen Daten Gaußsch sind: Wechselwirkung beginnt am Gatter, und der "
          "native Zustand ist gegen dessen Winkel blind. (2) P_viertel ist der selbstduale Punkt der "
          "Kompositionsfamilie. Dort gibt es eine einzige Skala; die freie Zahl G/R_vac^2 kann aus dem Kandidaten "
          "mit P_viertel grundsätzlich nicht als 10^-123 herauskommen, sondern nur über eine Anzahl (P_N). (3) Die "
          "Familie hat ein zweites Ende: die Durchgangs-Komposition B E B mit schwacher Phase eps. Sie ist vom "
          "Sandwich um zwei native Viertelschritte verschieden (Tausch mal inverse Uhr), ihre Anregung über dem "
          "Referenzzustand ist ein Dirac-Teilchen mit m a = eps, und ihre Vakuumskala ist exponentiell klein, "
          "m a = 4 e^(-(pi/2) cot(eps/2)). Das ist die Form der Ledger-Hierarchie e^(-1/alpha). (4) Der Preis: "
          "P_viertel entfällt oder gilt nur am Gitterende; neu sind P_durch, P_eps (eps = pi alpha, eine Zahl auf "
          "eine Zahl gelegt, Empfindlichkeit 137) und P_see. Der Vorfaktor ist um den Faktor 20 bis 8e4 nicht "
          "getroffen, die zweite Skala der Leiter gar nicht. (5) Für den Kandidaten heißt das: P_G lässt sich von "
          "einer freien Zahl 10^123 auf die Frage 'warum eps = g^2/4' und einen Vorfaktor der Ordnung 10^2 bis 10^4 "
          "umstellen. Das ist eine Umformulierung mit prüfbarem nächsten Schritt (eps aus der geeichten Uhr mit dem "
          "gelösten alpha), keine Herleitung. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z162_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z162.txt").write_text("\n".join(OUT) + "\n")
