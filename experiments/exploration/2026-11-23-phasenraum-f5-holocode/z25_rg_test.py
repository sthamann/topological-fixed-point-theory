"""Z25: RG-Test der H6-Lesart "Eichkopplungen an der Gitterskala induziert und gleich" (y1, z15, z21 S4/6.d).

H6 (y1): Die Eichkopplungen werden an der Gitterskala mu = 1/a durch Fermionschleifen induziert. Weil jede
Generation für SU(3), SU(2) und U(1)_Y (GUT-normiert) dieselbe Dynkin-Summe hat, sind sie dort gleich
(sin^2 theta_W = 3/8). z15: a/l_P liegt im Fenster 1,65 ... 14,3, im H6-Schema bei 10,8.

Gerechnet:
 0. Eingaben: PDG-Kopplungen bei M_Z wie v246/v249/v266, Planckeinheit aus z2 (l_P^2 = G), Fenster für 1/a,
    Higgsdublett als Import aus TFPT (PS.NCG.FLUCT.01).
 1. Kalibrierung: allgemeine Zwei-Schleifen-Formel (Machacek-Vaughn) reproduziert die SM-Matrix b_ij exakt,
    v246 (Spreizung 8,8 bei 2e16 in einer Schleife, 2-Schleifen-RK4), v266 (M_PS, M_GUT) und den
    2-Schleifen-Hybrid aus experiments/gauge-unification.
 2. Frage 1: reines SM, zwei Schleifen mit y_t, bis 1/a im ganzen Fenster.
 3. Frage 2: Zwischeninhalt aus der E8-Hülle {1, 10, 16, 45} (v247: sechs reelle 10, vier komplexe 16, eine
    reelle 45; kein 126) für Gleichheit genau an 1/a:
      A  SM -> PS mit vollständigen PS-Multipletts (13440 diskrete Inhalte, 1 und 2 Schleifen; dazu LP mit
         frei gewählten Schwellen oberhalb M_PS),
      B  SM-Eichgruppe bis 1/a (H6 ohne neue Eichbosonen) mit E8-Fragmenten auf freien Schwellen (LP),
      C  PS mit gespaltenen Multipletts (leichte SM-Fragmente unter M_PS; LP-Relaxation, ganzzahlige
         Aufzählung, 2-Schleifen-Bestätigung), Protonzerfall d = 6 wie v266 und durch skalare Tripletts.
 4. Frage 3: SO(10) vereinigt unterhalb von 1/a und läuft bis 1/a (mit den H6-Eltern (10,1)+(1,6)).
    Einheitliche Kopplung an 1/a gegen den y1-induzierten Wert (Formel aus y1 nachgerechnet und für den
    H6-Inhalt neu ausgewertet); Selbstkonsistenz von a/l_P mit der z15-Kante.
 5. Frage 4: welche Lesart fällt, welche überlebt, was wäre die Vorhersage.

Konventionen: alpha_i^-1 GUT-normiert, d alpha_i^-1/dt = -b_i/2pi - sum_j b_ij alpha_j/8pi^2 + c_i y_t^2/32pi^3,
t = ln mu. Schwellen einstufig (Stufenfunktion, Matching auf Baumniveau). y_t in führender Ordnung (eine
Schleife) unterhalb M_PS; oberhalb ist sein Zwei-Schleifen-Eichbeitrag < 0,01 in alpha^-1 und weggelassen.
"""
import itertools
import json
import time
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve, linprog

T0 = time.time()
HERE = Path(__file__).parent
OUT, CHECKS, BEFUNDE, ZAHLEN = [], {}, [], {}


def say(s=""):
    print(s)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(s):
    BEFUNDE.append(s)
    say(f"   [Befund] {s}")


def e(x, n=2):
    return f"{x:.{n}e}"


# ============================================================================ Gruppentheorie
REP = {  # (Dynkin-Index T, Casimir C2)
    "SU2": {1: (Fr(0), Fr(0)), 2: (Fr(1, 2), Fr(3, 4)), 3: (Fr(2), Fr(2))},
    "SU3": {1: (Fr(0), Fr(0)), 3: (Fr(1, 2), Fr(4, 3)), 8: (Fr(3), Fr(3))},
    "SU4": {1: (Fr(0), Fr(0)), 4: (Fr(1, 2), Fr(15, 8)), 6: (Fr(1), Fr(5, 2)), 15: (Fr(4), Fr(4))},
    "SO10": {1: (Fr(0), Fr(0)), 10: (Fr(1), Fr(9, 2)), 16: (Fr(2), Fr(45, 8)), 45: (Fr(8), Fr(8))},
}
CA = {"U1": Fr(0), "SU2": Fr(2), "SU3": Fr(3), "SU4": Fr(4), "SO10": Fr(8)}


def betas(groups, fields):
    """Ein- und Zwei-Schleifen-Koeffizienten (Machacek-Vaughn) für ein Produkt einfacher Faktoren und U(1)
    (GUT-normiert, Y in Standardnormierung). fields: (Darstellungen, Art, Anzahl); Art W = Weyl-Fermion,
    C = komplexer Skalar, R = reeller Skalar. b_ij = -34/3 C_A^2 d_ij + sum_W S_i [2 C_j + 10/3 C_A d_ij]
    + sum_C S_i [4 C_j + 2/3 C_A d_ij] (reell: halb)."""
    n = len(groups)
    b = [Fr(-11, 3) * CA[g] for g in groups]
    B = [[Fr(-34, 3) * CA[g] ** 2 if i == j else Fr(0) for j in range(n)] for i, g in enumerate(groups)]
    for reps, kind, mult in fields:
        TC, dims = [], []
        for g, r in zip(groups, reps):
            if g == "U1":
                y = Fr(r)
                TC.append((Fr(3, 5) * y * y, Fr(3, 5) * y * y))
                dims.append(1)
            else:
                TC.append(REP[g][r])
                dims.append(r)
        dtot = int(np.prod(dims))
        w1 = {"W": Fr(2, 3), "C": Fr(1, 3), "R": Fr(1, 6)}[kind]
        half = Fr(1, 2) if kind == "R" else Fr(1)
        for i in range(n):
            S_i = TC[i][0] * dtot / dims[i]
            b[i] += w1 * S_i * mult
            for j in range(n):
                if kind == "W":
                    term = 2 * S_i * TC[j][1] + (Fr(10, 3) * CA[groups[i]] * S_i if i == j else 0)
                else:
                    term = 4 * S_i * TC[j][1] + (Fr(2, 3) * CA[groups[i]] * S_i if i == j else 0)
                B[i][j] += half * term * mult
    return b, B


def fl(bB):
    b, B = bB
    return np.array([float(x) for x in b]), np.array([[float(x) for x in r] for r in B])


SMG, PSG = ["U1", "SU2", "SU3"], ["SU4", "SU2", "SU2"]
GEN = [((Fr(1, 6), 2, 3), "W"), ((Fr(-2, 3), 1, 3), "W"), ((Fr(1, 3), 1, 3), "W"),
       ((Fr(-1, 2), 2, 1), "W"), ((Fr(1), 1, 1), "W"), ((Fr(0), 1, 1), "W")]
SM_FERM = [(r, k, 3) for r, k in GEN]
HIGGS = [((Fr(1, 2), 2, 1), "C", 1)]
PS_FERM = [((4, 2, 1), "W", 3), ((4, 1, 2), "W", 3)]
# PS-Komponenten der E8-Hülle (v247): sechs reelle 10 -> (1,2,2)+(6,1,1); vier komplexe 16 -> (4,2,1)+(4b,1,2);
# eine reelle 45 -> (15,1,1)+(1,3,1)+(1,1,3)+(6,2,2). Anzahl je Komponente, Untergrenze (Brechung, Higgs).
PS_COMP = {"bd_r": ((1, 2, 2), "R", 6, 1), "t6_r": ((6, 1, 1), "R", 6, 0), "s412_c": ((4, 1, 2), "C", 4, 1),
           "s421_c": ((4, 2, 1), "C", 4, 0), "a15_r": ((15, 1, 1), "R", 1, 0), "a31_r": ((1, 3, 1), "R", 1, 0),
           "a13_r": ((1, 1, 3), "R", 1, 0), "a622_r": ((6, 2, 2), "R", 1, 0)}
PAR = list(PS_COMP)
# SM-Fragmente der E8-Hülle mit Höchstzahl (45: O8, T3, XY, Q, U, E; 10: H, D; 16: Q, U, E, D, L(=H))
FRAG = {"O8_r": ((0, 1, 8), "R", 1), "T3_r": ((0, 3, 1), "R", 1), "XY_c": ((Fr(-5, 6), 2, 3), "C", 1),
        "Q_c": ((Fr(1, 6), 2, 3), "C", 5), "U_c": ((Fr(-2, 3), 1, 3), "C", 5), "E_c": ((1, 1, 1), "C", 5),
        "H_c": ((Fr(1, 2), 2, 1), "C", 9), "D_c": ((Fr(-1, 3), 1, 3), "C", 10)}
FRAG_PARENT = {"bd_r": ["H_c"], "t6_r": ["D_c"], "s412_c": ["D_c", "U_c", "E_c"], "s421_c": ["Q_c", "H_c"],
               "a15_r": ["O8_r", "U_c"], "a31_r": ["T3_r"], "a13_r": ["E_c"], "a622_r": ["XY_c", "Q_c"]}
GEFAEHRLICH = ("D_c", "XY_c")      # Farbtripletts mit Diquark- und Leptoquark-Kopplung (d = 6 über Skalare)

bSM, BSM = fl(betas(SMG, SM_FERM + HIGGS))
bSM0, BSM0 = fl(betas(SMG, SM_FERM))
bPSF, BPSF = fl(betas(PSG, PS_FERM))
DB_PS = {k: fl(betas(PSG, PS_FERM + [(r, kind, 1)])) for k, (r, kind, _, _) in PS_COMP.items()}
DB_PS = {k: (v[0] - bPSF, v[1] - BPSF) for k, v in DB_PS.items()}
DB_FR = {k: fl(betas(SMG, SM_FERM + HIGGS + [(r, kind, 1)])) for k, (r, kind, _) in FRAG.items()}
DB_FR = {k: (v[0] - bSM, v[1] - BSM) for k, v in DB_FR.items()}
CT = np.array([17 / 10, 3 / 2, 2.0])
MM = np.array([[0, 0, 1], [0, 1, 0], [5 / 3, 0, -2 / 3]])    # (a4, a2L, a2R) aus (a1, a2, a3) bei M_PS

# ============================================================================ Teil 0: Eingaben
say("Z25: RG-Test von H6 — Gleichheit der Eichkopplungen an der Gitterskala 1/a")
say("\nTeil 0: Eingaben und Konventionen")
MZ = 91.1876
AINV = np.array([59.01, 29.59, 8.47])           # v246/v249/v266, GUT-normiert, PDG
DA3 = 0.065                                      # alpha_s = 0,1179(9) -> alpha_3^-1 +- 0,065
MT, YT_MT = 172.69, 0.9369                       # y_t(M_t) NNLO (Buttazzo et al. 2013), Eingabe
MPL = 1.2209e19                                  # 1/l_P mit l_P^2 = G (z2), wie z15 inv_lP
MBAR = 2.435323203e18
MS = (1 / (8 * np.pi)) ** 3.5 * MBAR             # Skalaronskala c3^(7/2) Mbar (v249/v266)
TAU_SK, TAU_HK, BAND = 2.4e34, 1.4e35, 3.0
ALP = {"untere Grenze 1,65": 1.65, "robust 3,32": 3.32, "H6-Schema 10,79": 10.79, "obere Grenze 14,27": 14.27}
MU_A = {k: MPL / v for k, v in ALP.items()}
WIN = (MPL / 14.27, MPL / 1.65)
say(f"   alpha_i^-1(M_Z) = {AINV.tolist()} (GUT-normiert, PDG; wie v246/v249/v266), alpha_s-Fehler -> +-{DA3} in alpha_3^-1.")
say(f"   y_t(M_t = {MT}) = {YT_MT} (Eingabe, führende Ordnung im Lauf). Skalaronskala M_s = c3^(7/2) Mbar = {e(MS)} GeV.")
say("   Planckeinheit: z2 setzt unter D eta = 1/4 je Fläche als Definition, 8 pi G = 2 pi/eta mit eta = kappa/a^2, also")
say("   l_P^2 = G = a^2/(4 kappa): nichtreduzierte Plancklänge, M_Pl = 1/l_P = 1,2209e19 GeV (z15 benutzt denselben Wert).")
for k, v in MU_A.items():
    say(f"   a/l_P = {ALP[k]:5.2f} ({k}): 1/a = {e(v)} GeV")
say(f"   Gegenprobe reduzierte Planckmasse (nicht gemeint): 1/a wäre um 1/sqrt(8 pi) = {1 / np.sqrt(8 * np.pi):.3f} kleiner.")
say("   Higgs: H6 hat kein EW-Higgsfeld (Phi ist SO(16)-Singulett mit Z4-Ladung 2). Das SM-Dublett (1,2,1/2) wird als")
say("   Eingabe aus TFPT importiert (PS.NCG.FLUCT.01: innere Fluktuationen -> genau ein Dublett, kein 126).")
ZAHLEN["eingaben"] = {"alpha_inv_MZ": AINV.tolist(), "M_Pl_GeV": MPL, "M_s_GeV": MS, "yt_Mt": YT_MT,
                      "inv_a_GeV": {k: v for k, v in MU_A.items()}, "fenster_inv_a_GeV": list(WIN)}

# ============================================================================ RG-Werkzeuge


def rhs(t, y, b, B, yuk):
    n = len(b)
    a = y[:n]
    al = 1 / a
    da = -b / (2 * np.pi) - B @ al / (8 * np.pi ** 2)
    if not yuk:
        return da
    yt = y[n]
    g2 = 4 * np.pi * al
    da = da + CT * yt ** 2 / (32 * np.pi ** 3)
    dyt = yt / (16 * np.pi ** 2) * (4.5 * yt ** 2 - 17 / 20 * g2[0] - 9 / 4 * g2[1] - 8 * g2[2])
    return np.concatenate([da, [dyt]])


def run(y0, t0, t1, b, B, yuk=False):
    if abs(t1 - t0) < 1e-14:
        return np.array(y0, float)
    s = solve_ivp(rhs, (t0, t1), np.array(y0, float), args=(b, B, yuk), method="DOP853", rtol=1e-11, atol=1e-12)
    return s.y[:, -1]


YT_MZ = YT_MT
for _ in range(8):
    YT_MZ *= YT_MT / run(np.r_[AINV, YT_MZ], 0, np.log(MT / MZ), bSM, BSM, True)[3]
Y0 = np.r_[AINV, YT_MZ]
say(f"   y_t(M_Z) = {YT_MZ:.4f} (aus y_t(M_t) zurückgelaufen).")


def sm_at(mu, loops=2, b=bSM, B=BSM, ain=AINV):
    T = np.log(mu / MZ)
    if loops == 1:
        return ain - b / (2 * np.pi) * T
    return run(np.r_[ain, YT_MZ], 0, T, b, B, True)[:3]


def sin2(a):
    return 0.6 / (a[0] / a[1] + 0.6)


def tau_p(M, ai):
    """d = 6 über Eichbosonen wie v266: 1e36 (M/1e16)^4 (alpha^-1/40)^2 Jahre, Band x3."""
    return 1e36 * (M / 1e16) ** 4 * (ai / 40) ** 2


# ============================================================================ Teil 1: Kalibrierung
say("\nTeil 1: Kalibrierung")
bF, BF = betas(SMG, SM_FERM + HIGGS)
B2_V246 = [[Fr(199, 50), Fr(27, 10), Fr(44, 5)], [Fr(9, 10), Fr(35, 6), Fr(12)], [Fr(11, 10), Fr(9, 2), Fr(-26)]]
check("Allgemeine Formel reproduziert SM b_i = (41/10, -19/6, -7) und die SM-Matrix b_ij von v246 exakt (Brüche)",
      bF == [Fr(41, 10), Fr(-19, 6), Fr(-7)] and BF == B2_V246)
b_min = fl(betas(PSG, PS_FERM + [((1, 2, 2), "C", 1), ((4, 1, 2), "C", 1)]))[0]
b_15c = fl(betas(PSG, PS_FERM + [((1, 2, 2), "C", 1), ((4, 1, 2), "C", 1), ((15, 1, 1), "C", 1)]))[0]
b_15r = fl(betas(PSG, PS_FERM + [((1, 2, 2), "C", 1), ((4, 1, 2), "C", 1), ((15, 1, 1), "R", 1)]))[0]
b_so10 = fl(betas(["SO10"], [((16,), "W", 3)]))[0][0]
check("PS-Einschleifen-b von v266 reproduziert: minimal (-31/3, -3, -7/3), +(15,1,1) (-9, -3, -7/3); SO(10)+3x16: -76/3",
      np.allclose(b_min, [-31 / 3, -3, -7 / 3]) and np.allclose(b_15c, [-9, -3, -7 / 3]) and abs(b_so10 + 76 / 3) < 1e-12)
say(f"   Hinweis: v266 zählt (15,1,1) komplex (Delta b4 = 4/3, gleich zwei reellen Adjungierten). Die E8-Hülle hat nur eine")
say(f"   reelle 45 (v247), also Delta b4 = 2/3: b = {np.round(b_15r, 4).tolist()}.")
spread_2e16 = float(np.ptp(sm_at(2e16, 1)))


def v246_rk4(t_end, dt=0.01):
    def deriv(ai):
        return -bSM / (2 * np.pi) - (BSM @ (1 / ai)) / (8 * np.pi ** 2)
    t, ai, best, bmu = 0.0, AINV.copy(), 1e9, MZ
    while t < t_end:
        k1 = deriv(ai); k2 = deriv(ai + 0.5 * dt * k1); k3 = deriv(ai + 0.5 * dt * k2); k4 = deriv(ai + dt * k3)
        ai = ai + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
        sp = float(ai.max() - ai.min())
        if sp < best:
            best, bmu = sp, MZ * np.exp(t)
    return best, bmu, float(np.ptp(ai))


rk_best, rk_mu, rk_2e16 = v246_rk4(np.log(2e16 / MZ))
dop_2e16 = float(np.ptp(run(AINV, 0, np.log(2e16 / MZ), bSM, BSM)))
say(f"   v246: Einschleifen-Spreizung bei 2e16 GeV = {spread_2e16:.3f}; RK4 zwei Schleifen: Minimum {rk_best:.3f} bei "
    f"{e(rk_mu)} GeV, bei 2e16 {rk_2e16:.3f}; DOP853 ohne y_t bei 2e16 {dop_2e16:.3f}")
check("v246 reproduziert: Spreizung 8,8 bei 2e16 (1 Schleife), 2 Schleifen Minimum 3,2 bei 1,7e14 und 8,7 bei 2e16",
      abs(spread_2e16 - 8.8) < 0.05 and abs(rk_best - 3.2) < 0.05 and abs(rk_mu / 1.7e14 - 1) < 0.03 and abs(rk_2e16 - 8.7) < 0.05)
check("Integrator: DOP853 und RK4 (v246) stimmen bei 2e16 auf 0,02 überein", abs(dop_2e16 - rk_2e16) < 0.02)


def ps1(bps, ain=AINV, bsm=bSM):
    """Einschleifen-Zweistufenlösung wie v249/v266: (M_PS, M_GUT, alpha_U^-1)."""
    a = ain
    Bs = bsm / (2 * np.pi)
    C = np.array(bps) / (2 * np.pi)
    lin = [(a[2], -Bs[2], -C[0]), (a[1], -Bs[1], -C[1]),
           ((5 / 3) * a[0] - (2 / 3) * a[2], -((5 / 3) * Bs[0] - (2 / 3) * Bs[2]), -C[2])]
    c4, c2l, c2r = lin
    A = np.array([[c4[1] - c2l[1], c4[2] - c2l[2]], [c4[1] - c2r[1], c4[2] - c2r[2]]])
    L1, L2 = np.linalg.solve(A, [-(c4[0] - c2l[0]), -(c4[0] - c2r[0])])
    return MZ * np.exp(L1), MZ * np.exp(L1 + L2), c4[0] + c4[1] * L1 + c4[2] * L2


v266 = {"minimal": ps1(b_min), "+(15,1,1) komplex (v266)": ps1(b_15c), "+(15,1,1) reell (E8)": ps1(b_15r)}
for k, (mp, mg, ag) in v266.items():
    say(f"   1 Schleife {k:26s}: M_PS = {e(mp)}, M_GUT = {e(mg)}, alpha_U^-1 = {ag:.2f}, tau_p = {e(tau_p(mg, ag))} a")
t15r = tau_p(v266["+(15,1,1) reell (E8)"][1], v266["+(15,1,1) reell (E8)"][2])
befund(f"v266 zählt (15,1,1) komplex, das sind zwei reelle Adjungierte; die E8-Hülle hat nur eine reelle 45. Damit sinkt M_GUT schon "
       f"in einer Schleife von 5,9e15 auf {e(v266['+(15,1,1) reell (E8)'][1], 1)} GeV und tau_p auf {e(t15r, 1)} a (< Super-K).")
check("v266 reproduziert: minimal M_PS = 4,2e13, M_GUT = 2,4e15; +(15,1,1) M_GUT = 5,9e15 (auf 3 %)",
      abs(v266["minimal"][0] / 4.2e13 - 1) < 0.03 and abs(v266["minimal"][1] / 2.4e15 - 1) < 0.03
      and abs(v266["+(15,1,1) komplex (v266)"][1] / 5.9e15 - 1) < 0.03)
# 2-Schleifen-Hybrid aus experiments/gauge-unification (2L SM ohne y_t, 1L PS, dortige PDG-Eingaben)
A_EXP = np.array([(3 / 5) * (1 - 0.23122) * 127.951, 0.23122 * 127.951, 1 / 0.1179])


def hybrid(bps):
    C = np.array(bps) / (2 * np.pi)

    def res(lnM):
        a = run(A_EXP, 0, lnM - np.log(MZ), bSM, BSM)
        a4, a2l, a2r = a[2], a[1], (5 / 3) * a[0] - (2 / 3) * a[2]
        t12, t13 = (a4 - a2l) / (C[0] - C[1]), (a4 - a2r) / (C[0] - C[2])
        return t12 - t13, t12
    lnM = brentq(lambda x: res(x)[0], np.log(1e12), np.log(1e15), xtol=1e-12)
    return np.exp(lnM), np.exp(lnM + res(lnM)[1])


hy_min, hy_15 = hybrid(b_min), hybrid(b_15c)
say(f"   Hybrid (2L SM, 1L PS) wie experiments/gauge-unification: minimal M_PS = {e(hy_min[0])}, M_GUT = {e(hy_min[1])};"
    f" +(15,1,1) M_GUT = {e(hy_15[1])} (dort 3,74e13/1,29e15 und 2,80e15)")
check("2-Schleifen-Hybrid aus results/pati_salam.json reproduziert (M_GUT 1,29e15 und 2,80e15 auf 1 %)",
      abs(hy_min[1] / 1.2906e15 - 1) < 0.01 and abs(hy_15[1] / 2.7989e15 - 1) < 0.01)
ZAHLEN["kalibrierung"] = {"v246_spreizung_1L_2e16": spread_2e16, "v246_2L_min": [rk_best, rk_mu], "v246_2L_2e16": rk_2e16,
                          "v266_1L": {k: list(v) for k, v in v266.items()}, "hybrid_2L": {"min": hy_min, "+15c": hy_15}}

# ============================================================================ Teil 2: Frage 1
say("\nTeil 2 (Frage 1): reines SM von M_Z bis 1/a, zwei Schleifen mit y_t")
q1 = {}
for k, mu in MU_A.items():
    a1l = sm_at(mu, 1)
    a2l = sm_at(mu, 2)
    a2n = sm_at(mu, 2, bSM0, BSM0)
    q1[k] = {"mu": mu, "a_1L": a1l.tolist(), "a_2L": a2l.tolist(), "D_1L": float(np.ptp(a1l)),
             "D_2L": float(np.ptp(a2l)), "D_2L_ohne_Higgs": float(np.ptp(a2n)), "sin2_2L": sin2(a2l)}
    say(f"   a/l_P = {ALP[k]:5.2f}, 1/a = {e(mu)}: alpha^-1 (2L) = {np.round(a2l, 2).tolist()}, Delta = {np.ptp(a2l):.2f}"
        f" (1L {np.ptp(a1l):.2f}; ohne Higgs-Import {np.ptp(a2n):.2f}); sin^2 theta_W(1/a) = {sin2(a2l):.4f} statt 0,375")
a_noy = run(AINV, 0, np.log(MU_A["H6-Schema 10,79"] / MZ), bSM, BSM)
dy = float(np.abs(a_noy - np.array(q1["H6-Schema 10,79"]["a_2L"])).max())
check(f"y_t in führender Ordnung genügt: Zwei-Schleifen-Yukawabeitrag an 1/a = {dy:.3f} < 0,1 in alpha^-1", dy < 0.1)
dense = solve_ivp(rhs, (0, np.log(1e19 / MZ)), Y0, args=(bSM, BSM, True), method="DOP853", rtol=1e-11, atol=1e-12,
                  dense_output=True)
cross = {}
for (i, j), lab in (((0, 1), "a1=a2"), ((0, 2), "a1=a3"), ((1, 2), "a2=a3")):
    t = brentq(lambda t: dense.sol(t)[i] - dense.sol(t)[j], 1, np.log(1e19 / MZ))
    cross[lab] = MZ * np.exp(t)
say("   Paarweise Kreuzungen (2L): " + ", ".join(f"{k} bei {e(v)} GeV" for k, v in cross.items()))
dmin = min(v["D_2L"] for v in q1.values())
say(f"   Kleinste Spreizung im Fenster: {dmin:.2f}. alpha_s-Fehler verschiebt um höchstens {DA3}; Zwei- gegen Einschleifen")
say("   ändert um < 0,8. Gitterschwellen von H6 (Eltern (10,1) = 5+5b, Verdoppler je ganze Generation) sind vollständige")
say("   SU(5)-Multipletts und verschieben alle drei gleich, also nicht die Spreizung.")
el = fl(betas(SMG, [((Fr(1, 2), 2, 1), "W", 1), ((Fr(-1, 2), 2, 1), "W", 1), ((Fr(-1, 3), 1, 3), "W", 1), ((Fr(1, 3), 1, 3), "W", 1)]))[0]
el = el - np.array([0, -22 / 3, -11.0])
check("H6-Eltern (10,1) unter SM: Delta b = (2/3, 2/3, 2/3), keine differentielle Schwelle an 1/a", np.allclose(el, 2 / 3))
need = 2 * np.pi * dmin
say(f"   Zum Schließen bräuchte es eine Schwellendifferenz Delta b x ln(Lambda/M) = 2 pi Delta = {need:.0f}.")
check("Frage 1: SM vereinigt nicht an 1/a — Spreizung (2L) > 1 im ganzen Fenster (negatives Ergebnis, Kill)",
      all(v["D_2L"] > 1.0 for v in q1.values()) and all(v["D_1L"] > 1.0 for v in q1.values()))
befund(f"Reines SM (mit importiertem Higgs) verfehlt die Gleichheit an 1/a im ganzen Fenster: Delta alpha^-1 = "
       f"{q1['untere Grenze 1,65']['D_2L']:.1f} / {q1['H6-Schema 10,79']['D_2L']:.1f} / {q1['obere Grenze 14,27']['D_2L']:.1f}"
       f" (a/l_P = 1,65 / 10,8 / 14,3, zwei Schleifen); sin^2 theta_W(1/a) = {q1['H6-Schema 10,79']['sin2_2L']:.3f} statt 0,375.")
ZAHLEN["frage1_sm"] = q1
ZAHLEN["frage1_kreuzungen_GeV"] = cross

# ============================================================================ vektorisierte 2-Schleifen-Zweistufe


def rk4_vec(f, y, nst):
    h = 1.0 / nst
    for _ in range(nst):
        k1 = f(y); k2 = f(y + 0.5 * h * k1); k3 = f(y + 0.5 * h * k2); k4 = f(y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return y


def sm_stage_vec(L1, nst=160):
    def f(y):
        al = 1 / y[:, :3]
        g2 = 4 * np.pi * al
        yt = y[:, 3]
        da = -bSM / (2 * np.pi) - al @ BSM.T / (8 * np.pi ** 2) + CT[None, :] * yt[:, None] ** 2 / (32 * np.pi ** 3)
        dy = yt / (16 * np.pi ** 2) * (4.5 * yt ** 2 - 0.85 * g2[:, 0] - 2.25 * g2[:, 1] - 8 * g2[:, 2])
        return np.c_[da, dy] * L1[:, None]
    return rk4_vec(f, np.tile(Y0, (len(L1), 1)), nst)


def gauge_stage_vec(A, L, b, B, nst=160):
    def f(a):
        return (-b / (2 * np.pi) - np.einsum("nij,nj->ni", B, 1 / a) / (8 * np.pi ** 2)) * L[:, None]
    return rk4_vec(f, A, nst)


def two_step_vec(bP, BP, L1, L2, iters=10):
    """Newton für (ln M_PS/M_Z, ln M_U/M_PS): a4 = a2L = a2R an M_U, SM 2L + y_t unten, PS 2L oben."""
    def F(L1, L2):
        A = sm_stage_vec(L1)[:, :3] @ MM.T
        A2 = gauge_stage_vec(A, L2, bP, BP)
        return np.c_[A2[:, 0] - A2[:, 1], A2[:, 0] - A2[:, 2]], A2
    eps = 1e-5
    for _ in range(iters):
        F0, _ = F(L1, L2)
        F1, _ = F(L1 + eps, L2)
        F2, _ = F(L1, L2 + eps)
        J = np.stack([(F1 - F0) / eps, (F2 - F0) / eps], axis=2)
        d = np.linalg.solve(J, -F0[:, :, None])[:, :, 0]
        d = np.clip(d, -3, 3)
        L1, L2 = L1 + d[:, 0], L2 + d[:, 1]
    F0, A2 = F(L1, L2)
    return L1, L2, A2[:, 0], np.abs(F0).max(axis=1)


def ps_bB(counts):
    """b, B der PS-Stufe für Anzahlen in der Reihenfolge PAR (Fermionen + Eichteil + Skalare)."""
    counts = np.atleast_2d(counts)
    b = bPSF[None, :] + counts @ np.array([DB_PS[k][0] for k in PAR])
    B = BPSF[None, :, :] + np.einsum("nk,kij->nij", counts, np.array([DB_PS[k][1] for k in PAR]))
    return b, B


# ============================================================================ Teil 3: Frage 2
say("\nTeil 3 (Frage 2): welcher E8-erlaubte Zwischeninhalt gibt Gleichheit genau an 1/a?")
say("   E8-Hülle (v247): SO(10)-Inhalt {1, 10, 16, 45}, kein 126. Höchstzahlen aus einer 248: sechs reelle 10 (aus (10,6)),")
say("   vier komplexe 16 (aus (16,4)+(16b,4b)), eine reelle 45. Fermionen fest: 3 x 16 leicht, Eltern (10,1)+(1,6) bei m/a.")
# --- 3A: diskreter Scan, vollständige PS-Multipletts
say("\n   3A  SM -> PS bei M_PS, vollständige PS-Multipletts (Überlebenshypothese), alle 13440 Inhalte:")
ranges = [range(lb, n + 1) for (_, _, n, lb) in PS_COMP.values()]
CNT = np.array(list(itertools.product(*ranges)), float)
bP_all, BP_all = ps_bB(CNT)
one = []
for bp in bP_all:
    try:
        one.append(ps1(bp))
    except np.linalg.LinAlgError:
        one.append((np.nan, np.nan, np.nan))
one = np.array(one)
ok1 = np.isfinite(one[:, 0]) & (one[:, 0] > 1e3) & (one[:, 1] > one[:, 0])
L1s = np.where(ok1, np.log(np.where(ok1, one[:, 0], 4e13) / MZ), np.log(4e13 / MZ))
L2s = np.where(ok1, np.log(np.where(ok1, one[:, 1] / one[:, 0], 50.0)), np.log(50.0))
L1v, L2v, aUv, resv = two_step_vec(bP_all, BP_all, L1s, L2s)
MPSv, MUv = MZ * np.exp(L1v), MZ * np.exp(L1v + L2v)
ok2 = ok1 & (resv < 1e-8) & (L2v > 0)
say(f"   1 Schleife: {ok1.sum()} gültige Lösungen, größtes M_U = {e(one[ok1, 1].max())} GeV; im Fenster: "
    f"{int(((one[:, 1] >= WIN[0]) & (one[:, 1] <= WIN[1]) & ok1).sum())}")
say(f"   2 Schleifen (vektorisiert, Newton): {ok2.sum()} konvergiert (Residuum < 1e-8), größtes M_U = {e(MUv[ok2].max())} GeV;"
    f" im Fenster: {int(((MUv >= WIN[0]) & (MUv <= WIN[1]) & ok2).sum())}")
# Gegenprobe gegen solve_ivp/fsolve für drei Inhalte


def two_step_ivp(sc_counts, guess):
    b, B = ps_bB(sc_counts)
    b, B = b[0], B[0]

    def F(v):
        y = run(Y0, 0, v[0], bSM, BSM, True)
        A2 = run(MM @ y[:3], 0, v[1], b, B)
        return [A2[0] - A2[1], A2[0] - A2[2]]
    v = fsolve(F, guess, xtol=1e-13)
    y = run(Y0, 0, v[0], bSM, BSM, True)
    return MZ * np.exp(v[0]), MZ * np.exp(v[0] + v[1]), run(MM @ y[:3], 0, v[1], b, B)[0]


def idx_of(cnt):
    return int(np.where((CNT == np.array(cnt, float)).all(1))[0][0])


C_MIN = [2, 0, 1, 0, 0, 0, 0, 0]          # komplexes Bidublett = zwei reelle, komplexe (4,1,2)
C_E8 = [2, 0, 1, 0, 1, 0, 0, 0]           # + eine reelle (15,1,1)
i_best = int(np.argmax(np.where(ok2, MUv, 0)))
cmp_ok = True
for cnt in (C_MIN, C_E8, CNT[i_best].tolist()):
    i = idx_of(cnt)
    ref = two_step_ivp(cnt, [L1v[i], L2v[i]])
    cmp_ok &= abs(ref[1] / MUv[i] - 1) < 2e-3 and abs(ref[2] - aUv[i]) < 0.01
check("Vektorisierter 2L-Löser stimmt mit solve_ivp/fsolve überein (M_U auf 0,2 %, alpha_U^-1 auf 0,01; drei Inhalte)", cmp_ok)
best_cnt = dict(zip(PAR, CNT[i_best].astype(int).tolist()))
say(f"   Inhalt mit größtem M_U (2L): {best_cnt}: M_PS = {e(MPSv[i_best])} (kappa = {MPSv[i_best] / MS:.2f}),"
    f" M_U = {e(MUv[i_best])}, alpha_U^-1 = {aUv[i_best]:.2f}")


def lp_ps_free(MPS, MU):
    """Volle PS-Multipletts, Basis (1 bd_r, 1 s412_c) bei M_PS, alle weiteren Komponenten auf freien Schwellen in [M_PS, M_U]."""
    L1, L2 = np.log(MPS / MZ), np.log(MU / MPS)
    a = AINV - bSM * L1 / (2 * np.pi)
    A0 = MM @ a - (bPSF + DB_PS["bd_r"][0] + DB_PS["s412_c"][0]) * L2 / (2 * np.pi)
    D = np.array([DB_PS[k][0] for k in PAR]) / (2 * np.pi)
    nmax = np.array([PS_COMP[k][2] - PS_COMP[k][3] for k in PAR], float)
    r = linprog(np.zeros(len(PAR)), A_eq=np.array([D[:, 0] - D[:, 1], D[:, 0] - D[:, 2]]),
                b_eq=[A0[0] - A0[1], A0[0] - A0[2]], bounds=[(0, n * L2) for n in nmax])
    return r.status == 0


def lp_ps_max(MPS):
    """Für festes M_PS ist alles linear in (L2, w): maximiere L2 = ln(M_U/M_PS)."""
    L1 = np.log(MPS / MZ)
    a = MM @ (AINV - bSM * L1 / (2 * np.pi))
    bb = bPSF + DB_PS["bd_r"][0] + DB_PS["s412_c"][0]
    D = np.array([DB_PS[k][0] for k in PAR])
    nmax = np.array([PS_COMP[k][2] - PS_COMP[k][3] for k in PAR], float)
    Aeq = np.array([np.r_[-(bb[0] - bb[1]), -(D[:, 0] - D[:, 1])], np.r_[-(bb[0] - bb[2]), -(D[:, 0] - D[:, 2])]]) / (2 * np.pi)
    Aub = np.c_[-nmax, np.eye(len(PAR))]
    r = linprog(np.r_[-1.0, np.zeros(len(PAR))], A_ub=Aub, b_ub=np.zeros(len(PAR)), A_eq=Aeq,
                b_eq=[-(a[0] - a[1]), -(a[0] - a[2])], bounds=[(0, np.log(1e20 / MPS))] + [(0, None)] * len(PAR))
    return MPS * np.exp(r.x[0]) if r.status == 0 else 0.0


grid_ps = np.geomspace(1e4, 1e17, 400)
muA = np.array([lp_ps_max(mp) for mp in grid_ps])
jA = int(np.argmax(muA))
fine = np.geomspace(grid_ps[max(jA - 1, 0)], grid_ps[min(jA + 1, len(grid_ps) - 1)], 200)
muAf = np.array([lp_ps_max(mp) for mp in fine])
mu_maxA, mps_maxA = float(muAf.max()), float(fine[int(np.argmax(muAf))])
feasA_ps = grid_ps[muA > 0]
say(f"   Freie Schwellen oberhalb M_PS (exaktes LP je M_PS, 1 Schleife): größtes M_U = {e(mu_maxA)} GeV bei M_PS = {e(mps_maxA)}"
    f" (kappa {mps_maxA / MS:.2f}); lösbar für M_PS = {e(feasA_ps.min())} ... {e(feasA_ps.max())}")
check("3A: Mit vollständigen E8-PS-Multipletts liegt kein M_U im Fenster (diskret 1L und 2L, freie Schwellen 1L) — negatives Ergebnis",
      ((one[:, 1] >= WIN[0]) & ok1).sum() == 0 and ((MUv >= WIN[0]) & ok2).sum() == 0 and mu_maxA < WIN[0])
befund(f"Mit vollständigen PS-Multipletts der E8-Hülle ist M_U höchstens {e(MUv[ok2].max(), 1)} GeV (2 Schleifen, 13440 Inhalte)"
       f" bzw. {e(mu_maxA, 1)} GeV (1 Schleife, freie Schwellen) — Faktor {WIN[0] / mu_maxA:.0f} unter dem Fenster.")
ZAHLEN["frage2_3A"] = {"anzahl_inhalte": int(len(CNT)), "max_MU_1L": float(one[ok1, 1].max()), "max_MU_2L": float(MUv[ok2].max()),
                       "max_MU_freie_schwellen_1L": float(mu_maxA), "bester_inhalt": best_cnt,
                       "bester_MPS_2L": float(MPSv[i_best]), "bester_MU_2L": float(MUv[i_best])}

# --- 3B: SM-Gruppe bis 1/a
say("\n   3B  SM-Eichgruppe bis 1/a (H6 ohne neue Eichbosonen), E8-Fragmente auf freien Schwellen >= M_low (LP, 1 Schleife):")
FK = list(FRAG)
DFR = np.array([DB_FR[k][0] for k in FK])
NFR = np.array([FRAG[k][2] for k in FK], float)


def lp_smfrag(MU, low, exclude=(), obj=None):
    T = np.log(MU / MZ)
    Lmax = np.log(MU / low)
    aMU = AINV - bSM * T / (2 * np.pi)
    Aeq = np.array([DFR[:, 0] - DFR[:, 1], DFR[:, 1] - DFR[:, 2]]) / (2 * np.pi)
    c = np.zeros(len(FK)) if obj is None else obj
    return linprog(c, A_eq=Aeq, b_eq=[aMU[0] - aMU[1], aMU[1] - aMU[2]],
                   bounds=[(0, 0 if k in exclude else n * Lmax) for k, n in zip(FK, NFR)])


def bisect(fun, lo, hi, n=50):
    for _ in range(n):
        m = np.sqrt(lo * hi)
        if fun(m):
            lo = m
        else:
            hi = m
    return lo


def lp_smfrag_max(low):
    """Linear in (T = ln M_U/M_Z, x): maximiere T bei Fragmenten zwischen low und M_U."""
    tl = np.log(low / MZ)
    Aeq = np.array([np.r_[-(bSM[0] - bSM[1]), -(DFR[:, 0] - DFR[:, 1])], np.r_[-(bSM[1] - bSM[2]), -(DFR[:, 1] - DFR[:, 2])]]) / (2 * np.pi)
    Aub = np.c_[-NFR, np.eye(len(FK))]
    r = linprog(np.r_[-1.0, np.zeros(len(FK))], A_ub=Aub, b_ub=-NFR * tl, A_eq=Aeq,
                b_eq=[-(AINV[0] - AINV[1]), -(AINV[1] - AINV[2])], bounds=[(tl, np.log(1e20 / MZ))] + [(0, None)] * len(FK))
    return MZ * np.exp(r.x[0]) if r.status == 0 else 0.0


muB = lp_smfrag_max(1e3)
lowB = bisect(lambda lo: lp_smfrag(WIN[0], lo).status == 0, 30, 1e5)
rB = lp_smfrag(WIN[0], 0.99 * lowB, obj=-np.ones(len(FK)))
LB = np.log(WIN[0] / (0.99 * lowB))
compB = {k: round(float(v / LB), 2) for k, v in zip(FK, rB.x) if v / LB > 5e-3}
safeB_any = lp_smfrag(WIN[0], 1e-3, exclude=GEFAEHRLICH).status == 0
say(f"   Größtes M_U mit Fragmenten >= 1 TeV: {e(muB)} GeV (Fensterrand {e(WIN[0])}).")
say(f"   Für M_U = Fensterrand müssen Fragmente bis {lowB:.0f} GeV hinunter; Lösung (Anzahl voll leicht): {compB}")
say(f"   Ohne Farbtripletts D, XY (protonsicher) erreichbar, selbst mit Fragmenten ab 1 MeV: {safeB_any}.")
check("3B: SM-Gruppe mit E8-Fragmenten erreicht das Fenster nicht mit Schwellen >= 1 TeV (exaktes LP-Maximum) — negatives Ergebnis",
      muB < WIN[0] and lowB < 1e3 and not safeB_any)
befund(f"Ohne neue Eichbosonen (H6 wie gebaut) reicht die Gleichheit mit E8-Fragmenten oberhalb 1 TeV nur bis {e(muB, 1)} GeV;"
       f" das Fenster verlangt ein Farboktett, Tripletts und Dubletts bei ~{lowB:.0f} GeV (LHC-ausgeschlossen).")
ZAHLEN["frage2_3B"] = {"max_MU_ab_1TeV": muB, "noetige_untergrenze_GeV": lowB, "zusammensetzung": compB,
                       "ohne_tripletts_erreichbar": bool(safeB_any)}

# --- 3C: PS mit gespaltenen Multipletts
say("\n   3C  SM -> (Fragmente) -> PS: leichte SM-Fragmente unter M_PS, ihre PS-Eltern ab M_PS (LP-Relaxation, 1 Schleife):")
PAIRS = [(p, k) for p in PAR for k in FRAG_PARENT[p]]
NZ, NX = len(PAR), len(PAIRS)


def lp_split(MPS, MU, low=1e3, exclude=(), obj=False):
    """Variablen: z_p (Eltern bei M_PS, kontinuierlich relaxiert), w_p (weitere Eltern oberhalb), x_kp (Fragmente
    unter M_PS, x = Anzahl x ln(M_PS/M_k), nur aus vorhandenen Eltern; ein H aus dem Basis-Bidublett ist der SM-Higgs)."""
    L1, L2, Lf = np.log(MPS / MZ), np.log(MU / MPS), np.log(MPS / low)
    if Lf <= 0 or L2 <= 0:
        return None
    a = AINV - bSM * L1 / (2 * np.pi)
    A0 = MM @ a - bPSF * L2 / (2 * np.pi)
    cols = [DB_PS[p][0] * L2 for p in PAR] + [DB_PS[p][0] for p in PAR] + [MM @ DB_FR[k][0] for _, k in PAIRS]
    C = np.array(cols).T / (2 * np.pi)
    n = 2 * NZ + NX
    Aub, bub = [], []
    for i, p in enumerate(PAR):
        row = np.zeros(n); row[NZ + i] = 1; row[i] = L2
        Aub.append(row); bub.append(PS_COMP[p][2] * L2)
    for j, (p, k) in enumerate(PAIRS):
        row = np.zeros(n); row[2 * NZ + j] = 1; row[PAR.index(p)] = -Lf
        Aub.append(row); bub.append(-Lf if (p, k) == ("bd_r", "H_c") else 0.0)
    bounds = [(PS_COMP[p][3], PS_COMP[p][2]) for p in PAR] + [(0, None)] * NZ
    bounds += [(0, 0) if k in exclude else (0, None) for _, k in PAIRS]
    c = np.zeros(n)
    if obj:
        c[2 * NZ:] = 1.0
    r = linprog(c, A_ub=np.array(Aub), b_ub=bub, A_eq=np.array([C[0] - C[1], C[0] - C[2]]), b_eq=[A0[0] - A0[1], A0[0] - A0[2]],
                bounds=bounds)
    return r


MU_GRID = np.geomspace(WIN[0], WIN[1], 20)
KBAND = np.geomspace(0.3 * MS, 3 * MS, 12)
KPIN = np.geomspace(1.0 * MS, 1.7 * MS, 6)


def any_feasible(mps_list, low=1e3, exclude=()):
    for mp in mps_list:
        for mu in MU_GRID:
            r = lp_split(mp, mu, low, exclude)
            if r is not None and r.status == 0:
                return True
    return False


feasC_pin = any_feasible(KPIN)
feasC_safe_pin = any_feasible(KPIN, exclude=GEFAEHRLICH)
safe_lo_ok = any_feasible([0.05 * MS], exclude=GEFAEHRLICH)
kmax_safe = bisect(lambda k: any_feasible([k * MS], exclude=GEFAEHRLICH), 0.05, 3.0, 30) if safe_lo_ok else float("nan")
low_pin = bisect(lambda lo: any_feasible(KPIN, low=lo), 1e2, 1e13, 40)
low_safe_pin = bisect(lambda lo: any_feasible(KPIN, low=lo, exclude=GEFAEHRLICH), 1e1, 1e8, 40) if feasC_safe_pin else float("nan")
r_ex = lp_split(MS, MU_A["H6-Schema 10,79"], obj=True)
x_ex = r_ex.x[2 * NZ:] / np.log(MS / 1e3)
needC = {f"{k} (aus {p})": round(float(v), 2) for (p, k), v in zip(PAIRS, x_ex) if v > 1e-6}
say(f"   kappa = M_PS/M_s in [1,0; 1,7] (PS.KAPPA.01): erreichbar {feasC_pin}; nur mit protonsicheren Fragmenten (ohne D, XY): "
    f"{feasC_safe_pin}. Protonsicher erreichbar nur bis kappa = {kmax_safe:.2f}, und nur mit leichten Fragmenten bis "
    f"{low_safe_pin:.0f} GeV hinunter (am Fensterrand: Oktett, ~5 (3b,1,-2/3), ~2 (3,2,1/6) an der Untergrenze).")
say(f"   Bei kappa in [1,0; 1,7] müssen leichte Fragmente bis {e(low_pin)} GeV hinunter; minimale leichte Last bei M_PS = M_s,")
say(f"   1/a(H6): {needC} (in Einheiten 'voll leicht ab 1 TeV').")
# ganzzahlige Aufzählung: Fragmente auf einer gemeinsamen Schwelle M_f, minimale Eltern + Zusatz t6/s421 bei M_PS
FT = list(FRAG)
RNG = {"O8_r": 2, "T3_r": 2, "XY_c": 2, "Q_c": 4, "U_c": 3, "E_c": 3, "H_c": 4, "D_c": 5}
CF = np.array(list(itertools.product(*[range(RNG[k]) for k in FT])), int)
O8, T3, XY, Q, U, E, H, Dd = CF.T
a15, a31, a622 = O8.copy(), T3.copy(), XY.copy()
s421 = np.maximum(0, Q - a622)
bd = 1 + np.maximum(0, H - s421)
s412 = np.maximum(1, U - a15)
a13 = (E > s412).astype(int)
s412 = np.maximum(s412, E - a13)
t6 = np.maximum(0, Dd - s412)
Z0 = np.stack([bd, t6, s412, s421, a15, a31, a13, a622], 1)
DBF = np.array([DB_FR[k][0] for k in FT])
DPP = np.array([DB_PS[p][0] for p in PAR])
MU_H6 = MU_A["H6-Schema 10,79"]
TU = np.log(MU_H6 / MZ)
sols = []
for te, se in itertools.product(range(7), range(5)):
    Z = Z0.copy(); Z[:, 1] += te; Z[:, 3] += se
    okZ = (Z[:, 0] <= 6) & (Z[:, 1] <= 6) & (Z[:, 2] <= 4) & (Z[:, 3] <= 4)
    bP = bPSF + Z @ DPP
    dF = CF @ DBF
    c0 = (MM @ AINV)[None, :] - bP * TU / (2 * np.pi)
    cL = -(MM @ bSM)[None, :] / (2 * np.pi) + bP / (2 * np.pi)
    cl = -(dF @ MM.T) / (2 * np.pi)
    A11, A12, A21, A22 = cL[:, 0] - cL[:, 1], cl[:, 0] - cl[:, 1], cL[:, 0] - cL[:, 2], cl[:, 0] - cl[:, 2]
    r1, r2 = -(c0[:, 0] - c0[:, 1]), -(c0[:, 0] - c0[:, 2])
    det = A11 * A22 - A12 * A21
    with np.errstate(all="ignore"):
        L1 = (r1 * A22 - A12 * r2) / det
        lf = (A11 * r2 - A21 * r1) / det
        MPS_, MF_ = MZ * np.exp(L1), MZ * np.exp(L1 - lf)
    ok = okZ & np.isfinite(L1) & (MPS_ > 1e3) & (MPS_ < MU_H6) & (lf >= 0) & (MF_ >= 1e3) & (CF.sum(1) > 0)
    for i in np.where(ok)[0]:
        sols.append((CF[i].tolist(), Z[i].tolist(), float(MPS_[i]), float(MF_[i])))
kap = np.array([s[2] / MS for s in sols])
safe = np.array([s[0][FT.index("D_c")] == 0 and s[0][FT.index("XY_c")] == 0 for s in sols])
nfr = np.array([sum(s[0]) for s in sols])
sel_pin = np.where((kap >= 1.0) & (kap <= 1.7))[0]
say(f"   Ganzzahlig (eine gemeinsame Fragmentschwelle >= 1 TeV, 1/a = {e(MU_H6)}): {len(sols)} Lösungen; kappa in [1; 1,7]: "
    f"{len(sel_pin)}; davon protonsicher {int(safe[sel_pin].sum())}; protonsicher bei kappa in [0,3; 3]: "
    f"{int((safe & (kap >= 0.3) & (kap <= 3)).sum())}.")
i_ex = sel_pin[int(np.argmax([sols[j][3] for j in sel_pin]))]
ex = sols[i_ex]
ex_frag = {k: n for k, n in zip(FT, ex[0]) if n}
ex_par = {p: n for p, n in zip(PAR, ex[1]) if n}
say(f"   Beispiel (größtes M_f bei kappa in [1; 1,7]): Fragmente {ex_frag} bei M_f = {e(ex[3])} GeV, PS-Eltern {ex_par} ab M_PS = {e(ex[2])}"
    f" (kappa = {ex[2] / MS:.2f}), M_U = 1/a = {e(MU_H6)}.")
# 2-Schleifen-Bestätigung des Beispiels: Unbekannte (ln M_PS, ln M_f), M_U fest
fr_fields = [(FRAG[k][0], FRAG[k][1], n) for k, n in ex_frag.items()]
bSF, BSF = fl(betas(SMG, SM_FERM + HIGGS + fr_fields))
bPX, BPX = ps_bB(np.array(ex[1], float))
bPX, BPX = bPX[0], BPX[0]


def F_ex(v, MU=MU_H6):
    lnps, lnf = v
    y = run(Y0, 0, lnf - np.log(MZ), bSM, BSM, True)
    y = run(y, lnf, lnps, bSF, BSF, True)
    A2 = run(MM @ y[:3], lnps, np.log(MU), bPX, BPX)
    return [A2[0] - A2[1], A2[0] - A2[2]], A2


v_ex = fsolve(lambda v: F_ex(v)[0], [np.log(ex[2]), np.log(ex[3])], xtol=1e-13)
res_ex, A_ex = F_ex(v_ex)
MPS_ex2, MF_ex2 = np.exp(v_ex)
say(f"   2 Schleifen: M_PS = {e(MPS_ex2)} (kappa = {MPS_ex2 / MS:.2f}), M_f = {e(MF_ex2)} GeV, alpha_U^-1(1/a) = {A_ex[0]:.2f},"
    f" Residuum {max(abs(x) for x in res_ex):.1e}")
tau_ex = tau_p(MU_H6, A_ex[0])
# Skalare Tripletts: Rate ~ (y y'/M_T^2)^2 statt (g^2/M_X^2)^2; y y' ~ y_u y_d (erste Generation, 10_H-artig)
yy = 1.2e-5 * 2.7e-5
gU2 = 4 * np.pi / A_ex[0]
MT_min = 1e16 * (TAU_SK / (1e36 * (A_ex[0] / 40) ** 2 * (gU2 / yy) ** 2)) ** 0.25
say(f"   Protonzerfall über Eichbosonen (d = 6, M_U = 1/a): tau_p = {e(tau_ex)} a — weit jenseits Hyper-K ({e(TAU_HK)} a).")
say(f"   Skalare Farbtripletts (3,1,-1/3) mit 10_H-artigen Kopplungen (y y' ~ y_u y_d = {yy:.1e}) verlangen M_T > {e(MT_min, 1)} GeV;")
say(f"   im Beispiel liegen {ex_frag.get('D_c', 0)} solche bei {e(MF_ex2, 1)} GeV: ohne zusätzliche Symmetrie (in H6 nicht hergeleitet) ausgeschlossen.")
check("3C: PS mit gespaltenen Multipletts erreicht 1/a (LP und ganzzahlig, 2L bestätigt mit M_f >= 1 TeV)",
      feasC_pin and len(sel_pin) > 0 and max(abs(x) for x in res_ex) < 1e-8 and MF_ex2 >= 1e3 and MPS_ex2 < MU_H6)
check("3C: protonsicher (ohne D, XY) bei kappa >= 1 nur am Bandrand (kappa_max < 1,2) mit Farbskalaren unter 3 TeV; ganzzahlig mit einer Schwelle keine sichere Lösung",
      kmax_safe < 1.2 and low_safe_pin < 3e3 and safe[sel_pin].sum() == 0)
befund(f"Gleichheit genau an 1/a geht mit E8-Inhalt nur über gespaltene Multipletts: bei M_PS auf der Skalaronskala "
       f"(kappa 1,0-1,7) mit einer Fragmentschwelle stets Farboktett + >= 3 Farbtripletts (3,1,-1/3) bei 1e3-1e5 GeV; "
       f"protonsicher nur bis kappa = {kmax_safe:.2f} mit ~8 Farbskalaren bei ~{low_safe_pin / 1e3:.0f} TeV.")
befund(f"In jeder Lesart mit Gleichheit genau an 1/a ist der Eich-Protonzerfall unbeobachtbar (tau_p ~ {e(tau_ex, 0)} a); Hyper-K sähe nichts.")
ZAHLEN["frage2_3C"] = {"lp_kappa_1_1.7_erreichbar": feasC_pin, "lp_protonsicher_kappa_1_1.7": feasC_safe_pin,
                       "kappa_max_protonsicher": kmax_safe, "protonsicher_fragmentskala_GeV": low_safe_pin, "noetige_fragmentskala_GeV": low_pin, "lp_minimale_last": needC,
                       "ganzzahlig_loesungen": len(sols), "ganzzahlig_kappa_1_1.7": int(len(sel_pin)),
                       "beispiel": {"fragmente": ex_frag, "eltern": ex_par, "MPS_1L": ex[2], "Mf_1L": ex[3], "MPS_2L": MPS_ex2,
                                    "Mf_2L": MF_ex2, "alpha_U_inv_2L": A_ex[0], "tau_p_eich_a": tau_ex,
                                    "MT_min_skalar_GeV": MT_min}}

# ============================================================================ Teil 4: Frage 3
say("\nTeil 4 (Frage 3): SO(10) vereinigt unterhalb, eine Kopplung läuft bis 1/a")
br = {"minimal (v266)": C_MIN, "+(15,1,1) komplex (v266, = 2 reelle)": None, "+(15,1,1) reell (E8, eine 45)": C_E8,
      "größtes M_U (3A)": CNT[i_best].tolist()}
q3 = {}
for k, cnt in br.items():
    if cnt is None:
        b, B = fl(betas(PSG, PS_FERM + [((1, 2, 2), "C", 1), ((4, 1, 2), "C", 1), ((15, 1, 1), "C", 1)]))

        def F(v, b=b, B=B):
            y = run(Y0, 0, v[0], bSM, BSM, True)
            A2 = run(MM @ y[:3], 0, v[1], b, B)
            return [A2[0] - A2[1], A2[0] - A2[2]]
        v = fsolve(F, [np.log(3.5e13 / MZ), np.log(80)], xtol=1e-13)
        y = run(Y0, 0, v[0], bSM, BSM, True)
        mp, mu, ag = MZ * np.exp(v[0]), MZ * np.exp(v[0] + v[1]), run(MM @ y[:3], 0, v[1], b, B)[0]
    else:
        i = idx_of(cnt)
        mp, mu, ag = MPSv[i], MUv[i], aUv[i]
    t = tau_p(mu, ag)
    q3[k] = {"M_PS": float(mp), "kappa": float(mp / MS), "M_GUT": float(mu), "alpha_U_inv": float(ag), "tau_p": float(t)}
    say(f"   2L {k:40s}: M_PS = {e(mp)} (kappa {mp / MS:.2f}), M_GUT = {e(mu)}, alpha_U^-1 = {ag:.2f}, tau_p = {e(t)} a"
        f" [{e(t / BAND, 1)} ... {e(t * BAND, 1)}] -> {'Super-K-sicher' if t / BAND > TAU_SK else ('nur im Band' if t * BAND > TAU_SK else 'ausgeschlossen')}")
das = [ps1(b_15r, ain=AINV + np.array([0, 0, d])) for d in (-DA3, DA3)]
rat = [tau_p(x[1], x[2]) / tau_p(v266["+(15,1,1) reell (E8)"][1], v266["+(15,1,1) reell (E8)"][2]) for x in das]
say(f"   alpha_s-Fehler (+-{DA3} in alpha_3^-1): M_GUT x {das[0][1] / v266['+(15,1,1) reell (E8)'][1]:.3f} / "
    f"{das[1][1] / v266['+(15,1,1) reell (E8)'][1]:.3f}, tau_p x {rat[0]:.2f} / {rat[1]:.2f} (1L, E8-Inhalt) — innerhalb des x3-Bands.")
check("v266-Inhalte bei vollen zwei Schleifen von Super-K ausgeschlossen (auch +(15,1,1) komplex am oberen Bandrand)",
      all(q3[k]["tau_p"] * BAND < TAU_SK for k in ("minimal (v266)", "+(15,1,1) komplex (v266, = 2 reelle)", "+(15,1,1) reell (E8, eine 45)")))
tau_all = tau_p(MUv, aUv)
kap_all = MPSv / MS
safe_c = ok2 & (tau_all > TAU_SK)
safe_b = ok2 & (tau_all * BAND > TAU_SK)
pin = (kap_all >= 1.0) & (kap_all <= 1.7)
hk = safe_c & (tau_all < TAU_HK)
say(f"   Scan 2L: Super-K-sicher (zentral) {safe_c.sum()} von {ok2.sum()}, mit Band x3 {safe_b.sum()}; davon kappa in [1; 1,7]:"
    f" {int((safe_c & pin).sum())} / {int((safe_b & pin).sum())}.")
if safe_c.any():
    say(f"   Sichere Inhalte: tau_p = {e(tau_all[safe_c].min())} ... {e(tau_all[safe_c].max())} a, kappa = "
        f"{kap_all[safe_c].min():.2f} ... {kap_all[safe_c].max():.2f}; innerhalb Hyper-K-Reichweite ({e(TAU_HK)}): {hk.sum()}.")
tp = tau_all[safe_b & pin]
say(f"   kappa in [1; 1,7] und Band-sicher: tau_p = {e(tp.min())} ... {e(tp.max())} a (oberer Bandrand {e(tp.max() * BAND)} a)"
    f" — Hyper-K entscheidet diese Gruppe vollständig.")
check("L3 mit kappa in [1; 1,7]: alle Band-sicheren E8-Inhalte liegen auch am oberen Bandrand unter der Hyper-K-Reichweite",
      len(tp) > 0 and tp.max() * BAND < TAU_HK)
# SO(10)-Stufe bis 1/a: Skalare in vollen SO(10)-Multipletts, H6-Eltern (10 Majorana-Weyl) ab m/a = 0,5/a
n10 = np.maximum(CNT[:, 0], CNT[:, 1])
n16 = np.maximum(CNT[:, 2], CNT[:, 3])
n45 = (CNT[:, 4:].sum(1) > 0).astype(float)
b10F, B10F = fl(betas(["SO10"], [((16,), "W", 3)]))
db10 = {k: tuple(x - y for x, y in zip(fl(betas(["SO10"], [((16,), "W", 3), ((d,), kind, 1)])), (b10F, B10F)))
        for k, d, kind in (("10", 10, "R"), ("16", 16, "C"), ("45", 45, "R"), ("par", 10, "W"))}
bS = b10F[0] + n10 * db10["10"][0][0] + n16 * db10["16"][0][0] + n45 * db10["45"][0][0]
BS = B10F[0, 0] + n10 * db10["10"][1][0, 0] + n16 * db10["16"][1][0, 0] + n45 * db10["45"][1][0, 0]
M_PAR = 0.5


def so10_at(mu_a, sel):
    L1 = np.log(M_PAR * mu_a / MUv[sel])
    L2 = np.full(sel.sum(), np.log(1 / M_PAR))
    a = gauge_stage_vec(aUv[sel][:, None], L1, bS[sel][:, None], BS[sel][:, None, None])
    return gauge_stage_vec(a, L2, (bS[sel] + db10["par"][0][0])[:, None], (BS[sel] + db10["par"][1][0, 0])[:, None, None])[:, 0]


aU_a = {k: so10_at(mu, safe_b) for k, mu in MU_A.items()}
for k in ("untere Grenze 1,65", "H6-Schema 10,79", "obere Grenze 14,27"):
    say(f"   alpha_U^-1(1/a) für a/l_P = {ALP[k]:5.2f}: {aU_a[k].min():.1f} ... {aU_a[k].max():.1f} (Band-sichere Inhalte),"
        f" g_U = {np.sqrt(4 * np.pi / aU_a[k].max()):.3f} ... {np.sqrt(4 * np.pi / aU_a[k].min()):.3f}")
say(f"   Eltern zwischen 0,5/a und 1/a: Delta alpha^-1 = (2/3)/(2 pi) ln 2 = {(2 / 3) / (2 * np.pi) * np.log(2):.3f}, für alle gleich.")

# --- y1: induzierte Kopplung nachrechnen
say("\n   Induzierte Kopplung nach y1 (naive Dirac-Fermionen, Fluss in einer Ebene, Masse M als Infrarotregler):")
s1m = np.array([[0, 1], [1, 0]]); s2m = np.array([[0, -1j], [1j, 0]])


def plane_operator(N, m):
    phi = 2 * np.pi * m / N ** 2
    idx = lambda x, y: (x % N) * N + (y % N)
    Tx = np.zeros((N * N, N * N), complex); Ty = np.zeros((N * N, N * N), complex)
    for x in range(N):
        for y in range(N):
            ux = np.exp(-1j * phi * N * y) if x == N - 1 else 1.0
            if x == N - 1:
                ux *= -1
            uy = np.exp(1j * phi * x) * (-1 if y == N - 1 else 1)
            Tx[idx(x, y), idx(x + 1, y)] = ux
            Ty[idx(x, y), idx(x, y + 1)] = uy
    Dx, Dy = (Tx - Tx.conj().T) / 2, (Ty - Ty.conj().T) / 2
    return np.kron(Dx, s1m) + np.kron(Dy, s2m)


def kappa_y1(N, Ms, m=1):
    mu1 = np.linalg.eigvalsh((lambda D: D.conj().T @ D)(plane_operator(N, m)))
    mu0 = np.linalg.eigvalsh((lambda D: D.conj().T @ D)(plane_operator(N, 0)))
    p = 2 * np.pi * (np.arange(N) + 0.5) / N
    s = (np.sin(p)[:, None] ** 2 + np.sin(p)[None, :] ** 2).ravel()
    phi = 2 * np.pi * m / N ** 2
    return np.array([(-np.sum(np.log(mu1[:, None] + s[None, :] + M * M)) + np.sum(np.log(mu0[:, None] + s[None, :] + M * M)))
                     / (0.5 * N ** 4 * phi ** 2) for M in Ms])


KM = [0.3, 0.45, 0.5, 0.6, 1.0]
k24 = kappa_y1(24, [0.3, 0.45, 0.6])
k32 = dict(zip(KM, kappa_y1(32, KM)))
check("y1 reproduziert: 1/g^2(M) = 0,40532 / 0,28405 / 0,20677 bei M = 0,3 / 0,45 / 0,6 (N = 24)",
      np.allclose(k24, [0.40532, 0.28405, 0.20677], atol=2e-5))
slope, icpt = np.polyfit(np.log(1 / np.array([0.3, 0.45, 0.6])), [k32[0.3], k32[0.45], k32[0.6]], 1)
say(f"   N = 32: kappa(M) = " + ", ".join(f"{m}: {k32[m]:.4f}" for m in KM) + f"; Anpassung 1/g^2 = {slope:.4f} ln(1/(aM)) + {icpt:.4f}")
say(f"   (Kontinuum 16 Geschmäcker: 16/(6 pi^2) = {16 / (6 * np.pi ** 2):.4f}). Formel für beliebigen Inhalt:")
say("      1/g_i^2(mu) = sum_f w_f(i) kappa(max(a mu, a m_f)),  w_f = Dirac-Gewicht x Dynkin-Index (Weyl 1/2, je naivem Feld 16 Geschmäcker).")
say("   H6-Inhalt: drei Generationen w = 3 je Faktor (gleiche Dynkin-Summe), Eltern (10,1) w = 1/2 bei m = 0,5; der importierte")
say("   Higgs bricht die Gleichheit minimal (komplexer Skalar, ein Geschmack: w = T/4/16).")
w_light, w_par = 3.0, 0.5
ind = {}
for lab, x in (("mu = 1/a", 1.0), ("mu = 0,5/a", 0.5)):
    k_naiv = w_light * k32[x] + w_par * k32[max(x, M_PAR)]
    ind[lab] = {"naiv_16": 4 * np.pi * k_naiv, "ein_geschmack": 4 * np.pi * k_naiv / 16,
                "g_naiv": 1 / np.sqrt(k_naiv), "g_ein": 4 / np.sqrt(k_naiv)}
    say(f"   {lab}: alpha_ind^-1 = {4 * np.pi * k_naiv:.2f} (naiv, alle 16 Geschmäcker leicht, Obergrenze) bzw. "
        f"{4 * np.pi * k_naiv / 16:.2f} (ein Geschmack, Wilson-artig); g_ind = {1 / np.sqrt(k_naiv):.2f} bzw. {4 / np.sqrt(k_naiv):.2f}")
wH = np.array([3 / 10, 1 / 2, 0]) / 4 / 16
kk = (w_light + w_par) * k32[1.0] + wH * k32[1.0]
say(f"   Mit Higgs: sin^2 theta_W,ind(1/a) = {0.6 / (kk[0] / kk[1] + 0.6):.5f} statt 0,375.")
# Obergrenze inkl. aller SO(10)-Skalaren mit vollem naivem Gewicht (großzügig)
w_s_max = float((n16 * 0.5 + n10 * 0.125 + n45 * 1.0)[safe_b].max()) if safe_b.any() else 0.0
ind_max = 4 * np.pi * (w_light + w_par + w_s_max) * k32[0.5]
aU_min = min(v.min() for v in aU_a.values())
sm_min = min(min(v["a_2L"]) for v in q1.values())
ind_sm = ind["mu = 0,5/a"]["naiv_16"]
say(f"   Großzügigste Obergrenze (naiv, mu = 0,5/a, alle SO(10)-Skalare mit 16-fachem Fermiongewicht T/4): alpha_ind^-1 <= {ind_max:.1f}")
say(f"   RG-Wert an 1/a: SO(10)-Lesart >= {aU_min:.1f} gegen <= {ind_max:.1f} (Verhältnis >= {aU_min / ind_max:.1f}); SM-Kopplungen"
    f" >= {sm_min:.1f} gegen <= {ind_sm:.1f} (Verhältnis >= {sm_min / ind_sm:.1f}); ein Geschmack: Verhältnis >= "
    f"{aU_min / ind['mu = 0,5/a']['ein_geschmack']:.0f}.")
mu_need = MZ * np.exp((AINV[0] - ind_sm) / (bSM[0] / (2 * np.pi)))
say(f"   Bedingung an a/l_P aus dem Wert: alpha_U^-1 wächst oberhalb M_GUT (asymptotisch frei, b = {bS[safe_b].min():.1f} ... {bS[safe_b].max():.1f}),")
say(f"   bleibt also über alpha_U^-1(M_GUT) ~ 45. Kein a/l_P schließt die Lücke; selbst alpha_1^-1 (SM, nicht asymptotisch frei)")
say(f"   fiele erst bei {e(mu_need, 0)} GeV auf {ind_sm:.0f}, also a/l_P ~ {MPL / mu_need:.0e}, weit unter der z15-Schranke 1,65.")
check("Frage 3: y1-induzierter Wert und RG-Wert an 1/a stimmen nirgends im Fenster überein (Verhältnis > 2 auch naiv, SO(10) und SM) — negatives Ergebnis",
      aU_min / ind_max > 2 and sm_min / ind_sm > 2)
befund(f"In der SO(10)-Lesart ist Gleichheit an 1/a automatisch, aber der Wert passt nicht: RG alpha_U^-1(1/a) = "
       f"{aU_min:.0f}-{max(v.max() for v in aU_a.values()):.0f}, y1-induziert höchstens {ind_max:.0f} (naiv) bzw. ~{ind['mu = 0,5/a']['ein_geschmack']:.1f} (ein Geschmack).")
befund(f"Die y1-Formel gibt für den H6-Inhalt g_ind(1/a) = {ind['mu = 1/a']['g_naiv']:.1f} (naiv) bis {ind['mu = 1/a']['g_ein']:.0f} "
       f"(ein Geschmack), nicht g = 0,55: Die induzierte Kopplung ist an der Gitterskala stark.")

# --- z15-Selbstkonsistenz von a/l_P
say("\n   Selbstkonsistenz mit z15 (a/l_P = 2 sqrt(kappa), Kante je Eichboson h_e - ln g, H6-Schema Wilson r = 1, schwer m = 0,5):")
Z15 = json.loads((HERE / "z15_ergebnis.json").read_text())
kw0, kw05 = Z15["kappa_dirac_wilson_r1"]["0.0"], Z15["kappa_dirac_wilson_r1"]["0.5"]
ks05, ks0 = Z15["kappa_skalar_massiv"]["0.5"], Z15["kappa_skalar_masselos"]
gez, C0g = Z15["kappa_photon_ez_destillierbar"], Z15["randfluss"]["innovationsvarianz_e"]


def H_dg(sig, d=0.0):
    K = int(np.ceil(12 * sig)) + 5
    n = np.arange(-K, K + 1)
    w = np.exp(-(n - d) ** 2 / (2 * sig * sig))
    p = w / w.sum()
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def kante(g):
    return float(np.mean([H_dg(np.sqrt(C0g) / g, d) for d in np.linspace(0, 1, 64, endpoint=False)]))


def kappa_tot(gauge, n_scal):
    return 24 * kw0 + 8 * kw05 + 2 * ks05 + n_scal * ks0 + sum(d * (gez + kante(g)) for d, g in gauge)


k_ref = kappa_tot([(1, 0.55), (3, 0.55), (8, 0.55)], 4)
check("z15 H6-Schema reproduziert: kappa = 29,126, a/l_P = 10,79", abs(k_ref - Z15["h6_schema"]["kappa"]) < 1e-3)


def fixpunkt(gfun, dims, n_scal, start=10.8):
    alp = start
    for _ in range(40):
        g = gfun(MPL / alp)
        alp_new = 2 * np.sqrt(kappa_tot(list(zip(dims, g)), n_scal))
        if abs(alp_new - alp) < 1e-6:
            break
        alp = alp_new
    return alp, g


def g_sm(mu):
    a = sm_at(mu)
    return np.sqrt(4 * np.pi / a)


alp_sm, g_sm_fp = fixpunkt(g_sm, [1, 3, 8], 4)
i_safe = np.where(safe_b)[0]
i_rep = i_safe[np.argmin(np.abs(kap_all[i_safe] - 1.3))] if len(i_safe) else i_best
n_sc10 = float(10 * n10[i_rep] + 32 * n16[i_rep] + 45 * n45[i_rep])


def g_so10(mu):
    sel = np.zeros(len(CNT), bool); sel[i_rep] = True
    return np.array([np.sqrt(4 * np.pi / so10_at(mu, sel)[0])])


alp_10, g10_fp = fixpunkt(g_so10, [45], n_sc10)
alp_10nab, _ = fixpunkt(g_so10, [25], n_sc10)
n_scPS = float(sum(c * {"R": 1, "C": 2}[PS_COMP[p][1]] * int(np.prod(PS_COMP[p][0])) for p, c in zip(PAR, ex[1])))
g_ps = np.sqrt(4 * np.pi / A_ex[0])
alp_ps, _ = fixpunkt(lambda mu: np.array([g_ps, g_ps, g_ps]), [15, 3, 3], n_scPS)
alp_ind = {lab: 2 * np.sqrt(kappa_tot([(1, g), (3, g), (8, g)], 4)) for lab, g in
           (("g_ind naiv (1/a)", ind["mu = 1/a"]["g_naiv"]), ("g_ind ein Geschmack (1/a)", ind["mu = 1/a"]["g_ein"]))}
say(f"   SM (12 Eichbosonen, g_i aus 2L-RG an 1/a): Fixpunkt a/l_P = {alp_sm:.2f}, g = {np.round(g_sm_fp, 3).tolist()} (z15 nahm 0,55).")
say(f"   PS-Beispiel 3C (21 Eichbosonen, g = {g_ps:.3f}, {n_scPS:.0f} reelle Skalare): a/l_P = {alp_ps:.2f}.")
say(f"   SO(10)-Lesart (45 Eichbosonen, Inhalt {dict(zip(PAR, CNT[i_rep].astype(int).tolist()))}, {n_sc10:.0f} reelle Skalare):"
    f" a/l_P = {alp_10:.2f} (nichtabelsch gezählt 25 Kopien: {alp_10nab:.2f}), g_U = {g10_fp[0]:.3f}.")
for k, v in alp_ind.items():
    say(f"   Mit {k} statt RG-g (SM-Inhalt): a/l_P = {v:.2f}.")
check("z15-Fixpunkt mit SM-RG-Kopplungen liegt im Fenster (Konsistenz der z15-Annahme g = 0,55)", 1.65 <= alp_sm <= 14.27)
check("SO(10)-Lesart: Fixpunkt a/l_P liegt über dem z15-Fenster (45 Eichkanten) — Bedingung an a/l_P",
      alp_10 > 14.27)
befund(f"Selbstkonsistentes a/l_P: SM {alp_sm:.1f}, PS-Beispiel {alp_ps:.1f}, SO(10) {alp_10:.1f} (nichtabelsch {alp_10nab:.1f});"
       f" mit y1-induziertem g fiele es auf {min(alp_ind.values()):.1f}-{max(alp_ind.values()):.1f}.")
ZAHLEN["frage3"] = {"zweige_2L": q3, "scan_supk_sicher_zentral": int(safe_c.sum()), "scan_supk_sicher_band": int(safe_b.sum()),
                    "scan_band_sicher_kappa_1_1.7": int((safe_b & pin).sum()),
                    "tau_sicher_bereich": [float(tau_all[safe_c].min()), float(tau_all[safe_c].max())] if safe_c.any() else None,
                    "alpha_U_inv_1a": {k: [float(v.min()), float(v.max())] for k, v in aU_a.items()},
                    "y1_kappa_N32": k32, "y1_fit": [float(slope), float(icpt)], "induziert": ind,
                    "induziert_obergrenze": ind_max, "verhaeltnis_min": aU_min / ind_max,
                    "a_lP_fixpunkt": {"SM": alp_sm, "PS_beispiel": alp_ps, "SO10": alp_10, "SO10_nichtabelsch": alp_10nab, **alp_ind}}

# ============================================================================ Teil 5: Frage 4
say("\nTeil 5 (Frage 4): Urteil über die Lesarten")
say(f"   L1  'induziert an 1/a, gleich, H6-Inhalt (SM + importierter Higgs)': fällt. Spreizung {dmin:.1f}-"
    f"{max(v['D_2L'] for v in q1.values()):.1f} an 1/a (2L); auch E8-Fragmente ohne neue Eichbosonen reichen nur bis {e(muB, 1)} GeV.")
say(f"   L2  'PS-Zwischenstufe, Gleichheit genau an 1/a': fällt mit vollständigen Multipletts (M_U <= {e(mu_maxA, 1)}). Überlebt nur als")
say(f"       Anpassung mit gespaltenen Multipletts (Farboktett + Farbtripletts bei TeV ... 1e5 GeV bei kappa ~ 1; protonsicher nur")
say(f"       bis kappa = {kmax_safe:.2f} mit ~8 Farbskalaren bei ~1 TeV), braucht 21 PS-Eichbosonen,")
say(f"       die der 3+2-Split von H6 nicht liefert, und gibt keinen sichtbaren Protonzerfall (tau_p ~ {e(tau_ex, 0)} a). Ihr Wert")
say(f"       alpha_U^-1(1/a) = {A_ex[0]:.1f} widerspricht der y1-Induktion (<= {ind_max:.0f}) ebenso wie in L3.")
say(f"   L3  'SO(10) unterhalb, eine Kopplung bis 1/a': Gleichheit automatisch, aber nicht induziert: RG alpha_U^-1(1/a) >= {aU_min:.0f}")
say(f"       gegen y1-induziert <= {ind_max:.0f}. Überlebt nur als nackte Gitterkopplung mit SO(10)-Symmetrie (45 Eichbosonen, nicht")
say(f"       aus dem 3+2-Split). Die v266-Inhalte sind bei zwei Schleifen ausgeschlossen; E8-Inhalte mit großem M_U überleben,")
say(f"       z. B. {best_cnt} mit tau_p = {e(q3['größtes M_U (3A)']['tau_p'], 1)} a (kappa {q3['größtes M_U (3A)']['kappa']:.2f}).")
say(f"       Mit kappa in [1; 1,7] (PS.KAPPA.01) bleiben {int((safe_b & pin).sum())} Inhalte im x3-Band ({int((safe_c & pin).sum())} zentral).")
say(f"   Vorhersage, falls L3 mit kappa in [1; 1,7]: tau_p = {e(tau_all[safe_b & pin].min(), 1)} ... {e(tau_all[safe_b & pin].max(), 1)} a"
    f" (Hyper-K-Reichweite {e(TAU_HK, 1)} a) und a/l_P = {alp_10nab:.0f}-{alp_10:.0f} aus z15 mit SO(10)-Eichkanten." if (safe_b & pin).any() else
    "   Keine Vorhersage im kappa-Band.")
dt = time.time() - T0
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)}). Laufzeit {dt:.0f} s.")

URTEIL = (f"H6s Lesart 'Eichkopplungen an der Gitterskala induziert und gleich' fällt in der gebauten Form: Mit SM-Inhalt und "
          f"importiertem Higgs ist die Spreizung an 1/a {dmin:.1f}-{max(v['D_2L'] for v in q1.values()):.1f} (2 Schleifen, a/l_P 1,65-14,3), "
          f"und ohne neue Eichbosonen erreichen E8-Fragmente oberhalb 1 TeV die Gleichheit nur bis {e(muB, 1)} GeV. Pati-Salam mit "
          f"vollständigen E8-Multipletts vereinigt höchstens bei {e(mu_maxA, 1)} GeV; Gleichheit genau an 1/a gelingt nur mit "
          f"gespaltenen Multipletts (bei kappa ~ 1 leichte Farbtripletts; protonsicher nur bis kappa = {kmax_safe:.2f} mit ~8 Farbskalaren bei ~1 TeV) und dann ohne "
          f"sichtbaren Protonzerfall. In der SO(10)-Lesart ist Gleichheit automatisch, aber der Wert widerspricht der y1-Induktion "
          f"(alpha_U^-1(1/a) >= {aU_min:.0f} gegen <= {ind_max:.0f}); sie überlebt nur als nackte SO(10)-Gitterkopplung, was 45 statt 12 "
          f"Eichbosonen verlangt (z15-Fixpunkt a/l_P ~ {alp_10:.0f}, außerhalb des Fensters), und die v266-Inhalte sind bei zwei Schleifen "
          f"von Super-K ausgeschlossen; mit kappa in [1; 1,7] überleben nur E8-Inhalte mit tau_p <= {e(tp.max(), 1)} a, die Hyper-K vollständig "
          f"entscheidet. sin^2 theta_W = 3/8 wird so von SO(10) geerbt, nicht von H6 induziert.")
json.dump({"checks": CHECKS, "befunde": BEFUNDE, "urteil": URTEIL, "zahlen": ZAHLEN},
          open(HERE / "z25_ergebnis.json", "w"), ensure_ascii=False, indent=1, default=float)
assert all(CHECKS.values()), [k for k, v in CHECKS.items() if not v]
