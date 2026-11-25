"""Z121: Kohärente Komposition -- ist 'Komposition = Viertelschritt der Rekombinationsenergie im viertelgedrehten
Transferrahmen' eine native Vorwärtsregel? Kette, Integrabilität, Kontinuum, drei Nähte (Folgeschritt zu z120).

Leitlinie 'Rekonstruktion statt Auswahl': Energie, Transfer und Uhr kommen aus den nativen Daten (z92: E^nat = L_0^osc + Q^2/2
je Quarter; z99: Gitter Z^2, su(2)_1-Transfer, E^nat einkörperig; z105: Maya-Diagramme, E^nat = Q^2/2 + |lam|; z81/z82:
R = i^(Q_a - Q_b)). Das Gatter von z120 ist Vergleichsziel, keine Prämisse. Firewall DIMENSION.UPLIFT.FIREWALL.01: Die Kette
ist ein Modell aus Nähten an geordneten Lagen, Kontinuum und 4D nur als Lesart. Gelöste Blöcke (Flavour, Massen, alpha, 16,
E8) bleiben unberührt, H6 wird nicht geöffnet, nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt),
L Lesart, N numerisch, O offen.

Bestand: z120 (U = B^+ e^(-i pi H_Q/2) B = P_T - i P_S, Wurzel aus SWAP, Leckfreiheit erzwingt theta = pi/2 = P_sym,
e^(-i pi H_Q/2) = -i R^(-1/2) auf N = 2, freie Energie trivial, cos eps = (1 + cos k)/2), z116/z119 (Transfer-SU(2),
Gamma(pi) = +1 auf dem Singulett, P_ident, ||[Q^2/2, E^+]|| = 1, KMS auf (1,1) = I/4), z106 (Klein-Kette W_12 = (-1)^(N_A N_B),
nicht beobachtbar), z108 (Kompositionssatz, I_AB|C, |P^2|Delta>|^2 = 32 Delta (Delta - 1), gamma = lambda^2/32), z99/z105
(Delta = 1 starr, Doppeltransfer ohne Box-Nachkommen), z110 (su(2)_1).

Teil 1  Die Phase aus der Rekombination. (a) Vier Lesarten der zwei Moden je Naht mit nativer Energie E^nat: Quarterpaar
        (a, b), zwei Stufen eines Fermions, zwei Spezies gleicher Stufe (z116/z119-Dublett), eta-Dublett (a-Teilchen,
        b-Loch); Nullmodenenergien Q = 0, 1, 2 mit Spininhalt, Uhr, Gatter je Lesart. (b) Ladungsgitter: g-Hülle (z99),
        U(2)_1 = u(1) + su(2)_1 im Fockfenster, spinblindes U(1)_1. (c) Winkel: Bei der Vierteldrehung ist jede Phase
        leckfrei (Spektralparameter u = tan(phi/2)); Leckfreiheit = Verträglichkeit [B^+ D B, D] = 0; [U, H_Q] = 0.
        (d) Rekombinationsregel unter P_null: i^(Q_A Q_B) als Wurzel der Klein-Kette, U^2 = -Gamma(pi), Niveau k gibt
        u = tan(pi/(4k)), Vorzeichen der Wurzel = Zeitrichtung. (e) P_ident gegen eine Kette von Nähten.
Teil 2  Die Kette F = (prod_gerade U)(prod_ungerade U). (a) Yang-Baxter (Zopf- und Yang-Form, Kontrollen); jedes
        SU(2)-invariante Zweiplatzgatter ist R(u). (b) Integrabilität: inhomogene Transfermatrizen, Lichtkegelfaktorisierung
        F = t(-1/2) t(1/2)^(-1), t(-1/2) t(1/2) = T^2, [t(u), t(v)] = [t(u), F] = 0 (L = 4..12), lokale Ladungen,
        Niveaustatistik gegen zwei Kontrollen, verborgene Halbschrittsymmetrie. (c) Einmagnonen-Dispersion gegen z120.
        (d) Kontinuum: Destri-de-Vega-Lesart bei Theta = 1/2 (selbstdual), Vakuum, Lücken, Verschränkung, v = 2 tanh(pi Theta).
Teil 3  Drei Nähte: verbundener Dreinahtteil (Chiralität mit festem Koeffizienten), Auswahlregeln (keine 1->2-Amplitude),
        Faktorisierung, radiale Struktur (keine Nachkommen), Vergleich mit lambda (gamma = lambda^2/32).
Teil 4  Bilanz und Urteil.
"""
import itertools
import json
import time
from pathlib import Path

import numpy as np
from scipy.linalg import expm, logm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(121)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def jw(n):
    I_, Zp, sm = np.eye(2), np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [0.0, 0.0]])
    out = []
    for j in range(n):
        m_ = [Zp]*j + [sm] + [I_]*(n - j - 1)
        op = m_[0]
        for M in m_[1:]:
            op = np.kron(op, M)
        out.append(op)
    return out


def nrm(A):
    return float(np.linalg.norm(A))


def bis_auf_phase(A, B):
    """min_phi ||A - e^(i phi) B||."""
    ph = np.vdot(B.ravel(), A.ravel())
    ph = ph/abs(ph) if abs(ph) > 1e-14 else 1.0
    return nrm(A - ph*B)


def schmidt_rang(U4, tol=1e-10):
    """Operator-Schmidt-Rang eines 4x4-Operators auf C^2 (x) C^2 (Basis |m>|m'>, Index 2 m + m')."""
    R = U4.reshape(2, 2, 2, 2).transpose(0, 2, 1, 3).reshape(4, 4)
    return int(np.sum(np.linalg.svd(R, compute_uv=False) > tol))


SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)
PSWAP = np.eye(4)[[0, 2, 1, 3]].astype(complex)                     # Spintausch auf C^2 (x) C^2
PT_, PS_ = (np.eye(4) + PSWAP)/2, (np.eye(4) - PSWAP)/2


def R_zopf(u):
    """Zopfform R(u) = (1 + i u P)/(1 + i u); u = 1 ist das Gatter von z120, u = oo der Tausch."""
    return (np.eye(4) + 1j*u*PSWAP)/(1 + 1j*u)


# =================================================================== Teil 1
say("Teil 1: Die Phase aus der Rekombination")
say("  (a) Lesarten der zwei Moden je Naht, native Energie E^nat (z92/z105)")
ALPHA = {"a": 0.25, "b": 0.75}


def e_nat_mode(q, x, typ):
    """E^nat = L_0 - (alpha - 1/2) Q einer Anregung: Teilchen bei x > 0 oder Loch bei x < 0 im Quarter q."""
    al = ALPHA[q]
    return (x - (al - 0.5)) if typ == "p" else (-x + (al - 0.5))


def l0_mode(x, typ):
    return x if typ == "p" else -x


LESARTEN = {
    "(a) Quarterpaar: Teilchen in a (1/4) und b (3/4)": [("a", 0.25, "p", "a"), ("b", 0.75, "p", "b")],
    "(b) zwei Stufen eines Fermions: a bei 1/4 und 5/4": [("a", 0.25, "p", "a"), ("a", 1.25, "p", "a")],
    "(c) zwei Spezies gleicher Stufe (Dublett z116/z119)": [("a", 0.25, "p", "a1"), ("a", 0.25, "p", "a2")],
    "(d) eta-Dublett/Fermipaar: a-Teilchen 1/4, b-Loch -1/4": [("a", 0.25, "p", "a"), ("b", -0.25, "h", "b")],
}
c4 = jw(4)
n4 = [x.T @ x for x in c4]
QA4, QB4 = n4[0] + n4[1], n4[2] + n4[3]
PHYS = [s_ for s_ in range(16) if np.isclose(QA4[s_, s_], 1) and np.isclose(QB4[s_, s_], 1)]
JP4 = c4[0].T @ c4[2] + c4[1].T @ c4[3]
JX4 = (JP4 + JP4.T)/2
B4 = expm(-1j*np.pi*JX4/2)
D4 = n4[0] @ n4[1] + n4[2] @ n4[3]                                   # Doppelbesetzungen (Doublonen)
HQ4 = (QA4 @ QA4 + QB4 @ QB4)/2
ZIEL = PT_ - 1j*PS_


def gatter(H, B=B4):
    return B.conj().T @ expm(-1j*np.pi*H/2) @ B


lesart_daten, ok_lesarten = {}, True
for name, moden in LESARTEN.items():
    eps = [e_nat_mode(q, x, t) for q, x, t, _ in moden]
    felder = sorted(set(f for *_, f in moden))
    # Einzelnaht: Zustände (n1, n2)
    zust = {}
    for n1, n2 in itertools.product((0, 1), repeat=2):
        occ = (n1, n2)
        E = sum(e*o for e, o in zip(eps, occ))
        Qf = {f: sum((1 if t == "p" else -1)*o for (q, x, t, ff), o in zip(moden, occ) if ff == f) for f in felder}
        E0 = sum(v*v/2 for v in Qf.values())
        L0 = sum(l0_mode(x, t)*o for (q, x, t, _), o in zip(moden, occ))
        Qa = sum((1 if t == "p" else -1)*o for (q, x, t, _), o in zip(moden, occ) if q == "a")
        Qb = sum((1 if t == "p" else -1)*o for (q, x, t, _), o in zip(moden, occ) if q == "b")
        zust[occ] = {"E_nat": E, "Nullmodenteil": E0, "|lam|": E - E0, "R = i^(Q_a - Q_b)": complex(1j**(Qa - Qb)),
                     "e^(2 pi i L_0)": complex(np.exp(2j*np.pi*L0))}
        ok_lesarten &= abs(E - E0 - round(E - E0)) < 1e-12 and E - E0 > -1e-12           # |lam| ganz, >= 0 (Maya)
        ok_lesarten &= abs(zust[occ]["R = i^(Q_a - Q_b)"] - zust[occ]["e^(2 pi i L_0)"]) < 1e-12  # Uhr z81
    dublett = abs(zust[(1, 0)]["E_nat"] - zust[(0, 1)]["E_nat"]) < 1e-12
    U_H = zust[(1, 1)]["E_nat"] - zust[(1, 0)]["E_nat"] - zust[(0, 1)]["E_nat"]
    U_H0 = zust[(1, 1)]["Nullmodenteil"] - zust[(1, 0)]["Nullmodenteil"] - zust[(0, 1)]["Nullmodenteil"]
    # Zwei Nähte: E^nat einkörperig, Transfer J^+ = sum_m c+_(A,m) c_(B,m) (Nullmode des su(2)_1-Stroms)
    H_nat = sum(eps[m]*(n4[m] + n4[2 + m]) for m in range(2))
    U_nat = gatter(H_nat)
    Up = U_nat[np.ix_(PHYS, PHYS)]
    H_null = sum(((sum(n4[2*s_ + m]*(1 if moden[m][2] == "p" else -1) for m in range(2) if moden[m][3] == f)) @
                  (sum(n4[2*s_ + m]*(1 if moden[m][2] == "p" else -1) for m in range(2) if moden[m][3] == f)))/2
                 for s_ in range(2) for f in felder)
    U_null = gatter(H_null)[np.ix_(PHYS, PHYS)]
    lesart_daten[name] = {
        "E^nat (n1, n2) = (0,0), (1,0), (0,1), (1,1)": [zust[o]["E_nat"] for o in ((0, 0), (1, 0), (0, 1), (1, 1))],
        "Nullmodenteil": [zust[o]["Nullmodenteil"] for o in ((0, 0), (1, 0), (0, 1), (1, 1))],
        "Dublett entartet": dublett, "Hubbard U_H = E(1,1) - E(1,0) - E(0,1)": U_H, "U_H des Nullmodenteils": U_H0,
        "||[E^nat, J^+]||": nrm(H_nat @ JP4 - JP4 @ H_nat),
        "||[Nullmodenteil, J^+]||_op": float(np.linalg.norm(H_null @ JP4 - JP4 @ H_null, 2)),
        "max |Matrixelement| [Nullmodenteil, J^+]": float(np.abs(H_null @ JP4 - JP4 @ H_null).max()),
        "Schmidt-Rang U(E^nat)": schmidt_rang(Up),
        "U(E^nat) reine Phase (Abstand)": bis_auf_phase(Up, np.eye(4)),
        "Schmidt-Rang U(Nullmodenteil)": schmidt_rang(U_null), "U(Nullmodenteil) = P_T - i P_S (Abstand)": bis_auf_phase(U_null, ZIEL)}
    say(f"   {name}: E^nat = {lesart_daten[name]['E^nat (n1, n2) = (0,0), (1,0), (0,1), (1,1)']}, Nullmodenteil "
        f"{lesart_daten[name]['Nullmodenteil']}, Dublett entartet {dublett}, U_H {U_H:g} (Nullmodenteil {U_H0:g}), "
        f"Schmidt-Rang U(E^nat) {schmidt_rang(Up)}, U(Nullmodenteil) {schmidt_rang(U_null)}")
ZAHLEN["lesarten"] = lesart_daten
LA, LB, LC, LD = (lesart_daten[k] for k in LESARTEN)
check("[E] Vier Lesarten der zwei Moden je Naht, jede mit der nativen Energie E^nat = L_0 - (alpha - 1/2) Q (z92/z105) und dem "
      "nativen Transfer J^+ = sum_m c+_(A,m) c_(B,m) (Nullmode des su(2)_1-Stroms, z99): E^nat ist einkörperig, die Zerlegung "
      "E^nat = Q^2/2 + |lam| ist in jeder Lesart ganzzahlig (Maya), und R = i^(Q_a - Q_b) = e^(2 pi i L_0) (z81) gilt auf allen "
      f"Zuständen ({ok_lesarten}). (a) Quarterpaar: E^nat = {LA['E^nat (n1, n2) = (0,0), (1,0), (0,1), (1,1)']}, Dublett "
      f"entartet, U_H = {LA['Hubbard U_H = E(1,1) - E(1,0) - E(0,1)']:g}; (b) zwei Stufen eines Fermions: "
      f"{LB['E^nat (n1, n2) = (0,0), (1,0), (0,1), (1,1)']}, nicht entartet, U_H = 0, aber der Nullmodenteil "
      f"{LB['Nullmodenteil']} ist spinblind mit U_H = {LB['U_H des Nullmodenteils']:g}; (c) zwei Spezies gleicher Stufe: "
      f"{LC['E^nat (n1, n2) = (0,0), (1,0), (0,1), (1,1)']}; (d) eta-Dublett: {LD['E^nat (n1, n2) = (0,0), (1,0), (0,1), (1,1)']}. "
      "In allen vier Lesarten vertauscht E^nat mit dem Transfer (max "
      f"{max(d['||[E^nat, J^+]||'] for d in lesart_daten.values()):.0e}), also ist U = B^+ e^(-i pi E^nat/2) B ein Produkt von "
      f"Einzelnahtphasen: Operator-Schmidt-Rang {[d['Schmidt-Rang U(E^nat)'] for d in lesart_daten.values()]}, in (a, c, d) "
      f"sogar eine reine Phase (Abstand {max(LA['U(E^nat) reine Phase (Abstand)'], LC['U(E^nat) reine Phase (Abstand)'], LD['U(E^nat) reine Phase (Abstand)']):.0e}). "
      "Nur der Nullmodenteil der Lesart (b), also E^nat ohne L_0^osc = |lam|, bricht die Transfer-SU(2) ([Q^2/2, J^+] = "
      f"J^+ (Q_A - Q_B + 1): größtes Matrixelement {LB['max |Matrixelement| [Nullmodenteil, J^+]']:g} wie z119, Operatornorm "
      f"{LB['||[Nullmodenteil, J^+]||_op']:.6f} = sqrt 2) und gibt das z120-Gatter: Schmidt-Rang "
      f"{LB['Schmidt-Rang U(Nullmodenteil)']}, Abstand zu P_T - i P_S {LB['U(Nullmodenteil) = P_T - i P_S (Abstand)']:.0e}; der "
      f"Nullmodenteil von (a, c, d) bleibt trivial (Rang {[LA['Schmidt-Rang U(Nullmodenteil)'], LC['Schmidt-Rang U(Nullmodenteil)'], LD['Schmidt-Rang U(Nullmodenteil)']]})",
      ok_lesarten and all(d["Schmidt-Rang U(E^nat)"] == 1 for d in lesart_daten.values())
      and max(d["||[E^nat, J^+]||"] for d in lesart_daten.values()) < 1e-12
      and max(LA["U(E^nat) reine Phase (Abstand)"], LC["U(E^nat) reine Phase (Abstand)"], LD["U(E^nat) reine Phase (Abstand)"]) < 1e-12
      and LA["Dublett entartet"] and LC["Dublett entartet"] and LD["Dublett entartet"] and not LB["Dublett entartet"]
      and all(abs(d["Hubbard U_H = E(1,1) - E(1,0) - E(0,1)"]) < 1e-12 for d in lesart_daten.values())
      and LB["Schmidt-Rang U(Nullmodenteil)"] == 4 and LB["U(Nullmodenteil) = P_T - i P_S (Abstand)"] < 1e-12
      and abs(LB["U_H des Nullmodenteils"] - 1) < 1e-12 and abs(LB["||[Nullmodenteil, J^+]||_op"] - np.sqrt(2)) < 1e-12
      and abs(LB["max |Matrixelement| [Nullmodenteil, J^+]"] - 1) < 1e-12
      and all(d["Schmidt-Rang U(Nullmodenteil)"] == 1 for d in (LA, LC, LD)))
befund("E: In keiner nativen Lesart ist das Gatter verschränkend. E^nat ist einkörperig (z99), der Transfer erhält die Modenmarke, "
       "also vertauscht E^nat mit der Transfer-SU(2), und U ist ein Produkt von Einzelnahtphasen. Die zwei Moden des z120-Modells "
       "sind am ehesten zwei Stufen eines Fermions (b): Dort ist der Ladungsterm Q^2/2 spinblind (Q zählt beide Stufen), aber "
       "L_0^osc = |lam| hebt die Entartung auf und macht E^nat additiv. Das z120-Gatter entsteht genau dann, wenn die Phase nur den "
       "Nullmodenteil sieht: Prämisse P_null ('die Komposition sieht nur das Ladungsgitter'). Die Uhr R = i^(Q_a - Q_b) ist "
       "oszillatorblind, weil e^(2 pi i n) = 1; ihr Viertelschritt i^n ist es nicht, P_null folgt also nicht aus der Uhr", True)

say("  (b) Ladungsgitter: g-Hülle, U(2)_1 im Fockfenster, spinblindes U(1)_1")
g_daten, ok_g = {}, True
for g in (-2, -1, 0, 1, 2):
    G = np.array([[1, g], [g, 1]], float)
    ev = np.linalg.eigvalsh(G)
    H_g = sum((n4[2*s_] @ n4[2*s_] + 2*g*n4[2*s_] @ n4[2*s_ + 1] + n4[2*s_ + 1] @ n4[2*s_ + 1])/2 for s_ in range(2))
    Ug = gatter(H_g)[np.ix_(PHYS, PHYS)]
    g_daten[str(g)] = {"Eigenwerte Gram": ev.tolist(), "psd": bool(ev.min() > -1e-12),
                       "U_H": float(g), "Abstand zu R(u=1) = P_T - i P_S": bis_auf_phase(Ug, R_zopf(1.0)),
                       "Abstand zu R(u=-1)": bis_auf_phase(Ug, R_zopf(-1.0)), "reine Phase": bis_auf_phase(Ug, np.eye(4))}
ok_g = (g_daten["0"]["reine Phase"] < 1e-12 and g_daten["1"]["Abstand zu R(u=1) = P_T - i P_S"] < 1e-12
        and g_daten["-1"]["Abstand zu R(u=-1)"] < 1e-12 and not g_daten["2"]["psd"] and not g_daten["-2"]["psd"]
        and all(g_daten[k]["psd"] for k in ("-1", "0", "1")) and abs(min(g_daten["1"]["Eigenwerte Gram"])) < 1e-12)
ZAHLEN["g_huelle"] = g_daten
# Fockfenster: zwei Spezies eines Quarters (U(2)_1) und eine Spezies (U(1)_1), Moden x in {-7/4, -3/4, 1/4, 5/4}
XS = [-1.75, -0.75, 0.25, 1.25]
c8 = jw(8)
n8 = [x.T @ x for x in c8]
vak8 = np.array([1 if x < 0 else 0 for x in XS]*2)
Enat8 = sum((XS[i % 4] + 0.25)*(n8[i] - vak8[i]*np.eye(256)) for i in range(8))
Q1 = sum(n8[i] - vak8[i]*np.eye(256) for i in range(4))
Q2 = sum(n8[4 + i] - vak8[4 + i]*np.eye(256) for i in range(4))
Sp = sum(c8[i].T @ c8[4 + i] for i in range(4))
Sz8 = (Q1 - Q2)/2
S2_8 = Sp @ Sp.T + Sz8 @ Sz8 - Sz8
dE, dQ = np.diag(Enat8), np.diag(Q1 + Q2)
fenster = {}
for Q, Es in ((0, (0.0, 1.0)), (1, (0.5, 1.5)), (2, (1.0, 2.0))):
    for E in Es:
        sel = [s_ for s_ in range(256) if abs(dQ[s_] - Q) < 1e-9 and abs(dE[s_] - E) < 1e-9]
        if not sel:
            continue
        cas = np.linalg.eigvalsh(S2_8[np.ix_(sel, sel)])
        spins = sorted(float(round((-1 + np.sqrt(1 + 4*max(c_, 0)))/2, 6)) for c_ in cas)
        fenster[f"Q = {Q}, E^nat = {E:g}"] = spins
# Sugawara: Primärzustände E = Q^2/4 + j(j+1)/3
sug = all(abs(min(dE[s_] for s_ in range(256) if abs(dQ[s_] - Q) < 1e-9) - (Q*Q/4 + j*(j + 1)/3)) < 1e-12
          for Q, j in ((0, 0), (1, 0.5), (2, 0)))
e1 = {Q: sorted(round(float(dE[s_]), 6) for s_ in range(256) if abs(np.diag(Q1)[s_] - Q) < 1e-9
                and abs(np.diag(Q2)[s_]) < 1e-9)[:2] for Q in (0, 1, 2)}
ZAHLEN["fockfenster"] = {"U(2)_1: Spins je (Q, E^nat)": fenster, "eine Spezies: zwei tiefste E^nat je Q": e1,
                         "Sugawara Q^2/4 + j(j+1)/3 auf den Primärzuständen": sug}
ok_fenster = (fenster["Q = 1, E^nat = 0.5"] == [0.5, 0.5] and fenster["Q = 2, E^nat = 1"] == [0.0]
              and fenster["Q = 2, E^nat = 2"] == [0.0, 1.0, 1.0, 1.0] and fenster["Q = 0, E^nat = 1"] == [0.0, 1.0, 1.0, 1.0]
              and sug and e1[1] == [0.5, 1.5] and e1[2][0] == 2.0)
check("[E] Ladungsgitter. g-Hülle (z99) je Naht für die zwei Komponenten, Gram [[1, g], [g, 1]]: positiv semidefinit nur für "
      f"g = -1, 0, 1 (g = +-2: Eigenwerte {g_daten['2']['Eigenwerte Gram']}); die Hubbard-Krümmung ist U_H = g. g = 0 (die "
      f"einzige lokale, nichtentartete Hülle) gibt eine reine Phase ({g_daten['0']['reine Phase']:.0e}), g = 1 genau das "
      f"z120-Gatter R(1) = P_T - i P_S ({g_daten['1']['Abstand zu R(u=1) = P_T - i P_S']:.0e}) mit dem Nullvektor (1, -1) "
      f"(Gram-Eigenwert {min(g_daten['1']['Eigenwerte Gram']):.0e}), g = -1 das inverse Gatter R(-1) "
      f"({g_daten['-1']['Abstand zu R(u=-1)']:.0e}). Fockfenster (Moden -7/4, -3/4, 1/4, 5/4 je Spezies, 256 Zustände): zwei "
      f"Spezies eines Quarters bilden U(2)_1, Spins {fenster}; die Primärzustände liegen bei Q^2/4 + j(j+1)/3 ({sug}): "
      "Q = 1 Dublett bei 1/2, Q = 2 Singulett bei 1, Triplett + Singulett bei 2. Also 'zwei Spezies' = 'Sugawara-Zerlegung "
      f"U(2)_1 = u(1) + su(2)_1' (dieselbe Theorie). Eine Spezies (U(1)_1) hat je Q nur einen tiefsten Zustand: Q = 1 bei {e1[1]}, "
      f"Q = 2 bei {e1[2]}, also Q^2/2 ohne Dublett",
      ok_g and ok_fenster)
befund("E: Native Ladungsgitter der Naht sind U(1)_1 (ein Fermion, kein Dublett) und U(2)_1 (zwei Spezies = Sugawara-Zerlegung; "
       "Dublett bei 1/2, Singulett bei 1, U_H = 0, Gatter trivial). Ein spinblindes U(1)_1 mit Dublett ist kein natives Gitter: "
       "Es ist die entartete Hülle g = 1 von z99 (Relativladung Nullvektor, Spinwechsel mit h = 0), die z99 als entartet "
       "ausgeschlossen hat. Optionen 2 und 3 der Aufgabe sind dieselbe und geben ein triviales Gatter; Option 1 ist bedingt "
       "(Prämisse P_null, gleichwertig g = 1)", True)

say("  (c) Die Winkel: Vierteldrehung erzwungen, Phase = Spektralparameter")
PPHYS = np.zeros((16, 16))
PPHYS[PHYS, PHYS] = 1
phis = np.linspace(0, np.pi, 13)
leck_phi, abst_phi = [], []
for ph in phis:
    Uf = gatter(2*ph/np.pi*HQ4)
    leck_phi.append(nrm((np.eye(16) - PPHYS) @ Uf @ PPHYS))
    u_ = np.tan(ph/2) if ph < np.pi - 1e-12 else None
    abst_phi.append(bis_auf_phase(Uf[np.ix_(PHYS, PHYS)], R_zopf(u_) if u_ is not None else PSWAP))
thetas = np.linspace(0, np.pi, 181)
phi_gen = 1.0
leckfrei, vertraeglich = [], []
for th in thetas:
    Bt = expm(-1j*th*JX4)
    Ut = Bt.conj().T @ expm(-1j*phi_gen*HQ4) @ Bt
    if nrm((np.eye(16) - PPHYS) @ Ut @ PPHYS) < 1e-10:
        leckfrei.append(round(float(th/np.pi), 4))
    Dt = Bt.conj().T @ D4 @ Bt
    if nrm(Dt @ D4 - D4 @ Dt) < 1e-10:
        vertraeglich.append(round(float(th/np.pi), 4))
U1 = gatter(HQ4)
komm_HQ = nrm(U1 @ HQ4 - HQ4 @ U1)
rho_kms = expm(-4*np.pi*HQ4)
rho_kms /= np.trace(rho_kms)
kms_inv = nrm(U1 @ rho_kms @ U1.conj().T - rho_kms)
rho11 = rho_kms[np.ix_(PHYS, PHYS)]/np.trace(rho_kms[np.ix_(PHYS, PHYS)])
ZAHLEN["winkel"] = {"Leckage bei theta = pi/2 über phi in [0, pi]": max(leck_phi),
                    "Abstand U(phi) zu R(tan(phi/2))": max(abst_phi), "leckfreie theta/pi bei phi = 1": leckfrei,
                    "theta/pi mit [B^+ D B, D] = 0": vertraeglich, "||[U, H_Q]||": komm_HQ,
                    "KMS(beta = 4 pi) invariant unter U": kms_inv, "KMS auf (1,1) minus I/4": nrm(rho11 - np.eye(4)/4)}
check(f"[E] Bei der Vierteldrehung theta = pi/2 ist das Gatter für jede Phase phi leckfrei (max Leckage {max(leck_phi):.0e}, 13 Werte "
      f"in [0, pi]) und gleich R(u) = (1 + i u P)/(1 + i u) mit u = tan(phi/2) (Abstand bis auf Phase {max(abst_phi):.0e}; phi = pi "
      "gibt den Tausch, u = oo): Die Phase ist der Spektralparameter der XXX-Familie, nicht durch Leckfreiheit festgelegt. Die "
      f"Drehung schon: bei der generischen Phase phi = 1 sind genau theta/pi = {leckfrei} leckfrei, und das sind genau die Winkel mit "
      f"[B^+ D B, D] = 0 ({vertraeglich}, D = Zahl der Doppelbesetzungen): Leckfreiheit heißt, die gedrehte Ladungsenergie ist mit "
      f"der nativen verträglich. Das Gatter vertauscht mit H_Q (||[U, H_Q]|| = {komm_HQ:.0e}) und lässt den KMS-Zustand bei "
      f"beta = 4 pi fest ({kms_inv:.0e}); auf (1,1) ist dieser I/4 (Abstand {nrm(rho11 - np.eye(4)/4):.0e}, z119): Das Gatter ist "
      "im Zustand unsichtbar, es wirkt nur auf reine Zustände und Dynamik",
      max(leck_phi) < 1e-12 and max(abst_phi) < 1e-12 and leckfrei == [0.0, 0.5, 1.0] and vertraeglich == [0.0, 0.5, 1.0]
      and komm_HQ < 1e-12 and kms_inv < 1e-14 and nrm(rho11 - np.eye(4)/4) < 1e-12)
befund("E: Von den zwei Winkeln ist nur die Drehung erzwungen (P_sym, z120), als Verträglichkeit [B^+ D B, D] = 0. Die Phase "
       "phi ist frei und durchläuft die ganze integrable Familie R(u), u = tan(phi/2). Der Wert u = 1 kommt allein aus dem "
       "Viertelschritt phi = pi/2 (Prämisse P_viertel)", True)

say("  (d) Rekombinationsregel unter P_null: Wurzel des Tauschs")
N4 = QA4 + QB4
sek2 = [s_ for s_ in range(16) if np.isclose(N4[s_, s_], 2)]
QAQB = QA4 @ QB4
P2 = np.zeros((16, 16))
P2[sek2, sek2] = 1
inter = nrm(P2 @ (expm(-1j*np.pi*HQ4/2) + expm(1j*np.pi/2*QAQB)) @ P2)              # e^(-i pi H_Q/2) = -i^(Q_A Q_B)
intra = nrm(P2 @ (expm(-1j*np.pi*HQ4/2) - np.exp(-1j*np.pi/2)*expm(-1j*np.pi/2*D4)) @ P2)  # = e^(-i pi/2) (-i)^D
W12 = np.diag([(-1)**int(round(QA4[s_, s_]*QB4[s_, s_])) for s_ in range(16)])
wurzel_kette = nrm(P2 @ (expm(1j*np.pi/2*QAQB) @ expm(1j*np.pi/2*QAQB) - W12) @ P2)
EPEM = JP4 + JP4.T
GAM = expm(-1j*np.pi*N4/2) @ expm(1j*np.pi/2*EPEM)
tausch_fehler = max(nrm(GAM @ c4[i] @ GAM.conj().T - c4[(i + 2) % 4]) for i in range(4))
gam_B = nrm(GAM - expm(-1j*np.pi*N4/2) @ np.linalg.matrix_power(B4.conj().T, 2))
Gp = GAM[np.ix_(PHYS, PHYS)]
U_quadrat = nrm(R_zopf(1.0) @ R_zopf(1.0) + Gp) + bis_auf_phase(U1[np.ix_(PHYS, PHYS)] @ U1[np.ix_(PHYS, PHYS)], Gp)
gedrehte_kette = bis_auf_phase((B4.conj().T @ W12 @ B4)[np.ix_(PHYS, PHYS)], Gp)
niveau = {}
for k in (1, 2, 3):
    Uk = (B4.conj().T @ expm(1j*np.pi/(2*k)*QAQB) @ B4)[np.ix_(PHYS, PHYS)]
    niveau[k] = bis_auf_phase(Uk, R_zopf(np.tan(np.pi/(4*k))))
zeit = bis_auf_phase(gatter(-HQ4)[np.ix_(PHYS, PHYS)], np.linalg.inv(ZIEL)) + bis_auf_phase(R_zopf(-1.0), np.linalg.inv(R_zopf(1.0)))
QAm, QBm = np.diag(QA4), np.diag(QB4)
Rd = np.array([1j**int(round(QAm[s_] - QBm[s_])) for s_ in sek2])
uhr = nrm(expm(-1j*np.pi/2*HQ4)[np.ix_(sek2, sek2)] - (-1j)*np.diag([1 if np.isclose(r_, 1) else -1j for r_ in Rd]))
ZAHLEN["rekombination"] = {"e^(-i pi H_Q/2) + i^(Q_A Q_B) auf N = 2": inter, "intra: e^(-i pi/2)(-i)^D": intra,
                           "(i^(Q_A Q_B))^2 - W_12 auf N = 2": wurzel_kette, "Gamma(pi) tauscht c_A <-> c_B": tausch_fehler,
                           "Gamma(pi) - e^(-i pi N/2) B^(-2)": gam_B, "U^2 + Gamma(pi) auf (1,1)": U_quadrat,
                           "B^+ W_12 B gegen Gamma(pi) auf (1,1)": gedrehte_kette,
                           "Niveau k: Abstand zu R(tan(pi/(4k)))": {str(k): v for k, v in niveau.items()},
                           "Zeitumkehr": zeit, "-i R^(-1/2) (z120)": uhr}
check("[E unter P_null] Auf N = 2 ist e^(-i pi H_Q/2) = -i^(Q_A Q_B) (Rekombination der Nähte A, B, Abweichung "
      f"{inter:.0e}) = e^(-i pi/2) (-i)^D (Rekombination der Komponenten jeder Naht, {intra:.0e}): beide Lesarten der "
      f"Rekombinationsenergie geben dasselbe Gatter. Der Viertelschritt i^(Q_A Q_B) ist eine Wurzel der Klein-Kette W_12 = "
      f"(-1)^(N_A N_B) (z106, {wurzel_kette:.0e}); der kanonische Tausch Gamma(pi) = e^(-i pi N/2) e^(i pi/2 (E^+ + E^-)) bildet "
      f"c_A <-> c_B ab ({tausch_fehler:.0e}) und ist e^(-i pi N/2) B^(-2) ({gam_B:.0e}): die Vierteldrehung ist die Wurzel des Tauschs. "
      f"Gedreht wird die Kette zum Tausch (B^+ W_12 B = Gamma(pi) auf (1,1) bis auf Phase, {gedrehte_kette:.0e}), und U^2 = "
      f"-Gamma(pi) auf (1,1) ({U_quadrat:.0e}): Das Gatter ist die Wurzel des kanonischen Nahttauschs. Ein Gitter auf Niveau k "
      f"(Delta H = Q_A Q_B/k) gibt R(tan(pi/(4k))) ({ {k: f'{v:.0e}' for k, v in niveau.items()} }): u = 1, sqrt 2 - 1, "
      f"tan(pi/12). Die andere Wurzel gibt R(-1) = R(1)^(-1) (Zeitumkehr, {zeit:.0e}); z120: e^(-i pi H_Q/2) = -i R^(-1/2) auf "
      f"N = 2 ({uhr:.0e})",
      inter < 1e-12 and intra < 1e-12 and wurzel_kette < 1e-12 and tausch_fehler < 1e-12 and gam_B < 1e-12
      and U_quadrat < 1e-12 and gedrehte_kette < 1e-12 and max(niveau.values()) < 1e-12 and zeit < 1e-12 and uhr < 1e-12)
befund("C: Unter P_null gilt die Regel exakt und ohne freien Parameter: Komposition = Viertelschritt der "
       "Rekombinationsenergie im viertelgedrehten Rahmen = Wurzel des kanonischen Nahttauschs. Prämissen: P_null (nur das "
       "Ladungsgitter, spinblind auf Niveau 1, gleichwertig g = 1), P_sym (Vierteldrehung, erzwungen als Verträglichkeit), "
       "P_viertel (Viertelschritt = Wurzel des Tauschs; legt u = 1 fest) und ein Vorzeichen der Wurzel (Zeitrichtung). "
       "Neu: Die gedrehte Klein-Kette ist der Tausch, und das Gatter ist seine Wurzel; das Niveau des Gitters bestimmt u", True)

say("  (e) P_ident gegen eine Kette")
c6 = jw(6)
n6 = [x.T @ x for x in c6]
evG = np.linalg.eigvalsh((Gp + Gp.conj().T)/2)
P_id2 = (np.eye(4) + Gp)/2
sing_auf_id = nrm(P_id2 - PS_)
U_auf_id = nrm(U1[np.ix_(PHYS, PHYS)] @ P_id2 + P_id2) + nrm(R_zopf(1.0) @ P_id2 + 1j*P_id2)


def gamma_paar(i, j, cs):
    ns = [x.T @ x for x in cs]
    Ep = sum(cs[2*i + s_].T @ cs[2*j + s_] for s_ in range(2))
    Nij = ns[2*i] + ns[2*i + 1] + ns[2*j] + ns[2*j + 1]
    return expm(-1j*np.pi*Nij/2) @ expm(1j*np.pi/2*(Ep + Ep.T))


PHYS3 = [s_ for s_ in range(64) if all(np.isclose(n6[2*i][s_, s_] + n6[2*i + 1][s_, s_], 1) for i in range(3))]
G_ab, G_bc = gamma_paar(0, 1, c6)[np.ix_(PHYS3, PHYS3)], gamma_paar(1, 2, c6)[np.ix_(PHYS3, PHYS3)]
M_id3 = (G_ab - np.eye(8)).conj().T @ (G_ab - np.eye(8)) + (G_bc - np.eye(8)).conj().T @ (G_bc - np.eye(8))
dim_id3 = int(np.sum(np.linalg.eigvalsh(M_id3) < 1e-10))
ZAHLEN["p_ident"] = {"Gamma(pi) auf (1,1): Eigenwerte": sorted(np.round(evG, 12).tolist()),
                     "P_ident-Projektor = Singulett": sing_auf_id, "U auf dem P_ident-Raum = -i": U_auf_id,
                     "dim P_ident-Raum, drei Nähte (1,1,1)": dim_id3}
check(f"[E] Gamma(pi) hat auf (1,1) die Eigenwerte {sorted(np.round(evG, 12).tolist())} (+1 auf dem Spin-Singulett, z119); "
      f"P_ident (z106) projiziert genau auf das Singulett ({sing_auf_id:.0e}), und dort ist das Gatter die reine Phase -1 "
      f"(z120-Normierung; R(1): -i; {U_auf_id:.0e}). Für drei Nähte (1,1,1) ist der gemeinsame +1-Raum von Gamma(AB) und Gamma(BC) {dim_id3}-dimensional "
      "(Pauli: kein total antisymmetrischer Zustand dreier Dubletts). Gälte P_ident zwischen Kettenplätzen, so wäre die Kette für "
      "L = 2 trivial und für L >= 3 leer",
      sorted(np.round(evG, 9).tolist()) == [-1.0, -1.0, -1.0, 1.0] and sing_auf_id < 1e-12 and U_auf_id < 1e-12 and dim_id3 == 0)
befund("E/C: P_ident (z91/z106: gleiche Nähte, symmetrisierte Komposition) und die Kette schließen sich aus. Die Kette braucht "
       "unterscheidbare, geordnete Nähte (Prämisse P_ord: eine Lagemarke je Naht, die die nativen Nahtdaten nicht tragen; bei "
       "z91 stecken die Lagen in der radialen Leiter jeder Naht). Das z120-Gatter ist also eine Regel zwischen Nähten "
       "verschiedener Lage, nie zwischen identifizierten", True)

# =================================================================== Teil 2
say("\nTeil 2: Die Kette F = (prod_gerade U)(prod_ungerade U)")


def swap_idx(L, x, y):
    idx = np.arange(2**L)
    bx, by = (idx >> x) & 1, (idx >> y) & 1
    return np.where(bx != by, idx ^ ((1 << x) | (1 << y)), idx)


def F_anw(V, L, u=1.0, sws=None):
    """F V mit Gattern R(u) = (1 + i u P)/(1 + i u): erst Bindungen (0,1), (2,3), ..., dann (1,2), ..., (L-1,0)."""
    a, b = 1/(1 + 1j*u), 1j*u/(1 + 1j*u)
    sws = sws or {x: swap_idx(L, x, (x + 1) % L) for x in range(L)}
    for x in range(0, L, 2):
        V = a*V + b*V[sws[x]]
    for x in range(1, L, 2):
        V = a*V + b*V[sws[x]]
    return V


def t_anw(V, L, u, xi, ableitung=False):
    """Inhomogene Transfermatrix t(u) = Tr_0 R_(0,L-1)(u - xi_(L-1)) ... R_(0,0)(u - xi_0), R(w) = w + i P_(0j)."""
    D = 2**L
    idx = np.arange(D)
    V = V.reshape(D, -1).astype(complex)
    out = np.zeros_like(V)
    for a0 in (0, 1):
        X = [np.zeros_like(V), np.zeros_like(V)]
        dX = [np.zeros_like(V), np.zeros_like(V)]
        X[a0] = V.copy()
        for j in range(L):
            w = u - xi[j]
            bj = ((idx >> j) & 1)[:, None]
            PX, PdX = [], []
            for a in (0, 1):
                src = idx ^ ((((idx >> j) & 1) ^ a) << j)
                PX.append(np.where(bj == 0, X[0][src], X[1][src]))
                PdX.append(np.where(bj == 0, dX[0][src], dX[1][src]))
            X, dX = [w*X[a] + 1j*PX[a] for a in (0, 1)], [X[a] + w*dX[a] + 1j*PdX[a] for a in (0, 1)]
        out += (dX if ableitung else X)[a0]
    return out


def T2_anw(V, L):
    D = 2**L
    idx = np.arange(D)
    img = ((idx << 2) | (idx >> (L - 2))) & (D - 1)
    out = np.zeros_like(V)
    out[img] = V
    return out


def sektor(L, nup):
    return np.array([s_ for s_ in range(2**L) if bin(s_).count("1") == nup])


def sektor_ops(L, nup):
    """Basis des S_z-Sektors, S^2 = 3L/4 + sum_(i<j) (P_ij - 1/2) und T^2 als Matrizen im Sektor."""
    D = 2**L
    b_ = sektor(L, nup)
    n_ = len(b_)
    pos_ = -np.ones(D, int)
    pos_[b_] = np.arange(n_)
    S2_ = (3*L/4 - L*(L - 1)/4)*np.eye(n_)
    for i in range(L):
        for j in range(i + 1, L):
            S2_[pos_[swap_idx(L, i, j)[b_]], np.arange(n_)] += 1
    img_ = ((b_ << 2) | (b_ >> (L - 2))) & (D - 1)
    T2_ = np.zeros((n_, n_))
    T2_[pos_[img_], np.arange(n_)] = 1
    return b_, pos_, S2_, T2_


def skalar_rest(A, B):
    c_ = np.vdot(B.ravel(), A.ravel())/np.vdot(B.ravel(), B.ravel())
    return nrm(A - c_*B)/nrm(A), c_


say("  (a) Yang-Baxter")
P12, P23 = np.kron(PSWAP, I2), np.kron(I2, PSWAP)
P13 = P12 @ P23 @ P12
ybe = []
for _ in range(5):
    u_, v_ = rng.normal(size=2)
    R12 = lambda w: (np.eye(8) + 1j*w*P12)/(1 + 1j*w)
    R23 = lambda w: (np.eye(8) + 1j*w*P23)/(1 + 1j*w)
    ybe.append(nrm(R12(u_) @ R23(u_ + v_) @ R12(v_) - R23(v_) @ R12(u_ + v_) @ R23(u_)))
    Y = lambda w, P: w*np.eye(8) + 1j*P
    ybe.append(nrm(Y(u_ - v_, P12) @ Y(u_, P13) @ Y(v_, P23) - Y(v_, P23) @ Y(u_, P13) @ Y(u_ - v_, P12)))
u_, v_ = 0.7, 1.3
kontr_add = nrm((u_*np.eye(8) + 1j*P12) @ ((u_ + v_)*np.eye(8) + 1j*P23) @ (v_*np.eye(8) + 1j*P12)
                - (v_*np.eye(8) + 1j*P23) @ ((u_ + v_)*np.eye(8) + 1j*P12) @ (u_*np.eye(8) + 1j*P23))
kontr_harm = nrm((u_*np.eye(8) + 1j*P12) @ ((u_*v_/(u_ + v_))*np.eye(8) + 1j*P23) @ (v_*np.eye(8) + 1j*P12)
                 - (v_*np.eye(8) + 1j*P23) @ ((u_*v_/(u_ + v_))*np.eye(8) + 1j*P12) @ (u_*np.eye(8) + 1j*P23))
ZZ = np.kron(SZ, SZ)
Rs = lambda w: R_zopf(w) @ expm(0.2j*ZZ)
kontr_zz = nrm(np.kron(Rs(u_), I2) @ np.kron(I2, Rs(u_ + v_)) @ np.kron(Rs(v_), I2)
               - np.kron(I2, Rs(v_)) @ np.kron(Rs(u_ + v_), I2) @ np.kron(I2, Rs(u_)))
gate_u1 = nrm(R_zopf(1.0) - ZIEL) + nrm((np.eye(4) + 1j*PSWAP)/(1 + 1j) - ZIEL)
dual = max(bis_auf_phase(R_zopf(w_) @ R_zopf(1/w_), PSWAP) for w_ in (0.3, 0.7, 1.9, 4.2)) + nrm(R_zopf(1.0) @ R_zopf(1.0) - PSWAP)
Stot = [np.kron(S_, I2) + np.kron(I2, S_) for S_ in (SX, SY, SZ)]
A_kom = np.vstack([np.kron(S_, np.eye(4)) - np.kron(np.eye(4), S_.T) for S_ in Stot])
dim_kom = 16 - np.linalg.matrix_rank(A_kom, tol=1e-10)
ZAHLEN["yang_baxter"] = {"max Fehler (Zopf additiv, Yang)": max(ybe), "Kontrolle u + iP additiv": kontr_add,
                         "u + iP harmonisch (1/w = 1/u + 1/v)": kontr_harm, "Kontrolle mit ZZ-Term": kontr_zz,
                         "Gatter = R(1)": gate_u1, "dim Kommutant SU(2) auf C^2 x C^2": int(dim_kom),
                         "Dualität R(u) R(1/u) ~ P, R(1)^2 = P": dual}
check(f"[E] Zopfform R(u) = (1 + i u P)/(1 + i u): R_12(u) R_23(u+v) R_12(v) = R_23(v) R_12(u+v) R_23(u), und Yangs Form "
      f"R(u) = u + i P: R_12(u-v) R_13(u) R_23(v) = R_23(v) R_13(u) R_12(u-v) (fünf Zufallspaare, max {max(ybe):.0e}). Das Gatter "
      f"ist R(1) ({gate_u1:.0e}); die in der Aufgabe genannte Form (u + i P)/(u + i) ist dieselbe Familie mit u -> 1/u: als Zopf "
      f"mit additivem Parameter verletzt sie die Gleichung ({kontr_add:.3f}), mit harmonischer Addition 1/w = 1/u + 1/v gilt sie "
      f"({kontr_harm:.0e}). Kontrolle: ein ZZ-Zusatz bricht Yang-Baxter ({kontr_zz:.3f}). Der Kommutant der SU(2) auf C^2 (x) C^2 "
      f"hat Dimension {dim_kom} (span{{1, P}}): Jedes SU(2)-invariante unitäre Zweiplatzgatter ist e^(i chi) R(u). Dualität: "
      f"R(u) R(1/u) ~ P, also R(1/u) ~ P R(u)^(-1); u <-> 1/u tauscht Identität und Tausch, und u = 1 ist ihr Fixpunkt, "
      f"R(1)^2 = P ({dual:.0e})",
      max(ybe) < 1e-12 and kontr_add > 0.1 and kontr_harm < 1e-12 and kontr_zz > 1e-3 and gate_u1 < 1e-12 and dim_kom == 2
      and dual < 1e-12)

say("  (b) Integrabilität")
lc, komm, kontr_t = {}, {}, {}
for L in (4, 6, 8, 10, 12):
    D = 2**L
    xi = np.array([-0.5 if j % 2 == 0 else 0.5 for j in range(L)])
    psi = rng.normal(size=(D, 2)) + 1j*rng.normal(size=(D, 2))
    tm, tp = t_anw(psi, L, -0.5, xi), t_anw(psi, L, 0.5, xi)
    rF, cF = skalar_rest(tm, F_anw(tp, L))
    rT, cT = skalar_rest(t_anw(tp, L, -0.5, xi), T2_anw(psi, L))
    lc[L] = (rF, rT, abs(abs(cF) - 1))
    u_, v_ = rng.normal(size=2) + 1j*rng.normal(size=2)
    tu = lambda V: t_anw(V, L, u_, xi)
    tv = lambda V: t_anw(V, L, v_, xi)
    komm[L] = (nrm(tu(tv(psi)) - tv(tu(psi)))/nrm(tu(tv(psi))), nrm(tu(F_anw(psi, L)) - F_anw(tu(psi), L))/nrm(tu(psi)))
    xi_f = np.array([-0.3 if j % 2 == 0 else 0.3 for j in range(L)])
    kontr_t[L] = nrm(t_anw(F_anw(psi, L), L, u_, xi_f) - F_anw(t_anw(psi, L, u_, xi_f), L))/nrm(t_anw(psi, L, u_, xi_f))
fam = {}
for tau in (0.37, 2.6):
    L = 8
    xi = np.array([-tau/2 if j % 2 == 0 else tau/2 for j in range(L)])
    psi = rng.normal(size=(2**L, 2)) + 1j*rng.normal(size=(2**L, 2))
    u_ = complex(rng.normal(), rng.normal())
    fam[str(tau)] = nrm(t_anw(F_anw(psi, L, u=tau), L, u_, xi) - F_anw(t_anw(psi, L, u_, xi), L, u=tau))/nrm(t_anw(psi, L, u_, xi))
ZAHLEN["transfermatrizen"] = {"ganze Familie R(tau), xi = -+tau/2: |[t(u), F]|": fam,
                              "Lichtkegel: Rest F t(1/2) gegen t(-1/2), Rest t(-1/2) t(1/2) gegen T^2, ||c_F| - 1|":
                                  {str(k): v for k, v in lc.items()},
                              "|[t(u), t(v)]|, |[t(u), F]|": {str(k): v for k, v in komm.items()},
                              "Kontrolle Inhomogenität 0.3: |[t(u), F]|": {str(k): v for k, v in kontr_t.items()}}
check("[E] Inhomogene Transfermatrix mit Yangs R(u) = u + i P und alternierenden Inhomogenitäten xi_j = -+1/2: "
      f"F = t(-1/2) t(1/2)^(-1) bis auf Phase und t(-1/2) t(1/2) = T^2 (Zweiplatztranslation) für L = 4..12 (max "
      f"{max(max(v) for v in lc.values()):.0e}, matrixfrei an Zufallsvektoren); [t(u), t(v)] = 0 und [t(u), F] = 0 für zufällige "
      f"komplexe u, v (max {max(max(v) for v in komm.values()):.0e}). Kontrolle: mit den falschen Inhomogenitäten -+0.3 vertauscht "
      f"t(u) nicht mit F ({min(kontr_t.values()):.2f} bis {max(kontr_t.values()):.2f}). Die Kette ist das Lichtkegelgitter von "
      "Destri-de Vega (rechte Schritte t(-1/2), linke t(1/2)^(-1)) bei Theta = 1/2; für die ganze Familie R(tau) gilt dasselbe mit "
      f"Theta = tau/2 (tau = 0.37, 2.6, L = 8: {max(fam.values()):.0e})",
      max(max(v) for v in lc.values()) < 1e-10 and max(max(v) for v in komm.values()) < 1e-10 and min(kontr_t.values()) > 1e-2
      and max(fam.values()) < 1e-10)

# lokale Ladungen: Anzahl im Raum SU(2)- und T^2-invarianter Operatoren der Reichweite <= 4, L = 8
L = 8
D = 2**L
SWS8 = {x: swap_idx(L, x, (x + 1) % L) for x in range(L)}
F8 = F_anw(np.eye(D, dtype=complex), L, sws=SWS8)


def perm_op(L, sites, perm):
    idx = np.arange(2**L)
    bits = [(idx >> s_) & 1 for s_ in sites]
    neu = idx.copy()
    for s_ in sites:
        neu &= ~(1 << s_)
    for i, s_ in enumerate(sites):
        neu |= bits[i] << sites[perm[i]]
    M = np.zeros((2**L, 2**L))
    M[neu, idx] = 1
    return M


H_nnn8 = sum(perm_op(L, [x, (x + 2) % L], [1, 0]) for x in range(L))
FA8 = F8 @ expm(1.0j*H_nnn8)
BASIS = [sum(perm_op(L, [(x + i) % L for i in range(4)], list(p_)) for x in range(0, L, 2))
         for p_ in itertools.permutations(range(4))]


def zahl_ladungen(Fm):
    A = [Fm @ B_ - B_ @ Fm for B_ in BASIS]
    GA = np.array([[np.vdot(a_, b_) for b_ in A] for a_ in A])
    GB = np.array([[np.vdot(a_, b_) for b_ in BASIS] for a_ in BASIS])
    nA = int(np.sum(np.linalg.eigvalsh(GA) < 1e-9*np.abs(GB).max()))
    nB = int(np.sum(np.linalg.eigvalsh(GB) < 1e-9*np.abs(GB).max()))
    return nA - nB, len(BASIS) - nB


n_int, dimraum = zahl_ladungen(F8)
n_ktrA, _ = zahl_ladungen(FA8)
xi8 = np.array([-0.5 if j % 2 == 0 else 0.5 for j in range(L)])
tp8 = t_anw(np.eye(D), L, -0.5, xi8)
dtp8 = t_anw(np.eye(D), L, -0.5, xi8, ableitung=True)
X8 = np.linalg.solve(tp8, dtp8)
Hb8 = (X8 - X8.conj().T)/2j
Ha8 = (X8 + X8.conj().T)/2
PA = np.array([I2, SX, SY, SZ])
buch = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
rows, cols, ps = buch[:L], buch[L:2*L], buch[2*L:3*L]
spec_ = rows + cols + "," + ",".join(f"{ps[k]}{cols[k]}{rows[k]}" for k in range(L)) + "->" + ps
koeff = np.einsum(spec_, Hb8.reshape([2]*(2*L)), *([PA]*L), optimize=True)/D
gew = {}
for multi in zip(*np.nonzero(np.abs(koeff) > 1e-12)):
    stellen = [L - 1 - k for k in range(L) if multi[k] != 0]                # Achse k <-> Platz L-1-k
    if not stellen:
        r_ = 0
    else:
        frei = sorted(set(range(L)) - set(stellen))
        luecke, lauf = 0, 0
        for x in list(range(L))*2:
            lauf = lauf + 1 if x in frei else 0
            luecke = max(luecke, min(lauf, L))
        r_ = L - luecke
    gew[r_] = float(gew.get(r_, 0.0) + abs(koeff[multi])**2)
ha_konst = nrm(Ha8 - np.trace(Ha8)/D*np.eye(D))
hb_F = nrm(Hb8 @ F8 - F8 @ Hb8)
hb_FA = nrm(Hb8 @ FA8 - FA8 @ Hb8)
ZAHLEN["lokale_ladungen"] = {"Raum SU(2)-, T^2-invariant, Reichweite <= 4 (Dimension)": dimraum,
                             "erhaltene Operatoren integrabel": n_int, "erhaltene Operatoren Kontrolle A": n_ktrA,
                             "Hb: Gewicht je Reichweite": {str(k): v for k, v in sorted(gew.items())},
                             "Ha - Konstante": ha_konst, "||[Hb, F]||": hb_F, "||[Hb, F_A]||": hb_FA}
check(f"[E] Lokale Erhaltungsgrößen (L = 8): Im {dimraum}-dimensionalen Raum der SU(2)- und T^2-invarianten Operatoren der "
      f"Reichweite <= 4 vertauschen {n_int} linear unabhängige mit F (Identität und zwei nichttriviale), mit der Kontrolle A "
      f"(F e^(i H_nnn), H_nnn = sum P_(x,x+2), SU(2)-invariant) nur {n_ktrA} (die Identität). Die Ladung aus der logarithmischen "
      f"Ableitung X = t(-1/2)^(-1) t'(-1/2): hermitescher Teil konstant ({ha_konst:.0e}), antihermitescher Teil Hb lokal mit "
      f"Gewichten je Reichweite { {k: round(v, 6) for k, v in sorted(gew.items())} } (nur Reichweite <= 3), [Hb, F] = {hb_F:.0e}, "
      f"[Hb, F_A] = {hb_FA:.2f}",
      n_int == 3 and n_ktrA == 1 and ha_konst < 1e-10 and max(gew) <= 3 and gew.get(3, 0) > 0.1 and hb_F < 1e-10 and hb_FA > 0.1)

# Niveaustatistik, L = 12
L = 12
D = 2**L
SWS12 = {x: swap_idx(L, x, (x + 1) % L) for x in range(L)}
bas, pos, S2s, T2s = sektor_ops(L, L//2)
nb = len(bas)
E12 = np.zeros((D, nb), complex)
E12[bas, np.arange(nb)] = 1
F12 = F_anw(E12, L, sws=SWS12)[bas]
Hn12 = np.zeros((nb, nb))
for x in range(L):
    Hn12[pos[swap_idx(L, x, (x + 2) % L)[bas]], np.arange(nb)] += 1
FA12 = F12 @ expm(1.0j*Hn12)
MZ = L//2


def bloecke(m):
    Pm = sum(np.exp(-2j*np.pi*m*r_/MZ)*np.linalg.matrix_power(T2s, r_) for r_ in range(MZ))/MZ
    w_, V_ = np.linalg.eigh((Pm + Pm.conj().T)/2)
    Bm = V_[:, w_ > 0.5]
    ws, Vs = np.linalg.eigh(Bm.conj().T @ S2s @ Bm)
    return {S: Bm @ Vs[:, np.abs(ws - S*(S + 1)) < 1e-6] for S in range(7) if np.any(np.abs(ws - S*(S + 1)) < 1e-6)}


def r_werte(ph):
    ph = np.sort(np.mod(ph, 2*np.pi))
    s_ = np.diff(np.concatenate([ph, [ph[0] + 2*np.pi]]))
    return np.minimum(s_, np.roll(s_, -1))/np.maximum(np.maximum(s_, np.roll(s_, -1)), 1e-300)


rv = {"integrabel": [], "Kontrolle A": []}
for m in (1, 2):
    for S, Bb in bloecke(m).items():
        if Bb.shape[1] < 8:
            continue
        for nm, Fm in (("integrabel", F12), ("Kontrolle A", FA12)):
            rv[nm].extend(r_werte(np.angle(np.linalg.eigvals(Bb.conj().T @ Fm @ Bb))))


def gatter2(V, x, y, G, L):
    idx = np.arange(2**L)
    bx, by = (idx >> x) & 1, (idx >> y) & 1
    basis_ = idx & ~((1 << x) | (1 << y))
    tl = bx + 2*by
    out = np.zeros_like(V)
    for l_ in range(4):
        out += G[tl, l_][:, None]*V[basis_ | ((l_ & 1) << x) | (((l_ >> 1) & 1) << y)]
    return out


def zufall_h():
    h = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
    return (h + h.conj().T)/2


UB_o, UB_e = ZIEL @ expm(0.2j*zufall_h()), ZIEL @ expm(0.2j*zufall_h())
idx12 = np.arange(D)
imgf = ((idx12 << 2) | (idx12 >> (L - 2))) & (D - 1)
gesehen, bahnen = np.zeros(D, bool), []
for s_ in range(D):
    if gesehen[s_]:
        continue
    b_ = [s_]
    t_ = imgf[s_]
    while t_ != s_:
        b_.append(t_)
        t_ = imgf[t_]
    gesehen[b_] = True
    bahnen.append(b_)
rB = []
for m in (1, 2):
    spalten = []
    for b_ in bahnen:
        if (m*len(b_)) % MZ:
            continue
        v_ = np.zeros(D, complex)
        for r_, o_ in enumerate(b_):
            v_[o_] = np.exp(2j*np.pi*m*r_/MZ)
        spalten.append(v_/np.linalg.norm(v_))
    Bm = np.array(spalten).T
    V_ = Bm.copy()
    for x in range(0, L, 2):
        V_ = gatter2(V_, x, (x + 1) % L, UB_o, L)
    for x in range(1, L, 2):
        V_ = gatter2(V_, x, (x + 1) % L, UB_e, L)
    rB.extend(r_werte(np.angle(np.linalg.eigvals(Bm.conj().T @ V_))))
rv["Kontrolle B"] = rB
# verborgene Halbschrittsymmetrie Q = U_gerade T bei gleichen Gattern
L = 8
Gh = np.linalg.qr(rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4)))[0]
Ue8 = np.eye(2**L, dtype=complex)
Uo8 = np.eye(2**L, dtype=complex)
for x in range(0, L, 2):
    Uo8 = gatter2(Uo8, x, x + 1, Gh, L)
for x in range(1, L, 2):
    Ue8 = gatter2(Ue8, x, (x + 1) % L, Gh, L)
idx8 = np.arange(2**L)
T1 = np.zeros((2**L, 2**L))
T1[((idx8 << 1) | (idx8 >> (L - 1))) & (2**L - 1), idx8] = 1
Fh = Ue8 @ Uo8
Qh = Ue8 @ T1
halbschritt = (nrm(Qh @ Fh - Fh @ Qh), nrm(Qh @ Qh - Fh @ T1 @ T1))
r_mittel = {k: float(np.mean(v)) for k, v in rv.items()}
ZAHLEN["niveaustatistik"] = {"<r>": r_mittel, "Zahl der r-Werte": {k: len(v) for k, v in rv.items()},
                             "Poisson": 2*np.log(2) - 1, "COE": 0.5307, "CUE": 0.5996,
                             "Halbschritt Q = U_g T: |[Q, F]|, |Q^2 - F T^2|": halbschritt}
check(f"[N] Niveaustatistik L = 12 (Quasienergien je Sektor S_z = 0, Impuls k = 2 pi m/6 mit m = 1, 2, Gesamtspin S): "
      f"<r> integrabel {r_mittel['integrabel']:.4f} ({len(rv['integrabel'])} Werte; Poisson 2 ln 2 - 1 = 0.3863), Kontrolle A "
      f"(F e^(i H_nnn), dieselben Symmetrien) {r_mittel['Kontrolle A']:.4f} (COE 0.5307), Kontrolle B (Gatter mit kleinem "
      f"zufälligem Zweiplatzterm, je Lage verschieden, nur Impuls erhalten) {r_mittel['Kontrolle B']:.4f} ({len(rB)} Werte; CUE "
      f"0.5996). Methodisch: Bei gleichen Gattern in beiden Lagen vertauscht der Halbschritt Q = U_gerade T mit F und Q^2 = F T^2 "
      f"(Haar-Gatter, {halbschritt[0]:.0e}, {halbschritt[1]:.0e}); er überlagert zwei Spektren, deshalb trägt Kontrolle B "
      "verschiedene Terme je Lage",
      abs(r_mittel["integrabel"] - 0.3863) < 0.04 and r_mittel["Kontrolle A"] > 0.49 and r_mittel["Kontrolle B"] > 0.55
      and max(halbschritt) < 1e-10)
befund("E/N: Die Kette ist integrabel: Lichtkegelfaktorisierung, vertauschende Transfermatrizen, lokale Ladungen und "
       "Poisson-Statistik; die SU(2)-symmetrische Kontrolle mit übernächster Nachbarkopplung ist chaotisch (COE). Weil jedes "
       "SU(2)-invariante Zweiplatzgatter R(u) ist, ist die Integrabilität bei gleichen Gattern in beiden Lagen automatisch", True)

say("  (c) Einmagnonen-Dispersion")
magnon, ok_mag = {}, True
for L in (4, 6, 8, 10, 12):
    D = 2**L
    b1 = sektor(L, 1)
    E1 = np.zeros((D, L), complex)
    E1[b1, np.arange(L)] = 1
    Fm = F_anw(E1, L)[b1]
    lam = np.linalg.eigvals(Fm)
    ks = 2*np.pi*np.arange(L//2)/(L//2)
    eps_ = np.arccos((1 - np.cos(ks))/2)
    erwartet = np.concatenate([-1j*np.exp(1j*eps_), -1j*np.exp(-1j*eps_)])
    d_ = max(min(abs(l_ - e_) for e_ in erwartet) for l_ in lam) + max(min(abs(l_ - e_) for l_ in lam) for e_ in erwartet)
    fm = F_anw(np.eye(D)[:, :1].astype(complex), L)[:, 0]
    magnon[str(L)] = {"Abstand zum Spektrum": float(d_), "F|FM> = |FM>": nrm(fm - np.eye(D)[:, 0])}
    ok_mag &= d_ < 1e-10 and nrm(fm - np.eye(D)[:, 0]) < 1e-12
kk = np.linspace(-np.pi, np.pi, 201)
bloch = 0.0
for k in kk:
    Vo = np.exp(-1j*np.pi/4)*np.array([[1, 1j], [1j, 1]])/np.sqrt(2)
    Ve = np.exp(-1j*np.pi/4)*np.array([[1, 1j*np.exp(-1j*k)], [1j*np.exp(1j*k), 1]])/np.sqrt(2)
    th_ = np.pi/4
    Ub = np.array([[np.cos(th_), -1j*np.sin(th_)], [-1j*np.sin(th_), np.cos(th_)]])
    Ub2 = np.array([[np.cos(th_), 1j*np.sin(th_)], [1j*np.sin(th_), np.cos(th_)]])
    Tk = np.diag([1, np.exp(1j*(k + np.pi))])
    z120 = (Tk.conj() @ Ub2 @ Tk) @ Ub
    a_ = np.sort_complex(np.linalg.eigvals(Ve @ Vo))
    b_ = np.sort_complex(-1j*np.linalg.eigvals(z120))
    bloch = max(bloch, min(np.abs(a_ - b_).max(), np.abs(a_ - b_[::-1]).max()))
km = 1e-3
lin = np.arccos((1 - np.cos(np.pi + km))/2)/km
gold = (np.pi/2 - np.arccos((1 - np.cos(km))/2))/km**2
ZAHLEN["magnon"] = {"L": magnon, "Bloch gegen z120 mit k -> k + pi und Faktor -i": bloch,
                    "eps/|k - pi| bei k -> pi": lin, "Goldstone: (pi/2 - eps)/k^2 bei k -> 0": gold}
check(f"[E] Ein umgeklappter Spin auf dem ferromagnetischen Zustand (F|FM> = |FM>): Die Eigenwerte von F sind genau "
      f"-i e^(+-i eps(k)) mit cos eps(k) = (1 - cos k)/2, k = 2 pi m/(L/2) (L = 4..12, max Abstand "
      f"{max(v['Abstand zum Spektrum'] for v in magnon.values()):.0e}). Das ist z120s cos eps = (1 + cos k)/2 nach k -> k + pi "
      f"(Orientierung der Hüpfphasen, Eichung) mal -i = (e^(-i pi/4))^2 (Bloch-Symbole auf 201 Impulsen gleich, {bloch:.0e}). "
      f"Bei k = pi berühren sich die Bänder linear, eps/|k - pi| = {lin:.6f} = 1/sqrt 2; bei k = 0 liegt die Goldstone-Mode "
      f"mit Eigenwert 1 und quadratischer Dispersion (pi/2 - eps)/k^2 = {gold:.4f} = 1/4",
      ok_mag and bloch < 1e-12 and abs(lin - 1/np.sqrt(2)) < 1e-4 and abs(gold - 0.25) < 1e-3)

say("  (d) Kontinuum: Destri-de-Vega-Lesart")


def S2_anw(V, L):
    out = (3*L/4)*V.copy()
    for i in range(L):
        for j in range(i + 1, L):
            out += V[swap_idx(L, i, j)] - 0.5*V
    return out


def verschr(psi, L):
    s_ = np.linalg.svd(psi.reshape(2**(L - L//2), 2**(L//2)), compute_uv=False)**2
    s_ = s_[s_ > 1e-15]
    return float(-(s_*np.log(s_)).sum())


def sektorbasis(L, nup, S):
    b_, _, S2_, T2_ = sektor_ops(L, nup)
    M_ = S2_ - S*(S + 1)*np.eye(len(b_))
    M_ = M_.conj().T @ M_ + 2*np.eye(len(b_)) - T2_ - T2_.conj().T
    w_, V_ = np.linalg.eigh(M_)
    return b_, V_[:, w_ < 1e-9]


def hb_sektor(L, nup, S, tau, BB=None):
    th = tau/2
    xi = np.array([-th if j % 2 == 0 else th for j in range(L)])
    b_, Bs = BB if BB is not None else sektorbasis(L, nup, S)
    EB = np.zeros((2**L, Bs.shape[1]), complex)
    EB[b_] = Bs
    X = np.linalg.solve(Bs.conj().T @ t_anw(EB, L, -th, xi)[b_], Bs.conj().T @ t_anw(EB, L, -th, xi, ableitung=True)[b_])
    w_, V_ = np.linalg.eigh((X - X.conj().T)/2j)
    return b_, Bs, w_, V_


def vakuum_daten(L, tau, basen):
    out = {}
    for nup, S in ((L//2, 0), (L//2 + 1, 1)):
        b_, Bs, w_, V_ = hb_sektor(L, nup, S, tau, basen[(nup, S)])
        psi = np.zeros(2**L, complex)
        psi[b_] = Bs @ V_[:, -1]
        Fpsi = F_anw(psi[:, None], L, u=tau)[:, 0]
        lam = np.vdot(psi, Fpsi)
        out[S] = (psi, lam, nrm(Fpsi - lam*psi), w_, V_, b_, Bs)
    return out


v_pred = lambda tau: 2*np.tanh(np.pi*tau/2)
dvd, basen_alle = {}, {}
for L in (6, 8, 10, 12):
    basen_alle[L] = {(nup, S): sektorbasis(L, nup, S) for nup, S in ((L//2, 0), (L//2 + 1, 1))}
    for tau in (0.5, 1.0):
        vd = vakuum_daten(L, tau, basen_alle[L])
        deps = -np.angle(vd[1][1]/vd[0][1])
        dvd[(L, tau)] = {"x_eff = -L deps/(2 pi v)": float(-L*deps/(2*np.pi*v_pred(tau))), "S(L/2)": verschr(vd[0][0], L),
                         "Residuen": max(vd[0][2], vd[1][2])}
        if tau == 1.0:
            # typische Eigenzustände desselben Sektors (S = 0, k = 0) als Gegenprobe zur Verschränkung
            psi0, lam0, r0, w0, V0, b0, B0 = vd[0]
            typ = []
            for k_ in range(V0.shape[1]):
                ps_ = np.zeros(2**L, complex)
                ps_[b0] = B0 @ V0[:, k_]
                typ.append(verschr(ps_, L))
            dvd[(L, tau)]["S typisch (Mittel im Sektor)"] = float(np.mean(typ))
# Heisenberg-Vergleich (H = sum S.S)
heis = {}
for L in (6, 8, 10, 12):
    D = 2**L
    swl = {x: swap_idx(L, x, (x + 1) % L) for x in range(L)}
    e_ = {}
    for nup in (L//2, L//2 + 1):
        b_ = sektor(L, nup)
        E_ = np.zeros((D, len(b_)))
        E_[b_, np.arange(len(b_))] = 1
        H_ = sum(0.5*E_[swl[x]] - 0.25*E_ for x in range(L))[b_]
        w_, V_ = np.linalg.eigh(H_)
        e_[nup] = (w_[0], V_[:, 0], b_)
    psi = np.zeros(D)
    psi[e_[L//2][2]] = e_[L//2][1]
    heis[L] = {"x_eff": float(L*(e_[L//2 + 1][0] - e_[L//2][0])/(2*np.pi*np.pi/2)), "S(L/2)": verschr(psi, L)}
# Nachverfolgung L = 8: vom Heisenberg-Grundzustand (tau = 0.02) stetig bis tau = 1
L = 8
prev, verfolgt = None, True
H8 = sum(0.5*np.eye(2**L)[swap_idx(L, x, (x + 1) % L)] - 0.25*np.eye(2**L) for x in range(L))
b8, B8 = basen_alle[8][(4, 0)]
for it, tau in enumerate(np.linspace(0.02, 1.0, 50)):
    _, _, w_, V_ = hb_sektor(L, 4, 0, tau, (b8, B8))
    if it == 0:
        hs = B8.conj().T @ H8[np.ix_(b8, b8)] @ B8
        g0 = np.linalg.eigh(hs)[1][:, 0]
        k_ = int(np.argmax(np.abs(V_.conj().T @ g0)**2))
        start_ueberlapp = float(np.abs(np.vdot(V_[:, k_], g0))**2)
    else:
        ov = np.abs(V_.conj().T @ prev)**2
        k_ = int(np.argmax(ov))
        verfolgt &= ov[k_] > 0.9
    prev = V_[:, k_]
verfolgt &= k_ == len(w_) - 1
Ls = np.array([6, 8, 10, 12])
Sv = np.array([dvd[(L, 1.0)]["S(L/2)"] for L in Ls])
Sh = np.array([heis[L]["S(L/2)"] for L in Ls])
St = np.array([dvd[(L, 1.0)]["S typisch (Mittel im Sektor)"] for L in Ls])
c_fit = 3*np.polyfit(np.log(Ls), Sv, 1)[0]
c_heis = 3*np.polyfit(np.log(Ls), Sh, 1)[0]
c_typ = 3*np.polyfit(np.log(Ls), St, 1)[0]
xe = {f"L={L}, tau={tau}": round(v["x_eff = -L deps/(2 pi v)"], 4) for (L, tau), v in dvd.items()}
xh = {f"L={L}": round(v["x_eff"], 4) for L, v in heis.items()}
masse = 4*np.exp(-np.pi/2)
ZAHLEN["kontinuum"] = {"Theta (Lichtkegel)": 0.5, "Selbstdualität u <-> 1/u": "u = 1 Fixpunkt",
                       "v(tau) = 2 tanh(pi tau/2), tau = 1": v_pred(1.0), "x_eff Kette": xe, "x_eff Heisenberg": xh,
                       "S(L/2) Vakuum": dict(zip(map(str, Ls), Sv.round(5).tolist())),
                       "S(L/2) Heisenberg": dict(zip(map(str, Ls), Sh.round(5).tolist())),
                       "S(L/2) typisch": dict(zip(map(str, Ls), St.round(5).tolist())),
                       "c_fit Vakuum, Heisenberg, typisch": [c_fit, c_heis, c_typ],
                       "Nachverfolgung L = 8 aus Heisenberg (stetig, endet im Maximum von Hb)": bool(verfolgt),
                       "Startüberlapp tau = 0.02": start_ueberlapp,
                       "DdV-Masse m a = 4 e^(-pi Theta) bei Theta = 1/2 (Lesart)": masse,
                       "Residuen der Floquet-Eigenzustände": max(v["Residuen"] for v in dvd.values())}
check(f"[N] Lichtkegelgitter bei Theta = 1/2 (u = 1, Fixpunkt der Dualität u <-> 1/u zwischen Identität und Tausch). Vakuum = "
      f"oberster Zustand der lokalen Ladung Hb im Sektor S = 0, k = 0, stetig aus dem Heisenberg-Grundzustand verfolgt (L = 8, "
      f"50 Schritte tau = 0.02..1, Startüberlapp {start_ueberlapp:.4f}, {verfolgt}); Anregung = oberster Zustand im Sektor S = 1, "
      f"k = 0; beide Floquet-Eigenzustände (Residuum {ZAHLEN['kontinuum']['Residuen der Floquet-Eigenzustände']:.0e}). Die "
      f"Floquet-Lücke skaliert wie 1/L: Delta eps = -2 pi v x_eff/L mit der Spinonengeschwindigkeit v = 2 tanh(pi tau/2) "
      f"(Bethe-Ansatz-Lesart des Lichtkegelgitters [L]; {v_pred(1.0):.4f} Plätze je Periode bei tau = 1, unter der "
      f"Lichtkegelgeschwindigkeit 2; Trotter-Grenze pi tau), x_eff = {xe} bei tau = 0.5 und 1 gegen Heisenberg {xh}: dieselbe "
      f"Zahl, x = 1/2 mit logarithmischer Korrektur, also bestätigt die Rechnung v(tau) [N]. Verschränkung S(L/2) des Vakuums {Sv.round(4).tolist()}, Anpassung (c/3) ln L: c = {c_fit:.3f} (Heisenberg "
      f"{c_heis:.3f}); typische Eigenzustände desselben Sektors {St.round(3).tolist()}, c_fit = {c_typ:.1f} (Volumengesetz). Eine "
      f"Masse m a = 4 e^(-pi/2) = {masse:.3f} (DdV-Formel, Korrelationslänge ~ 1.2 Plätze) würde L Delta eps linear wachsen lassen; "
      "beobachtet ist L Delta eps konstant",
      verfolgt and all(0.38 < v_ < 0.5 for v_ in xe.values()) and 0.75 < c_fit < 1.4 and 0.75 < c_heis < 1.2 and c_typ > 3
      and ZAHLEN["kontinuum"]["Residuen der Floquet-Eigenzustände"] < 1e-9)
befund("N/L: Bei u = 1 ist der Vakuumsektor in endlicher Größe kritisch und su(2)_1-artig: c nahe 1, x nahe 1/2 mit "
       "logarithmischer Korrektur wie bei Heisenberg, lineare Geschwindigkeit 2 tanh(pi/2). Bei festem Theta gibt es keinen "
       "Parameter, der eine Masse in Gittereinheiten gegen null schickt; der Skalenlimes des Vakuumsektors ist daher der "
       "IR-Fixpunkt, nach den Daten das su(2)_1-WZW-Modell mit marginaler Strom-Strom-Kopplung. Das chirale "
       "Gross-Neveu-Modell ist der Destri-de-Vega-Doppellimes Theta -> oo (Gatter -> Tausch), a -> 0 bei fester Masse "
       "m = (4/a) e^(-pi Theta) [L, Literatur]; er wird bei u = 1 nicht erreicht. Offen: das Vorzeichen der marginalen Kopplung "
       "(L <= 12 entscheidet es nicht) und eine Rechnung jenseits L = 12", True)

# =================================================================== Teil 3
say("\nTeil 3: Drei Nähte")
QA6, QB6, QC6 = (n6[2*i] + n6[2*i + 1] for i in range(3))


def nahtgatter(i, j, theta=np.pi/2):
    Ep = sum(c6[2*i + s_].T @ c6[2*j + s_] for s_ in range(2))
    Bij = expm(-1j*theta*(Ep + Ep.T)/2)
    Qi, Qj = n6[2*i] + n6[2*i + 1], n6[2*j] + n6[2*j + 1]
    return Bij.conj().T @ expm(-1j*np.pi*(Qi @ Qi + Qj @ Qj)/4) @ Bij


U_AB, U_BC = nahtgatter(0, 1), nahtgatter(1, 2)
W = U_BC @ U_AB
W_spin = np.kron(I2, R_zopf(1.0)) @ np.kron(R_zopf(1.0), I2)
spin_gleich = bis_auf_phase(W[np.ix_(PHYS3, PHYS3)], W_spin)
H_eff = 1j*logm(W_spin)
phasen_W = np.sort(np.angle(np.linalg.eigvals(W_spin)))
al_ = np.arccos(3/4)
phasen_soll = np.sort([-np.pi/2 - al_]*2 + [-np.pi/2 + al_]*2 + [0.0]*4)
phasen_ok = float(np.abs(phasen_W - phasen_soll).max())
PAULI3 = {}
for p_ in itertools.product(range(4), repeat=3):
    op = np.kron(np.kron(PA[p_[0]], PA[p_[1]]), PA[p_[2]])
    PAULI3[p_] = np.trace(op @ H_eff)/8
H3 = sum(c_*np.kron(np.kron(PA[p_[0]], PA[p_[1]]), PA[p_[2]]) for p_, c_ in PAULI3.items() if sum(x != 0 for x in p_) == 3)
eps3 = np.zeros((3, 3, 3))
for a_, b_, c_ in itertools.permutations(range(3)):
    eps3[a_, b_, c_] = np.linalg.det(np.eye(3)[[a_, b_, c_]])
CHI = sum(eps3[a_, b_, c_]*np.kron(np.kron(PA[a_ + 1], PA[b_ + 1]), PA[c_ + 1])
          for a_, b_, c_ in itertools.product(range(3), repeat=3))/8
c_chi = (np.trace(CHI @ H3)/np.trace(CHI @ CHI)).real
rest3 = nrm(H3 - c_chi*CHI)
c_formel = 4*np.arccos(3/4)/np.sqrt(7)
SdotS = lambda i, j: sum(np.kron(np.kron(*(PA[a_] if k == i or k == j else I2 for k in range(2))),
                                 PA[a_] if 2 in (i, j) else I2)/4 for a_ in (1, 2, 3))
j_paar = {nm: (np.trace(SdotS(i, j) @ H_eff)/np.trace(SdotS(i, j) @ SdotS(i, j))).real
          for nm, (i, j) in (("AB", (0, 1)), ("BC", (1, 2)), ("AC", (0, 2)))}
W_rev = np.kron(R_zopf(1.0), I2) @ np.kron(I2, R_zopf(1.0))
H_rev = 1j*logm(W_rev)
c_rev = (np.trace(CHI @ H_rev)/np.trace(CHI @ CHI)).real
c_von_u = {}
for uu in (0.25, 0.5, 1.0, 1.5):
    Hu = 1j*logm(np.kron(I2, R_zopf(uu)) @ np.kron(R_zopf(uu), I2))
    c_von_u[str(uu)] = float((np.trace(CHI @ Hu)/np.trace(CHI @ CHI)).real)
H_frei6 = sum(n_ for n_ in n6)/2


def freies_gatter(i, j):
    Ep = sum(c6[2*i + s_].T @ c6[2*j + s_] for s_ in range(2))
    Bij = expm(-1j*np.pi/2*(Ep + Ep.T)/2)
    return Bij.conj().T @ expm(-1j*np.pi*H_frei6/2) @ Bij


W_frei = freies_gatter(1, 2) @ freies_gatter(0, 1)
frei_phase = bis_auf_phase(W_frei[np.ix_(PHYS3, PHYS3)], np.eye(8))
TR = np.kron(np.kron(SY, SY), SY)
t_ungerade = nrm(TR @ CHI.conj() @ TR + CHI)
perm_ungerade = nrm(np.kron(PSWAP, I2) @ CHI @ np.kron(PSWAP, I2) + CHI)
ZAHLEN["drei_naehte_verbunden"] = {"W(Fock) = R_BC R_AB (Abstand)": spin_gleich, "c (Chiralität)": c_chi,
                                   "4 arccos(3/4)/sqrt 7": c_formel, "Rest des Dreikörperteils": rest3,
                                   "Paarkopplungen j_AB, j_BC, j_AC": j_paar, "c für U_AB U_BC": c_rev,
                                   "c(u)": c_von_u, "freie Energie: W reine Phase": frei_phase,
                                   "chi T-ungerade, ungerade unter Tausch": [t_ungerade, perm_ungerade]}
check(f"[E] Drei Nähte A, B, C (sechs Moden, Fock 64): W = U_BC U_AB ist auf (1,1,1) genau R_BC(1) R_AB(1) (Abstand "
      f"{spin_gleich:.0e}). Effektiver Erzeuger H = i log W, Pauli-Zerlegung (Möbius auf {{A, B, C}}): Der Dreikörperteil ist genau "
      f"c chi mit der Chiralität chi = s_A.(s_B x s_C) (Rest {rest3:.0e}), c = {c_chi:.6f} = 4 arccos(3/4)/sqrt 7 "
      f"({c_formel:.6f}; Eigenphasen von W: 0 auf S = 3/2, -pi/2 +- arccos(3/4) auf S = 1/2, {phasen_ok:.0e}): fest, ohne "
      f"freien Parameter; die "
      f"umgekehrte Reihenfolge U_AB U_BC gibt {c_rev:.6f}. Dazu eine induzierte "
      f"Paarkopplung j_AC = {j_paar['AC']:.6f} zwischen den nie direkt verbundenen Nähten (j_AB = {j_paar['AB']:.6f}, j_BC = "
      f"{j_paar['BC']:.6f}). Längs der Familie hängt c an u ({ {k: round(v, 5) for k, v in c_von_u.items()} }), u = 1 legt ihn "
      f"fest. chi ist zeitungerade und ungerade unter jedem Nahttausch ({t_ungerade:.0e}, {perm_ungerade:.0e}): ein "
      f"epsilon-(f-)Kanal, kein symmetrischer phi^3-Kanal. Kontrolle: mit der freien Energie ist W eine reine Phase ({frei_phase:.0e})",
      spin_gleich < 1e-12 and rest3 < 1e-12 and abs(c_chi - c_formel) < 1e-10 and abs(c_rev + c_formel) < 1e-10
      and phasen_ok < 1e-12
      and abs(j_paar["AC"]) > 1e-3 and frei_phase < 1e-12 and t_ungerade < 1e-12 and perm_ungerade < 1e-12
      and len(set(round(v, 6) for v in c_von_u.values())) == 4)

Nup6, Ndn6 = n6[0] + n6[2] + n6[4], n6[1] + n6[3] + n6[5]
D6 = n6[0] @ n6[1] + n6[2] @ n6[3] + n6[4] @ n6[5]
HQ6 = (QA6 @ QA6 + QB6 @ QB6 + QC6 @ QC6)/2
regeln = {nm: nrm(W @ X_ - X_ @ W) for nm, X_ in (("N_hoch", Nup6), ("N_runter", Ndn6), ("D", D6), ("H_Q", HQ6))}
W_th = nahtgatter(1, 2, np.pi/3) @ nahtgatter(0, 1, np.pi/3)
regel_kontr = nrm(W_th @ D6 - D6 @ W_th)
dd = np.diag(D6).round().astype(int)
amp_12 = max(abs(W[s_, t_]) for s_ in range(64) for t_ in range(64) if dd[s_] != dd[t_])
amp_12_kontr = max(abs(W_th[s_, t_]) for s_ in range(64) for t_ in range(64) if dd[s_] != dd[t_])
s20 = [s_ for s_ in range(16) if np.isclose(QA4[s_, s_], 2)]
amp_20_11 = max(abs(U1[t_, s_]) for s_ in s20 for t_ in PHYS)
zopf_def = (np.kron(R_zopf(1.0), I2) @ np.kron(I2, R_zopf(1.0)) @ np.kron(R_zopf(1.0), I2)
            - np.kron(I2, R_zopf(1.0)) @ np.kron(R_zopf(1.0), I2) @ np.kron(I2, R_zopf(1.0)))
p_, q_ = 1/(1 + 1j), 1j/(1 + 1j)
zopf_formel = nrm(zopf_def - p_*p_*q_*(P12 - P23))
drei_in_def = sum(abs(np.trace(np.kron(np.kron(PA[a_], PA[b_]), PA[c_]) @ zopf_def)/8)**2
                  for a_, b_, c_ in itertools.product(range(1, 4), repeat=3))
ZAHLEN["auswahlregeln"] = {"||[W, X]||": regeln, "Kontrolle theta = pi/3: ||[W, D]||": regel_kontr,
                           "max Amplitude zwischen verschiedenen D": amp_12, "dasselbe bei theta = pi/3": amp_12_kontr,
                           "zwei Nähte: (2,0) -> (1,1)": amp_20_11,
                           "Zopfdefekt R R R - R R R = p^2 q (P_AB - P_BC)": zopf_formel, "Dreikörpergewicht im Zopfdefekt": drei_in_def}
check(f"[E] Auswahlregeln im ganzen Fockraum dreier Nähte: W vertauscht mit N_hoch, N_runter, der Doppelbesetzungszahl D und "
      f"H_Q ({ {k: f'{v:.0e}' for k, v in regeln.items()} }); also keine Amplitude zwischen verschiedenen D (max {amp_12:.0e}), "
      f"insbesondere keine 1->2-Amplitude 'doppelt besetzte Naht -> zwei einfach besetzte' (zwei Nähte: (2,0) -> (1,1) "
      f"{amp_20_11:.0e}; Grund: J_x^2, J_y^2, J_z^2 vertauschen im Isospin 1). Kontrolle: bei theta = pi/3 ist [W, D] = "
      f"{regel_kontr:.2f} und die Amplitude {amp_12_kontr:.3f}. Faktorisierung: Die Zopfgleichung bei gleichen Parametern ist "
      f"um R R R - R R R = p^2 q (P_AB - P_BC) verletzt (p = 1/(1+i), q = i/(1+i), {zopf_formel:.0e}), ein reiner "
      f"Zweikörperdefekt (Dreikörpergewicht {drei_in_def:.0e}); mit additiven Rapiditäten verschwindet er (Teil 2a): "
      "Dreiteilchenstreuung = Produkt von Zweiteilchenstreuungen",
      all(v < 1e-12 for v in regeln.values()) and amp_12 < 1e-12 and amp_20_11 < 1e-12 and regel_kontr > 0.1
      and amp_12_kontr > 0.1 and zopf_formel < 1e-12 and drei_in_def < 1e-20)

# radiale Struktur: zwei native Nähte (ein Fermion je Naht), Moden x = -7/4, -3/4, 1/4, 5/4, Gatter aus Nullmoden
c8r = jw(8)
n8r = [x.T @ x for x in c8r]
vakr = np.array([1 if x < 0 else 0 for x in XS]*2)
Qr = [sum(n8r[4*i + k] - vakr[4*i + k]*np.eye(256) for k in range(4)) for i in range(2)]
Enat_r = sum((XS[k % 4] + 0.25)*(n8r[k] - vakr[k]*np.eye(256)) for k in range(8))
Lp = sum((XS[k] + 0.75)*c8r[4*i + k + 1].T @ c8r[4*i + k] for i in range(2) for k in range(3))   # sum (k+1) psi*_(k+1) psi_k
J0 = sum(c8r[k].T @ c8r[4 + k] for k in range(4))
J1 = sum(c8r[k + 1].T @ c8r[4 + k] for k in range(3))
Ur = (lambda Bm: Bm.conj().T @ expm(-1j*np.pi*(Qr[0] @ Qr[0] + Qr[1] @ Qr[1])/4) @ Bm)(expm(-1j*np.pi/4*(J0 + J0.T)))
Ur_k = (lambda Bm: Bm.conj().T @ expm(-1j*np.pi*(Qr[0] @ Qr[0] + Qr[1] @ Qr[1])/4) @ Bm)(expm(-1j*np.pi/4*(J1 + J1.T)))
radial = {"||[U, L^nat_1]||": nrm(Ur @ Lp - Lp @ Ur), "||[U, L^nat_-1]||": nrm(Ur @ Lp.T - Lp.T @ Ur),
          "||[U, E^nat]||": nrm(Ur @ Enat_r - Enat_r @ Ur), "Kontrolle J_1: ||[U_1, E^nat]||": nrm(Ur_k @ Enat_r - Enat_r @ Ur_k),
          "U nichttrivial (Abstand zur Phase)": bis_auf_phase(Ur, np.eye(256))}
w_l, V_l = np.linalg.eigh(Lp.T @ Lp)
prim = V_l[:, w_l < 1e-10]
radial["U erhält ker L^nat_1 (Primärzustände)"] = nrm(Lp @ Ur @ prim)
ZAHLEN["radial"] = radial
gamma_gatter = 0.0**2/32
check(f"[E] Radiale Struktur (z99/z83): zwei native Nähte (je ein Dirac-Fermion, Moden -7/4, -3/4, 1/4, 5/4, Fock 256), Gatter "
      f"aus den Nullmoden (Ladungsenergie Q_A^2/2 + Q_B^2/2 und Nullmode J_0 des Transferstroms, Vierteldrehung, Viertelschritt): "
      f"U vertauscht mit dem nativen Heber L^nat_(+-1) = sum (k+1) psi*_(k+1) psi_k (beide Nähte) und mit E^nat "
      f"({radial['||[U, L^nat_1]||']:.0e}, {radial['||[U, L^nat_-1]||']:.0e}, {radial['||[U, E^nat]||']:.0e}), ist aber nicht "
      f"trivial ({radial['U nichttrivial (Abstand zur Phase)']:.2f}); es bildet Primärzustände (ker L^nat_1) auf Primärzustände ab "
      f"({radial['U erhält ker L^nat_1 (Primärzustände)']:.0e}). Kontrolle: mit der verschobenen Mode J_1 vertauscht das Gatter "
      f"nicht mit E^nat ({radial['Kontrolle J_1: ||[U_1, E^nat]||']:.2f}). Also erzeugt die Komposition keine Box-Nachkommen, "
      f"Delta = 1 bleibt starr (z99, z105), und der phi^3-Term fehlt: lambda = 0, gamma = lambda^2/32 = {gamma_gatter:g}",
      radial["||[U, L^nat_1]||"] < 1e-12 and radial["||[U, L^nat_-1]||"] < 1e-12 and radial["||[U, E^nat]||"] < 1e-12
      and radial["Kontrolle J_1: ||[U_1, E^nat]||"] > 0.1 and radial["U nichttrivial (Abstand zur Phase)"] > 0.1
      and radial["U erhält ker L^nat_1 (Primärzustände)"] < 1e-12)
befund("E/L: Drei Nähte haben einen verbundenen Anteil mit festem Koeffizienten, c = +-4 arccos(3/4)/sqrt 7 vor der Chiralität "
       "(Vorzeichen = Reihenfolge der Kompositionen; dazu eine induzierte Kopplung A-C). Er ist zahlerhaltend, doublonenerhaltend, zeitungerade und antisymmetrisch: ein "
       "epsilon-Kanal wie die f-Struktur von z75/z99, kein phi^3. Eine verbundene 1->2-Amplitude gibt es nicht, weder in "
       "den Besetzungen noch radial (das Gatter vertauscht mit L^nat_(+-1)); Integrabilität verbietet Teilchenerzeugung ohnehin. "
       "Gegen z108: lambda = 0, gamma = 0; der Rekombinationsweg Box phi = lambda phi^3 bleibt offen", True)

# =================================================================== Teil 4
say("\nTeil 4: Bilanz")
BILANZ = [
    ("E", "Lesarten der zwei Moden je Naht (Quarterpaar, zwei Stufen, zwei Spezies, eta-Dublett): mit nativem E^nat ist das "
          "Gatter stets ein Produkt von Einzelnahtphasen (E^nat einkörperig, vertauscht mit dem Transfer)"),
    ("E", "Ladungsgitter: U(1)_1 (kein Dublett) und U(2)_1 = u(1) + su(2)_1 (Dublett 1/2, Singulett 1, U_H = 0) sind nativ; "
          "spinblindes U(1)_1 mit Dublett = entartete Hülle g = 1 (z99), U_H = g"),
    ("C", "z120-Gatter = Nullmodenteil der Lesart 'zwei Stufen eines Fermions' (Prämisse P_null), gleichwertig g = 1"),
    ("E", "Vierteldrehung = Verträglichkeit [B^+ D B, D] = 0 (nur theta = 0, pi/2, pi); die Phase phi ist frei: R(u), "
          "u = tan(phi/2); [U, H_Q] = 0, KMS-unsichtbar"),
    ("E unter P_null", "Viertelschritt der Rekombinationsenergie = Wurzel der Klein-Kette; gedrehte Kette = kanonischer Tausch; "
                       "U^2 = -Gamma(pi); Niveau k gibt u = tan(pi/(4k)); Wurzelvorzeichen = Zeitrichtung"),
    ("E", "P_ident und Kette schließen sich aus (L = 2 trivial, L >= 3 leer): die Kette braucht P_ord"),
    ("E", "Yang-Baxter, Lichtkegelfaktorisierung F = t(-1/2) t(1/2)^(-1), t t = T^2, vertauschende Transfermatrizen L = 4..12, "
          "zwei lokale Ladungen (Reichweite <= 3) gegen keine in der Kontrolle"),
    ("N", "Niveaustatistik L = 12: Poisson 0,386 gegen COE 0,535 (SU(2)-Kontrolle) und CUE 0,60 (Zufallskontrolle)"),
    ("E", "Einmagnonen-Dispersion -i e^(+-i eps), cos eps = (1 - cos k)/2 = z120 nach k -> k + pi"),
    ("N/L", "Kontinuum bei u = 1 (Theta = 1/2, selbstdual): Vakuumsektor kritisch, c ~ 1, x ~ 1/2 (log), v = 2 tanh(pi/2); "
            "Lesart su(2)_1-WZW; chirales Gross-Neveu nur im DdV-Doppellimes Theta -> oo"),
    ("E", "Drei Nähte: verbundener Teil c chi, c = +-4 arccos(3/4)/sqrt 7 (Vorzeichen = Reihenfolge), induziertes j_AC; "
          "epsilon-Kanal, zeitungerade"),
    ("E", "keine 1->2-Amplitude (N_s, D, H_Q erhalten; radial: [U, L^nat_(+-1)] = 0); lambda = 0, gamma = 0"),
    ("O", "eine native Begründung von P_null (Komposition sieht nur das Ladungsgitter) und von P_viertel (u = 1)"),
    ("O", "phi^3-Rekombination über drei Nähte: weiterhin ohne nativen Operator der Dimension 3"),
]
for st, txt in BILANZ:
    say(f"   [{st}] {txt}")
ZAHLEN["bilanz"] = [{"status": st, "aussage": txt} for st, txt in BILANZ]

urteil = ("Die Regel 'Komposition = Viertelschritt der Rekombinationsenergie im viertelgedrehten Transferrahmen' ist keine native "
          "Vorwärtsregel, sondern eine bedingte. Nativ ist E^nat einkörperig: In jeder Lesart der zwei Moden je Naht "
          "(Quarterpaar, zwei Stufen eines Fermions, zwei Spezies gleicher Stufe, eta-Dublett) vertauscht E^nat mit dem "
          "Transfer, und das Gatter ist ein Produkt von Einzelnahtphasen. Die nativen Gitter sind U(1)_1 ohne Dublett und "
          "U(2)_1 = u(1) + su(2)_1 mit U_H = 0; das spinblinde U(1)_1 mit Dublett ist die entartete Hülle g = 1 von z99. Das "
          "z120-Gatter ist genau der Nullmodenteil der Lesart 'zwei Stufen eines Fermions' (P_null). Unter P_null gilt die "
          "Regel exakt: i^(Q_A Q_B) ist eine Wurzel der Klein-Kette, die viertelgedrehte Kette ist der kanonische Tausch, und "
          "U ist seine Wurzel (U^2 = -Gamma(pi)). Erzwungen ist nur die Vierteldrehung (Verträglichkeit [B^+ D B, D] = 0); die "
          "Phase ist der Spektralparameter u = tan(phi/2), u = 1 kommt aus P_viertel, das Vorzeichen der Wurzel ist die "
          "Zeitrichtung. P_ident und Kette schließen sich aus; die Kette braucht unterscheidbare Nähte (P_ord). Die Kette ist "
          "die integrable XXX-Trotterisierung im selbstdualen Punkt u = 1, exakt das Destri-de-Vega-Lichtkegelgitter bei "
          "Theta = 1/2 (F = t(-1/2) t(1/2)^(-1), zwei lokale Ladungen, Poisson gegen COE/CUE, Einmagnonen-Dispersion = z120 "
          "geeicht). Ihr Vakuumsektor ist bei L <= 12 kritisch und su(2)_1-artig (c ~ 1, x ~ 1/2, v = 2 tanh(pi/2)); das "
          "chirale Gross-Neveu-Modell liegt erst im Doppellimes Theta -> oo. Drei Nähte tragen einen verbundenen Anteil mit "
          "festem Koeffizienten, c chi mit c = +-4 arccos(3/4)/sqrt 7 (Vorzeichen = Reihenfolge), aber im zeitungeraden "
          "epsilon-Kanal und ohne "
          "1->2-Amplitude: lambda = 0, gamma = 0, keine Box-Nachkommen. Für K15 heißt das: Die kohärente Komposition liefert eine "
          "parameterfreie, integrable Wechselwirkung (nichtabelsch, f-artig), aber keine Teilchenerzeugung und keinen "
          "phi^3-Vertex; die Szenarioachse 'Wechselwirkung' bekommt den Wert 'bedingt, integrabel, ohne Erzeugung'.")
say("\nUrteil\n   " + urteil)
say(f"\n   (Laufzeit {time.time() - T0:.1f} s)")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z121_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
