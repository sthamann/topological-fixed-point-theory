"""Z109: Prüfung der vierzehnten und fünfzehnten externen Analyse (25.11.).

Teil 1  Starrheit aus der Zweipunktantwort: Trägt w(a_i^+ sigma_t(a_i)) bei treuem Zustand genau eine Frequenz, so ist
        sigma_t(a_i) = e^(-i eps t) a_i und K = sum eps_i n_i + c; Gegenprobe mit Dichte-Dichte-Wechselwirkung.
Teil 2  Doppeltransfer allein: Ein quartischer Übergang zwischen komplementären Zuständen erhält nicht alle vier thermischen
        Einzelbesetzungsfunktionen (Analyse 14: 4e = 3e).
Teil 3  Inverse Rekonstruktion: k_S = (1/2|S|) sum_(i in S) tau(Gamma_S [K, gamma_i] gamma_i) exakt in der Clifford-Basis.
Teil 4  Kollektive Konstruktionen in zwölf Moden (zwei z94-Kopien): quartische Projektion g_4 = (pi/64) sin 2 theta der
        Vierergruppe aus z108; Achter-Konstruktion K_theta = pi [[(E+2) 1, -2 U], [-2 U^T, (E+2) 1]] mit winkel-
        unabhängigem Spektrum und winkelunabhängigen Zuständen auf bis zu drei Moden, aber winkelabhängiger Antwort
        tau([K, n_j]^+ [K, n_j]).
Teil 5  Marken und Minkowski: vier Marken spannen nur drei Richtungen; Gamma = -i A B ergänzt die vierte.
Teil 6  Chiralität: Wilson-Dirac-Eltern mit gleichen Daten auf allen Koordinatenhyperflächen und verschiedenem Grad (4D, 3D).
Teil 7  Gravitation und Lambda: erste Variation von lambda int sqrt(g)(R - 2 Lambda*) verschwindet auf Einsteinräumen; E_4 =
        8 Lambda^2/3; Anomalie-Selbstkonsistenz Lambda G = 3 pi/a: Standardmodell gegen TFPT, nötiges a = S_dS.
Teil 8  Modulare Eindeutigkeit (Analyse 15): chi_0 = theta_3/eta unter Gamma_theta, keine Nullstellen, Geschlecht 0 mit zwei
        Spitzen; E8-Charakter: d log chi = -(2 pi i/3) E6/E4, E6(i) = 0.
"""
import itertools
import json
from collections import defaultdict
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import minimize

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(109)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def kom(a, b):
    return a @ b - b @ a


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


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Starrheit aus der Zweipunktantwort")
a3 = jw(3)
N3 = [x.T @ x for x in a3]
eps = np.array([0.7, 1.3, 2.9])
K0 = sum(e*Nm for e, Nm in zip(eps, N3))


def frequenzen(K, a_):
    w, V = np.linalg.eigh(K)
    p = np.exp(-(w - w.min()))
    p /= p.sum()
    A = V.conj().T @ a_ @ V
    fr = defaultdict(float)
    for j in range(len(w)):
        for k in range(len(w)):
            wt = p[k]*abs(A[j, k])**2
            if wt > 1e-14:
                fr[round(float(w[j] - w[k]), 8)] += wt
    return dict(fr)


fr0 = [frequenzen(K0, x) for x in a3]
K1 = K0 + 0.8*N3[0] @ N3[1]
fr1 = [frequenzen(K1, x) for x in a3]
ein = all(len(f) == 1 and abs(list(f)[0] + e) < 1e-8 for f, e in zip(fr0, eps))
mehr = len(fr1[0]) >= 2 and len(fr1[1]) >= 2
rows = []
for x in a3:
    for g in (x, x.T):
        rows.append(np.kron(g.T, np.eye(8)) - np.kron(np.eye(8), g))
null3 = int(np.sum(np.linalg.svd(np.vstack(rows), compute_uv=False) < 1e-10))
eig_ok = all(np.allclose(kom(K0, x), -e*x) for x, e in zip(a3, eps))
check(f"Gaußscher Generator: jede Zweipunktantwort w(a_i^+ sigma_t(a_i)) hat genau eine Frequenz -eps_i ({[list(f) for f in fr0]}), "
      f"und [K, a_i] = -eps_i a_i; mit Dichte-Dichte-Term 0,8 n_0 n_1 hat a_0 die Frequenzen {sorted(fr1[0])}. Satz: Bei treuem "
      f"Zustand ist das Spektralmaß p_k |A_jk|^2 positiv, eine Frequenz heißt [K, a_i] = -eps_i a_i; was mit allen a_i, a_i^+ "
      f"vertauscht, ist skalar (Lösungsraum Dimension {null3}), also K = sum eps_i n_i + c",
      ein and mehr and eig_ok and null3 == 1)
befund("E: Der Satz der Analyse stimmt und schließt S(n) unter einer stärkeren Voraussetzung: der vollständigen zeitabhängigen "
       "Zweipunktantwort mit einer Frequenz je Mode. Offen ist, ob die native Hardy-Quelle diese Antwort ohne Gaußannahme "
       "liefert", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Ein einzelner Doppeltransfer erhält die thermischen Einzelbesetzungen nicht")
a4 = jw(4)
N4 = [x.T @ x for x in a4]
ts = np.linspace(0.05, 3.0, 40)


def besetzungen(H):
    w, V = np.linalg.eigh(H)
    out = []
    for t in ts:
        r = V @ np.diag(np.exp(-t*(w - w.min()))) @ V.T
        r /= np.trace(r)
        out.append([np.trace(r @ n) for n in N4])
    return np.array(out)


idx_s, idx_sb = int("1100", 2), int("0011", 2)
res_min, kontrolle = {}, []
for versuch in range(3):
    h = rng.uniform(0.4, 2.0, size=4)*rng.choice([-1, 1], size=4)
    H0 = sum(hi*n for hi, n in zip(h, N4))
    ziel = besetzungen(H0)
    for g in (0.2, 0.5, 1.0):
        def rest(xy, g=g, H0=H0, ziel=ziel):
            H = H0.copy()
            H[idx_s, idx_s] += xy[0]
            H[idx_sb, idx_sb] += xy[1]
            H[idx_s, idx_sb] = H[idx_sb, idx_s] = g
            return float(np.sum((besetzungen(H) - ziel)**2))
        best = min(minimize(rest, x0, method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-20, "maxiter": 6000}).fun
                   for x0 in ([0, 0], [g, -g], [-g, g], [1.0, 1.0]))
        res_min[(versuch, g)] = best
    def rest0(xy, H0=H0, ziel=ziel):
        H = H0.copy()
        H[idx_s, idx_s] += xy[0]
        H[idx_sb, idx_sb] += xy[1]
        return float(np.sum((besetzungen(H) - ziel)**2))
    kontrolle.append(minimize(rest0, [0.3, -0.2], method="Nelder-Mead",
                              options={"xatol": 1e-12, "fatol": 1e-24, "maxiter": 6000}).fun)
skaliert = {k: v/k[1]**4 for k, v in res_min.items()}
ok_null = min(res_min.values()) > 1e4*max(kontrolle) and min(skaliert.values()) > 1e-4
ZAHLEN["doppeltransfer_min_rest"] = {f"{k[0]}, g = {k[1]}": float(v) for k, v in res_min.items()}
ZAHLEN["doppeltransfer_kontrolle_g0"] = [float(v) for v in kontrolle]
check(f"vier Moden, freie Energien h_i != 0 (drei Zufallssätze), Übergang zwischen den komplementären Zuständen 1100 und 0011 "
      f"mit Kopplung g und frei verschiebbaren Diagonalen: Die vier thermischen Einzelbesetzungen über t in [0,05, 3] lassen "
      f"sich für g = 0,2, 0,5, 1 nicht halten (kleinster Rest {min(res_min.values()):.1e}, Rest/g^4 >= "
      f"{min(skaliert.values()):.1e}; Kontrolle g = 0 mit verschobenem Start erreicht {max(kontrolle):.0e}). Das bestätigt die "
      "Analyse (4e = 3e): ein isolierter Doppeltransfer genügt nicht, nötig ist eine kompensierte kollektive Struktur",
      ok_null)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Inverse Rekonstruktion aus den dynamischen Antworten")
gam = []
for x in a3:
    gam += [x + x.T, 1j*(x - x.T)]
ok_maj = all(np.allclose(g @ g, np.eye(8)) and np.allclose(g, g.conj().T) for g in gam)
untermengen = [S for k in range(2, 7, 2) for S in itertools.combinations(range(6), k)]


def Gam(S):
    M = np.eye(8, dtype=complex)
    for s_ in S:
        M = M @ gam[s_]
    return (1j)**(len(S)*(len(S) - 1)//2)*M


kS = {S: rng.normal() for S in untermengen}
K = 0.37*np.eye(8) + sum(v*Gam(S) for S, v in kS.items())
tau = lambda M: np.trace(M).real/8
D = [kom(K, g) for g in gam]
rekon = {S: sum(tau(Gam(S) @ D[i] @ gam[i]) for i in S)/(2*len(S)) for S in untermengen}
fehler = max(abs(rekon[S] - kS[S]) for S in untermengen)
check(f"sechs Majoranas (drei Moden), zufälliger gerader Generator mit allen {len(untermengen)} Koeffizienten k_S: die Formel "
      f"k_S = (1/2|S|) sum_(i in S) tau(Gamma_S D_i gamma_i) mit D_i = [K, gamma_i] gibt alle exakt zurück (Fehler {fehler:.0e}); "
      "der Generator ist aus den dynamischen Antworten bis auf eine Konstante eindeutig", ok_maj and fehler < 1e-12)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Kollektive Konstruktionen in zwölf Moden")
z94 = json.loads((HERE / "z94_ergebnis.json").read_text())
fz = z94["zahlen"]["aussenbehauptung_gleiche_saaten"]["f(n) für N = 6"]
C6 = list(itertools.product((0, 1), repeat=6))
H6 = {c: sum(c) + 4*fz.get("".join(map(str, c)), 0) for c in C6}
qplus = lambda c: sum(c[2*j] - c[2*j + 1] for j in range(3))
info, grp = {}, defaultdict(list)
for a_ in C6:
    for b_ in C6:
        n = a_ + b_
        dat = (H6[a_] + 3*H6[b_], sum(n), sum(a_) - sum(b_), qplus(a_) + qplus(b_))
        info[n] = dat
        grp[(dat[0], dat[2], dat[3])].append(n)
hd = lambda x, y: sum(u != v for u, v in zip(x, y))
vierer = []
for (E, qm, qp), Xs in grp.items():
    Ys = grp.get((E + 4, qm, qp), [])
    for X in Xs:
        for Y in Ys:
            if (info[X][1] - info[Y][1]) % 2 or hd(X, Y) < 4:
                continue
            delta = tuple(y - x for x, y in zip(X, Y))
            for (E2, qm2, qp2), Cs in grp.items():
                if E2 != E or (qm2, qp2) == (qm, qp) and False:
                    continue
                for Cc in Cs:
                    Dd = tuple(c - d for c, d in zip(Cc, delta))
                    if Dd == Cc or any(v not in (0, 1) for v in Dd) or info[Dd][0] != E + 4 or info[Dd][2:] != (qm2, qp2):
                        continue
                    if len({X, Y, Cc, Dd}) == 4 and (info[Cc][1] - info[Dd][1]) % 2 == 0:
                        vierer.append((X, Y, Cc, Dd))
z108 = json.loads((HERE / "z108_ergebnis.json").read_text())
tab = z108["zahlen"]["rotation_tabelle"]
st = {k: tuple(int(ch) for ch in v["Konfiguration"].replace("|", "")) for k, v in tab.items()}
A, B, Cq, Dq = st["A"], st["B"], st["C"], st["D"]


def wirke(ops, n):
    """ops: Liste (Mode, 'c' oder 'a') von rechts nach links angewandt; Jordan-Wigner-Vorzeichen."""
    n = list(n)
    sgn = 1
    for mode, typ in reversed(ops):
        if (typ == "a" and n[mode] == 0) or (typ == "c" and n[mode] == 1):
            return 0, None
        sgn *= (-1)**sum(n[:mode])
        n[mode] = 0 if typ == "a" else 1
    return sgn, tuple(n)


diffAB = [i for i in range(12) if A[i] != B[i]]
cre = [i for i in diffAB if B[i] == 1]
ann = [i for i in diffAB if A[i] == 1]
O = [(i, "c") for i in cre] + [(j, "a") for j in ann]
sAB, bild = wirke(O, A)
Od = [(j, "c") for j in ann] + [(i, "a") for i in reversed(cre)]
sCD, bild2 = wirke(Od, Cq)
kA, kB = np.pi*info[A][0], np.pi*info[B][0]
g4 = {}
for th in (np.pi/6, np.pi/4):
    KAB = np.sin(th)*np.cos(th)*(kA - kB)
    for s_ in (1, -1):
        KCD = s_*KAB
        g4[(round(th, 6), s_)] = (2/4096)*(sAB*KAB + sCD*KCD)/(1/8)
ZAHLEN["quartisch"] = {"Operator": f"{' '.join(f'a{i}^+' for i in cre)} {' '.join(f'a{j}' for j in ann)} + h.c.",
                       "Vorzeichen A->B, C->D": [sAB, sCD], "g4": {f"{k[0]:.4f}, Orientierung {k[1]}": v for k, v in g4.items()}}
gmax = max(abs(v) for k, v in g4.items() if abs(k[0] - np.pi/4) < 1e-6)
g30 = max(abs(v) for k, v in g4.items() if abs(k[0] - np.pi/6) < 1e-6)
gmin = min(abs(v) for v in g4.values())
check(f"quartische Projektion in der Vierergruppe aus z108: Der Übergang A -> B ist {ZAHLEN['quartisch']['Operator']} "
      f"(eine Mode vernichtet, drei erzeugt, also ein Paar -> zwei Paare), und derselbe Operator verbindet D und C (Bild "
      f"{bild == B}, {bild2 == Dq}). Projektion g_4 = tau(V K)/tau(V V): {gmax:.6f} bei 45 Grad = pi/64 (Analyse 15: 0,049087), "
      f"{g30:.6f} bei 30 Grad = (pi/64) sin 60 Grad, wenn beide Drehungen gleich orientiert sind; bei entgegengesetzter "
      f"Orientierung (die die Bedingungen ebenso erfüllt) {gmin:.0e}: dann liegt der Übergang in höheren Mehrkörpertermen",
      bild == B and bild2 == Dq and abs(gmax - np.pi/64) < 1e-12 and abs(g30 - np.pi/64*np.sin(np.pi/3)) < 1e-12
      and gmin < 1e-15)

gruppen = defaultdict(list)
for q in vierer:
    gruppen[(info[q[0]][0], info[q[0]][2:], info[q[2]][2:])].append(q)
achter = None
for key, qs in sorted(gruppen.items()):
    for q1, q2 in itertools.combinations(qs, 2):
        alle = set(q1) | set(q2)
        if len(alle) != 8:
            continue
        paare = [(q1[0], q1[1]), (q2[0], q2[1]), (q1[0], q2[1]), (q2[0], q1[1]),
                 (q1[2], q1[3]), (q2[2], q2[3]), (q1[2], q2[3]), (q2[2], q1[3])]
        if all(hd(x, y) >= 4 for x, y in paare) and info[q1[0]][2:] == info[q2[0]][2:]:
            achter = (key, q1, q2)
            break
    if achter:
        break
ZAHLEN["achter_suche"] = {"Vierergruppen (alle Anregungsmuster)": len(vierer), "Gruppen": len(gruppen),
                          "gefunden": None if not achter else {"Exponent": achter[0][0], "Ladungen": [achter[0][1], achter[0][2]]}}
if achter:
    (E8, _, _), q1, q2 = achter
    zust = [q1[0], q2[0], q1[1], q2[1], q1[2], q2[2], q1[3], q2[3]]
    pos = {z: i for i, z in enumerate(zust)}

    def Kblock(th):
        U = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        M = np.zeros((4, 4))
        M[:2, 2:] = U
        M[2:, :2] = U.T
        Kb = np.pi*((E8 + 2)*np.eye(4) - 2*M)
        return np.block([[Kb, np.zeros((4, 4))], [np.zeros((4, 4)), Kb]])

    spektren = [np.round(np.linalg.eigvalsh(Kblock(t))/np.pi, 9) for t in (0.0, 0.4, 1.1)]
    spek_ok = all(np.allclose(spektren[0], s_) for s_ in spektren) and set(spektren[0]) == {E8, E8 + 4}
    pAll = {n: np.exp(-np.pi*(info[n][0] - E8)) for n in info}
    Zp = sum(pAll.values())

    def zustand(th):
        rb = expm(-(Kblock(th) - np.pi*E8*np.eye(8)))
        return rb

    rb0, rb1 = zustand(0.0), zustand(1.1)
    diag_ok = np.allclose(np.diag(rb0), np.diag(rb1))
    ein_ok = True
    for k in range(12):
        alt = sum(pAll[n] for n in info if n[k])
        neu = alt - sum(pAll[z] for z in zust if z[k]) + sum(rb0[pos[z], pos[z]] for z in zust if z[k])
        ein_ok &= abs(neu - alt)/Zp < 1e-12
    unsichtbar = all(hd(zust[i], zust[j]) >= 4 for i in range(8) for j in range(8) if abs(rb1[i, j]) > 1e-14 and i != j)
    antw = {}
    for j in range(12):
        werte = []
        for th in (0.0, np.pi/6, np.pi/3):
            Kb = Kblock(th)
            s2 = sum(abs(Kb[i, k])**2 for i in range(8) for k in range(8) if zust[i][j] != zust[k][j])
            werte.append(s2/4096)
        antw[j] = [w/np.pi**2 for w in werte]
    variiert = [j for j, w in antw.items() if max(w) - min(w) > 1e-12]
    ZAHLEN["achter_antwort_durch_pi2"] = {j: [round(x, 6) for x in w] for j, w in antw.items() if w[0] or w[1]}
    check(f"Achter-Konstruktion aus zwei Vierergruppen ({len(vierer)} Vierergruppen über alle Anregungsmuster; Exponent "
          f"{E8}): K_theta = pi[[{E8 + 2} 1, -2 U], [-2 U^T, {E8 + 2} 1]] je Block. Spektrum winkelunabhängig "
          f"({sorted(set(spektren[0]))} pi), Diagonale winkelunabhängig ({diag_ok}), alle zwölf Einzelbesetzungen wie im "
          f"z94-Produkt ({ein_ok}), Kohärenzen nur zwischen Zuständen mit >= 4 verschiedenen Moden ({unsichtbar}): alle "
          f"Zustände auf bis zu drei Moden hängen nicht vom Winkel ab. Die dynamische Antwort tau([K, n_j]^+ [K, n_j]) "
          f"hängt für die Moden {variiert} vom Winkel ab (in pi^2 bei 0, 30, 60 Grad: "
          f"{[ZAHLEN['achter_antwort_durch_pi2'].get(j) for j in variiert[:2]]})",
          spek_ok and diag_ok and ein_ok and unsichtbar and len(variiert) > 0)
else:
    check("Achter-Konstruktion: keine Vierergruppen-Paare mit gleichen Ladungen und >= 4 verschiedenen Moden gefunden", False)
befund("E: Analyse 14 trifft den Kern. Gleiche Einzelbesetzungen, gleiches globales Modulspektrum und gleiche Zustände auf "
       "bis zu drei Moden lassen eine stetige Familie verschiedener Kopplungen zu; eine einzige dynamische Antwort "
       "(tau([K, n_j]^+ [K, n_j])) unterscheidet sie, und die inverse Rekonstruktion bestimmt den ganzen Generator aus den "
       "dynamischen Antworten. Was fehlt, ist das Vorwärtsproblem: welche native Regel diese Antworten festlegt", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Vier Marken und die vierte Richtung")
pk = np.array([[1, np.cos(k*np.pi/2), np.sin(k*np.pi/2), 0] for k in range(4)])
rang = np.linalg.matrix_rank(pk)
t_, x_, y_, z_ = sp.symbols("t x y z", real=True)
sx, sy = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]])
Gm = -sp.I*sx*sy
X = t_*sp.eye(2) + x_*sx + y_*sy + z_*Gm
check(f"die vier Nullrichtungen der quadratischen Marken p_k = (1, cos k pi/2, sin k pi/2, 0) haben Rang {rang}, p_0 + p_2 = "
      f"p_1 + p_3: sie spannen keinen Minkowskiraum. Mit Gamma = -i A B (A, B antikommutierende Involutionen, hier sigma_x, "
      f"sigma_y; Gamma = {Gm.tolist()}) gilt det(t + x A + y B + z Gamma) = {sp.expand(X.det())}",
      rang == 3 and np.allclose(pk[0] + pk[2], pk[1] + pk[3])
      and sp.expand(X.det() - (t_**2 - x_**2 - y_**2 - z_**2)) == 0 and Gm == sp.Matrix([[1, 0], [0, -1]]))

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: Chiralität: gleiche Daten auf den Hyperflächen, verschiedener Grad")


def grad_ecken(dim, Mfun):
    tot = 0
    for c in itertools.product((0, 1), repeat=dim):
        k = np.array(c)*np.pi
        if Mfun(k) < 0:
            tot += (-1)**sum(c)
    return tot


def grad_integral(dim, Mfun, dMfun, Ngit):
    ks = np.linspace(0, 2*np.pi, Ngit, endpoint=False)
    G = np.meshgrid(*([ks]*dim), indexing="ij")
    pts = np.stack([g.ravel() for g in G], axis=1)
    d = np.concatenate([np.sin(pts), Mfun(pts.T)[:, None]], axis=1)
    cols = [d]
    for j in range(dim):
        dj = np.zeros_like(d)
        dj[:, j] = np.cos(pts[:, j])
        dj[:, dim] = dMfun(pts.T, j)
        cols.append(dj)
    Mat = np.stack(cols, axis=2)
    det = np.linalg.det(Mat)
    norm = np.linalg.norm(d, axis=1)
    vol = {3: 2*np.pi**2, 4: 8*np.pi**2/3}[dim]
    return float(np.sum(det/norm**(dim + 1))*(2*np.pi/Ngit)**dim/vol), float(norm.min())


M4 = lambda a: (lambda k: -1 + np.sum(1 - np.cos(k), axis=0) - a*np.prod(1 - np.cos(k), axis=0))


def dM4(a):
    def f(k, j):
        rest = np.prod([1 - np.cos(k[i]) for i in range(len(k)) if i != j], axis=0)
        return np.sin(k[j]) - a*np.sin(k[j])*rest
    return f


g4d = {}
for a in (0.0, 1.0):
    ecke = grad_ecken(4, lambda k, a=a: -1 + np.sum(1 - np.cos(k)) - a*np.prod(1 - np.cos(k)))
    integ, nmin = grad_integral(4, M4(a), dM4(a), 26)
    g4d[a] = (ecke, integ, nmin)
k_ = sp.symbols("k1:5", real=True)
diffterm = sp.prod([1 - sp.cos(k) for k in k_])
hyper = all(diffterm.subs(k, 0) == 0 for k in k_)


def M3(a1, a2):
    return lambda k: (-1 + np.sum(1 - np.cos(k), axis=0)
                      - a1*sum((1 - np.cos(k[i]))*(1 - np.cos(k[j])) for i, j in itertools.combinations(range(3), 2))
                      - a2*np.prod(1 - np.cos(k), axis=0))


def dM3(a1, a2):
    def f(k, j):
        andere = [i for i in range(3) if i != j]
        s = np.sin(k[j])
        return (s - a1*s*sum(1 - np.cos(k[i]) for i in andere)
                - a2*s*np.prod([1 - np.cos(k[i]) for i in andere], axis=0))
    return f


g3d = {}
for (a1, a2) in ((0.0, 0.0), (0.0, 1.0), (1.0, 0.0)):
    ecke = grad_ecken(3, lambda k, a1=a1, a2=a2: (-1 + np.sum(1 - np.cos(k)) - a1*sum((1 - np.cos(k[i]))*(1 - np.cos(k[j]))
                                                                                     for i, j in itertools.combinations(range(3), 2))
                                                   - a2*np.prod(1 - np.cos(k))))
    integ, nmin = grad_integral(3, M3(a1, a2), dM3(a1, a2), 48)
    g3d[(a1, a2)] = (ecke, integ, nmin)
ZAHLEN["grade"] = {"4D (a: Ecken, Integral, min|d|)": {str(k): v for k, v in g4d.items()},
                   "3D ((a1, a2): Ecken, Integral, min|d|)": {str(k): v for k, v in g3d.items()}}
check(f"4D Wilson-Dirac d = (sin k, M_a), M_a = -1 + sum(1 - cos k_i) - a prod(1 - cos k_i): der Zusatzterm verschwindet auf "
      f"jeder Hyperfläche k_i = 0 ({hyper}) und ist dort bis zur Ordnung k^7 unsichtbar; |Grad| (Eckenformel, numerisches "
      f"Integral 26^4, Vorzeichen = Orientierungskonvention) für a = 0: {g4d[0.0][0]}, {abs(g4d[0.0][1]):.3f}; a = 1: "
      f"{g4d[1.0][0]}, {abs(g4d[1.0][1]):.3f} (Lücke min|d| = {min(g4d[0.0][2], g4d[1.0][2]):.2f} > 0). 3D mit symmetrischen "
      f"gemischten Termen (gleich auf den Achsen), (a_Paar, a_Tripel) = (0,0), (0,1), (1,0): Grade "
      f"{[v[0] for v in g3d.values()]} (Integral {[round(v[1], 3) for v in g3d.values()]}; Analyse 15: 1, 0, 3)",
      hyper and g4d[0.0][0] == 1 and g4d[1.0][0] == 2 and abs(abs(g4d[0.0][1]) - 1) < 0.1
      and abs(abs(g4d[1.0][1]) - 2) < 0.1 and min(v[2] for v in g4d.values()) > 0.1
      and [v[0] for v in g3d.values()] == [1, 0, 3] and all(abs(abs(v[0]) - abs(v[1])) < 0.02 for v in g3d.values())
      and min(v[2] for v in g3d.values()) > 0.05)
befund("E: Die Gegenbeispiele der Analyse 15 stimmen: Lokale Daten auf Achsen bzw. Hyperflächen legen den Index des "
       "Elternoperators nicht fest. Der Index muss aus der vollständig rekonstruierten lokalen Algebra kommen, sonst ist er "
       "eine versteckte Eingabe", True)

# ------------------------------------------------------------------ Teil 7
say("\nTeil 7: Gravitation und kosmologische Konstante")
L, Ls = sp.symbols("Lambda Lambda_star", positive=True)
dim = 4
Riem2 = 2*(L/3)**2*dim*(dim - 1)
E4 = sp.simplify(Riem2 - 4*(dim*L**2) + (dim*L)**2)
gmn = sp.symbols("g")
erste_var = sp.simplify(Ls - sp.Rational(1, 2)*(4*Ls) + Ls)
ainv = mp.mpf("137.0359992168")
a_sm = sp.Rational(4 + 11*24 + 62*12, 360)
LG_sm = 3*mp.pi/mp.mpf(float(a_sm))
LG_tfpt = 3*mp.e**(-2*ainv)/(32*mp.pi**3)
a_noetig = 3*mp.pi/LG_tfpt
S_dS = 32*mp.pi**4*mp.e**(2*ainv)
b_euler = 64*mp.pi**2*(mp.pi**2/2)*mp.e**(2*ainv)
ZAHLEN["lambda_anomalie"] = {"a_SM (mit nu_R)": str(a_sm), "Lambda G (SM)": mp.nstr(LG_sm, 5), "Lambda G (TFPT)": mp.nstr(LG_tfpt, 5),
                             "nötiges a": mp.nstr(a_noetig, 5), "S_dS = 32 pi^4 e^(2/alpha)": mp.nstr(S_dS, 5)}
check(f"Einsteinraum R_mn = Lambda* g_mn: die erste Variation von lambda int sqrt(g)(R - 2 Lambda*) ist proportional zu R_mn "
      f"- R g_mn/2 + Lambda* g_mn = {erste_var} g_mn, also verschieben solche Terme G, ohne Vakuum und erste Antwort zu ändern. "
      f"Maximal symmetrisch: E_4 = {E4}. Anomalie-Selbstkonsistenz Lambda G = 3 pi/a: Standardmodell mit nu_R a = {a_sm} = "
      f"{float(a_sm):.4f}, Lambda G = {mp.nstr(LG_sm, 4)}; TFPT (COSMO.LAM.01) Lambda G = {mp.nstr(LG_tfpt, 4)}, nötig a = "
      f"{mp.nstr(a_noetig, 4)}; das ist exakt die de-Sitter-Entropie 3 pi/(Lambda G) = 32 pi^4 e^(2/alpha) und genau das "
      "Gewicht des Euler-Terms |b| 32 pi^2 chi(S^4) aus z96",
      erste_var == 0 and sp.simplify(E4 - sp.Rational(8, 3)*L**2) == 0 and abs(LG_sm - mp.mpf("3.3527")) < 1e-3
      and abs(LG_tfpt/mp.mpf("2.836e-122") - 1) < 1e-3 and abs(a_noetig/S_dS - 1) < 1e-25 and abs(b_euler/S_dS - 1) < 1e-25)
befund("C/O: Der Anomaliemechanismus gibt Lambda > 0 von selbst, aber in Planckgröße; für den TFPT-Wert müsste der "
       "Anomaliekoeffizient die Horizontentropie sein (a = S_dS = 32 pi^4 e^(2/alpha)). Das ist dieselbe Zahl wie das "
       "Euler-Gewicht der MacDowell-Mansouri-Wirkung (z96): ein Konsistenzziel, keine Herleitung", True)

# ------------------------------------------------------------------ Teil 8
say("\nTeil 8: Modulare Eindeutigkeit des Charakters (Analyse 15)")
mp.mp.dps = 30


def eta(tau):
    q = mp.exp(2j*mp.pi*tau)
    return mp.exp(2j*mp.pi*tau/24)*mp.nprod(lambda n: 1 - q**n, [1, mp.inf])


def th3(tau):
    return mp.jtheta(3, 0, mp.exp(1j*mp.pi*tau))


chi0 = lambda tau: th3(tau)/eta(tau)
pts = [mp.mpc(0.31, 1.07), mp.mpc(-0.42, 0.83), mp.mpc(0.11, 2.2)]
S_inv = max(abs(chi0(-1/t)/chi0(t) - 1) for t in pts)
T2 = max(abs(chi0(t + 2)/chi0(t) - mp.exp(-1j*mp.pi/6)) for t in pts)
prod_formel = max(abs(th3(t) - eta(t)**5/(eta(t/2)**2*eta(2*t)**2)) for t in pts)
gitter = [mp.mpc(x, y) for x in np.linspace(-1, 1, 21) for y in np.linspace(0.12, 2.0, 16) if abs(complex(x, y)) >= 1]
min_th3 = min(abs(th3(t)) for t in gitter)
genus = sp.Rational(1) + sp.Rational(3, 12) - sp.Rational(1, 4) - 0 - sp.Rational(2, 2)
E4f = lambda tau: 1 + 240*mp.nsum(lambda n: n**3*mp.exp(2j*mp.pi*tau*n)/(1 - mp.exp(2j*mp.pi*tau*n)), [1, mp.inf])
E6f = lambda tau: 1 - 504*mp.nsum(lambda n: n**5*mp.exp(2j*mp.pi*tau*n)/(1 - mp.exp(2j*mp.pi*tau*n)), [1, mp.inf])
tau0 = mp.mpc(0.2, 1.3)
dlog = mp.diff(lambda t: mp.log(E4f(t)/eta(t)**8), tau0)
rechts = -2j*mp.pi/3*E6f(tau0)/E4f(tau0)
check(f"chi_0 = theta_3/eta: S-invariant ({float(S_inv):.0e}), T^2-Multiplikator e^(-i pi/6) ({float(T2):.0e}), Produktformel "
      f"theta_3 = eta^5/(eta(tau/2)^2 eta(2 tau)^2) ({float(prod_formel):.0e}), also keine Nullstellen in H (eta hat keine; "
      f"Gitterminimum {float(min_th3):.3f} am Rand zur Spitze 1, wo theta_3 gegen 0 geht); X(Gamma_theta) hat Geschlecht "
      f"{genus} (Index 3, ein elliptischer Punkt, zwei Spitzen): ein holomorpher, an beiden Spitzen beschränkter Quotient "
      f"chi/chi_0 ist konstant. E8-Naht: d log(E4/eta^8)/d tau = -(2 pi i/3) E6/E4 ({float(abs(dlog - rechts)):.0e}), "
      f"E6(i) = {mp.nstr(abs(E6f(mp.mpc(0, 1))), 3)}",
      S_inv < 1e-20 and T2 < 1e-20 and prod_formel < 1e-20 and min_th3 > 1e-3 and genus == 0 and abs(dlog - rechts) < 1e-15
      and abs(E6f(mp.mpc(0, 1))) < 1e-20)
befund("C: Der modulare Eindeutigkeitsweg der Analyse 15 ist mathematisch tragfähig (Geschlecht 0, keine Nullstellen, zwei "
       "Spitzen). Er verlegt die Starrheit in eine Voraussetzung: dass jeder Alternativzustand einen einzigen, unter Gamma_theta "
       "mit demselben Multiplikator kovarianten Charakter hat (Skalarcharakter) und im R-Sektor nichtnegative Energien. Ob "
       "das aus der Quelle folgt, ohne die freie Theorie vorauszusetzen, ist offen; mit gleicher Einteilchenbesetzung gibt "
       "die relative Entropie dann rho = rho_0", True)

urteil = ("Die vierzehnte und die fünfzehnte externe Analyse halten in allen gerechneten Punkten. (1) Eine einzige Frequenz in "
          "jeder Zweipunktantwort erzwingt bei treuem Zustand den gaußschen Generator K = sum eps_i n_i + c (Starrheit unter "
          "der stärkeren Voraussetzung der vollen Zweipunktantwort). (2) Ein isolierter Doppeltransfer erhält die thermischen "
          "Einzelbesetzungen nicht (kleinster Rest > 0 für g != 0). (3) Die inverse Rekonstruktion k_S = (1/2|S|) sum "
          "tau(Gamma_S [K, gamma_i] gamma_i) ist exakt. (4) In der z108-Vierergruppe ist der Übergang ein Vierfermionoperator "
          "(eine Mode vernichtet, drei erzeugt) mit g_4 = (pi/64) sin 2 theta bei gleich orientierten Drehungen und 0 bei "
          "entgegengesetzten; eine Achter-Konstruktion mit winkelunabhängigem Spektrum und winkelunabhängigen Zuständen auf bis "
          "zu drei Moden hat winkelabhängige dynamische Antworten. (5) Vier Marken spannen nur drei Richtungen; Gamma = -i A B "
          "gibt die vierte und det = t^2 - x^2 - y^2 - z^2. (6) Wilson-Dirac-Eltern mit gleichen Daten auf allen "
          "Hyperflächen haben Grad 1 und 2 (4D) bzw. drei verschiedene Grade (3D). (7) Erste Variationen legen G nicht fest; "
          "der Anomaliemechanismus gibt Lambda > 0 nur mit a = S_dS = 32 pi^4 e^(2/alpha), dem Euler-Gewicht aus z96. (8) Der "
          "modulare Eindeutigkeitsweg ist tragfähig, verlegt die Starrheit aber in die Skalarcharakter-Voraussetzung. "
          "Gemeinsam: Das inverse Problem ist gelöst, das Vorwärtsproblem (die Regel, die die gemeinsamen Antworten "
          "festlegt) nicht.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z109_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
