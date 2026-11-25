"""Z92: Folgefrage F4 'Materiespin verkleben und die fehlende Drei' -- zehnte und elfte externe Analyse geprüft.

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: DIMENSION.UPLIFT.FIREWALL.01 (jede 3+1-Aussage ist Kinematik der radialen
Einpaarschicht oder Lesart), GRAV.NONCIRCULAR.01 (kein Schluss modular => Raumzeit). Gelöste Bausteine (16, E8, Flavour) werden
benutzt, nicht neu geöffnet; eine richtige Dimension ist keine richtige Darstellung. Bestand: z82 (Träger e_m = z^m/(z^4 - 1),
R = diag(1, i, -1, -i, 1), R^2- und V4-Schnitt, das Original erreicht nur die 2, z_X), z83 (radiale Standardform, CCR-Schritte,
Helizität h = (m - mquer)/2, nativer Heber L^nat), z81 (R = i^(Q_a - Q_b), R^2 = (-1)^F), z79 (eta-Dublett), z34 (D4 auf dem
Träger, V4-Schnitt).

Teil 1  Die fehlende Drei (zehnte Analyse): Sym^2(W)/(uv)^2 mit W = span(u^3 v, u v^3) exakt, Eindeutigkeit des Divisors (ggT),
        D4-Äquivarianz bis auf den Charakter von (uv)^2. Darstellungstheorie: Äquivarianzalgebra in g_SM (R^2- und V4-Schnitt),
        Hyperladung, Schwachtriplett in Lambda^* C^5, SU(5)-äquivariante Abbildungen, E8 aus D5 + A3 (Wurzeln, Wurzelpaare).
        Native Paaroperatoren des Originals: Pauli, Ad_R-Charaktere, Verdrehung durch das niedrigste R = -1-Boson, U(1)^2-Ladungen.
Teil 2  Verklebung Spin x_Z2 Z4 im radialen Bild: 2 pi-Drehung auf der vollen Standardform (Helizitätssektoren), R auf der Leiter,
        native radiale Energie E^nat auf allen Ladungssektoren, Kovarianz der Vertexoperatoren, Familie U^osc c(Q) möglicher
        Fortsetzungen: was genau scheitert.
Teil 3  Elfte Analyse: Verzweigungsfelder der Z2-Permutationsorbifaltung an den Marken (h = c/16), Vierteldrehung und Uhr auf
        ihnen, radiale Helizität, Blätterwechsel auf der Standardform.
Teil 4  Bilanz und Urteil für F4.
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import sympy as sp
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(92)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def nullraum(A, tol=1e-9):
    _, s, vh = np.linalg.svd(A)
    rang = int(np.sum(s > tol*max(1.0, s[0] if len(s) else 1.0)))
    return vh[rang:].conj().T


# =================================================================== Teil 1
say("Teil 1: Die fehlende Drei -- Sym^2(W)/(uv)^2")

# --- 1.1 exakte Konstruktion
u, v, z = sp.symbols("u v z")
W_MON = [u**3*v, u*v**3]                                          # e_3, e_1 (e_m = u^m v^(4-m)/(u^4 - v^4) = z^m/(z^4 - 1))
SYM2 = [sp.expand(a*b) for a, b in itertools.combinations_with_replacement(W_MON, 2)]
GGT = sp.gcd_list(SYM2)
QUOT = [sp.cancel(p/GGT) for p in SYM2]
m_von = lambda mon: sp.degree(mon, u)                             # Index m des Trägermonoms u^m v^(4-m)
fix_ok = sp.expand(sp.I*u*v - u*v - (sp.I - 1)*u*v) == 0          # (u, v) -> (i u, v) fix genau auf uv = 0: Punkte 0, unendlich
quarten = all(sp.Poly(q, u, v).is_homogeneous and sp.Poly(q, u, v).total_degree() == 4 for q in QUOT)
bild = sorted(m_von(q) for q in QUOT)
gl = lambda mon: sp.I**sp.degree(mon, u)                          # kanonische Linearisierung: Rückzug z -> i z, R e_m = i^m e_m
om = sp.exp(sp.I*sp.pi/4)
slz = lambda mon: sp.simplify(om**(sp.degree(mon, u) - sp.degree(mon, v)))   # SL(2)-Hebung diag(om, om^-1)
ok_ch = True
for lift in (gl, slz):
    for p, q in zip(SYM2, QUOT):
        ok_ch &= sp.simplify(lift(p)/lift(GGT) - lift(q)) == 0
char_gl = [gl(q) for q in sorted(QUOT, key=m_von)]
twist = {"kanonisch (z -> i z, z82)": str(gl(GGT)), "SL(2)-Hebung": str(slz(GGT)),
         "Markenschnitt u^4 - v^4 in der SL(2)-Hebung": str(sp.simplify(sp.expand((om*u)**4 - (v/om)**4)/(u**4 - v**4)))}
E_MOD = [z**m/(z**4 - 1) for m in range(5)]
rat_ok = all(sp.cancel(E_MOD[i]*E_MOD[j]/E_MOD[2] - E_MOD[i + j - 2]) == 0 for i, j in ((1, 1), (1, 3), (3, 3)))
ZAHLEN["sym2_W"] = {"Sym^2 W": [str(p) for p in SYM2], "ggT": str(GGT), "Quotient": [str(q) for q in QUOT],
                    "Charaktere (e_0, e_2, e_4), kanonisch": [str(c) for c in char_gl], "Charakter von (uv)^2": twist}
check("Sym^2(W) mit W = span(u^3 v, u v^3) = {e_3, e_1} ist span(u^6 v^2, u^4 v^4, u^2 v^6); ihr ggT ist genau (uv)^2, und uv = 0 "
      "sind die Fixpunkte 0, unendlich der Vierteldrehung. Der Quotient ist span(u^4, u^2 v^2, v^4) = {e_4, e_2, e_0}, also "
      "genau der R^2 = +1-Eigenraum, als Funktionen e_1^2/e_2 = e_0, e_1 e_3/e_2 = e_2, e_3^2/e_2 = e_4 (exakt). Charaktere "
      "(e_0, e_2, e_4) = (+1, -1, +1); die Teilungsabbildung ist in beiden Hebungen äquivariant. Jeder quartische Teiler, "
      "der Sym^2 W in den Träger bringt, teilt den ggT, ist also c (uv)^2: der Divisor ist erzwungen, nicht gewählt",
      fix_ok and quarten and bild == [0, 2, 4] and ok_ch and [str(c) for c in char_gl] == ["1", "-1", "1"] and rat_ok
      and sp.Poly(GGT, u, v).total_degree() == 4 and sp.expand(GGT - u**2*v**2) == 0)

# --- 1.2 D4-Äquivarianz bis auf den Charakter von (uv)^2
D4 = {f"r^{k}": (lambda w, k=k: sp.I**k*w) for k in range(4)}
D4.update({f"r^{k} s": (lambda w, k=k: sp.I**k/w) for k in range(4)})


def rueckzug(g):
    M = sp.zeros(5, 5)
    for k in range(5):
        poly = sp.Poly(sp.cancel(E_MOD[k].subs(z, g(z))*(z**4 - 1)), z)
        assert poly.degree() <= 4
        for m in range(5):
            M[m, k] = poly.coeff_monomial(z**m)
    return M


def koord(f):
    poly = sp.Poly(sp.cancel(f*(z**4 - 1)), z)
    assert poly.degree() <= 4
    return sp.Matrix([poly.coeff_monomial(z**m) for m in range(5)])


def funk(c):
    return sum(c[m]*E_MOD[m] for m in range(5))


ok_d4, chi, spur = True, {}, {}
for nm, g in D4.items():
    G = rueckzug(g)
    chi[nm] = sp.simplify(sp.cancel(E_MOD[2].subs(z, g(z))/E_MOD[2]))
    ok_d4 &= chi[nm].free_symbols == set()
    for i, j in ((1, 1), (1, 3), (3, 3)):
        links = koord(funk(G[:, i])*funk(G[:, j])/E_MOD[2])
        rechts = chi[nm]*G*koord(E_MOD[i]*E_MOD[j]/E_MOD[2])
        ok_d4 &= sp.simplify(links - rechts) == sp.zeros(5, 1)
    Gw = G.extract([1, 3], [1, 3])
    s2 = sp.simplify((Gw.trace()**2 + (Gw*Gw).trace())/2)          # tr Sym^2 g = (tr(g)^2 + tr(g^2))/2
    t3 = sp.simplify(G.extract([0, 2, 4], [0, 2, 4]).trace())
    spur[nm] = (str(s2), str(t3), str(chi[nm]))
    ok_d4 &= sp.simplify(s2 - chi[nm]*t3) == 0
mult_sym2 = sorted(str(sp.simplify(c)) for c in (sp.I**2, sp.I**4, sp.I**6))
mult_drei = sorted(str(sp.I**m) for m in (0, 2, 4))
ZAHLEN["D4_spuren (Sym^2 W, Drei, chi)"] = spur
check("D4 (Vierteldrehung r, Spiegelung s: z -> 1/z, Rückzug auf dem Träger, alle 8 Elemente exakt): die Teilungsabbildung "
      "T(a, b) = a b/e_2 erfüllt T(g a, g b) = chi(g) g T(a, b) mit dem Charakter chi von e_2 = (uv)^2 (chi(r) = chi(s) = -1), "
      "und tr_(Sym^2 W)(g) = chi(g) tr_3(g) für alle g. Als Z4-Moduln: Sym^2 W hat die Charaktere {-1, +1, -1}, die Drei "
      "{+1, -1, +1}; gleich erst nach der Verdrehung mit chi (kanonische Linearisierung; in der SL(2)-Hebung trägt statt "
      "(uv)^2 der Markenschnitt u^4 - v^4 das Vorzeichen)",
      ok_d4 and chi["r^1"] == -1 and chi["r^0 s"] == -1 and mult_sym2 == ["-1", "-1", "1"] and mult_drei == ["-1", "1", "1"])

# --- 1.3 Darstellungstheorie in g_SM: Äquivarianzalgebra, Hyperladung, Schwachtriplett
E5 = np.eye(5)


def T_num(a, b):
    out = np.zeros(5, complex)
    for i in (1, 3):
        for j in (1, 3):
            out[i + j - 2] += a[i]*b[j]
    return out


def gsm_basis(Wb, Cb):
    """komplexifizierte g_SM: sl(3) auf der Farbe, sl(2) auf dem Dublett, Y = -1/3 P_Farbe + 1/2 P_schwach (letztes Element)."""
    B = [Cb @ (np.eye(3)[:, [i]] @ np.eye(3)[[j], :]) @ Cb.conj().T for i in range(3) for j in range(3) if i != j]
    B += [Cb @ np.diag(d) @ Cb.conj().T for d in ([1, -1, 0], [0, 1, -1])]
    B += [Wb @ Mw @ Wb.conj().T for Mw in (np.array([[0, 1], [0, 0]]), np.array([[0, 0], [1, 0]]), np.diag([1, -1]))]
    B.append(-1/3*Cb @ Cb.conj().T + 1/2*Wb @ Wb.conj().T)
    return B


PW = E5[:, [1, 3]] @ E5[[1, 3], :]
P3 = E5[:, [0, 2, 4]] @ E5[[0, 2, 4], :]


def aequivarianz(B):
    zeilen = []
    for X in B:
        r = [((np.eye(5) - PW) @ X @ PW).ravel(), ((np.eye(5) - P3) @ X @ P3).ravel()]
        for a, b in itertools.product((1, 3), repeat=2):
            r.append(X @ T_num(E5[a], E5[b]) - T_num(PW @ X @ E5[a], E5[b]) - T_num(E5[a], PW @ X @ E5[b]))
        zeilen.append(np.concatenate(r))
    return nullraum(np.array(zeilen).T)


SCHNITTE = {"R^2-Schnitt": (E5[:, [1, 3]], E5[:, [0, 2, 4]]),
            "V4-Schnitt (z34)": (np.array([(E5[0] - E5[4])/np.sqrt(2), E5[2]]).T,
                                 np.array([E5[1], E5[3], (E5[0] + E5[4])/np.sqrt(2)]).T)}
aeq = {}
for nm, (Wb, Cb) in SCHNITTE.items():
    N = aequivarianz(gsm_basis(Wb, Cb))
    schwach_rang = np.linalg.matrix_rank(N[8:11, :], tol=1e-9) if N.shape[1] else 0
    aeq[nm] = {"dim": int(N.shape[1]), "Y-Anteil max": float(np.abs(N[11, :]).max()) if N.shape[1] else 0.0,
               "Rang des sl(2)-Anteils": int(schwach_rang)}
Wv, Cv = SCHNITTE["V4-Schnitt (z34)"]
P3v = P3
schnitt_v4 = {"mit dem Dublett": int(np.linalg.matrix_rank(Wv.conj().T @ P3v, tol=1e-9)),
              "mit der Farbe": int(np.linalg.matrix_rank(Cv.conj().T @ P3v, tol=1e-9))}
ZAHLEN["aequivarianzalgebra_in_g_SM (dim 12)"] = aeq
ZAHLEN["V4: Bild {e_0, e_2, e_4} geschnitten"] = schnitt_v4
Y_quelle, Y_ziel, Y_e2 = Fr(1, 2) + Fr(1, 2), Fr(-1, 3), Fr(-1, 3)
ZAHLEN["hyperladung"] = {"Y(Sym^2 W)": str(Y_quelle), "Y(Drei)": str(Y_ziel), "nötige Y des Divisors": str(Y_quelle - Y_ziel),
                         "Y(e_2) = Y((uv)^2-Schnitt)": str(Y_e2)}
check("[negativ, E] Darstellung: von der 12-dimensionalen g_SM bleibt im R^2-Schnitt eine 3-dimensionale Äquivarianzalgebra "
      "(die diagonale sl(2): Schwachdrehungen, auf der Drei als Sym^2 = so(3) gesperrt), ohne Y-Anteil; im V4-Schnitt bleibt "
      "1 Dimension, und das Bild {e_0, e_2, e_4} zerfällt dort in Dublett (2) und eine Farbrichtung (1). Hyperladung: "
      "Y(Sym^2 W) = 1 statt -1/3; die Teilung bräuchte einen Divisor mit Y = 4/3, der Schnitt (uv)^2 = e_2 hat Y = -1/3. Die "
      "Konstruktion ist (1, 3)_1 -> (3, 1)_(-1/3): keine G_SM-äquivariante Abbildung",
      aeq["R^2-Schnitt"]["dim"] == 3 and aeq["R^2-Schnitt"]["Y-Anteil max"] < 1e-9 and aeq["R^2-Schnitt"]["Rang des sl(2)-Anteils"] == 3
      and aeq["V4-Schnitt (z34)"]["dim"] == 1 and schnitt_v4 == {"mit dem Dublett": 2, "mit der Farbe": 1}
      and Y_quelle - Y_ziel == Fr(4, 3) and Y_e2 != Fr(4, 3))


def car(n):
    """Jordan-Wigner: Vernichter a_0..a_(n-1) auf (C^2)^n, Zustand |n_0 ... n_(n-1)>, Bit j = Besetzung von Mode j."""
    dim = 2**n
    ops = []
    for j in range(n):
        A = np.zeros((dim, dim))
        for s in range(dim):
            if (s >> j) & 1:
                A[s ^ (1 << j), s] = (-1)**bin(s & ((1 << j) - 1)).count("1")
        ops.append(A)
    return ops


a5 = car(5)
SIG = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
Jw = [sum(a5[i].T @ a5[j]*(S[p, q]/2) for p, i in enumerate((1, 3)) for q, j in enumerate((1, 3))) for S in SIG]
cas = sum(J @ J for J in Jw)
ew_cas = sorted(set(np.round(np.linalg.eigvalsh(cas).real, 9)))
anti = [a5[1].T @ a5[1].T, a5[1].T @ a5[3].T + a5[3].T @ a5[1].T, a5[3].T @ a5[3].T]
ZAHLEN["schwache_Casimir_auf_Lambda*C^5"] = [float(x) for x in ew_cas]
check("Schwachtriplett: auf dem ganzen Fockraum Lambda^* C^5 (32 Zustände, darin die 16) hat die schwache su(2) nur die Casimir-"
      "werte {0, 3/4}: kein j = 1, also kommt Sym^2 W = (1, 3)_1 in der Materie nicht vor; die symmetrischen Paare "
      "a_1^+ a_1^+, a_1^+ a_3^+ + a_3^+ a_1^+, a_3^+ a_3^+ sind null (Pauli), nur Lambda^2 W (Singulett) bleibt",
      ew_cas == [0.0, 0.75] and all(np.abs(X).max() < 1e-12 for X in anti))

# --- 1.4 SU(5) und E8 (D5 + A3 über Z4 verklebt)
gen5 = [E5[:, [i]] @ E5[[j], :] for i in range(5) for j in range(5) if i != j] + \
       [np.diag(E5[i] - E5[i + 1]) for i in range(4)]
spalten = []
for i in range(5):
    for j, k in itertools.combinations_with_replacement(range(5), 2):
        Bt = np.zeros((5, 5, 5))
        Bt[i, j, k] = Bt[i, k, j] = 1
        spalten.append(np.concatenate([(np.einsum('il,ljk->ijk', X, Bt) - np.einsum('ilk,lj->ijk', Bt, X)
                                        - np.einsum('ijl,lk->ijk', Bt, X)).ravel() for X in gen5]))
dim_su5 = nullraum(np.array(spalten).T).shape[1]
wurzeln = []
for i, j in itertools.combinations(range(8), 2):                  # D8-Wurzeln: D5 (40) + A3 = D3 (12) + (10, 6) (60)
    for si, sj in itertools.product((1, -1), repeat=2):
        w = np.zeros(8)
        w[i], w[j] = si, sj
        wurzeln.append(w)
for s in itertools.product((0.5, -0.5), repeat=8):                # (16, 4) + (16quer, 4quer): 128
    if sum(1 for x in s if x < 0) % 2 == 0:
        wurzeln.append(np.array(s))
WZ = np.array(wurzeln)
n5 = np.round((WZ[:, :5]**2).sum(1), 9)
n3 = np.round((WZ[:, 5:]**2).sum(1), 9)
klassen = {f"D5-Norm {a}, A3-Norm {b}": int(np.sum((n5 == a) & (n3 == b))) for a, b in sorted(set(zip(n5, n3)))}
gram = WZ @ WZ.T
FARBE, SCHWACH = [0, 2, 4], [1, 3]
ziele = []
for c in FARBE:
    for w1, w2 in itertools.combinations_with_replacement(SCHWACH, 2):
        t = np.zeros(5)
        t[c] += 1
        t[w1] -= 1
        t[w2] -= 1
        ziele.append(t)
einzeln = sum(1 for t in ziele for r in WZ if np.allclose(r[:5], t))
S2 = WZ[:, None, :] + WZ[None, :, :]
paare, paare_a3null, a3_normen = 0, 0, set()
for t in ziele:
    treffer = np.all(np.abs(S2[:, :, :5] - t) < 1e-9, axis=2)
    paare += int(treffer.sum())
    a3n = np.round((S2[:, :, 5:]**2).sum(2), 9)
    paare_a3null += int((treffer & (a3n < 1e-9)).sum())
    a3_normen |= set(a3n[treffer].tolist())
leptoquark = sum(1 for r in WZ if np.allclose(r[5:], 0) and sum(r[c] for c in FARBE) == 1 and sum(r[w] for w in SCHWACH) == -1
                 and np.count_nonzero(r) == 2)
ZAHLEN["E8"] = {"Wurzeln": int(len(WZ)), "Klassen (D5-Norm, A3-Norm)": klassen, "su5_aequivariante_Sym2->5": int(dim_su5),
                "einzelne Wurzeln mit D5-Teil e_c - e_w - e_w'": einzeln, "geordnete Wurzelpaare mit diesem D5-Teil": paare,
                "davon A3-neutral": paare_a3null, "A3-Normen dieser Paare": sorted(a3_normen),
                "A3-neutrale Leptoquark-Wurzeln e_c - e_w (lineare 2 -> 3)": leptoquark}
check("[negativ, E] SU(5) und E8: es gibt keine su(5)-äquivariante symmetrische Abbildung 5 x 5 -> 5 (Nullraum 0 von 75). E8 aus "
      "D5 + A3 (D8-Wurzeln und gerade Spinoren, Aufspaltung 5 + 3): 240 Wurzeln der Norm 2 in den Klassen (2, 0) 40, (0, 2) 12, "
      "(1, 1) 60, (5/4, 3/4) 128 (Kleber 5/4 + 3/4 = 2); keine Wurzel hat den D5-Teil e_c - e_w - e_w' oder e_c - 2 e_w (Normen 3, "
      f"5 > 2). Von den {paare} geordneten Wurzelpaaren mit diesem D5-Teil ist keines A3-neutral (A3-Norm {sorted(a3_normen)}): "
      f"zwei Dublett-Einheiten werden in E8 nur unter Abgabe eines A3-Gewichts (aus 4 x 4 = 6 + 10: über einen (10, 6)-"
      f"Erzeuger oder zwei Spinorerzeuger) zu einer Farbeinheit; "
      f"Kontrolle: {leptoquark} A3-neutrale Leptoquark-Wurzeln e_c - e_w realisieren die lineare Abbildung 2 -> 3",
      dim_su5 == 0 and len(WZ) == 240 and np.allclose(np.diag(gram), 2) and einzeln == 0 and paare > 0 and paare_a3null == 0
      and klassen == {"D5-Norm 0.0, A3-Norm 2.0": 12, "D5-Norm 1.0, A3-Norm 1.0": 60, "D5-Norm 1.25, A3-Norm 0.75": 128,
                      "D5-Norm 2.0, A3-Norm 0.0": 40} and leptoquark == 6)

# --- 1.5 native Paaroperatoren des Originals: zwei Fermipaare (x, -x), x = 1/4, 5/4
# Moden: 0: psi_(1/4), 1: psi_(5/4) (Quarter a, Charakter i), 2: psi_(-1/4), 3: psi_(-5/4) (Quarter b, Charakter -i)
c4 = car(4)
cd = [X.T for X in c4]
CHAR = [1j, 1j, -1j, -1j]
Rdiag = np.diag([np.prod([CHAR[j] if (s >> j) & 1 else 1 for j in range(4)]) for s in range(16)])


def adR(X):
    """Ad_R-Charakter von X, falls X ein R-Eigenoperator ist, sonst None."""
    Y = Rdiag @ X @ np.linalg.inv(Rdiag)
    for ch in (1, 1j, -1, -1j):
        if np.allclose(Y, ch*X):
            return ch
    return None


paar_sym_einzeln = [cd[0] @ cd[0], cd[0] @ cd[2] + cd[2] @ cd[0], cd[2] @ cd[2]]  # ein Fermipaar: symmetrischer Teil
T_pp = cd[0] @ cd[1] - cd[1] @ cd[0]
T_pm = cd[0] @ cd[3] - cd[1] @ cd[2]
T_mm = cd[2] @ cd[3] - cd[3] @ cd[2]
T_mp = cd[2] @ cd[1] - cd[3] @ cd[0]
trip = [T_pp, T_pm, T_mm]
adj = [cd[0] @ c4[2], cd[2] @ c4[0], cd[0] @ c4[0] - cd[2] @ c4[2]]
D = cd[0] @ c4[2]                                                 # psi*_(1/4) psi_(-1/4): eta-Umklapp am niedrigsten Paar
PAR = np.diag([(-1)**bin(s).count("1") for s in range(16)])
verdreht = [X @ D.conj().T for X in trip]
ch_trip = [adR(X) for X in trip]
ch_adj = [adR(X) for X in adj]
ch_verdreht = [adR(X) for X in verdreht]
gerade = all(np.allclose(PAR @ X @ PAR, X) for X in trip + verdreht + [D])
E_lad = lambda al, Q: Fr(Q*Q, 2) + (al - Fr(1, 2))*Q
sekt_m1 = sorted((E_lad(Fr(1, 4), Qa) + E_lad(Fr(3, 4), Qb), Qa, Qb) for Qa in range(-4, 5) for Qb in range(-4, 5)
                 if (Qa + Qb) % 2 == 0 and (Qa - Qb) % 4 == 2)
ladungen = {"Dublett": [(1, 0), (0, 1)], "Sym^2-Paare": [(2, 0), (1, 1), (0, 2)], "D": (1, -1)}
verdreht_lad = [(qa - 1, qb + 1) for qa, qb in ladungen["Sym^2-Paare"]]
konstant = [(p, q) for p in range(-3, 4) for q in range(-3, 4) if (p, q) != (0, 0)
            and len({p*qa + q*qb for qa, qb in verdreht_lad}) == 1]
verhaeltnisse = sorted({Fr(p*verdreht_lad[0][0] + q*verdreht_lad[0][1], p*1 + q*0) for p, q in konstant if p != 0})
ZAHLEN["native_paare"] = {"Ad_R Sym^2-Tripel": [str(c) for c in ch_trip], "Ad_R eta-Adjungierte": [str(c) for c in ch_adj],
                          "Ad_R verdreht (mal D^+)": [str(c) for c in ch_verdreht],
                          "niedrigster R = -1-Bosonsektor (L0, Q_a, Q_b)": [str(x) for x in sekt_m1[0]],
                          "U(1)^2-Ladungen der verdrehten Drei": verdreht_lad,
                          "Y(Drei)/Y(Dublett) der einzigen konstanten U(1)": [str(x) for x in verhaeltnisse]}
check("[E] native Paaroperatoren (CAR zweier Fermipaare (x, -x), x = 1/4, 5/4): auf EINEM Paar verschwindet der symmetrische "
      "Teil (Pauli), Sym^2 des eta-Dubletts braucht zwei Bahnstufen (innen symmetrisch, Bahn antisymmetrisch, T_+- = T_-+). Das "
      "Sym^2-Tripel hat die Ad_R-Charaktere (-1, +1, -1), die eta-Adjungierte (-1, -1, +1): als Multimenge beide {+1, -1, -1}, "
      "die Drei des Trägers {+1, +1, -1}. Die Verdrehung durch das eindeutige niedrigste R = -1-Boson D = psi*_(1/4) psi_(-1/4) "
      "(Sektor (1, -1) bei L_0 = 1/2, k = 2 mod 4) gibt (+1, -1, +1), aber alle diese Operatoren sind gerade (bosonisch)",
      all(np.abs(X).max() < 1e-12 for X in paar_sym_einzeln) and all(np.abs(X).max() > 0.5 for X in trip)
      and np.allclose(T_pm, T_mp) and ch_trip == [-1, 1, -1] and ch_adj == [-1, -1, 1] and ch_verdreht == [1, -1, 1]
      and gerade and adR(D) == -1 and (sekt_m1[0][0], sekt_m1[0][1], sekt_m1[0][2]) == (Fr(1, 2), 1, -1)
      and sekt_m1[1][0] > sekt_m1[0][0])
check("[negativ, E] Ladungen der verdrehten Drei: die U(1)^2-Ladungen (Q_a, Q_b) sind (1, 1), (0, 2), (-1, 3); konstant auf der "
      "Drei ist nur Q_a + Q_b (= 2), auf dem Dublett 1: Verhältnis Y(Drei)/Y(Dublett) = 2 statt -2/3 wie in der 16. Im Original "
      "wirkt innen nur der Torus U(1) x U(1) (z82), keine SU(3): Farbe ist nativ nicht identifizierbar, nur Z4- und U(1)-Ladungen",
      verdreht_lad == [(1, 1), (0, 2), (-1, 3)] and all(p == q for p, q in konstant) and verhaeltnisse == [Fr(2)])
befund("Teil 1: Die Konstruktion der zehnten Analyse stimmt exakt (Quotient = R^2 = +1-Eigenraum, Charaktere (+1, -1, +1), "
       "Divisor = ggT, erzwungen). Ihre Bedeutung ist aber nur Z4/D4-darstellungstheoretisch: Sym^2 eines Dubletts ist ein "
       "Schwachtriplett (1, 3)_1, die Teilung durch (uv)^2 ist mit Y unverträglich, g_SM bleibt nur als diagonale sl(2) erhalten, "
       "SU(5) hat keine solche Abbildung, und E8 erreicht sie nur unter Abgabe eines A3-Gewichts. Im Original sind die passenden "
       "Paaroperatoren bosonisch und tragen das falsche Y-Verhältnis. Die Drei ist algebraisch, nicht physikalisch", True)

# =================================================================== Teil 2
say("\nTeil 2: Verklebung Spin x_Z2 Z4 im radialen Bild")

# --- 2.1 volle Standardform (z83, CCR-Schritte): 2 pi-Drehung und Helizitätssektoren
AMP_CCR = (lambda m, k: np.sqrt(k + 1), lambda m, k: np.sqrt(m - k))
SCHRITT = (lambda k: k + 1, lambda k: k)                          # S_1 Teilchen hoch, S_2 Loch runter (z83)


def standardform(M):
    """(+)_(m, mquer <= M) H_m (x) H-quer_mquer wie z83 (voll=True): P[a,b] = S_a (x) S-quer_b, A = Bilineare links, B rechts."""
    lab = [(m, mb, k, kb) for m in range(1, M + 1) for mb in range(1, M + 1) for k in range(m) for kb in range(mb)]
    ix = {s: i for i, s in enumerate(lab)}
    n = len(lab)
    P = {}
    for a in range(2):
        for b in range(2):
            r, c, w = [], [], []
            for (m, mb, k, kb), i in ix.items():
                if m < M and mb < M:
                    r.append(ix[(m + 1, mb + 1, SCHRITT[a](k), SCHRITT[b](kb))])
                    c.append(i)
                    w.append(AMP_CCR[a](m, k)*AMP_CCR[b](mb, kb))
            P[(a, b)] = sps.csr_matrix((w, (r, c)), shape=(n, n)).toarray()

    def bilinear(c_, a_, seite):
        r, cc, w = [], [], []
        for (m, mb, k, kb), i in ix.items():
            mm, kk = (m, k) if seite == 0 else (mb, kb)
            k1 = kk - (1 if a_ == 0 else 0)
            if mm == 1 or not 0 <= k1 <= mm - 2:
                continue
            k2 = SCHRITT[c_](k1)
            r.append(ix[(m, mb, k2, kb) if seite == 0 else (m, mb, k, k2)])
            cc.append(i)
            w.append(AMP_CCR[a_](mm - 1, k1)*AMP_CCR[c_](mm - 1, k1))
        return sps.csr_matrix((w, (r, cc)), shape=(n, n)).toarray()

    A = {(c_, a_): bilinear(c_, a_, 0) for c_ in range(2) for a_ in range(2)}
    B = {(c_, a_): bilinear(c_, a_, 1) for c_ in range(2) for a_ in range(2)}
    return lab, ix, P, A, B, np.diag([0.5*(s[0] + s[1]) for s in lab])


M2 = 7
LAB, IXS, PS, AS, BS, ES = standardform(M2)
DG = {(c_, a_): AS[(c_, a_)] - BS[(a_, c_)] for c_ in range(2) for a_ in range(2)}   # räumliche Drehungen (z83)
JX, JY, JZ = (DG[(0, 1)] + DG[(1, 0)])/2, (DG[(0, 1)] - DG[(1, 0)])/2j, (DG[(0, 0)] - DG[(1, 1)])/2
achse = rng.normal(size=3)
achse /= np.linalg.norm(achse)
JN = achse[0]*JX + achse[1]*JY + achse[2]*JZ
MV, MBV = np.array([s[0] for s in LAB]), np.array([s[1] for s in LAB])
U2 = np.zeros((len(LAB), len(LAB)), complex)
f2, f4, spins_halb = 0.0, 0.0, True
CJ = JX @ JX + JY @ JY + JZ @ JZ
for m in range(1, M2 + 1):
    for mb in range(1, M2 + 1):
        sel = np.where((MV == m) & (MBV == mb))[0]
        blk = expm(2j*np.pi*JN[np.ix_(sel, sel)])
        U2[np.ix_(sel, sel)] = blk
        f2 = max(f2, np.abs(blk - (-1)**(m - mb)*np.eye(len(sel))).max())
        f4 = max(f4, np.abs(blk @ blk - np.eye(len(sel))).max())
        if (m - mb) % 2:
            jj = (-1 + np.sqrt(1 + 4*np.linalg.eigvalsh(CJ[np.ix_(sel, sel)].real)))/2
            spins_halb &= np.allclose((2*jj) % 2, 1, atol=1e-6)
hel = {}
for zwei_h in range(-3, 4):
    bloecke = sorted((float((m + mb)/2), m*mb) for m in range(1, M2 + 1) for mb in range(1, M2 + 1) if m - mb == zwei_h)
    hel[str(Fr(zwei_h, 2))] = {"U_2pi": (-1)**zwei_h, "E_min": bloecke[0][0], "Dimensionen je E": [d for _, d in bloecke[:4]]}
weyl = hel["1/2"]["Dimensionen je E"] == [n*(n + 1) for n in range(1, 5)]
ZAHLEN["helizitaetssektoren (M = 7)"] = hel
ZAHLEN["U_2pi_fehler"] = {"U_2pi - (-1)^(m - mquer)": f2, "U_4pi - 1": f4}
check("[E] volle radiale Standardform (z83, CCR, M = 7, 784 Zustände): die 2 pi-Drehung um eine zufällige Achse, U_2pi = "
      f"exp(2 pi i n.J) mit den räumlichen Drehungen A_ca - B_ac, ist auf jedem Block (m, mquer) gleich (-1)^(m - mquer) = "
      f"(-1)^(2h) = e^(2 pi i E) (Fehler {f2:.0e}), U_4pi = 1. Halbzahlige Helizität h = +-1/2, +-3/2: U_2pi = -1, halbzahlige "
      "Drehimpulse, E_min = |h| + 1; der Sektor h = 1/2 hat die Weyl-Zählung n(n+1) = 2, 6, 12, 20 bei E = n + 1/2",
      f2 < 1e-12 and f4 < 1e-12 and spins_halb and weyl and hel["1/2"]["E_min"] == 1.5 and hel["3/2"]["E_min"] == 2.5)

# --- 2.2 R auf der Leiter
R_leiter = np.exp(2j*np.pi*MV)                                    # R = e^{2 pi i L_0}, L_0 = m auf den Paarschalen (z83 Teil 5)
Rq_leiter = np.exp(2j*np.pi*MBV)
U2_diag = np.real(np.diag(U2))
halb = (MV - MBV) % 2 == 1
falsch = np.abs(U2_diag - np.real(R_leiter**2)) > 0.5
anteil_falsch = float(falsch.mean())
ZAHLEN["leiter: Anteil U_2pi != R^2"] = anteil_falsch
check("[negativ, E] auf der ganzen Leiter ist L_0 = m und L-quer_0 = mquer ganz, also R = R-quer = 1 und R^2 = +1; U_2pi = R^2 "
      f"gilt genau auf den ganzzahligen Helizitäten und scheitert auf allen halbzahligen ({anteil_falsch:.1%} der Zustände): "
      "die radialen Spinoren sind neutrale Paare, R-Bosonen. Die wörtliche Verklebung auf der Leiter selbst ist falsch",
      np.allclose(R_leiter, 1) and np.allclose(Rq_leiter, 1) and np.array_equal(falsch, halb) and 0.4 < anteil_falsch < 0.6)

# --- 2.3 native radiale Energie auf allen Ladungssektoren (Jordan-Wigner, exakte Brüche)
KMIN, KMAX = -9, 9


def bit(k):
    return 1 << (k - KMIN)


def nunter(s, k):
    return bin(s & ((1 << (k - KMIN)) - 1)).count("1")


def bil(a, b, vec, coef=1):
    """coef psi*_a psi_b auf {Bitmaske: Koeffizient} (wie z83)."""
    out = {}
    for s, cc in vec.items():
        if not (s >> (b - KMIN)) & 1:
            continue
        sg = -1 if nunter(s, b) % 2 else 1
        s1 = s ^ bit(b)
        if (s1 >> (a - KMIN)) & 1:
            continue
        if nunter(s1, a) % 2:
            sg = -sg
        s2 = s1 | bit(a)
        out[s2] = out.get(s2, 0) + coef*sg*cc
    return {k: x for k, x in out.items() if x != 0}


def erzeuge(a, vec):
    out = {}
    for s, cc in vec.items():
        if (s >> (a - KMIN)) & 1:
            continue
        out[s | bit(a)] = out.get(s | bit(a), 0) + (-1 if nunter(s, a) % 2 else 1)*cc
    return out


def vadd(p, q, f=1):
    w = dict(p)
    for k, cc in q.items():
        w[k] = w.get(k, 0) + f*cc
    return {k: cc for k, cc in w.items() if cc != 0}


def opsum(terme):
    def f(vec):
        out = {}
        for coef, a, b in terme:
            if KMIN <= a <= KMAX and KMIN <= b <= KMAX and coef != 0:
                out = vadd(out, bil(a, b, vec, coef))
        return out
    return f


LNAT = {-1: opsum([(Fr(k + 1), k + 1, k) for k in range(KMIN, KMAX + 1)]),
        1: opsum([(Fr(k), k - 1, k) for k in range(KMIN, KMAX + 1)])}


def partitionen(n, maxpart=None):
    maxpart = n if maxpart is None else maxpart
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitionen(n - k, k):
            yield (k,) + rest


def maya(Q, lam):
    s = 0
    for i in range(len(lam)):
        s |= bit(lam[i] - (i + 1) + Q)
    for y in range(KMIN, Q - len(lam)):
        s |= bit(y)
    return s


def ladung(s):
    return sum(1 for k in range(0, KMAX + 1) if (s >> (k - KMIN)) & 1) - sum(1 for k in range(KMIN, 0) if not (s >> (k - KMIN)) & 1)


def L0(s, al):
    return sum((k + al) for k in range(0, KMAX + 1) if (s >> (k - KMIN)) & 1) - \
        sum((k + al) for k in range(KMIN, 0) if not (s >> (k - KMIN)) & 1)


ZUST = [(Q, lam, maya(Q, lam)) for Q in range(-2, 3) for n in range(4) for lam in partitionen(n)]
ok_enat, ok_l0 = True, True
ENAT = {}
for Q, lam, s in ZUST:
    kom = vadd(LNAT[1](LNAT[-1]({s: 1})), LNAT[-1](LNAT[1]({s: 1})), -1)
    soll = {s: 2*sum(lam) + Q*Q} if 2*sum(lam) + Q*Q else {}
    ok_enat &= kom == soll
    ENAT[s] = Fr(2*sum(lam) + Q*Q, 2)
    for al in (Fr(1, 4), Fr(3, 4)):
        ok_l0 &= ladung(s) == Q and L0(s, al) == sum(lam) + E_lad(al, Q) and ENAT[s] == L0(s, al) - (al - Fr(1, 2))*Q
ein_teilchen = all(LNAT[-1]({maya(1, (p,) if p else ()): 1}) == {maya(1, (p + 1,)): p + 1} for p in range(0, 5))
glied = all(abs(np.exp(2j*np.pi*float(n + Fr(Q*Q, 2))) - (-1)**Q) < 1e-12 for Q in range(-4, 5) for n in range(3))
statistik = all(abs(np.exp(2j*np.pi*float(na + Fr(Qa*Qa, 2) + nb + Fr(Qb*Qb, 2)))
                    - np.exp(4j*np.pi*float(na + E_lad(Fr(1, 4), Qa) + nb + E_lad(Fr(3, 4), Qb)))) < 1e-12
                and abs(np.exp(4j*np.pi*float(na + E_lad(Fr(1, 4), Qa) + nb + E_lad(Fr(3, 4), Qb))) - (-1)**(Qa + Qb)) < 1e-12
                for Qa in range(-3, 4) for Qb in range(-3, 4) for na in range(3) for nb in range(3))
check(f"[E] native radiale sl(2,R) auf allen Ladungssektoren ({len(ZUST)} Maya-Zustände, |Q| <= 2, Stufe <= 3, Fenster 19 Moden): "
      "E^nat = [L^nat_1, L^nat_-1]/2 ist diagonal mit E^nat = L_0^osc + Q^2/2 = L_0 - (alpha - 1/2) Q für beide Quarter; das "
      "Einzelteilchen psi*_p Omega ist die Leiter D_1/2 (E^nat = p + 1/2, Amplitude p + 1, ein Zustand je Stufe, ohne "
      "Schrittart-Index). Also e^(2 pi i E^nat) = (-1)^Q, und für beide Quarter e^(2 pi i (E^nat_a + E^nat_b)) = (-1)^(Q_a + Q_b) "
      "= R^2: die radiale 2 pi-Phase des nativen Hebers ist die Statistik, aber die halbe Zahl stammt aus Q^2/2 (Ladungsgitter), "
      "nicht aus einem Spinor-K-Typ", ok_enat and ok_l0 and ein_teilchen and glied and statistik)

# --- 2.4 Vertexoperatoren und die Familie U^osc c(Q)
ok_ad, ok_enat_ad, n_ad = True, True, 0
for al in (Fr(1, 4), Fr(3, 4)):
    for Q, lam, s in ZUST:
        for r in range(-3, 4):
            w = erzeuge(r, {s: 1})
            if not w:
                continue
            t = next(iter(w))
            d_osc = (L0(t, al) - E_lad(al, Q + 1)) - (L0(s, al) - E_lad(al, Q))
            ok_ad &= d_osc == r - Q                                     # Ad_U psi*_x = (-1)^(r - Q) psi*_x, U = (-1)^(L_0^osc)
            ok_enat_ad &= (L0(t, al) - (al - Fr(1, 2))*(Q + 1)) - (L0(s, al) - (al - Fr(1, 2))*Q) == r + Fr(1, 2)
            n_ad += 1


def c_fam(kappa, Q):
    """Kovarianz: c(Q+1)/c(Q) = kappa (-1)^Q, c(0) = 1."""
    c = 1 + 0j
    if Q >= 0:
        for q in range(Q):
            c *= kappa*(-1)**q
    else:
        for q in range(-1, Q - 1, -1):
            c /= kappa*(-1)**q
    return c


QS = range(-3, 4)
fam = {}
for k8 in range(8):
    kap = np.exp(2j*np.pi*k8/8)
    cs = [c_fam(kap, Q) for Q in QS]
    fam[f"kappa = e^(2 pi i {k8}/8)"] = {
        "Involution": all(abs(c*c - 1) < 1e-12 for c in cs),
        "V = R^2 auf den Ladungsgrundzuständen": all(abs(c - (-1)**Q) < 1e-12 for c, Q in zip(cs, QS)),
        "V^2 = R^2": all(abs(c*c - (-1)**Q) < 1e-12 for c, Q in zip(cs, QS)),
        "c(Q), Q = -3..3": [complex(np.round(c, 12)) for c in cs]}
beides = [k for k, d in fam.items() if d["Involution"] and d["V = R^2 auf den Ladungsgrundzuständen"]]
invol = [k for k, d in fam.items() if d["Involution"]]
quadrat = [k for k, d in fam.items() if d["V^2 = R^2"]]
ist_enat = all(abs(c_fam(1j, Q) - 1j**(Q*Q)) < 1e-12 for Q in QS) and all(
    abs((-1)**sum(lam)*c_fam(1j, Q) - np.exp(1j*np.pi*float(ENAT[s]))) < 1e-12 for Q, lam, s in ZUST)
ZAHLEN["familie_U_osc_c(Q)"] = {"involutiv": invol, "V^2 = R^2": quadrat, "involutiv und V = R^2": beides,
                                "Kovarianz geprüft an (Zustand, Mode)": n_ad}
say(f"   Familie V = U^osc c(Q): involutiv {invol}; V^2 = R^2 {quadrat}; beides (Verklebung) {beides}")
check(f"[negativ, E] Vertexoperatoren ({n_ad} Paare aus Zustand und Mode, beide Quarter): mit der radialen 2 pi-Drehung auf dem "
      "Oszillatorfaktor, U^osc = (-1)^(L_0^osc) (z83 R-neutral, universell nach z82), gilt Ad_U psi*_x = (-1)^(r - Q) psi*_x, "
      "das Vorzeichen hängt vom Ladungssektor ab: psi*_x ist nicht kovariant, und auf dem tiefsten Fermion (Ladung 1, "
      "Oszillatorvakuum) ist U = +1, R^2 = -1. Jede Fortsetzung V = U^osc c(Q) ist genau dann kovariant, wenn c(Q+1)/c(Q) = "
      "kappa (-1)^Q, also c(Q) = kappa^Q (-1)^(Q(Q-1)/2) (der Faktor ist der konforme Spin Q^2/2 des Ladungsvertex). Über alle "
      "kappa (8. Einheitswurzeln): involutiv nur kappa = +-1, und dann verfehlt V das R^2 = (-1)^Q bei Q = 1 oder Q = 2; "
      "V^2 = R^2 nur für kappa = +-i, und kappa = i ist genau e^(i pi E^nat) (Ad = i(-1)^r, kovariant), aber kein Element der "
      "Ordnung 2. Keine Fortsetzung ist zugleich Involution, kovariant und gleich R^2",
      ok_ad and ok_enat_ad and not beides and sorted(invol) == ["kappa = e^(2 pi i 0/8)", "kappa = e^(2 pi i 4/8)"]
      and sorted(quadrat) == ["kappa = e^(2 pi i 2/8)", "kappa = e^(2 pi i 6/8)"] and ist_enat)
befund("Teil 2: Eine einzige Lorentzwirkung, deren -1 gleich R^2 ist und die die Vertexoperatoren trägt, gibt es in der nativen "
       "Familie nicht. Genau scheitert dreierlei: (a) die halbzahligen Helizitäten der radialen Standardform sind neutrale Paare "
       "mit R = 1; (b) auf dem Oszillatorfaktor ist die 2 pi-Drehung R-blind, die Vertexoperatoren bekommen ein vom "
       "Ladungssektor abhängiges Vorzeichen; (c) die einzige kovariante Korrektur ist der quadratische Ladungsspin Q^2/2: sie "
       "macht aus der 2 pi-Drehung ein Element mit Quadrat R^2 (eine halbe Drehung im Spinorsinn), nicht eines mit 2 pi = R^2. "
       "Positiv und neu: e^(2 pi i E^nat) = R^2 = (-1)^F exakt; die Statistik sitzt in der radialen Energie des Ladungsgitters", True)

# =================================================================== Teil 3
say("\nTeil 3: Elfte Analyse -- Verzweigungsfelder der Z2-Permutationsorbifaltung an den Marken")

# --- 3.1 Gewicht h = c/16 und Ising-Koset
cS, NS = sp.symbols("c N", positive=True)
h_tw = cS/24*(NS - 1/NS)                                          # zyklische Orbifaltung Z_N: Grundenergie des Twistsektors
c_E8 = lambda k: sp.Rational(248*k, k + 30)                       # c(G_k) = k dim G/(k + h^v), h^v(E8) = 30
kac = lambda r, s_, m: sp.Rational((r*(m + 1) - s_*m)**2 - 1, 4*m*(m + 1))
h8 = h_tw.subs({NS: 2, cS: 8})
ZAHLEN["verzweigungsfelder"] = {"h_twist(N = 2)": str(sp.simplify(h_tw.subs(NS, 2))), "c = 8": str(h8),
                                "c = 1 (ein Quarter)": str(h_tw.subs({NS: 2, cS: 1})),
                                "c(E8_1 x E8_1) - c(E8_2)": str(2*c_E8(1) - c_E8(2)),
                                "Ising h_(2,1), h_(1,2)": [str(kac(2, 1, 3)), str(kac(1, 2, 3))],
                                "tiefstes Fermion des Originals L_0": str(E_lad(Fr(1, 4), 1))}
check("[E] Gewicht: die Z2-Permutationsorbifaltung hat h_twist = (c/24)(N - 1/N) = c/16, mit c = 8 (Naht (E8)_1) also h = 1/2; "
      "Koset (E8)_1 x (E8)_1 = (E8)_2 x Ising (16 = 31/2 + 1/2), und h = 1/2 ist das Gewicht des Ising-Fermions h_(2,1) "
      "(gleiche Chiralitäten). Zum Vergleich: ein Quarter (c = 1) gäbe 1/16 = h_(1,2), das tiefste Fermion des Originals hat "
      "L_0 = 1/4 (Beobachtung, keine Zuordnung)",
      sp.simplify(h_tw.subs(NS, 2) - cS/16) == 0 and h8 == sp.Rational(1, 2) and c_E8(1) == 8 and c_E8(2) == sp.Rational(31, 2)
      and 2*c_E8(1) - c_E8(2) == sp.Rational(1, 2) and kac(2, 1, 3) == h8 and kac(1, 2, 3) == sp.Rational(1, 16))

# --- 3.2 Vierteldrehung und Uhr auf den Verzweigungsfeldern
C4 = np.roll(np.eye(4), 1, axis=0)                                # Marke j -> j + 1 (z -> i z permutiert mu_4 zyklisch)
vt = {}
for nm, spin in (("chiral (1/2, 0)", 0.5), ("gleiche Chiralitäten (1/2, 1/2)", 0.0)):
    Rt = np.exp(1j*np.pi/2*spin)*C4                               # sigma_j -> (dw/dz)^h (dw/dz-quer)^h-quer sigma_(j+1), dw/dz = i
    F_ = np.exp(2j*np.pi*spin)                                    # (-1)^F = e^{2 pi i s}
    ew2 = np.linalg.eigvals(Rt @ Rt)
    R4 = np.linalg.matrix_power(Rt, 4)
    vt[nm] = {"R^4": complex(np.round(R4[0, 0], 12)), "R^4 skalar": bool(np.allclose(R4, R4[0, 0]*np.eye(4))),
              "Eigenwerte R^2": sorted({str(complex(np.round(e, 9))) for e in ew2}),
              "(-1)^F": complex(np.round(F_, 12)), "R^2 = (-1)^F auf": int(np.sum(np.abs(ew2 - F_) < 1e-9)),
              "R^4 = (-1)^F": bool(np.allclose(R4, F_*np.eye(4)))}
uhr_halb = complex(np.round(np.exp(2j*np.pi*0.5), 12))           # Uhr R = e^{2 pi i L_0} auf einem Zustand mit L_0 = 1/2
ZAHLEN["verzweigungsfelder_R"] = vt
ZAHLEN["Uhr auf L_0 = 1/2"] = {"R": uhr_halb, "R^2": complex(np.round(uhr_halb**2, 12))}
ch, gl_ = vt["chiral (1/2, 0)"], vt["gleiche Chiralitäten (1/2, 1/2)"]
check("[negativ, E] die Vierteldrehung permutiert die vier Marken zyklisch und dreht jede lokal um pi/2: R sigma_j = i^s sigma_(j+1). "
      "Chiral (h, h-quer) = (1/2, 0): R^4 = -1, R hat Ordnung 8, R^2 hat die Eigenwerte +-i; R^2 = (-1)^F gilt auf 0 von 4 "
      "Fouriermoden, stattdessen R^4 = (-1)^F (Spin x_Z2 Z8, nicht Z4), und die Uhr des Originals hat R^4 = 1. Gleiche "
      "Chiralitäten (1/2, 1/2): Spin 0, R^4 = 1, R^2 = (+1, -1, +1, -1) auf den Fouriermoden, (-1)^F = +1: Treffer 2 von 4. Die "
      "Uhr R = e^(2 pi i L_0) gibt einem Zustand mit L_0 = 1/2 den Wert -1 und R^2 = +1: den Charakter des Divisors (uv)^2 und "
      "des Umklapp-Bosons D aus Teil 1, also R^2-gerade",
      ch["R^4 skalar"] and abs(ch["R^4"] + 1) < 1e-12 and ch["R^2 = (-1)^F auf"] == 0 and ch["R^4 = (-1)^F"]
      and gl_["R^4 skalar"] and abs(gl_["R^4"] - 1) < 1e-12 and gl_["R^2 = (-1)^F auf"] == 2
      and abs(uhr_halb + 1) < 1e-12 and abs(uhr_halb**2 - 1) < 1e-12)

# --- 3.3 radiales Bild: Blätterwechsel und Helizität
perm = np.array([IXS[(s[1], s[0], s[3], s[2])] for s in LAB])
FL = np.zeros((len(LAB), len(LAB)))
FL[perm, np.arange(len(LAB))] = 1
HH = np.diag([0.5*(s[0] - s[1]) for s in LAB])
flip = {"F E F - E": float(np.abs(FL @ ES @ FL - ES).max()), "F h F + h": float(np.abs(FL @ HH @ FL + HH).max()),
        "F P_ab F - P_ba": float(max(np.abs(FL @ PS[(a, b)] @ FL - PS[(b, a)]).max() for a in range(2) for b in range(2))),
        "F J F - (-Jx, Jy, -Jz)": float(max(np.abs(FL @ JX @ FL + JX).max(), np.abs(FL @ JY @ FL - JY).max(),
                                            np.abs(FL @ JZ @ FL + JZ).max())),
        "F U_2pi F - U_2pi": float(np.abs(FL @ U2 @ FL - U2).max())}
spek = {(Fr(m + mb, 2), Fr(m - mb, 2)) for m in range(1, M2 + 1) for mb in range(1, M2 + 1)}
spek_ok = all((2*h).denominator == 1 and E - abs(h) - 1 >= 0 and (E - abs(h) - 1).denominator == 1 for E, h in spek)
tw = {"(1/2, 0)": (Fr(1, 4), Fr(1, 4)), "(1/2, 1/2)": (Fr(1, 2), Fr(0))}
tw_raus = all(p not in spek and (p[0] < abs(p[1]) + 1) for p in tw.values())
ZAHLEN["blaetterwechsel"] = flip
ZAHLEN["verzweigungsfeld_radial (E, h)"] = {k: [str(x) for x in p] for k, p in tw.items()}
check("[negativ, E] radiales Bild mit den Blättern = Original und Tomita-Partner (J Delta J = Delta^-1): der Blätterwechsel F "
      "(H_m (x) H-quer_mquer -> H_mquer (x) H-quer_m) erhält E, kehrt h um, vertauscht P_ab -> P_ba und wirkt auf die Drehungen "
      "als (-Jx, Jy, -Jz): er ist die räumliche Spiegelung y -> -y (orientierungsumkehrend, Helizitätsumkehr) und vertauscht "
      "mit U_2pi. Ein Verzweigungsfeld (L_0, L-quer_0) = (1/2, 0) hätte (E, h) = (1/4, 1/4), eines mit (1/2, 1/2) (E, h) = "
      "(1/2, 0); beide liegen außerhalb des Spektrums der Standardform (2h ganz, E - |h| - 1 ganz >= 0): keine halbzahlige "
      "Helizität, kein masseloser K-Typ",
      max(flip.values()) < 1e-12 and spek_ok and tw_raus)
befund("Teil 3: Die Gewichtsrechnung der elften Analyse stimmt (c/16 = 1/2 bei c = 8, Ising-Fermion des Kosets). Mit der "
       "Verklebungsfrage ist sie nicht verträglich: An den Marken gibt die Vierteldrehung R^4 = (-1)^F (Z8) statt R^2 = (-1)^F, "
       "die Uhr gibt einem Gewicht 1/2 den Charakter -1 (R^2 = +1, wie der Divisor (uv)^2), im radialen Bild hat ein "
       "Verzweigungsfeld h = 1/4 oder 0, und der Blätterwechsel zwischen den Tomita-konjugierten Blättern ist eine räumliche "
       "Spiegelung (Helizitätsumkehr), keine 2 pi-Drehung. Die Warnung der Analyse trifft: hier ist der Blätterwechsel gerade "
       "die orientierungsumkehrende Spiegelung", True)

# =================================================================== Teil 4
say("\nTeil 4: Bilanz für F4")
STATUS = {
    "Drei = Sym^2(W)/(uv)^2": ("E", "exakt; Divisor = ggT (erzwungen); Charaktere (+1, -1, +1); D4-äquivariant bis auf chi((uv)^2)"),
    "Drei als Farbtriplett mit Y = -1/3": ("E (negativ, gescopt)", "Sym^2 W = (1, 3)_1; g_SM-Äquivarianz dim 3 (R^2) bzw. 1 (V4) von "
                                                                    "12, ohne Y; Divisor bräuchte Y = 4/3; SU(5): 0; E8 nur mit "
                                                                    "A3-Gewicht; kein j = 1 in Lambda^* C^5"),
    "native Paaroperatoren (k = 2 mod 4)": ("E (negativ)", "Charaktere erst nach Verdrehung durch das Umklapp-Boson D (L_0 = 1/2); "
                                                           "gerade, bosonisch; Y-Verhältnis 2 statt -2/3"),
    "2 pi-Drehung der radialen Standardform": ("E", "(-1)^(m - mquer) = e^(2 pi i E); halbzahlige Helizität -1; Weyl-Zählung"),
    "radiale Spinoren und R": ("E (negativ)", "R = 1 auf der ganzen Leiter; U_2pi = -1 != R^2 = +1 auf allen halbzahligen "
                                              "Helizitätszuständen"),
    "native radiale Energie": ("E", "E^nat = L_0^osc + Q^2/2; e^(2 pi i E^nat) = (-1)^F = R^2 auf allen Ladungssektoren"),
    "Verklebung Spin x_Z2 Z4": ("E (negativ, in der Familie U^osc c(Q))", "kovariant heißt c = kappa^Q (-1)^(Q(Q-1)/2); involutiv "
                                                                         "(kappa = +-1) verfehlt R^2, V^2 = R^2 (kappa = +-i, "
                                                                         "e^(i pi E^nat)) ist keine Involution"),
    "elfte Analyse (Verzweigungsfelder)": ("E (negativ, gescopt)", "h = c/16 = 1/2 richtig; R^4 = (-1)^F bzw. Uhr R = -1; radial "
                                                                   "h = 1/4 oder 0; Blätterwechsel = Spiegelung"),
    "Verklebung mit neuer Struktur": ("O", "ein Spinorindex auf dem Ladungsfaktor oder eine Regel (-1)^(2h) = R^2 für Ladung x "
                                           "Paar; nativ nicht vorhanden (eta-SU(2) im KMS-Zustand nicht implementierbar, z82)"),
}
for k, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k}: {t_}")
ZAHLEN["status"] = STATUS
ZAHLEN["F4_urteil_kurz"] = {"Drei": "algebraisch (Z4/D4), nicht physikalisch", "Verklebung": "nein (nativ); bedingt nur mit "
                            "neuer Struktur", "elfte Analyse": "Gewicht richtig, Spin-Z4 und Helizität nicht"}
urteil = ("F4 ergibt zwei gescopte Negative und einen positiven Baustein. Die fehlende Drei: Die Konstruktion der zehnten Analyse "
          "stimmt exakt. Sym^2(W) mit W = span(u^3 v, u v^3) hat den ggT (uv)^2, den Fixpunktdivisor der Vierteldrehung im "
          "Quadrat; der Quotient ist span(u^4, u^2 v^2, v^4) = {e_4, e_2, e_0}, genau der R^2 = +1-Eigenraum, mit den Charakteren "
          "(+1, -1, +1). Der Divisor ist erzwungen, und die Abbildung ist D4-äquivariant bis auf den Charakter von (uv)^2. Ihre "
          "Bedeutung ist aber nur Z4/D4-darstellungstheoretisch. Sym^2 des Dubletts ist das Schwachtriplett (1, 3)_1, das in "
          "Lambda^* C^5 gar nicht vorkommt (Pauli). Von g_SM bleibt im R^2-Schnitt nur die diagonale sl(2) mit dim 3 von 12 "
          "äquivariant, im V4-Schnitt dim 1, ohne Hyperladung. Die Teilung bräuchte einen Divisor mit Y = 4/3, (uv)^2 = e_2 hat "
          f"-1/3. SU(5) hat keine solche Abbildung; in E8 aus D5 + A3 realisieren {paare} Wurzelpaare den D5-Schritt, aber keines "
          "A3-neutral. Im Original tragen nur gerade, bosonische Paaroperatoren die Charaktere der Drei, und zwar erst nach "
          "Verdrehung durch das niedrigste R = -1-Boson; ihr einziges konstantes U(1) gibt das Verhältnis 2 statt -2/3. Die "
          "Drei ist also algebraisch, nicht die physikalische Farbe. Die Verklebung Spin x_Z2 Z4: Auf der radialen Standardform "
          "ist die 2 pi-Drehung (-1)^(2h) = e^(2 pi i E), mit Weyl-Zählung im Sektor h = 1/2. Diese Spinoren sind aber neutrale "
          f"Paare mit R = 1, sodass U_2pi = R^2 auf allen halbzahligen Helizitäten scheitert ({anteil_falsch:.0%} der Zustände). "
          "Neu und positiv: "
          "Die native radiale Energie ist E^nat = L_0^osc + Q^2/2 auf allen Ladungssektoren, also e^(2 pi i E^nat) = (-1)^F = R^2 "
          "exakt. Die Statistik sitzt damit in der radialen Energie des Ladungsgitters, nicht in einem Spinor-K-Typ (das "
          "Einzelteilchen ist eine D_1/2-Leiter ohne Schrittart-Index). Eine Lorentzwirkung, deren -1 gleich R^2 ist und die "
          "die Vertexoperatoren trägt, gibt es in der nativen Familie U^osc c(Q) nicht. Kovarianz erzwingt c(Q) = kappa^Q "
          "(-1)^(Q(Q-1)/2), den konformen Spin Q^2/2 des Ladungsvertex. Die Involutionen (kappa = +-1) verfehlen R^2 bei Q = 1 "
          "oder 2; mit V^2 = R^2 (kappa = +-i, e^(i pi E^nat)) ist V keine Involution, also eine halbe statt einer ganzen Drehung. "
          "Urteil Verklebung: nein (nativ), bedingt nur mit neuer Struktur, etwa einem Spinorindex auf dem Ladungsfaktor oder "
          "einer Regel (-1)^(2h) = R^2 für Ladung x Paar. Die elfte Analyse: h = c/16 = 1/2 bei c = 8 stimmt (Ising-Fermion des "
          "Kosets). An den Marken gibt die Vierteldrehung aber R^4 = (-1)^F (Z8) statt R^2, die Uhr gibt dem Gewicht 1/2 den "
          "Charakter -1 (R^2 = +1), radial hat ein Verzweigungsfeld h = 1/4 oder 0, und der Blätterwechsel zwischen den "
          "Tomita-Blättern ist die räumliche Spiegelung. Mit der Verklebung ist sie nicht verträglich. Nichts ist befördert; "
          "3+1-Aussagen bleiben Kinematik der Einpaarschicht.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z92_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
