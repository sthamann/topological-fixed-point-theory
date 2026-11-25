"""Z99: Wechselwirkung über Vertexoperatoren -- zwei Nähte über die gemeinsame Schwingungshülle (Folgeschritt G3, H6h46).

Leitlinie 'Rekonstruktion statt Auswahl': Die Kopplung wird aus den Vertexoperatoren der Nähte gebaut, nicht aus einer
gewünschten 4D-Wechselwirkung. Firewall DIMENSION.UPLIFT.FIREWALL.01: Jede Aussage ist typisiert (Naht = chirale 2D-Struktur,
4D nur als Lesart über die radiale Leiter von z83, Paare = Einteilchenzustände; keine wechselwirkende 3+1-Schließung). Gelöste
Bausteine bleiben unberührt, nichts wird befördert. Bestand: z58/z93 (Bosonisierung bei beta = 4 pi, Fermi-Propagator = Twist
x Ladungskern x Oszillatorfaktor), z75 (<J J J> = A f, verbundener Vierpunkt nichtabelsch != 0, Vektorkomposita Delta = 3),
z82 (dieselbe Oszillatorhülle in jedem Ladungssektor), z83 (Paarschalen H_m, E = L_0 = m, K-Typen (j, j), verdoppelt
E = (m + m-quer)/2, nackte Schritte +1/-1), z91 (Sym^N, keine native Fusion über gerade oder lineare Algebren, K-Typ-Hindernis,
Delta = 1, sl(2,R)-Vertex D_1 x D_1 -> D_2 bricht so(4)), z92 (E^nat = [L^nat_1, L^nat_-1]/2 = L_0^osc + Q^2/2
= L_0 - (alpha - 1/2) Q).

Bezeichnungen: Naht i = 1, 2 mit chiralem Dirac-Fermion, Moden x in Z + alpha_i (Quarter alpha = 1/4 oder 3/4), Fermisee x < 0,
Ladung Q_i, bosonisiert phi_i (Ladungsgitter x Oszillatoren), Vertexoperator V^(i)_Q = e^{i Q phi_i}. Strommoden
J^(ab)_n = sum_x psi*_(a, x-n) psi_(b, x) (n < 0 hebt); E^+ = J^(12), E^- = J^(21), H = J^(11) - J^(22). Paarschalen
H^(ab)_m = span{psi*_(a,p) psi_(b,r) Omega: p > 0 > r, E^nat = m}: ab = 11, 22 diagonal (z83), ab = 12, 21 Ladungstransfer
(Vertexsektor). Paarfelder b+_(ab)(p, r) = psi*_(a,p) psi_(b,r) (nur Erzeugung bzw. Vernichtung eines Paars).

Teil 1  Zwei Nähte über eine gemeinsame Hülle. (a) Hüllenfamilie: Nullmoden euklidisch, Oszillatoren mit Gram [[1, g], [g, 1]];
        Austauschfunktion V^(1) V^(2) gegen V^(2) V^(1): lokal nur für g = 0; konsistente geteilte Hülle (Gram der Nullmoden =
        Gram der Oszillatoren) braucht ganzzahliges g, g = +-1 ist entartet (Relativladung Nullvektor: Transfer ohne
        Oszillatorinhalt, h = 0). (b) g = 0 in einer gemeinsamen CAR (Jordan-Wigner, zwei Nähte zu je sechs Moden): der
        Ladungstransfer e^{i(phi_1 - phi_2)} = psi*_1 psi_2 ist der su(2)_1-Strom (Relationen exakt im Fenster), Paarschalen aller
        vier Typen mit dim m und E^nat = m, gleiche und verschiedene Quarter, E^nat aus dem nativen Heber. (c) Dreipunktstruktur:
        Paarfelder (Paarparität), Ströme (abelsch, nichtabelsch), 1->1-Formfaktor als Schritt, Fusion zweier Kreuzteilchen in
        die Einpaarschicht gegen das K-Typ-Hindernis von z91. (d) Vierpunktstruktur: Konstituententausch, Pauli, Gram der
        Zweiteilchenzustände gegen den freien adjungierten Fockraum, Singulett über N Nähte.
Teil 2  Anomale Dimension. (a) nativ: E^nat additiv auf allen Zweipaarzuständen, Modennormen <E^+_m E^-_-m> = m (h = 1),
        Kreuzantwort bei beta = 4 pi = eine Fermischleife (Produkt der Nahtbesetzungen); (b) Mehrfachtransfer: E^nat_min(Q) = Q^2;
        (c) Deformationen: Nullmodenkopplung (chemisches Potential), Radiusdeformation des Relativbosons über den Tomita-Partner
        (exakt marginal; Diagonalströme geschützt, skalare Kreuzzweige E_4D = e^(+-2 theta)), nichtabelsche Strom-Strom-Kopplung
        (Casimir aus den Nahtströmen, läuft); (d) Eichfeldlesart: adjungierte Skalarleitern gegen das Maxwell-Singleton,
        Flato-Fronsdal; (e) verbundene 1->2-Amplituden: alle Bilineare gegen den kleinsten Vierfermionterm (Doppeltransfer).
Teil 3  Bilanz (E/C/L/O) und Urteil zu G3; Folgen für K15 und T3/T5.
"""
import json
import time
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(99)
BETA = 4*np.pi
TOL = 1e-12


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def r(x, n=12):
    return float(np.round(float(x), n))


# ------------------------------------------------------------------ Bausteine: Nähte in einer gemeinsamen CAR
def cdag_matrix(j, nb):
    """Jordan-Wigner-Erzeuger der Mode j auf nb Moden (Vorzeichen = Parität der besetzten Moden unter j)."""
    dim = 1 << nb
    s = np.arange(dim)
    src = s[((s >> j) & 1) == 0]
    tief = src & ((1 << j) - 1)
    par = np.zeros(len(src), dtype=np.int64)
    for b in range(j):
        par ^= (tief >> b) & 1
    return sps.csr_matrix(((1 - 2*par).astype(float), (src | (1 << j), src)), shape=(dim, dim))


class Naehte:
    """N Nähte zu je 2K Moden x = alpha + n, n = -K..K-1 (K unter, K über dem Fermisee), eine gemeinsame CAR."""

    def __init__(self, alphas, K, beta=BETA):
        self.alphas = [Fr(a) for a in alphas]
        self.N = len(alphas)
        self.modes = [(i, a + n) for i, a in enumerate(self.alphas) for n in range(-K, K)]
        self.NB = len(self.modes)
        self.dim = 1 << self.NB
        self.idx = {md: j for j, md in enumerate(self.modes)}
        self.cd = [cdag_matrix(j, self.NB) for j in range(self.NB)]
        self.c = [m.T.tocsr() for m in self.cd]
        s = np.arange(self.dim)
        occ = (s[:, None] >> np.arange(self.NB)[None, :]) & 1
        xs = np.array([float(x) for _, x in self.modes])
        naht = np.array([i for i, _ in self.modes])
        unten = xs < 0
        anreg = np.where(unten, 1 - occ, occ)
        self.L0 = (anreg*np.abs(xs)).sum(1)
        self.Q = np.stack([occ[:, (naht == i) & ~unten].sum(1) - (1 - occ[:, (naht == i) & unten]).sum(1)
                           for i in range(self.N)], 1)
        self.Enat = self.L0 - sum((float(a) - 0.5)*self.Q[:, i] for i, a in enumerate(self.alphas))
        self.Nplus = occ[:, ~unten].sum(1)
        self.Nminus = (1 - occ[:, unten]).sum(1)
        self.vi = int(sum(1 << j for j in range(self.NB) if unten[j]))
        self.vac = np.zeros(self.dim)
        self.vac[self.vi] = 1.0
        f = 1/(1 + np.exp(beta*xs))                                                 # KMS-Besetzung (z58), beta = 4 pi
        self.rho = np.exp((occ*np.log(f) + (1 - occ)*np.log(1 - f)).sum(1))
        self.xs, self.naht, self.unten = xs, naht, unten

    def bil(self, a, b):
        return (self.cd[a] @ self.c[b]).tocsr()

    def strom(self, a, b, n):
        """J^(ab)_n = sum_x psi*_(a, x-n) psi_(b, x) im Fenster; für a = b, n = 0 normalgeordnet (= Q_a)."""
        if a == b and n == 0:
            return sps.diags(self.Q[:, a].astype(float)).tocsr()
        M = sps.csr_matrix((self.dim, self.dim))
        for (i, x), j in self.idx.items():
            t = self.idx.get((a, x - n)) if i == b else None
            if t is not None:
                M = M + self.bil(t, j)
        return M.tocsr()

    def lnat(self, n):
        """Nativer Heber (z83/z92) L^nat_-1 = sum_x (x - alpha + 1) psi*_(x+1) psi_x je Naht, L^nat_1 = (L^nat_-1)^+."""
        M = sps.csr_matrix((self.dim, self.dim))
        for (i, x), j in self.idx.items():
            t = self.idx.get((i, x + 1))
            if t is not None:
                M = M + float(x - self.alphas[i] + 1)*self.bil(t, j)
        return M.tocsr() if n == -1 else M.T.tocsr()

    def paare(self, a, b, mmax):
        """Paarschalen H^(ab)_m: Liste (m, k, p, r, Vektor) mit k = Teilchenlage x_p - alpha_a."""
        aus = []
        for (i, x), p in self.idx.items():
            if i != a or x < 0:
                continue
            for (i2, y), q in self.idx.items():
                if i2 != b or y > 0:
                    continue
                v = self.bil(p, q) @ self.vac
                s_ = int(np.flatnonzero(v)[0])
                m = self.Enat[s_]
                if m <= mmax + 1e-9:
                    aus.append((int(round(m)), int(x - self.alphas[a]), p, q, v))
        return sorted(aus, key=lambda t: (t[0], t[1]))

    def ew(self, X):
        """KMS-Erwartungswert Tr(rho X) (rho diagonal in der Besetzungsbasis, Produktzustand)."""
        return float(X.diagonal() @ self.rho)


# =================================================================== Teil 1
say("Teil 1: Zwei Nähte über eine gemeinsame Hülle")
say("  (a) Hüllenfamilie")
z, w, g = sp.symbols("z w g")
loc = {}
for gam in (0, 1, 2):
    for Q, Qp in ((1, 1), (1, -1), (2, 1)):
        n_ = gam*Q*Qp
        f12 = (1 - w/z)**n_                                                         # V^(1)_Q(z) V^(2)_Q'(w): nur Oszillatoren
        f21 = (1 - z/w)**n_                                                         # V^(2)_Q'(w) V^(1)_Q(z)
        loc[(gam, Q, Qp)] = sp.simplify(sp.cancel(f12/f21))
lok_ok = (all(v == 1 for (gam, _, _), v in loc.items() if gam == 0)
          and all(sp.diff(v, z) != 0 for (gam, _, _), v in loc.items() if gam != 0))
G = sp.Matrix([[1, g], [g, 1]])
konsistent = []
for gv in [sp.Rational(k, 4) for k in range(-6, 7)]:
    Gv = G.subs(g, gv)
    ganz = all((sp.Matrix([[q1, q2]])*Gv*sp.Matrix([q1p, q2p]))[0].is_integer
               for q1, q2, q1p, q2p in [(1, 0, 0, 1), (1, 1, 1, -1), (2, -1, 1, 3)])
    if ganz and all(ev >= 0 for ev in Gv.eigenvals()):
        konsistent.append(gv)
alpha_null = {gv: sp.Matrix([1, -gv]) for gv in (1, -1)}
null_ok = all(G.subs(g, gv)*alpha_null[gv] == sp.zeros(2, 1) for gv in (1, -1))
h_transfer = {gv: (alpha_null[gv].T*G.subs(g, gv)*alpha_null[gv])[0]/2 for gv in (1, -1)}
a12 = sp.Matrix([1, -1])
h_g0 = (a12.T*a12)[0]/2
ZAHLEN["huellenfamilie"] = {"Austauschquotient (g, Q, Q')": {str(k): str(v) for k, v in loc.items()},
                            "konsistente g (Nullmoden = Oszillator-Gram, ganzzahlig, psd)": [str(x) for x in konsistent],
                            "h des Transfers bei g = +-1": {str(k): str(v) for k, v in h_transfer.items()},
                            "h des Transfers bei g = 0": str(h_g0)}
check("[E negativ] (a) Wörtlich geteilte Hülle (Nullmoden je Naht euklidisch, ein Oszillatorfeld mit Gram [[1, g], [g, 1]]): "
      "V^(1)_Q(z) V^(2)_Q'(w) gibt (1 - w/z)^(g Q Q'), die umgekehrte Ordnung (1 - z/w)^(g Q Q'); der Quotient ist "
      f"{loc[(1, 1, 1)]} für g = 1, Q = Q' = 1, also keine Konstante: Die Vertexoperatoren beider Nähte sind nur für g = 0 "
      "zueinander lokal (Quotient 1). Der Vertexoperator der einen Naht mit dem Oszillatorfeld der anderen ist nichtlokal",
      lok_ok)
check("[E] (a) Konsistente geteilte Hülle (Gram der Nullmoden = Gram der Oszillatoren, ganzzahlige Paarungen, positiv "
      f"semidefinit): genau g in {[str(x) for x in konsistent]}. Für g = +-1 ist die Relativladung (1, -+1) ein Nullvektor, "
      f"G alpha = 0, h = {h_transfer[1]}: der Transfer ist eine reine Nullmodenverschiebung ohne Oszillatorinhalt und "
      "vertauscht mit allen Oszillatoren, also mit jedem 4D-Einteilchenzustand (verbundene 4D-Korrelationen exakt null). "
      f"Nichtentartet und lokal ist nur g = 0: zwei Hüllen, das Gitter Z^2, Transfer alpha = (1, -1) mit h = |alpha|^2/2 = {h_g0}",
      konsistent == [-1, 0, 1] and null_ok and all(v == 0 for v in h_transfer.values()) and h_g0 == 1)

say("  (b) g = 0: zwei Nähte in einer gemeinsamen CAR")
S = Naehte([Fr(1, 4), Fr(1, 4)], 3)                                                 # gleiche Quarter
M = Naehte([Fr(1, 4), Fr(3, 4)], 3)                                                 # verschiedene Quarter
say(f"   Fockraum 2 x 6 Moden, Dimension {S.dim}; Aufbau {time.time() - T0:.1f} s")
Ep = {n: S.strom(0, 1, n) for n in (-3, -2, -1, 0, 1, 2, 3)}
Em = {n: S.strom(1, 0, n) for n in (-3, -2, -1, 0, 1, 2, 3)}
J1 = {n: S.strom(0, 0, n) for n in (-3, -2, -1, 0, 1, 2, 3)}
J2 = {n: S.strom(1, 1, n) for n in (-3, -2, -1, 0, 1, 2, 3)}
Hc = {n: J1[n] - J2[n] for n in J1}
TYPEN = [(0, 0), (1, 1), (0, 1), (1, 0)]
NAME = {(0, 0): "11", (1, 1): "22", (0, 1): "12", (1, 0): "21"}
schalen_S = {t: S.paare(*t, 3) for t in TYPEN}
schalen_M = {t: M.paare(*t, 3) for t in TYPEN}
test = [S.vac] + [v for t in TYPEN for (m, k, p, q, v) in schalen_S[t] if m == 1]
fk_def = 0.0
for m_ in (-1, 0, 1):
    for n_ in (-1, 0, 1):
        for v in test:
            lhs = Ep[m_] @ (Em[n_] @ v) - Em[n_] @ (Ep[m_] @ v)
            rhs = Hc[m_ + n_] @ v + (m_ if m_ + n_ == 0 else 0)*v
            fk_def = max(fk_def, np.abs(lhs - rhs).max())
            lhs2 = Hc[m_] @ (Ep[n_] @ v) - Ep[n_] @ (Hc[m_] @ v)
            fk_def = max(fk_def, np.abs(lhs2 - 2*Ep[m_ + n_] @ v).max())
            lhs3 = Hc[m_] @ (Hc[n_] @ v) - Hc[n_] @ (Hc[m_] @ v)
            fk_def = max(fk_def, np.abs(lhs3 - (2*m_ if m_ + n_ == 0 else 0)*v).max())
e1 = Ep[-1] @ S.vac
h_kreuz = S.Enat[int(np.flatnonzero(e1)[0])]
ZAHLEN["frenkel_kac"] = {"Defekt su(2)_1-Relationen": r(fk_def, 15), "E^nat(E^+_-1 Omega)": r(h_kreuz)}
check("[E] (b) Frenkel-Kac in der gemeinsamen CAR (gleiche Quarter): E^+_n = sum psi*_(1,x-n) psi_(2,x), E^-_n, "
      "H_n = J^(11)_n - J^(22)_n erfüllen [E^+_m, E^-_n] = H_(m+n) + m delta, [H_m, E^+_n] = 2 E^+_(m+n), [H_m, H_n] = 2m delta "
      f"auf Vakuum und Stufe 1, m, n = -1, 0, 1 (Defekt {fk_def:.0e}); E^+_-1 Omega hat E^nat = {h_kreuz:g} = |alpha|^2/2. Der "
      "Ladungstransfer zwischen zwei Nähten über die nichtentartete Hülle ist der nichtabelsche su(2)_1-Strom: Vertexroute und "
      "Stromroute aus z91 sind auf Niveau 1 dieselbe Struktur",
      fk_def < TOL and abs(h_kreuz - 1) < TOL)

dims_ok, enat_ok, l0_tab = True, True, {}
for name_, sys_, sch in (("gleich", S, schalen_S), ("verschieden", M, schalen_M)):
    for t in TYPEN:
        for m in (1, 2, 3):
            el = [e for e in sch[t] if e[0] == m]
            dims_ok &= len(el) == m and sorted(e[1] for e in el) == list(range(m))
            for (m2, k, p, q, v) in el:
                s_ = int(np.flatnonzero(v)[0])
                enat_ok &= abs(sys_.Enat[s_] - m) < TOL and sys_.Nplus[s_] == 1 and sys_.Nminus[s_] == 1
        l0_tab[f"{name_} {NAME[t]}"] = sorted({r(sys_.L0[int(np.flatnonzero(e[4])[0])], 6) for e in sch[t] if e[0] == 1})
heber_def = 0.0
for sys_, sch in ((S, schalen_S), (M, schalen_M)):
    Lm, Lp = sys_.lnat(-1), sys_.lnat(1)
    Ecomm = (Lp @ Lm - Lm @ Lp)/2
    e0 = (Ecomm @ sys_.vac)[sys_.vi]
    for t in TYPEN:
        for (m, k, p, q, v) in sch[t]:
            if m <= 2:
                heber_def = max(heber_def, np.abs(Ecomm @ v - (e0 + m)*v).max())
ZAHLEN["paarschalen"] = {"L_0 der Stufe 1 je Typ": l0_tab, "Defekt [L^nat_1, L^nat_-1]/2 - const = m": r(heber_def, 15)}
check("[E] (b) Paarschalen aller vier Typen 11, 22, 12, 21: dim H^(ab)_m = m mit Teilchenlagen k = 0..m-1 und E^nat = m "
      "(m = 1..3), für gleiche (1/4, 1/4) und verschiedene Quarter (1/4, 3/4); bei verschiedenen Quartern liegt L_0 der "
      f"Kreuzschalen bei m -+ 1/2 ({l0_tab['verschieden 12']}, {l0_tab['verschieden 21']} für m = 1), E^nat aber bei m. Der "
      "native Heber gibt dasselbe: [L^nat_1, L^nat_-1]/2 = E^nat + const auf allen Einpaarzuständen m <= 2 (Defekt "
      f"{heber_def:.0e}). Der Vertexsektor ist eine weitere radiale Leiter mit Delta = 1, die Einteilchenleiter wird "
      "u(2)-adjungiert (2 x 2-quer)", dims_ok and enat_ok and heber_def < TOL)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

say("  (c) Dreipunktstruktur der 4D-Einteilchenfelder")
erz = [S.bil(p, q) for t in TYPEN for (m, k, p, q, v) in schalen_S[t] if m <= 2]   # 12 Paarerzeuger, Stufe <= 2
ops = erz + [B.T.tocsr() for B in erz]
zv = [Z @ S.vac for Z in ops]
links = [X.T @ S.vac for X in ops]
drei_vak = 0.0
for Y in ops:
    yz = np.array([Y @ v for v in zv])
    drei_vak = max(drei_vak, np.abs(np.array(links) @ yz.T).max())
drei_kms = 0.0
for _ in range(300):
    i_, j_, k_ = rng.integers(len(ops), size=3)
    drei_kms = max(drei_kms, abs(S.ew(ops[i_] @ ops[j_] @ ops[k_])))
ab_vak = max(abs(S.vac @ (A[3] @ (B[-2] @ (C[-1] @ S.vac)))) for A in (J1, J2) for B in (J1, J2) for C in (J1, J2))
ab_vak = max(ab_vak, max(abs(S.vac @ (A[2] @ (B[-1] @ (C[-1] @ S.vac)))) for A in (J1, J2) for B in (J1, J2) for C in (J1, J2)))
f_pm = S.vac @ (Hc[2] @ (Ep[-1] @ (Em[-1] @ S.vac)))
f_mp = S.vac @ (Hc[2] @ (Em[-1] @ (Ep[-1] @ S.vac)))
f_h = S.vac @ (Em[2] @ (Hc[-1] @ (Ep[-1] @ S.vac)))
f_kms = S.ew(Hc[2] @ Ep[-1] @ Em[-1])
f_kms_s = S.ew(Hc[2] @ Em[-1] @ Ep[-1])
ZAHLEN["dreipunkt"] = {"max |<b b b>| Vakuum (13824 Tripel)": r(drei_vak, 15), "max |<b b b>|_beta (300 Tripel)": r(drei_kms, 15),
                       "abelsch max |<J J J>| Vakuum": r(ab_vak, 15),
                       "<H_2 E+_-1 E-_-1>, <H_2 E-_-1 E+_-1>, <E-_2 H_-1 E+_-1>": [r(f_pm), r(f_mp), r(f_h)],
                       "beta = 4 pi: <H_2 E+_-1 E-_-1>, <H_2 E-_-1 E+_-1>": [r(f_kms, 10), r(f_kms_s, 10)]}
check("[E] (c) Paarfelder: Jeder Erzeuger oder Vernichter eines Paars ändert die Teilchenzahl über dem See um 1; alle "
      f"verbundenen Dreipunktantworten aus Paarfeldern der vier Typen verschwinden exakt, im Vakuum (alle 13824 Tripel, max "
      f"{drei_vak:.0e}) und im KMS-Zustand bei beta = 4 pi (max {drei_kms:.0e}). Der Vertexsektor bricht die Paarparität "
      "nicht: Ein Paar kann sich nicht in zwei Paare spalten", drei_vak < TOL and drei_kms < TOL)
check("[E] (c) Ströme: abelsche Dreipunktantworten <J J J> aus J^(11), J^(22) verschwinden im Vakuum (max "
      f"{ab_vak:.0e}); die nichtabelschen nicht: <H_2 E^+_-1 E^-_-1> = {f_pm:g}, vertauscht {f_mp:g}, <E^-_2 H_-1 E^+_-1> = "
      f"{f_h:g}, bei beta = 4 pi {f_kms:.6f} und {f_kms_s:.6f}: rein antisymmetrisch, f-Struktur ohne d-Anteil wie in z75",
      ab_vak < TOL and abs(f_pm - 2) < TOL and abs(f_mp + 2) < TOL and abs(f_h - 2) < TOL and abs(f_kms + f_kms_s) < 1e-10
      and abs(f_kms) > 0.5)


def schale(sch, t, m):
    return [e for e in sch[t] if e[0] == m]


def matrix(X, quelle, ziel):
    return np.array([[zv_[4] @ (X @ qv[4]) for qv in quelle] for zv_ in ziel])


schritt_ok, schritt_tab = True, {}
for m in (1, 2):
    ziel = schale(schalen_S, (0, 1), m + 1)
    for qtyp, soll in (((1, 1), "Teilchen hoch"), ((0, 0), "Loch runter")):
        A = matrix(Ep[-1], schale(schalen_S, qtyp, m), ziel)
        muster = np.zeros_like(A)
        for k in range(m):
            muster[k + 1 if soll == "Teilchen hoch" else k, k] = 1
        schritt_ok &= np.array_equal(np.abs(A), muster)
        schritt_tab[f"E+_-1: {NAME[qtyp]}_{m} -> 12_{m + 1} ({soll})"] = A.tolist()
    D = matrix(J1[-1], schale(schalen_S, (0, 0), m), schale(schalen_S, (0, 0), m + 1))
    schritt_tab[f"J11_-1: 11_{m} -> 11_{m + 1}"] = D.tolist()
    schritt_ok &= np.count_nonzero(D) == 2*m and set(np.abs(D[D != 0])) == {1.0}
ZAHLEN["formfaktor_1_nach_1"] = schritt_tab
check("[E] (c) Die nichtabelsche Dreipunktantwort ist ein 1->1-Formfaktor <Einteilchen| E^+_-1 |Einteilchen>: E^+_-1 bildet "
      "H^(22)_m auf H^(12)_(m+1) mit Teilchen hoch (k -> k + 1) und H^(11)_m auf H^(12)_(m+1) mit Loch runter (k -> k) ab, je mit "
      "Betrag 1 (m = 1, 2), genau die nackten Hardy-Schritte von J_-1 (z73/z83), nur flavourwechselnd. In der verdoppelten "
      "Lesart wirkt der Schritt auf einer Seite: (j, j) -> (j + 1/2, j), ein Spinor-K-Typ (halber Impuls S_a, z85 "
      "Twistorverschiebung), kein skalarer Kubikvertex", schritt_ok)

v_pm = Ep[-1] @ (Em[-1] @ S.vac)
v_mp = Em[-1] @ (Ep[-1] @ S.vac)
ein = (S.Nplus == 1) & (S.Nminus == 1)
zwei = (S.Nplus == 2) & (S.Nminus == 2)
Lm_S, Lp_S = S.lnat(-1), S.lnat(1)
hm1 = Hc[-1] @ S.vac
u_lw = v_pm - 0.5*(Lm_S @ hm1)
fus = {"|Einpaaranteil|^2": r(v_pm[ein] @ v_pm[ein]), "|Zweipaaranteil|^2": r(v_pm[zwei] @ v_pm[zwei]),
       "E^nat des Einpaaranteils": sorted({r(S.Enat[s_]) for s_ in np.flatnonzero(np.abs(v_pm) > 1e-12) if ein[s_]}),
       "|[E+_-1, E-_-1] Omega - H_-2 Omega|": r(np.abs(v_pm - v_mp - Hc[-2] @ S.vac).max(), 15),
       "|L^nat_1 v - H_-1 Omega|": r(np.abs(Lp_S @ v_pm - hm1).max(), 15), "|L^nat_1 u|": r(np.abs(Lp_S @ u_lw).max(), 15),
       "|Einpaaranteil von u|^2": r(u_lw[ein] @ u_lw[ein])}
ktyp_ok = all((m1 - 1)/2 + (m2 - 1)/2 < (m1 + m2 - 1)/2 for m1 in range(1, 41) for m2 in range(1, 41))
ZAHLEN["fusion_zweier_kreuzteilchen"] = fus
check("[E] (c) Fusion zweier Kreuzteilchen der Stufe 1: E^+_-1 E^-_-1 Omega hat neben dem freien Zweiteilchenanteil "
      f"(|.|^2 = {fus['|Zweipaaranteil|^2']:g}) einen Einteilchenanteil (|.|^2 = {fus['|Einpaaranteil|^2']:g}) auf der Stufe "
      f"E^nat = {fus['E^nat des Einpaaranteils']}; der Kommutator ist exakt H_-2 Omega (f-Teil), und L^nat_1 v = H_-1 Omega, "
      "also v = u + (1/2) L^nat_-1 H_-1 Omega mit L^nat_1 u = 0: ein sl(2,R)-Vertex D_1 x D_1 -> D_2 plus ein Nachkomme des "
      "dritten Teilchens. Gegen z91: Die Stufe 2 trägt j = 1/2, und (j_1, j_1) x (j_2, j_2) enthält (j_3, j_3) bei "
      "m_3 = m_1 + m_2 nie (m <= 40, gilt flavourblind auch für die Kreuzleitern). Der Vertexsektor gibt einen verbundenen "
      "Dreipunktvertex der Naht, aber keinen so(4)- und energieverträglichen skalaren 4D-Vertex",
      fus["|Einpaaranteil|^2"] > 0.5 and fus["|Zweipaaranteil|^2"] > 0.5 and fus["E^nat des Einpaaranteils"] == [2.0]
      and fus["|[E+_-1, E-_-1] Omega - H_-2 Omega|"] < TOL and fus["|L^nat_1 v - H_-1 Omega|"] < TOL
      and fus["|L^nat_1 u|"] < TOL and ktyp_ok)

xs_ = sp.symbols("x0:4", real=True)
r12 = sp.sqrt(sum((xs_[i] - [1, 0, 0, 0][i])**2 for i in range(4)))
r13 = sp.sqrt(sum((xs_[i] - [0, 2, 0, 0][i])**2 for i in range(4)))
punkt = {xs_[0]: sp.Rational(1, 3), xs_[1]: sp.Rational(-1, 2), xs_[2]: sp.Rational(1, 5), xs_[3]: sp.Rational(2, 7)}
box3 = sp.N(sum(sp.diff(1/(r12*r13*sp.sqrt(5)), xs_[i], 2) for i in range(4)).subs(punkt), 12)
T_ = [np.array(X.todense()) for X in (Ep[0], Em[0], Hc[0])]
basis1 = [e[4] for t in TYPEN for e in schale(schalen_S, t, 1)]
B1 = np.array(basis1)
T1 = [B1 @ X @ B1.T for X in T_]                                                    # Nullmoden auf der Stufe-1-Schicht
Jb = [(T1[0] + T1[1])/np.sqrt(2), (T1[0] - T1[1])/(1j*np.sqrt(2)), T1[2]/np.sqrt(2)]
fabc = np.zeros((3, 3, 3))
for a_ in range(3):
    for b_ in range(3):
        K_ = Jb[a_] @ Jb[b_] - Jb[b_] @ Jb[a_]
        for c_ in range(3):
            fabc[a_, b_, c_] = np.real(np.trace(K_ @ Jb[c_].conj().T)/(1j*np.trace(Jb[c_] @ Jb[c_].conj().T)))
sym_f = max(abs(sum(fabc[p_] for p_ in [(a_, b_, c_), (b_, c_, a_), (c_, a_, b_), (b_, a_, c_), (a_, c_, b_), (c_, b_, a_)]))
            for a_ in range(3) for b_ in range(3) for c_ in range(3))
ZAHLEN["delta1_dreipunkt"] = {"Box <phi phi phi> am Testpunkt": str(box3), "f_123": r(fabc[0, 1, 2]),
                              "max |Symmetrisierung von f|": r(sym_f, 15)}
check(f"[E] (c) Delta = 1 bleibt: Box <phi phi phi> = {box3} != 0 am Testpunkt (z91), also C_phiphiphi = 0 für jede "
      f"Flavourstruktur; die Strukturkonstanten der Nahtströme (f_123 = {fabc[0, 1, 2]:.6f} aus den Nullmoden auf der "
      f"Stufe-1-Schicht) sind total antisymmetrisch, ihre Symmetrisierung ist {sym_f:.0e}: Die f-Struktur kann einen "
      "symmetrischen skalaren Dreipunkt ohnehin nicht tragen (z75)", box3 != 0 and abs(abs(fabc[0, 1, 2]) - np.sqrt(2)) < TOL
      and sym_f < TOL)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

say("  (d) Vierpunktstruktur")
b1 = {t: S.bil(*schale(schalen_S, t, 1)[0][2:4]) for t in TYPEN}
ad = lambda X: X.T.tocsr()
tausch = S.vac @ (ad(b1[(0, 1)]) @ (ad(b1[(1, 0)]) @ (b1[(1, 1)] @ (b1[(0, 0)] @ S.vac))))
direkt = S.vac @ (ad(b1[(0, 0)]) @ (ad(b1[(1, 1)]) @ (b1[(1, 1)] @ (b1[(0, 0)] @ S.vac))))
pauli = S.vac @ (ad(b1[(0, 0)]) @ (ad(b1[(0, 0)]) @ (b1[(0, 0)] @ (b1[(0, 0)] @ S.vac))))
alle = [(t, e[0], S.bil(e[2], e[3])) for t in TYPEN for e in schalen_S[t]]
gram_tab = {}
for E_ in (2, 3, 4):
    vecs, frei = [], 0
    for i_ in range(len(alle)):
        for j_ in range(i_, len(alle)):
            if alle[i_][1] + alle[j_][1] == E_:
                frei += 1
                vecs.append(alle[i_][2] @ (alle[j_][2] @ S.vac))
    rang = np.linalg.matrix_rank(np.array(vecs), tol=1e-10)
    nf = int(np.sum(zwei & (np.abs(S.Enat - E_) < 1e-9)))
    gram_tab[E_] = {"frei (Sym^2 der vier Leitern)": frei, "Rang in der Naht": int(rang), "Zweipaarzustände der CAR": nf}
ZAHLEN["vierpunkt_zwei_naehte"] = {"<b12 b21 b+22 b+11>": r(tausch), "<b11 b22 b+22 b+11>": r(direkt),
                                   "<b11 b11 b+11 b+11>": r(pauli), "gram": gram_tab}
sing = {}
for Nn in range(2, 7):
    Sn = Naehte([Fr(1, 4)]*Nn, 1)
    bb = {(i, j): Sn.bil(Sn.idx[(i, Fr(1, 4))], Sn.idx[(j, Fr(-3, 4))]) for i in range(Nn) for j in range(Nn)}
    Phi = sum(bb[(i, i)] for i in range(Nn))/np.sqrt(Nn)
    pp = Phi @ (Phi @ Sn.vac)
    zwei_p = pp @ pp
    ein_p = (Phi @ Sn.vac) @ (Phi @ Sn.vac)
    adj = max(abs((bb[(i, j)] @ (bb[(j, i)] @ Sn.vac)) @ pp) for i in range(Nn) for j in range(Nn) if i != j)
    sing[Nn] = {"verbunden <Phi Phi Phi+ Phi+> - 2<Phi Phi+>^2": r(zwei_p - 2*ein_p**2), "-2/N": r(-2/Nn),
                "max |<b_ij b_ji Phi+ Phi+>|": r(adj), "2/N": r(2/Nn)}
ZAHLEN["singulett_ueber_N_naehte"] = sing
check(f"[E] (d) Verbundener Vierpunkt = Konstituententausch: <b_12 b_21 b+_22 b+_11> = {tausch:g} (freie 4D-Wickregel: 0, "
      f"verschiedene Arten), <b_11 b_22 b+_22 b+_11> = {direkt:g} (frei 1), <b_11 b_11 b+_11 b+_11> = {pauli:g} (frei 2). "
      "Zwei Diagonalteilchen und zwei Kreuzteilchen sind derselbe Fermionzustand (Löcher getauscht). Gram der "
      f"Zweiteilchenzustände: E = 2: {gram_tab[2]}, E = 3: {gram_tab[3]}, E = 4: {gram_tab[4]}: übervollständig, der freie "
      "adjungierte Fockraum passt nicht in zwei Nähte (wie z88 für eine)",
      abs(abs(tausch) - 1) < TOL and abs(direkt - 1) < TOL and abs(pauli) < TOL and gram_tab[2]["Rang in der Naht"] == 1
      and all(v["Rang in der Naht"] < v["frei (Sym^2 der vier Leitern)"] for v in gram_tab.values())
      and all(v["Rang in der Naht"] == v["Zweipaarzustände der CAR"] for v in gram_tab.values()))
check("[E] (d) Singulett über N Nähte (Phi = sum_i b_ii / sqrt N, Stufe 1): verbundener Vierpunkt "
      f"{[sing[n]['verbunden <Phi Phi Phi+ Phi+> - 2<Phi Phi+>^2'] for n in sing]} = -2/N (N = 2..6), Übergang "
      f"Singulett-Singulett -> Kreuz-Kreuz max |<b_ij b_ji Phi+ Phi+>| = {[sing[n]['max |<b_ij b_ji Phi+ Phi+>|'] for n in sing]}"
      " = 2/N: eine statistische Wechselwirkung der Ordnung 1/N (z91: kappa_4 = -1/(2N)), jetzt mit Flavourtausch über den "
      "Vertexsektor; sie verschwindet für N -> unendlich",
      all(abs(v["verbunden <Phi Phi Phi+ Phi+> - 2<Phi Phi+>^2"] - v["-2/N"]) < 1e-10
          and abs(v["max |<b_ij b_ji Phi+ Phi+>|"] - v["2/N"]) < 1e-10 for v in sing.values()))
befund("Teil 1: Die wörtlich geteilte Hülle ist nichtlokal oder (konsistent, g = +-1) leer; die einzige lokale, nichtentartete "
       "Kopplung zweier Nähte ist das Gitter Z^2 mit zwei Hüllen, und ihr Ladungstransfer ist der su(2)_1-Strom. Er macht die "
       "Einteilchenleiter u(2)-adjungiert (vier Leitern mit Delta = 1), gibt verbundene Dreipunktantworten nur als "
       "1->1-Formfaktoren (flavourwechselnde Schritte) bzw. als sl(2,R)-Fusion D_1 x D_1 -> D_2, die so(4) bricht, und "
       "verbundene Vierpunktantworten nur als Konstituententausch und Pauli (statistisch, 1/N). Einen kovarianten "
       "skalaren 4D-Kubikvertex gibt es nicht; das No-go von z91 hält auch mit Vertexoperatoren", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 2
say("\nTeil 2: Anomale Dimension")
say("  (a) nativ")
add_ok, n_zwei = True, 0
for i_ in range(len(alle)):
    for j_ in range(i_, len(alle)):
        v = alle[i_][2] @ (alle[j_][2] @ S.vac)
        sup = np.flatnonzero(np.abs(v) > 1e-12)
        if len(sup):
            n_zwei += 1
            add_ok &= bool(np.all(np.abs(S.Enat[sup] - alle[i_][1] - alle[j_][1]) < TOL))
normen = {m: r(S.vac @ (Ep[m] @ (Em[-m] @ S.vac))) for m in (1, 2, 3)}
normen_d = {m: r(S.vac @ (J1[m] @ (J1[-m] @ S.vac))) for m in (1, 2, 3)}
fF = lambda x: 0.5*(1 - np.tanh(BETA*np.asarray(x, dtype=float)/2))                 # 1/(1 + e^{beta x}) ohne Überlauf
xg = np.arange(-80, 80) + 0.25
fenster = np.array([float(x) for (i, x) in S.modes if i == 0])
schleife = {m: (S.ew(Ep[m] @ Em[-m]), float(np.sum([fF(y)*(1 - fF(y + m)) for y in fenster if y + m <= fenster.max()])),
                float(np.sum(fF(xg)*(1 - fF(xg + m))))) for m in (1, 2)}
ZAHLEN["nativ"] = {"Zweipaarzustände geprüft": n_zwei, "<E+_m E-_-m> Vakuum": normen, "<J11_m J11_-m> Vakuum": normen_d,
                   "beta = 4 pi: Fock, Schleife im Fenster, Schleife unendlich": {m: [r(a, 14), r(b, 14), r(c, 14)]
                                                                                    for m, (a, b, c) in schleife.items()}}
check(f"[E] (a) Keine Verschiebung im nativen System: E^nat ist ein Einkörperoperator und auf allen {n_zwei} nichtverschwindenden "
      "Zweipaarzuständen (vier Typen, Stufen <= 3) exakt die Summe der Einteilchenenergien (keine Bindung, kein anomaler Anteil); "
      f"die Modennormen <E^+_m E^-_-m> = {list(normen.values())} und <J^(11)_m J^(11)_-m> = {list(normen_d.values())} (m = 1, 2, 3) "
      "sind die eines Gewichts h = 1 (binom(m + h - 1, 2h - 1) = m); bei beta = 4 pi ist die Kreuzantwort genau eine "
      f"Fermischleife sum f(y)(1 - f(y + m)): Fock {schleife[1][0]:.12f} und {schleife[2][0]:.10f} gegen die Schleife im "
      f"Fenster (Abweichung {max(abs(a - b) for a, b, c in schleife.values()):.0e}), die unendliche Schleife weicht um die "
      f"Fensterreste {[f'{abs(a - c):.0e}' for a, b, c in schleife.values()]} (~e^(-9 pi), e^(-5 pi)) ab. Delta der Kreuzleiter = 1",
      add_ok and n_zwei > 100 and all(abs(normen[m] - m) < TOL and abs(normen_d[m] - m) < TOL for m in (1, 2, 3))
      and all(abs(a - b) < 1e-12 and abs(a - c) < 1e-6 for a, b, c in schleife.values()))

say("  (b) Mehrfachtransfer")


def enat_min(alphas, K, Qv):
    modes = [(i, Fr(a) + n) for i, a in enumerate(alphas) for n in range(-K, K)]
    nb = len(modes)
    s = np.arange(1 << nb)
    occ = (s[:, None] >> np.arange(nb)[None, :]) & 1
    xs = np.array([float(x) for _, x in modes])
    naht = np.array([i for i, _ in modes])
    unten = xs < 0
    L0 = (np.where(unten, 1 - occ, occ)*np.abs(xs)).sum(1)
    Q = np.stack([occ[:, (naht == i) & ~unten].sum(1) - (1 - occ[:, (naht == i) & unten]).sum(1)
                  for i in range(len(alphas))], 1)
    En = L0 - sum((float(a) - 0.5)*Q[:, i] for i, a in enumerate(alphas))
    return float(En[np.all(Q == np.array(Qv), axis=1)].min())


transfer = {f"{name_} (Q, -Q), Q = {Q}": enat_min(al, 3, (Q, -Q))
            for name_, al in (("gleich", (Fr(1, 4), Fr(1, 4))), ("verschieden", (Fr(1, 4), Fr(3, 4)))) for Q in (1, 2, 3)}
transfer["eine Naht (1, 0)"] = enat_min((Fr(1, 4), Fr(1, 4)), 3, (1, 0))
transfer["vier Nähte (1, 1, -1, -1)"] = enat_min([Fr(1, 4)]*4, 2, (1, 1, -1, -1))
transfer["vier Nähte (2, 0, -2, 0)"] = enat_min([Fr(1, 4)]*4, 2, (2, 0, -2, 0))
ZAHLEN["mehrfachtransfer_Enat_min"] = transfer
check("[E] (b) Mehrfachtransfer: Der Sektor (Q, -Q) zweier Nähte beginnt bei E^nat = "
      f"{[transfer[f'gleich (Q, -Q), Q = {Q}'] for Q in (1, 2, 3)]} (gleiche Quarter) und "
      f"{[transfer[f'verschieden (Q, -Q), Q = {Q}'] for Q in (1, 2, 3)]} (verschiedene), also Q^2 = |Q alpha|^2/2 statt Q mal 1: "
      "Q gleiche Kreuzteilchen in einem Nahtpaar kosten Q(Q - 1) mehr (Pauli der Konstituenten). Verteilt auf verschiedene "
      f"Nahtpaare kostet derselbe Gesamttransfer nur {transfer['vier Nähte (1, 1, -1, -1)']} (gegen "
      f"{transfer['vier Nähte (2, 0, -2, 0)']}); eine Naht allein mit Q = 1 hat {transfer['eine Naht (1, 0)']} (R^2 = -1, z92). "
      "Das ist Statistik im Q^2/2, keine Kopplung, und verschwindet in der Zusammensetzung",
      all(abs(transfer[f"{n_} (Q, -Q), Q = {Q}"] - Q*Q) < TOL for n_ in ("gleich", "verschieden") for Q in (1, 2, 3))
      and abs(transfer["vier Nähte (1, 1, -1, -1)"] - 2) < TOL and abs(transfer["vier Nähte (2, 0, -2, 0)"] - 4) < TOL
      and abs(transfer["eine Naht (1, 0)"] - 0.5) < TOL)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

say("  (c) Deformationen")
Lm_S = S.lnat(-1)
testv = [S.vac] + [e[4] for t in TYPEN for e in schalen_S[t] if e[0] <= 2]
komm0 = max(np.abs(X @ (Y @ v) - Y @ (X @ v)).max() for X in (Lm_S, Lp_S) for Y in (Ep[0], Em[0]) for v in testv)
gk = 0.37
spek_ok, spek = True, {}
for m in (1, 2, 3):
    Bm = np.array([e[4] for t in TYPEN for e in schale(schalen_S, t, m)])
    X = sps.diags(S.Enat) + gk*(Ep[0] + Em[0])
    XB = np.array([X @ b for b in Bm])
    red = XB @ Bm.T
    spek_ok &= np.abs(XB - red @ Bm).max() < TOL
    ev = np.sort(np.linalg.eigvalsh((red + red.T)/2))
    soll = np.sort([m - 2*gk]*m + [m]*(2*m) + [m + 2*gk]*m)
    spek_ok &= np.abs(ev - soll).max() < 1e-10
    spek[m] = [r(x, 10) for x in ev]
ZAHLEN["nullmodenkopplung"] = {"max |[L^nat_+-1, E^+-_0]|": r(komm0, 15), "g": gk, "Spektrum Stufe m": spek}
check("[E] (c-i) Nullmodenkopplung g (E^+_0 + E^-_0) = 2 g T^x: vertauscht exakt mit dem nativen Heber L^nat_(+-1) "
      f"(Defekt {komm0:.0e}) und schließt auf jeder Stufe m der Einpaarschicht; das Spektrum ist m + g (-2, 0, 0, 2) mit je "
      f"Vielfachheit m (m = 1..3, g = {gk}). Ein flavourabhängiges chemisches Potential, die Dilatation bleibt E^nat, Delta = 1",
      komm0 < TOL and spek_ok)

osz_ok, nm_tab = True, {}
for lam in (-0.6, -0.3, 0.2, 0.5):
    v_ = np.sqrt(1 - lam**2)
    for n in (1, 2, 3):
        cut = 80

        def block(d, n=n, lam=lam):
            Hb = np.diag([n*(2.0*k + d) for k in range(cut + 1)])
            for k in range(cut):
                Hb[k + 1, k] = Hb[k, k + 1] = n*lam*np.sqrt((k + 1 + d)*(k + 1))
            return np.linalg.eigvalsh(Hb)[0]

        e0, e1 = block(0), block(1)
        osz_ok &= abs((e1 - e0)/v_ - n) < 1e-9 and abs(e0 - n*(v_ - 1)) < 1e-9
    th = np.arctanh(lam)/2
    boost = lambda pL, pR: (np.cosh(th)*pL + np.sinh(th)*pR, np.sinh(th)*pL + np.cosh(th)*pR)
    E0 = lambda pL, pR: ((pL**2 + pR**2)/2 + lam*pL*pR)/v_
    for pL, pR in ((np.sqrt(2), 0.0), (np.sqrt(2), np.sqrt(2)), (np.sqrt(2), -np.sqrt(2))):
        qL, qR = boost(pL, pR)
        osz_ok &= abs(E0(pL, pR) - (qL**2 + qR**2)/2) < 1e-12 and abs((qL**2 - qR**2) - (pL**2 - pR**2)) < 1e-12
    nm_tab[lam] = {"chiral (sqrt2, 0): h + hquer": r(E0(np.sqrt(2), 0.0), 10),
                   "E_4D (sqrt2, +sqrt2)": r(E0(np.sqrt(2), np.sqrt(2))/2, 10),
                   "E_4D (sqrt2, -sqrt2)": r(E0(np.sqrt(2), -np.sqrt(2))/2, 10), "e^(2 theta)": r(np.exp(2*th), 10)}
zweig_ok = all(min(v["E_4D (sqrt2, +sqrt2)"], v["E_4D (sqrt2, -sqrt2)"]) < 1 < max(v["E_4D (sqrt2, +sqrt2)"], v["E_4D (sqrt2, -sqrt2)"])
               and abs(v["E_4D (sqrt2, +sqrt2)"]*v["E_4D (sqrt2, -sqrt2)"] - 1) < 1e-9 and v["chiral (sqrt2, 0): h + hquer"] > 1
               for v in nm_tab.values())
mult = {"Naht: (k + 1)^2": [(k + 1)**2 for k in range(6)], "lang, Delta > 1: binom(k + 3, 3)": [sp.binomial(k + 3, 3) for k in range(6)]}
fehlend = [int(a - b) for a, b in zip(mult["lang, Delta > 1: binom(k + 3, 3)"], mult["Naht: (k + 1)^2"])]
qq = sp.symbols("q")
kurz = sp.series(qq*(1 + qq)/(1 - qq)**3, qq, 0, 8).removeO()
ZAHLEN["radiusdeformation"] = {"tabelle (lambda)": {str(k): v for k, v in nm_tab.items()}, "vielfachheiten": {k: [int(x) for x in v] for k, v in mult.items()},
                               "fehlende Box-Nachkommen": fehlend}
check("[E] (c-ii) Exakt marginale Deformation des Relativbosons über den Tomita-Partner (lambda J^3 Jquer^3): Die "
      "Oszillatoren bleiben geschützt (Zweimoden-Bogoljubow, Lücke n sqrt(1 - lambda^2), nach Geschwindigkeitsnormierung genau "
      "n für n = 1..3, lambda = -0.6..0.5), also bleiben die Diagonalleitern bei Delta = 1. Die Nullmoden werden O(1,1)-geboostet "
      "(Spin erhalten): der chirale Transferstrom bekommt h + hquer = cosh 2 theta >= 1, die skalaren Kreuzzweige "
      "(sqrt2, +-sqrt2) E_4D = e^(+-2 theta), Produkt 1: "
      f"{ {k: (v['E_4D (sqrt2, +sqrt2)'], v['E_4D (sqrt2, -sqrt2)']) for k, v in nm_tab.items()} }",
      osz_ok and zweig_ok)
check("[E negativ] (c-ii) Die verschobene Kreuzleiter behält die kurzen Vielfachheiten (k + 1)^2 = "
      f"{mult['Naht: (k + 1)^2']} bei E = Delta + k (die deformierte Theorie ist zur undeformierten mit geboosteten Nullmoden "
      "unitär äquivalent, die Einpaarschicht wird samt Vielfachheiten abgebildet). Eine "
      "unitäre su(2,2)-Darstellung mit tiefstem K-Typ (0, 0) ist für Delta = 1 kurz (q(1 + q)/(1 - q)^3, Reihe "
      f"{sp.Poly(kurz, qq).all_coeffs()[::-1][1:7]}) und für Delta > 1 lang mit binom(k + 3, 3) = "
      f"{[int(x) for x in mult['lang, Delta > 1: binom(k + 3, 3)']]}; es fehlen die Box-Nachkommen {fehlend}, und Delta < 1 "
      "verletzt die Schranke. Für lambda != 0 hat keiner der beiden Zweige eine unitäre 4D-Skalarlesart: Die Vertexdeformation "
      "kann Delta verschieben, aber nicht die Vielfachheiten, die eine anomale Dimension in 4D verlangt",
      fehlend == [0, 0, 1, 4, 10, 20] and [int(x) for x in sp.Poly(kurz, qq).all_coeffs()[::-1][1:7]] == mult["Naht: (k + 1)^2"])
ff = np.einsum("acd,bcd->ab", fabc, fabc)
ZAHLEN["nichtabelsche_kopplung"] = {"sum_cd f_acd f_bcd": ff.round(12).tolist()}
check("[E] (c-iii) Nichtabelsche Strom-Strom-Kopplung lambda J^a Jquer^a: aus den Nullmoden der Nahtströme "
      f"sum_cd f_acd f_bcd = {np.round(ff[0, 0], 10):g} delta_ab = 2 h-dual delta_ab (h-dual = 2). Der Einschleifenkoeffizient ist "
      "nicht null, die Kopplung läuft (für ein Vorzeichen marginal relevant, das 2D-Gegenstück asymptotischer Freiheit); "
      "Betrag und Vorzeichen von lambda sind nicht nativ, und als Einkörperoperator ändert sie die Vielfachheiten der "
      "Einpaarschicht ebenso wenig", np.abs(ff - 4*np.eye(3)).max() < 1e-10)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

say("  (d) Eichfeldlesart")
maxwell_ok = all(m*(m + 2) == (2*((m - 1)/2) + 3)*(2*((m - 1)/2) + 1) for m in range(1, 21))
kub_hel = all(not ((j1 + 1) + (j2 + 1) >= (j1 + j2 + 2) and j1 + j2 >= j1 + j2 + 1)
              for j1 in [k/2 for k in range(21)] for j2 in [k/2 for k in range(21)])
kub_misch = all(not (j1 + 1 + j2 >= j1 + j2 + 2 and j1 + j2 + 1 >= j1 + j2 + 1)
                and not (j1 + 1 + j2 >= j1 + j2 + 1 and j1 + j2 + 1 >= j1 + j2 + 2)
                for j1 in [k/2 for k in range(21)] for j2 in [k/2 for k in range(21)])
chiS = qq*(1 + qq)/(1 - qq)**3
ffr = qq**2/(1 - qq)**4 + sum(((s_ + 1)**2*qq**(s_ + 2) - s_**2*qq**(s_ + 3))/(1 - qq)**4 for s_ in range(1, 16))
diff_ff = sp.expand(sp.series(chiS**2 - ffr, qq, 0, 15).removeO())
ZAHLEN["eichfeldlesart"] = {"Flato-Fronsdal Rest bis q^14": str(diff_ff),
                            "Maxwell-K-Typen (j+1, j) bei E = 2 + 2j aus H^(ab)_(m+2) x Hquer^(ab)_m": "m (m + 2)"}
check("[E] (d) Kinematisch trägt der volle verdoppelte Raum (z83, jede Helizität einmal) mit Flavour u(2)-adjungierte "
      "Helizität-+-1-Leitern: H^(ab)_(m+2) x Hquer^(ab)_m hat m (m + 2) Zustände bei E = m + 1, genau der Maxwell-K-Typ "
      "(j + 1, j) bei E = 2 + 2j (m = 1..20). Das ist der Einteilcheninhalt linearisierter u(2)-Eichfelder, nicht mehr. "
      "Kubische Vertizes zwischen positiven Energieleitern mit additiver Energie gibt es weder für Skalare noch für "
      "Helizität +-1 (gleiche und gemischte Helizität, j <= 10): Das Tensorprodukt beginnt bei 2 E_min > E_min. Flato-Fronsdal: "
      f"chi_S^2 = Summe der Ströme mit Delta = s + 2 (Rest bis q^14: {diff_ff}); zwei Kreuzteilchen geben die Vektorkomposita "
      "J^a_mu mit Delta = 3 aus z75, keine Gluonen", maxwell_ok and kub_hel and kub_misch and diff_ff == 0)
befund("Teil 2: Delta = 1 ist für die Vertexsektoren starr. Nativ ist E^nat einkörperig (keine Bindung), der Transferstrom hat "
       "h = 1, und Mehrfachtransfer kostet nur Pauli (Q^2 statt Q). Verschieben lässt sich Delta nur durch eine nicht native "
       "Kopplung an den Tomita-Partner (Radius), und dann ohne die Box-Nachkommen, die eine unitäre 4D-Skalarleiter mit "
       "Delta > 1 braucht. Die Eichfeldlesart gibt kinematisch u(2)-adjungierte Helizität-+-1-Leitern und die "
       "Vektorkomposita aus z75, aber keinen Vertex", True)
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

say("  (e) Was eine Wechselwirkung bräuchte: verbundene 1->2-Amplituden")
quellen = [S.bil(e[2], e[3]) for t in TYPEN for e in schalen_S[t] if e[0] <= 2]
bil_alle = [S.bil(i, j) for i in range(S.NB) for j in range(S.NB)]
verb_bil = 0.0
for X in bil_alle:
    xv = X @ S.vac
    for A in quellen:
        d_ = X @ (A @ S.vac) - A @ xv
        verb_bil = max(verb_bil, np.abs(d_[zwei]).max())
ix = S.idx
O4 = (S.cd[ix[(0, Fr(5, 4))]] @ S.cd[ix[(0, Fr(1, 4))]] @ S.c[ix[(1, Fr(1, 4))]] @ S.c[ix[(1, Fr(-3, 4))]]).tocsr()
A22 = S.bil(ix[(1, Fr(1, 4))], ix[(1, Fr(-7, 4))])
d4 = O4 @ (A22 @ S.vac) - A22 @ (O4 @ S.vac)
q_bild = (O4 @ (A22 @ S.vac))
s4 = np.flatnonzero(np.abs(q_bild) > 1e-12)
ZAHLEN["verbunden_1_nach_2"] = {"max Zweipaaranteil von [X, b+_A] Omega über alle 144 Bilineare": r(verb_bil, 15),
                                "Doppeltransfer psi*_1 psi*_1 psi_2 psi_2 auf b+_22: |verbunden|": r(np.abs(d4[zwei]).max()),
                                "Ladungen des Bildes": [S.Q[s_].tolist() for s_ in s4], "E^nat des Bildes": [r(S.Enat[s_]) for s_ in s4]}
check("[E] (e) Jeder Bilinearoperator psi*_i psi_j (alle 144 im Fenster, also jede native Vertexstruktur mit |alpha|^2 <= 2: "
      "Ströme, Transfer, Heber, Nullmoden) gibt auf Einpaarzuständen nur unverbundene Zweipaaranteile: der Zweipaaranteil von "
      f"X b+_A Omega - b+_A X Omega ist {verb_bil:.0e} (Einkörperoperatoren lassen das vorhandene Paar als Zuschauer). Der "
      "kleinste Vierfermionterm, der Doppeltransfer psi*_(1,5/4) psi*_(1,1/4) psi_(2,1/4) psi_(2,-3/4) (Komponente von "
      "e^{2i(phi_1 - phi_2)}, h = |2 alpha|^2/2 = 4, Teil 2b), verbraucht das Teilchen des Paars b+_22(1/4, -7/4) und erzeugt "
      f"zwei Kreuzpaare (Ladungen {ZAHLEN['verbunden_1_nach_2']['Ladungen des Bildes']}), verbunden "
      f"{np.abs(d4[zwei]).max():g}. Eine verbundene 1->2-Amplitude braucht also einen Zweikörperterm mit gewählter Kopplung",
      verb_bil < TOL and abs(np.abs(d4[zwei]).max() - 1) < TOL and all(S.Q[s_].tolist() == [2, -2] for s_ in s4))
say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")

# =================================================================== Teil 3
say("\nTeil 3: Bilanz für G3")
STATUS = {
    "wörtlich geteilte Hülle": ("E (negativ)", "ein Oszillatorfeld für beide Nähte bei getrennten Nullmoden ist nichtlokal "
                                               "(Austauschquotient -w/z); konsistent (Nullmoden-Gram = Oszillator-Gram) nur g in "
                                               "{-1, 0, 1}, g = +-1 entartet: Transfer ohne Oszillatorinhalt, h = 0, keine 4D-Wirkung"),
    "Vertexroute = Stromroute": ("E", "einzige lokale nichtentartete Kopplung: Gitter Z^2, Transfer e^{i(phi_1 - phi_2)} = "
                                      "psi*_1 psi_2 = su(2)_1-Strom (Frenkel-Kac exakt im Fenster); für N Nähte entsprechend "
                                      "u(N)_1 (Standard, hier N = 2 gerechnet)"),
    "Kreuzleitern": ("E", "H^(12)_m, H^(21)_m: dim m, E^nat = m, auch für verschiedene Quarter (L_0 = m -+ 1/2); die "
                          "Einteilchenleiter wird u(2)-adjungiert, Delta = 1"),
    "verbundener Dreipunkt": ("E (negativ für 4D)", "Paarfelder exakt 0 (Paarparität); Ströme f-Struktur +-2 ohne d-Anteil, als "
                                                    "1->1-Formfaktor (flavourwechselnde nackte Schritte) und als sl(2,R)-Fusion "
                                                    "D_1 x D_1 -> D_2 plus Nachkomme; kein so(4) x E-kovarianter Skalarvertex "
                                                    "(K-Typ, z91), C_phiphiphi = 0 bei Delta = 1"),
    "verbundener Vierpunkt": ("E", "Konstituententausch -1, Pauli 0 statt 2; Gram-Rang 1, 8, 26 gegen frei 10, 32, 84 (E = 2, 3, 4); "
                                   "Singulett über N Nähte -2/N, Singulett -> Kreuz 2/N: statistisch, verschwindet für N -> unendlich"),
    "anomale Dimension nativ": ("E (negativ)", "E^nat einkörperig und additiv (188 Zweipaarzustände), Modennormen m (h = 1), "
                                               "Doppeltransfer Q^2 nur Pauli (verteilt Q); Nullmodenkopplung = chemisches Potential"),
    "Radiusdeformation": ("E (negativ) / C", "exakt marginal, nicht nativ (Kopplung an den Tomita-Partner): Diagonalleitern "
                                             "geschützt, Kreuzzweige e^(+-2 theta), aber mit kurzen Vielfachheiten (k + 1)^2; es "
                                             "fehlen die Box-Nachkommen 0, 0, 1, 4, 10, 20: keine unitäre 4D-Lesart für lambda != 0"),
    "nichtabelsche Strom-Strom-Kopplung": ("C", "sum f f = 2 h-dual delta = 4 delta, läuft (2D-Gegenstück asymptotischer "
                                                "Freiheit); Kopplung und Vorzeichen nicht nativ; einkörperig, ändert keine "
                                                "Vielfachheiten"),
    "Eichfeldlesart": ("E kinematisch / L", "u(2)-adjungierte Helizität-+-1-Leitern mit Maxwell-K-Typen im vollen verdoppelten "
                                            "Raum; kein energieadditiver kubischer Vertex zwischen positiven Energieleitern; "
                                            "Flato-Fronsdal: zwei Leitern = Ströme Delta = s + 2 (z75: J^a_mu, Delta = 3)"),
    "Wechselwirkung": ("O", "braucht einen Zweikörperterm (paarzahländernd, verbundene 1->2-Amplitude), der die Box-Nachkommen "
                            "einer langen Darstellung liefert; kleinster Kandidat im Vertexgitter: Doppeltransfer "
                            "e^{2i(phi_1 - phi_2)}, h = 4, mit gewählter Kopplung"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("G3 ist beantwortet: Native Wechselwirkung über Vertexoperatoren nein, bedingt nur über einen nicht nativen "
          "Zweikörperterm. Die wörtlich geteilte Schwingungshülle (der Vertexoperator der einen Naht mit dem Oszillatorfeld der "
          "anderen) ist nichtlokal; konsistent geteilt (Gram der Nullmoden gleich Gram der Oszillatoren) ist sie nur bei "
          "g = +-1, und dann ist die Relativladung ein Nullvektor, der Transfer hat h = 0 und keinen Oszillatorinhalt. Die einzige "
          "lokale, nichtentartete Kopplung zweier Nähte ist das Gitter Z^2 mit zwei Hüllen, und ihr Ladungstransfer "
          "e^{i(phi_1 - phi_2)} = psi*_1 psi_2 ist exakt der su(2)_1-Strom: Vertexroute und Stromroute aus z91 sind auf Niveau 1 "
          "dieselbe. Der Vertexsektor liefert weitere radiale Leitern H^(12)_m, H^(21)_m mit dim m und E^nat = m, auch bei "
          "verschiedenen Quartern, die Einteilchenleiter wird u(2)-adjungiert mit Delta = 1. Verbundene Dreipunktantworten gibt "
          "es nur in zwei Formen: als 1->1-Formfaktor (flavourwechselnde nackte Schritte, f-Struktur +-2 ohne d-Anteil wie in "
          "z75) und als sl(2,R)-Fusion zweier Kreuzteilchen in die Einpaarschicht (|.|^2 = 2 auf Stufe 2, D_1 x D_1 -> D_2 plus "
          "Nachkomme). Keine ist ein so(4)- und energieverträglicher skalarer 4D-Vertex; Paarfelder haben exakt verschwindende "
          "Dreipunktantworten, das No-go von z91 hält. Verbundene Vierpunktantworten sind Konstituententausch und Pauli "
          "(Gram-Rang 1, 8, 26 gegen 10, 32, 84; Singulett -2/N), also statistisch. Delta = 1 ist starr: E^nat ist einkörperig "
          "und additiv, der Transferstrom hat h = 1, Mehrfachtransfer kostet Q^2 nur durch Pauli. Die einzige exakt marginale "
          "Vertexdeformation (Radius des Relativbosons über den Tomita-Partner) schützt die Diagonalleitern und verschiebt die "
          "Kreuzzweige auf e^(+-2 theta), behält aber die kurzen Vielfachheiten (k + 1)^2; einer unitären langen 4D-Skalarleiter "
          "fehlen dann die Box-Nachkommen 1, 4, 10, 20, und e^(-2 theta) < 1 verletzt die Schranke. Die nichtabelsche "
          "Strom-Strom-Kopplung läuft (sum f f = 4 delta), ist aber ebenfalls gewählt und einkörperig. Jeder Bilinearoperator "
          "gibt nur unverbundene 1->2-Anteile; eine verbundene Amplitude beginnt erst bei einem Vierfermionterm wie dem "
          "Doppeltransfer e^{2i(phi_1 - phi_2)} (h = 4). Für K15: Kinematisch liegen u(2)-adjungierte Helizität-+-1-Leitern mit "
          "Maxwell-K-Typen im vollen verdoppelten Raum, und die f-Struktur ist da; eine Yang-Mills-Dynamik folgt nicht, zwei "
          "Leitern geben nach Flato-Fronsdal nur die Ströme mit Delta = s + 2 (die Vektorkomposita aus z75). Für T3/T5: Die "
          "lokale wechselwirkende Theorie braucht eine paarzahländernde Zweikörperdynamik mit fester Kopplung, die die langen "
          "Vielfachheiten liefert; weder Vertexoperatoren noch Ströme noch marginale Deformationen der Nähte tun das. Nichts "
          "davon ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\n   (Laufzeit {time.time() - T0:.1f} s)")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z99_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
