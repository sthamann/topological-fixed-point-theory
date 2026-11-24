"""Z67: Lücke (E) der Vogelperspektive -- Familiengewichte gamma und TM1-Brechung aus dem Spin-2-Zustand.

Leitlinie 'Rekonstruktion statt Auswahl': keine neue Auswahlregel. Der Zustand Psi im Fünferträger H0(O(4)) (H6-Marken =
Majorana-Sterne von (|2> - |-2>)/sqrt2, Quartik 1 - z^4, J = 0, tau = i) wirkt über die einzige SL(2)-äquivariante Kontraktion
Psi_AB^CD phi_CD (zweite Überschiebung (Psi, phi)^(2)) auf die Familien H0(K(D)) = H0(O(2)) als symmetrische, spurfreie
3x3-Matrix Q (Weylmatrix auf Lambda^2_+). Frage: Liefert Q die Familiengewichte (gamma, TM1-Brechung, Massenmuster)?

Teil 1  Q exakt. Allgemeine Quartik: tr Q = 0, tr Q^2 = 1152 I, det Q = 27648 J, SL(2)-Kovarianz. H6-Zustand: Eigenwerte
        (24i, -24i, 0). Büschel a + b z^4 (mu4-halbinvariant, J = 0 überall) und (1 - z^4) + t z^2 (V4-invariantes Dublett).
Teil 2  H6-Familienbasis (omega = P dz/(z^4 - 1), Rückzug, Charaktere {-1, i, -i}). Q antivertauscht mit der Vierteldrehung
        und vertauscht mit V4: Die Nullfamilie ist durch Symmetrie die Familie mit Charakter -1; Residuen: die drei
        Eigenfamilien sind die drei Paarungen der vier Marken.
Teil 3  Wörterbuch zum z57-Raum, VOR jedem Datenvergleich festgelegt: Residuen (H1 der punktierten Sphäre = A3) und
        Gitterisometrie A3 = D3 auf die Koordinaten 5, 6, 7 des E8-Bilds der Clocks (alle 48 geprüft). Q ist stets diagonal in
        den Koordinaten und null auf der mu4-Achse; J-verträglich ist die Achse e5 (Richtung von F_J). Q in der zyklischen
        Dreierbasis; Beziehungen zu F_C, F_J, u.
Teil 4  Vorhersagen mit festem Protokoll: (a) gamma, (b) TM1-Brechung gegen die eingefrorenen Werte, (c) Massenmuster.
        Alle Lesarten und Wörterbücher werden gezählt (Look-elsewhere).
Teil 5  Bilanz (E/C/L/O).
"""
import itertools
import json
import random
import re
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import brentq, fsolve
from sympy.physics.quantum.cg import CG

HERE = Path(__file__).parent
REPO = HERE.parents[2]
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
ALTERNATIVEN = []


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def alternative(lesart, woerterbuch, ergebnis):
    ALTERNATIVEN.append({"lesart": lesart, "woerterbuch": woerterbuch, "ergebnis": ergebnis})


x, y, z = sp.symbols("x y z")
I = sp.I
QUAD = [y**2, x*y, x**2]                      # Familienpolynome P = 1, z, z^2 in omega = P dz/(z^4 - 1)


def tv2(f, g):
    """Zweite Überschiebung (f, g)^(2): die SL(2)-äquivariante Kontraktion Psi_AB^CD phi_CD."""
    return sp.expand(sp.diff(f, x, 2)*sp.diff(g, y, 2) - 2*sp.diff(f, x, y)*sp.diff(g, x, y)
                     + sp.diff(f, y, 2)*sp.diff(g, x, 2))


def weyl(Psi):
    M = sp.zeros(3, 3)
    for j, bj in enumerate(QUAD):
        r = sp.Poly(tv2(Psi, bj), x, y)
        for i, bi in enumerate(QUAD):
            M[i, j] = r.coeff_monomial(bi)
    return M


def binform(pz, n):
    return sp.expand(y**n*pz.subs(z, x/y))


def zero(M):
    return all(sp.simplify(v) == 0 for v in M)


# ------------------------------------------------------------------ Teil 1
say("Teil 1: die Weylmatrix Q exakt")
a_, b_, c_, d_, e_ = sp.symbols("a b c d e")
gen = a_*x**4 + 4*b_*x**3*y + 6*c_*x**2*y**2 + 4*d_*x*y**3 + e_*y**4
Qg = weyl(gen)
Iinv = a_*e_ - 4*b_*d_ + 3*c_**2
Jinv = a_*c_*e_ + 2*b_*c_*d_ - a_*d_**2 - b_**2*e_ - c_**3
check("allgemeine Quartik: tr Q = 0, tr Q^2 = 1152 I, det Q = 27648 J (exakte Polynomidentitäten). Also gilt bei J = 0 stets "
      "Eigenwerte (lambda, -lambda, 0)", sp.expand(Qg.trace()) == 0 and sp.expand((Qg**2).trace() - 1152*Iinv) == 0
      and sp.expand(Qg.det() - 27648*Jinv) == 0)
Bf = sp.Matrix([[0, 0, 1], [0, -sp.Rational(1, 2), 0], [1, 0, 0]])     # invariante Bilinearform (polarisierte Diskriminante)
random.seed(67)


def rq():
    return sp.Rational(random.randint(-7, 7), random.randint(1, 5))


p_, q_, r_ = rq(), rq(), rq()
while p_ == 0:
    p_ = rq()
s_ = (1 + q_*r_)/p_
sub = {x: p_*x + q_*y, y: r_*x + s_*y}
F4 = sum(rq()*x**k*y**(4 - k) for k in range(5))
G2 = sum(rq()*x**k*y**(2 - k) for k in range(3))
lhs = tv2(sp.expand(F4.subs(sub, simultaneous=True)), sp.expand(G2.subs(sub, simultaneous=True)))
rhs = sp.expand(tv2(F4, G2).subs(sub, simultaneous=True))
Tg = sp.zeros(3, 3)
for j, bj in enumerate(QUAD):
    pol = sp.Poly(sp.expand(bj.subs(sub, simultaneous=True)), x, y)
    for i, bi in enumerate(QUAD):
        Tg[i, j] = pol.coeff_monomial(bi)
check("Q ist symmetrisch bezüglich der invarianten Form (Q^T B = B Q, allgemein) und SL(2)-kovariant: (g.F, g.G)^(2) = "
      "g.(F, G)^(2) und T_g^T B T_g = B für ein zufälliges g in SL(2, Q)",
      zero(Qg.T*Bf - Bf*Qg) and sp.expand(lhs - rhs) == 0 and p_*s_ - q_*r_ == 1 and zero(Tg.T*Bf*Tg - Bf))

psi_z = 1 - z**4
Psi6 = binform(psi_z, 4)
coef = [sp.Poly(Psi6, x, y).coeff_monomial(x**k*y**(4 - k)) for k in range(5)]
comp = [coef[k]/sp.sqrt(sp.binomial(4, k)) for k in range(5)]          # Komponente bei m = k - 2
norm = sp.sqrt(sum(abs(cc)**2 for cc in comp))
comp = [sp.simplify(cc/norm) for cc in comp]
stars = set(sp.solve(psi_z, z))
Q6 = weyl(Psi6)
I6 = Iinv.subs({a_: coef[4], b_: 0, c_: 0, d_: 0, e_: coef[0]})
J6 = Jinv.subs({a_: coef[4], b_: 0, c_: 0, d_: 0, e_: coef[0]})
say(f"   H6-Zustand 1 - z^4 = y^4 - x^4: Komponenten m = -2..2: {comp}; Sterne {sorted(stars, key=str)}; I = {I6}, J = {J6}")
say(f"   Q(H6) in der Basis P = (1, z, z^2): {Q6.tolist()}, Eigenwerte {Q6.eigenvals()}")
check("H6: Psi = -(|2> - |-2>)/sqrt2 (derselbe Strahl), Sterne = mu4, I = -1, J = 0; Q hat genau die Eigenwerte (24i, -24i, 0)",
      comp[0] == -comp[4] and all(cc == 0 for cc in comp[1:4]) and stars == {1, -1, I, -I} and I6 == -1 and J6 == 0
      and Q6.eigenvals() == {24*I: 1, -24*I: 1, 0: 1})
ZAHLEN["Q_H6_Basis_1_z_z2"] = str(Q6.tolist())
# Normierungen: Weyl Q/24, Einheitszustand, Clebsch-Gordan
QW = Q6/24
cM = {m: comp[m + 2] for m in range(-2, 3)}
Qcg = sp.zeros(3, 3)
for i, m in enumerate((-1, 0, 1)):
    for j, mp in enumerate((-1, 0, 1)):
        Qcg[i, j] = sp.nsimplify(sum(cM[M]*CG(2, M, 1, mp, 1, m).doit() for M in range(-2, 3)))
Tm = sp.diag(1, sp.sqrt(2), 1)                                        # |1,m> -> sqrt(binom(2,k)) z^k, k = m + 1
QcgP = sp.simplify(Tm*Qcg*Tm.inv())
mu = sp.simplify(QcgP[0, 2]/Q6[0, 2])
ev_unit = sorted([complex(v) for v in (QW/norm).eigenvals()], key=lambda v: v.imag)
ev_cg = sorted([complex(v) for v in Qcg.eigenvals()], key=lambda v: v.imag)
say(f"   Normierungen: Weyl Q/24 hat tr = 2I ({sp.simplify((QW**2).trace())}); Einheitszustand: Eigenwerte "
    f"{[round(v.imag, 6) for v in ev_unit]} i; Clebsch-Gordan-Kopplung (2 x 1 -> 1): {[round(v.imag, 6) for v in ev_cg]} i "
    f"(= +-sqrt(3/10)), proportional zu Q mit Faktor {mu}")
check("drei Normierungen (Überschiebung, Weyl, Clebsch-Gordan) geben dieselbe Matrix bis auf einen Faktor: Nur Verhältnisse und "
      "Eigenvektoren von Q sind zustandsbestimmt, die absolute Skala nicht", zero(QcgP - mu*Q6) and mu != 0
      and sp.simplify((QW**2).trace() - 2*I6) == 0 and abs(abs(ev_cg[0]) - np.sqrt(0.3)) < 1e-12)
ZAHLEN["lambda_Normierungen"] = {"Ueberschiebung": "24i", "Weyl Q/24": "i", "Einheitszustand Weyl": "i/sqrt2",
                                 "Clebsch-Gordan": "i sqrt(3/10)"}

Ap, Bp, sq = sp.symbols("A B s", complex=True)
QP = weyl(binform(Ap + Bp*z**4, 4))
lam = sp.symbols("lambda")
cpP = sp.factor((QP - lam*sp.eye(3)).det())
JP = Jinv.subs({a_: Bp, b_: 0, c_: 0, d_: 0, e_: Ap})
QPs = QP.subs(Bp, Ap*sq**2)                                           # B = A s^2, also s = sqrt(B/A) ohne Zweigwahl
vpm = [sp.Matrix([1, 0, sq]), sp.Matrix([1, 0, -sq])]
rr_, ph_ = sp.symbols("r varphi", positive=True)
s_abs2 = sp.simplify(abs(rr_*sp.exp(I*ph_))**2)                       # |s|^2 für s = r e^(i varphi)
wgt = [sp.simplify(1/(1 + s_abs2)), sp.simplify(s_abs2/(1 + s_abs2))]
eq_w = sp.solve(sp.Eq(wgt[0], wgt[1]), rr_)
say(f"   Büschel A + B z^4: Q = {QP.tolist()}, charakteristisches Polynom {cpP}, J = {JP}; mit B = A s^2, s = r e^(i varphi): "
    f"Eigenvektoren 1 +- s z^2 mit Gewichten {wgt} auf den Familien (1, z^2); gleich genau für r = {eq_w}")
check("mu4-Büschel A + B z^4: J = 0 für alle A, B; Q z = 0 für alle A, B (Nullfamilie immer P = z); +-lambda = +-24 A s auf "
      "1 +- s z^2 (s^2 = B/A) -- gleiche Gewichte 1/2 auf den Familien 1 und z^2 genau für |A| = |B| (H6: Marken auf dem Einheitskreis)",
      JP == 0 and zero(QP*sp.Matrix([0, 1, 0])) and sp.expand(cpP + lam*(lam**2 - 576*Ap*Bp)) == 0
      and zero(QPs*vpm[0] - 24*Ap*sq*vpm[0]) and zero(QPs*vpm[1] + 24*Ap*sq*vpm[1]) and eq_w == [1])
t = sp.symbols("t")
QV = weyl(binform(1 - z**4 + t*z**2, 4))
V4b = [sp.Matrix([1, 0, I]), sp.Matrix([0, 1, 0]), sp.Matrix([1, 0, -I])]        # 1 + i z^2, z, 1 - i z^2
evV = [4*t + 24*I, -8*t, 4*t - 24*I]
JV = sp.factor(Jinv.subs({a_: -1, b_: 0, c_: t/6, d_: 0, e_: 1}))
tz = set(sp.solve(JV, t))
say(f"   V4-Dublett (1 - z^4) + t z^2: Q = {QV.tolist()}; Eigenwerte auf (1 + i z^2, z, 1 - i z^2): {evV}; J = {JV}")
zero_fam = {}
for tv in tz:
    zero_fam[str(tv)] = [k for k in range(3) if sp.simplify(evV[k].subs(t, tv)) == 0]
rad = sorted({round(float(abs(complex(rt))), 9) for rt in sp.solve(1 - z**4 + 6*I*z**2, z)})
say(f"   J = 0 genau für t in {sorted(tz, key=str)}; Nullfamilie je t (0 = 1 + i z^2, 1 = z, 2 = 1 - i z^2): {zero_fam}; "
    f"Sternradien bei t = 6i: {rad}")
check("V4-Büschel: Q ist für jedes t diagonal in derselben Basis (1 + i z^2, z, 1 - i z^2) mit Eigenwerten (4t + 24i, -8t, "
      "4t - 24i); J = -t(t^2 + 36)/216 verschwindet genau bei t = 0, +-6i. Bei t = 0 ist die Nullfamilie z; bei t = +-6i wandert "
      "sie zu 1 -+ i z^2, und die Sterne liegen nicht mehr auf dem Einheitskreis (nicht die H6-Marken)",
      all(zero(QV*V4b[k] - evV[k]*V4b[k]) for k in range(3)) and sp.simplify(JV + t*(t**2 + 36)/216) == 0
      and tz == {0, 6*I, -6*I} and zero_fam["0"] == [1] and zero_fam[str(6*I)] == [2] and zero_fam[str(-6*I)] == [0]
      and np.allclose(rad, [np.sqrt(2) - 1, np.sqrt(2) + 1]))
befund("Q ist exakt und zustandsbestimmt bis auf einen Faktor. Bei J = 0 (tau = i) hat es die Eigenwerte (lambda, -lambda, 0); "
       "welche Familie die Null trägt, legt die Lage des Marken-Quadrats fest", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: H6-Familienbasis, Symmetrie, Residuen")


def mob(g, zz):
    al, be, ga, de = g
    return (al*zz + be)/(ga*zz + de)


GRP = {"r": (I, 0, 0, 1), "r2": (-1, 0, 0, 1), "s": (0, 1, 1, 0), "rs": (0, I, 1, 0), "r3s": (0, -I, 1, 0)}


def fam_rep(g):           # Rückzug von omega = P dz/(z^4 - 1)
    M = sp.zeros(3, 3)
    gz = mob(g, z)
    dg = sp.diff(gz, z)
    for k in range(3):
        Pt = sp.Poly(sp.expand(sp.cancel(gz**k*dg*(z**4 - 1)/(gz**4 - 1))), z)
        assert Pt.degree() <= 2
        for i in range(3):
            M[i, k] = Pt.coeff_monomial(z**i)
    return M


def car_rep(g):           # Rückzug von f = Q/(z^4 - 1), deg Q <= 4
    M = sp.zeros(5, 5)
    gz = mob(g, z)
    for k in range(5):
        Pt = sp.Poly(sp.expand(sp.cancel(gz**k*(z**4 - 1)/(gz**4 - 1))), z)
        for i in range(5):
            M[i, k] = Pt.coeff_monomial(z**i)
    return M


RF = {n: fam_rep(g) for n, g in GRP.items()}
RC = {n: car_rep(g) for n, g in GRP.items()}
psi_vec = sp.Matrix([1, 0, 0, 0, -1])
chi = {}
for n in GRP:
    A1, A2 = Q6*RF[n], RF[n]*Q6
    k0 = next(i for i in range(9) if A2[i] != 0)
    cval = sp.simplify(A1[k0]/A2[k0])
    chi[n] = cval if zero(A1 - cval*A2) else None
say(f"   Vierteldrehung r auf (1, z, z^2): {RF['r'].tolist()}; Charakter der Vertauschung Q rho(g) = chi rho(g) Q: {chi}")
check("Familiencharaktere aus dem Rückzug P -> i P(iz): (i, -1, -i) auf (1, z, z^2), also {-1, i, -i}; der Zustand 1 - z^4 (die "
      "konstante Trägerfunktion) ist unter ganz D4 invariant",
      RF["r"] == sp.diag(I, -1, -I) and all(zero(RC[n]*psi_vec - psi_vec) for n in GRP))
check("Q antivertauscht mit der Vierteldrehung (chi(r) = -1) und vertauscht mit V4 = {r^2, rs, r^3 s} (chi = 1): Q bildet die "
      "Familie i auf -i ab, und die Familie mit Charakter -1 muss auf den Eigenraum +1 gehen, den es nicht gibt -- Q z = 0 ist "
      "durch Symmetrie erzwungen", chi["r"] == -1 and chi["r2"] == 1 and chi["rs"] == 1 and chi["r3s"] == 1
      and zero(Q6*sp.Matrix([0, 1, 0])))
v4ch = []
for v in V4b:
    ch = []
    for n in ("r2", "rs"):
        w_ = RF[n]*v
        k0 = next(i for i in range(3) if v[i] != 0)
        ch.append(sp.simplify(w_[k0]/v[k0]))
        assert zero(w_ - ch[-1]*v)
    v4ch.append(tuple(ch))
evQ = [sp.simplify((Q6*v)[next(i for i in range(3) if v[i] != 0)]/v[next(i for i in range(3) if v[i] != 0)]) for v in V4b]
say(f"   V4-Charaktere (r^2, rs) von (1 + i z^2, z, 1 - i z^2): {v4ch}; Q-Eigenwerte darauf: {evQ}")
check("die Eigenfamilien von Q sind genau die V4-Charakterbasis (drei verschiedene nichttriviale Charaktere = 'Farbe = Familien "
      "als V4-Darstellung', Z34): Nullfamilie z, +-24i auf 1 +- i z^2, also je gleiche Gewichte 1/2 auf den Familien i und -i",
      len(set(v4ch)) == 3 and (1, 1) not in v4ch and all(zero(Q6*V4b[k] - evQ[k]*V4b[k]) for k in range(3))
      and evQ == [24*I, 0, -24*I])
Pchar = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])                  # Spalten: z (-1), 1 (i), z^2 (-i)
Qchar = Pchar.T*Q6*Pchar
say(f"   Q in der H6-Familienbasis (Charaktere -1, i, -i): {Qchar.tolist()}; hermitesch normiert Q/(24i): "
    f"{(Qchar/(24*I)).tolist()}")
ZAHLEN["Q_H6_Charakterbasis_(-1,i,-i)"] = str(Qchar.tolist())
ZAHLEN["Gewichte_Eigenfamilien_in_Charakterbasis"] = {"0": [1, 0, 0], "+lambda": [0, 0.5, 0.5], "-lambda": [0, 0.5, 0.5]}
MARKS = [sp.Integer(1), I, sp.Integer(-1), -I]
Res = sp.Matrix(4, 3, lambda n, k: MARKS[n]**(k + 1)/4)                # Res_zeta(P dz/(z^4-1)) = zeta P(zeta)/4
pat = []
for v in V4b:
    rv = Res*v
    k0 = next(i for i in range(4) if rv[i] != 0)
    pat.append(tuple(sp.simplify(rv[i]/rv[k0]) for i in range(4)))
say(f"   Residuenmuster der Eigenfamilien auf den Marken (1, i, -1, -i): {pat}")
check("Residuenabbildung: Bild = summenfreie Funktionen auf den Marken (Rang 3, Residuensatz). Die Eigenfamilien von Q gehen "
      "auf die drei reellen Paarungsmuster; die Nullfamilie auf (1, -1, 1, -1), die Antipodenpaarung {1, -1} | {i, -i}",
      Res.rank() == 3 and zero(sp.ones(1, 4)*Res) and pat[1] == (1, -1, 1, -1)
      and sorted([pat[0], pat[2]]) == sorted([(1, 1, -1, -1), (1, -1, -1, 1)]))
befund("Die H6-Familienbasis (Charaktere) ist nicht die Eigenbasis von Q: Q mischt die Familien i und -i maximal und lässt die "
       "Familie -1 als Nullfamilie stehen. Die Eigenbasis ist die reelle V4-Basis der drei Markenpaarungen", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Wörterbuch zum z57-Raum (festgelegt vor jedem Datenvergleich)")
S = sp.zeros(4, 4)                         # Rückzug auf Residuen: (S res)(zeta) = res(i zeta)
for n, zt in enumerate(MARKS):
    S[n, MARKS.index(sp.expand(I*zt))] = 1
check("die Vierteldrehung wirkt auf den Residuen als zyklische Verschiebung der Marken, verträglich mit dem Rückzug der "
      "Familien (S Res = Res rho(r))", zero(S*Res - Res*RF["r"]))
Qres = Res*Q6*(Res.H*Res).inv()*Res.H
verts = sp.Matrix([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])
Jb = sp.Matrix([[0, -1], [1, 0]])          # J auf (e6, e7): e6 -> e7, e7 -> -e6
klassen, w1_maps, lattice_ok, chosen = {}, [], True, None
for perm in itertools.permutations(range(4)):
    for sg in (1, -1):
        phi = sp.Matrix(4, 3, lambda n, k: sg*verts[perm[n], k]/2)
        ok = phi.T*phi == sp.eye(3)
        for n1, n2 in itertools.combinations(range(4), 2):
            vv = phi.T*(sp.eye(4)[:, n1] - sp.eye(4)[:, n2])
            ok &= sorted(abs(v) for v in vv) == [0, 1, 1]
        lattice_ok &= ok
        Mr = phi.T*S*phi
        axis = [k for k in range(3) if Mr[k, k] == -1]
        Qz = sp.simplify(phi.T*Qres*phi)
        diag = tuple(sp.simplify(Qz[k, k]/(24*I)) for k in range(3))
        isdiag = zero(Qz - sp.diag(*[Qz[k, k] for k in range(3)]))
        sense = None
        if axis == [0]:
            sense = "J" if Mr[1:, 1:] == Jb else ("J^-1" if Mr[1:, 1:] == -Jb else None)
        key = (axis[0] if len(axis) == 1 else None, diag, sense, isdiag)
        klassen[key] = klassen.get(key, 0) + 1
        if axis == [0] and sense == "J" and diag == (0, 1, -1) and chosen is None:
            chosen = (perm, sg, phi)
say("   alle 48 Gitterisometrien A3 (summenfrei auf den Marken) -> D3 (Koordinaten e5, e6, e7):")
for k_, v_ in sorted(klassen.items(), key=str):
    say(f"     mu4-Achse e{5 + k_[0]}, Q/(24i) = diag{k_[1]}, Drehsinn auf (e6, e7): {k_[2]}, diagonal: {k_[3]} -- {v_} Abbildungen")
check("alle 48 Abbildungen sind Gitterisometrien (Wurzeln e_i - e_j -> +-e_a +- e_b). Für jede ist Q diagonal in den "
      "Koordinaten (e5, e6, e7), null auf der mu4-Achse und +-24i auf den beiden anderen; je Achse 16 Abbildungen",
      lattice_ok and all(k_[3] for k_ in klassen) and sum(klassen.values()) == 48
      and all(sorted(k_[1], key=str) == sorted((0, 1, -1), key=str) and k_[1][k_[0]] == 0 for k_ in klassen)
      and all(sum(v_ for k_, v_ in klassen.items() if k_[0] == ax) == 16 for ax in range(3)))
nJ = {k_[1]: v_ for k_, v_ in klassen.items() if k_[0] == 0 and k_[2] == "J"}
check("J-verträglich (mu4-Achse = e5, die Richtung von F_J, und Drehsinn wie J auf (e6, e7)) sind 8 Abbildungen; sie geben beide "
      "Vorzeichen von Q (je 4). Das Vorzeichen ist also nicht festgelegt -- wie beim Strahl Psi = -Psi",
      sum(nJ.values()) == 8 and nJ.get((0, 1, -1)) == 4 and nJ.get((0, -1, 1)) == 4)
perm_c, sg_c, PHI = chosen
QH = sp.diag(0, 1, -1)                      # Q/(24i) im Wörterbuch W1; das Vorzeichen geht in kappa ein
say(f"   Festgelegt: W1 = Residuen + Isometrie (Marken 1, i, -1, -i -> Ecken {[list(verts[k, :]*sg_c) for k in perm_c]}), "
    f"Q_H = Q/(24i) = diag(0, 1, -1) in (e5, e6, e7); Vorzeichen frei (kappa reell, beide Vorzeichen)")
omega = sp.Rational(-1, 2) + sp.sqrt(3)*I/2
Fz = sp.Matrix(3, 3, lambda j, k: omega**(j*k)/sp.sqrt(3))
Qcyc = sp.simplify(Fz.H*QH*Fz)
S3 = sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
say(f"   Q_H in der zyklischen Dreierbasis (z57, e = f_0): {Qcyc.tolist()}")
check("in der zyklischen Dreierbasis ist Q_H = (i/sqrt3)(S - S^T) mit der zyklischen Verschiebung S: nullte Diagonale, rein "
      "imaginär-antisymmetrisch (hermitesch)", zero(Qcyc - I/sp.sqrt(3)*(S3 - S3.T)))
ZAHLEN["Q_H_zyklische_Basis"] = str(Qcyc.tolist())
FC = sp.Matrix([[13, 1, -1], [1, 5, -5], [-1, -5, 5]])/16
FJ = sp.diag(1, 0, 0)
u = sp.Matrix([0, 1, 1])/sp.sqrt(2)
wv = sp.Matrix([0, 1, -1])/sp.sqrt(2)
e5 = sp.Matrix([1, 0, 0])
f0 = Fz[:, 0]
TBM = sp.Matrix.hstack(u, e5, wv)
check("Beziehungen im z57-Raum: Q_H F_J = F_J Q_H = 0, tr(Q_H F_J) = tr(Q_H F_C) = 0, [Q_H, F_C] != 0; Q_H u = w; w hat keinen "
      "e-Anteil (|<f_0|w>|^2 = 0), e5 ist trimaximal (1/3), u trägt 2/3; in der TBM-Basis (u, e5, w) hat Q_H nur den "
      "Eintrag u-w (Diagonale null); Q_H^2 = I - F_J",
      zero(QH*FJ) and zero(FJ*QH) and (QH*FJ).trace() == 0 and sp.simplify((QH*FC).trace()) == 0 and not zero(QH*FC - FC*QH)
      and zero(QH*u - wv) and sp.simplify(abs((f0.H*wv)[0])**2) == 0 and sp.simplify(abs((f0.H*e5)[0])**2) == sp.Rational(1, 3)
      and sp.simplify(abs((f0.H*u)[0])**2) == sp.Rational(2, 3) and zero(TBM.T*QH*TBM - sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]))
      and zero(QH**2 - (sp.eye(3) - FJ)))
al, be = sp.symbols("alpha beta")
Qstate = sp.simplify(PHI.T*Res*weyl(binform(al*(1 - z**4) + be*z**2, 4))*(Res.H*Res).inv()*Res.H*PHI)
check("V4-Dublett unter W1: Q(alpha(1 - z^4) + beta z^2) = 24i alpha Q_H + beta (4 I - 12 F_J). Das F_J-Gewicht steckt genau im "
      "z^2-Anteil des Zustands", zero(Qstate - (24*I*al*QH + be*(4*sp.eye(3) - 12*FJ))))
phiM2 = sp.Matrix.hstack(sp.Matrix([1, -1, 1, -1])/2, sp.Matrix([1, 1, -1, -1])/2, sp.Matrix([1, -1, -1, 1])/2) * \
    sp.Matrix.hstack(e5, u, wv).T                                     # Muster -> (e5, u, w): nicht gittertreu
QM2 = sp.simplify(phiM2.T*Qres*phiM2/(24*I))
root_img = phiM2.T*(sp.eye(4)[:, 0] - sp.eye(4)[:, 1])
check("Gegenprobe M2 (Muster auf e5, u, w gedreht, um 45 Grad gegen W1): Q/(24i) = |u><u| - |w><w|, erhält u (TM1) und ist ein "
      "reiner gamma-Term; M2 ist aber keine Gitterabbildung (Wurzel e_1 - e_i -> nicht ganzzahlig)",
      zero(QM2 - (u*u.T - wv*wv.T)) and any(not sp.sympify(v).is_integer for v in root_img))
befund("Unter dem gittertreuen, J-verträglichen Wörterbuch W1 ist Q orthogonal zur J-Clock (Q F_J = 0) und koppelt nur die "
       "TM1-Spalte u an w: Es ist die TM2-Richtung (die trimaximale Spalte e5 bleibt unberührt). Das Quadrat Q^2 = I - F_J und "
       "der z^2-Anteil des Zustands sind reine gamma-Terme", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Vorhersagen mit festem Protokoll")
PROTOKOLL = {
    "Woerterbuch": "W1: Residuen H0(K(D)) -> summenfrei(Marken) = A3, Isometrie A3 = D3 auf (e5, e6, e7), mu4-Achse = e5 "
                   "(J-verträglich); Q_H = Q/(24i) = diag(0, 1, -1); Vorzeichen frei (kappa reell)",
    "Flavourbasis": "zyklische Dreierbasis aus z57 (e = f_0, mu = f_1, tau = f_2), bedingt wie in z57",
    "Familienoperator": "H = F_C + gamma F_J + kappa Q_H (F_C, F_J aus z57)",
    "Etiketten": "nu3 = Spalte mit kleinstem |U_e|^2, nu1 = übrige mit größtem |U_e|^2 (datenartig, ohne Eigenwerte)",
    "Lesarten": {"R0": "Beitrag von Q zu gamma", "R1": "gamma := Q-Gewicht auf der F_J-Richtung (= 0), kappa an s13 angepasst",
                 "R1'": "gamma = 1/8 (Analyse), kappa an s13 angepasst", "R2": "gamma und kappa frei (reine Anpassung)",
                 "R3": "Q selbst als Massenoperator", "R4": "phasenfreies Gewicht kappa' Q^2 >= 0"},
    "Daten": "eingefroren sin^2 theta12 = 1/3 - phi0/2, sin^2 theta13 = phi0 e^(-5/6), theta23 = 45 Grad und delta = 240 Grad "
             "zugewiesen; NuFIT 6.1 (NO) sin^2 theta12 = 0.3088 +- 0.0067; sin^2 theta13 (NuFIT 6.0) aus predictions_frozen.json",
    "Reihenfolge": "Wörterbuch und Lesarten aus der Struktur (Teil 3) vor dem Datenvergleich; Prototypen liefen vorher, daher "
                   "werden alle Lesarten und Wörterbücher gezählt",
}
for k_, v_ in PROTOKOLL.items():
    say(f"   {k_}: {v_}")
frozen = json.loads((REPO / "verification/predictions_frozen.json").read_text())
pr = {p["id"]: p for p in frozen["predictions"]}
tex = {p["id"]: p for p in frozen["assigned_texture_values"]}
s12f, s13f = float(pr["SIN2_THETA12_SEED"]["frozen_value"]), float(pr["SIN2_THETA13"]["frozen_value"])
dlt, s23a = float(tex["DELTA_CP_NU_DEG"]["frozen_value"]), float(tex["SIN2_THETA23"]["frozen_value"])
m13 = re.search(r"NuFIT 6\.0: ([0-9.]+) \+- ([0-9.]+)", pr["SIN2_THETA13"]["experiment"])
nf13 = (float(m13.group(1)), float(m13.group(2)))
NF61 = (0.3088, 0.0067)
FCn, FJn, QHn = (np.array(M.evalf(), float) for M in (FC, FJ, QH))
Fn = np.array(Fz.evalf(), complex)
un, wn, e5n = (np.array(v.evalf(), float).ravel() for v in (u, wv, e5))


def mix(H):
    w, V = np.linalg.eigh(H)
    U = Fn.conj().T@V
    A = abs(U)**2
    i3 = int(np.argmin(A[0]))
    rest = [i for i in range(3) if i != i3]
    i1 = max(rest, key=lambda i: A[0, i])
    i2 = [i for i in rest if i != i1][0]
    s13 = A[0, i3]
    Uo = U[:, [i1, i2, i3]]
    s12, s23 = A[0, i2]/(1 - s13), A[1, i3]/(1 - s13)
    Jcp = np.imag(Uo[0, 0]*Uo[1, 1]*np.conj(Uo[0, 1])*np.conj(Uo[1, 0]))
    jmax = np.sqrt(s12*(1 - s12)*s23*(1 - s23))*(1 - s13)*np.sqrt(s13)
    return dict(s12=s12, s13=s13, s23=s23, col1=A[:, i1], J=Jcp, sind=Jcp/jmax, ev=w[[i1, i2, i3]], V=V[:, [i1, i2, i3]])


def roots(fun, lo=-1.0, hi=1.0, n=4001):
    ks = np.linspace(lo, hi, n)
    v = [fun(k) for k in ks]
    out = []
    for i in range(n - 1):
        if v[i]*v[i + 1] < 0:
            r = brentq(fun, ks[i], ks[i + 1])
            if abs(fun(r)) < 1e-12:
                out.append(r)
    return out


def sig(s12):
    return (s12 - NF61[0])/NF61[1]


tm1 = (1 - 3*s13f)/(3*(1 - s13f))
say(f"\n   Daten: eingefroren s12 = {s12f:.6f}, s13 = {s13f:.6f}, s23 = {s23a} (zugewiesen), delta = {dlt:g} Grad (zugewiesen); "
    f"TM1 bei diesem s13: s12 = {tm1:.6f}; NuFIT 6.1: {NF61[0]} +- {NF61[1]}; NuFIT 6.0 s13 = {nf13[0]} +- {nf13[1]}")

say("\n   (a) gamma")
th0 = 0.5*np.arctan(2*np.sqrt(2)/3)
m0 = mix(FCn)
say(f"   R0: Q_H hat auf e5 (F_J-Richtung) das Gewicht {QH[0, 0]}; reines TM1 mit gamma = 0: theta = {np.degrees(th0):.3f} Grad, "
    f"s13 = {np.sin(th0)**2/3:.6f} (eingefroren {s13f:.6f})")
check("R0: Q liefert bei tau = i (J = 0) kein gamma: Q_H F_J = 0. Reines TM1 mit gamma = 0 gibt theta = 21.66 Grad und "
      "s13 = 0.0454, etwa das Doppelte des eingefrorenen Werts", QH[0, 0] == 0 and abs(m0["s13"] - np.sin(th0)**2/3) < 1e-12
      and m0["s13"] > 1.9*s13f)
alternative("R0", "W1", f"gamma_Q = 0; TM1 mit gamma = 0: s13 = {m0['s13']:.5f}")
th_s = np.arcsin(np.sqrt(3*s13f))
g_star = (2*np.sqrt(2)/np.tan(2*th_s) - 3)/16
g_mirror = -3/8 - g_star
g15 = (2*np.sqrt(6) - 3)/16
ms_, mm_ = mix(FCn + g_star*FJn), mix(FCn + g_mirror*FJn)
say(f"   reines TM1 bräuchte für das eingefrorene s13 gamma* = {g_star:.5f} (15 Grad: {g15:.5f}, Analyse 1/8); Spiegelast "
    f"gamma -> -3/8 - gamma (theta -> 90 Grad - theta, nu2 und nu3 vertauscht): {g_mirror:.5f}")
check("Gegenprobe: gamma* reproduziert im reinen Clock-Operator das eingefrorene s13, und der Spiegelwert -3/8 - gamma* gibt "
      "dieselben Beträge (z57 liest die Eigenwerte nicht): gamma ist in z57 nur bis auf diese Spiegelung bestimmt",
      abs(ms_["s13"] - s13f) < 1e-12 and abs(mm_["s13"] - s13f) < 1e-12 and abs(ms_["s12"] - mm_["s12"]) < 1e-12
      and abs(ms_["s12"] - tm1) < 1e-12)
ZAHLEN["gamma_stern"], ZAHLEN["gamma_spiegel"] = round(g_star, 6), round(g_mirror, 6)
r4 = roots(lambda k: mix(FCn + k*(np.eye(3) - FJn))["s13"] - s13f, 0.0, 1.0)
kq = [0.0, 0.05, 0.1, 0.15]
s13q = [mix(FCn + kk*(np.eye(3) - FJn))["s13"] for kk in kq]
say(f"   R4: kappa' Q_H^2 = kappa' (I - F_J), also gamma_eff = -kappa'. Direkter Ast (kappa' < 3/16): s13 für kappa' = {kq}: "
    f"{np.round(s13q, 5).tolist()}; Lösungen kappa' >= 0 für das eingefrorene s13: {[round(float(k), 6) for k in r4]}, "
    f"s12 = {[round(float(mix(FCn + k*(np.eye(3) - FJn))['s12']), 6) for k in r4]}")
check("R4: Das phasenfreie Gewicht kappa' Q^2 (kappa' >= 0) ist ein reiner gamma-Term mit gamma_eff = -kappa' und erhält TM1. "
      "Auf dem direkten Ast wächst s13 von 0.0454 weg; die einzige Lösung liegt auf dem Spiegelast bei kappa' = 3/8 + gamma* "
      "und gibt genau TM1 (s12 = 0.317563): angepasst, keine Brechung, keine Vorhersage",
      all(np.diff(s13q) > 0) and s13q[0] > 0.045 and len(r4) == 1 and abs(r4[0] + g_mirror) < 1e-6
      and abs(mix(FCn + r4[0]*(np.eye(3) - FJn))["s12"] - tm1) < 1e-9)
alternative("R4", "W1", f"nur Spiegelast kappa' = {r4[0]:.5f}, s12 = TM1")
sreal = sp.symbols("s", real=True)
Js = sp.factor(sp.expand(JV.subs(t, I*sreal)))
sz = set(sp.solve(sp.expand(Js), sreal))
say(f"   z^2-Beimischung t = i s (hermitesch, gamma_eff = -kappa s/2): J = {Js}, null genau für s in {sorted(sz)}")
check("Ein gamma != 0 aus dem Zustand verlangt einen z^2-Anteil s mit J != 0 (also tau != i); die einzigen weiteren Nullstellen "
      "s = +-6 sind andere Quadrate mit verschobener Nullfamilie (nicht die H6-Marken)", sz == {0, 6, -6})
befund("(a) Q legt gamma nicht fest: linear ist gamma_Q = 0 (Q F_J = 0), das Quadrat ist ein gamma-Term freier Größe (passt "
       "nur angepasst und nur auf dem Spiegelast), und ein z^2-Anteil des Zustands bräuchte J != 0. Weder 1/8 noch "
       "(2 sqrt6 - 3)/16 folgt aus dem Spin-2-Zustand", True)

say("\n   (b) TM1-Brechung")
c12, s12r, c13, s13r = np.sqrt(1 - s12f), np.sqrt(s12f), np.sqrt(1 - s13f), np.sqrt(s13f)
c23 = s23 = np.sqrt(0.5)
ed = np.exp(1j*np.radians(dlt))
Ufr = np.array([[c12*c13, s12r*c13, s13r/ed],
                [-s12r*c23 - c12*s23*s13r*ed, c12*c23 - s12r*s23*s13r*ed, s23*c13],
                [s12r*s23 - c12*c23*s13r*ed, -c12*s23 - s12r*c23*s13r*ed, c23*c13]])
col_fr = abs(Ufr[:, 0])**2
col_mt = np.array([col_fr[0], (1 - col_fr[0])/2, (1 - col_fr[0])/2])
dUe1 = col_fr[0] - 2/3
ds12 = s12f - tm1
say(f"   eingefrorene erste Spalte (delta = 240 Grad): {np.round(col_fr, 6).tolist()}, mit mu-tau-Reflexion: "
    f"{np.round(col_mt, 6).tolist()}; TM1 (2/3, 1/6, 1/6); nötig: Delta|U_e1|^2 = {dUe1:+.6f}, Delta s12 = {ds12:+.6f}")
check("Die eingefrorenen Werte verlangen |U_e1|^2 = c12^2 c13^2 = 0.677233 > 2/3 (Delta +0.010566) und s12 um 0.010816 unter "
      "TM1; die mu-tau-Asymmetrie der Spalte kommt allein vom zugewiesenen delta",
      np.allclose(Ufr@Ufr.conj().T, np.eye(3)) and abs(col_fr[0] - (1 - s12f)*(1 - s13f)) < 1e-14 and dUe1 > 0.01
      and ds12 < -0.01 and abs(col_fr[1] - col_fr[2]) > 0.05)
ZAHLEN["eingefroren_erste_Spalte_delta240"] = np.round(col_fr, 6).tolist()
ZAHLEN["noetig_dUe1"], ZAHLEN["noetig_ds12"] = round(dUe1, 6), round(ds12, 6)
r1 = roots(lambda k: mix(FCn + k*QHn)["s13"] - s13f)
r1_list = [(round(float(k), 6), round(float(mix(FCn + k*QHn)["s12"]), 6)) for k in r1]
viable = [k for k in r1 if abs(sig(mix(FCn + k*QHn)["s12"])) < 3]
say(f"   R1 (H = F_C + kappa Q_H): alle Lösungen kappa in [-1, 1] für s13 = eingefroren: {r1_list}")
k1 = viable[0]
M1 = mix(FCn + k1*QHn)
comp1 = [float(M1["V"][:, 0]@vv) for vv in (un, e5n, wn)]
frac = (M1["s12"] - tm1)/ds12
say(f"   R1, einzige Lösung im NuFIT-3-sigma-Bereich: kappa = {k1:.6f}, s12 = {M1['s12']:.6f} ({sig(M1['s12']):+.2f} sigma), "
    f"|U_e1|^2 = {M1['col1'][0]:.6f}, erste Spalte {np.round(M1['col1'], 5).tolist()}, s23 = {M1['s23']:.6f}, "
    f"sin delta = {M1['sind']:+.6f}; nu1 in (u, e5, w) = {np.round(comp1, 5).tolist()}; Anteil an der nötigen Brechung "
    f"{100*frac:.1f} %")
check("R1: von vier Lösungen liegt genau eine im NuFIT-3-sigma-Bereich (kappa = 0.0499). Sie bricht TM1 in die richtige Richtung, "
      "aber nur zu etwa 36 %: s12 = 0.3136 liegt zwischen eingefroren (0.3067) und TM1 (0.3176)",
      len(r1) == 4 and len(viable) == 1 and abs(k1 - 0.049871) < 1e-5 and s12f < M1["s12"] < tm1 and 0.3 < frac < 0.4)
check("R1 ist mu-tau-reflexionssymmetrisch (s23 = 1/2, |sin delta| = 1): Das folgt aus der Reellheit von H in den Koordinaten, "
      "nicht aus Q; delta = 240 Grad (zugewiesen) ist so nicht erreichbar",
      abs(M1["s23"] - 0.5) < 1e-12 and abs(abs(M1["sind"]) - 1) < 1e-9 and abs(M1["col1"][1] - M1["col1"][2]) < 1e-12)
band = []
for s13t in (nf13[0] - nf13[1], nf13[0], nf13[0] + nf13[1]):
    kk = brentq(lambda k: mix(FCn + k*QHn)["s13"] - s13t, 0, 0.15)
    band.append(mix(FCn + kk*QHn)["s12"])
say(f"   R1 mit s13 aus NuFIT 6.0 ({nf13[0]} -+ {nf13[1]}): s12 = {np.round(band, 6).tolist()}")
ZAHLEN["R1"] = {"kappa": round(k1, 6), "s12": round(M1["s12"], 6), "sigma_NuFIT61": round(sig(M1["s12"]), 2),
                "Ue1_quadrat": round(M1["col1"][0], 6), "Anteil_Brechung": round(frac, 3), "nu1_in_u_e5_w": np.round(comp1, 5).tolist(),
                "s12_bei_NuFIT60_s13": np.round(band, 6).tolist(), "alle_Loesungen": r1_list}
alternative("R1", "W1 (Achse e5, beide Vorzeichen)", f"4 Lösungen, 1 im 3-sigma-Bereich: s12 = {M1['s12']:.4f}")
r1p = roots(lambda k: mix(FCn + FJn/8 + k*QHn)["s13"] - s13f)
r1p_list = [(round(float(k), 6), round(float(mix(FCn + FJn/8 + k*QHn)["s12"]), 6)) for k in r1p]
kp = min(r1p, key=abs)
say(f"   R1' (gamma = 1/8 fest): Lösungen {r1p_list}; nächste kappa = {kp:.6f}, s12 = {mix(FCn + FJn/8 + kp*QHn)['s12']:.6f}")
check("R1': Mit gamma = 1/8 braucht es nur kappa = -0.0038, und s12 = 0.3182 bleibt praktisch TM1 -- Q bricht dort nichts Nennenswertes",
      abs(kp) < 0.01 and abs(mix(FCn + FJn/8 + kp*QHn)["s12"] - 0.318156) < 1e-5)
alternative("R1'", "W1", f"kappa = {kp:.4f}, s12 = {mix(FCn + FJn/8 + kp*QHn)['s12']:.4f} (TM1-nah)")


def f2(pv):
    m = mix(FCn + pv[0]*FJn + pv[1]*QHn)
    return [m["s12"] - s12f, m["s13"] - s13f]


sols = set()
for g0 in np.linspace(-2, 1, 31):
    for k0 in np.linspace(-1, 1, 41):
        sv, info, ier, msg = fsolve(f2, [g0, k0], full_output=True)
        if ier == 1 and max(abs(np.array(f2(sv)))) < 1e-11 and abs(sv[0]) < 5 and abs(sv[1]) < 5:
            sols.add((round(float(sv[0]), 6), round(float(sv[1]), 6)))
sols = sorted(sols, key=lambda p: abs(p[1]))
dist = {p: min(abs(p[0] - 0.125), abs(p[0] + 0.5)) for p in sols}
gbest = min(sols, key=lambda p: dist[p])
say(f"   R2 (gamma, kappa frei, an s12 und s13 eingefroren angepasst): {len(sols)} Lösungen {sols}; am nächsten an 1/8 oder "
    f"dessen Spiegelwert -1/2: {gbest} (Abstand {dist[gbest]:.3f})")
check("R2 passt exakt, aber nicht eindeutig (mehrere Lösungen), und jede braucht kappa != 0 neben einem verschobenen gamma: "
      "Das ist eine Anpassung mit zwei Parametern an zwei Zahlen, keine Vorhersage",
      len(sols) >= 2 and all(abs(p[1]) > 0.01 for p in sols) and all(dist[p] > 0.05 for p in sols))
ZAHLEN["R2_Loesungen"] = sols
alternative("R2", "W1", f"{len(sols)} exakte Anpassungen, nächst an 1/8 bzw. -1/2: {gbest}")
for nm, Qa in (("W1-Achse e6", np.diag([1., 0, -1])), ("W1-Achse e7", np.diag([1., -1, 0])),
               ("M2 (nicht gittertreu)", np.array(QM2.evalf(), float))):
    rr = roots(lambda k: mix(FCn + k*Qa)["s13"] - s13f)
    vals = [(round(float(k), 5), round(float(mix(FCn + k*Qa)["s12"]), 5)) for k in rr]
    vi = [v for v in vals if abs(sig(v[1])) < 3]
    say(f"   Alternative R1 mit {nm}: Lösungen {vals}; im 3-sigma-Bereich: {vi}")
    ZAHLEN[f"R1_{nm}"] = vals
    alternative("R1", nm, f"im 3-sigma-Bereich: {vi}")
    if nm.startswith("M2"):
        check("M2: Q ist dort genau der gamma-Term (Lösungen kappa = gamma* und Spiegelwert -3/8 - gamma*), beide mit s12 = TM1 = "
              "0.317563 -- keine Brechung", len(vi) == 2 and all(abs(v[1] - tm1) < 1e-5 for v in vi)
              and sorted(abs(v[0] - g) < 1e-4 for v in vi for g in (g_star, g_mirror)).count(True) == 2)
befund("(b) Unter W1 bricht Q TM1 in der TM2-Richtung (u -> w). Mit gamma = 0 aus (a) und einem angepassten kappa entsteht die "
       "Summenregel s12 = 0.3136 bei s13 = 0.023108: richtige Richtung, ein Drittel der eingefrorenen Brechung. Mit freiem "
       "gamma ist alles Anpassung", True)

say("\n   (c) Massenmuster der Nullfamilie gegen +-lambda")
need = {"NO (m1 = 0)": col_fr[0], "IO (m3 = 0)": s13f}


def zero_state_ok(col):
    return {"NO": abs(col[0] - need["NO (m1 = 0)"]) < 0.05, "IO": col[0] < 0.05 and abs(col[1] - col[2]) < 0.3}


tests = []
for perm in itertools.permutations(range(3)):                  # W0: Charakterfamilien (z, 1, z^2) -> (e, mu, tau)
    col = np.zeros(3)
    col[perm[0]] = 1.0                                         # Nullfamilie z ist eine reine Flavourrichtung
    tests.append(("W0 " + str(perm), col, zero_state_ok(col)))
for ax in range(3):
    ev_ = np.eye(3)[ax]
    col = abs(Fn.conj().T@ev_)**2
    tests.append((f"W1 Achse e{5 + ax}", col, zero_state_ok(col)))
colM2 = abs(Fn.conj().T@e5n)**2
tests.append(("M2", colM2, zero_state_ok(colM2)))
for nm, col, ok in tests:
    say(f"   R3 {nm}: Nullfamilie hat Flavourgewichte {np.round(col, 4).tolist()}; NO {ok['NO']}, IO {ok['IO']}")
    alternative("R3", nm, f"Nullfamilie {np.round(col, 3).tolist()}, NO {ok['NO']}, IO {ok['IO']}")
r3_aus = not any(ok["NO"] or ok["IO"] for _, _, ok in tests) and len(tests) == 10
check("R3: Als Massenoperator passt Q in keinem der 10 Wörterbücher: Die Nullfamilie ist entweder eine reine Flavourrichtung "
      "(W0) oder trimaximal 1/3 (jede Gitterabbildung, M2); gebraucht wären |U_e|^2 = 0.677 (NO, m1 = 0) oder 0.023 (IO, m3 = 0)",
      r3_aus)
Tson = sp.diag(1, sp.sqrt(2), 1)                                      # SU(2)-orthonormale Koordinaten -> P-Koeffizienten
Mmaj = sp.simplify(Tson.T*Bf*Q6*Tson)                                  # Majorana-Lesart: Psi(phi, phi') = B(phi, Q phi')
sv_maj = np.linalg.svd(np.array(Mmaj.evalf(), complex), compute_uv=False)
sv_op = np.sort(abs(np.array([complex(v) for v in Q6.eigenvals()])))[::-1]
say(f"   Majorana-Lesart M = B Q (symmetrisch) in SU(2)-orthonormalen Koordinaten: {Mmaj.tolist()}, Takagi-Werte "
    f"{np.round(sv_maj, 9).tolist()}; Operator-Lesart |Eigenwerte| {np.round(sv_op, 9).tolist()}")
check("in beiden Lesarten (Operator Q, Majorana-Form B Q) sind die Massen (24, 24, 0): Das Paar ist entartet (pseudo-Dirac, Typ "
      "L_e - L_mu - L_tau, umgekehrte Ordnung mit m3 = 0), Delta m^2 im Paar ist exakt null; der Nullzustand ist in beiden die "
      "Familie z, die Mischung im entarteten Paar ist unbestimmt",
      Mmaj == Mmaj.T and np.allclose(sv_maj, [24, 24, 0]) and np.allclose(sv_op, [24, 24, 0]) and Mmaj[1, :] == sp.zeros(1, 3))
evR1 = np.sort(abs(M1["ev"]))
ratio = (evR1[1]**2 - evR1[0]**2)/(evR1[2]**2 - evR1[0]**2)
say(f"   Nebenbefund: Die Eigenwerte von H (R1) als Massen gelesen gäben Delta m^2_klein/Delta m^2_groß = {ratio:.3f}, weit weg "
    f"von der beobachteten Kleinheit (Größenordnung 10^-2); z57 liest sie nicht als Massen")
befund("(c) Q gibt kein brauchbares Massenmuster: In W1 ist die Nullfamilie e5 die trimaximale TBM-Spalte nu2, die nie masselos "
       "sein kann (m2 > m1 >= 0), und Q hat in der TBM-Basis die Diagonale null (keine Aufspaltung erster Ordnung)", r3_aus)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Bilanz")
ZAHLEN["alternativen_anzahl"] = len(ALTERNATIVEN)
ZAHLEN["normierungen_geprueft"] = 3
ZAHLEN["gitterabbildungen_geprueft"] = 48
say(f"   Look-elsewhere: {len(ALTERNATIVEN)} Konfigurationen aus Lesart und Wörterbuch ausgewertet, dazu 3 Normierungen und 48 "
    f"Gitterabbildungen (6 verschiedene Q, J-verträglich 8 mit beiden Vorzeichen)")
for a in ALTERNATIVEN:
    say(f"     {a['lesart']} / {a['woerterbuch']}: {a['ergebnis']}")
tafel = [
    ("E", "Q exakt; tr Q^2 = 1152 I, det Q = 27648 J; bei J = 0 Eigenwerte (lambda, -lambda, 0); die Nullfamilie ist durch die "
          "Vierteldrehung erzwungen (Charakter -1, P = z); die Eigenfamilien sind die V4-Basis = die drei Markenpaarungen"),
    ("E unter W1", "Q diagonal in den E8-Familienkoordinaten, null auf der mu4-Achse e5; Q F_J = 0, Q^2 = I - F_J, Q u = w; in der "
                   "zyklischen Basis Q = (i/sqrt3)(S - S^T)"),
    ("C", "(a) linear gamma = 0; Q^2 = I - F_J ist ein gamma-Term freier Größe (passt nur angepasst, auf dem Spiegelast); ein "
          "z^2-Anteil bräuchte J != 0: Q legt gamma nicht fest"),
    ("L", "(b) R1: gamma = 0, kappa an s13 angepasst: s12 = 0.3136 (NuFIT 6.1 +0.7 sigma, TM1 +1.3 sigma, eingefroren -0.3 sigma); "
          "richtige Richtung, 36 % der eingefrorenen Brechung"),
    ("E negativ", "(c) Q als Massenoperator in allen 10 Wörterbüchern ausgeschlossen; entartetes Paar"),
    ("O", "relative Skala kappa zwischen Weylmatrix und Coxeter-Clock (keine Normregel), Phase, das Wörterbuch selbst (z57-Raum "
          "= H6-Familien; zyklische Leptonbasis), delta = 240 Grad"),
]
for k_, v_ in tafel:
    say(f"   {k_}: {v_}")
ZAHLEN["bilanz"] = [list(r) for r in tafel]
urteil = ("Die Weylmatrix des H6-Spin-2-Zustands auf den Familien ist exakt: Bei J = 0 hat sie die Eigenwerte (lambda, -lambda, 0), "
          "die Nullfamilie ist durch die Vierteldrehung erzwungen (Charakter -1), und die Eigenfamilien sind die drei "
          "Markenpaarungen (V4-Basis). Unter jeder Gitterabbildung A3 = D3 ist Q diagonal in den E8-Familienkoordinaten; im "
          "J-verträglichen Wörterbuch gilt Q F_J = 0, Q^2 = I - F_J und Q u = w. Damit legt der Zustand gamma nicht fest (linear "
          "0, quadratisch ein gamma-Term freier Größe, ein z^2-Anteil bräuchte J != 0); linear gibt er nur eine TM2-artige "
          "Kopplung der TM1-Spalte. "
          "Mit gamma = 0 und einem angepassten Faktor kappa = 0.0499 folgt die Summenregel s12 = 0.3136 bei s13 = 0.023108: "
          "zwischen TM1 (0.3176) und dem eingefrorenen Wert (0.3067), NuFIT 6.1 +0.7 sigma, aber nur 36 % der eingefrorenen "
          "Brechung und eine Lesart unter 18 ausgewerteten. Als Massenoperator ist Q ausgeschlossen. Nichts davon reproduziert die "
          "eingefrorenen Werte ohne Anpassung; die Idee (E) gibt die Familiengewichte nicht.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z67_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "protokoll": PROTOKOLL, "alternativen": ALTERNATIVEN,
     "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
