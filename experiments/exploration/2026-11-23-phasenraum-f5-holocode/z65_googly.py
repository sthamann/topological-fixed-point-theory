"""Z65: Lücke (C) -- die nicht-selbstduale Hälfte der Gravitation (googly-Problem) und die zwei Blätter der Naht.

Frage: Können TFPTs zwei Blätter (holomorphes Blatt und Konjugierten-/Purifikationsblatt, Tomita-Konjugation J) die beiden
Helizitäten des Gravitons tragen, und wie weit reicht das?

Konventionen wie z61 (Penrose-Rindler): Twistor Z = (omega^A, pi_A'), dualer Twistor W = (lambda_A, mu^A'),
Z.W = omega^A lambda_A + pi_A' mu^A'. Inzidenz omega = i x pi (Gerade L_x in PT, Blatt 1), mu = -i x^T lambda (Gerade L~_x in
PT*, Blatt 2). eps_01 = eps^01 = 1. Koordinaten x^{AA'} mit Index a = 2A + A', Metrik g_ab = eps_AB eps_A'B' (ds^2 = 2 det dx).
Psi~_{A'B'C'D'} heißt wie in z61 'selbstdual' (Helizität +2), Psi_{ABCD} 'antiselbstdual' (Helizität -2).
Statuskürzel: E exakt, C bedingt, L Lesart, O offen.

Teil 1  Linear, zwei Blätter: Grad -6 auf PT gibt Psi~, Grad -6 auf PT* gibt Psi (Residuen = geschlossene Formen, Bianchi
        exakt). Grad +2 auf PT gibt ebenfalls Psi: linear trägt schon ein Blatt beide Helizitäten.
Teil 2  Eine Metrik mit beiden Hälften: h = h+ + h- aus je einem Elementarzustand pro Blatt (Spinheben des Grad-(-2)-Keims
        mit o^A = B^A bzw. o~^A' = D^A'). R^(1)_ab = 0 als Polynomidentität, Weylteile = Penrose-Transformierte beider Blätter,
        allgemein über Fouriersymbole. Jeder Elementarzustand allein ist eine exakte nichtlineare Kerr-Schild-Lösung einer
        Chiralität (Reihe in eps bis Ordnung 10, exakt an zwei komplexen Punkten).
Teil 3  Nichtlinear: Kreuzterm R^(2)[h+, h-] != 0 (nicht kollinear, auch J-reell); kollinear (pp-Welle) überlagern sich beide
        Hälften exakt (allgemeine Funktionen, symbolisch).
Teil 4  Ambitwistoren: L_x x L~_x liegt in Z.W = 0; Q_x und Q_y treffen sich genau bei Nullabstand in einem Punkt (der
        Nullgeodäte); Klein-Korrespondenz; J (Tomita, x -> x^+) ist die Realstruktur PT -> PT*; Realteil von A = PN;
        J-Fixpunkte von Q_x = Himmelssphäre = reine Zustände (z62).
Teil 5  J und Helizität: J (Realstruktur) bildet h+ auf seinen Realitätspartner ab (anderer Tubus), J mit Spiegelung
        (CPT, Bisognano-Wichmann-Form) auf ein h- gleicher Frequenz (Helizität -2 auf Blatt 2); J-reelle Störungen tragen
        stets beide Hälften; alpha/beta-Hindernis: h+ hindert Blatt 1 und deformiert Blatt 2, h- umgekehrt
        (googly-Dualität), linear und nichtlinear.
Teil 6  Bilanz: bekannt (Penrose, Ward, Witten, Isenberg-Yasskin-Green, LeBrun, Baston-Mason, Mason-Skinner, Adamo-Mason),
        was TFPT bräuchte, was offen bleibt.
"""
import itertools
import json
import random
from pathlib import Path

import numpy as np
import sympy as sp
from sympy import QQ_I

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")
    return bool(cond)


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ Grundlagen
I = sp.I
X = sp.Matrix(2, 2, sp.symbols("x00 x01 x10 x11"))
XS = [X[0, 0], X[0, 1], X[1, 0], X[1, 1]]                  # Koordinatenindex a = 2A + A'
EPS = ((0, 1), (-1, 0))                                    # eps_AB = eps^AB numerisch, eps_01 = 1
IDX = [divmod(a, 2) for a in range(4)]                     # a -> (A, A')
SIG = [2*Ap + A for (A, Ap) in IDX]                        # Konjugation tauscht ungestrichen <-> gestrichen
GINV = [[EPS[A][C]*EPS[Ap][Cp] for (C, Cp) in IDX] for (A, Ap) in IDX]   # g^ab = g_ab = eps eps
z = sp.Symbol("z")


def lo(k):                                                 # kappa_B = kappa^A eps_AB
    return [sp.expand(sum(k[A]*EPS[A][B] for A in range(2))) for B in range(2)]


def up(k):                                                 # kappa^A = eps^AB kappa_B
    return [sp.expand(sum(EPS[A][B]*k[B] for B in range(2))) for A in range(2)]


def dot2(u, v):
    return sp.expand(u[0]*v[0] + u[1]*v[1])


random.seed(65)


def rq(n=9, m=5):
    return sp.Rational(random.randint(-n, n), random.randint(1, m))


def cq():
    while True:
        w = rq() + I*rq()
        if w != 0:
            return w


def herm():
    t, s, u, v = (rq() for _ in range(4))
    return sp.Matrix([[t + s, u - I*v], [u + I*v, t - s]])


def ebene(lam, Y):            # dualer Twistor (lam_A, mu^A') mit mu^A' = -i lam_A Y^AA': die Ebene enthält L_Y
    return list(lam), [sp.expand(-I*sum(lam[A]*Y[A, k] for A in range(2))) for k in range(2)]


def punkt(pil, Y):            # Twistor (omega^A, pi_A') auf L_Y: omega^A = i Y^AA' pi_A'
    return [sp.expand(I*sum(Y[A, k]*pil[k] for k in range(2))) for A in range(2)], list(pil)


def auf_L(dt):                # A.Z auf L_x = a^A'(x) pi_A' mit a^A' = i lam_A x^AA' + mu^A'
    lam, mu = dt
    return [sp.expand(I*sum(lam[A]*X[A, k] for A in range(2)) + mu[k]) for k in range(2)]


def auf_Lt(tw):               # C.W auf L~_x = c^A(x) lambda_A mit c^A = C^A - i x^AA' C_A'
    om, pil = tw
    return [sp.expand(om[A] - I*sum(X[A, k]*pil[k] for k in range(2))) for A in range(2)]


def residue(numer, la, lb, order):
    """Residuum bei der Nullstelle von la.(1,z) von numer/((la.(1,z))^order (lb.(1,z))) (wie z61)."""
    zp = -la[0]/la[1]
    g = numer/((lb[0] + lb[1]*z)*la[1]**order)
    return sp.diff(g, z, order - 1).subs(z, zp)/sp.factorial(order - 1)


def P(e):
    return sp.Poly(sp.expand(e), *XS, domain=QQ_I)


def ev(poly, pt):
    return QQ_I.from_sympy(sp.expand(poly(*pt)))


def qi(e):
    return QQ_I.from_sympy(sp.expand(e))


def cconj(w):
    return QQ_I(w.x, -w.y)


class Keim:
    """Skalarer Keim chi = 1/Q (Q quadratisch); d_idx chi = N_idx / Q^(|idx|+1) mit Polynom N_idx (exakt, QQ_I)."""

    def __init__(self, Q):
        self.Q = P(Q)
        self.dQ = [self.Q.diff(v) for v in XS]
        self.cache = {(): P(1)}

    def N(self, idx):
        idx = tuple(sorted(idx))
        if idx not in self.cache:
            prev, last, m = self.N(idx[:-1]), idx[-1], len(idx)
            self.cache[idx] = prev.diff(XS[last])*self.Q - m*prev*self.dQ[last]
        return self.cache[idx]


def h_num(art, keim, s_lo, s_up, a, b, extra=()):
    """Zähler von d_extra h_ab über Q^(3 + |extra|).
    art '+': h+_ab = o_A o_B o^C o^D d_(C A') d_(D B') chi   (Blatt 1, Spinheben mit o^A)
    art '-': h-_ab = o~_A' o~_B' o~^C' o~^D' d_(A C') d_(B D') chi~   (Blatt 2, Spinheben mit o~^A')"""
    (A, Ap), (B, Bp) = IDX[a], IDX[b]
    tot = P(0)
    for C in range(2):
        for D in range(2):
            if art == "+":
                pref, idx = s_lo[A]*s_lo[B]*s_up[C]*s_up[D], (2*C + Ap, 2*D + Bp)
            else:
                pref, idx = s_lo[Ap]*s_lo[Bp]*s_up[C]*s_up[D], (2*A + C, 2*B + D)
            if pref != 0:
                tot += keim.N(idx + tuple(extra))*P(pref)
    return tot


def riemann_lin(art, keim, s_lo, s_up):
    """Linearisierter Riemann R_abcd = 1/2 (d_b d_c h_ad + d_a d_d h_bc - d_a d_c h_bd - d_b d_d h_ac), Zähler über Q^5."""
    cache = {}

    def hh(p, q, e, f):
        key = (min(p, q), max(p, q), min(e, f), max(e, f))
        if key not in cache:
            cache[key] = h_num(art, keim, s_lo, s_up, key[0], key[1], (key[2], key[3]))
        return cache[key]
    half = P(sp.Rational(1, 2))
    return {(a, b, c, d): (hh(a, d, b, c) + hh(b, c, a, d) - hh(b, d, a, c) - hh(a, c, b, d))*half
            for a, b, c, d in itertools.product(range(4), repeat=4)}


def weyl_aus_riemann(R, scale=None):
    """Psi_ABCD = 1/4 R_{AA'BB'CC'DD'} eps^A'B' eps^C'D' und Psi~_A'B'C'D' = 1/4 R_{AA'BB'CC'DD'} eps^AB eps^CD.
    R als dict (a,b,c,d) -> Wert; die freien Indizes heißen hier (F1..F4), die summierten (s1..s4)."""
    Psi, Psit = {}, {}
    for F in itertools.product(range(2), repeat=4):
        s_u = s_g = None
        for s in itertools.product(range(2), repeat=4):
            e = EPS[s[0]][s[1]]*EPS[s[2]][s[3]]
            if not e:
                continue
            t_u = R[tuple(2*F[i] + s[i] for i in range(4))]*e     # frei ungestrichen, summiert gestrichen -> Psi
            t_g = R[tuple(2*s[i] + F[i] for i in range(4))]*e     # frei gestrichen, summiert ungestrichen -> Psi~
            s_u = t_u if s_u is None else s_u + t_u
            s_g = t_g if s_g is None else s_g + t_g
        Psi[F] = s_u*scale if scale is not None else s_u
        Psit[F] = s_g*scale if scale is not None else s_g
    return Psi, Psit


def komp(n):                                               # Komponente mit n Einsen: (0,..,0,1,..,1)
    return tuple([0]*(4 - n) + [1]*n)


# ------------------------------------------------------------------ Daten der beiden Blätter
Y1 = I*sp.eye(2) + herm()        # Zentrum des Elementarzustands auf Blatt 1 (Im Y1 = 1 > 0: regulär auf dem reellen Schnitt)
Y2 = I*sp.eye(2) + herm()        # Blatt 2, gleicher Tubus (gleiche Frequenz wie Blatt 1)
A_dt, B_dt = ebene([cq(), cq()], Y1), ebene([cq(), cq()], Y1)
C_tw, D_tw = punkt([cq(), cq()], Y2), punkt([cq(), cq()], Y2)
a_, b_ = auf_L(A_dt), auf_L(B_dt)
c_, d_ = auf_Lt(C_tw), auf_Lt(D_tw)
b_lo, d_lo = lo(b_), lo(d_)
Q1 = dot2(a_, b_lo)                                        # a^A' b_A'
Q2 = dot2(c_, d_lo)                                        # c^A d_A
o_up = up(B_dt[0])                                         # o^A = B^A, also o^A B_A = 0
o_lo = lo(o_up)
ot_up = up(D_tw[1])                                        # o~^A' = D^A', also o~^A' D_A' = 0
ot_lo = lo(ot_up)
oA = dot2(o_up, A_dt[0])
otC = dot2(ot_up, C_tw[1])
ZAHLEN["daten"] = {"Y1": str(Y1), "Y2": str(Y2), "o.A": str(oA), "o~.C": str(otC)}

# ------------------------------------------------------------------ Teil 1
say("Teil 1: Zwei Blätter, linear (Penrose-Transformation auf PT und PT*)")
# Die Residuen hängen von x nur über die Koeffizienten der Linearformen a.pi, b.pi (bzw. c.lambda, d.lambda) ab. Die
# Identitäten werden daher mit generischen Symbolen (g_a, g_b, g_p) exakt geprüft und dann x-abhängig eingesetzt.
GA, GB, GP = list(sp.symbols("ga0 ga1")), list(sp.symbols("gb0 gb1")), list(sp.symbols("gp0 gp1"))
GB_lo = lo(GB)
GAB = dot2(GA, GB_lo)
r1 = sp.cancel(Q1/(X - Y1).det())
r2 = sp.cancel(Q2/(X - Y2).det())
check("E: Keime 1/((A.Z)(B.Z)) auf PT und 1/((C.W)(D.W)) auf PT* geben chi = 1/(a.b), chi~ = 1/(c.d) mit a.b = "
      f"({r1}) det(x - Y1), c.d = ({r2}) det(x - Y2): Elementarzustände mit Zentrum im Tubus Im Y = 1 (regulär auf dem "
      "reellen Schnitt); Residuum = 1/(a.b) für generische Linearformen",
      r1.free_symbols == set() and r2.free_symbols == set() and r1 != 0 and r2 != 0 and
      sp.cancel(residue(sp.Integer(1), GA, GB, 1) - 1/GAB) == 0)

ok = all(sp.cancel(residue(z**n, GA, GB, 5) - GB_lo[0]**(4 - n)*GB_lo[1]**n/GAB**5) == 0 for n in range(5))
check("E: Grad -6 auf beiden Blättern: Res[pi^4/((a.pi)^5 (b.pi))] = b_A'b_B'b_C'b_D'/(a.b)^5 für generische a, b (alle "
      "fünf Komponenten); eingesetzt gibt das Psi~_PT = b^4/(a.b)^5 (Blatt 1) und Psi_PT* = d^4/(c.d)^5 (Blatt 2)", ok)


def feldgl(Nk, keimQ, gestrichen):
    """Bianchi im Vakuum für Psi = N_n/Q^5 (N_n Zähler der Komponente mit n Einsen); True, wenn alle Gleichungen exakt 0."""
    Qp = P(keimQ)
    ok_ = True
    for B in range(2):
        for m in range(4):
            if gestrichen:     # nabla^(AA') Psi~_A'... : d_(B,1') Psi~_(0',rest) - d_(B,0') Psi~_(1',rest)
                v1, v0 = XS[2*B + 1], XS[2*B]
            else:              # nabla^(AA') Psi_A... : d_(1,B') Psi_(0,rest) - d_(0,B') Psi_(1,rest)
                v1, v0 = XS[2 + B], XS[B]
            N0, N1 = Nk[m], Nk[m + 1]
            num = (N0.diff(v1)*Qp - 5*N0*Qp.diff(v1)) - (N1.diff(v0)*Qp - 5*N1*Qp.diff(v0))
            ok_ &= num.is_zero
    return ok_


N_pt = [P(b_lo[0]**(4 - n)*b_lo[1]**n) for n in range(5)]
N_ptd = [P(d_lo[0]**(4 - n)*d_lo[1]**n) for n in range(5)]
ok_pt, ok_ptd = feldgl(N_pt, Q1, True), feldgl(N_ptd, Q2, False)
check("E: nabla^(AA') Psi~_A'B'C'D' = 0 (Blatt 1, Helizität +2) und nabla^(AA') Psi_ABCD = 0 (Blatt 2, Helizität -2) als "
      "Polynomidentitäten (je 8 Gleichungen): beide Blätter geben exakt die linearisierte Vakuum-Bianchi-Gleichung",
      ok_pt and ok_ptd)

# Grad +2 auf PT: Psi_ABCD = Res[d^4 F/d omega^4], F = (P.Z)^4/((A.Z)(B.Z))
sPs, sAs, sBs = sp.symbols("sP sA sB")
P_dt = ([cq(), cq()], [cq(), cq()])
p_ = auf_L(P_dt)


def Dw(e, A):
    return P_dt[0][A]*sp.diff(e, sPs) + A_dt[0][A]*sp.diff(e, sAs) + B_dt[0][A]*sp.diff(e, sBs)


def br(u, v):                                              # <u v> = u^A' v_A'
    return dot2(u, lo(v))


def res_formel(m, j, k, aa, bb, pp):
    """Res_(a-Pol) (p.pi)^m/((a.pi)^j (b.pi)^k), m = j + k - 2: (-1)^(j-1) C(j+k-2, j-1) <ap>^(k-1) <bp>^(j-1)/<ab>^(j+k-1);
    k = 0 gibt 0 (nur ein Pol). Rückgabe: (Zähler, Potenz von <ab>)."""
    if k == 0 or j <= 0:
        return sp.Integer(0), 0
    return (-1)**(j - 1)*sp.binomial(j + k - 2, j - 1)*br(aa, pp)**(k - 1)*br(bb, pp)**(j - 1), j + k - 1


ok_rf = True
for j in range(1, 6):
    for k in range(1, 7 - j):
        m = j + k - 2
        g_ = (GP[0] + GP[1]*z)**m/((GB[0] + GB[1]*z)**k*GA[1]**j)
        r_ = sp.diff(g_, z, j - 1).subs(z, -GA[0]/GA[1])/sp.factorial(j - 1)
        nf, pw = res_formel(m, j, k, GA, GB, GP)
        ok_rf &= sp.cancel(r_ - nf/br(GA, GB)**pw) == 0
N_g2 = []
for n in range(5):
    e = sPs**4/(sAs*sBs)
    for A in komp(n):
        e = Dw(e, A)
    tot = P(0)
    for term in sp.Add.make_args(sp.expand(e)):
        pw_ = term.as_powers_dict()
        j, k, m = -pw_.get(sAs, 0), -pw_.get(sBs, 0), pw_.get(sPs, 0)
        nf, pw = res_formel(m, j, k, a_, b_, p_)
        if pw:
            tot += P(term.subs({sPs: 1, sAs: 1, sBs: 1})*nf)*P(Q1)**(5 - pw)
    N_g2.append(tot)                                       # Psi_(n Einsen) = N/(a.b)^5
ok_g2 = ok_rf and feldgl(N_g2, Q1, False) and any(not N.is_zero for N in N_g2)
check("E: Grad +2 auf PT (F = (P.Z)^4/((A.Z)(B.Z)), Psi_ABCD = Res d^4F/d omega^4; Residuenformel für alle 15 Polordnungen "
      "symbolisch geprüft) gibt Psi = N/(a.b)^5 mit nabla^(AA') Psi_ABCD = 0 exakt und Psi != 0: dasselbe Blatt trägt "
      "linear auch die Helizität -2", ok_g2)
lin_beide = CHECKS[list(CHECKS)[-1]] and ok_pt
befund("E: Linear trägt schon ein Blatt beide Helizitäten (Grad -6 -> Psi~, Grad +2 -> Psi); das googly-Problem ist kein "
       "lineares Problem. Penroses Deformationen von PT sind die Grad-+2-Klassen H1(PT, O(2)), also Psi (Helizität -2), "
       "nicht das Grad-(-6)-Feld Psi~ aus z61", lin_beide)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Eine Metrik mit beiden Hälften, je ein Elementarzustand pro Blatt")
K1, K2 = Keim(Q1), Keim(Q2)
ok = sp.expand(dot2(o_up, B_dt[0])) == 0 and sp.expand(dot2(ot_up, D_tw[1])) == 0 and oA != 0 and otC != 0
for n in range(5):
    pr = komp(n)
    phi1 = P(0)
    phi2 = P(0)
    for A, B, C, D in itertools.product(range(2), repeat=4):
        w1 = o_up[A]*o_up[B]*o_up[C]*o_up[D]
        if w1 != 0:
            phi1 += K1.N((2*A + pr[0], 2*B + pr[1], 2*C + pr[2], 2*D + pr[3]))*P(w1)
        w2 = ot_up[A]*ot_up[B]*ot_up[C]*ot_up[D]
        if w2 != 0:
            phi2 += K2.N((2*pr[0] + A, 2*pr[1] + B, 2*pr[2] + C, 2*pr[3] + D))*P(w2)
    ok &= (phi1 - P(24*oA**4)*N_pt[n]).is_zero and (phi2 - P(24*otC**4)*N_ptd[n]).is_zero
check("E: Spinheben ist Penrose: (o^A d_AA')^4 chi = 24 (o.A)^4 Psi~_PT und (o~^A' d_AA')^4 chi~ = 24 (o~.C)^4 Psi_PT* "
      "exakt, mit o^A = B^A, o~^A' = D^A' (o.B = 0, o~.D = 0): der Grad-(-6)-Elementarzustand ist das Spinheben des "
      "Grad-(-2)-Keims", ok)

Rp = riemann_lin("+", K1, o_lo, o_up)
Rm = riemann_lin("-", K2, ot_lo, ot_up)
ric_p = [sum((Rp[(a, b, c, d)]*P(GINV[a][c]) for a in range(4) for c in range(4) if GINV[a][c]), P(0))
         for b in range(4) for d in range(4)]
ric_m = [sum((Rm[(a, b, c, d)]*P(GINV[a][c]) for a in range(4) for c in range(4) if GINV[a][c]), P(0))
         for b in range(4) for d in range(4)]
hp0 = {(a, b): h_num("+", K1, o_lo, o_up, a, b) for a in range(4) for b in range(4)}
hm0 = {(a, b): h_num("-", K2, ot_lo, ot_up, a, b) for a in range(4) for b in range(4)}
sym_ok = all((hp0[(a, b)] - hp0[(b, a)]).is_zero and (hm0[(a, b)] - hm0[(b, a)]).is_zero for a in range(4) for b in range(4))
tr_ok = sum((hp0[(a, b)]*P(GINV[a][b]) for a in range(4) for b in range(4) if GINV[a][b]), P(0)).is_zero and \
    sum((hm0[(a, b)]*P(GINV[a][b]) for a in range(4) for b in range(4) if GINV[a][b]), P(0)).is_zero
div_ok = True
for b in range(4):
    dp = sum((h_num("+", K1, o_lo, o_up, a, b, (c,))*P(GINV[a][c]) for a in range(4) for c in range(4) if GINV[a][c]), P(0))
    dm = sum((h_num("-", K2, ot_lo, ot_up, a, b, (c,))*P(GINV[a][c]) for a in range(4) for c in range(4) if GINV[a][c]),
             P(0))
    div_ok &= dp.is_zero and dm.is_zero
nonzero_h = any(not v.is_zero for v in hp0.values()) and any(not v.is_zero for v in hm0.values())
check("E: h = h+ + h- (h+ = o_A o_B o^C o^D d_CA' d_DB' chi, h- = o~_A' o~_B' o~^C' o~^D' d_AC' d_BD' chi~) ist symmetrisch, "
      "spurfrei, transversal, != 0 und löst R^(1)_ab = 0 exakt (alle 16 Komponenten als Polynomidentität über Q^5)",
      sym_ok and tr_ok and div_ok and nonzero_h and all(r.is_zero for r in ric_p) and all(r.is_zero for r in ric_m))

Psi_p, Psit_p = weyl_aus_riemann(Rp, scale=P(sp.Rational(1, 4)))
Psi_m, Psit_m = weyl_aus_riemann(Rm, scale=P(sp.Rational(1, 4)))
ok_w = any(not v.is_zero for v in Psit_p.values()) and any(not v.is_zero for v in Psi_m.values())
for key in itertools.product(range(2), repeat=4):
    n = sum(key)
    ok_w &= (Psit_p[key] - P(-12*oA**4)*N_pt[n]).is_zero and Psi_p[key].is_zero
    ok_w &= (Psi_m[key] - P(-12*otC**4)*N_ptd[n]).is_zero and Psit_m[key].is_zero
check("E: Weylteile von h: Psi~[h] = -12 (o.A)^4 Psi~_PT (nur aus Blatt 1), Psi[h] = -12 (o~.C)^4 Psi_PT* (nur aus Blatt 2), "
      "beide != 0, exakt in allen 16 Komponenten: h enthält beide Helizitäten, jede aus ihrem Blatt", ok_w)
ZAHLEN["weyl_konstanten"] = {"Psi~[h]/Psi~_PT": str(sp.expand(-12*oA**4)), "Psi[h]/Psi_PT*": str(sp.expand(-12*otC**4))}

# allgemein: Fouriersymbole (chi = e^{i k.x}, beliebige konstante o, o~)
Kk = sp.Matrix(2, 2, sp.symbols("k00 k01 k10 k11"))
KS = [Kk[0, 0], Kk[0, 1], Kk[1, 0], Kk[1, 1]]
os_up, ots_up = list(sp.symbols("o0 o1")), list(sp.symbols("p0 p1"))
os_lo, ots_lo = lo(os_up), lo(ots_up)
ka0, ka1, kt0, kt1 = sp.symbols("kappa0 kappa1 kappat0 kappat1")
kap, kapt = [ka0, ka1], [kt0, kt1]


def symbol_h(art):
    h = {}
    for a in range(4):
        for b in range(4):
            (A, Ap), (B, Bp) = IDX[a], IDX[b]
            if art == "+":
                h[(a, b)] = os_lo[A]*os_lo[B]*(I*sum(os_up[C]*KS[2*C + Ap] for C in range(2)))*(
                    I*sum(os_up[D]*KS[2*D + Bp] for D in range(2)))
            else:
                h[(a, b)] = ots_lo[Ap]*ots_lo[Bp]*(I*sum(ots_up[C]*KS[2*A + C] for C in range(2)))*(
                    I*sum(ots_up[D]*KS[2*B + D] for D in range(2)))
    R = {}
    for a, b, c, d in itertools.product(range(4), repeat=4):
        R[(a, b, c, d)] = sp.Rational(1, 2)*(-KS[b]*KS[c]*h[(a, d)] - KS[a]*KS[d]*h[(b, c)] + KS[a]*KS[c]*h[(b, d)]
                                             + KS[b]*KS[d]*h[(a, c)])
    return h, R


def mod_det(e):                       # Rest modulo det k (k00 = k01 k10/k11 einsetzen; det k ist irreduzibel)
    r = sp.together(sp.expand(e).subs(KS[0], KS[1]*KS[2]/KS[3]))
    return sp.expand(sp.numer(r)) == 0


ok_s = True
for art in "+-":
    hs, Rs = symbol_h(art)
    for b in range(4):
        for d in range(4):
            ok_s &= mod_det(sum(GINV[a][c]*Rs[(a, b, c, d)] for a in range(4) for c in range(4)))
    Ps, Pts = weyl_aus_riemann(Rs, scale=sp.Rational(1, 4))
    kegel = {KS[2*A + Ap]: kap[A]*kapt[Ap] for (A, Ap) in IDX}
    ok_n = dot2(os_up, kap)
    ot_n = dot2(ots_up, kapt)
    for key in itertools.product(range(2), repeat=4):
        if art == "+":
            ok_s &= mod_det(Ps[key])
            ok_s &= sp.expand(Pts[key].subs(kegel) + sp.Rational(1, 2)*ok_n**4*sp.prod([kapt[i] for i in key])) == 0
        else:
            ok_s &= mod_det(Pts[key])
            ok_s &= sp.expand(Ps[key].subs(kegel) + sp.Rational(1, 2)*ot_n**4*sp.prod([kap[i] for i in key])) == 0
check("E: allgemein (Fouriersymbole, beliebige Lösungen chi, chi~ von Box = 0, beliebige konstante o, o~): R^(1)[h+ + h-] "
      "= 0 und Psi[h+] = 0 = Psi~[h-] modulo det k; auf dem Kegel k = kappa kappa~ ist Psi~[h+] = -1/2 (o.kappa)^4 kappa~^4 "
      "und Psi[h-] = -1/2 (o~.kappa~)^4 kappa^4: die Helizitäten trennen sich sauber nach Blättern", ok_s)

# ---- exakte nichtlineare Krümmung an einem Punkt: Reihe in eps mit Koeffizienten in QQ_I
ZERO, ONE = QQ_I.zero, QQ_I.one
HALB, VIERTEL = qi(sp.Rational(1, 2)), qi(sp.Rational(1, 4))
ETA = np.array([[qi(GINV[a][b]) for b in range(4)] for a in range(4)], dtype=object)
ID4 = np.array([[ONE if a == b else ZERO for b in range(4)] for a in range(4)], dtype=object)


def nzarr(T):
    return any(bool(v) for v in np.asarray(T).flat)


class Reihe:
    """Abgeschnittene Potenzreihe in eps mit Tensorkoeffizienten (exakt über QQ_I)."""

    def __init__(self, terme, nmax):
        self.n = nmax
        self.t = {k: v for k, v in terme.items() if k <= nmax and nzarr(v)}

    def __add__(self, o):
        t = dict(self.t)
        for k, v in o.t.items():
            t[k] = t[k] + v if k in t else v
        return Reihe(t, min(self.n, o.n))

    def __neg__(self):
        return Reihe({k: -v for k, v in self.t.items()}, self.n)

    def __sub__(self, o):
        return self + (-o)

    def mal(self, c):
        return Reihe({k: v*c for k, v in self.t.items()}, self.n)

    def ein(self, spec, o):
        t, n = {}, min(self.n, o.n)
        for i, A in self.t.items():
            for j, B in o.t.items():
                if i + j <= n:
                    r = np.einsum(spec, A, B)
                    t[i + j] = t[i + j] + r if (i + j) in t else r
        return Reihe(t, n)

    def perm(self, spec):
        return Reihe({k: np.einsum(spec, v) for k, v in self.t.items()}, self.n)


def kruemmung(h0, h1, h2, nmax, E=None):
    """Exakte Krümmung von g = eta + eps h an einem Punkt (h, d_e h, d_e d_f h dort gegeben), als Reihe in eps bis nmax.
    R^a_bcd = d_c Gam^a_db - d_d Gam^a_cb + Gam^a_ce Gam^e_db - Gam^a_de Gam^e_cb, Ric_bd = R^a_bad, R_abcd = g_ae R^e_bcd.
    Mit Rahmen E[mu, a] zusätzlich die Rahmenkomponenten R_abcd (für die Weylspinoren)."""
    g = Reihe({0: ETA, 1: h0}, nmax)
    dg, ddg = Reihe({1: h1}, nmax), Reihe({1: h2}, nmax)
    M, T, terms = ETA.dot(h0), ETA.copy(), {0: ETA.copy()}
    for k in range(1, nmax + 1):                            # g^-1 = sum (-eps eta h)^k eta (eta^-1 = eta)
        T = -M.dot(T)
        terms[k] = T
    gi = Reihe(terms, nmax)
    dgi = -(gi.ein("ac,ecd->ead", dg).ein("ead,db->eab", gi))
    S = dg.perm("bdc->dbc") + dg.perm("cdb->dbc") - dg
    dS = ddg.perm("ebdc->edbc") + ddg.perm("ecdb->edbc") - ddg
    Gam = gi.ein("ad,dbc->abc", S).mal(HALB)
    dGam = (dgi.ein("ead,dbc->eabc", S) + gi.ein("ad,edbc->eabc", dS)).mal(HALB)
    Rup = (dGam.perm("cadb->abcd") - dGam.perm("dacb->abcd") + Gam.ein("ace,edb->abcd", Gam)
           - Gam.ein("ade,ecb->abcd", Gam))
    out = {"Ric": Rup.perm("abad->bd"), "R": g.ein("ae,ebcd->abcd", Rup)}
    if E is not None:
        T1 = out["R"].ein("mnrs,sd->mnrd", E)
        T2 = T1.ein("mnrd,rc->mncd", E)
        T3 = T2.ein("mncd,nb->mbcd", E)
        out["Rfr"] = T3.ein("mbcd,ma->abcd", E)
    return out


def weyl_np(Rarr):
    d = {key: Rarr[key] for key in itertools.product(range(4), repeat=4)}
    return weyl_aus_riemann(d, scale=VIERTEL)


def punktwerte(art, keim, s_lo, s_up, pt):
    """h, d_e h, d_e d_f h am Punkt pt (exakt, QQ_I)."""
    Qv = ev(keim.Q, pt)
    h0, h1, h2 = (np.empty((4, 4), dtype=object), np.empty((4, 4, 4), dtype=object),
                  np.empty((4, 4, 4, 4), dtype=object))
    for a in range(4):
        for b in range(4):
            h0[a, b] = ev(h_num(art, keim, s_lo, s_up, a, b), pt)/Qv**3
            for e in range(4):
                h1[e, a, b] = ev(h_num(art, keim, s_lo, s_up, a, b, (e,)), pt)/Qv**4
                for f in range(4):
                    h2[e, f, a, b] = ev(h_num(art, keim, s_lo, s_up, a, b, (e, f)), pt)/Qv**5
    return h0, h1, h2


def rpunkt():
    return tuple(rq() + I*rq() for _ in range(4))


x0, x1 = rpunkt(), rpunkt()
PW = {"+": punktwerte("+", K1, o_lo, o_up, x0), "-": punktwerte("-", K2, ot_lo, ot_up, x0)}
PW1 = {"+": punktwerte("+", K1, o_lo, o_up, x1), "-": punktwerte("-", K2, ot_lo, ot_up, x1)}
ks_ok, ks_info = True, {}
for (art, keim, s_lo, s_up, lvec, Vfak, Qx), (xp, PWp) in itertools.product(
        (("+", K1, o_lo, o_up, (o_lo, b_lo), -2*oA**2, Q1), ("-", K2, ot_lo, ot_up, (d_lo, ot_lo), -2*otC**2, Q2)),
        ((x0, PW), (x1, PW1))):
    # Kerr-Schild-Form: h+ = V l l mit l_AA' = o_A b_A'(x), h- = V~ l~ l~ mit l~_AA' = d_A(x) o~_A'
    l_sym = [sp.expand(lvec[0][A]*lvec[1][Ap]) for (A, Ap) in IDX]
    l_up_sym = [sp.expand(sum(GINV[a][b]*l_sym[b] for b in range(4))) for a in range(4)]
    geod = all(sp.expand(sum(l_up_sym[a]*sp.diff(l_sym[b], XS[a]) for a in range(4))) == 0 for b in range(4))
    null = sp.expand(sum(l_up_sym[a]*l_sym[a] for a in range(4))) == 0
    x0_, (h0, h1, h2) = xp, PWp[art]
    Qv = ev(P(Qx), x0_)
    lv = np.array([qi(sp.expand(l_sym[a].subs(dict(zip(XS, x0_))))) for a in range(4)], dtype=object)
    lu = ETA.dot(lv)
    V = qi(Vfak)/Qv**3
    ks_form = all(not (h0[a, b] - V*lv[a]*lv[b]) for a in range(4) for b in range(4))
    E = Reihe({0: ID4, 1: -np.outer(lu, lv)*(V*HALB)}, 10)          # E^mu_a = delta - eps V l^mu l_a / 2
    K = kruemmung(h0, h1, h2, 10, E)
    ric0 = not K["Ric"].t
    Psi_o, Psit_o = {}, {}
    for k_, arr in K["Rfr"].t.items():
        ps, pts = weyl_np(arr)
        if any(bool(v) for v in ps.values()):
            Psi_o[k_] = ps
        if any(bool(v) for v in pts.values()):
            Psit_o[k_] = pts
    if art == "+":
        ok_art = geod and null and ks_form and ric0 and not Psi_o and set(Psit_o) == {1}
        ks_info["h+"] = {"Ordnungen Psi~": sorted(Psit_o), "Ordnungen Psi": sorted(Psi_o)}
        lin_ok = all(not (Psit_o[1][key] - qi(-12*oA**4)*ev(N_pt[sum(key)], x0_)/Qv**5) for key in Psit_o[1])
    else:
        ok_art = geod and null and ks_form and ric0 and not Psit_o and set(Psi_o) == {1}
        ks_info["h-"] = {"Ordnungen Psi": sorted(Psi_o), "Ordnungen Psi~": sorted(Psit_o)}
        lin_ok = all(not (Psi_o[1][key] - qi(-12*otC**4)*ev(N_ptd[sum(key)], x0_)/Qv**5) for key in Psi_o[1])
    ks_ok &= ok_art and lin_ok
ZAHLEN["kerr_schild"] = ks_info
check("E: jeder Elementarzustand allein ist exakt nichtlinear: h+ = V l l mit l_AA' = o_A b_A'(x) (null, geodätisch), "
      "g = eta + eps h+ hat Ric = 0 in allen Ordnungen von eps (bis 10, exakt an zwei komplexen Punkten), Psi = 0 in allen "
      "Ordnungen und Psi~ nur in Ordnung 1 (= linearer Wert); h- spiegelbildlich (Psi~ = 0, Psi nur linear): "
      "Kerr-Schild/Plebanski-Gravitonen je einer Chiralität (Gürses-Gürsey: geodätisches l macht R^a_b linear)", ks_ok)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Nichtlinear: Überlagerung der beiden Hälften")
# Werkzeugprobe: Reihenkrümmung gegen direkte symbolische Rechnung für ein generisches (nicht Vakuum-) polynomiales h
epsS = sp.Symbol("eps")
random.seed(651)
hW = sp.zeros(4, 4)
for a in range(4):
    for b in range(a, 4):
        hW[a, b] = hW[b, a] = sum(rq(3, 3)*mon for mon in (1, XS[0], XS[1]*XS[2], XS[3]**2, I*XS[2]))
etaS = sp.Matrix(4, 4, lambda a, b: GINV[a][b])
GW = etaS + epsS*hW
GWi = etaS - epsS*etaS*hW*etaS                              # g^-1 bis eps^1 genügt für Ric bis eps^2
GamW = [[[sp.expand(sum(GWi[a, d]*(sp.diff(GW[d, c], XS[b]) + sp.diff(GW[d, b], XS[c]) - sp.diff(GW[b, c], XS[d]))
                        for d in range(4))/2) for c in range(4)] for b in range(4)] for a in range(4)]
ptW = (sp.Rational(1, 2) + I, sp.Rational(-2, 3), I/3, sp.Rational(5, 4) - I)
subW = dict(zip(XS, ptW))
h0W = np.array([[qi(hW[a, b].subs(subW)) for b in range(4)] for a in range(4)], dtype=object)
h1W = np.array([[[qi(sp.diff(hW[a, b], XS[e]).subs(subW)) for b in range(4)] for a in range(4)] for e in range(4)],
               dtype=object)
h2W = np.array([[[[qi(sp.diff(hW[a, b], XS[e], XS[f]).subs(subW)) for b in range(4)] for a in range(4)] for f in range(4)]
                for e in range(4)], dtype=object)
KW = kruemmung(h0W, h1W, h2W, 2)
ok_w = True
for b in range(4):
    for d in range(4):
        e_ = sum(sp.diff(GamW[a][d][b], XS[a]) - sp.diff(GamW[a][a][b], XS[d])
                 + sum(GamW[a][a][k]*GamW[k][d][b] - GamW[a][d][k]*GamW[k][a][b] for k in range(4)) for a in range(4))
        e_ = sp.expand(e_.subs(subW))
        for k in range(3):
            ok_w &= qi(e_.coeff(epsS, k)) == KW["Ric"].t.get(k, np.full((4, 4), ZERO, dtype=object))[b, d]
check("E (Werkzeug): die exakte Reihenkrümmung stimmt für ein generisches polynomiales h (kein Vakuum) in den Ordnungen "
      "eps^0, eps^1, eps^2 mit der direkten symbolischen Ricci-Rechnung überein",
      ok_w and 1 in KW["Ric"].t and 2 in KW["Ric"].t)
hp_x0, hm_x0 = PW["+"], PW["-"]
Ksum = kruemmung(hp_x0[0] + hm_x0[0], hp_x0[1] + hm_x0[1], hp_x0[2] + hm_x0[2], 2)
R1sum = Ksum["Ric"].t.get(1)
R2sum = Ksum["Ric"].t.get(2)
R2max = max(abs(complex(float(w.x), float(w.y))) for w in R2sum.flat) if R2sum is not None else 0.0
ZAHLEN["kreuzterm_R2_max"] = R2max
check("E: g = eta + eps (h+ + h-): Ric = 0 in Ordnung eps, aber Ric^(2) = Ric^(2)[h+, h-] != 0 (Kreuzterm, exakt am Punkt; "
      f"max |Ric^(2)_ab| = {R2max:.3g}); die Selbstterme verschwinden (Kerr-Schild): die beiden Hälften wechselwirken, "
      "ihre Summe ist keine Lösung", R1sum is None and R2sum is not None and nzarr(R2sum))

u_, ze_, zt_, v_ = XS
ff, gg = sp.Function("f")(u_, ze_), sp.Function("g")(u_, zt_)
Hpp = ff + gg
Gpp = sp.Matrix(4, 4, lambda a, b: GINV[a][b])
Gpp[0, 0] += Hpp
Gpi = Gpp.inv()
Gam_s = [[[sp.expand(sum(Gpi[a, d]*(sp.diff(Gpp[d, c], XS[b]) + sp.diff(Gpp[d, b], XS[c]) - sp.diff(Gpp[b, c], XS[d]))
                         for d in range(4))/2) for c in range(4)] for b in range(4)] for a in range(4)]
Rup_s = {}
for a, b, c, d in itertools.product(range(4), repeat=4):
    Rup_s[(a, b, c, d)] = sp.expand(sp.diff(Gam_s[a][d][b], XS[c]) - sp.diff(Gam_s[a][c][b], XS[d])
                                    + sum(Gam_s[a][c][e]*Gam_s[e][d][b] - Gam_s[a][d][e]*Gam_s[e][c][b] for e in range(4)))
ric_pp = [sp.simplify(sum(Rup_s[(a, b, a, d)] for a in range(4))) for b in range(4) for d in range(4)]
Rlo_s = {k: sp.expand(sum(Gpp[k[0], e]*Rup_s[(e,) + k[1:]] for e in range(4))) for k in Rup_s}
e_v = sp.eye(4)
e_v[3, 0] = Hpp/2                                          # e^(11') = dv + H du/2, sonst Koordinatenformen
Ev = e_v.inv()                                             # Ev[mu, a]
fr_ok = sp.simplify(e_v.T*sp.Matrix(4, 4, lambda a, b: GINV[a][b])*e_v - Gpp) == sp.zeros(4, 4)
nzE = {a: [(m, Ev[m, a]) for m in range(4) if Ev[m, a] != 0] for a in range(4)}
Rfr_s = {}
for a, b, c, d in itertools.product(range(4), repeat=4):
    Rfr_s[(a, b, c, d)] = sp.expand(sum(wm*wn*wr*ws*Rlo_s[(m, n, r, s)] for (m, wm) in nzE[a] for (n, wn) in nzE[b]
                                        for (r, wr) in nzE[c] for (s, ws) in nzE[d]))
Psi_pp, Psit_pp = weyl_aus_riemann(Rfr_s, scale=sp.Rational(1, 4))
hpp = {(a, b): (Hpp if (a, b) == (0, 0) else 0) for a in range(4) for b in range(4)}
Rlin_pp = {(a, b, c, d): sp.Rational(1, 2)*(sp.diff(hpp[(a, d)], XS[b], XS[c]) + sp.diff(hpp[(b, c)], XS[a], XS[d])
                                            - sp.diff(hpp[(b, d)], XS[a], XS[c]) - sp.diff(hpp[(a, c)], XS[b], XS[d]))
           for a, b, c, d in itertools.product(range(4), repeat=4)}
lin_gleich = all(sp.simplify(Rfr_s[k] - Rlin_pp[k]) == 0 for k in Rfr_s)
psi_nur_f = all(not Psi_pp[k].has(gg) for k in Psi_pp) and any(Psi_pp[k].has(ff) for k in Psi_pp)
psit_nur_g = all(not Psit_pp[k].has(ff) for k in Psit_pp) and any(Psit_pp[k].has(gg) for k in Psit_pp)
psi_nur_g = all(not Psi_pp[k].has(ff) for k in Psi_pp) and any(Psi_pp[k].has(gg) for k in Psi_pp)
psit_nur_f = all(not Psit_pp[k].has(gg) for k in Psit_pp) and any(Psit_pp[k].has(ff) for k in Psit_pp)
getrennt = (psi_nur_f and psit_nur_g) or (psi_nur_g and psit_nur_f)
traeger = {"Psi": "f(u, zeta)" if psi_nur_f else "g(u, zeta~)", "Psi~": "g(u, zeta~)" if psit_nur_g else "f(u, zeta)"}
ZAHLEN["pp_welle"] = {"Metrik": "ds^2 = 2 du dv - 2 dzeta dzeta~ + (f(u,zeta) + g(u,zeta~)) du^2", **traeger,
                      "Psi_nichtnull": {str(k): str(v) for k, v in Psi_pp.items() if v != 0},
                      "Psi~_nichtnull": {str(k): str(v) for k, v in Psit_pp.items() if v != 0}}
check("E: kollinear (pp-Welle ds^2 = 2 du dv - 2 dzeta dzeta~ + (f(u,zeta) + g(u,zeta~)) du^2, f, g beliebig): volle "
      f"nichtlineare Ric = 0 exakt, Riemann = linearer Riemann, Psi trägt nur {traeger['Psi']}, Psi~ nur {traeger['Psi~']}, "
      "beide != 0: hier überlagern sich die beiden Hälften exakt (J-reell für g = conj f, zeta~ = conj zeta)",
      fr_ok and all(r == 0 for r in ric_pp) and lin_gleich and getrennt)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Ambitwistoren -- die beiden Blätter als Quadrik Z.W = 0")
pis, lams = list(sp.symbols("pi0 pi1")), list(sp.symbols("la0 la1"))
om_s = [I*sum(X[A, k]*pis[k] for k in range(2)) for A in range(2)]
mu_s = [-I*sum(X[A, k]*lams[A] for A in range(2)) for k in range(2)]
check("E: Z.W = omega^A lambda_A + pi_A' mu^A' verschwindet identisch auf L_x x L~_x (symbolisch, alle x, pi, lambda): "
      "Q_x = L_x x L~_x = CP1 x CP1 liegt in der Ambitwistorquadrik A = {Z.W = 0} (dim A = 3 + 3 - 1 = 5 = Raum der "
      "komplexen Nullgeodäten 4 + 2 - 1)",
      sp.expand(dot2(om_s, lams) + dot2(pis, mu_s)) == 0 and 3 + 3 - 1 == 4 + 2 - 1)
rng = np.random.default_rng(65)


def crand(*s):
    return rng.normal(size=s) + 1j*rng.normal(size=s)


def Lz(x, pi_):
    return np.concatenate([1j*x @ pi_, pi_])


def Ltz(x, lam):
    return np.concatenate([lam, -1j*x.T @ lam])


def zw(Z, W):
    return Z[:2] @ W[:2] + Z[2:] @ W[2:]


def kern(Mx):
    _, s, vh = np.linalg.svd(Mx)
    return vh[-1].conj(), s


ok_null, ok_nn = True, True
for _ in range(200):
    x = crand(2, 2)
    al, be = crand(2), crand(2)
    y = x + np.outer(al, be)
    p0, s1 = kern(x - y)
    l0, s2 = kern((x - y).T)
    Z0, W0 = Lz(x, p0), Ltz(x, l0)
    ok_null &= s1[0] > 1e-6 and s1[1] < 1e-10 and s2[1] < 1e-10
    ok_null &= np.allclose(Z0, Lz(y, p0)) and np.allclose(W0, Ltz(y, l0)) and abs(zw(Z0, W0)) < 1e-10
    for t in rng.normal(size=3) + 1j*rng.normal(size=3):
        zt = x + t*(y - x)
        ok_null &= np.allclose(Z0, Lz(zt, p0)) and np.allclose(W0, Ltz(zt, l0))
    yn = crand(2, 2)
    ok_nn &= np.linalg.svd(x - yn, compute_uv=False).min() > 1e-6
check("E: Q_x und Q_y schneiden sich genau bei Nullabstand, und dann in genau einem Punkt (Z0, W0) von A (200 Nullpaare: "
      "Kerne eindimensional; 200 Zufallspaare: leer); (Z0, W0) liegt auf Q_z für alle z der Nullgeraden durch x, y: ein "
      "Punkt von A = eine Nullgeodäte", ok_null and ok_nn)


def inz_matrix(pi_, lam):          # Inzidenzsystem für vec(x) = (x00, x01, x10, x11)
    return np.array([[1j*pi_[0], 1j*pi_[1], 0, 0], [0, 0, 1j*pi_[0], 1j*pi_[1]],
                     [-1j*lam[0], 0, -1j*lam[1], 0], [0, -1j*lam[0], 0, -1j*lam[1]]])


ok_kl = True
for _ in range(200):
    xt, pi_, lam = crand(2, 2), crand(2), crand(2)
    Z, W = Lz(xt, pi_), Ltz(xt, lam)
    Mi = inz_matrix(pi_, lam)
    rhs = np.concatenate([Z[:2], W[2:]])
    sv = np.linalg.svd(Mi, compute_uv=False)
    sol = np.linalg.lstsq(Mi, rhs, rcond=None)[0]
    v, _ = kern(Mi)
    vm = v.reshape(2, 2)
    ok_kl &= sv[2] > 1e-8 and sv[3] < 1e-10 and np.allclose(Mi @ sol, rhs) and abs(zw(Z, W)) < 1e-10
    ok_kl &= abs(np.linalg.det(vm)) < 1e-10 and np.allclose(vm @ pi_, 0) and np.allclose(lam @ vm, 0)
    W2 = np.concatenate([lam, crand(2)])
    sol2 = np.linalg.lstsq(Mi, np.concatenate([Z[:2], W2[2:]]), rcond=None)[0]
    res2 = np.linalg.norm(Mi @ sol2 - np.concatenate([Z[:2], W2[2:]]))
    ok_kl &= (abs(zw(Z, W2)) > 1e-6) == (res2 > 1e-8)
check("E: Klein-Korrespondenz: das Inzidenzsystem omega = i x pi, mu = -i x^T lambda hat stets Rang 3; es ist genau dann "
      "lösbar, wenn Z.W = 0, und dann ist die Lösungsmenge eine komplexe Nullgerade mit Richtung lambda^A pi^A' (200 + 200 "
      "Stichproben)", ok_kl)

Hs = np.block([[np.zeros((2, 2)), np.eye(2)], [np.eye(2), np.zeros((2, 2))]])
ok_J = np.allclose(np.sort(np.linalg.eigvalsh(Hs)), [-1, -1, 1, 1])
for _ in range(200):
    x, pi_ = crand(2, 2), crand(2)
    Z = Lz(x, pi_)
    JZ = np.concatenate([pi_.conj(), Z[:2].conj()])          # J(omega, pi) = (lambda = conj pi, mu = conj omega)
    ok_J &= np.allclose(JZ, Ltz(x.conj().T, pi_.conj()))
    sig = zw(Z, JZ)
    imh = (x - x.conj().T)/(2j)
    ok_J &= abs(sig.imag) < 1e-10 and abs(sig - (-2*pi_.conj() @ imh @ pi_)) < 1e-9 and abs(sig - Z.conj() @ Hs @ Z) < 1e-9
    xh = (x + x.conj().T)/2
    ok_J &= abs(zw(Lz(xh, pi_), np.concatenate([pi_.conj(), (1j*xh @ pi_).conj()]))) < 1e-10
    A_ = crand(2, 2)
    xp = xh + 1j*(A_ @ A_.conj().T + 0.1*np.eye(2))
    sp_ = zw(Lz(xp, pi_), np.concatenate([pi_.conj(), (1j*xp @ pi_).conj()])).real
    sm_ = zw(Lz(xp.conj().T, pi_), np.concatenate([pi_.conj(), (1j*xp.conj().T @ pi_).conj()])).real
    ok_J &= sp_ < 0 < sm_
check("E: Die Twistor-Realstruktur J: PT -> PT*, (omega, pi) -> (conj pi, conj omega), induziert auf den Geraden "
      "x -> x^+, also genau die Tomita-Konjugation der Standardform (z62): J(L_x) = L~_(x^+) (200 Stichproben); "
      "Sigma(Z) = Z.JZ ist hermitesch mit Signatur (2,2) und auf L_x gleich -2 pi^+ (Im x) pi: x reell -> L_x in PN, "
      "Im x > 0 -> L_x in PT-, Im x < 0 -> PT+; J vertauscht die beiden Hardy-Hälften", ok_J)

ok_R = True
for _ in range(200):
    xh = crand(2, 2)
    xh = (xh + xh.conj().T)/2
    pi_ = crand(2)
    Z = Lz(xh, pi_)
    JZ = np.concatenate([pi_.conj(), Z[:2].conj()])
    ok_R &= np.allclose(JZ, Ltz(xh, pi_.conj())) and abs(zw(Z, JZ)) < 1e-10
    v, _ = kern(inz_matrix(pi_, pi_.conj()))
    w = np.array([[0, 1], [-1, 0]]) @ pi_.conj()            # w = eps conj(pi)
    rho = np.outer(w, w.conj())
    vm = v.reshape(2, 2)
    c_ = np.vdot(rho.ravel(), vm.ravel())/np.vdot(rho.ravel(), rho.ravel())
    ok_R &= np.allclose(vm, c_*rho) and np.allclose(rho, rho.conj().T) and np.linalg.eigvalsh(rho).min() > -1e-12
    ok_R &= abs(np.linalg.det(rho)) < 1e-12 and np.trace(rho).real > 0
check("E: Realteil: (Z, JZ) liegt in A genau für Z in PN; für reelles x ist {(Z, JZ): Z in L_x} die J-Fixmenge von Q_x, eine "
      "Sphäre S^2. Der zugehörige Lichtstrahl hat die Richtung rho = w w^+ mit w = eps conj(pi): ein reiner Zustand des "
      "Dubletts (200 Stichproben). Die Naht als Himmelssphäre (z62) ist also der J-reelle Teil von Q_x", ok_R)
befund("C (unter T und V' = conj V): Die Tomita-Paarung liefert die reelle Form der Ambitwistorquadrik: Z.JZ = 0 ist PN, A ist "
       "ihre Komplexifizierung (zweites Blatt holomorph unabhängig statt W = JZ), und die Naht ist die J-Fixmenge von Q_x. "
       "Die Quadrik selbst ist die kanonische Dualitätspaarung PT x PT* -> C; sie braucht keine neue Wahl, wohl aber, dass "
       "das zweite Blatt als holomorpher Dualraum (nicht nur als J-Bild) zählt",
       all(CHECKS[k] for k in list(CHECKS)[-5:]))

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: J und Helizität -- Realitätspartner, CPT-Bild und das Hindernis der Blätter")


def kj(e):
    return sp.expand(sp.conjugate(e))


CJ = ([kj(w) for w in A_dt[1]], [kj(w) for w in A_dt[0]])  # J(A) als Twistor: omega^A = conj mu^A', pi_A' = conj lambda_A
DJ = ([kj(w) for w in B_dt[1]], [kj(w) for w in B_dt[0]])
cJ, dJ = auf_Lt(CJ), auf_Lt(DJ)
QJ = dot2(cJ, lo(dJ))
otJ_up = up(DJ[1])
otJ_lo = lo(otJ_up)
rJ = sp.cancel(QJ/(X - Y1.H).det())
KJ = Keim(QJ)
RmJ = riemann_lin("-", KJ, otJ_lo, otJ_up)
Psi_mJ, Psit_mJ = weyl_aus_riemann(RmJ, scale=P(sp.Rational(1, 4)))


def hpunkt():
    t, s, u, v = (rq() for _ in range(4))
    return (t + s, u - I*v, u + I*v, t - s)


ok_j = all(sp.expand(otJ_up[k] - kj(o_up[k])) == 0 for k in range(2)) and rJ.free_symbols == set() and rJ != 0
ok_j &= all(v.is_zero for v in Psit_mJ.values())
ok_real, nz_p, nz_m = True, False, False
for _ in range(3):
    pt = hpunkt()
    q1, qj = ev(P(Q1), pt), ev(P(QJ), pt)
    for a in range(4):
        for b in range(4):
            ok_j &= ev(h_num("-", KJ, otJ_lo, otJ_up, a, b), pt)/qj**3 == cconj(ev(hp0[(SIG[a], SIG[b])], pt)/q1**3)
    for F in itertools.product(range(2), repeat=4):
        psJ, pst = ev(Psi_mJ[F], pt)/qj**5, ev(Psit_p[F], pt)/q1**5
        ok_real &= psJ == cconj(pst)
        nz_p |= bool(pst)
        nz_m |= bool(psJ)
check("E: J (reine Konjugation = Realstruktur) bildet die Daten von Blatt 1 auf Daten von Blatt 2 ab (C = conj A, "
      f"D = conj B, o~ = conj o): Jh+ ist ein h- mit Zentrum Y1^+ (c.d = ({rJ}) det(x - Y1^+), anderer Tubus: "
      "entgegengesetzte Frequenz, der Realitätspartner von h+); an drei reellen Punkten gilt (Jh+)_ab = "
      "conj h+_(sigma a)(sigma b) exakt, Psi~[Jh+] = 0", ok_j)

# CPT-Bild (Bisognano-Wichmann-Form von J: Konjugation mal Spiegelung x -> -x): auf Twistoren (omega, pi) -> (-omega, pi)
CP_ = ([-kj(w) for w in A_dt[1]], [kj(w) for w in A_dt[0]])
DP_ = ([-kj(w) for w in B_dt[1]], [kj(w) for w in B_dt[0]])
cP, dP = auf_Lt(CP_), auf_Lt(DP_)
QP = dot2(cP, lo(dP))
rP = sp.cancel(QP/(X + Y1.H).det())
KP = Keim(QP)
otP_up = up(DP_[1])
otP_lo = lo(otP_up)
RmP = riemann_lin("-", KP, otP_lo, otP_up)
Psi_mP, Psit_mP = weyl_aus_riemann(RmP, scale=P(sp.Rational(1, 4)))
imY1, imYP = (Y1 - Y1.H)/(2*I), ((-Y1.H) - (-Y1.H).H)/(2*I)
ok_cpt = rP.free_symbols == set() and rP != 0 and all(v.is_zero for v in Psit_mP.values())
ok_cpt &= sp.simplify(imY1 - sp.eye(2)) == sp.zeros(2, 2) and sp.simplify(imYP - imY1) == sp.zeros(2, 2)
nz_cpt = False
for _ in range(3):
    pt = hpunkt()
    mpt = tuple(-w for w in pt)
    q1m, qp = ev(P(Q1), mpt), ev(P(QP), pt)
    for a in range(4):
        for b in range(4):
            ok_cpt &= ev(h_num("-", KP, otP_lo, otP_up, a, b), pt)/qp**3 == cconj(ev(hp0[(SIG[a], SIG[b])], mpt)/q1m**3)
    for F in itertools.product(range(2), repeat=4):
        vP = ev(Psi_mP[F], pt)/qp**5
        ok_cpt &= vP == cconj(ev(Psit_p[F], mpt)/q1m**5)
        nz_cpt |= bool(vP)
check("E: CPT-Bild (Konjugation mal Spiegelung x -> -x, auf Twistoren (omega, pi) -> (-omega, pi), dann J): "
      f"(Theta J h+)(x) = conj h+(-x) ist ein h- mit Zentrum -Y1^+ (c.d = ({rP}) det(x + Y1^+)) im selben Tubus wie Y1 "
      "(Im = 1, gleiche Frequenz) und reinem Psi = conj Psi~[h+](-x) != 0 (drei reelle Punkte, exakt): CPT bildet die "
      "Helizität +2 von Blatt 1 auf eine Helizität -2 von Blatt 2 ab", ok_cpt and nz_cpt)
pt_r = hpunkt()
hpR, hmR = punktwerte("+", K1, o_lo, o_up, pt_r), punktwerte("-", KJ, otJ_lo, otJ_up, pt_r)
KR = kruemmung(hpR[0] + hmR[0], hpR[1] + hmR[1], hpR[2] + hmR[2], 2)
R2R = KR["Ric"].t.get(2)
reell = R2R is not None and nzarr(R2R) and all(cconj(R2R[b, d]) == R2R[SIG[b], SIG[d]] for b in range(4) for d in range(4))
check("E: J-reelle Störung h_R = h+ + Jh+: an reellen Punkten Psi[h_R] = conj Psi~[h_R] != 0 (alle 16 Komponenten, exakt): in "
      "Lorentzsignatur gibt es kein reelles chirales Graviton, J-Realität erzwingt beide Hälften; der Kreuzterm Ric^(2) ist "
      "!= 0 und reell (conj Ric^(2)_ab = Ric^(2)_(sigma a)(sigma b))",
      ok_real and nz_p and nz_m and KR["Ric"].t.get(1) is None and reell)

# Hindernis erster Ordnung: R(m, m) für alpha-Bivektoren m = eps^AB pi^A' pi^B' (Blatt 1) und beta-Bivektoren
# m~ = lambda^A lambda^B eps^A'B' (Blatt 2); allgemein R(m_pi, m_pi) = 4 Psi~(pi^4), R(m~, m~) = 4 Psi(lambda^4).
q1x, q2x = ev(P(Q1), x0), ev(P(Q2), x0)
Rpv = {k: ev(v, x0)/q1x**5 for k, v in Rp.items()}
Rmv = {k: ev(v, x0)/q2x**5 for k, v in Rm.items()}
Pstp = {F: ev(v, x0)/q1x**5 for F, v in Psit_p.items()}
Psm = {F: ev(v, x0)/q2x**5 for F, v in Psi_m.items()}


def rmm(Rv, m):
    tot = ZERO
    for a, b, c, d in itertools.product(range(4), repeat=4):
        if m[a][b] and m[c][d] and Rv[(a, b, c, d)]:
            tot += Rv[(a, b, c, d)]*m[a][b]*m[c][d]
    return tot


def viert(T, s):
    tot = ZERO
    for F in itertools.product(range(2), repeat=4):
        tot += T[F]*s[F[0]]*s[F[1]]*s[F[2]]*s[F[3]]
    return tot


ok_h = True
VIER = qi(4)
for _ in range(5):
    pu = [qi(cq()), qi(cq())]
    lu = [qi(cq()), qi(cq())]
    m_a = [[qi(EPS[IDX[a][0]][IDX[b][0]])*pu[IDX[a][1]]*pu[IDX[b][1]] for b in range(4)] for a in range(4)]
    m_b = [[lu[IDX[a][0]]*lu[IDX[b][0]]*qi(EPS[IDX[a][1]][IDX[b][1]]) for b in range(4)] for a in range(4)]
    ap, am = rmm(Rpv, m_a), rmm(Rmv, m_a)
    bp, bm = rmm(Rpv, m_b), rmm(Rmv, m_b)
    ok_h &= ap == VIER*viert(Pstp, pu) and bool(ap) and not am
    ok_h &= bm == VIER*viert(Psm, lu) and bool(bm) and not bp
check("E: Hindernis erster Ordnung (5 zufällige pi, lambda, exakt): h+ (Grad -6 von Blatt 1) gibt R(m_pi, m_pi) = "
      "4 Psi~(pi^4) != 0 für alpha-Ebenen (Blatt 1 gehindert) und R(m~, m~) = 0 für beta-Ebenen (Blatt 2 bleibt); h- "
      "umgekehrt", ok_h)
googly = ok_h and ks_ok
befund("E: googly-Dualität. Das Grad-(-6)-Feld eines Blatts zerstört dieses Blatt und deformiert das andere: h+ (Psi~, das "
       "Feld aus z61) hat keine alpha-Flächen, aber beta-Flächen, auch nichtlinear (Kerr-Schild: Psi = 0 in allen "
       "Ordnungen); es ist ein nichtlineares Graviton von PT* (Blatt 2), nicht von PT. Penroses Deformation von PT trägt die "
       "Gegenhelizität. Je eine Helizität lebt also nichtlinear auf je einem Blatt, über Kreuz zur linearen Zuordnung",
       googly)
befund("E: Für beide Helizitäten zugleich überlebt kein Blatt: Kreuzterm != 0 (Teil 3), und eine J-reelle nichtflache "
       "Störung hat Psi und Psi~ != 0, also weder alpha- noch beta-Flächen. Übrig bleibt nur A (Nullgeodäten, Realteil PN); "
       "einzig der kollineare pp-Sektor überlagert sich exakt",
       CHECKS[[k for k in CHECKS if "Kreuzterm, exakt" in k][0]] and CHECKS[[k for k in CHECKS if "J-reelle" in k][0]]
       and CHECKS[[k for k in CHECKS if "kollinear (pp-Welle" in k][0]])

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: Bilanz -- was bekannt ist, was TFPT bräuchte")
bekannt = [
    ("selbstdual, nichtlinear [E, Literatur]",
     "Penrose 1976 (nichtlineares Graviton): Deformationen von PT mit Faserung über CP1 und O(2)-wertiger symplektischer Form "
     "<-> Vakuummetriken einer Chiralität (hier: Psi~ = 0, Psi != 0); Ward 1977 ebenso für selbstduale Yang-Mills-Felder"),
    ("googly [O, Literatur]",
     "Grad-(-6)-Felder auf PT (z61) haben keine Deformationsbeschreibung auf PT; ein gekrümmter Twistorraum für beide "
     "Helizitäten ist nicht bekannt (Penrose 2015, palatial twistor theory, nicht abgeschlossen)"),
    ("volles Yang-Mills [E, Literatur]",
     "Witten 1978, Isenberg-Yasskin-Green 1978: jeder Zusammenhang <-> Bündel auf A, trivial auf Q_x; Yang-Mills <=> "
     "Fortsetzung in die formale Umgebung dritter Ordnung von A in PT x PT*; N = 3 Super-Yang-Mills ohne Fortsetzung "
     "(Integrabilität auf Super-Lichtstrahlen)"),
    ("konforme Struktur [E, Literatur]",
     "LeBrun 1983: A (Raum der komplexen Nullgeodäten mit Kontaktstruktur) kodiert jede komplexe konforme Struktur, ohne "
     "Feldgleichungen"),
    ("Einstein, linearisiert [E, Literatur]",
     "Baston-Mason 1987 (laut Abstract): Deformationen von PN setzen sich stets bis zur vierten formalen Umgebung fort, bis "
     "zur fünften genau "
     "bei Bach = 0, bis zur sechsten bei konform-Einstein (algebraisch allgemeine Weylkrümmung)"),
    ("Einstein, voll gekrümmt [E, Literatur]",
     "LeBrun 1991 (Thickenings and conformal gravity): Twistorkorrespondenz für Bach- und Eastwood-Dighton-flache "
     "Raumzeiten (algebraisch allgemein: konform Einstein); Lambda = 0 und die Skala brauchen Zusatzdaten"),
    ("störungstheoretisch [C, Literatur]",
     "Mason 2005, Mason-Skinner 2008, Adamo-Mason 2013: Twistorwirkungen; die zweite Helizität ist ein Feld vom Gewicht -6 "
     "auf demselben PT, das Nicht-Selbstduale eine Störung um den selbstdualen Sektor (keine direkte klassische "
     "Äquivalenz bewiesen). Ambitwistor-Strings: volle Baumamplituden (Mason-Skinner 2014); auf gekrümmtem Hintergrund "
     "ist Anomaliefreiheit äquivalent zu den Supergravitationsgleichungen (Adamo-Casali-Skinner 2015)"),
]
for k_, v_ in bekannt:
    say(f"   {k_}: {v_}")
braucht = [
    ("zweites Blatt = PT* [C]", "unter T und V' = conj V ist das J-Bild von Blatt 1 der duale Twistorraum (J(L_x) = L~_(x^+))"),
    ("Quadrik [E/C]", "Z.JZ = 0 ist PN, die Naht ist die J-Fixmenge von Q_x; A ist die Komplexifizierung. Gebraucht wird, "
                      "dass Blatt 2 holomorph unabhängige Koordinaten trägt (Lesart: unabhängige Anregungen der Kommutante im "
                      "TFD) [L, nicht hergeleitet]"),
    ("lineare Gravitation [E]", "beide Helizitäten aus je einem Blatt, h = h+ + h- mit R^(1) = 0; linear genügt aber schon "
                                "ein Blatt (Grad -6 und +2)"),
    ("nichtlinear, eine Helizität [E/C]", "Deformation des jeweils anderen Blatts (googly-Dualität); Kerr-Schild-"
                                          "Elementarzustände exakt, allgemein Penrose"),
    ("nichtlinear, beide Helizitäten [O]",
     "weder PT noch PT* überlebt; nötig ist eine Deformation von A mit Daten in der formalen Umgebung (Ordnung 3 für "
     "Yang-Mills; 5 bzw. 6 für Bach bzw. Einstein linearisiert; voll gekrümmt LeBrun) und eine Skala (Lambda = 0, G). "
     "TFPT hat dafür keine Struktur; ein Kandidat wären Korrelationen der Blätter jenseits der J-Diagonale (Z.W != 0), "
     "nicht gezeigt"),
]
for k_, v_ in braucht:
    say(f"   TFPT braucht -- {k_}: {v_}")
ZAHLEN["bekannt"] = dict(bekannt)
ZAHLEN["tfpt_braucht"] = dict(braucht)
ZAHLEN["status"] = {
    "E": ["Grad -6 auf PT -> Psi~, auf PT* -> Psi, Grad +2 auf PT -> Psi (Bianchi exakt)",
          "h = h+ + h- mit R^(1) = 0, Weylteile = Penrose-Transformierte beider Blätter; allgemein über Symbole",
          "Kerr-Schild-Elementarzustände exakt nichtlinear je einer Chiralität; Kreuzterm != 0; pp-Sektor exakt",
          "Ambitwistorquadrik, Nullabstand, Klein-Korrespondenz, J = Realstruktur, Naht = J-Fixmenge von Q_x",
          "J-Realität erzwingt beide Hälften; googly-Dualität linear und nichtlinear (Kerr-Schild)"],
    "C": ["Blatt 2 = PT* und Tomita-Paarung = reelle Form der Quadrik, unter T und V' = conj V"],
    "L": ["TFPTs J als Bisognano-Wichmann-Konjugation (CPT mal Drehung) lesen: dann bildet es Helizität +2 von Blatt 1 "
          "auf Helizität -2 von Blatt 2 gleicher Frequenz ab (die Mathematik dazu ist E, Teil 5); unabhängige "
          "Anregungen beider Blätter als TFD-Kommutante"],
    "O": ["nichtlineare nicht-selbstduale Gravitation: Deformation von A mit formaler Umgebung und Skala, in TFPT nicht "
          "vorhanden; das googly-Problem ist hier nicht gelöst"],
}
befund("O: Das nichtlineare nicht-selbstduale Problem ist hier nicht gelöst. Die zwei Blätter tragen die zwei Helizitäten "
       "linear vollständig und nichtlinear je einzeln (über Kreuz); für beide zugleich bleibt nur A mit Fortsetzungsdaten, "
       "die TFPT nicht hat", googly and CHECKS[[k for k in CHECKS if "Kreuzterm, exakt" in k][0]])

urteil = ("Die zwei Blätter der Naht können die zwei Helizitäten tragen, aber nur linear und je einzeln: Grad -6 auf PT (Blatt "
          "1) gibt Psi~, Grad -6 auf PT* (Blatt 2) gibt Psi, und h = h+ + h- aus je einem Elementarzustand löst R^(1) = 0 "
          "exakt mit beiden Weylhälften. Linear genügt allerdings schon ein Blatt (Grad +2 gibt Psi), und nichtlinear gilt "
          "eine googly-Dualität: das Grad-(-6)-Feld eines Blatts (z61s Psi~) ist ein nichtlineares Graviton des anderen "
          "Blatts (Kerr-Schild exakt), Penroses Deformation von PT trägt die Gegenhelizität. Beide zugleich zerstören beide "
          "Blätter (Kreuzterm != 0; J-reelle Störungen tragen stets beide Hälften); es bleibt die Ambitwistorquadrik A. "
          "TFPTs Tomita-Konjugation ist genau die Realstruktur PT -> PT* (J(L_x) = L~_(x^+)), Z.JZ = 0 ist PN und die Naht "
          "die J-Fixmenge von Q_x; J allein liefert den Realitätspartner, J mit Spiegelung (CPT) eine Helizität -2 auf "
          "Blatt 2 gleicher Frequenz; die volle Gravitation verlangt aber Daten in der formalen Umgebung von A (Baston-Mason, "
          "LeBrun) und eine Skala, die TFPT nicht hat. Das nichtlineare nicht-selbstduale Problem ist nicht gelöst.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z65_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
