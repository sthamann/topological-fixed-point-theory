"""Z73: Ist die 4D-Struktur in der Teilchen-Loch-Struktur des Originals nativ? (Prüfung einer externen Analyse, 24.11.)

Behauptung der Analyse: Die m Teilchen-Loch-Paare der Strommode m tragen einen Spin j = (m-1)/2; die Tomita-Verdopplung gibt
(j) x (j-quer) = K-Typ (j, j); die Summe ist der Einteilchencharakter des masselosen konform gekoppelten 4D-Skalars, also der
SU(2,2)-Leiterdarstellung. Original: quasifreie Quarterkovarianz c(x) = 1/(1 + e^{4 pi x}), x in Z + alpha, alpha = 1/4, 3/4,
KMS bei beta = 4 pi, Strom J_m = sum_x psi*_x psi_{x+m} (z58). Fockvakuum Omega = Fermisee (Hardy-Schnitt r < 0).

Teil 1  Bilokale Auslese D_m(gamma) = sum_r (c_r - c_{r+m}) e^{i gamma (2r+m+1)}: gamma ist (bis auf die gemeinsame Phase e^{i gamma})
        der halbe Spreizwinkel des punktgespaltenen Bilinears psi*(Theta + gamma) psi(Theta - gamma); D_m ist dessen Kubo-Antwort
        auf J_m. Faktorisierung D_m = D_1 sin(m gamma)/sin(gamma) (Teleskopsumme, für jedes Besetzungsprofil), die Identität
        d^2_gamma W(tau, 0) = -(1/3)(d^2_tau - 1) W(tau, 0) für W = sum_m (D_m/D_1) g_m(tau), und ihre exakte 2D-Form: sin(gamma) W
        ist harmonisch in (tau, gamma) mit Dirichlet-Rand, d_gamma(sin(gamma) W) ist die native Nahtwinkelantwort W2 aus z71.
Teil 2  Natives su(2) auf dem Paarraum H_m = span{psi*_{r+m} psi_r Omega} (exakt im Fockraum, Jordan-Wigner, Fenster von 22
        Moden): No-go für Bilineare (jede Kompression ist diagonal), native Cartan-Größe, Stromhüpfen (Hardy-Hälften von
        J_{-1}, J_1 und P J_{-1} J_1 P), Möbiushüpfen (P L_{-1} L_1 P im Twistsektor), Obstruktion für jede feste
        Polynomgewichtung, Dicke-su(2) der m disjunkten Paare und die Schwinger-Gewichte sqrt(n_1 + 1) sqrt(n_2), die su(2)
        erst herstellen.
Teil 3  Tomita-Verdopplung: sum m^2 q^m = q(1+q)/(1-q)^3 = Einteilchencharakter des konformen Skalars auf R x S^3 (K-Typen (j, j)
        bei Energie 2j + 1, voller Toruscharakter), native verdoppelte Auslese, Standardform von End(C^m), Kompression des
        KMS-Zustands auf H_m (exakt der Spurzustand) und die Kubo-Gewichte der Paare (exakte Verhältnisse in q = e^{-pi}).
Teil 4  SU(2,2) auf V + V-quer* (Hermitesche Form der Signatur (2,2)), Oszillatorquantisierung Z = (omega, pi),
        P_{AA'} = pi-quer_A pi_{A'}: [P, P] = 0, P^0 >= 0, P^2 = 0, E = P^0 + K^0, K-Typen im Helizität-0-Sektor; Vergleich mit
        den nativen Hebern (Toeplitz statt CCR).
Teil 5  Bilanz (E/C/L/O) und Urteil.
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.sparse as sps
import sympy as sp

HERE = Path(__file__).parent
mp.mp.dps = 50
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(73)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


B4 = 4*mp.pi
ALPHAS = (mp.mpf(1)/4, mp.mpf(3)/4)


def occ(x, beta=B4):
    return 1/(1 + mp.exp(beta*x))


def D_all(g, alpha, beta, M):
    """native bilokale Auslese D_1..D_M bei Winkel g (direkt aus c, ohne Faktorisierung)"""
    R = int(mp.ceil(150/beta)) + 2
    cv = {k: occ(k + alpha, beta) for k in range(-R - M - 1, R + M + 2)}
    ks = range(-R - M, R + 1)
    ph = {k: mp.expj(g*(2*(k + alpha) + 1)) for k in ks}
    return [mp.expj(g*m)*mp.fsum((cv[k] - cv[k + m])*ph[k] for k in ks) for m in range(1, M + 1)]


def chi(m, g):
    return mp.sin(m*g)/mp.sin(g) if g != 0 else mp.mpf(m)


# ================================================================== Teil 1
say("Teil 1: Bilokale Auslese, Faktorisierung und die Identität für W")
xs = sp.symbols("x")
KS = 4
csym = {k: sp.Symbol(f"c{k}") for k in range(-KS, KS)}


def cval(k, prof):
    if prof == "stufe":
        return sp.Integer(1 if k < 0 else 0)
    if k < -KS:
        return sp.Integer(1)
    if k >= KS:
        return sp.Integer(0)
    return csym[k]


def Dsym(m, prof):                                    # Faktor x^(2 alpha) abgespalten, x = e^{i gamma}
    return sp.expand(sum((cval(k, prof) - cval(k + m, prof))*xs**(2*k + m + 1) for k in range(-KS - m - 2, KS + 2)))


ok_sym = True
for m in range(1, 7):
    ok_sym &= sp.expand(Dsym(m, "frei")*(xs - 1/xs) - Dsym(1, "frei")*(xs**m - xs**(-m))) == 0
D1_stufe = Dsym(1, "stufe")
ok_stufe = D1_stufe == 1 and all(sp.expand(Dsym(m, "stufe") - sum(xs**(m - 1 - 2*k) for k in range(m))) == 0
                                 for m in range(1, 7))
check("[E] D_m(gamma) sin(gamma) = D_1(gamma) sin(m gamma) gilt exakt für ein BELIEBIGES Besetzungsprofil (16 freie Symbole "
      "c_k, m = 1..6, sympy): die Faktorisierung ist eine Teleskopsumme, keine Eigenschaft des Zustands; bei T = 0 (Stufe) ist "
      "D_1 = e^{2 i alpha gamma} und die m Paare tragen je Gewicht 1 mit w - alpha = -j..j", ok_sym and ok_stufe)

res = mp.mpf(0)
for beta in (mp.mpf(1), 2*mp.pi, B4, 8*mp.pi):
    for al in ALPHAS:
        for g in (mp.mpf("0.3"), mp.mpf("1.1"), mp.mpf("2.0"), mp.mpc("0.7", "0.2")):
            Ds = D_all(g, al, beta, 10)
            for m in range(1, 11):
                res = max(res, abs(Ds[m - 1]*mp.sin(g) - Ds[0]*mp.sin(m*g)))
ZAHLEN["faktorisierung_max_residuum"] = mp.nstr(res, 5)
check(f"[E] numerisch für das Original: beide Sektoren, beta in {{1, 2 pi, 4 pi, 8 pi}}, m = 1..10, vier Winkel (einer komplex), "
      f"50 Stellen: max |D_m sin(gamma) - D_1 sin(m gamma)| = {mp.nstr(res, 3)}", res < mp.mpf(10)**-40)


def gm(m, tau):
    return (mp.exp(-m*tau) + mp.exp(-m*(B4 - tau)))/(1 - mp.exp(-B4*m))


def W(tau, g, al=ALPHAS[0], norm=True):
    M = int(110/min(tau, B4 - tau)) + 5
    Ds = D_all(g, al, B4, M)
    return mp.fsum((Ds[m - 1]/Ds[0] if norm else Ds[m - 1])*gm(m, tau) for m in range(1, M + 1))


def Gcur(tau):                                        # z58: thermische Stromantwort
    M = int(110/min(tau, B4 - tau)) + 5
    return mp.fsum(m/(1 - mp.exp(-B4*m))*mp.exp(-m*tau) + m/(mp.exp(B4*m) - 1)*mp.exp(m*tau) for m in range(1, M + 1))


dev_n, dev_u, dev_line = mp.mpf(0), mp.mpf(1)*10**6, mp.mpf(0)
for al in ALPHAS:
    for tau in (mp.mpf("1.3"), mp.mpf("2.9"), mp.mpf("6.1")):
        lhs = mp.diff(lambda g: W(tau, g, al), 0, 2)
        rhs = -(mp.diff(lambda t: W(t, 0, al), tau, 2) - W(tau, 0, al))/3
        dev_n = max(dev_n, abs(lhs - rhs)/abs(rhs))
        dev_line = max(dev_line, abs(W(tau, 0, al) - Gcur(tau))/Gcur(tau))
        lhs2 = mp.diff(lambda g: W(tau, g, al, False), 0, 2)
        rhs2 = -(mp.diff(lambda t: W(t, 0, al, False), tau, 2) - W(tau, 0, al, False))/3
        dev_u = min(dev_u, abs(lhs2 - rhs2)/abs(rhs2))
ZAHLEN["W_identitaet"] = {"normiert_max_rel": mp.nstr(dev_n, 5), "linie_gleich_stromantwort_max_rel": mp.nstr(dev_line, 5),
                          "unnormiert_min_rel_abweichung": mp.nstr(dev_u, 5)}
check("[E] W(tau, gamma) = sum_m (D_m/D_1)(gamma) g_m(tau) mit den nativen D_m: W(tau, 0) ist die Stromantwort aus z58, und "
      f"d^2_gamma W(tau, 0) = -(1/3)(d^2_tau - 1) W(tau, 0) gilt (beide Sektoren, drei tau, rel. {mp.nstr(dev_n, 3)}); ohne "
      f"Normierung durch D_1 (Phase e^(2 i alpha gamma), thermische Verschmierung) verfehlt sie die Identität (rel. >= "
      f"{mp.nstr(dev_u, 3)})", dev_n < mp.mpf(10)**-35 and dev_line < mp.mpf(10)**-35 and dev_u > mp.mpf("0.1"))


def W2(tau, th):                                      # z71: native Antwort längs des Nahtwinkels
    M = int(110/min(tau, B4 - tau)) + 5
    return mp.fsum(m*mp.cos(m*th)*gm(m, tau) for m in range(1, M + 1))


al = ALPHAS[0]
d2u, ds3, dW2, d2W = mp.mpf(0), mp.mpf(0), mp.mpf(0), mp.mpf(10)**6
for tau, g in ((mp.mpf("1.3"), mp.mpf("0.4")), (mp.mpf("2.9"), mp.mpf("1.1")), (mp.mpf("6.1"), mp.mpf("2.0"))):
    def u(t, gg):
        return mp.sin(gg)*W(t, gg, al)
    lap_u = mp.diff(u, (tau, g), (2, 0)) + mp.diff(u, (tau, g), (0, 2))
    Wv = W(tau, g, al)
    s3 = (mp.diff(lambda t: W(t, g, al), tau, 2) + mp.diff(lambda gg: W(tau, gg, al), g, 2)
          + 2*mp.cot(g)*mp.diff(lambda gg: W(tau, gg, al), g, 1) - Wv)
    d_u = mp.diff(lambda gg: u(tau, gg), g, 1)
    lap_W = mp.diff(lambda t: W(t, g, al), tau, 2) + mp.diff(lambda gg: W(tau, gg, al), g, 2)
    d2u = max(d2u, abs(lap_u)/abs(Wv))
    ds3 = max(ds3, abs(s3)/abs(Wv))
    dW2 = max(dW2, abs(d_u - W2(tau, g))/abs(W2(tau, g)))
    d2W = min(d2W, abs(lap_W)/abs(Wv))
ZAHLEN["zweiD_aequivalenz"] = {"laplace_sin_W": mp.nstr(d2u, 5), "S3_zonal": mp.nstr(ds3, 5),
                               "d_gamma(sin W) - W2": mp.nstr(dW2, 5), "laplace_W_selbst_min": mp.nstr(d2W, 5)}
check("[E] dieselben Daten zweidimensional: u = sin(gamma) W ist harmonisch in (tau, gamma) (Dirichlet bei gamma = 0, pi), "
      "W erfüllt die S^3-Zonalgleichung (d^2_tau + d^2_gamma + 2 cot(gamma) d_gamma - 1) W = 0, und d_gamma(sin(gamma) W) ist "
      "exakt die native Nahtwinkelantwort W2(tau, gamma) aus z71; W selbst ist nicht 2D-harmonisch (drei Punkte, 50 Stellen)",
      d2u < mp.mpf(10)**-30 and ds3 < mp.mpf(10)**-30 and dW2 < mp.mpf(10)**-30 and d2W > mp.mpf("0.01"))
befund("Teil 1 (E): Die Faktorisierung ist kinematisch (für jedes Profil). Die Identität gilt nur für das durch D_1 normierte W. "
       "Die S^3-Zonalform ist dieselbe Information wie eine 2D-Dirichlet-Antwort im Spreizwinkel, W = (1/sin(gamma)) "
       "int_0^gamma W2 dtheta. Die Linien- und Bilokaldaten entscheiden also nicht zwischen 2D und 4D: Nur Klassenfunktionen "
       "(ein Winkel von S^3) sind nativ", True)

# ================================================================== Teil 2
say("\nTeil 2: Natives su(2) auf dem Paarraum H_m (exakt im Fockraum)")
KMIN, KMAX = -11, 10
MODES = range(KMIN, KMAX + 1)


def bit(k):
    return 1 << (k - KMIN)


VAC = sum(bit(k) for k in range(KMIN, 0))


def nbelow(s, k):
    return bin(s & ((1 << (k - KMIN)) - 1)).count("1")


def bil(a, b, vec, coef=1):
    """coef psi*_a psi_b auf einem Fockvektor {Bitmaske: Koeffizient} (Jordan-Wigner, aufsteigende Ordnung)"""
    out = {}
    for s, cc in vec.items():
        if not (s >> (b - KMIN)) & 1:
            continue
        sg = -1 if nbelow(s, b) % 2 else 1
        s1 = s ^ bit(b)
        if (s1 >> (a - KMIN)) & 1:
            continue
        if nbelow(s1, a) % 2:
            sg = -sg
        s2 = s1 | bit(a)
        out[s2] = out.get(s2, 0) + coef*sg*cc
    return out


def vadd(u, v, f=1):
    w = dict(u)
    for k, cc in v.items():
        w[k] = w.get(k, 0) + f*cc
    return {k: cc for k, cc in w.items() if cc != 0}


def opsum(terms):
    def f(vec):
        out = {}
        for coef, a, b in terms:
            if KMIN <= a <= KMAX and KMIN <= b <= KMAX:
                out = vadd(out, bil(a, b, vec, coef))
        return out
    return f


def comp(*fs):
    def g(vec):
        for f in reversed(fs):
            vec = f(vec)
        return vec
    return g


def inner(u, v):
    return sum(cc*v.get(s, 0) for s, cc in u.items())


def Jop(n):                                           # J_n = sum psi*_x psi_{x+n}
    return opsum([(1, k, k + n) for k in MODES])


def Lop(n, alpha):                                    # Möbius/Virasoro, Gewicht 1/2, Twistsektor: s = alpha + k
    return opsum([(alpha + k - Fr(n, 2), k - n, k) for k in MODES])


Xp = opsum([(1, k + 1, k) for k in MODES if k >= 0])            # Hardy-Hälfte von J_{-1}: Teilchen hoch
Xpd = opsum([(1, k, k + 1) for k in MODES if k >= 0])           # adjungiert: Teilchen runter
Xh = opsum([(1, k, k - 1) for k in MODES if k <= -1])           # Hardy-Hälfte von J_{-1}: Loch tiefer
Xhd = opsum([(1, k - 1, k) for k in MODES if k <= -1])          # adjungiert (Teil von J_1): Loch flacher


def basis(m):                                         # b_k = psi*_{k+m} psi_k Omega, k = -m..-1 (Gewicht w - alpha = k + (m+1)/2)
    return [bil(kr + m, kr, {VAC: 1}) for kr in range(-m, 0)]


def compress(f, Bs):
    return [[inner(Bs[i], f(Bs[j])) for j in range(len(Bs))] for i in range(len(Bs))]


def leak(f, Bs):
    tot = 0
    for b in Bs:
        v = f(b)
        for bb in Bs:
            v = vadd(v, bb, -inner(bb, v))
        tot += sum(abs(complex(cc)) for cc in v.values())
    return tot


def L0_of(state, alpha):
    return sum((alpha + k)*(((state >> (k - KMIN)) & 1) - (1 if k < 0 else 0)) for k in MODES)


ok_pair = True
for m in range(1, 8):
    Bs = basis(m)
    v = Jop(-m)({VAC: 1})
    co = [inner(b, v) for b in Bs]
    rest = v
    for b, cc in zip(Bs, co):
        rest = vadd(rest, b, -cc)
    gram = [[inner(Bs[i], Bs[j]) for j in range(m)] for i in range(m)]
    ok_pair &= (not rest and all(abs(cc) == 1 for cc in co) and gram == [[int(i == j) for j in range(m)] for i in range(m)]
                and all(L0_of(next(iter(b)), al_) == m for b in Bs for al_ in (Fr(1, 4), Fr(3, 4))))
check("[E] J_{-m} Omega = sum von genau m orthonormalen Ein-Paar-Zuständen psi*_{r+m} psi_r Omega (Loch r in (-m, 0)), alle mit "
      "L_0 = m (m = 1..7, beide Sektoren): H_m ist m-dimensional", ok_pair)

nogo, diag_seen = True, set()
for m in range(1, 7):
    Bs = basis(m)
    for a in MODES:
        for b in MODES:
            Mx = compress(opsum([(1, a, b)]), Bs)
            if a != b:
                nogo &= all(Mx[i][j] == 0 for i in range(m) for j in range(m))
            else:
                diag_seen |= {Mx[i][i] for i in range(m)}
check("[E] No-go für Bilineare: P_m psi*_a psi_b P_m = 0 für alle a != b im Fenster (22 Moden, m = 1..6); nur die "
      "Besetzungszahlen wirken (diagonal), Paarbildner psi psi, psi* psi* wechseln die Ladung. Also komprimiert JEDE "
      "Bilinearform (J_n, L_n, W-Moden, beliebige Gewichte) auf H_m diagonal und erzeugt nur eine abelsche Algebra -- kein "
      "Heber oder Senker", nogo and diag_seen == {0, 1})


def Hcartan(vec, alpha):
    """H = 1/2 (L_0^(Teilchen) - L_0^(Loch)) + 1/2 - alpha, diagonal (nativ mit Hardy-Schnitt)"""
    out = {}
    for s, cc in vec.items():
        lp = sum(alpha + k for k in MODES if k >= 0 and (s >> (k - KMIN)) & 1)
        lh = sum(-(alpha + k) for k in MODES if k < 0 and not (s >> (k - KMIN)) & 1)
        out[s] = cc*(Fr(1, 2)*(lp - lh) + Fr(1, 2) - alpha)
    return out


ok_cart = True
for m in range(1, 8):
    Bs = basis(m)
    for al_ in (Fr(1, 4), Fr(3, 4)):
        Hm = compress(lambda v: Hcartan(v, al_), Bs)
        wt = [Fr(kr) + Fr(m + 1, 2) for kr in range(-m, 0)]
        bil_w = [(2*(al_ + kr) + m + 1)/2 - al_ for kr in range(-m, 0)]
        ok_cart &= all(Hm[i][j] == (wt[i] if i == j else 0) for i in range(m) for j in range(m)) and wt == bil_w
check("[E] native Cartan-Größe: H = (L_0^Teilchen - L_0^Loch)/2 + 1/2 - alpha ist auf H_m diagonal mit den Gewichten -j..j, "
      "genau die bilokalen Gewichte (2r + m + 1)/2 - alpha (m = 1..7, beide Sektoren)", ok_cart)


def su2_ratios(Mx, m):
    """|e_w|^2 / ((j - w)(j + w + 1)) für den w-hebenden Teil; konstant <=> Abschluss zu su(2) bis auf Skala"""
    j = Fr(m - 1, 2)
    out = []
    for i in range(m - 1):
        w = -j + i
        e = Mx[i + 1][i]
        out.append(sp.nsimplify(sp.Rational(str(Fr(e))) if not isinstance(e, sp.Basic) else e)**2
                   / sp.Rational(str((j - w)*(j + w + 1))))
    return out


def closes(rat):
    return len(set(rat)) <= 1 and all(r != 0 for r in rat)


tab_cur, tab_moeb, ok_cur_struct, amp_ok = {}, {}, True, True
for m in range(2, 8):
    Bs = basis(m)
    Ep = comp(Xp, Xhd)
    Mx = compress(Ep, Bs)
    QJ = compress(comp(Jop(-1), Jop(1)), Bs)
    ok_cur_struct &= leak(Ep, Bs) == 0 and all(abs(Mx[i + 1][i]) == 1 for i in range(m - 1))
    ok_cur_struct &= all(QJ[i][k] == Mx[i][k] + Mx[k][i] for i in range(m) for k in range(m) if i != k)
    ok_cur_struct &= all(Mx[i][k] == 0 for i in range(m) for k in range(m) if i != k + 1)
    Mn = np.array(Mx, dtype=float)
    ok_cur_struct &= np.allclose(Mn @ Mn.T - Mn.T @ Mn, np.diag([-1.0] + [0.0]*(m - 2) + [1.0]))     # P_oben - P_unten
    rat = su2_ratios(Mx, m)
    tab_cur[m] = {"verhaeltnisse": [str(r) for r in rat], "schliesst": closes(rat)}
    for al_ in (Fr(1, 4), Fr(3, 4)):
        QL = compress(comp(Lop(-1, al_), Lop(1, al_)), Bs)
        amp_ok &= all(abs(QL[i + 1][i]) == abs((al_ - m + i + Fr(1, 2))*(al_ + i + Fr(1, 2))) for i in range(m - 1))
        ratL = su2_ratios(QL, m)
        tab_moeb[(m, str(al_))] = {"verhaeltnisse": [str(r) for r in ratL], "schliesst": closes(ratL)}
ZAHLEN["stromhuepfen_su2_test"] = tab_cur
ZAHLEN["moebiushuepfen_su2_test"] = {f"m={k[0]}, alpha={k[1]}": v for k, v in tab_moeb.items()}
say(f"   Stromhüpfen |e_w|^2/((j-w)(j+w+1)): {{m: {', '.join(str(m) + ': ' + str(v['verhaeltnisse']) for m, v in tab_cur.items())}}}")
check("[E] Stromhüpfen: E_+ = X_p X_h^+ (Produkt der Hardy-Hälften von J_{-1} und J_1) lässt H_m exakt invariant, hebt w um 1 "
      "mit |Koeffizient| = 1 und ist genau der hebende Teil von P_m J_{-1} J_1 P_m; es schließt (bis auf eine Skala) zu su(2) "
      "mit Spin (m-1)/2 genau für m = 2, 3 und für kein m = 4..7 ([E_+, E_-] = P_oben - P_unten statt 2H)",
      ok_cur_struct and [m for m, v in tab_cur.items() if v["schliesst"]] == [2, 3])
check("[E] Möbiushüpfen: P_m L_{-1} L_1 P_m (Virasoro-Gewicht 1/2 im Quartersektor) hat die Hüpfamplituden (r + 1/2)(p + 1/2) "
      "(Loch r, Teilchen p vor dem Schritt) und schließt nur für m = 2 zu su(2) (m = 2..7, beide Sektoren)",
      amp_ok and sorted({k[0] for k, v in tab_moeb.items() if v["schliesst"]}) == [2])

hh0, hh1, hh2, cc_ = sp.symbols("h0 h1 h2 c")
eqs6 = [sp.expand((hh0 + hh1*k + hh2*k**2)**2 - cc_*(k + 1)*(5 - k)) for k in range(5)]
sol6 = sp.solve(eqs6, [hh0, hh1, hh2, cc_], dict=True)
bb_, gg_ = sp.symbols("b g")
A4 = [(1 + bb_*k)*(1 + gg_*(3 - k)) for k in range(3)]
w4 = [(k + 1)*(3 - k) for k in range(3)]
sol4 = sp.solve([sp.expand(A4[k]**2*w4[0] - A4[0]**2*w4[k]) for k in (1, 2)], [bb_, gg_], dict=True)
sol4_real = [s for s in sol4 if all(sp.im(sp.N(v)) == 0 for v in s.values())]
ZAHLEN["polynomgewichte"] = {"m=6, Grad<=2: Loesungen": [str(s) for s in sol6],
                             "m=4, affin x affin: reelle Abstimmungen": [{str(k): float(sp.N(v)) for k, v in s.items()}
                                                                        for s in sol4_real]}
check("[E] allgemein: Ein Hüpfen mit Amplitude h(n_1) = Polynom vom Grad <= 2 in der Paarlage (umfasst alle Produkte affiner "
      "Teilchen- und Lochgewichte, also Strom, Möbius und Mischungen, mal jeder Funktion von L_0) schließt bei m = 6 nie zu "
      f"su(2): einzige Lösung h = 0 (sympy, {sol6}). Einzelne Stufen lassen sich abstimmen (m = 4: {len(sol4_real)} reelle, "
      "irrationale Gewichtungen), das ist stufenweise Auswahl, keine Rekonstruktion",
      sol6 == [{cc_: 0, hh0: 0, hh1: 0, hh2: 0}] and len(sol4_real) > 0
      and all(any(abs(float(sp.N(v)) - round(float(sp.N(v)))) > 1e-6 for v in s.values()) for s in sol4_real))

ok_dicke = True
dicke = {}
for m in range(2, 7):
    # Unterraum der m disjunkten Paare (k, k+m), k = -m..-1: 2^m Zustände durch Umklappen aus Omega
    states = []
    for flips in itertools.product((0, 1), repeat=m):
        s = VAC
        for f, kr in zip(flips, range(-m, 0)):
            if f:
                s = (s ^ bit(kr)) | bit(kr + m)
        states.append(s)
    Tp = opsum([(1, kr + m, kr) for kr in range(-m, 0)])
    Tm = opsum([(1, kr, kr + m) for kr in range(-m, 0)])

    def Tz(vec, m=m):
        return {s: cc*Fr(sum(((s >> (kr + m - KMIN)) & 1) - ((s >> (kr - KMIN)) & 1) for kr in range(-m, 0)), 2)
                for s, cc in vec.items()}
    Sb = [{s: 1} for s in states]
    C1 = compress(lambda v: vadd(comp(Tp, Tm)(v), comp(Tm, Tp)(v), -1), Sb)
    C2 = compress(lambda v: {k: 2*c_ for k, c_ in Tz(v).items()}, Sb)
    ok_dicke &= C1 == C2
    Bs = basis(m)
    cas = compress(lambda v: vadd(vadd(comp(Tz, Tz)(v), {k: c_*Fr(1, 2) for k, c_ in comp(Tp, Tm)(v).items()}),
                                  {k: c_*Fr(1, 2) for k, c_ in comp(Tm, Tp)(v).items()}), Bs)
    ev = sorted(np.round(np.linalg.eigvalsh(np.array(cas, dtype=float)), 9))
    want = sorted([(m/2 - 1)*(m/2)]*(m - 1) + [(m/2)*(m/2 + 1)])
    dicke[m] = ev
    ok_dicke &= np.allclose(ev, want) and leak(Tp, Bs) != 0
ZAHLEN["dicke_casimir_auf_H_m"] = {m: [float(x) for x in v] for m, v in dicke.items()}
check("[E] Dicke-su(2) der m disjunkten Paare (Fensterteil von J_{-m}, J_m, nativ mit Hardy-Schnitt): schließt exakt, "
      "[T_+, T_-] = 2 T_z auf allen 2^m Zuständen, wechselt aber die Paarzahl; H_m ist ihr Ein-Flip-Gewichtsraum und zerfällt in "
      "Spin m/2 (einmal) und Spin m/2 - 1 ((m-1)-mal), nicht in Spin (m-1)/2 (m = 2..6)", ok_dicke)

ok_schw = True
for m in range(2, 8):
    Bs = basis(m)
    Xps = opsum([(sp.sqrt(k + 1), k + 1, k) for k in MODES if k >= 0])
    Xhds = opsum([(sp.sqrt(-k - 1), k, k + 1) for k in MODES if k <= -2])
    Ep = comp(Xps, Xhds)
    Em_ = comp(opsum([(sp.sqrt(-k - 1), k + 1, k) for k in MODES if k <= -2]),
               opsum([(sp.sqrt(k + 1), k, k + 1) for k in MODES if k >= 0]))
    Mp = sp.Matrix(compress(Ep, Bs))
    Mm = sp.Matrix(compress(Em_, Bs))
    Hm = sp.diag(*[sp.Rational(2*kr + m + 1, 2) for kr in range(-m, 0)])
    ok_schw &= leak(Ep, Bs) == 0 and sp.simplify(Mp*Mm - Mm*Mp - 2*Hm) == sp.zeros(m, m)
    ok_schw &= sp.simplify(Hm*Mp - Mp*Hm - Mp) == sp.zeros(m, m) and Mm == Mp.T
    ok_schw &= sp.simplify((Hm**2 + (Mp*Mm + Mm*Mp)/2) - sp.Rational((m - 1)*(m + 1), 4)*sp.eye(m)) == sp.zeros(m, m)
check("[L] Erst mit eingesetzten Schwinger-Gewichten, E_+ = (sum sqrt(n_1+1) psi*_{k+1} psi_k)(sum sqrt(n_2) psi*_k psi_{k+1}) "
      "(Teilchen- und Lochleiter als Oszillatorleitern gelesen), schließt das Hüpfen für alle m = 2..7 exakt zu su(2) mit "
      "Casimir j(j+1), j = (m-1)/2", ok_schw)
befund("Teil 2: Native schließende Operatoren gibt es nicht. (E) Bilineare sind auf H_m diagonal; die natürlichen quartischen "
       "Heber (Strom, Möbius) schließen nur für m <= 3 bzw. m = 2 (wie beim gleichförmigen Hüpfen ohne perfekten "
       "Zustandstransfer ab N = 4). Allgemein scheitert jede feste Polynomgewichtung vom Grad d ab m = 2d + 2, denn h^2 = "
       "c (n+1)(m-1-n) an m - 1 > 2d Punkten erzwingt eine Polynomidentität mit einfachen Nullstellen rechts. (E) Die einzige "
       "native su(2) aus Bilinearen ist die Dicke-su(2) der Paarbesetzung, und die steht quer zu H_m. (L) su(2) mit Spin "
       "(m-1)/2 verlangt die Gewichte sqrt((j-w)(j+w+1)) = sqrt(n_1+1) sqrt(n_2), also eine Oszillatorlesart der Teilchen- und "
       "Lochleitern, die das Original nicht liefert: das su(2) ist eine Lesart", True)

# ================================================================== Teil 3
say("\nTeil 3: Tomita-Verdopplung, 4D-Charakter und die thermischen Gewichte")
q = sp.symbols("q")
ser = sp.series(q*(1 + q)/(1 - q)**3, q, 0, 41).removeO()
ok_ser = all(ser.coeff(q, m) == m**2 for m in range(1, 41)) and ser.coeff(q, 0) == 0
check("[E] sum_{m>=1} m^2 q^m = q(1+q)/(1-q)^3 (Koeffizienten bis q^40 exakt)", ok_ser)

u_, v_ = sp.symbols("u v")
X4 = [u_*v_, u_/v_, v_/u_, 1/(u_*v_)]


def hcomp(n):
    if n < 0:
        return sp.Integer(0)
    return sum(sp.prod([X4[i] for i in c]) for c in itertools.combinations_with_replacement(range(4), n))


def chis(n, z):
    return sum(z**(n - 2*k) for k in range(n + 1))


ok_kt = True
x1, x2, x3, x4 = sp.symbols("x1:5", real=True)
r2 = x1**2 + x2**2 + x3**2 + x4**2


def lap4(f):
    return sum(sp.diff(f, xx, 2) for xx in (x1, x2, x3, x4))


harm = {1: x1, 2: x1**2 - x2**2 + 2*x3*x4, 3: x1*x2*x3 + x4**3 - 3*x4*x1**2, 4: (x1 + sp.I*x2)**4}
for n in range(0, 7):
    ok_kt &= sp.simplify(sp.expand(hcomp(n) - hcomp(n - 2) - chis(n, u_)*chis(n, v_))) == 0
    ok_kt &= sp.binomial(n + 3, 3) - (sp.binomial(n + 1, 3) if n >= 2 else 0) == (n + 1)**2
    ok_kt &= sp.sqrt(n*(n + 2) + 1) == n + 1
for n, f in harm.items():                             # Delta_{S^3} Y = -n(n+2) Y über die Gradnull-Fortsetzung f/r^n
    g = f/r2**sp.Rational(n, 2)
    ok_kt &= lap4(f) == 0 and sp.simplify(lap4(g) + n*(n + 2)*f/r2**sp.Rational(n + 2, 2)) == 0
qq = mp.mpf("0.3")
for a, b in ((mp.mpf("0.4"), mp.mpf("1.3")), (mp.mpf("2.2"), mp.mpf("0.9"))):
    xsn = [mp.expj(a + b), mp.expj(a - b), mp.expj(b - a), mp.expj(-a - b)]
    lhs = mp.fsum(chi(n + 1, a)*chi(n + 1, b)*qq**(n + 1) for n in range(0, 200))
    rhs = qq*(1 - qq**2)/mp.fprod(1 - qq*xx for xx in xsn)
    ok_kt &= abs(lhs - rhs) < mp.mpf(10)**-40
check("[E] 4D-Skalar auf R x S^3: harmonische Polynome vom Grad n auf R^4 haben den Toruscharakter chi_{n/2}(a) chi_{n/2}(b) "
      "(n = 0..6, sympy), Dimension (n+1)^2, Delta_{S^3} = -n(n+2) und konforme Energie sqrt(n(n+2) + 1) = n + 1: K-Typen (j, j) "
      "bei 2j + 1; Erzeugende sum chi chi q^(n+1) = q(1 - q^2)/prod(1 - q x_i)", ok_kt)

dev_t = mp.mpf(0)
for al_ in ALPHAS:
    for a, b in ((mp.mpf("0.4"), mp.mpf("1.3")), (mp.mpf("2.2"), mp.mpf("0.9"))):
        Da, Db = D_all(a, al_, B4, 120), D_all(b, al_, B4, 120)
        xsn = [mp.expj(a + b), mp.expj(a - b), mp.expj(b - a), mp.expj(-a - b)]
        lhs = mp.fsum(qq**m*(Da[m - 1]/Da[0])*(Db[m - 1]/Db[0]) for m in range(1, 121))
        rhs = qq*(1 - qq**2)/mp.fprod(1 - qq*xx for xx in xsn)
        dev_t = max(dev_t, abs(lhs - rhs))
ZAHLEN["torus_charakter_abweichung"] = mp.nstr(dev_t, 5)
check("[E] die verdoppelte native Auslese des Originals (beta = 4 pi, beide Sektoren), sum_m q^m (D_m/D_1)(a) (D_m/D_1)(b), ist "
      f"exakt der volle Toruscharakter des 4D-Skalars q(1 - q^2)/prod(1 - q x_i) (Abweichung {mp.nstr(dev_t, 3)}): "
      "Übereinstimmung auf Torusebene (U(1) x U(1))", dev_t < mp.mpf(10)**-35)


def spin_mats(m):
    j = (m - 1)/2
    ws = [-j + i for i in range(m)]
    Jp = np.zeros((m, m))
    for i in range(m - 1):
        Jp[i + 1, i] = np.sqrt((j - ws[i])*(j + ws[i] + 1))
    return [(Jp + Jp.T)/2, (Jp - Jp.T)/(2j), np.diag(ws)], j


ok_std = True
for m in range(2, 6):
    T, j = spin_mats(m)
    wr = rng.uniform(0.2, 1.0, m)
    for rho in (np.eye(m)/m, np.diag(wr/wr.sum())):
        xi = np.sqrt(rho)                                                              # rho diagonal
        ok_std &= np.allclose(xi.conj().T, xi)                                         # J xi = xi
        for _ in range(5):
            A = rng.normal(size=(m, m)) + 1j*rng.normal(size=(m, m))
            S_A = (np.sqrt(rho) @ (A @ xi) @ np.linalg.inv(np.sqrt(rho))).conj().T     # J Delta^(1/2) (A xi), J x = x^+
            ok_std &= np.allclose(S_A, A.conj().T @ xi)                                # = A^+ xi (Tomita S)
    Lc = sum(np.kron(t, np.eye(m)) @ np.kron(t, np.eye(m)) for t in T)                 # vec(x) zeilenweise: x -> t x
    Rc = sum(np.kron(np.eye(m), t.T) @ np.kron(np.eye(m), t.T) for t in T)             # x -> x t
    ok_std &= np.allclose(Lc, j*(j + 1)*np.eye(m*m)) and np.allclose(Rc, j*(j + 1)*np.eye(m*m))
    ok_std &= all(np.allclose((np.kron(t, np.eye(m)) - np.kron(np.eye(m), t.T)) @ np.eye(m).ravel(), 0) for t in T)
check("[E; su(2) aus der Lesart] Standardform eines treuen Zustands auf End(C^m) (Spurzustand und Nicht-Spurzustand, m = 2..5): "
      "m^2-dimensional, S = J Delta^(1/2) mit J x = x^+, J bildet die Links- auf die Rechtswirkung ab; mit dem su(2) der Lesart "
      "haben Links- und Rechtscasimir j(j+1), also (j) x (j-quer) = K-Typ (j, j), und die Einheit ist der Diagonalsinglett "
      "(Zonalvektor)", ok_std)

tq = sp.symbols("t", positive=True)                  # t = e^{-beta/4}: e^{beta x} = t^{-4x}, ganzzahlige Potenzen


def c_t(x4):                                          # x4 = 4x (ganzzahlig)
    return 1/(1 + tq**(-x4))


ok_kms = True
for m in range(1, 7):
    for a4 in (1, 3):                                 # alpha = a4/4
        wts = []
        for kr in range(-m, 0):
            r4, p4 = 4*kr + a4, 4*(kr + m) + a4
            wts.append(sp.simplify((1 - c_t(r4))*c_t(p4)/(c_t(r4)*(1 - c_t(p4)))))
        ok_kms &= all(sp.simplify(wv - tq**(4*m)) == 0 for wv in wts)
check("[E] die Kompression des Originalzustands (quasifrei, KMS) auf H_m ist exakt der Spurzustand: alle m Paarzustände haben "
      "relativ zum Vakuum das Gewicht e^{-beta m} (sympy, beliebiges beta, m = 1..6, beide Sektoren); der TFD-Anteil in "
      "H_m x H_m-quer ist der maximal verschränkte Vektor, Delta = 1 auf End(H_m): innerhalb einer Stufe keine Modulzeit", ok_kms)

qe = mp.exp(-mp.pi)
kubo = {}
ok_kubo = True
for al_ in ALPHAS:
    for m in range(1, 7):
        R = 40
        ks = range(-R - m, R)
        K = {k: occ(k + al_)*(1 - occ(k + al_ + m)) for k in ks}
        tot = mp.fsum(K.values())
        ok_kubo &= abs(tot - m/(1 - mp.exp(-B4*m))) < mp.mpf(10)**-40
        p = {k: K[k]/tot for k in ks}
        win = list(range(-m, 0))
        td = (mp.fsum(abs(p[k] - mp.mpf(1)/m) for k in win) + mp.fsum(p[k] for k in ks if k not in win))/2
        wv = [K[k] for k in win]
        kubo[(float(al_), m)] = {"fenster_gewichte": [mp.nstr(x, 8) for x in wv],
                                 "ausserhalb": mp.nstr(mp.fsum(K[k] for k in ks if k not in win), 8),
                                 "spurabstand": mp.nstr(td, 6)}
        rand = wv[0] if al_ == ALPHAS[0] else wv[-1]           # Fensterkante an der Fermikante (Hardy-Marke)
        if m >= 2:
            ok_kubo &= abs(rand - 1/(1 + qe)) < mp.mpf("2e-4")
        ok_kubo &= abs(mp.fsum(K[k] for k in ks if k not in win) - qe/(1 + qe)) < mp.mpf("2e-4")
q_ = sp.symbols("q", positive=True)
edge = {m: f"(1+q^3)(1+q^{4*m - 3})/((1+q)(1+q^{4*m - 1})) = "
           + mp.nstr((1 + qe**3)*(1 + qe**(4*m - 3))/((1 + qe)*(1 + qe**(4*m - 1))), 10) for m in (2, 3, 4)}
ZAHLEN["kubo_gewichte"] = {f"alpha={k[0]}, m={k[1]}": v for k, v in kubo.items()}
ZAHLEN["randverhaeltnis_unten_oben_alpha1/4 (q = e^-pi)"] = {m: str(v) for m, v in edge.items()}
ZAHLEN["rho1_gewichte"] = [mp.nstr(1/(1 + qe), 8), mp.nstr(qe/(1 + qe), 8)]
say(f"   Spurabstand der Kubo-Gewichte zum Spurzustand (alpha = 1/4): "
    f"{[kubo[(0.25, m)]['spurabstand'] for m in range(1, 7)]}")
check("[E] die Kubo-Gewichte c_r(1 - c_{r+m}) der Paare (die Gewichte von D_m) sind der Spurzustand bis auf Kanteneffekte: die "
      "Paarkante an der Fermikante trägt 1/(1 + e^-pi) = 95,8576 %, außerhalb des Fensters liegen e^-pi/(1 + e^-pi) = 4,1424 % "
      "(die rho_1-Gewichte aus z63), Summe exakt m/(1 - e^{-4 pi m}); Spurabstand ~ 0,04/m", ok_kubo)
befund("Teil 3 (E): Die Zählung stimmt: m^2 = dim End(H_m), sum m^2 q^m ist der Einteilchencharakter des 4D-Skalars, und die "
       "verdoppelte Auslese gibt sogar den vollen Toruscharakter. Der Originalzustand auf H_m ist exakt der Spurzustand (alle "
       "Paare entartet in L_0); seine Tomita-Verdopplung hat Delta = 1. Die Einheit von End(H_m) entspricht in der su(2)-Lesart "
       "der Zonalfunktion chi_j, also gerade der bilokalen Auslese. (j, j) als SU(2) x SU(2)-Darstellung setzt das su(2) aus "
       "Teil 2 voraus (L)", True)

# ================================================================== Teil 4
say("\nTeil 4: SU(2,2) auf V + V-quer* und Vergleich mit dem nativen Inhalt")
Sig = np.block([[np.zeros((2, 2)), np.eye(2)], [np.eye(2), np.zeros((2, 2))]])
sig_ev = np.sort(np.linalg.eigvalsh(Sig))
basisX = []
for i in range(4):
    for k in range(4):
        for ph in (1, 1j):
            E_ = np.zeros((4, 4), dtype=complex)
            E_[i, k] = ph
            basisX.append(E_)
rows = []
for E_ in basisX:
    cond = E_.conj().T @ Sig + Sig @ E_
    rows.append(np.concatenate([cond.real.ravel(), cond.imag.ravel(), [np.trace(E_).real, np.trace(E_).imag]]))
dim_su22 = 32 - np.linalg.matrix_rank(np.array(rows).T, tol=1e-10)
ok_tw = np.allclose(sig_ev, [-1, -1, 1, 1]) and dim_su22 == 15
for _ in range(20):
    xh = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
    xh = (xh + xh.conj().T)/2
    yh = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
    yh = (yh + yh.conj().T)/2
    Tx = np.block([[np.zeros((2, 2)), 1j*xh], [np.zeros((2, 2)), np.zeros((2, 2))]])
    Ty = np.block([[np.zeros((2, 2)), 1j*yh], [np.zeros((2, 2)), np.zeros((2, 2))]])
    ok_tw &= np.allclose(Tx.conj().T @ Sig + Sig @ Tx, 0) and np.allclose(Tx @ Ty - Ty @ Tx, 0)
    pi0 = rng.normal(size=2) + 1j*rng.normal(size=2)
    Zx = np.concatenate([1j*xh @ pi0, pi0])                     # Punkt der Geraden L_x: omega = i x pi
    ok_tw &= abs(Zx.conj() @ Sig @ Zx) < 1e-12
    om0 = rng.normal(size=2) + 1j*rng.normal(size=2)
    Zs = np.concatenate([om0, np.zeros(2)])                     # Naht P(V) = {pi = 0}
    ok_tw &= abs(Zs.conj() @ Sig @ Zs) < 1e-12 and np.allclose(Tx @ Zs, 0)
check("[E] Twistorraum V + V-quer* mit Sigma((omega, pi)) = omega.pi-quer + pi.omega-quer: Signatur (2,2), su(2,2) hat Dimension "
      "15, die Translationen omega -> omega + i x pi (x in Herm(2)) sind abelsch und liegen in su(2,2), reelle Punkte geben "
      "Nullgeraden; die Naht P(V) = {pi = 0} ist eine Nullgerade (Punkt im Unendlichen) und translationsfest", ok_tw)

Lt = 7
dd = Lt + 1
a1 = sps.diags(np.sqrt(np.arange(1, dd)), 1, format="csr")
I1 = sps.identity(dd, format="csr")


def emb(op, pos):
    mats = [I1]*4
    mats[pos] = op
    out = mats[0]
    for Mm in mats[1:]:
        out = sps.kron(out, Mm, format="csr")
    return out


Aop = [emb(a1, 0), emb(a1, 1)]
Bop = [emb(a1, 2), emb(a1, 3)]
Ad = [x.T.tocsr() for x in Aop]
Bd = [x.T.tocsr() for x in Bop]
dim = dd**4
occn = np.array(np.unravel_index(np.arange(dim), (dd, dd, dd, dd))).T
Na, Nb = occn[:, 0] + occn[:, 1], occn[:, 2] + occn[:, 3]
s2 = 1/np.sqrt(2)
omega = [s2*(Aop[i] - Bd[i]) for i in range(2)]
pi_ = [s2*(Aop[i] + Bd[i]) for i in range(2)]
pib = [s2*(Ad[i] + Bop[i]) for i in range(2)]
omb = [s2*(Ad[i] - Bop[i]) for i in range(2)]
Z, Zb = omega + pi_, pib + omb
in1 = np.where(occn.max(axis=1) <= Lt - 1)[0]
in2 = np.where(occn.max(axis=1) <= Lt - 2)[0]


def mx(Mx, idx):
    sub = Mx[:, idx]
    return 0.0 if sub.nnz == 0 else float(np.abs(sub.data).max())


ok_ccr = True
for a_ in range(4):
    for b_ in range(4):
        ok_ccr &= mx(Z[a_] @ Zb[b_] - Zb[b_] @ Z[a_] - (1.0 if a_ == b_ else 0.0)*sps.identity(dim, format="csr"), in1) < 1e-12
        ok_ccr &= mx(Z[a_] @ Z[b_] - Z[b_] @ Z[a_], in1) < 1e-12 and mx(Zb[a_] @ Zb[b_] - Zb[b_] @ Zb[a_], in1) < 1e-12
G =[[Z[a_] @ Zb[b_] for b_ in range(4)] for a_ in range(4)]
ok_u22 = True
zero = sps.csr_matrix((dim, dim))
for a_, b_, c_, d_ in itertools.product(range(4), repeat=4):
    Cm = G[a_][b_] @ G[c_][d_] - G[c_][d_] @ G[a_][b_]
    Tg = (G[c_][b_] if a_ == d_ else zero) - (G[a_][d_] if b_ == c_ else zero)
    ok_u22 &= mx(Cm - Tg, in2) < 1e-10
check("[E] Oszillatorquantisierung Z = (omega, pi) = ((a - b^+)/sqrt2, (a + b^+)/sqrt2), Z-quer = Sigma Z^+: [Z^alpha, Z-quer_beta] "
      "= delta, alle übrigen Kommutatoren 0, und die 16 Bilineare Z^alpha Z-quer_beta schließen zu u(2,2) (alle 256 Kommutatoren, "
      "abgeschnittener Fockraum 8^4, Innenbereich)", ok_ccr and ok_u22)

P = [[pib[A_] @ pi_[Ap] for Ap in range(2)] for A_ in range(2)]
Kc = [[omega[A_] @ omb[Ap] for Ap in range(2)] for A_ in range(2)]
pauli = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
Pmu = [sum(0.5*s[Ap, A_]*P[A_][Ap] for A_ in range(2) for Ap in range(2)) for s in pauli]
h0in = np.where((Na == Nb) & (occn.max(axis=1) <= Lt - 2))[0]
ok_pp = all(mx(Pmu[i] @ Pmu[k] - Pmu[k] @ Pmu[i], h0in) < 1e-11 for i in range(4) for k in range(4))
detP = mx(P[0][0] @ P[1][1] - P[0][1] @ P[1][0], h0in)
P2 = mx(Pmu[0] @ Pmu[0] - Pmu[1] @ Pmu[1] - Pmu[2] @ Pmu[2] - Pmu[3] @ Pmu[3], h0in)
sub1 = np.where((Na == Nb) & (occn.max(axis=1) <= Lt - 1))[0]
P0c = Pmu[0][sub1][:, sub1].toarray()
p0min = float(np.linalg.eigvalsh((P0c + P0c.conj().T)/2).min())
Hel = sps.diags((Na - Nb).astype(float))
hel_ok = max(mx(P[i][k] @ Hel - Hel @ P[i][k], np.arange(dim)) for i in range(2) for k in range(2)) < 1e-12
Econf = sps.diags(0.5*(Na + Nb) + 1.0)
E_ok = mx(Pmu[0] + 0.5*(Kc[0][0] + Kc[1][1]) - Econf, in1) < 1e-12
EP = mx(Econf @ Pmu[0] - Pmu[0] @ Econf, in1)
ZAHLEN["P_tests"] = {"det_P": detP, "P^2": P2, "P0_min_eig_kompression": p0min, "[E, P0]_max": EP}
check(f"[E] P_(AA') = pi-quer_A pi_(A'): [P_mu, P_nu] = 0, det P = P^2 = 0 ({detP:.1e}, {P2:.1e}), P^0 = (1/2) sum pi_A^+ pi_A "
      f">= 0 (Kompression im Helizität-0-Sektor, kleinster Eigenwert {p0min:.3f}), P erhält die Helizität (N_a - N_b)/2, und "
      f"E = P^0 + K^0 = (N_a + N_b)/2 + 1 ist der konforme Hamiltonoperator; [E, P^0] != 0",
      ok_pp and detP < 1e-11 and P2 < 1e-11 and p0min > -1e-12 and hel_ok and E_ok and EP > 0.1)

JL = [0.5*sum(s[i, k]*Ad[i] @ Aop[k] for i in range(2) for k in range(2)) for s in pauli[1:]]
JR = [0.5*sum(s[i, k]*Bd[i] @ Bop[k] for i in range(2) for k in range(2)) for s in pauli[1:]]
CL, CR = sum(x @ x for x in JL), sum(x @ x for x in JR)
ktab, ok_k = {}, True
for n in range(0, Lt + 1):
    idx = np.where((Na == n) & (Nb == n))[0]
    jj = n/2
    okn = (len(idx) == (n + 1)**2 and np.allclose(CL[idx][:, idx].toarray(), jj*(jj + 1)*np.eye(len(idx)))
           and np.allclose(CR[idx][:, idx].toarray(), jj*(jj + 1)*np.eye(len(idx)))
           and np.allclose(Econf.diagonal()[idx], n + 1))
    ktab[n] = (len(idx), bool(okn))
    ok_k &= okn
ZAHLEN["K_typen_helizitaet0"] = {n: {"dim": v[0], "(j,j) bei E = n+1": v[1]} for n, v in ktab.items()}
check("[E] K-Typen im Helizität-0-Sektor (N_a = N_b = n <= 7): E = n + 1 mit Vielfachheit (n+1)^2, Casimir links = rechts = "
      "j(j+1) mit j = n/2, also (j, j) bei 2j + 1: Charakter sum (n+1)^2 q^(n+1) = q(1+q)/(1-q)^3 wie in Teil 3", ok_k)

# native Heber auf dem Paarraum und ihr Tomita-Partner
Mlev = 7
pbasis = [(m, kr) for m in range(1, Mlev + 1) for kr in range(-m, 0)]
pvec = {pk: bil(pk[1] + pk[0], pk[1], {VAC: 1}) for pk in pbasis}
NP = len(pbasis)


def pmat(f):
    Mx = np.zeros((NP, NP))
    for j_, pj in enumerate(pbasis):
        img = f(pvec[pj])
        for i_, pi_k in enumerate(pbasis):
            Mx[i_, j_] = inner(pvec[pi_k], img)
    return Mx


XP, XH = pmat(Xp), pmat(Xh)
lev = np.array([pk[0] for pk in pbasis])
n1 = np.array([pk[0] + pk[1] for pk in pbasis])                  # Teilchenstufe k_p
inner1 = lev <= Mlev - 1
inner2 = lev <= Mlev - 2
P0n = np.diag((n1 == 0).astype(float))
iso = (np.allclose((XP.T @ XP)[np.ix_(inner1, inner1)], np.eye(inner1.sum()))
       and np.allclose((XP @ XP.T)[np.ix_(inner1, inner1)], (np.eye(NP) - P0n)[np.ix_(inner1, inner1)])
       and np.allclose((XP @ XH - XH @ XP)[:, inner2], 0))
Xs = [XP, XH]
R = [[np.kron(Xs[A_], Xs[Ap]) for Ap in range(2)] for A_ in range(2)]
in2d = np.kron(inner2, inner2).astype(bool)
in1d = np.kron(inner1, inner1).astype(bool)
eqlev = np.kron(lev[:, None], np.ones(NP)).ravel() == np.kron(np.ones(NP), lev[None, :]).ravel()
ok_R = all(np.allclose((R[a][b] @ R[c][d] - R[c][d] @ R[a][b])[:, in2d], 0)
           for a, b, c, d in itertools.product(range(2), repeat=4))
ok_R &= np.allclose((R[0][0] @ R[1][1] - R[0][1] @ R[1][0])[:, in2d], 0)
Cn = (R[0][0] @ R[0][0].T - R[0][0].T @ R[0][0])
tgtn = -(np.kron(P0n, np.eye(NP)) + np.kron(np.eye(NP), P0n) - np.kron(P0n, P0n))
sel = in1d & eqlev
ok_R &= np.allclose(Cn[np.ix_(sel, sel)], tgtn[np.ix_(sel, sel)])
evn = sorted(set(np.round(np.diag(Cn)[sel], 9)))
osc = -(occn[:, 0] + occn[:, 2] + 1.0)
Cosc = (Ad[0] @ Bd[0]) @ (Aop[0] @ Bop[0]) - (Aop[0] @ Bop[0]) @ (Ad[0] @ Bd[0])
ok_osc = mx(Cosc - sps.diags(osc, format="csr"), in1) < 1e-12
evo = sorted(set(osc[np.where((Na == Nb) & (occn.max(axis=1) <= Lt - 1))[0]]))
ZAHLEN["vergleich_hebekommutator"] = {"nativ_[R,R^+]_eigenwerte": [float(x) for x in evn],
                                      "oszillator_[E+,E-]_eigenwerte": [float(x) for x in evo]}
check("[E] native Heber X_p (Teilchen hoch), X_h (Loch tiefer) sind die Hardy-Hälften von J_{-1}: vertauschende Isometrien mit "
      "X^+ X = 1, X X^+ = 1 - P_0 (Toeplitz, q = 0). R_(AA') = X_A (x) X-quer_(A') auf dem verdoppelten Paarraum vertauschen, "
      "det R = 0, und R hebt m und m-quer um 1 wie a_A^+ b_(A')^+; aber [R_11', R_11'^+] = -(P_0 (x) 1 + 1 (x) P_0 - P_0 (x) P_0) "
      f"hat die Eigenwerte {[float(x) for x in evn]}, der Oszillator dagegen -(N_a1 + N_b1 + 1) mit {[float(x) for x in evo][:4]}...: "
      "die nativen Heber schließen nicht zu su(2,2)", iso and ok_R and ok_osc and evn == [-1.0, 0.0] and max(evo) <= -1)
befund("Teil 4 (E): Die 4D-Leiterdarstellung ist auf dem verdoppelten Paarraum realisierbar: gleiche Graduierung, gleiche "
       "Vielfachheiten. Mit Naht = {pi = 0} ist P bilinear im Partnerdublett pi aus V-quer* (nach Inversion vertauschen P "
       "und K). Die Translationen brauchen dafür die CCR-Normierung "
       "a^+ = X sqrt(N+1), und die Heber des Originals sind Toeplitz-Isometrien (q = 0). In der 4D-Lesart ist L_0 des "
       "Originals der konforme Hamiltonoperator E = P^0 + K^0 (ganzzahlig, periodisch, verträglich mit z70); P^0 selbst "
       "verlangt die nicht nativen nichtkompakten Erzeuger", True)

# ================================================================== Teil 5
say("\nTeil 5: Bilanz")
ZAHLEN["status"] = {
    "E": ["Paarraum H_m (dim m), J_{-m} Omega = Summe der m Paare, L_0 = m",
          "bilokale Faktorisierung D_m = D_1 chi_j (Teleskopsumme, für jedes Profil; nicht zustandsspezifisch)",
          "normiertes W erfüllt die S^3-Identität; äquivalent: sin(gamma) W ist 2D-harmonisch, d_gamma(sin W) = W2 (z71)",
          "native U(1)-Gewichte w - alpha = -j..j (Cartan) und voller Toruscharakter der verdoppelten Auslese",
          "No-go: jede Bilinearform komprimiert auf H_m diagonal",
          "Strom- bzw. Möbiushüpfen schließt nur für m <= 3 bzw. m = 2; keine feste Polynomgewichtung schließt für alle m "
          "(Grad <= 2: nie bei m = 6); Dicke-su(2) steht quer zu H_m",
          "KMS-Kompression auf H_m = Spurzustand exakt; Kubo-Gewichte = Spurzustand bis auf rho_1-Kanten (4,1424 %)",
          "SU(2,2) auf V + V-quer*: [P, P] = 0, P^0 >= 0, P^2 = 0, K-Typen (j, j) bei 2j + 1, Charakter = sum m^2 q^m",
          "native Heber: Toeplitz (X^+X = 1), [R, R^+] = Randprojektoren statt -(N + N + 1)"],
    "C": ["4D-Lesart unter T' (Naht = P(V) = {pi = 0} in P(V + V-quer*), Gerade des Punkts im Unendlichen): Translationen "
          "bilinear im Partnerdublett (nach Inversion vertauschen P und K)"],
    "L": ["su(2) mit Spin (m-1)/2 auf H_m (Schwinger-/Oszillatorgewichte eingesetzt)",
          "(j, j) als SU(2) x SU(2)-K-Typ und die S^3-Gruppenwirkung",
          "su(2,2)-Translationen mit CCR-Normierung auf dem verdoppelten Paarraum"],
    "O": ["ein Prinzip des Originals, das die Oszillatornormierung sqrt(n+1) der Teilchen- und Lochleitern erzwingt",
          "Nichtzonale S^3-Daten (zwei weitere Winkel) aus Antworten des Originals"],
}
for k_, v_ in ZAHLEN["status"].items():
    for item in v_:
        say(f"   {k_}: {item}")
befund("Urteil: Die 4D-Struktur ist im Original auf der Ebene von Zählung und Torus nativ (E), nicht auf der Ebene der Gruppe. "
       "Das su(2) auf dem Paarraum und damit SU(2) x SU(2) und SU(2,2) sind Lesarten (L). Die masselose 4D-Darstellung ist "
       "mit dem Original verträglich, aber nicht aus ihm rekonstruiert", True)

urteil = ("Die externe Behauptung hält auf der Ebene von Zählung und Torus, nicht auf der Ebene der Gruppe. Exakt und nativ "
          "sind: der m-dimensionale Paarraum H_m von J_{-m} Omega mit den bilokalen Gewichten -j..j; die Faktorisierung "
          "D_m = D_1 sin(m gamma)/sin(gamma), die aber eine Teleskopsumme ist und für jedes Besetzungsprofil gilt; die "
          "S^3-Identität des normierten W, die exakt einer 2D-Dirichlet-Antwort (sin(gamma) W harmonisch, d_gamma(sin W) = W2 aus "
          "z71) gleichwertig ist; die Kompression des Originalzustands auf H_m, die exakt der Spurzustand ist; und der volle "
          "Toruscharakter: sum m^2 q^m und die verdoppelte Auslese sind der Einteilchencharakter des konformen 4D-Skalars. "
          "Nicht nativ ist das su(2): Jede Bilinearform komprimiert auf H_m diagonal, das native Strom- und Möbiushüpfen "
          "schließt nur für m <= 3 bzw. m = 2, keine feste Polynomgewichtung schließt für alle m, und die Dicke-su(2) der "
          "Paare steht quer zu H_m. Spin (m-1)/2 für alle m verlangt die Oszillatorgewichte sqrt(n_1+1) sqrt(n_2). Ebenso sind die nativen Heber Toeplitz-Isometrien (q = 0) "
          "statt CCR-Oszillatoren, und ihre Kommutatoren sind Randprojektoren statt u(2) + u(2). Die SU(2,2)-Leiterdarstellung "
          "auf V + V-quer* ist korrekt ([P, P] = 0, P^0 >= 0, P^2 = 0, K-Typen (j, j) bei 2j + 1) und auf dem verdoppelten "
          "Paarraum realisierbar. Sie ist mit dem Original verträglich, aber nicht aus ihm rekonstruiert; die fehlende Zutat "
          "ist die Oszillatornormierung der Teilchen- und Lochleitern.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z73_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
