"""Z95: Prüfung der elften externen Analyse 'Kompositionsgesetze und Double Cover' (25.11.).

Teil 1  Klassische Bitgesetze: 16 binäre Verknüpfungen, 6 assoziativ und kommutativ, 2 Gruppen (XOR, XNOR); von den 24
        reversiblen Zweibitoperationen vertragen sich mit der Uhr R (x) R nur Identität und Tausch.
Teil 2  Quantengesetze: unitäre Verknüpfungen, die R (x) R erhalten und tauschsymmetrisch sind, bilden bis auf Phase eine
        dreiparametrige Klasse (Uhr, Hüpfen, n_1 n_2); die Klasse U^4 = 1 hat 64 Gatter, davon 16 gaußsch (Matchgate);
        Yang-Baxter (Zopfrelation) wird gezählt.
Teil 3  Kollektive Nichtgaußizität identischer Saaten: Zustände rho(n) ~ y^(|n| + 4 f(|n|)), y = e^-pi, erhalten die Uhr
        exakt; gesucht sind nichttriviale f mit exakt erhaltenen lokalen Randverteilungen, für N = 2..8 Saaten.
Teil 4  Relative Paritätsphase W_12 = (-1)^(n_1 n_2) und die Uhr R~ = (R (x) R) W.
Teil 5  E8-Komposition: c((E8)_2) = 31/2, relativer Sektor c = 1/2, T_rel = (T_1 + T_2 - Summe J_1 J_2)/32,
        c_rel(N) = 8N(N-1)/(N+30); E8 + E8(-1): Diagonale isotrop vom Rang 8.
Teil 6  Double Cover: J R J = R^-1 auf dem Quadrattorus, Purifikation der Saat, Verbindungsoperator B, Twistgewicht c/16.
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Klassische Bitgesetze")
ops = list(itertools.product((0, 1), repeat=4))                   # Tabelle f(00), f(01), f(10), f(11)
f_ = lambda t, a, b: t[2*a + b]
assoz_komm = [t for t in ops if all(f_(t, a, b) == f_(t, b, a) for a in (0, 1) for b in (0, 1))
              and all(f_(t, f_(t, a, b), c) == f_(t, a, f_(t, b, c)) for a in (0, 1) for b in (0, 1) for c in (0, 1))]
gruppen = [t for t in assoz_komm if any(all(f_(t, e, a) == a for a in (0, 1)) for e in (0, 1))
           and all(len({f_(t, a, b) for b in (0, 1)}) == 2 for a in (0, 1))]
RR = np.diag([1, 1j, 1j, -1])                                     # R (x) R auf 00, 01, 10, 11
perm_ok = 0
for p in itertools.permutations(range(4)):
    P = np.zeros((4, 4))
    for i, j in enumerate(p):
        P[j, i] = 1
    perm_ok += np.allclose(P @ RR, RR @ P)
check(f"16 Verknüpfungen, {len(assoz_komm)} assoziativ und kommutativ, {len(gruppen)} Gruppen (XOR, XNOR); von 24 reversiblen "
      f"Zweibitoperationen vertauschen {perm_ok} mit R (x) R (Identität, Tausch)",
      len(ops) == 16 and len(assoz_komm) == 6 and len(gruppen) == 2 and perm_ok == 2)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Quantengesetze mit erhaltener Uhr und Tauschsymmetrie")
SW = np.eye(4)[[0, 2, 1, 3]]
basis = []
for k in range(16):
    for ph in (1, 1j):
        H = np.zeros(16, complex)
        H[k] = ph
        basis.append(H.reshape(4, 4))
rows = []
for H in basis:
    cond = np.concatenate([(H - H.conj().T).ravel(), (H @ RR - RR @ H).ravel(), (H @ SW - SW @ H).ravel()])
    rows.append(np.concatenate([cond.real, cond.imag]))
u, s_, vh = np.linalg.svd(np.array(rows).T)
dim = int(np.sum(s_ < 1e-9)) + max(0, len(basis) - len(s_))
check(f"hermitesche Erzeuger, die mit R (x) R und dem Tausch vertauschen: Dimension {dim} (= 3 + globale Phase: Uhr n_1 + n_2, "
      "Hüpfen, n_1 n_2), also die Blockform diag(a, [[b, c], [c, b]], d)", dim == 4)


def gatter(a, e1, e2, d):
    b, c = (e1 + e2)/2, (e1 - e2)/2
    return np.array([[a, 0, 0, 0], [0, b, c, 0], [0, c, b, 0], [0, 0, 0, d]])


mu4 = [1, 1j, -1, -1j]
klassen = {}
for a, e1, e2, d in itertools.product(range(4), repeat=4):
    key = min(tuple((x - s) % 4 for x in (a, e1, e2, d)) for s in range(4))
    klassen[key] = gatter(*[mu4[x] for x in key])
gauss = {k for k in klassen if (k[0] + k[3]) % 4 == (k[1] + k[2]) % 4}   # det A = det B (Matchgate)
I2 = np.eye(2)


def zopf(U):
    U12, U23 = np.kron(U, I2), np.kron(I2, U)
    return np.allclose(U12 @ U23 @ U12, U23 @ U12 @ U23, atol=1e-12)


yb = {k for k, U in klassen.items() if zopf(U)}
ZAHLEN["U4_klasse"] = {"gatter": len(klassen), "gaussisch": len(gauss), "zopfrelation": len(yb),
                       "zopf_nichtgaussisch": len(yb - gauss)}
say(f"   {ZAHLEN['U4_klasse']}")
check(f"U^4 = 1 bis auf Phase: {len(klassen)} Gatter, davon {len(gauss)} gaußsch (Matchgate det A = det B) und "
      f"{len(klassen) - len(gauss)} nichtgaußsch (die Analyse: 64, 16, 48)", len(klassen) == 64 and len(gauss) == 16)
check(f"Zopfrelation (Yang-Baxter in Zopfform) erfüllen {len(yb)} Gatter, davon {len(yb - gauss)} nichtgaußsch (die Analyse "
      "nennt 17 und 12)", len(yb) == 17 and len(yb - gauss) == 12)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Kollektive Nichtgaußizität identischer Saaten (Randverteilungen und Uhr exakt)")
Y = sp.Symbol("y", positive=True)


def loesungen(N, rng=range(-2, 3)):
    """f: {0..N} -> Z (f(0) = 0), Gewichte y^(k + 4 f(k)) mal C(N, k); Randbedingung als Multimengen-Identität."""
    sols = []
    for vals in itertools.product(rng, repeat=N):
        f = (0,) + vals
        if all(v == 0 for v in vals):
            continue
        links, rechts = {}, {}
        for k in range(N + 1):
            if k >= 1:
                e = k + 4*f[k]
                links[e] = links.get(e, 0) + sp.binomial(N - 1, k - 1)
            if k <= N - 1:
                e = k + 4*f[k] + 1
                rechts[e] = rechts.get(e, 0) + sp.binomial(N - 1, k)
        if {e: m for e, m in links.items() if m} == {e: m for e, m in rechts.items() if m}:
            sols.append(f)
    return sols


tab_sym = {N: loesungen(N) for N in range(2, 9)}
ZAHLEN["symmetrische_f_loesungen"] = {N: len(s_) for N, s_ in tab_sym.items()}
check("mit symmetrischen f (nur von der Teilchenzahl abhängig, |f| <= 2) gibt es für N = 2..8 keine nichttriviale Lösung: die "
      "kollektive Struktur braucht eine Zerlegung, die die Moden unterscheidet", all(len(s_) == 0 for s_ in tab_sym.values()))

from scipy.optimize import milp, LinearConstraint, Bounds


def allgemein(N, V=(-1, 0, 1)):
    """ganzzahlige Suche: f: {0,1}^N -> V, f(0) = 0, nichttrivial, Randverteilungen exakt (Multimengen)."""
    zust = list(itertools.product((0, 1), repeat=N))
    nV = len(V)
    nv = len(zust)*nV
    idx = lambda a, v: nV*a + V.index(v)
    A, lo, hi = [], [], []
    for a in range(len(zust)):                                   # eine Wahl je Zustand
        r = np.zeros(nv)
        r[nV*a:nV*a + nV] = 1
        A.append(r); lo.append(1); hi.append(1)
    r = np.zeros(nv); r[idx(0, 0)] = 1
    A.append(r); lo.append(1); hi.append(1)                      # f(0) = 0
    for i in range(N):
        for E in range(4*min(V) - 1, N + 4*max(V) + 2):
            r = np.zeros(nv)
            for a, n in enumerate(zust):
                k = sum(n)
                for v in V:
                    if n[i] == 1 and k + 4*v == E:
                        r[idx(a, v)] += 1
                    if n[i] == 0 and k + 4*v + 1 == E:
                        r[idx(a, v)] -= 1
            if np.any(r):
                A.append(r); lo.append(0); hi.append(0)
    r = np.zeros(nv)
    for a in range(len(zust)):
        for v in V:
            if v != 0:
                r[idx(a, v)] = 1
    A.append(r); lo.append(1); hi.append(np.inf)                 # nichttrivial
    res = milp(c=r, constraints=LinearConstraint(np.array(A), lo, hi), integrality=np.ones(nv), bounds=Bounds(0, 1))
    if res.status != 0:
        return None
    x = np.round(res.x).astype(int)
    return {zust[a]: [v for v in V if x[idx(a, v)]][0] for a in range(len(zust))}


eng = {N: allgemein(N) for N in range(2, 8)}
weit_V = tuple(range(-2, 7))
weit = {N: allgemein(N, weit_V) for N in range(2, 7)}
ZAHLEN["diagonal_eng_f_-1_bis_1"] = {N: (s_ is not None) for N, s_ in eng.items()}
ZAHLEN["diagonal_weit_f_-2_bis_6"] = {N: (s_ is not None) for N, s_ in weit.items()}
ok6 = False
if weit[6] is not None:
    f6 = weit[6]
    zust6 = list(f6)
    w6 = {n: Y**(sum(n) + 4*f6[n]) for n in zust6}
    Z6 = sum(w6.values())
    rand6 = all(sp.simplify(sum(w6[n] for n in zust6 if n[i])/Z6 - Y/(1 + Y)) == 0 for i in range(6))
    yv6 = float(sp.exp(-sp.pi))
    pv6 = {n: float(w6[n].subs(Y, yv6)) for n in zust6}
    Zv6 = sum(pv6.values())
    p1_6 = yv6/(1 + yv6)
    paare6 = [sum(pv6[n] for n in zust6 if n[i] and n[j])/Zv6 - p1_6**2 for i, j in itertools.combinations(range(6), 2)]
    S6 = -sum(v/Zv6*np.log(v/Zv6) for v in pv6.values())
    S6p = 6*(-(p1_6*np.log(p1_6) + (1 - p1_6)*np.log(1 - p1_6)))
    ok6 = rand6 and S6 < S6p - 1e-12 and sum(abs(c_) > 1e-12 for c_ in paare6) == 15
    ZAHLEN["sechs_saaten"] = {"f_bereich": [min(f6.values()), max(f6.values())],
                              "verschobene_zustaende": sum(1 for v in f6.values() if v), "S": S6, "S_Produkt": S6p}
    say(f"   sechs Saaten: {ZAHLEN['sechs_saaten']}")
check("Zustände diagonal in der Besetzungsbasis (exakte ganzzahlige Suche, Vakuum fest): mit f in {-1, 0, 1} keine Lösung für "
      "N = 2..7, mit f in [-2, 6] keine für N = 2..5, aber eine exakte für N = 6: alle Randverteilungen symbolisch gleich der "
      "Saat, Uhr exakt, alle 15 Paarkorrelationen ungleich null, Entropie kleiner als das Produkt. Die Schwelle N = 6 der "
      "Analyse ist bestätigt (Korrektur H6h46: der erste Suchraum war zu eng, nötig sind Verschiebungen bis f = 5, vgl. z94)",
      all(s_ is None for s_ in eng.values()) and all(weit[N] is None for N in range(2, 6)) and ok6)

from math import comb


def gleichmaessig(N, F=(-2, -1, 1, 2)):
    """teilchenzahlerhaltende Zustände mit gleichmäßigen Projektoren (jede Mode gleich besetzt), Dimensionen d_(k,f)."""
    keys = [(k, f) for k in range(1, N + 1) for f in F]
    rows = {}
    for j, (k, f) in enumerate(keys):
        for e, c in ((k + 4*f, k), (k + 4*f + 1, k - N), (k, -k), (k + 1, -(k - N))):
            rows.setdefault(e, np.zeros(len(keys)))[j] += c
    A = [r for r in rows.values() if np.any(r)]
    lo, hi = [0]*len(A), [0]*len(A)
    for k in range(1, N + 1):
        r = np.array([1.0 if kk == k else 0.0 for kk, _ in keys])
        A.append(r); lo.append(0); hi.append(comb(N, k))
    A.append(np.ones(len(keys))); lo.append(1); hi.append(np.inf)
    res = milp(c=np.ones(len(keys)), constraints=LinearConstraint(np.array(A), lo, hi),
               integrality=np.ones(len(keys)), bounds=Bounds(0, np.inf))
    if res.status != 0:
        return None
    d = np.round(res.x).astype(int)
    return {keys[j]: int(d[j]) for j in range(len(keys)) if d[j]}


tab_g = {N: gleichmaessig(N) for N in range(2, 10)}
erstes = min((N for N, s_ in tab_g.items() if s_ is not None), default=None)
ok_ex = False
if erstes:
    N, d = erstes, tab_g[erstes]
    stufen = {k: [(comb(N, k) - sum(v for (kk, f), v in d.items() if kk == k), 0)]
                 + [(v, f) for (kk, f), v in d.items() if kk == k] for k in range(N + 1)}
    Z = sum(m*Y**(k + 4*f) for k, lst in stufen.items() for m, f in lst)
    rand = sp.simplify(sum(sp.Rational(k, N)*m*Y**(k + 4*f) for k, lst in stufen.items() for m, f in lst)/Z - Y/(1 + Y))
    yv = float(sp.exp(-sp.pi))
    wv = [(m, k, float(yv**(k + 4*f))) for k, lst in stufen.items() for m, f in lst if m]
    Zv = sum(m*w for m, _, w in wv)
    S = -sum(m*(w/Zv)*np.log(w/Zv) for m, _, w in wv)
    p1 = yv/(1 + yv)
    S_prod = N*(-(p1*np.log(p1) + (1 - p1)*np.log(1 - p1)))
    ok_ex = rand == 0 and S < S_prod - 1e-12
    ZAHLEN["gleichmaessig_beispiel"] = {"N": N, "verschobene_dimensionen": {f"k={k},f={f}": v for (k, f), v in d.items()},
                                        "S": S, "S_Produkt": S_prod}
    say(f"   gleichmäßige Projektoren: erstes N = {erstes}, {ZAHLEN['gleichmaessig_beispiel']}")
lp = {}
for N in range(2, 8):
    zust = list(itertools.product((0, 1), repeat=N))
    nv = len(zust)*3
    idx = lambda a, v: 3*a + (v + 1)
    A, lo, hi = [], [], []
    for a in range(len(zust)):
        r = np.zeros(nv); r[3*a:3*a + 3] = 1; A.append(r); lo.append(1); hi.append(1)
    r = np.zeros(nv); r[idx(0, 0)] = 1; A.append(r); lo.append(1); hi.append(1)
    for i in range(N):
        for E in range(-5, N + 6):
            r = np.zeros(nv)
            for a, n in enumerate(zust):
                k = sum(n)
                for v in (-1, 0, 1):
                    if n[i] == 1 and k + 4*v == E:
                        r[idx(a, v)] += 1
                    if n[i] == 0 and k + 4*v + 1 == E:
                        r[idx(a, v)] -= 1
            if np.any(r):
                A.append(r); lo.append(0); hi.append(0)
    r = np.zeros(nv)
    for a in range(len(zust)):
        r[idx(a, -1)] = r[idx(a, 1)] = 1
    A.append(r); lo.append(1); hi.append(np.inf)
    res = milp(c=-r, constraints=LinearConstraint(np.array(A), lo, hi), integrality=np.zeros(nv), bounds=Bounds(0, 1))
    lp[N] = res.status == 0
ZAHLEN["lineare_lockerung_zulaessig"] = lp
check(f"teilchenzahlerhaltend mit gleichmäßigen Projektoren (realisierbar: jeder Sektor hat eine Orthonormalbasis aus "
      f"zyklischen Fourierzuständen mit gleicher Besetzung aller Moden): keine Lösung für N <= {erstes - 1 if erstes else '?'}, "
      f"exakte Lösung bei N = {erstes} (Randverteilungen symbolisch gleich der Saat, Uhr exakt, Entropie kleiner als das "
      "Produkt); die lineare Lockerung (notwendige Bedingung für beliebige teilchenzahlerhaltende Zustände) ist ab N = 5 "
      "zulässig, Realisierbarkeit dort nicht gezeigt", erstes == 9 and ok_ex and not lp[4] and lp[5])
befund("Die Analyse hält auch in der Zahl: Identische Saaten können ab N = 6 kollektiv nichtgaußsch sein, mit exakt freien "
       "Randverteilungen und exakter Uhr (Exponentenkollisionen modulo 4; diagonal mit Verschiebungen bis 5, gleichmäßig ab 9). "
       "Für das Original (verschiedene m, jedes ungerade m einmal) bleibt die Starrheit in allen gerechneten Fällen (z94): Die "
       "Nichtentartung trägt sie, nicht die Ladungsauswahl. Nativ ist ein solcher Zustand nicht (z91)", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Relative Paritätsphase")
W = np.diag([1, 1, 1, -1])
SWf = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, -1]])
Rt = RR @ W
check("W_12 = (-1)^(n_1 n_2) ist genau der Unterschied zwischen fermionischem und gewöhnlichem Tausch (SW_f = SW W); "
      "R~ = (R (x) R) W erfüllt R~^2 = (-1)^(n_1 + n_2) und R~^4 = 1",
      np.allclose(SWf, SW @ W) and np.allclose(Rt @ Rt, np.diag([1, -1, -1, 1])) and np.allclose(np.linalg.matrix_power(Rt, 4), np.eye(4)))

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: E8-Komposition")
k, h, dimg = sp.symbols("k h dimg", positive=True)
c_sug = lambda lev: sp.Rational(248)*lev/(lev + 30)
Ns = sp.Symbol("N", positive=True, integer=True)
c_rel_N = sp.simplify(8*Ns - 248*Ns/(Ns + 30) - 8*Ns*(Ns - 1)/(Ns + 30))
# T_i = (1/62) J_i^2, T_diag = (1/64)(J_1 + J_2)^2: T_1 + T_2 - T_diag = (1/32)(T_1 + T_2 - J_1.J_2), wenn J^2 = 62 T
a1 = sp.Rational(1, 62) - sp.Rational(1, 64)
ok_T = a1*62 == sp.Rational(1, 32) and sp.Rational(2, 64) == sp.Rational(1, 32)
check(f"c((E8)_1) = {c_sug(1)}, c((E8)_2) = {c_sug(2)}, relativer Sektor 16 - 31/2 = {16 - c_sug(2)} (Ising); "
      "T_rel = (T_1 + T_2 - Summe J_1^a J_2^a)/32 aus den Sugawara-Normierungen; c_rel(N) = 8N(N-1)/(N+30)",
      c_sug(1) == 8 and c_sug(2) == sp.Rational(31, 2) and 16 - c_sug(2) == sp.Rational(1, 2) and ok_T and c_rel_N == 0)
# E8-Wurzeln
wurzeln = []
for i, j in itertools.combinations(range(8), 2):
    for si in (1, -1):
        for sj in (1, -1):
            v = [0]*8
            v[i], v[j] = si, sj
            wurzeln.append(tuple(Fr(x) for x in v))
for signs in itertools.product((1, -1), repeat=8):
    if signs.count(-1) % 2 == 0:
        wurzeln.append(tuple(Fr(s, 2) for s in signs))
norm = lambda v: sum(x*x for x in v)
ok_null = len(wurzeln) == 240 and all(norm(v) == 2 for v in wurzeln) and all(norm(v) - norm(v) == 0 for v in wurzeln)
G = np.array([[float(sum(a*b for a, b in zip(u_, v_))) for v_ in wurzeln[:8]] for u_ in wurzeln[:8]])
basis8 = np.array([[float(x) for x in v] for v in wurzeln])
rang = np.linalg.matrix_rank(basis8)
check(f"E8: {len(wurzeln)} Wurzeln der Norm 2; in E8 + E8(-1) (Signatur (8,8)) hat jede Diagonale (v, v) die Norm 0, die "
      f"Diagonale ist isotrop vom Rang {rang} = 8, also maximal", ok_null and rang == 8)

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: Double Cover, Purifikation, Twistgewicht")
rng = np.random.default_rng(95)
mod = lambda z: complex(z.real % 1, z.imag % 1)
Rg = lambda z: 1j*z + 0.5
Rinv = lambda z: -1j*(z - 0.5)
Jg = lambda z: np.conj(z) + 0.5j
ok_t = True
for _ in range(256):
    z = complex(*rng.random(2))
    d_ = mod(Jg(Rg(Jg(z)))) - mod(Rinv(z))
    ok_t &= min(abs(d_ - (p + 1j*q_)) for p in (-1, 0, 1) for q_ in (-1, 0, 1)) < 1e-12
fix = [y_ for y_ in np.linspace(0, 1, 401)[:-1] if abs((2*y_ - 0.5) - round(2*y_ - 0.5)) < 1e-12]
yv = np.exp(-np.pi)
xx = 2*np.exp(-np.pi/2)/(1 + yv)
check(f"auf dem Quadrattorus gilt J R J = R^-1 (256 Punkte modulo Z + iZ), J hat die Fixkreise Im z = {fix}; Purifikation "
      f"(|00> + e^-pi/2 |11>)/sqrt(1 + e^-pi): <X (x) X> = {xx:.6f}; Twistgewicht c/24 (n - 1/n) bei n = 2, c = 8: "
      f"{sp.Rational(8, 24)*(2 - sp.Rational(1, 2))}",
      ok_t and fix == [0.25, 0.75] and abs(xx - 0.398537) < 5e-7 and sp.Rational(8, 24)*(2 - sp.Rational(1, 2)) == sp.Rational(1, 2))
befund("Die elfte Analyse hält in den gerechneten Punkten; die physikalischen Lesarten (E8-Coset als Kompositionsgesetz, Marken "
       "als Twistoperatoren, Double Cover als konjugiertes Paar) sind Kandidaten und gehen an F1, F3, F4 und F6", True)

urteil = ("Die elfte externe Analyse hält in den gerechneten Punkten. Von 16 Bitgesetzen sind 6 assoziativ und kommutativ, 2 sind "
          "Gruppen; von 24 reversiblen Zweibitoperationen vertragen sich nur Identität und Tausch mit der Uhr. Uhr- und "
          "tauscherhaltende Quantengesetze bilden eine dreiparametrige Klasse (Uhr, Hüpfen, n_1 n_2), mit U^4 = 1 gibt es 64, "
          "davon 16 gaußsch. Identische Saaten können ab N = 6 kollektiv nichtgaußsch sein (diagonal, Verschiebungen bis 5; die "
          "Schwelle der Analyse ist bestätigt, Korrektur H6h46); für das Original trägt die Nichtentartung die Starrheit (z94). "
          "E8-Komposition: relativer Sektor c = 1/2, "
          "T_rel = (T_1 + T_2 - Summe J_1 J_2)/32; E8 + E8(-1) hat eine isotrope Diagonale vom Rang 8. J R J = R^-1 auf dem "
          "Torus, Purifikation <X X> = 0,398537, Twistgewicht 1/2. Die physikalischen Lesarten bleiben Kandidaten.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z95_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
