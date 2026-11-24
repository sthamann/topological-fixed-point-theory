"""Z75: Nichtabelsche Mehrpunktantworten der Naht und der Spin-2-Kanal (nächste Hürde der externen Analyse).

Statuskürzel: E exakt, C bedingt, L Lesart, O offen. Original wie z58: c(x) = 1/(1 + e^{4 pi x}), x = r + alpha.

Teil 1  Nichtabelsche Drei- und Vierpunktantworten aus Wick (Fermionschleifen, gegen Fockraum geprüft): su(N) c u(8) c
        so(16)_1 c E8_1 im Originalzustand (acht komplexe Fermionen im Quarterrahmen alpha = 1/4, beta = 4 pi), im
        NS-Zustand und im Vakuum. <J J J> = f F, kein d-Anteil, Ward-Identität; verbundene Vierpunktantwort != 0
        (nichtabelsch), = 0 (abelsch). Vielfachheiten: Nullvektor (E^alpha_-1)^2|0> = 0, Stufe-2-Zählungen.
Teil 2  Externe Behauptung 'keine freien 4D-Felder': Delta = 1 in 4D erzwingt Box phi = 0 und damit verschwindende
        Dreipunktfunktion; Vertexoperatoren brauchen den 2D-Logarithmus; E8-Bilineare brauchen die 248 (Niveau 30).
        Was eine 4D-Lesart sein müsste: zusammengesetzte Felder.
Teil 3  Spin 2: Einteilchencharakter von T(z), Paarraum, Kern der linearisierten Weyl-Gleichung (Polynomlösungen) und
        4D-Charaktere; Linientest (Delta); Sym^2 des Skalarsingletons (Flato-Fronsdal in 4D).
Teil 4  Araki-Woods-Dubleton der externen Analyse (d_i, e_i aus zwei Quartermoden und Tomita-Partnern): Helizitätssektoren
        D = N_d - N_e, Charaktervergleich, was davon Dynamik des Originals ist.
Teil 5  Weinberg-Witten-Typisierung, Bilanz, was für Wechselwirkung und nichtlineare Gravitation bleibt.
"""
import functools
import itertools
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import sympy as sp
from sympy.polys.matrices import DomainMatrix

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


PI4 = 4*np.pi


def occ_factory(beta):
    """Besetzung c(x) der originalen Kovarianz bei inverser Temperatur beta (beta = inf: Vakuum)."""
    if np.isinf(beta):
        return lambda x: 1.0 if x < 0 else 0.0
    return lambda x: 1.0/(1.0 + np.exp(np.clip(beta*x, -700, 700)))


def key(x):
    return int(round(4*x))


# ------------------------------------------------------------------ Fockraum (Prüfstand für die Schleifenformel)
def fock_ops(nmodes):
    I2 = sps.identity(2, format="csr")
    Zp = sps.diags([1.0, -1.0], format="csr")
    sm = sps.csr_matrix(np.array([[0.0, 1.0], [0.0, 0.0]]))
    ops = []
    for j in range(nmodes):
        mats = [Zp]*j + [sm] + [I2]*(nmodes - j - 1)
        op = mats[0]
        for M in mats[1:]:
            op = sps.kron(op, M, format="csr")
        ops.append(op)
    return ops


class Fermi:
    """F komplexe Fermionen mit Moden xs; J(T)_m = sum_x c^+_x T c_{x+m} senkt die Energie um m."""

    def __init__(self, F, xs, occf):
        self.F, self.xs, self.M = F, list(xs), len(xs)
        self.c = fock_ops(F*self.M)
        self.dim = self.c[0].shape[0]
        self.num = [np.real((op.T @ op).diagonal()) for op in self.c]
        rho = np.ones(self.dim)
        for k, nn in enumerate(self.num):
            p = occf(self.xs[k % self.M])
            rho = rho*(p*nn + (1 - p)*(1 - nn))
        self.rho = rho
        self.idx = {key(x): j for j, x in enumerate(self.xs)}

    def J(self, T, m):
        op = sps.csr_matrix((self.dim, self.dim), dtype=complex)
        for j, x in enumerate(self.xs):
            k = self.idx.get(key(x + m))
            if k is None:
                continue
            for f in range(self.F):
                for g in range(self.F):
                    if T[f, g] != 0:
                        op = op + T[f, g]*(self.c[f*self.M + j].T @ self.c[g*self.M + k])
        return op.tocsr()

    def ev(self, op):
        return complex(np.dot(self.rho, op.diagonal()))

    def vacuum(self):
        bits = [1 if self.xs[k % self.M] < 0 else 0 for k in range(self.F*self.M)]
        n = len(bits)
        v = np.zeros(self.dim)
        v[sum(b << (n - 1 - k) for k, b in enumerate(bits))] = 1.0
        return v


# ------------------------------------------------------------------ Schleifenformel (exakte Wick-Kontraktion)
def loop_sum(ms, sigma, occf, xs):
    """Eine gerichtete Fermionschleife sigma (sigma[0] = 0) durch die Ströme mit Moden ms (Operatorordnung 0..n-1).
    Kontraktion c (Strom sigma_k) mit c^+ (Strom sigma_k+1): 1 - c(y), wenn sigma_k links steht, sonst c(y).
    Vorzeichen = (-1)^(Kreuzungen) der Paarung im Originalprodukt c^+_0 c_0 c^+_1 c_1 ..."""
    n = len(ms)
    assert abs(sum(ms)) < 1e-12, "geschlossene Schleife verlangt sum m = 0 (L_0-invarianter Zustand)"
    pairs = []
    for k in range(n):
        a, b = 2*sigma[k] + 1, 2*sigma[(k + 1) % n]
        pairs.append((min(a, b), max(a, b)))
    cross = sum(1 for (a, b), (c, d) in itertools.combinations(pairs, 2) if a < c < b < d or c < a < d < b)
    keys = {key(x) for x in xs}
    tot = 0.0
    for x in xs:
        y, w, ok = x, 1.0, True
        for k in range(n):
            y = y + ms[sigma[k]]
            if key(y) not in keys:
                ok = False
                break
            w *= (1 - occf(y)) if sigma[k] < sigma[(k + 1) % n] else occf(y)
        if ok:
            tot += w
    return (-1)**cross*tot


def corr_conn(Ts, ms, occf, xs):
    """Verbundene n-Punkt-Antwort <J(T_1)_{m_1} ... J(T_n)_{m_n}>_c als Summe der (n-1)! gerichteten Schleifen."""
    n = len(ms)
    val = 0j
    for perm in itertools.permutations(range(1, n)):
        sigma = (0,) + perm
        tr = np.trace(functools.reduce(np.matmul, [Ts[s] for s in sigma]))
        if abs(tr) < 1e-15:
            continue
        val += tr*loop_sum(ms, sigma, occf, xs)
    return val


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Nichtabelsche Drei- und Vierpunktantworten der Nahtströme (Wick)")
say("   Original auf der Naht: 16 Majorana = 8 komplexe Fermionen; die Quarter alpha = 3/4 sind ladungskonjugierte")
say("   alpha = 1/4-Fermionen, also alle acht im Quarterrahmen x in Z + 1/4 mit c(x) = 1/(1 + e^{4 pi x}). Der Zustand ist")
say("   U(8)-invariant (die Vierteldrehung bricht so(16) -> u(8)); gerechnet für su(N) c u(8) c so(16)_1 c E8_1.")

# 1.1 Schleifenformel gegen den exakten Fockraum (gleiche endliche Modenmenge, generisches beta, zufällige Matrizen)
rng = np.random.default_rng(75)
xs_small = [0.25 + r for r in range(-3, 3)]
occ_t = occ_factory(1.3)
sysF = Fermi(2, xs_small, occ_t)
Tr4 = [rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)) for _ in range(4)]
Jop = {}


def Jb(i, m):
    if (i, m) not in Jop:
        Jop[(i, m)] = sysF.J(Tr4[i], m)
    return Jop[(i, m)]


err = 0.0
for ms in [(1, -1), (2, -2)]:
    bf = sysF.ev(Jb(0, ms[0]) @ Jb(1, ms[1]))
    lf = corr_conn([Tr4[0], Tr4[1]], ms, occ_t, xs_small)
    err = max(err, abs(bf - lf)/max(1.0, abs(bf)))
for ms in [(1, 1, -2), (2, -1, -1), (-1, 2, -1), (1, -2, 1)]:
    bf = sysF.ev(Jb(0, ms[0]) @ Jb(1, ms[1]) @ Jb(2, ms[2]))
    lf = corr_conn(Tr4[:3], ms, occ_t, xs_small)
    err = max(err, abs(bf - lf)/max(1.0, abs(bf)))
for ms in [(1, -1, 1, -1), (1, 1, -1, -1), (2, -1, 1, -2), (1, -2, 2, -1)]:
    ops = [Jb(i, m) for i, m in enumerate(ms)]
    full = sysF.ev(ops[0] @ ops[1] @ ops[2] @ ops[3])
    disc = (sysF.ev(ops[0] @ ops[1])*sysF.ev(ops[2] @ ops[3]) + sysF.ev(ops[0] @ ops[2])*sysF.ev(ops[1] @ ops[3])
            + sysF.ev(ops[0] @ ops[3])*sysF.ev(ops[1] @ ops[2]))
    lf = corr_conn(Tr4, ms, occ_t, xs_small)
    err = max(err, abs(full - disc - lf)/max(1.0, abs(full)))
ZAHLEN["schleifenformel_gegen_fockraum_maxfehler"] = err
check("Schleifenformel (Kreuzungsvorzeichen, c bzw. 1 - c je Kontraktion, Farbspur) = exakter Fockraum (2 Flavours x 6 "
      "Quartermoden, beta = 1,3, zufällige komplexe Matrizen) für 2-, 3- und verbundene 4-Punkt-Antworten (Fehler < 1e-10)",
      err < 1e-10)

# 1.2 Dreipunktantwort: f- und d-Anteil, Ward-Identität, abelsche Kontrolle
R = 40
ZUST = {
    "Original (alpha = 1/4, beta = 4 pi)": ([0.25 + r for r in range(-R, R)], occ_factory(PI4)),
    "NS (alpha = 1/2, beta = 4 pi, so(16)-invariant)": ([0.5 + r for r in range(-R, R)], occ_factory(PI4)),
    "Vakuum (alpha = 1/4)": ([0.25 + r for r in range(-R, R)], occ_factory(np.inf)),
    "Vakuum (NS)": ([0.5 + r for r in range(-R, R)], occ_factory(np.inf)),
    "heiß (alpha = 1/4, beta = 1)": ([0.25 + r for r in range(-R, R)], occ_factory(1.0)),
}
gm = [np.array(M, dtype=complex) for M in (
    [[0, 1, 0], [1, 0, 0], [0, 0, 0]], [[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], [[1, 0, 0], [0, -1, 0], [0, 0, 0]],
    [[0, 0, 1], [0, 0, 0], [1, 0, 0]], [[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
    [[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], np.diag([1, 1, -2])/np.sqrt(3))]
T3 = [g/2 for g in gm]
fT = np.array([[[np.real(-2j*np.trace((a @ b - b @ a) @ c)) for c in T3] for b in T3] for a in T3])
dT = np.array([[[np.real(2*np.trace((a @ b + b @ a) @ c)) for c in T3] for b in T3] for a in T3])
tr3 = np.array([[[np.trace(a @ b @ c) for c in T3] for b in T3] for a in T3])
tr3r = np.transpose(tr3, (0, 2, 1))                                   # Tr(T^a T^c T^b)
MS3 = [(1, 1, -2), (2, -1, -1), (1, -3, 2), (3, -1, -2)]
ok_fd, ok_ab3, f_coef = True, True, {}
for name, (xs, oc) in ZUST.items():
    rows = []
    for ms in MS3:
        L1, L2 = loop_sum(ms, (0, 1, 2), oc, xs), loop_sum(ms, (0, 2, 1), oc, xs)
        C = tr3*L1 + tr3r*L2
        A = np.sum(C*fT)/np.sum(fT*fT)
        B = np.sum(C*dT)/np.sum(dT*dT)
        res = np.abs(C - A*fT - B*dT).max()
        ok_fd &= res < 1e-12 and abs(B) < 1e-12 and (abs(A) > 0.1 or ms == (1, -3, 2))
        ok_ab3 &= abs(L1 + L2) < 1e-12                                # U(1): <J J J> = L1 + L2
        rows.append((ms, complex(A)))
    f_coef[name] = {str(ms): f"{A.imag:.12f} i" for ms, A in rows}
ZAHLEN["dreipunkt_f_koeffizient_su3"] = f_coef
say(f"   f-Koeffizient von <J^a_m1 J^b_m2 J^c_m3> (su(3), Tr T^a T^b = delta/2): {f_coef}")
check("in allen Zuständen (Original, NS, Vakuum, heiß) ist <J^a J^b J^c> = A(m) f^{abc} exakt (Residuum 0; A != 0 außer "
      "für (1, -3, 2), wo J_2 im Vakuum vernichtet); der "
      "d^{abc}-Anteil (S1 - S2)/4 verschwindet, obwohl su(3) ein d-Symbol hat", ok_fd)
check("abelsche Kontrolle: die U(1)-Dreipunktantwort L1 + L2 verschwindet in allen Zuständen exakt (gaußsche Oszillatoren, "
      "z58), also trägt allein die nichtabelsche Struktur die Dreipunktantwort", ok_ab3)
Xr = [sum(rng.normal()*t for t in T3) for _ in range(3)]
ok_ward = True
for name, (xs, oc) in ZUST.items():
    for (m1, m2, m3) in MS3:
        lhs = corr_conn(Xr, (m1, m2, m3), oc, xs) - corr_conn([Xr[1], Xr[0], Xr[2]], (m2, m1, m3), oc, xs)
        rhs = corr_conn([Xr[0] @ Xr[1] - Xr[1] @ Xr[0], Xr[2]], (m1 + m2, m3), oc, xs)
        ok_ward &= abs(lhs - rhs) < 1e-12*max(1, abs(rhs)) and (abs(rhs) > 1e-3 or (m1, m2, m3) == (1, -3, 2))
check("Ward-Identität <[J^X_m1, J^Y_m2] J^Z_m3> = <J^[X,Y]_{m1+m2} J^Z_m3> für zufällige su(3)-Elemente in allen Zuständen: "
      "der Kommutator der Ströme ist ein Operator, keine Zahl (freie Felder hätten eine c-Zahl)", ok_ward)
A16 = [rng.normal(size=(16, 16)) for _ in range(3)]
A16 = [a - a.T for a in A16]
check("so(16) hat im Vektor (16 Majorana) kein d-Symbol: Tr(A{B, C}) = 0 für antisymmetrische A, B, C; die NS-Antwort ist "
      "so(16)-invariant, also ganz f^{abc} F (für E8 ebenso: Casimir-Grade 2, 8, 12, ..., kein Grad 3)",
      abs(np.trace(A16[0] @ (A16[1] @ A16[2] + A16[2] @ A16[1]))) < 1e-9)

# 1.3 Verbundene Vierpunktantwort
S2m = [np.array(M, dtype=complex)/2 for M in ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]])]
ok_vac4, ok_ab4, vier = True, True, {}
for name, (xs, oc) in ZUST.items():
    for ms in [(1, 1, -1, -1), (1, -1, 1, -1), (2, -1, -1, 0 + 0)][:2]:
        ab = sum(loop_sum(ms, (0,) + p, oc, xs) for p in itertools.permutations(range(1, 4)))
        ok_ab4 &= abs(ab) < 1e-12
    if name.startswith("Vakuum"):
        for _ in range(3):
            Ys = [sum(rng.normal()*t for t in T3) for _ in range(4)]
            got = corr_conn(Ys, (1, 1, -1, -1), oc, xs)
            Yb = Ys[1] @ Ys[2] - Ys[2] @ Ys[1]
            ref = np.trace(Ys[0] @ (Yb @ Ys[3] - Ys[3] @ Yb))
            ok_vac4 &= abs(got - ref) < 1e-12 and abs(ref) > 1e-3
    vier[name] = {f"{a}{b}{c}{d}": f"{corr_conn([S2m[a], S2m[b], S2m[c], S2m[d]], (1, 1, -1, -1), oc, xs).real:.10f}"
                  for (a, b, c, d) in [(0, 0, 0, 0), (0, 1, 0, 1), (0, 0, 1, 1), (0, 1, 1, 0)]}
ZAHLEN["vierpunkt_verbunden_su2_m(1,1,-1,-1)"] = vier
say(f"   verbundene Vierpunktantwort su(2) (Indizes 0,1,2 = x,y,z): {vier}")
check("abelsch verschwindet die verbundene Vierpunktantwort exakt (Summe der sechs Schleifen = 0, alle Zustände, zwei "
      "Modenmuster); im Vakuum ist sie nichtabelsch genau Tr(X[[Y, Z], W]) (Ward, zufällige su(3)-Elemente), also != 0",
      ok_ab4 and ok_vac4)
v4o = vier["Original (alpha = 1/4, beta = 4 pi)"]
check("im Originalzustand ist die verbundene su(2)-Vierpunktantwort für (x,y,x,y) und (x,x,y,y) von null verschieden, für "
      "einen einzelnen Generator (x,x,x,x) null", abs(float(v4o["0101"])) > 1e-3 and abs(float(v4o["0011"])) > 1e-3
      and abs(float(v4o["0000"])) < 1e-12)

# 1.4 Vielfachheiten: Nullvektor und Stufenzählung
xs_ns = [0.5 + r for r in range(-3, 3)]
sysV = Fermi(2, xs_ns, occ_factory(np.inf))
vac = sysV.vacuum()
Eop = sysV.J(np.array([[0, 1], [0, 0]], dtype=complex), -1)
Hop = sysV.J(np.diag([1, -1]).astype(complex)/np.sqrt(2), -1)
nE1, nE2 = np.linalg.norm(Eop @ vac)**2, np.linalg.norm(Eop @ (Eop @ vac))**2
nH2 = np.linalg.norm(Hop @ (Hop @ vac))**2
ZAHLEN["nullvektor"] = {"|E_-1|0>|^2": nE1, "|E_-1^2|0>|^2": nE2, "|H_-1^2|0>|^2": nH2}
check("Niveau 1: der Wurzelstrom hat |E_-1|0>|^2 = 1, aber (E_-1)^2|0> = 0 (Nullvektor, lokales Feld :E E: = 0); der "
      "Cartanstrom hat |(H_-1)^2|0>|^2 = 2 wie ein freies Boson. Ein freies Feld hätte in keiner Richtung einen Nullvektor",
      abs(nE1 - 1) < 1e-12 and nE2 < 1e-24 and abs(nH2 - 2) < 1e-12)
t = sp.Symbol("t")
ferm = sp.expand(sp.prod([(1 + t**(2*k + 1))**16 for k in range(4)]))
ferm_m = sp.expand(sp.prod([(1 - t**(2*k + 1))**16 for k in range(4)]))
so16_g2 = int(sp.Poly((ferm + ferm_m)/2, t).coeff_monomial(t**4))
so16_g1 = int(sp.Poly((ferm + ferm_m)/2, t).coeff_monomial(t**2))
qq = sp.Symbol("q")
E4s = 1 + 240*sum(sp.divisor_sigma(n, 3)*qq**n for n in range(1, 4))
eta8 = sp.series(sp.prod([(1 - qq**n)**(-8) for n in range(1, 4)]), qq, 0, 4).removeO()
e8 = sp.Poly(sp.expand(E4s*eta8), qq)
e8_g1, e8_g2 = int(e8.coeff_monomial(qq)), int(e8.coeff_monomial(qq**2))


def frei(d):
    return d + d*(d + 1)//2


ZAHLEN["stufe_2"] = {"so16_1": [so16_g1, so16_g2, frei(120)], "E8_1": [e8_g1, e8_g2, frei(248)]}
check(f"Stufenzählung des Vakuummoduls: so(16)_1 hat {so16_g1} Ströme und {so16_g2} Zustände auf Stufe 2 (16 NS-Majorana, "
      f"gerade), E8_1 {e8_g1} und {e8_g2} (E4/eta^8); freie Felder über denselben Strömen gäben {frei(120)} bzw. {frei(248)}",
      so16_g1 == 120 and so16_g2 == 2076 and e8_g1 == 248 and e8_g2 == 4124)
befund("E: Die nichtabelschen Antworten der Naht sind nicht gaußsch: <J J J> = A f^{abc} != 0, verbundene Vierpunktantwort != 0, "
       "Kommutator operatorwertig, Nullvektoren auf Niveau 1. Die abelschen (Cartan-)Ströme allein sind gaußsch (z58)", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Können die Ströme freie 4D-Felder sein? (externe Behauptung)")
y = sp.symbols("y0:4", real=True)
x2, x3 = (0, 0, 0, 0), (3, 1, 0, 2)
r12 = sp.sqrt(sum((y[i] - x2[i])**2 for i in range(4)))
r13 = sp.sqrt(sum((y[i] - x3[i])**2 for i in range(4)))
r23 = sp.sqrt(sum((x2[i] - x3[i])**2 for i in range(4)))
lap = lambda f: sum(sp.diff(f, yi, 2) for yi in y)  # noqa: E731
pt = dict(zip(y, (1, 2, -1, 1)))
box2 = sp.simplify(lap(1/r12**2))
box3 = sp.N(lap(1/(r12*r13*r23)).subs(pt), 20)
ZAHLEN["box_dreipunkt_delta1"] = str(box3)
check("4D: Box(1/x12^2) = 0 (Delta = 1 ist frei), aber Box_1 der einzigen konformen Dreipunktstruktur 1/(|x12||x13||x23|) "
      f"ist {sp.N(box3, 8)} != 0: eine Dreipunktfunktion von Delta = 1-Skalaren muss verschwinden", box2 == 0 and abs(box3) > 1e-4)
z1, z2, z3 = sp.symbols("z1 z2 z3")
S2d = 1/((z1 - z2)*(z1 - z3)*(z2 - z3))
Dl = sp.Symbol("Delta", positive=True)
R12, R13, R23 = sp.symbols("R12 R13 R23", positive=True)
S4d = 1/(R12*R13*R23)**Dl
S4d_sw = S4d.subs({R13: R23, R23: R13}, simultaneous=True)                      # x1 <-> x2: R12 fest, R13 <-> R23
check("Bose-Symmetrie: in 2D ist 1/(z12 z13 z23) antisymmetrisch unter 1 <-> 2, verträgt also f^{abc}; für 4D-Skalare "
      "gleicher Dimension (jedes Delta) ist die einzige konforme Struktur 1/(r12 r13 r23)^Delta symmetrisch, also muss ihr "
      "f^{abc}-Anteil verschwinden",
      sp.simplify(S2d + S2d.subs({z1: z2, z2: z1}, simultaneous=True)) == 0 and sp.simplify(S4d - S4d_sw) == 0)
xr, qv = sp.symbols("r q", positive=True)
d2 = sp.simplify(xr*sp.diff(sp.log(xr**(-2*qv**2)), xr))
d4 = sp.simplify(xr*sp.diff(qv**2/(4*sp.pi**2*xr**2), xr))
check("Vertexoperatoren: in 2D ist d log<V V>/d log r = -2 q^2 konstant (Gewicht q^2/2 = 1 für |q|^2 = 2), in 4D ist "
      "d log exp(q^2/(4 pi^2 r^2))/d log r = -q^2/(2 pi^2 r^2) ortsabhängig: e^{i q phi} hat in 4D keine Skalendimension",
      sp.diff(d2, xr) == 0 and sp.diff(d4, xr) != 0)
roots = [v for v in itertools.product([0, 1, -1], repeat=8) if sum(abs(c) for c in v) == 2]
roots = [np.array(v, dtype=float) for v in roots]
roots += [np.array(s)/2 for s in itertools.product([1, -1], repeat=8) if s.count(-1) % 2 == 0]
Rm = np.array(roots)
idx_adj = Rm.T @ Rm
check("E8-Wurzeln: 112 + 128 = 240, alle mit Norm 2; sum_alpha alpha_i alpha_j = 60 delta_ij, der Vektor von so(16) gibt "
      "2 delta_ij: adjungierte Majorana-Fermionen tragen Niveau 60/2 = 30 = h^v(E8), nicht 1",
      len(roots) == 240 and np.allclose((Rm**2).sum(1), 2) and np.allclose(idx_adj, 60*np.eye(8)))


def linie(n, D):
    return math.comb(n + D - 1, 2*D - 1)


rat = [linie(n, 2)/linie(n, 1) for n in range(1, 6)]
check("Linienantwort einer 4D-Primären der Dimension Delta ist C(n + Delta - 1, 2 Delta - 1): Delta = 1 gibt n (= 2D-Strom, z71), "
      "ein Bilinear :Phi^+ T Phi: (Delta = 2) gibt n(n^2 - 1)/6; das Verhältnis ist nicht konstant",
      [linie(n, 1) for n in range(1, 6)] == [1, 2, 3, 4, 5] and len(set(rat)) == 5)
befund("E: Die externe Behauptung hält und wird schärfer. Jede skalare 4D-Lesart der nichtabelschen Ströme (die K-Typen "
       "(j, j) der Paar-Verdopplungs-Abbildung, die Linienfortsetzung z71) ist ausgeschlossen, frei oder nicht: Bose-Symmetrie "
       "tötet den f^{abc}-Anteil skalarer Dreipunktfunktionen für jedes Delta, und Delta = 1 erzwingt zusätzlich Box phi = 0. "
       "Eine 4D-Lesart muss Spin tragen: Vektorströme J^a_mu (Delta = 3), deren Dreipunktfunktion f-Strukturen hat", True)
befund("E/L: Was eine 4D-Lesart sein müsste: zusammengesetzte Vektorströme. so(16) c E8: J^ab_mu als Bilineare eines freien "
       "16-Pletts (etwa psi-bar gamma_mu psi oder Phi d_mu Phi, Delta = 3; die z71-Linienidentität mit Delta = 1 geht verloren). Die 128 Spinorströme von E8 haben keine freie "
       "4D-Lesart: Sie sind Vertexoperatoren (brauchen den 2D-Logarithmus), und als Bilineare bräuchte E8 Felder in der 248 "
       "mit Niveau 30 statt 1. E8_1 bleibt eine Algebra der 2D-Naht (Defekt) oder ein globaler Symmetrieinhalt mit "
       "zusammengesetzten Strömen einer wechselwirkenden Theorie", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Spin 2 -- Virasoro-Antwort der Naht gegen die masselosen Spin-2-K-Typen von SU(2,2)")
xs1 = [0.5 + r for r in range(-5, 5)]
sys1 = Fermi(1, xs1, occ_factory(np.inf))
vac1 = sys1.vacuum()
one = np.array([[1.0 + 0j]])


def Lop(n):
    op = sps.csr_matrix((sys1.dim, sys1.dim), dtype=complex)
    for j, x in enumerate(xs1):
        k = sys1.idx.get(key(x + n))
        if k is not None:
            op = op + (x + n/2)*(sys1.c[j].T @ sys1.c[k])
    return op.tocsr()


Lm = {n: Lop(-n) for n in range(1, 5)}
Jm = {n: sys1.J(one, -n) for n in range(1, 5)}
v = {n: Lm[n] @ vac1 for n in range(1, 5)}
w = {n: Jm[n] @ vac1 for n in range(1, 5)}
ok_gram = abs(np.linalg.norm(v[1])) < 1e-14 and all(abs(np.vdot(v[n], v[n]) - (n**3 - n)/12) < 1e-12 and
                                                     abs(np.vdot(w[n], w[n]) - n) < 1e-12 for n in range(2, 5))
check("Vakuummodul von T (ein komplexes Fermion, c = 1): L_-1|0> = 0, |L_-n|0>|^2 = c(n^3 - n)/12, |J_-n|0>|^2 = n; der "
      "Einteilchenraum von T hat genau einen Zustand je n >= 2: chi_T(q) = q^2/(1 - q)", ok_gram)


def pair_vec(xp, yh):
    jp, jh = sys1.idx[key(xp)], sys1.idx[key(yh)]
    return sys1.c[jp].T @ (sys1.c[jh] @ vac1)


ok_pair, ok_qp = True, True
Pbasis = {}
for n in range(1, 5):
    P = np.array([pair_vec(xp, xp - n) for xp in [0.5 + i for i in range(n)]]).T
    Pbasis[n] = P
    coef_v, *_ = np.linalg.lstsq(P, v[n], rcond=None)
    coef_w, *_ = np.linalg.lstsq(P, w[n], rcond=None)
    ok_pair &= (np.linalg.matrix_rank(P) == n and np.linalg.norm(P @ coef_v - v[n]) < 1e-12
                and np.linalg.norm(P @ coef_w - w[n]) < 1e-12 and abs(np.vdot(w[n], v[n])) < 1e-12
                and (n == 1 or np.linalg.matrix_rank(np.array([v[n], w[n]]).T) == 2))
L1op = Lop(1)
for k in range(1, 5):
    img = L1op @ Pbasis[k]
    kern = Pbasis[k].shape[1] - np.linalg.matrix_rank(img, tol=1e-10)
    ok_qp &= kern == 1
check("Paarraum P_n (Teilchen + Loch, Energie n) hat dim n und enthält J_-n|0> und L_-n|0> (orthogonal, Rang 2 für n >= 2): "
      "Strom und Spannungstensor teilen sich dieselben n Paare", ok_pair)
check("in jedem P_k (k = 1..4) gibt es genau einen neuen Quasiprimären (Kern von L_1 eindimensional): J, T, W_3, W_4; die "
      "'m Paare der Strommode m' der externen Behauptung sind der ganze W_{1+inf}-Turm, nicht der Strom allein", ok_qp)


def kernel_dim(s2, k):
    """Polynomlösungen vom Grad k der masselosen Spin-(s2/2)-Gleichung in 4D (komplex); s2 = 0: Box phi = 0."""
    mons = list(itertools.combinations_with_replacement(range(4), k))
    mid = {m: i for i, m in enumerate(mons)}

    def deriv(mon, a):                                   # d/dx^a eines Monoms (Tupel von Variablenindizes)
        if a not in mon:
            return None, 0
        lst = list(mon)
        c = lst.count(a)
        lst.remove(a)
        return tuple(lst), c
    if s2 == 0:
        low = list(itertools.combinations_with_replacement(range(4), k - 2)) if k >= 2 else []
        lid = {m: i for i, m in enumerate(low)}
        rows = [[0]*len(mons) for _ in low]
        for j, mon in enumerate(mons):
            for (a1, a2, sg) in ((0, 3, 1), (1, 2, -1)):              # Box = 2(d00 d11 - d01 d10), a = 2A + A'
                m1, c1 = deriv(mon, a1)
                if m1 is None:
                    continue
                m2, c2 = deriv(m1, a2)
                if m2 is None:
                    continue
                rows[lid[m2]][j] += sg*c1*c2
        nunk = len(mons)
    else:
        low = list(itertools.combinations_with_replacement(range(4), k - 1)) if k >= 1 else []
        lid = {m: i for i, m in enumerate(low)}
        nunk = (s2 + 1)*len(mons)
        rows = []
        for Ep in (0, 1):
            for tt in range(s2):
                block = [[0]*nunk for _ in low]
                for j, mon in enumerate(mons):
                    m1, c1 = deriv(mon, 2 + Ep)                    # d_{1E'} phi_t
                    if m1 is not None:
                        block[lid[m1]][tt*len(mons) + j] += c1
                    m0, c0 = deriv(mon, 0 + Ep)                    # - d_{0E'} phi_{t+1}
                    if m0 is not None:
                        block[lid[m0]][(tt + 1)*len(mons) + j] -= c0
                rows += block
    if not rows:
        return nunk
    rk = DomainMatrix([[sp.QQ(x) for x in r] for r in rows], (len(rows), nunk), sp.QQ).rank()
    return nunk - rk


kern = {s2: [kernel_dim(s2, k) for k in range(5)] for s2 in (0, 2, 4)}
ZAHLEN["4D_polynomloesungen_grad_0..4"] = {"Skalar": kern[0], "Maxwell (1,0)": kern[2], "Weyl (2,0)": kern[4]}
check("4D-Kern der Feldgleichungen (Polynomlösungen vom Grad k = 0..4, exakt über Q): Skalar (k+1)^2, Maxwell (k+3)(k+1), "
      f"linearisierte Weyl-/Bianchi-Gleichung (k+5)(k+1); gefunden {kern}",
      kern[0] == [(k + 1)**2 for k in range(5)] and kern[2] == [(k + 3)*(k + 1) for k in range(5)]
      and kern[4] == [(k + 5)*(k + 1) for k in range(5)])
# Einteilchencharaktere auf R x S^3: Energie n = Delta + k, je Helizität (n + s)(n - s)
d4 = {n: (2*(n*n - 4) if n >= 3 else 0) for n in range(1, 13)}
d4_kern = {3 + k: 2*kern[4][k] for k in range(5)}
dT = {n: (1 if n >= 2 else 0) for n in range(1, 13)}
dP = {n: (n*n if n >= 2 else 0) for n in range(1, 13)}
dMax = {n: (2*(n*n - 1) if n >= 2 else 0) for n in range(1, 13)}
qs = sp.Symbol("q")
ser = sp.series(2*(5*qs**3 - 3*qs**4)/(1 - qs)**3, qs, 0, 13).removeO()
ok_char = all(ser.coeff(qs, n) == d4[n] for n in range(1, 13)) and all(d4[n] == d4_kern[n] for n in d4_kern)
ZAHLEN["charaktere_d_n_n=1..12"] = {"T_einteilchen": list(dT.values()), "T_paar_verdopplung": list(dP.values()),
                                   "4D_spin2_beide_helizitaeten": list(d4.values()), "4D_Maxwell": list(dMax.values())}
say(f"   d_n (n = 1..12): T einteilchen {list(dT.values())}")
say(f"                    T über Paare x Tomita {list(dP.values())}")
say(f"                    4D Spin 2 (+-2)       {list(d4.values())}")
check("4D-Spin-2-Einteilchencharakter 2(5q^3 - 3q^4)/(1 - q)^3 = sum 2(n^2 - 4) q^n stimmt mit den Kernzahlen überein; "
      "T ergibt q^2/(1 - q) (wörtlich) bzw. sum n^2 q^n (Paar-Verdopplungs-Abbildung): keine Stufe n stimmt überein "
      "(n^2 = 2(n^2 - 4) hieße n^2 = 8), niedrigste Energie 2 gegen 3",
      ok_char and all(dT[n] != d4[n] and dP[n] != d4[n] for n in range(2, 13)) and dT[2] == 1 and d4[2] == 0)
befund("E: SU(2)xSU(2)-Inhalt. Die Paar-Verdopplungs-Abbildung gibt bei Energie n = 2j + 1 nur (j, j), weil die "
       "Tomita-Verdopplung einer Mode ihre Energie erhält (beide Hälften vom Grad 2j). Spin 2 braucht (j+1, j-1) + (j-1, j+1), "
       "also Hälften der Energien n + 2 und n - 2; deren Kreuzantwort verschwindet im L_0-invarianten Originalzustand", True)
gramT = {n: (n**3 - n)/12 for n in range(2, 9)}
r2 = [gramT[n]/linie(n, 2) for n in range(2, 9)]
r3 = [gramT[n]/linie(n, 3) for n in range(3, 9)]
check("Linientest: die T-Antwort c(n^3 - n)/12 ist genau die Linienantwort einer Delta = 2-Primären (Verhältnis konstant "
      "c/2), nicht die des Weyl-Tensors Delta = 3 (C(n+2, 5), Verhältnis nicht konstant)",
      max(r2) - min(r2) < 1e-12 and abs(r2[0] - 0.5) < 1e-12 and len(set(np.round(r3, 12))) == len(r3))


# Sym^2 des Skalarsingletons (Flato-Fronsdal in 4D), K-Typen (2E, 2j1, 2j2)
def prod2(a, b):
    return list(range(abs(a - b), a + b + 1, 2))


def verma(E2, a, b, Emax2):
    ch = Counter()
    k = 0
    while E2 + 2*k <= Emax2:
        for l_ in range(k//2 + 1):
            c = k - 2*l_
            for A in prod2(a, c):
                for B in prod2(b, c):
                    ch[(E2 + 2*k, A, B)] += 1
        k += 1
    return ch


Emax2 = 2*12
single = [(2*(1 + k), k, k) for k in range(12)]
sym2 = Counter()
for i, u in enumerate(single):
    for v_ in single[i:]:
        E2 = u[0] + v_[0]
        if E2 > Emax2:
            continue
        if u == v_:
            for A in range(2*u[1], -1, -4):
                for B in range(2*u[2], -1, -4):
                    sym2[(E2, A, B)] += 1
            for A in range(2*u[1] - 2, -1, -4):
                for B in range(2*u[2] - 2, -1, -4):
                    sym2[(E2, A, B)] += 1
        else:
            for A in prod2(u[1], v_[1]):
                for B in prod2(u[2], v_[2]):
                    sym2[(E2, A, B)] += 1
rest, found, ok_peel = Counter(sym2), [], True
while True:
    keys_ = [k_ for k_, m_ in rest.items() if m_ != 0]
    if not keys_:
        break
    e2 = min(k_[0] for k_ in keys_)
    for k_ in sorted(kk for kk in keys_ if kk[0] == e2):
        mult = rest[k_]
        if mult < 0:
            ok_peel = False
        E2, a, b = k_
        if a == b and a >= 1 and E2 == 2*(a + 2):
            typ, ch = f"Strom s = {a}", verma(E2, a, b, Emax2)
            ch.subtract(verma(E2 + 2, a - 1, b - 1, Emax2))
        else:
            typ, ch = ("Skalar Delta = 2" if (a, b) == (0, 0) and E2 == 4 else f"lang {k_}"), verma(E2, a, b, Emax2)
        for kk, mm in ch.items():
            rest[kk] -= mult*mm
        found.append((k_, typ, mult))
    if len(found) > 40:
        ok_peel = False
        break
ZAHLEN["sym2_skalar_zerlegung"] = [f"{t_} bei E = {k_[0]//2}, Vielfachheit {m_}" for k_, t_, m_ in found]
say(f"   Sym^2(Skalar) bis E = 12: {ZAHLEN['sym2_skalar_zerlegung']}")
expect = [((4, 0, 0), "Skalar Delta = 2", 1)] + [((2*(s + 2), s, s), f"Strom s = {s}", 1) for s in (2, 4, 6, 8, 10)]
check("Sym^2 des 4D-Skalarsingletons (Sugawara-Bild T = :J J:/2 unter J -> phi) zerfällt exakt in phi^2 (Delta = 2, passend "
      "zum Linientest) und die erhaltenen Ströme gerader Spins, darunter T_mu nu (Delta = 4); keine masselose Helizität-2-"
      "Darstellung (2,0)/(0,2) bei Delta = 3", ok_peel and found == expect)
befund("E: Der Spannungstensor der Naht bildet sich nicht auf das 4D-Graviton ab. Wörtlich hat sein Einteilchenraum den "
       "Charakter q^2/(1 - q); die Paar-Verdopplungs-Abbildung, die den Strom auf (j, j) mit m^2 schickt, schickt T auf "
       "dieselben K-Typen (j, j), weil J und T denselben Paarraum teilen; die Linienantwort ist die von Delta = 2. Gegen "
       "2(n^2 - 4) bei (j+1, j-1) + (j-1, j+1) ab n = 3 stimmt keine Stufe", True)
befund("E/L: Im Sugawara-Bild landet T im Zweiteilchensektor des Skalars: phi^2 (Delta = 2) und der 4D-Spannungstensor "
       "T_mu nu als zusammengesetzter Strom (Delta = 4). Das ist 4D-Materie mit Energietensor, kein Graviton", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Araki-Woods-Dubleton der externen Analyse (zwei Quartermoden + Tomita-Partner)")
say("   a_i = J^(i)_m/sqrt(m) (Quarter i = 1, 2), b_i = J a_i J; d_i = (a_i - q b_i^+)/sqrt(1-q^2), e_i = (b_i - q a_i^+)/sqrt(1-q^2)")
nmax = 40
Aop = np.diag(np.sqrt(np.arange(1, nmax)), 1)
Id = np.eye(nmax)
a_, b_ = np.kron(Aop, Id), np.kron(Id, Aop)
Na, Nb = a_.T @ a_, b_.T @ b_
low = [i*nmax + j for i in range(16) for j in range(16)]
ok_aw, ok_kms, ok_D, ok_E, ok_comp = True, True, True, True, True
bcoef = {}
for qv_ in (0.3, float(np.exp(-2*np.pi))):
    s_ = np.sqrt(1 - qv_**2)
    d_ = (a_ - qv_*b_.T)/s_
    e_ = (b_ - qv_*a_.T)/s_
    Om = np.zeros(nmax*nmax)
    for n in range(nmax):
        Om[n*nmax + n] = s_*qv_**n
    ok_aw &= np.linalg.norm(d_ @ Om) < 1e-12 and np.linalg.norm(e_ @ Om) < 1e-12
    for X, Y, val in ((d_, d_.T, 1), (e_, e_.T, 1), (d_, e_, 0), (d_, e_.T, 0)):
        Cm = X @ Y - Y @ X
        ok_aw &= np.abs(Cm[np.ix_(low, low)] - val*np.eye(len(low))).max() < 1e-10
    ratio = (Om @ a_ @ a_.T @ Om)/(Om @ a_.T @ a_ @ Om)
    Dh = np.diag(qv_**(np.diag(Na) - np.diag(Nb)))                  # Delta^{1/2} = q^{N_a - N_b}
    ok_kms &= abs(ratio - qv_**-2) < 1e-9*ratio and abs(np.linalg.norm(Dh @ a_ @ Om) - np.linalg.norm(a_.T @ Om)) < 1e-10
    Nd, Ne = d_.T @ d_, e_.T @ e_
    ok_D &= np.abs((Na - Nb - Nd + Ne)[np.ix_(low, low)]).max() < 1e-9
    Cn = (Nd + Ne) @ a_ - a_ @ (Nd + Ne)
    target = -((1 + qv_**2)*a_ - 2*qv_*b_.T)/(1 - qv_**2)
    ok_E &= np.abs((Cn - target)[np.ix_(low, low)]).max() < 1e-9
    bcoef[f"q = {qv_:.6g}"] = 2*qv_/(1 - qv_**2)
    psi = np.linalg.matrix_power(d_.T, 4) @ Om
    w_low = sum(abs(psi[i*nmax + j])**2 for i in range(4) for j in range(nmax))
    ok_comp &= w_low < 1e-20*np.vdot(psi, psi).real and np.vdot(psi, psi).real > 1
ZAHLEN["AW_b_dagger_koeffizient_in_[N_d+N_e, a]"] = bcoef
check("Araki-Woods (q = 0,3 und q = e^{-2 pi}, 40 Niveaus): d Omega = e Omega = 0 für den Thermofeldvektor, "
      "[d, d^+] = [e, e^+] = 1, [d, e] = [d, e^+] = 0: d_i, e_i sind freie Oszillatoren mit Omega als Fockvakuum", ok_aw)
check("KMS-Anschluss an das Original: <a a^+>/<a^+ a> = q^-2 = e^{4 pi m} für q = e^{-2 pi m} (z58), und "
      "|Delta^{1/2} a Omega| = |a^+ Omega| mit Delta^{1/2} = q^{N_a - N_b}: der Modulgenerator ist K = 4 pi m (N_a - N_b)", ok_kms)
check("N_a - N_b = N_d - N_e als Operatoridentität: die Helizitätszahl D der externen Analyse ist (bis 4 pi m) genau der "
      "Modulgenerator des Originals, nicht seine Energie", ok_D)
check("die konforme Energie E = 1 + (N_d + N_e)/2 erhält die Algebra des Originals nicht: [N_d + N_e, a] = -((1+q^2) a - "
      "2q b^+)/(1 - q^2) enthält den Kommutanten-Operator b^+ (Koeffizient 2q/(1-q^2) != 0)", ok_E)
check("(d^+)^4 Omega liegt ganz bei N_a >= 4: die Zustände mit D = 4 sind Komposita aus mindestens vier Stromquanten "
      "(Gewicht bei N_a <= 3 exakt null)", ok_comp)
cnt = Counter()
for occs in itertools.product(range(13), repeat=4):
    nd, ne = occs[0] + occs[1], occs[2] + occs[3]
    if nd + ne <= 12:
        cnt[(nd - ne, nd + ne)] += 1
okD = True
for D, s2_, kl in ((0, 0, kern[0]), (2, 2, kern[2]), (-2, 2, kern[2]), (4, 4, kern[4]), (-4, 4, kern[4])):
    for k in range(5):
        okD &= cnt[(D, abs(D) + 2*k)] == kl[k]
emin = {D: min(tot for (DD, tot) in cnt if DD == D)/2 + 1 for D in range(-8, 9)}
check("Fockraum der vier Oszillatoren (d_1, d_2, e_1, e_2): der Sektor D hat bei E = 1 + |D|/2 + k genau die Vielfachheit der "
      "4D-Kernzahlen (D = 0: (k+1)^2, D = +-2: (k+3)(k+1), D = +-4: (k+5)(k+1)) und die niedrigste Energie 1 + |D|/2: "
      "Helizität +-2 bei E = 3, exakt der 4D-Spin-2-Charakter", okD and all(emin[D] == 1 + abs(D)/2 for D in emin))
ok_shift = np.linalg.norm(Lop(-1) @ w[1] - w[2]) < 1e-12
check("[L_-1, J_-1] = J_-2 im Fockraum: die Virasoro-Moden verschieben die Strommode m; T wirkt nicht innerhalb eines "
      "Dubletons fester Mode und ist nicht der Erzeuger seines D = +-4-Sektors", ok_shift)
befund("E (Mathematik): Die externe Behauptung zum Dubleton hält: Mit Omega als Fockvakuum der d, e ist der Oszillatorraum "
       "die Summe aller masselosen Einteilchenräume, je Helizität h = D/2 einmal, Helizität +-2 ab E = 3 mit dem exakten "
       "4D-Charakter. Das ist die Oszillator- bzw. Twistorquantisierung von SU(2,2) (Günaydin-Saclioglu, Penrose)", True)
befund("E/L: Was davon Dynamik des Originals ist: D ist der Modulgenerator (Zeit des Zustands), nicht die Energie; E erhält "
       "die Algebra des Originals nicht (mischt mit dem Kommutanten); die Helizität-2-Zustände sind Vierquantenkomposita; "
       "jede Mode m liefert eine volle Kopie aller Helizitäten (q = e^{-2 pi m} ändert den Inhalt nicht), dazu alle "
       "halbzahligen Helizitäten aus bosonischen Oszillatoren; die SU(2)_L-Drehung der beiden Quarter ist nur auf den "
       "Oszillatoren eine Symmetrie (L). Kinematik im Sinn von DIMENSION.UPLIFT.FIREWALL.01", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Weinberg-Witten-Typisierung und Bilanz")
G = np.zeros((4, 4), dtype=complex)
G[1, 2], G[2, 1] = -1j, 1j                                              # Drehung um z auf (t, x, y, z)
evv = np.linalg.eigvalsh(G)
symb = [(i, j) for i in range(4) for j in range(i, 4)]
GS = np.zeros((10, 10), dtype=complex)
for c_, (i, j) in enumerate(symb):
    for k in range(4):
        for (p_, r_) in ((k, j), (i, k)):
            coeff = G[k, i] if (p_, r_) == (k, j) else G[k, j]
            key2 = tuple(sorted((p_, r_)))
            GS[symb.index(key2), c_] += coeff
evs = np.linalg.eigvals(GS)
check("Drehladungen um die Impulsachse: Vektor max |sigma| = 1, symmetrischer Tensor max |sigma| = 2. Weinberg-Witten: "
      "Vorwärts-Matrixelemente <p', h|O|p, h> eines kovarianten Operators verlangen 2|h| <= max |sigma|, also |h| <= 1/2 "
      "für J^mu und |h| <= 1 für T^mu nu; der Dubleton-Sektor D = +-4 hat 2|h| = 4",
      abs(np.max(np.abs(evv)) - 1) < 1e-12 and abs(np.max(np.abs(evs)) - 2) < 1e-12)
befund("Typisierung (GRAV.WEINBERG.WITTEN.01): Sugawara-/Paarweg: keine Teilchen mit |h| > 1, T_mu nu ist ein "
       "zusammengesetzter Strom (Delta = 4), WW schweigt, kein Graviton. Dubletonweg: Helizität +-2 existiert, kann aber an "
       "keinen kovarianten erhaltenen T^mu nu koppeln; umgangene Prämisse ist die Existenz eines lokalen kovarianten "
       "T^mu nu mit nichtverschwindenden Vorwärts-Matrixelementen. Der Dubleton hat nur globale SU(2,2)-Erzeuger, und "
       "seine Energie E ist nicht einmal eine Automorphie des Originals. Graviton als Eichfeld bräuchte eine "
       "Eichstruktur, die die Naht nicht liefert", True)
tafel = [
    ("nichtabelsche Dreipunktantwort", "E: A(m) f^{abc} in Original, NS, Vakuum; kein d-Anteil; Ward exakt"),
    ("verbundene Vierpunktantwort", "E: != 0 nichtabelsch, = 0 abelsch; Vakuum = Tr(X[[Y,Z],W])"),
    ("freie 4D-Felder", "E: ausgeschlossen, ebenso jede skalare Lesart (Bose-Symmetrie; Delta = 1 zusätzlich Box phi = 0)"),
    ("4D-Lesart der Ströme", "E/L: nur Vektorkomposita J^a_mu; so(16) als Bilineare, die 128 Spinorströme von E8 gar nicht"),
    ("T -> Spin 2 wie J -> Skalar", "E negativ: T teilt den Paarraum mit J, gibt (j, j), Linienantwort Delta = 2"),
    ("4D-Spin-2-Charakter", "E: 2(5q^3 - 3q^4)/(1-q)^3, aus dem Kern der linearisierten Weyl-Gleichung"),
    ("Araki-Woods-Dubleton", "E Mathematik: D = +-4 gibt Helizität +-2 bei E = 3; L als Physik (D = Modulgenerator)"),
    ("Weinberg-Witten", "keine kovariante T^mu nu-Kopplung für D = +-4; umgangene Prämisse benannt"),
    ("Wechselwirkung, nichtlineare Gravitation", "O: braucht lokale Eichstruktur, kubische Ecke (z65), Zeit (z71)"),
]
for k_, v_ in tafel:
    say(f"   {k_}: {v_}")
ZAHLEN["bilanz"] = dict(tafel)
urteil = ("Die nichtabelschen Antworten der Naht sind exakt gerechnet und nicht gaußsch: <J J J> = A(m) f^{abc} im "
          "Originalzustand, im NS-Zustand und im Vakuum, ohne d-Anteil, mit exakter Ward-Identität; die verbundene "
          "Vierpunktantwort ist nichtabelsch von null verschieden und abelsch null; auf Niveau 1 gibt es Nullvektoren "
          "(E_-1)^2|0> = 0. Damit hält die externe Behauptung in schärferer Form: Keine skalare 4D-Lesart (frei oder nicht, "
          "jedes Delta) trägt f^{abc}, eine 4D-Lesart der Ströme muss aus Vektorkomposita J^a_mu bestehen, und die 128 Spinorströme von E8 haben "
          "gar keine freie 4D-Lesart. Der Spannungstensor der Naht geht nicht wie der Strom auf das Graviton: Sein "
          "Einteilchencharakter ist q^2/(1 - q), er teilt den Paarraum mit dem Strom (also dieselben K-Typen (j, j)), seine "
          "Linienantwort ist die von Delta = 2, und im Sugawara-Bild landet er bei phi^2 und dem zusammengesetzten "
          "4D-Spannungstensor. Gegen den 4D-Spin-2-Charakter 2(5q^3 - 3q^4)/(1 - q)^3 stimmt keine Stufe. Der "
          "Araki-Woods-Dubleton der externen Analyse gibt Helizität +-2 bei E = 3 mit exaktem Charakter, aber als "
          "Vierquantenkomposita, mit dem Modulgenerator als Helizität und einer Energie, die nicht die Dynamik des "
          "Originals ist; nach Weinberg-Witten können diese Zustände an keinen kovarianten T^mu nu koppeln. Offen bleiben "
          "Wechselwirkung und nichtlineare Gravitation: eine lokale Eichstruktur für den D = +-4-Sektor, die kubische "
          "Gravitonecke (Daten um die Ambitwistorquadrik, z65) und die physikalische Zeit (z71).")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z75_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
