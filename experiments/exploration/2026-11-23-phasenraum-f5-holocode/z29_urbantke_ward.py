"""Z29: Gravitation aus der Familienstruktur? Urbantke-Metrik, Wardidentität, positiver Kegel.

Prüft die Gravitationsrechnungen der vierten externen Analyse (24.11.):
Teil 1  Invariante symmetrische Formen auf Lambda^2 F: Vierteldrehung gegen volle Familiengruppe S4.
Teil 2  Urbantke-Metrik und Plebanski-Simplizität für die äquivarianten Zweiformen
        B^i = a dt^dx^i + (b/2) eps^ijk dx^j^dx^k; Lorentzsche Realitätsbedingung.
Teil 3  Positiver Kegel von Herm(2): Jede quadratische Form, die auf allen Rang-1-Richtungen verschwindet,
        ist proportional zur Minkowskiform.
Teil 4  Linearisierte Diffeomorphismen: Wardidentität q^mu K_{mu nu, rho sigma} = 0 für den allgemeinen
        kovarianten quadratischen Operator bis zu zwei Ableitungen (euklidisch, d = 4).
Teil 5  Nichtminimale Kopplung xi: diffeomorphismusinvariant und doch G-verschiebend (a_1 = (1/6 - xi) R).
Teil 6  Longitudinaler Test auf dem Z27-Gitter: reine Koordinatenänderung z' = z + (eps/k) sin kz
        (flache Metrik g_zz = (1 + eps cos kz)^2). Kovariant wäre die Antwort null.
"""
import itertools
import json
from pathlib import Path

import numpy as np
import sympy as sp

rng = np.random.default_rng(29)
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ Teil 1
say("Teil 1: invariante Formen auf Lambda^2 F")
R90 = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.0]])
C3 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0.0]])
grp = {tuple(np.eye(3).ravel())}
frontier = [np.eye(3)]
while frontier:
    new = []
    for g in frontier:
        for s_ in (R90, C3):
            h = np.round(s_ @ g, 12)
            if tuple(h.ravel()) not in grp:
                grp.add(tuple(h.ravel())); new.append(h)
    frontier = new
G24 = [np.array(t).reshape(3, 3) for t in grp]
def inv_forms(gens):
    Bs = [np.outer(np.eye(3)[i], np.eye(3)[j]) + np.outer(np.eye(3)[j], np.eye(3)[i]) for i in range(3) for j in range(i, 3)]
    rows = np.array([np.concatenate([(g.T @ B @ g - B).ravel() for g in gens]) for B in Bs]).T
    return len(Bs) - np.linalg.matrix_rank(rows)
check("eigentliche Würfelgruppe hat 24 Elemente, alle mit det = 1", len(G24) == 24 and all(np.isclose(np.linalg.det(g), 1) for g in G24))
check("mu4 (Vierteldrehung): 2 invariante symmetrische Formen", inv_forms([R90]) == 2)
check("S4 (24 Elemente): 1 invariante Form, also delta_ij", inv_forms(G24) == 1)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Urbantke-Metrik und Simplizität")
a, b = sp.symbols("a b")
eps3 = lambda i, j, k: sp.LeviCivita(i, j, k)
eps4 = lambda *x: sp.LeviCivita(*x)
B = []
for i in range(3):
    M = sp.zeros(4, 4)
    M[0, i + 1], M[i + 1, 0] = a, -a
    for j in range(3):
        for k in range(3):
            M[j + 1, k + 1] += b * eps3(i, j, k)
    B.append(M)
gt = sp.zeros(4, 4)
idx = range(4)
for mu in idx:
    for nu in idx:
        s_ = 0
        for i, j, k in itertools.permutations(range(3)):
            e = eps3(i, j, k)
            for al, be, ga, de in itertools.permutations(range(4)):
                s_ += e * B[i][mu, al] * B[j][be, ga] * B[k][de, nu] * eps4(al, be, ga, de)
        gt[mu, nu] = sp.expand(s_)
c0 = sp.simplify(gt[0, 0] / a**3)
check("Urbantke: g~ = c diag(a^3, a b^2, a b^2, a b^2)", all(sp.simplify(gt[m, n] - (c0 * (a**3 if m == 0 else a * b**2) if m == n else 0)) == 0
                                                          for m in idx for n in idx))
say(f"   Normierung c = {c0}; g~ = {c0} diag(a^3, a b^2, a b^2, a b^2)")
wedge = lambda X, Y: sp.expand(sum(X[m, n] * Y[r, s] * eps4(m, n, r, s) for m, n, r, s in itertools.permutations(range(4))) / 4)
W = sp.Matrix(3, 3, lambda i, j: wedge(B[i], B[j]))
check("Simplizität: B^i ^ B^j = 2ab delta^ij vol (bis auf Normierung der Volumenform)",
      all(sp.simplify(W[i, j] - (W[0, 0] if i == j else 0)) == 0 for i in range(3) for j in range(3))
      and sp.simplify(W[0, 0] / (a * b)) != 0)
say(f"   B^i ^ B^j = ({sp.factor(W[0, 0])}) delta^ij")
lor = [sp.simplify((gt[m, m] / gt[1, 1]).subs({a: sp.I, b: 1})) for m in range(4)]
euk = [sp.simplify((gt[m, m] / gt[1, 1]).subs({a: 1, b: 1})) for m in range(4)]
check("a = i, b = 1: Signatur (-,+,+,+)", lor == [-1, 1, 1, 1])
check("a, b reell: euklidisch (+,+,+,+)", euk == [1, 1, 1, 1])
ar, ai, br, bi = sp.symbols("ar ai br bi", real=True)
Bc = [M.subs({a: ar + sp.I * ai, b: br + sp.I * bi}) for M in B]
Wmix = sp.expand(wedge(Bc[0], Bc[0].conjugate()))
check("Realitätsbedingung B^i ^ conj(B^j) = 0 <=> Re(a conj b) = 0", sp.simplify(Wmix - Wmix.subs({}) ) == 0
      and sp.simplify(Wmix / (ar * br + ai * bi)) .is_number)
say(f"   B^1 ^ conj(B^1) = {sp.factor(Wmix)}")
befund("Familiengruppe S4 erzwingt die isotrope Simplizität; die Lorentzsche Phase folgt nicht aus S4", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: positiver Kegel von Herm(2)")
sig = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
coords = lambda X: np.array([0.5 * np.trace(s @ X).real for s in sig])
pts = []
for _ in range(40):
    v = rng.normal(size=2) + 1j * rng.normal(size=2)
    pts.append(coords(np.outer(v, v.conj())))           # Rang-1, positiv
mono = [(i, j) for i in range(4) for j in range(i, 4)]
A = np.array([[p[i] * p[j] for i, j in mono] for p in pts])
ns = np.linalg.svd(A)[2][-(10 - np.linalg.matrix_rank(A)):]
check("quadratische Formen auf allen Rang-1-Richtungen null: Lösungsraum 1-dimensional", ns.shape[0] == 1)
Q = np.zeros((4, 4))
for (i, j), c in zip(mono, ns[0]):
    Q[i, j] += c / (1 if i == j else 2); Q[j, i] = Q[i, j]
Q /= Q[0, 0]
check("... und proportional zu diag(1, -1, -1, -1)", np.allclose(Q, np.diag([1, -1, -1, -1]), atol=1e-8))
check("Positivität X >= 0 <=> t >= |x| (Zukunftskegel)", all(
    (np.linalg.eigvalsh(sum(c * s for c, s in zip(v, sig))).min() >= -1e-12) == (v[0] >= np.linalg.norm(v[1:]) - 1e-12)
    for v in rng.normal(size=(400, 4))))

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Wardidentität für den quadratischen Operator")
q = sp.symbols("q0:4")
q2 = sum(x**2 for x in q)
d = lambda i, j: 1 if i == j else 0
cs = sp.symbols("c1:6")
m1, m2 = sp.symbols("m1 m2")
def K(mu, nu, rh, si):
    I = sp.Rational(1, 2) * (d(mu, rh) * d(nu, si) + d(mu, si) * d(nu, rh))
    return (cs[0] * q2 * I + cs[1] * q2 * d(mu, nu) * d(rh, si)
            + cs[2] * (q[mu] * q[rh] * d(nu, si) + q[mu] * q[si] * d(nu, rh) + q[nu] * q[rh] * d(mu, si) + q[nu] * q[si] * d(mu, rh))
            + cs[3] * (q[mu] * q[nu] * d(rh, si) + q[rh] * q[si] * d(mu, nu))
            + cs[4] * q[mu] * q[nu] * q[rh] * q[si] / q2 + m1 * I + m2 * d(mu, nu) * d(rh, si))
eqs = []
for nu, rh, si in itertools.product(range(4), repeat=3):
    expr = sp.together(sp.expand(sum(q[mu] * K(mu, nu, rh, si) for mu in range(4)) * q2))
    num = sp.numer(expr)
    eqs += sp.Poly(sp.expand(num), *q).coeffs()
sol = sp.solve(eqs, list(cs) + [m1, m2], dict=True)
say(f"   Lösungen: {sol}")
free = [s for s in (list(cs) + [m1, m2]) if s not in sol[0]] if sol else []
local = {k: v for k, v in sol[0].items()} if sol else {}
check("Wardidentität erzwingt m1 = m2 = 0 (keine impulsunabhängigen Terme)", sol and local.get(m1, None) == 0 and local.get(m2, None) == 0)
# lokale Lösung (c5 = 0): genau eine Richtung
sol_loc = sp.solve(eqs + [cs[4]], list(cs) + [m1, m2], dict=True)
free_loc = [s for s in cs[:4] if s not in sol_loc[0]]
check("lokal (c5 = 0): genau ein freier Gesamtkoeffizient (linearisiertes Einstein-Hilbert, Fierz-Pauli masselos)",
      len(free_loc) == 1)
fp = {s: sp.simplify(sol_loc[0].get(s, s) / free_loc[0]) for s in cs[:4]}
say(f"   Verhältnisse c1 : c2 : c3 : c4 = {[fp[s] for s in cs[:4]]}")
befund("Die Wardidentität legt die Spin-2-Struktur fest, aber nicht den Gesamtkoeffizienten (also nicht G)", len(free_loc) == 1)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: nichtminimale Kopplung")
xi = sp.symbols("xi")
a1 = sp.Rational(1, 6) - xi
say(f"   Wärmekern a_1 = ({a1}) R für (-Box + m^2 + xi R): induziertes 1/(16 pi G) proportional zu {a1}")
check("xi R phi^2 ist diffeomorphismusinvariant und verschiebt G (a_1 hängt von xi ab)", sp.diff(a1, xi) != 0)
z27 = json.loads(Path(__file__).with_name("z27_ergebnis.json").read_text())
xi_star = z27["zahlen"]["xi_stern_skalar"]
say(f"   Z27: R = 1,42 des Skalars entspricht xi = {xi_star:.4f}; dieser Teil der Kopplungsmehrdeutigkeit ist kovariant")
check("Z27-Äquivalent xi* liegt nahe 0 (kovarianter, Ward-verträglicher Anteil der Mehrdeutigkeit)", abs(xi_star) < 0.05)

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: longitudinaler Test (reine Koordinatenänderung) auf dem Z27-Gitter, Skalar, Kopplung K")
phat = lambda p: 2 * np.sin(p / 2)
def c_long(m, N=96):
    p = (np.arange(N) + 0.5) * 2 * np.pi / N - np.pi
    px, py, pz = np.meshgrid(p, p, p, indexing="ij")
    w = np.sqrt(m**2 + phat(px)**2 + phat(py)**2 + phat(pz)**2)
    return np.mean(3 * phat(pz)**2 / w - phat(pz)**4 / w**3) / 8
cl = {m: c_long(m) for m in (0.0, 0.5, 1.0)}
cl_fine = c_long(0.0, 160)
say(f"   analytisch C_long = (1/8) <3 p_z^2/w - p_z^4/w^3>: m = 0: {cl[0.0]:.6f} (N = 96), {cl_fine:.6f} (N = 160); "
    f"m = 0,5: {cl[0.5]:.6f}; m = 1: {cl[1.0]:.6f}")
def e2_lattice(Lz, Lp=24, m=0.0, eps=1e-3, j=1):
    k = 2 * np.pi * j / Lz
    zs = np.arange(Lz)
    tot = []
    for e in (eps, -eps, 0.0):
        sq = 1 + e * np.cos(k * zs)                       # sqrt(h) an den Plätzen
        wl = 1 / (1 + e * np.cos(k * (zs + 0.5)))          # sqrt(h) h^zz an der Kantenmitte (Regel K)
        Lap = np.zeros((Lz, Lz))
        for n_ in range(Lz):
            n2 = (n_ + 1) % Lz
            Lap[n_, n_] += wl[n_]; Lap[n2, n2] += wl[n_]; Lap[n_, n2] -= wl[n_]; Lap[n2, n_] -= wl[n_]
        E = 0.0
        p = 2 * np.pi * np.arange(Lp) / Lp
        for px in p:
            for py in p:
                perp = phat(px)**2 + phat(py)**2 + m**2
                Kmat = Lap + np.diag(sq * perp)
                S = np.diag(1 / np.sqrt(sq))
                ev = np.linalg.eigvalsh(S @ Kmat @ S)
                E += 0.5 * np.sqrt(np.clip(ev, 0, None)).sum()
        tot.append(E / (Lz * Lp * Lp))
    return (tot[0] + tot[1] - 2 * tot[2]) / (2 * eps**2)      # Koeffizient von eps^2 (halbe zweite Ableitung)
lat = {Lz: e2_lattice(Lz) for Lz in (16, 32, 64)}
lat05 = e2_lattice(64, m=0.5)
say(f"   Gitter m = 0,5 (L_z = 64): {lat05:.6f} gegen analytisch {cl[0.5]:.6f} (Wert der externen Analyse 0,247893)")
say(f"   Gitter (L_perp = 24, Modulation j = 1): Koeffizient von eps^2 je Platz {dict((k, round(v, 6)) for k, v in lat.items())}")
ext = lat[64] + (lat[64] - lat[32]) / 3                   # Richardson in k^2
say(f"   Extrapolation k -> 0 (Richardson): {ext:.6f}")
ZAHLEN["C_long"] = {"m0_N96": cl[0.0], "m0_N160": cl_fine, "m05": cl[0.5], "m1": cl[1.0], "gitter": lat, "extrapoliert": ext, "gitter_m05": lat05}
check("Gitterantwort auf reine Koordinatenänderung ist positiv und konvergiert gegen C_long (auf 2 %)",
      all(v > 0 for v in lat.values()) and abs(ext - cl_fine) / cl_fine < 0.02 and abs(lat05 - cl[0.5]) / cl[0.5] < 0.02)
check("C_long > 0 für m = 0; 0,5; 1 (keine kovariante Regel K)", all(v > 0 for v in cl.values()))
befund("Festes Gitter verletzt die longitudinale Wardidentität schon bei k -> 0 (Abschneideanisotropie)", True)
say("   Lesart: Bei festen Koordinatenkanten ändert die Koordinatenänderung den eigentlichen Gitterabstand; die")
say("   Vakuumenergie hängt von ihm ab. Kein kovarianter Gegenterm (sqrt(g) Lambda ist invariant) hebt das auf;")
say("   es braucht nichtkovariante Gegenterme oder Kanten, die die Metrik selbst definiert.")

urteil = ("Die Familiengruppe S4 lässt auf Lambda^2 F nur die isotrope Form delta_ij zu (mu4 allein zwei), also die "
          "Plebanski-Simplizität B^i ^ B^j ~ delta^ij. Für B^i = a dt^dx^i + (b/2) eps dx^dx ist die Urbantke-Metrik "
          "c diag(a^3, ab^2, ab^2, ab^2): reell euklidisch, mit a = i lorentzsch; die Realitätsbedingung "
          "B^i ^ conj B^j = 0 heißt Re(a conj b) = 0 und folgt nicht aus S4. Der positive Kegel von Herm(2) legt die "
          "Minkowskiform eindeutig fest (eine quadratische Form verschwindet auf allen Rang-1-Richtungen nur als Vielfaches "
          "von det). Die Wardidentität erzwingt verschwindende Massenterme und die Fierz-Pauli-Struktur mit genau einem "
          "freien Gesamtkoeffizienten; xi R phi^2 ist kovariant und verschiebt G, also legt die Wardidentität G nicht fest. "
          "Auf dem festen Z27-Gitter verletzt schon eine reine Koordinatenänderung die Wardidentität: Die Vakuumantwort "
          "konvergiert gegen C_long = (1/8)<3 p_z^2/w - p_z^4/w^3> > 0. Für Gravitation braucht H6 also Kanten, die die "
          "Metrik selbst definiert, und eine eigene Regel für den Gesamtkoeffizienten.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
Path(__file__).with_name("z29_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
