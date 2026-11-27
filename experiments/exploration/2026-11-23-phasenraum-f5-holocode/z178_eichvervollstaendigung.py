"""Z178: F1, die Naht als notwendige Eichvervollständigung. Schließt die Lücke 'Superauswahl ist noch keine Eichung'.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Bausteine bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, N numerisch, C bedingt, L Lesart, O offen.
Bestand: z170 (27 Spiegelvervollständigungen der 16 Uhrsektoren, A = Z4 x Z4, q = (5 x^2 + 3 y^2)/8), z176 (die Eichsumme
der Paaruhr über E8 (x) E8-quer ist die diagonale Vervollständigung; K19 und K20 folgen aus P_std und A* nur unter der
Lesart 'Superauswahl der Paaruhr = Eichsumme'; einseitige Uhr: Anomalie der Ordnung 2), z174/z177 (Modultest: gerade
Nahtzahl, normale Statistik nur für zwei Paare, N = 1..12).

Teil A  Sektorebene, exakt: die uhrneutrale Algebra N = (E8 (x) E8-quer)^U, ihre Sektoren N^perp/N, alle Vervollständigungen
        über N, ihre Fixpunktalgebren und Ladungen, die Eindeutigkeit der ladungsfreien Erweiterung, die Modulgleichung.
Teil B  Standardform, endliches Modell: GNS-Raum der superselektierten Algebra = U-invarianter Teilraum; die Projektion
        E(X) = (1/4) Summe U^k X U^-k ist die Jones-Projektion; Gegenmodell 'Superauswahl ohne Eichung'.
Teil C  Assoziator und Anomalie: expliziter abelscher 3-Kozykel der 16 Sektoren (Pentagon, zwei Hexagone), Anomalie der
        einseitigen Uhr als Invariante, Additivität (Nahtzahl gerade für jedes N), kommutative Algebra auf den
        Vervollständigungen.
Teil D  Dieselbe Summe als Spinstruktursumme: unabhängige GSO-Summen der beiden Deckhälften (10 | 6 reelle Eltern).
Teil E  Nahtzahl über den Testbereich hinaus: Modultest aus z177 für N bis 24.
Teil F  Bilanz: was aus P_std, A* und der GNS-Regel folgt, und was Annahme bleibt.
"""
import contextlib
import io
import json
import sys
import time
from fractions import Fraction as Fr
from itertools import product
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[2] / "research" / "originlab"))
import clock_sector_search as cs  # noqa: E402  (nur lesend)

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T_START = time.time()


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def zeit():
    return f"(Laufzeit bisher {time.time() - T_START:.0f} s)"


# =================================================================== Teil A
say("Teil A: die uhrneutrale Algebra und ihre Vervollständigungen (exakt)")
A = cs.A
G = cs.G
ZERO = ((0, 0), (0, 0))
COMP = [frozenset(s) for s in cs.isotropic_subgroups() if len(s) == 16]
ISO = [frozenset(s) for s in cs.isotropic_subgroups()]


def glue(s):
    return [(k % 4, (s * k) % 4) for k in range(4)]


def neutral(s, t):
    """Uhrneutrale Algebra: Sektoren (k, s k) (x) (k, t k) von E8 (x) E8-quer mit Paaruhr U = i^(k - k')."""
    return frozenset(((k % 4, (s * k) % 4), (k % 4, (t * k) % 4)) for k in range(4))


def e8e8(s, t):
    return frozenset((a, b) for a in glue(s) for b in glue(t))


def perp(H):
    return frozenset(g for g in G if all(cs.B_lr(g, h) == 0 for h in H))


def typ(L):
    left = {g[0] for g in L if g[1] == (0, 0)}
    right = {g[1] for g in L if g[0] == (0, 0)}
    spin = {cs.Q_lr(g) for g in L}
    art = "Graph" if len(left) == 1 and len(right) == 1 else ("faktorisiert" if len(left) == 4 and len(right) == 4
                                                                else f"gemischt {len(left)}x{len(right)}")
    return art, "bosonisch" if spin == {Fr(0)} else "Spin"


BITS = [(s, t) for s in (1, 3) for t in (1, 3)]


def flux(a, s):
    """Uhrfluss eines Sektors a relativ zum E8-Kleber {(k, s k)}: Monodromie b(a, (1, s)) = flux/4."""
    return (a[0] - s * a[1]) % 4


tab_a1 = {}
for s, t in BITS:
    N = neutral(s, t)
    ueber = [L for L in COMP if N <= L]
    unt = e8e8(s, t)
    rows = []
    for L in ueber:
        u = L & unt
        ladungen = sorted({(g[0][0] - g[1][0]) % 4 for g in u})
        rows.append({"typ": typ(L), "unverdrillt": len(u), "Ladungen_U": ladungen, "ladungsfrei": ladungen == [0],
                     "Flüsse": sorted({flux(g[0], s) for g in L}),
                     "Ladung bei Fluss 1": sorted({(g[0][0] - g[1][0]) % 4 for g in L if flux(g[0], s) == 1}),
                     "Spin bei Fluss 1": sorted({str(cs.Q_lr(g)) for g in L if flux(g[0], s) == 1})})
    tab_a1[(s, t)] = rows
alle_iso = all(all(cs.Q_lr(g) == 0 for g in neutral(s, t)) for s, t in BITS)
fuenf = all(len(v) == 5 for v in tab_a1.values())
groessen = all(sorted(r["unverdrillt"] for r in v if r["typ"][1] == "bosonisch") == [4, 8, 16]
               and sorted(r["unverdrillt"] for r in v if r["typ"][1] == "Spin") == [4, 8] for v in tab_a1.values())
frei_typen = {str(k): sorted(r["typ"] for r in v if r["ladungsfrei"]) for k, v in tab_a1.items()}
frei_ok = all(v == [("Graph", "Spin"), ("Graph", "bosonisch")] for v in frei_typen.values())
fluss1 = all({(r["typ"][1], tuple(r["Ladung bei Fluss 1"]), tuple(r["Spin bei Fluss 1"])) for r in v if r["ladungsfrei"]}
             == {("bosonisch", (0,), ("0",)), ("Spin", (2,), ("1/2",))} for v in tab_a1.values())
typen = {str(k): sorted((r["typ"], r["unverdrillt"], tuple(r["Ladungen_U"]), tuple(r["Flüsse"])) for r in v)
         for k, v in tab_a1.items()}
ZAHLEN["A1"] = {"Vervollständigungen über N je Spiegelbit": {str(k): len(v) for k, v in tab_a1.items()}, "Typen": typen}
check("A1 [E] Die uhrneutrale Algebra N = (E8 (x) E8-quer)^U (Sektoren (k, s k) (x) (k, t k), Paaruhr U = i^(k - k')) ist für "
      "alle vier Spiegelbits (s, t) eine isotrope Untergruppe der Ordnung 4. Unter den 27 Vervollständigungen enthalten "
      "genau fünf die Algebra N. Drei sind bosonisch: E8 (x) E8-quer selbst (unverdrillter Teil 16 Sektoren, Ladungen 0..3, "
      "kein Uhrfluss), eine Zwischenform (8 Sektoren, Ladungen 0 und 2, Flüsse 0 und 2: nur U^2 = (-1)^(F+F-quer) geeicht) "
      "und ein Graph (unverdrillter Teil = N, alle vier Flüsse). Zwei sind Spin-Formen (halbzahliger Spin im verdrillten "
      "Teil). Kein Feld mit Paaruhr-Ladung im unverdrillten Teil tragen genau zwei, beide Graphen mit allen vier Flüssen: "
      "der bosonische (Fluss 1 mit Ladung 0, Spin ganz) und eine Spin-Form (Fluss 1 mit Ladung 2, Spin 1/2)",
      alle_iso and fuenf and groessen and frei_ok and fluss1)

# A2: Sektoren von N und ihre Lagrange-Untergruppen
s0, t0 = 1, 1
N0 = neutral(s0, t0)
P0 = perp(N0)
cosets = sorted({frozenset(cs.g_add(g, h) for h in N0) for g in P0}, key=lambda c: sorted(c))
qv = {c: cs.Q_lr(next(iter(c))) for c in cosets}
assert all(len({cs.Q_lr(g) for g in c}) == 1 for c in cosets)
cidx = {c: i for i, c in enumerate(cosets)}
rep = [sorted(c)[0] for c in cosets]


def coset_of(g):
    return cidx[frozenset(cs.g_add(g, h) for h in N0)]


# Erzeuger: e = Ladung (unverdrillt, U = i), m = Fluss (uhrverdrillt, Spin 0)
unt0 = e8e8(s0, t0)
e_kand = [c for c in cosets if all(g in unt0 for g in c) and (next(iter(c))[0][0] - next(iter(c))[1][0]) % 4 == 1]
m_kand = [c for c in cosets if qv[c] == 0 and not any(g in unt0 for g in c)]
basis = None
for ce in e_kand:
    for cm in m_kand:
        ge, gm = sorted(ce)[0], sorted(cm)[0]
        if cs.B_lr(ge, gm) == Fr(1, 4) and cs.Q_lr(ge) == 0:
            basis = (ge, gm)
            break
    if basis:
        break
ge, gm = basis


def mult(g, n):
    out = ZERO
    for _ in range(n % 4):
        out = cs.g_add(out, g)
    return out


dz4 = {}
for x in range(4):
    for y in range(4):
        dz4[(x, y)] = coset_of(cs.g_add(mult(gm, x), mult(ge, y)))
hyper = len(set(dz4.values())) == 16 and all(qv[cosets[dz4[(x, y)]]] == Fr(x * y, 4) % 1 for x in range(4) for y in range(4))
lag16 = [L for L in ISO if len(L) == 16 and N0 <= L]
lag16_bos = [L for L in lag16 if typ(L)[1] == "bosonisch"]
einseitig = frozenset(((0, 0), (k % 4, (t0 * k) % 4)) for k in range(4))
lag_eins = [L for L in COMP if einseitig <= L]
eins_typen = sorted(typ(L) for L in lag_eins)
eins_fluss = [sorted({flux(g[0], 1) for g in L}) for L in lag_eins]
ZAHLEN["A2"] = {"Sektoren von N": len(cosets), "Basis (Ladung e, Fluss m)": [str(ge), str(gm)],
                "einseitig: Vervollständigungen": eins_typen, "einseitig: Flüsse": eins_fluss}
check("A2 [E] Sektoren: N^perp/N hat 16 Elemente und ist die hyperbolische Form q(x, y) = x y/4 auf Z4 x Z4 (x = Fluss, "
      "y = Ladung der Paaruhr), also das Doppel D(Z4) ohne Anomalie. Seine Lagrange-Untergruppen mit q = 0 sind genau drei: "
      "{(0, y)} (alle Ladungen, kein Fluss: E8 (x) E8-quer), {(x, 0)} (alle Flüsse, keine Ladung: die Eichsumme) und "
      "{(2a, 2b)}; dazu kommen die zwei Spin-Formen <(1, 2)> und <(2, 1)> (q = 1/2). Die Eichsumme ist das einzige "
      "bosonische Komplement der Ladungen: {(x, c x)} hat q = 0 nur für c = 0. Gegenprobe einseitige Uhr "
      "(N_1 = E8^R (x) E8-quer): über N_1 liegen drei Vervollständigungen, alle faktorisiert, und keine enthält einen "
      "Sektor mit ungeradem Uhrfluss. Einseitig lässt sich R^2 = (-1)^F eichen, R nicht (Anomalie der Ordnung 2)",
      len(cosets) == 16 and hyper and len(lag16) == 5 and len(lag16_bos) == 3
      and [c for c in range(4) if all(Fr(c * x * x, 4) % 1 == 0 for x in range(4))] == [0]
      and len(lag_eins) == 3 and all(t[0] == "faktorisiert" for t in eins_typen)
      and not any(f % 2 for fl in eins_fluss for f in fl))

# A3: Modulgleichung auf den 16 Sektoren von N
zeta = sp.I
Smat = sp.Matrix(16, 16, lambda i, j: sp.Rational(1, 4) * zeta ** int(-4 * cs.B_lr(rep[i], rep[j]) % 4))
Tmat = sp.diag(*[zeta ** int(4 * qv[c] % 4) for c in cosets])
s_ok = sp.simplify(Smat * Smat.H - sp.eye(16)) == sp.zeros(16) and sp.simplify((Smat * Tmat) ** 3 - Smat ** 2) == sp.zeros(16)
inv = (Smat - sp.eye(16)).col_join(Tmat - sp.eye(16)).nullspace()
ind = {name: sp.Matrix([1 if i in {dz4[k] for k in keys} else 0 for i in range(16)])
       for name, keys in (("Ladungen", [(0, y) for y in range(4)]), ("Flüsse", [(x, 0) for x in range(4)]),
                          ("gerade", [(2 * a, 2 * b) for a in range(2) for b in range(2)]))}
spann = sp.Matrix.hstack(*ind.values()).rank() == 3 and sp.Matrix.hstack(*inv).rank() == 3 \
    and sp.Matrix.hstack(*inv, *ind.values()).rank() == 3
# Zwang: Vakuum einmal, kein unverdrillter Sektor mit Ladung
a_, b_, c_ = sp.symbols("a b c")
n = a_ * ind["Ladungen"] + b_ * ind["Flüsse"] + c_ * ind["gerade"]
lsg = sp.solve([n[dz4[(0, 0)]] - 1] + [n[dz4[(0, y)]] for y in (1, 2, 3)], [a_, b_, c_], dict=True)
inv2 = (Smat - sp.eye(16)).col_join(Tmat ** 2 - sp.eye(16)).nullspace()
ind5 = [sp.Matrix([1 if frozenset(cs.g_add(rep[i], h) for h in N0) <= L else 0 for i in range(16)]) for L in lag16]
spann5 = sp.Matrix.hstack(*inv2, *ind5).rank() == len(inv2) == sp.Matrix.hstack(*ind5).rank()
frei5 = [k for k, v in enumerate(ind5) if v[dz4[(0, 0)]] == 1 and all(v[dz4[(0, y)]] == 0 for y in (1, 2, 3))]
ZAHLEN["A3"] = {"Dimension des invarianten Raums": len(inv), "Lösung": str(lsg), "Dimension nur S und T^2": len(inv2),
                "ladungsfreie Lagrange-Lösungen ohne T": len(frei5)}
check("A3 [E] Modulgleichung: Auf den 16 Sektoren von N wirken S = (1/4) e^(-2 pi i b) und T = e^(2 pi i q) (S unitär, "
      "(S T)^3 = S^2). Der Raum der S- und T-invarianten Vektoren hat die Dimension 3 und wird von den drei "
      "Lagrange-Untergruppen aufgespannt. Mit 'Vakuum genau einmal' und 'kein unverdrillter Sektor mit Paaruhr-Ladung' "
      "bleibt genau eine Lösung: jeder Fluss x genau einmal, mit Ladung 0. Die Vielfachheiten der uhrverdrillten Sektoren "
      "sind damit nicht wählbar. Ohne T (nur S und T^2, also mit fester Spinstruktur) hat der invariante Raum die Dimension "
      f"{len(inv2)}, aufgespannt von allen fünf Vervollständigungen, und es gibt zwei ladungsfreie: Die Spin-Form ist das "
      "Gegenmodell, das erst T ausschließt", s_ok and len(inv) == 3 and spann and lsg == [{a_: 0, b_: 1, c_: 0}]
      and len(inv2) == 5 and spann5 and len(frei5) == 2)

# A4: Fixpunktalgebren aller 27 Vervollständigungen
def neutral_gleich(s, t):
    """Neutrale Algebra der gleichsinnigen Paaruhr U' = i^(k + k') (Ladungskonjugation auf einer Seite)."""
    return frozenset(((k % 4, (s * k) % 4), ((-k) % 4, (-t * k) % 4)) for k in range(4))


fix_tab = []
for L in sorted(COMP, key=lambda s: sorted(s)):
    unter = [H for H in ISO if H <= L and len(H) == 4 and all(cs.Q_lr(g) == 0 for g in H)]
    cyc = [H for H in unter if any(len({mult(g, k) for k in range(4)}) == 4 for g in H)]      # zyklisch Z4
    fix_tab.append({"typ": typ(L), "Index-4-Unteralgebren": len(unter), "zyklisch": len(cyc),
                    "N von U": [str(b) for b in BITS if neutral(*b) <= L],
                    "N von U'": [str(b) for b in BITS if neutral_gleich(*b) <= L]})
je_bit = {str(b): sum(neutral(*b) <= L for L in COMP) for b in BITS}
n_enth = sum(bool(r["N von U"]) for r in fix_tab)
bos_graphen = [r for r in fix_tab if r["typ"] == ("Graph", "bosonisch")]
ZAHLEN["A4"] = {"Tabelle": fix_tab, "Vervollständigungen, die ein N enthalten": n_enth, "je Spiegelbit": je_bit}
check("A4 [E] Fixpunktalgebren aller 27 Vervollständigungen: Jede hat die duale Symmetriegruppe der Ordnung 16; ihre "
      "bosonischen Index-4-Unteralgebren sind die Untergruppen der Ordnung 4 mit q = 0. Jede uhrneutrale Algebra N(s, t) "
      f"liegt in genau fünf der 27 (zusammen {n_enth}). Von den vier bosonischen Graphen (jeder Sektor einmal, z170) sind "
      "zwei die Eichsummen der Paaruhr U = i^(k - k') für je zwei Spiegelbits; die anderen zwei sind die Eichsummen der "
      "gleichsinnigen Uhr U' = i^(k + k'), also dasselbe mit Ladungskonjugation auf einer Seite. Jeder bosonische Graph "
      "ist damit Eichsumme, und keine zweite bosonische Vervollständigung mit allen 16 Sektoren teilt mit ihm eine "
      "neutrale Algebra",
      set(je_bit.values()) == {5} and len(bos_graphen) == 4
      and sorted((len(r["N von U"]), len(r["N von U'"])) for r in bos_graphen) == [(0, 2), (0, 2), (2, 0), (2, 0)])

# A5: Zustand bei tau = i: Gewichte
mp.mp.dps = 40
tau = mp.mpc(0, 1)
qn = mp.exp(mp.pi * 1j * tau)
th2, th3, th4 = mp.jtheta(2, 0, qn), mp.jtheta(3, 0, qn), mp.jtheta(4, 0, qn)
eta = mp.exp(2j * mp.pi * tau / 24) * mp.qp(mp.exp(2j * mp.pi * tau))


def so_char(nn):
    return [(th3 ** nn + th4 ** nn) / (2 * eta ** nn), th2 ** nn / (2 * eta ** nn), (th3 ** nn - th4 ** nn) / (2 * eta ** nn),
            th2 ** nn / (2 * eta ** nn)]                                    # Klassen 0 (o), 1 (s), 2 (v), 3 (c)


c10, c6 = so_char(5), so_char(3)
chi = {(x, y): c10[x] * c6[y] for x in range(4) for y in range(4)}
z_e8 = sum(chi[a] for a in glue(1))
z_diag = sum(abs(v) ** 2 for v in chi.values())
z_unt = sum(abs(chi[a]) ** 2 for a in glue(1))
gew_verdr = sum(abs(chi[a]) ** 2 for a in A if flux(a, 1) % 2) / z_diag
gew_f2 = sum(abs(chi[a]) ** 2 for a in A if flux(a, 1) == 2) / z_diag
ziel = (12 * mp.sqrt(2) - 13) / 14
ZAHLEN["A5"] = {"Z_E8(i)": mp.nstr(z_e8.real, 20), "Z_diag(i)": mp.nstr(z_diag, 25), "Sp über H^U": mp.nstr(z_unt, 20),
                "Gewicht ungerader Fluss": mp.nstr(gew_verdr, 20), "Gewicht Fluss 2": mp.nstr(gew_f2, 20)}
check("A5 [N, 40 Stellen] Zustand bei tau = i: Die E8-Naht hat Z = 12; der GNS-Raum H^U der neutralen Algebra allein hat "
      f"die Spur {mp.nstr(z_unt, 12)} und ist nicht S-invariant; die Eichvervollständigung hat Z = 80 + 48 sqrt(2) = "
      f"{mp.nstr(z_diag, 14)} und trägt in den Sektoren mit ungeradem Uhrfluss das Gewicht (12 sqrt(2) - 13)/14 = "
      f"{mp.nstr(gew_verdr, 12)} (z176); jeder der drei Flüsse 1, 2, 3 trägt die Hälfte davon",
      abs(z_e8 - 12) < mp.mpf(10) ** -30 and abs(z_diag - (80 + 48 * mp.sqrt(2))) < mp.mpf(10) ** -28
      and abs(gew_verdr - ziel) < mp.mpf(10) ** -30 and abs(2 * gew_f2 - ziel) < mp.mpf(10) ** -30)
befund("Teil A: Über der uhrneutralen Algebra gibt es genau fünf Vervollständigungen, drei bosonische und zwei Spin-Formen. "
       "Genau eine ist S- und T-invariant und ohne Feld mit Paaruhr-Ladung. Sie ist die Eichsumme; ihre Vielfachheiten "
       "sind durch S und T festgelegt", True)

# =================================================================== Teil B
say("\nTeil B: Standardform der superselektierten Naht (endliches Modell)")
rng = np.random.default_rng(178)


def modell(dims, beta_pi=True):
    """Naht mit Uhr R = i^Q, Q-Eigenräume der Dimensionen dims[k]; rho uhrinvariant, sonst generisch."""
    d = sum(dims)
    Q = np.concatenate([[k] * m for k, m in enumerate(dims)]).astype(int)
    R = np.diag(1j ** Q)
    blocks = []
    for m in dims:
        X = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))
        blocks.append(X @ X.conj().T + 0.3 * np.eye(m))
    rho = np.zeros((d, d), dtype=complex)
    o = 0
    for Bk in blocks:
        m = Bk.shape[0]
        rho[o:o + m, o:o + m] = Bk
        o += m
    return d, Q, R, rho / np.trace(rho).real


def wurzel(M):
    w, v = np.linalg.eigh((M + M.conj().T) / 2)
    return (v * np.sqrt(np.clip(w, 0, None))) @ v.conj().T


b_ok, b_zahlen = True, {}
for dims in ((2, 1, 2, 1), (3, 2, 2, 1), (1, 3, 1, 2)):
    d, Q, R, rho = modell(dims)
    Om = wurzel(rho)                                    # Standardvektor als Matrix; H (x) H-quer = Matrizen
    Uop = lambda X: R @ X @ R.conj().T                  # Paaruhr R (x) R-quer^-1
    Eop = lambda X: sum(np.linalg.matrix_power(R, k) @ X @ np.linalg.matrix_power(R, k).conj().T for k in range(4)) / 4
    # (1) Omega ist U-invariant; die U-Ladung von |k><l| ist k - l
    inv_om = np.allclose(Uop(Om), Om)
    # (2) Jones-Projektion: der Projektor auf den Abschluss von M^R Omega ist P_U = (1/4) Summe U^k
    mats = []
    for a in range(d):
        for b in range(d):
            if Q[a] == Q[b]:
                Eab = np.zeros((d, d), dtype=complex)
                Eab[a, b] = 1
                mats.append((Eab @ Om).reshape(-1))
    span_dim = np.linalg.matrix_rank(np.array(mats), tol=1e-10)
    # Projektor P_U als Matrix auf d^2
    PU = np.zeros((d * d, d * d), dtype=complex)
    for a in range(d):
        for b in range(d):
            Eab = np.zeros((d, d), dtype=complex)
            Eab[a, b] = 1
            PU[:, a * d + b] = Eop(Eab).reshape(-1)
    Qm, _ = np.linalg.qr(np.array(mats).T)
    Pspan = Qm[:, :span_dim] @ Qm[:, :span_dim].conj().T
    jones = np.allclose(Pspan, PU, atol=1e-9) and span_dim == sum(m * m for m in dims)
    # (3) e x e = E(x) e für zufälliges x (x wirkt von links)
    X = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    Lx = np.kron(X, np.eye(d))
    LE = np.kron(Eop(X), np.eye(d))
    bim = np.allclose(PU @ Lx @ PU, LE @ PU, atol=1e-9)
    # (4) Kegelvektoren: xi_phi = sigma^(1/2); U-invariant genau für uhrinvariante Zustände
    sig = modell(dims)[3]
    Y = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    sig2 = Y @ Y.conj().T
    sig2 /= np.trace(sig2).real
    kegel = np.allclose(Uop(wurzel(sig)), wurzel(sig)) and not np.allclose(Uop(wurzel(sig2)), wurzel(sig2))
    gel = np.allclose(wurzel(Eop(sig2)), Eop(wurzel(Eop(sig2))))
    # (5) Der U-invariante Teil des Kegels spannt H^U auf
    pos = []
    for _ in range(4 * d * d):
        s_ = modell(dims)[3]
        pos.append(wurzel(s_).reshape(-1))
    spann_kegel = np.linalg.matrix_rank(np.array(pos), tol=1e-9) == sum(m * m for m in dims)
    # (6) GNS der vollen neutralen Algebra B(H (x) H-quer)^U auf Omega
    erz = []
    for _ in range(3 * d * d):
        Z1 = rng.normal(size=(d * d, d * d)) + 1j * rng.normal(size=(d * d, d * d))
        Umat = np.kron(R, R.conj())
        Zn = sum(np.linalg.matrix_power(Umat, k) @ Z1 @ np.linalg.matrix_power(Umat, k).conj().T for k in range(4)) / 4
        erz.append(Zn @ Om.reshape(-1))
    gns = np.linalg.matrix_rank(np.array(erz), tol=1e-8) == sum(m * m for m in dims)
    # (7) im Quarterzustand ist die Paaruhr ein Element des Modulflusses
    rho_q = np.diag(np.exp(-np.pi * Q.astype(float)))
    rho_q /= np.trace(rho_q)
    w = np.diag(rho_q).astype(complex)
    mod = np.allclose((w[:, None] / w[None, :]) ** (-0.5j), 1j ** (Q[:, None] - Q[None, :]))   # Delta^(-i/2) = U
    b_zahlen[str(dims)] = {"dim H (x) H-quer": d * d, "dim H^U": int(span_dim)}
    b_ok = b_ok and inv_om and jones and bim and kegel and gel and spann_kegel and gns and mod
ZAHLEN["B1"] = b_zahlen
check("B1 [E, drei Modelle] Standardform der Naht (H (x) H-quer, Omega = vec(rho^(1/2)), Paaruhr U = R (x) R-quer^-1). "
      "(1) U Omega = Omega. (2) Der Abschluss von M^R Omega, also der GNS-Raum der uhrneutralen Algebra, ist genau der "
      "U-invariante Teilraum H^U = (+)_k H_k (x) H_k-quer; die Projektion darauf ist P_U = (1/4) Summe U^k (Jones-Projektion "
      "der Inklusion M^R in M). (3) P_U x P_U = E(x) P_U mit E(X) = (1/4) Summe R^k X R^-k: Die bedingte Erwartung ist die "
      "Kompression auf H^U. (4) Ein Kegelvektor ist genau dann U-invariant, wenn sein Zustand uhrinvariant ist, und die "
      "U-invarianten Kegelvektoren spannen H^U auf. (5) Auch die volle Algebra der Paaruhr-neutralen Operatoren auf "
      "H (x) H-quer erzeugt aus Omega genau H^U. (6) Im Quarterzustand rho ~ e^(-pi Q) ist U = Delta^(-i/2), also ein "
      "Element des Modulflusses (z130, z174): E existiert kanonisch und erhält den Zustand (Takesaki)", b_ok)

# Gegenmodell: Superauswahl ohne Eichung
unt = e8e8(1, 1)
N = neutral(1, 1)
modelle = []
for L in COMP:
    u = L & unt
    modelle.append({"bos": typ(L)[1] == "bosonisch", "enthält_N": N <= L,
                    "ladungsfrei": all((g[0][0] - g[1][0]) % 4 == 0 for g in u), "GNS": u == N})
zahl = lambda f: [sum(1 for m in modelle if f(m)), sum(1 for m in modelle if f(m) and m["bos"])]
gm = {"S-invariant [alle, davon auch T-invariant]": zahl(lambda m: True),
      "kein lokales Feld mit Paaruhr-Ladung": zahl(lambda m: m["ladungsfrei"]),
      "enthält die neutrale Algebra": zahl(lambda m: m["enthält_N"]),
      "unverdrillter Teil = GNS-Raum H^U": zahl(lambda m: m["GNS"])}
ZAHLEN["B2"] = gm
check("B2 [E] Gegenmodelle auf der Sektorebene (Spiegelbit (1, 1); Zahlen: alle S-invarianten, davon auch T-invariant): "
      f"Modulinvarianz allein lässt {gm['S-invariant [alle, davon auch T-invariant]']} Vervollständigungen zu. 'Kein lokales "
      f"Feld mit Paaruhr-Ladung' erfüllen {gm['kein lokales Feld mit Paaruhr-Ladung']}, 'enthält die neutrale Algebra' "
      f"{gm['enthält die neutrale Algebra']}, 'unverdrillter Teil = GNS-Raum H^U' {gm['unverdrillter Teil = GNS-Raum H^U']}. "
      "Zwei Gegenmodelle sind damit benannt. (a) 'E8 (x) E8-quer mit Auswahlregel' erfüllt A* und die Superauswahl, sein "
      "Zustandsraum ist aber nicht der GNS-Raum des Zustands über den Observablen: Omega ist dort nicht zyklisch (B1). "
      "(b) Die Spin-Form erfüllt Superauswahl, GNS-Regel und S, aber nicht T (Fluss 1 mit Spin 1/2). Die Lücke zwischen "
      "Superauswahl und Eichung schließen also die GNS-Regel des Bestands und die volle Basiswechselsymmetrie aus A*",
      gm["S-invariant [alle, davon auch T-invariant]"] == [27, 10] and gm["enthält die neutrale Algebra"] == [5, 3]
      and gm["unverdrillter Teil = GNS-Raum H^U"] == [2, 1])
befund("Teil B: Superauswahl der Paaruhr plus GNS-Regel gibt den Gauß-Zwang U = 1 auf dem Zustandsraum als Satz; mit der "
       "S-Invarianz aus A* folgen die uhrverdrillten Sektoren mit festen Vielfachheiten", True)

# =================================================================== Teil C
say("\nTeil C: Assoziator, Hexagone und Anomalie")
M5, M3 = 5, 3


def carry(b, c):
    return (b + c) // 4


def om_exp(a, b, c):
    """Exponent von omega in Einheiten 1/8: omega = (-1)^(5 a_x [b_x + c_x >= 4] + 3 a_y [b_y + c_y >= 4])."""
    return (4 * (M5 * a[0] * carry(b[0], c[0]) + M3 * a[1] * carry(b[1], c[1]))) % 8


def c_exp(a, b):
    return (M5 * a[0] * b[0] + M3 * a[1] * b[1]) % 8


def add(a, b):
    return ((a[0] + b[0]) % 4, (a[1] + b[1]) % 4)


pent = all((om_exp(b, c, d) + om_exp(a, add(b, c), d) + om_exp(a, b, c) - om_exp(add(a, b), c, d) - om_exp(a, b, add(c, d))) % 8 == 0
           for a in A for b in A for c in A for d in A)
hex1 = all((om_exp(b, c, a) + c_exp(a, add(b, c)) + om_exp(a, b, c) - c_exp(a, c) - om_exp(b, a, c) - c_exp(a, b)) % 8 == 0
           for a in A for b in A for c in A)
hex2 = all((-om_exp(c, a, b) + c_exp(add(a, b), c) - om_exp(a, b, c) - c_exp(a, c) + om_exp(a, c, b) - c_exp(b, c)) % 8 == 0
           for a in A for b in A for c in A)
qform = all(Fr(c_exp(a, a), 8) % 1 == cs.q_form(a) for a in A)
nichttriv = sorted({om_exp(a, b, c) for a in A for b in A for c in A})


def anomalie(a, om=om_exp, addf=add, order=4):
    """Invariante der zyklischen Gruppe <a>: prod_j omega(a, j a, a) = exp(2 pi i k/order); Rückgabe k."""
    e, cur = 0, (0, 0) if isinstance(a[0], int) else None
    cur = tuple(0 for _ in a) if isinstance(a[0], int) else tuple((0, 0) for _ in a)
    for _ in range(order):
        e += om(a, cur, a)
        cur = addf(cur, a)
    return Fr(e % 8, 8) * order % order


hebungen = [add((1, 0), (k, k)) for k in range(4)]
anom = {str(h): int(anomalie(h)) for h in hebungen}
q_heb = {str(h): str(cs.q_form(h)) for h in hebungen}


def om_n(nl, nr):
    """Kozykel auf nl holomorphen und nr antiholomorphen Kopien, eingeschränkt auf die diagonale Uhr."""
    def om(a, b, c):
        return ((nl - nr) * om_exp(a, b, c)) % 8
    return om


add_tab = {(nl, nr): int(anomalie((1, 0), om=om_n(nl, nr))) for nl in range(0, 7) for nr in range(0, 7)}
add_ok = all(k == (2 * (nl - nr)) % 4 for (nl, nr), k in add_tab.items())
ZAHLEN["C1"] = {"Anomalie der vier Hebungen (in Z4)": anom, "q der Hebungen": q_heb,
                "Anomalie(N_L, N_R) = 2 (N_L - N_R) mod 4": add_ok}
check("C1 [E] Expliziter abelscher 3-Kozykel der 16 Uhrsektoren: c(a, b) = e^(2 pi i (5 a_x b_x + 3 a_y b_y)/8), "
      "omega(a, b, c) = (-1)^(5 a_x u_x + 3 a_y u_y) mit dem Übertrag u von b + c. Pentagon (65536 Gleichungen) und beide "
      "Hexagone (je 4096) gelten exakt, c(a, a) = e^(2 pi i q(a)), und omega ist nicht trivial. Die Anomalie der "
      "einseitigen Uhr ist die Invariante prod_j omega(a, j a, a) = e^(2 pi i k/4): k = 2 für alle vier Hebungen "
      "(1, 0) + (k, k) der Uhr (q = 5/8, 7/8, 1/8, 3/8). Sie ist additiv: Für N_L holomorphe und N_R antiholomorphe Nähte mit "
      "diagonaler Uhr ist k = 2 (N_L - N_R) mod 4 (geprüft für 0..6). Die Summe über die einseitige Uhr existiert also genau "
      "für N_L = N_R mod 2: mit Spiegelpartner immer, über N Nähte ohne Partner genau für gerades N. Das ist der Beweis der "
      "Geradzahligkeit für jedes N, den der Modultest (z174, z177) für N = 1..12 numerisch zeigte",
      pent and hex1 and hex2 and qform and nichttriv == [0, 4] and set(anom.values()) == {2} and add_ok)


# kommutative Algebra auf einer Untergruppe H von A (+) A-quer: kappa mit d kappa = omega_LR, kappa(a,b)/kappa(b,a) = c_LR
def lr_om(g, h, k):
    return (om_exp(g[0], h[0], k[0]) - om_exp(g[1], h[1], k[1])) % 8


def lr_c(g, h):
    return (c_exp(g[0], h[0]) - c_exp(g[1], h[1])) % 8


def smith_loesbar(M, r, K):
    """Lösbarkeit von M x = r über Z/2^K, exakt über die Smith-Form (Zeilen- und Spaltenoperationen).
    Rückgabe (lösbar, größte Pivotbewertung)."""
    mod = 1 << K
    M = np.array(M, dtype=np.int64) % mod
    r = np.array(r, dtype=np.int64) % mod
    zl = np.ones(M.shape[0], bool)
    sl = np.ones(M.shape[1], bool)
    vmax = 0
    while True:
        zi, si = np.flatnonzero(zl), np.flatnonzero(sl)
        sub = M[np.ix_(zi, si)]
        if not sub.any():
            break
        low = sub & -sub
        low[low == 0] = mod
        a, b = np.unravel_index(np.argmin(low), low.shape)
        v = int(low[a, b]).bit_length() - 1
        i, j = zi[a], si[b]
        uinv = pow(int(M[i, j]) >> v, -1, mod)
        M[i] = (M[i] * uinv) % mod
        r[i] = (r[i] * uinv) % mod
        f = M[:, j] >> v
        f[i] = 0
        f[~zl] = 0
        M = (M - np.outer(f, M[i])) % mod
        r = (r - f * r[i]) % mod
        M[i, :] = 0
        M[i, j] = 1 << v
        if r[i] % (1 << v):
            return False, max(vmax, v)
        vmax = max(vmax, v)
        zl[i] = False
        sl[j] = False
    return bool(not r[zl].any()), vmax


def algebra(H, K=12, c_vorzeichen=1):
    """Gibt es kappa: H x H -> U(1) mit d kappa = omega (assoziativ) und kappa(a, b)/kappa(b, a) = c(a, b) (kommutativ)?"""
    Hs = sorted(H)
    n = len(Hs)
    pos = {g: i for i, g in enumerate(Hs)}
    f = (1 << K) // 8
    M = np.zeros((n ** 3 + n * n, n * n), dtype=np.int64)
    r = np.zeros(n ** 3 + n * n, dtype=np.int64)
    z = 0
    for g in Hs:
        for h in Hs:
            for k in Hs:
                M[z, pos[g] * n + pos[h]] += 1
                M[z, pos[cs.g_add(g, h)] * n + pos[k]] += 1
                M[z, pos[h] * n + pos[k]] -= 1
                M[z, pos[g] * n + pos[cs.g_add(h, k)]] -= 1
                r[z] = f * lr_om(g, h, k)
                z += 1
    for g in Hs:
        for h in Hs:
            M[z, pos[g] * n + pos[h]] += 1
            M[z, pos[h] * n + pos[g]] -= 1
            r[z] = c_vorzeichen * f * lr_c(g, h)
            z += 1
    return smith_loesbar(M, r, K)


alg, vmax = {}, 0
N = neutral(1, 1)
for name, H in [("N", N)] + [(f"über N: {typ(L)}, unverdrillt {len(L & e8e8(1, 1))}", L) for L in COMP if N <= L] + \
        [("einseitige Uhr <(1,0)> (q = 5/8)", frozenset(mult(((1, 0), (0, 0)), k) for k in range(4))),
         ("<(2,0)> (Vektor, q = 1/2)", frozenset(mult(((2, 0), (0, 0)), k) for k in range(2))),
         ("E8-Kleber einseitig", frozenset((a, (0, 0)) for a in glue(1)))]:
    ok, v = algebra(H)
    ok2, _ = algebra(H, c_vorzeichen=-1)
    alg[name] = [ok, ok2, all(cs.Q_lr(g) == 0 for g in H)]
    vmax = max(vmax, v)
for k_, v_ in alg.items():
    say(f"      {k_}: kommutative Algebra {'ja' if v_[0] else 'nein'}; q = 0 auf H: {'ja' if v_[2] else 'nein'}")
ZAHLEN["C2"] = {"Algebra [lösbar, lösbar mit c^-1, q = 0]": alg, "größte Pivotbewertung (Modul 2^12)": vmax}
check("C2 [E] Assoziativität: Eine Erweiterung um eine Untergruppe H trägt eine assoziative und kommutative Multiplikation "
      "genau dann, wenn es kappa mit d kappa = omega und kappa(a, b)/kappa(b, a) = c(a, b) gibt (lineares System, exakt "
      "über die Smith-Form modulo 2^12 entschieden). In allen neun Fällen ist es genau dann lösbar, wenn q = 0 auf H ist: "
      "lösbar für die neutrale Algebra N, für den E8-Kleber und für die drei bosonischen Vervollständigungen über N (volle "
      "Ordnung 16); nicht lösbar für die zwei Spin-Formen, für die einseitige Uhr <(1, 0)> (q = 5/8: Anomalie) und für "
      "<(2, 0)> (q = 1/2: fermionisch). Die Aufhebung der Anomalie und die Assoziativität der lokalen Fusion sind dieselbe "
      "Bedingung",
      all(v_[0] == v_[2] and v_[0] == v_[1] for v_ in alg.values()) and sum(v_[0] for v_ in alg.values()) == 5
      and len(alg) == 9 and vmax <= 8)
befund("Teil C: Die Anomalie der einseitigen Uhr ist die Assoziatorklasse 2 in Z4. Sie hebt sich genau mit Spiegelpartner "
       "oder bei gerader Nahtzahl auf, exakt und für jedes N", True)

# =================================================================== Teil D
say("\nTeil D: Die Eichsumme als unabhängige Spinstruktursumme der beiden Deckhälften")
mp.mp.dps = 40


def thetas(tau_):
    qq = mp.exp(mp.pi * 1j * tau_)
    e_ = mp.exp(2j * mp.pi * tau_ / 24) * mp.qp(mp.exp(2j * mp.pi * tau_))
    return mp.jtheta(2, 0, qq), mp.jtheta(3, 0, qq), mp.jtheta(4, 0, qq), e_


def z_summen(tau_):
    t2, t3, t4, e_ = thetas(tau_)
    gso = lambda nn: (abs(t2) ** (2 * nn) + abs(t3) ** (2 * nn) + abs(t4) ** (2 * nn)) / (2 * abs(e_) ** (2 * nn))
    chi10 = [(t3 ** 5 + t4 ** 5) / (2 * e_ ** 5), t2 ** 5 / (2 * e_ ** 5), (t3 ** 5 - t4 ** 5) / (2 * e_ ** 5), t2 ** 5 / (2 * e_ ** 5)]
    chi6 = [(t3 ** 3 + t4 ** 3) / (2 * e_ ** 3), t2 ** 3 / (2 * e_ ** 3), (t3 ** 3 - t4 ** 3) / (2 * e_ ** 3), t2 ** 3 / (2 * e_ ** 3)]
    diag = sum(abs(chi10[x] * chi6[y]) ** 2 for x in range(4) for y in range(4))
    e8 = abs((t2 ** 8 + t3 ** 8 + t4 ** 8) / (2 * e_ ** 8)) ** 2
    return {"getrennt je Hälfte": gso(5) * gso(3), "gemeinsam": gso(8), "chiral je Seite (E8 x E8-quer)": e8, "diagonal": diag}


pts = [mp.mpc('0.137', '1.093'), mp.mpc('-0.211', '0.947'), mp.mpc('0.389', '1.211')]
ident = max(abs(z_summen(t)["getrennt je Hälfte"] - z_summen(t)["diagonal"]) / z_summen(t)["diagonal"] for t in pts)
inv_max = mp.mpf(0)
for t in pts:
    z0 = z_summen(t)
    for t2_ in (-1 / t, t + 1):
        z1 = z_summen(t2_)
        inv_max = max(inv_max, max(abs(z0[k] - z1[k]) / abs(z0[k]) for k in z0))
verschieden = len({mp.nstr(v, 12) for v in z_summen(pts[0]).values()}) == 3
ZAHLEN["D1"] = {"Identität (rel. Fehler)": mp.nstr(ident, 3), "Modulinvarianz (rel. Fehler)": mp.nstr(inv_max, 3),
                "Werte bei tau = i": {k: mp.nstr(v, 15) for k, v in z_summen(mp.mpc(0, 1)).items()}}
check("D1 [E/N, 40 Stellen] Die 16 reellen Eltern mit Spiegelpartner als freie Fermionen: Die Summe über die "
      "Spinstrukturen, getrennt für die zwei Deckhälften (10 und 6 Eltern, links und rechts gemeinsam), ist "
      "[(|th2|^10 + |th3|^10 + |th4|^10)/(2 |eta|^10)] [(|th2|^6 + |th3|^6 + |th4|^6)/(2 |eta|^6)] und stimmt mit der "
      f"diagonalen Vervollständigung Summe_a |chi_a|^2 der 16 Uhrsektoren überein (rel. Fehler {mp.nstr(ident, 2)}). "
      "Die gemeinsame Summe über alle 16 Eltern gibt so(16)_1 diagonal, die chirale Summe je Seite E8 (x) E8-quer; alle "
      f"drei sind modulinvariant (rel. Fehler {mp.nstr(inv_max, 2)}) und voneinander verschieden. Die Eichsumme der "
      "Paaruhr ist also dieselbe Operation wie in H6 der Schritt 'A* gibt die GSO-Summe', angewandt auf die Deckhälften "
      "einzeln: Die Paaruhr wirkt auf die Eltern beider Seiten als Deckinvolution (+1 auf 10, -1 auf 6)",
      ident < mp.mpf(10) ** -30 and inv_max < mp.mpf(10) ** -30 and verschieden)
befund("Teil D: Eichung der Paaruhr = Unabhängigkeit von der Spinstruktur je Deckhälfte. A* allein wählt zwischen den drei "
       "Summen nicht; die Auswahl trifft die Superauswahl der Paaruhr (Teil A, B)", True)

# =================================================================== Teil E
say("\nTeil E: Nahtzahl im Modultest für N bis 24 (Methode z177)")
SRC177 = (HERE / "z177_nahtpaar_kandidaten_nullkegel.py").read_text()
kopf = SRC177.split("GW = {k: gitter(B[k]) for k in B}")[0]
ns = {"__file__": str(HERE / "z177_nahtpaar_kandidaten_nullkegel.py"), "__name__": "z177_teil"}
with contextlib.redirect_stdout(io.StringIO()):
    exec(kopf, ns)
gitter, moden, comb, konj, EPS = ns["gitter"], ns["moden"], ns["comb"], ns["konj"], ns["EPS"]
Bn, YS = ns["B"], ns["YS"]
NS2 = tuple(range(2, 25, 2))


def tafel_mp(GWx):
    """Wie tafel in z177, aber mit mpmath-Zahlen (die normierten Minima überschreiten für N >= 18 den Gleitkommabereich)."""
    out = {}
    for N in NS2:
        M = {k: moden(Bn[k], GWx[k], N) for k in Bn}
        m0 = M["1"]
        z = {}
        for j in range(4):
            m_ = comb([m0, M["G"], M["X"], konj(M["X"])], [1, (-1) ** j, (-1j) ** j, 1j ** j])
            z[j] = min(m_[Y][s_].real / abs(m0[Y][0]) for Y in YS for s_ in m_[Y])
        out[N] = z
    return out


GW80 = {k: gitter(Bn[k], 80) for k in Bn}
GW60 = {k: gitter(Bn[k], 60) for k in Bn}
mp.mp.dps = 80
T80, T60 = tafel_mp(GW80), tafel_mp(GW60)
relabw = max(abs(T80[N][j] - T60[N][j]) / max(1, abs(T80[N][j])) for N in NS2 for j in range(4))
pos = {N: [j for j in range(4) if T80[N][j] > EPS] for N in NS2}
stabil = all((T80[N][j] > EPS) == (T60[N][j] > EPS) for N in NS2 for j in range(4))
normal = [N for N in NS2 if 2 in pos[N]]
vertauscht = {N: [j for j in pos[N] if j % 2] for N in NS2 if any(j % 2 for j in pos[N])}
for N in NS2:
    say(f"      N = {N}: " + "; ".join(f"X = i^{j}: {mp.nstr(T80[N][j], 3)}" for j in range(4)))
mp.mp.dps = 40
ZAHLEN["E1"] = {"positive Eigenwerte i^j je N": {str(N): v for N, v in pos.items()}, "normale Statistik": normal,
                "vertauschte Statistik": {str(k): v for k, v in vertauscht.items()},
                "Vorzeichen stabil (60 gegen 80 Stellen)": stabil, "rel. Abweichung 60 gegen 80 Stellen": mp.nstr(relabw, 3),
                "Tafel": {str(N): {f"i^{j}": mp.nstr(T80[N][j], 6) for j in range(4)} for N in NS2}}
check(f"E1 [N, 80 Stellen] Modultest der Z4-Verklebung für gerade N = 2..24 Paare: normale Statistik (Eigenwert -1) besteht "
      f"nur für N = {normal}; vertauschte Statistik besteht für {vertauscht}; für N = 8..24 besteht kein Eigenwert, und "
      "der Sektor des Skalarvakuums besteht nie. Gitterwerte bei 60 und 80 Stellen geben dieselben Vorzeichen und "
      f"dieselben Minima (relative Abweichung {mp.nstr(relabw, 2)}). Der Testbereich aus z177 (N <= 12) ist damit "
      "verdoppelt; ein Beweis für jedes N ist das nicht",
      normal == [2] and stabil and relabw < mp.mpf(10) ** -20 and all(pos[N] == [] for N in NS2 if N >= 8)
      and all(0 not in pos[N] for N in NS2))
befund("Teil E: Geradzahligkeit ist für jedes N bewiesen (C1). 'Genau zwei Paare' bleibt eine numerische Aussage, jetzt "
       "für N <= 24", True)
say(f"   {zeit()}")

# =================================================================== Teil F
say("\nTeil F: Bilanz")
kette = [
    "[Annahme] P_std: die Paaruhr ist superselektiert (jede Observable vertauscht mit U)",
    "[Bestand] GNS-Regel: der Zustandsraum ist der GNS-Raum des Zustands von A* über den Observablen",
    "[Satz, B1] GNS-Raum = H^U, die Vakuumdarstellung der uhrneutralen Algebra N; Gauß-Zwang U = 1 auf allen Zuständen",
    "[A*] Basiswechsel sind Symmetrien: die Zustandssumme ist S- und T-invariant; Sp über H^U allein ist es nicht (A5)",
    "[Satz, A1 bis A3] Über N gibt es genau eine S- und T-invariante Vervollständigung ohne Feld mit Paaruhr-Ladung: jeder "
    "Fluss einmal. Das ist die diagonale Vervollständigung mit allen 16 Uhrsektoren. Ohne T bliebe eine Spin-Form",
    "[Folge] Der Zustand von A* enthält die uhrverdrillten Sektoren mit dem Gewicht (12 sqrt(2) - 13)/14 (K19)",
    "[Folge] Der einseitige Modulator rho des Zustands gibt T, mit der jetzt hergeleiteten Fusion S, die Uhr als "
    "Monodromie, den Defekt und den Blatttausch S^2 = C (z176 C4); der Partner ist die Purifikation (K20)",
]
for k in kette:
    say("   " + k)
check("F1 [E] Kette: Aus P_std, der GNS-Regel und A* folgt die diagonale Vervollständigung mit allen 16 Uhrsektoren als "
      "Satz (A1 bis A3, B1, B2). Die Lesart 'Superauswahl = Eichsumme' aus z176 ist damit durch eine Herleitung ersetzt; "
      "eine Auswahl unter den 27 findet nicht statt. Eingesteckt bleibt P_std",
      all(CHECKS.values()))
befund("K19 folgt aus P_std, der GNS-Regel und A* (kein eigener offener Punkt mehr, solange P_std gezählt wird)", True)
befund("K20 folgt unter denselben Voraussetzungen und der Identifikation der Uhr (P_zentral): Elternobjekt ist der "
       "einseitige Modulator rho des Zustands von A*; die Teilung Naht | Partner ist geometrisch (z180)", True)
befund("P_std ist hergeleitet", False)
befund("P_paar ist hergeleitet (genau zwei Paare für jedes N)", False)
befund("Komplettlösung", False)

urteil = (
    "F1 ist gerechnet. (1) Über der uhrneutralen Algebra N = (E8 (x) E8-quer)^U liegen genau fünf der 27 "
    "Vervollständigungen. Drei sind bosonisch: E8 (x) E8-quer (alle Ladungen der Paaruhr, kein Uhrfluss), eine Zwischenform "
    "(nur die Fermionparität geeicht) und die Eichsumme (alle Flüsse, keine Ladung); zwei sind Spin-Formen. Die Sektoren "
    "von N sind das Doppel D(Z4); die Eichsumme ist das einzige bosonische Komplement der Ladungen. (2) Die Lücke zwischen Superauswahl und Eichung ist ein Satz der "
    "Standardform: Der GNS-Raum des Zustands über den Paaruhr-neutralen Observablen ist genau der U-invariante Teilraum, "
    "und die Projektion E(X) = (1/4) Summe U^k X U^-k ist die Jones-Projektion. Das Gegenmodell 'E8 (x) E8-quer mit "
    "Auswahlregel' scheitert an der GNS-Regel des Bestands, nicht an A*. (3) Die S- und T-Invarianz aus A* legt dann die "
    "uhrverdrillten Sektoren mit Vielfachheit fest: Der invariante Raum ist dreidimensional, und 'Vakuum einmal, keine "
    "Ladung' hat genau eine Lösung. Ohne T bleibt als zweites Modell eine Spin-Form (Fluss 1 mit Ladung 2 und Spin 1/2); "
    "sie ist als Gegenmodell verzeichnet. K19 und K20 folgen damit aus P_std, der GNS-Regel und A* ohne Lesart. (4) Die Anomalie "
    "der einseitigen Uhr ist die Assoziatorklasse 2 in Z4 eines expliziten Kozykels (Pentagon und Hexagone exakt); sie ist "
    "additiv, also ist die Nahtzahl ohne Partner für jedes N gerade. Assoziativität der Fusion und Aufhebung der Anomalie "
    "sind dieselbe Bedingung. (5) Die Eichsumme ist die Spinstruktursumme der beiden Deckhälften getrennt, dieselbe "
    "Operation wie 'A* gibt GSO' in H6. (6) 'Genau zwei Paare' besteht numerisch jetzt für N <= 24; bewiesen ist nur die "
    "Geradzahligkeit. Offen bleiben P_std selbst und der Rest von P_paar. Keine Komplettlösung, nichts befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z178_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil},
                                                    indent=1, ensure_ascii=False, default=str))
(HERE / "ausgabe_z178.txt").write_text("\n".join(x for x in OUT if "Laufzeit bisher" not in x) + "\n")
say(f"{zeit()}")
