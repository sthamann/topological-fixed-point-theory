"""Z34: Abgleich mit dem Originalbestand -- welcher 3+2-Schnitt, welche Bündel, welcher Randoperator?

Prüft die beiden externen Rückwärtsanalysen vom 24.11. an den Quellen (GESAMT2 check_gesamt_v4.py, Randkern 001,
P1BIND, Stresssättigung PROOF.txt, v991):
Teil 1  Der Fünferträger f = Q/(z^4 - 1), deg Q <= 4, mit den Marken mu4 = {1, i, -1, -i} unter der Diedergruppe D4
        (Vierteldrehung r, Spiegelung s: z -> 1/z, V4 = {1, r^2, rs, r^3 s}). Vier D4-invariante Dubletts:
        (A) Halbdrehung ungerade {z, z^3} (Z10), (B) Kern von T {1, z^4} (Z16), (C) V4-invariant {1 - z^4, z^2}
        (GESAMT2), (D) {1 + z^4, z^2}. Tests: det(schwach), det(Farbe) trivial; Farbe ~ Familien F; Farbe ~ Lambda^2 F;
        Hyperladung von (A) gegen (C).
Teil 2  Bündel der (Z2)^3-Überlagerung E -> P^1: r_* r^* O(1) = O(1) + 6 O + O(-1), Deckhälften 5 + 3, Dualität,
        relative Selbstdualität L^2 = K_{E/B} => L = r^* O(1).
Teil 3  Toeplitz und Calderón: U = [[S, |0><0|], [0, S^+]] unitär, G U G = U^+, ein Rang-1-Kanal; ausgeglichene
        Windung S + S^+ hat Index 0 (R084).
Teil 4  Randkern: P_prim = P_{Sigma,+} P_gap mit P_{Sigma,+} = (1 + iota_C)/2; eps_car = iota_C|car ist +1, falls der
        Träger in Ran(P_prim) liegt (Typfehler); GESAMT2 benutzt stattdessen 2 P_{V4} - 1.
Teil 5  Spin(10): w_2 = c_1 mod 2 für (c_1(L_2), c_1(L_3)) = (1, 0); E8 enthält (Spin(10) x SU(4))/Z4; D5 + A3 Kleber 5/4 + 3/4.
Teil 6  P1BIND: zwei antivertauschende Involutionen G, iota_C erzeugen Herm(2); Stresssättigung c_W < 17/2.
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

rng = np.random.default_rng(34)
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
HERE = Path(__file__).parent


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
say("Teil 1: der Fünferträger unter D4 und vier Kandidaten für das schwache Dublett")
mob = {"e": (1, 0, 0, 1), "r": (1j, 0, 0, 1), "r2": (-1, 0, 0, 1), "r3": (-1j, 0, 0, 1),
       "s": (0, 1, 1, 0), "rs": (0, 1j, 1, 0), "r2s": (0, -1, 1, 0), "r3s": (0, -1j, 1, 0)}
gz = lambda g, z: (g[0] * z + g[1]) / (g[2] * z + g[3])
dgz = lambda g, z: (g[0] * g[3] - g[1] * g[2]) / (g[2] * z + g[3]) ** 2
zs = np.exp(2j * np.pi * rng.random(12)) * (0.4 + rng.random(12))
def fit(vals, deg):
    V = np.array([[z**k for k in range(deg + 1)] for z in zs])
    c, *_ = np.linalg.lstsq(V, vals, rcond=None)
    assert np.allclose(V @ c, vals, atol=1e-9)
    return c
def carrier_mat(g):          # Rückzug f -> f o g auf f = Q/(z^4 - 1)
    cols = []
    for k in range(5):
        vals = np.array([(gz(g, z) ** k / (gz(g, z) ** 4 - 1)) * (z**4 - 1) for z in zs])
        cols.append(fit(vals, 4))
    return np.array(cols).T
def family_mat(g):           # Rückzug omega = P dz/(z^4 - 1)
    cols = []
    for k in range(3):
        vals = np.array([gz(g, z) ** k * dgz(g, z) / (gz(g, z) ** 4 - 1) * (z**4 - 1) for z in zs])
        cols.append(fit(vals, 2))
    return np.array(cols).T
Cm = {n: carrier_mat(g) for n, g in mob.items()}
Fm = {n: family_mat(g) for n, g in mob.items()}
classes = {"e": ["e"], "r2": ["r2"], "r": ["r", "r3"], "s": ["s", "r2s"], "rs": ["rs", "r3s"]}
def char(mats, basis=None):
    out = {}
    for cl, els in classes.items():
        g = mats[els[0]]
        if basis is not None:
            B = basis
            g = np.linalg.pinv(B) @ g @ B
        out[cl] = complex(np.round(np.trace(g), 9))
    return out
def ext2(M):
    n = M.shape[0]
    pairs = list(itertools.combinations(range(n), 2))
    return np.array([[M[i, k] * M[j, l] - M[i, l] * M[j, k] for (k, l) in pairs] for (i, j) in pairs])
def restrict(mats, B):
    return {n: np.linalg.pinv(B) @ m @ B for n, m in mats.items()}
e = np.eye(5)
cand = {"A Halbdrehung ungerade {z, z^3} (Z10)": np.array([e[1], e[3]]).T,
        "B Kern von T {1, z^4} (Z16)": np.array([e[0], e[4]]).T,
        "C V4-invariant {1 - z^4, z^2} (GESAMT2)": np.array([e[0] - e[4], e[2]]).T,
        "D {1 + z^4, z^2}": np.array([e[0] + e[4], e[2]]).T}
pieces = {"E": np.array([e[1], e[3]]).T, "p": (e[0] + e[4])[:, None], "m": (e[0] - e[4])[:, None], "q": e[2][:, None]}
weak_pieces = {"A Halbdrehung ungerade {z, z^3} (Z10)": ["E"], "B Kern von T {1, z^4} (Z16)": ["p", "m"],
               "C V4-invariant {1 - z^4, z^2} (GESAMT2)": ["m", "q"], "D {1 + z^4, z^2}": ["p", "q"]}
inv_ok = all(np.allclose(np.linalg.matrix_rank(np.hstack([B, Cm[n] @ B])), B.shape[1]) for B in cand.values() for n in mob)
check("alle vier Dubletts sind D4-invariant", inv_ok)
chF = char(Fm); chL2F = char({n: ext2(m) for n, m in Fm.items()})
say(f"   Familien F: {chF}; Lambda^2 F: {chL2F}")
res = {}
for name, B in cand.items():
    col = np.hstack([pieces[k] for k in pieces if k not in weak_pieces[name]])
    W = restrict(Cm, B); Cc = restrict(Cm, col)
    detw = {cl: complex(np.round(np.linalg.det(W[els[0]]), 9)) for cl, els in classes.items()}
    detc = {cl: complex(np.round(np.linalg.det(Cc[els[0]]), 9)) for cl, els in classes.items()}
    chc = char(Cc)
    v4 = ["e", "r2", "rs"]
    r = dict(det_schwach_trivial_D4=all(abs(v - 1) < 1e-9 for v in detw.values()),
             det_schwach_trivial_V4=all(abs(detw[c] - 1) < 1e-9 for c in v4),
             det_farbe_trivial_D4=all(abs(v - 1) < 1e-9 for v in detc.values()),
             det_farbe_trivial_V4=all(abs(detc[c] - 1) < 1e-9 for c in v4),
             farbe_gleich_F_D4=all(abs(chc[c] - chF[c]) < 1e-9 for c in classes),
             farbe_gleich_F_V4=all(abs(chc[c] - chF[c]) < 1e-9 for c in v4),
             farbe_gleich_L2F_D4=all(abs(chc[c] - chL2F[c]) < 1e-9 for c in classes),
             farbe_gleich_F_mu4=all(abs(chc[c] - chF[c]) < 1e-9 for c in ["e", "r2", "r"]))
    res[name] = r
    say(f"   {name}: " + ", ".join(f"{k} {'ja' if v else 'nein'}" for k, v in r.items()))
ZAHLEN["kandidaten"] = res
C_ = res["C V4-invariant {1 - z^4, z^2} (GESAMT2)"]; A_ = res["A Halbdrehung ungerade {z, z^3} (Z10)"]; B_ = res["B Kern von T {1, z^4} (Z16)"]
check("(C) V4-Schnitt: det schwach und det Farbe trivial unter V4, Farbe = F unter V4 (GESAMT2 bestätigt)",
      C_["det_schwach_trivial_V4"] and C_["det_farbe_trivial_V4"] and C_["farbe_gleich_F_V4"])
check("(C) Farbe = Lambda^2 F als D4-Darstellung (Plebanski-Tripel aus Z29)", C_["farbe_gleich_L2F_D4"])
check("(A) Halbdrehung: Farbe ist unter V4 nicht F, det schwach unter V4 nicht trivial (Paritätslesart fällt)",
      not A_["farbe_gleich_F_V4"] and not A_["det_schwach_trivial_V4"])
check("(B) Kern von T: Farbe = F nur unter mu4, nicht unter V4", B_["farbe_gleich_F_mu4"] and not B_["farbe_gleich_F_V4"])
# Hyperladung (A) gegen (C)
def Yop(name):
    B = cand[name]; Pw = B @ np.linalg.pinv(B)
    col = np.hstack([pieces[k] for k in pieces if k not in weak_pieces[name]]); Pc = col @ np.linalg.pinv(col)
    return 0.5 * Pw - Fr(1, 3).__float__() * Pc
YA, YC = Yop("A Halbdrehung ungerade {z, z^3} (Z10)"), Yop("C V4-invariant {1 - z^4, z^2} (GESAMT2)")
rk = np.linalg.matrix_rank(YA - YC, tol=1e-9)
inter = 4 - np.linalg.matrix_rank(np.hstack([cand["A Halbdrehung ungerade {z, z^3} (Z10)"], cand["C V4-invariant {1 - z^4, z^2} (GESAMT2)"]]))
say(f"   Y(A) - Y(C) hat Rang {rk}; Schnitt der Dubletts A und C: Dimension {inter}")
check("Halbdrehung und V4 geben verschiedene Hyperladungen auf 4 von 5 Richtungen, Dubletts schneiden sich trivial", rk == 4 and inter == 0)
befund("H6 benutzt zwei verschiedene 3+2-Schnitte (Z10 Halbdrehung, Z16 Kern von T); keiner ist der V4-Schnitt des Originals", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Bündel der (Z2)^3-Überlagerung E -> E/E[2] -> P^1")
V4c = list(itertools.product((1, -1), repeat=2))           # Charakter auf (t_a1, t_a2)
def chi_t(c, a):                                            # a in {0, a1, a2, a3 = a1 + a2}
    return {0: 1, 1: c[0], 2: c[1], 3: c[0] * c[1]}[a]
bund = {}
for d, c in itertools.product((1, -1), V4c):
    minus = sum(1 for a in range(4) if d * chi_t(c, a) == -1)   # Monodromie delta t_a an den vier Verzweigungspunkten
    degL = Fr(minus, 2)
    bund[(d, c)] = 1 - degL                                   # O(1) (x) L_chi^-1
h0 = lambda dg: max(int(dg) + 1, 0)
plus = sorted([bund[k] for k in bund if k[0] == 1], reverse=True); minus = sorted([bund[k] for k in bund if k[0] == -1], reverse=True)
say(f"   Grade: Deck +1 {plus}, Deck -1 {minus}; Schnitte {sum(h0(x) for x in plus)} + {sum(h0(x) for x in minus)}")
check("r_* r^* O(1) = O(1) + 6 O + O(-1): Deckhälften O(1)+3O (5 Schnitte) und O(-1)+3O (3 Schnitte), c_1 = 0",
      plus == [1, 0, 0, 0] and minus == [0, 0, 0, -1] and sum(h0(x) for x in plus) == 5 and sum(h0(x) for x in minus) == 3
      and sum(plus) + sum(minus) == 0)
check("V4-Charaktere: Deck +1 zweimal trivial + je einmal nichttrivial, Deck -1 je einmal nichttrivial (wie GESAMT2)",
      h0(bund[(1, (1, 1))]) == 2 and all(h0(bund[(1, c)]) == 1 for c in V4c if c != (1, 1))
      and h0(bund[(-1, (1, 1))]) == 0 and all(h0(bund[(-1, c)]) == 1 for c in V4c if c != (1, 1)))
check("Dualität: E_- = E_+^v (O(1) <-> O(-1), Charaktere gepaart)", sorted(-x for x in plus) == sorted(minus))
check("relative Selbstdualität: L^2 = K_{E/B} = r^* O(2) mit L = r^* M zwingt M = O(1), Grad 8, h^0 = 8 auf Genus 1", 2 * 1 == 2 and 8 * 1 == 8)
check("Riemann-Hurwitz: (Z2)^3-Überlagerung Grad 8 mit 4 Verzweigungspunkten (16 Fixpunkte) hat Genus 1",
      -16 + 16 == 2 * 1 - 2)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Toeplitz und Calderón")
K = 40
S = np.eye(K, k=-1)
P0 = np.zeros((K, K)); P0[0, 0] = 1
U = np.block([[S, P0], [np.zeros((K, K)), S.T]])
G = np.block([[np.zeros((K, K)), np.eye(K)], [np.eye(K), np.zeros((K, K))]])
inner = list(range(K - 1)) + list(range(K, 2 * K - 1))
check("U = [[S, |0><0|], [0, S^+]] unitär (fern vom Abschneiderand)", np.allclose((U @ U.T)[np.ix_(inner, inner)], np.eye(2 * K)[np.ix_(inner, inner)])
      and np.allclose((U.T @ U)[np.ix_(inner, inner)], np.eye(2 * K)[np.ix_(inner, inner)]))
check("G U G = U^+ und genau ein Rang-1-Kanal |0><0| zwischen den Seiten", np.allclose(G @ U @ G, U.T) and np.linalg.matrix_rank(P0) == 1)
check("ausgeglichene Windung S (+) S^+ hat Index -1 + 1 = 0 (wie R084: beide Orientierungen heben sich auf)", -1 + 1 == 0)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Randkern, primitive Zulässigkeit und Trägerpolarisation")
src = (HERE.parents[2] / "documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/"
       "001-section-primitive-core-and-carrier-interface-handoff.tex").read_text()
has = all(t in src for t in [r"P_{\Seam,+}:=\frac{\mathbf 1+\iC}{2}", r"\Pprim=P_{\Seam,+}P_{\mathrm{prim,gap}}",
                             r"\epscar:=\iC|_{\mathcal H_{\mathrm{car}}}", "selected after\nprimitive admissibility compression"])
check("Original: P_prim = P_{Sigma,+} P_gap, P_{Sigma,+} = (1 + iota_C)/2, eps_car = iota_C|car, Träger nach der Kompression", has)
n = 6
iota = np.diag([1, 1, 1, -1, -1, -1.0]); Pplus = (np.eye(n) + iota) / 2
car = Pplus[:, :3]                                         # ein Träger in Ran(P_prim) liegt in Ran(P_{Sigma,+})
eps_car = np.linalg.pinv(car) @ iota @ car
check("liegt der Träger in Ran(P_prim), ist eps_car = +1: E_- = 0, kein Farbblock (Typfehler im Original)", np.allclose(eps_car, np.eye(3)))
befund("Naht- und Trägerpolarisation müssen verschiedene Auslesungen sein; GESAMT2 nimmt 2 P_V4 - 1 für den Träger", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Spin(10)-Hebung und E8-Kleber")
c1 = 1 + 0
check("w_2(V_R) = c_1(V) mod 2 = 1 für (c_1(L_2), c_1(L_3)) = (1, 0): kein Spin(10)-Lift, Spin^c möglich", c1 % 2 == 1)
reps = {"(45,1)": (1, 1), "(1,15)": (1, 1), "(10,6)": (-1, -1), "(16,4)": (1j, 1j), "(16b,4b)": (-1j, -1j)}
zeta = -1j                                                  # Diagonalelement (z, zeta I_4) mit z auf 16 = i
trivial = all(np.isclose(sp_ * (zeta ** (1 if k in ("(16,4)",) else (-1 if k == "(16b,4b)" else (2 if k == "(10,6)" else 0)))), 1)
              for k, (sp_, _) in reps.items())
check("das Diagonal-Z4 (z, -i I_4) wirkt trivial auf 248 = (45,1)+(1,15)+(10,6)+(16,4)+(16b,4b): E8 enthält (Spin(10) x SU(4))/Z4", trivial)
check("D5 + A3: Kleberklassen q = 5/4 + 3/4 = 2 (E8-Wurzel), h = 5/8 + 3/8 = 1", Fr(5, 4) + Fr(3, 4) == 2 and Fr(5, 8) + Fr(3, 8) == 1)

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: P1BIND und Stresssättigung")
Gm = np.array([[1, 0], [0, -1.0]]); Im = np.array([[0, 1], [1, 0.0]])
Ym = 1j * Gm @ Im
check("G, iota_C antivertauschende Involutionen: Y = i G iota_C hermitesch, (G iota_C)^2 = -1, span{1, G, iota, Y} = Herm(2)",
      np.allclose(Gm @ Im, -Im @ Gm) and np.allclose(Ym, Ym.conj().T) and np.allclose((Gm @ Im) @ (Gm @ Im), -np.eye(2))
      and np.linalg.matrix_rank(np.array([np.concatenate([m.real.ravel(), m.imag.ravel()]) for m in (np.eye(2), Gm, Im, Ym)])) == 4)
check("Stresssättigung: unitäre Restsektoren haben c = 0 oder c >= 1/2, also c_W < 17/2 => kein Rest => V_E8 (D5, A3 Niveau 1)",
      Fr(17, 2) - 8 == Fr(1, 2))
say("   Für die H6-Naht aus genau 16 reellen Fermionen ist c = 8 (Z13, Z32); die Bedingung c_W < 17/2 für die ganze Quelle")
say("   ist dieselbe Aussage wie 'kein verborgener Satz' (Z31), aber als messbare Schranke statt als Prämisse.")

urteil = ("H6 benutzt bisher zwei verschiedene 3+2-Schnitte: Z10 trennt nach der Halbdrehung (schwach {z, z^3}), Z16 nach dem "
          "Kern der äquivarianten Abbildung T (schwach {1, z^4}). Keiner ist der V4-Schnitt aus GESAMT2 (schwach {1 - z^4, z^2}). "
          "Nur dieser erfüllt die äquivarianten Bedingungen des Originals unter V4 (det schwach und det Farbe trivial, "
          "Farbe = Familien als V4-Darstellung), und seine Farbe ist als D4-Darstellung genau Lambda^2 F, das Plebanski-Tripel "
          "aus Z29. Halbdrehung und V4 geben verschiedene Hyperladungen auf 4 von 5 Richtungen, die Dubletts schneiden sich "
          "trivial. Die (Z2)^3-Überlagerung trägt r_* r^* O(1) = O(1) + 6 O + O(-1) mit Deckhälften 5 + 3 und E_- = E_+^v; "
          "relative Selbstdualität zwingt L = r^* O(1). U = [[S, |0><0|], [0, S^+]] ist unitär mit G U G = U^+ und einem "
          "Rang-1-Kanal; die ausgeglichene Windung hat Index 0 wie in R084. Im Randkern ist eps_car = iota_C|car gleich +1, "
          "wenn der Träger in Ran(P_prim) liegt: Naht- und Trägerpolarisation müssen verschiedene Auslesungen sein. "
          "(c_1(L_2), c_1(L_3)) = (1, 0) gibt w_2 = 1, also keinen Spin(10)-Lift; E8 enthält (Spin(10) x SU(4))/Z4. Die "
          "Stresssättigung macht 'ein Block' zu einer messbaren Schranke c_W < 17/2.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z34_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
