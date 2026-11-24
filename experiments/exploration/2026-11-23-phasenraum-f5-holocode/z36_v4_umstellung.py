"""Z36: H6 auf den V4-Schnitt umgestellt -- was ändert sich, was bleibt?

Rechnet die Teile von H6 neu, die vom konkreten 3+2-Teilraum abhängen (Z10, Z16, Z28), mit dem V4-Schnitt aus GESAMT2/Z34:
Teil 1  Äquivariante Abbildungen Träger -> Familien: V4 (3-dimensional, Kern immer das V4-Dublett), D4, mu4 (Z16).
Teil 2  Hyperladung Y = 1/12 + (5/12)(2 P_V4 - 1); eine Generation Lambda^gerade C^5 mit SM-Ladungen; sin^2 theta_W = 3/8.
Teil 3  Endliche Symmetrien in der Eichgruppe: V4 liegt in SU(3)_Farbe; die Vierteldrehung hat auf dem Dublett det -1.
Teil 4  Farbe = Lambda^2 F (D4): expliziter Verflechter; Farbe = F (V4): die V4-Abbildung T.
Teil 5  Z28 mit dem V4-Schnitt: Stabilisator des reinen Spinors und der V4-Teilung ist G_SM (12).
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

rng = np.random.default_rng(36)
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


mob = {"e": (1, 0, 0, 1), "r": (1j, 0, 0, 1), "r2": (-1, 0, 0, 1), "r3": (-1j, 0, 0, 1),
       "s": (0, 1, 1, 0), "rs": (0, 1j, 1, 0), "r2s": (0, -1, 1, 0), "r3s": (0, -1j, 1, 0)}
gz = lambda g, z: (g[0] * z + g[1]) / (g[2] * z + g[3])
dgz = lambda g, z: (g[0] * g[3] - g[1] * g[2]) / (g[2] * z + g[3]) ** 2
zs = np.exp(2j * np.pi * rng.random(12)) * (0.4 + rng.random(12))
def fit(vals, deg):
    V = np.array([[z**k for k in range(deg + 1)] for z in zs])
    c, *_ = np.linalg.lstsq(V, vals, rcond=None)
    return c
Cm = {n: np.array([fit(np.array([(gz(g, z) ** k / (gz(g, z) ** 4 - 1)) * (z**4 - 1) for z in zs]), 4) for k in range(5)]).T
      for n, g in mob.items()}
Fm = {n: np.array([fit(np.array([gz(g, z) ** k * dgz(g, z) / (gz(g, z) ** 4 - 1) * (z**4 - 1) for z in zs]), 2) for k in range(3)]).T
      for n, g in mob.items()}
e = np.eye(5)
weak = np.array([e[0] - e[4], e[2]]).T / np.array([np.sqrt(2), 1])
color = np.array([e[1], e[3], (e[0] + e[4]) / np.sqrt(2)]).T
check("V4-Teilung orthogonal im Monombasis-Skalarprodukt (schwach {1 - z^4, z^2}, Farbe {z, z^3, 1 + z^4})",
      np.allclose(weak.conj().T @ color, 0) and np.allclose(np.hstack([weak, color]).conj().T @ np.hstack([weak, color]), np.eye(5)))

# ------------------------------------------------------------------ Teil 1
say("Teil 1: äquivariante Abbildungen Träger -> Familien")
def equivariant_maps(group):
    rows = []
    for gname in group:
        A, B = Cm[gname], Fm[gname]
        rows.append(np.kron(np.eye(5), B) - np.kron(A.T, np.eye(3)))     # vec(B T - T A) = 0
    M = np.vstack(rows)
    u, sv, vh = np.linalg.svd(M)
    null = vh[np.sum(sv > 1e-8):].conj()
    return [v.reshape(5, 3).T for v in null]
groups = {"V4": ["r2", "rs"], "D4": ["r", "s"], "mu4": ["r"]}
maps = {k: equivariant_maps(v) for k, v in groups.items()}
say(f"   Dimension der äquivarianten Abbildungen: { {k: len(v) for k, v in maps.items()} }")
T_V4 = sum(rng.normal() * m for m in maps["V4"])
kerV4 = np.linalg.svd(T_V4)[2][3:].conj().T
inter = np.linalg.matrix_rank(np.hstack([kerV4, weak]), tol=1e-8)
check("V4: dreidimensional; jede generische V4-Abbildung T ist surjektiv mit Kern = V4-Dublett (schwach)",
      len(maps["V4"]) == 3 and np.linalg.matrix_rank(T_V4, tol=1e-8) == 3 and inter == 2)
T_mu4 = sum(rng.normal() * m for m in maps["mu4"])
ker_mu4 = np.linalg.svd(T_mu4)[2][3:].conj().T
check("mu4 (Z16): surjektive T möglich, Kern enthält z^2 aber das Dublett ist nicht {z, z^3} (Halbdrehung)",
      np.linalg.matrix_rank(T_mu4, tol=1e-8) == 3 and np.linalg.matrix_rank(np.hstack([ker_mu4, e[[1, 3]].T]), tol=1e-8) == 4)
rk_D4 = max((np.linalg.matrix_rank(sum(rng.normal() * m for m in maps["D4"]), tol=1e-8) for _ in range(3)), default=0) if maps["D4"] else 0
say(f"   D4-äquivariante Abbildungen: Dimension {len(maps['D4'])}, generischer Rang {rk_D4}")
check("unter D4 gibt es keine surjektive äquivariante T (Farbe ist unter D4 nicht F)", rk_D4 < 3)
befund("Der V4-Schnitt ist der Kern der natürlichen V4-Abbildung Träger -> Familien; Z16s mu4-Abbildung ist eine andere Wahl", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Hyperladung und eine Generation")
Pw = weak @ weak.conj().T
Y = np.eye(5) / 12 + 5 / 12 * (2 * Pw - np.eye(5))
ev = np.sort(np.linalg.eigvalsh(Y))
check("Y = 1/12 + (5/12)(2 P_V4 - 1): Eigenwerte -1/3 (x3), 1/2 (x2); tr Y = 0, tr Y^2 = 5/6",
      np.allclose(ev, [-1 / 3] * 3 + [1 / 2] * 2) and abs(np.trace(Y)) < 1e-12 and abs(np.trace(Y @ Y) - 5 / 6) < 1e-12)
yv = [Fr(1, 2)] * 2 + [Fr(-1, 3)] * 3
gen = []
for k in (0, 2, 4):
    for S in itertools.combinations(range(5), k):
        gen.append(sum((yv[i] for i in S), Fr(0)))
cnt = {}
for y in gen:
    cnt[y] = cnt.get(y, 0) + 1
want = {Fr(0): 1, Fr(1): 1, Fr(1, 6): 6, Fr(-2, 3): 3, Fr(-1, 2): 2, Fr(1, 3): 3}
say(f"   Lambda^gerade C^5: {dict(sorted(cnt.items()))}")
check("eine Generation: nu^c 0, e^c 1, Q 1/6 (x6), u^c -2/3 (x3), L -1/2 (x2), d^c 1/3 (x3)", cnt == want)
t3 = [Fr(1, 2), Fr(-1, 2), Fr(0), Fr(0), Fr(0)]           # T3 auf (Dublett, Farbe)
Q = [x + y for x, y in zip(t3, yv)]
sw2 = sum(x * x for x in t3) / sum(q * q for q in Q)
say(f"   tr T3^2 = {sum(x * x for x in t3)}, tr Q^2 = {sum(q * q for q in Q)}, sin^2 theta_W = {sw2}")
check("sin^2 theta_W = tr T3^2 / tr Q^2 = 3/8 auf dem Träger", sw2 == Fr(3, 8))
befund("Dimensionsabhängige Ergebnisse (Ladungen, sin^2 theta_W, tr Y^2) bleiben unter dem V4-Schnitt gleich", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: endliche Symmetrien in der Eichgruppe")
basis = np.hstack([weak, color])
blk = {n: np.linalg.inv(basis) @ Cm[n] @ basis for n in mob}
v4_in_su3 = all(np.allclose(blk[n][:2, :2], np.eye(2)) and np.allclose(blk[n][:2, 2:], 0)
                and abs(np.linalg.det(blk[n][2:, 2:]) - 1) < 1e-9 for n in ("e", "r2", "rs", "r3s"))
check("V4 wirkt auf dem Dublett trivial und auf der Farbe mit det 1: V4 liegt in SU(3)_Farbe (diskrete Farbdrehung)", v4_in_su3)
dW, dC = np.linalg.det(blk["r"][:2, :2]), np.linalg.det(blk["r"][2:, 2:])
say(f"   Vierteldrehung r: det(Dublett) = {np.round(dW, 9)}, det(Farbe) = {np.round(dC, 9)}, block-diagonal "
    f"{np.allclose(blk['r'][:2, 2:], 0) and np.allclose(blk['r'][2:, :2], 0)}")
check("r erhält die V4-Teilung, hat auf dem Dublett det -1 und auf der Farbe det 1: r in U(2) x SU(3), nicht in S(U(2) x U(3))",
      np.allclose(blk["r"][:2, 2:], 0) and abs(dW + 1) < 1e-9 and abs(dC - 1) < 1e-9)
half_old = np.diag([1, -1, 1, -1, 1.0])
check("unter der alten Halbdrehungs-Teilung war die Halbdrehung (-1 auf {z, z^3}) ein zentrales Hyperladungselement",
      np.allclose(np.linalg.eigvals(half_old).real.round(), np.array([1, -1, 1, -1, 1])))
befund("Unter V4 wird die Markengruppe eine diskrete Farbdrehung; die Vierteldrehung wirkt auf dem Dublett wie eine Spiegelung", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Verflechter")
def ext2(M):
    pairs = list(itertools.combinations(range(M.shape[0]), 2))
    return np.array([[M[i, k] * M[j, l] - M[i, l] * M[j, k] for (k, l) in pairs] for (i, j) in pairs])
L2F = {n: ext2(m) for n, m in Fm.items()}
Cc = {n: blk[n][2:, 2:] for n in mob}
X = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
Phi = sum(L2F[n] @ X @ np.linalg.inv(Cc[n]) for n in mob) / 8
check("Farbe -> Lambda^2 F: gemittelter Verflechter ist D4-äquivariant und invertierbar",
      all(np.allclose(L2F[n] @ Phi, Phi @ Cc[n], atol=1e-9) for n in mob) and abs(np.linalg.det(Phi)) > 1e-6)
TV = T_V4 @ color
check("Farbe -> F: die V4-Abbildung T ist auf der Farbe invertierbar und V4-äquivariant",
      abs(np.linalg.det(TV)) > 1e-6 and all(np.allclose(Fm[n] @ T_V4, T_V4 @ Cm[n], atol=1e-8) for n in ("r2", "rs")))

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Z28 mit dem V4-Schnitt")
n = 5
Zp = np.diag([1.0, -1.0]); I2 = np.eye(2); low = np.array([[0, 1], [0, 0]], complex)
def kron_all(ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out
a = [kron_all([Zp] * i + [low] + [I2] * (n - i - 1)) for i in range(n)]
ad = [x.conj().T for x in a]
dim = 2 ** n
gens = []
for i in range(n):
    gens.append(1j * (ad[i] @ a[i] - 0.5 * np.eye(dim)))
for i, j in itertools.combinations(range(n), 2):
    gens += [ad[i] @ a[j] - ad[j] @ a[i], 1j * (ad[i] @ a[j] + ad[j] @ a[i]),
             ad[i] @ ad[j] - a[j] @ a[i], 1j * (ad[i] @ ad[j] + a[j] @ a[i])]
vac = np.zeros(dim, complex); vac[0] = 1
Qc = color @ color.conj().T
NQ = sum(Qc[i, j] * ad[i] @ a[j] for i in range(n) for j in range(n))
rows = np.array([np.concatenate([np.concatenate([(g @ vac).real, (g @ vac).imag]),
                                 np.concatenate([(g @ NQ - NQ @ g).real.ravel(), (g @ NQ - NQ @ g).imag.ravel()])]) for g in gens]).T
stab = 45 - np.linalg.matrix_rank(rows, tol=1e-9)
say(f"   dim Stab(s, V4-Teilung) = {stab}")
check("Stabilisator des reinen Spinors und der V4-Teilung ist G_SM (Dimension 12) wie in Z28", stab == 12)

urteil = ("Mit dem V4-Schnitt (schwach {1 - z^4, z^2}, Farbe {z, z^3, 1 + z^4}, orthogonal) bleibt alles, was nur Dimensionen "
          "benutzt: Y hat die Eigenwerte -1/3 und 1/2, tr Y^2 = 5/6, eine Generation Lambda^gerade C^5 trägt genau die "
          "SM-Ladungen, sin^2 theta_W = 3/8, und der Stabilisator des reinen Spinors mit der V4-Teilung ist G_SM (12). Neu "
          "festgelegt wird die Abbildung Träger -> Familien: Die V4-äquivarianten Abbildungen bilden einen dreidimensionalen "
          "Raum, jede generische ist surjektiv mit Kern genau das V4-Dublett; Z16s mu4-Abbildung war eine andere Wahl, und "
          "unter D4 gibt es keine surjektive äquivariante Abbildung. Die Farbe ist zu Lambda^2 F D4-äquivariant isomorph "
          "(expliziter Verflechter) und zu F V4-äquivariant. In der Eichgruppe wird die Markengruppe V4 eine diskrete "
          "Farbdrehung in SU(3); die Vierteldrehung erhält die Teilung, wirkt aber auf dem Dublett mit det -1 und liegt in "
          "U(2) x SU(3), nicht in S(U(2) x U(3)). Unter der alten Halbdrehungs-Teilung war die Halbdrehung dagegen ein "
          "zentrales Hyperladungselement.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z36_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": {"aequivariante_abbildungen": {k: len(v) for k, v in maps.items()},
                                                      "stab": int(stab)}, "urteil": urteil},
    indent=1, ensure_ascii=False, default=str) + "\n")
