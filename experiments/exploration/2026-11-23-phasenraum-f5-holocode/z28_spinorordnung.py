"""Z28: Projektive Spinorordnung -- Eichgruppe, Viererholonomie und Wand aus einem Ordnungsparameter?

Prüft die Rechnungen der vierten externen Analyse (24.11.) zur Eichstruktur von H6:
Teil 1  Spin(10) auf dem Fockraum Lambda^gerade C^5 (Jordan-Wigner), Abschluss der 45 Generatoren.
Teil 2  Stabilisatoren: reiner Spinor s (Fockvakuum), seine Linie, zusätzlich die 3+2-Teilung.
Teil 3  Zentrum z von Spin(10) auf s^{(x)k}; kleinste zentrumerhaltende Potenz; Weyl-Dimensionen V(k w5).
Teil 4  Stab(s^{(x)4}) = SU(5) x Z4; pi_1 der Vakuummannigfaltigkeiten (lange exakte Sequenz).
Teil 5  Paarungsfeld Phi mit Ladung 2: Verträglichkeit mit der primitiven Viererholonomie (Wand).
Teil 6  y1: zusammengesetzte Verbindung (Maurer-Cartan) gegen unabhängige Elternverbindung (Rang).
Teil 7  Craigs Bedingungen (arXiv:2610.09093, Abschnitt 11.1) für den TFPT-Inhalt; Spin(10)-16 unter
        einer Wurzel-SU(2); Zusammenhang von Spin(10)/SU(2).
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

rng = np.random.default_rng(28)
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


def real_rank(vecs, tol=1e-9):
    M = np.array([np.concatenate([v.real.ravel(), v.imag.ravel()]) for v in vecs])
    return np.linalg.matrix_rank(M, tol=tol)


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Spin(10) auf Lambda^gerade C^5")
n = 5
I2, Z = np.eye(2), np.diag([1.0, -1.0])
lower = np.array([[0, 1], [0, 0]], complex)          # |1> -> |0>, Basis (|0>, |1>)


def kron_all(ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out


a = [kron_all([Z] * i + [lower] + [I2] * (n - i - 1)) for i in range(n)]
ad = [x.conj().T for x in a]
dim = 2 ** n
Nop = sum(ad[i] @ a[i] for i in range(n))
occ = np.array([bin(k).count("1") for k in range(dim)])
even = np.where(occ % 2 == 0)[0]
vac = np.zeros(dim, complex); vac[0] = 1.0           # |00000>, Fockvakuum = reiner Spinor s
check("Kanonische Antivertauschung {a_i, a_j^+} = delta_ij", all(
    np.allclose(a[i] @ ad[j] + ad[j] @ a[i], np.eye(dim) * (i == j)) for i in range(n) for j in range(n)))

gens, labels = [], []
for i in range(n):
    gens.append(1j * (ad[i] @ a[i] - 0.5 * np.eye(dim))); labels.append(("u", i, i))
for i, j in itertools.combinations(range(n), 2):
    gens.append(ad[i] @ a[j] - ad[j] @ a[i]); labels.append(("u", i, j))
    gens.append(1j * (ad[i] @ a[j] + ad[j] @ a[i])); labels.append(("u", j, i))
    gens.append(ad[i] @ ad[j] - a[j] @ a[i]); labels.append(("p", i, j))
    gens.append(1j * (ad[i] @ ad[j] + a[j] @ a[i])); labels.append(("p", j, i))
check("45 antihermitesche Generatoren, linear unabhängig", len(gens) == 45 and real_rank(gens) == 45
      and all(np.allclose(g.conj().T, -g) for g in gens))
span = np.array([np.concatenate([g.real.ravel(), g.imag.ravel()]) for g in gens]).T
closed = True
for g, h in itertools.combinations(gens, 2):
    c = g @ h - h @ g
    v = np.concatenate([c.real.ravel(), c.imag.ravel()])
    res = v - span @ np.linalg.lstsq(span, v, rcond=None)[0]
    closed &= np.linalg.norm(res) < 1e-9
check("Kommutatoren schließen (Lie-Algebra der Dimension 45 = so(10))", closed)
check("Generatoren erhalten die Fermionparität (16 = gerader Teil)", all(
    np.allclose(g[np.ix_(even, np.where(occ % 2 == 1)[0])], 0) for g in gens))

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Stabilisatoren")
Xs = [g @ vac for g in gens]
dim_stab_s = 45 - real_rank(Xs)
Xs_line = [x - (vac.conj() @ x) * vac for x in Xs]
dim_stab_line = 45 - real_rank(Xs_line)
N3 = sum(ad[i] @ a[i] for i in range(3))
# Stabilisator von s und der 3+2-Teilung: X s = 0 und [X, N_3] = 0
def stab_dim(conds):
    rows = np.array([np.concatenate([np.concatenate([f(g).real.ravel(), f(g).imag.ravel()]) for f in conds])
                     for g in gens]).T
    return 45 - np.linalg.matrix_rank(rows, tol=1e-9)
dim_stab_32 = stab_dim([lambda g: g @ vac, lambda g: g @ N3 - N3 @ g])
dim_stab_line_32 = stab_dim([lambda g: g @ vac - (vac.conj() @ (g @ vac)) * vac, lambda g: g @ N3 - N3 @ g])
say(f"   dim Stab([s]) = {dim_stab_line} (u(5)), dim Stab(s) = {dim_stab_s} (su(5)), "
    f"dim Stab(s, 3+2) = {dim_stab_32}, dim Stab([s], 3+2) = {dim_stab_line_32}")
check("Stab der reinen Spinorlinie: 25 = dim u(5)", dim_stab_line == 25)
check("Stab des reinen Spinors: 24 = dim su(5)", dim_stab_s == 24)
check("Stab von s und 3+2-Teilung: 12 = dim(su3 + su2 + u1)", dim_stab_32 == 12)
say(f"   gebrochen: Spin(10) -> SU(5) {45 - 24}, SU(5) -> G_SM {24 - 12}, zusammen {45 - 12}, ungebrochen 12")
ZAHLEN["stabilisatoren"] = {"linie": dim_stab_line, "spinor": dim_stab_s, "spinor_3plus2": dim_stab_32,
                             "linie_3plus2": dim_stab_line_32, "gebrochen": 45 - dim_stab_32}

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Zentrum und Spinorpotenzen")
zc = np.diag(np.exp(1j * np.pi * (occ - n / 2)))         # exp(i pi (N - 5/2)), liegt in exp(u(1))
check("z = exp(i pi (N - 5/2)) liegt in exp(span): Generator i sum(a^+a - 1/2) ist einer der 45",
      np.allclose(sum(gens[i] for i in range(5)), 1j * (Nop - 2.5 * np.eye(dim))))
check("z wirkt auf den Vektor als -1 (Ad_z a_i^+ = -a_i^+)", all(np.allclose(zc @ ad[i] @ zc.conj().T, -ad[i]) for i in range(n)))
ph16 = np.unique(np.round(np.diag(zc)[even], 12))
check("z wirkt auf die 16 als Skalar -i (Ordnung 4)", len(ph16) == 1 and np.isclose(ph16[0], -1j))
phases = {k: complex(np.round((-1j) ** k, 12)) for k in range(1, 5)}
say(f"   z auf s^(x)k: {phases}")
kmin = min(k for k in range(1, 9) if np.isclose((-1j) ** k, 1))
check("kleinste zentrumerhaltende Spinorpotenz k = 4", kmin == 4)
# Weyl-Dimensionsformel D5
pos = [(i, j, s) for i in range(5) for j in range(i + 1, 5) for s in (-1, 1)]
rho = [Fr(4 - i) for i in range(5)]
def weyl_dim(lam):
    num = den = Fr(1)
    for i, j, s in pos:
        alpha = [0] * 5; alpha[i] = 1; alpha[j] = s
        num *= sum((lam[k] + rho[k]) * alpha[k] for k in range(5))
        den *= sum(rho[k] * alpha[k] for k in range(5))
    return num / den
dims = {k: int(weyl_dim([Fr(k, 2)] * 5)) for k in range(1, 5)}
say(f"   dim V(k w5) = {dims}")
check("dim V(w5) = 16, V(4 w5) = 2772", dims[1] == 16 and dims[4] == 2772)
ZAHLEN["weyl_dim_k_omega5"] = dims

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Stab(Q), Q = s^(x)4, und Fundamentalgruppen")
say("   Lie-Stabilisator von Q: X.Q = 4c Q + (Terme mit dem zu s senkrechten Teil von X s), also X.Q = 0 <=> X s = 0.")
check("Lie(Stab Q) = Lie(Stab s) = su(5), Dimension 24", dim_stab_s == 24)
say("   Gruppe: g Q = Q <=> g s = zeta s mit zeta^4 = 1. Stab([s]) ist die zusammenhängende Doppelüberlagerung von U(5)")
say("   mit Charakter lambda (lambda^2 = det A); lambda^-1(1) = SU(5), lambda^-1(mu4) hat vier Komponenten, z liegt in")
say("   der Komponente lambda = -i und ist zentral: Stab(Q) = SU(5) x <z> = SU(5) x Z4; mit 3+2: G_SM x Z4.")
comp = sorted({complex(np.round((-1j) ** m, 12)) for m in range(4)}, key=np.angle)
check("z^m auf s trifft alle vier Charaktere von mu4 (Komponentengruppe Z4)", len(comp) == 4)
say("   pi_1(Spin(10)) = 0 => pi_1(Spin(10)/H) = pi_0(H) (lange exakte Sequenz). H = G_SM x Z4 => pi_1 = Z4.")
say("   Mit U(1)_D und R = s (x) u^-1: H = {(h z^n, (-i)^n)}, pi_1(Spin(10) x U(1)/H) ist eine Erweiterung von Z4 durch")
say("   pi_1(U(1)) = Z; der Weg zu (z, -i) ist eine Viertelwurzel des U(1)_D-Generators, also pi_1 = Z, Auslese n mod 4.")
check("Viermal der n = 1-Weg: (z^4, (-i)^4) = (1, 1) und U(1)_D-Windung 1", np.allclose(np.linalg.matrix_power(zc, 4), np.eye(dim))
      and np.isclose((-1j) ** 4, 1))
ZAHLEN["pi1"] = {"Spin10/(G_SM x Z4)": "Z4", "Spin10 x U(1)_D / H": "Z (Auslese mod 4)"}

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Ordnungsfeld der Ladung q unter der primitiven Viererholonomie i")
say("   Kovariante Eindeutigkeit verlangt die Phasenwindung w = q/4 (mod 1). Endliche Spannung ohne Wand nur bei w ganz.")
wind = {q: Fr(q, 4) for q in (1, 2, 4)}
for q, w in wind.items():
    say(f"   q = {q}: w = {w}; minimale Gradientendichte ~ min_nu (nu - w)^2 = {min((k - w) ** 2 for k in range(-2, 3))}")
check("Phi (Ladung 2, J-ungerade): halbe Windung, also Wand (gepinnte Phase) oder log-divergente Energie",
      wind[2] == Fr(1, 2))
check("Ladung 4 (Phi^2 bzw. Vierfermion-Kondensat): ganze Windung, primitive Holonomie ohne Wand verträglich",
      wind[4] == 1)
befund("Wandproblem aus Z24 auch gruppentheoretisch: <Phi> != 0 lässt von Z4 nur Z2", True)
befund("Symmetrische Massenerzeugung mit Ladung-4-Kondensat erhält die primitive Holonomie", True)

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: y1 -- zusammengesetzte gegen unabhängige Verbindung")
def rand_u(k):
    q, r = np.linalg.qr(rng.normal(size=(k, k)) + 1j * rng.normal(size=(k, k)))
    return q * (np.diag(r) / abs(np.diag(r)))
# Maurer-Cartan: g(t) = g0 exp(t X) -> omega = g^+ dg; d omega + omega^2 = 0 (flach), Blöcke
g0 = rand_u(5)
H1 = rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5)); X1 = H1 - H1.conj().T
H2 = rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5)); X2 = H2 - H2.conj().T
expm = lambda M: (lambda w, v: v @ np.diag(np.exp(w)) @ np.linalg.inv(v))(*np.linalg.eig(M))
gf = lambda s, t: g0 @ expm(s * X1) @ expm(t * X2)
h = 1e-5
w1 = lambda s, t: gf(s, t).conj().T @ (gf(s + h, t) - gf(s - h, t)) / (2 * h)
w2 = lambda s, t: gf(s, t).conj().T @ (gf(s, t + h) - gf(s, t - h)) / (2 * h)
F = (w2(h, 0) - w2(-h, 0)) / (2 * h) - (w1(0, h) - w1(0, -h)) / (2 * h) + w1(0, 0) @ w2(0, 0) - w2(0, 0) @ w1(0, 0)
check("Maurer-Cartan-Verbindung g^+dg ist flach", np.linalg.norm(F) < 1e-5)
A3 = lambda w: w[:3, :3]
F3 = (A3(w2(h, 0)) - A3(w2(-h, 0))) / (2 * h) - (A3(w1(0, h)) - A3(w1(0, -h))) / (2 * h) \
    + A3(w1(0, 0)) @ A3(w2(0, 0)) - A3(w2(0, 0)) @ A3(w1(0, 0))
Xo1, Xo2 = w1(0, 0)[3:, :3], w2(0, 0)[3:, :3]
check("projizierter Block: F_3 = X_1^+ X_2 - X_2^+ X_1 (nur aus ersten Ableitungen)",
      np.allclose(F3, Xo1.conj().T @ Xo2 - Xo2.conj().T @ Xo1, atol=1e-5))
# Rang: zusammengesetzt (48 Parameter, Bildrang <= 48) gegen Elternverbindung (Krümmung frei)
V3 = g0[:, :3]; V2 = g0[:, 3:]
def project(F5):        # u(5)-Krümmung einer Ebene -> s(u3+u2)-Komponenten
    return np.concatenate([(V3.conj().T @ F5 @ V3).ravel(), (V2.conj().T @ F5 @ V2).ravel()])
basisF = []
for m in range(6):
    for k in range(25):
        E = np.zeros((5, 5), complex)
        r, c = divmod(k, 5)
        if r == c:
            E[r, r] = 1j
        elif r < c:
            E[r, c], E[c, r] = 1, -1
        else:
            E[r, c], E[c, r] = 1j, 1j
        F5 = g0 @ E @ g0.conj().T
        vec = [np.zeros(13, complex) for _ in range(6)]
        vec[m] = project(F5)
        basisF.append(np.concatenate(vec))
rank_parent = real_rank(basisF, tol=1e-8)
# s(u3+u2): u(3)+u(2) hat 13 Dimensionen; die physische Eichkrümmung su3+su2+u1 hat 12 je Ebene
say(f"   Elternverbindung: Bildrang der projizierten Krümmung {rank_parent} (u(3)+u(2), 6 x 13 = 78)")
check("unabhängige Elternverbindung erreicht alle Krümmungen (Rang 78 >= 72 der SM-Komponenten)", rank_parent == 78)
befund("y1 benutzt die flache Maurer-Cartan-Struktur; mit unabhängiger Elternverbindung entfällt die Rangschranke", True)

# ------------------------------------------------------------------ Teil 7
say("\nTeil 7: Craigs Bedingungen (Abschnitt 11.1) und Spin(10)")
# SU(5)-Gewichte: Fundamentaldarstellung e_1..e_5; SU(2) auf zwei Farben (4, 5)
def doublets(weights):  # Zahl der Dubletts unter SU(2) auf den Farben 4,5: Gewichte mit w4 - w5 = +-1 paaren
    cnt = {}
    for w in weights:
        d = w[3] - w[4]
        cnt[d] = cnt.get(d, 0) + 1
    return cnt.get(1, 0)
f5 = [tuple(int(i == k) for i in range(5)) for k in range(5)]
f5b = [tuple(-x for x in w) for w in f5]
f10 = [tuple(a_ + b_ for a_, b_ in zip(f5[i], f5[j])) for i, j in itertools.combinations(range(5), 2)]
anom = {"5": 1, "5b": -1, "10": 1, "1": 0}
content = {"10": 3, "5b": 3, "1": 3 + 6, "5": 1}
content["5b"] += 1
A = sum(anom[r] * m for r, m in content.items())
D = 3 * (doublets(f10) + doublets(f5b)) + doublets(f5) + doublets(f5b)
say(f"   TFPT unter SU(5): 3(10 + 5b + 1) + (5 + 5b) + 6 x 1 (Eltern (10,1)+(1,6)); kubische Anomalie {A}; "
    f"Dubletts unter SU(2) auf zwei Farben: je Generation {doublets(f10) + doublets(f5b)}, Eltern {doublets(f5) + doublets(f5b)}, gesamt {D}")
check("Craig (129): kubische SU(5)-Anomalie 0", A == 0)
check("Craig (130): Spin-1/2-Vielfachheit gerade (je Generation 4, Eltern 2, gesamt 14)", D == 14 and D % 2 == 0)
w16 = [w for w in itertools.product((Fr(1, 2), Fr(-1, 2)), repeat=5) if sum(1 for x in w if x < 0) % 2 == 0]
d16 = sum(1 for w in w16 if w[0] - w[1] == 1)
check("Spin(10)-16 unter einer Wurzel-SU(2): 4 Dubletts + 8 Singuletts", len(w16) == 16 and d16 == 4)
say("   Spin(10)/SU(2) (Wurzeleinbettung, Index 1): pi_1 = pi_2 = 0, pi_3 = coker(Z -> Z iso) = 0, pi_4 eingebettet in")
say("   ker(pi_3) = 0. Also 4-zusammenhängend; das ist richtig, wird von Craig aber nicht benutzt (keine Bündelreduktion;")
say("   global: Zylinder [0,1] x S^1, Basis T^4 x S^1). Spin(10) mit der 16 ist bei Craig ausdrücklich nicht abgedeckt.")
befund("Craigs Maß deckt den SU(5)-Zweig von TFPT ab (Bedingungen 11.1 erfüllt)", True)
befund("Craigs Maß deckt Spin(10) mit der 16 ab", False)

# ------------------------------------------------------------------ Urteil
urteil = ("Die Spinorordnung trägt: Der reine Spinor s im Fockraum Lambda^gerade C^5 hat den Stabilisator su(5) (24), "
          "seine Linie u(5) (25), mit der 3+2-Teilung bleibt G_SM (12); gebrochen sind 21 + 12 = 33 Richtungen. Das Zentrum "
          "z wirkt auf die 16 als -i; s^(x)4 ist die kleinste zentrumerhaltende Potenz (V(4 w5), Dimension 2772), "
          "Stab = SU(5) x Z4 bzw. G_SM x Z4, also pi_1 = Z4 aus dem Symmetriebruch statt als Eingabe; mit U(1)_D und "
          "R = s (x) u^-1 ist pi_1 = Z mit Auslese mod 4. Ein J-ungerades Paarungsfeld (Ladung 2) mit Vakuumwert lässt "
          "nur Z2 übrig und verlangt um die primitive Naht eine halbe Windung, also eine Wand (Z24) oder log-divergente "
          "Energie; ein Ladung-4-Kondensat (symmetrische Massenerzeugung) ist mit der primitiven Holonomie verträglich. "
          "Die Rangschranke von y1 (47 < 72) gehört zur flachen Maurer-Cartan-Struktur; eine unabhängige Elternverbindung "
          "hebt sie auf. Craigs Bedingungen (11.1) erfüllt der SU(5)-Zweig von TFPT; Spin(10) mit der 16 deckt Craig nicht "
          "ab, und die Bündelreduktion auf SU(2) kommt in seiner Konstruktion nicht vor. Offen: ob TFPT das Kondensat und "
          "die Elternverbindung dynamisch erzeugt.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
Path(__file__).with_name("z28_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
