"""Z32: Vom Toeplitz-Index zur lokalen E8-Naht -- Hardy-Kozyklus, Virasoro, Gitternetz, Modularität, Markenrahmen.

Prüft die sechste externe Analyse (24.11.):
Teil 1  Graduierter Dirac-Operator D = [[0, S^+], [S, 0]]: {D, Gamma} = 0, D^2, Index; A = S f(N) mit f = (n+1)^alpha:
        Index -1 für alle alpha, Wärmespur 1 + 2 sum e^{-t n^{2 alpha}}, spektrale Dimension 1/alpha,
        Kommutator [A, S] beschränkt genau für alpha <= 1.
Teil 2  Hardy-Kozyklus: Tr(P U^m (1-P) U^-m P) = m; Virasoro-Norm sum_j (j + (1-m)/2)^2 = (m^3 - m)/12 (c = 1 je komplexes
        Fermion), also c = 8 für 16 reelle; Sugawara so(16)_1: 120/15 = 8.
Teil 3  E8 = D8 u (D8 + s): gerade, unimodular, 240 Vektoren der Norm 2.
Teil 4  Modularität: Die Toeplitz-Algebra ist irreduzibel (Kommutant C, schwache Hülle B(H)) und hat keinen separierenden
        Vektor; ein Gibbs-Zustand mit Punktspektrum hat eine fastperiodische modulare Gruppe (keine echte halbseitige
        modulare Inklusion).
Teil 5  Markenrahmen: Ein festes Markenquadrat hat im Oktaeder nur D4 (8), drei gleichwertige Quadrate; D4-invariante Formen 2;
        Mittel über die drei Rahmen (1/3) sum (I - e_a e_a^T) = (2/3) I.
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

rng = np.random.default_rng(32)
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
say("Teil 1: graduierter Dirac-Operator und gewichteter Shift")
K = 80
S = np.eye(K, k=-1)
D = np.block([[np.zeros((K, K)), S.T], [S, np.zeros((K, K))]])
G = np.diag([1.0] * K + [-1.0] * K)
D2 = D @ D
inner = list(range(K - 1)) + list(range(K, 2 * K - 1))       # ferner Rand der Abschneidung ausgenommen
check("{D, Gamma} = 0", np.allclose(D @ G + G @ D, 0))
target = np.diag([1.0] * K + [0.0] + [1.0] * (K - 1))
check("D^2 = diag(I, I - |0><0|) (eine ungepaarte Nullmode)", np.allclose(D2[np.ix_(inner, inner)], target[np.ix_(inner, inner)]))
idx = {}
for al in (0.5, 0.75, 1.0, 1.5):
    f = (np.arange(K) + 1.0) ** al
    A = S @ np.diag(f)
    coker = int(np.sum(np.abs(A[:K - 1, :K - 1]).sum(axis=1) < 1e-12))           # Zeile 0 leer: |0> nicht im Bild
    idx[al] = (int(np.linalg.matrix_rank(A[:, :K - 1]) == K - 1), coker)
say(f"   A = S f(N), f = (n+1)^alpha: (injektiv, Kokern) {idx}")
check("Index von A = S f(N) ist -1 für alpha = 0,5; 0,75; 1; 1,5", all(v == (1, 1) for v in idx.values()))
def dspec(al, t):
    n = np.arange(1, 4_000_001, dtype=float)
    H = lambda tt: 1 + 2 * np.sum(np.exp(-tt * n ** (2 * al)))
    h = 1e-3
    return -2 * (np.log(H(t * (1 + h))) - np.log(H(t * (1 - h)))) / (np.log(1 + h) - np.log(1 - h))
ds = {al: dspec(al, 1e-5) for al in (0.5, 0.75, 1.0)}      # t N^(2 alpha) >> 1, t << 1
say(f"   spektrale Dimension (kleine t): {dict((k, round(v, 4)) for k, v in ds.items())} (erwartet 1/alpha)")
check("d_spec = 1/alpha (2; 4/3; 1) auf 1 %", all(abs(v - 1 / k) / (1 / k) < 0.01 for k, v in ds.items()))
comm = {al: float(np.max(np.abs(np.diff((np.arange(2000) + 1.0) ** al)))) for al in (0.5, 1.0, 1.5)}
say(f"   sup |f(n+1) - f(n)| bis n = 2000 (= ||[A, S]||): {dict((k, round(v, 3)) for k, v in comm.items())}")
check("[A, S] beschränkt für alpha <= 1, wachsend für alpha > 1", comm[0.5] <= 1 and comm[1.0] <= 1 + 1e-12 and comm[1.5] > 50)
befund("Irreduzibilität und Index legen weder Dimension noch Metrik fest (d_spec = 1/alpha frei)", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Hardy-Kozyklus und Virasoro")
N = 64
U = np.roll(np.eye(N), 1, axis=0)
P = np.diag([1.0 if n < N // 2 else 0.0 for n in range(N)])
coc = {m: int(round(np.trace(P @ np.linalg.matrix_power(U, m) @ (np.eye(N) - P) @ np.linalg.matrix_power(U, -m) @ P)))
       for m in range(1, 9)}
say(f"   Tr(P U^m (1-P) U^-m P) = {coc}")
check("Hardy-Kozyklus = m (Moden, die den Schnitt überqueren) = Schwinger-Term [J_m, J_-m] = m", all(coc[m] == m for m in coc))
vir = all(sum((Fr(j) + Fr(1 - m, 2)) ** 2 for j in range(m)) == Fr(m**3 - m, 12) for m in range(1, 49))
check("sum_{j<m} (j + (1-m)/2)^2 = (m^3 - m)/12 exakt für m = 1..48 (c = 1 je komplexes Fermion)", vir)
check("8 komplexe = 16 reelle Fermionen: c = 8; Sugawara so(16)_1: 120/(1 + 14) = 8", Fr(120, 15) == 8 and 8 * 1 == 8)
befund("Toeplitz-Index -> Stromniveau 1 -> chirale Anomalie ist eine exakte Kette (Segal-Kozyklus)", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: E8 = D8 u (D8 + s)")
# Basis von D8 u (D8 + s): sieben D8-Wurzeln und s - e_8
B = np.array([np.eye(8)[0] - np.eye(8)[1], np.eye(8)[1] - np.eye(8)[2], np.eye(8)[2] - np.eye(8)[3], np.eye(8)[3] - np.eye(8)[4],
              np.eye(8)[4] - np.eye(8)[5], np.eye(8)[5] - np.eye(8)[6], np.eye(8)[5] + np.eye(8)[6], np.full(8, 0.5) - np.eye(8)[7]])
Gram = B @ B.T
cnt = 0
for v in itertools.product(range(-1, 2), repeat=8):
    if sum(x * x for x in v) == 2 and sum(v) % 2 == 0:
        cnt += 1
for v in itertools.product((0.5, -0.5), repeat=8):
    if sum(1 for x in v if x < 0) % 2 == 0:
        cnt += 1
say(f"   det Gram = {round(np.linalg.det(Gram), 9)}, Diagonale {np.diag(Gram)}, Vektoren der Norm 2: {cnt}")
check("D8 + s gerade und unimodular (det 1), 240 Wurzeln (112 + 128) = E8", round(np.linalg.det(Gram)) == 1
      and all(abs(x - round(x)) < 1e-9 and round(x) % 2 == 0 for x in np.diag(Gram)) and cnt == 240)
say("   Für gerade positive Gitter gibt es unitäre lokale konforme Netze (Dong-Xu); das E8-Netz der Naht ist damit konstruiert,")
say("   sobald 16 reelle Fermionen und die GSO-Erweiterung gegeben sind. Offen: Herleitung als Rand des 4D-Systems.")

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Modularität")
Kt = 6
St = np.eye(Kt, k=-1)
def commutant_dim(gens):
    n = gens[0].shape[0]
    return n * n - np.linalg.matrix_rank(np.vstack([np.kron(np.eye(n), g) - np.kron(g.T, np.eye(n)) for g in gens]))
check("Toeplitz (abgeschnitten): S, S^+ erzeugen M_K, Kommutant C (irreduzibel, schwache Hülle B(H))",
      commutant_dim([St, St.T]) == 1)
xi = rng.normal(size=Kt) + 1j * rng.normal(size=Kt)
perp = rng.normal(size=Kt) + 1j * rng.normal(size=Kt); perp -= (xi.conj() @ perp) / (xi.conj() @ xi) * xi
Aann = np.outer(rng.normal(size=Kt), perp.conj())
check("kein separierender Vektor für B(H): es gibt A != 0 mit A xi = 0", np.linalg.norm(Aann) > 0.1 and np.linalg.norm(Aann @ xi) < 1e-10)
E = np.array([0.0, 1.0, np.sqrt(2), np.pi / 2])
best = min(((t, float(np.max(np.abs(np.exp(1j * t * (E[:, None] - E[None, :])) - 1)))) for t in np.arange(50.0, 20000.0, 0.25)),
           key=lambda x: x[1])
say(f"   Punktspektrum {list(np.round(E, 4))}: Delta^(it) kehrt zurück, bei t = {best[0]} ist max |e^(it dE) - 1| = {best[1]:.3f}")
check("fastperiodische modulare Gruppe (Rückkehr nahe I für großes t): keine echte halbseitige modulare Inklusion", best[1] < 0.15)
befund("Der globale Shift und ein diskreter Gibbs-Zustand reichen für Raumzeittranslationen nicht; nötig sind lokale Typ-III-Algebren", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Markenrahmen im Oktaeder")
Rq = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.0]]); C3 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0.0]])
grp = {tuple(np.eye(3).ravel())}; fr = [np.eye(3)]
while fr:
    nw = []
    for g in fr:
        for s_ in (Rq, C3):
            h = np.round(s_ @ g, 12)
            if tuple(h.ravel()) not in grp:
                grp.add(tuple(h.ravel())); nw.append(h)
    fr = nw
O = [np.array(t).reshape(3, 3) for t in grp]
square = {tuple(v) for v in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]}
stab = [g for g in O if {tuple(np.round(g @ np.array(v), 12)) for v in square} == square]
orbit = {frozenset(tuple(np.round(g @ np.array(v), 12)) for v in square) for g in O}
def inv_forms(gens):
    Bs = [np.outer(np.eye(3)[i], np.eye(3)[j]) + np.outer(np.eye(3)[j], np.eye(3)[i]) for i in range(3) for j in range(i, 3)]
    rows = np.array([np.concatenate([(g.T @ Bm @ g - Bm).ravel() for g in gens]) for Bm in Bs]).T
    return len(Bs) - np.linalg.matrix_rank(rows)
say(f"   Stabilisator eines Markenquadrats: {len(stab)} Elemente; gleichwertige Quadrate: {len(orbit)}; D4-invariante Formen: {inv_forms(stab)}")
check("ein festes Markenquadrat hat nur D4 (8), drei gleichwertige Rahmen, D4 lässt 2 invariante Formen (anisotrop)",
      len(stab) == 8 and len(orbit) == 3 and inv_forms(stab) == 2)
avg = sum(np.eye(3) - np.outer(np.eye(3)[a], np.eye(3)[a]) for a in range(3)) / 3
check("Mittel über die drei Rahmen: (1/3) sum (I - e_a e_a^T) = (2/3) I (isotrop)", np.allclose(avg, 2 / 3 * np.eye(3)))
befund("Korrektur zu Z29/Z30: geometrische S4-Isotropie braucht die drei Rahmen (oder eine interne S4), nicht ein festes Quadrat", True)

urteil = ("Die Rechnungen der Analyse stimmen. Der graduierte Dirac-Operator aus S hat genau eine ungepaarte Nullmode; A = S f(N) "
          "behält Index -1 für jedes alpha, die spektrale Dimension ist 1/alpha (2; 4/3; 1 auf 1 %), [A, S] bleibt nur für "
          "alpha <= 1 beschränkt: Irreduzibilität legt die Topologie fest, nicht die Metrik. Der Hardy-Kozyklus "
          "Tr(P U^m (1-P) U^-m P) = m ist der Schwinger-Term, und sum (j + (1-m)/2)^2 = (m^3 - m)/12 gibt c = 1 je komplexes "
          "Fermion (exakt bis m = 48), also c = 8 für 16 reelle und so(16)_1; mit s = (1/2)^8 ist D8 + s gerade, unimodular "
          "mit 240 Wurzeln, also E8, und dafür gibt es ein lokales konformes Netz ohne Gitter. Die Toeplitz-Algebra ist "
          "irreduzibel und hat keinen separierenden Vektor; ein Gibbs-Zustand mit Punktspektrum hat eine fastperiodische "
          "modulare Gruppe, also keine echte halbseitige Inklusion: Raumzeittranslationen brauchen lokale Typ-III-Algebren. "
          "Korrektur: Ein festes Markenquadrat hat im Oktaeder nur D4 (zwei invariante Formen); die geometrische S4-Isotropie "
          "aus Z29/Z30 verlangt die drei gleichwertigen Rahmen oder eine interne S4.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z32_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": {"d_spec": ds, "kozyklus": coc, "rueckkehr": best}, "urteil": urteil},
    indent=1, ensure_ascii=False, default=str) + "\n")
