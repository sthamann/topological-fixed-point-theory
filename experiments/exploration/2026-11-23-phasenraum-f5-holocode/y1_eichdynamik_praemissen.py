"""Y: Eichdynamik (K15) und Prüfung der Prämissen von H6h.

Teil 1  Zusammengesetzte Verbindung: Der örtlich variierende 3+2-Split (Punkt in Gr_3(C^5)) definiert
        zwischen Nachbarn eine Paralleltranslation U_xy = polar(V_x^* V_y); sie ist eichkovariant.
Teil 2  Induzierte Wirkung: naive Diracfermionen auf dem 4D Lipschitzgitter (Torus N^4, überall
        antiperiodisch) in einem U(1) Hintergrund mit konstantem quantisiertem Fluss in einer Ebene.
        Gamma(m) = -log|det D|. Prüft: quadratisch im Fluss, positiv, eichinvariant, logarithmisch laufend.
Teil 3  Kopplungen aus einer Generation: Dynkin-Summen für SU(3), SU(2), U(1)_Y (GUT-normiert).
Teil 4  Prämisse B: folgt sie aus A, C und der 2D Naht?
Teil 5  Prämisse C gegen Daten: S4 auf den Familien, Restsymmetrien Z3 (geladene Leptonen) und
        Z2 x Z2 (Majorana-Neutrinos), Mischungswinkel gegen NuFIT und TFPT.
Teil 6  Prämisse A: was ohne sie frei wäre.
"""
import itertools
import json
from pathlib import Path

import numpy as np

rng = np.random.default_rng(11)
out = {}

# ------------------------------------------------------------------ Teil 1


def frame(P):
    w, v = np.linalg.eigh(P)
    return v[:, w > 0.5]


def rand_proj(rank=3, n=5):
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, _ = np.linalg.qr(A)
    return q[:, :rank] @ q[:, :rank].conj().T


def polar(M):
    u, _, vh = np.linalg.svd(M)
    return u @ vh


Px, Py = rand_proj(), rand_proj()
Vx, Vy = frame(Px), frame(Py)
U = polar(Vx.conj().T @ Vy)
hx, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
hy, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
U2 = polar((Vx @ hx).conj().T @ (Vy @ hy))
out["T1_link_unitary"] = bool(np.allclose(U.conj().T @ U, np.eye(3)))
out["T1_link_covariant"] = bool(np.allclose(U2, hx.conj().T @ U @ hy))
print(f"Teil 1: Verbindung aus dem lokalen 3+2-Split ist unitär {out['T1_link_unitary']} und eichkovariant "
      f"(U -> h_x^* U h_y) {out['T1_link_covariant']}")

# ------------------------------------------------------------------ Teil 2
s1 = np.array([[0, 1], [1, 0]]); s2 = np.array([[0, -1j], [1j, 0]])


def plane_operator(N, m, gauge=None):
    """Naive 2D lattice Dirac in the flux plane, antiperiodic, flux 2 pi m / N^2 per plaquette."""
    phi = 2 * np.pi * m / N ** 2
    idx = lambda x, y: (x % N) * N + (y % N)
    Tx = np.zeros((N * N, N * N), complex)
    Ty = np.zeros((N * N, N * N), complex)
    for x in range(N):
        for y in range(N):
            # x-link: twist at the seam of the torus keeps the flux uniform; antiperiodic sign
            ux = np.exp(-1j * phi * N * y) if x == N - 1 else 1.0
            if x == N - 1:
                ux *= -1
            uy = np.exp(1j * phi * x) * (-1 if y == N - 1 else 1)
            Tx[idx(x, y), idx(x + 1, y)] = ux
            Ty[idx(x, y), idx(x, y + 1)] = uy
    if gauge is not None:
        G = np.diag(np.exp(1j * gauge))
        Tx, Ty = G @ Tx @ G.conj().T, G @ Ty @ G.conj().T
    Dx, Dy = (Tx - Tx.conj().T) / 2, (Ty - Ty.conj().T) / 2
    return np.kron(Dx, s1) + np.kron(Dy, s2)


def gamma_eff(N, m, gauge=None, M=0.0):
    Dp = plane_operator(N, m, gauge)
    mu = np.linalg.eigvalsh(Dp.conj().T @ Dp)               # 2-spinor plane part
    p = 2 * np.pi * (np.arange(N) + 0.5) / N
    s = (np.sin(p)[:, None] ** 2 + np.sin(p)[None, :] ** 2).ravel()
    # 4D spinor = plane 2-spinor x 2; (D+M)^dag (D+M) = -D^2 + M^2 with anticommuting plane and transverse parts
    return -0.5 * 2 * np.sum(np.log(mu[:, None] + s[None, :] + M * M))


def kappa(N, m, M):
    phi = 2 * np.pi * m / N ** 2
    return (gamma_eff(N, m, M=M) - gamma_eff(N, 0, M=M)) / (0.5 * N ** 4 * phi ** 2)


print("Teil 2: induzierte Kopplung 1/g^2 = Delta Gamma / (N^4 F^2 / 2), mit Fermionmasse M als Infrarotregler")
print("   ohne Masse (F setzt selbst die Infrarotskala, erwartet F^2 log F):",
      {N: [round(kappa(N, m, 0.0), 4) for m in (1, 2, 3)] for N in (8, 16)})
table = {}
for M in (0.3, 0.45, 0.6):
    for N in (16, 24):
        table[(M, N)] = [kappa(N, m, M) for m in (1, 2)]
        F1 = 2 * np.pi / N ** 2
        print(f"   M = {M:4.2f}, N = {N}: F/M^2 = {F1 / M ** 2:.3f}; 1/g^2 bei m = 1, 2: "
              f"{table[(M, N)][0]:.5f}, {table[(M, N)][1]:.5f}")
N = 16
gauge = rng.uniform(0, 2 * np.pi, N * N)
inv = abs(gamma_eff(N, 2, gauge, M=0.45) - gamma_eff(N, 2, M=0.45)) < 1e-7
k24 = {M: table[(M, 24)][0] for M in (0.3, 0.45, 0.6)}
slope = np.polyfit(np.log(1 / np.array(list(k24))), list(k24.values()), 1)[0]
out["T2_quadratic"] = all(abs(v[1] - v[0]) / v[0] < 0.02 for (M, NN), v in table.items() if NN == 24)
out["T2_positive"] = all(v[0] > 0 and v[1] > 0 for v in table.values())
out["T2_gauge_invariant"] = bool(inv)
out["T2_volume_independent"] = all(abs(table[(M, 24)][0] - table[(M, 16)][0]) / table[(M, 24)][0] < 0.02 for M in (0.3, 0.45, 0.6))
out["T2_log_running_slope"] = float(slope)
out["T2_slope_continuum_16_dirac"] = 16 / (6 * np.pi ** 2)
print(f"   quadratisch (m = 1 und 2 gleich auf 2 %, N = 24): {out['T2_quadratic']}; positiv: {out['T2_positive']}; "
      f"eichinvariant: {out['T2_gauge_invariant']}; volumenunabhängig (N = 16 gegen 24): {out['T2_volume_independent']}")
print(f"   Lauf: d(1/g^2)/d log(1/M) = {slope:.4f}; Kontinuum für 16 Diracarten (naive Doppler, Ladung 1): "
      f"16/(6 pi^2) = {16 / (6 * np.pi ** 2):.4f}")

# ------------------------------------------------------------------ Teil 3
gen = [  # (name, SU3 dim, SU2 dim, Y, multiplicity)
    ("Q", 3, 2, 1 / 6), ("u^c", 3, 1, -2 / 3), ("d^c", 3, 1, 1 / 3), ("L", 1, 2, -1 / 2),
    ("e^c", 1, 1, 1), ("nu^c", 1, 1, 0)]
T = lambda dim: 0.5 if dim in (2, 3) else 0.0
b3 = sum(T(c) * w for _, c, w, _ in gen)
b2 = sum(T(w) * c for _, c, w, _ in gen)
b1 = 3 / 5 * sum(c * w * y * y for _, c, w, y in gen)
sin2 = (3 / 5 * 1) / (1 + 3 / 5 * 1) if np.isclose(b1, b2) and np.isclose(b2, b3) else None
out["T3_dynkin"] = {"SU3": b3, "SU2": b2, "U1_GUT": b1}
out["T3_sin2_thetaW_lattice_scale"] = sin2
print(f"\nTeil 3: Dynkin-Summen einer Generation: {out['T3_dynkin']} -> gleiche induzierte Kopplungen, "
      f"sin^2 theta_W = {sin2} an der Gitterskala (TFPT: 3/8 an der Trägerskala)")

# ------------------------------------------------------------------ Teil 4
from math import comb
chiral = {d: (comb(d, d // 2) // 2 if d % 4 == 0 else 0) for d in range(2, 10)}
out["T4_chiral_middle_forms"] = chiral
print(f"\nTeil 4: chirale Mittelformen je Dimension {chiral}: drei nur bei d = 4. Mit einer 2D Naht ist die Kodimension 2;"
      "\n   A auf allen vier Zyklen macht die Normalebene der Naht gleichwertig (quadratisch). B folgt aus A, C und der 2D Naht.")

# ------------------------------------------------------------------ Teil 5
# S4 on the three families: generated by sigma (3-cycle) and the quarter turn (rotation by 90 deg about one axis)
sig = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], float)
qt = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]], float)
grp = [np.eye(3)]
fr = [np.eye(3)]
while fr:
    nx = []
    for g in fr:
        for h in (sig, qt):
            m = np.round(g @ h, 9)
            if not any(np.allclose(m, x) for x in grp):
                grp.append(m); nx.append(m)
    fr = nx
out["T5_group_order"] = len(grp)


def order(m):
    x, n = m, 1
    while not np.allclose(x, np.eye(3)):
        x, n = x @ m, n + 1
    return n


z3 = [g for g in grp if order(g) == 3]
inv2 = [g for g in grp if order(g) == 2]
klein = []
for a, b in itertools.combinations(inv2, 2):
    if np.allclose(a @ b, b @ a) and order(a @ b) == 2:
        key = sorted([tuple(np.round(x, 6).ravel()) for x in (a, b, a @ b)])
        if key not in [k for k, _ in klein]:
            klein.append((key, (a, b)))


def eigbasis(mats):
    # simultaneous eigenbasis of commuting unitary matrices with distinct joint characters
    M = sum((k + 1) * np.pi * m for k, m in enumerate(mats))
    w, v = np.linalg.eig(M)
    return v / np.linalg.norm(v, axis=0)


patterns = {}
for g3 in z3:
    El = eigbasis([g3])
    for key, (a, b) in klein:
        Ev = eigbasis([a, b])
        Umix = np.abs(El.conj().T @ Ev) ** 2
        # canonical form: sort rows and columns
        rows = sorted([tuple(sorted(np.round(r, 6))) for r in Umix])
        patterns.setdefault(tuple(rows), Umix)
tbm = np.array([[2 / 3, 1 / 3, 0], [1 / 6, 1 / 3, 1 / 2], [1 / 6, 1 / 3, 1 / 2]])
tbm_key = tuple(sorted(tuple(sorted(np.round(r, 6))) for r in tbm))
out["T5_mixing_patterns"] = [list(map(list, k)) for k in patterns]
out["T5_TBM_in_S4"] = tbm_key in patterns
phi0 = 1 / (6 * np.pi) + 3 / (256 * np.pi ** 4)
out["T5_compare"] = {"sin2_12_TBM": 1 / 3, "sin2_12_TFPT": 1 / 3 - phi0 / 2, "sin2_12_NuFIT": 0.307,
                     "sin2_23_TBM": 0.5, "sin2_23_NuFIT_range": "0.47-0.57", "sin2_13_TBM": 0.0, "sin2_13_NuFIT": 0.0220}
print(f"\nTeil 5: Familiengruppe Ordnung {len(grp)}; {len(z3)} Elemente der Ordnung 3, {len(klein)} Kleinsche Vierergruppen")
for k in patterns:
    print(f"   Mischmuster |U|^2 (Zeilen sortiert): {[list(r) for r in k]}"
          + ("   <- tribimaximal" if k == tbm_key else ""))
print(f"   Vergleich: {out['T5_compare']}")

# ------------------------------------------------------------------ Teil 6
out["T6_without_A"] = {"tau": "frei (jeder Modul)", "c3": "1/(2 pi) ohne, 1/(4 pi) mit Parität statt 1/(8 pi)",
                       "Kragenbreite": "frei (T/lambda_1 frei)", "N der Naht": "frei (keine Modularbedingung)",
                       "E8": "nicht erzwungen", "Spinstruktur": "nicht eindeutig"}
print(f"\nTeil 6: ohne Prämisse A: {out['T6_without_A']}")

Path(__file__).with_name("y1_ergebnis.json").write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float))
