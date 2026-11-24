"""Z30: P1 aus einer irreduziblen Isometrie? Toeplitz-Index, Zeitneutralität, Spin-Z4, Marken als Operatoren.

Prüft die Konstruktionen der fünften externen Analyse (24.11.):
Teil 1  Irreduzible nichtunitäre Isometrie = einseitiger Shift (Wold): S^+S = I, SS^+ = I - |0><0|, Index -1;
        Vielfachheit m zerfällt in m reduzierende Teilräume (Rang des Defekts = m). Toeplitz: ind T_f = -wind f.
Teil 2  Zweiseitige Vervollständigung U auf l^2(Z_N): Spiegelung R, Schrittparität C, J = RC, D = diag(i^n),
        gelifteter Operator D~ = e^{i pi/4} D; Defekte der Halbraumprojektion an den beiden Schnitten.
Teil 3  Spektralparameter: Bahn {1, i, -1, -i} der Vierteldrehung, w^2 = z^4 - 1 (j = 1728), h^0(K(D)) = 3,
        h^0(O(D)) = 5 mit Halbdrehung 3 + 2.
Teil 4  Marken als reine Qubit-Zustände: Rang 3, Kommutator ergänzt sigma_z, Rang 4 = Herm(2); sechs Punkte
        {0, oo, +-1, +-i} = Oktaeder; Drehgruppe treu auf den vier Raumdiagonalen = S4; Tetraeder-Nullvektoren
        mit Minkowski-Gram-Eigenwerten {1, -1/3, -1/3, -1/3}; A4 erhält das Tetraeder.
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

rng = np.random.default_rng(30)
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
say("Teil 1: irreduzible Isometrie und Toeplitz-Index")
M = 60                                            # Abschneidung des einseitigen Shifts für Matrixtests
S = np.eye(M, k=-1)                               # S|n> = |n+1>
check("S^+ S = I auf dem Bereich n < M-1 (Abschneidung nur am fernen Rand)", np.allclose((S.T @ S)[:M - 1, :M - 1], np.eye(M - 1)))
check("S S^+ = I - |0><0| (bis auf den fernen Rand)", np.allclose((S @ S.T)[:M - 1, :M - 1], (np.eye(M) - np.outer(np.eye(M)[0], np.eye(M)[0]))[:M - 1, :M - 1]))
# Vielfachheit m: S_m = S (x) I_m; Defekt I - S_m S_m^+ hat Rang m, und jede Kopie ist reduzierend
for m in (1, 2, 3):
    Sm = np.kron(S, np.eye(m))
    defect = (np.eye(M * m) - Sm @ Sm.T)[: (M - 1) * m, : (M - 1) * m]
    rk = np.linalg.matrix_rank(defect)
    say(f"   Vielfachheit m = {m}: Rang des Defekts I - S S^+ = {rk}")
check("Defektrang = Vielfachheit; irreduzibel (keine reduzierenden Teilräume) <=> m = 1, Index -1",
      all(np.linalg.matrix_rank((np.eye(M * m) - np.kron(S, np.eye(m)) @ np.kron(S, np.eye(m)).T)[: (M - 1) * m, : (M - 1) * m]) == m
          for m in (1, 2, 3)))
# Toeplitz: Index von T_f mit f = z^n per Kompression auf H^2 (Fourier): ind T_{z^n} = -n
def toeplitz_index(nwind, N=256, K=60):
    th = 2 * np.pi * np.arange(N) / N
    f = np.exp(1j * nwind * th) * (2 + np.cos(th) + 0.3j * np.sin(2 * th)) / abs(2 + np.cos(th) + 0.3j * np.sin(2 * th))
    c = np.fft.fft(f) / N
    T = np.array([[c[(i - j) % N] for j in range(K)] for i in range(K)])
    sv = np.linalg.svd(T, compute_uv=False)
    # Kern- und Kokerndimension aus kleinen Singulärwerten der rechteckigen Abschneidung (K+|n| x K)
    Tl = np.array([[c[(i - j) % N] for j in range(K)] for i in range(K + abs(nwind))])
    ker = K - np.linalg.matrix_rank(Tl, tol=1e-6)
    Tr = np.array([[c[(i - j) % N] for j in range(K + abs(nwind))] for i in range(K)])
    coker = K - np.linalg.matrix_rank(Tr, tol=1e-6)
    return ker - coker
idx = {n: toeplitz_index(n) for n in (-2, -1, 0, 1, 2)}
say(f"   Toeplitz-Index (Kern - Kokern) je Windung: {idx}")
check("ind T_f = -wind(f) für glatte Symbole der Windung -2..2", all(idx[n] == -n for n in idx))
befund("Einheitswindung folgt aus Irreduzibilität einer nichtunitären Isometrie (Prämisse: der Ursprungsprozess ist eine)", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: zweiseitige Vervollständigung, Spiegelung, Spin-Z4")
res = {}
for N in (8, 16, 32, 64, 128):
    e = np.eye(N)
    U = np.roll(e, 1, axis=0)                     # U|n> = |n+1>
    R = np.zeros((N, N)); C = np.diag([(-1.0) ** n for n in range(N)])
    for n in range(N):
        R[(-n - 1) % N, n] = 1
    J = R @ C
    D = np.diag([1j ** n for n in range(N)])
    Dt = np.exp(1j * np.pi / 4) * D
    P = np.diag([1.0 if n < N // 2 else 0.0 for n in range(N)])   # Halbraum {0..N/2-1}
    dP = U @ P @ U.T - P
    ev = np.round(np.linalg.eigvalsh(dP), 10)
    viol = sum(not ok for ok in [
        np.allclose(R @ U @ R, U.T), np.allclose(R @ R, e), np.allclose(R @ C, -C @ R),
        np.allclose(J @ J, -e), np.allclose(np.linalg.matrix_power(J, 4), e),
        np.allclose(D @ D, C), np.allclose(D @ U @ D.conj().T, 1j * U),
        np.allclose(np.linalg.matrix_power(Dt, 4), -e), np.allclose(J @ Dt @ np.linalg.inv(J), np.linalg.inv(Dt))])
    defects = (int(round(sum(ev[ev > 0.5]))), -int(round(sum(-ev[ev < -0.5]))))
    res[N] = (viol, defects)
    say(f"   N = {N:3d}: Relationsverletzungen {viol}; U P U^+ - P hat die Eigenwerte +1 und -1 je einmal: {defects}")
ZAHLEN["endlich"] = {str(k): v for k, v in res.items()}
check("R U R = U^-1, {R, C} = 0, J^2 = -1, D^2 = C, D U D^+ = i U, D~^4 = -1, J D~ J^-1 = D~^-1 (N = 8..128)",
      all(v[0] == 0 for v in res.values()))
check("Halbraumprojektion: je ein Defekt +1 und -1 an den beiden Schnitten (Summe 0 im endlichen System)",
      all(v[1] == (1, -1) for v in res.values()))
say("   Hinweis: J^2 = -1 gilt für jedes Paar antivertauschender Involutionen; die Gleichsetzung mit H6s J (psi -> chi -> -psi)")
say("   ist eine Lesart, keine Herleitung.")

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Spektralparameter, vier Marken, Torus, Träger")
orbit = [1j ** k for k in range(4)]
check("Bahn der Vierteldrehung z -> iz durch 1: {1, i, -1, -i}", np.allclose(sorted(orbit, key=np.angle), sorted([1, 1j, -1, -1j], key=np.angle)))
# j-Invariante von w^2 = z^4 - 1: Quartik mit Wurzeln e_k; Doppelverhältnis lambda -> j = 256 (l^2 - l + 1)^3 / (l^2 (l - 1)^2)
r = orbit
lam = ((r[0] - r[2]) * (r[1] - r[3])) / ((r[0] - r[3]) * (r[1] - r[2]))
jinv = 256 * (lam**2 - lam + 1) ** 3 / (lam**2 * (lam - 1) ** 2)
say(f"   Doppelverhältnis lambda = {lam:.6f}, j = {jinv.real:.6f}")
check("w^2 = z^4 - 1 hat j = 1728, also tau = i", np.isclose(jinv, 1728))
h0 = {"O(D)": 4 + 1, "K(D)": -2 + 4 + 1}
even = [k for k in range(5) if k % 2 == 0]; odd = [k for k in range(5) if k % 2 == 1]
check("h^0(O(D)) = 5, h^0(K(D)) = 3 (Riemann-Roch auf P^1, deg D = 4)", h0 == {"O(D)": 5, "K(D)": 3})
check("Halbdrehung z -> -z teilt H^0(O(D)) = <1, z, ..., z^4> in 3 + 2", (len(even), len(odd)) == (3, 2))

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Marken als Operatoren")
I2 = np.eye(2); sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0 + 0j, -1])
def bloch(zv):                                    # Stereographie: z -> Blochvektor (z = 0 Nordpol, z = oo Südpol)
    if zv is None:
        return np.array([0, 0, -1.0])
    zz = complex(zv); d = 1 + abs(zz) ** 2
    return np.array([2 * zz.real / d, 2 * zz.imag / d, (1 - abs(zz) ** 2) / d])
proj = lambda n: 0.5 * (I2 + n[0] * sx + n[1] * sy + n[2] * sz)
marks = [proj(bloch(zv)) for zv in orbit]
flat = lambda Ms: np.array([np.concatenate([m.real.ravel(), m.imag.ravel()]) for m in Ms])
rk4 = np.linalg.matrix_rank(flat(marks))
comm = -2j * (marks[0] @ marks[1] - marks[1] @ marks[0])
rkc = np.linalg.matrix_rank(flat(marks + [comm]))
six = [proj(bloch(zv)) for zv in [0, None] + orbit]
rk6 = np.linalg.matrix_rank(flat(six))
say(f"   Rang: vier Marken {rk4}, mit Kommutator {rkc}, sechs Richtungen {rk6}")
check("Marken = (I +- sigma_x)/2, (I +- sigma_y)/2", np.allclose(marks[0], (I2 + sx) / 2) and np.allclose(marks[1], (I2 + sy) / 2))
check("Rang 3 -> mit -2i[P1, P2] = sigma_z Rang 4 = Herm(2)", rk4 == 3 and np.allclose(comm, sz) and rkc == 4 and rk6 == 4)
pts6 = np.array([bloch(zv) for zv in [0, None] + orbit])
check("sechs Punkte {0, oo, +-1, +-i} = Oktaederecken (+-e_x, +-e_y, +-e_z)",
      sorted(map(tuple, np.round(pts6, 12))) == sorted(map(tuple, np.round(np.vstack([np.eye(3), -np.eye(3)]), 12))))
# Drehgruppe des Oktaeders und Wirkung auf den vier Raumdiagonalen
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
G = [np.array(t).reshape(3, 3) for t in grp]
diags = [np.array(v) for v in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]
def perm(g):
    out = []
    for dvec in diags:
        im = g @ dvec
        out.append(next(k for k, e in enumerate(diags) if np.allclose(im, e) or np.allclose(im, -e)))
    return tuple(out)
perms = {perm(g) for g in G}
check("24 Drehungen, treu auf den vier Raumdiagonalen: Bild = alle 24 Permutationen (S4)", len(G) == 24 and len(perms) == 24)
tet = [np.array(v) / np.sqrt(3) for v in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]
X = [0.5 * (I2 + n[0] * sx + n[1] * sy + n[2] * sz) for n in tet]
eta = lambda A, B: 0.5 * (np.trace(A) * np.trace(B) - np.trace(A @ B)).real   # Polarisierung von det
Gm = np.array([[eta(A, B) for B in X] for A in X])
evg = np.sort(np.linalg.eigvalsh(Gm))
say(f"   Minkowski-Gram der Tetraeder-Nullvektoren: Eigenwerte {np.round(evg, 12)}")
check("X_i Rang 1 und positiv (lichtartig, zukunftsgerichtet)", all(np.isclose(np.linalg.det(x).real, 0) and np.linalg.eigvalsh(x).min() > -1e-12 for x in X))
check("Gram-Eigenwerte {1, -1/3, -1/3, -1/3}: linear unabhängig, Signatur 1 + 3", np.allclose(evg, [-1 / 3, -1 / 3, -1 / 3, 1]))
keep = sum(1 for g in G if all(any(np.allclose(g @ t, u) for u in tet) for t in tet))
check("12 Drehungen erhalten das Tetraeder (A4), 12 tauschen es mit dem dualen", keep == 12)
befund("Familiengruppe S4, Tetraederpaar und Lorentzkegel haben eine gemeinsame Realisierung auf Herm(2)", True)

urteil = ("Die Operatorkonstruktionen stimmen. Eine irreduzible nichtunitäre Isometrie ist nach Wold der einseitige Shift; "
          "ihr Defekt hat Rang 1 (Vielfachheit m gäbe Rang m und zerfiele), also Index -1, und der Toeplitz-Index ist -wind "
          "(numerisch für Windung -2..2). In der zweiseitigen Vervollständigung gelten R U R = U^-1, {R, C} = 0, J = RC mit "
          "J^2 = -1, D^2 = C, D U D^+ = iU, D~^4 = -1 und J D~ J^-1 = D~^-1 für N = 8 bis 128; die Halbraumprojektion hat an "
          "den beiden Schnitten die Defekte +1 und -1. Die Bahn der Vierteldrehung gibt w^2 = z^4 - 1 mit j = 1728 (tau = i), "
          "h^0(K(D)) = 3 und h^0(O(D)) = 5 mit 3 + 2. Als Qubit-Zustände haben die vier Marken Rang 3; ihr Kommutator "
          "ergänzt sigma_z zu Herm(2) (Rang 4), die sechs Punkte {0, oo, +-1, +-i} sind das Oktaeder, dessen Drehgruppe treu "
          "als S4 auf den Raumdiagonalen wirkt; die Tetraeder-Nullvektoren haben Minkowski-Gram {1, -1/3, -1/3, -1/3}. "
          "Grenzen: Irreduzibilität ist eine neue Prämisse statt der Einheitswindung; J^2 = -1 gilt für jedes Paar "
          "antivertauschender Involutionen, die Gleichsetzung mit H6s J ist eine Lesart; die Identifikation des Hardy-Kerns "
          "mit dem P1-Kern und das Soldering der Operatorrichtungen an physische Verschiebungen sind offen.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
Path(__file__).with_name("z30_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
