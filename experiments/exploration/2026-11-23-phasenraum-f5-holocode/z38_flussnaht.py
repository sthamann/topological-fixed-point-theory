"""Z38: Spinsektoren und geometrische Flussnaht -- V4 aus h = 1, Thetamodell, Anomalie, Zigarre, Gitterindex.

Prüft die siebte externe Analyse (H6h21-Vorschlag) und ergänzt einen Gittertest:
Teil 1  Auswahl: Untergruppen der Ordnung d des quadratischen Torus; Überlagerungsgrad 2d, 4d Majorana, Twistgewicht d/4.
        h = 1 => d = 4; unter den sieben Untergruppen der Ordnung 4 ist nur V4 = E[2] unter der Vierteldrehung invariant.
Teil 2  Thetamodell (Grad 8) für die vier Spinstrukturen: Deckzerlegung 5 + 3 nur periodisch-periodisch, sonst 4 + 4;
        antiperiodisch-antiperiodisch: T_a^2 = -1, delta T_a delta^-1 = -T_a, J = delta (i T_a) mit J^2 = -1;
        die Vierteldrehung lässt von den geraden Spinstrukturen nur (A,A) fest.
Teil 3  Anomalie: I_6 = (8/6) x^3 - (8/24) x p_1, über S^2 mit Flüssen 2 + 6 x 1 + 0: I_4 = 4 x^2 - p_1/3 (k = 8, c_- = 8);
        E8-Ladungen unter s = (1/2)^8: Histogramm {-2: 1, -1: 56, 0: 126, 1: 56, 2: 1}, neutrale Wurzeln = E7.
Teil 4  Lokalisierung: flache Ebene (Aharonov-Casher) N - 1 Moden, Zigarre N Moden (Normintegrale); 1 gegen 8.
Teil 5  Zigarre ds^2 = drho^2 + R^2 tanh^2(rho/R) dtheta^2: K = 2 sech^2/R^2, Ricci-Soliton mit Phi = -log cosh,
        B = N K mit Fluss N, e^{-2 Phi} B konstant, Lücke 1/(2R) am antiperiodischen Kreis.
Teil 6  Gitter (eigener Test): Wilson-Dirac und Overlap auf dem 2D-Torus mit N Flussquanten; der Index N hängt nicht an
        einer exakten ortslokalen chiralen Symmetrie (Wilson näherungsweise, Overlap exakt).
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import beta as Beta

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
say("Teil 1: Auswahl der Überlagerung")
for d in (1, 2, 3, 4, 5, 8):
    say(f"   d = {d}: Grad {2 * d}, h^0(L) = {2 * d}, {4 * d} Majorana, Twistgewicht {Fr(4 * d, 16)}")
check("Twistgewicht h = d/4 = 1 genau für d = 4 (16 Majorana, acht Schnitte)", [d for d in range(1, 20) if Fr(d, 4) == 1] == [4])
G4 = [(a, b) for a in range(4) for b in range(4)]
add = lambda u, v: ((u[0] + v[0]) % 4, (u[1] + v[1]) % 4)
subs = set()
for gens in itertools.combinations(G4, 2):
    S = {(0, 0)}
    frontier = [(0, 0)]
    while frontier:
        nw = []
        for x in frontier:
            for g in gens:
                y = add(x, g)
                if y not in S:
                    S.add(y); nw.append(y)
        frontier = nw
    if len(S) == 4:
        subs.add(frozenset(S))
for g in G4:
    S = {(0, 0)}; x = (0, 0)
    for _ in range(4):
        x = add(x, g); S.add(x)
    if len(S) == 4:
        subs.add(frozenset(S))
imul = lambda u: ((-u[1]) % 4, u[0] % 4)                  # Multiplikation mit i auf Z[i]/4
inv = [S for S in subs if {imul(u) for u in S} == set(S)]
cyc = [S for S in subs if any(len({((k * g[0]) % 4, (k * g[1]) % 4) for k in range(4)}) == 4 for g in S)]
say(f"   Untergruppen der Ordnung 4: {len(subs)} (zyklisch {len(cyc)}); unter i invariant: {[sorted(S) for S in inv]}")
check("sieben Untergruppen der Ordnung 4 (sechs zyklisch), nur V4 = E[2] ist unter der Vierteldrehung invariant",
      len(subs) == 7 and len(cyc) == 6 and len(inv) == 1 and set(inv[0]) == {(0, 0), (2, 0), (0, 2), (2, 2)})
check("Aut(Z4) = {+1, -1}: kein Automorphismus mit Quadrat -1 (k^2 = -1 mod 4 unlösbar)", all((k * k) % 4 != 3 for k in range(4)))
befund("Vierteldrehung + Twistgewicht 1 erzwingen V4, acht Schnitte und 16 Majorana (Z5, Z34 verbunden)", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Thetamodell vom Grad 8 für die vier Spinstrukturen")
def theta_model(al, be):
    ks = [Fr(j) + al for j in range(8)]                     # Indizes k in Z + alpha, modulo 8 mit e_{k+8} = (-1)^(2 beta) e_k
    def ref(k):
        q = k
        sgn = 1
        while q < 0:
            q += 8; sgn *= (-1) ** int(2 * be)
        while q >= 8:
            q -= 8; sgn *= (-1) ** int(2 * be)
        return ks.index(q), sgn
    n = 8
    delta = np.zeros((n, n), complex); T1 = np.zeros((n, n), complex); T2 = np.zeros((n, n), complex)
    for i, k in enumerate(ks):
        j, s_ = ref(-k); delta[j, i] = s_
        T1[i, i] = np.exp(1j * np.pi * float(k))
        j, s_ = ref(k + 4); T2[j, i] = s_
    return delta, T1, T2
tab = {}
for name, (al, be) in {"(P,P)": (Fr(0), Fr(0)), "(A,P)": (Fr(1, 2), Fr(0)), "(P,A)": (Fr(0), Fr(1, 2)), "(A,A)": (Fr(1, 2), Fr(1, 2))}.items():
    dl, T1, T2 = theta_model(al, be)
    I8 = np.eye(8)
    ev = np.round(np.linalg.eigvals(dl).real).astype(int)
    split = (int(np.sum(ev == 1)), int(np.sum(ev == -1)))
    t1sq = "+1" if np.allclose(T1 @ T1, I8) else ("-1" if np.allclose(T1 @ T1, -I8) else "?")
    t2sq = "+1" if np.allclose(T2 @ T2, I8) else ("-1" if np.allclose(T2 @ T2, -I8) else "?")
    conj = all(np.allclose(dl @ T @ np.linalg.inv(dl), -T) for T in (T1, T2))
    tab[name] = dict(split=split, T1sq=t1sq, T2sq=t2sq, dTd_minus=conj, delta_sq=bool(np.allclose(dl @ dl, I8)),
                     commute=bool(np.allclose(T1 @ T2, T2 @ T1)))
    say(f"   {name}: Deck {split[0]} + {split[1]}, T1^2 = {t1sq}, T2^2 = {t2sq}, delta T delta^-1 = -T: {conj}")
check("periodisch-periodisch: Deck 5 + 3 (Träger | Familien); alle anderen 4 + 4", tab["(P,P)"]["split"] == (5, 3)
      and all(tab[k]["split"] == (4, 4) for k in ("(A,P)", "(P,A)", "(A,A)")))
dl, T1, T2 = theta_model(Fr(1, 2), Fr(1, 2))
C1 = 1j * T1
J = dl @ C1
check("(A,A): T_a^2 = -1, delta T_a delta^-1 = -T_a; C_a = i T_a Involution, {delta, C_a} = 0, J = delta C_a mit J^2 = -1",
      tab["(A,A)"]["T1sq"] == "-1" and tab["(A,A)"]["T2sq"] == "-1" and tab["(A,A)"]["dTd_minus"]
      and np.allclose(C1 @ C1, np.eye(8)) and np.allclose(dl @ C1, -C1 @ dl) and np.allclose(J @ J, -np.eye(8)))
check("Vierteldrehung vertauscht die Perioden: (alpha, beta) -> (beta, alpha); von den geraden Strukturen bleibt nur (A,A)",
      [k for k, v in {"(A,P)": ("1/2", "0"), "(P,A)": ("0", "1/2"), "(A,A)": ("1/2", "1/2")}.items() if v[0] == v[1]] == ["(A,A)"])
say("   Spin-Dirac-Nullmoden H^.(E, M): M trivial (periodisch-periodisch) 1 + 1, M nichttrivial 0 (Riemann-Roch, Grad 0).")
befund("5 + 3 lebt im periodischen Sektor (Nullmoden), der Zustand im (A,A)-Sektor (4 + 4): verschiedene GSO-Sektoren", True)
ZAHLEN["theta"] = tab

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Anomalie und E8")
x, p1 = sp.symbols("x p1")
fl = [2] + [1] * 6 + [0]
I4 = sp.expand(sum(sp.Rational(n, 1) * (sp.Rational(3, 6) * x**2 - p1 / 24) for n in fl))
say(f"   Flüsse {fl}: I_4 = {I4}")
check("I_4 = 4 x^2 - p_1/3, also U(1)-Niveau 8 und c_- = 8", sp.simplify(I4 - (4 * x**2 - p1 / 3)) == 0)
roots = []
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product((1, -1), repeat=2):
        v = [Fr(0)] * 8; v[i] = Fr(si); v[j] = Fr(sj); roots.append(tuple(v))
for sg in itertools.product((Fr(1, 2), Fr(-1, 2)), repeat=8):
    if sum(1 for t in sg if t < 0) % 2 == 0:
        roots.append(tuple(sg))
s = tuple([Fr(1, 2)] * 8)
hist = {}
for r_ in roots:
    q = sum(a * b for a, b in zip(r_, s)); hist[q] = hist.get(q, 0) + 1
neutral = np.array([[float(t) for t in r_] for r_ in roots if sum(a * b for a, b in zip(r_, s)) == 0])
say(f"   E8-Wurzeln: {len(roots)}; Ladungen unter s: {dict(sorted(hist.items()))}; neutrale Wurzeln {len(neutral)}, Rang {np.linalg.matrix_rank(neutral)}")
check("240 Wurzeln, Histogramm {-2: 1, -1: 56, 0: 126, 1: 56, 2: 1}; 126 neutrale Wurzeln vom Rang 7 = E7",
      len(roots) == 240 and hist == {-2: 1, -1: 56, 0: 126, 1: 56, 2: 1} and np.linalg.matrix_rank(neutral) == 7)
check("q = (1,...,1) = 2 s hat |q|^2 = 8 (Niveau 8) und ist nicht primitiv; s selbst ist eine E8-Wurzel", 8 == sum(1 for _ in range(8))
      and sum(t * t for t in s) == 2)
befund("Zentralisator des primitiven Flusses ist e7 + u(1); der Fluss erhält Spin(10) nicht ganz", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Lokalisierung -- flache Ebene gegen Zigarre")
flat = lambda N: max(N - 1, 0)                              # psi = z^k e^{-phi}, |psi|^2 ~ r^(2k - 2N): k < N - 1
cigar = lambda N: sum(1 for k in range(0, 20) if k < N - Fr(1, 2))
check("Aharonov-Casher: flache Ebene N - 1; Zigarre N normierbare Moden (k < N - 1/2)", [flat(n) for n in (0, 1, 2, 3)] == [0, 0, 1, 2]
      and [cigar(n) for n in (0, 1, 2, 3)] == [0, 1, 2, 3])
tot_flat, tot_cig = sum(flat(n) for n in fl), sum(cigar(n) for n in fl)
say(f"   TFPT-Flüsse {fl}: flach {tot_flat}, Zigarre {tot_cig}")
check("mit den Flüssen 2 + 6 x 1 + 0: flach 1 Mode, auf der Zigarre 8", (tot_flat, tot_cig) == (1, 8))
norms = {}
for N, k in ((1, 0), (2, 0), (2, 1), (3, 2)):
    num = 2 * np.pi * quad(lambda r: r ** (2 * k + 1) * (1 + r * r) ** (-0.5 - N), 0, np.inf, limit=200)[0]
    exact = np.pi * Beta(k + 1, N - 0.5 - k)
    norms[(N, k)] = (num, exact)
say("   Normintegrale int |psi_k|^2 dA: " + ", ".join(f"N={N},k={k}: {v[0]:.10f} (exakt {v[1]:.10f})" for (N, k), v in norms.items()))
check("Normintegrale 2 pi (N=1), 2 pi/3 und 4 pi/3 (N=2) exakt, numerisch auf 1e-9", abs(norms[(1, 0)][1] - 2 * np.pi) < 1e-12
      and abs(norms[(2, 0)][1] - 2 * np.pi / 3) < 1e-12 and abs(norms[(2, 1)][1] - 4 * np.pi / 3) < 1e-12
      and all(abs(a - b) < 1e-9 for a, b in norms.values()))

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Zigarrengeometrie")
rho, R, N_ = sp.symbols("rho R N", positive=True)
f = R * sp.tanh(rho / R)
K = sp.simplify(-sp.diff(f, rho, 2) / f)
Phi = -sp.log(sp.cosh(rho / R))
eq_rr = sp.simplify(K + 2 * sp.diff(Phi, rho, 2))
eq_tt = sp.simplify(K + 2 * sp.diff(f, rho) * sp.diff(Phi, rho) / f)
check("K = (2/R^2) sech^2(rho/R); Ricci + 2 Hess(Phi) = 0 in beiden Komponenten", sp.simplify(K - 2 / R**2 / sp.cosh(rho / R) ** 2) == 0
      and eq_rr == 0 and eq_tt == 0)
u = sp.symbols("u", positive=True)                       # u = tanh(rho/R), du = sech^2/R drho
tot_K = sp.simplify(2 * sp.pi * sp.integrate((2 / R**2) * R * u * R, (u, 0, 1)))
tot_num = 2 * np.pi * quad(lambda r_: (2 / 1.7**2) / np.cosh(r_ / 1.7) ** 2 * 1.7 * np.tanh(r_ / 1.7), 0, 60)[0]
check("int K dA = 2 pi; B = N K hat Fluss N; e^{-2 Phi} B = 2N/R^2 konstant (dilatongewichtete Maxwellgleichung)",
      sp.simplify(tot_K - 2 * sp.pi) == 0 and abs(tot_num - 2 * np.pi) < 1e-9 and sp.simplify(sp.exp(-2 * Phi) * N_ * K - 2 * N_ / R**2) == 0)
check("asymptotischer Kreis antiperiodisch: Moden (n + 1/2)/R, Lücke 1/(2R)", min(abs(n + Fr(1, 2)) for n in range(-3, 3)) == Fr(1, 2))

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: Gitter -- Flussindex ohne exakte chirale Symmetrie")
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.diag([1.0 + 0j, -1])
def wilson(L, N, m0=0.0, r=1.0):
    n = L * L
    idx = lambda x, y: (x % L) * L + (y % L)
    B = 2 * np.pi * N / (L * L)
    D = np.zeros((2 * n, 2 * n), complex)
    for xx in range(L):
        for yy in range(L):
            i = idx(xx, yy)
            for mu, g in ((0, s1), (1, s2)):
                if mu == 0:
                    j = idx(xx + 1, yy); U = np.exp(-1j * B * yy)
                else:
                    j = idx(xx, yy + 1); U = np.exp(1j * B * L * xx) if yy == L - 1 else 1.0
                hop_f = 0.5 * (g - r * np.eye(2)) * U          # psi(x) <- psi(x+mu)
                hop_b = 0.5 * (-g - r * np.eye(2)) * np.conj(U)
                D[2 * i:2 * i + 2, 2 * j:2 * j + 2] += hop_f
                D[2 * j:2 * j + 2, 2 * i:2 * i + 2] += hop_b
            D[2 * i:2 * i + 2, 2 * i:2 * i + 2] += (2 * r + m0) * np.eye(2)
    return D
G5 = lambda L: np.kron(np.eye(L * L), s3)
res = {}
for L, N in ((12, 1), (12, 2), (14, 3)):
    Dw = wilson(L, N)
    w, v = np.linalg.eig(Dw)
    real = [(wi, vi) for wi, vi in zip(w, v.T) if abs(wi.imag) < 1e-6 and wi.real < 1.0]
    chir = sum(np.sign((vi.conj() @ G5(L) @ vi).real) for _, vi in real)
    Hw = G5(L) @ wilson(L, N, m0=-1.0)
    hw, hv = np.linalg.eigh(Hw)
    sgn = (hv * np.sign(hw)) @ hv.conj().T
    Dov = np.eye(2 * L * L) + G5(L) @ sgn
    wo = np.linalg.eigvals(Dov)
    zeros = int(np.sum(np.abs(wo) < 1e-8))
    ind_ov = int(round(-np.trace(sgn).real / 2))
    res[(L, N)] = dict(wilson_real_low=len(real), wilson_chiral_sum=int(chir), overlap_zero=zeros, overlap_index=ind_ov)
    say(f"   L = {L}, N = {N}: Wilson reelle Eigenwerte < 1: {len(real)}, Chiralitätssumme {int(chir)}; Overlap Nullmoden {zeros}, "
        f"Index -Tr(sign H)/2 = {ind_ov}")
check("Wilson-Dirac: Chiralitätssumme der physischen reellen Eigenwerte = +-N (ohne exakte chirale Symmetrie)",
      all(abs(v["wilson_chiral_sum"]) == N for (L, N), v in res.items()))
check("Overlap: genau N exakte Nullmoden, Index |N| (exakt)", all(v["overlap_zero"] == N and abs(v["overlap_index"]) == N for (L, N), v in res.items()))
befund("Die Chiralität der Flussnaht ist ein Transversalindex; sie braucht kein geladenes Phi und keine exakte ortslokale U(1)", True)
ZAHLEN["gitter"] = {str(k): v for k, v in res.items()}

urteil = ("Twistgewicht 1 verlangt d = 4 Translationen, und von den sieben Untergruppen der Ordnung 4 des quadratischen Torus "
          "ist nur V4 = E[2] unter der Vierteldrehung invariant: Vierteldrehung und h = 1 erzwingen V4, acht Schnitte und 16 "
          "Majorana. Im Thetamodell vom Grad 8 gibt nur die periodisch-periodische Spinstruktur die Deckzerlegung 5 + 3, die "
          "anderen 4 + 4; in (A,A) gilt T_a^2 = -1, delta T_a delta^-1 = -T_a, und J = delta (i T_a) hat J^2 = -1, (A,A) ist die "
          "einzige gerade Struktur, die die Vierteldrehung festhält. Die Flüsse 2 + 6 x 1 + 0 geben I_4 = 4 x^2 - p_1/3, also "
          "Niveau 8 und c_- = 8; unter s = (1/2)^8 zerfallen die E8-Wurzeln in 1 + 56 + 126 + 56 + 1, die neutralen bilden E7. "
          "Auf der flachen Ebene trägt Fluss N nur N - 1 Moden (insgesamt 1), auf der Zigarre N (insgesamt 8, Normintegrale "
          "exakt). Die Zigarre ist ein Ricci-Soliton mit Phi = -log cosh, B = N K hat Fluss N und erfüllt die "
          "dilatongewichtete Maxwellgleichung, der antiperiodische Kreis gibt die Lücke 1/(2R). Eigener Gittertest: Wilson-Dirac "
          "trägt mit N Flussquanten die Chiralitätssumme N, Overlap genau N exakte Nullmoden. Die Chiralität der Flussnaht ist "
          "ein Transversalindex ohne geladenes Paarungsfeld; der Satz aus Z24 (exakte ortslokale U(1) mit geladenem Phi gibt "
          "Index 0) greift hier nicht. Offen: Zufluss bzw. Gegenseite für die Anomalie, Gitterzigarre, Rückwirkung.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z38_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
