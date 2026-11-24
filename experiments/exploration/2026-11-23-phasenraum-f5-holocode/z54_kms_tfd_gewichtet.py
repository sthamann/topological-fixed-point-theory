"""Z54: Prüfung der externen Analyse zu H6h27 (Zweig research/h6h27-source-dynamics-20261124) -- Zustand und Operator.

Teil 1  KMS: Für einen Oszillator mit alpha_t(a) = e^{-it} a erzwingt die volle KMS-Bedingung p_(n+1) = e^{-beta} p_n, also den
        gaußschen Zustand (1-r) r^N. Gleiche mittlere Besetzung allein genügt nicht (Gegenbeispiel rho_alt). Präzisierung
        zu z46/z48: Der Mittelwertvergleich wählt q nur innerhalb der geometrischen Familie; KMS rekonstruiert sie.
Teil 2  Polarverschiebung: J_-1 wirkt auf der Besetzungsleiter als a^+ (Faktoren sqrt(n+1)), nicht als Isometrie; der
        isometrische Schritt ist V = a^+ (N+1)^(-1/2). Korrektur zu z48 ('J_-1 wirkt als S^4' gilt nur auf den
        Einteilchen-Marken).
Teil 3  Thermofeldoperator: d_L, d_R vernichten den Thermofeldzustand; H = d_L^+ d_L + d_R^+ d_R ist im Paarsektor
        proportional zu A^+ (N+1) A mit (A f)_n = q f_n - f_(n+1): diskrete Leiter Delta n statt Band, Randantwort
        (1-r) Sum r^n/(x + Delta n).
Teil 4  Gewichteter Stern: Kantengewicht W = n+1 in der Tiefe; der symmetrische Sektor ist genau diese Leiter, die übrigen
        Sektoren liegen dazwischen, unter jeder Schwelle bleiben endlich viele Moden; die Nullmode bleibt universell.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "research" / "originlab"))
from process_rules import metric_star  # noqa: E402

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ Teil 1
say("Teil 1: KMS statt Mittelwert")
nmax = 80
beta = 1.1
r = np.exp(-beta)
p = (1 - r)*r**np.arange(nmax)
eps = (1 - r)*r/4
alt = p.copy()
alt[0] += eps
alt[1] -= 2*eps
alt[2] += eps
n = np.arange(nmax)
kms = lambda pp: np.max(np.abs(pp[1:] - r*pp[:-1]))
wick = lambda pp: abs(np.sum(n*(n - 1)*pp) - 2*np.sum(n*pp)**2)          # <a+ a+ a a> = 2 <a+ a>^2
say(f"   thermisch: KMS-Defekt {kms(p):.1e}, Wick-Defekt {wick(p):.1e}; rho_alt: Spur {alt.sum():.12f}, min p {alt.min():.4f}, "
    f"Mittel {np.sum(n*alt):.12f} (thermisch {np.sum(n*p):.12f}), KMS-Defekt {kms(alt):.3f}, Wick-Defekt {wick(alt):.3f}")
check("die volle KMS-Bedingung (p_(n+1) = e^-beta p_n auf allen Matrixeinheiten) legt den gaußschen Zustand (1-r) r^N fest",
      kms(p) < 1e-15 and wick(p) < 1e-12)
check("rho_alt = rho_beta + eps (P0 - 2 P1 + P2) ist positiv, normiert, hat dieselbe mittlere Besetzung und verletzt KMS und Wick",
      abs(alt.sum() - 1) < 1e-12 and alt.min() > 0 and abs(np.sum(n*alt) - np.sum(n*p)) < 1e-12 and kms(alt) > 1e-3
      and wick(alt) > 1e-3)
befund("Präzisierung zu z46/z48: Der Mittelwertvergleich wählt q nur innerhalb der geometrischen Familie; die volle KMS-"
       "Bedingung des Stromoszillators rekonstruiert den Zustand und gibt dasselbe q = e^-beta/2 (bei beta = 4 pi: e^-2pi)", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Polarverschiebung")
N = 40
ad = np.diag(np.sqrt(np.arange(1, N)), -1)                                   # a^+ |n> = sqrt(n+1) |n+1>
Om = np.zeros(N)
Om[0] = 1
nrm = np.linalg.norm(ad @ ad @ Om)**2
V = ad @ np.diag(1/np.sqrt(np.arange(1, N + 1)))
iso = np.max(np.abs((V.T @ V)[:N - 1, :N - 1] - np.eye(N - 1)))
say(f"   |J_-1^2 Omega|^2 = {nrm:.6f}; V = J_-1 (N+1)^(-1/2): |V^+ V - I| = {iso:.1e} (ohne Abschneiderand)")
check("J_-1 ist auf der Besetzungsleiter keine Isometrie (|J_-1^2 Omega|^2 = 2); die Polarphase V ist der isometrische Schritt",
      abs(nrm - 2) < 1e-12 and iso < 1e-12)
befund("Korrektur zu z48: 'J_-1 wirkt als S^4' gilt für die Einteilchen-Marken (r -> r+1, k -> k+4); auf der "
       "Besetzungsleiter des Thermofeldzustands ist der isometrische Schritt V = J_-1 (N+1)^(-1/2)", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Thermofeldoperator")
q = 0.5
r = q*q
M = 30
a = np.diag(np.sqrt(np.arange(1, M)), 1)
I = np.eye(M)
aL, aR = np.kron(a, I), np.kron(I, a)
dL = (aL - q*aR.T)/np.sqrt(1 - r)
dR = (aR - q*aL.T)/np.sqrt(1 - r)
idx = lambda nl, nr: nl*M + nr
eta = np.zeros(M*M)
for k in range(M):
    eta[idx(k, k)] = np.sqrt(1 - r)*q**k
inner = [idx(k, k) for k in range(M - 3)]
check("d_L und d_R vernichten den Thermofeldzustand sqrt(1-q^2) Sum q^n |n,n> (bis auf den Abschneiderand)",
      np.max(np.abs((dL @ eta)[inner])) < 1e-12 and np.max(np.abs((dR @ eta)[inner])) < 1e-12)
Hq = dL.T @ dL + dR.T @ dR
pair = [idx(k, k) for k in range(M)]
Hp = Hq[np.ix_(pair, pair)]
off = np.max(np.abs(Hq[np.ix_(pair, [i for i in range(M*M) if i not in set(pair)])]))
Afl = np.zeros((M - 1, M))
for k in range(M - 1):
    Afl[k, k], Afl[k, k + 1] = q, -1
K = Afl.T @ np.diag(np.arange(1, M)) @ Afl
c = Hp[0, 0]/K[0, 0]
dev = np.max(np.abs(Hp[:M - 2, :M - 2] - c*K[:M - 2, :M - 2]))
say(f"   Paarsektor: Kopplung nach außen {off:.1e}; H_Paar = {c:.6f} A^+ (N+1) A (Abweichung {dev:.1e}); 2/(1-q^2) = {2/(1 - r):.6f}")
check("H = d_L^+ d_L + d_R^+ d_R erhält den Paarsektor und ist dort 2/(1-q^2) mal A^+ (N+1) A mit (A f)_n = q f_n - f_(n+1)",
      off < 1e-12 and dev < 1e-10 and abs(c - 2/(1 - r)) < 1e-12)
L = 200
Al = np.zeros((L - 1, L))
for k in range(L - 1):
    Al[k, k], Al[k, k + 1] = q, -1
lam = 1.0
KL = lam*Al.T @ np.diag(np.arange(1, L)) @ Al
ev = np.sort(np.linalg.eigvalsh(KL))[:8]
Delta = lam*(1 - r)
xs = np.array([0.2, 1.0, 3.0])
Cnum = [np.linalg.solve(x*np.eye(L) + KL, np.eye(L)[:, 0])[0] for x in xs]
Cform = [(1 - r)*sum(r**k/(x + Delta*k) for k in range(400)) for x in xs]
say(f"   K_lambda (lambda = 1, q = 1/2): tiefste Eigenwerte {np.round(ev, 8).tolist()}, Delta = {Delta}")
check("K_lambda = lambda A^+ (N+1) A hat das Spektrum Delta n mit Delta = lambda (1 - q^2): diskrete Leiter statt Band",
      np.allclose(ev, Delta*np.arange(8), atol=1e-9))
check("die Randantwort ist (1-r) Sum r^n/(x + Delta n): positives Spektralmaß ausdrücklich", np.allclose(Cnum, Cform, atol=1e-10))
befund("Der Zustand legt die Dynamik nicht fest: A^+ A (Band) und A^+ (N+1) A (Leiter) haben denselben geometrischen "
       "Grundzustand; erst 'quadratisch in den Strommoden' wählt die Leiter", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Gewichteter Stern")


def star_L(qv, m, weighted):
    st = metric_star(qv, m)
    W = np.array([n + 1 if weighted else 1 for j in range(4) for n in range(m)], float)
    return st, st["B"].T @ np.diag(W) @ st["B"]


q = 0.5
r = q*q
st, Lw = star_L(q, 60, True)
check("die Nullmode bleibt eta_v = M_v/M_eff (universell, eta_v/M_v = 1/M_eff)", np.max(np.abs(Lw @ st["eta"])) < 1e-12)
nodes = st["nodes"]
ix = {v: i for i, v in enumerate(nodes)}
sym = np.zeros((len(nodes), 61))
sym[ix["c"], 0] = 1
for n_ in range(60):
    for j in range(4):
        sym[ix[(j, n_)], n_ + 1] = 0.5
Ls = sym.T @ Lw @ sym
evs = np.sort(np.linalg.eigvalsh(Ls))[:6]
evall = np.sort(np.linalg.eigvalsh(Lw))
check("der symmetrische Sektor (Zentrum plus Armmittel) ist genau die Leiter Delta n, Delta = 1 - q^2",
      np.allclose(evs, (1 - r)*np.arange(6), atol=1e-8))
ladder = (1 - r)*np.arange(40)
others = [e for e in evall[:120] if np.min(np.abs(e - ladder)) > 1e-6]
inter = all(any(ladder[k] < e < ladder[k + 1] for k in range(39)) for e in others if e < ladder[-1])
say(f"   tiefste Nicht-Leiter-Eigenwerte: {np.round(others[:6], 4).tolist()}")
check("die übrigen Sektoren liegen zwischen benachbarten Leiterstufen (verschränkt)", inter and len(others) > 3)
cnt = {}
for m in (20, 40, 80):
    _, Lf = star_L(q, m, False)
    _, Lwm = star_L(q, m, True)
    cnt[m] = (int(np.sum(np.linalg.eigvalsh(Lf) < 2.0)), int(np.sum(np.linalg.eigvalsh(Lwm) < 2.0)))
say(f"   Moden unter der Schwelle 2 lambda (flach, gewichtet): {cnt}")
check("unter einer festen Schwelle wächst die Modenzahl beim flachen Stern mit der Tiefe, beim gewichteten nicht",
      cnt[80][0] > cnt[20][0] and cnt[80][1] == cnt[20][1])
befund("Der gewichtete Stern ist der natürliche Kandidat, wenn die Dynamik quadratisch in den Strommoden ist: diskretes "
       "Spektrum, endlich viele leichte Moden, gleiche universelle Nullmode", True)

urteil = ("Die externe Analyse hält in den nachgerechneten Punkten. Die volle KMS-Bedingung des Stromoszillators legt den "
          "gaußschen Zustand fest; gleiche mittlere Besetzung allein genügt nicht (rho_alt ist positiv, normiert, hat dasselbe "
          "Mittel und verletzt KMS und Wick). Korrektur zu z48: J_-1 ist auf der Besetzungsleiter keine Isometrie "
          "(|J_-1^2 Omega|^2 = 2), der isometrische Schritt ist V = J_-1 (N+1)^(-1/2). Der Thermofeldoperator d_L^+ d_L + "
          "d_R^+ d_R ist im Paarsektor 2/(1-q^2) A^+ (N+1) A: Spektrum Delta n, Randantwort (1-r) Sum r^n/(x + Delta n). Derselbe "
          "Zustand trägt also zwei Dynamiken (Band oder Leiter); quadratisch in den Strommoden wählt die Leiter. Im "
          "gewichteten Stern ist der symmetrische Sektor genau diese Leiter, die übrigen Sektoren liegen dazwischen, unter "
          "jeder Schwelle bleiben endlich viele Moden, und die Nullmode bleibt universell.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z54_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
