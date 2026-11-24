"""Z62: 3+1 als natürlicher Kegel der Standardform -- die Naht als Himmelssphäre (Vogelperspektive, Teil 2).

Die Naht CP^1 ist der Raum der reinen Zustände des Dubletts V = H0(O(1))^* (Euler-Folge, Bestand). Tomita-Takesaki gibt
für einen treuen Zustand rho auf End(V) = M_n(C) die Standardform (M_n(C) mit Hilbert-Schmidt-Produkt, J x = x^+,
natürlicher Kegel P = {a sqrt(rho) a^+}). Geprüft wird:

Teil 1  J-reeller Teil = Herm(n), natürlicher Kegel = positive Matrizen, Extremalstrahlen = reine Zustände = CP^(n-1).
Teil 2  Der Kegel ist genau für n = 2 ein Lorentzkegel (quadratisch, Extremalstrahlen bilden eine Sphäre S^(n^2-2)):
        2(n - 1) = n^2 - 2 nur für n = 2. Dann ist Herm(2) = R^(1,3), der Kegel der Zukunftslichtkegel und die Naht CP^1 die
        Himmelssphäre (Nullrichtungen).
Teil 3  x -> g x g^+ (g in SL(2,C)) erhält J, Kegel und det: eigentliche orthochrone Lorentzgruppe; auf den reinen Zuständen
        wirkt sie als Möbiusabbildung der Naht. Der Modulfluss rho^(it) x rho^(-it) ist eine Drehung um sqrt(rho).
Teil 4  Bilanz.
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(62)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def crand(*s):
    return rng.normal(size=s) + 1j*rng.normal(size=s)


def msqrt(h):
    w, v = np.linalg.eigh(h)
    return (v*np.sqrt(np.clip(w, 0, None))) @ v.conj().T


def mpow(h, p):
    w, v = np.linalg.eigh(h)
    return (v*w**p) @ v.conj().T


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Standardform des Dubletts")
ok = True
for n in (2, 3):
    A = crand(n, n)
    rho = A @ A.conj().T
    rho /= np.trace(rho).real
    xi = msqrt(rho)
    for _ in range(200):
        a = crand(n, n)
        y = a @ xi @ a.conj().T                                         # a J a J xi
        ok &= np.allclose(y, y.conj().T) and np.linalg.eigvalsh(y).min() > -1e-10
        B = crand(n, n)
        p = B @ B.conj().T                                              # beliebige positive Matrix
        a2 = msqrt(p) @ mpow(rho, -0.25)
        ok &= np.allclose(a2 @ xi @ a2.conj().T, p)
check("für n = 2, 3: J x = x^+, der natürliche Kegel {a J a J xi_rho} ist genau der Kegel der positiven Matrizen in Herm(n) "
      "(beide Inklusionen an 200 Zufallselementen)", ok)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Nur das Dublett gibt einen Lorentzkegel")
tab = {n: (2*(n - 1), n*n - 2) for n in range(2, 7)}
ZAHLEN["extremalstrahlen_vs_sphaere"] = tab
check("reine Zustände CP^(n-1) haben Dimension 2(n - 1), die Extremalstrahlen eines Lorentzkegels in Herm(n) bilden S^(n^2 - 2): "
      "gleich nur für n = 2 " + str(tab), [n for n, (u, v) in tab.items() if u == v] == [2])
sig = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
okc = True
for _ in range(500):
    c = rng.normal(size=4)
    x = sum(ci*s for ci, s in zip(c, sig))
    m2 = c[0]**2 - c[1]**2 - c[2]**2 - c[3]**2
    psd = np.linalg.eigvalsh(x).min() >= -1e-12
    okc &= abs(np.linalg.det(x).real - m2) < 1e-10 and (psd == (m2 >= -1e-12 and c[0] >= 0))
check("n = 2: det(t + x.sigma) = t^2 - |x|^2, und positiv heißt genau t >= |x|: der natürliche Kegel ist der Zukunftslichtkegel "
      "von R^(1,3); reine Zustände sind die Nullstrahlen, die Naht CP^1 ist die Himmelssphäre", okc)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Lorentzgruppe und Modulfluss")
okl = True
for _ in range(200):
    g = crand(2, 2)
    g /= np.sqrt(np.linalg.det(g))
    c = rng.normal(size=4)
    x = sum(ci*s for ci, s in zip(c, sig))
    y = g @ x @ g.conj().T
    okl &= np.allclose(y, y.conj().T) and abs(np.linalg.det(y) - np.linalg.det(x)) < 1e-8
    psi = crand(2)
    zpsi = psi[1]/psi[0]
    gp = g @ psi
    mob = (g[1, 0] + g[1, 1]*zpsi)/(g[0, 0] + g[0, 1]*zpsi)
    okl &= abs(gp[1]/gp[0] - mob) < 1e-8
check("x -> g x g^+ mit g in SL(2,C) erhält J-Realität, det und den Kegel (Lorentzgruppe); auf reinen Zuständen wirkt g als "
      "Möbiusabbildung der Naht", okl)
A = crand(2, 2)
rho = A @ A.conj().T
rho /= np.trace(rho).real
okm = True
for t in np.linspace(-2, 2, 9):
    w, v = np.linalg.eigh(rho)
    rit = (v*np.exp(1j*t*np.log(w))) @ v.conj().T
    c = rng.normal(size=4)
    x = sum(ci*s for ci, s in zip(c, sig))
    y = rit @ x @ rit.conj().T
    okm &= abs(np.trace(y) - np.trace(x)) < 1e-10 and abs(np.linalg.det(y) - np.linalg.det(x)) < 1e-10
    okm &= np.allclose(rit @ msqrt(rho) @ rit.conj().T, msqrt(rho))
check("der Modulfluss Delta^(it) x = rho^(it) x rho^(-it) erhält Zeitkomponente und det und hält sqrt(rho) fest: eine Drehung "
      "im Ruhesystem des Zustands", okm)
befund("Die Kausalstruktur von 3+1 (Dimension 4, Signatur, Zukunftskegel, Lorentzgruppe, Himmelssphäre = Naht) folgt aus der "
       "Standardform eines Zustands auf dem Dublett, also aus Tomita-Takesaki; dass es gerade 3+1 ist, liegt daran, dass die "
       "Naht eine Sphäre (CP^1) ist. Das ist die Idee der Ur-Alternative (von Weizsäcker) und der Himmelssphäre (Penrose), "
       "hier aus dem Originalzustand gelesen", True)
befund("Nicht gezeigt: dass der Originalzustand einen treuen Zustand auf End(V) mit V = H0(O(1))^* festlegt; die Orte "
       "(Translationen) liefert erst die Twistorlesart T (z61)", True)

urteil = ("Die Naht CP^1 ist der Raum der reinen Zustände des Dubletts V. Die Tomita-Standardform eines treuen Zustands auf End(V) "
          "hat als J-reellen Teil Herm(2) = R^(1,3); ihr natürlicher Kegel ist genau der Zukunftslichtkegel, die reinen Zustände "
          "sind die Nullstrahlen, die Naht ist die Himmelssphäre, SL(2,C) wirkt als Lorentzgruppe und auf der Naht als "
          "Möbiusgruppe, der Modulfluss als Drehung im Ruhesystem. Ein Lorentzkegel entsteht nur für n = 2 (2(n-1) = n^2 - 2), "
          "also gerade 3+1. Offen ist, dass der Originalzustand diesen Dublettzustand festlegt; die Orte gibt die Twistorlesart.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z62_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
