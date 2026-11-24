"""Z71: Prüfung der externen Analyse 'Weiterentwicklung von H6h33' (24.11.): lokale Modulzeit, 4D-Fortsetzung, Ursprung.

Teil 1  Lokale Modulzeit: Dieselbe originale Kovarianz, eingeschränkt auf einen Bogen der Naht, hat ein lokales
        Modulspektrum, das mit feinerer Auflösung dicht wird (kein pi-Gitter); global bleibt es 4 pi x.
Teil 2  4D-Fortsetzung: W(tau, 0) der S^3-Fortsetzung ist exakt die thermische Stromantwort (z58). Dieselbe Linienantwort hat
        auch die 2D-Fortsetzung cos(m theta) längs der Naht; die native Querantwort ist die 2D-Antwort längs des Nahtwinkels.
        Die vorgeschlagene Identität d^2_gamma W = -(1/3)(d^2_tau - 1) W gilt für S^3, nicht für den nativen Nahtwinkel (dort
        gilt d^2_theta W = -d^2_tau W). Die Linienantwort erzwingt die Dimension also nicht.
Teil 3  Ursprungskandidat: Zwei binäre Beobachtungen mit Überlapp p geben vier Marken mit harmonischem Kreuzverhältnis
        (tau = i) genau bei p = 1/2 (maximal inkompatibel, MUB). Die Paarung nach Beobachtung ist die Antipodenpaarung, also das
        kanonische eta aus z70; das gleich gewichtete Gemisch ist der Spurzustand (keine Modulzeit).
"""
import json
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).parent
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
say("Teil 1: Lokale Modulzeit auf einem Nahtbogen")
alpha = 0.25


def local_spectrum(M, frac=0.5):
    xs = np.arange(-M//2, M//2) + alpha
    c = 1/(1 + np.exp(np.clip(4*np.pi*xs, -700, 700)))
    th = 2*np.pi*np.arange(M)/M
    U = np.exp(1j*np.outer(th, xs))/np.sqrt(M)                  # unitäre DFT zwischen Moden und Gitterpunkten
    C = (U*c) @ U.conj().T
    glob = np.sort(np.linalg.eigvalsh(C))
    A = np.arange(int(frac*M))
    nu = np.linalg.eigvalsh(C[np.ix_(A, A)])
    nu = nu[(nu > 1e-13) & (nu < 1 - 1e-13)]
    eps = np.sort(np.log((1 - nu)/nu))
    return glob, c, eps


res = {}
for M in (128, 256, 512):
    glob, c, eps = local_spectrum(M)
    win = eps[np.abs(eps) < 4*np.pi]
    res[M] = (len(win), float(np.median(np.diff(win))) if len(win) > 2 else np.inf,
              float(np.mean(np.abs(win/np.pi - np.round(win/np.pi)) < 0.05)), np.allclose(np.sort(c), glob, atol=1e-10))
ZAHLEN["lokales_spektrum"] = {M: {"anzahl_in_|eps|<4pi": r[0], "median_abstand": r[1], "anteil_nahe_pi_Z": r[2]}
                              for M, r in res.items()}
say(f"   {ZAHLEN['lokales_spektrum']}")
check("global ist das Kovarianzspektrum exakt c(x), also Modulenergien 4 pi x (Abstand 4 pi) -- auf dem halben Bogen wächst "
      "die Zahl lokaler Modulenergien in |eps| < 4 pi mit der Auflösung (128, 256, 512 Punkte), der mittlere Abstand sinkt, und "
      "die Werte liegen nicht auf dem pi-Gitter: das lokale Spektrum wird dicht",
      all(r[3] for r in res.values()) and res[128][0] < res[256][0] < res[512][0]
      and res[512][1] < res[128][1] and max(r[2] for r in res.values()) < 0.3)
befund("Die externe Analyse hält in diesem Punkt: Die lokale Modulzeit desselben Zustands ist nicht periodisch; damit ist der "
       "spektrale Widerspruch zur periodischen globalen Uhr aufgelöst. Dass sie die physikalische Zeit ist, verlangt noch lokale "
       "Translationsbeziehungen (halbseitige Inklusionen)", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Die 4D-Fortsetzung und die native Querantwort")
B = 4*np.pi
ms = np.arange(1, 60)


def W4(tau, g):
    s = np.sin(ms*g)/np.sin(g) if g != 0 else ms.astype(float)
    return np.sum(s*(np.exp(-ms*tau) + np.exp(-ms*(B - tau)))/(1 - np.exp(-B*ms)))


def W2(tau, th):                       # native 2D-Stromantwort längs des Nahtwinkels (Realteil der holomorphen Antwort)
    return np.sum(ms*np.cos(ms*th)*(np.exp(-ms*tau) + np.exp(-ms*(B - tau)))/(1 - np.exp(-B*ms)))


def Gcur(tau):                         # z58: <J_m J_-m> e^{-m tau} + <J_-m J_m> e^{m tau}
    up = ms/(1 - np.exp(-B*ms))
    dn = ms/(np.exp(B*ms) - 1)
    return np.sum(up*np.exp(-ms*tau) + dn*np.exp(ms*tau))


taus = [0.7, 1.3, 2.9, 6.1]
ok_line = all(abs(W4(t, 0) - Gcur(t)) < 1e-10*abs(Gcur(t)) and abs(W2(t, 0) - Gcur(t)) < 1e-10*abs(Gcur(t)) for t in taus)
check("W4(tau, 0) = W2(tau, 0) = thermische Stromantwort aus z58 (vier tau-Werte): die Linienantwort unterscheidet S^3 und die "
      "Naht nicht", ok_line)
m, g, t = sp.symbols("m gamma tau", positive=True)
d2_s3 = sp.limit(sp.diff(sp.sin(m*g)/sp.sin(g), g, 2), g, 0)
d2_s1 = sp.limit(sp.diff(m*sp.cos(m*g), g, 2), g, 0)
id_s3 = sp.simplify(d2_s3 - (-sp.Rational(1, 3)*(m**2 - 1)*m))
id_s1_vs_4d = sp.simplify(d2_s1 - (-sp.Rational(1, 3)*(m**2 - 1)*m))
id_s1_2d = sp.simplify(d2_s1 - (-m**2*m))
check("je Mode gilt für S^3: d^2_gamma (sin m gamma/sin gamma)|_0 = -(1/3)(m^2 - 1) m, also die vorgeschlagene Identität "
      "d^2_gamma W = -(1/3)(d^2_tau - 1) W; für den nativen Nahtwinkel gilt stattdessen d^2_theta W = -d^2_tau W (2D-Welle), die "
      "4D-Identität verfehlt", id_s3 == 0 and id_s1_vs_4d != 0 and id_s1_2d == 0)
befund("Die native Querantwort der Quelle ist zweidimensional (Nahtwinkel). Die 4D-Fortsetzung braucht eine Winkelvariable, die "
       "das Original nicht hat; sie ist eine konsistente Fortsetzung, aber nicht erzwungen (wie die Analyse selbst einräumt). "
       "Erzwungen wäre sie, wenn T' die Querrichtungen als Lagen der Naht liefert", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Ursprung aus zwei binären Beobachtungen")


def bloch_z(psi):
    return psi[1]/psi[0] if abs(psi[0]) > 1e-14 else np.inf


def cross(z1, z2, z3, z4):
    return ((z1 - z3)*(z2 - z4))/((z1 - z4)*(z2 - z3))


rows = []
for p in (0.1, 0.3, 0.5, 0.7, 0.9):
    a, ap = np.array([1, 0j]), np.array([0j, 1])
    b = np.array([np.sqrt(p), np.sqrt(1 - p)*np.exp(0.4j)])
    bp = np.array([-np.sqrt(1 - p)*np.exp(-0.4j), np.sqrt(p)])
    zs = [0.0, None, bloch_z(b), bloch_z(bp)]
    lam = (zs[2] - zs[0])/(zs[3] - zs[0])                        # Kreuzverhältnis mit z2 = unendlich
    mix = (np.outer(a, a.conj()) + np.outer(ap, ap.conj()) + np.outer(b, b.conj()) + np.outer(bp, bp.conj()))/4
    rows.append((p, complex(lam), np.allclose(mix, np.eye(2)/2)))
ZAHLEN["kreuzverhaeltnis_p"] = [(p, round(l.real, 6), round(l.imag, 6)) for p, l, _ in rows]
harm = [p for p, l, _ in rows if abs(l + 1) < 1e-9]
check("zwei binäre Beobachtungen mit Überlapp p: die vier Ausgänge haben das Kreuzverhältnis -(1-p)/p (reell, ein Großkreis), "
      "harmonisch (tau = i) genau bei p = 1/2, also maximal inkompatibel; das gleichgewichtete Gemisch ist immer I/2",
      harm == [0.5] and all(abs(l + (1 - p)/p) < 1e-9 for p, l, _ in rows) and all(m_ for *_, m_ in rows))
befund("Die Paarung nach Beobachtung ({a, a-perp | b, b-perp}) ist die Antipodenpaarung, also das kanonische eta aus z70; die "
       "Clifford-Drehung a -> b -> a-perp -> b-perp ist die Vierteldrehung. Der Ursprungskandidat erklärt Marken, tau = i und eta; "
       "der Zustand bleibt offen: Das gleichgewichtete Gemisch ist der Spurzustand ohne Modulzeit, der Originalzustand braucht "
       "das Gewicht rho_m (z63)", True)

urteil = ("Die zweite externe Analyse hält in ihren prüfbaren Kernpunkten, mit einer Einschränkung. Die lokale Modulzeit desselben "
          "Originalzustands auf einem Nahtbogen hat ein dicht werdendes Spektrum ohne pi-Gitter, also löst sie den spektralen "
          "Widerspruch zur periodischen globalen Uhr. Die S^3-Fortsetzung W(tau, gamma) reproduziert die Stromantwort exakt auf der "
          "Linie gamma = 0, aber dieselbe Linienantwort hat die native 2D-Fortsetzung längs des Nahtwinkels, und für diesen gilt "
          "die 2D-Wellenidentität statt der vorgeschlagenen 4D-Identität: Die Dimension ist durch die Quelle nicht erzwungen. Der "
          "Ursprungskandidat 'zwei maximal inkompatible binäre Beobachtungen' gibt tau = i genau bei p = 1/2 und als Paarung das "
          "kanonische eta; der Zustand bleibt offen.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z71_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
