"""H6: Verklebung als verschränkter Zustand zweier Z4 Sektorsysteme, modulare Wirkung als Dynamik."""
import itertools, numpy as np
from fractions import Fraction as F
from common import Checks
c = Checks("h06_bell")
hD = [F(0), F(5, 8), F(1, 2), F(5, 8)]; hA = [F(a * (4 - a), 8) for a in range(4)]
c.ok("Sektorgewichte: D5 hat 0, 5/8, 1/2, 5/8 und A3 hat 0, 3/8, 1/2, 3/8", hA == [0, F(3, 8), F(1, 2), F(3, 8)])
c.ok("Kitaevs Stufen: Wirbelgewicht nu/16 gibt 5/8 bei nu = 10 und 3/8 bei nu = 6; 10 + 6 = 16", F(10, 16) == hD[1] and F(6, 16) == hA[1])
w = np.exp(2j * np.pi / 4)
SD = np.array([[w ** (a * b) for b in range(4)] for a in range(4)]) / 2          # Bicharakter 2 h(a+b) - ... = ab/4 (mod 1) für D5
SA = SD.conj()                                                                   # -ab/4 für A3
TD = np.diag([np.exp(2j * np.pi * float(h)) for h in hD]); TA = np.diag([np.exp(2j * np.pi * float(h)) for h in hA])
bich = lambda h, a, b: (h[(a + b) % 4] - h[a] - h[b]) % 1
c.ok("die S Matrizen sind die diskreten Fouriertransformationen zu den Bicharakteren +ab/4 (D5) und -ab/4 (A3): zueinander konjugiert",
     all(bich(hD, a, b) == F(a * b, 4) % 1 and bich(hA, a, b) == F(-a * b, 4) % 1 for a in range(4) for b in range(4)))
gm = lambda h: sum(np.exp(2j * np.pi * float(x)) for x in h) / 2
c.ok("Gauß Milgram: Summe der Phasen / 2 = exp(2 pi i c/8) mit c = 5 und c = 3; zusammen c = 8",
     abs(gm(hD) - np.exp(2j * np.pi * 5 / 8)) < 1e-12 and abs(gm(hA) - np.exp(2j * np.pi * 3 / 8)) < 1e-12)
c.ok("die Diagonale {(a,a)} ist isotrop: h_D(a) + h_A(a) ist ganzzahlig für alle a (Lagrange Untergruppe der Ordnung 4 = Wurzel aus 16)",
     all((hD[a] + hA[a]) % 1 == 0 for a in range(4)), [str(hD[a] + hA[a]) for a in range(4)])
v = np.zeros(16, complex)
for a in range(4): v[4 * a + a] = 0.5
S2 = np.kron(SD, SA); T2 = np.kron(TD, TA)
c.ok("der Klebevektor (1/2) Summe |a,a> ist fest unter S und T: das ist der E8 Abschluss als Zustand", np.allclose(S2 @ v, v) and np.allclose(T2 @ v, v))
rho = (v.reshape(4, 4) @ v.reshape(4, 4).conj().T)
c.ok("er ist maximal verschränkt: reduzierter Zustand auf jeder Seite ist 1/4", np.allclose(rho, np.eye(4) / 4))
inv = np.linalg.svd(np.vstack([S2 - np.eye(16), T2 - np.eye(16)]))[1]
v2 = np.zeros(16, complex)
for a in range(4): v2[4 * a + (-a) % 4] = 0.5
nfix = int(np.sum(inv < 1e-9))
c.ok("NICHT eindeutig: der Fixraum ist zweidimensional. Neben der Diagonale ist auch die Gegendiagonale (a, -a) fest (Ladungskonjugation auf einer Seite); beide geben E8",
     nfix == 2 and np.allclose(S2 @ v2, v2) and np.allclose(T2 @ v2, v2) and all((hD[a] + hA[(-a) % 4]) % 1 == 0 for a in range(4)), nfix)
# Dynamik? die von S, T erzeugte Gruppe ist endlich und wirkt als Produkt
def closure(gens, cap=5000):
    key = lambda M: tuple(np.round(M / (M.flat[np.flatnonzero(np.abs(M) > 1e-9)[0]] / abs(M.flat[np.flatnonzero(np.abs(M) > 1e-9)[0]])), 6).flatten())
    seen = {key(np.eye(len(gens[0]))): np.eye(len(gens[0]), dtype=complex)}; front = list(seen.values())
    while front and len(seen) < cap:
        nxt = []
        for M in front:
            for g in gens:
                N = g @ M; k = key(N)
                if k not in seen: seen[k] = N; nxt.append(N)
        front = nxt
    return len(seen)
gD = closure([SD, TD])
c.ok("die modulare Wirkung auf einer Seite erzeugt eine endliche Gruppe (projektiv)", gD < 5000, gD)
c.ok("DYNAMIK FÄLLT: S und T wirken als Produkt U_D (x) U_A. Die reduzierte Entwicklung jeder Seite ist unitär, alle Eigenwerte haben Betrag 1; ein Spektrum {1, (2/3)^6, (1/3)^6} ist unmöglich",
     np.allclose(np.abs(np.linalg.eigvals(SD)), 1) and np.allclose(np.abs(np.linalg.eigvals(TD)), 1))
c.ok("ABBRUCHKRITERIUM GREIFT: der reduzierte Zustand ist 1/4 auf vier Sektoren; die Referenz der vorhandenen Rückabbildung ist I/5 auf fünf Plätzen. Verschiedene Räume", rho.shape[0] == 4 != 5)
c.done({"gruppenordnung_projektiv": gD,
        "urteil": "Struktur exakt: der E8 Abschluss ist ein maximal verschränkter fester Zustand der sechzehn Sektoren (Kitaev 10 und 6); es gibt genau zwei, Diagonale und Gegendiagonale. Als Dynamik ausgeschlossen: endliche Gruppe, Produktwirkung, keine Relaxation. Das eigene Abbruchkriterium greift."})
