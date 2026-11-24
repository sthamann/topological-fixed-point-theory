"""Z77: Prüfung der fünften externen Analyse 'Der vollständige relative Antwortkanal' (24.11.).

Teil 1  [J_-1, B_(r,m)] = B_(r,m+1) - B_(r-1,m+1) für B_(r,m) = psi*_(r+m) psi_r gilt als Operatoridentität auf dem vollen
        endlichen Fockraum (die Mehrpaaranteile kürzen sich im Kommutator).
Teil 2  Thermischer Zweipfadkanal: Die Formel der Analyse p_m(r) = Summe_j (j+1) w(r+j)/((m+1) Summe_j w(r+j)) ist der
        thermische Mittelwert der Vakuumgewichte über die Kantenlage und geht im Vakuum in z76 über; die einfache normierte
        Zustandserhaltung W_(m+1)/(m+1) = [p W_m(r) + (1 - p(r+1)) W_m(r+1)]/m erfüllt sie für m >= 2 nicht.
Teil 3  Vakuum gegen thermisch: Im Vakuum gibt der Petz-Rückkanal exakt l(l+1); die thermischen Werte der Analyse sind
        ohne ihre Konditionierung nicht reproduziert.
"""
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sps

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
say("Teil 1: Kommutatoridentität auf dem vollen Fockraum")


def fock(n):
    I2 = sps.identity(2, format="csr")
    Zp = sps.diags([1.0, -1.0], format="csr")
    sm = sps.csr_matrix(np.array([[0.0, 1.0], [0.0, 0.0]]))
    ops = []
    for j in range(n):
        mats = [Zp]*j + [sm] + [I2]*(n - j - 1)
        op = mats[0]
        for M in mats[1:]:
            op = sps.kron(op, M, format="csr")
        ops.append(op)
    return ops


n = 10
cs = fock(n)
Jm1 = sum(cs[s + 1].T @ cs[s] for s in range(n - 1))


def B(r, m):
    return cs[r + m].T @ cs[r]


err = 0.0
for m in range(1, 5):
    for r in range(1, n - m - 1):
        lhs = Jm1 @ B(r, m) - B(r, m) @ Jm1
        rhs = B(r, m + 1) - B(r - 1, m + 1)
        err = max(err, abs(lhs - rhs).max())
check(f"[J_-1, B_(r,m)] = B_(r,m+1) - B_(r-1,m+1) auf dem vollen Fockraum von {n} Moden (alle inneren r, m = 1..4), "
      f"maximaler Fehler {err:.1e}", err < 1e-12)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Thermischer Zweipfadkanal")


def thermal_ok(beta, alpha, m):
    """Relatives Residuum der einfachen normierten Zustandserhaltung mit der Formel der Analyse (nur Paare mit W > 1e-12)."""
    c = lambda x: 1/(1 + np.exp(np.clip(beta*(x + alpha), -700, 700)))
    w = lambda r: c(r) - c(r + 1)
    Wm = lambda mm, r: c(r) - c(r + mm)
    p = lambda r: sum((j + 1)*w(r + j) for j in range(m))/((m + 1)*sum(w(r + j) for j in range(m)))
    res = []
    for r in range(-m - 3, 4):
        if min(Wm(m, r), Wm(m, r + 1), Wm(m + 1, r)) < 1e-12:
            continue
        lhs = Wm(m + 1, r)/(m + 1)
        rhs = (p(r)*Wm(m, r) + (1 - p(r + 1))*Wm(m, r + 1))/m
        res.append(abs(lhs - rhs)/lhs)
    return max(res) if res else 0.0, 0, 1


res_max = 0.0
for beta in (2.0, 4*np.pi):
    for alpha in (0.25, 0.75):
        for m in range(2, 6):
            e, pmin, pmax = thermal_ok(beta, alpha, m)
            res_max = max(res_max, e)
ZAHLEN["residuum_einfache_erhaltung"] = res_max
check(f"die Formel p_m(r) = Summe (j+1) w(r+j)/((m+1) Summe w(r+j)) erfüllt die einfache normierte Zustandserhaltung "
      f"W_(m+1)/(m+1) = [p W_m(r) + (1 - p(r+1)) W_m(r+1)]/m NICHT (m = 2..5, relatives Residuum bis {res_max:.2f}); sie ist der thermische "
      f"Mittelwert der Vakuumgewichte (j+1)/(m+1) über die Lage der Fermikante in der Paarschale", res_max > 1e-4)
beta, m = 400.0, 3
cz = lambda x: 1/(1 + np.exp(np.clip(beta*(x + 0.25), -700, 700)))
vac_ok = True
for r in range(-m, 0):
    j_edge = -r - 1                                        # Kante zwischen r+j = -1 und 0 (x = -3/4 | 1/4)
    pr = sum((j + 1)*(cz(r + j) - cz(r + j + 1)) for j in range(m))/((m + 1)*(cz(r) - cz(r + m)))
    vac_ok &= abs(pr - (j_edge + 1)/(m + 1)) < 1e-9
check("im Vakuumgrenzfall (beta = 400) wird die Formel genau zu den Gewichten (k+1)/(m+1) aus z76, mit k = Kantenlage", vac_ok)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Vakuum gegen thermisch (zentrale Paare, Petz-Rückkanal)")


def petz_spectrum(weights_m, weights_m1, p):
    """Kanal Phi: diag-Zustand auf m Paaren -> m+1 Paare mit Pfadamplituden sqrt(p), sqrt(1-p); Petz bzgl. weights_m."""
    m = len(weights_m)
    Ap = np.zeros((m + 1, m))
    Ah = np.zeros((m + 1, m))
    for k in range(m):
        Ap[k + 1, k] = np.sqrt(p[k])
        Ah[k, k] = np.sqrt(1 - p[k])
    sig = np.diag(weights_m)
    Phi = lambda X: Ap @ X @ Ap.T + Ah @ X @ Ah.T
    img = Phi(sig)
    isq = np.diag(1/np.sqrt(np.diag(img)))
    sq = np.diag(np.sqrt(weights_m))
    Rec = lambda Y: sq @ (Ap.T @ (isq @ Y @ isq) @ Ap + Ah.T @ (isq @ Y @ isq) @ Ah) @ sq
    S = np.zeros((m*m, m*m))
    for a in range(m*m):
        E = np.zeros(m*m)
        E[a] = 1
        S[:, a] = Rec(Phi(E.reshape(m, m))).reshape(-1)
    return np.sort(np.linalg.eigvals(m*(m + 1)*(np.eye(m*m) - S)).real), img


m = 2
vac_w = np.ones(m)/m
p_vac = [(k + 1)/(m + 1) for k in range(m)]
ev_vac, _ = petz_spectrum(vac_w, None, p_vac)
c = lambda x: 1/(1 + np.exp(np.clip(4*np.pi*(x + 0.25), -700, 700)))
rows = [r for r in range(-m, 0)]                          # zentrale Paare (r, r+m) über der Fermikante
w_th = np.array([c(r) - c(r + m) for r in rows])
w_th = w_th/w_th.sum()
p_th = [sum((j + 1)*(c(r + j) - c(r + j + 1)) for j in range(m))/((m + 1)*(c(r) - c(r + m))) for r in rows]
ev_th, _ = petz_spectrum(w_th, None, p_th)
ZAHLEN["petz_m2"] = {"vakuum": np.round(ev_vac, 6).tolist(), "thermisch_beta4pi": np.round(ev_th, 6).tolist()}
say(f"   Vakuum {np.round(ev_vac, 4).tolist()}, thermisch (beta = 4 pi, zentrale Paare) {np.round(ev_th, 4).tolist()}")
check("im Vakuum ist das skalierte Petz-Spektrum exakt (0, 2, 2, 2) = l(l+1)", np.allclose(ev_vac, [0, 2, 2, 2], atol=1e-9))
befund("Thermisch (beta = 4 pi, zentrale Paare, Gewichte nach der Formel der Analyse) ergibt diese Auslesung " +
       str(np.round(ev_th, 4).tolist()) + "; die Werte der Analyse (0; 1,9864; 1,9864; 2,0281) sind ohne ihre Konditionierung nicht "
       "reproduziert, die Aufspaltung im thermischen Fall ist also nur qualitativ belegt", True)
befund("Damit ist die offene Frage 'Vakuum gegen thermisch' (H6h35 (c)) für die Winkelstruktur präzisiert: Die Sphäre ist "
       "exakt im Vakuumsektor und im thermischen Original leicht aufgespalten; der Unterschied ist das thermische Gedächtnis "
       "der Auslesung. Die Young-Gitter-Lesart der höheren Paare (Plancherel) und die Wahl des Instruments (kohärent gegen "
       "gemessen) sind nicht nachgerechnet", True)

urteil = ("Die fünfte externe Analyse hält nur teilweise. [J_-1, B_(r,m)] = B_(r,m+1) - B_(r-1,m+1) gilt als Operatoridentität "
          "auf dem vollen Fockraum. Ihre thermischen Übergangswahrscheinlichkeiten p_m(r) = Summe (j+1) w(r+j)/((m+1) Summe w(r+j)) "
          "sind der thermische Mittelwert der Vakuumgewichte über die Kantenlage und gehen im Vakuum in z76 über; die einfache "
          "normierte Zustandserhaltung zwischen den Schalen erfüllen sie nicht, ihre Eindeutigkeit hängt an einer Konditionierung, "
          "die ohne das Paket nicht nachvollziehbar ist. Im Vakuum ist das Petz-Spektrum exakt l(l+1); die thermische Aufspaltung "
          "(0; 1,9864; 1,9864; 2,0281) ist nicht reproduziert. Die Sphäre bleibt eine Vakuumgröße.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z77_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
