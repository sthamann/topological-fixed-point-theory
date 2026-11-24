"""Z76: Prüfung der vierten externen Analyse 'Weiterentwicklung von H6h34' (24.11.): Ladungsuhr, Zweipfadkanal, Winkeloperator.

Teil 1  K/pi = Q_a - Q_b (mod 4): Die modulare Vierteldrehung ist die Differenz der Quarterladungen (vgl. z74: R ist
        Ladungsdrehung mal Parität). Damit sind rho_1 (geladene Feldalgebra) und rho_4 (neutrale Algebra C[z^4]) zwei
        Auflösungen desselben Prozesses.
Teil 2  Zweipfadkanal zwischen den Paarschalen H_m -> H_(m+1) mit den beiden nativen Elementarschritten (Teilchen hoch, Loch
        runter): Die Forderung, den Originalzustand (auf H_m der Spurzustand, z73) zu erhalten, legt die Gewichte eindeutig fest:
        p_k = (k+1)/(m+1), q_k = (m-k)/(m+1). Ihre Wurzeln sind die Schwinger-Oszillatorgewichte (CCR), also genau die
        Normierung, die z73 vermisst hat.
Teil 3  Petz-Rückkanal R: m(m+1)(I - R Phi) = Summe [J_a,[J_a, .]] (Spin (m-1)/2), Eigenwerte l(l+1), l = 0..m-1, mit
        Vielfachheit 2l+1: der Laplace der quantisierten Sphäre.
Teil 4  Zusammengesetzte Kanäle H_m -> H_(m+n) haben Krausrang n+1 (statt 2^n Pfadfolgen).
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
say("Teil 1: Ladungsuhr")
Qr = range(-12, 13)
ok1 = all((sp.Rational(4)*((sp.Rational(1, 2)*qa**2 - sp.Rational(1, 4)*qa) + (sp.Rational(1, 2)*qb**2 + sp.Rational(1, 4)*qb))
           - (qa - qb)) % 4 == 0 for qa in Qr for qb in Qr)
check("K/pi = 4(E_a(Q_a) + E_b(Q_b)) + 4 N_osc mit E_alpha(Q) = Q^2/2 + (alpha - 1/2)Q ist kongruent Q_a - Q_b mod 4 "
      "(alle |Q| <= 12): die modulare Vierteldrehung misst die Ladungsdifferenz der beiden Quarter", ok1)
befund("rho_1 (geladene Feldalgebra, Lücke pi) und rho_4 (neutrale Algebra C[z^4] = C[z]^Z4, Lücke 4 pi) sind zwei "
       "Auflösungen desselben Prozesses, keine konkurrierenden Lösungen (deckt sich mit z74); die Zuordnung der vier Marken zur "
       "Quotientengeometrie bleibt offen", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Zustandserhaltender Zweipfadkanal H_m -> H_(m+1)")


def weights_from_state(m):
    """Unbekannte p_k (k = 0..m-1), q_k = 1 - p_k; Forderung Phi(I/m) = I/(m+1) (Diagonale j = 0..m)."""
    p = sp.symbols(f"p0:{m}")
    eqs = []
    for j in range(m + 1):
        lhs = (p[j - 1] if j >= 1 else 0) + ((1 - p[j]) if j <= m - 1 else 0)
        eqs.append(sp.Eq(lhs/m, sp.Rational(1, m + 1)))
    sol = sp.solve(eqs, p, dict=True)
    return sol


ok2 = True
for m in range(2, 8):
    sol = weights_from_state(m)
    ok2 &= len(sol) == 1 and all(sol[0][sp.Symbol(f"p{k}")] == sp.Rational(k + 1, m + 1) for k in range(m))
check("für m = 2..7 legt die Forderung 'Phi bildet den Originalzustand auf H_m (Spurzustand) auf den auf H_(m+1) ab' die "
      "Gewichte eindeutig fest: p_k = (k+1)/(m+1), q_k = (m-k)/(m+1)", ok2)


def kraus(m):
    Ap = np.zeros((m + 1, m))
    Ah = np.zeros((m + 1, m))
    for k in range(m):
        Ap[k + 1, k] = np.sqrt((k + 1)/(m + 1))          # Teilchen eine Stufe hoch
        Ah[k, k] = np.sqrt((m - k)/(m + 1))              # Loch eine Stufe runter
    return [Ap, Ah]


ok_ccr = True
for m in range(2, 8):
    Ap, Ah = kraus(m)
    # Schwinger: Spin j = (m-1)/2 mit n1 = k, n2 = m-1-k; a1^+ -> sqrt(n1+1), a2^+ -> sqrt(n2+1) = sqrt(m-k)
    ok_ccr &= np.allclose(np.diag(Ap, -1)*np.sqrt(m + 1), np.sqrt(np.arange(1, m + 1)))
    ok_ccr &= np.allclose(np.diag(Ah[:m, :m])*np.sqrt(m + 1), np.sqrt(np.arange(m, 0, -1)))
    ok_ccr &= np.allclose(sum(A.T @ A for A in (Ap, Ah)), np.eye(m))
check("die Krausamplituden sind sqrt(k+1)/sqrt(m+1) und sqrt(m-k)/sqrt(m+1), genau die Schwinger-Oszillatorgewichte "
      "(CCR) der Teilchen- und Lochleiter, und der Kanal ist spurerhaltend", ok_ccr)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Petz-Rückkanal und Winkeloperator")


def superop(fn, n):
    S = np.zeros((n*n, n*n))
    for a in range(n*n):
        E = np.zeros(n*n)
        E[a] = 1
        S[:, a] = fn(E.reshape(n, n)).reshape(-1)
    return S


def spin_ops(j):
    ms = np.arange(j, -j - 1, -1)
    n = len(ms)
    Jz = np.diag(ms)
    Jp = np.zeros((n, n))
    for i in range(1, n):
        Jp[i - 1, i] = np.sqrt(j*(j + 1) - ms[i]*(ms[i] + 1))
    Jx, Jy = (Jp + Jp.T)/2, (Jp - Jp.T)/(2j)
    return Jx, Jy, Jz


ok3 = True
spectra = {}
for m in range(2, 8):
    Ks = kraus(m)
    Phi = lambda X, Ks=Ks: sum(A @ X @ A.T for A in Ks)
    Rec = lambda Y, Ks=Ks, m=m: (m + 1)/m*sum(A.T @ Y @ A for A in Ks)      # Petz mit Referenz I/m, Phi(I/m) = I/(m+1)
    RPhi = superop(lambda X: Rec(Phi(X)), m)
    L = m*(m + 1)*(np.eye(m*m) - RPhi)
    Js = spin_ops((m - 1)/2)
    Cas = superop(lambda X, Js=Js: sum(J @ (J @ X - X @ J) - (J @ X - X @ J) @ J for J in Js), m).real
    ev = np.sort(np.linalg.eigvals(L).real)
    expect = np.sort(np.concatenate([[l*(l + 1)]*(2*l + 1) for l in range(m)]))
    ok3 &= np.allclose(L, Cas, atol=1e-10) and np.allclose(ev, expect, atol=1e-9)
    spectra[m] = sorted(set(np.round(ev, 6).tolist()))
ZAHLEN["winkelspektren"] = spectra
check("m(m+1)(I - R Phi) ist exakt der Casimir-Superoperator Summe [J_a,[J_a, .]] zu Spin (m-1)/2 (m = 2..7), Eigenwerte "
      "l(l+1) mit Vielfachheit 2l+1, l = 0..m-1 (m = 5: 0, 2, 6, 12, 20 auf 1+3+5+7+9 = 25 Richtungen)", ok3)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Zusammengesetzte Kanäle")


def compose_rank(m, n):
    Ks = [np.eye(m)]
    dim = m
    for s in range(n):
        Kn = kraus(dim)
        Ks = [B @ A for A in Ks for B in Kn]
        dim += 1
    choi = sum(np.outer(A.reshape(-1), A.reshape(-1)) for A in Ks)
    return np.linalg.matrix_rank(choi, tol=1e-10), len(Ks)


ranks = {n: compose_rank(2, n) for n in (1, 2, 3, 4, 6)}
ZAHLEN["krausrang"] = ranks
check("n Schritte H_2 -> H_(2+n): 2^n Pfadfolgen, Krausrang n+1 (n = 1, 2, 3, 4, 6): die Pfadinformation verdichtet sich auf die "
      "Zahl der Teilchenschritte", all(r == n + 1 and k == 2**n for n, (r, k) in ranks.items()))
befund("Die Oszillatornormierung (offener Punkt (b) aus H6h35) folgt bedingt: Mit den beiden nativen Elementarschritten als "
       "Krausstruktur erzwingt die Erhaltung des Originalzustands die Schwinger-CCR-Gewichte, und Kanal mit Petz-Rückkanal gibt "
       "die su(2)-Casimir-Struktur der Sphäre. Einschränkung (wie die Analyse selbst): Der native Stromoperator koppelt auch an "
       "weitere Paare; der Kanal ist eine zustandsbalancierte Auslesung, nicht die ganze Dynamik", True)
befund("Die vierte Analyse nennt außerdem die Klich-Determinante für Mehrpunktantworten (Standard) und den Einstein-Dreipunkt"
       "vertex mit Konstruierbarkeit der Baumgravitation (Standard unter lokaler unitärer 4D-Streuung); beides ist bedingt, G "
       "folgt nicht", True)

urteil = ("Die vierte externe Analyse hält in den geprüften Punkten. K/pi ist kongruent Q_a - Q_b mod 4, also sind rho_1 und rho_4 "
          "zwei Auflösungen desselben Prozesses. Der Zweipfadkanal zwischen den Paarschalen mit Teilchen- und Lochschritt ist durch "
          "die Erhaltung des Originalzustands eindeutig: p_k = (k+1)/(m+1), q_k = (m-k)/(m+1). Seine Krausamplituden sind genau die "
          "Schwinger-CCR-Gewichte, und Kanal mit Petz-Rückkanal gibt exakt den Casimir m(m+1)(I - R Phi) = Summe [J_a,[J_a, .]] mit "
          "l(l+1). Damit ist die in H6h35 fehlende Oszillatornormierung bedingt durch den Zustand erzwungen; offen bleibt, dass "
          "der Kanal nur eine zustandsbalancierte Auslesung der vollen Dynamik ist.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z76_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
