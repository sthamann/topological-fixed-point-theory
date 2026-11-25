"""Z80: Prüfung der siebten und achten externen Analyse (25.11.) gegen den Bestand z65, z71, z74-z78.

Original wie z76/z77: c(x) = 1/(1 + e^(beta x)) auf x = r + alpha, beta = 4 pi; w(r) = c(r) - c(r+1),
W_m(r) = c(r) - c(r+m) = Summe_(j<m) w(r+j). Paar (r, r+m) = psi*_(r+m) psi_r: Loch bei r, Teilchen bei r+m.

Teil 1  (Analyse 7) Die Zuordnung in z77 war vertauscht. Mit u_m(r) = Summe_j (m-j) w(r+j)/((m+1) W_m(r)) für 'Teilchen hoch'
        und 1 - u_m(r) für 'Loch runter' gilt u(r) W_m(r) + (1 - u(r+1)) W_m(r+1) = m/(m+1) W_(m+1)(r) exakt (symbolisch für
        allgemeines w, numerisch für jede Temperatur). Die Formel der fünften Analyse ist genau das Lochgewicht 1 - u, der
        Vakuumgrenzfall ist z76 mit richtigen Etiketten, die in [0, 1] beschränkte Lösung ist eindeutig (homogen: C/W_m(r)).
Teil 2  z77 Teil 3 neu gerechnet (Petz-Rückkanal auf den zentralen Paaren, m = 2, beta = 4 pi, alpha = 1/4, gleiche
        Abschneidung und Normierung) mit der richtigen Zuordnung; Vergleich mit (0; 1,9864; 1,9864; 2,0281) und Varianten.
Teil 3  (Analyse 8 §1) Mischung verschobener Vakua: c(x) = Summe_e w(e) [x <= e]; W_m und u_m sind Mischungsmittel der
        z76-Größen, Petz je Kantenkomponente gibt l(l+1). Gegenbeispiel: der quasifreie Zustand ist nicht die Mischung. Der
        behauptete Anteil 1/(m+1) 'verfehlter' Übergangsamplituden von J_-1 auf dem Jordan-Wigner-Fockraum.
Teil 4  (Analysen 7 §3, 8 §2) Relative Zustände aus z78: U(1)-invariante reelle Quadriken in c sind genau c^+ X c, X in Herm(2);
        <X|Y> = 1/det(I - (i/2)(Y - X)) auf dem Fockraum; s^2 = 4(1 - Re 1/K) = det(Y - X), K = 4/(r^2 - (t + 2i)^2); das
        ist der Kern nu = 1 aus z74 an den Röhrenpunkten X + iI.
Teil 5  (Analyse 8 §4) Die Wurzel des Kerns ist nicht positiv: E[(det Q)^2] = -3/16 bei nu = 1/2, Gramtest; Wallach-Menge.
Teil 6  (Analyse 7 §5) Zusammengesetzter Spin 2 eines freien masselosen 4D-Skalars: rho_2(s) ~ s^2, stetig, kein Pol.
Teil 7  (Analyse 8 §3) phi_ABCD = c_A c_B c_C c_D e^(i p.x), p = c c^+, löst die masselose Spin-2-Gleichung (Helizität 2).
Teil 8  (Analyse 7 §6) A^2 = B^2 = I, irreduzibel, Vorzeichensymmetrie A -> -A erzwingt AB + BA = 0: Cl_2(C) = M_2(C),
        Kreuzverhältnis -1 (tau = i), Überlapp p = 1/2; Vergleich mit z71 Teil 3.
Teil 9  (Analyse 8 §6) Flavourauslesung |U_ai|^2 = Tr(P_a Q_i), J = Im Tr(P_a Q_i P_b Q_j): Standard (Jarlskog); TFPT-STAND §5.
Teil 10 (Analyse 8 §5) Jacobson-Kette: nicht gerechnet, nur Befund (Firewall GRAV.NONCIRCULAR.01).
"""
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

HERE = Path(__file__).parent
ROOT = HERE.parents[2]
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
mp.mp.dps = 50
rng = np.random.default_rng(80)

I2 = np.eye(2)
SIG = [np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]),
       np.diag([1.0 + 0j, -1.0])]
ETA = np.diag([1.0, -1.0, -1.0, -1.0])


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def herm(v):
    return sum(v[k]*SIG[k] for k in range(4))


def vec(H):
    return np.array([np.trace(H @ s).real/2 for s in SIG])


def w_mp(beta, alpha, r):
    """w(r) = c(r) - c(r+1) ohne Auslöschung: e^(b x)(e^b - 1)/((1 + e^(b x))(1 + e^(b (x+1)))), x = r + alpha."""
    x = r + alpha
    return mp.exp(beta*x)*mp.expm1(beta)/((1 + mp.exp(beta*x))*(1 + mp.exp(beta*(x + 1))))


def c_mp(beta, alpha, r):
    return 1/(1 + mp.exp(beta*(r + alpha)))


def fock(n):
    """Jordan-Wigner wie z77: Index 0 je Faktor = leer, Index 1 = besetzt."""
    I2s = sps.identity(2, format="csr")
    Zp = sps.diags([1.0, -1.0], format="csr")
    sm = sps.csr_matrix(np.array([[0.0, 1.0], [0.0, 0.0]]))
    ops = []
    for j in range(n):
        mats = [Zp]*j + [sm] + [I2s]*(n - j - 1)
        op = mats[0]
        for M in mats[1:]:
            op = sps.kron(op, M, format="csr")
        ops.append(op)
    return ops


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Thermischer Zweipfadkanal mit richtiger Zuordnung (Analyse 7)")
ok_sym, ok_alt, ok_komp = True, True, True
for m in range(2, 11):
    w = sp.symbols(f"w0:{m + 1}", positive=True)                  # w_j = w(r+j), j = 0..m
    Wr, Wr1, Wn = sum(w[:m]), sum(w[1:]), sum(w)                   # W_m(r), W_m(r+1), W_(m+1)(r)
    u_r = sum((m - j)*w[j] for j in range(m))/((m + 1)*Wr)
    u_r1 = sum((m - j)*w[j + 1] for j in range(m))/((m + 1)*Wr1)
    ok_sym &= sp.cancel(u_r*Wr + (1 - u_r1)*Wr1 - sp.Rational(m, m + 1)*Wn) == 0
    p_r = sum((j + 1)*w[j] for j in range(m))/((m + 1)*Wr)        # Formel der fünften Analyse
    p_r1 = sum((j + 1)*w[j + 1] for j in range(m))/((m + 1)*Wr1)
    rest = sp.expand(sp.cancel(p_r*Wr + (1 - p_r1)*Wr1 - sp.Rational(m, m + 1)*Wn))
    ok_alt &= rest != 0 and sp.expand(rest - ((1 - m)*(w[0] + w[m]) + 2*sum(w[1:m]))/(m + 1)) == 0
    ok_komp &= sp.cancel(p_r + u_r - 1) == 0
check("symbolisch für allgemeines w (Symbole w_0..w_m, m = 2..10): mit u_m(r) = Summe_j (m-j) w(r+j)/((m+1) W_m(r)) für "
      "'Teilchen hoch' und 1 - u für 'Loch runter' gilt u(r) W_m(r) + (1 - u(r+1)) W_m(r+1) = m/(m+1) W_(m+1)(r) exakt", ok_sym)
check("die Formel der fünften Analyse p = Summe (j+1) w(r+j)/((m+1) W_m(r)) ist identisch 1 - u, also genau das Gewicht "
      "'Loch runter'; als 'Teilchen hoch' eingesetzt (z77) bleibt das Residuum [(1-m)(w_0 + w_m) + 2 Summe w_1..w_(m-1)]/(m+1) "
      "!= 0", ok_komp and ok_alt)

BETAS = {"2": mp.mpf(2), "4 pi": 4*mp.pi, "60": mp.mpf(60)}
ALPHAS = {"1/4": mp.mpf(1)/4, "1/2": mp.mpf(1)/2, "3/4": mp.mpf(3)/4}
res_max, u_band = mp.mpf(0), True
for B in BETAS.values():
    for A in ALPHAS.values():
        wv = {r: w_mp(B, A, r) for r in range(-30, 31)}
        for m in range(2, 13):
            Wf = lambda mm, r, wv=wv: mp.fsum(wv[r + j] for j in range(mm))
            uf = lambda r, wv=wv, m=m, Wf=Wf: mp.fsum((m - j)*wv[r + j] for j in range(m))/((m + 1)*Wf(m, r))
            for r in range(-m - 8, 9):
                lhs = uf(r)*Wf(m, r) + (1 - uf(r + 1))*Wf(m, r + 1)
                rhs = mp.mpf(m)/(m + 1)*Wf(m + 1, r)
                res_max = max(res_max, abs(lhs - rhs)/rhs)
                u_band &= mp.mpf(1)/(m + 1) - mp.mpf(10)**-45 <= uf(r) <= mp.mpf(m)/(m + 1) + mp.mpf(10)**-45
ZAHLEN["residuum_neu_max"] = float(res_max)


def z77_residuum(beta, alpha, m, zuordnung):
    """Wörtlich die Schleife aus z77 Teil 2 (Gleitkomma); zuordnung 'alt' = Formel als Teilchen hoch, 'neu' = 1 - Formel."""
    c = lambda x: 1/(1 + np.exp(np.clip(beta*(x + alpha), -700, 700)))
    w_ = lambda r: c(r) - c(r + 1)
    Wm = lambda mm, r: c(r) - c(r + mm)
    p = lambda r: sum((j + 1)*w_(r + j) for j in range(m))/((m + 1)*sum(w_(r + j) for j in range(m)))
    q = p if zuordnung == "alt" else (lambda r: 1 - p(r))
    res = []
    for r in range(-m - 3, 4):
        if min(Wm(m, r), Wm(m, r + 1), Wm(m + 1, r)) < 1e-12:
            continue
        lhs = Wm(m + 1, r)/(m + 1)
        rhs = (q(r)*Wm(m, r) + (1 - q(r + 1))*Wm(m, r + 1))/m
        res.append(abs(lhs - rhs)/lhs)
    return max(res) if res else 0.0


res_alt = max(z77_residuum(b, a, m, "alt") for b in (2.0, 4*np.pi) for a in (0.25, 0.75) for m in range(2, 6))
res_neu_f = max(z77_residuum(b, a, m, "neu") for b in (2.0, 4*np.pi) for a in (0.25, 0.75) for m in range(2, 6))
ZAHLEN["residuum_z77_alt"] = res_alt
ZAHLEN["residuum_z77_schleife_neu"] = res_neu_f
check(f"numerisch (mpmath, 50 Stellen; beta = 2, 4 pi, 60; alpha = 1/4, 1/2, 3/4; m = 2..12; alle r mit |r| <= m + 8): "
      f"maximales relatives Residuum {float(res_max):.1e}, und u liegt stets in [1/(m+1), m/(m+1)] c [0, 1]. Zum Vergleich in "
      f"der wörtlichen z77-Schleife: alte Zuordnung Residuum {res_alt:.4f} (z77: 0,94), richtige Zuordnung {res_neu_f:.1e}",
      res_max < mp.mpf(10)**-40 and u_band and abs(res_alt - 0.93774) < 1e-4 and res_neu_f < 1e-12)

ok_vak = True
for m in range(2, 9):
    for r in range(-m, 0):
        wv = {s: sp.Integer(1 if s == -1 else 0) for s in range(r, r + m + 1)}   # Vakuum: Kante zwischen -1 und 0
        u = sum((m - j)*wv[r + j] for j in range(m))/((m + 1)*sum(wv[r + j] for j in range(m)))
        k = r + m                                                               # Teilchenlage (Etikett von z76)
        ok_vak &= u == sp.Rational(k + 1, m + 1) and 1 - u == sp.Rational(m - k, m + 1)
ok_etikett, verschieden = True, 0
for m in range(2, 7):
    cz = lambda x: 1/(1 + np.exp(np.clip(400*(x + 0.25), -700, 700)))
    for r in range(-m, 0):
        pr = sum((j + 1)*(cz(r + j) - cz(r + j + 1)) for j in range(m))/((m + 1)*(cz(r) - cz(r + m)))
        k = r + m
        ok_etikett &= abs(pr - (m - k)/(m + 1)) < 1e-9
        if 2*k != m - 1:
            ok_etikett &= abs(pr - (k + 1)/(m + 1)) > 1e-3
            verschieden += 1
check("Vakuumgrenzfall (exakt, m = 2..8): u -> (k+1)/(m+1) und 1 - u -> (m-k)/(m+1) mit k = r + m (Teilchenlage, Etikett von "
      "z76). Die z77-Formel geht bei beta = 400 in (m-k)/(m+1) über, das Lochgewicht von z76; z77 hat mit k = Kantenlage "
      f"-r-1 = m-1-k verglichen, also gespiegelt ({verschieden} Paare mit (k+1)/(m+1) != (m-k)/(m+1), m = 2..6)",
      ok_vak and ok_etikett)

with mp.workdps(260):
    B_, A_, m = 4*mp.pi, mp.mpf(1)/4, 3
    wv = {r: w_mp(B_, A_, r) for r in range(-36, 37)}
    Wf = lambda mm, r: mp.fsum(wv[r + j] for j in range(mm))
    uf = lambda r: mp.fsum((m - j)*wv[r + j] for j in range(m))/((m + 1)*Wf(m, r))
    r0 = -2

    def lauf(delta):
        """Die Rekursion vorwärts und rückwärts von u(r0) + delta aus integrieren (sie ist erster Ordnung)."""
        v = {r0: uf(r0) + delta}
        for r in range(r0, 30):
            v[r + 1] = 1 - (mp.mpf(m)/(m + 1)*Wf(m + 1, r) - v[r]*Wf(m, r))/Wf(m, r + 1)
        for r in range(r0 - 1, -31, -1):
            v[r] = (mp.mpf(m)/(m + 1)*Wf(m + 1, r) - (1 - v[r + 1])*Wf(m, r + 1))/Wf(m, r)
        return v

    v0 = lauf(0)
    dev0 = max(abs(v0[r] - uf(r)) for r in range(-30, 31))          # Rundung 1e-260, verstärkt um W(r0)/W(r) <= 1e164
    ok_hom, raus = dev0 < mp.mpf(10)**-80, {}
    for delta in (mp.mpf(10)**-12, -mp.mpf(10)**-12):
        vd = lauf(delta)
        form = max(abs((vd[r] - uf(r))*Wf(m, r)/Wf(m, r0) - delta)/abs(delta) for r in range(-30, 31))
        out = [r for r in range(-30, 31) if not (0 <= vd[r] <= 1)]
        links, rechts = max([r for r in out if r < r0], default=None), min([r for r in out if r > r0], default=None)
        ok_hom &= form < mp.mpf(10)**-100 and links is not None and rechts is not None
        raus[mp.nstr(delta, 2)] = [links, rechts]
    inv_w = float(mp.log10(1/Wf(m, 30)))
ZAHLEN["eindeutigkeit"] = {"abweichung_delta_0": float(dev0), "verlassen_von_[0,1]_bei_r": raus, "log10_1/W_3(30)": inv_w}
check(f"Eindeutigkeit (beta = 4 pi, alpha = 1/4, m = 3, 260 Stellen): die Rekursion ist erster Ordnung, jede Lösung ist "
      f"u + C/W_m(r) (C = delta W_m(r0), r0 = -2, exakt bestätigt); 1/W_3 wächst bis 10^{inv_w:.0f} bei r = 30, eine Störung "
      f"delta = +-1e-12 verlässt [0, 1] auf beiden Seiten (innerste Austrittsstellen r = {raus}): die beschränkte Lösung "
      "in [0, 1] ist eindeutig u", ok_hom)
befund("Das Negativ von z77/H6h37 ('hält nicht wie behauptet', Residuum 0,94) war unsere Fehlzuordnung: Die Formel der fünften "
       "Analyse ist das Gewicht 'Loch runter', z77 hat sie als 'Teilchen hoch' eingesetzt und den Vakuumvergleich mit dem "
       "gespiegelten Etikett k = Kantenlage geführt. Korrigiert ist der thermische Zweipfadkanal bei jeder Temperatur exakt "
       "zustandserhaltend, eindeutig und geht im Vakuum in z76 über", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Petz-Rückkanal auf den zentralen Paaren, neu mit richtiger Zuordnung (z77 Teil 3)")


def petz_spektrum(gew, p, m, ref=None):
    """Wie z77: Kanal auf len(gew) Paaren -> len(gew)+1 Paare, Teilchen hoch sqrt(p), Loch runter sqrt(1-p); Petz bzgl. gew.
    ref = None: Petz mit dem Bild Phi(sigma) (z77); sonst mit dem angegebenen Zielzustand. Skalierung m(m+1)."""
    n = len(gew)
    Ap, Ah = np.zeros((n + 1, n)), np.zeros((n + 1, n))
    for k in range(n):
        Ap[k + 1, k] = np.sqrt(p[k])
        Ah[k, k] = np.sqrt(1 - p[k])
    Phi = lambda X: Ap @ X @ Ap.T + Ah @ X @ Ah.T
    img = Phi(np.diag(gew)) if ref is None else np.diag(ref)
    isq = np.diag(1/np.sqrt(np.diag(img)))
    sq = np.diag(np.sqrt(gew))
    Rec = lambda Y: sq @ (Ap.T @ (isq @ Y @ isq) @ Ap + Ah.T @ (isq @ Y @ isq) @ Ah) @ sq
    S = np.zeros((n*n, n*n))
    for a in range(n*n):
        E = np.zeros(n*n)
        E[a] = 1
        S[:, a] = Rec(Phi(E.reshape(n, n))).reshape(-1)
    return np.sort(np.linalg.eigvals(m*(m + 1)*(np.eye(n*n) - S)).real)


def thermisch(beta, alpha, rows, m=2):
    c = lambda x: 1/(1 + np.exp(np.clip(beta*(x + alpha), -700, 700)))
    W = np.array([c(r) - c(r + m) for r in rows])
    p_alt = np.array([sum((j + 1)*(c(r + j) - c(r + j + 1)) for j in range(m))/((m + 1)*(c(r) - c(r + m))) for r in rows])
    ziel = np.array([c(r) - c(r + m + 1) for r in range(rows[0] - 1, rows[-1] + 1)])*m/(m + 1)
    return W/W.sum(), 1 - p_alt, p_alt, ziel/W.sum()


ANALYSE = np.array([0.0, 1.9864, 1.9864, 2.0281])
Z77_ALT = np.array([0.0, 1.9802, 1.9802, 3.9589])
gew, u_th, p_alt, ziel = thermisch(4*np.pi, 0.25, [-2, -1])
ev_alt = petz_spektrum(gew, p_alt, 2) + 0.0
ev_neu = np.round(petz_spektrum(gew, u_th, 2), 12) + 0.0
ZAHLEN["petz_m2_beta4pi_alpha1/4"] = {"z77_alte_zuordnung": np.round(ev_alt, 6).tolist(),
                                      "richtige_zuordnung": np.round(ev_neu, 6).tolist(), "analyse": ANALYSE.tolist()}
say(f"   alte Zuordnung {np.round(ev_alt, 4).tolist()}, richtige Zuordnung {np.round(ev_neu, 4).tolist()}, "
    f"Analyse {ANALYSE.tolist()}")
check("mit der alten Zuordnung reproduziert die Rechnung z77 Teil 3 wörtlich (0; 1,9802; 1,9802; 3,9589): gleiche Konstruktion",
      np.allclose(ev_alt, Z77_ALT, atol=1e-4))
check(f"mit der richtigen Zuordnung (u für Teilchen hoch) ergibt dieselbe Auslesung {np.round(ev_neu, 4).tolist()}: Null, ein "
      "Dublett knapp unter 2 und ein Singulett knapp über 2, alle innerhalb 0,05 von l(l+1) = (0, 2, 2, 2); das alte Singulett "
      "3,96 war ein Artefakt der Vertauschung",
      np.max(np.abs(ev_neu - np.array([0, 2, 2, 2]))) < 0.05 and abs(ev_neu[1] - ev_neu[2]) < 1e-9 and abs(ev_alt[3] - 2) > 1)
varianten = {}
for name, (b, a, rows, ref) in {
        "alpha = 1/4, Zeilen -2..-1 (z77)": (4*np.pi, 0.25, [-2, -1], None),
        "alpha = 3/4, Zeilen -2..-1": (4*np.pi, 0.75, [-2, -1], None),
        "alpha = 1/2, Zeilen -2..-1": (4*np.pi, 0.5, [-2, -1], None),
        "alpha = 1/4, Zeilen -3..-1": (4*np.pi, 0.25, [-3, -2, -1], None),
        "alpha = 1/4, Zeilen -2..0": (4*np.pi, 0.25, [-2, -1, 0], None),
        "alpha = 1/4, Zeilen -3..0": (4*np.pi, 0.25, [-3, -2, -1, 0], None),
        "alpha = 1/4, Petz bzgl. exaktem Zielzustand W_3": (4*np.pi, 0.25, [-2, -1], "ziel"),
        "alpha = 1/4, Spurreferenz I/2, thermisches u": (4*np.pi, 0.25, [-2, -1], "spur")}.items():
    g_, u_, _, z_ = thermisch(b, a, rows)
    if ref == "spur":
        ev = petz_spektrum(np.ones(len(rows))/len(rows), u_, 2)
    else:
        ev = petz_spektrum(g_, u_, 2, ref=z_ if ref == "ziel" else None)
    varianten[name] = ((np.round(ev[:4], 4) + 0.0).tolist(), float(np.max(np.abs(ev[:4] - ANALYSE))))
ZAHLEN["petz_varianten_(niedrigste_vier,_max_abweichung)"] = varianten
for k_, v_ in varianten.items():
    say(f"   {k_}: {v_[0]}, max. Abweichung von der Analyse {v_[1]:.4f}")
dmin = min(v_[1] for v_ in varianten.values())
f_ = lambda bb: petz_spektrum(*thermisch(bb, 0.25, [-2, -1])[:2], 2)[1] - 1.9864
b_fit = brentq(f_, 5, 60)
ev_fit = petz_spektrum(*thermisch(b_fit, 0.25, [-2, -1])[:2], 2)
ZAHLEN["beta_angepasst_an_dublett"] = {"beta/pi": b_fit/np.pi, "spektrum": np.round(ev_fit, 4).tolist()}
check(f"die Werte der Analyse (0; 1,9864; 1,9864; 2,0281) sind auch mit richtiger Zuordnung in keiner natürlichen Variante "
      f"reproduziert (kleinste maximale Abweichung {dmin:.4f}); selbst mit beta angepasst ans Dublett (beta = {b_fit/np.pi:.3f} pi) "
      f"wird das Singulett {ev_fit[3]:.4f} statt 2,0281", dmin > 1e-3 and abs(ev_fit[3] - 2.0281) > 1e-3)
ev34 = varianten["alpha = 3/4, Zeilen -2..-1"][0]
check("alpha = 3/4 gibt exakt dasselbe Spektrum wie alpha = 1/4 (Teilchen-Loch-Spiegelung x -> -x der Kovarianz)",
      np.allclose(ev34, varianten["alpha = 1/4, Zeilen -2..-1 (z77)"][0], atol=1e-4))
befund("Qualitativ stimmt das korrigierte Spektrum mit der Analyse überein (Null, Dublett knapp unter 2, Singulett knapp über 2, "
       "Summe nahe 6); die vierstelligen Werte hängen an einer nicht dokumentierten Konditionierung. Breitere Abschneidungen "
       "sind keine Sphärenauslesung mehr (zusätzliche Moden, niedrigere Werte)", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Mischung verschobener Vakua (Analyse 8 §1) und der Anteil 1/(m+1)")
dev_mix = mp.mpf(0)
with mp.workdps(120):                                       # c direkt: Differenzen nahe 1 brauchen Stellen
    for B in (mp.mpf(2), 4*mp.pi):
        for A in (mp.mpf(1)/4, mp.mpf(3)/4):
            wv = {e: w_mp(B, A, e) for e in range(-30, 400)}
            for x in range(-8, 9):
                cx = c_mp(B, A, x)
                dev_mix = max(dev_mix, abs(mp.fsum(wv[e] for e in range(x, 400)) - cx))    # c(x) = Summe_(e >= x) w(e)
            for m in range(2, 7):
                for r in range(-6, 5):
                    Wd = c_mp(B, A, r) - c_mp(B, A, r + m)                                   # direkt
                    schale = [e for e in range(-30, 30) if r <= e <= r + m - 1]             # Kanten, deren Schale r enthält
                    Wmix = mp.fsum(wv[e] for e in schale)
                    u_mix = mp.fsum(wv[e]*mp.mpf(r + m - e)/(m + 1) for e in schale)/Wmix   # z76: (k_e + 1)/(m+1)
                    u_dir = mp.fsum((m - j)*(c_mp(B, A, r + j) - c_mp(B, A, r + j + 1)) for j in range(m))/((m + 1)*Wd)
                    # Ausgang: Mischung der verschobenen z76-Kanäle auf Ziel r' = r gegen m/(m+1) W_(m+1)(r)
                    aus = mp.fsum(wv[e]*mp.mpf(m)/(m + 1) for e in range(-30, 30) if r <= e <= r + m)
                    Wn = c_mp(B, A, r) - c_mp(B, A, r + m + 1)
                    dev_mix = max(dev_mix, abs(Wmix - Wd), abs(u_mix - u_dir), abs(aus - mp.mpf(m)/(m + 1)*Wn))
ZAHLEN["mischung_max_abweichung"] = float(dev_mix)
check(f"exakt bis {float(dev_mix):.0e} (beta = 2, 4 pi; alpha = 1/4, 3/4; m = 2..6): c(x) = Summe_e w(e) [x <= e] (Fermiseen "
      "mit Kante zwischen e und e+1), W_m(r) = Mischung der Vakuumpaargewichte, u_m(r) = Mittel der verschobenen z76-Gewichte "
      "(r + m - e)/(m+1) unter P(e | r), und der Ausgang der verschobenen z76-Kanäle ist m/(m+1) W_(m+1): die thermische "
      "Zustandserhaltung aus Teil 1 ist die w-Mittelung der Vakuumidentität", dev_mix < mp.mpf(10)**-40)


def kraus76(m):
    Ap, Ah = np.zeros((m + 1, m)), np.zeros((m + 1, m))
    for k in range(m):
        Ap[k + 1, k] = np.sqrt((k + 1)/(m + 1))
        Ah[k, k] = np.sqrt((m - k)/(m + 1))
    return Ap, Ah


def block_petz(m, gew_e):
    """Kantenaufgelöster Kanal: direkte Summe der z76-Kanäle, Referenz Summe_e w(e) I/m, auf blockdiagonalen Matrizen."""
    E = len(gew_e)
    Ap, Ah = kraus76(m)
    KA, KH = sla.block_diag(*[Ap]*E), sla.block_diag(*[Ah]*E)
    sg = np.concatenate([np.full(m, g/m) for g in gew_e])
    Phi = lambda X: KA @ X @ KA.T + KH @ X @ KH.T
    img = np.diag(Phi(np.diag(sg)))
    isq, sq = np.diag(1/np.sqrt(img)), np.diag(np.sqrt(sg))
    Rec = lambda Y: sq @ (KA.T @ isq @ Y @ isq @ KA + KH.T @ isq @ Y @ isq @ KH) @ sq
    basis = [(e*m + i, e*m + j) for e in range(E) for i in range(m) for j in range(m)]
    S = np.zeros((len(basis), len(basis)))
    leck = 0.0
    for a, (i, j) in enumerate(basis):
        X = np.zeros((E*m, E*m))
        X[i, j] = 1
        Y = Rec(Phi(X))
        S[:, a] = [Y[k, l] for (k, l) in basis]
        leck = max(leck, np.abs(Y).sum() - np.abs(S[:, a]).sum())
    ev = np.sort(np.linalg.eigvals(m*(m + 1)*(np.eye(len(basis)) - S)).real)
    expect = np.sort(np.concatenate([[l*(l + 1)]*(2*l + 1)*E for l in range(m)]))
    return np.max(np.abs(ev - expect)), leck


ok_blk = True
for beta in (2.0, 4*np.pi):
    wl = [float(w_mp(mp.mpf(beta), mp.mpf(1)/4, e)) for e in range(-3, 2)]
    for m in (2, 3, 4):
        for gw in (wl, rng.random(5)):
            d_, l_ = block_petz(m, gw)
            ok_blk &= d_ < 1e-9 and l_ < 1e-12
check("Petz je Kantenkomponente (direkte Summe der verschobenen z76-Kanäle, fünf Kanten, Gewichte w(e) bei beta = 2, 4 pi und "
      "zufällig; m = 2, 3, 4): m(m+1)(I - R Phi) hat exakt l(l+1) mit Vielfachheit 5(2l+1), unabhängig von Gewichten und "
      "Temperatur; kein Übersprechen zwischen den Kanten", ok_blk)

n2, s0 = 10, 5
cs2 = fock(n2)
occ = np.array([(cs2[s].T @ cs2[s]).diagonal() for s in range(n2)]).astype(int)        # occ[s, Index]
beta_g, alpha_g = 2.0, 0.25
cg = lambda r: 1/(1 + np.exp(beta_g*(r + alpha_g)))
cvec = np.array([cg(s - s0) for s in range(n2)])
P_qf = np.prod(np.where(occ == 1, cvec[:, None], 1 - cvec[:, None]), axis=0)          # quasifrei: Produkt in Besetzungsbasis
P_mix = np.zeros(2**n2)
for e in range(-1, n2):                                                                 # Fermisee: Plätze 0..e besetzt
    idx = np.where((occ == (np.arange(n2)[:, None] <= e).astype(int)).all(axis=0))[0][0]
    P_mix[idx] = (1 - cvec[0]) if e == -1 else (cvec[n2 - 1] if e == n2 - 1 else cvec[e] - cvec[e + 1])
ew = lambda P, O: float(P @ O.diagonal())
nops = [cs2[s].T @ cs2[s] for s in range(n2)]
ok_ge, cov_mix_max = abs(P_mix.sum() - 1) < 1e-14 and abs(P_qf.sum() - 1) < 1e-14, 0.0
for x in range(n2):
    ok_ge &= abs(ew(P_mix, nops[x]) - cvec[x]) < 1e-14 and abs(ew(P_qf, nops[x]) - cvec[x]) < 1e-14
    for y in range(x + 1, n2):
        cov_m = ew(P_mix, nops[x] @ nops[y]) - cvec[x]*cvec[y]
        cov_q = ew(P_qf, nops[x] @ nops[y]) - cvec[x]*cvec[y]
        ok_ge &= abs(cov_m - cvec[y]*(1 - cvec[x])) < 1e-14 and abs(cov_q) < 1e-14
        cov_mix_max = max(cov_mix_max, cov_m)
for m in (1, 2, 3):
    for r in range(0, n2 - m):
        Bo = cs2[r + m].T @ cs2[r]
        BdB, Kom = Bo.T @ Bo, Bo.T @ Bo - Bo @ Bo.T
        ok_ge &= abs(ew(P_mix, BdB) - (cvec[r] - cvec[r + m])) < 1e-14
        ok_ge &= abs(ew(P_qf, BdB) - cvec[r]*(1 - cvec[r + m])) < 1e-14
        ok_ge &= abs(ew(P_mix, Kom) - (cvec[r] - cvec[r + m])) < 1e-14 and abs(ew(P_qf, Kom) - (cvec[r] - cvec[r + m])) < 1e-14
ZAHLEN["kovarianz_mischung_max_(beta=2)"] = cov_mix_max
check(f"Gegenbeispiel auf dem Fockraum (10 Moden, beta = 2, alpha = 1/4): Mischung und quasifreier Zustand haben dieselben "
      f"c(x) und dieselben Kommutatorgewichte w([B^+, B]) = W_m(r), aber für x < y ist Cov(n_x, n_y) = c(y)(1 - c(x)) in der "
      f"Mischung (bis {cov_mix_max:.3f}) und 0 im quasifreien Zustand; w(B^+ B) ist W_m(r) gegen c(r)(1 - c(r+m)): der volle "
      "thermische Zustand ist nicht die Mischung", ok_ge and cov_mix_max > 0.05)

n = 14
cs = fock(n)
dim = 2**n
Jm1 = sum(cs[s + 1].T @ cs[s] for s in range(n - 1))
Nop = sum(cs[s].T @ cs[s] for s in range(n))
leer = np.zeros(dim)
leer[0] = 1.0


def see(e):
    v = leer.copy()
    for s in range(e, -1, -1):
        v = cs[s].T @ v
    return v


def Bop(r, m):
    return cs[r + m].T @ cs[r]


muster, ok_fock = {}, True
for m in range(2, 6):
    pro_kante = []
    for e in (5, 6, 7):
        F = see(e)
        tot = multi = 0.0
        an = {}
        for r in range(e - m + 1, e + 1):
            v = Bop(r, m) @ F
            ok_fock &= abs(np.linalg.norm(v) - 1) < 1e-12 and abs(v @ (Nop @ v) - (e + 1)) < 1e-12
            Jv = Jm1 @ v
            proj = 0.0
            for rp in range(e - m, e + 1):
                amp = (Bop(rp, m + 1) @ F) @ Jv
                an[rp - e] = an.get(rp - e, 0.0) + amp**2
                proj += amp**2
                soll = 1.0 if rp == r else (-1.0 if rp == r - 1 else 0.0)                 # Teilchen hoch +1, Loch runter -1
                ok_fock &= abs(amp - soll) < 1e-12
            nrm = Jv @ Jv
            innen = 1.0 if (r != e and r + m != e + 1) else 0.0                          # Paarbildung an der Kante
            ok_fock &= abs(nrm - proj - innen) < 1e-12
            tot += nrm
            multi += nrm - proj
        pro_kante.append((round(tot, 10), round(multi, 10), tuple(round(an[t], 10) for t in sorted(an))))
    ok_fock &= len(set(pro_kante)) == 1
    muster[m] = pro_kante[0]
check("Fockraum (14 Moden, Kanten e = 5, 6, 7; m = 2..5): J_-1 psi*_(r+m) psi_r |FS_e> = (Teilchen hoch, Amplitude +1) - (Loch "
      "runter, -1) + (neues Paar an der Kante, Norm 1 für innere Paare, 0 an den Schalenenden); Ladung e+1 je Kante "
      "(Superauswahl: verschiedene Kanten sind verschiedene Ladungssektoren); das Muster ist translationsinvariant", ok_fock)

anteile = {}
ok_a, ok_b, ok_c = True, True, True
for m in range(2, 6):
    tot, multi, an = muster[m]
    f_a = multi/tot
    f_b = sum(2 - a for a in an)/(2*(m + 1))
    ok_a &= abs(f_a - (m - 2)/(3*m - 2)) < 1e-12 and (abs(f_a - 1/(m + 1)) > 1e-3 if m != 4 else True)
    ok_b &= abs(f_b - 1/(m + 1)) < 1e-12
    f_bt, f_ct_neu, f_ct_alt = [], [], []
    for beta in (2.0, 4*np.pi):
        for alpha in (0.25, 0.75):
            wf = {e: float(w_mp(mp.mpf(beta), mp.mpf(alpha), e)) for e in range(-50, 51)}
            ank = {t: an[t + m] for t in range(-m, 1)}
            ankunft = sum(wf[e]*ank[rp - e] for rp in range(-40, 41) for e in range(rp, rp + m + 1))
            bedarf = sum(2*wf[e] for rp in range(-40, 41) for e in range(rp, rp + m + 1))
            f_bt.append(1 - ankunft/bedarf)
            Wm_ = {r: sum(wf[r + j] for j in range(m)) for r in range(-45, 45)}
            Wn_ = {r: sum(wf[r + j] for j in range(m + 1)) for r in range(-45, 45)}
            u_ = {r: sum((m - j)*wf[r + j] for j in range(m))/((m + 1)*Wm_[r]) for r in range(-45, 45)}
            for zu, liste in ((u_, f_ct_neu), ({r: 1 - u_[r] for r in u_}, f_ct_alt)):
                aus = sum(zu[r]*Wm_[r] + (1 - zu[r + 1])*Wm_[r + 1] for r in range(-40, 41))
                liste.append(1 - aus/sum(Wn_[r] for r in range(-40, 41)))
    ok_b &= max(abs(x - 1/(m + 1)) for x in f_bt) < 1e-12
    ok_c &= max(abs(x - 1/(m + 1)) for x in f_ct_neu + f_ct_alt) < 1e-12
    Ap, Ah = kraus76(m)
    D = np.zeros((m + 1, m))
    for k in range(m):
        D[k + 1, k], D[k, k] = 1.0, -1.0
    f_d = 1 - (np.sum(Ap*D)**2/np.sum(Ap*Ap) + np.sum(Ah*D)**2/np.sum(Ah*Ah))/np.sum(D*D)
    anteile[m] = {"(a)_mehrpaaranteil": round(float(f_a), 6) + 0.0, "(b)_pfaddefizit_schalenenden": round(float(f_b), 6),
                  "(b)_thermisch_ueber_alle_kanten": round(float(max(f_bt)), 12), "(c)_gewichtsbilanz_neu/alt":
                  [round(float(max(f_ct_neu)), 12), round(float(max(f_ct_alt)), 12)],
                  "(d)_HS_ausserhalb_kraus": round(float(f_d), 6), "1/(m+1)": round(1/(m + 1), 6)}
    ok_a &= abs(f_d - 1/(m + 1)) > 1e-3
ZAHLEN["anteil_1/(m+1)_definitionen"] = anteile
for k_, v_ in anteile.items():
    say(f"   m = {k_}: {v_}")
check("Anteil 1/(m+1), Zähldefinitionen: (b) die nativen Zweipfadübergänge erreichen die beiden Schalenenden von H_(m+1) nur "
      "über einen Pfad; das fehlende Quadratgewicht ist je Kante 2 von 2(m+1), thermisch über alle Kanten summiert exakt "
      "1/(m+1) (beta = 2, 4 pi; alpha = 1/4, 3/4); (c) die Gewichtsbilanz m gegen m+1 gibt 1/(m+1) für jeden spurerhaltenden "
      "Kanal, auch für die falsche z77-Zuordnung", ok_b and ok_c)
check("Anteil 1/(m+1), dynamische Definitionen nicht: (a) der Mehrpaaranteil von |J_-1 Paar>^2 (Paarbildung an der Kante) ist "
      "(m-2)/(3m-2) (nur bei m = 4 zufällig 1/5), (d) der Hilbert-Schmidt-Anteil des nativen J_-1 außerhalb der Kraus-Hülle "
      "ist ein anderer Wert", ok_a)
befund("Die Mischungslesart der Analyse 8 hält für die Einteilchendaten exakt (c, W_m, u_m, Kanalausgang) und kantenaufgelöst "
       "gibt Petz l(l+1) bei jeder Temperatur. Sie ist aber eine Darstellung der Einteilchendaten, nicht des Zustands: Der "
       "quasifreie Zustand hat dieselben c und Kommutatorgewichte, aber keine Kantenkorrelationen. Die über die Kante "
       "vergröberte Auslesung (Teil 2) ist nicht exakt l(l+1)", True)
befund("Die Zahl 1/(m+1) ist ohne Definition der Analyse nur als Zählidentität reproduziert (2m native Pfade auf 2(m+1) "
       "Plätze, bzw. Gewicht m gegen m+1); als dynamische Aussage über J_-1 (Mehrpaaranteil) gilt sie nicht", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Relative Zustände aus z78 und ihr Überlapp (Analysen 7 §3, 8 §2)")
J4 = np.block([[0*I2, -I2], [I2, 0*I2]])
symb = []
for i in range(4):
    for j in range(i, 4):
        S = np.zeros((4, 4))
        S[i, j] = S[j, i] = 1
        symb.append(S)
dim_inv = 10 - np.linalg.matrix_rank(np.array([(S @ J4 - J4 @ S).reshape(-1) for S in symb]).T)
realform = lambda X: np.block([[X.real, -X.imag], [X.imag, X.real]])
rank_h = np.linalg.matrix_rank(np.array([realform(s).reshape(-1) for s in SIG]).T)
ok_q = dim_inv == 4 and rank_h == 4
ok_q &= all(np.allclose(realform(s), realform(s).T) and np.allclose(realform(s) @ J4, J4 @ realform(s)) for s in SIG)
for _ in range(20):
    cc = rng.normal(size=2) + 1j*rng.normal(size=2)
    X = herm(rng.normal(size=4))
    vv = np.concatenate([cc.real, cc.imag])
    ok_q &= abs(vv @ realform(X) @ vv - (cc.conj() @ X @ cc).real) < 1e-12
tt, xx, yy, zz = sp.symbols("t x y z", real=True)
Zs = sp.Matrix([[tt + zz, xx - sp.I*yy], [xx + sp.I*yy, tt - zz]])
ok_q &= sp.expand(Zs.det() - (tt**2 - xx**2 - yy**2 - zz**2)) == 0
check(f"die U(1)-invarianten reellen Quadriken auf C^2 = R^4 (Kommutante der komplexen Struktur in den symmetrischen 4x4) "
      f"haben Dimension {dim_inv} und sind genau Re-Formen c^+ X c mit X in Herm(2) (Rang {rank_h}); "
      "det(t + x.sigma) = t^2 - x^2 - y^2 - z^2", ok_q)


def fock4(nn):
    a = sps.diags(np.sqrt(np.arange(1, nn)), 1, format="csr")
    Id = sps.identity(nn, format="csr")
    ops = []
    for k in range(4):
        mats = [Id]*4
        mats[k] = a
        op = mats[0]
        for M in mats[1:]:
            op = sps.kron(op, M, format="csr")
        ops.append(op)
    return ops                              # a1, a2, b1, b2


nf = 12
a1_, a2_, b1_, b2_ = fock4(nf)
aa, bb = [a1_, a2_], [b1_, b2_]
err_ov, err_s2, ovs = 0.0, 0.0, []
for qv in (0.2, float(np.exp(-2*np.pi)), 0.35):
    om = np.zeros([nf]*4)
    for i in range(nf):
        for j in range(nf):
            om[i, j, i, j] = (1 - qv**2)*qv**(i + j)
    om = om.reshape(-1).astype(complex)
    kap = np.sqrt((1 - qv)/(1 + qv))
    ct = [kap*(aa[i] + bb[i].T) for i in range(2)]
    PX = lambda X: 0.5*sum(ct[i].conj().T @ ct[j]*X[i, j] for i in range(2) for j in range(2))
    for _ in range(4):
        X, Y = herm(rng.normal(size=4)*0.3), herm(rng.normal(size=4)*0.3)
        sx, sy = spla.expm_multiply(1j*PX(X), om), spla.expm_multiply(1j*PX(Y), om)
        ov = np.vdot(sx, sy)
        Kf = 1/np.linalg.det(I2 - 0.5j*(Y - X))
        err_ov = max(err_ov, abs(ov - Kf))
        err_s2 = max(err_s2, abs(4*(1 - (1/ov).real) - np.linalg.det(Y - X).real))
        ovs.append(abs(ov))
ZAHLEN["relativ_ueberlapp_max_fehler"] = err_ov
ZAHLEN["relativ_s2_max_fehler"] = err_s2
check(f"auf dem Fockraum (vier Moden a_1, a_2, b_1, b_2, je 12 Stufen; q = 0,2, e^-2pi, 0,35; je vier zufällige X, Y) gilt "
      f"<X|Y> = <Omega| e^(-i P(X)) e^(i P(Y)) |Omega> = 1/det(I - (i/2)(Y - X)) mit maximalem Fehler {err_ov:.1e}, und "
      f"4(1 - Re 1/<X|Y>) = det(Y - X) bis {err_s2:.1e}: unabhängig von q (c ist standard-komplex-gaußsch)",
      err_ov < 1e-6 and err_s2 < 1e-5)
Ks = 1/(sp.eye(2) - sp.I/2*Zs).det()
r2s = xx**2 + yy**2 + zz**2
ok_k = sp.simplify(Ks - 4/(r2s - (tt + 2*sp.I)**2)) == 0
ok_k &= sp.simplify(4*(1 - sp.re(sp.expand(1/Ks))) - Zs.det()) == 0
ok_cay, normmax = True, 0.0
for _ in range(50):
    X, Y = herm(rng.normal(size=4)), herm(rng.normal(size=4))
    K = 1/np.linalg.det(I2 - 0.5j*(Y - X))
    ZX, ZY = np.linalg.inv(X + 2j*I2) @ X, np.linalg.inv(Y + 2j*I2) @ Y          # Cayley (z74) der Röhrenpunkte X + iI
    normmax = max(normmax, np.linalg.norm(ZX, 2), np.linalg.norm(ZY, 2))
    rhs = 16/np.linalg.det(I2 - ZY @ ZX.conj().T)/(np.linalg.det(Y + 2j*I2)*np.conj(np.linalg.det(X + 2j*I2)))
    ok_cay &= abs(K - rhs) < 1e-12*abs(K)
check("symbolisch: K = 1/det(I - (i/2)Z) = 4/(r^2 - (t + 2i)^2), also 16 pi^2 mal die masselose 4D-Wightmanfunktion "
      "1/(4 pi^2 (r^2 - t^2)) bei der imaginären Zeitverschiebung t -> t + 2i, und s^2 = 4(1 - Re 1/K) = det(Z) = t^2 - r^2; "
      f"numerisch (50 Paare) ist K = 16 det(I - Z_Y Z_X^+)^-1/(det(Y + 2i) conj det(X + 2i)) mit Z_X = Cayley(X + iI), "
      f"|Z| <= {normmax:.3f} < 1: der Kern nu = 1 aus z74 an den Röhrenpunkten X + iI", ok_k and ok_cay and normmax < 1)
befund("Neu als Verbindung: Die Überlappe der relativen Zustände aus z78 sind exakt der Kern nu = 1 det(I - Z W^+)^-1 aus z74 in "
       "Röhrenkoordinaten (Cayley der Punkte X + iI, bis auf den Automorphiefaktor), und der Minkowskiabstand det(Y - X) ist "
       "aus dem Überlapp auslesbar. Offen bleibt wie in z78, dass der Quarterindex das Dublett der Naht trägt und die "
       "Lorentzwirkung darauf nativ ist", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Die Wurzel des Kerns ist nicht positiv (Analyse 8 §4)")
nu = sp.symbols("nu")
zs_ = sp.symbols("z0:4", real=True)
detm = (1 - sp.I*zs_[0]/2)**2 + (zs_[1]**2 + zs_[2]**2 + zs_[3]**2)/4         # det(I - i Z/2), Z = z.sigma
box = lambda f: sp.diff(f, zs_[0], 2) - sum(sp.diff(f, zs_[k], 2) for k in (1, 2, 3))
phi_nu = detm**(-nu)
b1 = box(phi_nu)
nullp = {zz_: 0 for zz_ in zs_}
E1 = sp.factor(sp.simplify(-b1.subs(nullp)))
E2 = sp.factor(sp.simplify(box(b1).subs(nullp)))
ZAHLEN["momente"] = {"E[det Q]": str(E1), "E[(det Q)^2]": str(E2)}
check(f"aus E[exp(i tr(QZ)/2)] = det(I - iZ/2)^-nu folgt (Reihenentwicklung, det Q = q_0^2 - |q|^2 als Box auf der "
      f"charakteristischen Funktion) E[det Q] = {E1}, E[(det Q)^2] = {E2}; bei nu = 1/2 also {E1.subs(nu, sp.Rational(1, 2))} "
      f"und {E2.subs(nu, sp.Rational(1, 2))} < 0: kein positives Maß; bei nu = 1 beide 0 (Q = c c^+ hat Rang eins)",
      E2.subs(nu, sp.Rational(1, 2)) == sp.Rational(-3, 16) and E2.subs(nu, 1) == 0 and E1.subs(nu, 1) == 0)


def Knu(Z, nuv):
    lam = np.linalg.eigvalsh(Z)
    return np.prod((1 - 0.5j*lam)**(-nuv))                       # Hauptzweig je Eigenwert, stetig


h = 0.05
punkte, koeff = [np.zeros(4)], [-4/h**2]
for mu, et in enumerate([1, -1, -1, -1]):
    for sgn in (1, -1):
        v = np.zeros(4)
        v[mu] = sgn*h
        punkte.append(v)
        koeff.append(-et/h**2)
koeff = np.array(koeff)
form = {}
for nuv in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
    G = np.array([[Knu(herm(punkte[j] - punkte[i]), nuv) for j in range(9)] for i in range(9)])
    form[nuv] = float((koeff @ G @ koeff).real)
pz = rng.normal(size=(60, 4))
mineig = {}
for nuv in (0.5, 1.0):
    G = np.array([[Knu(herm(pz[j] - pz[i]), nuv) for j in range(60)] for i in range(60)])
    ew_ = np.linalg.eigvalsh((G + G.conj().T)/2)
    mineig[nuv] = (float(ew_.min()), float(ew_.max()))
ZAHLEN["gram_schablone_h0.05"] = form
ZAHLEN["gram_zufall_60_min_max"] = {str(k_): v_ for k_, v_ in mineig.items()}
check(f"Gramtest: die Neun-Punkt-Schablone der Form (det Q)^2 (h = 0,05) gibt {form[0.5]:.4f} -> -3/16 bei nu = 1/2 und ist "
      f"negativ für nu = 1/4, 1/2, 3/4, nichtnegativ für nu = 1, 3/2, 2; 60 Zufallspunkte: kleinster Eigenwert "
      f"{mineig[0.5][0]:.3f} für K^(1/2), {mineig[1.0][0]:.1e} >= 0 für K",
      all(form[v_] < 0 for v_ in (0.25, 0.5, 0.75)) and all(form[v_] >= 0 for v_ in (1.0, 1.5, 2.0))
      and abs(form[0.5] + 3/16) < 0.01 and mineig[0.5][0] < -1e-3 and mineig[1.0][0] > -1e-10*mineig[1.0][1])
treffer = {f: (HERE / f).read_text().count("Wallach") for f in ("z40_gemeinsamer_ursprung.py", "z41_altbestand_spin2.py",
                                                               "z74_tprime_m_kern.py")}
ZAHLEN["wallach_im_bestand"] = treffer
befund(f"Das ist die Wallach-Menge der 2x2-Matrixkugel (nu in {{0, 1}} oder nu > 1, also {{0}} u [1, oo)); nu = 1 ist der "
       f"kleinste nichttriviale Punkt. Schon im Bestand: 'Wallach' in z40 ({treffer['z40_gemeinsamer_ursprung.py']}x, "
       f"Ikosaeder-Test negativ bei nu = 1/2), z41 ({treffer['z41_altbestand_spin2.py']}x, Wallach-Menge = 4D-"
       f"Unitaritätsschranke) und z74 ({treffer['z74_tprime_m_kern.py']}x, (nu)_m >= 0 nur auf der Wallach-Menge); das "
       "Moment -3/16 ist eine neue Darstellung, kein neuer Befund", all(v_ > 0 for v_ in treffer.values()))

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: Zusammengesetzter Spin 2 eines freien masselosen 4D-Skalars (Analyse 7 §5)")


def boost(v, nvec):
    nvec = np.asarray(nvec, dtype=float)/np.linalg.norm(nvec)
    g = 1/np.sqrt(1 - v**2)
    L = np.eye(4)
    L[0, 0] = g
    L[0, 1:] = L[1:, 0] = g*v*nvec
    L[1:, 1:] += (g - 1)*np.outer(nvec, nvec)
    return L


def proj2(p):
    pl = ETA @ p
    th = ETA - np.outer(pl, pl)/(p @ ETA @ p)
    P2 = 0.5*(np.einsum("ac,bd->abcd", th, th) + np.einsum("ad,bc->abcd", th, th)) - np.einsum("ab,cd->abcd", th, th)/3
    return P2, th


def rho_spin(s, v=0.0, xi=0.0, nth=24, nph=48):
    """Zweiteilchen-Phasenraum (masselos, Integral dPhi_2 = 1/(8 pi)) von M Pi M; M = <p1 p2|T_mu nu|0>, T mit Verbesserung xi."""
    L = boost(v, (0.3, -0.5, 0.8))
    xg, wg = leggauss(nth)
    p = L @ np.array([np.sqrt(s), 0, 0, 0])
    P2, th = proj2(p)
    pl = ETA @ p
    a2 = a0 = 0.0
    for ct_, wt in zip(xg, wg):
        st = np.sqrt(1 - ct_**2)
        for f in 2*np.pi*np.arange(nph)/nph:
            nn = np.array([st*np.cos(f), st*np.sin(f), ct_])
            p1, p2 = L @ (np.sqrt(s)/2*np.array([1, *nn])), L @ (np.sqrt(s)/2*np.array([1, *(-nn)]))
            M = np.outer(ETA @ p1, ETA @ p2) + np.outer(ETA @ p2, ETA @ p1) - ETA*(p1 @ ETA @ p2)
            M = M - 2*xi*(np.outer(pl, pl) - ETA*(p @ ETA @ p))
            Mu = ETA @ M @ ETA
            a2 += wt/nph*np.einsum("ab,abcd,cd->", Mu, P2, Mu)
            a0 += wt/nph*np.einsum("ab,ab->", Mu, th)**2/3
    return a2/2/5/(8*np.pi), a0/2/(8*np.pi)


pt = np.array([1.3, 0.2, -0.4, 0.5])
P2t, tht = proj2(pt)
P2m = P2t.reshape(16, 16)
Gm = np.kron(ETA, ETA)
ok_pr = np.allclose(P2m @ Gm @ P2m, P2m) and abs(np.trace(P2m @ Gm) - 5) < 1e-12
ok_pr &= np.allclose(np.einsum("abcd,c->abd", P2t, ETA @ ETA @ pt), 0)
ok_pr &= np.allclose(np.einsum("abcd,cd->ab", P2t, np.linalg.inv(ETA)), 0)
konst, ok_rho = 1/(240*np.pi), True
werte = {}
for s in (0.01, 0.1, 1.0, 10.0):
    r2, r0 = rho_spin(s)
    r2b, _ = rho_spin(s, v=0.7)
    r2i, r0i = rho_spin(s, xi=1/6)
    werte[s] = (r2/s**2, r2b/s**2, r0/s**2, r0i/s**2)
    ok_rho &= max(abs(x/konst - 1) for x in (r2/s**2, r2b/s**2, r2i/s**2)) < 1e-10
    ok_rho &= abs(r0/s**2*24*np.pi - 1) < 1e-10 and abs(r0i) < 1e-12*s**2
xg, wg = leggauss(8)
atom = {}
for eps in (1.0, 0.1, 0.01):
    ss = eps*(xg + 1)/2
    atom[eps] = float(eps/2*sum(wt*rho_spin(s_, nth=8, nph=12)[0] for s_, wt in zip(ss, wg)))
ok_atom = all(abs(atom[e_]/(e_**3/(720*np.pi)) - 1) < 1e-8 for e_ in atom)
ZAHLEN["spin2_dichte_rho2/s^2_(ruhe,_boost,_spin0_minimal,_spin0_verbessert)"] = {str(k_): v_ for k_, v_ in werte.items()}
ZAHLEN["spin2_gewicht_in_[0,eps]"] = {str(k_): v_ for k_, v_ in atom.items()}
check("Spin-2-Projektor Pi^(2)(p) ist idempotent, hat Spur 5 und annulliert p_mu und eta_mu nu (Verbesserungsterm)", ok_pr)
check("Zweiteilchen-Phasenraum (s = 0,01..10, Ruhesystem und Boost v = 0,7): die Spin-2-Spektraldichte von <T T> ist "
      "rho_2(s) = s^2/(240 pi) theta(s), Lorentz-invariant und unabhängig von der Verbesserung; Spin 0: s^2/(24 pi) für das "
      "minimale T, 0 für das verbesserte (xi = 1/6)", ok_rho)
check(f"kein isolierter masseloser Pol: das Spektralgewicht in [0, eps] ist eps^3/(720 pi) (eps = 1; 0,1; 0,01: "
      f"{', '.join(f'{v_:.2e}' for v_ in atom.values())}) und geht gegen 0, ohne delta(s)-Anteil", ok_atom)
befund("Deckt sich mit z75 (T(z) der Naht gibt kein Graviton; im Sugawara-Bild landet T bei phi^2 und dem zusammengesetzten "
       "T_mu nu mit Delta = 4, kein (2,0) bei Delta = 3) und mit v1 H6h21 ('kein zusammengesetzter masseloser Spin-2-Pol'): "
       "zusammengesetzte Spin-2-Zustände eines freien Feldes bilden ein Kontinuum s^2, ein Teilchen fehlt; Weinberg-Witten "
       "verbietet ohnehin ein masseloses Spin-2-Teilchen mit kovariantem erhaltenem T_mu nu. Nichts Neues", True)

# ------------------------------------------------------------------ Teil 7
say("\nTeil 7: Spin-2-Ebene-Wellen aus c (Analyse 8 §3)")


def feld(cv, x):
    p = vec(np.outer(cv, cv.conj()))
    return np.einsum("a,b,c,d->abcd", cv, cv, cv, cv)*np.exp(1j*(p @ ETA @ x)), p


def weyl_res(cv, x, ops):
    ph, p = feld(cv, x)
    res = sum(np.einsum("xa,abcd->xbcd", ops[mu], 1j*(ETA @ p)[mu]*ph) for mu in range(4))
    return np.abs(res).max()/(np.abs(ph).max()*np.abs(p).max())


D_rechts = [SIG[0], SIG[1], SIG[2], SIG[3]]          # d_0 + sigma.grad: auf e^(ip.x) gibt das i adj(P), adj(P) c = 0
D_links = [SIG[0], -SIG[1], -SIG[2], -SIG[3]]        # Gegenoperator: i P, P c = |c|^2 c != 0
res_r, res_l, err_h, err_lg, err_lor = 0.0, 1e9, 0.0, 0.0, 0.0
for _ in range(100):
    cv = rng.normal(size=2) + 1j*rng.normal(size=2)
    x = rng.normal(size=4)
    res_r = max(res_r, weyl_res(cv, x, D_rechts))
    res_l = min(res_l, weyl_res(cv, x, D_links))
    th_ = rng.uniform(0, 2*np.pi)
    ph0, p0 = feld(cv, x)
    ph1, p1 = feld(np.exp(1j*th_)*cv, x)
    err_h = max(err_h, np.abs(ph1 - np.exp(4j*th_)*ph0).max()/np.abs(ph0).max() + np.abs(p1 - p0).max())
    chi = rng.uniform(0, 2*np.pi)
    nhat = p0[1:]/np.linalg.norm(p0[1:])
    Mlg = sla.expm(-0.5j*chi*sum(nhat[k]*SIG[k + 1] for k in range(3)))
    ph2, p2 = feld(Mlg @ cv, x)
    err_lg = max(err_lg, np.abs(ph2 - np.exp(-2j*chi)*ph0).max()/np.abs(ph0).max() + np.abs(p2 - p0).max())
    Gc = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
    M = Gc/np.sqrt(np.linalg.det(Gc))
    x2 = vec(M @ herm(x) @ M.conj().T)
    ph3, _ = feld(M @ cv, x2)
    rhs = np.einsum("ae,bf,cg,dh,efgh->abcd", M, M, M, M, ph0)
    err_lor = max(err_lor, np.abs(ph3 - rhs).max()/np.abs(rhs).max())
ZAHLEN["spin2_ebene_welle"] = {"weyl_residuum": res_r, "gegenoperator_min": res_l, "phase": err_h, "kleine_gruppe": err_lg,
                               "lorentz": err_lor}
check(f"phi_ABCD = c_A c_B c_C c_D e^(i p.x) mit p = c c^+ (null) löst die masselose Spin-2-Gleichung (d_0 + sigma.grad) phi = 0 "
      f"(100 zufällige c, x; Residuum {res_r:.1e}; der Gegenoperator gibt >= {res_l:.2f}, der Test trennt)",
      res_r < 1e-13 and res_l > 0.1)
check(f"Helizität 2: c -> e^(i theta) c lässt p fest und gibt phi -> e^(4 i theta) phi; die Drehung um p um chi (kleine "
      f"Gruppe) gibt phi -> e^(-2 i chi) phi, also |h| = 2 (Fehler {max(err_h, err_lg):.1e}); Lorentzkovarianz "
      f"phi_(Mc)(Lambda x) = M^(x4) phi_c(x) für 100 zufällige M in SL(2,C) (Fehler {err_lor:.1e})",
      err_h < 1e-12 and err_lg < 1e-12 and err_lor < 1e-10)
befund("Standard-Spinor-Helizität (Penrose; masselose Felder als Produkte von Nullspinoren), linear. z65 hat bereits beide "
       "Helizitäten linear aus beiden Blättern (Penrose-Transformation, Bianchi exakt) und das nichtlineare googly-Problem; die "
       "Ebene Welle ist die Impulsform derselben linearen Lösungen. Lokale freie Felder über modulare Lokalisierung "
       "(Brunetti-Guido-Longo) sind ebenfalls Standard und keine TFPT-Dynamik. Nichts Neues", True)

# ------------------------------------------------------------------ Teil 8
say("\nTeil 8: Clifford aus der Vorzeichensymmetrie (Analyse 7 §6)")
A = SIG[3]
thetas = np.sort(np.append(np.linspace(0.05, np.pi - 0.05, 60), np.pi/2))    # 60 Gitterwinkel (ohne Mitte) plus pi/2
impl, ok_schur = [], True
for th_ in thetas:
    Bm = np.cos(th_)*SIG[3] + np.sin(th_)*SIG[1]
    ok_schur &= np.allclose(A @ Bm + Bm @ A, 2*np.cos(th_)*I2)
    Lm = []
    for i in range(2):
        for j in range(2):
            Eu = np.zeros((2, 2), dtype=complex)
            Eu[i, j] = 1
            Lm.append(np.concatenate([(Eu @ A + A @ Eu).ravel(), (Eu @ Bm - Bm @ Eu).ravel()]))
    sv = np.linalg.svd(np.array(Lm).T, compute_uv=False)
    if np.sum(sv < 1e-10) > 0:
        impl.append(float(th_))
ok_1d = True
for a_ in (1.0, -1.0):                                          # 1D-Darstellungen A = a, B = b: u a + a u = 0, u b - b u = 0
    for b_ in (1.0, -1.0):
        ok_1d &= np.linalg.matrix_rank(np.array([[2*a_], [b_ - b_]])) == 1           # nur u = 0: nie implementierbar
check("irreduzible Darstellungen von <A, B | A^2 = B^2 = I>: in 2D ist {A, B} = 2 cos(theta) I zentral (Schur); die "
      "Automorphie A -> -A, B -> B ist in derselben Darstellung genau dann unitär implementierbar, wenn theta = pi/2 "
      f"(Intertwinergleichungen, {len(thetas)} Winkel: lösbar nur bei {[round(t_/np.pi, 4) for t_ in impl]} pi), in 1D nie: "
      "also AB + BA = 0",
      ok_schur and ok_1d and len(impl) == 1 and abs(impl[0] - np.pi/2) < 1e-12)
Bc = SIG[1]
ok_cl = np.linalg.matrix_rank(np.array([M.ravel() for M in (I2, A, Bc, A @ Bc)])) == 4
ok_cl &= np.allclose((A @ Bc) @ (A @ Bc), -I2)


def bloch_z(psi):
    return psi[1]/psi[0] if abs(psi[0]) > 1e-12 else np.inf


kv = {}
for th_ in (np.pi/3, np.pi/2, 2*np.pi/3):
    Bm = np.cos(th_)*SIG[3] + np.sin(th_)*SIG[1]
    ea, va = np.linalg.eigh(A)
    eb, vb = np.linalg.eigh(Bm)
    za = [bloch_z(va[:, 1]), bloch_z(va[:, 0])]                  # +1 -> 0, -1 -> unendlich
    zb = [bloch_z(vb[:, 1]), bloch_z(vb[:, 0])]
    lam = (zb[0] - za[0])/(zb[1] - za[0])                         # Kreuzverhältnis mit z2 = unendlich (wie z71)
    pp = abs(np.vdot(va[:, 1], vb[:, 1]))**2
    kv[round(th_/np.pi, 4)] = (complex(np.round(lam, 12)), round(pp, 12))
    ok_cl &= abs(lam + (1 - pp)/pp) < 1e-12
lam_h, p_h = kv[0.5]
jinv = 256*(lam_h**2 - lam_h + 1)**3/(lam_h**2*(lam_h - 1)**2)
ok_cl &= abs(lam_h + 1) < 1e-12 and abs(p_h - 0.5) < 1e-12 and abs(jinv - 1728) < 1e-9
ZAHLEN["kreuzverhaeltnis_ueberlapp"] = {str(k_): [str(v_[0]), v_[1]] for k_, v_ in kv.items()}
check(f"bei AB + BA = 0: I, A, B, AB spannen M_2(C) (Cl_2(C) = M_2(C)), (AB)^2 = -I; die vier Marken (Eigenvektoren von A "
      f"und B) haben das Kreuzverhältnis -(1-p)/p mit p = cos^2(theta/2), also -1 (harmonisch, j = {jinv.real:.0f}, tau = i) "
      "und Überlapp p = 1/2", ok_cl)
rhos_inv = []
for name, gens in (("beide", (A, Bc)), ("nur_A->-A", (Bc,))):
    Lm = []
    for i in range(2):
        for j in range(2):
            Eu = np.zeros((2, 2), dtype=complex)
            Eu[i, j] = 1
            Lm.append(np.concatenate([(g @ Eu @ g - Eu).ravel() for g in gens]))
    sv = np.linalg.svd(np.array(Lm).T, compute_uv=False)
    rhos_inv.append(int(np.sum(sv < 1e-10)))
bwert, tmod = 0.4, 1.0
rho_b = (I2 + bwert*Bc)/2
sig_t = sla.expm(1j*tmod*sla.logm(rho_b)) @ A @ sla.expm(-1j*tmod*sla.logm(rho_b))
ok_st = rhos_inv == [1, 2] and np.linalg.norm(sig_t - A) > 0.1
check(f"Zustände: invariant unter beiden Vorzeichenspiegelungen (A -> -A und B -> -B, implementiert durch B bzw. A) ist nur "
      f"I/2 (Lösungsraum Dimension {rhos_inv[0]}), deren Modulfluss trivial ist; invariant nur unter A -> -A bleibt "
      f"rho = (I + bB)/2 (Dimension {rhos_inv[1]}) mit nichttrivialem Modulfluss (b = 0,4: |sigma_t(A) - A| = "
      f"{np.linalg.norm(sig_t - A):.2f})", ok_st)
befund("Gegenüber z71 Teil 3 (zwei binäre Beobachtungen, harmonisch genau bei p = 1/2, dort gesetzt) folgt p = 1/2 hier aus der "
       "Vorzeichensymmetrie des Gesetzes, sofern sie in derselben irreduziblen Darstellung implementiert ist; diese Prämisse "
       "ist nicht hergeleitet. Der Zustand darf die volle Vorzeichengruppe nicht teilen (sonst Spurzustand ohne Modulzeit, wie "
       "das gleichgewichtete Gemisch in z71); er braucht das Gewicht rho_m (z63)", True)

# ------------------------------------------------------------------ Teil 9
say("\nTeil 9: Flavourauslesung (Analyse 8 §6)")


def haar_u(nn):
    Zc = (rng.normal(size=(nn, nn)) + 1j*rng.normal(size=(nn, nn)))/np.sqrt(2)
    Qm, Rm = np.linalg.qr(Zc)
    return Qm*(np.diag(Rm)/np.abs(np.diag(Rm)))


ok_fl = True
for _ in range(20):
    U = haar_u(3)
    P = [np.diag(np.eye(3)[a]) for a in range(3)]
    Q = [np.outer(U[:, i], U[:, i].conj()) for i in range(3)]
    ok_fl &= all(abs(np.trace(P[a] @ Q[i]).real - abs(U[a, i])**2) < 1e-13 for a in range(3) for i in range(3))
    Js = []
    for a in range(3):
        for b in range(3):
            for i in range(3):
                for j in range(3):
                    if a != b and i != j:
                        Jt = np.trace(P[a] @ Q[i] @ P[b] @ Q[j]).imag
                        ok_fl &= abs(Jt - (U[a, i]*U[b, j]*np.conj(U[a, j])*np.conj(U[b, i])).imag) < 1e-13
                        Js.append(abs(Jt))
    ok_fl &= max(Js) - min(Js) < 1e-13
    D1, D2 = np.diag(np.exp(1j*rng.uniform(0, 2*np.pi, 3))), np.diag(np.exp(1j*rng.uniform(0, 2*np.pi, 3)))
    U2 = D1 @ U @ D2
    ok_fl &= abs(abs((U2[0, 1]*U2[1, 2]*np.conj(U2[0, 2])*np.conj(U2[1, 1])).imag) - Js[0]) < 1e-13
t12, t23, t13, dl = 0.227, 0.042, 0.0037, 1.2
s12, s23, s13, c12, c23, c13 = np.sin(t12), np.sin(t23), np.sin(t13), np.cos(t12), np.cos(t23), np.cos(t13)
R23 = np.array([[1, 0, 0], [0, c23, s23], [0, -s23, c23]])
U13 = np.array([[c13, 0, s13*np.exp(-1j*dl)], [0, 1, 0], [-s13*np.exp(1j*dl), 0, c13]])
R12 = np.array([[c12, s12, 0], [-s12, c12, 0], [0, 0, 1]])
Vp = R23 @ U13 @ R12
P = [np.diag(np.eye(3)[a]) for a in range(3)]
Q = [np.outer(Vp[:, i], Vp[:, i].conj()) for i in range(3)]
J_proj = np.trace(P[0] @ Q[1] @ P[1] @ Q[2]).imag
J_pdg = c12*s12*c23*s23*c13**2*s13*np.sin(dl)
ok_fl &= abs(J_proj - J_pdg) < 1e-15
check(f"für 20 Haar-zufällige U in U(3): |U_ai|^2 = Tr(P_a Q_i), Im Tr(P_a Q_i P_b Q_j) = Im(U_ai U_bj U_aj* U_bi*), |J| gleich "
      f"für alle a != b, i != j und invariant unter Umphasung; Standardparametrisierung: J = c12 s12 c23 s23 c13^2 s13 sin delta "
      f"= {J_pdg:.3e} wird exakt ausgelesen (Jarlskog)", ok_fl)
stand = (ROOT / "docs/research/TFPT-STAND.md").read_text()
abschnitt5 = stand.split("## 5. Was fehlt")[1].split("## 6.")[0]
rest_zeile = [z_.strip().lstrip("- ").replace("**", "").rstrip(".") for z_ in abschnitt5.splitlines() if "Restflavor" in z_]
ZAHLEN["tfpt_stand_5_restflavor"] = rest_zeile
befund(f"Standard (Jarlskog 1985): die Identitäten gelten für jedes unitäre U und legen nichts fest. Flavour, Massen, alpha, die 16 "
       f"und E8 sind gelöste Bausteine (TFPT-STAND §3) und werden nicht neu geöffnet. Beim Flavour führt TFPT-STAND §5 als offen "
       f"nur '{rest_zeile[0] if rest_zeile else '(Zeile nicht gefunden)'}'; keine dieser Größen folgt aus den "
       "Auslesungsidentitäten, am Stand ändert sich nichts", bool(rest_zeile))

# ------------------------------------------------------------------ Teil 10
say("\nTeil 10: Jacobson-Kette (Analyse 8 §5), nicht gerechnet")
befund("Die Kette 'Verschränkungsgleichgewicht => Einstein' (Jacobson 1995/2015) setzt das Gleichgewicht der Verschränkungs"
       "entropie kleiner Kugeln und die Flächennormierung (Entropiedichte 1/(4 G)) voraus; beides ist hier nicht hergeleitet "
       "(GRAV.ENTROPY.EQUILIBRIUM.01 bleibt bedingt, die Flächennormierung steht in TFPT-STAND §5). Firewall "
       "GRAV.NONCIRCULAR.01: kein Schluss 'modular => Bisognano-Wichmann'; BW gilt erst in einem lokalen Poincaré-kovarianten "
       "Netz mit Vakuum, nie umgekehrt aus einem Modulfluss", True)

urteil = ("Die siebte Analyse hat in ihrem Hauptpunkt recht: Die thermische Zuordnung in z77 war vertauscht. Mit "
          "u_m(r) = Summe (m-j) w(r+j)/((m+1) W_m(r)) für 'Teilchen hoch' und 1 - u für 'Loch runter' gilt "
          "u W_m(r) + (1 - u(r+1)) W_m(r+1) = m/(m+1) W_(m+1)(r) exakt für jedes w (symbolisch m = 2..10, numerisch bei jeder "
          "Temperatur), die Formel der fünften Analyse ist genau das Lochgewicht 1 - u, der Vakuumgrenzfall ist z76 und die "
          "beschränkte Lösung ist eindeutig. Das Negativ von z77/H6h37 (Residuum 0,94, 'hält nicht wie behauptet') war unsere "
          "Fehlzuordnung; korrigiert ist der thermische Zweipfadkanal bei jeder Temperatur exakt zustandserhaltend. Der Petz-"
          "Rückkanal auf den zentralen Paaren gibt korrigiert (0; 1,9794; 1,9794; 2,0398) statt (0; 1,9802; 1,9802; 3,9589), "
          "nahe l(l+1); die Werte der Analyse (0; 1,9864; 1,9864; 2,0281) sind nicht exakt reproduziert. Die Mischung "
          "verschobener Vakua (Analyse 8) gilt exakt für die Einteilchendaten, und kantenaufgelöst gibt Petz l(l+1) bei jeder "
          "Temperatur; der quasifreie Zustand ist aber nicht diese Mischung (Cov(n_x, n_y) = c(y)(1 - c(x)) gegen 0), und der "
          "Anteil 1/(m+1) ergibt sich nur als Zählidentität, nicht als Mehrpaaranteil von J_-1 ((m-2)/(3m-2)). Die relativen "
          "Zustände aus z78 haben den Überlapp 1/det(I - (i/2)(Y - X)) = 4/(r^2 - (t + 2i)^2), daraus s^2 = det(Y - X); das ist "
          "der Kern nu = 1 aus z74 an den Röhrenpunkten X + iI. Die Wurzel des Kerns ist nicht positiv (E[(det Q)^2] = -3/16), "
          "das ist die bekannte Wallach-Menge. Zusammengesetzter Spin 2 eines freien Skalars ist ein Kontinuum s^2/(240 pi) ohne "
          "Pol, wie z75. Die Spin-2-Ebene-Welle ist Standard-Spinor-Helizität, linear, nichts Neues gegenüber z65. Die "
          "Vorzeichensymmetrie A -> -A erzwingt in einer irreduziblen Darstellung AB + BA = 0, also Cl_2(C) = M_2(C), tau = i "
          "und p = 1/2; gegenüber z71 folgt p = 1/2 damit aus einer Symmetrie statt gesetzt zu sein, bedingt auf diese Prämisse. "
          "Die Flavourauslesung ist Standard und ändert TFPT-STAND §5 nicht; die Jacobson-Kette bleibt bedingt.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z80_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
