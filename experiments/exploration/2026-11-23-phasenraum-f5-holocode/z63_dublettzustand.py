"""Z63: Lücke (A) -- legt die originale Kovarianz einen treuen Zustand auf dem Dublett der Naht fest?

Leitlinie 'Rekonstruktion statt Auswahl': kein neues Auswahlprinzip, nur Antworten und Symmetrien des Originals.
Original: c(x) = 1/(1 + e^{4 pi x}) auf x in Z + alpha, alpha = 1/4, 3/4 (Hardy-Marken k = 4x ungerade), KMS bei
beta = 4 pi, rekonstruierter Einteilchengenerator h = 4 pi x (z58). Wörterbuch des Auftrags (H6, z48): Multiplikation mit
z auf der Naht ist der Übertragungsschritt S^4 = J_-1, also H0(O(1)) = span{1, z} <-> zwei benachbarte Moden {x, x+1}
eines Sektors (m = 4 Hardy-Marken je Einheitswindung). z62 braucht einen treuen Zustand auf End(V), V = H0(O(1))^*.

Teil 1  Kandidaten auf W_x = span{e_x, e_(x+1)}: (i) C_W/Tr und (1 - C_W)/Tr, (ii) N=1-Ecke des Originalzustands,
        (iv) Hardy-Produkt mit Gewichten aus dem KMS-Zustand, gerade Ecke {|00>, |11>}: Treue, Basispunkt- und
        Sektorunabhängigkeit, effektives beta.
Teil 2  Eindeutigkeit: Die KMS-Bedingung für h_W bei beta = 1 hat genau eine Lösung; die N=1-Ecke liegt im
        Zentralisator, ihr Modulfluss ist die Einschränkung des Originalflusses (der von C_W/Tr nicht).
Teil 3  Symmetrien: R = Delta^(-i/2) ist auf W skalar; Spiegel x -> -x allein kehrt den Zustand um, mit
        Teilchen-Loch-Tausch (c(-x) = 1 - c(x)) erhält er ihn; der Spurzustand ist symmetrisch, aber nur KMS bei beta = 0.
Teil 4  Weitere Dubletts und Wörterbuch: (iii) Fermipaar k = -1, +1 (Lücke 2 pi), Teilchen/Loch und Sektordublett
        (Spurzustand); (v) Fock-Ecke {Omega, Twistzustand} (Lücke pi, R = Vierteldrehung); (vi) Stromleiter {|0>, |4>}
        aus eta_q und Ein-Quant-Ecke benachbarter Strommoden (Lücke 4 pi). Allgemein: m Marken je Einheitswindung geben
        Lücke pi m; R hat auf P(V_m) Ordnung 4/ggT(m, 4); Sym^m rho_1 auf dem R-invarianten Dublett ist rho_m. No-go: auf
        jedem fermionischen Einteilchendublett hat R projektiv Ordnung <= 2.
Teil 5  Standardform des kanonischen Zustands: Blochlänge, Rapidität von sqrt(rho), det sqrt(rho), Ruhesystem,
        Modulfluss als Drehung der Naht, Achse = Polachse des Szegő-Kreises, Aberration, Schmidt-Verhältnis = q.
Teil 6  Bilanz und Urteil.
"""
import json
from itertools import product
from math import gcd
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla

HERE = Path(__file__).parent
mp.mp.dps = 40
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(63)


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


PI4 = 4*np.pi
ALPHAS = (0.25, 0.75)


def c(x):
    return 1/(1 + np.exp(PI4*x))


def cb(x):                                            # 1 - c(x) ohne Auslöschung
    return 1/(1 + np.exp(-PI4*x))


def gibbs(gap):
    w = np.array([1.0, np.exp(-gap)])
    return np.diag(w/w.sum())


def beta_eff(rho, gap=PI4):                           # in Einheiten des Originals (beta = 1 <-> 4 pi je x-Einheit)
    p = np.real(np.diag(rho))
    return float(np.log(p[0]/p[1])/gap)


def kms_matrix(H, beta=1.0):
    """Zeilen: Tr(rho (A e^{-beta H} B e^{beta H} - B A)) = 0 für alle Matrixeinheiten A, B (rho zeilenweise)."""
    n = H.shape[0]
    E, Ei = sla.expm(-beta*H), sla.expm(beta*H)
    rows = []
    for i, j, k, l in product(range(n), repeat=4):
        A = np.zeros((n, n))
        A[i, j] = 1
        B = np.zeros((n, n))
        B[k, l] = 1
        rows.append((A @ E @ B @ Ei - B @ A).T.ravel())
    return np.array(rows, dtype=complex)


def kms_loesungen(H, beta=1.0):
    L = kms_matrix(H, beta)
    _, s, vh = np.linalg.svd(L)
    rank = int(np.sum(s > 1e-10*s.max()))
    n = H.shape[0]
    return [v.conj().reshape(n, n) for v in vh[rank:]]


def proj_ordnung(U, nmax=8):
    for n in range(1, nmax + 1):
        P = np.linalg.matrix_power(U, n)
        if np.allclose(P, np.trace(P)/P.shape[0]*np.eye(P.shape[0]), atol=1e-10):
            return n
    return None


# Zweimoden-Fockraum (Jordan-Wigner), Basis |n1 n2>, |1> = besetzt
a1 = np.array([[0.0, 1.0], [0.0, 0.0]])
A1, A2 = np.kron(a1, np.eye(2)), np.kron(np.diag([1.0, -1.0]), a1)
VAC = np.eye(4)[0]
F = np.column_stack([A1.T @ VAC, A2.T @ VAC])          # Einteilchenbasis (e_1, e_2) im Fockraum = N=1-Sektor


def fock2(x1, x2):
    return np.kron(np.diag([cb(x1), c(x1)]), np.diag([cb(x2), c(x2)]))


def ecke(rhoF, basis=F):
    blk = basis.T @ rhoF @ basis
    return blk/np.trace(blk)


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Kandidaten auf W_x = span{e_x, e_(x+1)} (Schnitte 1, z von O(1) im Wörterbuch z = S^4)")
punkte = [(al, r + al) for al in ALPHAS for r in range(-5, 5)]
RHO4 = gibbs(PI4)
ok_cov, ok_dg = True, True
kand = {k: [] for k in ("(i) C_W/Tr", "(i') (1 - C_W)/Tr", "(ii) N=1-Ecke", "(vii) gerade Ecke {|00>, |11>}")}
for al, x in punkte:
    rF = fock2(x, x + 1)
    C2 = np.array([[np.trace(rF @ Ai.T @ Aj) for Aj in (A1, A2)] for Ai in (A1, A2)])
    ok_cov &= np.allclose(C2, np.diag([c(x), c(x + 1)]), rtol=1e-12, atol=1e-300)
    Am = crand(2, 2)
    dG = sum(Am[i, j]*(Ai.T @ Aj) for i, Ai in enumerate((A1, A2)) for j, Aj in enumerate((A1, A2)))
    ok_dg &= np.allclose(F.T @ dG @ F, Am)
    kand["(i) C_W/Tr"].append(np.diag([c(x), c(x + 1)])/(c(x) + c(x + 1)))
    kand["(i') (1 - C_W)/Tr"].append(np.diag([cb(x), cb(x + 1)])/(cb(x) + cb(x + 1)))
    kand["(ii) N=1-Ecke"].append(ecke(rF))
    ev = np.array([cb(x)*cb(x + 1), c(x)*c(x + 1)])
    kand["(vii) gerade Ecke {|00>, |11>}"].append(np.diag(ev/ev.sum()))
check("der Zweimoden-Fockzustand des Originals reproduziert <a_i^+ a_j> = diag(c(x), c(x+1)) (20 Paare, beide Sektoren); "
      "A -> dGamma(A) ist auf dem N=1-Sektor die Identität, also End(W) = N=1-Ecke (ohne Wahl)", ok_cov and ok_dg)

tab1 = {}
for k, lst in kand.items():
    b = np.array([beta_eff(r) for r in lst])
    treu = all(np.linalg.eigvalsh(r).min() > 0 for r in lst)
    tab1[k] = {"treu": bool(treu), "beta_eff_min": float(b.min()), "beta_eff_max": float(b.max()),
               "basispunktfrei": bool(b.max() - b.min() < 1e-9)}
    say(f"   {k}: treu {treu}, beta_eff in [{b.min():.6f}, {b.max():.6f}] über x in [-5+alpha, 4+alpha], beide Sektoren")
ZAHLEN["kandidaten_W"] = tab1
dev_ii = max(np.abs(r - RHO4).max()/RHO4[1, 1] for r in kand["(ii) N=1-Ecke"])
check("(ii) die N=1-Ecke ist für alle 20 Paare (beide Seiten der Fermikante, beide Sektoren) dieselbe treue Dichtematrix "
      f"diag(1, e^(-4 pi))/(1 + e^(-4 pi)) (rel. Abweichung {dev_ii:.1e}): basispunkt- und sektorfrei, beta_eff = 1",
      dev_ii < 1e-9 and tab1["(ii) N=1-Ecke"]["basispunktfrei"] and tab1["(ii) N=1-Ecke"]["treu"])
b_h = {al: beta_eff(np.diag([c(al), c(al + 1)])) for al in ALPHAS}
b_f = {al: beta_eff(np.diag([c(al - 1), c(al)])) for al in ALPHAS}
ZAHLEN["(i)_beta_eff_hardy_basispunkt"] = b_h
ZAHLEN["(i)_beta_eff_ueber_fermikante"] = b_f
check("(i) C_W/Tr hängt vom Basispunkt ab (beta_eff von ~0 tief im Fermisee bis ~1 weit oben) und am Hardy-Basispunkt "
      f"x = alpha vom Sektor: beta_eff = {b_h[0.25]:.6f} (alpha = 1/4) gegen {b_h[0.75]:.6f} (alpha = 3/4); über der "
      f"Fermikante {b_f[0.25]:.4f} bzw. {b_f[0.75]:.4f}: kein Zustand auf H0(O(1)) ohne Wahl",
      not tab1["(i) C_W/Tr"]["basispunktfrei"] and tab1["(i) C_W/Tr"]["beta_eff_min"] < 0.01
      and tab1["(i) C_W/Tr"]["beta_eff_max"] > 0.99 and abs(b_h[0.25] - b_h[0.75]) > 1e-4)
check("(i') und (vii) sind ebenso basispunktabhängig (gerade Ecke: beta_eff = 2x + 1)",
      not tab1["(i') (1 - C_W)/Tr"]["basispunktfrei"] and not tab1["(vii) gerade Ecke {|00>, |11>}"]["basispunktfrei"]
      and all(abs(beta_eff(r) - (2*x + 1)) < 1e-9 for r, (al, x) in zip(kand["(vii) gerade Ecke {|00>, |11>}"], punkte)))

say("   (iv) Hardy-Produkt <s, s'> = sum_x w(x) conj(s_x) s'_x auf {x, x+1}, Gewicht aus dem KMS-Zustand:")
gew = {
    "1 (reines Hardy-Produkt)": lambda x: 1.0 + 0*x,
    "c (Fermi)": c,
    "1 - c (Loch)": cb,
    "c(1 - c) (Fluktuation)": lambda x: c(x)*cb(x),
    "e^{-h} (Boltzmann = Modularoperator)": lambda x: np.exp(-PI4*x),
    "e^{+h} (gespiegelt)": lambda x: np.exp(PI4*x),
    "e^{-h/2} (halber Modularoperator)": lambda x: np.exp(-2*np.pi*x),
}
tab_gew = {}
for k, w in gew.items():
    b = np.array([np.log(w(x)/w(x + 1))/PI4 for al, x in punkte])
    frei = bool(b.max() - b.min() < 1e-9)
    tab_gew[k] = {"basispunktfrei": frei, "beta_eff": round(float(b.mean()), 12) if frei else
                  [float(b.min()), float(b.max())], "kms_beta1": bool(frei and abs(b.mean() - 1) < 1e-9)}
    say(f"      {k}: basispunktfrei {frei}, beta_eff {tab_gew[k]['beta_eff']}")
ZAHLEN["kandidaten_gewichte_iv"] = tab_gew
frei_set = sorted(k for k, v in tab_gew.items() if v["basispunktfrei"])
kms_set = [k for k, v in tab_gew.items() if v["kms_beta1"]]
check("(iv) basispunktfrei sind genau die Exponentialgewichte e^{-s h} (s = 0, 1, -1, 1/2), also KMS-Zustände bei beta = s; "
      "die Temperatur des Originals (beta = 1) lässt genau eines: e^{-h}, und das ist (ii)",
      frei_set == sorted(["1 (reines Hardy-Produkt)", "e^{-h} (Boltzmann = Modularoperator)", "e^{+h} (gespiegelt)",
                          "e^{-h/2} (halber Modularoperator)"])
      and kms_set == ["e^{-h} (Boltzmann = Modularoperator)"])

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Eindeutigkeit (KMS) und Verträglichkeit mit dem Modulfluss des Originals")
HW = PI4*np.diag([0.0, 1.0])
sol = kms_loesungen(HW)
rho_kms = sol[0]/np.trace(sol[0]) if len(sol) == 1 else np.zeros((2, 2))
same_L = np.allclose(kms_matrix(HW), kms_matrix(HW + PI4*1.75*np.eye(2)))
ok_rand = True
for n in (2, 3):
    for _ in range(5):
        X = crand(n, n)
        Hr = (X + X.conj().T)/2
        s_ = kms_loesungen(Hr)
        ok_rand &= len(s_) == 1 and np.allclose(s_[0]/np.trace(s_[0]), sla.expm(-Hr)/np.trace(sla.expm(-Hr)))
check("die KMS-Bedingung (beta = 1) für den eingeschränkten Generator h_W = 4 pi diag(x, x+1) hat auf End(W) genau eine "
      "Lösung, rho = diag(1, e^(-4 pi))/(1 + e^(-4 pi)); die Konstante 4 pi x fällt aus der Bedingung heraus "
      "(Basispunkt); ebenso eindeutig für 10 zufällige Generatoren (n = 2, 3)",
      len(sol) == 1 and np.allclose(rho_kms, RHO4) and same_L and ok_rand)
x0 = 0.25
rF = fock2(x0, x0 + 1)
P1 = F @ F.T
rc = ecke(rF)
ld, lc = np.log(np.diag(rF)), np.log(np.real(np.diag(rc)))
ok_fl = True
for t in (0.013, 0.1, 0.37):
    UF, Uc = np.diag(np.exp(1j*t*ld)), np.diag(np.exp(1j*t*lc))
    for _ in range(10):
        Am = crand(2, 2)
        ok_fl &= np.allclose(F.T @ UF @ (F @ Am @ F.T) @ UF.conj().T @ F, Uc @ Am @ Uc.conj().T)
check("die N=1-Projektion liegt im Zentralisator des Originalzustands ([rho_F, P_1] = 0), und der Modulfluss der Ecke ist "
      "die Einschränkung des Originalflusses (Pedersen-Takesaki), Drehrate 4 pi", np.allclose(rF @ P1, P1 @ rF) and ok_fl)
rate_h = {al: float(np.log(c(al)/c(al + 1))) for al in ALPHAS}
rate_f = {al: float(np.log(c(al - 1)/c(al))) for al in ALPHAS}
ZAHLEN["drehrate_C_W"] = {"hardy_basispunkt": rate_h, "ueber_fermikante": rate_f, "original": PI4}
check("C_W/Tr hat dagegen die Drehrate log(c_x/c_(x+1)): am Hardy-Basispunkt "
      f"{rate_h[0.25]:.4f} bzw. {rate_h[0.75]:.6f}, über der Fermikante {rate_f[0.25]:.4f} bzw. {rate_f[0.75]:.4f} statt "
      f"4 pi = {PI4:.6f}: sein Modulfluss ist nicht die Zeit des Originals",
      all(abs(v - PI4) > 1e-5 for v in list(rate_h.values()) + list(rate_f.values())))

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Symmetrien des Originals")
modes = [r + al for al in ALPHAS for r in range(-10, 10)]
ok_R = all(abs(np.exp(1j*PI4*x/2) - 1j**int(round(4*x))) < 1e-9 for x in modes)
ok_Rs = all(abs(1j**int(round(4*x)) - 1j**int(round(4*(x + 1)))) < 1e-12 for x in modes)
check("die Vierteldrehung R (Charakter i^k, k = 4x) ist der Modulfluss bei t = -1/2: e^{i h/2} = i^k auf allen Moden; auf "
      "W_x = {x, x+1} ist R skalar (i^k = i^(k+4)) und stellt keine Bedingung", ok_R and ok_Rs)
units = [np.outer(np.eye(2)[i], np.eye(2)[j]) for i in range(2) for j in range(2)]
dim_comm = 4 - np.linalg.matrix_rank(np.array([(HW @ E - E @ HW).ravel() for E in units]).T)
check("modulflussinvariante Zustände auf End(W) bilden nur die Familie diag(p, 1 - p) (Kommutante von h_W hat Dimension 2): "
      "der Fluss legt die Achse fest, die Länge legt erst KMS fest", dim_comm == 2)
SW = np.array([[0.0, 1.0], [1.0, 0.0]])
UPH = np.zeros((4, 4))
for n1, n2 in product((0, 1), repeat=2):
    UPH[2*(1 - n1) + (1 - n2), 2*n1 + n2] = 1            # Teilchen-Loch: Besetzungsmuster (n1, n2) -> (1-n1, 1-n2)
ph_N1 = F.T @ UPH @ F
ok_c = all(abs(c(-x) - cb(x)) < 1e-15 for x in modes)
ok_sp, ok_spph = True, True
for al, x in punkte:
    rW = ecke(fock2(x, x + 1))                           # Basis (e_x, e_(x+1))
    rWs = ecke(fock2(-x, -x - 1))                        # Spiegelpaar in Basis (M e_x, M e_(x+1)) = (e_-x, e_(-x-1))
    ok_sp &= np.allclose(rWs, SW @ rW @ SW) and not np.allclose(rWs, rW)
    ok_spph &= np.allclose(ph_N1 @ rW @ ph_N1.T, rWs)
check("Spiegel x -> -x (Sektor alpha -> 1 - alpha) allein bildet den Dublettzustand auf seine Blochumkehr ab (w -> 1/w, "
      "Hardy-Pol <-> Gegenpol: Pfeil der Modulzeit); zusammen mit dem Teilchen-Loch-Tausch (c(-x) = 1 - c(x), auf N=1 der "
      "Tausch) ist er eine Symmetrie und erhält den Zustand", ok_c and ok_sp and ok_spph
      and np.allclose(np.abs(ph_N1), SW))
spur = np.eye(2)/2
check("der Spurzustand ist unter allen diesen Abbildungen invariant und treu, aber KMS nur bei beta = 0 (kein Modulfluss "
      "auf End(V)); die Temperatur des Originals schließt ihn aus",
      np.allclose(SW @ spur @ SW, spur) and len(kms_loesungen(HW, 0.0)) == 1
      and np.allclose(kms_loesungen(HW, 0.0)[0]/np.trace(kms_loesungen(HW, 0.0)[0]), spur)
      and not np.allclose(rho_kms, spur))

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Weitere Dubletts und die Abhängigkeit vom Wörterbuch")
rP = ecke(fock2(-0.25, 0.25))                            # Basis (e_-1/4 [alpha = 3/4, k = -1], e_+1/4 [alpha = 1/4, k = +1])
RP = np.diag([1j**-1, 1j**1])
check("(iii) Fermipaar k = -1 (x = -1/4, alpha = 3/4) und k = +1 (x = +1/4, alpha = 1/4): N=1-Ecke = Gibbs mit Lücke "
      "2 pi; R wirkt darauf als Halbdrehung (projektive Ordnung 2); kein Bild von {1, z} unter z = S^4 (Abstand 1/2, zwei "
      "Sektoren), wohl aber unter z = S^2", np.allclose(rP, gibbs(2*np.pi)) and proj_ordnung(RP) == 2)
ptl = np.diag([c(0.25), cb(-0.25)])
ptl = ptl/np.trace(ptl)
th = [mp.nsum(lambda Q: mp.exp(-4*mp.pi*(Q**2/2 + (al - mp.mpf(1)/2)*Q)), [-mp.inf, mp.inf])
      for al in (mp.mpf(1)/4, mp.mpf(3)/4)]
check("(iii') Teilchen bei k = +1 und Loch bei k = -1 sind gleich besetzt (c(1/4) = 1 - c(-1/4)), das Sektordublett hat "
      "gleiche Thetagewichte (Q -> -Q): beides gibt den Spurzustand (KMS nur bei beta = 0)",
      np.allclose(ptl, spur, atol=1e-15) and abs(th[0] - th[1]) < mp.mpf(10)**-30)


def fock_sektor(alpha, rmin=-4, rmax=3):
    mo = [r + alpha for r in range(rmin, rmax + 1)]
    zst = []
    for occ in product((0, 1), repeat=len(mo)):
        dl0 = sum(x for x, n in zip(mo, occ) if x > 0 and n) + sum(-x for x, n in zip(mo, occ) if x < 0 and not n)
        Q = sum(n for x, n in zip(mo, occ) if x > 0) - sum(1 - n for x, n in zip(mo, occ) if x < 0)
        w = np.prod([c(x) if n else cb(x) for x, n in zip(mo, occ)])
        zst.append((int(round(4*dl0)), Q, w))
    w0 = [w for M, Q, w in zst if M == 0][0]
    return [(M, Q, w/w0) for M, Q, w in zst]


def p_part(n):                                          # Partitionszahl p(n) durch dynamische Programmierung
    if n < 0:
        return 0
    tab = [1] + [0]*n
    for part in range(1, n + 1):
        for s in range(part, n + 1):
            tab[s] += tab[s - part]
    return tab[n]


ok_boltz, ok_mult, mults = True, True, {}
ecken_v = {}
for al in ALPHAS:
    zst = fock_sektor(al)
    ok_boltz &= all(abs(w*np.exp(np.pi*M) - 1) < 1e-9 for M, Q, w in zst if M <= 16)
    mult = [sum(1 for M, Q, w in zst if M == j) for j in range(9)]
    mb = []
    for j in range(9):
        s_ = 0
        for Q in range(-6, 7):
            Nq = j/4 - (Q*Q/2 + (al - 0.5)*Q)
            if abs(Nq - round(Nq)) < 1e-9:
                s_ += p_part(int(round(Nq)))
        mb.append(s_)
    ok_mult &= mult == mb
    mults[al] = mult
    w1 = [w for M, Q, w in zst if M == 1]
    w4 = [w for M, Q, w in zst if M == 4]
    ecken_v[al] = (len(w1), w1[0], len(w4), w4[0], [Q for M, Q, w in zst if M == 1][0])
ZAHLEN["fock_markenvielfachheit_M0_bis_8"] = {str(k): v for k, v in mults.items()}
check("Fockzustand jedes Sektors (256 Besetzungsmuster): Gewicht relativ zu Omega ist genau e^{-pi M} mit der Marke "
      "M = 4 Delta L0 (KMS des vollen Modulgenerators); Vielfachheiten der Marken 0..8 = Bosonisierung (Jacobi) "
      f"{mults[0.25]}: Marke 2 ist in einem Sektor leer", ok_boltz and ok_mult and mults[0.25][2] == 0)
check("(v) Fock-Ecke {Omega, Twistzustand} (Marken 0 und 1, je einfach; alpha = 1/4: Teilchen bei x = 1/4, alpha = 3/4: "
      "Loch bei x = -1/4): Gewichtsverhältnis e^(-pi) in beiden Sektoren, also Gibbs mit Lücke pi; {Omega, J_-1 Omega} "
      "(Marken 0, 4) gibt e^(-4 pi)",
      all(v[0] == 1 and abs(v[1]*np.exp(np.pi) - 1) < 1e-9 and v[2] == 1 and abs(v[3]*np.exp(PI4) - 1) < 1e-9
          for v in ecken_v.values()))
R1 = 40


def strom(m, al=mp.mpf(1)/4):                          # omega(J_-m J_m) = sum c(x + m)(1 - c(x)), 1 - c ohne Auslöschung
    return mp.fsum(1/(1 + mp.exp(4*mp.pi*(r + al + m)))/(1 + mp.exp(-4*mp.pi*(r + al))) for r in range(-R1, R1))


om = {m: strom(m) for m in range(1, 6)}
q2 = om[1]/(1 + om[1])                                  # z48: x = omega/(1 + omega) = q^2
leiter = np.diag([1.0, float(q2)])/(1 + float(q2))
bm = {m: (om[m]/m)/(1 + om[m]/m) for m in om}
check("(vi) Stromleiter: q^2 = omega/(1 + omega) aus der gemessenen Stromnorm (z48) gibt auf {|0>, |4>} = {1, z} die Ecke "
      "diag(1, q^2)/(1 + q^2) = rho_4; die Ein-Quant-Ecke benachbarter Strommoden {m, m+1} (Besetzungen aus den "
      "Fermisummen, m = 1..4) hat ebenfalls das Verhältnis e^(-4 pi)",
      np.allclose(leiter, RHO4, rtol=1e-12) and abs(q2*mp.exp(4*mp.pi) - 1) < mp.mpf(10)**-25
      and all(abs(bm[m + 1]/bm[m]*mp.exp(4*mp.pi) - 1) < mp.mpf(10)**-25 for m in range(1, 5)))

say("   Allgemein: m Hardy-Marken je Einheitswindung (z = S^m) geben auf V_m = H0(O(1)) die Lücke pi m:")
ordn, bahnen = {}, {}
for m in (1, 2, 4):
    U = np.diag([1, 1j**m])
    ordn[m] = proj_ordnung(U)
    pts = sorted({(round(float(np.real(1j**(m*n))), 9), round(float(np.imag(1j**(m*n))), 9)) for n in range(8)})
    bahnen[m] = pts
    say(f"      m = {m}: Lücke {m} pi, R auf P(V_m) projektive Ordnung {ordn[m]}, R-Bahn von [1:1]: {len(pts)} Punkt(e)")
ZAHLEN["R_ordnung_auf_P(V_m)"] = ordn
check("R wirkt auf P(V_m) mit Ordnung 4/ggT(m, 4) = 4, 2, 1 für m = 1, 2, 4; vier getrennte Marken {1, i, -1, -i} als "
      "R-Bahn gibt es nur für m ungerade (Twistornaht); im Wörterbuch z = S^4 ist R auf der Naht trivial und die vier "
      "Marken (zeta^4 = 1) fallen auf den einen Punkt w = 1",
      [ordn[m] for m in (1, 2, 4)] == [4 // gcd(m, 4) for m in (1, 2, 4)] == [4, 2, 1] and len(bahnen[1]) == 4
      and len(bahnen[4]) == 1)


def sym_zustand(rho, m):
    big = rho
    for _ in range(m - 1):
        big = np.kron(big, rho)
    B = np.zeros((2**m, m + 1))
    for idx in range(2**m):
        B[idx, bin(idx).count("1")] = 1
    B = B/np.linalg.norm(B, axis=0)
    S = B.T @ big @ B
    return S/np.trace(S)


RHO1, RHO2 = gibbs(np.pi), gibbs(2*np.pi)
S4, S2 = sym_zustand(RHO1, 4), sym_zustand(RHO1, 2)
e04 = np.diag([S4[0, 0], S4[4, 4]])/(S4[0, 0] + S4[4, 4])
e02 = np.diag([S2[0, 0], S2[2, 2]])/(S2[0, 0] + S2[2, 2])
Rsym = np.array([1j**j for j in range(5)])
grad0 = [j for j in range(5) if abs(Rsym[j] - 1) < 1e-12]
psi_mark = np.array([1.0, -1.0])/np.sqrt(2)              # 1 - zeta^4 im R-invarianten Dublett {zeta^0, zeta^4}
ew_mark = float(psi_mark @ RHO4 @ psi_mark)
check("Sym^4 rho_1 auf dem Träger H0(O(4)) eingeschränkt auf das R-invariante Dublett {1, zeta^4} = H0(O_w(1)) ist rho_4 "
      "(ebenso Sym^2 rho_1 auf {1, zeta^2} = rho_2): die Zustände der verschiedenen Nähte sind verträglich; der "
      f"H6-Markenzustand 1 - zeta^4 liegt in diesem Dublett auf dem Äquator der Zustandsachse (Erwartung {ew_mark:.3f})",
      np.allclose(e04, RHO4) and np.allclose(e02, RHO2) and grad0 == [0, 4] and abs(ew_mark - 0.5) < 1e-12)
marks = np.array([1, 1j, -1, -1j])
stab = [k/16 for k in range(32) if np.allclose(np.sort_complex(np.round(np.exp(-1j*np.pi*k/16)*marks, 9)),
                                               np.sort_complex(np.round(marks, 9)))]
check("die vier Marken brechen die kontinuierliche Modulrotation der Naht (zeta -> e^{-i pi t} zeta) genau auf <R>: "
      f"invariant nur für t in {stab} (Periode 2)", stab == [0.0, 0.5, 1.0, 1.5])


def sektorladung(besetzt_rel):
    """besetzt_rel: dict Mode -> +1 (Teilchen über dem See) / -1 (Loch im See); Rückgabe (Q_1/4, Q_3/4, Parität)."""
    q = {0.25: 0, 0.75: 0}
    for x, s in besetzt_rel.items():
        q[0.25 if abs((x - 0.25) % 1) < 1e-9 else 0.75] += s
    return (q[0.25], q[0.75], sum(abs(s) for s in besetzt_rel.values()) % 2)


real = {
    "m = 4, {e_x, e_x+1}, x = 1/4": ({0.25: 1}, {1.25: 1}),
    "m = 4, {Omega, J_-1 Omega}": ({}, {0.25: 1, -0.75: -1}),
    "m = 2, Fermipaar k = -1, +1": ({}, {0.25: 1, -0.25: -1}),
    "m = 1, {Omega, Twistzustand}": ({}, {0.25: 1}),
}
neutral = {}
for k, (s0, s1) in real.items():
    q0, q1 = sektorladung(s0), sektorladung(s1)
    neutral[k] = (q0[0] == q1[0] and q0[1] == q1[1] and q0[2] == q1[2])
ZAHLEN["End(V)_in_neutraler_Antwortalgebra"] = neutral
check("Ladungsbuchhaltung der Matrixeinheiten |b><a| von End(V): neutral je Quartersektor und gerade (also im Bereich der "
      "Stromantworten) nur für m = 4; m = 2 verschiebt Ladung zwischen den Sektoren, m = 1 ist ungerade und geladen",
      [neutral[k] for k in real] == [True, True, False, False])
ks = [k for k in range(-41, 42) if k % 2]
Rd = np.array([1j**k for k in ks])
Ei = [i for i, k in enumerate(ks) if k % 4 == 1]
Emi = [i for i, k in enumerate(ks) if k % 4 == 3]
ok_ng = np.allclose(Rd**2, -1)
ords = set()
for _ in range(200):
    vecs = []
    for _ in range(2):
        v = np.zeros(len(ks), complex)
        idx = Ei if rng.random() < 0.5 else Emi
        v[idx] = crand(len(idx))
        vecs.append(v)
    Qm, _ = np.linalg.qr(np.column_stack(vecs))
    RW = Qm.conj().T @ (Rd[:, None]*Qm)
    ok_ng &= np.allclose(RW @ RW, -np.eye(2))
    ords.add(proj_ordnung(RW))
check("No-go: auf dem fermionischen Einteilchenraum (ungerade Marken) ist R^2 = -1; jedes R-invariante Dublett (200 "
      f"zufällige) hat projektive R-Ordnung in {sorted(ords)}: die Vierteldrehung der vier Marken ist auf keinem Dublett "
      "aus zwei Originalmoden R", ok_ng and ords <= {1, 2})

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Standardform des kanonischen Dublettzustands (sqrt(rho) in Herm(2) = R^(1,3))")


def standardform(gap):
    rho = gibbs(gap)
    p0, p1 = rho[0, 0], rho[1, 1]
    sq = np.diag(np.sqrt([p0, p1]))
    t_, z_ = (np.sqrt(p0) + np.sqrt(p1))/2, (np.sqrt(p0) - np.sqrt(p1))/2
    det_sq = float(np.sqrt(p0*p1))
    g = np.diag([(p1/p0)**0.125, (p0/p1)**0.125])
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    tt = 0.01
    U = np.diag([p0**(1j*tt), p1**(1j*tt)])
    winkel = float(np.angle((U @ X @ U.conj().T)[0, 1]))/tt
    tP = 2*np.pi/gap
    UP = np.diag([p0**(1j*tP), p1**(1j*tP)])
    wb = np.exp(1j*np.linspace(0, 2*np.pi, 7)[:-1])
    zb = (1 - np.abs(np.exp(gap/4)*wb)**2)/(1 + np.abs(np.exp(gap/4)*wb)**2)
    S = -float(p0*np.log(p0) + p1*np.log(p1))
    return {
        "rho": [float(p0), float(p1)], "bloch": float(p0 - p1), "reinheit": float(p0**2 + p1**2),
        "sqrt_rho_t_z": [float(t_), float(z_)], "rapiditaet": float(np.arctanh(z_/t_)), "gamma": float(np.cosh(np.arctanh(z_/t_))),
        "det_sqrt_rho": det_sq, "eigenlaenge": float(np.sqrt(det_sq)), "entropie": S,
        "ruhe_ok": bool(abs(np.linalg.det(g) - 1) < 1e-12 and np.allclose(g @ sq @ g.T, np.sqrt(det_sq)*np.eye(2))),
        "cosh_rel": float(np.trace(sq)/(2*np.sqrt(det_sq))), "modulwinkel_je_einheit": winkel, "periode": tP,
        "periode_ok": bool(np.allclose(UP @ X @ UP.conj().T, X)), "aberration_z": [float(v) for v in zb],
        "aberration_rad": float(np.arccos(-zb[0])), "schmidt": float(np.sqrt(p1/p0)),
    }


SF = {m: standardform(np.pi*m) for m in (4, 1, 2)}
ok_sf = True
for m, d in SF.items():
    gap = np.pi*m
    ok_sf &= abs(d["bloch"] - np.tanh(gap/2)) < 1e-12 and abs(d["rapiditaet"] - gap/4) < 1e-9
    ok_sf &= abs(d["det_sqrt_rho"] - 1/(2*np.cosh(gap/2))) < 1e-14 and d["ruhe_ok"]
    ok_sf &= abs(d["cosh_rel"] - np.cosh(gap/4)) < 1e-9 and abs(d["modulwinkel_je_einheit"] - gap) < 1e-9
    ok_sf &= d["periode_ok"] and np.allclose(d["aberration_z"], -np.tanh(gap/4))
    ok_sf &= abs(d["aberration_rad"] - 2*np.arctan(np.exp(-gap/4))) < 1e-12 and abs(d["schmidt"] - np.exp(-gap/2)) < 1e-12
ZAHLEN["standardform"] = {f"m={m}": d for m, d in SF.items()}
d4 = SF[4]
say(f"   m = 4 (Wörterbuch z = S^4): rho = diag({d4['rho'][0]:.10f}, {d4['rho'][1]:.6e}), Blochlänge "
    f"{d4['bloch']:.10f} = tanh 2 pi, Reinheit {d4['reinheit']:.10f}")
say(f"      sqrt(rho) = t + z sigma_3 mit (t, z) = ({d4['sqrt_rho_t_z'][0]:.8f}, {d4['sqrt_rho_t_z'][1]:.8f}); Rapidität "
    f"{d4['rapiditaet']:.10f} = pi, gamma = cosh pi = {d4['gamma']:.6f}")
say(f"      det sqrt(rho) = {d4['det_sqrt_rho']:.6e} = 1/(2 cosh 2 pi), Minkowski-Länge {d4['eigenlaenge']:.6f}, "
    f"Entropie {d4['entropie']:.3e}")
say(f"      Modulfluss: Drehung um die Achse mit {d4['modulwinkel_je_einheit']:.6f} = 4 pi je Moduleinheit, Periode "
    f"{d4['periode']:.4f}; Aberration: die Naht liegt im Ruhesystem {np.degrees(d4['aberration_rad']):.3f} Grad vom Gegenpol")
for m in (1, 2):
    d = SF[m]
    say(f"   m = {m}: Blochlänge {d['bloch']:.6f}, Rapidität {d['rapiditaet']:.6f} = {m} pi/4, gamma {d['gamma']:.4f}, "
        f"det sqrt(rho) {d['det_sqrt_rho']:.6f}, Modulwinkel {d['modulwinkel_je_einheit']:.6f} je Einheit, Periode "
        f"{d['periode']:.3f}")
check("Standardform für m = 4, 1, 2: Blochlänge tanh(pi m/2), Rapidität von sqrt(rho) gegen das Hardy-System pi m/4, "
      "det sqrt(rho) = 1/(2 cosh(pi m/2)), Boost g in SL(2,C) ins Ruhesystem, Modulfluss = Drehung um die Achse mit "
      "pi m je Moduleinheit (Periode 2/m), Aberration 2 arctan e^{-pi m/4}, Schmidt-Verhältnis e^{-pi m/2}", ok_sf)
check("für m = 4 ist das Schmidt-Verhältnis von sqrt(rho) genau q = e^(-2 pi) aus z46, und die Purifikation "
      "(|1>|1> + q|z>|z>)/sqrt(1 + q^2) ist der auf zwei Glieder gekürzte Übertragungs-Nullvektor eta_q",
      abs(d4["schmidt"] - float(mp.sqrt(q2))) < 1e-15
      and np.allclose(np.array([1.0, float(mp.sqrt(q2))])/np.sqrt(1 + float(q2)), np.sqrt(d4["rho"])))
ok_moeb = True
for s in (0.3, 1.0, 2.5):
    t = -s/PI4
    U = np.diag([RHO4[0, 0]**(1j*t), RHO4[1, 1]**(1j*t)])
    for w in (0.5 + 0.2j, 2.0 - 1j, np.exp(0.7j)):
        v = U @ np.array([1.0, w])
        ok_moeb &= abs(v[1]/v[0] - np.exp(1j*s)*w) < 1e-10
th_ = np.linspace(0, 2*np.pi, 64, endpoint=False)
ok_sz = True
for _ in range(5):
    a_, b_ = crand(2)
    f = a_ + b_*np.exp(1j*th_)
    ok_sz &= abs(np.mean(f) - a_) < 1e-12 and abs(np.mean(np.exp(-1j*th_)*f) - b_) < 1e-12
check("der Modulfluss bewegt die Punkte der Naht P(V) als Drehung w -> e^{is} w (s = x-Zeit, t = -s/(4 pi)): die Zeit des "
      "Originals ist die Drehung der Naht um ihre Polachse; der dominante reine Zustand ist der Schnitt 1 = <1, .> = "
      "Auswertung bei w = 0 (Szegő-Reproduktion), der Szegő-Kreis |w| = 1 ist der Äquator", ok_moeb and ok_sz)
befund("Der kanonische Dublettzustand ist fast rein (1 - Blochlänge = 7e-6): sein Beobachter bewegt sich gegen das "
       "Hardy-System mit Rapidität pi (gamma = 11,6) zum Hardy-Pol w = 0; in seinem Ruhesystem liegt der Szegő-Kreis "
       "innerhalb von 5 Grad um den Gegenpol", True)

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: Bilanz")
tafel = [
    ("(i) C_W/Tr, (1 - C_W)/Tr", "treu, aber basispunkt- und sektorabhängig, nicht KMS: nicht kanonisch"),
    ("(ii) N=1-Ecke {x, x+1}", "treu, basispunkt- und sektorfrei, einzige KMS-Lösung: rho_4 (kanonisch für m = 4)"),
    ("(iii) Fermipaar k = -1, +1", "treu, KMS mit Lücke 2 pi (rho_2), gehört zum Wörterbuch z = S^2, Sektorladung"),
    ("(iii') Teilchen/Loch, Sektordublett", "Spurzustand: treu, symmetrisch, KMS nur bei beta = 0"),
    ("(iv) Hardy-Produkt mit Gewicht", "basispunktfrei nur e^{-s h}; KMS bei beta = 1 nur e^{-h} = (ii)"),
    ("(v) Fock-Ecke {Omega, Twistzustand}", "treu, KMS mit Lücke pi (rho_1), Twistornaht m = 1, ungerade und geladen"),
    ("(vi) Stromleiter {|0>, |4>}, Strommoden", "rho_4, gleich (ii), Schmidt-Verhältnis q = e^(-2 pi)"),
    ("(vii) gerade Ecke {|00>, |11>}", "basispunktabhängig (beta_eff = 2x + 1)"),
]
for k_, v_ in tafel:
    say(f"   {k_}: {v_}")
ZAHLEN["bilanz"] = dict(tafel)
ZAHLEN["alternativen_gezaehlt"] = {"kandidaten_auf_W": 4, "gewichte_iv": len(gew), "weitere_dubletts": 5,
                                   "woerterbuecher_m": [1, 2, 4], "summe_zustandskonstruktionen": 4 + len(gew) + 5}
ZAHLEN["status"] = {
    "E": ["für jede Naht mit m Marken je Einheitswindung genau ein treuer KMS-Zustand auf End(H0(O(1))), "
          "rho_m = diag(1, e^{-pi m})/(1 + e^{-pi m}), basispunkt- und sektorfrei, gleich der Ecke des Originalzustands",
          "m = 4: rho_4, Blochlänge tanh 2 pi, Rapidität pi, det sqrt(rho) = 1/(2 cosh 2 pi), Modulfluss 4 pi je Einheit",
          "R = Delta^(-i/2); auf P(V_m) Ordnung 4/ggT(m, 4); No-go für Vierteldrehung auf Originaldubletts",
          "rho_4 = Sym^4 rho_1 auf dem R-invarianten Dublett"],
    "C": ["dass die 3+1-Naht die Strom-/Übertragungsnaht z = S^4 ist (Wörterbuch des Auftrags; gestützt dadurch, dass nur "
          "m = 4 im neutralen Bereich der Stromantworten liegt)"],
    "L": ["Beobachter des Originals mit Rapidität pi gegen das Hardy-System; Zeit = Drehung der Naht um die Polachse"],
    "O": ["welche Naht die Himmelssphäre ist: z = S^4 (rho_4, R trivial auf P(V)) oder die Twistornaht mit vier getrennten "
          "Marken z = S (rho_1, ungerade Feldalgebra); hängt mit (B) zusammen"],
}
befund("Ohne Wörterbuch bleibt genau eine diskrete Freiheit: die Zahl m der Hardy-Marken je Einheitswindung. Sie steckt "
       "nicht im Zustand, sondern in der Zuordnung der Naht; ein freier kontinuierlicher Parameter bleibt nicht", True)

urteil = ("Lücke (A) ist bedingt gelöst. Exakt: Für jede Naht, deren Einheitswindung m Hardy-Marken umfasst, legt das "
          "Original auf End(H0(O(1))) genau einen treuen Zustand fest, der KMS für den eingeschränkten rekonstruierten "
          "Generator ist: rho_m = diag(1, e^{-pi m})/(1 + e^{-pi m}) in der Hardy-Basis (1, z). Er ist basispunkt- und "
          "sektorfrei und gleich der Ecke des Originalzustands; R = Delta^(-i/2) und Spiegel mit Teilchen-Loch-Tausch "
          "erhalten ihn, der Spiegel allein kehrt ihn um. Die Einteilchenkovarianz selbst (C/Tr C) ist dagegen basispunkt- "
          "und sektorabhängig und nicht KMS. Im Wörterbuch des Auftrags (z = S^4, m = 4) ist rho = diag(1, e^{-4 pi})/"
          "(1 + e^{-4 pi}): Blochlänge tanh 2 pi, sqrt(rho) hat gegen das Hardy-System die Rapidität pi (gamma = cosh pi), "
          "det sqrt(rho) = 1/(2 cosh 2 pi), der Modulfluss ist die Drehung der Naht um die Polachse (4 pi je "
          "Moduleinheit, 1 je Einheit x-Zeit), und das Schmidt-Verhältnis ist q = e^{-2 pi}. Bedingt ist das Wörterbuch: "
          "Die Twistornaht mit vier getrennten Marken verlangt m ungerade (R ist dort die Vierteldrehung von P(V)); das "
          "gibt rho_1 mit Lücke pi auf {Omega, Twistzustand}, außerhalb der neutralen Stromantworten. Beide Zustände "
          "sind verträglich: rho_4 ist Sym^4 rho_1 auf dem R-invarianten Dublett. Welche Naht die Himmelssphäre ist, "
          "bleibt offen und hängt mit (B) zusammen.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z63_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
