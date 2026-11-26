"""Z174: Schritt 1 der Kernstrategie aus z173. Was gibt die Standardform des Nahtzustands für Statistik und Komposition?

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Bausteine bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, N numerisch, C bedingt, L Lesart, O offen.
Bestand: z113 (Modultest der 2 pi-Drehung: pro Naht nur U = +1 oder G = -1; die Z4-Verklebung wurde dort nicht geprüft),
z117 (Verklebung Z_SU (x) R = +1 auf dem Spinormodul, ein Bit), z122 (Z_SU = i^(2h-2), Sektorcharakter X = R i^(-2h)),
z125 (Gatterfamilie U = B_theta^+ e^(-i phi H_Q) B_theta, P_univ), z130 (Standardform einer Naht mit Tomita-Partner, P_univ
nativ in der Partnerlesart), z144 (P_uhr: Phase als Funktion der Uhrwerte genau für phi in (pi/2) Z), z173 (Abstand zur
Komplettlösung: 2 Klauseln, 1 freie Zahl, 7 Prämissen, dazu K19 und K20).

Teil A  Standardform Naht + Tomita-Partner im endlichen Modell (K = 2 und 3 Moden je Seite, KMS bei beta = 4 pi):
        Standardvektor, Transfer-SU(2), Gatterfamilie auf dem Standardvektor, Klein-Identität, Uhrbahnen der beiden nativen
        Standardvektoren (bosonisch, fermionisch), Wechselwirkung, Symmetrien des natürlichen Kegels.
Teil B  Modultest der Z4-Verklebung nach der Methode von z113 (nichtchirale Naht = Original (x) Partner, N Kopien):
        Bahnen, Anomalie, Positivität für X = R i^(-2h), den Spiegel X' = R i^(2h) und die J-kovariante Fassung
        X_p = R R-quer^-1 i^(-2h); robuste Auswertung der langen Nähte über Normalformen.
Teil C  Eintrag ins Suchprogramm (OriginLab reconstruction_search) und neuer Abstand zur Komplettlösung.
"""
import json
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[2] / "research" / "originlab"))
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T_START = time.time()


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def zeit():
    return f"(Laufzeit bisher {time.time() - T_START:.0f} s)"


# =================================================================== Teil A
say("Teil A: Standardform einer Naht mit ihrem Tomita-Partner (endliches Modell)")


def jw(n):
    I_, Zp, sm = np.eye(2), np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [0.0, 0.0]])
    out = []
    for j in range(n):
        m_ = [Zp]*j + [sm] + [I_]*(n - j - 1)
        op = m_[0]
        for x in m_[1:]:
            op = np.kron(op, x)
        out.append(op)
    return out


def nrm(X):
    return float(np.abs(X).max()) if np.size(X) else 0.0


def modell(K, kfun):
    """K Moden je Seite, Jordan-Wigner-Ordnung A_0..A_(K-1), B_0..B_(K-1); Zustand rho = e^(-kfun(Q))/Z auf der Naht."""
    C = jw(2*K)
    d = 2**K
    NUM = [np.diag(x.T @ x).round().astype(int) for x in C]
    QA, QB = sum(NUM[:K]), sum(NUM[K:])
    JP = sum(C[c].T @ C[K + c] for c in range(K))
    JM = JP.T
    qa = np.array([bin(i).count("1") for i in range(d)])
    p = np.exp(-kfun(qa))
    p = p/p.sum()
    Om = np.diag(np.sqrt(p)).reshape(-1)                   # vec(rho^(1/2)) in Zeilenordnung: |i>_A (x) |j>_B
    F = np.zeros((d*d, d*d))                               # Flip der Seiten: vec(Y) -> vec(Y^T)
    for i in range(d):
        for j in range(d):
            F[j*d + i, i*d + j] = 1
    return {"C": C, "K": K, "d": d, "QA": QA, "QB": QB, "N": QA + QB, "JP": JP, "JM": JM, "JX": (JP + JM)/2,
            "JY": (JP - JM)/(2j), "JZ": np.diag((QA - QB)/2.0), "p": p, "Om": Om, "F": F, "NUM": NUM,
            "HQ": np.diag((QA**2 + QB**2)/2.0), "HM": np.diag((QA**2 - QB**2)/2.0),
            "RA": np.diag(1j**(QA % 4)), "pair": np.diag(1j**((QA - QB) % 4))}


KMS = lambda q: np.pi*q                                    # K = 4 pi * Q/4: Quarterenergie 1/4 bei beta = 4 pi
KMS2 = lambda q: 2*np.pi*q**2 - np.pi*q                    # K = 4 pi E_lad(1/4, Q) (z113, z130)
MOD = {K: modell(K, KMS) for K in (2, 3)}
MOD2 = {K: modell(K, KMS2) for K in (2, 3)}

# --- A1 Standardform aus der Definition
a1 = {}
for K, m in MOD.items():
    d = m["d"]
    p = m["p"]
    Dl = np.array([p[i]/p[j] for i in range(d) for j in range(d)])             # Delta vec(Y) = vec(rho Y rho^-1)
    rest = 0.0
    for a in range(d):
        for b in range(d):
            x = np.zeros((d, d))
            x[a, b] = 1
            xO = np.kron(x, np.eye(d)) @ m["Om"]
            xdO = np.kron(x.T, np.eye(d)) @ m["Om"]
            rest = max(rest, nrm(m["F"] @ (np.sqrt(Dl)*xO).conj() - xdO))      # J Delta^(1/2) x Omega = x^+ Omega, J = F K
    traeger = nrm(m["Om"][m["QA"] != m["QB"]])
    moden = max(nrm(m["Om"][m["NUM"][c] != m["NUM"][K + c]]) for c in range(K))
    paaruhr = nrm(np.diag(Dl**(-0.5j)) - m["pair"])                             # Delta^(it) bei t = -1/2
    a1[K] = (rest, traeger, moden, paaruhr, nrm(m["F"] @ m["Om"].conj() - m["Om"]))
ZAHLEN["A1_reste (S, Träger Q_A != Q_B, Moden, Paaruhr, J Omega)"] = {str(K): v for K, v in a1.items()}
check("A1 [E] Standardform aus der Definition (K = 2: 16, K = 3: 64 Zustände): Omega = vec(rho^(1/2)), J = Seitentausch mit "
      "Konjugation, Delta = rho (x) rho^-1; S = J Delta^(1/2) erfüllt S x Omega = x^+ Omega auf allen Matrixeinheiten. Omega "
      "trägt nur n_A = n_B je Mode (also Q_A = Q_B), J Omega = Omega, und Delta^(-i/2) ist die Paaruhr i^(Q_A - Q_B) (z130)",
      all(max(v) < 1e-12 for v in a1.values()))

# --- A2 Transfer-Singulett
rng = np.random.default_rng(174)
a2 = {}
for K, m in MOD.items():
    r = max(nrm(m["JP"] @ m["Om"]), nrm(m["JM"] @ m["Om"]),
            max(nrm(expm(-1j*th*m["JX"]) @ m["Om"] - m["Om"]) for th in np.linspace(0, 2*np.pi, 25)))
    d = m["d"]
    X = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
    Uz, _ = np.linalg.qr(X)                                                     # Kontrolle: Zustand nicht diagonal in der Besetzung
    rho = Uz @ np.diag(m["p"]) @ Uz.conj().T
    w, V = np.linalg.eigh(rho)
    Omz = (V @ np.diag(np.sqrt(w)) @ V.conj().T).reshape(-1)
    a2[K] = (r, nrm(m["JP"] @ Omz), nrm(MOD2[K]["JP"] @ MOD2[K]["Om"]))
ZAHLEN["A2 (Rest, Kontrolle gedrehter Zustand, zweiter Zustand)"] = {str(K): v for K, v in a2.items()}
check("A2 [E] Der Standardvektor ist ein Transfer-Singulett: J^+ Omega = J^- Omega = 0, und der Strahlteiler "
      "B_theta = e^(-i theta J_x) lässt Omega für jeden Winkel fest (K = 2, 3; ebenso für K = 4 pi E_lad). Der Winkel theta "
      "der Komposition ist auf dem Standardvektor unsichtbar. Kontrolle: Ein in der Besetzungsbasis nicht diagonaler "
      f"Zustand ist kein Singulett (|J^+ Omega| = {a2[2][1]:.3f}, {a2[3][1]:.3f})",
      all(v[0] < 1e-12 and v[2] < 1e-12 and v[1] > 1e-2 for v in a2.values()))


def gatter(m, th, phi, gesetz="HQ"):
    B = expm(-1j*th*m["JX"])
    return B.conj().T @ expm(-1j*phi*m[gesetz]) @ B


# --- A3 Gatterfamilie auf dem Standardvektor
a3 = {}
for K, m in MOD.items():
    rp = rm = 0.0
    for th in np.linspace(0, np.pi, 13):
        for phi in np.linspace(0, 2*np.pi, 9):
            rp = max(rp, nrm(gatter(m, th, phi) @ m["Om"] - expm(-1j*phi*m["HQ"]) @ m["Om"]))
            rm = max(rm, nrm(gatter(m, th, phi, "HM") @ m["Om"] - m["Om"]))
    a3[K] = (rp, rm)
ZAHLEN["A3_reste (+ Gesetz, - Gesetz)"] = {str(K): v for K, v in a3.items()}
check("A3 [E] Gatterfamilie U = B_theta^+ e^(-i phi H) B_theta auf Omega (13 Winkel, 9 Phasen): Mit dem nahtsymmetrischen "
      "Gesetz H_Q = (Q_A^2 + Q_B^2)/2 gilt U Omega = e^(-i phi H_Q) Omega unabhängig von theta, auf dem Paarsektor also die "
      "Phase e^(-i phi Q^2). Mit dem modularen Gesetz (Q_A^2 - Q_B^2)/2 bleibt Omega für alle theta, phi fest",
      all(max(v) < 1e-12 for v in a3.values()))

# --- A4 Klein-Identität und fermionischer Standardvektor
klein = all(abs(np.exp(-1j*np.pi*Q*Q/2) - (1 + 1j*(-1)**Q)/(1 + 1j)) < 1e-12
            and abs(np.exp(-1j*np.pi*Q*Q/2) - 1j**((-Q) % 4)*(-1)**(Q*(Q - 1)//2)) < 1e-12 for Q in range(-8, 9))
a4 = {}
OMG = {}
for K, m in MOD.items():
    q = np.exp(-np.pi/2)
    Og = np.zeros(m["d"]**2)
    Og[0] = 1
    for c in range(K):
        Og = Og + q*(m["C"][c].T @ (m["C"][K + c].T @ Og))                      # prod_c (1 + q a_c^+ b_c^+) |0>
    Og = Og/np.linalg.norm(Og)
    OMG[K] = Og
    sgn = np.diag((-1.0)**(m["QB"]*(m["QB"] - 1)//2))
    a4[K] = (nrm(Og - sgn @ m["Om"]),
             nrm(gatter(m, 0.7, np.pi/2) @ m["Om"] - np.linalg.inv(m["RA"]) @ Og),
             nrm(gatter(m, 0.7, -np.pi/2) @ m["Om"] - m["RA"] @ Og),
             nrm(gatter(m, 0.7, np.pi) @ m["Om"] - m["RA"] @ m["RA"] @ m["Om"]))
ZAHLEN["A4_reste"] = {str(K): v for K, v in a4.items()}
check("A4 [E] Klein-Identität: e^(-i pi Q^2/2) = (1 + i (-1)^Q)/(1 + i) = i^(-Q) (-1)^(Q(Q-1)/2) für alle ganzen Q (|Q| <= 8): "
      "Die Viertelphase des Gatters ist auf dem Paarsektor genau die Klein-Verdrehung Z der gradierten Standardform. Der "
      "fermionische Standardvektor Omega_g = prod_c (1 + q a_c^+ b_c^+)|0> ist (-1)^(Q(Q-1)/2) Omega (z122, z130). Daher "
      "U(theta, +pi/2) Omega = R^-1 Omega_g, U(theta, -pi/2) Omega = R Omega_g und U(theta, pi) Omega = R^2 Omega",
      klein and all(max(v) < 1e-12 for v in a4.values()))

# --- A5 Uhrbahnen der beiden nativen Standardvektoren
def in_bahn(v, w, RA):
    return [k for k in range(4) if abs(abs(np.vdot(np.linalg.matrix_power(RA, k) @ w, v)) - 1) < 1e-10]


a5 = {}
for name, mods in (("K = pi Q", MOD), ("K = 4 pi E_lad", MOD2)):
    for K, m in mods.items():
        Og = np.diag((-1.0)**(m["QB"]*(m["QB"] - 1)//2)) @ m["Om"]
        bos, fer = {}, {}
        for kk in range(24):
            v = expm(-1j*kk*np.pi/12*m["HQ"]) @ m["Om"]
            if in_bahn(v, m["Om"], m["RA"]):
                bos[kk] = in_bahn(v, m["Om"], m["RA"])
            if in_bahn(v, Og, m["RA"]):
                fer[kk] = in_bahn(v, Og, m["RA"])
        a5[f"{name}, K = {K}"] = {"bosonisch": bos, "fermionisch": fer}
ZAHLEN["A5_treffer (phi = k pi/12 -> Uhrpotenz)"] = {k: {a: {str(x): y for x, y in b.items()} for a, b in v.items()}
                                                     for k, v in a5.items()}
check("A5 [E, entscheidend] Uhrbahn der Standardvektoren: Unter den 24 Phasen phi = k pi/12 liegt U(theta, phi) Omega genau "
      "für phi = 0 und pi in der Uhrbahn {R^j Omega} des bosonischen Standardvektors und genau für phi = pi/2 und 3 pi/2 in "
      "der Uhrbahn {R^j Omega_g} des fermionischen; sonst in keiner (K = 2, 3; beide KMS-Zustände). Das ist die Menge "
      "phi in (pi/2) Z aus z144 (P_uhr), hier als Eigenschaft des Zustandsvektors",
      all(v == {"bosonisch": {0: [0], 12: [2]}, "fermionisch": {6: [3], 18: [1]}} for v in a5.values()))


def bahnphasen(k, l, qmax=6, mmax=3):
    """exakt: Phase -phi Q^2 - 2 chi m (phi = k/24, chi = l/24 Umdrehungen) auf dem Paarsektor mit Oszillatorniveau m je Seite."""
    ph = {(Q, mm): (Fr(-k*Q*Q, 24) - Fr(2*l*mm, 24)) % 1 for Q in range(qmax + 1) for mm in range(mmax + 1)}
    bos = [j for j in range(4) if all(v == Fr(j*Q, 4) % 1 for (Q, mm), v in ph.items())]
    fer = [j for j in range(4) if all(v == (Fr(j*Q, 4) + Fr(Q*(Q - 1), 4)) % 1 for (Q, mm), v in ph.items())]
    return bos, fer


gitter = {(k, l): bahnphasen(k, l) for k in range(24) for l in range(24)}
treffer = sorted((k, l) for (k, l), (b, f) in gitter.items() if b or f)
ZAHLEN["A6_gitter_treffer (phi = k pi/12, chi = l pi/12)"] = [list(x) for x in treffer]
check("A6 [E] Mit Oszillatoranteil (Phasengesetz e^(-i phi H_Q - i chi L_0^osc) je Seite, exakte Phasenarithmetik, Q <= 6, "
      "Niveau m <= 3, Gitter 24 x 24): Das Bild liegt genau dann in einer der beiden Uhrbahnen, wenn phi in (pi/2) Z und "
      "chi in pi Z. Der Rest chi = pi ist (-1)^(L_0^osc + L-quer_0^osc), auf niveaugleichen Zuständen die 2 pi-Drehung "
      "U_2pi aus z113, die Omega festlässt. Bis auf U_2pi ist das Kriterium also P_uhr (z144: phi in (pi/2) Z, chi = 0)",
      treffer == sorted((k, l) for k in (0, 6, 12, 18) for l in (0, 12)))


def linearitaet(m, U):
    """Normalisiert U die Transfer-su(2) in jedem Sektor fester Gesamtladung? Rest der besten linearen Darstellung."""
    worst = 0.0
    for N in sorted(set(m["N"].tolist())):
        idx = np.where(m["N"] == N)[0]
        Js = [J[np.ix_(idx, idx)] for J in (m["JX"], m["JY"], m["JZ"])] + [np.eye(len(idx))]
        A = np.array([J.reshape(-1) for J in Js]).T
        Us = U[np.ix_(idx, idx)]
        for J in Js[:3]:
            y = (Us @ J @ Us.conj().T).reshape(-1)
            x = np.linalg.lstsq(A, y, rcond=None)[0]
            worst = max(worst, float(np.linalg.norm(A @ x - y)))
    return worst


a7 = {}
for K, m in MOD.items():
    a7[K] = {kk: (round(linearitaet(m, gatter(m, np.pi/2, kk*np.pi/4)), 6),
                  round(linearitaet(m, gatter(m, np.pi/2, kk*np.pi/4, "HM")), 6)) for kk in range(8)}
ZAHLEN["A7_linearitaetsrest (phi = k pi/4: + Gesetz, - Gesetz)"] = {str(K): {str(k): v for k, v in t.items()}
                                                                 for K, t in a7.items()}
check("A7 [E] Wechselwirkung (z144: ein Gatter ist einkörperig, wenn es die Transfer-su(2) je Gesamtladung normalisiert): "
      "Im symmetrischen Gesetz ist das Gatter bei theta = pi/2 genau für phi in {0, pi} einkörperig, im modularen Gesetz "
      "für jede Phase. Jedes Mitglied beider Familien, das Omega festlässt, ist damit wechselwirkungsfrei; die "
      "wechselwirkenden Gatter bewegen Omega, und unter ihnen landen nur phi = +-pi/2 auf einem nativen Standardvektor",
      all(all((t[kk][0] < 1e-9) == (kk in (0, 4)) and t[kk][1] < 1e-9 for kk in range(8)) for t in a7.values()))


def kegel(m, U, n=200):
    """kleinster Eigenwert und Hermitizitätsrest der Bilder positiver Rang-1-Matrizen (natürlicher Kegel = PSD-Matrizen)."""
    d = m["d"]
    r_ = np.random.default_rng(1740 + d)
    worst = herm = 0.0
    for _ in range(n):
        v = r_.normal(size=d) + 1j*r_.normal(size=d)
        W = (U @ np.outer(v, v.conj()).reshape(-1)).reshape(d, d)
        herm = max(herm, nrm(W - W.conj().T))
        worst = min(worst, float(np.linalg.eigvalsh((W + W.conj().T)/2).min()))
    return worst, herm


a8 = {}
for K, m in MOD.items():
    Bpi = expm(-1j*np.pi*m["JX"])
    kozyklus = np.diag(1j**(m["N"] % 4)*(-1.0)**(m["QA"]*m["QB"]))
    x = rng.normal(size=(m["d"], m["d"]))
    tausch = nrm(m["F"] @ np.kron(x, np.eye(m["d"])) @ m["F"] - np.kron(np.eye(m["d"]), x))
    kand = {"Seitentausch F": m["F"], "Modulfluss": np.diag(np.exp(-0.37j*np.pi*(m["QA"] - m["QB"]))), "Paaruhr": m["pair"],
            "Gatter (pi/2, pi/2)": gatter(m, np.pi/2, np.pi/2), "Gatter (pi/2, pi)": gatter(m, np.pi/2, np.pi),
            "Klein Z_A": np.diag(np.exp(-0.5j*np.pi*m["QA"]**2))}
    a8[K] = {"F = i^N (-1)^(Q_A Q_B) B_pi": nrm(m["F"] - kozyklus @ Bpi), "F M F = M'": tausch,
             "kegel": {k: kegel(m, U) for k, U in kand.items()},
             "fest": {k: nrm(U @ m["Om"] - m["Om"]) for k, U in kand.items()}}
ZAHLEN["A8_kegel (min Eigenwert, Hermitizitätsrest)"] = {str(K): {"F_rest": v["F = i^N (-1)^(Q_A Q_B) B_pi"],
                                                                "kegel": v["kegel"], "Omega fest": v["fest"]}
                                                       for K, v in a8.items()}
erh = ("Seitentausch F", "Modulfluss", "Paaruhr")
check("A8 [E] Symmetrien der Standardform (unitär, erhalten den natürlichen Kegel und vertauschen mit J): Seitentausch F "
      "(F M F = M'), Modulfluss und Paaruhr erhalten den Kegel (200 Stichproben je Modell) und lassen Omega fest. Exakt ist "
      "F = i^(Q_A + Q_B) (-1)^(Q_A Q_B) B_pi: der halbe Transferumlauf mal Produktuhr mal Klein-Kozyklus (z106). Das "
      "Viertelgatter, der einkörperige Tausch bei phi = pi und die Klein-Verdrehung erhalten den Kegel nicht: Die einzige "
      "nicht faktorisierte Symmetrie der Standardform ist der Seitentausch, und er ist einkörperig",
      all(v["F = i^N (-1)^(Q_A Q_B) B_pi"] < 1e-12 and v["F M F = M'"] < 1e-12
          and all(v["kegel"][k][0] > -1e-10 and v["kegel"][k][1] < 1e-10 and v["fest"][k] < 1e-12 for k in erh)
          and all(v["kegel"][k][0] < -1e-2 or v["kegel"][k][1] > 1e-2 for k in v["kegel"] if k not in erh)
          for v in a8.values()))
befund("Teil A: In der Standardform (Naht + Tomita-Partner) ist P_univ nativ (z130) und P_uhr gleichwertig dazu, dass die "
       "Komposition den Standardvektor in die Uhrbahn eines nativen Standardvektors abbildet. Beides zusammen ist eine "
       "Prämisse P_std (Komposition als Morphismus der gradierten Standardform), nicht zwei. Hergeleitet aus A* ist P_std "
       "nicht; P_ord (Kette aus mehr als zwei Nähten) berührt das Paar nicht", True)
say(f"   {zeit()}")

# =================================================================== Teil B
say("\nTeil B: Modultest der Z4-Verklebung (Methode z113, robuste Auswertung über Normalformen)")
mp.mp.dps = 50
I = mp.mpc(0, 1)
PI = mp.pi
TP = [mp.mpc('0.137', '1.093'), mp.mpc('-0.211', '0.947'), mp.mpc('0.389', '1.211')]
NULL = mp.mpf(10)**-30
GEN = {"S": (0, -1, 1, 0), "T": (1, 1, 0, 1)}
EINS = (1, 0, 0, 1)
QU, HA = Fr(1, 4), Fr(1, 2)


def R_(x):
    return mp.mpf(x.numerator)/x.denominator if isinstance(x, Fr) else mp.mpf(x)


def theta(a, b, tau):
    """theta[a;b](tau) = sum_n exp(i pi tau (n+a)^2 + 2 pi i (n+a) b)."""
    a, b = R_(a), R_(b)
    nmax = int(mp.sqrt((mp.mp.dps*2.31 + 12)/(PI*tau.imag))) + 3
    return mp.fsum(mp.exp(I*PI*tau*(n + a)**2 + 2*I*PI*(n + a)*b) for n in range(-nmax, nmax + 1))


def eta(tau):
    kmax = int(mp.sqrt((mp.mp.dps*2.31 + 12)/(3*PI*tau.imag))) + 3
    return mp.exp(2*I*PI*tau/24)*mp.fsum((-1)**k*mp.exp(I*PI*tau*k*(3*k - 1)) for k in range(-kmax, kmax + 1))


def act(g, tau):
    return (g[0]*tau + g[1])/(g[2]*tau + g[3])


def mul(g, h):
    return (g[0]*h[0] + g[1]*h[2], g[0]*h[1] + g[1]*h[3], g[2]*h[0] + g[3]*h[2], g[2]*h[1] + g[3]*h[3])


def wurzel(c, nmax=4608):
    if abs(abs(c) - 1) > 1e-12:
        return None
    r = Fr(float(mp.arg(c)/(2*PI))).limit_denominator(nmax)
    return r % 1 if abs(mp.exp(2*I*PI*R_(r)) - c) < 1e-12 else None


def werte(F, g):
    return [F(act(g, t)) for t in TP]


def proportional(v, w):
    if any(abs(x) < NULL for x in v + w):
        return mp.mpc(1) if all(abs(x) < NULL for x in v + w) else None
    c = v[0]/w[0]
    return c if all(abs(x/y - c) < mp.mpf(10)**-20*abs(c) for x, y in zip(v[1:], w[1:])) else None


def bahn(F, maxn=600):
    """SL(2,Z)-Bahn {F o g}: Vertreter, Werte, Permutationen und Phasen (Einheitswurzeln) von S und T (wie z113)."""
    reps, vals, eltern = [EINS], [werte(F, EINS)], {0: None}
    perm, ph = {"S": {}, "T": {}}, {"S": {}, "T": {}}
    k = 0
    while k < len(reps):
        for gn, gen in GEN.items():
            g2 = mul(reps[k], gen)
            v = werte(F, g2)
            for j, w in enumerate(vals):
                c = proportional(v, w)
                if c is not None:
                    perm[gn][k], ph[gn][k] = j, wurzel(c)
                    break
            else:
                reps.append(g2)
                vals.append(v)
                eltern[len(reps) - 1] = (k, gn)
                perm[gn][k], ph[gn][k] = len(reps) - 1, Fr(0)
                assert len(reps) <= maxn
        k += 1
    assert all(x is not None for d in ph.values() for x in d.values())
    return {"F": F, "reps": reps, "vals": vals, "perm": perm, "ph": ph, "n": len(reps), "eltern": eltern}


def lokalisiere(B, H):
    v = werte(H, EINS)
    for j, w in enumerate(B["vals"]):
        c = proportional(v, w)
        if c is not None:
            return j, wurzel(c)
    return None


def konsistenz(B, chi):
    """Invariante Kombination sum_j e^(2 pi i V_j) F_j mit Charakter chi; bad = Widersprüche (Anomalie)."""
    V = [None]*B["n"]
    V[0] = Fr(0)
    stack = [0]
    while stack:
        i = stack.pop()
        for gn in ("S", "T"):
            j, r = B["perm"][gn][i], B["ph"][gn][i]
            if V[j] is None:
                V[j] = (V[i] + r - chi[gn]) % 1
                stack.append(j)
    bad = {(V[i] + B["ph"][gn][i] - chi[gn] - V[B["perm"][gn][i]]) % 1 for i in range(B["n"]) for gn in ("S", "T")}
    return V, bad - {Fr(0)}


def charaktere(B):
    out = []
    for k in range(24):
        for j in range(4):
            chi = {"T": Fr(k, 24), "S": Fr(j, 4)}
            V, bad = konsistenz(B, chi)
            if not bad:
                out.append((chi, V))
    return out


def chistr(ch):
    return [f"T {c['T']}, S {c['S']}" for c, _ in ch]


def L_lad(u, t, tau):
    """sum_Q q^((Q+u)^2/2) e^(2 pi i t Q): Spur von e^(2 pi i t Q) über die Ladungen eines Quarters (z113)."""
    return mp.exp(-2*I*PI*R_(t)*R_(u))*theta(u, t, tau)


def osz(r, tau):
    """Tr e^(2 pi i r L_0^osc) q^(L_0^osc - 1/24) = e^(2 pi i r/24)/eta(tau + r) für eine chirale Bose-Leiter."""
    return mp.exp(2*I*PI*R_(r)/24)/eta(tau + R_(r))


def hol(par):
    ta, tb, ra, rb = par
    return lambda tau: L_lad(-QU, ta, tau)*L_lad(QU, tb, tau)*osz(ra, tau)*osz(rb, tau)


def naht(name, N=1):
    fl, fr = hol(EINSCHUB[name][0]), hol(EINSCHUB[name][1])
    return lambda tau: (fl(tau)*mp.conj(fr(tau)))**N


# Einschübe (links: Original; rechts: Partner, komplex konjugiert): Ladung e^(2 pi i (ta Q_a + tb Q_b)), Oszillatoren
# e^(2 pi i r L_0^osc) je Leiter. R = i^(Q_a - Q_b), i^(-2h) = i^(-(L_0^osc - L-quer_0^osc)) (z113: (-1)^(2h) = U_2pi).
Q4 = QU
EINSCHUB = {
    "1": ((0, 0, 0, 0), (0, 0, 0, 0)),
    "X = R i^(-2h)": ((Q4, -Q4, -Q4, -Q4), (0, 0, -Q4, -Q4)),
    "X' = R i^(2h)": ((Q4, -Q4, Q4, Q4), (0, 0, Q4, Q4)),
    "G = U R^2": ((HA, HA, HA, HA), (0, 0, HA, HA)),
    "X_p = R Rq^-1 i^(-2h)": ((Q4, -Q4, -Q4, -Q4), (Q4, -Q4, -Q4, -Q4)),
    "X_p^3": ((-Q4, Q4, Q4, Q4), (-Q4, Q4, Q4, Q4)),
    "i^(-2h)": ((0, 0, -Q4, -Q4), (0, 0, -Q4, -Q4)),
    "G_tot = U (-1)^(F+Fq)": ((HA, HA, HA, HA), (HA, HA, HA, HA)),
}


def xgcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x, y = xgcd(b, a % b)
    return g, y, x - (a//b)*y


def normalform(par):
    """Thetas [a;b] und Etas eta((A tau + B)/D) eines holomorphen Faktors."""
    ta, tb, ra, rb = par
    et = tuple((Fr(r).denominator, Fr(r).numerator % Fr(r).denominator, Fr(r).denominator) for r in (ra, rb))
    return (((-QU) % 1, Fr(ta) % 1), (QU % 1, Fr(tb) % 1)), et


def schritt(st, gen):
    """Normalform von N o gen bis auf eine Konstante: theta[a;b] -> theta[a; a+b+1/2] (T), theta[b; -a] (S);
    eta(M tau) mit M = [[A, B], [0, D]] -> Hermite-Normalform von M gen."""
    th, et = st
    if gen == "T":
        return tuple((a, (a + b + HA) % 1) for a, b in th), tuple((A, (B + A) % D, D) for A, B, D in et)
    et2 = []
    for A, B, D in et:
        g, x, _ = xgcd(B, D)
        D2 = A*D//g
        et2.append((g, (-x*A) % D2, D2))
    return tuple((b % 1, (-a) % 1) for a, b in th), tuple(et2)


def nwert(st, tau):
    th, et = st
    v = mp.mpc(1)
    for a, b in th:
        v *= theta(a, b, tau)
    for A, B, D in et:
        v /= eta((A*tau + B)/D)
    return v


def tabelle(name):
    """Bahn der Spur mit Einschub, Normalformen N_j und Konstanten c_j mit F o g_j = c_j N_j (an zwei Punkten geprüft)."""
    F = naht(name)
    B = bahn(F)
    forms = {0: (normalform(EINSCHUB[name][0]), normalform(EINSCHUB[name][1]))}
    for j in range(1, B["n"]):
        k, gn = B["eltern"][j]
        forms[j] = (schritt(forms[k][0], gn), schritt(forms[k][1], gn))
    consts, err = [], mp.mpf(0)
    for j in range(B["n"]):
        c, d = B["reps"][j][2], B["reps"][j][3]
        pts = TP[:2] if c == 0 else [mp.mpc(-mp.mpf(d)/c, mp.mpf(1)/abs(c)),
                                     mp.mpc(-mp.mpf(d)/c + mp.mpf('0.3')/abs(c), mp.mpf(1)/abs(c))]
        cs = [F(act(B["reps"][j], t))/(nwert(forms[j][0], t)*mp.conj(nwert(forms[j][1], t))) for t in pts]
        err = max(err, abs(cs[0] - cs[1])/abs(cs[0]))
        consts.append(cs[0])
    B["forms"], B["consts"] = forms, consts
    return B, float(err)


# --- B1 Bausteine und Normalformen
bau = []
for r in (-QU, QU, HA):
    q = mp.exp(2*I*PI*(TP[0] + R_(r)))
    prod = mp.mpc(1)
    for n in range(1, 400):
        prod *= 1 - q**n
    bau.append(abs(osz(r, TP[0]) - mp.exp(-2*I*PI*TP[0]/24)/prod))
TAB, nf_err = {}, {}
for name in EINSCHUB:
    TAB[name], nf_err[name] = tabelle(name)
bahnen = {name: B["n"] for name, B in TAB.items()}
ZAHLEN["B1"] = {"Oszillatorspur_Fehler": float(max(bau)), "Normalform_Konstanzfehler": nf_err, "Bahnen": bahnen,
                "max |c| der Vertreter": max(abs(g[2]) for B in TAB.values() for g in B["reps"])}
check("B1 [E] Bausteine (mpmath, 50 Stellen): Tr e^(2 pi i r L_0^osc) q^(L_0^osc - 1/24) = e^(2 pi i r/24)/eta(tau + r) gegen das "
      "Produkt für r = +-1/4, 1/2. Bahnen der Spuren mit Einschub: ohne Einschub 6, G = U R^2 und G_tot je 24, die "
      f"Z4-Einschübe X, X', X_p, X_p^3 und i^(-2h) je 96 Elemente ({bahnen}). Robuste Auswertung: Jedes Bahnelement ist "
      "c_j mal eine Normalform aus Thetafunktionen und eta((A tau + B)/D) (lange Nähte, A D = 16 bzw. 4); die Konstanz von "
      f"c_j ist an zwei Punkten je Element geprüft (max. Fehler {max(nf_err.values()):.0e})",
      max(bau) < mp.mpf(10)**-40 and max(nf_err.values()) < 1e-40
      and bahnen == {"1": 6, "X = R i^(-2h)": 96, "X' = R i^(2h)": 96, "G = U R^2": 24, "X_p = R Rq^-1 i^(-2h)": 96,
                     "X_p^3": 96, "i^(-2h)": 96, "G_tot = U (-1)^(F+Fq)": 24})
say(f"   {zeit()}")

# --- B2 Anomalie
NS = (1, 2, 3, 4)
BN = {(name, N): (TAB[name] if N == 1 else bahn(naht(name, N))) for name in EINSCHUB for N in NS
      if name not in ("X_p^3", "i^(-2h)")}
CH = {key: charaktere(B) for key, B in BN.items()}
anom = {name: {N: chistr(CH[(name, N)]) for N in NS} for name in ("X = R i^(-2h)", "X' = R i^(2h)", "X_p = R Rq^-1 i^(-2h)")}
enth = {}
for other in ("X_p^3", "i^(-2h)"):
    j, r = lokalisiere(TAB["X_p = R Rq^-1 i^(-2h)"], naht(other))
    V = [V for c, V in CH[("X_p = R Rq^-1 i^(-2h)", 1)] if c["T"] == 0 and c["S"] == 0][0]
    enth[other] = str((V[j] - r) % 1)
ZAHLEN["B2_charaktere"] = {k: {str(N): v for N, v in t.items()} for k, t in anom.items()}
ZAHLEN["B2_Xp_Bahn_enthaelt (rel. Koeffizient)"] = enth
zwei = ["T 0, S 0", "T 1/2, S 1/2"]
check("B2 [E, entscheidend] Anomalie: Die einseitige Verklebung X = R i^(-2h) (Sektorcharakter aus z122) und ihr Spiegel "
      "X' = R i^(2h) haben für N = 1 und N = 3 keinen widerspruchsfreien Charakter: Über ein einzelnes Paar Naht + Partner "
      "lässt sich die Verklebung nicht modulinvariant summieren. Für N = 2 und 4 ist die Bahnsumme frei (Charaktere "
      f"{zwei}). Die J-kovariante Fassung X_p = R R-quer^-1 i^(-2h) ist für jedes N frei; ihre Bahn enthält X_p^3 und "
      f"i^(-2h) mit Koeffizient {enth} (auf niveaugleichen Zuständen ist X_p = i^(-2h))",
      all(anom[x][1] == [] and anom[x][3] == [] and anom[x][2] == zwei and anom[x][4] == zwei
          for x in ("X = R i^(-2h)", "X' = R i^(2h)"))
      and all(anom["X_p = R Rq^-1 i^(-2h)"][N] == zwei for N in NS) and enth == {"X_p^3": "0", "i^(-2h)": "0"})
say(f"   {zeit()}")

# --- B3 Positivität
YS = (mp.mpf('0.35'), mp.mpf('0.6'), mp.mpf('1.0'))
KF = 24
GRID = [(Y, k) for Y in YS for k in range(KF)]


def gitterwerte(B):
    mp.mp.dps = 30
    out = [[B["consts"][j]*nwert(B["forms"][j][0], mp.mpc(mp.mpf(k)/KF, Y))*mp.conj(nwert(B["forms"][j][1], mp.mpc(mp.mpf(k)/KF, Y)))
            for (Y, k) in GRID] for j in range(B["n"])]
    mp.mp.dps = 50
    return out


GW = {name: gitterwerte(TAB[name]) for name in EINSCHUB if name not in ("X_p^3", "i^(-2h)")}
say(f"   Gitterwerte berechnet {zeit()}")


def moden(name, N):
    """Spinmoden |s| <= 5 der Bahnsumme von F^N mit trivialem Charakter (Koeffizient 1 für die Spur selbst)."""
    B1, B = TAB[name], BN[(name, N)]
    ch = [V for c, V in CH[(name, N)] if c["T"] == 0 and c["S"] == 0]
    if not ch:
        return None
    F1 = naht(name)
    tot = [mp.mpc(0)]*len(GRID)
    for k in range(B["n"]):
        g = B["reps"][k]
        j, r = (k, Fr(0)) if N == 1 else lokalisiere(B1, lambda tau, g=g: F1(act(g, tau)))
        ph = mp.exp(2*I*PI*(R_(ch[0][k]) + N*R_(r)))
        tot = [t + ph*x**N for t, x in zip(tot, GW[name][j])]
    return {Y: {s: mp.fsum(x*mp.exp(-2*I*PI*s*k/KF) for k, x in enumerate(tot[iy*KF:(iy + 1)*KF]))/KF
                for s in range(-5, 6)} for iy, Y in enumerate(YS)}


def comb(ms, cs):
    return {Y: {s: sum(c*m_[Y][s] for m_, c in zip(ms, cs)) for s in ms[0][Y]} for Y in YS}


def konj(m_):
    """Moden des inversen Einschubs: Z_(X^-1)(tau) = conj(Z_X(-conj tau)), also m_s -> conj(m_s)."""
    return {Y: {s: mp.conj(m_[Y][s]) for s in m_[Y]} for Y in YS}


def minimum(m_, m0):
    """kleinster normierter Spinmode (negativ: keine Hilbertraumspur) und größter Imaginärteil relativ zum größten Betrag."""
    gross = max(abs(m_[Y][s]) for Y in YS for s in m_[Y])
    return (float(min(m_[Y][s].real/abs(m0[Y][0]) for Y in YS for s in m_[Y])),
            float(max(abs(m_[Y][s].imag) for Y in YS for s in m_[Y])/gross))


pos = {}
for N in NS:
    M = {name: moden(name, N) for name in GW}
    m0 = M["1"]
    zeile = {"Z_1": minimum(m0, m0), "Z_1 + Z_G": minimum(comb([m0, M["G = U R^2"]], [1, 1]), m0),
             "Z_1 - Z_G": minimum(comb([m0, M["G = U R^2"]], [1, -1]), m0)}
    for xn in ("X = R i^(-2h)", "X' = R i^(2h)"):
        if M[xn]:
            for j in range(4):                                                  # Projektion auf den Eigenwert i^j
                tot = comb([m0, M["G = U R^2"], M[xn], konj(M[xn])], [1, (-1)**j, (-I)**j, I**j])
                zeile[f"{xn[:2].strip()} = i^{j}"] = minimum(tot, m0)
    for b in (1, -1):
        tot = comb([m0, M["G_tot = U (-1)^(F+Fq)"], M["X_p = R Rq^-1 i^(-2h)"]], [1, 1, b])
        zeile[f"X_p = {b:+d}"] = minimum(tot, m0)
    pos[N] = zeile
    say(f"      N = {N}: " + "; ".join(f"{k}: {v[0]:.2e}" for k, v in zeile.items()))
    say(f"   {zeit()}")
ZAHLEN["B3_positivitaet (min normierter Spinmode, rel. Imaginärteil; Y = 0.35, 0.6, 1.0, |s| <= 5)"] = \
    {str(N): z for N, z in pos.items()}
EPS = -1e-8
ok = lambda x: x[0] > EPS
reell = max(v[1] for z in pos.values() for v in z.values())
z113 = json.loads((HERE / "z113_ergebnis.json").read_text())["zahlen"]
ref = z113["positivitaet (min normierter Spinmode, Y = 0.35, 0.6, 1.0, |s| <= 5)"]
kontrolle = max(abs(pos[N]["Z_1 + Z_G"][0] - ref[str(N)]["G = U R^2"]["beta = +1"])
                / max(1.0, abs(ref[str(N)]["G = U R^2"]["beta = +1"])) for N in NS)
ZAHLEN["B3_kontrolle_gegen_z113"] = kontrolle
check("B3 [N] Kontrolle gegen z113: Z_1 + Z_G (Projektion G = +1) reproduziert für N = 1..4 die dortigen kleinsten "
      f"Spinmoden (N = 1: {pos[1]['Z_1 + Z_G'][0]:.2f}, relative Abweichung {kontrolle:.0e}); alle Summen sind reell "
      f"(größter relativer Imaginärteil {reell:.0e})", kontrolle < 1e-6 and reell < 1e-12)
eins = {N: sorted(k for k, v in pos[N].items() if k.startswith("X =") and ok(v)) for N in (2, 4)}
spiegel = {N: sorted(k for k, v in pos[N].items() if k.startswith("X' =") and ok(v)) for N in (2, 4)}
ZAHLEN["B4_positive_Eigenwerte"] = {"X": {str(N): v for N, v in eins.items()}, "X'": {str(N): v for N, v in spiegel.items()}}
check("B4 [N, entscheidend] Positivität der verklebten Zustandssumme (Projektion auf X = i^j, alle Spinmoden |s| <= 5 bei "
      f"Y = 0,35, 0,6, 1,0): Für zwei Paare (N = 2, diagonal X_1 X_2) sind genau {eins[2]} Hilbertraumspuren, im Spiegel "
      f"{spiegel[2]} (dieselben Eigenwerte -1 und -i); X = +1 und X' = +1 (die Sektoren mit dem Skalarvakuum R = 1, h = 0) "
      f"geben negative Vielfachheiten. Für N = 4 bleibt nur {eins[4]} bzw. {spiegel[4]} (ungerades j, also G = -1). Mit "
      "normaler Statistik (G = X^2 = +1) ist die Verklebung genau über zwei Paare summierbar und dort genau mit dem "
      "Eigenwert -1",
      eins[2] == ["X = i^2", "X = i^3"] and spiegel[2] == ["X' = i^2", "X' = i^3"]
      and eins[4] == ["X = i^1"] and spiegel[4] == ["X' = i^1"])
check("B5 [N] J-kovariante Fassung: X_p = +1 ist für N = 1..4 eine Hilbertraumspur, X_p = -1 nie. Auf niveaugleichen "
      "Zuständen heißt X_p = +1 aber 2h = 0 mod 4: Die J-kovariante Verklebung ist modulverträglich, weil sie im "
      "unverdrehten Sektor keine Spinoren enthält (die Z4-Fortsetzung von U = +1 aus z113)",
      all(ok(pos[N]["X_p = +1"]) and not ok(pos[N]["X_p = -1"]) for N in NS))

# --- B6 Spintabelle der Verklebung (exakt)
tab_glue, tab_spiegel = {}, {}
for zweih in (-2, -1, 0, 1, 2):                                                 # 2h
    tab_glue[zweih] = [T for T in range(4) if (T + zweih) % 4 == 2]             # X' = i^T i^(2h) = -1
    tab_spiegel[zweih] = [T for T in range(4) if (T - zweih) % 4 == 2]          # X  = i^T i^(-2h) = -1
zsu = all(((2 - zweih) % 4) == tab_glue[zweih][0] for zweih in tab_glue)        # R = Z_SU^-1, Z_SU = i^(2h-2)
ZAHLEN["B6_spintabelle (2h -> Uhrexponent T, R = i^T)"] = {"X' = -1": {str(k): v for k, v in tab_glue.items()},
                                                           "X = -1": {str(k): v for k, v in tab_spiegel.items()}}
check("B6 [E] Was der Eigenwert -1 bedeutet: X' = R i^(2h) = -1 ist R = Z_SU^-1 mit Z_SU = i^(2h-2) (z97, z122), die "
      "Verklebung Z_SU (x) R = +1 aus z117, auf allen Zuständen. Sie gibt die Tabelle h = 0: R = -1 (Skalare, die (10,6)); "
      "h = +-1: R = 1 (Vektoren, die Ströme (45,1) + (1,15)); h = +1/2: R = i (die 16); h = -1/2: R = -i (die 16-quer). "
      "X = -1 ist der Spiegel (16 bei h = -1/2), das Bit P_IF. Mit X^2 = G = U R^2 = +1 gilt die normale Statistik",
      zsu and tab_glue == {-2: [0], -1: [3], 0: [2], 1: [1], 2: [0]}
      and tab_spiegel == {-2: [0], -1: [1], 0: [2], 1: [3], 2: [0]})
befund("Teil B: Die Z4-Verklebung besteht den Modultest von z113 nicht pro Naht (anomal auf einem Paar), aber diagonal über "
       "zwei Paare, und dort genau im Sektor der physikalischen Spintabelle (Vektoren uhrneutral, Skalare R = -1, "
       "Weyl-Fermionen R = +-i) bis auf das Spiegelbit; der Sektor mit dem Skalarvakuum fällt durch. Über vier Paare als "
       "Ganzes geht nur vertauschte Statistik. Neue benannte Prämisse: P_paar (die Statistikregel wird über ein Paar von "
       "Nähten summiert)", True)
say(f"   {zeit()}")

# =================================================================== Teil C
say("\nTeil C: Eintrag ins Suchprogramm und Abstand zur Komplettlösung")
import reconstruction_search as rs                                              # noqa: E402

(HERE / "z174_ergebnis.json").write_text(json.dumps(                            # die Belege, auf die das Suchprogramm zugreift
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": "vorläufig (Teil C läuft)"},
    indent=1, ensure_ascii=False, default=str))
dist = rs.complete_solution_distance()
objs, meta = rs.scenario_candidates()
PAAR = "Z4 glue summed over a pair of seams (P_glue, P_paar)"
alt = [o for o in objs if o["scenario"]["statistics"] == "Z4 glue on the spinor module (P_glue)"
       and {"admissible", "coherent"} <= o["clauses"]]
neu = [o for o in objs if o["scenario"]["statistics"] == PAAR and {"admissible", "coherent"} <= o["clauses"]]
ZAHLEN["C_abstand"] = {"distance": dist["distance"], "open_total": dist["open_total"], "by_question": dist["by_question"],
                       "best": dist["best"], "open_items": sorted(dist["open_items"])}
check("C1 [E] Suchprogramm (OriginLab, reconstruction_search): neuer Achsenwert 'Z4-Verklebung, über ein Paar von Nähten "
      "summiert' mit Beleg B2, B4 und Prämissen P_IF (Konvention) und P_paar; der alte Wert (pro Naht) trägt die Klausel "
      f"modular_statistics weiterhin nicht ({sum('modular_statistics' in o['clauses'] for o in alt)} von {len(alt)}), der "
      f"neue trägt sie ({sum('modular_statistics' in o['clauses'] for o in neu)} von {len(neu)})",
      alt and neu and not any("modular_statistics" in o["clauses"] for o in alt)
      and all("modular_statistics" in o["clauses"] for o in neu))
impl = {(i["replaces"], tuple(i["by"])): i["holds"] for i in meta["implications"]}
check("C2 [E] Implikationen: P_univ -> P_std (z130: nativ in der Partnerlesart) und P_uhr -> P_std (A5, A6) sind "
      "eingetragen und gelten; P_std ist als Prämisse benannt",
      impl.get(("P_univ", ("P_std",))) and impl.get(("P_uhr", ("P_std",))) and "P_std" in rs.PREMISE_MEANING
      and "P_paar" in rs.PREMISE_MEANING)
check("C3 [E] Abstand zur Komplettlösung nach der Runde: 1 fehlende Klausel (graviton), 1 freie Zahl (P_G), 7 Prämissen "
      f"über A* hinaus ({sorted(x for x in dist['open_items'] if x.startswith(('P_', 'I_')))}), dazu K19 und K20: zehn "
      "offene Punkte statt elf. Vorher (z173): 2 Klauseln, 1 freie Zahl, 7 Prämissen",
      dist["distance"] == {"missing_clauses": 1, "free_numbers": 1, "premises_beyond_axiom": 7}
      and dist["open_total"] == 10 and not dist["complete"]
      and sorted(dist["open_items"]) == ["I_step", "P_G", "P_curv", "P_lift", "P_ord", "P_paar", "P_std", "graviton"])
befund("Teil C: Keine Komplettlösung. Die Statistikfrage ist von 'keine Variante besteht' zu 'besteht unter P_paar' "
       "geworden, die Kompositionsfrage von drei Prämissen zu zwei (P_std, P_ord)", True)

urteil = (
    "Schritt 1 ist gerechnet, mit zwei positiven und zwei negativen Befunden. (1) Standardform: Der Standardvektor von Naht "
    "und Tomita-Partner ist ein Transfer-Singulett; der Winkel der Komposition ist auf ihm unsichtbar, und die Gatterfamilie "
    "wirkt nur über ihr Phasengesetz. Die Bilder des Standardvektors liegen genau für phi in (pi/2) Z in der Uhrbahn eines "
    "nativen Standardvektors: für 0 und pi beim bosonischen (kein Gatter, einkörperiger Tausch), für +-pi/2 beim "
    "fermionischen (die Wechselwirkung). Die Viertelphase ist exakt die Klein-Verdrehung der gradierten Standardform. Damit "
    "fallen P_univ (z130) und P_uhr in eine Prämisse P_std zusammen: Die Komposition ist ein Morphismus der gradierten "
    "Standardform. Negativ: Kein Mitglied der Familie, das den Standardvektor festlässt, wechselwirkt, und die einzige nicht "
    "faktorisierte Symmetrie der Standardform ist der einkörperige Seitentausch. Die Wechselwirkung liegt also nicht in der "
    "Standardform einer Naht, sie braucht P_std als Zusatz; P_ord bleibt unberührt. (2) Statistik: Der in z113 offene "
    "Modultest der Z4-Verklebung ist gerechnet. Pro Naht ist sie anomal. Diagonal über zwei Paare ist sie modulinvariant "
    "und positiv, und zwar genau im Sektor R = Z_SU^-1 mit normaler Statistik: Vektoren uhrneutral, Skalare bei R = -1, "
    "Weyl-Fermionen bei R = +-i, bis auf das Spiegelbit. Der Sektor mit dem Skalarvakuum fällt durch. Die J-kovariante "
    "Fassung besteht immer, enthält aber keine Spinoren. Die diagonale Vervollständigung einer einzelnen Naht trägt die "
    "Spinorstatistik damit nicht; nötig ist die Summe über ein Paar (P_paar). Abstand zur Komplettlösung: 1 Klausel, 1 "
    "freie Zahl, 7 Prämissen, dazu K19 und K20, zehn Punkte statt elf. Keine Komplettlösung.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z174_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil},
                                                    indent=1, ensure_ascii=False, default=str))
(HERE / "ausgabe_z174.txt").write_text("\n".join(x for x in OUT if "Laufzeit bisher" not in x) + "\n")
say(f"{zeit()}")
