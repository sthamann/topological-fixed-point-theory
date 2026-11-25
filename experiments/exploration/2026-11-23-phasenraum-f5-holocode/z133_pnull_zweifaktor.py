"""Z133: P_null in der Zwei-Faktoren-Architektur von z82 -- welches Transfergatter respektiert die native Struktur?

Folgeprüfung zu z120/z121 (Empfehlung V1). Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour,
Massen, alpha, 16, E8) bleiben unberührt, nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt),
L Lesart, N numerisch, O offen.

Bestand: z82 (innen Ladungsgitter + Uhr, außen Oszillatoren; Hülle universell, KMS-Zustand ein Produkt, Sugawara
L_n = L_n^osc + mu_Q J_n), z92 (E^nat = L_0^osc + Q^2/2, e^(2 pi i E^nat) = R^2), z106 (Klein-Kette, Gamma(pi)), z119
(Q^2/2 bricht die Transfer-SU(2); relativer Transferrahmen über die Nahtuhr beobachtbar), z120/z121 (U = B^+ e^(-i phi H_Q) B,
B = e^(-i pi J_x/2); nativ E^nat einkörperig, Gatter trivial; unter P_null U = P_T - i P_S, U^2 = -Gamma(pi), u = tan(phi/2);
[U, L^nat_(+-1)] = 0), z125 (P_univ; Ladungsenergie auf N = 3 nicht erhalten).

Modell: zwei Nähte A, B, je ein chirales Dirac-Fermion des Quarters alpha = 1/4 (Lesart 'zwei Stufen eines Fermions', z121 b),
Moden y in {-M, ..., M-1} (x = y + 1/4), See y < 0 gefüllt. Behalten werden alle Zustände mit E^nat_A + E^nat_B <= 4; jede
Anregung am Fensterrand kostet >= M - 1/2 = 4.5, also ist der behaltene Raum exakt der native (kein Abschneidefehler), und er
ist unter Transfer-Nullmode, E^nat, H_Q, Ladungen und Uhr abgeschlossen. Zwei-Faktoren-Zerlegung je Naht über das
Maya-Diagramm: innen Q (Gitter, Q^2/2, Uhr i^Q), außen die Partition lambda (Oszillatoren, L_0^osc = |lambda|).

Teil 1  Innen/außen: Maya-Zerlegung, E^nat = Q^2/2 + |lambda|, Uhr oszillatorblind; die drei Gatter aus E^nat, H_Q = Q^2/2 und
        L_0^osc; U_in U_out = U_full, U_full = e^(-i pi E^nat/2) ein Produkt von Einzelnahtphasen.
Teil 2  Fünf Strukturtests für U_full, U_in, U_out und das rein innere Gatter U_rein (Wurzel des Ladungstauschs, Oszillatoren
        bleiben): (a) Faktorzerlegung (Energieaustausch innen/außen, Innen-außen-Verschränkung, Q_A = Q_B-Sektor),
        (b) Uhr (global, relativ, je Naht) und Tomita-J (reelle Struktur), (c) KMS bei beta = 4 pi, (d) native sl(2)
        (L^nat_(+-1), E^nat), (e) Lokalität je Nahtpaar (drei Nähte). Satz: rein innen + [., L^nat_(+-1)] = 0 => diagonal.
Teil 3  P_viertel: Kandidaten Uhrordnung, Viertelschritt, Modulperiode beta = 4 pi, Selbstdualität, U^2 = -Gamma(pi).
Teil 4  Bilanz und Urteil (P_null, P_viertel, P_ord).
"""
import itertools
import json
import os
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(133)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def nrm(A):
    return float(np.linalg.norm(A))


def bis_auf_phase(A, B):
    """min_phi ||A - e^(i phi) B||."""
    ph = np.vdot(B.ravel(), A.ravel())
    ph = ph/abs(ph) if abs(ph) > 1e-14 else 1.0
    return nrm(A - ph*B)


def komm(A, B):
    return A @ B - B @ A


# =================================================================== Fockfenster
M = 5                                   # Moden je Naht: y = -M .. M-1
NM = 2*M                                # Moden je Naht
EMAX = 4.0


def modus(s, y):
    return s*NM + (y + M)


def naht_daten(occ):
    """Ladung Q, Energie E^nat und Partition lambda einer Naht aus ihrer Besetzung (Liste über y = -M .. M-1)."""
    Q = sum(occ[y + M] for y in range(0, M)) - sum(1 - occ[y + M] for y in range(-M, 0))
    E = sum((y + 0.5)*occ[y + M] for y in range(0, M)) + sum((-y - 0.5)*(1 - occ[y + M]) for y in range(-M, 0))
    besetzt = sorted((y for y in range(-M, M) if occ[y + M]), reverse=True)
    lam = [besetzt[i] + (i + 1) - Q for i in range(len(besetzt))]
    lam = tuple(p for p in lam if p != 0)
    return Q, E, lam


def fenster(nahte, emax):
    """Alle Besetzungen von nahte Nähten mit sum E^nat <= emax (Bitmaske über nahte*NM Moden, Bit k = Modus k)."""
    je = []
    for bits in itertools.product((0, 1), repeat=NM):
        Q, E, lam = naht_daten(bits)
        if E <= emax + 1e-9:
            je.append((bits, Q, E, lam))
    zust = []
    for kombi in itertools.product(je, repeat=nahte):
        if sum(k[2] for k in kombi) <= emax + 1e-9:
            maske = 0
            for s, k in enumerate(kombi):
                for j, b in enumerate(k[0]):
                    maske |= b << (s*NM + j)
            zust.append((maske, tuple((k[1], k[3]) for k in kombi), sum(k[2] for k in kombi)))
    zust.sort(key=lambda z: (round(z[2]*4), z[0]))
    return zust


def bilinear(zust, idx, i, j, koeff=1.0):
    """koeff * c+_i c_j (Jordan-Wigner nach Modusindex) auf dem behaltenen Raum; Bilder außerhalb werden gemeldet."""
    n = len(zust)
    A = np.zeros((n, n))
    raus = 0
    for c, (m, _, _) in enumerate(zust):
        if not (m >> j) & 1:
            continue
        m2 = m ^ (1 << j)
        s = bin(m2 & ((1 << j) - 1)).count("1")
        if (m2 >> i) & 1:
            continue
        s += bin(m2 & ((1 << i) - 1)).count("1")
        m3 = m2 | (1 << i)
        r = idx.get(m3)
        if r is None:
            raus += 1
            continue
        A[r, c] += koeff*(-1)**s
    return A, raus


ZU = fenster(2, EMAX)
IDX = {m: k for k, (m, _, _) in enumerate(ZU)}
DIM = len(ZU)
QA = np.diag([float(z[1][0][0]) for z in ZU])
QB = np.diag([float(z[1][1][0]) for z in ZU])
LA = np.diag([float(sum(z[1][0][1])) for z in ZU])          # |lambda_A| = L_0^osc der Naht A
LB = np.diag([float(sum(z[1][1][1])) for z in ZU])
ENAT = np.diag([z[2] for z in ZU])
HQ = (QA @ QA + QB @ QB)/2
L0OSC = LA + LB
QT = QA + QB
dQT = np.diag(QT).round().astype(int)
dE = np.diag(ENAT)
I_ = np.eye(DIM)


def transfer(shift=0):
    """J^+ = sum_y c+_(A,y+shift) c_(B,y) (shift = 0: Nullmode des su(2)_1-Stroms)."""
    Jp = np.zeros((DIM, DIM))
    raus = 0
    for y in range(-M, M):
        if -M <= y + shift < M:
            A_, r_ = bilinear(ZU, IDX, modus(0, y + shift), modus(1, y))
            Jp += A_
            raus += r_
    return Jp, raus


JP, raus_J = transfer()
JX, JY = (JP + JP.T)/2, (JP - JP.T)/(2j)
JZ = (QA - QB)/2


def L_heber(s):
    """L^nat_+ = sum_y (y + 1) c+_(y+1) c_y der Naht s (z92/z121, Amplitude p + 1)."""
    Lp = np.zeros((DIM, DIM))
    for y in range(-M, M - 1):
        A_, _ = bilinear(ZU, IDX, modus(s, y + 1), modus(s, y), float(y + 1))
        Lp += A_
    return Lp


LP = L_heber(0) + L_heber(1)
LM = LP.T
TIEF = dE <= EMAX - 1 + 1e-9                                       # Spalten, deren Bild unter L^nat_+ im Raum liegt
LABEL = {(z[1][0][0], z[1][0][1], z[1][1][0], z[1][1][1]): k for k, z in enumerate(ZU)}
INNEN = sorted(set((z[1][0][0], z[1][1][0]) for z in ZU))
AUSSEN = sorted(set((z[1][0][1], z[1][1][1]) for z in ZU))
I_POS = {q: k for k, q in enumerate(INNEN)}
A_POS = {l: k for k, l in enumerate(AUSSEN)}
EA = np.diag([z[1][0][0]**2/2 + sum(z[1][0][1]) for z in ZU])
EB = ENAT - EA
SEKTOREN = sorted(set((int(dQT[k]), round(float(dE[k])*2)/2) for k in range(DIM)))


def sektor(qt, e):
    return [k for k in range(DIM) if dQT[k] == qt and abs(dE[k] - e) < 1e-9]


def innen_aussen_rang(v, tol=1e-10):
    """Schmidt-Rang eines Zustands über den Schnitt innen (Q_A, Q_B) | außen (lambda_A, lambda_B)."""
    Mt = np.zeros((len(INNEN), len(AUSSEN)), complex)
    for k, z in enumerate(ZU):
        Mt[I_POS[(z[1][0][0], z[1][1][0])], A_POS[(z[1][0][1], z[1][1][1])]] = v[k]
    sv = np.linalg.svd(Mt, compute_uv=False)
    return int(np.sum(sv > tol*max(sv.max(), 1e-300)))


def zustand(qa, la, qb, lb):
    v = np.zeros(DIM, complex)
    v[LABEL[(qa, tuple(la), qb, tuple(lb))]] = 1
    return v


def gatter(H, phi=np.pi/2, chi=0.0, theta=np.pi/2):
    """U = B^+ e^(-i phi H) B mit B = e^(-i theta (cos chi J_x + sin chi J_y))."""
    Bm = expm(-1j*theta*(np.cos(chi)*JX + np.sin(chi)*JY))
    return Bm.conj().T @ expm(-1j*phi*H) @ Bm


# =================================================================== Teil 1
say("Teil 1: Innen/außen im nativen Fockfenster")
p_zahl = [1, 1, 2, 3, 5]
erw = 0
for qa in range(-3, 4):
    for qb in range(-3, 4):
        rest = EMAX - (qa*qa + qb*qb)/2
        if rest < -1e-9:
            continue
        erw += sum(p_zahl[na]*p_zahl[nb] for na in range(5) for nb in range(5) if na + nb <= rest + 1e-9)
maya_fehler = max(abs(z[2] - sum(q*q/2 + sum(l) for q, l in z[1])) for z in ZU)
RA = np.diag([1j**int(z[1][0][0]) for z in ZU])
RB = np.diag([1j**int(z[1][1][0]) for z in ZU])
uhr_L0 = nrm(expm(2j*np.pi*(EA - QA/4)) - RA) + nrm(expm(2j*np.pi*(EB - QB/4)) - RB)
uhr_blind = nrm(expm(2j*np.pi*(L0OSC)) - I_)
r2_stat = nrm(expm(2j*np.pi*ENAT) - RA @ RA @ RB @ RB)
LPC = komm(LP, LM)
sl2 = nrm((LPC + 2*ENAT)[:, TIEF])
ZAHLEN["fenster"] = {"M (Moden je Naht 2M)": M, "E_max": EMAX, "Zustände": DIM, "erwartet aus Partitionen": erw,
                     "Transfer verlässt den Raum (Bilder)": raus_J, "Maya E^nat - Q^2/2 - |lambda|": maya_fehler,
                     "Uhr e^(2 pi i L_0) = i^Q": uhr_L0, "e^(2 pi i L_0^osc) = 1": uhr_blind,
                     "e^(2 pi i E^nat) = R^2 (z92)": r2_stat, "[L^nat_+, L^nat_-] + 2 E^nat (E <= 3)": sl2,
                     "Innenetiketten": len(INNEN), "Außenetiketten": len(AUSSEN), "Sektoren (Q_tot, E)": len(SEKTOREN)}
check(f"[E] Natives Fockfenster zweier Nähte (je ein Dirac-Fermion, alpha = 1/4, Moden y = -{M}..{M - 1}, alle Zustände mit "
      f"E^nat <= {EMAX:g}; Randanregungen kosten >= {M - 0.5}): {DIM} Zustände = Zählung aus Partitionen ({erw}); der Transfer "
      f"J^+ = sum_y c+_(A,y) c_(B,y) verlässt den Raum nie ({raus_J} Bilder). Maya-Zerlegung je Naht (innen Q, außen lambda): "
      f"E^nat = Q^2/2 + |lambda| exakt ({maya_fehler:.0e}), {len(INNEN)} Innen- und {len(AUSSEN)} Außenetiketten. Die Uhr "
      f"e^(2 pi i L_0) = i^Q je Naht ({uhr_L0:.0e}) ist oszillatorblind (e^(2 pi i L_0^osc) = 1, {uhr_blind:.0e}), "
      f"e^(2 pi i E^nat) = R_A^2 R_B^2 (z92, {r2_stat:.0e}); die native sl(2): [L^nat_+, L^nat_-] = -2 E^nat ({sl2:.0e})",
      DIM == erw and raus_J == 0 and maya_fehler < 1e-12 and uhr_L0 < 1e-10 and uhr_blind < 1e-10 and r2_stat < 1e-10
      and sl2 < 1e-10)

B0 = expm(-1j*np.pi/2*JX)
U_FULL = gatter(ENAT)
U_IN = gatter(HQ)
U_OUT = gatter(L0OSC)
full_diag = nrm(U_FULL - expm(-1j*np.pi/2*ENAT))
produkt = nrm(U_IN @ U_OUT - U_FULL)
inv_beziehung = nrm(U_OUT - U_FULL @ U_IN.conj().T)
einkoerper = nrm(U_FULL - expm(-1j*np.pi/2*EA) @ expm(-1j*np.pi/2*EB))
transfer_in = nrm(komm(U_IN, QA))
transfer_out = nrm(komm(U_OUT, QA))
ZAHLEN["drei_gatter"] = {"U_full - e^(-i pi E^nat/2)": full_diag, "U_in U_out - U_full": produkt,
                         "U_out - U_full U_in^+": inv_beziehung, "U_full - e^(-i pi E_A/2) e^(-i pi E_B/2)": einkoerper,
                         "||[U_in, Q_A]|| (Transfer)": transfer_in, "||[U_out, Q_A]||": transfer_out}
check(f"[E] Drei Gatter U_H = B^+ e^(-i pi H/2) B (B = e^(-i pi J_x/2)) für H = E^nat, H_Q = Q^2/2 (innen), L_0^osc (außen). "
      f"Weil der Transfer mit E^nat vertauscht, ist U_full = e^(-i pi E^nat/2) ({full_diag:.0e}) = Produkt der "
      f"Einzelnahtphasen ({einkoerper:.0e}): trivial. Weil H_Q und L_0^osc vertauschen und E^nat = H_Q + L_0^osc, gilt exakt "
      f"U_in U_out = U_full ({produkt:.0e}), also U_out = U_full U_in^(-1) ({inv_beziehung:.0e}): Innen- und Außengatter sind "
      f"bis auf die freie Einzelnahtphase zueinander invers (die zwei Wurzelvorzeichen von z121 = Zeitrichtungen). Beide "
      f"übertragen Ladung ([U, Q_A]: {transfer_in:.2f}, {transfer_out:.2f})",
      full_diag < 1e-10 and produkt < 1e-10 and inv_beziehung < 1e-10 and einkoerper < 1e-10 and transfer_in > 0.1
      and transfer_out > 0.1)
befund("E: Innen gegen außen ist keine Wahl zwischen trivial und nichttrivial, sondern zwischen zwei zueinander inversen "
       "Gattern: U_out = U_full U_in^(-1). Jeder Test, der freie Einzelnahtphasen und Zeitumkehr nicht sieht, behandelt P_null "
       "und 'P_osc' gleich; unterscheiden kann nur ein Test, der auf die freie Phase e^(-i pi E^nat/2) reagiert", True)

# rein inneres Gatter: Wurzel des Ladungstauschs Q_A <-> Q_B bei festen Oszillatoren (Vorzeichen +1 bzw. Klein (-1)^(Q_A Q_B))
S_IN = np.zeros((DIM, DIM))
S_IN_K = np.zeros((DIM, DIM))
for k, z in enumerate(ZU):
    (qa, la), (qb, lb) = z[1]
    r = LABEL[(qb, la, qa, lb)]
    S_IN[r, k] = 1
    S_IN_K[r, k] = (-1)**(qa*qb)
PH_QT = np.diag(np.exp(-1j*np.pi*np.diag(QT)**2/8))
U_REIN = PH_QT @ ((1 - 1j)*I_ + (1 + 1j)*S_IN)/2
U_REIN_K = PH_QT @ ((1 - 1j)*I_ + (1 + 1j)*S_IN_K)/2
GATTER = {"U_full": U_FULL, "U_in": U_IN, "U_out": U_OUT, "U_rein": U_REIN}
ZAHLEN["gatter_unitaer"] = {k: nrm(G @ G.conj().T - I_) for k, G in list(GATTER.items()) + [("U_rein (Klein)", U_REIN_K)]}

# =================================================================== Teil 2
say("\nTeil 2: Welches Gatter respektiert die native Struktur?")
say("  (a) Faktorzerlegung innen/außen")
aust = {}
for nm, G in GATTER.items():
    je_sektor = {}
    for qt, e in SEKTOREN:
        s_ = sektor(qt, e)
        v_ = nrm(komm(G, HQ)[np.ix_(s_, s_)])
        if v_ > 1e-10:
            je_sektor[f"Q_tot={qt}, E={e:g}"] = round(v_, 6)
    aust[nm] = {"||[U, H_Q]||": nrm(komm(G, HQ)), "Sektoren mit Austausch": je_sektor}
# Beispiel: Transferspin 3/2 bei Q_tot = 1, E = 5/2, Eingang |2, -1; leer, leer>
v32 = U_IN @ zustand(2, (), -1, ())
hq_vorher = 9/4 + 1/4
hq_nachher = float(np.real(np.vdot(v32, HQ @ v32)))
osc_nachher = float(np.real(np.vdot(v32, L0OSC @ v32)))
ZAHLEN["energieaustausch"] = {**{k: v for k, v in aust.items()},
                              "U_in |2,-1; leer> : <H_Q> vorher, nachher, <L_0^osc> nachher": [hq_vorher, hq_nachher, osc_nachher]}
ungerade_nur = all(int(k.split(",")[0].split("=")[1]) % 2 for k in aust["U_in"]["Sektoren mit Austausch"])
check("[E] Energieaufteilung E^nat = H_Q + L_0^osc: U_full und U_rein erhalten H_Q (und damit L_0^osc) exakt "
      f"({aust['U_full']['||[U, H_Q]||']:.0e}, {aust['U_rein']['||[U, H_Q]||']:.0e}); U_in und U_out erhalten E^nat, "
      f"tauschen aber Energie zwischen Ladungsgitter und Oszillatoren aus ({aust['U_in']['||[U, H_Q]||']:.3f}), genau in den "
      f"Sektoren ungerader Gesamtladung mit Transferspin >= 3/2: {aust['U_in']['Sektoren mit Austausch']}. Beispiel: U_in auf "
      f"|Q_A = 2, Q_B = -1; Oszillatorvakuum> (H_Q = {hq_vorher:g}, L_0^osc = 0) gibt <H_Q> = {hq_nachher:.4f}, <L_0^osc> = "
      f"{osc_nachher:.4f}: das aus der Gitterenergie gebaute Gatter pumpt Gitterenergie in die Oszillatoren (z125: H_Q auf N = 3 "
      "nicht erhalten; hier als Austausch innen -> außen bei erhaltenem E^nat)",
      aust["U_full"]["||[U, H_Q]||"] < 1e-10 and aust["U_rein"]["||[U, H_Q]||"] < 1e-10 and aust["U_in"]["||[U, H_Q]||"] > 0.1
      and aust["U_out"]["||[U, H_Q]||"] > 0.1 and ungerade_nur and len(aust["U_in"]["Sektoren mit Austausch"]) > 0
      and abs(hq_nachher + osc_nachher - hq_vorher) < 1e-10 and osc_nachher > 0.1)

raenge = {}
for nm, G in GATTER.items():
    raenge[nm] = max(innen_aussen_rang(G[:, k]) for k in range(DIM))
bsp = U_IN @ zustand(2, (1,), 0, ())
bsp_rang = innen_aussen_rang(bsp)
bsp_gewichte = sorted((round(float(abs(bsp[k])**2), 6), ZU[k][1]) for k in range(DIM) if abs(bsp[k]) > 1e-9)
# Kontrolle: der native Heber L^nat_+ selbst verschränkt Gitter und Oszillatoren (Kopplung Q J_(-1), z82)
v0 = (zustand(1, (), 0, ()) + zustand(2, (), 0, ()))/np.sqrt(2)
v_sl2 = v0 + 0.3*(LP @ v0)
rang_sl2 = innen_aussen_rang(v_sl2)
lin_Q = []
for q in (0, 1, 2):
    a_, b_ = LABEL[(q, (), 0, ())], LABEL[(q, (1,), 0, ())]
    lin_Q.append(float(LP[b_, a_]))
ZAHLEN["innen_aussen_verschraenkung"] = {"max Schmidt-Rang über Basis-Eingänge": raenge,
                                          "U_in |2,0; (1), leer>: Rang, Gewichte": [bsp_rang, [[w, str(l)] for w, l in bsp_gewichte]],
                                          "Kontrolle (1 + 0.3 L^nat_+)(|1,0> + |2,0>) x Vakuum: Rang": rang_sl2,
                                          "<Q, (1)| L^nat_+ |Q, leer> für Q = 0, 1, 2": lin_Q}
check(f"[E] Innen-außen-Verschränkung aus Produkteingängen |Q_A, Q_B> x |lambda_A, lambda_B> (größter Schmidt-Rang): {raenge}. "
      f"U_in auf |2, 0; (1), leer> gibt Rang {bsp_rang} ({bsp_gewichte}): Ladungsgitter und Oszillatoren werden verschränkt. "
      f"Kontrolle, dass das strenge Produkt V (x) W keine native Forderung ist: schon der native Heber koppelt beide "
      f"(<Q,(1)|L^nat_+|Q,leer> = {lin_Q} = Q, die Kopplung Q J_(-1) von z82), und (1 + 0.3 L^nat_+)(|1,0> + |2,0>) x Vakuum hat "
      f"Rang {rang_sl2}",
      raenge["U_full"] == 1 and raenge["U_rein"] == 1 and raenge["U_in"] >= 2 and raenge["U_out"] >= 2 and bsp_rang == 2
      and rang_sl2 == 2 and np.allclose(lin_Q, [0, 1, 2]))

# Gerade Gesamtladung: U_in = e^(-i pi Q_tot^2/8) ((1 - i) + (1 + i) e^(i pi J_y))/2 (Wurzel des vollen Nahttauschs)
GER = [k for k in range(DIM) if dQT[k] % 2 == 0]
UNG = [k for k in range(DIM) if dQT[k] % 2 == 1]
RY = expm(1j*np.pi*JY)
form_ger = nrm((U_IN - PH_QT @ ((1 - 1j)*I_ + (1 + 1j)*RY)/2)[np.ix_(GER, GER)])
form_ung = nrm((U_IN - np.exp(-1j*np.pi/8)*PH_QT @ (expm(0.5j*np.pi*JY) + expm(-0.5j*np.pi*JY))/np.sqrt(2))[np.ix_(UNG, UNG)])
# Tausch = Tausch innen (x) Tausch außen: |e^(i pi J_y)| ist die Permutation (Q_A, l_A; Q_B, l_B) -> (Q_B, l_B; Q_A, l_A)
S_VOLL = np.zeros((DIM, DIM))
for k, z in enumerate(ZU):
    (qa, la), (qb, lb) = z[1]
    S_VOLL[LABEL[(qb, lb, qa, la)], k] = 1
tausch_voll = nrm(np.abs(RY[np.ix_(GER, GER)]) - S_VOLL[np.ix_(GER, GER)])
DIAG = [k for k in range(DIM) if ZU[k][1][0][0] == ZU[k][1][1][0]]
leck_diag = nrm(np.delete(U_IN[:, DIAG], DIAG, axis=0))
# auf Q_A = Q_B: Wurzel des Oszillatortauschs lambda_A <-> lambda_B (mit Vorzeichen), innen die Eins
S_OSZ = (np.abs(RY) > 0.5)*np.sign(RY.real)
wurzel_osz = nrm((U_IN - PH_QT @ ((1 - 1j)*I_ + (1 + 1j)*S_OSZ)/2)[np.ix_(DIAG, DIAG)])
osz_permutation = nrm(np.abs(RY[np.ix_(DIAG, DIAG)]) - S_VOLL[np.ix_(DIAG, DIAG)])
ZAHLEN["struktur_U_in"] = {"gerade Q_tot: U_in - Wurzel von e^(i pi J_y)": form_ger,
                           "ungerade Q_tot: U_in - e^(-i pi/8) sqrt2 cos(pi J_y/2)": form_ung,
                           "|e^(i pi J_y)| = voller Nahttausch (innen x außen)": tausch_voll,
                           "Leckage aus Q_A = Q_B": leck_diag, "Q_A = Q_B: U_in - Wurzel des Oszillatortauschs": wurzel_osz,
                           "Q_A = Q_B: |e^(i pi J_y)| = Tausch lambda_A <-> lambda_B": osz_permutation}
check(f"[E] Struktur des Innengatters auf dem ganzen nativen Raum: U_in = e^(-i pi Q_tot^2/8) e^(-i pi J_y^2/2). Gerade "
      f"Gesamtladung (ganzer Transferspin): U_in = e^(-i pi Q_tot^2/8) ((1 - i) + (1 + i) e^(i pi J_y))/2 ({form_ger:.0e}), und "
      f"e^(i pi J_y) ist bis auf Vorzeichen der volle Nahttausch (Q_A, lambda_A) <-> (Q_B, lambda_B) ({tausch_voll:.0e}): das "
      f"Gatter ist die Wurzel des Tauschs von Gitter UND Oszillatoren. Ungerade Gesamtladung: U_in = e^(-i pi Q_tot^2/8) "
      f"e^(-i pi/8) sqrt 2 cos(pi J_y/2) ({form_ung:.0e}), Summe der Vierteldrehungen um y (Spin 1/2 trivial, ab 3/2 "
      f"Gitter-Oszillator-Mischung). Auf Q_A = Q_B (z120s physikalischer Sektor, alle Oszillatorinhalte) ist U_in leckfrei "
      f"({leck_diag:.0e}) und gleich 1_innen (x) Wurzel des Oszillatortauschs lambda_A <-> lambda_B ({wurzel_osz:.0e}, "
      f"Permutation {osz_permutation:.0e}): dort wirkt das 'innere' Gatter nur auf den Außenfaktor",
      form_ger < 1e-10 and form_ung < 1e-10 and tausch_voll < 1e-10 and leck_diag < 1e-10 and wurzel_osz < 1e-10
      and osz_permutation < 1e-10)
befund("E: Das Innengatter respektiert die Zwei-Faktoren-Zerlegung nicht. Der native Transfer ist die Nullmode eines "
       "Vertexoperators (Frenkel-Kac: die Wurzeln J_z = +-1 sind Gittervektoren, das Cartan-Mitglied J_z = 0 des Tripletts ist "
       "eine Oszillatoranregung), also mischt jede Vierteldrehung Gitter und Oszillatoren. Ergebnis: bei ungerader "
       "Gesamtladung Energieaustausch Gitter <-> Oszillatoren, bei gerader Verschränkung beider Faktoren (Wurzel des vollen "
       "Tauschs), und im physikalischen Sektor Q_A = Q_B wirkt es als 1 (x) sqrt(SWAP) auf die OSZILLATOREN. Das z120-Gatter "
       "'sieht' das Gitter nur in den virtuellen Zwischenzuständen (2,0), (0,2); seine Wirkung liegt im Außenfaktor. Die "
       "Faktorzerlegung respektieren nur U_full (trivial) und das rein innere U_rein", True)
w_fk = JP @ zustand(-1, (), 1, ())
fk = [(ZU[k][1], round(float(abs(w_fk[k])**2), 6)) for k in range(DIM) if abs(w_fk[k]) > 1e-9]
ZAHLEN["frenkel_kac"] = {"J^+ |Q_A=-1, Q_B=1; Vakuum>": [[str(a), b] for a, b in fk]}
check(f"[E] Frenkel-Kac-Kontrolle: J^+ auf |Q_A = -1, Q_B = 1; Oszillatorvakuum> (J_z = -1, E^nat = 1) gibt {fk}: das "
      "J_z = 0-Mitglied des Transfer-Tripletts liegt bei Q_A = Q_B = 0 mit einer Oszillatoranregung |lambda| = 1 (relativer "
      "Oszillator), nicht im Gitter. Die Transfer-SU(2) ist keine Symmetrie des Ladungsgitters allein",
      len(fk) == 2 and all(a[0][0] == 0 and a[1][0] == 0 and sum(a[0][1]) + sum(a[1][1]) == 1 for a, _ in fk)
      and abs(sum(b for _, b in fk) - 2) < 1e-10)

say("  (b) Uhr und Tomita-J")
R_GLOB, R_REL = RA @ RB, RA @ RB.conj()
uhr = {}
for nm, G in GATTER.items():
    uhr[nm] = {"global R_A R_B": nrm(komm(G, R_GLOB)), "relativ R_A R_B^-1": nrm(komm(G, R_REL)),
               "je Naht R_A": nrm(komm(G, RA)), "je Naht R_A^2 = (-1)^Q_A": nrm(komm(G, RA @ RA)),
               "R_A auf Q_A = Q_B": nrm((RA @ G @ RA.conj().T - G)[np.ix_(DIAG, DIAG)])}
U_CHI = {c: gatter(HQ, chi=c) for c in (0.0, 0.3, np.pi/2)}
rahmen = nrm(RA @ U_IN @ RA.conj().T - U_CHI[np.pi/2])
ZAHLEN["uhr"] = uhr | {"R_A U_in(chi = 0) R_A^+ - U_in(chi = pi/2)": rahmen}
rein_ung = nrm(komm(U_REIN, R_REL)[np.ix_(UNG, UNG)])
rein_ger = nrm(komm(U_REIN, R_REL)[np.ix_(GER, GER)])
ZAHLEN["uhr"]["U_rein: relative Uhr auf gerader / ungerader Gesamtladung"] = [rein_ger, rein_ung]
nat3 = ("U_full", "U_in", "U_out")
check(f"[E] Uhr: U_full, U_in und U_out vertauschen mit der globalen Uhr R_A R_B, der relativen Uhr R_A R_B^-1 = e^(i pi J_z) "
      f"(z125-Paaruhr) und der Nahtparität R_A^2 (P_univ) (max "
      f"{max(max(uhr[k]['global R_A R_B'], uhr[k]['relativ R_A R_B^-1'], uhr[k]['je Naht R_A^2 = (-1)^Q_A']) for k in nat3):.0e}). "
      f"Das rein innere U_rein verletzt relative Uhr und Nahtparität ({uhr['U_rein']['relativ R_A R_B^-1']:.3f}; gerade "
      f"Gesamtladung {rein_ger:.0e}, ungerade {rein_ung:.3f}): Der Ladungstausch bei festen Oszillatoren vertauscht die "
      f"Nahtparitäten. Mit der Nahtuhr R_A = i^Q_A vertauscht nur U_full ({uhr['U_full']['je Naht R_A']:.0e}); U_in, U_out "
      f"({uhr['U_in']['je Naht R_A']:.3f}) nicht: R_A dreht den Transferrahmen um eine Viertel (R_A U_in(chi) R_A^+ = "
      f"U_in(chi + pi/2), {rahmen:.0e}); auf Q_A = Q_B ist U_in uhrinvariant ({uhr['U_in']['R_A auf Q_A = Q_B']:.0e}). Das "
      "Innengatter braucht also außerhalb des physikalischen Sektors einen gewählten Transferrahmen (z119: über die Nahtuhr "
      "beobachtbar)",
      all(uhr[k]["global R_A R_B"] < 1e-10 and uhr[k]["relativ R_A R_B^-1"] < 1e-10 and uhr[k]["je Naht R_A^2 = (-1)^Q_A"] < 1e-10
          for k in nat3)
      and uhr["U_rein"]["global R_A R_B"] < 1e-10 and rein_ger < 1e-10 and rein_ung > 0.1
      and uhr["U_full"]["je Naht R_A"] < 1e-10 and uhr["U_in"]["je Naht R_A"] > 0.1
      and uhr["U_out"]["je Naht R_A"] > 0.1 and rahmen < 1e-10 and uhr["U_in"]["R_A auf Q_A = Q_B"] < 1e-10)

tom = {nm: nrm(G - G.T) for nm, G in GATTER.items()}
tom_chi = {f"{c:.4f}": nrm(G - G.T) for c, G in U_CHI.items()}
zeitumkehr = {nm: nrm(G.conj() - G.conj().T) for nm, G in GATTER.items()}
ZAHLEN["tomita"] = {"J U J - U^-1 (J = K, reelle Besetzungsbasis)": tom, "Transferrahmen chi": tom_chi}
check(f"[L/E] Tomita-J als reelle Struktur der KMS-Zustände in der reellen Besetzungsbasis (J = K; J Delta^(1/2) bildet U auf "
      f"U^+ ab genau dann, wenn U^T = U): J U J = U^-1 gilt für alle vier Gatter ({ {k: f'{v:.0e}' for k, v in tom.items()} }): "
      f"jedes ist zeitumkehrkovariant (die Zeitumkehr ist die Inverse, das zweite Wurzelvorzeichen). Kontrolle: im "
      f"Transferrahmen chi = 0.3 verletzt ({tom_chi['0.3000']:.3f}), chi = pi/2 (Uhrbild) erfüllt ({tom_chi['1.5708']:.0e}): "
      "J und Uhr zusammen lassen genau die Uhrbahn chi in (pi/2) Z der Transferrahmen zu",
      all(v < 1e-10 for v in tom.values()) and tom_chi["0.3000"] > 0.1 and tom_chi["1.5708"] < 1e-10
      and tom_chi["0.0000"] < 1e-10)

say("  (c) KMS bei beta = 4 pi")
RHO = expm(-4*np.pi*ENAT)
RHO /= np.trace(RHO)
kms = {nm: nrm(komm(G, RHO))/nrm(RHO) for nm, G in GATTER.items()}
Xz = rng.normal(size=(DIM, DIM))
Xz = (Xz + Xz.T)/2
U_ZUF = expm(-0.3j*Xz)
kms_kontr = nrm(komm(U_ZUF, RHO))/nrm(RHO)
kms_theta = nrm(komm(gatter(HQ, theta=np.pi/3), RHO))/nrm(RHO)
rho_in = np.diag([sum(RHO[k, k] for k in range(DIM) if (ZU[k][1][0][0], ZU[k][1][1][0]) == q) for q in INNEN]).real
rho_aus = np.diag([sum(RHO[k, k] for k in range(DIM) if (ZU[k][1][0][1], ZU[k][1][1][1]) == l) for l in AUSSEN]).real
produkt_kms = max(abs(RHO[k, k] - rho_in[I_POS[(z[1][0][0], z[1][1][0])], I_POS[(z[1][0][0], z[1][1][0])]]
                      * rho_aus[A_POS[(z[1][0][1], z[1][1][1])], A_POS[(z[1][0][1], z[1][1][1])]])
                  for k, z in enumerate(ZU))/RHO.max()
ZAHLEN["kms"] = {"||[U, rho_4pi]||/||rho||": kms, "Kontrolle zufällig": kms_kontr, "theta = pi/3 (undicht)": kms_theta,
                 "rho = rho_innen x rho_außen (relativ, im Fenster)": float(produkt_kms)}
check(f"[E] KMS bei beta = 4 pi (Gibbs von E^nat): alle vier Gatter lassen den Zustand fest ({ {k: f'{v:.0e}' for k, v in kms.items()} }), "
      f"ebenso das undichte theta = pi/3 ({kms_theta:.0e}); Kontrolle e^(-0.3 i X) mit zufälligem X: {kms_kontr:.3f}. Der "
      f"Zustand ist im Fenster ein Produkt innen x außen bis auf die Fensterabschneidung (relativ {produkt_kms:.0e}, z82). Der "
      "KMS-Test trennt die Kandidaten nicht: jede E^nat-erhaltende Komposition besteht ihn, und das Gatter ist im Zustand "
      "unsichtbar (z121)",
      all(v < 1e-12 for v in kms.values()) and kms_theta < 1e-12 and kms_kontr > 1e-3 and produkt_kms < 1e-8)

say("  (d) Native sl(2): L^nat_(+-1), E^nat")


def ad_phase(G, X):
    """G X G^+ = c X auf den Spalten mit E <= E_max - 1: (c, relativer Rest)."""
    Y = (G @ X @ G.conj().T)[:, TIEF]
    Xt = X[:, TIEF]
    c_ = np.vdot(Xt.ravel(), Y.ravel())/np.vdot(Xt.ravel(), Xt.ravel())
    return complex(np.round(c_, 12)), nrm(Y - c_*Xt)/nrm(Xt)


sl2_d = {}
for nm, G in GATTER.items():
    cp, rp = ad_phase(G, LP)
    cm, rm = ad_phase(G, LM)
    sl2_d[nm] = {"[U, E^nat]": nrm(komm(G, ENAT)), "Ad L_+ = c L_+: c": str(cp), "Rest +": rp, "Ad L_- : c": str(cm),
                 "Rest -": rm}
fam_b = {}
for b in range(5):
    cp, rp = ad_phase(gatter(HQ + b*L0OSC), LP)
    fam_b[str(b)] = {"c": str(cp), "Rest": rp}
ZAHLEN["sl2"] = {"Gatter": sl2_d, "Familie e^(-i pi (H_Q + b L_0^osc)/2)": fam_b}
streng = [nm for nm, d in sl2_d.items() if d["Rest +"] < 1e-10 and d["Rest -"] < 1e-10
          and abs(complex(d["Ad L_+ = c L_+: c"]) - 1) < 1e-10
          and d["[U, E^nat]"] < 1e-10]
kovariant = [nm for nm, d in sl2_d.items() if d["Rest +"] < 1e-10 and d["Rest -"] < 1e-10 and d["[U, E^nat]"] < 1e-10]
check(f"[E] Native sl(2) (L^nat_(+-1), E^nat; Spalten E <= 3): alle vier Gatter erhalten E^nat. Streng mit L^nat_(+-1) "
      f"vertauscht nur {streng}; Ad_U L^nat_+ = c L^nat_+ mit c = {[sl2_d[k]['Ad L_+ = c L_+: c'] for k in GATTER]} für "
      f"(U_full, U_in, U_out, U_rein): U_full und U_out drehen die sl(2) um die Viertelphase -i (sie enthalten die freie "
      f"Konformzeit e^(-i pi L_0^osc/2)), U_rein zerstört sie (Rest {sl2_d['U_rein']['Rest +']:.3f}: der Ladungstausch bei festen "
      f"Oszillatoren verletzt die Kopplung Q_s J_(-1,s)). Familie e^(-i pi (H_Q + b L_0^osc)/2), b = 0..4: c = "
      f"{[v['c'] for v in fam_b.values()]}: strenge Vertauschung genau für b in 4Z, also genau wenn die Phase die Oszillatoren "
      "nicht sieht (e^(-2 pi i L_0^osc) = 1)",
      streng == ["U_in"] and set(kovariant) == {"U_full", "U_in", "U_out"} and sl2_d["U_rein"]["Rest +"] > 0.1
      and all((abs(complex(fam_b[str(b)]["c"]) - 1) < 1e-10) == (b % 4 == 0) for b in range(5))
      and abs(complex(sl2_d["U_full"]["Ad L_+ = c L_+: c"]) + 1j) < 1e-10
      and all(v["Rest"] < 1e-10 for v in fam_b.values()))

# Satz: rein innen (V (x) 1) und [V, L^nat_(+-1)] = 0 => [V, Q_A] = [V, Q_B] = 0
mix = {}
for t in np.linspace(0, np.pi, 7):
    Vt = np.cos(t)*I_ - 1j*np.sin(t)*S_IN
    mix[f"{t/np.pi:.4f}"] = round(nrm(komm(Vt, LP)[:, TIEF]), 6)
KLEIN4 = np.diag([1j**int(z[1][0][0]*z[1][1][0]) for z in ZU])
klein_sl2 = nrm(komm(KLEIN4, LP)) + nrm(komm(KLEIN4, ENAT))
klein_ein = nrm(KLEIN4 - expm(1j*np.pi*(QT @ QT - QA @ QA - QB @ QB)/4))
ZAHLEN["satz_rein_innen"] = {"V_t = cos t - i sin t S_innen: ||[V_t, L^nat_+]|| (t/pi)": mix,
                             "i^(Q_A Q_B): ||[., L^nat_+]|| + ||[., E^nat]||": klein_sl2,
                             "i^(Q_A Q_B) = e^(i pi (Q_tot^2 - Q_A^2 - Q_B^2)/4)": klein_ein}
check(f"[E] Satz (rein innen + native sl(2) => diagonal): L^nat_+ = sum_s (L^osc_(+,s) + Q_s J_(-1,s)) mit "
      f"Q-unabhängigen Außenoperatoren (Matrixelement = Q, Teil a); für V (x) 1 folgt [V (x) 1, L^nat_+] = sum_s [V, Q_s] (x) "
      f"J_(-1,s), also [V, Q_A] = [V, Q_B] = 0. Probe: die Wurzelfamilie des Ladungstauschs V_t = cos t - i sin t S_innen "
      f"vertauscht nur bei t in pi Z ({mix}). Diagonal bleibt etwa die Z4-Wechselphase i^(Q_A Q_B) (vertauscht mit L^nat, "
      f"E^nat: {klein_sl2:.0e}); sie ist e^(i pi Q_tot^2/4) mal Einzelnahtphasen ({klein_ein:.0e}), bei superselektierter "
      "Gesamtladung also einkörperig: Ein rein inneres, sl(2)-verträgliches Gatter wechselwirkt nicht",
      all((v < 1e-10) == (abs(float(k) - round(float(k))) < 1e-9) for k, v in mix.items()) and klein_sl2 < 1e-10
      and klein_ein < 1e-10)
befund("E: Trilemma der Zwei-Faktoren-Architektur. Faktorzerlegung (a), Uhr/Nahtparität (b), strenge sl(2)-Verträglichkeit (d) "
       "und Nichttrivialität sind nicht zugleich erfüllbar: U_full erfüllt a, b, c, e, (d kovariant), ist aber trivial; U_rein "
       "erfüllt a, verletzt aber Nahtparität und sl(2); U_in ist nichttrivial und erfüllt b (Paaruhr, Parität), c, d streng, "
       "e, verletzt aber a (Energieaustausch und Verschränkung Gitter <-> Oszillatoren). Die Architektur selbst (Teil a, "
       "Kontrolle) ist kein strenges Produkt: Der native Heber koppelt Q an J_(-1). P_null folgt also nicht aus der "
       "Faktorzerlegung; es ist genau die strenge Vertauschung mit der nativen sl(2)", True)

say("  (e) Lokalität je Nahtpaar (drei Nähte)")


def jw(n):
    Zp, sm = np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [0.0, 0.0]])
    out = []
    for j in range(n):
        op = np.array([[1.0]])
        for M_ in [Zp]*j + [sm] + [np.eye(2)]*(n - j - 1):
            op = np.kron(op, M_)
        out.append(op)
    return out


C9 = jw(9)                                        # drei Nähte x Moden y = -1, 0, 1 (x = -3/4, 1/4, 5/4)
N9 = [c.T @ c for c in C9]
Q9 = [sum(N9[3*s + k] for k in range(3)) - np.eye(512) for s in range(3)]   # See: y = -1 besetzt
E9 = [sum((y + 0.5)*N9[3*s + y + 1] for y in (0, 1)) + 0.5*(np.eye(512) - N9[3*s]) for s in range(3)]


def paar_gatter(i, j, H):
    Jp = sum(C9[3*i + k].T @ C9[3*j + k] for k in range(3))
    Bm = expm(-1j*np.pi/4*(Jp + Jp.T))
    return Bm.conj().T @ expm(-1j*np.pi/2*H) @ Bm


U_AB = paar_gatter(0, 1, (Q9[0] @ Q9[0] + Q9[1] @ Q9[1])/2)
U_AB_full = paar_gatter(0, 1, E9[0] + E9[1])
U_AB_kontr = paar_gatter(0, 1, (Q9[0] @ Q9[0] + Q9[1] @ Q9[1])/2 + 0.5*Q9[0] @ Q9[2])
C_OPS = [C9[6 + k] for k in range(3)] + [C9[6 + k].T @ C9[6 + l] for k in range(3) for l in range(3)]
lok = {nm: max(nrm(komm(G, X)) for X in C_OPS) for nm, G in (("U_in", U_AB), ("U_full", U_AB_full), ("Kontrolle Q_A Q_C", U_AB_kontr))}
U_BC = paar_gatter(1, 2, (Q9[1] @ Q9[1] + Q9[2] @ Q9[2])/2)
ordnung = nrm(komm(U_AB, U_BC))
U_BC_full = paar_gatter(1, 2, E9[1] + E9[2])
ordnung_full = nrm(komm(U_AB_full, U_BC_full))
ZAHLEN["lokalitaet"] = {"max ||[U_AB, Operator von C]||": lok, "||[U_AB, U_BC]|| (U_in)": ordnung,
                        "||[U_AB, U_BC]|| (U_full)": ordnung_full}
check(f"[E] Lokalität (drei Nähte, je drei Moden, Fock 512): U_AB vertauscht mit allen Operatoren der Naht C (c_C mit "
      f"Jordan-Wigner-Kette und alle Bilinearen; U_in {lok['U_in']:.0e}, U_full {lok['U_full']:.0e}); Kontrolle mit einem "
      f"Q_A Q_C-Term: {lok['Kontrolle Q_A Q_C']:.3f}. Benachbarte Innengatter vertauschen nicht ([U_AB, U_BC] = {ordnung:.3f}; "
      f"U_full: {ordnung_full:.0e}): Die Reihenfolge der Kompositionen zählt, die Kette braucht P_ord (z121)",
      lok["U_in"] < 1e-10 and lok["U_full"] < 1e-10 and lok["Kontrolle Q_A Q_C"] > 0.1 and ordnung > 0.1
      and ordnung_full < 1e-10)

# =================================================================== Teil 3
say("\nTeil 3: P_viertel -- ist phi = pi/2 (u = 1) nativ festgelegt?")
N_ABS = QT + 2*M*I_                                         # absolute Teilchenzahl im Fenster (See: 2M Moden besetzt)
GAM_X = expm(-1j*np.pi/2*N_ABS) @ expm(1j*np.pi*JX)         # kanonischer Tausch z106/z121 (Transferphase reell)
GAM_Y = expm(-1j*np.pi/2*N_ABS) @ expm(1j*np.pi*JY)         # derselbe Tausch im uhrgedrehten Rahmen
SEK_G = [(qt, e) for qt, e in SEKTOREN if qt % 2 == 0]
SEK_U = [(qt, e) for qt, e in SEKTOREN if qt % 2 == 1]
PHIS = [k*np.pi/12 for k in range(1, 24)]


def je_sektor_max(A, B, sek):
    return max(bis_auf_phase(A[np.ix_(sektor(*s), sektor(*s))], B[np.ix_(sektor(*s), sektor(*s))])
               / max(nrm(B[np.ix_(sektor(*s), sektor(*s))]), 1e-300) for s in sek)


def ist_skalar(A, sek):
    return max(bis_auf_phase(A[np.ix_(sektor(*s), sektor(*s))], np.eye(len(sektor(*s)))) for s in sek)


wurzel_y, wurzel_x, wurzel_u, ordnung, selbstdual, uhrpotenz, beta_sel = [], [], [], {}, [], [], []
RHOS = {b: np.diag(np.exp(-b*dE)) for b in (4*np.pi, 2*np.pi, 1.0)}
CHARS = [np.diag(np.exp(1j*kap*np.diag(QA - QB))) for kap in np.arange(8)*np.pi/4]
for ph in PHIS:
    U2 = gatter(HQ, phi=ph) @ gatter(HQ, phi=ph)
    lab = f"{ph/np.pi:.4f}"
    if je_sektor_max(U2, GAM_Y, SEK_G) < 1e-10:
        wurzel_y.append(lab)
    if je_sektor_max(U2, GAM_X, SEK_G) < 1e-10:
        wurzel_x.append(lab)
    if je_sektor_max(U2, GAM_Y, SEK_U) < 1e-10 or je_sektor_max(U2, GAM_X, SEK_U) < 1e-10:
        wurzel_u.append(lab)
    Uk, k_min = I_.copy(), None
    U1 = gatter(HQ, phi=ph)
    for k in range(1, 25):
        Uk = Uk @ U1
        if nrm(komm(Uk, QA)) < 1e-9:
            k_min = k
            break
    ordnung[lab] = k_min
    if ist_skalar(gatter(HQ, phi=ph) @ gatter(HQ, phi=np.pi - ph).conj().T, SEKTOREN) < 1e-10:
        selbstdual.append(lab)
    P = expm(-1j*ph*HQ)
    if all(any(bis_auf_phase(P[np.ix_(s_, s_)], C[np.ix_(s_, s_)]) < 1e-10 for C in CHARS)
           for s_ in [[k for k in range(DIM) if dQT[k] == q] for q in sorted(set(dQT.tolist()))]):
        uhrpotenz.append(lab)
    if all(nrm(komm(U1, R_)) < 1e-10 for R_ in RHOS.values()):
        beta_sel.append(lab)
diag_x = je_sektor_max(U_IN @ U_IN, GAM_X, [(2, 1.0), (2, 2.0), (2, 3.0)])
DIAG11 = [k for k in DIAG if dQT[k] == 2]
auf_11 = nrm((U_IN @ U_IN + GAM_X)[np.ix_(DIAG11, DIAG11)])
skalar_u = ist_skalar(U_IN @ U_IN, SEK_U)
ZAHLEN["p_viertel"] = {"phi/pi mit U^2 ~ Gamma_y je geradem Sektor": wurzel_y,
                       "phi/pi mit U^2 ~ Gamma_x je geradem Sektor": wurzel_x,
                       "phi/pi mit U^2 ~ Gamma auf ungeraden Sektoren": wurzel_u,
                       "U(pi/2)^2 + Gamma_x auf Q_A = Q_B = 1 (alle Oszillatoren)": auf_11,
                       "U(pi/2)^2 je ungeradem Sektor skalar (Abstand)": skalar_u,
                       "kleinste Ordnung k mit U(phi)^k diagonal (k <= 24)": ordnung,
                       "selbstdual U(phi) ~ U(pi - phi) je Sektor": selbstdual,
                       "e^(-i phi H_Q) = Uhrcharakter i^(kappa (Q_A - Q_B)) je Q_tot": uhrpotenz,
                       "[U(phi), rho_beta] = 0 für beta = 4 pi, 2 pi, 1": len(beta_sel)}
vier = [k for k, v in ordnung.items() if v == 4]
check(f"[E] Kandidaten für P_viertel auf dem ganzen nativen Raum, Phasen phi = k pi/12 (k = 1..23): (1) Modulperiode: "
      f"U(phi) lässt den KMS-Zustand für beta = 4 pi, 2 pi, 1 bei allen {len(beta_sel)} Phasen fest: beta wählt phi nicht. (2) "
      f"Uhr-Viertelschritt: e^(-i phi H_Q) ist je Gesamtladung ein Uhrcharakter i^(kappa (Q_A - Q_B)) nur für phi/pi = "
      f"{uhrpotenz} (der Tausch u = oo; Charakter je Sektor frei), nie für pi/2: die Gatterphase ist quadratisch (Gaußsumme), die Uhr linear; z120s "
      f"e^(-i pi H_Q/2) = -i R^(-1/2) ist eine Aussage auf N = 2 (z120). (3) Uhrordnung: kleinste k mit U(phi)^k diagonal "
      f"{ordnung}; Ordnung 4 = Ordnung der Uhr genau für phi/pi = {vier}. (4) Selbstdualität U(phi) ~ U(pi - phi) je Sektor: "
      f"{selbstdual}. (5) Wurzel des Tauschs: U^2 ~ Gamma je geradem Sektor für phi/pi = {wurzel_y} (Tausch im uhrgedrehten "
      f"Rahmen Gamma_y), mit dem reellen Gamma_x von z121 nur {wurzel_x}; auf Q_A = Q_B = 1 mit allen Oszillatoren gilt "
      f"U^2 = -Gamma_x ({auf_11:.0e}, z121). Auf ungerader Gesamtladung ist U(pi/2)^2 je Sektor skalar ({skalar_u:.0e}) und "
      f"keine Phase macht U^2 zum Tausch ({wurzel_u})",
      len(beta_sel) == len(PHIS) and uhrpotenz == ["1.0000"] and vier == ["0.5000", "1.5000"]
      and selbstdual == ["0.5000", "1.5000"] and wurzel_y == ["0.5000", "1.5000"] and wurzel_x == []
      and auf_11 < 1e-10 and skalar_u < 1e-10 and wurzel_u == [] and diag_x > 1e-3)
befund("C: P_viertel ist nicht nativ festgelegt. Modulperiode und Uhr-Viertelschritt wählen phi nicht; Uhrordnung 4, "
       "Selbstdualität und 'U^2 = Tausch' wählen alle dasselbe Paar phi = pi/2, 3 pi/2 (u = +-1, die zwei Zeitrichtungen), "
       "sind aber gleichwertige Neuformulierungen derselben Prämisse, keine Herleitung. Die Wurzelform gilt nur bei gerader "
       "Gesamtladung, und auf dem ganzen geraden Sektor nur mit dem Tausch im uhrgedrehten Rahmen (Gamma_y); das reelle "
       "Gamma_x von z121 trifft nur Q_A = Q_B. Bei ungerader Gesamtladung ist U^2 eine Sektorphase, keine Wurzel eines Tauschs", True)

# =================================================================== Teil 4
say("\nTeil 4: Bilanz")
BILANZ = [
    ("E", "Zwei-Faktoren-Fenster exakt nativ (242 Zustände, E^nat <= 4): innen Q (Q^2/2, Uhr i^Q), außen lambda (L_0^osc); "
          "U_full = e^(-i pi E^nat/2) trivial; U_in U_out = U_full, Innen- und Außengatter zueinander invers"),
    ("E", "(a) U_in verletzt die Faktorzerlegung: Energieaustausch Gitter <-> Oszillatoren bei ungerader Gesamtladung "
          "(Transferspin >= 3/2), Innen-außen-Verschränkung bei gerader (Wurzel des vollen Nahttauschs), auf Q_A = Q_B "
          "1 (x) sqrt(SWAP) auf den Oszillatoren; Grund: Transfer-Nullmode = Vertexoperator (Frenkel-Kac)"),
    ("E", "(a) Kontrolle: der native Heber L^nat_+ = L^osc_+ + Q J_(-1) koppelt selbst Gitter und Oszillatoren; das strenge "
          "Produkt ist keine native Forderung"),
    ("E", "(b) U_full, U_in, U_out: global-, Paaruhr- und nahtparitätsinvariant; Nahtuhr R_A dreht den Transferrahmen "
          "(U_in nur auf Q_A = Q_B invariant); U_rein verletzt Paaruhr und Nahtparität; J = K: alle J U J = U^-1, Rahmen "
          "chi in (pi/2) Z"),
    ("E", "(c) KMS beta = 4 pi: alle E^nat-erhaltenden Gatter invariant; trennt nicht"),
    ("E", "(d) strenge sl(2)-Vertauschung nur U_in; U_full, U_out drehen L^nat_(+-1) um -i; U_rein bricht sl(2); Familie "
          "e^(-i pi (H_Q + b L_0^osc)/2) vertauscht genau für b in 4Z"),
    ("E", "Satz: rein innen + sl(2) => diagonal in (Q_A, Q_B) => bei superselektierter Gesamtladung einkörperig"),
    ("E", "(e) Lokalität je Nahtpaar; [U_AB, U_BC] != 0: P_ord bleibt nötig"),
    ("C", "P_null ist innerhalb der Familie gleichwertig mit P_sl2 'die Komposition vertauscht streng mit der nativen "
          "sl(2) (enthält keine freie Konformzeit)'; aus der Faktorzerlegung folgt es nicht"),
    ("C", "P_viertel: Uhrordnung 4, Selbstdualität, U^2 = Tausch (gerade Gesamtladung, Gamma_y) wählen gleichwertig "
          "phi = pi/2, 3 pi/2; Modulperiode und Uhrschritt wählen nichts; ungerade Gesamtladung: U^2 Sektorphase"),
    ("O", "eine native Begründung von P_sl2 (etwa: Komposition als Morphismus der radialen Darstellung) und von P_viertel"),
]
for st, txt in BILANZ:
    say(f"   [{st}] {txt}")
ZAHLEN["bilanz"] = [{"status": st, "aussage": txt} for st, txt in BILANZ]

urteil = ("P_null folgt nicht aus der Zwei-Faktoren-Architektur von z82, und das Innengatter ist nicht das einzige, das sie "
          "respektiert; die Erwartung 'nur das Innengatter besteht' ist widerlegt. Im exakt nativen Fockfenster zweier Nähte "
          "(242 Zustände) sind Innen- und Außengatter zueinander invers bis auf die triviale freie Phase (U_in U_out = U_full "
          "= e^(-i pi E^nat/2)). Die Faktorzerlegung respektieren nur das triviale U_full und das rein innere U_rein (Wurzel "
          "des Ladungstauschs bei festen Oszillatoren); U_rein verletzt aber Paaruhr, Nahtparität (P_univ) und die native "
          "sl(2), und allgemein gilt: rein innen und sl(2)-verträglich heißt diagonal im Ladungsgitter, also bei "
          "superselektierter Gesamtladung wechselwirkungsfrei. Das z120-Innengatter U_in verletzt die Zerlegung: Weil die "
          "Transfer-Nullmode ein Vertexoperator ist (das Cartan-Mitglied des Tripletts ist eine Oszillatoranregung), "
          "tauscht es bei ungerader Gesamtladung Energie zwischen Gitter und Oszillatoren aus, verschränkt bei gerader "
          "beide Faktoren (Wurzel des vollen Nahttauschs) und wirkt im physikalischen Sektor Q_A = Q_B als 1 (x) sqrt(SWAP) "
          "auf den Oszillatoren: Das Gitter liefert nur die Phase der virtuellen Zwischenzustände, die Wirkung liegt außen. "
          "Was U_in auszeichnet, ist die strenge Vertauschung mit der nativen sl(2) (L^nat_(+-1), E^nat): Die Phasenfamilie "
          "e^(-i pi (H_Q + b L_0^osc)/2) vertauscht genau für b in 4Z, U_full und U_out drehen die sl(2) um -i. P_null ist "
          "damit gleichwertig mit der Prämisse P_sl2 ('die Komposition enthält keine freie Konformzeit'), bedingt, nicht "
          "hergeleitet. Uhr (global, Paaruhr, Parität), Tomita-J (J U J = U^-1), KMS bei beta = 4 pi und Lokalität je "
          "Nahtpaar bestehen U_full, U_in und U_out gleichermaßen; die Nahtuhr R_A dreht den Transferrahmen des Innengatters "
          "(invariant nur auf Q_A = Q_B). P_viertel ist nicht nativ: Modulperiode und Uhr-Viertelschritt wählen phi nicht, "
          "Uhrordnung 4, Selbstdualität und 'U^2 = Tausch' wählen gleichwertig phi = pi/2, 3 pi/2 (u = +-1, Zeitrichtung). "
          "Die Wurzelform U^2 = -Gamma(pi) gilt nur bei gerader Gesamtladung (auf dem ganzen geraden Sektor mit dem Tausch "
          "im uhrgedrehten Rahmen, mit dem reellen Gamma(pi) nur auf Q_A = Q_B); bei ungerader ist U^2 eine Sektorphase. "
          "P_ord bleibt nötig ([U_AB, U_BC] != 0, Lokalität je Paar exakt). Für K15: Der Szenariowert 'lokales kohärentes "
          "Transfergatter' bleibt bedingt, mit den Prämissen P_sl2 (= P_null in der Familie), P_viertel und P_ord, gültig im "
          "Paarsektor gerader Gesamtladung, und er ist als Wechselwirkung der Außenfaktoren zu lesen (sqrt(SWAP) der "
          "Oszillatoren zweier Nähte, ausgelöst durch die Gitterenergie), nicht als rein innere Regel. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z133_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z133.txt").write_text("\n".join(OUT) + "\n")
