"""Z135: Yang-Mills aus den Schleifen des Transfernetzes? Plaketten, Normierung beta_lat, Kontinuumskopplung und Dynamik
(dritter Weg zu K15 neben der Twistorroute z110/z124 und der kohärenten Komposition z121; Folgeschritt zu z120 Teil 5).

Leitlinie 'Rekonstruktion statt Auswahl': Die Gatter kommen aus z120: die Vierteldrehung B = e^(-i pi J_x/2) (auf einer
Anregung der Transfer mit Mischwinkel pi/4, z120 Teil 4) und das Austauschgatter U = B^+ e^(-i pi H_Q/2) B (auf einer
Anregung je Naht P_T - i P_S, z120/z121). Die Eichgruppe wirkt auf den Nahtetiketten (U(1)_Y, SU(2)_L, SU(3) auf der 16;
SU(2) auf dem Dublett; zufälliges U(d) als generischer Fall); Linkvariablen sind Paralleltransporte in den Transferströmen
J^+ = sum c+_A U c_B. Die Schleifenantwort W(g) = (d^2 - |Tr g|^2)/(4 d^4) der achtzehnten Analyse ist Vergleichsziel,
keine Prämisse. Firewall DIMENSION.UPLIFT.FIREWALL.01: 2D- und 3D-Anordnungen der Nähte sind Konstruktionen, Kontinuum
und 4D nur Lesart. H6 und die gelösten Blöcke (alpha^-1 = 137.03599921684, die 16, ihre Hyperladungen, E8) bleiben fest
und dienen nur als Vergleich; nichts wird angepasst, nichts befördert. Statuskürzel: E exakt, C bedingt (Prämisse
genannt), L Lesart, N numerisch, O offen.

Teil 1  Plaketten aus dem Gatternetz. (a) Eichkovarianz der z120-Gatter im Fockraum, Einzelbindung = reine Eichung,
        Einschränkungen auf eine Anregung (B) und auf eine Anregung je Naht (-i R(1) mit Etikettentransport).
        (b) 2D-Netz der Vierteldrehungen in vier Teilschichten (Ordnung +x, +y, -x, -y wie das RLBL-Floquetnetz):
        Rückkehrblock cos^4 + sin^4 U_p, jede Plakette genau einmal je Periode; Spurantwort = fundamentale Wilson-Funktion.
        (c) Elementarer Austauschring (vier Nähte): Spurantwort = W(g) exakt, |Tr g|^2 aus zwei gegenläufigen Zyklen.
        (d) W = adjungierte Wilson-Wirkung (Darstellung V x V-quer): zentrumsblind, quadratisch gleich der fundamentalen.
        (e) Gitter: Austauschnetz nicht additiv (Schleifengas), Einzelanregungsnetz additiv. (f) 3D-Anordnung mit
        Ebenenteilperioden: bei theta = pi/4 längere Schleifen gleichen Gewichts. (g) Orientierungsungerader Teil auf der
        16: Tr X = Tr X^3 = 0 für so(10), erste Ordnung X^5 = Pfaffsche. (h) 1D-Kette: keine Plakette.
Teil 2  Kopplung: naiver Kontinuumslimes 1/g^2 = kappa_s T_R xi = kappa_t T_R/xi (anisotrop numerisch geprüft); Spuren
        über die 16 exakt: 1/g_i^2 = lambda k_i mit k_i = 2 Tr_16 T_i^2 = 4 für SU(3), SU(2), U(1)_Y (GUT), also
        sin^2 theta_W = 3/8 in jeder Lesart; Kontrollen (Träger-5, unvollständige Etiketten); drei Lesarten von lambda
        (Netzspur, W, Niveau mit P_norm) und Vergleich mit alpha und dem z25-Fenster, ohne Anpassung.
Teil 3  Dynamik: Z4-Quantenlinks (Uhr als Linkgruppe): Netz statisch, Gauß-Gesetz, Besetzungen eingefroren;
        Reflexionspositivität über die Charakterentwicklung; Vorzeichen der Antwort je Lesart; Lichtkegel des 2D-Netzes;
        was für echte Yang-Mills-Dynamik fehlt.
Teil 4  Bilanz, Urteil, Folgen für K15 und die Szenarioachse Wechselwirkung.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import itertools  # noqa: E402
import json  # noqa: E402
from fractions import Fraction as Fr  # noqa: E402
from functools import reduce  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.linalg import expm  # noqa: E402
from scipy.special import iv  # noqa: E402

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(135)
ALPHA_INV = 137.03599921684                      # gelöster Block, nur Vergleich


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def e0(x):
    return f"{float(x):.0e}"


def jw(n):
    I_, Zp, sm = np.eye(2), np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [0.0, 0.0]])
    out = []
    for j in range(n):
        m_ = [Zp]*j + [sm] + [I_]*(n - j - 1)
        op = m_[0]
        for M in m_[1:]:
            op = np.kron(op, M)
        out.append(op)
    return out


def haar(d):
    Z = (rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d)))/np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q*(np.diag(R)/abs(np.diag(R)))


def herm(d, scale=1.0):
    H = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
    return scale*(H + H.conj().T)/2


def W(g):
    d = g.shape[0]
    return float((d**2 - abs(np.trace(g))**2)/(4*d**4))


def floquet(nsite, d, layers, links, th, signs=None):
    """Einteilchen-Floquetoperator: Schichten aus Gattern (ia, ib, key); psi_A <- c psi_A - i s U psi_B,
    psi_B <- c psi_B - i s U^+ psi_A (U transportiert von B nach A)."""
    c, s = np.cos(th), np.sin(th)
    D = nsite*d
    F = np.eye(D, dtype=complex)
    for li, lay in enumerate(layers):
        sg = 1 if signs is None else signs[li]
        G = np.eye(D, dtype=complex)
        for ia, ib, key in lay:
            U = links[key]
            sa, sb = slice(ia*d, ia*d + d), slice(ib*d, ib*d + d)
            G[sa, sa] = c*np.eye(d)
            G[sb, sb] = c*np.eye(d)
            G[sa, sb] = -sg*1j*s*U
            G[sb, sa] = -sg*1j*s*U.conj().T
        F = G @ F
    return F


def huepfpfad(x, layers, links, d):
    """Pfad, der in jeder Schicht springt: Endpunkt, Transport, besuchte Orte."""
    pos, M, besucht = x, np.eye(d, dtype=complex), [x]
    for lay in layers:
        for ia, ib, key in lay:
            if pos == ia:
                M, pos = links[key].conj().T @ M, ib
                break
            if pos == ib:
                M, pos = links[key] @ M, ia
                break
        besucht.append(pos)
    return pos, M, besucht


def perm_trace(n, d, gates, links, a_st, a_sw, nur_voll=False):
    """Tr ueber (C^d)^(x n) von prod (a_st 1 + a_sw P_e(U)), P_e(U)|a>_p|b>_q = |U b>_p |U^+ a>_q; exakte
    Permutationsentwicklung (jede Teilmenge der Täusche ist eine Permutation mit Transporten)."""
    total = 0j
    I_ = np.eye(d, dtype=complex)
    for S in itertools.product((0, 1), repeat=len(gates)):
        if nur_voll and sum(S) < len(gates):
            continue
        pos, M, at, w = list(range(n)), [I_]*n, list(range(n)), 1.0 + 0j
        for (p, q, key), sw in zip(gates, S):
            if sw:
                U = links[key]
                lp, lq = at[p], at[q]
                M[lp], M[lq] = U.conj().T @ M[lp], U @ M[lq]
                pos[lp], pos[lq] = q, p
                at[p], at[q] = lq, lp
                w *= a_sw
            else:
                w *= a_st
        seen, tr = [False]*n, 1.0 + 0j
        for i in range(n):
            if seen[i]:
                continue
            j, P = i, I_
            while not seen[j]:
                seen[j] = True
                P = M[j] @ P
                j = pos[j]
            tr *= np.trace(P)
        total += w*tr
    return total


def zyklen_voll(n, gates):
    """Zyklen der Permutation 'alle Gatter tauschen' mit den Orten, die jedes Etikett durchläuft."""
    pos, at, wege = list(range(n)), list(range(n)), [[i] for i in range(n)]
    for p, q, _ in gates:
        lp, lq = at[p], at[q]
        pos[lp], pos[lq] = q, p
        at[p], at[q] = lq, lp
        wege[lp].append(q)
        wege[lq].append(p)
    seen, zyk = [False]*n, []
    for i in range(n):
        if seen[i]:
            continue
        j, z, weg = i, [], [i]
        while not seen[j]:
            seen[j] = True
            z.append(j)
            weg += wege[j][1:]
            j = pos[j]
        zyk.append((z, weg))
    return zyk


def dense_exchange(n, d, gates, links, a_st, a_sw):
    D = d**n
    F = np.eye(D, dtype=complex)
    for p, q, key in gates:
        U = links[key]
        T = np.zeros((D, D), dtype=complex)
        for idx_ in itertools.product(range(d), repeat=n):
            col = np.ravel_multi_index(idx_, (d,)*n)
            a, b = idx_[p], idx_[q]
            for ap in range(d):
                for bq in range(d):
                    amp = U[ap, b]*np.conj(U[a, bq])
                    if abs(amp) < 1e-300:
                        continue
                    new = list(idx_)
                    new[p], new[q] = ap, bq
                    T[np.ravel_multi_index(new, (d,)*n), col] += amp
        F = (a_st*np.eye(D) + a_sw*T) @ F
    return F


def pfaff(A):
    A = A.astype(complex).copy()
    n, res = A.shape[0], 1.0 + 0j
    for k in range(0, n - 1, 2):
        kp = k + 1 + int(np.argmax(abs(A[k, k + 1:])))
        if kp != k + 1:
            A[[k + 1, kp], :] = A[[kp, k + 1], :]
            A[:, [k + 1, kp]] = A[:, [kp, k + 1]]
            res *= -1
        if A[k, k + 1] == 0:
            return 0.0
        res *= A[k, k + 1]
        if k + 2 < n:
            tau = A[k, k + 2:]/A[k, k + 1]
            A[k + 2:, k + 2:] += np.outer(A[k + 1, k + 2:], tau) - np.outer(tau, A[k + 1, k + 2:])
    return res


TH = np.pi/4                                       # Vierteldrehung: Mischwinkel pi/4 (z120, P_sym)
CC, SS = np.cos(TH), np.sin(TH)
A_ST, A_SW = 1/(1 + 1j), 1j/(1 + 1j)               # R(1) = (1 + i P)/(1 + i)

# ------------------------------------------------------------------ die 16 (Lambda^gerade C^5, z82)
c5 = jw(5)
PAR = np.array([bin(i).count("1") % 2 for i in range(32)])
EV, OD = np.where(PAR == 0)[0], np.where(PAR == 1)[0]


def on16(X):
    return sum(c5[i].T @ (X[i, j]*c5[j]) for i in range(5) for j in range(5) if X[i, j] != 0)[np.ix_(EV, EV)]


Y5 = np.diag([-1/3]*3 + [1/2]*2)
Y16 = on16(Y5).real
GM = np.zeros((8, 3, 3), complex)
GM[0][0, 1] = GM[0][1, 0] = 1
GM[1][0, 1], GM[1][1, 0] = -1j, 1j
GM[2][0, 0], GM[2][1, 1] = 1, -1
GM[3][0, 2] = GM[3][2, 0] = 1
GM[4][0, 2], GM[4][2, 0] = -1j, 1j
GM[5][1, 2] = GM[5][2, 1] = 1
GM[6][1, 2], GM[6][2, 1] = -1j, 1j
GM[7] = np.diag([1, 1, -2])/np.sqrt(3)
PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0]).astype(complex)]
T3_16, T2_16 = [], []
for a in range(8):
    X_ = np.zeros((5, 5), complex)
    X_[:3, :3] = GM[a]/2
    T3_16.append(on16(X_))
for a in range(3):
    X_ = np.zeros((5, 5), complex)
    X_[3:, 3:] = PAULI[a]/2
    T2_16.append(on16(X_))


def su2L(theta):
    return expm(1j*sum(t*T for t, T in zip(theta, T2_16)))


def u1Y(phi):
    return np.diag(np.exp(1j*phi*np.diag(Y16)))


def su2():
    g = haar(2)
    return g/np.sqrt(np.linalg.det(g))


# =================================================================== Teil 1
say("Teil 1: Plaketten aus dem Gatternetz")
say("  (a) Eichkovarianz der z120-Gatter, Einzelbindung, Einschränkungen")
c4 = jw(4)
cA, cB = c4[:2], c4[2:]
n4 = [x.T @ x for x in c4]
QA, QB = n4[0] + n4[1], n4[2] + n4[3]
HQ = (QA @ QA + QB @ QB)/2


def jplus(U):
    return sum(cA[m].T @ (U[m, mp]*cB[mp]) for m in range(2) for mp in range(2))


def gatter(U):
    J = jplus(U)
    Bm = expm(-1j*np.pi*(J + J.conj().T)/4)
    return Bm.conj().T @ expm(-1j*np.pi*HQ/2) @ Bm


def gamma2(H, modes):
    return expm(1j*sum(modes[m].T @ (H[m, mp]*modes[mp]) for m in range(2) for mp in range(2)))


HU, HhA, HhB = herm(2), herm(2), herm(2)
Ulink, hA, hB = expm(1j*HU), expm(1j*HhA), expm(1j*HhB)
GAB = gamma2(HhA, cA) @ gamma2(HhB, cB)
kov = np.linalg.norm(GAB @ gatter(Ulink) @ GAB.conj().T - gatter(hA @ Ulink @ hB.conj().T))
kov_falsch = np.linalg.norm(GAB @ gatter(Ulink) @ GAB.conj().T - gatter(hB @ Ulink @ hA.conj().T))
GBd = gamma2(-HU, cB)
rein = np.linalg.norm(gatter(Ulink) - GBd @ gatter(np.eye(2)) @ GBd.conj().T)
vac16 = np.zeros(16)
vac16[0] = 1
BAS = np.array([cA[m].T @ cB[mp].T @ vac16 for m in range(2) for mp in range(2)]).T


def Pswap(U):
    d = U.shape[0]
    T = np.zeros((d*d, d*d), dtype=complex)
    for a in range(d):
        for b in range(d):
            T[:, a*d + b] = np.kron(U[:, b], U.conj().T[:, a])
    return T


G4 = BAS.conj().T @ gatter(Ulink) @ BAS
R1U = (np.eye(4) + 1j*Pswap(Ulink))/(1 + 1j)
phys_abw = np.linalg.norm(G4 - (-1j)*R1U)
phys_falsch = np.linalg.norm(G4 - (-1j)*(np.eye(4) + 1j*Pswap(Ulink.conj().T))/(1 + 1j))
ein = [cA[m].T @ vac16 for m in range(2)] + [cB[m].T @ vac16 for m in range(2)]
E1 = np.array(ein).T
G1 = E1.conj().T @ gatter(Ulink) @ E1
ein_abw = np.linalg.norm(G1 - np.exp(-1j*np.pi/4)*np.eye(4))
ZAHLEN["fock_gatter"] = {"Kovarianz": float(kov), "Kontrolle falsch transformiert": float(kov_falsch),
                         "Einzelbindung reine Eichung": float(rein), "(1,1): Abstand zu -i R(1)(U)": float(phys_abw),
                         "(1,1): Kontrolle U^+ statt U": float(phys_falsch), "eine Anregung: Abstand zu e^(-i pi/4)": float(ein_abw)}
check(f"[E] Link U im Transferstrom J^+ = sum c+_(A,m) U_(m m') c_(B,m') (Fockraum zweier Nähte, z120): Gamma(h_A + h_B) "
      f"U(U) Gamma^+ = U(h_A U h_B^+) (Abweichung {e0(kov)}; mit h_B U h_A^+ verfehlt um {kov_falsch:.3f}); eine "
      f"Einzelbindung ist reine Eichung, U(U) = Gamma_B(U^+) U(1) Gamma_B(U^+)^+ ({e0(rein)}). Auf dem physikalischen Raum "
      f"(eine Anregung je Naht) ist das Gatter -i R(1)(U) = -i (1 + i P(U))/(1 + i) mit dem Etikettentausch samt Transport "
      f"P(U)|a>_A|b>_B = |U b>_A |U^+ a>_B ({e0(phys_abw)}; mit U^+ statt U {phys_falsch:.3f}); auf einer Anregung allein "
      f"bleibt es auch mit Link die reine Phase e^(-i pi/4) ({e0(ein_abw)}). Der Einteilchentransfer ist also die "
      f"Vierteldrehung B selbst (Mischwinkel pi/4, z120 Teil 4), der Etikettenaustausch das z120/z121-Gatter",
      kov < 1e-12 and kov_falsch > 0.1 and rein < 1e-12 and phys_abw < 1e-12 and phys_falsch > 0.1 and ein_abw < 1e-12)

say("  (b) 2D-Netz der Vierteldrehungen: Plaketten als Elementarzyklen")
BS2 = [(1, 0), (0, 1), (-1, 0), (0, -1)]
L2 = 4
SITES2 = [(x, y) for x in range(L2) for y in range(L2)]
IDX2 = {s_: i for i, s_ in enumerate(SITES2)}
A2 = [s_ for s_ in SITES2 if (s_[0] + s_[1]) % 2 == 0]
BONDS2 = {}
for j, b in enumerate(BS2):
    for a in A2:
        BONDS2[(j, a)] = (IDX2[a], IDX2[((a[0] + b[0]) % L2, (a[1] + b[1]) % L2)])


def schichten2(order=(0, 1, 2, 3)):
    return [[(BONDS2[(j, a)][0], BONDS2[(j, a)][1], (j, a)) for a in A2] for j in order]


LAY2 = schichten2()
NS2 = len(SITES2)


def torus_analyse(links, d):
    F = floquet(NS2, d, LAY2, links, TH)
    abw, plak, flaechen, spur_hol = 0.0, [], [], 0j
    for x in range(NS2):
        end, M, bes = huepfpfad(x, LAY2, links, d)
        blk = F[x*d:(x + 1)*d, x*d:(x + 1)*d]
        abw = max(abw, np.linalg.norm(blk - (CC**4*np.eye(d) + SS**4*M)), 0.0 if end == x else 1.0)
        spur_hol += np.trace(M)
        plak.append(tuple(sorted(set(bes))))
        # orientierte Fläche aus den Sprungvektoren (entfaltet)
        pts, cur = [(0, 0)], (0, 0)
        for u, v in zip(bes[:-1], bes[1:]):
            du = ((SITES2[v][0] - SITES2[u][0] + 1) % L2) - 1, ((SITES2[v][1] - SITES2[u][1] + 1) % L2) - 1
            cur = (cur[0] + du[0], cur[1] + du[1])
            pts.append(cur)
        flaechen.append(0.5*sum(pts[i][0]*pts[i + 1][1] - pts[i + 1][0]*pts[i][1] for i in range(len(pts) - 1)))
    spur_abw = abs(np.trace(F) - (NS2*d*CC**4 + SS**4*spur_hol))
    return F, abw, len(set(plak)), sorted(set(flaechen)), spur_abw


ergebnisse = {}
fall = {"U(16) zufällig": (16, lambda: haar(16)), "SU(2)_L c 16": (16, lambda: su2L(rng.normal(size=3))),
        "U(1)_Y c 16": (16, lambda: u1Y(rng.normal())), "SU(2) auf dem Dublett": (2, lambda: su2())}
for name, (d, gen) in fall.items():
    links = {k: gen() for k in BONDS2}
    F, abw, nplak, fl, sabw = torus_analyse(links, d)
    unit = np.linalg.norm(F.conj().T @ F - np.eye(NS2*d))
    ergebnisse[name] = {"Rückkehrblock": float(abw), "Plaketten": nplak, "orientierte Flächen": fl, "Spur": float(sabw),
                        "Unitarität": float(unit)}
d = 2
links = {k: haar(d) for k in BONDS2}
hx = [haar(d) for _ in range(NS2)]
links_h = {k: hx[BONDS2[k][0]] @ links[k] @ hx[BONDS2[k][1]].conj().T for k in BONDS2}
Hbl = np.zeros((NS2*d, NS2*d), complex)
for x in range(NS2):
    Hbl[x*d:(x + 1)*d, x*d:(x + 1)*d] = hx[x]
F_ = floquet(NS2, d, LAY2, links, TH)
kov2 = np.linalg.norm(floquet(NS2, d, LAY2, links_h, TH) - Hbl @ F_ @ Hbl.conj().T)
eins = {k: np.eye(d, dtype=complex) for k in BONDS2}
F1 = floquet(NS2, d, LAY2, eins, TH)
lay_k = schichten2((0, 2, 1, 3))
Fk, Fk1 = floquet(NS2, d, lay_k, links, TH), floquet(NS2, d, lay_k, eins, TH)
kontr_ord = abs(np.trace(Fk) - np.trace(Fk1))
kontr_blk = max(np.linalg.norm(Fk[x*d:(x + 1)*d, x*d:(x + 1)*d] - CC**4*np.eye(d)) for x in range(NS2))
links_pg = {k: hx[BONDS2[k][0]] @ hx[BONDS2[k][1]].conj().T for k in BONDS2}
kontr_pg = abs(np.trace(floquet(NS2, d, LAY2, links_pg, TH)) - np.trace(F1))
wilson = SS**4*sum(np.trace(np.eye(d) - huepfpfad(x, LAY2, links, d)[1]).real for x in range(NS2))
wil_abw = abs((np.trace(F1) - np.trace(F_)).real - wilson)
ZAHLEN["netz_2d"] = {"Fälle": ergebnisse, "Kovarianz F[U^h] = H F H^+": float(kov2),
                     "Kontrolle Schrittfolge (+x, -x, +y, -y): Spurantwort": float(kontr_ord),
                     "Kontrolle Schrittfolge: Rückkehrblock - cos^4": float(kontr_blk),
                     "Kontrolle reine Eichlinks: Spurantwort": float(kontr_pg),
                     "Re Tr(F[1] - F[U]) - sin^4 sum Re Tr(1 - U_p)": float(wil_abw), "kappa = sin^4(pi/4)": SS**4}
ok_faelle = all(v["Rückkehrblock"] < 1e-12 and v["Plaketten"] == NS2 and v["orientierte Flächen"] == [-1.0]
                and v["Spur"] < 1e-10 and v["Unitarität"] < 1e-10 for v in ergebnisse.values())
check(f"[E] 2D-Netz der Vierteldrehungen (Torus 4 x 4, bipartit, vier Teilschichten in der Folge +x, +y, -x, -y wie das "
      f"RLBL-Floquetnetz, Link in jedem Sprung): Der Rückkehrblock ist <x|F|x> = cos^4 theta 1 + sin^4 theta U_p(x) "
      f"(max Abweichung {e0(max(v['Rückkehrblock'] for v in ergebnisse.values()))} für {', '.join(ergebnisse)}), "
      f"der Pfad mit vier Sprüngen umläuft eine Plakette, jede der 16 Plaketten genau einmal je Periode und alle mit "
      f"derselben Orientierung (orientierte Flächen {ergebnisse['U(16) zufällig']['orientierte Flächen']}). Also Tr F[U] = "
      f"N d cos^4 theta + sin^4 theta sum_p Tr U_p (max {e0(max(v['Spur'] for v in ergebnisse.values()))}) und Re Tr(F[1] - "
      f"F[U]) = sin^4 theta sum_p Re Tr(1 - U_p) ({e0(wil_abw)}): genau die fundamentale Wilson-Funktion, additiv, mit "
      f"kappa = sin^4(pi/4) = 1/4 je Periode. Eichkovarianz F[U^h] = H F[U] H^+ ({e0(kov2)}). Kontrollen: Die Schrittfolge "
      f"+x, -x, +y, -y schließt keine Schleife (Spurantwort {e0(kontr_ord)}, Rückkehrblock = cos^4, {e0(kontr_blk)}); "
      f"reine Eichlinks h_x h_y^+ geben keine Antwort ({e0(kontr_pg)})",
      ok_faelle and kov2 < 1e-12 and kontr_ord < 1e-12 and kontr_blk < 1e-12 and kontr_pg < 1e-12 and wil_abw < 1e-10)

say("  (c) Elementarer Austauschring: die Schleifenantwort W(g) der achtzehnten Analyse")
RING = [(0, 1, "e01"), (1, 2, "e12"), (2, 3, "e23"), (3, 0, "e30")]


def ring_hol(links):
    return (links["e30"].conj().T @ links["e23"].conj().T @ links["e12"].conj().T @ links["e01"].conj().T)


ring_dicht = {}
for d in (2, 3):
    links = {k: haar(d) for _, _, k in RING}
    eins = {k: np.eye(d, dtype=complex) for k in links}
    Fd = dense_exchange(4, d, RING, links, A_ST, A_SW)
    tP, t1 = perm_trace(4, d, RING, links, A_ST, A_SW), perm_trace(4, d, RING, eins, A_ST, A_SW)
    g = ring_hol(links)
    ring_dicht[d] = {"dicht gegen Permutation": float(abs(np.trace(Fd) - tP)),
                     "(Tr F[U] - Tr F[1])/d^4 - W": float(abs((tP - t1)/d**4 - W(g))), "W": W(g),
                     "Tr F[1]": [float(t1.real), float(t1.imag)]}
ring16 = {}
eins16 = {k: np.eye(16, dtype=complex) for _, _, k in RING}
t1_16 = perm_trace(4, 16, RING, eins16, A_ST, A_SW)
for name, links in [("U(16) zufällig", {k: haar(16) for _, _, k in RING}),
                    ("SU(2)_L c 16", {k: su2L(rng.normal(size=3)) for _, _, k in RING})]:
    tP = perm_trace(4, 16, RING, links, A_ST, A_SW)
    ring16[name] = float(abs((tP - t1_16)/16**4 - W(ring_hol(links))))
phi_werte = (1e-3, 0.3, 1.0)
y_antwort = {}
for ph in phi_werte:
    links = {"e01": u1Y(ph), "e12": np.eye(16, dtype=complex), "e23": np.eye(16, dtype=complex),
             "e30": np.eye(16, dtype=complex)}
    tP = perm_trace(4, 16, RING, links, A_ST, A_SW)
    y_antwort[ph] = ((tP - t1_16)/16**4, W(u1Y(ph)))
y_klein = y_antwort[1e-3][0].real/1e-6
zyk = zyklen_voll(4, RING)


def windung(weg):
    """Windungszahl eines geschlossenen Orts-Wegs auf dem Ring 0-1-2-3."""
    w = 0
    for u, v in zip(weg[:-1], weg[1:]):
        w += {1: 1, 3: -1, 0: 0}[(v - u) % 4]
    return w//4


windungen = sorted(windung(weg) for _, weg in zyk)
RING_alt = [(0, 1, "e01"), (2, 3, "e23"), (1, 2, "e12"), (3, 0, "e30")]
links = {k: haar(2) for _, _, k in RING}
eins = {k: np.eye(2, dtype=complex) for k in links}
alt_abw = abs((perm_trace(4, 2, RING_alt, links, A_ST, A_SW) - perm_trace(4, 2, RING_alt, eins, A_ST, A_SW))/16
              - W(ring_hol(links)))
lay_ring = [[(p, q, k)] for p, q, k in RING]
bring = {}
for d in (2, 16):
    links = {k: haar(d) for _, _, k in RING}
    bring[d] = float(abs(np.trace(floquet(4, d, lay_ring, links, TH)) - (4*CC**2*d + SS**4*np.trace(ring_hol(links)))))
tf1_formel = {d: -(d**4 - 5*d**2 + 4j*(d**3 - d))/4 for d in (2, 3, 16)}
tf1_abw = max(abs(ring_dicht[2]["Tr F[1]"][0] + 1j*ring_dicht[2]["Tr F[1]"][1] - tf1_formel[2]),
              abs(ring_dicht[3]["Tr F[1]"][0] + 1j*ring_dicht[3]["Tr F[1]"][1] - tf1_formel[3]), abs(t1_16 - tf1_formel[16]))
ZAHLEN["austauschring"] = {"d = 2, 3": {str(k): v for k, v in ring_dicht.items()}, "d = 16: |Antwort/d^4 - W|": ring16,
                           "U(1)_Y: (Antwort/d^4, W) je phi": {str(k): [float(v[0].real), float(v[0].imag), v[1]]
                                                                for k, v in y_antwort.items()},
                           "U(1)_Y: Antwort/(d^4 phi^2) bei phi = 1e-3": y_klein, "5/24576": 5/24576,
                           "Windungen der Zyklen (alle Täusche)": windungen,
                           "Zyklen (alle Täusche)": [z for z, _ in zyk], "andere Gatterfolge: |Antwort - W|": float(alt_abw),
                           "Einzelanregungsring: |Tr F - (4 cos^2 d + sin^4 Tr g)|": bring,
                           "Tr F[1] = -(d^4 - 5 d^2 + 4 i (d^3 - d))/4: max Abweichung": float(tf1_abw)}
ok_y = all(abs(v[0] - v[1]) < 1e-12 for v in y_antwort.values())
check(f"[E] Elementarer Austauschring (vier Nähte in zyklischer Folge, Gatter -i R(1)(U), eine Anregung je Naht, "
      f"Etiketten in C^d): (Tr F[U] - Tr F[1])/d^4 = (d^2 - |Tr g|^2)/(4 d^4) = W(g) exakt (d = 2, 3 dicht gegen "
      f"Permutationsentwicklung {e0(max(v['dicht gegen Permutation'] for v in ring_dicht.values()))}, Abweichung von W "
      f"{e0(max(v['(Tr F[U] - Tr F[1])/d^4 - W'] for v in ring_dicht.values()))}; d = 16: U(16) {e0(ring16['U(16) zufällig'])}, "
      f"SU(2)_L {e0(ring16['SU(2)_L c 16'])}; U(1)_Y bei phi = 0.001, 0.3, 1 gleich W ({ok_y}), Antwort/phi^2 = "
      f"{y_klein:.10f} = 5/24576 = {5/24576:.10f}). Grund: Nur der Term mit allen vier Täuschen sieht die Holonomie (jede "
      f"echte Teilmenge ist ein Baum), und seine Permutation {[z for z, _ in zyk]} hat zwei Zyklen mit den Windungen "
      f"{windungen}: ein Etikett läuft um, die drei anderen laufen gegen es um, also Tr g Tr g^-1 = |Tr g|^2 = "
      f"chi_(V x V-quer)(g); der Vorfaktor ist (i/(1 + i))^4 = -1/4, und 1/d^4 normiert auf den gemischten "
      f"Etikettenzustand der vier Nähte (Tr F[1] = -(d^4 - 5 d^2 + 4 i (d^3 - d))/4, {e0(tf1_abw)}). Andere Gatterfolge "
      f"(01, 23, 12, 30): wieder W ({e0(alt_abw)}). Der Einzelanregungsring gibt dagegen Tr F = 4 cos^2 theta d + sin^4 "
      f"theta Tr g ({e0(max(bring.values()))}): fundamental",
      all(v["dicht gegen Permutation"] < 1e-12 and v["(Tr F[U] - Tr F[1])/d^4 - W"] < 1e-14 for v in ring_dicht.values())
      and max(ring16.values()) < 1e-14 and ok_y and abs(y_klein - 5/24576) < 1e-8 and windungen == [-1, 1]
      and alt_abw < 1e-14 and max(bring.values()) < 1e-12 and tf1_abw < 1e-8)
befund("E: Die Schleifenantwort W(g) der achtzehnten Analyse ist rekonstruiert: Sie ist die Einperioden-Spurantwort des "
       "elementaren Austauschrings mit dem nativen Gatter, normiert auf den Etikettenzustand. Ihre Form |Tr g|^2 kommt "
       "daher, dass bei voller Besetzung jedes Umlaufen eines Etiketts einen Gegenumlauf der anderen erzwingt. Das "
       "Einzelanregungsnetz gibt dagegen die fundamentale Spur Tr g. Beide sind Plakettenfunktionen; sie unterscheiden "
       "sich in Darstellung, Normierung und Vorzeichen", True)

say("  (d) W ist die adjungierte, nicht die fundamentale Wilson-Wirkung")
wid, zen, quad, quart = 0.0, [], [], []
for d in (2, 4, 16):
    g = haar(d)
    wid = max(wid, abs(W(g) - (np.trace(np.eye(d*d) - np.kron(g, g.conj())).real)/(4*d**4)))
    z = np.exp(0.7j)
    zen.append((abs(W(z*g) - W(g)), abs(np.trace(np.eye(d) - z*g).real - np.trace(np.eye(d) - g).real)))
    A = herm(d)
    A -= np.trace(A)/d*np.eye(d)
    r = []
    for eps in (1e-2, 5e-3):
        gg = expm(1j*eps*A)
        fund = np.trace(np.eye(d) - gg).real/(2*d**3)
        r.append((W(gg)/fund - 1)/eps**2)
    quad.append(abs(W(expm(1j*1e-3*A))/(np.trace(np.eye(d) - expm(1j*1e-3*A)).real/(2*d**3)) - 1))
    quart.append(r + [float(-np.trace(A @ A).real/(4*d))])
w_su2 = W(-np.eye(2))
f_su2 = np.trace(np.eye(2) + np.eye(2)).real
quart_stabil = all(abs(r[1] - r[2]) < 1e-3*abs(r[2]) and abs(r[2]) > 1e-3 for r in quart)
ZAHLEN["W_gegen_Wilson"] = {"W - Re Tr_(V x V-quer)(1 - g x g-quer)/(4 d^4)": float(wid),
                            "Zentrum: |W(zg) - W(g)|, |Delta Re Tr(1 - g)|": [[float(a), float(b)] for a, b in zen],
                            "SU(2): W(-1), Re Tr(1 - (-1))": [w_su2, float(f_su2)],
                            "quadratisch: W/(Re Tr(1 - g)/(2 d^3)) - 1 bei eps = 1e-3": [float(x) for x in quad],
                            "quartisch: (W/fund - 1)/eps^2 bei eps = 0.01, 0.005; -Tr A^2/(4 d)": [[float(x) for x in r]
                                                                                                   for r in quart]}
check(f"[E] W(g) = Re Tr_(V x V-quer)(1 - g x g-quer)/(4 d^4) identisch ({e0(wid)}; d = 2, 4, 16): eine Wilson-Wirkung in "
      f"der Darstellung V x V-quer = adjungiert + Singulett, nicht in der fundamentalen. Folgen: zentrumsblind, W(e^(0.7 i) "
      f"g) = W(g) ({e0(max(a for a, _ in zen))}), während die fundamentale Funktion sich um bis zu {max(b for _, b in zen):.3f} "
      f"ändert; auf dem Dublett ist W(-1) = {w_su2:g} gegen Re Tr(1 - (-1)) = {f_su2:g}: W sieht SO(3), nicht SU(2). "
      f"Quadratisch gilt W = Re Tr_V(1 - g)/(2 d^3) + O(A^4) für Tr A = 0 (Abweichung bei eps = 1e-3: "
      f"{e0(max(quad))}), die quartischen Terme unterscheiden sich: (W/fund - 1)/eps^2 = {[round(float(r[1]), 4) for r in quart]} "
      f"= -Tr A^2/(4 d) = {[round(float(r[2]), 4) for r in quart]} (d = 2, 4, 16)",
      wid < 1e-14 and max(a for a, _ in zen) < 1e-14 and min(b for _, b in zen) > 1e-3 and w_su2 == 0 and f_su2 == 4
      and max(quad) < 1e-5 and quart_stabil)

say("  (e) Gitter: Einzelanregungsnetz additiv, Austauschnetz ein Schleifengas")
SITES3 = [(x, y) for x in range(3) for y in range(3)]
IDX3 = {s_: i for i, s_ in enumerate(SITES3)}
A3 = [s_ for s_ in SITES3 if (s_[0] + s_[1]) % 2 == 0]
GATES3 = []
for j, b in enumerate(BS2):
    for a in A3:
        bb = (a[0] + b[0], a[1] + b[1])
        if bb in IDX3:
            GATES3.append((IDX3[a], IDX3[bb], (j, a)))
KEYS3 = [k for _, _, k in GATES3]


def kante3(u, v):
    for ia, ib, k in GATES3:
        if {ia, ib} == {IDX3[u], IDX3[v]}:
            return k


COT = {(0, 0): kante3((0, 0), (1, 0)), (1, 0): kante3((2, 0), (2, 1)), (0, 1): kante3((0, 1), (0, 2)),
       (1, 1): kante3((1, 2), (2, 2))}
LAY3 = [[(ia, ib, k) for ia, ib, k in GATES3 if k[0] == j] for j in range(4)]


def links3(hols, d):
    L_ = {k: np.eye(d, dtype=complex) for k in KEYS3}
    for p, g in hols.items():
        L_[COT[p]] = g
    return L_


d = 2
hcmp = {(0, 0): haar(d)}
dicht3 = abs(np.trace(dense_exchange(9, d, GATES3, links3(hcmp, d), A_ST, A_SW))
             - perm_trace(9, d, GATES3, links3(hcmp, d), A_ST, A_SW))
f00 = perm_trace(9, d, GATES3, links3({}, d), A_ST, A_SW)
zeilen, werte, zentrum3 = [], [], 0.0
for t in range(6):
    g = haar(d)
    f = perm_trace(9, d, GATES3, links3({(1, 1): g}, d), A_ST, A_SW)
    fz = perm_trace(9, d, GATES3, links3({(1, 1): np.exp(0.7j)*g}, d), A_ST, A_SW)
    zentrum3 = max(zentrum3, abs(f - fz))
    tg = np.trace(g)
    zeilen.append([abs(tg)**2 - d**2, tg - d, np.conj(tg) - d])
    werte.append(f - f00)
koef, *_ = np.linalg.lstsq(np.array(zeilen), np.array(werte), rcond=None)
fit_rest = float(np.linalg.norm(np.array(zeilen) @ koef - np.array(werte)))
g1, g2 = haar(d), haar(d)
gemischt = {}
for p1, p2 in [((0, 0), (1, 1)), ((0, 0), (1, 0))]:
    for name, fun in [("Einzelanregung", lambda Ls: np.trace(floquet(9, d, LAY3, Ls, TH))),
                      ("Austausch", lambda Ls: perm_trace(9, d, GATES3, Ls, A_ST, A_SW))]:
        m_ = (fun(links3({p1: g1, p2: g2}, d)) - fun(links3({p1: g1}, d)) - fun(links3({p2: g2}, d))
              + fun(links3({}, d)))
        gemischt[f"{name} {p1}-{p2}"] = float(abs(m_))
einzel = float(abs(perm_trace(9, d, GATES3, links3({(0, 0): g1}, d), A_ST, A_SW) - f00))
PLAK3 = {p: {IDX3[(p[0], p[1])], IDX3[(p[0] + 1, p[1])], IDX3[(p[0], p[1] + 1)], IDX3[(p[0] + 1, p[1] + 1)]}
         for p in COT}
zyk3 = []
for z_, weg in zyklen_voll(9, GATES3):
    pts = [SITES3[w] for w in weg]
    flaeche = 0.5*sum(pts[i][0]*pts[i + 1][1] - pts[i + 1][0]*pts[i][1] for i in range(len(pts) - 1))
    zyk3.append((len(z_), flaeche, [p for p, c in PLAK3.items() if set(weg) == c]))
zyk3.sort()
voll_ok = ([z_[0] for z_ in zyk3] == [1, 1, 1, 1, 5] and all(z_[1] == -1.0 and len(z_[2]) == 1 for z_ in zyk3[:4])
           and sorted(z_[2][0] for z_ in zyk3[:4]) == sorted(COT) and zyk3[4][1] == 4.0)
ZAHLEN["gitter_3x3"] = {"dicht gegen Permutation (512 Zustände)": float(dicht3),
                        "Einplakettenantwort = c (|Tr g|^2 - d^2): c": [float(koef[0].real), float(koef[0].imag)],
                        "Anteile Tr g, Tr g*": [float(abs(koef[1])), float(abs(koef[2]))], "Fitrest": fit_rest,
                        "zentrumsblind": float(zentrum3), "gemischte Differenzen": gemischt,
                        "Einplakettenantwort |f(g1) - f(1)|": einzel,
                        "alle tauschen: (Zykluslänge, orientierte Fläche, Plakette)": [[a, b, [str(x) for x in c]]
                                                                                        for a, b, c in zyk3]}
check(f"[E] 3 x 3-Gitter (neun Nähte, vier Plaketten, zwölf Gatter in derselben Teilschichtfolge, d = 2; dichter Fockraum "
      f"512 gegen Permutationsentwicklung {e0(dicht3)}): Eine einzelne Plakette mit Holonomie g gibt im Austauschnetz "
      f"exakt c (|Tr g|^2 - d^2) mit c = {complex(np.round(koef[0], 6))} (Fitrest {e0(fit_rest)}, Anteile Tr g und Tr g* "
      f"{e0(max(abs(koef[1]), abs(koef[2])))}; zentrumsblind {e0(zentrum3)}), also wieder V x V-quer, aber mit komplexem, "
      f"geometrieabhängigem Koeffizienten. Zwei Plaketten sind nicht additiv: gemischte Differenz f(g1, g2) - f(g1, 1) - "
      f"f(1, g2) + f(1, 1) = {gemischt['Austausch (0, 0)-(1, 1)']:.4f} (diagonal) und {gemischt['Austausch (0, 0)-(1, 0)']:.4f} "
      f"(benachbart) gegen die Einplakettenantwort {einzel:.4f}: ein Schleifengas mit Produkten von Plaketten. Im Term "
      f"'alle tauschen' umläuft je ein Etikett jede der vier Plaketten (Zykluslängen {[z_[0] for z_ in zyk3]}, Flächen "
      f"{[z_[1] for z_ in zyk3]}), und die fünf Randetiketten umlaufen gemeinsam das ganze Gitter gegenläufig: vier "
      f"Plakettenspuren mal eine Randschleife, Gesamtwindung null ({voll_ok}). Kontrolle: Das Einzelanregungsnetz auf "
      f"demselben Gitter ist additiv "
      f"({e0(max(gemischt['Einzelanregung (0, 0)-(1, 1)'], gemischt['Einzelanregung (0, 0)-(1, 0)']))})",
      voll_ok and dicht3 < 1e-10 and fit_rest < 1e-10 and max(abs(koef[1]), abs(koef[2])) < 1e-10 and zentrum3 < 1e-10
      and abs(koef[0].imag) > 0.1 and min(gemischt["Austausch (0, 0)-(1, 1)"], gemischt["Austausch (0, 0)-(1, 0)"]) > 0.05
      and max(gemischt["Einzelanregung (0, 0)-(1, 1)"], gemischt["Einzelanregung (0, 0)-(1, 0)"]) < 1e-12)
befund("E (gescopt negativ): Sum_p W(U_p) ist nicht das Gitterfunktional des Austauschnetzes. W ist exakt nur für den "
       "isolierten Ring; auf dem Gitter koppeln die Plaketten (Schleifengas, Produktterme gleicher Ordnung, komplexe "
       "Koeffizienten). Eine Wilson-Wirkung als Summe über Plaketten liefert nur das Einzelanregungsnetz, und zwar in der "
       "fundamentalen Darstellung", True)

say("  (f) 3D-Anordnung: Ebenenteilperioden xy, yz, zx")
L3 = 6                                             # L = 4 erlaubte Windungen um den Torus in einer Periode
SITES3D = [(x, y, z) for x in range(L3) for y in range(L3) for z in range(L3)]
IDX3D = {s_: i for i, s_ in enumerate(SITES3D)}
A3D = [s_ for s_ in SITES3D if sum(s_) % 2 == 0]
EV3 = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1)}


def add3(p, e, sg):
    return tuple((p[i] + sg*EV3[e][i]) % L3 for i in range(3))


LAY3D = []
for u, v in [("x", "y"), ("y", "z"), ("z", "x")]:
    for e, sg in [(u, 1), (v, 1), (u, -1), (v, -1)]:
        LAY3D.append([(IDX3D[a], IDX3D[add3(a, e, sg)], (a, e, sg)) for a in A3D])
KEYS3D = sorted({k for lay in LAY3D for _, _, k in lay})
KI3 = {k: i for i, k in enumerate(KEYS3D)}
LAYIDX3 = [(np.array([g[0] for g in lay]), np.array([g[1] for g in lay]), np.array([KI3[g[2]] for g in lay]))
           for lay in LAY3D]


def spur3d(phasen, th):
    """Einperiodenspur für U(1)-Links (d = 1), schichtweise als Zeilenoperationen."""
    c, s = np.cos(th), np.sin(th)
    vec = np.zeros(len(KEYS3D))
    for k, v in phasen.items():
        vec[KI3[k]] = v
    U = np.exp(1j*vec)
    F = np.eye(len(SITES3D), dtype=complex)
    for IA, IB, KK in LAYIDX3:
        Fa, Fb, u = F[IA].copy(), F[IB].copy(), U[KK][:, None]
        F[IA] = c*Fa - 1j*s*u*Fb
        F[IB] = c*Fb - 1j*s*np.conj(u)*Fa
    return np.trace(F)


test_links = {k: np.array([[np.exp(1j*rng.normal())]]) for k in KEYS3D}
schnell = abs(spur3d({k: np.angle(v[0, 0]) for k, v in test_links.items()}, TH)
              - np.trace(floquet(len(SITES3D), 1, LAY3D, test_links, TH)))


def plaketten3d(k):
    """Alle Plaketten (als Eckenmengen), die den Link k enthalten."""
    a, e, sg = k
    b = add3(a, e, sg)
    out = set()
    for f in "xyz":
        if f == e:
            continue
        for sf in (1, -1):
            out.add(frozenset({a, b, add3(a, f, sf), add3(b, f, sf)}))
    return out


k1 = (A3D[0], "x", 1)
P_k1 = plaketten3d(k1)
verh3d, n_ohne = {}, 0
for th in (TH, np.pi/16):
    t00 = spur3d({}, th)
    t10 = spur3d({k1: 0.9}, th)
    mit, ohne, n_ohne = 0.0, 0.0, 0
    for k2 in KEYS3D:
        if k2 == k1:
            continue
        m_ = abs(spur3d({k1: 0.9, k2: 1.3}, th) - t10 - spur3d({k2: 1.3}, th) + t00)
        if P_k1 & plaketten3d(k2):
            mit = max(mit, m_)
        else:
            ohne = max(ohne, m_)
            n_ohne += m_ > 1e-12
    verh3d[round(th/np.pi, 4)] = (mit, ohne, n_ohne)
# 2D-Gegenstück: auf dem 4 x 4-Torus ohne gemeinsame Plakette exakt null


def spur2d(phasen):
    L_ = {k: np.array([[np.exp(1j*phasen.get(k, 0.0))]]) for k in BONDS2}
    return np.trace(floquet(NS2, 1, LAY2, L_, TH))


def plaketten_der_kante(k):
    out = set()
    for x in range(NS2):
        bes = huepfpfad(x, LAY2, {kk: np.eye(1) for kk in BONDS2}, 1)[2]
        kanten = set()
        for u_, v_ in zip(bes[:-1], bes[1:]):
            for kk, (ia, ib) in BONDS2.items():
                if {ia, ib} == {u_, v_}:
                    kanten.add(kk)
        if k in kanten:
            out.add(x)
    return out


kk1 = (0, A2[0])
ohne2d = 0.0
p_k1 = plaketten_der_kante(kk1)
t00, t10 = spur2d({}), spur2d({kk1: 0.9})
for kk2 in BONDS2:
    if kk2 == kk1 or plaketten_der_kante(kk2) & p_k1:
        continue
    ohne2d = max(ohne2d, abs(spur2d({kk1: 0.9, kk2: 1.3}) - t10 - spur2d({kk2: 1.3}) + t00))
r_nat = verh3d[0.25][1]/verh3d[0.25][0]
r_klein = verh3d[0.0625][1]/verh3d[0.0625][0]
s4c8, s8c4 = SS**4*CC**8, SS**8*CC**4
ZAHLEN["netz_3d"] = {"Teilperioden": "xy, yz, zx (je +a, +b, -a, -b), Torus 6^3, d = 1",
                     "Zeilenform gegen dichte Floquetmatrix": float(schnell),
                     "max gemischte Differenz (mit, ohne gemeinsame Plakette, Zahl ohne != 0) je theta/pi":
                         {str(k): [float(v[0]), float(v[1]), int(v[2])] for k, v in verh3d.items()},
                     "Verhältnis ohne/mit bei theta = pi/4": float(r_nat), "bei theta = pi/16": float(r_klein),
                     "sin^4 cos^8 (Plakette), sin^8 cos^4 (Plakettenpaar) bei pi/4": [s4c8, s8c4],
                     "2D-Torus: ohne gemeinsame Plakette": float(ohne2d)}
check(f"[E] 3D-Anordnung (Torus 6^3, Teilperioden xy, yz, zx mit je +a, +b, -a, -b, jede Richtung zweimal je Periode; "
      f"Zeilenform gegen dichte Matrix {e0(schnell)}): Für {verh3d[0.25][2]} Partnerlinks ohne gemeinsame Plakette mit dem "
      f"Link k1 ist die gemischte Differenz nicht null; ihr Maximum relativ zu Partnern mit gemeinsamer Plakette ist "
      f"{r_nat:.3f} bei theta = pi/4 und {r_klein:.4f} bei theta = pi/16 (fällt wie sin^4). Ein Pfad umläuft in einer Periode "
      f"Plaketten verschiedener Ebenen nacheinander (Gewicht sin^8 cos^4 gegen sin^4 cos^8 der Einzelplakette, bei pi/4 "
      f"beide {s4c8:.6f}): Bei der nativen Vierteldrehung ist die 3D-Antwort keine reine Plakettenwirkung, nur im Grenzfall "
      f"kleiner Winkel. Kontrolle 2D-Torus: ohne gemeinsame Plakette exakt null ({e0(ohne2d)})",
      schnell < 1e-10 and r_nat > 0.2 and r_klein < 0.01 and verh3d[0.25][2] > 0 and abs(s4c8 - s8c4) < 1e-15
      and ohne2d < 1e-12)

say("  (g) Orientierungsungerader Teil der Netzspur auf der 16")
GAM = []
for j in range(5):
    GAM.append(c5[j] + c5[j].T)
    GAM.append(1j*(c5[j].T - c5[j]))
cliff = max(np.linalg.norm(GAM[a] @ GAM[b] + GAM[b] @ GAM[a] - 2*(a == b)*np.eye(32)) for a in range(10) for b in range(10))
so10 = []
for t in range(4):
    om = rng.normal(size=(10, 10))
    om = om - om.T
    X = (1j/4)*sum(om[a, b]*GAM[a] @ GAM[b] for a in range(10) for b in range(10))
    X16, X16q = X[np.ix_(EV, EV)], X[np.ix_(OD, OD)]
    pf = pfaff(om)
    so10.append({"herm": float(np.linalg.norm(X - X.conj().T)), "Tr X": float(abs(np.trace(X16))),
                 "Tr X^3": float(abs(np.trace(X16 @ X16 @ X16))),
                 "Tr X^5/Pf (16)": float((np.trace(np.linalg.matrix_power(X16, 5))/pf).real),
                 "Tr X^5/Pf (16-quer)": float((np.trace(np.linalg.matrix_power(X16q, 5))/pf).real),
                 "Pf^2 - det": float(abs(pf**2 - np.linalg.det(om))/abs(np.linalg.det(om)))})
yq = [Fr(0)] + [Fr(1)] + [Fr(-2, 3)]*3 + [Fr(1, 6)]*6 + [Fr(1, 3)]*3 + [Fr(-1, 2)]*2
trY = {k: sum(y**k for y in yq) for k in (1, 2, 3, 5)}
im_y = np.sum(np.sin(1e-2*np.diag(Y16)))/1e-2**5
Hsu = herm(16)
Hsu -= np.trace(Hsu)/16*np.eye(16)
su16_3 = abs(np.trace(Hsu @ Hsu @ Hsu))
ZAHLEN["chiral"] = {"Clifford": float(cliff), "so(10) Stichproben": so10, "Tr Y^k (k = 1, 2, 3, 5)": {str(k): str(v) for k, v in trY.items()},
                    "Im Tr_16 e^(i phi Y)/phi^5 bei phi = 0.01": float(im_y), "Tr Y^5/120": float(trY[5]/120),
                    "Kontrolle su(16): |Tr X^3|": float(su16_3)}
c5vals = [s_["Tr X^5/Pf (16)"] for s_ in so10]
check(f"[E] Der orientierungsungerade Teil der Netzspur ist sin^4 theta sum_p Im Tr U_p (alle Plaketten gleich orientiert). "
      f"Auf der 16 verschwindet er für jedes Spin(10)-Feld bis zur vierten Ordnung: Tr_16 X = Tr_16 X^3 = 0 (so(10) hat "
      f"keine kubische Invariante; Stichproben max {e0(max(max(s_['Tr X'], s_['Tr X^3']) for s_ in so10))}, Clifford "
      f"{e0(cliff)}), und die fünfte Ordnung ist die Pfaffsche: Tr_16 X^5 = {np.mean(c5vals):.6f} Pf(omega), Tr_16-quer "
      f"X^5 = {np.mean([s_['Tr X^5/Pf (16-quer)'] for s_ in so10]):.6f} Pf(omega) (Streuung {e0(np.ptp(c5vals))}, "
      f"Pf^2 = det bis {e0(max(s_['Pf^2 - det'] for s_ in so10))}): das Vorzeichen ist die Chiralität. Hyperladung exakt: Tr Y "
      f"= {trY[1]}, Tr Y^3 = {trY[3]} (Anomaliefreiheit), Tr Y^5 = {trY[5]}; Im Tr_16 e^(i phi Y)/phi^5 = {im_y:.6f} -> "
      f"Tr Y^5/120 = {float(trY[5]/120):.6f}. Im Kontinuum ist der Term O(a^10), irrelevant. Kontrolle: ein zufälliges "
      f"su(16)-Element hat |Tr X^3| = {su16_3:.2f}",
      cliff < 1e-12 and max(max(s_["Tr X"], s_["Tr X^3"]) for s_ in so10) < 1e-10 and np.ptp(c5vals) < 1e-8
      and abs(np.mean(c5vals)) > 1 and all(abs(s_["Tr X^5/Pf (16)"] + s_["Tr X^5/Pf (16-quer)"]) < 1e-8 for s_ in so10)
      and trY[1] == 0 and trY[3] == 0 and trY[5] == Fr(5, 9) and abs(im_y - 1/216) < 1e-4 and su16_3 > 1)

say("  (h) 1D-Kette (z120 Teil 4): keine Plakette")
L1 = 8
LAY1 = [[(i, i + 1, ("l", i)) for i in range(0, L1, 2)], [(i, (i + 1) % L1, ("l", i)) for i in range(1, L1, 2)]]
links1 = {("l", i): haar(2) for i in range(L1)}
eins1 = {k: np.eye(2, dtype=complex) for k in links1}
F1d, F1d1 = (floquet(L1, 2, LAY1, links1, TH, signs=[1, -1]), floquet(L1, 2, LAY1, eins1, TH, signs=[1, -1]))
kette = {n: float(abs(np.trace(np.linalg.matrix_power(F1d, n)) - np.trace(np.linalg.matrix_power(F1d1, n))))
         for n in (1, 2, 3, 4)}
hz = [haar(2) for _ in range(L1)]
links_pg1 = {("l", i): hz[i] @ hz[(i + 1) % L1].conj().T for i in range(L1)}
kette_pg = float(abs(np.trace(np.linalg.matrix_power(floquet(L1, 2, LAY1, links_pg1, TH, signs=[1, -1]), 4))
                     - np.trace(np.linalg.matrix_power(F1d1, 4))))


def zwei_perioden(a1, a2):
    """Tr(F[a2] F[a1]) mit flachen U(1)-Links e^(i a_t) auf allen Kanten (gleichförmiges A_x, Periode t)."""
    Fs = [floquet(L1, 1, LAY1, {("l", i): np.array([[np.exp(1j*a)]]) for i in range(L1)}, TH, signs=[1, -1])
          for a in (a1, a2)]
    return np.trace(Fs[1] @ Fs[0])


el_statisch = abs(zwei_perioden(0.8, 0.8) - zwei_perioden(0.0, 0.0))
el_fluss = abs(zwei_perioden(0.0, 0.8) - zwei_perioden(0.0, 0.0))
el_differenz = abs(zwei_perioden(0.5, 1.3) - zwei_perioden(0.0, 0.8))
ZAHLEN["kette_1d"] = {"|Tr F^n[U] - Tr F^n[1]|": {str(k): v for k, v in kette.items()},
                      "n = 4, Polyakov-Schleife trivial (reine Eichung)": kette_pg,
                      "zwei Perioden, flache Links: statisch a = a' = 0.8": float(el_statisch),
                      "zwei Perioden: a = 0 -> 0.8 (elektrischer Fluss)": float(el_fluss),
                      "zwei Perioden: nur a' - a zählt": float(el_differenz)}
check(f"[E] Kette der alternierenden Transfers (z120 Teil 4, Ring mit {L1} Nähten, zufällige U(2)-Links): Tr F^n hängt für "
      f"n = 1, 2, 3 nicht von den Links ab ({e0(max(kette[1], kette[2], kette[3]))}); erst bei n = L/2 = 4 erscheint die "
      f"Polyakov-Windung ({kette[4]:.4f}), und sie verschwindet bei reiner Eichung ({e0(kette_pg)}). In 1+1D hat das Netz "
      f"keine Plakette. Die Materie reagiert aber auf zeitartige Plaketten, wenn die Links sich von Periode zu Periode "
      f"ändern: Mit flachen Links e^(i a_t) ist Tr(F[a_2] F[a_1]) bei a_1 = a_2 = 0.8 gleich dem Wert ohne Feld "
      f"({e0(el_statisch)}), bei a_1 = 0, a_2 = 0.8 nicht ({el_fluss:.4f}), und es zählt nur a_2 - a_1 = elektrischer Fluss "
      f"({e0(el_differenz)}). Die Zeitabhängigkeit der Links liefert das Netz nicht (Teil 3a)",
      max(kette[1], kette[2], kette[3]) < 1e-12 and kette[4] > 1e-3 and kette_pg < 1e-12 and el_statisch < 1e-12
      and el_fluss > 1e-3 and el_differenz < 1e-12)

# =================================================================== Teil 2
say("\nTeil 2: Die Kontinuumskopplung")
say("  (a) Naiver Kontinuumslimes, anisotrop")


def kontinuum(Ns, Nt, kap_s, kap_t, eps=0.5, ladungen=np.diag(Y16)):
    a_s, a_t = 1/Ns, 1/Nt
    T = float(np.sum(ladungen**2))
    x2 = np.arange(Ns)*a_s
    th_m = a_s*eps*(np.sin(2*np.pi*x2) - np.sin(2*np.pi*(x2 + a_s)))
    S_mag = kap_s*Ns**2*Nt*np.sum(1 - np.cos(np.outer(th_m, ladungen)))
    t = np.arange(Nt)*a_t
    th_e = a_s*eps*(np.sin(2*np.pi*t) - np.sin(2*np.pi*(t + a_t)))
    S_el = kap_t*Ns**3*np.sum(1 - np.cos(np.outer(th_e, ladungen)))
    xi = a_s/a_t
    int_F2 = 2*np.pi**2*eps**2
    return S_mag/(kap_s*T*xi/2*int_F2), S_el/(kap_t*T/(2*xi)*int_F2)


konti = {}
for Ns, Nt in [(32, 64), (64, 128), (128, 256), (128, 64)]:
    konti[f"N_s = {Ns}, N_t = {Nt}"] = kontinuum(Ns, Nt, 0.25, 0.7)
abw_k = [float(max(abs(v[0] - 1), abs(v[1] - 1))) for v in konti.values()]
ordnung = np.log(abw_k[0]/abw_k[1])/np.log(2)
ZAHLEN["kontinuum"] = {"S_lat/S_kont (magnetisch, elektrisch)": {k: [float(v[0]), float(v[1])] for k, v in konti.items()},
                       "Konvergenzordnung": float(ordnung)}
check(f"[E] Naiver Kontinuumslimes für S = kappa_s sum_(räuml. p) Re Tr_R(1 - U_p) + kappa_t sum_(zeitart. p) Re Tr_R(1 - "
      f"U_p) auf einem anisotropen 4D-Gitter (a_s, a_t, xi = a_s/a_t; R = Hyperladung der 16, T_R = Tr Y^2 = 10/3; glatte "
      f"magnetische und elektrische Testfelder): S_lat/S_kont -> 1 mit 1/g^2 = kappa_s T_R xi (magnetisch) und "
      f"kappa_t T_R/xi (elektrisch); Abweichungen {[round(x, 6) for x in abw_k]} für (N_s, N_t) = {list(konti)}, Ordnung "
      f"{ordnung:.2f} (a^2). Isotrop (xi = 1): 1/g^2 = kappa T_R; Lorentzinvarianz verlangt kappa_s xi = kappa_t/xi, also "
      f"g^2 = 1/(T_R sqrt(kappa_s kappa_t))",
      max(abw_k[2:]) < 1e-3 and abw_k[0] > abw_k[1] > abw_k[2] and abs(ordnung - 2) < 0.2)

say("  (b) Verhältnisse der Faktoren: Spuren über die 16")
tr3 = np.array([[np.trace(A_ @ B_).real for B_ in T3_16] for A_ in T3_16])
tr2 = np.array([[np.trace(A_ @ B_).real for B_ in T2_16] for A_ in T2_16])
trYY = float(np.trace(Y16 @ Y16))
misch = max(abs(np.trace(A_ @ Y16)) for A_ in T3_16 + T2_16)
k_i = {"SU(3)": Fr(2)*Fr(2), "SU(2)": Fr(2)*Fr(2), "U(1)_Y (GUT, T_1 = sqrt(3/5) Y)": Fr(2)*Fr(3, 5)*trY[2]}


def sin2(trT3, trYY_):
    return trT3/(trT3 + trYY_)


sin2_16 = sin2(Fr(2), trY[2])
adj = []
for X_ in [T3_16[7], T2_16[2], Y16]:
    adj.append((np.trace(np.linalg.matrix_power(np.kron(X_, np.eye(16)) - np.kron(np.eye(16), X_.conj()), 2)).real,
                2*16*np.trace(X_ @ X_).real))
sin2_W = (adj[1][0])/(adj[1][0] + adj[2][0])
PHI0 = 1e-3
schleife = {}
for name, X_ in [("SU(3) lambda_3/2", T3_16[2]), ("SU(3) lambda_8/2", T3_16[7]), ("SU(2) T_3", T2_16[2]), ("U(1) Y", Y16)]:
    gg = expm(1j*PHI0*X_)
    schleife[name] = (W(gg)/PHI0**2, SS**4*np.trace(np.eye(16) - gg).real/PHI0**2)
SIN2_MZ = 0.23122                                   # PDG, MS-quer bei M_Z: nur Vergleich, nicht verwendet
traeger = {"Tr_5 Y^2": Fr(5, 6), "Tr_5 T_3^2": Fr(1, 2)}
sin2_5 = sin2(traeger["Tr_5 T_3^2"], traeger["Tr_5 Y^2"])
dubl_Y = 6*Fr(1, 36) + 2*Fr(1, 4)
sin2_dubl = sin2(Fr(2), dubl_Y)
farb_Y = 6*Fr(1, 36) + 3*Fr(4, 9) + 3*Fr(1, 9)
ZAHLEN["verhaeltnisse"] = {"Tr_16 T_a T_b SU(3) (Diagonale)": [float(x) for x in np.diag(tr3)],
                           "Tr_16 T_a T_b SU(2) (Diagonale)": [float(x) for x in np.diag(tr2)], "Tr_16 Y^2": trYY,
                           "Nebendiagonale max": float(max(np.abs(tr3 - np.diag(np.diag(tr3))).max(), np.abs(tr2 - np.diag(np.diag(tr2))).max())),
                           "Mischspuren Tr(T Y)": float(misch), "Niveaus k_i = 2 Tr_16 T_i^2": {k: str(v) for k, v in k_i.items()},
                           "sin^2 theta_W (16)": str(sin2_16), "W-Lesart: T_(V x V-quer) = 2 d T_V": [[float(a), float(b)] for a, b in adj],
                           "sin^2 theta_W (W-Lesart)": float(sin2_W), "Träger-5": str(sin2_5), "nur Dubletts (Q, L)": str(sin2_dubl),
                           "nur Farbträger (Q, u^c, d^c): Tr Y^2": str(farb_Y),
                           "je Schleife/phi^2 (W, Netzspur kappa = 1/4)": {k: [float(a), float(b)] for k, v in schleife.items()
                                                                           for a, b in [v]},
                           "sin^2 theta_W(M_Z) gemessen (nur Vergleich)": SIN2_MZ}
check(f"[E] Spuren über die 16 (Lambda^gerade C^5, z82): Tr_16(T_a T_b) = 2 delta_ab für SU(3) und SU(2) (Diagonalen "
      f"{sorted({float(x) for x in np.round(np.diag(tr3), 12)})}, {sorted({float(x) for x in np.round(np.diag(tr2), 12)})}, Nebendiagonale "
      f"{e0(ZAHLEN['verhaeltnisse']['Nebendiagonale max'])}, Mischspuren mit Y {e0(misch)}), Tr_16 Y^2 = {trY[2]} = "
      f"{trYY:.12f}: Tr Y^2 : Tr T_3^2 = 10/3 : 2 = 5/3, die SU(5)-Normierung. Jede Spurantwort über die vollständige 16 gibt "
      f"1/g_i^2 = lambda k_i mit dem induzierten Niveau k_i = 2 Tr_16 T_i^2 = {', '.join(str(v) for v in k_i.values())} "
      f"(SU(3), SU(2), U(1)_Y GUT-normiert) aus den Niveau-1-Nahtströmen, also g_3 = g_2 = g_1 und sin^2 theta_W = "
      f"Tr T_3^2/(Tr T_3^2 + Tr Y^2) = {sin2_16}, unabhängig von theta. In der W-Lesart ist T_(V x V-quer) = 2 d T_V für "
      f"spurfreie Erzeuger ({[(round(float(a), 6), round(float(b), 6)) for a, b in adj]}), dasselbe Verhältnis ({sin2_W:.12f}). "
      f"Je Schleife und phi^2: W gibt {', '.join(f'{k} {v[0]:.10f}' for k, v in schleife.items())} (= 1/8192 für die "
      f"nichtabelschen Erzeuger, 5/24576 für Y), die Netzspur {', '.join(f'{v[1]:.6f}' for v in schleife.values())} (= 1/4 "
      f"bzw. 5/12); das Verhältnis Y : T ist in beiden {schleife['U(1) Y'][0]/schleife['SU(2) T_3'][0]:.6f} = 5/3. "
      f"Kontrollen: der Träger 5 (vollständiges SU(5)-Multiplett) gibt wieder {sin2_5}; nur die Dubletts Q, L als Etiketten "
      f"geben {sin2_dubl}: Das Verhältnis 3/8 hängt an der Vollständigkeit der Etiketten, nicht an einer Wahl. Vergleich: "
      f"gemessen sin^2 theta_W(M_Z) = {SIN2_MZ}; 3/8 gilt nur an der Anpassungsskala, der Weg dorthin ist das Laufen (z25)",
      np.allclose(np.diag(tr3), 2) and np.allclose(np.diag(tr2), 2) and ZAHLEN["verhaeltnisse"]["Nebendiagonale max"] < 1e-12
      and misch < 1e-12 and trY[2] == Fr(10, 3) and abs(trYY - 10/3) < 1e-12 and all(v == 4 for v in k_i.values())
      and sin2_16 == Fr(3, 8) and all(abs(a - b) < 1e-9 for a, b in adj) and abs(sin2_W - 3/8) < 1e-12
      and sin2_5 == Fr(3, 8) and sin2_dubl == Fr(3, 4)
      and all(abs(v[0] - 1/8192) < 1e-9 and abs(v[1] - 0.25) < 1e-6 for k, v in schleife.items() if k != "U(1) Y")
      and abs(schleife["U(1) Y"][0] - 5/24576) < 1e-9 and abs(schleife["U(1) Y"][1] - 5/12) < 1e-6)

say("  (c) Absolutwert: drei native Lesarten von lambda (nur Vergleich, keine Anpassung)")
lesarten = {"R1 Netzspur je Periode (kappa = sin^4 theta, P_exp, xi = 1)": SS**4/2,
            "R2 W je Etikettenzustand (kappa = 1/(4 d^4) auf V x V-quer, d = 16)": 1/(4*16**3),
            "R3 Niveau k = 1 mit P_norm (z110: g_eff^2 = 1/k)": 1.0}
absolut = {}
for name, lam in lesarten.items():
    inv_g2 = lam*4
    aU = 4*np.pi*inv_g2
    # drei Generationen: linear in der Spur (R1, R3); für W hängt es am Etikettenraum (d = 48 gemeinsam oder 3 x 16)
    drei = [3*aU] if not name.startswith("R2") else [4*np.pi*6/(2*48**3), 3*aU]
    absolut[name] = {"lambda": lam, "1/g_U^2": inv_g2,
                     "beta_lat (Wilson, xi = 1): SU(3) 6/g^2, SU(2) 4/g^2, U(1)_Y 1/g_Y^2": [6*inv_g2, 4*inv_g2, 5/3*inv_g2],
                     "alpha_U^-1": aU, "alpha_U^-1 drei Generationen": drei,
                     "alpha_em^-1 (sin^2 = 3/8, ohne Laufen)": 8/3*aU, "alpha^-1/alpha_em^-1": ALPHA_INV/(8/3*aU)}
werte_g2 = [v["1/g_U^2"] for v in absolut.values()]
werte_em = [v["alpha_em^-1 (sin^2 = 3/8, ohne Laufen)"] for v in absolut.values()]
spanne = max(werte_g2)/min(werte_g2)
treffer = [n for n, v in absolut.items() if abs(v["alpha^-1/alpha_em^-1"] - 1) < 1e-2]
naechst = min(absolut.values(), key=lambda v: abs(v["alpha^-1/alpha_em^-1"] - 1))
z25 = [n for n, v in absolut.items() if any(60 <= a <= 74 for a in [v["alpha_U^-1"]] + v["alpha_U^-1 drei Generationen"])]
ZAHLEN["absolutwert"] = {"Lesarten": absolut, "Spanne 1/g^2 (max/min)": spanne, "Treffer alpha^-1 innerhalb 1 %": treffer,
                         "im z25-Fenster alpha_U^-1(1/a) = 60 ... 74": z25, "alpha^-1 (gelöst, nur Vergleich)": ALPHA_INV}
for name, v in absolut.items():
    bl = v["beta_lat (Wilson, xi = 1): SU(3) 6/g^2, SU(2) 4/g^2, U(1)_Y 1/g_Y^2"]
    say(f"   {name}: 1/g_U^2 = {v['1/g_U^2']:.6g}, beta_lat (SU(3), SU(2), U(1)_Y) = "
        f"({', '.join(f'{x:.6g}' for x in bl)}), alpha_U^-1 = {v['alpha_U^-1']:.6g} (drei Generationen"
        f"{' (d = 48 gemeinsam bzw. 3 x d = 16)' if name.startswith('R2') else ''} "
        f"{', '.join(f'{a:.4g}' for a in v['alpha_U^-1 drei Generationen'])}), alpha_em^-1 = "
        f"{v['alpha_em^-1 (sin^2 = 3/8, ohne Laufen)']:.6g}, alpha^-1/alpha_em^-1 = {v['alpha^-1/alpha_em^-1']:.6g}")
check(f"[C/L] Absolutwert nicht nativ: Mit denselben nativen Daten (d = 16, theta = pi/4, k = 1, also k_i = 4) geben die drei "
      f"Lesarten 1/g_U^2 = lambda k_i = {', '.join(f'{x:.6g}' for x in werte_g2)} (Spanne Faktor {spanne:.0f}); "
      f"alpha_em^-1 an der Anpassungsskala (sin^2 theta_W = 3/8, ohne Laufen) = {', '.join(f'{x:.4f}' for x in werte_em)} "
      f"gegen alpha^-1 = {ALPHA_INV}: kein Treffer innerhalb 1 %. Am nächsten liegt R3 mit 128 pi/3 = "
      f"{naechst['alpha_em^-1 (sin^2 = 3/8, ohne Laufen)']:.4f} ({100*(naechst['alpha^-1/alpha_em^-1'] - 1):.1f} % "
      f"daneben), aber als Wert vor jedem Laufen und für eine Generation (drei Generationen: Faktor 3); das ist keine "
      f"Übereinstimmung. Das z25-Fenster alpha_U^-1(1/a) = 60 ... 74 trifft keine Lesart, auch nicht mit drei Generationen "
      f"(R1: {absolut[list(absolut)[0]]['alpha_U^-1']:.2f} bzw. {absolut[list(absolut)[0]]['alpha_U^-1 drei Generationen'][0]:.2f}, "
      f"R3: {naechst['alpha_U^-1']:.2f} bzw. {naechst['alpha_U^-1 drei Generationen'][0]:.2f}); R1 ist wie die "
      f"y1-Induktion zu stark",
      spanne > 1000 and not treffer and not z25)
befund("C: Nativ ist das Verhältnis (k_3 = k_2 = k_1 = 4, sin^2 theta_W = 3/8 an der Anpassungsskala), nicht der Betrag. "
       "lambda hängt an der Lesart der Antwort als Wirkung (P_exp: Spur je Periode, je Etikettenzustand oder Stromnormierung "
       "P_norm) und an der Anisotropie xi; diese Wahl verschiebt 1/g^2 über vier Größenordnungen. Mit z25 überlebt "
       "'gleich an 1/a' nur als nackte SO(10)-Gitterkopplung", True)

# =================================================================== Teil 3
say("\nTeil 3: Dynamik")
say("  (a) Z4-Quantenlinks: das Netz misst den Fluss, es bewegt ihn nicht")
NL, NSR = 4, 4
DL = 4**NL
CLOCK = np.diag([1, 1j, -1, -1j])
SHIFT = np.roll(np.eye(4), 1, axis=0)
RINGP = [(0, 1), (1, 2), (2, 3), (3, 0)]


def op_link(e, O):
    mats = [np.eye(4)]*NL
    mats[e] = O
    return np.kron(reduce(np.kron, mats), np.eye(NSR))


def ctrl(e):
    p, q = RINGP[e]
    out = np.zeros((DL*NSR, DL*NSR), dtype=complex)
    for m in range(4):
        u = 1j**m
        G = np.eye(NSR, dtype=complex)
        G[p, p] = G[q, q] = CC
        G[p, q], G[q, p] = -1j*SS*u, -1j*SS*np.conj(u)
        proj = [np.eye(4)]*NL
        P_ = np.zeros((4, 4))
        P_[m, m] = 1
        proj[e] = P_
        out += np.kron(reduce(np.kron, proj), G)
    return out


FQ = reduce(lambda A_, B_: B_ @ A_, [ctrl(e) for e in range(NL)], np.eye(DL*NSR, dtype=complex))


def gauss(x):
    G = np.eye(DL*NSR, dtype=complex)
    for e, (p, q) in enumerate(RINGP):
        if p == x:
            G = G @ op_link(e, SHIFT)
        if q == x:
            G = G @ op_link(e, SHIFT.T)
    ph = np.eye(NSR, dtype=complex)
    ph[x, x] = 1j
    return G @ np.kron(np.eye(DL), ph)


kU = max(np.linalg.norm(FQ @ op_link(e, CLOCK) - op_link(e, CLOCK) @ FQ) for e in range(NL))
kX = max(np.linalg.norm(FQ @ op_link(e, SHIFT) - op_link(e, SHIFT) @ FQ) for e in range(NL))
kG = max(np.linalg.norm(FQ @ gauss(x) - gauss(x) @ FQ) for x in range(NSR))
psi_l = rng.normal(size=DL) + 1j*rng.normal(size=DL)
psi_l /= np.linalg.norm(psi_l)
psi0 = np.kron(psi_l, np.eye(NSR)[0])
HE = sum(op_link(e, SHIFT) + op_link(e, SHIFT.T) for e in range(NL))
UE = expm(-1j*0.3*HE)
ps, pe = psi0.copy(), psi0.copy()
for _ in range(12):
    ps = FQ @ ps
    pe = UE @ (FQ @ pe)
bes_s = (np.abs(ps.reshape(DL, NSR))**2).sum(1)
bes_e = (np.abs(pe.reshape(DL, NSR))**2).sum(1)
rho_l = ps.reshape(DL, NSR) @ ps.reshape(DL, NSR).conj().T
rein_l = float(np.trace(rho_l @ rho_l).real)
ZAHLEN["quantenlinks_z4"] = {"||[F, U_e]||": float(kU), "||[F, X_e]||": float(kX), "||[F, G_x]||": float(kG),
                             "Besetzungsänderung statisch (12 Perioden)": float(np.abs(bes_s - np.abs(psi_l)**2).max()),
                             "mit elektrischem Term e^(-0.3 i sum(X + X^+))": float(np.abs(bes_e - np.abs(psi_l)**2).max()),
                             "Reinheit des reduzierten Linkzustands nach 12 Perioden": rein_l}
check(f"[E] Z4-Quantenlinks (die Uhr Z4 als Linkgruppe; Ring aus vier Nähten mit einer Anregung, Linkraum C^4 je Kante, "
      f"1024 Zustände, Netz als gesteuertes Gatter): F vertauscht mit allen Linkoperatoren U_e ({e0(kU)}) und mit dem "
      f"Gauß-Gesetz G_x = (Linkverschiebungen an x) x (Ladungsphase an x) ({e0(kG)}), nicht mit den Verschiebungen X_e "
      f"({kX:.2f}). Die Besetzungen der Linkkonfigurationen bleiben über 12 Perioden exakt "
      f"({e0(ZAHLEN['quantenlinks_z4']['Besetzungsänderung statisch (12 Perioden)'])}), nur die Kohärenzen zerfallen "
      f"(Reinheit {rein_l:.4f}): Das Netz misst den Fluss, es bewegt ihn nicht. Kontrolle: mit einem elektrischen Term "
      f"e^(-0.3 i sum (X_e + X_e^+)) ändern sich die Besetzungen um bis zu "
      f"{ZAHLEN['quantenlinks_z4']['mit elektrischem Term e^(-0.3 i sum(X + X^+))']:.4f}",
      kU < 1e-12 and kG < 1e-12 and kX > 1 and ZAHLEN["quantenlinks_z4"]["Besetzungsänderung statisch (12 Perioden)"] < 1e-12
      and ZAHLEN["quantenlinks_z4"]["mit elektrischem Term e^(-0.3 i sum(X + X^+))"] > 1e-3 and rein_l < 0.99)

say("  (b) Reflexionspositivität über die Charakterentwicklung der Plakettengewichte")
MQ = 4096
PHI = 2*np.pi*np.arange(MQ)/MQ


def fourier(w, nmax=6):
    return [float(np.real(np.mean(w*np.exp(-1j*n*PHI)))) for n in range(-nmax, nmax + 1)]


xq, wq = np.polynomial.legendre.leggauss(1200)
THQ = (xq + 1)*np.pi/2
MASS = (2/np.pi)*np.sin(THQ)**2*wq*np.pi/2


def charakter(w, j2max=8):
    return [float(np.sum(MASS*w*np.sin((j2 + 1)*THQ)/np.sin(THQ))) for j2 in range(j2max + 1)]


kap = SS**4
rp = {"U(1) fundamental kappa = +1/4": fourier(np.exp(kap*np.cos(PHI))),
      "U(1) fundamental kappa = -1/4": fourier(np.exp(-kap*np.cos(PHI))),
      "SU(2) fundamental kappa = +1/4": charakter(np.exp(kap*2*np.cos(THQ))),
      "SU(2) fundamental kappa = -1/4": charakter(np.exp(-kap*2*np.cos(THQ))),
      "SU(2) adjungiert e^(+|Tr U|^2/4)": charakter(np.exp(kap*(2*np.cos(THQ))**2)),
      "SU(2) W-Gewicht e^(-W)": charakter(np.exp(-(4 - (2*np.cos(THQ))**2)/(4*16)))}
bessel = max(abs(a - iv(abs(n), kap)) for a, n in zip(rp["U(1) fundamental kappa = +1/4"], range(-6, 7)))
pos = {k: min(v) for k, v in rp.items()}
halb = max(abs(x) for k in ("SU(2) adjungiert e^(+|Tr U|^2/4)", "SU(2) W-Gewicht e^(-W)") for x in rp[k][1::2])
ganz = min(x for k in ("SU(2) adjungiert e^(+|Tr U|^2/4)", "SU(2) W-Gewicht e^(-W)") for x in rp[k][0::2])
sigma = -np.log(rp["SU(2) fundamental kappa = +1/4"][1]/(2*rp["SU(2) fundamental kappa = +1/4"][0]))
ZAHLEN["reflexionspositivitaet"] = {"Koeffizienten": rp, "min je Gewicht": pos, "Bessel I_n(1/4)": float(bessel),
                                    "halbzahlige Koeffizienten adjungiert/W": float(halb),
                                    "Stringspannung 1+1D, fundamental, kappa = 1/4": float(sigma)}
check(f"[E] Osterwalder-Seiler-Kriterium: Die Transfermatrix der Linkvariablen (temporale Eichung) ist positiv, wenn das "
      f"zeitartige Plakettengewicht positiv definit ist, also nichtnegative Charakterkoeffizienten hat (für die Faltung auf "
      f"einer Kante ist das auch notwendig). U(1) mit e^(kappa cos): "
      f"Koeffizienten = I_n(kappa) ({e0(bessel)}), min {pos['U(1) fundamental kappa = +1/4']:.2e} > 0 für kappa = +1/4, "
      f"{pos['U(1) fundamental kappa = -1/4']:.3f} < 0 für -1/4; SU(2) fundamental min {pos['SU(2) fundamental kappa = +1/4']:.2e} "
      f"gegen {pos['SU(2) fundamental kappa = -1/4']:.3f}; adjungiert e^(|Tr U|^2/4) und W-Gewicht e^(-W): ganzzahlige "
      f"Koeffizienten min {ganz:.2e} > 0. Positiv ist also jedes Gewicht e^(kappa chi) mit kappa > 0 und einer positiven "
      f"Summe von Charakteren chi, fundamental wie adjungiert. Beim adjungierten und beim W-Gewicht sind alle halbzahligen Koeffizienten null ({e0(halb)}): "
      f"Fundamentale (Dublett-)Wilson-Schleifen verschwinden, die Theorie sieht nur SO(3). Fundamental bei kappa = 1/4 ist "
      f"die 1+1D-Stringspannung {sigma:.3f} (starke Kopplung)",
      bessel < 1e-12 and pos["U(1) fundamental kappa = +1/4"] > 0 and pos["U(1) fundamental kappa = -1/4"] < -1e-3
      and pos["SU(2) fundamental kappa = +1/4"] > 0 and pos["SU(2) fundamental kappa = -1/4"] < -1e-3
      and ganz > 0 and halb < 1e-12)

say("  (c) Vorzeichen der Antwort je Lesart")
vorz = {}
for d in (2, 3, 16):
    g = expm(1j*0.2*herm(d))
    links = {"e01": g, "e12": np.eye(d, dtype=complex), "e23": np.eye(d, dtype=complex), "e30": np.eye(d, dtype=complex)}
    eins = {k: np.eye(d, dtype=complex) for k in links}
    tX, tX1 = perm_trace(4, d, RING, links, A_ST, A_SW), perm_trace(4, d, RING, eins, A_ST, A_SW)
    tB, tB1 = (np.trace(floquet(4, d, lay_ring, links, TH)), np.trace(floquet(4, d, lay_ring, eins, TH)))
    vorz[d] = {"Austausch: Re Delta Tr F": float((tX - tX1).real), "Austausch: Delta |Tr F|": float(abs(tX) - abs(tX1)),
               "Einzelanregung: Re Delta Tr F": float((tB - tB1).real), "Einzelanregung: Delta |Tr F|": float(abs(tB) - abs(tB1))}
ZAHLEN["vorzeichen"] = {str(k): v for k, v in vorz.items()}
check(f"[E] Vorzeichen (Fluss auf dem Ring, kleine Holonomie): Das Einzelanregungsnetz ist in beiden Lesarten diamagnetisch, "
      f"Re Delta Tr F < 0 und Delta |Tr F| < 0 für d = 2, 3, 16 "
      f"({[(round(v['Einzelanregung: Re Delta Tr F'], 5), round(v['Einzelanregung: Delta |Tr F|'], 5)) for v in vorz.values()]}): "
      f"kappa > 0, reflexionspositiv. Der Austauschring gewinnt dagegen Spur, Re Delta Tr F = +d^4 W > 0 für alle d "
      f"({[round(v['Austausch: Re Delta Tr F'], 5) for v in vorz.values()]}), und der Betrag |Tr F| steigt bei d = 2, fällt "
      f"aber bei d = 3, 16 ({[round(v['Austausch: Delta |Tr F|'], 5) for v in vorz.values()]}), weil Re Tr F[1] = -(d^4 - 5 "
      f"d^2)/4 bei d^2 = 5 das Vorzeichen wechselt. Das Vorzeichen einer Wirkung aus W ist also eine Lesart (P_exp), nicht "
      f"nativ; in der linearen Lesart wäre kappa < 0",
      all(v["Einzelanregung: Re Delta Tr F"] < 0 and v["Einzelanregung: Delta |Tr F|"] < 0 and v["Austausch: Re Delta Tr F"] > 0
          for v in vorz.values()) and vorz[2]["Austausch: Delta |Tr F|"] > 0 and vorz[3]["Austausch: Delta |Tr F|"] < 0
      and vorz[16]["Austausch: Delta |Tr F|"] < 0)

say("  (d) Lichtkegel des 2D-Netzes")


def bloch(kx, ky, th):
    c, s = np.cos(th), np.sin(th)
    F = np.eye(2, dtype=complex)
    for b in BS2:
        ph = np.exp(1j*(kx*b[0] + ky*b[1]))
        F = np.array([[c, -1j*s*ph], [-1j*s*np.conj(ph), c]]) @ F
    return F


KX, KY = np.meshgrid(np.linspace(-np.pi, np.pi, 161), np.linspace(-np.pi, np.pi, 161), indexing="ij")


def quasienergie(th):
    """eps(k) >= 0 auf dem Gitter; det F = 1, also cos eps = Re Tr F/2."""
    c, s = np.cos(th), np.sin(th)
    F = np.broadcast_to(np.eye(2, dtype=complex), KX.shape + (2, 2)).copy()
    for b in BS2:
        ph = np.exp(1j*(KX*b[0] + KY*b[1]))
        St = np.zeros(KX.shape + (2, 2), dtype=complex)
        St[..., 0, 0] = St[..., 1, 1] = c
        St[..., 0, 1], St[..., 1, 0] = -1j*s*ph, -1j*s*np.conj(ph)
        F = St @ F
    return np.arccos(np.clip(np.real(F[..., 0, 0] + F[..., 1, 1])/2, -1, 1))


luecke_pi, linie0 = {}, {}
for j in range(1, 8):
    th = j*np.pi/16
    eps_k = quasienergie(th)
    luecke_pi[f"{j}/16"] = float(np.min(np.pi - eps_k))
    m0 = 0.0
    for kx in np.linspace(0, np.pi, 9):
        ev = np.angle(np.linalg.eigvals(bloch(kx, np.pi - kx, th)))
        m0 = max(m0, float(min(abs(ev))))
    linie0[f"{j}/16"] = m0
kegel = []
for ang in (0, np.pi/8, np.pi/4, np.pi/2, 3*np.pi/4):
    ev = np.angle(np.linalg.eigvals(bloch(1e-4*np.cos(ang), 1e-4*np.sin(ang), TH)))
    kegel.append(float(np.mean(np.pi - np.abs(ev))/1e-4))
quer = float(np.mean(np.abs(np.angle(np.linalg.eigvals(bloch(np.pi/2 + 1e-4/np.sqrt(2), np.pi/2 + 1e-4/np.sqrt(2), TH)))))/1e-4)
ZAHLEN["lichtkegel_2d"] = {"min Abstand zu eps = pi je theta/pi": luecke_pi, "eps = 0 auf |k_x| + |k_y| = pi: max |eps|": linie0,
                           "Kegel bei eps = pi, k = 0: Steigung je Richtung": kegel, "eps = 0-Linie: Quersteigung": quer}
andere = [v for k, v in luecke_pi.items() if k != "4/16"]
check(f"[N/L] Das 2D-Netz hat bei der nativen Vierteldrehung zwei verschiedene masselose Sektoren: Die Lücke bei eps = pi "
      f"schließt sich nur bei theta = pi/4 (min Abstand zu pi bei theta/pi = 1/16 ... 7/16, Gitter 161^2: "
      f"{[round(v, 4) for v in luecke_pi.values()]}) in einem isotropen Kegel bei k = 0 mit Steigung "
      f"{[round(x, 6) for x in kegel]} = sqrt 2 Plätze je Periode; bei eps = 0 liegt für jedes theta die genestete Linie "
      f"|k_x| + |k_y| = pi (max |eps| {e0(max(linie0.values()))}; Quersteigung bei theta = pi/4 {quer:.6f}), eine Fermilinie "
      f"des Quadratgitters, nicht lorentzartig. Ein einziger nativer Lichtkegel, der xi festlegen könnte, existiert nicht",
      luecke_pi["4/16"] < 1e-6 and min(andere) > 0.1 and np.allclose(kegel, np.sqrt(2), atol=1e-6)
      and max(linie0.values()) < 1e-10 and abs(quer - np.sqrt(2)) < 1e-5)
befund("E/C: Das Netz liefert eine eichinvariante, statische Funktion der Links (Fluss wird gemessen, nicht bewegt). Eine "
       "reflexionspositive Transfermatrix mit Links als Variablen braucht vier Zutaten, die das Netz nicht hat: (1) P_link: "
       "Links als eigene Freiheitsgrade mit Haarmaß; (2) P_E: einen elektrischen Term bzw. zeitartige Plaketten mit "
       "kappa_t > 0. Je Periode hat das Netz nur räumliche Zyklen (in 1+1D gar keine); auf zeitabhängige Links reagiert die "
       "Materie über Mehrperiodenschleifen (Teil 1h), erzeugt sie aber nicht. Mit P_link würde ein elektrischer Term "
       "höchstens induziert (Sakharov-Weg von y1), und dessen Wert an 1/a ist nach z25 zu stark; (3) P_exp: die Regel "
       "Antwort -> Wirkung mit Vorzeichen und Koeffizient (für W lesartabhängig); (4) xi bzw. P_Lorentz, denn g^2 = "
       "1/(T_R sqrt(kappa_s kappa_t)) und das Netz hat keinen einzigen Lichtkegel. In 3D kommt (5) hinzu: Bei theta = pi/4 "
       "ist die Antwort keine reine Plakettenwirkung. Die native Reflexion (z84) ist die thermische Klasse; eine "
       "Gitterspiegelebene für Linkvariablen gehört nicht dazu", True)

# =================================================================== Teil 4
say("\nTeil 4: Bilanz")
BILANZ = [
    ("E", "Link im Transferstrom: z120-Gatter eichkovariant, Einzelbindung reine Eichung; auf einer Anregung je Naht "
          "-i R(1)(U) mit Etikettentransport, auf einer Anregung die Vierteldrehung B"),
    ("E", "2D-Netz der Vierteldrehungen: <x|F|x> = cos^4 + sin^4 U_p, jede Plakette einmal je Periode, gleich orientiert; "
          "Re Tr(F[1] - F[U]) = sin^4 theta sum_p Re Tr(1 - U_p): fundamentale Wilson-Funktion, kappa = 1/4, additiv"),
    ("E", "Austauschring: (Tr F[U] - Tr F[1])/d^4 = W(g) exakt; |Tr g|^2 aus Umlauf und Gegenumlauf; 5/24576 phi^2 für Y; "
          "Normierung 1/(4 d^4) = (i/(1+i))^4 mal Etikettennormierung"),
    ("E", "W = adjungierte Wilson-Wirkung (V x V-quer): zentrumsblind, auf dem Dublett SO(3); quadratisch = fundamentale "
          "mit 1/(2 d^3)"),
    ("E (negativ)", "Austauschnetz auf dem Gitter: Schleifengas, nicht additiv, komplexe Koeffizienten; sum_p W(U_p) ist "
                    "nicht sein Funktional"),
    ("E (negativ, gescopt)", "3D mit Ebenenteilperioden: bei theta = pi/4 längere Schleifen gleichen Gewichts; reine "
                             "Plakettenwirkung nur für kleine Winkel"),
    ("E", "Orientierungsungerader Teil auf der 16: Tr X = Tr X^3 = 0 für so(10), Tr X^5 = -+60 Pf; Tr Y^5 = 5/9; im Kontinuum "
          "irrelevant"),
    ("E (negativ)", "1D-Kette: keine Plakette, nur Polyakov-Windungen; zeitartige Plaketten sieht die Materie erst über zwei "
                    "Perioden mit zeitabhängigen Links (nur a_2 - a_1 zählt), die das Netz nicht erzeugt"),
    ("E", "Kontinuum: 1/g^2 = kappa_s T_R xi = kappa_t T_R/xi; Lorentz: g^2 = 1/(T_R sqrt(kappa_s kappa_t))"),
    ("E", "Verhältnis nativ: k_i = 2 Tr_16 T_i^2 = 4 für SU(3), SU(2), U(1)_Y (GUT), sin^2 theta_W = 3/8 in jeder Lesart; "
          "an die Vollständigkeit der Etiketten gebunden (Dubletts allein: 3/4)"),
    ("C/L (negativ)", "Absolutwert: lambda je Lesart 1/8 (Netzspur), 1/16384 (W), 1 (P_norm): 1/g_U^2 = 1/2, 1/4096, 4; "
                      "kein Treffer bei alpha^-1, keiner im z25-Fenster"),
    ("E", "Dynamik: Netz statisch (misst den Fluss), Gauß-Gesetz exakt; RP für kappa > 0 fundamental wie adjungiert; "
          "adjungiert ohne halbzahlige Charaktere"),
    ("E", "Vorzeichen: Einzelanregungsnetz diamagnetisch; Austauschring lesartabhängig (d = 2 gegen d = 3, 16)"),
    ("N/L", "Lichtkegel: Kegel bei eps = pi nur bei theta = pi/4 (Steigung sqrt 2), Fermilinie bei eps = 0 für alle theta"),
    ("O", "native Herkunft von P_link, P_E, P_exp und xi; eine 3D-Anordnung mit reiner Plakettenantwort"),
]
for st, txt in BILANZ:
    say(f"   [{st}] {txt}")
ZAHLEN["bilanz"] = [{"status": st, "aussage": txt} for st, txt in BILANZ]

urteil = ("Yang-Mills mit fester Kopplung ist aus den Schleifen des Transfernetzes nicht nativ; der Weg ist bedingt, mit "
          "gescopten Negativen. Nativ und exakt ist die Plakettenstruktur. Mit dem Link im Transferstrom ist jedes "
          "z120-Gatter eichkovariant. Im 2D-Netz der Vierteldrehungen (vier Teilschichten +x, +y, -x, -y) ist der "
          "Rückkehrblock cos^4 theta + sin^4 theta U_p, jede Plakette wird einmal je Periode umlaufen, und Re Tr(F[1] - F[U]) "
          "= (1/4) sum_p Re Tr(1 - U_p) ist genau die fundamentale Wilson-Funktion. Die Schleifenantwort W(g) = "
          "(d^2 - |Tr g|^2)/(4 d^4) der achtzehnten Analyse ist rekonstruiert: Sie ist die Einperioden-Spurantwort des "
          "elementaren Austauschrings mit dem nativen Gatter -i R(1), je Etikettenzustand normiert. Bei voller Besetzung "
          "erzwingt jeder Umlauf eines Etiketts einen Gegenumlauf, daher |Tr g|^2. W ist damit eine adjungierte "
          "Wilson-Wirkung: zentrumsblind, auf dem Dublett SO(3), quadratisch gleich der fundamentalen mit dem Faktor "
          "1/(2 d^3) und für die Hyperladung der 16 genau 5/24576 phi^2. Auf dem Gitter ist das Austauschnetz aber ein "
          "nicht additives Schleifengas, und in 3D trägt die Vierteldrehung längere Schleifen gleichen Gewichts. Der "
          "orientierungsungerade Teil verschwindet auf der 16 bis zur fünften Ordnung (so(10) ohne kubische Invariante, "
          "dann die Pfaffsche). Nativ ist auch das Verhältnis der Kopplungen: Jede Spurantwort über die vollständige 16 "
          "gibt 1/g_i^2 = lambda k_i mit k_i = 2 Tr_16 T_i^2 = 4 für alle drei Faktoren, also sin^2 theta_W = 3/8 an der "
          "Anpassungsskala; das ist Tr Y^2 : Tr T_3^2 = 10/3 : 2, die SU(5)-Normierung. Nicht nativ ist der Betrag: "
          "lambda = 1/8 (Netzspur je Periode), 1/16384 (W) oder 1 (Niveau mit P_norm) gibt 1/g_U^2 = 1/2, 1/4096, 4. "
          "Keine Lesart trifft alpha^-1 = 137,036 (R3 gibt 128 pi/3 = 134,04 vor jedem Laufen und für eine Generation, "
          "keine Übereinstimmung), keine das z25-Fenster 60 ... 74. Dynamik gibt es nicht: Das Netz vertauscht mit allen "
          "Linkoperatoren, es misst den Fluss, bewegt ihn aber nicht. Je Periode hat es nur räumliche Zyklen, in 1+1D gar "
          "keine; auf zeitabhängige Links reagiert die Materie über Mehrperiodenschleifen, aber die Zeitabhängigkeit erzeugt "
          "das Netz nicht. Eine reflexionspositive Transfermatrix entsteht für kappa > 0 (fundamental "
          "wie adjungiert), aber nur mit den Prämissen P_link (Links mit Haarmaß), P_E (elektrischer Term, kappa_t > 0), "
          "P_exp (Antwort -> Wirkung mit Vorzeichen; für W lesartabhängig) und xi bzw. P_Lorentz (g^2 = 1/(T_R sqrt(kappa_s "
          "kappa_t)); das Netz hat keinen einzigen Lichtkegel). Für K15 heißt das: bleibt L, präzisiert. Die "
          "Plakettenstruktur und das Verhältnis 3/8 sind nativ, die Kopplung ist frei, die Dynamik fehlt. Die Szenarioachse "
          "'Wechselwirkung' erhält den Wert 'Gitter-Yang-Mills aus Transferschleifen: bedingt auf P_link + P_E + P_exp, "
          "Verhältnis nativ, Betrag frei'; mit z25 nur als nackte SO(10)-Gitterkopplung datenverträglich. Nichts ist "
          "befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z135_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n", encoding="utf-8")
(HERE / "ausgabe_z135.txt").write_text("\n".join(OUT) + "\n", encoding="utf-8")
