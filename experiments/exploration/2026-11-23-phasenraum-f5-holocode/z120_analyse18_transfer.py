"""Z120: Prüfung der achtzehnten externen Analyse 'Fortsetzung der fundamentalen Lösung' (25.11.).

Teil 1  Eindeutigkeit aus einer Energieantwort: Zustandssummen sind ganzzahlige Polynome in y = e^-pi; gleiche mittlere
        Modulenergie y Z'/Z heißt Z' Z_0 - Z Z_0' = 0 bei transzendentem y, also Z = Z_0 und mit gleichen
        Einteilchenbesetzungen rho = rho_0. Kontrolle am z94-Gegenbeispiel (0,2485429930 gegen 0,0903035897).
Teil 2  Globale Mittelung über N Nähte: Gewichte w_(N,k) aus der Transfergruppe; reduzierte Paarkorrelation von zwei
        gegen drei Nähte: (1 - r)/(5 + r) gegen (1 - r)/(5 + 3 r), r = e^-4pi, keine konsistente Familie.
Teil 3  Lokales kohärentes Transfergatter U = B^+ e^(-i pi H_Q/2) B, B = e^(-i pi J_x/2): auf dem physikalischen
        Unterraum U = P_Triplett - i P_Singulett (inverse Wurzel aus SWAP), Ordnung 4, verschränkend, nicht gaußsch.
Teil 4  Alternierende Transfers: cos eps(k) = (1 + cos k)/2, eps ~ |k|/sqrt 2; Wilson-Fortsetzung A_q(k) = (S^+ - w(k))
        (1 - q S^4): ein Weylpunkt am Rand, sieben Reglermoden am fernen Rand, Volumenlücke >= 1/sqrt 2.
Teil 5  Schleifenantworten: W(g) = (d^2 - |Tr g|^2)/(4 d^4) ~ Tr A^2/(4 d^3); Hyperladung auf der 16: 5/24576 phi^2;
        orientierte Antwort d_s d_t W_Gamma = -Tr(Gamma A B)/(4 d^3), im Spinor eps-artig.
"""
import itertools
import json
from collections import Counter
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(120)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


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


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Eindeutigkeit aus einer Energieantwort")
z94 = json.loads((HERE / "z94_ergebnis.json").read_text())
fz = z94["zahlen"]["aussenbehauptung_gleiche_saaten"]["f(n) für N = 6"]
C6 = list(itertools.product((0, 1), repeat=6))
H6 = {c: sum(c) + 4*fz.get("".join(map(str, c)), 0) for c in C6}
mp.mp.dps = 50
yv = mp.exp(-mp.pi)
Zk = mp.fsum(yv**H6[c] for c in C6)
Ek = mp.fsum(H6[c]*yv**H6[c] for c in C6)/Zk
Eg = 6*yv/(1 + yv)
rand_k = [mp.fsum(yv**H6[c] for c in C6 if c[k])/Zk for k in range(6)]
y = sp.symbols("y")
Zs = sum(y**H6[c] for c in C6)
Z0 = (1 + y)**6
Pd = sp.expand(sp.diff(Zs, y)*Z0 - Zs*sp.diff(Z0, y))
koeff_int = all(c_.is_integer for c_ in sp.Poly(Pd, y).coeffs())
ZAHLEN["energie"] = {"<E> Gauß": mp.nstr(Eg, 12), "<E> kollektiv": mp.nstr(Ek, 12),
                     "Grad des Differenzpolynoms": sp.degree(Pd, y)}
check(f"z94-Zustand aus sechs Saaten: Einzelbesetzungen gleich (max Abweichung {float(max(abs(r - yv/(1 + yv)) for r in rand_k)):.0e}), "
      f"aber mittlere Modulenergie {mp.nstr(Ek, 10)} gegen {mp.nstr(Eg, 10)} (Gauß); das Differenzpolynom Z' Z_0 - Z Z_0' hat "
      f"ganzzahlige Koeffizienten ({koeff_int}, Grad {sp.degree(Pd, y)}) und ist nicht null. Satz: Stimmt die Energie bei "
      "y = e^-pi (transzendent, Gelfond) überein, so verschwindet dieses ganzzahlige Polynom identisch, (Z/Z_0)' = 0, Z = Z_0 "
      "(gleiche Zustandszahl), und mit gleichen Einteilchenbesetzungen gibt die relative Entropie rho = rho_0",
      float(max(abs(r - yv/(1 + yv)) for r in rand_k)) < 1e-40 and abs(Eg - mp.mpf("0.2485429930")) < 1e-10
      and abs(Ek - mp.mpf("0.0903035897")) < 1e-10 and koeff_int and Pd != 0)
befund("C: Mit einer einzigen zusätzlichen skalaren Antwort, der mittleren Modulenergie, wird S(n) für jedes n trivial "
       "(Transzendenz statt Kombinatorik; vgl. z107: Starrheit ist arithmetisch). Bedingung: Die native Spannungsantwort "
       "(A*, <T> = 0, z107) muss die tatsächliche Modulenergie festlegen; im unendlichen System braucht es einen "
       "kontrollierten Grenzwert", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Globale Mittelung über N Nähte")
r = np.exp(-4*np.pi)


def nahtmodell(N):
    c = jw(2*N)
    num = [x.T @ x for x in c]
    Q = [num[2*i] + num[2*i + 1] for i in range(N)]
    Ntot = sum(Q)
    Sx = sum(c[2*i].T @ c[2*i + 1] + c[2*i + 1].T @ c[2*i] for i in range(N))/2
    Sy = sum(-1j*(c[2*i].T @ c[2*i + 1] - c[2*i + 1].T @ c[2*i]) for i in range(N))/2
    Sz = sum(num[2*i] - num[2*i + 1] for i in range(N))/2
    S2 = Sx @ Sx + Sy @ Sy + Sz @ Sz
    sek = np.where(np.isclose(np.diag(Ntot), N))[0]
    phys = [s_ for s_ in sek if all(np.isclose(Q[i][s_, s_], 1) for i in range(N))]
    dE = np.array([sum(Q[i][s_, s_]**2 for i in range(N))/2 - N/2 for s_ in range(2**(2*N))])
    S2s = S2[np.ix_(sek, sek)]
    ev, V = np.linalg.eigh(S2s)
    gew = {}
    for J2 in sorted(set(np.round(ev, 6))):
        P = V[:, np.isclose(ev, J2)]
        D = np.diag(np.exp(-4*np.pi*dE[sek]))
        gew[J2] = np.trace(P.conj().T @ D @ P).real/P.shape[1]
    S2p = S2[np.ix_(phys, phys)]
    evp, Vp = np.linalg.eigh(S2p)
    rho = sum(gew[J2]*Vp[:, np.isclose(evp, J2)] @ Vp[:, np.isclose(evp, J2)].conj().T for J2 in gew if np.any(np.isclose(evp, J2)))
    rho = rho/np.trace(rho)
    return rho, gew, phys, c


def paarkorrelation(rho, N):
    # Physikalischer Raum: je Naht ein Spin; Basis aus den Fockzuständen, Spin der Naht i = n_(2i) - n_(2i+1)
    dim = rho.shape[0]
    return rho, dim


rho2, gew2, phys2, c2m = nahtmodell(2)
rho3, gew3, phys3, c3m = nahtmodell(3)


def sz_naht(c, phys, i):
    n_up, n_dn = c[2*i].T @ c[2*i], c[2*i + 1].T @ c[2*i + 1]
    return ((n_up - n_dn)/2)[np.ix_(phys, phys)]


def sx_naht(c, phys, i):
    return ((c[2*i].T @ c[2*i + 1] + c[2*i + 1].T @ c[2*i])/2)[np.ix_(phys, phys)]


C2 = 4*np.trace(rho2 @ sz_naht(c2m, phys2, 0) @ sz_naht(c2m, phys2, 1)).real
C3 = 4*np.trace(rho3 @ sz_naht(c3m, phys3, 0) @ sz_naht(c3m, phys3, 1)).real
C2x = 4*np.trace(rho2 @ sx_naht(c2m, phys2, 0) @ sx_naht(c2m, phys2, 1)).real
w_formel = lambda N, k: sum(sp.binomial(k, d)*sp.binomial(N - k + 1, d)*sp.exp(-4*sp.pi*d) for d in range(k + 1))/sp.binomial(N + 1, k)
w3 = sorted(gew3.values())
w3_formel = sorted(float(w_formel(3, k)) for k in (0, 1))
ZAHLEN["mittelung"] = {"Gewichte N = 2": {str(k): v for k, v in gew2.items()}, "Gewichte N = 3": {str(k): v for k, v in gew3.items()},
                       "C2": C2, "C3": C3, "C2 Formel": (1 - r)/(5 + r), "C3 Formel": (1 - r)/(5 + 3*r)}
check(f"N Nähte mit je zwei Moden, je eine Anregung, Mittelung über die Transfergruppe (Schur: Spur/Dimension je Spinsektor): "
      f"Gewichte N = 3 {[round(float(x), 10) for x in w3]} = Formel w_(N,k) {[round(x, 10) for x in w3_formel]}; Paarkorrelation "
      f"<sigma_A sigma_B> zweier Nähte allein {C2:.10f} = (1 - r)/(5 + r) (isotrop: x-Komponente {C2x:.10f}), als Teil von drei "
      f"Nähten {C3:.10f} = (1 - r)/(5 + 3 r): die reduzierten Zustände bilden keine konsistente Familie (Differenz "
      f"{C2 - C3:.2e} ~ r)",
      np.allclose(w3, w3_formel, atol=1e-12) and abs(C2 - (1 - r)/(5 + r)) < 1e-12 and abs(C3 - (1 - r)/(5 + 3*r)) < 1e-12
      and abs(C2 - C2x) < 1e-12 and C2 - C3 > 0)
befund("E: Die Korrektur der Analyse 18 stimmt: Die globale Gruppenmittelung gibt einen effektiven Generator, aber keine "
       "konsistente Familie lokaler Zustände; bei r = 0 ist sie eine Mischung ungekoppelter Produkte. Lokale Wechselwirkung "
       "verlangt eine geordnete, kohärente Zusammensetzung", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Lokales kohärentes Transfergatter")
c = jw(4)
num = [x.T @ x for x in c]
QA, QB = num[0] + num[1], num[2] + num[3]
Jp = c[0].T @ c[2] + c[1].T @ c[3]
Jx = (Jp + Jp.T)/2
HQ = (QA @ QA + QB @ QB)/2
Bm = expm(-1j*np.pi*Jx/2)
U = Bm.conj().T @ expm(-1j*np.pi*HQ/2) @ Bm
phys = [s_ for s_ in range(16) if np.isclose(QA[s_, s_], 1) and np.isclose(QB[s_, s_], 1)]
Up = U[np.ix_(phys, phys)]
leck = np.linalg.norm(U[:, phys] - np.eye(16)[:, phys] @ Up)
Sx = (c[0].T @ c[1] + c[1].T @ c[0] + c[2].T @ c[3] + c[3].T @ c[2])/2
Sy = (-1j*(c[0].T @ c[1] - c[1].T @ c[0]) - 1j*(c[2].T @ c[3] - c[3].T @ c[2]))/2
Sz = (num[0] - num[1] + num[2] - num[3])/2
S2 = (Sx @ Sx + Sy @ Sy + Sz @ Sz)[np.ix_(phys, phys)]
Pt = S2/2
Ps = np.eye(4) - Pt
ziel = Pt - 1j*Ps
gleich_bis_phase = min(np.linalg.norm(Up - ph*ziel) for ph in (1, -1, 1j, -1j))
ordnung4 = np.allclose(np.linalg.matrix_power(Up, 4), np.linalg.matrix_power(Up, 4)[0, 0]*np.eye(4))
wurzel_swap = np.allclose(ziel @ ziel, Pt - Ps)
cr = [x.T for x in c]
vak = np.eye(16)[:, 0]
u1 = np.array([[vak @ c[k_] @ U @ cr[i] @ vak for i in range(4)] for k_ in range(4)])
abw = 0.0
for i, j in itertools.combinations(range(4), 2):
    psi = cr[i] @ cr[j] @ vak
    a_i = sum(u1[k_, i]*cr[k_] for k_ in range(4))
    a_j = sum(u1[k_, j]*cr[k_] for k_ in range(4))
    abw = max(abw, np.linalg.norm(U @ psi - a_i @ a_j @ vak))
schmidt = np.linalg.svd(Up[:, phys.index(9)].reshape(2, 2), compute_uv=False)
ZAHLEN["gatter"] = {"Leckage aus dem physikalischen Raum": float(leck), "Abstand zu P_T - i P_S bis auf Phase": float(gleich_bis_phase),
                    "Einteilchenblock": str(np.round(u1, 6).tolist()), "Abweichung von der gaußschen Fortsetzung": float(abw),
                    "Schmidt-Werte von U|A hoch, B runter>": schmidt.tolist()}
check(f"zwei Nähte (je zwei Moden), H_Q = (Q_A^2 + Q_B^2)/2, J_x = (J^+ + J^-)/2 aus den Transferströmen: U = B^+ e^(-i pi H_Q/2) B "
      f"lässt den Raum mit je einer Anregung invariant (Leckage {leck:.0e}) und ist dort bis auf eine globale Phase P_Triplett - "
      f"i P_Singulett (Abstand {gleich_bis_phase:.0e}), eine Wurzel aus SWAP ({wurzel_swap}), Ordnung 4 ({ordnung4}), maximal "
      f"verschränkend (Schmidt-Werte {np.round(schmidt, 6).tolist()}). Nicht gaußsch: Der Einteilchenblock ist die reine Phase "
      f"e^(-i pi/4) (Abweichung {np.linalg.norm(u1 - np.exp(-1j*np.pi/4)*np.eye(4)):.0e}), die gaußsche Fortsetzung auf zwei "
      f"Teilchen weicht um {abw:.4f} ab (= |1 - i|/sqrt 2: Singulettphase -1 statt -i)",
      leck < 1e-12 and gleich_bis_phase < 1e-12 and wurzel_swap and ordnung4 and np.allclose(schmidt, [2**-0.5]*2)
      and np.linalg.norm(u1 - np.exp(-1j*np.pi/4)*np.eye(4)) < 1e-12 and abw > 0.5)
befund("C: Das Gatter ist aus Ladungsenergie Q^2/2 und Transferstrom mit Viertelwinkeln gebaut. Anders als die globale "
       "Mittelung erhält es die Reihenfolge der Ereignisse; ob seine Winkel frei sind, prüft Teil 3b", True)

say("\nTeil 3b: Die Winkel des Gatters sind nicht frei")
Pphys = np.zeros((16, 16))
Pphys[phys, phys] = 1


def gatter(th, ph, H=HQ):
    Bt = expm(-1j*th*Jx)
    return Bt.conj().T @ expm(-1j*ph*H) @ Bt


thetas = np.linspace(0, np.pi, 181)
leckagen = [np.linalg.norm((np.eye(16) - Pphys) @ gatter(t, np.pi/2) @ Pphys) for t in thetas]
null_th = [round(float(t/np.pi), 4) for t, l in zip(thetas, leckagen) if l < 1e-10]
trivial_th = [round(float(t/np.pi), 4) for t in thetas if np.linalg.norm((np.eye(16) - Pphys) @ gatter(t, np.pi/2) @ Pphys) < 1e-10
              and np.linalg.norm(gatter(t, np.pi/2)[np.ix_(phys, phys)] - gatter(t, np.pi/2)[phys[0], phys[0]]*np.eye(4)) < 1e-10]
ZAHLEN["winkel"] = {"leckfreie Drehwinkel theta/pi": null_th, "davon trivial": trivial_th,
                    "maximale Leckage": float(max(leckagen))}
check(f"Transferdrehung B_theta = e^(-i theta J_x) bei Quadratphase pi/2: Der physikalische Raum (eine Anregung je Naht) bleibt "
      f"nur für theta/pi = {null_th} invariant (sonst Leckage bis {max(leckagen):.3f} in Doppelbesetzungen); bei {trivial_th} ist "
      "das Gatter eine reine Phase. Der einzige nichttriviale leckfreie Winkel ist die Vierteldrehung theta = pi/2, also "
      "A -> B -> -A innerlich: genau P_sym (z94)",
      null_th == [0.0, 0.5, 1.0] and trivial_th == [0.0, 1.0])
QAm, QBm = np.diag(QA), np.diag(QB)
N2 = [s_ for s_ in range(16) if np.isclose(QAm[s_] + QBm[s_], 2)]
Rd = np.array([1j**int(round(QAm[s_] - QBm[s_])) for s_ in N2])
R_inv_wurzel = np.diag([1 if np.isclose(r_, 1) else -1j for r_ in Rd])
T_N2 = expm(-1j*np.pi/2*HQ)[np.ix_(N2, N2)]
wurzel_uhr = np.linalg.norm(T_N2 - (-1j)*R_inv_wurzel)
H_wechsel = np.diag(QAm*QBm)
wechsel_N2 = np.linalg.norm(expm(1j*np.pi/2*H_wechsel)[np.ix_(N2, N2)] + T_N2)
U_wechsel = gatter(np.pi/2, -np.pi/2, H_wechsel)[np.ix_(phys, phys)]
abst_wechsel = min(np.linalg.norm(U_wechsel - ph*ziel) for ph in np.exp(1j*np.pi*np.arange(8)/4))
H_frei = sum(n_ @ n_ for n_ in num)/2
U_frei = gatter(np.pi/2, np.pi/2, H_frei)[np.ix_(phys, phys)]
frei_trivial = np.linalg.norm(U_frei - U_frei[0, 0]*np.eye(4))
phis = np.linspace(0, np.pi, 9)
rel_phase = [complex(np.round((gatter(np.pi/2, p)[np.ix_(phys, phys)] @ np.linalg.eigh(Ps)[1][:, -1]
                               ).conj() @ np.linalg.eigh(Ps)[1][:, -1], 6)) for p in phis]
ZAHLEN["wurzel_der_uhr"] = {"|e^(-i pi H_Q/2) + i R^(-1/2)| auf N = 2": float(wurzel_uhr),
                            "|i^(Q_A Q_B) + e^(-i pi H_Q/2)| auf N = 2": float(wechsel_N2),
                            "Abstand Z4-Wechselgatter zu P_T - i P_S": float(abst_wechsel),
                            "freie Energie: Abweichung von reiner Phase": float(frei_trivial)}
check(f"Auf dem Sektor N = 2 ist die Quadratphase die inverse Wurzel der Uhr: e^(-i pi H_Q/2) = -i R^(-1/2) mit R = "
      f"i^(Q_A - Q_B) (Abweichung {wurzel_uhr:.0e}), also U = -i (B^+ R B)^(-1/2): das Gatter ist die Wurzel der um eine "
      f"Vierteldrehung gedrehten Uhr. Gleichwertig ist die bilineare Z4-Wechselphase zweier Nähte: Auf N = 2 gilt i^(Q_A Q_B) "
      f"= -e^(-i pi H_Q/2) (Abweichung {wechsel_N2:.0e}), also U ~ B^+ i^(Q_A Q_B) B = P_T - i P_S (Abstand {abst_wechsel:.0e}). "
      f"Kontrolle: Mit der freien Energie sum n_s^2/2 (gaußsch, ohne Doppelbesetzungsterm) ist das Gatter eine reine Phase "
      f"(Abweichung {frei_trivial:.0e}); die Wechselwirkung sitzt im spinblinden Term Q^2 = n_hoch + n_runter + 2 n_hoch n_runter, "
      f"gleichwertig im Kreuzterm Q_A Q_B",
      wurzel_uhr < 1e-12 and wechsel_N2 < 1e-12 and abst_wechsel < 1e-12 and frei_trivial < 1e-12)
befund("C: Die Winkel des Gatters sind keine neuen Zahlen. Die Drehung ist durch Leckfreiheit auf die Vierteldrehung "
       "festgelegt, und das ist die Ursprungsprämisse P_sym. Die Phase ist die Wurzel der nativen Uhr (Ordnung 4 wie R). "
       "Offen bleiben drei Punkte: ob die spinblinde Ladungsenergie (nicht die freie) die native Antwort der Naht ist, "
       "ob die Uhr bei der Komposition als Wurzel wirkt, und das Vorzeichen der Wurzel", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Alternierende Transfers, lineare Dispersion, Wilson-Fortsetzung")
ks = np.linspace(-np.pi, np.pi, 401)
th = np.pi/4
Ub = np.array([[np.cos(th), -1j*np.sin(th)], [-1j*np.sin(th), np.cos(th)]])
Ub2 = np.array([[np.cos(th), 1j*np.sin(th)], [1j*np.sin(th), np.cos(th)]])
abw_disp = 0.0
for k in ks:
    T = np.array([[1, 0], [0, np.exp(1j*k)]])
    UF = (T.conj() @ Ub2 @ T) @ Ub
    ph = np.angle(np.linalg.eigvals(UF))
    abw_disp = max(abw_disp, abs(np.cos(ph[0]) - (1 + np.cos(k))/2), abs(np.cos(ph[1]) - (1 + np.cos(k))/2))
kk = 1e-3
eps_klein = np.arccos((1 + np.cos(kk))/2)
qv = np.exp(-2*np.pi)
Nh = 65
S = np.zeros((Nh, Nh))
for n in range(Nh - 1):
    S[n + 1, n] = 1
ecken = list(itertools.product((0, np.pi), repeat=3))
rand_moden, gehoben, sv_ecken = 0, 0, []
for kv in ecken:
    w = sum(1 - np.cos(x) for x in kv)
    A = (S.T - w*np.eye(Nh)) @ (np.eye(Nh) - qv*np.linalg.matrix_power(S, 4))
    uu, sv, vh = np.linalg.svd(A)
    sv_ecken.append(float(sv[-1]))
    if sv[-1] < 1e-8 and np.sum(np.abs(vh[-1][:8])**2) > 0.99 and np.sum(np.abs(uu[:, -1][-8:])**2) > 0.99:
        rand_moden += 1
    elif sv[-1] > (w - 1)*(1 - qv) - 1e-9 and w > 1:
        gehoben += 1
nullmodus = np.array([qv**(n//4) if n % 4 == 0 else 0.0 for n in range(Nh)])
nullmodus /= np.linalg.norm(nullmodus)
A0 = S.T @ (np.eye(Nh) - qv*np.linalg.matrix_power(S, 4))
gitter = np.linspace(-np.pi, np.pi, 41)
K1, K2, K3 = np.meshgrid(gitter, gitter, gitter, indexing="ij")
wg = 3 - np.cos(K1) - np.cos(K2) - np.cos(K3)
h2 = (np.sin(K1)**2 + np.sin(K2)**2 + np.sin(K3)**2)/2
luecke = float(np.sqrt(h2 + ((1 - qv)*(1 - wg))**2).min())
ZAHLEN["dispersion_wilson"] = {"eps(k)/k bei k = 1e-3": eps_klein/kk, "kleinste Singulärwerte an den Ecken": sv_ecken,
                               "Volumenlücke": luecke}
check(f"alternierende Transfergatter (Mischwinkel pi/4, entgegengesetzte Orientierung): cos eps(k) = (1 + cos k)/2 exakt "
      f"(Abweichung {abw_disp:.0e}), eps(k)/k -> {eps_klein/kk:.6f} = 1/sqrt 2: lineare Dispersion. Wilson-Fortsetzung A_q(k) "
      f"= (S^+ - w(k))(1 - q S^4) mit w = sum(1 - cos k_j), 65 Hardy-Zustände: an den acht Ecken genau {rand_moden} Nullmodus "
      f"(k = 0, am Rand, = sqrt(1 - q^2) sum q^n |4n>, Residuum {np.linalg.norm(A0 @ nullmodus):.0e}; sein Kokern sitzt am "
      f"Abschneiderand), an den {gehoben} anderen Ecken (w = 2, 4, 6) ist der kleinste Singulärwert >= (w - 1)(1 - q): "
      f"{[round(x, 4) for x in sv_ecken]}; Volumenlücke min sqrt(|h|^2 + (1 - q)^2 (1 - w)^2) = {luecke:.6f} = 1/sqrt 2 "
      f"(erreicht bei w = 1, |h|^2 = 1/2; Gitter 41^3)",
      abw_disp < 1e-12 and abs(eps_klein/kk - 1/np.sqrt(2)) < 1e-4 and rand_moden == 1 and gehoben == 7
      and np.linalg.norm(A0 @ nullmodus) < 1e-12 and abs(luecke - 1/np.sqrt(2)) < 1e-9)
befund("C/L: Lineare Dispersion und Wilson-Term folgen aus der Impulsantwort des Transfergatters, und genau ein Weylpunkt "
       "bleibt am ursprünglichen Rand (die sieben anderen Ecken haben Index 0). Die dreidimensionale Anordnung der Gatter ist "
       "eine Konstruktion; native Entstehung und wechselwirkender Kontinuumslimes sind offen", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Schleifenantworten: Yang-Mills- und orientierte Kontraktion")
W = lambda g: (g.shape[0]**2 - abs(np.trace(g))**2)/(4*g.shape[0]**4)
ok_ym = True
for d in (2, 4, 16):
    A = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
    A = (A + A.conj().T)/2
    A -= np.trace(A)/d*np.eye(d)
    e_ = 1e-4
    ok_ym &= abs(W(expm(1j*e_*A))/(e_**2*np.trace(A @ A).real/(4*d**3)) - 1) < 1e-6
y_lad = [sp.Rational(-1, 3)]*3 + [sp.Rational(1, 2)]*2
Yv = [float(sum(y_lad[i] for i in S_)) for k in (0, 2, 4) for S_ in itertools.combinations(range(5), k)]
Yd = np.diag(Yv)
phi = 1e-3
ym16 = W(expm(1j*phi*Yd))/phi**2
s0, sx_, sy_, sz_ = np.eye(2), np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0]).astype(complex)
g = [np.kron(sx_, s0), np.kron(sy_, s0), np.kron(sz_, sx_), np.kron(sz_, sy_)]
G5 = g[0] @ g[1] @ g[2] @ g[3]
Gam = G5 if np.allclose(G5 @ G5, np.eye(4)) else 1j*G5
Aa, Bb = g[0] @ g[1]/2*1j, g[2] @ g[3]/2*1j
WG = lambda gm: np.trace(Gam @ gm)*np.trace(np.linalg.inv(gm))/(4*4**4)
h_ = 1e-3
gst = lambda s_, t_: WG(expm(1j*s_*Aa) @ expm(1j*t_*Bb))
mix = (gst(h_, h_) - gst(h_, -h_) - gst(-h_, h_) + gst(-h_, -h_))/(4*h_**2)
erwartet = -np.trace(Gam @ Aa @ Bb)/(4*4**3)
ok_or = abs(mix - erwartet) < 1e-6*max(1, abs(erwartet)) and abs(erwartet) > 1e-6
sym = abs(np.trace(Gam @ Aa @ Aa)) < 1e-12
check(f"W(g) = (d^2 - |Tr g|^2)/(4 d^4) gibt Tr A^2/(4 d^3) + O(A^4) für d = 2, 4, 16 ({ok_ym}); Hyperladung auf der 16 "
      f"(Tr Y^2 = 10/3): W(e^(i phi Y))/phi^2 = {ym16:.10f} = 5/24576 = {5/24576:.10f}; orientierte Antwort "
      f"W_Gamma = Tr(Gamma g) Tr(g^-1)/(4 d^4): gemischte Ableitung {np.round(mix, 8)} = -Tr(Gamma A B)/(4 d^3) = "
      f"{np.round(erwartet, 8)} für A ~ gamma_01, B ~ gamma_23 (eps-artig; Tr(Gamma A A) = 0: {sym})",
      ok_ym and abs(ym16 - 5/24576) < 1e-9 and ok_or and sym)
befund("E/C: Die Schleifenidentitäten stimmen als Algebra: dieselbe Spurantwort trägt die Yang-Mills-Kontraktion Tr A^2 und, mit "
       "der Orientierung Gamma, die antisymmetrische Kontraktion Tr(Gamma A B). Dass diese Schleifen das native Funktional "
       "bilden, ihre Normierung im Kontinuum und G sind offen", True)

urteil = ("Die achtzehnte externe Analyse hält in allen gerechneten Punkten. (1) Eine zusätzliche skalare Antwort, die mittlere "
          "Modulenergie, macht die Starrheit für jedes n trivial: ganzzahlige Zustandssummen in y = e^-pi, Transzendenz, Z = Z_0, "
          "rho = rho_0; am z94-Gegenbeispiel unterscheidet die Energie (0,0903 gegen 0,2485) bei gleichen Einzelbesetzungen. "
          "Bedingung: A* muss die Modulenergie festlegen. (2) Die globale Mittelung über N Nähte gibt die Gewichte w_(N,k), aber "
          "keine konsistente Familie: Paarkorrelation (1 - r)/(5 + r) gegen (1 - r)/(5 + 3 r). (3) Das lokale Gatter U = B^+ "
          "e^(-i pi H_Q/2) B ist auf dem physikalischen Raum P_T - i P_S (eine Wurzel aus SWAP), Ordnung 4, maximal "
          "verschränkend, nicht gaußsch. Seine Winkel sind nicht frei: Leckfreiheit erzwingt die Vierteldrehung (= P_sym, "
          "z94), und die Quadratphase ist auf N = 2 die inverse Wurzel der Uhr, U = -i (B^+ R B)^(-1/2), gleichwertig die "
          "bilineare Z4-Wechselphase i^(Q_A Q_B) zweier Nähte. Mit der freien Energie wäre das Gatter trivial; die "
          "Wechselwirkung sitzt im spinblinden Ladungsterm bzw. im Kreuzterm Q_A Q_B. (4) Alternierende Gatter geben "
          "cos eps = (1 + cos k)/2, also "
          "lineare Dispersion; die Wilson-Fortsetzung A_q(k) = (S^+ - w)(1 - q S^4) lässt genau einen Weylpunkt am Rand, "
          "Volumenlücke >= 1/sqrt 2. (5) Die Schleifenantwort trägt Tr A^2 (Yang-Mills, 5/24576 phi^2 für die Hyperladung der "
          "16) und mit Gamma die orientierte Kontraktion Tr(Gamma A B). Offen: native Kompositionsregel (Winkel), "
          "Kontinuumslimes, Identifikation der Schleifen mit dem Funktional, G und Lambda.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z120_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
