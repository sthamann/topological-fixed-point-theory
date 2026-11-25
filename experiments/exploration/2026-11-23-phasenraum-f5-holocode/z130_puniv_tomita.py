"""Z130: Ist P_univ strukturell? Native Zwischennahtoperatoren, Einzeltransfer Naht -> Tomita-Partner, Folgen für den Ursprung.

Folgeprüfung zu z125 (P_sym aus nichttrivialer leckfreier Komposition unter der Superauswahlprämisse P_univ). Leitlinie
'Rekonstruktion statt Auswahl'; H6 und die gelösten Bausteine bleiben unberührt, nichts wird befördert. Statuskürzel:
E exakt, C bedingt, L Lesart, O offen.
Bestand: z125 (Leckformel |cos sin||1 - e^(-i phi)|, Z4-Superauswahltabelle, Tomita-Paaruhr Delta^(-i/2) ~ i^(Q_A - Q_B) des
Original-Partner-Paars, P_univ), z99 (einzige lokale Kopplung zweier Nähte: su(2)_1-Strom psi*_A psi_B), z105 (Doppeltransfer,
Tomita-kombiniert, uhrinvariant), z106 (gemeinsame CAR = Super-Tensorprodukt, Gamma(pi) in der U(2) der Ströme, gerade Daten
wählen die Komposition nicht), z113 (Modulinvarianz der nichtchiralen Naht erzwingt R = R-quer und (-1)^(F+F-quer) = +1;
e^(2 pi i L_0) = R), z116/z119 (virtuelle Transfers, Twirl, Q^2/2 bricht die Transfer-SU(2)), z120 (Gatter
U = B^+ e^(-i pi H_Q/2) B), z122 (R_phys = R R-quer^-1, J R_phys J = R_phys, a_c Omega = q b_c^+ Omega), z93 (Ladungssektor
ist Teil des Zustands), z94/z107 (P_clock: Delta^(-i/2) = R folgt aus A*; A* allein gibt P_sym nicht).

Teil 1  Klassifikation der nativen Zwischennahtoperatoren nach (-1)^Q_A, (-1)^Q_B und nach Paar- und Produktuhr: Ströme E^+-_n,
        Cartan H_n, Vertex psi*_A psi_B, Kreuzpaar psi*_A psi*_B, Doppeltransfer C^+, virtuelle Transfers (Superaustausch, Twirl),
        Sandwich U_theta, Transferdrehung B_theta, kanonischer Tausch Gamma(pi); J-Bilder im verdoppelten Modell (zwei Nähte mit
        je eigenem Partner, 256 Zustände).
Teil 2  Standardform einer Naht mit ihrem Tomita-Partner (GNS des KMS-Zustands bei beta = 4 pi; S aus der Definition, Delta, J):
        Modulgruppe, J, M' = JMJ und gradierter Kommutant, KMS-Relationen, Einzeltransfer und Kreuzpaar, 2 pi-Drehung und
        Niveaugleichheit (z113) als physikalischer Raum, Gleichwertigkeit mit P_univ, Besetzung der Sektoren.
Teil 3  Folgen: reduzierter Ursprung, Verträglichkeit mit z107 und z94, Lesart zweier gleichchiraler Nähte, K19.
Teil 4  Bilanz, Prämissen, Urteil.
"""
import json
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.stats import unitary_group

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(130)
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


def nrm(X):
    return float(np.linalg.norm(X))


def kom(X, Y):
    return nrm(X @ Y - Y @ X)


def klein(x, schwelle=1e-12):
    """Ausgabe kleiner Reste ohne Rundungsrauschen."""
    return f"< {schwelle:.0e}" if abs(x) < schwelle else f"{x:.2e}"


def zahl(x, st=6, null=1e-12):
    """JSON-Wert: Reste unter der Schwelle null als 0, sonst auf st signifikante Stellen."""
    return 0.0 if abs(x) < null else float(f"{x:.{st}g}")


def j_marken(th):
    """wie z94/z125: j-Invariante der vier Marken {1, -1, e^(i th), -e^(i th)}."""
    z = [1, -1, np.exp(1j*th), -np.exp(1j*th)]
    lam = (z[0] - z[2])*(z[1] - z[3])/((z[0] - z[3])*(z[1] - z[2]))
    return float((256*(1 - lam + lam**2)**3/(lam**2*(1 - lam)**2)).real)


# ------------------------------------------------------------------ gemeinsames Modell (16 Zustände)
# Moden in Jordan-Wigner-Reihenfolge A0, A1, B0, B1 (z125: zwei Seiten mit je zwei Moden). In Teil 1 ist B eine zweite Naht,
# in Teil 2 der Tomita-Partner von A (Purifikation); die Operatoren sind dieselben, verschieden ist die Lesart.
C4 = jw(4)
NUM = [np.diag(x.T @ x).round().astype(int) for x in C4]
QA, QB = NUM[0] + NUM[1], NUM[2] + NUM[3]
NT = QA + QB
GA, GB = np.diag((-1.0)**QA), np.diag((-1.0)**QB)
GAM = GA @ GB
R_PAAR = np.diag(1j**((QA - QB) % 4))                     # i^(Q_A - Q_B)
R_PROD = np.diag(1j**((QA + QB) % 4))                     # i^(Q_A + Q_B)
JP = C4[0].T @ C4[2] + C4[1].T @ C4[3]                    # J^+ = E^+_0 (B -> A)
JM = JP.T
JX, JY = (JP + JM)/2, (JP - JM)/(2j)
JZ = np.diag((QA - QB)/2.0)
NOP = np.diag(NT.astype(float))
HQ = np.diag((QA**2 + QB**2)/2.0)
I16 = np.eye(16)
GAMMA_PI = expm(-1j*np.pi/2*NOP) @ expm(1j*np.pi/2*(JP + JM))   # kanonischer Tausch (z106, z119)
TH181 = np.linspace(0, np.pi, 181)


def gatter(th, phi=np.pi/2):
    """z120/z125: U = B_theta^+ e^(-i phi H_Q) B_theta mit B_theta = e^(-i theta J_x)."""
    B = expm(-1j*th*JX)
    return B.conj().T @ expm(-1j*phi*HQ) @ B


def paritaetsklassen(X, ga=GA, gb=GB):
    """Anteile |X_(sA,sB)|^2/|X|^2 der Klassen ga X ga = sA X, gb X gb = sB X."""
    out = {}
    for sA in (1, -1):
        for sB in (1, -1):
            Y = (X + sA*ga @ X @ ga + sB*gb @ X @ gb + sA*sB*(ga @ gb) @ X @ (gb @ ga))/4
            out[(sA, sB)] = nrm(Y)**2/nrm(X)**2
    return out


def uhrklassen(X, Ru):
    """Anteile der Klassen k = 0..3 mit Ru X Ru^-1 = i^k X."""
    out = {}
    for k in range(4):
        Y = sum((1j)**(-k*m)*np.linalg.matrix_power(Ru, m) @ X @ np.linalg.matrix_power(Ru.conj().T, m) for m in range(4))/4
        out[k] = nrm(Y)**2/nrm(X)**2
    return out


def klassen_text(fr, tol=1e-12):
    reine = [k for k, v in fr.items() if v > tol]
    lab = {(1, 1): "(+,+)", (1, -1): "(+,-)", (-1, 1): "(-,+)", (-1, -1): "(-,-)"}
    if len(reine) == 1:
        return f"rein {lab[reine[0]]}" if isinstance(reine[0], tuple) else f"rein k = {reine[0]}"
    return " + ".join(f"{lab[k] if isinstance(k, tuple) else 'k = ' + str(k)} {fr[k]:.3f}" for k in reine)


def kommutant_projektor(gens):
    """HS-Orthogonalprojektor auf den Kommutanten {gens}' (Zeilen-vec: vec(GX - XG) = (G (x) 1 - 1 (x) G^T) vec X); für die
    unitäre Darstellung einer kompakten Gruppe ist das genau das Haar-Mittel (Twirl)."""
    d = gens[0].shape[0]
    M = np.vstack([np.kron(G, np.eye(d)) - np.kron(np.eye(d), G.T) for G in gens])
    _, s, Vh = np.linalg.svd(M)
    rang = int(np.sum(s > 1e-10))
    V = Vh[rang:].conj().T
    return V @ V.conj().T, V.shape[1]


def blockpartition(diags):
    """gemeinsame Eigenräume einer Familie diagonaler Operatoren (deterministisch sortiert)."""
    schluessel = {}
    for i in range(len(diags[0])):
        k = tuple((round(float(np.real(d[i])), 9), round(float(np.imag(d[i])), 9)) for d in diags)
        schluessel.setdefault(k, []).append(i)
    return sorted(tuple(v) for v in schluessel.values())


# =================================================================== Teil 1
say("Teil 1: Klassifikation der nativen Zwischennahtoperatoren")
# 1.1 Operatorliste
P11 = np.diag(((QA == 1) & (QB == 1)).astype(float))
PD = np.diag((((QA == 2) & (QB == 0)) | ((QA == 0) & (QB == 2))).astype(float))
HOP = JP + JM
K_VIRT = -P11 @ HOP @ PD @ HOP @ P11                       # zweite Ordnung (H_Q - E_0 = 1 auf (2,0), (0,2)), z119 (a-iv)
P_U2, DIM_U2 = kommutant_projektor([JX, JY, JZ, NOP])


def twirl(X, P=P_U2):
    return (P @ X.reshape(-1)).reshape(X.shape)


BOLTZ = expm(-4*np.pi*HQ)
TWIRL = twirl(BOLTZ)
OPS = {
    "Strom E^+_0 (z99, Vertex e^(i(phi_A - phi_B)))": JP.astype(complex),
    "Strom E^+_1": (C4[0].T @ C4[3]).astype(complex),
    "Strom E^+_-1": (C4[1].T @ C4[2]).astype(complex),
    "Strom E^-_0": JM.astype(complex),
    "Cartan H_0": np.diag(QA - QB).astype(complex),
    "Cartan H_1": (C4[0].T @ C4[1] - C4[2].T @ C4[3]).astype(complex),
    "Kreuzpaar psi*_A psi*_B (Vertex e^(i(phi_A + phi_B)))": (C4[0].T @ C4[2].T).astype(complex),
    "Doppeltransfer C^+ (z105)": (C4[0].T @ C4[1].T @ C4[3] @ C4[2]).astype(complex),
    "Superaustausch (virtuell, z119)": K_VIRT.astype(complex),
    "Twirl e^(-4 pi H_Q) über die Transfer-U(2) (z116)": TWIRL,
    "Sandwich U_pi/2 (z120)": gatter(np.pi/2),
    "Sandwich U_pi/3": gatter(np.pi/3),
    "Transferdrehung B_pi/3": expm(-1j*np.pi/3*JX),
    "kanonischer Tausch Gamma(pi) (z106)": GAMMA_PI,
}
KLASS = {}
for name, X in OPS.items():
    fr = paritaetsklassen(X)
    KLASS[name] = {"Parität": klassen_text(fr), "Paaruhr": klassen_text(uhrklassen(X, R_PAAR)),
                   "Produktuhr": klassen_text(uhrklassen(X, R_PROD)), "_fr": fr}
    say(f"      {name}: Parität {KLASS[name]['Parität']}, Paaruhr {KLASS[name]['Paaruhr']}, "
        f"Produktuhr {KLASS[name]['Produktuhr']}")
tausch_fehl = max(nrm(GAMMA_PI @ C4[s] @ GAMMA_PI.conj().T - C4[2 + s]) + nrm(GAMMA_PI @ C4[2 + s] @ GAMMA_PI.conj().T - C4[s])
                  for s in range(2))
tausch_par = nrm(GAMMA_PI @ GA @ GAMMA_PI.conj().T - GB)
erw = {"Strom E^+_0 (z99, Vertex e^(i(phi_A - phi_B)))": ("rein (-,-)", "rein k = 2", "rein k = 0"),
       "Strom E^+_1": ("rein (-,-)", "rein k = 2", "rein k = 0"), "Strom E^+_-1": ("rein (-,-)", "rein k = 2", "rein k = 0"),
       "Strom E^-_0": ("rein (-,-)", "rein k = 2", "rein k = 0"), "Cartan H_0": ("rein (+,+)", "rein k = 0", "rein k = 0"),
       "Cartan H_1": ("rein (+,+)", "rein k = 0", "rein k = 0"),
       "Kreuzpaar psi*_A psi*_B (Vertex e^(i(phi_A + phi_B)))": ("rein (-,-)", "rein k = 0", "rein k = 2"),
       "Doppeltransfer C^+ (z105)": ("rein (+,+)", "rein k = 0", "rein k = 0"),
       "Superaustausch (virtuell, z119)": ("rein (+,+)", "rein k = 0", "rein k = 0"),
       "Twirl e^(-4 pi H_Q) über die Transfer-U(2) (z116)": ("rein (+,+)", "rein k = 0", "rein k = 0"),
       "Sandwich U_pi/2 (z120)": ("rein (+,+)", "rein k = 0", "rein k = 0")}
gemischt = ["Sandwich U_pi/3", "Transferdrehung B_pi/3", "kanonischer Tausch Gamma(pi) (z106)"]
klass_ok = all((KLASS[k]["Parität"], KLASS[k]["Paaruhr"], KLASS[k]["Produktuhr"]) == v for k, v in erw.items()) and all(
    KLASS[k]["_fr"][(1, 1)] > 1e-3 and KLASS[k]["_fr"][(-1, -1)] > 1e-3 and KLASS[k]["Produktuhr"] == "rein k = 0"
    for k in gemischt)
ZAHLEN["klassifikation"] = {k: {kk: vv for kk, vv in v.items() if kk != "_fr"} for k, v in KLASS.items()}
ZAHLEN["klassifikation_tausch"] = {"Gamma(pi) c_A <-> c_B (Fehler)": zahl(tausch_fehl), "Gamma(pi) (-1)^Q_A Gamma(pi)^+ - (-1)^Q_B":
                                   zahl(tausch_par)}
check(f"[E] Klassifikation (zwei Nähte, je zwei Moden, 16 Zustände) nach den Nahtparitäten (-1)^Q_A, (-1)^Q_B und nach Paaruhr "
      f"i^(Q_A - Q_B) bzw. Produktuhr i^(Q_A + Q_B) (Klasse k: Ad Uhr = i^k): z99s Ströme E^+-_n (n = 0, +-1; Vertex "
      f"e^(+-i(phi_A - phi_B))) und damit die Transferströme J^+- sind rein (-,-), also Einzelfermiontransfers, Paaruhr k = 2, "
      f"Produktuhr k = 0; das Kreuzpaar psi*_A psi*_B ist rein (-,-) mit Paaruhr k = 0 und Produktuhr k = 2; rein (+,+) und "
      f"uhrneutral sind die Cartanströme H_0, H_1, z105s Doppeltransfer, die virtuellen Transfers (Superaustausch zweiter "
      f"Ordnung, Twirl von e^(-4 pi H_Q)) und das Sandwich U_pi/2; gemischt sind U_pi/3 "
      f"({KLASS['Sandwich U_pi/3']['Parität']}), B_pi/3 ({KLASS['Transferdrehung B_pi/3']['Parität']}) und der Tausch "
      f"Gamma(pi) ({KLASS['kanonischer Tausch Gamma(pi) (z106)']['Parität']}; er vertauscht c_A <-> c_B, Fehler "
      f"{klein(tausch_fehl)}, und (-1)^Q_A <-> (-1)^Q_B, Rest {klein(tausch_par)})",
      klass_ok and tausch_fehl < 1e-12 and tausch_par < 1e-12)

# 1.2 Warum die virtuellen Transfers gerade sind: die Nahtparität liegt in der Transfer-U(2)
jz_rel = nrm(JZ - (JP @ JM - JM @ JP)/2)
ga_u2 = nrm(GA - expm(1j*np.pi*(NOP/2 + JZ)))
X_zuf = rng.normal(size=(16, 16)) + 1j*rng.normal(size=(16, 16))
X_zuf = X_zuf + X_zuf.conj().T
tw_zuf = twirl(X_zuf)
P_JX, _ = kommutant_projektor([JX])
tw_x = twirl(X_zuf, P_JX)
i11 = np.where((QA == 1) & (QB == 1))[0]
ew11 = np.sort(np.linalg.eigvalsh(TWIRL[np.ix_(i11, i11)].real))/np.exp(-4*np.pi)
w_S = float(ew11[0])
z116 = json.loads((HERE / "z116_ergebnis.json").read_text())
w_S_z116 = float(z116["zahlen"]["austausch"]["w_Singulett exakt"])
ZAHLEN["transfer_u2"] = {"|J_z - [J^+, J^-]/2|": zahl(jz_rel), "|(-1)^Q_A - e^(i pi (N/2 + J_z))|": zahl(ga_u2),
                         "dim Kommutant der Transfer-U(2)": DIM_U2, "w_Singulett (Twirl, relativ e^-4pi)": zahl(w_S, 12),
                         "w_Singulett z116": w_S_z116, "|[Twirl(X), (-1)^Q_A]| zufällig": zahl(kom(tw_zuf, GA)),
                         "|[Twirl_Jx(X), (-1)^Q_A]| (nur e^(-i theta J_x))": zahl(kom(tw_x, GA))}
check(f"[E] Die Nahtparität ist ein Element der Transfer-U(2): (-1)^Q_A = e^(i pi (N/2 + J_z)) mit J_z = [J^+, J^-]/2 (Reste "
      f"{klein(ga_u2)}, {klein(jz_rel)}), die halbe Drehung um die Nahtachse (z125). Jedes Haar-Mittel über die U(2) (Projektion "
      f"auf den Kommutanten, Dimension {DIM_U2}) ist daher nahtgerade: Twirl eines zufälligen hermiteschen Operators "
      f"|[., (-1)^Q_A]| = {klein(kom(tw_zuf, GA))}; z116s Gewicht w_S = (1 + 2 e^-4pi)/3 = {w_S:.10f} (z116: {w_S_z116:.10f}). "
      f"Kontrolle: Mittel nur über die Transferdrehung e^(-i theta J_x) (enthält die Halbdrehung nicht) gibt "
      f"{kom(tw_x, GA):.3f} != 0. Die Geradheit der virtuellen Transfers ist also strukturell (Gruppenmittel), nicht dynamisch",
      ga_u2 < 1e-12 and jz_rel < 1e-12 and kom(tw_zuf, GA) < 1e-10 and kom(TWIRL, JP) < 1e-12 and kom(TWIRL, GA) < 1e-12
      and abs(w_S - (1 + 2*np.exp(-4*np.pi))/3) < 1e-12 and abs(w_S - w_S_z116) < 1e-10 and kom(tw_x, GA) > 0.1)

# 1.3 J-Bilder: zwei Nähte mit je eigenem Partner (verdoppeltes Modell, Lesart zweier gleichchiraler Nähte)
d = 16
p_ab = np.exp(-np.pi*(QA + QB).astype(float))
p_ab /= p_ab.sum()
OMM = np.diag(np.sqrt(p_ab))                                    # Omega = vec(rho^(1/2)), Produktzustand rho_A (x) rho_B
OM2 = OMM.reshape(-1)
Tm = np.zeros((d*d, d*d), dtype=complex)
Tdm = np.zeros((d*d, d*d), dtype=complex)
for k_, (i, j) in enumerate((i, j) for i in range(d) for j in range(d)):
    v = np.zeros((d, d), dtype=complex)
    v[i, :] = OMM[j, :]
    Tm[:, k_] = v.reshape(-1)
    v = np.zeros((d, d), dtype=complex)
    v[j, :] = OMM[i, :]
    Tdm[:, k_] = v.reshape(-1)
L2 = Tdm @ np.linalg.inv(Tm).conj()                             # S v = L conj(v)
D2 = (L2.conj().T @ L2).conj()
w2, V2 = np.linalg.eigh((D2 + D2.conj().T)/2)
LJ2 = L2 @ (V2 @ np.diag(w2**-0.5) @ V2.conj().T).conj()        # J = L_J K


def emb(x):
    return np.kron(x, I16)


def J2(X):
    return LJ2 @ X.conj() @ LJ2.conj().T


GA2, GB2, GAq, GBq = emb(GA), emb(GB), np.kron(I16, GA), np.kron(I16, GB)
DA, DB = GA2 @ GAq, GB2 @ GBq                                   # kanonische Parität der verdoppelten Nähte
LA = np.diag(np.kron(1j**(QA % 4), 1j**((-QA) % 4)))            # i^(Q_A - Q_A-quer), Niveaugleichheit je Naht
LB = np.diag(np.kron(1j**(QB % 4), 1j**((-QB) % 4)))
LTOT, GTOT = LA @ LB, GA2 @ GB2 @ GAq @ GBq
j_par = max(nrm(J2(GA2) - GAq), nrm(J2(GB2) - GBq))
j_omega = nrm(LJ2 @ OM2.conj() - OM2)
J_LISTE = ["Strom E^+_0 (z99, Vertex e^(i(phi_A - phi_B)))", "Cartan H_0", "Doppeltransfer C^+ (z105)", "Sandwich U_pi/2 (z120)",
           "Sandwich U_pi/3", "Transferdrehung B_pi/3", "kanonischer Tausch Gamma(pi) (z106)"]
j_tab, j_klasse_ok = {}, True
for name in J_LISTE:
    x = OPS[name]
    X1 = emb(x)
    XJ = J2(X1)
    fr_o = paritaetsklassen(X1, GA2, GB2)
    fr_p = paritaetsklassen(XJ, GAq, GBq)
    j_klasse_ok &= max(abs(fr_o[k] - fr_p[k]) for k in fr_o) < 1e-12
    XX = X1 @ XJ
    homogen = max(KLASS[name]["_fr"].values()) > 1 - 1e-12
    j_tab[name] = {"homogen": homogen, "J-gepaart: [., (-1)^(F_A + F_A-quer)]": zahl(kom(XX, DA)),
                   "J-gepaart: [., Niveau je Naht L_A]": zahl(kom(XX, LA)), "J-gepaart: [., L_A L_B]": zahl(kom(XX, LTOT)),
                   "x (x) 1: [., (-1)^(F_A + F_A-quer)]": zahl(kom(X1, DA))}
homogen_ok = all((v["J-gepaart: [., (-1)^(F_A + F_A-quer)]"] == 0.0) == v["homogen"] for v in j_tab.values())
scan = {"L_A L_B (Gesamtdrehung) x (x) 1": [], "L_A L_B (Gesamtdrehung) J-gepaart": [], "(-1)^(F_A + F_A-quer) x (x) 1": [],
        "(-1)^(F_A + F_A-quer) J-gepaart": [], "L_A je Naht x (x) 1": [], "L_A je Naht J-gepaart": []}
TH37 = np.linspace(0, np.pi, 37)
for th in TH37:
    X1 = emb(gatter(th))
    XJ = X1 @ J2(X1)
    for art, X in (("x (x) 1", X1), ("J-gepaart", XJ)):
        if kom(X, LTOT) < 1e-9 and kom(X, GTOT) < 1e-9:
            scan[f"L_A L_B (Gesamtdrehung) {art}"].append(round(float(th/np.pi), 4))
        if kom(X, DA) < 1e-9 and kom(X, DB) < 1e-9:
            scan[f"(-1)^(F_A + F_A-quer) {art}"].append(round(float(th/np.pi), 4))
        if kom(X, LA) < 1e-9 and kom(X, LB) < 1e-9:
            scan[f"L_A je Naht {art}"].append(round(float(th/np.pi), 4))
ZAHLEN["j_bilder_verdoppelt"] = {"|J (-1)^Q_A J - (-1)^Q_A-quer|": zahl(j_par), "|J Omega - Omega|": zahl(j_omega),
                                 "Operatoren": j_tab,
                                 "Gatterscan theta/pi (37 Winkel)": {k: (v if len(v) < 6 else f"{len(v)} von 37")
                                                                     for k, v in scan.items()}}
check(f"[E] J-Bilder (zwei Nähte mit je eigenem Partner, Produkt-KMS-Zustand, 256 Zustände, J aus S = J Delta^(1/2)): J bildet "
      f"(-1)^Q_A auf (-1)^Q_A-quer ab (Rest {klein(j_par)}), J Omega = Omega ({klein(j_omega)}), und die Paritätsklasse jedes "
      f"Operators auf dieselbe Klasse der Partner ({j_klasse_ok}). Gepaart mit seinem J-Bild, x J x J, vertauscht ein Operator "
      f"genau dann mit der kanonischen Parität (-1)^(F_A + F_A-quer) jeder verdoppelten Naht, wenn er homogen ist "
      f"({homogen_ok}: Strom, Cartan, Doppeltransfer, U_pi/2 ja; U_pi/3, B_pi/3, Gamma(pi) nein). Gatterscan: die "
      f"Gesamtdrehung L_A L_B = i^(Q_A - Q_A-quer + Q_B - Q_B-quer) und die Gesamtparität fordern nichts "
      f"({scan['L_A L_B (Gesamtdrehung) x (x) 1'] if len(scan['L_A L_B (Gesamtdrehung) x (x) 1']) < 6 else '37 von 37'} "
      f"Winkel, auch J-gepaart), die kanonische Parität je verdoppelter Naht erzwingt theta/pi = "
      f"{scan['(-1)^(F_A + F_A-quer) x (x) 1']} (J-gepaart {scan['(-1)^(F_A + F_A-quer) J-gepaart']}), die Niveaugleichheit "
      f"je Naht lässt nur {scan['L_A je Naht x (x) 1']} (J-gepaart {scan['L_A je Naht J-gepaart']}): nur triviale Gatter. "
      f"J-gepaart erhalten Strom und Doppeltransfer L_A (Reste {klein(j_tab[J_LISTE[0]]['J-gepaart: [., Niveau je Naht L_A]'])}, "
      f"{klein(j_tab['Doppeltransfer C^+ (z105)']['J-gepaart: [., Niveau je Naht L_A]'])}), U_pi/2 nicht "
      f"({j_tab['Sandwich U_pi/2 (z120)']['J-gepaart: [., Niveau je Naht L_A]']:.3f}: es mischt (2,0) und (0,2))",
      j_par < 1e-12 and j_omega < 1e-12 and j_klasse_ok and homogen_ok
      and len(scan["L_A L_B (Gesamtdrehung) x (x) 1"]) == 37 and len(scan["L_A L_B (Gesamtdrehung) J-gepaart"]) == 37
      and scan["(-1)^(F_A + F_A-quer) x (x) 1"] == [0.0, 0.5, 1.0] and scan["(-1)^(F_A + F_A-quer) J-gepaart"] == [0.0, 0.5, 1.0]
      and scan["L_A je Naht x (x) 1"] == [0.0, 1.0] and scan["L_A je Naht J-gepaart"] == [0.0, 1.0]
      and j_tab[J_LISTE[0]]["J-gepaart: [., Niveau je Naht L_A]"] == 0.0
      and j_tab["Doppeltransfer C^+ (z105)"]["J-gepaart: [., Niveau je Naht L_A]"] == 0.0
      and j_tab["Sandwich U_pi/2 (z120)"]["J-gepaart: [., Niveau je Naht L_A]"] > 0.1)
befund("E: Native Einzelfermiontransfers gibt es als Erzeuger: z99s su(2)_1-Ströme E^+-_n (der Vertex e^(i(phi_A - phi_B))) und "
       "die Transferströme J^+- sind in jeder Nahtparität ungerade, für zwei gleichchirale Nähte aber uhrneutral (Produktuhr "
       "k = 0, h = 1) und damit lokal-physikalisch. Die vollendeten Kompositionen und Effektivterme (z105-Doppeltransfer, "
       "z116/z119-Twirl und Superaustausch, z120-Sandwich bei pi/2) sind nahtgerade; beim Twirl, weil die Nahtparität selbst "
       "die Halbdrehung der Transfer-U(2) ist. V5s wörtliche Prüfung ('jeder native Zwischennahtoperator ist je Naht gerade') "
       "fällt auf der Ebene der Erzeuger negativ aus; mit dem J-Bild gepaart gilt sie in der kanonischen Parität der "
       "verdoppelten Naht genau für homogene Operatoren", True)

# =================================================================== Teil 2
say("\nTeil 2: Standardform einer Naht mit ihrem Tomita-Partner (A = Original, B = Partner)")
CA = jw(2)                                                       # Moden der Naht auf ihrem 4-dim. Raum (Index 2 n_0 + n_1)
QA4 = np.array([0, 1, 1, 2])
N0, N1 = np.array([0, 0, 1, 1]), np.array([0, 1, 0, 1])
BX = []
for i in range(4):
    for j in range(4):
        m_ = np.zeros((4, 4), dtype=complex)
        m_[i, j] = 1
        BX.append(np.kron(m_, np.eye(4)))
M_ERZ = [C4[0], C4[1], C4[0].T, C4[1].T]                         # Erzeuger von M = B(H_A) (x) 1


def standardform(rhoA):
    """Omega = vec(rho^(1/2)) (Zeilenordnung), S aus der Definition S x Omega = x^+ Omega auf allen 16 Matrixeinheiten,
    Delta = S^* S, J = S Delta^(-1/2) = L_J K."""
    w, V = np.linalg.eigh(rhoA)
    Om = (V @ np.diag(np.sqrt(w)) @ V.conj().T).reshape(-1)
    T = np.array([x @ Om for x in BX]).T
    Td = np.array([x.conj().T @ Om for x in BX]).T
    L = Td @ np.linalg.inv(T).conj()
    Dl = (L.conj().T @ L).conj()
    wd, Vd = np.linalg.eigh((Dl + Dl.conj().T)/2)
    LJ = L @ (Vd @ np.diag(wd**-0.5) @ Vd.conj().T).conj()
    return {"Om": Om, "Delta": Dl, "LJ": LJ, "wd": wd, "Vd": Vd, "rho": rhoA}


def Jop(sf, X):
    return sf["LJ"] @ X.conj() @ np.linalg.inv(sf["LJ"])


def Jvec(sf, v):
    return sf["LJ"] @ v.conj()


def Dit(sf, t):
    return sf["Vd"] @ np.diag(sf["wd"]**(1j*t)) @ sf["Vd"].conj().T


def rho_aus(K):
    r = expm(-K)
    return r/np.trace(r).real


GAMMA0 = CA[0] + CA[0].T                                         # Majorana der Naht (ungerade)
KVAR = {"pi Q (z125)": np.diag(np.pi*QA4).astype(complex),
        "4 pi E_lad(1/4, Q) = 2 pi Q^2 - pi Q": np.diag(2*np.pi*QA4**2 - np.pi*QA4).astype(complex),
        "4 pi Q^2/2 (ohne Viertelverschiebung)": np.diag(2*np.pi*QA4**2).astype(complex),
        "zwei Viertel (eps = pi, 3 pi)": np.diag(np.pi*N0 + 3*np.pi*N1).astype(complex),
        "paritätsbrechend pi Q + 0,7 gamma_0 (Kontrolle)": np.diag(np.pi*QA4).astype(complex) + 0.7*GAMMA0}
SF = {k: standardform(rho_aus(K)) for k, K in KVAR.items()}
GERADE = list(KVAR)[:4]

# 2.1 Tomita-Grunddaten
tom = {}
for k, sf in SF.items():
    rho = sf["rho"]
    ref = np.kron(rho, np.linalg.inv(rho).T)
    tom[k] = {"|Delta - rho (x) rho^-1^T|/|.|": zahl(nrm(sf["Delta"] - ref)/nrm(ref), 2, 1e-10),
              "|J^2 - 1|": zahl(nrm(sf["LJ"] @ sf["LJ"].conj() - I16), 2, 1e-10),
              "|J Omega - Omega|": zahl(nrm(Jvec(sf, sf["Om"]) - sf["Om"]), 2, 1e-10),
              "max |[J x J, y]| (x, y Erzeuger von M)": zahl(max(kom(Jop(sf, x), y) for x in M_ERZ for y in M_ERZ), 2, 1e-10),
              "|J M J| != 0": zahl(min(nrm(Jop(sf, x)) for x in M_ERZ))}
ZAHLEN["tomita_grunddaten"] = tom
check(f"[E] Standardform für fünf Zustände der Naht (zwei Moden, KMS bei beta = 4 pi; K = pi Q wie z125, K = 4 pi E_lad(1/4, Q) "
      f"= 2 pi Q^2 - pi Q mit Viertelverschiebung (z113), K = 4 pi Q^2/2 ohne sie (z119), zwei Viertel mit eps = pi, 3 pi, "
      f"und paritätsbrechend als Kontrolle): S aus der Definition auf allen 16 Matrixeinheiten, Delta = rho (x) rho^-1 (relativ "
      f"{klein(max(v['|Delta - rho (x) rho^-1^T|/|.|'] for v in tom.values()), 1e-10)}), J^2 = 1, J Omega = Omega ("
      f"{klein(max(v['|J^2 - 1|'] + v['|J Omega - Omega|'] for v in tom.values()), 1e-10)}), J M J vertauscht mit M ("
      f"{klein(max(v['max |[J x J, y]| (x, y Erzeuger von M)'] for v in tom.values()), 1e-10)}; die Reste stammen aus "
      f"Gewichten bis e^(-8 pi)): M' = J M J",
      all(v["|Delta - rho (x) rho^-1^T|/|.|"] < 1e-6 and v["|J^2 - 1|"] < 1e-10 and v["|J Omega - Omega|"] < 1e-10
          and v["max |[J x J, y]| (x, y Erzeuger von M)"] < 1e-10 and v["|J M J| != 0"] > 0.5 for v in tom.values()))

# 2.2 Modulgruppe
TS = [-1.3, -0.5, 0.25, 0.7, 2.1]
mod = {}
for k, sf in SF.items():
    mod[k] = {"max_t |[Delta^it, (-1)^Q_A]|": zahl(max(kom(Dit(sf, t), GA) for t in TS)),
              "max_t |[Delta^it, (-1)^Q_B]|": zahl(max(kom(Dit(sf, t), GB) for t in TS)),
              "max_t |[Delta^it, Q_A]|": zahl(max(kom(Dit(sf, t), np.diag(QA.astype(float))) for t in TS))}
ZAHLEN["modulgruppe"] = mod
kontrolle = mod["paritätsbrechend pi Q + 0,7 gamma_0 (Kontrolle)"]
check(f"[E] Modulgruppe: Für die vier geraden Zustände vertauscht Delta^it (fünf Zeiten t) mit (-1)^Q_A und (-1)^Q_B und erhält "
      f"Q_A und Q_B einzeln (max {klein(max(max(v.values()) for k, v in mod.items() if k in GERADE))}): sie transportiert keine "
      f"Ladung zwischen Naht und Partner. Kontrolle: Der paritätsbrechende Zustand (gamma_0 = a_0 + a_0^+ im Modulhamilton) "
      f"bricht das (|[Delta^it, (-1)^Q_A]| = {kontrolle['max_t |[Delta^it, (-1)^Q_A]|']:.3f}): die Geradheit der Modulgruppe "
      f"ist die Univalenz des Zustands, nicht der Algebra",
      all(max(v.values()) < 1e-10 for k, v in mod.items() if k in GERADE) and kontrolle["max_t |[Delta^it, (-1)^Q_A]|"] > 0.1)

# 2.3 Delta^(-i/2) als 2 pi-Drehung
R_ZWEI = np.diag(1j**(((NUM[0] - NUM[1]) - (NUM[2] - NUM[3])) % 4))
UHREN = {"Paaruhr i^(Q_A - Q_B)": R_PAAR, "Produktuhr i^(Q_A + Q_B)": R_PROD, "Gesamtparität": GAM,
         "zwei Viertel i^((Q_a - Q_b) - (Q_a-quer - Q_b-quer))": R_ZWEI}


def phasenspanne(D, R):
    if nrm(D - np.diag(np.diag(D))) > 1e-10:
        return 99.0
    return float(np.ptp(np.angle(np.diag(D)/np.diag(R))))


dmi = {k: {u: zahl(phasenspanne(Dit(SF[k], -0.5), R), 4) for u, R in UHREN.items()} for k in GERADE}
treffer = {k: sorted(u for u, s in v.items() if s < 1e-10) for k, v in dmi.items()}
e_lad_ok = all((Fraction(2*Q*Q - Q, 4) - Fraction(Q, 4)).denominator == 1
               and (Fraction(2*Q*Q + Q, 4) + Fraction(Q, 4)).denominator == 1
               and (Fraction(Q*Q, 2) - Fraction(Q, 2)).denominator == 1 for Q in range(-8, 9))
ZAHLEN["delta_-i/2"] = {"Phasenspanne gegen Uhr": dmi, "Treffer": treffer,
                        "e^(2 pi i E_lad(1/4,Q)) = i^Q, e^(2 pi i E_lad(3/4,Q)) = i^-Q, e^(2 pi i Q^2/2) = (-1)^Q (|Q| <= 8)": e_lad_ok}
check(f"[E] Delta^(-i/2) ist die 2 pi-Drehung der verdoppelten Naht: für K = pi Q und K = 4 pi E_lad(1/4, Q) genau die Paaruhr "
      f"i^(Q_A - Q_B) = R R-quer^-1 (Treffer {treffer['pi Q (z125)']}, {treffer['4 pi E_lad(1/4, Q) = 2 pi Q^2 - pi Q']}); "
      f"exakt gilt e^(2 pi i E_lad(1/4, Q)) = i^Q und e^(2 pi i E_lad(3/4, Q)) = i^-Q ({e_lad_ok}, |Q| <= 8), also "
      f"e^(2 pi i (L_0 - L-quer_0)) = R R-quer^-1 auf dem Ladungsgitter (z113) und P_clock bei beta = 4 pi. Mit zwei Vierteln "
      f"(eps = pi, 3 pi) ist es R R-quer^-1 mit R = i^(Q_a - Q_b) ({treffer['zwei Viertel (eps = pi, 3 pi)']}). Kontrolle: "
      f"ohne Viertelverschiebung (4 pi Q^2/2, e^(2 pi i Q^2/2) = (-1)^Q) ist Delta^(-i/2) nur die Gesamtparität "
      f"({treffer['4 pi Q^2/2 (ohne Viertelverschiebung)']}); die Produktuhr trifft nie",
      treffer["pi Q (z125)"] == ["Paaruhr i^(Q_A - Q_B)"] and treffer["4 pi E_lad(1/4, Q) = 2 pi Q^2 - pi Q"] == ["Paaruhr i^(Q_A - Q_B)"]
      and treffer["4 pi Q^2/2 (ohne Viertelverschiebung)"] == ["Gesamtparität"]
      and treffer["zwei Viertel (eps = pi, 3 pi)"] == ["zwei Viertel i^((Q_a - Q_b) - (Q_a-quer - Q_b-quer))"] and e_lad_ok)

# 2.4 J: Paritäten, Uhren, Klassen
SF0 = SF["pi Q (z125)"]
SF1 = SF["4 pi E_lad(1/4, Q) = 2 pi Q^2 - pi Q"]
j_uhr = {}
for k, sf in (("pi Q", SF0), ("4 pi E_lad", SF1)):
    j_uhr[k] = {"|J (-1)^Q_A J - (-1)^Q_B|": zahl(nrm(Jop(sf, GA) - GB)),
                "|J Paaruhr J - Paaruhr|": zahl(nrm(Jop(sf, R_PAAR) - R_PAAR)),
                "|J Produktuhr J - Produktuhr^-1|": zahl(nrm(Jop(sf, R_PROD) - R_PROD.conj().T)),
                "|J Produktuhr J - Produktuhr|": zahl(nrm(Jop(sf, R_PROD) - R_PROD))}
kl_map_ok, uhr_map_ok = True, True
for sA in (1, -1):
    for sB in (1, -1):
        X = rng.normal(size=(16, 16)) + 1j*rng.normal(size=(16, 16))
        X = (X + sA*GA @ X @ GA + sB*GB @ X @ GB + sA*sB*GAM @ X @ GAM)/4
        kl_map_ok &= paritaetsklassen(Jop(SF0, X))[(sB, sA)] > 1 - 1e-12
for kk in range(4):
    X = rng.normal(size=(16, 16)) + 1j*rng.normal(size=(16, 16))
    X = sum((1j)**(-kk*m)*np.linalg.matrix_power(R_PAAR, m) @ X @ np.linalg.matrix_power(R_PAAR.conj().T, m) for m in range(4))/4
    uhr_map_ok &= uhrklassen(Jop(SF0, X), R_PAAR)[(-kk) % 4] > 1 - 1e-12
ZAHLEN["J_uhren"] = {"Reste": j_uhr, "Klasse (sA, sB) -> (sB, sA)": bool(kl_map_ok), "Paaruhrklasse k -> -k": bool(uhr_map_ok)}
check(f"[E] J vertauscht die Seitenparitäten (J (-1)^Q_A J = (-1)^Q_B, Rest {klein(j_uhr['pi Q']['|J (-1)^Q_A J - (-1)^Q_B|'])}), "
      f"bildet die Klasse (sA, sB) auf (sB, sA) ab ({kl_map_ok}) und die Paaruhrklasse k auf -k ({uhr_map_ok}); es erhält die "
      f"Paaruhr (Rest {klein(j_uhr['pi Q']['|J Paaruhr J - Paaruhr|'])}, wie z122: J R_phys J = R_phys) und invertiert die "
      f"Produktuhr (Rest {klein(j_uhr['pi Q']['|J Produktuhr J - Produktuhr^-1|'])}, Abstand zur Produktuhr selbst "
      f"{j_uhr['pi Q']['|J Produktuhr J - Produktuhr|']:.3f}): Der Tomita-Partner ist gegenchiral, J ist ein Tausch, kein "
      "Transfer, und aus Einzeltransfers macht es wieder Einzeltransfers (k = 2 -> 2)",
      all(v["|J (-1)^Q_A J - (-1)^Q_B|"] == 0.0 and v["|J Paaruhr J - Paaruhr|"] == 0.0 and v["|J Produktuhr J - Produktuhr^-1|"]
          == 0.0 and v["|J Produktuhr J - Produktuhr|"] > 1 for v in j_uhr.values()) and kl_map_ok and uhr_map_ok)

# 2.5 Kommutant, gradierter Kommutant, Einzeltransfer als Operator
ZK = (I16 + 1j*GAM)/(1 + 1j)                                     # Klein-Verdrehung
NAHTPAAR = C4[0].T @ C4[1].T                                     # a_0^+ a_1^+: nahtgerade, Ladung +2 auf der Naht allein
KREUZPAAR = C4[0].T @ C4[2].T                                    # a_0^+ b_0^+: Kreuzpaar Naht/Partner
enth = {"Kreuzpaar: [., Paaruhr]": zahl(kom(KREUZPAAR, R_PAAR)), "Kreuzpaar: [., (-1)^Q_A]": zahl(kom(KREUZPAAR, GA), 4),
        "Nahtpaar: [., (-1)^Q_A] + [., (-1)^Q_B]": zahl(kom(NAHTPAAR, GA) + kom(NAHTPAAR, GB)),
        "Nahtpaar: [., Paaruhr]": zahl(kom(NAHTPAAR, R_PAAR), 4)}
komm = {}
for c_ in range(2):
    a, b = C4[c_], C4[2 + c_]                                    # a_c Naht, b_c gradierter Partner (JW)
    f = Jop(SF0, a)                                              # f_c = J a_c J in M'
    th_ = ZK.conj().T @ f @ ZK                                   # Theta a_c Theta^*, Theta = Z^* J
    komm[c_] = {"max |[f_c, a_d]|": zahl(max(kom(f, C4[e]) for e in range(2))),
                "|f_c - (-1)^Q_A b_c|": zahl(nrm(f - GA @ b)),
                "max |{Theta a_c Theta^*, a_d}|": zahl(max(nrm(th_ @ C4[e] + C4[e] @ th_) for e in range(2))),
                "|Theta a_c Theta^* + i (-1)^Q_B b_c|": zahl(nrm(th_ + 1j*GB @ b)),
                "max |{b_c, a_d}|": zahl(max(nrm(b @ C4[e] + C4[e] @ b) for e in range(2)))}
z2 = nrm(ZK @ ZK - GAM)
T_EIN = C4[2].T @ C4[0]                                          # b_0^+ a_0: ein Fermion Naht -> Partner
m_ = GA @ C4[0]
mp_ = (GA @ C4[2]).conj().T
t_mm = nrm(T_EIN - mp_ @ m_)
mp_in = max(kom(mp_, x) for x in M_ERZ)
m_in = max(kom(m_, Jop(SF0, x)) for x in M_ERZ)
t_kl = paritaetsklassen(T_EIN)
parts = {"globale Univalenz: Kommutant der Gesamtparität": blockpartition([np.diag(GAM)]),
         "Kommutant der Paaruhr": blockpartition([np.diag(R_PAAR)]),
         "Kommutant beider Nahtparitäten (= M_ev v M'_ev)": blockpartition([np.diag(GA), np.diag(GB)]),
         "Kommutant von N und Nahtparität": blockpartition([NT, np.diag(GA)]),
         "Kommutant von N und Paaruhr": blockpartition([NT, np.diag(R_PAAR)])}
dims = {k: sum(len(b)**2 for b in v) for k, v in parts.items()}
d_glob, d_paar, d_naht, d_nA, d_nR = dims.values()
n_gleich = parts["Kommutant von N und Nahtparität"] == parts["Kommutant von N und Paaruhr"]
ZAHLEN["kommutant"] = {"Moden": komm, "|Z^2 - Gamma|": zahl(z2), "|T - m' m|": zahl(t_mm), "m' in M' (Rest)": zahl(mp_in),
                       "m in M (Rest)": zahl(m_in), "Klasse T": klassen_text(t_kl), "Dimensionen": dims, "Enthaltensein": enth,
                       "mit N: dieselbe Blockzerlegung": n_gleich}
check(f"[E] Kommutantstruktur: f_c = J a_c J liegt im gewöhnlichen Kommutanten (vertauscht mit M, Rest "
      f"{klein(max(v['max |[f_c, a_d]|'] for v in komm.values()))}) und ist der Klein-verdrehte Partner (-1)^Q_A b_c (Rest "
      f"{klein(max(v['|f_c - (-1)^Q_A b_c|'] for v in komm.values()))}); Theta = Z^* J mit Z = (1 + i Gamma)/(1 + i), Z^2 = Gamma "
      f"({klein(z2)}), bildet a_c auf den gradierten Partner -i (-1)^Q_B b_c ab, der mit M antivertauscht (Koszul, Reste "
      f"{klein(max(v['max |{Theta a_c Theta^*, a_d}|'] + v['|Theta a_c Theta^* + i (-1)^Q_B b_c|'] for v in komm.values()))}): "
      f"Z vertauscht M' und den gradierten Kommutanten. Der Einzeltransfer T = b_0^+ a_0 (Naht -> Partner, {klassen_text(t_kl)}) "
      f"ist exakt m' m mit m = (-1)^Q_A a_0 in M und m' = ((-1)^Q_A b_0)^+ in M' (Reste {klein(t_mm)}, {klein(mp_in)}, "
      f"{klein(m_in)}): Er liegt in (M v M')^gerade. Dimensionen der Kommutanten: globale Univalenz {d_glob}, Paaruhr "
      f"{d_paar}, beide Nahtparitäten (M_ev v M'_ev = 8 x 8) {d_naht}; Paaruhr und Nahtparitäten enthalten einander nicht "
      f"(das Kreuzpaar a_0^+ b_0^+ vertauscht mit der Paaruhr, Rest {klein(enth['Kreuzpaar: [., Paaruhr]'])}, nicht mit "
      f"(-1)^Q_A, {enth['Kreuzpaar: [., (-1)^Q_A]']:.3f}; das Nahtpaar a_0^+ a_1^+ umgekehrt, Reste "
      f"{klein(enth['Nahtpaar: [., (-1)^Q_A] + [., (-1)^Q_B]'])} und {enth['Nahtpaar: [., Paaruhr]']:.3f}), mit N fallen sie zusammen "
      f"({d_nA} = {d_nR}, dieselbe Blockzerlegung: {n_gleich}). Kommutant, Koszul-Vorzeichen und Verdrehung schließen den "
      "Einzeltransfer also nicht aus (gescopt negativ für diesen Mechanismus)",
      all(v["max |[f_c, a_d]|"] == 0.0 and v["|f_c - (-1)^Q_A b_c|"] == 0.0 and v["max |{Theta a_c Theta^*, a_d}|"] == 0.0
          and v["|Theta a_c Theta^* + i (-1)^Q_B b_c|"] == 0.0 and v["max |{b_c, a_d}|"] == 0.0 for v in komm.values())
      and z2 < 1e-12 and t_mm < 1e-12 and mp_in < 1e-12 and m_in < 1e-12 and t_kl[(-1, -1)] > 1 - 1e-12
      and [d_glob, d_paar, d_naht, d_nA, d_nR] == [128, 72, 64, 38, 38] and n_gleich
      and enth["Kreuzpaar: [., Paaruhr]"] == 0.0 and enth["Kreuzpaar: [., (-1)^Q_A]"] > 0.5
      and enth["Nahtpaar: [., (-1)^Q_A] + [., (-1)^Q_B]"] == 0.0 and enth["Nahtpaar: [., Paaruhr]"] > 0.5)

# 2.6 KMS-Relationen
kms = {}
Om0 = SF0["Om"]
q0 = np.exp(-np.pi/2)
U_PAAR = np.diag((-1.0)**(QB*(QB - 1)//2))                       # Om_g = U Om: Paarprodukt-Purifikation (gleicher Zustand auf M)
Omg = U_PAAR @ Om0
for c_ in range(2):
    a, b, f = C4[c_], C4[2 + c_], Jop(SF0, C4[c_])
    kms[f"Mode {c_}"] = {"|a Om - q f^+ Om|": zahl(nrm(a @ Om0 - q0*f.conj().T @ Om0)),
                         "|a Om_g - q b^+ Om_g|": zahl(nrm(a @ Omg - q0*b.conj().T @ Omg)),
                         "min_q |a Om - q b^+ Om| (natürlicher Kegel, gradiert)": zahl(
                             nrm(a @ Om0 - (np.vdot(b.conj().T @ Om0, a @ Om0)/nrm(b.conj().T @ Om0)**2)*(b.conj().T @ Om0)))}
gleicher_zustand = max(abs(np.vdot(Omg, x @ Omg) - np.vdot(Om0, x @ Om0)) for x in BX)
tomita_form = {}
for k in GERADE:
    sf = SF[k]
    Dh = sf["Vd"] @ np.diag(sf["wd"]**0.5) @ sf["Vd"].conj().T
    tomita_form[k] = zahl(max(nrm(Jvec(sf, Dh @ (C4[c_] @ sf["Om"])) - C4[c_].T @ sf["Om"]) for c_ in range(2)))
Om1 = SF1["Om"]
v_e = Jop(SF1, C4[0]).conj().T @ Om1
skalar_elad = nrm(C4[0] @ Om1 - (np.vdot(v_e, C4[0] @ Om1)/nrm(v_e)**2)*v_e)     # bestes skalares q
f0, f1 = Jop(SF0, C4[0]), Jop(SF0, C4[1])
T_F = f1.conj().T @ C4[0]                                        # Transfer a_0 -> Partnermode 1 (M'-Fermion)
t_auf_om = nrm(T_F @ Om0 - q0*f1.conj().T @ f0.conj().T @ Om0)
theta_b = np.arctan(q0)
G_BOG = sum(theta_b*(C4[c_].T @ Jop(SF0, C4[c_]).conj().T - Jop(SF0, C4[c_]) @ C4[c_]) for c_ in range(2))
vak = np.zeros(16)
vak[0] = 1
bog = nrm(expm(G_BOG) @ vak - Om0)
bog_kl, bog_uhr = paritaetsklassen(G_BOG), uhrklassen(G_BOG, R_PAAR)
om_ga = float(np.vdot(Om0, GA @ Om0).real)                       # <(-1)^Q_A> = p_gerade - p_ungerade
om_ga_soll = ((1 - np.exp(-np.pi))/(1 + np.exp(-np.pi)))**2
ZAHLEN["kms"] = {"Moden": kms, "gleicher Zustand auf M (Om_g gegen Om)": zahl(gleicher_zustand),
                 "Tomita-Form J Delta^(1/2) a Om = a^+ Om (max)": tomita_form,
                 "skalare Form bei 4 pi E_lad (bester Rest)": zahl(skalar_elad), "|T Om - q f_1^+ f_0^+ Om|": zahl(t_auf_om),
                 "|e^G vac - Om| (thermischer Bogoljubow)": zahl(bog), "Klasse G": klassen_text(bog_kl),
                 "Paaruhrklasse G": klassen_text(bog_uhr), "<Omega, (-1)^Q_A Omega>": zahl(om_ga, 10)}
check(f"[E] KMS-Relationen (K = pi Q): a_c Omega = q f_c^+ Omega mit q = e^(-pi/2) für den Vektor im natürlichen Kegel und den "
      f"Kommutantpartner f_c = J a_c J (Rest {klein(max(v['|a Om - q f^+ Om|'] for v in kms.values()))}); in z122s gradierter "
      f"Form a_c Omega_g = q b_c^+ Omega_g für die Paarprodukt-Purifikation Omega_g = (-1)^(Q_B(Q_B-1)/2) Omega (Rest "
      f"{klein(max(v['|a Om_g - q b^+ Om_g|'] for v in kms.values()))}; gleicher Zustand auf M, {klein(gleicher_zustand)}), "
      f"nicht dagegen mit b_c und Omega zugleich (Rest {min(v['min_q |a Om - q b^+ Om| (natürlicher Kegel, gradiert)'] for v in kms.values()):.3f}). "
      f"Mit der Viertelverschiebung (4 pi E_lad, nicht quasifrei) gilt keine skalare Form (bester Rest {skalar_elad:.3f}), "
      f"die Tomita-Form J Delta^(1/2) a_c Omega = a_c^+ Omega aber für alle vier geraden Zustände (max "
      f"{klein(max(tomita_form.values()))}). Folgen: Der Transfer auf Omega ist eine Partner-Paarerzeugung (T Omega = q f_1^+ "
      f"f_0^+ Omega, Rest {klein(t_auf_om)}), und Omega selbst ist e^G |0> mit dem thermischen Bogoljubow-Erzeuger G = sum "
      f"theta (a^+ f^+ - f a), tan theta = q (Rest {klein(bog)}): G ist {klassen_text(bog_kl)}, aber Paaruhr {klassen_text(bog_uhr)}, "
      f"und Omega überlagert beide Nahtparitäten kohärent (<(-1)^Q_A> = {om_ga:.6f} = ((1 - e^-pi)/(1 + e^-pi))^2). Die "
      "KMS-Struktur koppelt Naht und Partner durch Kreuzpaare, nicht durch Transfers; eine strikte Nahtparitäts-Superauswahl "
      "verböte schon die Purifikation",
      all(v["|a Om - q f^+ Om|"] == 0.0 and v["|a Om_g - q b^+ Om_g|"] == 0.0
          and v["min_q |a Om - q b^+ Om| (natürlicher Kegel, gradiert)"] > 1e-3 for v in kms.values())
      and gleicher_zustand < 1e-12 and skalar_elad > 1e-3 and max(tomita_form.values()) < 1e-10 and t_auf_om < 1e-12
      and bog < 1e-12 and bog_kl[(-1, -1)] > 1 - 1e-12 and bog_uhr[0] > 1 - 1e-12 and abs(om_ga - om_ga_soll) < 1e-12
      and om_ga < 0.99)

# 2.7 Niveaugleichheit (z113) als physikalischer Raum
z113 = json.loads((HERE / "z113_ergebnis.json").read_text())
z113_check = [v for k, v in z113["checks"].items() if k.startswith("[E] Nichtchirale Naht, Bahn der Spur ohne Einschub")]
z113_koeff = z113["zahlen"]["nichtchirale_naht"]["1"]["1"]["enthält (rel. Koeffizient)"]
z113_proj = z113["zahlen"]["urteil_kurz"]["native Projektion"]
P_NIV = np.diag((((QA - QB) % 4) == 0).astype(float))
P_GSO = np.diag(((NT % 2) == 0).astype(float))
P_PHYS = P_NIV @ P_GSO
Q_PHYS = I16 - P_PHYS


def leck(X):
    return nrm(Q_PHYS @ X @ P_PHYS)


om_phys = {k: zahl(nrm(SF[k]["Om"] - P_PHYS @ SF[k]["Om"])) for k in GERADE}
erhalt = {"Delta^it (fünf t, beide Zustände)": max(leck(Dit(sf, t)) for sf in (SF0, SF1) for t in TS),
          "J (beide Zustände)": max(nrm(Jop(sf, P_PHYS) - P_PHYS) for sf in (SF0, SF1)),
          "Bogoljubow-Erzeuger G": leck(G_BOG), "Kreuzpaar a_0^+ f_0^+": leck(C4[0].T @ f0.conj().T),
          "Tausch Gamma(pi)": leck(GAMMA_PI), "Gatter U_pi/2": leck(gatter(np.pi/2)), "Transferdrehung B_pi": leck(expm(-1j*np.pi*JX)),
          "J-gepaartes Nahtpaar x J x J": leck(NAHTPAAR @ Jop(SF0, NAHTPAAR))}
verletzt = {"Einzeltransfer T = b_0^+ a_0": leck(T_EIN), "Strom E^+_0 (Naht <- Partner)": leck(JP),
            "Gatter U_pi/3": leck(gatter(np.pi/3)), "Transferdrehung B_pi/2": leck(expm(-1j*np.pi/2*JX)),
            "Nahtpaar a_0^+ a_1^+ allein (nahtgerade)": leck(NAHTPAAR)}
t_innen = nrm(P_PHYS @ T_EIN @ P_PHYS)
t_aussen = nrm(T_EIN @ P_PHYS)
ZAHLEN["niveaugleichheit"] = {"z113 Prüfung bestanden": z113_check, "z113 Koeffizienten": z113_koeff, "z113 Projektion": z113_proj,
                              "dim P_phys": int(round(np.trace(P_PHYS))), "|Omega - P_phys Omega|": om_phys,
                              "Leck erhaltender Objekte": {k: zahl(v) for k, v in erhalt.items()},
                              "Leck verletzender Objekte": {k: zahl(v, 4) for k, v in verletzt.items()},
                              "|P T P|, |T P|": [zahl(t_innen), zahl(t_aussen)]}
check(f"[E] Physikalischer Raum der verdoppelten Naht: z113 (Prüfung bestanden: {z113_check == [True]}) gibt als native "
      f"Projektion '{z113_proj}' mit festen Koeffizienten {z113_koeff}; hier P_phys = (i^(Q_A - Q_B) = 1) und (Gesamtparität = "
      f"1), Dimension {int(round(np.trace(P_PHYS)))} von 16. Omega liegt darin (max {klein(max(om_phys.values()))}, alle vier "
      f"geraden Zustände). Ihn erhalten Delta^it, J, der Bogoljubow-Erzeuger, das Kreuzpaar a^+ f^+, der Tausch Gamma(pi), "
      f"U_pi/2, B_pi und das mit seinem J-Bild gepaarte Nahtpaar a_0^+ a_1^+ J a_0^+ a_1^+ J (max Leck "
      f"{klein(max(erhalt.values()))}). Der Einzeltransfer Naht -> Partner hat dort kein Matrixelement: |P T P| = "
      f"{klein(t_innen)}, |T P| = {t_aussen:.3f}; er dreht die 2 pi-Drehung um -1 (Spin 1/2 mod 1). Ebenso verlassen den Raum "
      f"der Strom E^+_0 ({verletzt['Strom E^+_0 (Naht <- Partner)']:.3f}), U_pi/3 ({verletzt['Gatter U_pi/3']:.6f}), B_pi/2 "
      f"({verletzt['Transferdrehung B_pi/2']:.3f}) und das nahtgerade Nahtpaar allein "
      f"({verletzt['Nahtpaar a_0^+ a_1^+ allein (nahtgerade)']:.3f}): Der physikalische Raum ist nicht die Nahtparitätsklasse",
      z113_check == [True] and z113_koeff == {"1": "0", "(-1)^(F+Fq)": "0", "R Rq^-1": "0"}
      and z113_proj == "R = R-quer und (-1)^(F+F-quer) = +1" and max(om_phys.values()) < 1e-12 and max(erhalt.values()) < 1e-10
      and t_innen < 1e-12 and t_aussen > 0.5 and min(verletzt.values()) > 0.1)

# 2.8 Gleichwertigkeit mit P_univ auf den physikalischen Sektoren
ident = 0.0
for N in range(5):
    PN = np.diag((NT == N).astype(float))
    soll = (I16 + (-1)**(N//2)*GA)/2 @ PN if N % 2 == 0 else 0*PN
    ident = max(ident, nrm(P_PHYS @ PN - soll))
P_EVEN = np.diag(((NT % 2) == 0).astype(float))
P_PROD1 = np.diag((((QA + QB) % 4) == 0).astype(float))
lf = {"P_phys erhalten": [], "[U, (-1)^Q_A] = 0 auf geradem N": [], "[U, Paaruhr] = 0 auf geradem N": [],
      "Produktuhr-Projektion erhalten (Kontrolle)": []}
for th in TH181:
    U = gatter(th)
    tt = round(float(th/np.pi), 4)
    if leck(U) < 1e-10:
        lf["P_phys erhalten"].append(tt)
    if nrm((U @ GA - GA @ U) @ P_EVEN) < 1e-10:
        lf["[U, (-1)^Q_A] = 0 auf geradem N"].append(tt)
    if nrm((U @ R_PAAR - R_PAAR @ U) @ P_EVEN) < 1e-10:
        lf["[U, Paaruhr] = 0 auf geradem N"].append(tt)
    if nrm((I16 - P_PROD1) @ U @ P_PROD1) < 1e-10:
        lf["Produktuhr-Projektion erhalten (Kontrolle)"].append(tt)


def zufall_unitaer(bloecke):
    U = np.zeros((16, 16), dtype=complex)
    for b in bloecke:
        b = list(b)
        u = unitary_group.rvs(len(b), random_state=rng) if len(b) > 1 else np.array([[np.exp(2j*np.pi*rng.uniform())]])
        U[np.ix_(b, b)] = u
    return U


B_N = blockpartition([NT])
B_NP = blockpartition([NT, np.diag(P_PHYS)])
zuf = {"generisch": [], "physikalisch": []}
for art, bl in (("generisch", B_N), ("physikalisch", B_NP)):
    for _ in range(20):
        U = zufall_unitaer(bl)
        zuf[art].append((leck(U) < 1e-10, nrm((U @ GA - GA @ U) @ P_EVEN) < 1e-10, nrm((U @ GA - GA @ U) @ (I16 - P_EVEN))))
aequi = all(a == b for v in zuf.values() for a, b, _ in v)
z125 = json.loads((HERE / "z125_ergebnis.json").read_text())
sa = z125["zahlen"]["superauswahl_z4"]
z125_ok = "(1, 3)" in sa["erzwingt pi/2"] and "(1, 1)" in sa["leer"]
ZAHLEN["gleichwertigkeit"] = {"|P_phys P_N - (1 + (-1)^(N/2) (-1)^Q_A)/2 P_N|": zahl(ident),
                              "leckfreie Winkel theta/pi (181)": {k: (v if len(v) < 6 else f"{len(v)} von 181") for k, v in lf.items()},
                              "Zufall generisch (P_phys erhalten, Nahtparität gerades N)": [sum(x[0] for x in zuf["generisch"]),
                                                                                           sum(x[1] for x in zuf["generisch"])],
                              "Zufall physikalisch": [sum(x[0] for x in zuf["physikalisch"]), sum(x[1] for x in zuf["physikalisch"])],
                              "physikalisch: max |[U, (-1)^Q_A]| auf ungeradem N": zahl(max(x[2] for x in zuf["physikalisch"]), 4),
                              "z125 Tabelle: Paaruhr (1,3) erzwingt, Produktuhr (1,1) leer": z125_ok}
check(f"[E] Gleichwertigkeit: Exakt ist P_phys P_N = (1 + (-1)^(N/2) (-1)^Q_A)/2 P_N für gerades N und 0 für ungerades (Rest "
      f"{klein(ident)}): Für teilchenzahlerhaltende Kompositionen heißt 'physikalischen Raum erhalten' genau 'Nahtparität auf "
      f"den geraden Sektoren erhalten'. z125s Familie (181 Winkel): P_phys erhalten bei theta/pi = {lf['P_phys erhalten']}, "
      f"Nahtparität bei {lf['[U, (-1)^Q_A] = 0 auf geradem N']}, Paaruhr bei {lf['[U, Paaruhr] = 0 auf geradem N']}; die "
      f"Produktuhr-Projektion als Kontrolle bei {len(lf['Produktuhr-Projektion erhalten (Kontrolle)'])} von 181 (fordert nichts; "
      f"z125-Tabelle: Paaruhr erzwingt, Produktuhr leer: {z125_ok}). Zufällige N-erhaltende Unitäre: generisch erhalten "
      f"{sum(x[0] for x in zuf['generisch'])}/20 den Raum und {sum(x[1] for x in zuf['generisch'])}/20 die Nahtparität, "
      f"blockphysikalische {sum(x[0] for x in zuf['physikalisch'])}/20 und {sum(x[1] for x in zuf['physikalisch'])}/20 "
      f"(Äquivalenz in allen 40: {aequi}); auf den ungeraden, unphysikalischen Sektoren bleibt die Nahtparität frei (max "
      f"{max(x[2] for x in zuf['physikalisch']):.3f})",
      ident < 1e-12 and lf["P_phys erhalten"] == [0.0, 0.5, 1.0] and lf["[U, (-1)^Q_A] = 0 auf geradem N"] == [0.0, 0.5, 1.0]
      and lf["[U, Paaruhr] = 0 auf geradem N"] == [0.0, 0.5, 1.0] and len(lf["Produktuhr-Projektion erhalten (Kontrolle)"]) == 181
      and z125_ok and aequi and sum(x[0] for x in zuf["generisch"]) == 0 and sum(x[0] for x in zuf["physikalisch"]) == 20
      and max(x[2] for x in zuf["physikalisch"]) > 0.1)

# 2.9 Besetzung der Sektoren: TFD gegen Produkt zweier unabhängiger Nähte
tfd_aus = {k: zahl(nrm(np.diag((QA != QB).astype(float)) @ SF[k]["Om"])**2) for k in GERADE}
p_prod = np.exp(-4*np.pi*(QA**2 + QB**2)/2.0)
s20 = int(np.where((QA == 2) & (QB == 0))[0][0])
s11 = int(np.where((QA == 1) & (QB == 1))[0][0])
verh = float(p_prod[s20]/p_prod[s11])
ZAHLEN["besetzung"] = {"TFD: Gewicht auf Q_A != Q_B": tfd_aus, "Produkt-KMS (H_Q, beta = 4 pi): (2,0)/(1,1) je Zustand": zahl(verh, 6)}
check(f"[E] Besetzung: Der KMS-Vektor der Standardform besetzt keinen Sektor mit Q_A != Q_B (Gewicht {klein(max(tfd_aus.values()))} "
      f"für alle vier geraden Zustände), der Produkt-KMS-Zustand zweier unabhängiger Nähte (H_Q, beta = 4 pi) dagegen den "
      f"Doppelsektor mit {verh:.3e} = e^(-4 pi) je Zustand relativ zu (1,1) (z125 aus z93). z125s Forderung 'Leckfreiheit auf "
      "allen Sektoren' gehört zur Lesart zweier unabhängiger Nähte; in der Original-Partner-Lesart, in der die Paaruhr "
      "Delta^(-i/2) ist, liegt der native Zustand im physikalischen Raum und der Doppelsektor außerhalb",
      max(tfd_aus.values()) < 1e-24 and abs(verh - np.exp(-4*np.pi)) < 1e-15)
befund("E/C: In der Standardform einer Naht mit ihrem Tomita-Partner erzeugt kein Tomita-Objekt einen physikalischen "
       "Einzeltransfer. Modulgruppe (gerade für gerade Zustände), J (Tausch der Seiten, gegenchiral), M' = JMJ und der "
       "gradierte Kommutant Z M' Z^* lassen den Einzeltransfer als Operator zu (m' m in (M v M')^gerade), die KMS-Relationen "
       "koppeln Naht und Partner durch Kreuzpaare a^+ f^+ (paaruhrneutral, nahtparitätsungerade). Was ausschließt, ist die "
       "2 pi-Drehung der verdoppelten Naht: Delta^(-i/2) = i^(Q_A - Q_B) = R R-quer^-1 = e^(2 pi i (L_0 - L-quer_0)) (mit der "
       "Viertelverschiebung bei beta = 4 pi), und z113s native Projektion macht sie auf physikalischen Zuständen zu 1. Der "
       "Einzeltransfer dreht sie um -1 und hat im physikalischen Raum kein Matrixelement. Für N-erhaltende Kompositionen ist "
       "das exakt P_univ auf den geraden Sektoren. Die native Form von P_univ ist die Niveaugleichheit (Paaruhr); die "
       "Nahtparitätsform ist ihr N-erhaltender Spezialfall und ginge als strikte Superauswahl zu weit (sie verböte das "
       "KMS-Kreuzpaar und die Purifikation selbst)", True)

# =================================================================== Teil 3
say("\nTeil 3: Folgen")
# 3.1 Kette zu P_sym
i1 = np.where(NT == 1)[0]
A1 = (-1j*R_PAAR)[np.ix_(i1, i1)]
anti = {}
for th in (np.pi/3, np.pi/2, 2*np.pi/3):
    Bt = expm(-1j*th*JX[np.ix_(i1, i1)])
    B1 = Bt.conj().T @ A1 @ Bt
    anti[round(float(th/np.pi), 4)] = zahl(nrm(A1 @ B1 + B1 @ A1 - 2*np.cos(th)*np.eye(len(i1))))
anti_pi2 = nrm(A1 @ (expm(-1j*np.pi/2*JX[np.ix_(i1, i1)]).conj().T @ A1 @ expm(-1j*np.pi/2*JX[np.ix_(i1, i1)]))
               + (expm(-1j*np.pi/2*JX[np.ix_(i1, i1)]).conj().T @ A1 @ expm(-1j*np.pi/2*JX[np.ix_(i1, i1)])) @ A1)
gso_n1 = nrm(P_PHYS[np.ix_(i1, i1)])
ZAHLEN["kette_psym"] = {"|{A, B} - 2 cos theta| (theta/pi)": anti, "|{A, B}| bei pi/2": zahl(anti_pi2),
                        "|P_phys auf N = 1|": zahl(gso_n1), "j(pi/2)": zahl(j_marken(np.pi/2), 10)}
check(f"[E] Kette: Erhaltung des physikalischen Raums erzwingt in z125s Familie theta = pi/2 (Teil 2), und auf dem "
      f"Einanregungssektor ist A = -iR, B = B_theta^+ A B_theta mit {{A, B}} = 2 cos theta (Reste {anti}), bei pi/2 also "
      f"{{A, B}} = 0 ({klein(anti_pi2)}): P_sym mit tau = i, j = {j_marken(np.pi/2):.1f}. Der Einanregungssektor selbst ist "
      f"gesamtungerade und liegt außerhalb von P_phys (|P_phys| dort {klein(gso_n1)}): z94s Binäroperationen wirken auf dem "
      "nackten Fermionsektor, entschieden wird auf dem physikalischen Sektor N = 2 (dieselbe Transferdrehung)",
      max(anti.values()) < 1e-12 and anti_pi2 < 1e-12 and gso_n1 < 1e-12 and abs(j_marken(np.pi/2) - 1728) < 1e-6)

# 3.2 Verträglichkeit mit z107
z107 = json.loads((HERE / "z107_ergebnis.json").read_text())
j_z107 = float(z107["zahlen"]["psym_minimalitaet"]["j(Marken) theta = pi/3"])
leck_pi3 = leck(gatter(np.pi/3))
ohne = max(leck(expm(-1j*phi*HQ)) for phi in (np.pi/5, np.pi/2, 2.0))
ZAHLEN["z107"] = {"j(pi/3) z107": j_z107, "j(pi/3) hier": zahl(j_marken(np.pi/3), 10), "Leck theta = pi/3": zahl(leck_pi3, 10),
                  "Leck ohne Komposition (theta = 0, drei phi)": zahl(ohne)}
check(f"[E] Verträglichkeit mit z107: Der spannungsfreie Fall theta = pi/3 (A* gilt, j = {j_z107:.1f} laut z107, hier "
      f"{j_marken(np.pi/3):.1f}) verlässt den physikalischen Raum (Leck {leck_pi3:.6f} = |cos sin||1 - e^(-i pi/2)|), und ohne "
      f"Komposition (theta = 0) bleibt jedes Phasengesetz im Raum (Leck {klein(ohne)}): A* allein legt theta nicht fest, P_sym "
      "braucht weiterhin die Existenz einer nichttrivialen Komposition",
      abs(j_z107 - j_marken(np.pi/3)) < 1e-6 and abs(leck_pi3 - np.cos(np.pi/3)*np.sin(np.pi/3)*np.sqrt(2)) < 1e-12 and ohne < 1e-12)

# 3.3 Verträglichkeit mit z94
phi_tab = {}
for phi in (np.pi/5, np.pi/3, np.pi/2, 2.0, np.pi):
    phi_tab[f"{phi:.4f}"] = [round(float(th/np.pi), 4) for th in TH181 if leck(gatter(th, phi)) < 1e-10]
pclock = {k: any(u.startswith(("Paaruhr", "zwei Viertel")) for u in v) for k, v in treffer.items()}   # Delta^(-i/2) = R R-quer^-1
ZAHLEN["z94"] = {"leckfrei je phi": phi_tab, "P_clock (Delta^(-i/2) = R R-quer^-1) je Zustand": pclock}
check(f"[E] Verträglichkeit mit z94: Die leckfreie Menge ist für phi = pi/5, pi/3, pi/2, 2, pi dieselbe "
      f"({sorted({tuple(v) for v in phi_tab.values()})}), und P_phys ist durch R R-quer^-1 definiert, also beta-blind: Die "
      f"Komposition legt theta fest, nicht beta. P_clock (Delta^(-i/2) = R R-quer^-1 mit der Uhr R des jeweiligen "
      f"Ladungsgitters) gilt umgekehrt unabhängig von theta genau für die Zustände mit Viertelverschiebung "
      f"({ {k: v for k, v in pclock.items()} }). P_sym und P_clock bleiben "
      "unabhängig; P_univ in der Niveauform braucht P_clock nicht (es benutzt R R-quer^-1 selbst), nur seine Gleichsetzung mit "
      "Delta^(-i/2)",
      all(v == [0.0, 0.5, 1.0] for v in phi_tab.values()) and pclock["pi Q (z125)"]
      and pclock["4 pi E_lad(1/4, Q) = 2 pi Q^2 - pi Q"] and pclock["zwei Viertel (eps = pi, 3 pi)"]
      and not pclock["4 pi Q^2/2 (ohne Viertelverschiebung)"])

# 3.4 Lesart zweier gleichchiraler Nähte und K19
gleich = {"Strom E^+_0: Produktuhrklasse": KLASS["Strom E^+_0 (z99, Vertex e^(i(phi_A - phi_B)))"]["Produktuhr"],
          "Gatterscan Gesamtdrehung": scan["L_A L_B (Gesamtdrehung) x (x) 1"],
          "Gatterscan kanonische Parität je verdoppelter Naht": scan["(-1)^(F_A + F_A-quer) x (x) 1"],
          "Gatterscan Niveau je Naht": scan["L_A je Naht x (x) 1"]}
ZAHLEN["gleichchiral"] = {k: (v if not isinstance(v, list) or len(v) < 6 else f"{len(v)} von 37") for k, v in gleich.items()}
check(f"[C] Lesart zweier gleichchiraler Nähte (z99s Rahmen, zwei Hüllen, beide holomorph): die 2 pi-Drehung des Komposits ist "
      f"die Produktuhr, z99s Strom ist dort uhrneutral ({gleich['Strom E^+_0: Produktuhrklasse']}, h = 1, lokal), und die native "
      f"Gesamtdrehung fordert nichts (37 von 37 Winkeln, Teil 1.3). P_univ ist dort gleichwertig zur Univalenz je verdoppelter "
      f"Naht, (-1)^(F_A + F_A-quer) auch unter Kopplung superselektiert (theta/pi = "
      f"{gleich['Gatterscan kanonische Parität je verdoppelter Naht']}); die stärkere Niveaugleichheit je Naht ließe nur "
      f"{gleich['Gatterscan Niveau je Naht']} zu. In dieser Lesart bleibt P_univ eine eigene Prämisse (gescopt negativ). K19: "
      "In der Original-Partner-Lesart ist der zweite Kompositionspartner der GNS-Partner jedes treuen Zustands (Teil 2 rechnet "
      "ihn aus dem Zustand einer Naht); gebraucht wird eine Naht mit KMS-Zustand, nicht zwei. Ihre Existenz bleibt L",
      gleich["Strom E^+_0: Produktuhrklasse"] == "rein k = 0" and len(gleich["Gatterscan Gesamtdrehung"]) == 37
      and gleich["Gatterscan kanonische Parität je verdoppelter Naht"] == [0.0, 0.5, 1.0] and gleich["Gatterscan Niveau je Naht"] == [0.0, 1.0])
befund("C: P_univ ist in der Original-Partner-Lesart, in der z125 die Paaruhr als Delta^(-i/2) liest, keine eigene Prämisse: Es "
       "ist die Niveaugleichheit (Univalenz der 2 pi-Drehung) der verdoppelten Naht, die z113 aus der Modulinvarianz nativ "
       "gibt, angewandt auf N-erhaltende Kompositionen. Der reduzierte Ursprung ist A* + Existenz einer nichttrivialen "
       "kohärenten Komposition der Naht mit ihrem Tomita-Partner + Rahmen (Hardy-Naht, Viertelzählung, CAR, Modulinvarianz "
       "der verdoppelten Naht). Für zwei gleichchirale Nähte ist die 2 pi-Drehung die Produktuhr, P_univ bleibt dort eine "
       "Prämisse (Univalenz je verdoppelter Naht)", True)

# =================================================================== Bilanz
say("\nTeil 4: Bilanz")
STATUS = {
    "Klassifikation der Erzeuger": ("E", "z99-Ströme E^+-_n, J^+- und Vertex e^(i(phi_A - phi_B)) rein (-,-) = Einzeltransfers; "
                                         "Kreuzpaar (-,-) paaruhrneutral; Cartan, Doppeltransfer, virtuelle Transfers, U_pi/2 "
                                         "rein (+,+); U_theta, B_theta, Gamma(pi) gemischt"),
    "virtuelle Transfers": ("E", "nahtgerade, weil (-1)^Q_A = e^(i pi (N/2 + J_z)) in der Transfer-U(2) liegt (Twirl strukturell)"),
    "J-Paarung (verdoppelt)": ("E", "x J x J erhält (-1)^(F_A + F_A-quer) genau für homogene x; Gesamtdrehung fordert nichts; "
                                    "Niveau je Naht verbietet jede Komposition"),
    "Tomita-Objekte": ("E", "Delta^it gerade (Univalenz des Zustands), J Tausch und gegenchiral, J a J = (-1)^Q_A b, Theta = Z^* J "
                            "gradiert; Einzeltransfer = m' m in (M v M')^gerade: Kommutant/Koszul schließen nicht aus"),
    "KMS": ("E", "a Omega = q f^+ Omega (Kegel) bzw. a Omega_g = q b^+ Omega_g (gradiert); Omega = e^G|0> mit Kreuzpaar-"
                 "Erzeuger G, paaruhrneutral, nahtparitätsungerade"),
    "2 pi-Drehung": ("E", "Delta^(-i/2) = i^(Q_A - Q_B) = R R-quer^-1 = e^(2 pi i (L_0 - L-quer_0)) mit Viertelverschiebung bei "
                          "beta = 4 pi (P_clock); ohne sie nur (-1)^F"),
    "Einzeltransfer Naht -> Partner": ("E", "kein Matrixelement im physikalischen Raum (z113: R = R-quer, (-1)^(F+F-quer) = +1)"),
    "P_univ (Original-Partner-Lesart)": ("E unter z113 und der Lesart", "= Erhaltung des physikalischen Raums für N-erhaltende "
                                                                        "Kompositionen; native Form: Niveaugleichheit"),
    "P_univ (zwei gleichchirale Nähte)": ("C (gescopt negativ)", "Gesamtdrehung = Produktuhr fordert nichts; minimales "
                                                                 "Zusatzdatum: Univalenz je verdoppelter Naht"),
    "Verträglichkeit z107, z94": ("E", "theta = pi/3 spannungsfrei, aber undicht; leckfreie Menge phi- und beta-blind; P_clock "
                                       "unabhängig von theta"),
    "K19 Existenz der Naht": ("L", "Kompositionspartner = GNS-Partner; eine Naht mit treuem KMS-Zustand genügt, ihre Existenz "
                                   "bleibt unhergeleitet"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
PRAEMISSEN = {
    "A* (Nahtprämisse)": "(a) S-Invarianz des Inhalts, (b) <T> = 0; gibt P_tn/P_clock (z107) und S(n) (z125)",
    "Existenz einer nichttrivialen kohärenten Komposition": "die Naht tauscht mit ihrem Tomita-Partner Anregungen kohärent aus "
                                                            "(U auf dem physikalischen Raum nicht diagonal)",
    "Lesart P_TP (benannt)": "Kompositionspartner ist der Tomita-Partner (gegenchiral); dieselbe Lesart, in der z125 die "
                             "Paaruhr als Delta^(-i/2) liest",
    "Rahmen": "Hardy-Naht, Viertelzählung, CAR, Modulinvarianz der verdoppelten Naht (z113: R = R-quer, (-1)^(F+F-quer) = +1)",
    "entfällt": "P_univ als eigene Superauswahlprämisse (in der Lesart P_TP)",
}
ZAHLEN["praemissen"] = PRAEMISSEN
check("[C] Auftragskriterien: Klassifikation aller genannten nativen Zwischennahtoperatoren mit J-Bildern (Teil 1), Test der "
      "Tomita-Quellen (Modulgruppe, J, M' = JMJ, KMS) auf einen Einzeltransfer mit Kontrollen (Teil 2), Folgen für Ursprung, "
      "z107, z94 und K19 sowie die Gegenlesart (Teil 3)",
      all(v for k, v in CHECKS.items()))

urteil = ("P_univ ist in der Original-Partner-Lesart nativ, für zwei gleichchirale Nähte nicht. (1) Klassifikation: Die nativen "
          "Erzeuger der Zwischennahtkopplung sind Einzelfermiontransfers. z99s Ströme E^+-_n (der Vertex e^(+-i(phi_A - "
          "phi_B))) und die Transferströme J^+- sind rein (-,-) in den Nahtparitäten, das Kreuzpaar psi*_A psi*_B ebenso. Rein "
          "(+,+) sind die Cartanströme, z105s Doppeltransfer, die virtuellen Transfers (Superaustausch und Twirl; beim Twirl, "
          "weil (-1)^Q_A = e^(i pi (N/2 + J_z)) selbst in der Transfer-U(2) liegt) und das Sandwich U_pi/2; U_theta für theta "
          "außerhalb {0, pi/2, pi}, B_theta und Gamma(pi) sind gemischt. V5s wörtliche Prüfung fällt damit auf der Ebene der "
          "Erzeuger negativ aus. Mit dem J-Bild gepaart (x J x J, verdoppeltes Modell mit zwei Partnern) erhält ein Operator die "
          "kanonische Parität (-1)^(F_A + F_A-quer) jeder verdoppelten Naht genau dann, wenn er homogen ist. (2) Standardform: "
          "In der GNS-Darstellung des KMS-Zustands bei beta = 4 pi (S aus der Definition, Delta, J) erzeugt kein Tomita-Objekt "
          "einen physikalischen Einzeltransfer, aber nicht aus dem vermuteten Grund. Die Modulgruppe ist für gerade Zustände je "
          "Seite gerade und ladungserhaltend (ein paritätsbrechender Zustand bricht das). J vertauscht die Seitenparitäten, "
          "erhält die Paaruhr und invertiert die Produktuhr: Der Partner ist gegenchiral. J a_c J = (-1)^Q_A b_c liegt im "
          "gewöhnlichen Kommutanten, Theta = Z^* J bildet a_c auf den gradierten Partner ab, Z vertauscht beide. Kommutant, "
          "Koszul-Vorzeichen und Verdrehung lassen den Einzeltransfer als Operator zu: T = b^+ a = m' m liegt in "
          "(M v M')^gerade (128 gegen 64 Dimensionen). Die KMS-Relation a_c Omega = q f_c^+ Omega (natürlicher Kegel; in z122s "
          "gradierter Form für die Paarprodukt-Purifikation) macht den Transfer auf Omega zu einer Partner-Paarerzeugung, und "
          "Omega selbst entsteht durch den Kreuzpaar-Erzeuger G = sum theta (a^+ f^+ - f a): nahtparitätsungerade, aber "
          "paaruhrneutral. Ausschlaggebend ist die 2 pi-Drehung der verdoppelten Naht. Mit der Viertelverschiebung bei beta = "
          "4 pi ist Delta^(-i/2) = i^(Q_A - Q_B) = R R-quer^-1 = e^(2 pi i (L_0 - L-quer_0)) (P_clock; ohne Verschiebung nur "
          "(-1)^F). z113s native Projektion (Modulinvarianz der nichtchiralen Naht: R = R-quer, (-1)^(F+F-quer) = +1) macht "
          "sie auf physikalischen Zuständen zu 1. Omega liegt in diesem Raum, alle Tomita-Objekte erhalten ihn, und der "
          "Einzeltransfer Naht -> Partner hat dort kein Matrixelement (|P T P| = 0, |T P| = 1: Spin 1/2 mod 1). Für "
          "N-erhaltende Kompositionen ist die Erhaltung des Raums exakt die Nahtparität auf den geraden Sektoren "
          "(P_phys P_N = (1 + (-1)^(N/2) (-1)^Q_A)/2 P_N); z125s Familie ist genau bei theta/pi in {0, 0.5, 1} physikalisch, "
          "die Produktuhr als Kontrolle fordert nichts. Die native Form von P_univ ist die Niveaugleichheit (Paaruhr); die "
          "Nahtparitätsform ist ihr N-erhaltender Spezialfall und ginge als strikte Superauswahl zu weit (sie verböte das "
          "KMS-Kreuzpaar und damit die Purifikation, die beide Nahtparitäten kohärent überlagert). z125s Forderung 'Leckfreiheit auf allen Sektoren' stammt aus der Lesart zweier unabhängiger Nähte "
          "(Doppelsektor e^(-4 pi)); der KMS-Vektor der Standardform besetzt Q_A != Q_B nicht. (3) Folgen: In der "
          "Original-Partner-Lesart (P_TP, dieselbe, in der z125 die Paaruhr als Delta^(-i/2) liest) reduziert sich der "
          "Ursprung auf A* + Existenz einer nichttrivialen kohärenten Komposition der Naht mit ihrem Tomita-Partner + Rahmen "
          "(Hardy-Naht, Viertelzählung, CAR, Modulinvarianz der verdoppelten Naht); P_univ entfällt als eigene Prämisse. "
          "Verträglich mit z107 (theta = pi/3 spannungsfrei, j = 3905.8, aber undicht 0.612; ohne Komposition keine "
          "theta-Bedingung) und z94 (die leckfreie Menge ist phi- und beta-blind, P_clock gilt unabhängig von theta und nur mit "
          "Viertelverschiebung). K19: Der zweite Kompositionspartner ist der GNS-Partner; gebraucht wird eine Naht mit treuem "
          "KMS-Zustand, nicht zwei, ihre Existenz bleibt L. Für zwei gleichchirale Nähte (z99s Rahmen) ist die 2 pi-Drehung die "
          "Produktuhr, Univalenz fordert dann nichts, z99s Ströme sind physikalisch (h = 1), und P_univ ist dort gleichwertig zur "
          "Univalenz je verdoppelter Naht ((-1)^(F_A + F_A-quer) auch unter Kopplung superselektiert), dem minimalen "
          "Zusatzdatum; die Niveaugleichheit je Naht verböte jede Komposition. Nichts befördert, Ledger unverändert.")
ZAHLEN["laufzeit unter 3 min"] = time.time() - T_START < 180
say("\nUrteil\n   " + urteil)
say(f"\nLaufzeit unter 3 min: {ZAHLEN['laufzeit unter 3 min']}")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z130_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z130.txt").write_text("\n".join(OUT) + "\n", encoding="utf-8")
