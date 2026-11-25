"""Z87: Ursprung der Saat (Strategie S6): vom Qubit zum Original.

Auftrag S6 (GESAMTSTAND-UND-STRATEGIEN.md): Folgt die ganze Einteilchenstruktur des Originals aus der Saat rho_1 mit nativen
Regeln? Ist die zweite Quantisierung die einzige Vervollständigung mit Uhr und KMS? Was fixiert die Lücke pi der Saat?
Leitlinie 'Rekonstruktion statt Auswahl': keine neuen Auswahlprinzipien; jeder Schritt ist gerechnete Folge oder benannte
Prämisse. Gelöste Bausteine werden nicht geöffnet, nichts wird befördert.
Bestand: z81 (R = i^(Q_a - Q_b), R^2 = (-1)^F, Paarkovarianz rho_m), z63 (rho_m je Windung, Sym^m, Fock-Ecke V_1), z74
(Uhrenquotient m | 4), z80 Teil 3 (Mischung verschobener Vakua) und Teil 8 (Cl_2, tau = i, Spurzustand), z71 Teil 3,
z46 (q = e^(-2 pi)), z48 (k = 4x, J_-1 = S^4), z58 (h = 4 pi x, gaußsche Antworten), BOUNDARY-SOURCE (Calderón-Halbannulus:
C_A = (1 + e^(2 pi B_cov))^-1, B_cov = 2(B_1/4 + B_3/4)), ELLIPTIC-CURRENT-SOURCE (Gaußquantisierung als Prämisse).

Teil 1  Saatfamilie: Die Paarkovarianz auf (x, -x), x = m/4, ist exakt rho_1^m/Tr; die ganze Kovarianz folgt aus rho_1.
        Die Quartermoden sind die NS-Moden der Überlagerung (x = B_cov/2), R ist ihre Halbdrehung mit R^2 = -1, ungerade m
        heißt NS. Alternativen (nichtlineare Dispersion, gerade m, andere beta, andere Zwiste) gegen fünf native Kriterien:
        KMS bezüglich der Drehuhr, Spin-Z4, Uhrenidentität Delta^(-i/2) = e^(2 pi i L_0), Hardy/Szegő-Drehkovarianz,
        Uhrenquotient (z74).
Teil 2  Statistik: Die Z4-Graduierung allein trennt CAR und CCR nicht (beide Fockräume tragen R^2 = (-1)^N, R^4 = 1).
        Was CAR erzwingt: binäre Saat (Cl_2 = CAR eines Modus), Basisunabhängigkeit der Saat (Paulionen scheitern),
        KMS-Positivität bei negativen Modulenergien (für jeden Zustand), Lokalität mit Saatnorm auf der Überlagerung
        (Spin-Statistik).
Teil 3  Eindeutigkeit der Vervollständigung: KMS für die freie Dynamik ist eindeutig. Ohne Dynamik gibt es ungleiche
        Vervollständigungen mit denselben Einteilchendaten (Mischung verschobener Vakua, z80; treue Diagonaldeformation).
        Die native Uhrenidentität Delta^(-i/2) = R (z79/z81) schließt beide aus und erzwingt, mit Ladungssuperauswahl, in allen
        gerechneten Abschneidungen den quasifreien Zustand: Gitterspektrum, exakte Ladungsidentitäten über Z[y] (e^pi
        transzendent), Subadditivität.
Teil 4  Die Saat selbst: Cl_2 aus der Vorzeichensymmetrie (z80); gerade Zustände (1 + b Gamma)/2; ihr Modulfluss dreht die
        beiden Binäroperationen; Viertelschritt bei t = -1/2 genau für Lücke g in pi(4Z + 1). Kette tau = i -> q = e^(-2 pi)
        -> beta_angle = 2 pi -> beta = 4 pi -> Lücke pi über den Calderón-Halbannulus; Status jedes Schritts.
Teil 5  Bilanz und Urteil für S6: Prämissenzählung, Erfolgs- und Abbruchkriterien, K19.
"""
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.linalg as sla
import sympy as sp

HERE = Path(__file__).parent
ROOT = HERE.parents[2]
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
mp.mp.dps = 40
rng = np.random.default_rng(87)
PI4 = 4*np.pi
YN = np.exp(-np.pi)                                                    # y = e^(-pi), Saatgewicht


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def c_orig(x):
    return 1/(1 + np.exp(PI4*x))


def proj_ordnung(U, nmax=12):
    P = np.eye(U.shape[0], dtype=complex)
    for k in range(1, nmax + 1):
        P = P @ U
        if np.allclose(P, P[0, 0]*np.eye(U.shape[0]), atol=1e-10):
            return k
    return None


def jw(n):
    """Jordan-Wigner-Vernichter, Index 0 je Faktor = leer, 1 = besetzt."""
    I2, Zp, sm = np.eye(2), np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [0.0, 0.0]])
    ops = []
    for j in range(n):
        mats = [Zp]*j + [sm] + [I2]*(n - j - 1)
        op = mats[0]
        for M in mats[1:]:
            op = np.kron(op, M)
        ops.append(op)
    return ops


def entropie(rho):
    w = np.linalg.eigvalsh((rho + rho.conj().T)/2)
    w = w[w > 1e-300]
    return float(-np.sum(w*np.log(w)))


# =================================================================== Teil 1
say("Teil 1: Die Saatfamilie aus rho_1 und nativen Regeln")

# --- 1.1 Paarkovarianz = rho_1^m/Tr (exakt) und die ganze Kovarianz aus rho_1
ok11 = True
r1 = sp.diag(1, sp.exp(-sp.pi))
for m in range(1, 17, 2):
    x = sp.Rational(m, 4)
    paar = sp.diag(1/(1 + sp.exp(-4*sp.pi*x)), 1/(1 + sp.exp(4*sp.pi*x)))   # diag(c(-x), c(x))
    pot = r1**m
    pot = pot/pot.trace()
    ok11 &= sp.simplify(paar - pot) == sp.zeros(2, 2) and sp.simplify(paar.trace() - 1) == 0


def c_aus_saat(x):
    """Nur aus rho_1: Besetzungsgewicht e^(-pi |4x|) relativ zum leeren Zustand; Löcher über c(-x) = 1 - c(x)."""
    w = YN**abs(int(round(4*x)))
    return w/(1 + w) if x > 0 else 1/(1 + w)


xs_all = [r + a for r in range(-12, 12) for a in (0.25, 0.75)]
dev12 = max(abs(c_aus_saat(x) - c_orig(x))/max(c_orig(x), 1e-300) for x in xs_all)
check("auf jedem Fermipaar (x, -x), x = m/4, m = 1, 3, .., 15 ist die Paarkovarianz diag(c(-x), c(x)) exakt rho_1^m/Tr(rho_1^m) "
      f"mit Spur 1 (symbolisch); die ganze Kovarianz c(x) auf 48 Moden folgt allein aus rho_1 (rel. Abweichung {dev12:.0e})",
      ok11 and dev12 < 1e-12)


def sym_ecke(rho, m):
    big = rho
    for _ in range(m - 1):
        big = np.kron(big, rho)
    B = np.zeros((2**m, m + 1))
    for idx in range(2**m):
        B[idx, bin(idx).count("1")] = 1
    B = B/np.linalg.norm(B, axis=0)
    S = B.T @ big @ B
    e = np.diag([S[0, 0], S[m, m]])
    return e/np.trace(e)


rho1 = np.diag([1.0, YN])/(1 + YN)
gibbs = lambda g: np.diag([1.0, np.exp(-g)])/(1 + np.exp(-g))
check("Sym^m rho_1, eingeschränkt auf die äußere Ecke {|0..0>, |1..1>}, ist rho_m für m = 1..7 (z63 hatte m = 2, 4): Potenz und "
      "symmetrische Potenz der Saat geben dieselbe Familie",
      all(np.allclose(sym_ecke(rho1, m), gibbs(np.pi*m), rtol=1e-12, atol=0) for m in range(1, 8)))

# --- 1.2 NS-Überlagerung
ok13 = True
for n in range(-14, 14):
    bc = n + 0.5                                                        # B_cov-Eigenwert (NS)
    x = bc/2                                                            # Quarterenergie
    ok13 &= abs(1/(1 + np.exp(2*np.pi*bc)) - c_orig(x)) < 1e-15        # beta_angle = 2 pi <-> beta = 4 pi
    ok13 &= abs(np.exp(1j*np.pi*bc) - 1j**int(round(4*x))) < 1e-12     # V = T_pi = R
    ok13 &= abs(np.exp(2j*np.pi*bc) + 1) < 1e-12                        # V^2 = -1
quarter = sorted({round(((n + 0.5)/2) % 1, 9) for n in range(-14, 14)})
ramond = [(n, n/2, int(4*n/2) % 2, np.exp(2j*np.pi*n)) for n in range(-3, 4)]
ok_R = all(m4 == 0 and abs(v - 1) < 1e-12 for n, x, m4, v in ramond) and abs(c_orig(0.0) - 0.5) < 1e-15
check("die Quartermoden sind die NS-Moden der Überlagerung: x = B_cov/2 mit B_cov in Z + 1/2 läuft genau über die Zwiste "
      f"{quarter}, c_cov = 1/(1 + e^(2 pi B_cov)) = c(x), und die Halbdrehung V = e^(i pi B_cov) (BOUNDARY-SOURCE: V = T_pi) ist "
      "R = i^(4x) mit V^2 = -1; Ramond (B_cov ganz) gäbe gerade 4x, V^2 = +1 und eine halb besetzte Nullmode",
      ok13 and quarter == [0.25, 0.75] and ok_R)
befund("Ungerade m ist dieselbe Aussage wie NS auf der Überlagerung und wie R^2 = -1 auf Fermionen: die Saatfamilie ist das "
       "NS-Fermion der quadratischen Überlagerung, aufgeteilt nach den Eigenwerten i, -i der Halbdrehung", True)

# --- 1.3 Alternativen gegen fünf native Kriterien


def moden(alphas, R=8):
    return sorted({r + a for a in alphas for r in range(-R, R)})


def kriterien(h, xs):
    hv = np.array([h(x) for x in xs])
    nicht_null = [x for x in xs if abs(x) > 1e-12]
    nullmode = any(abs(x) < 1e-12 for x in xs)
    bet = np.array([h(x)/x for x in nicht_null])
    K = (not nullmode) and np.ptp(bet) < 1e-9*abs(bet).max()
    Rv = np.exp(1j*hv/2)
    Z = bool(np.allclose(Rv**4, 1, atol=1e-9) and np.allclose(Rv**2, -1, atol=1e-9))
    C = bool(np.allclose(Rv, np.exp(2j*np.pi*np.array(xs)), atol=1e-9))
    sek = defaultdict(list)
    for x in xs:
        sek[round(x % 1, 9)].append(x)
    H = all(np.ptp([h(x + 1) - h(x) for x in v]) < 1e-9 for v in sek.values())
    xmin = min(x for x in xs if h(x) > 1e-12)
    meff = h(xmin)/np.pi
    ganz = abs(meff - round(meff)) < 1e-9 and round(meff) in (1, 2, 4)
    ordn = proj_ordnung(np.diag([1, np.exp(1j*h(xmin)/2)]))
    U = bool(ganz and ordn == 4)
    return {"K": bool(K), "Z": Z, "C": C, "H": bool(H), "U": U, "m_eff": round(float(meff), 6), "Eckordnung": ordn}


lin = lambda b: (lambda x: b*np.pi*x)
alternativen = {
    "Original 4 pi x, alpha 1/4, 3/4": (lin(4), moden((0.25, 0.75)), set()),
    "nichtlinear 4 pi x (1 + x^2/10)": (lambda x: PI4*x*(1 + 0.1*x*x), moden((0.25, 0.75)), set("KZCHU")),
    "gitterkubisch pi (4x)^3": (lambda x: np.pi*(4*x)**3, moden((0.25, 0.75)), set("KH")),
    "massiv 4 pi sign(x) sqrt(x^2 + 1/4)": (lambda x: PI4*np.sign(x)*np.sqrt(x*x + 0.25), moden((0.25, 0.75)), set("KZCHU")),
    "gerade m: alpha = 1/2": (lin(4), moden((0.5,)), set("ZU")),
    "gerade m: alpha = 0 (Nullmode)": (lin(4), moden((0.0,)), set("KZU")),
    "Zwist alpha = 1/8, 7/8": (lin(4), moden((0.125, 0.875)), set("ZU")),
    "Zwist alpha = 1/3, 2/3": (lin(4), moden((1/3, 2/3)), set("ZU")),
    "beta = 2 pi": (lin(2), moden((0.25, 0.75)), set("ZCU")),
    "beta = 8 pi": (lin(8), moden((0.25, 0.75)), set("ZCU")),
    "beta = 12 pi": (lin(12), moden((0.25, 0.75)), set("CU")),
    "beta = 20 pi": (lin(20), moden((0.25, 0.75)), set("U")),
    "beta = 4,04 pi": (lin(4.04), moden((0.25, 0.75)), set("ZCU")),
}
tab1, ok_tab = {}, True
for name, (h, xs, soll) in alternativen.items():
    kr = kriterien(h, xs)
    fail = {k for k in "KZCHU" if not kr[k]}
    tab1[name] = {**kr, "scheitert_an": "".join(sorted(fail)) or "-"}
    ok_tab &= fail == soll
    say(f"      {name:38s} K {kr['K']!s:5} Z {kr['Z']!s:5} C {kr['C']!s:5} H {kr['H']!s:5} U {kr['U']!s:5} "
        f"m_eff {kr['m_eff']}, Eckordnung {kr['Eckordnung']}")
ZAHLEN["teil1_alternativen"] = tab1
ZAHLEN["teil1_kriterien"] = {"K": "jedes Fermipaar KMS bezüglich derselben Drehuhr L_0 (h(x)/x konstant, keine Nullmode)",
                             "Z": "Spin-Z4: R = e^(ih/2) mit R^4 = 1, R^2 = -1 auf allen Moden",
                             "C": "Uhrenidentität Delta^(-i/2) = e^(2 pi i L_0) auf allen Moden",
                             "H": "Hardy/Szegő: [h, S] = omega S (Modulfluss = Drehung der Naht)",
                             "U": "Uhrenquotient: m_eff = h(x_min)/pi teilt 4, Fock-Ecke mit projektiver R-Ordnung 4"}
check("nur das Original erfüllt alle fünf Kriterien; jede der 12 Alternativen scheitert an genau der erwarteten Menge: nichtlinear "
      "und massiv an allen, gerade m und andere Zwiste an Z und U (alpha = 0 zusätzlich an K: Nullmode ohne Zeit), beta = 2, 8, "
      "4,04 pi an Z, C, U, beta = 12 pi an C, U (Delta^(-i/2) = R^-1), beta = 20 pi nur am Uhrenquotienten (m_eff = 5), "
      "gitterkubisch nur an K und H", ok_tab)
befund("Die Saatregel rho_m = rho_1^m/Tr braucht keine freie Wahl über die benannten Kriterien hinaus: lineare Dispersion ist "
       "die Drehkovarianz (H, Hardy/Szegő-Rahmen), ungerade m ist Spin-Z4 (Z), beta = 4 pi ist Uhrenidentität (C = Prämisse "
       "P_clock, Teil 4) plus Uhrenquotient (U, z74) bzw. tau = i (Teil 4: Im tau = 5 gäbe beta = 20 pi). Die Z4-Bedingungen "
       "allein lassen die gitterkubische Dispersion und beta = 20 pi zu", True)

# =================================================================== Teil 2
say("\nTeil 2: Statistik -- erzwingt der Spin-Z4 die CAR?")
sx, sy_, sz = (np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]), np.diag([1.0 + 0j, -1.0]))
I2 = np.eye(2)
A, B = sx, sy_
psi = (A + 1j*B)/2
Gam = -1j*A @ B
check("Saatebene: die zwei Binäroperationen A, B (A^2 = B^2 = 1, AB + BA = 0, z80) sind zwei Majoranas eines Fermimodus: "
      "psi = (A + iB)/2 hat psi^2 = 0, {psi, psi^+} = 1, und Gamma = -iAB = 1 - 2 psi^+ psi ist (-1)^F; Cl_2 = M_2 = CAR(C)",
      np.allclose(psi @ psi, 0) and np.allclose(psi @ psi.conj().T + psi.conj().T @ psi, I2)
      and np.allclose(Gam, I2 - 2*psi.conj().T @ psi) and np.allclose(A @ B + B @ A, 0))


def boson(nc):
    return np.diag(np.sqrt(np.arange(1, nc)), 1)


spread = {nc: float(np.linalg.eigvalsh(boson(nc) + boson(nc).T).max()) for nc in (4, 8, 16, 32)}
ZAHLEN["ccr_feld_spektrum_max"] = spread
check("ein CCR-Feld a + a^+ ist keine Binäroperation: sein Spektrum wächst mit der Abschneidung (größter Eigenwert "
      f"{', '.join(f'{v:.2f}' for v in spread.values())} bei 4, 8, 16, 32 Stufen), A^2 = 1 ist unmöglich",
      all(b > a for a, b in zip(list(spread.values()), list(spread.values())[1:])) and spread[32] > 9)

car4 = jw(4)
sm, Zp = np.array([[0.0, 1.0], [0.0, 0.0]]), np.diag([1.0, -1.0])
pau = [np.kron(np.kron(np.kron(sm, I2), I2), I2), np.kron(np.kron(np.kron(I2, sm), I2), I2),
       np.kron(np.kron(np.kron(I2, I2), sm), I2), np.kron(np.kron(np.kron(I2, I2), I2), sm)]
dev_bu = {}
for name, ops in (("CAR", car4), ("Paulionen", pau)):
    w2 = wa = 0.0
    for _ in range(30):
        f = rng.normal(size=4) + 1j*rng.normal(size=4)
        f /= np.linalg.norm(f)
        P = sum(fk*o for fk, o in zip(f, ops))
        w2 = max(w2, np.abs(P @ P).max())
        wa = max(wa, np.abs(P @ P.conj().T + P.conj().T @ P - np.eye(16)).max())
    dev_bu[name] = (float(w2), float(wa))
ZAHLEN["basisunabhaengigkeit_psi(f)"] = dev_bu
check("Basisunabhängigkeit der Saat (jede normierte Einteilchenrichtung f trägt dieselbe Cl_2, psi(f)^2 = 0, {psi(f), psi(f)^+} = 1; "
      f"30 zufällige f in C^4): CAR erfüllt sie (Abweichung {max(dev_bu['CAR']):.0e}), Paulionen mit vertauschenden Moden nicht "
      f"(Abweichung {max(dev_bu['Paulionen']):.2f})", max(dev_bu["CAR"]) < 1e-14 and min(dev_bu["Paulionen"]) > 0.1)

# Z4-Graduierung auf beiden Fockräumen
X4 = [0.25, -0.75, 0.75, -0.25]
occ4 = np.array([np.diag(c.T @ c) for c in car4]).T.round().astype(int)
R_car = np.diag(np.exp(2j*np.pi*(occ4 @ np.array(X4))))
par_car = np.diag((-1.0)**occ4.sum(axis=1))
nc = 10
a = boson(nc)
a1, a2 = np.kron(a, np.eye(nc)), np.kron(np.eye(nc), a)
n1, n2 = np.diag(a1.T @ a1).round().astype(int), np.diag(a2.T @ a2).round().astype(int)
R_ccr = np.diag(np.exp(2j*np.pi*(0.25*n1 - 0.25*n2)))
par_ccr = np.diag((-1.0)**(n1 + n2))
ok_grad = (np.allclose(R_car @ R_car, par_car) and np.allclose(np.linalg.matrix_power(R_car, 4), np.eye(16))
           and np.allclose(R_ccr @ R_ccr, par_ccr) and np.allclose(np.linalg.matrix_power(R_ccr, 4), np.eye(nc*nc)))
grad_car = np.allclose(car4[0] @ car4[2] + car4[2] @ car4[0], 0)
grad_ccr = np.allclose(a1 @ a2 - a2 @ a1, 0) and not np.allclose(a1 @ a2 + a2 @ a1, 0)
check("die Z4-Graduierung R = Gamma(R_1) mit R_1 = e^(2 pi i x) existiert auf beiden Fockräumen: CAR (4 Moden) und CCR (Moden "
      "x = 1/4, -1/4, 10 Stufen) haben R^2 = (-1)^N und R^4 = 1; aber nur bei CAR ist R^2 die Statistikgraduierung (ungerade "
      "Erzeuger antivertauschen), bei CCR vertauschen sie", ok_grad and grad_car and grad_ccr)
befund("R^2 = (-1)^F allein wählt die Statistik nicht, wenn (-1)^F nur als Zahlparität gelesen wird; als Statistikgraduierung "
       "gelesen ist es genau die Spin-Statistik-Beziehung der Überlagerung und verlangt CAR", True)

hh, nn = sp.symbols("h n", real=True)
nB = sp.solve(sp.Eq(nn, sp.exp(-hh)*(nn + 1)), nn)[0]
nF = sp.solve(sp.Eq(nn, sp.exp(-hh)*(1 - nn)), nn)[0]
paarB, paarF = sp.simplify(nB + nB.subs(hh, -hh)), sp.simplify(nF + nF.subs(hh, -hh))
neg = [x for x in xs_all if x < 0]
nB_neg = [float(nB.subs(hh, sp.Float(PI4*x, 30))) for x in neg]
Zb = {nc_: float(np.sum(np.exp(np.pi*np.arange(nc_)))) for nc_ in (6, 10, 14)}
ZAHLEN["ccr_kms"] = {"<a^+a> bei x = -1/4": float(nB.subs(hh, -sp.pi)), "Gibbsnorm h = -pi, Stufen 6/10/14": Zb,
                     "Paarspur CAR": str(paarF), "Paarspur CCR": str(paarB)}
check("KMS-Positivität (für jeden Zustand, nicht nur quasifrei): aus KMS und CCR folgt <a^+a> = 1/(e^h - 1), aus KMS und CAR "
      f"1/(e^h + 1); auf allen {len(neg)} Moden mit x < 0 des Originals (h < 0) wäre der CCR-Wert negativ (x = -1/4: "
      f"{float(nB.subs(hh, -sp.pi)):.4f}), die Gibbsnorm eines Bosemodus mit h = -pi divergiert mit der Abschneidung; die "
      f"Paarspur ist {paarF} (CAR, Zustand rho_m) gegen {paarB} (CCR, kein Zustand)",
      sp.simplify(nB - 1/(sp.exp(hh) - 1)) == 0 and sp.simplify(nF - 1/(sp.exp(hh) + 1)) == 0 and paarF == 1 and paarB == -1
      and max(nB_neg) < 0 and Zb[14]/Zb[10] > 1e5)

ker = {}
for u in (0.5, 1.0, 2.0, 3.0):
    eps = mp.mpf("1e-7")
    S = mp.exp((1j*u - eps)/2)/(1 - mp.exp(1j*u - eps))                # Summe_(n >= 0) e^((iu - eps)(n + 1/2))
    KF, KB = 2*mp.re(S), 2*mp.im(S)                                     # Summe e^(iru), Summe sign(r) e^(iru) / i
    ker[u] = (float(KF), float(KB), float(1/mp.sin(u/2)))
ZAHLEN["spin_statistik_kerne"] = {str(u): {"CAR {psi(0), psi(u)}": v[0], "CCR [phi(0), phi(u)]/i": v[1], "1/sin(u/2)": v[2]}
                                  for u, v in ker.items()}
check("Spin-Statistik auf der Überlagerung (reelles Feld, NS-Moden r in Z + 1/2, Saatnorm 1 je Mode, Abel-Regularisierung): "
      "der CAR-Antikommutator ist lokal (Kern -> 0 bei u = 0,5, 1, 2, 3), der CCR-Kommutator mit derselben Norm und positiver "
      "Energie ist i/sin(u/2), also nichtlokal; lokale CCR verlangt die Norm r je Mode (Spin 1), also keine binäre Saat",
      all(abs(v[0]) < 1e-5 and abs(v[1] - v[2]) < 1e-5*v[2] for v in ker.values()))
befund("Die Statistik ist erzwungen: CAR folgt aus der binären Saat (Cl_2 = CAR eines Modus), aus der Basisunabhängigkeit der "
       "Saat, aus der KMS-Positivität bei den negativen Modulenergien des Originals und aus der Lokalität mit Saatnorm. R^2 = "
       "(-1)^F ist dann die Spin-Statistik-Beziehung der NS-Überlagerung, kein unabhängiger Wahlschritt", True)

# =================================================================== Teil 3
say("\nTeil 3: Eindeutigkeit der Vervollständigung")
SEK4 = [0, 0, 1, 1]                                                    # Quarter a = {1/4, -3/4}, b = {3/4, -1/4}
C4 = car4
N4 = [c.T @ c for c in C4]
dim4 = 16
sea4 = (np.array(X4) < 0).astype(float)
h4 = PI4*np.array(X4)
K0d = (occ4 - sea4) @ h4
rho_qf = np.diag(np.exp(-K0d))
rho_qf /= np.trace(rho_qf)
c4 = c_orig(np.array(X4))
Qa = np.diag(((occ4 - sea4)*(np.array(SEK4) == 0)).sum(axis=1))
Qb = np.diag(((occ4 - sea4)*(np.array(SEK4) == 1)).sum(axis=1))
R4 = np.diag(1j**(np.diag(Qa) - np.diag(Qb)).round().astype(int))


def kov(rho):
    return np.array([[np.trace(rho @ C4[i].T @ C4[j]) for j in range(4)] for i in range(4)])


cov0 = kov(rho_qf)
check("Prüfsystem: zwei Fermipaare (m = 1, 3), 16 Fockzustände; rho_qf = e^(-K_0)/Z ist das Produkt der Besetzungen c(x), "
      "und die Uhr R = e^(iK_0/2) ist genau i^(Q_a - Q_b) (z81)",
      np.allclose(np.diag(cov0).real, c4, rtol=1e-12) and np.allclose(cov0 - np.diag(np.diag(cov0)), 0)
      and np.allclose(np.diag(np.exp(1j*K0d/2)), R4))


# --- 3.1 KMS für die freie Dynamik ist eindeutig
def kms_matrix(H, beta=1.0):
    n = H.shape[0]
    E, Ei = sla.expm(-beta*H), sla.expm(beta*H)
    rows = []
    for i, j, k, l in itertools.product(range(n), repeat=4):
        Am = np.zeros((n, n))
        Am[i, j] = 1
        Bm = np.zeros((n, n))
        Bm[k, l] = 1
        rows.append((Am @ E @ Bm @ Ei - Bm @ Am).T.ravel())
    return np.array(rows, dtype=complex)


X3 = [0.25, -0.25, 0.75]
C3 = jw(3)
occ3 = np.array([np.diag(c.T @ c) for c in C3]).T.round().astype(int)
K3 = np.diag((occ3 - (np.array(X3) < 0)) @ (PI4*np.array(X3)))
_, s3, vh3 = np.linalg.svd(kms_matrix(K3))
nul3 = vh3[int(np.sum(s3 > 1e-10*s3.max())):]
rho3 = nul3[0].conj().reshape(8, 8) if len(nul3) == 1 else np.zeros((8, 8))
rho3 = rho3/np.trace(rho3)
prod3 = np.diag(np.prod(np.where(occ3 == 1, c_orig(np.array(X3)), 1 - c_orig(np.array(X3))), axis=1))
check("KMS (beta = 1) für die freie Dynamik e^(itK_0), K_0 = dGamma(4 pi x), auf dem ganzen Fockraum (3 Moden, alle 4096 "
      "Matrixeinheitenpaare): genau eine Lösung, und sie ist der quasifreie Produktzustand mit c(x)",
      len(nul3) == 1 and np.allclose(rho3, prod3, atol=1e-12))

# --- 3.2 Mischung verschobener Vakua (z80) bei beta = 4 pi
Xs = [r + 0.25 for r in range(-3, 4)]                                  # eine Quarterkette
Cs = jw(len(Xs))
occs = np.array([np.diag(c.T @ c) for c in Cs]).T.round().astype(int)
cs = c_orig(np.array(Xs))
cbs = c_orig(-np.array(Xs))                                            # 1 - c(x) = c(-x), ohne Auslöschung
P_mix = np.zeros(2**len(Xs))
for e in range(-1, len(Xs)):
    muster = (np.arange(len(Xs)) <= e).astype(int)
    idx = int(np.where((occs == muster).all(axis=1))[0][0])
    P_mix[idx] = cbs[0] if e == -1 else (cs[-1] if e == len(Xs) - 1 else cs[e] - cs[e + 1])
P_qfs = np.prod(np.where(occs == 1, cs, cbs), axis=1)
nb = lambda P, j: float(P @ occs[:, j])
nn2 = lambda P, i, j: float(P @ (occs[:, i]*occs[:, j]))
leer2 = lambda P, i, j: float(P @ ((1 - occs[:, i])*(1 - occs[:, j])))
i_, j_ = 2, 3                                                          # x = -3/4, y = 1/4 (über der Fermikante)
hx, hy = PI4*Xs[i_], PI4*Xs[j_]
kms_qf = abs(np.exp(-(hx + hy))*leer2(P_qfs, i_, j_) - nn2(P_qfs, i_, j_))/nn2(P_qfs, i_, j_)
kms_mix = abs(np.exp(-(hx + hy))*leer2(P_mix, i_, j_) - nn2(P_mix, i_, j_))/nn2(P_mix, i_, j_)
cov_mix = nn2(P_mix, i_, j_) - cs[i_]*cs[j_]
ZAHLEN["mischung_beta4pi"] = {"rang": int(np.sum(P_mix > 0)), "dim": 2**len(Xs), "KMS_rel_abweichung_qf": kms_qf,
                             "KMS_rel_abweichung_mix": kms_mix, "Cov(n_-3/4, n_1/4)": cov_mix}
check(f"Mischung verschobener Vakua (z80) bei beta = 4 pi auf einer Quarterkette (7 Moden): dieselben c(x) "
      f"(Abweichung {max(abs(nb(P_mix, j) - cs[j]) for j in range(len(Xs))):.0e}), aber Rang {int(np.sum(P_mix > 0))} von "
      f"{2**len(Xs)} (nicht treu, kein Modulfluss), Cov(n_-3/4, n_1/4) = c(1/4)(1 - c(-3/4)) = {cov_mix:.2e} statt 0 und "
      f"Verletzung der Zweiteilchen-KMS-Relation um {kms_mix:.1%} (quasifrei {kms_qf:.0e})",
      max(abs(nb(P_mix, j) - cs[j]) for j in range(len(Xs))) < 1e-14 and int(np.sum(P_mix > 0)) == len(Xs) + 1
      and abs(cov_mix - cs[j_]*cbs[i_]) < 1e-15 and kms_mix > 0.01 and kms_qf < 1e-12)

# --- 3.3 treue Diagonaldeformation
lam = 0.5
Dd = (np.diag(N4[0]) - c4[0])*(np.diag(N4[3]) - c4[3])                 # (n_1/4 - c)(n_-1/4 - c'), zwei Quarter
rho_l = rho_qf @ np.diag(1 + lam*Dd)
covl = kov(rho_l)
treu = np.linalg.eigvalsh(rho_l).min() > 0
komm = all(np.allclose(rho_l @ O, O @ rho_l) for O in (Qa, Qb, np.diag(K0d), R4))
cov_nn = float(np.trace(rho_l @ N4[0] @ N4[3]).real - c4[0]*c4[3])
soll_nn = lam*c4[0]*(1 - c4[0])*c4[3]*(1 - c4[3])
psi_a, n_c = C4[0], N4[3]
lhs = np.exp(-h4[0])*np.trace(rho_l @ psi_a @ n_c @ psi_a.T).real      # KMS: w(A sigma_i(B)) = w(BA), A = psi, B = n psi^+
rhs = np.trace(rho_l @ n_c @ psi_a.T @ psi_a).real
uhr_l = sla.expm(-0.5j*sla.logm(rho_l)) @ np.linalg.inv(R4)
uhr_q = sla.expm(-0.5j*sla.logm(rho_qf)) @ np.linalg.inv(R4)
check(f"treue Diagonaldeformation rho_lambda = rho_qf (1 + lambda D), D = (n_1/4 - c)(n_-1/4 - c'), lambda = {lam}: dieselbe "
      f"volle Kovarianz (Abweichung {np.abs(covl - cov0).max():.0e}), treu, vertauscht mit Q_a, Q_b, K_0 und R (stationär "
      f"unter der freien Uhr); aber Cov(n_1/4, n_-1/4) = {cov_nn:.3e} = lambda c(1-c)c'(1-c') (Wick verletzt) und die KMS-"
      f"Relation für psi_1/4, n_-1/4 psi_1/4^+ gilt nicht ({lhs:.4e} gegen {rhs:.4e})",
      np.abs(covl - cov0).max() < 1e-15 and treu and komm and abs(cov_nn - soll_nn) < 1e-15 and abs(lhs - rhs) > 1e-4*abs(rhs))
check("die Uhrenidentität Delta^(-i/2) = R (z79/z81) für den eigenen Modulfluss: rho_qf^(-i/2) R^-1 ist ein Vielfaches der Eins, "
      "rho_lambda^(-i/2) R^-1 nicht (die Deformation verstimmt die Uhr)",
      np.allclose(uhr_q, uhr_q[0, 0]*np.eye(dim4), atol=1e-10) and not np.allclose(uhr_l, uhr_l[0, 0]*np.eye(dim4), atol=1e-6))


def takesaki_rest(rho, op, t=0.3):
    L = sla.logm(rho)
    st = sla.expm(1j*t*L) @ op @ sla.expm(-1j*t*L)
    coef = np.trace(op.conj().T @ st)/np.trace(op.conj().T @ op)
    return float(np.linalg.norm(st - coef*op)/np.linalg.norm(op))


tk_q, tk_l = takesaki_rest(rho_qf, C4[0]), takesaki_rest(rho_l, C4[0])
S_qf, S_l = entropie(rho_qf), entropie(rho_l)
S_marg = float(sum(-(c*np.log(c) + (1 - c)*np.log(1 - c)) for c in c4))
ZAHLEN["entropie"] = {"S(rho_qf)": S_qf, "Summe Einmodenentropien": S_marg, "S(rho_lambda)": S_l}
check(f"Takesaki-Form der Gaußprämisse: unter dem Modulfluss des Ganzen bleibt die Saatalgebra des Modus 1/4 invariant für "
      f"rho_qf (Rest {tk_q:.0e}), nicht für rho_lambda (Rest {tk_l:.3f}); Entropie: S(rho_qf) = Summe der Einmodenentropien "
      f"= {S_qf:.6f}, S(rho_lambda) = {S_l:.6f} kleiner (Subadditivität)",
      tk_q < 1e-12 and tk_l > 1e-3 and abs(S_qf - S_marg) < 1e-12 and S_l < S_qf - 1e-6)

# --- 3.4 Starrheit unter der Uhrenidentität
yS = sp.symbols("y")


def aufbau(ms):
    X = sorted([sp.Rational(m, 4) for m in ms] + [-sp.Rational(m, 4) for m in ms])
    sek = [0 if (x % 1) == sp.Rational(1, 4) else 1 for x in X]
    blocks = defaultdict(list)
    for p in itertools.product((0, 1), repeat=len(X)):
        Em = sum(int(4*x)*(pk - (1 if x < 0 else 0)) for x, pk in zip(X, p))
        qa = sum((pk - (1 if x < 0 else 0)) for x, pk, s in zip(X, p, sek) if s == 0)
        qb = sum((pk - (1 if x < 0 else 0)) for x, pk, s in zip(X, p, sek) if s == 1)
        blocks[(qa, qb)].append(Em)
    for k in blocks:
        blocks[k].sort()
    cy = lambda x: yS**int(4*x)/(1 + yS**int(4*x)) if x > 0 else 1/(1 + yS**int(-4*x))
    Qa_ = sp.together(sum(cy(x) - (1 if x < 0 else 0) for x, s in zip(X, sek) if s == 0))
    Qb_ = sp.together(sum(cy(x) - (1 if x < 0 else 0) for x, s in zip(X, sek) if s == 1))
    Na, Da = sp.fraction(sp.factor(Qa_))
    Nb, Db = sp.fraction(sp.factor(Qb_))
    Dg = sp.lcm(Da, Db)
    Na, Nb, Dg = sp.expand(Na*sp.cancel(Dg/Da)), sp.expand(Nb*sp.cancel(Dg/Db)), sp.expand(Dg)
    return dict(blocks), Dg, Na, Nb


def koeff(expr, L=32):
    arr = np.zeros(L, dtype=np.int64)
    for (e,), cf in sp.Poly(sp.expand(expr), yS).terms():
        arr[e] += int(cf)
    return arr


def ladungsvektor(Dv, Nav, Nbv, key, W):
    qa, qb = key
    return np.concatenate([qa*np.convolve(W, Dv) - np.convolve(W, Nav), qb*np.convolve(W, Dv) - np.convolve(W, Nbv)])


def suche_zwei_paare(ms, fenster=range(-8, 25), OFF=16, L=72):
    blocks, Dg, Na, Nb = aufbau(ms)
    Dv, Nav, Nbv = koeff(Dg), koeff(Na), koeff(Nb)
    keys = sorted(blocks)

    def beitrag(key, exps):
        W = np.zeros(L, dtype=np.int64)
        for e in exps:
            W[e + OFF] += 1
        return ladungsvektor(Dv, Nav, Nbv, key, W)
    orig_ok = not sum(beitrag(k, blocks[k]) for k in keys).any()
    cands = {k: [list(c) for c in itertools.combinations_with_replacement(
        [e for e in fenster if e % 4 == (k[0] - k[1]) % 4], len(blocks[k]))] for k in keys}
    pre = {k: [beitrag(k, c) for c in cands[k]] for k in keys}
    g1, g2, p1, p2 = [], [], 1, 1
    for k in sorted(keys, key=lambda k: -len(cands[k])):
        if p1 <= p2:
            g1.append(k)
            p1 *= len(cands[k])
        else:
            g2.append(k)
            p2 *= len(cands[k])
    H = defaultdict(list)
    for idx in itertools.product(*[range(len(cands[k])) for k in g2]):
        H[sum(pre[k][j] for k, j in zip(g2, idx)).tobytes()].append(idx)
    loes = []
    for idx in itertools.product(*[range(len(cands[k])) for k in g1]):
        v = sum(pre[k][j] for k, j in zip(g1, idx))
        for idx2 in H.get((-v).tobytes(), []):
            s = {k: cands[k][j] for k, j in zip(g1, idx)}
            s.update({k: cands[k][j] for k, j in zip(g2, idx2)})
            loes.append(s)

    def verschiebung(s):
        d = {e for k in keys for e in (np.array(s[k]) - np.array(blocks[k])).tolist()}
        return d.pop() if len(d) == 1 else None
    vs = [verschiebung(s) for s in loes]
    return {"D": str(Dg), "Kandidaten": p1*p2, "Loesungen": len(loes),
            "globale_Verschiebungen": sorted(int(v) for v in vs if v is not None),
            "nichttrivial": sum(v is None for v in vs), "original_erfuellt": orig_ok}


starr = {}
for ms in ((1, 3), (1, 5), (3, 5), (1, 7)):
    starr[str(ms)] = suche_zwei_paare(ms)
    r = starr[str(ms)]
    say(f"      Paare m = {ms}: D = {r['D']}, {r['Kandidaten']:.2e} Blockspektren, Lösungen {r['Loesungen']}, davon globale "
        f"Verschiebungen y^4k mit 4k in {r['globale_Verschiebungen']}, nichttrivial {r['nichttrivial']}")
ZAHLEN["starrheit_zwei_paare"] = starr
check("Starrheit (zwei Fermipaare, vier Abschneidungen m = (1,3), (1,5), (3,5), (1,7), Exponenten -8..24): Unter der "
      "Uhrenidentität hat K' auf jedem Ladungsblock ein Spektrum in pi (Q_a - Q_b) + 4 pi Z, die Blockgewichte sind also "
      "Summen von y-Potenzen; da e^pi transzendent ist (Gelfond), müssen die Ladungsmittel <Q_a>, <Q_b> als Identität über Z[y] "
      "stimmen. Alle Lösungen sind globale Verschiebungen (derselbe Zustand): blockisospektral",
      all(r["nichttrivial"] == 0 and r["original_erfuellt"] and 0 in r["globale_Verschiebungen"] for r in starr.values()))

# drei Paare: kleine Umordnungen
blocks3, D3, Na3, Nb3 = aufbau((1, 3, 5))
Dv3, Nav3, Nbv3 = koeff(D3), koeff(Na3), koeff(Nb3)
OFF3, L3 = 16, 80
zuege = []
for key, es in blocks3.items():
    for e in sorted(set(es)):
        for k in (-3, -2, -1, 1, 2, 3):
            W = np.zeros(L3, dtype=np.int64)
            W[e + 4*k + OFF3] += 1
            W[e + OFF3] -= 1
            zuege.append(((key, e, k), ladungsvektor(Dv3, Nav3, Nbv3, key, W), W))


def gueltig(kombi):
    cnt = Counter((z[0][0], z[0][1]) for z in kombi)
    return all(n_ <= blocks3[key].count(e) for (key, e), n_ in cnt.items())


def klasse(kombi):
    neu = {k: Counter(v) for k, v in blocks3.items()}
    for (key, e, k), _, _ in kombi:
        neu[key][e] -= 1
        neu[key][e + 4*k] += 1
    if any(min(c.values()) < 0 for c in neu.values()):
        return None
    if all(+neu[k] == Counter(blocks3[k]) for k in blocks3):
        return "trivial"
    return "T=0" if not sum(z[2] for z in kombi).any() else "T!=0"


paare3 = defaultdict(list)
for i in range(len(zuege)):
    for j in range(i, len(zuege)):
        paare3[(zuege[i][1] + zuege[j][1]).tobytes()].append((i, j))
einzel3 = defaultdict(list)
for i, z in enumerate(zuege):
    einzel3[z[1].tobytes()].append(i)
gefunden = set()
for i, z in enumerate(zuege):
    if not z[1].any():
        gefunden.add((i,))
for kb, lst in paare3.items():
    v = np.frombuffer(kb, dtype=np.int64)
    if not v.any():
        gefunden |= {tuple(sorted(p)) for p in lst}
    for i3 in einzel3.get((-v).tobytes(), []):
        gefunden |= {tuple(sorted(p + (i3,))) for p in lst}
    for p2 in paare3.get((-v).tobytes(), []):
        gefunden |= {tuple(sorted(p + p2)) for p in lst}
klassen = Counter()
beispiel = None
for kombi_idx in gefunden:
    kombi = [zuege[i] for i in kombi_idx]
    if not gueltig(kombi):
        continue
    kl = klasse(kombi)
    if kl is None:
        continue
    klassen[kl] += 1
    if kl == "T=0" and beispiel is None and len(kombi) == 4:
        beispiel = [z[0] for z in kombi]
ZAHLEN["starrheit_drei_paare_kleine_umordnungen"] = {"Zuege": len(zuege), "Klassen": dict(klassen),
                                                    "Beispiel_T=0": str(beispiel)}
check(f"drei Fermipaare (m = 1, 3, 5; 16 Ladungsblöcke, 64 Zustände), alle Umordnungen aus bis zu vier Eigenwertzügen um 4k, "
      f"|k| <= 3 ({len(zuege)} Einzelzüge): Lösungen der Ladungsidentitäten {dict(klassen)}; jede nichttriviale Lösung hat T = "
      "Summe der Gewichtsänderungen = 0, ist also global isospektral; keine Lösung mit T != 0",
      klassen.get("T!=0", 0) == 0 and klassen.get("T=0", 0) > 0)
say(f"      Beispiel einer global isospektralen Ladungsverschiebung: {beispiel}")

# Subadditivitätsschritt numerisch: isospektrale Drehungen im Ladungskommutanten ändern die Besetzungen
bl = defaultdict(list)
for s_ in range(dim4):
    bl[(int(round(Qa[s_, s_].real)), int(round(Qb[s_, s_].real)))].append(s_)
abw = []
for _ in range(20):
    U = np.zeros((dim4, dim4), dtype=complex)
    for idx in bl.values():
        Zc = rng.normal(size=(len(idx), len(idx))) + 1j*rng.normal(size=(len(idx), len(idx)))
        Qm, Rm = np.linalg.qr(Zc)
        U[np.ix_(idx, idx)] = Qm*(np.diag(Rm)/np.abs(np.diag(Rm)))
    rU = U @ rho_qf @ U.conj().T
    abw.append(float(np.abs(np.diag(kov(rU)).real - c4).max()))
check("Subadditivitätsschritt: S(rho') = S(rho_qf) zusammen mit denselben Einmodenbesetzungen erzwingt rho' = rho_qf; 20 "
      f"zufällige isospektrale Drehungen im Ladungskommutanten verschieben die Besetzungen stets (kleinste Abweichung "
      f"{min(abw):.1e})", min(abw) > 1e-8)

quelle = (ROOT / "experiments/theory-contracts/native-arc-charge-source-interface-20261026/ELLIPTIC-CURRENT-SOURCE.txt")
text = quelle.read_text() if quelle.exists() else ""
zeile = "The free Gaussian quantization is an EXPLICIT existing premise."
befund(f"Die gaußsche Struktur aus z58 (Wick, Jacobi-Tripelprodukt) ist eine Ausgabe des quasifreien Zustands; die Quelle "
       f"nennt sie ausdrücklich Prämisse ('{zeile[:58]}...'). Hier: KMS für die freie Dynamik erzwingt sie (Standard); ohne "
       "Dynamik gibt es ungleiche Vervollständigungen; die native Uhrenidentität Delta^(-i/2) = R erzwingt sie zusammen mit den "
       "Einteilchendaten und der Ladungssuperauswahl in allen gerechneten Abschneidungen. Allgemein (beliebig viele Paare, "
       "beliebige Exponenten) ist das nicht bewiesen", zeile in text)

# =================================================================== Teil 4
say("\nTeil 4: Die Saat selbst -- was fixiert die Lücke pi?")
thetas = np.sort(np.append(np.linspace(0.05, np.pi - 0.05, 60), np.pi/2))
impl = []
for th_ in thetas:
    Bm = np.cos(th_)*sz + np.sin(th_)*sx
    Lm = []
    for i in range(2):
        for j in range(2):
            Eu = np.zeros((2, 2), dtype=complex)
            Eu[i, j] = 1
            Lm.append(np.concatenate([(Eu @ sz + sz @ Eu).ravel(), (Eu @ Bm - Bm @ Eu).ravel()]))
    if np.sum(np.linalg.svd(np.array(Lm).T, compute_uv=False) < 1e-10) > 0:
        impl.append(float(th_))
check(f"Vorzeichensymmetrie (z80 Teil 8, nachgerechnet): A -> -A ist in derselben irreduziblen 2D-Darstellung nur bei theta = "
      f"pi/2 implementierbar ({len(thetas)} Winkel, lösbar bei {[round(t_/np.pi, 4) for t_ in impl]} pi): AB + BA = 0",
      len(impl) == 1 and abs(impl[0] - np.pi/2) < 1e-12)


def modfluss(rho, X, t):
    Lg = sla.logm(rho)
    return sla.expm(1j*t*Lg) @ X @ sla.expm(-1j*t*Lg)


komm_G = [M for M in (I2, sx, sy_, sz) if np.allclose(M @ Gam, Gam @ M)]


def modfluss_g(g, X, t):
    """rho_g = e^(g Gamma/2)/Z, also rho_g^(it) X rho_g^(-it) = e^(itg Gamma/2) X e^(-itg Gamma/2), ohne logm."""
    U = sla.expm(0.5j*t*g*Gam)
    return U @ X @ U.conj().T


lsg, tab4 = [], {}
for g in np.linspace(0.25, 12, 48)*np.pi:
    sA, sB = modfluss_g(g, A, -0.5), modfluss_g(g, B, -0.5)
    if np.allclose(sA, B, atol=1e-9) and np.allclose(sB, -A, atol=1e-9):
        lsg.append(float(round(g/np.pi, 6)))
for g in (0.5, 1.0, 2.0, 3.0):
    sA = modfluss_g(g*np.pi, A, -0.5)
    tab4[f"g = {g} pi"] = [round(float(np.trace(sA @ P).real/2), 6) for P in (A, B)]
ok_rot = True
for g in (0.37, 1.0, 2.6):
    rho_g = (I2 + np.tanh(g/2)*Gam)/2
    for t in (0.1, 0.7):
        ok_rot &= np.allclose(modfluss(rho_g, A, t), np.cos(g*t)*A - np.sin(g*t)*B, atol=1e-10)
        ok_rot &= np.allclose(modfluss(rho_g, A, t), modfluss_g(g, A, t), atol=1e-10)
ZAHLEN["modfluss_sigma_-1/2(A)_in_(A,B)"] = tab4
check("gerade Zustände (vertauschen mit Gamma = (-1)^F, also mit der gemeinsamen Spiegelung A, B -> -A, -B) sind genau "
      f"(1 + b Gamma)/2 (Kommutant der Dimension {len(komm_G)}); ihr Modulfluss dreht die Binäroperationen: sigma_t(A) = "
      "cos(gt) A - sin(gt) B mit g = log((1+b)/(1-b)); der Viertelschritt sigma_-1/2: A -> B -> -A tritt im Raster g in (0, 12 pi] "
      f"genau bei g = {lsg} pi auf, also g in pi(4Z + 1)", len(komm_G) == 2 and ok_rot and lsg == [1.0, 5.0, 9.0])
impl_r = sla.expm(-0.5j*sla.logm((I2 + np.tanh(np.pi/2)*Gam)/2))
rho1_q = (I2 + np.tanh(np.pi/2)*Gam)/2
check("die kleinste Lösung g = pi gibt rho_1 = diag(1, e^(-pi))/(1 + e^(-pi)), und rho_1^(-i/2) ist proportional zu diag(1, i), "
      "der Vierteldrehung R auf der Fock-Ecke V_1 = {Omega, psi^+ Omega} (z79); ohne Uhrenbedingung ist jedes b in (-1, 1) ein "
      "treuer gerader Zustand mit Modulzeit: das Qubit allein hat keine Skala",
      np.allclose(np.diag(rho1_q).real, [1/(1 + YN), YN/(1 + YN)]) and abs(impl_r[1, 1]/impl_r[0, 0] - 1j) < 1e-12
      and all(np.linalg.eigvalsh((I2 + b*Gam)/2).min() > 0 for b in (0.1, 0.5, 0.99)))

marken = []
for M in (A, B):
    w_, v_ = np.linalg.eigh(M)
    marken += [complex(v_[1, k]/v_[0, k]) for k in range(2)]
kv = sorted({round(((z1 - z3)*(z2 - z4)/((z1 - z4)*(z2 - z3))).real, 9) for z1, z2, z3, z4 in itertools.permutations(marken)})
lam_l = mp.mpf(1)/2
tau = 1j*mp.ellipk(1 - lam_l)/mp.ellipk(lam_l)


def eisen(tau_, k):
    q_ = mp.exp(2j*mp.pi*tau_)
    sig = lambda n: mp.fsum(d**(k - 1) for d in range(1, n + 1) if n % d == 0)
    c_ = {4: 240, 6: -504}[k]
    return 1 + c_*mp.fsum(sig(n)*q_**n for n in range(1, 80))


E4, E6 = eisen(tau, 4), eisen(tau, 6)
jinv = 1728*E4**3/(E4**3 - E6**2)
e6_achse = {float(y_): float(mp.re(eisen(1j*y_, 6))) for y_ in (0.5, 0.8, 1.0, 1.25, 2.0)}
ZAHLEN["torus"] = {"marken": [str(np.round(z, 9)) for z in marken], "kreuzverhaeltnisse": kv, "tau": str(mp.nstr(tau, 20)),
                   "E4(i)": float(mp.re(E4)), "E6(i)": float(abs(E6)), "j": float(mp.re(jinv)), "E6(iy)": e6_achse}
check(f"die vier Marken (Eigenvektoren von A und B) sind {{1, -1, i, -i}} auf dem Äquator, senkrecht zur Energieachse Gamma; "
      f"Kreuzverhältnisse {kv}; die verzweigte Doppelüberlagerung hat tau = i K(1/2)/K(1/2) = {mp.nstr(tau, 12)}, E6(i) = "
      f"{float(abs(E6)):.0e}, E4(i) = {float(mp.re(E4)):.6f}, j = {float(mp.re(jinv)):.6f}; auf der Achse tau = iy wechselt "
      "E6 nur bei y = 1 das Vorzeichen (H6-Route <T> = 0)",
      sorted(round(abs(z), 9) for z in marken) == [1.0]*4 and kv == [-1.0, 0.5, 2.0] and abs(tau - 1j) < 1e-30
      and abs(E6) < 1e-30 and abs(jinv - 1728) < 1e-20 and e6_achse[0.5] < 0 and e6_achse[0.8] < 0 and e6_achse[1.25] > 0
      and e6_achse[2.0] > 0 and abs(e6_achse[1.0]) < 1e-25)

q_nom = mp.exp(2j*mp.pi*tau)
qJ = mp.exp(1j*mp.pi*tau)
om = 1/(mp.exp(4*mp.pi) - 1)                                           # z46: omega(J_-1 J_1)
check(f"Nomen: q = e^(2 pi i tau) = {mp.nstr(mp.re(q_nom), 12)} = e^(-2 pi) (z46), e^(i pi tau) = {mp.nstr(mp.re(qJ), 12)} = "
      "e^(-pi) ist das Saatgewicht, das Stromquant J_-1 (Hardy-Marke + 4, z48) hat q^2 = omega/(1 + omega) = e^(-4 pi) mit "
      "omega = 1/(e^(4 pi) - 1) (z46): eine Zahl, drei Lesarten",
      abs(q_nom - mp.exp(-2*mp.pi)) < 1e-30 and abs(qJ - mp.exp(-mp.pi)) < 1e-30
      and abs(om/(1 + om) - mp.exp(-4*mp.pi)) < mp.mpf(10)**-35 and abs(rho1_q[1, 1].real - float(qJ.real/(1 + qJ.real))) < 1e-15)

scan = {}
for y_ in (0.5, 1, 2, 3, 5, 9):
    passt = all(abs(1/(1 + mp.exp(2*(mp.mpf(y_)/2)*2*mp.pi*(n + mp.mpf(1)/2))) - 1/(1 + mp.exp(2*mp.pi*(n + mp.mpf(1)/2))))
                < 1e-30 for n in range(-6, 6))
    xs_q = [r + a for r in range(-5, 5) for a in (0.25, 0.75)]
    uhr = all(abs(np.exp(1j*PI4*y_*x/2) - np.exp(2j*np.pi*x)) < 1e-9 for x in xs_q)
    scan[y_] = {"Original": bool(passt), "Saatluecke/pi": y_, "Uhrenidentitaet": bool(uhr)}
ZAHLEN["calderon_scan_Im_tau"] = {str(k): v for k, v in scan.items()}
check("Calderón-Halbannulus über tau = iy (Breite W = y/2, NS-Rand B_Sigma = 2 pi (n + 1/2), C_A = (1 + e^(2W B_Sigma))^-1): das "
      "Original entsteht genau bei y = 1; die Saatlücke ist pi y, beta_angle = 2 pi y, beta = 4 pi y; die Uhrenidentität "
      f"gilt genau für y in 4Z + 1 ({[k for k, v in scan.items() if v['Uhrenidentitaet']]}): tau = i wählt y = 1",
      [k for k, v in scan.items() if v["Original"]] == [1] and [k for k, v in scan.items() if v["Uhrenidentitaet"]] == [1, 5, 9])

kette = {
    "1 zwei Binäroperationen mit Vorzeichensymmetrie -> Cl_2 = M_2 = CAR(1)": ("C", "Prämisse P_sym: Symmetrie in derselben "
                                                                               "irreduziblen Darstellung (z80)"),
    "2 Marken = Quadrat, Kreuzverhältnis -1, j = 1728, tau = i": ("E", "gerechnet aus 1"),
    "3 Doppelüberlagerung = Torus C/(Z + iZ), zwei Blätter": ("E/C", "Beziehung E, physische Lesart H6 K1"),
    "4 Zeit = zweiter Zyklus, KMS für die Drehung bei beta = L": ("C", "Prämisse P_clock/P_time; in H6 unter A* hergeleitet "
                                                                   "(K6); äquivalent Delta^(-i/2) = R"),
    "5 NS längs der Naht": ("E", "= R^2 = -1 auf Fermionen = gemeinsame Spiegelung (-1)^F des Qubits; H6 K4"),
    "6 Lücke pi = 2 pi Im(tau) x 1/2, beta = 4 pi, q = e^(-2 pi)": ("E", "aus 2, 4, 5 gerechnet"),
    "7 auf dem Qubit: Lücke pi <-> sigma_-1/2 = Viertelschritt, minimal": ("E", "g in pi(4Z+1), minimal = Uhrenquotient m = 1"),
}
ZAHLEN["kette_saat"] = kette
for k_, (s_, t_) in kette.items():
    say(f"   [{s_}] {k_}: {t_}")
befund("Die Lücke pi ist keine eigene Zahl: sie ist tau = i (P_sym) unter der Uhrenbedingung (Zeit = zweiter Torus-Zyklus, "
       "äquivalent Delta^(-i/2) = R) mit NS-Rand. Das Qubit allein hat keine Skala; die Uhrenbedingung wählt pi(4Z + 1), "
       "tau = i bzw. der Uhrenquotient wählt pi", True)
befund("Lesart (L) zu K19: Die Naht ist die Bahn des Modulflusses der Saat durch die vier Marken, der Äquator des Saat-Qubits "
       "(Szegő-Kreis, z63 Teil 5); sie existiert, sobald der Zustand die Einzelspiegelungen bricht (sonst Spurzustand ohne Zeit, "
       "z80). Das ist eine Lesart, kein Beweis der Existenz der Naht", True)

# =================================================================== Teil 5
say("\nTeil 5: Bilanz für S6")
STATUS = {
    "S6(a) Saatfamilie aus rho_1": ("E", "rho_m = rho_1^m/Tr exakt; lineare Dispersion = Drehkovarianz (Szegő), ungerade m = "
                                         "Spin-Z4 = NS-Überlagerung, beta = Uhrenidentität + Uhrenquotient; 12 Alternativen "
                                         "scheitern je an benannten Kriterien"),
    "Statistik CAR": ("E", "binäre Saat, Basisunabhängigkeit, KMS-Positivität (jeder Zustand), Lokalität mit Saatnorm; R^2 = "
                           "(-1)^F allein wählt nicht"),
    "S6(b) Eindeutigkeit mit freier Dynamik": ("E", "KMS für dGamma(h) eindeutig = quasifrei"),
    "S6(b) Eindeutigkeit ohne Dynamik": ("E (negativ)", "ungleiche Vervollständigungen mit denselben Einteilchendaten: "
                                                       "Mischung verschobener Vakua (nicht treu), treue Diagonaldeformation"),
    "S6(b) Eindeutigkeit mit Uhrenidentität": ("E (gescopt) / O", "starr in allen gerechneten Abschneidungen (zwei Paare "
                                                                  "vollständig im Fenster, drei Paare bis vier Züge); "
                                                                  "allgemein offen"),
    "S6(c) Lücke pi": ("E unter P_sym, P_clock", "tau = i, beta_angle = 2 pi, beta = 4 pi, q = e^(-2 pi), Lücke pi"),
    "K19 Existenz der Naht": ("L", "Naht = Modulbahn der Marken auf dem Saat-Qubit; Existenz nicht hergeleitet"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
PRAEMISSEN = {
    "P_sym": "zwei Binäroperationen, Vorzeichensymmetrie des Gesetzes in derselben irreduziblen Darstellung (z80) -> Cl_2, tau = i",
    "P_clock": "der Modulfluss des Zustands ist bei t = -1/2 die geometrische Vierteldrehung (Delta^(-i/2) = e^(2 pi i L_0) = R); "
               "äquivalent: Zeit = zweiter Zyklus des quadratischen Torus, KMS bei beta = L (H6 K6 unter A*)",
    "Rahmen": "Hardy/Szegő-Naht mit Drehung L_0 und Viertelzählung k = 4x (P1-Daten, H6 K1/K3), Ladungssuperauswahl je Quarter",
    "P_free (bedingt entbehrlich)": "quasifreie Vervollständigung; in allen gerechneten Abschneidungen Folge von P_clock, "
                                    "allgemein offen",
}
ZAHLEN["praemissen"] = PRAEMISSEN
ZAHLEN["praemissen_anzahl"] = {"benannt": 2, "bedingt (offen, in Abschneidungen Folge)": 1, "Rahmen": 1}
check("S6-Kriterien: Erfolg 'Original = Saat + Spin-Z4 + zweite Quantisierung' erfüllt für Familie und Statistik; "
      "Abbruch 'Familie braucht freie Wahl' nicht erfüllt; Abbruch 'mehrere ungleiche Vervollständigungen' nur ohne "
      "Uhrenidentität erfüllt, mit ihr in allen gerechneten Abschneidungen nicht",
      ok_tab and CHECKS[[k for k in CHECKS if k.startswith("Starrheit")][0]] and klassen.get("T!=0", 0) == 0)

urteil = ("S6 trägt in gescopter Form. (a) Die Saatfamilie ist ohne freie Wahl: Die Paarkovarianz auf (x, -x), x = m/4, ist exakt "
          "rho_1^m/Tr, die ganze Kovarianz folgt aus rho_1; die Quartermoden sind die NS-Moden der quadratischen Überlagerung "
          "(x = B_cov/2), R ist deren Halbdrehung mit R^2 = -1, ungerade m heißt NS. Von zwölf Alternativen scheitert jede an "
          "benannten nativen Kriterien; die Z4-Bedingungen allein ließen die gitterkubische Dispersion pi(4x)^3 und beta = 20 pi "
          "zu, erst die Szegő-Drehkovarianz und der Uhrenquotient (bzw. tau = i) schließen sie aus. (b) CAR ist erzwungen, aber nicht durch "
          "R^2 = (-1)^F allein (auch der CCR-Fockraum trägt R^2 = (-1)^N, R^4 = 1): durch die binäre Saat (Cl_2 = CAR eines "
          "Modus), die Basisunabhängigkeit der Saat, die KMS-Positivität bei negativen Modulenergien (CCR gäbe <a^+a> = -1,045 bei "
          "x = -1/4, Paarspur -1) und die Lokalität mit Saatnorm (CCR-Kern i/sin(u/2)). Die quasifreie Vervollständigung ist "
          "eindeutig, wenn der Zustand KMS für die freie Dynamik ist; ohne Dynamik gibt es ungleiche Vervollständigungen mit "
          "denselben Einteilchendaten (Mischung verschobener Vakua, nicht treu; treue Diagonaldeformation, stationär und "
          "R-invariant). Die native Uhrenidentität Delta^(-i/2) = R schließt beide aus (die Mischung hat keinen Modulfluss, die "
          "Deformation verstimmt die Uhr) und erzwingt mit Ladungssuperauswahl den "
          "quasifreien Zustand in allen gerechneten Abschneidungen (vier Zwei-Paar-Abschneidungen vollständig im "
          "Exponentenfenster, drei Paare bis vier Eigenwertzüge): Gitterspektrum, exakte Ladungsidentitäten über Z[y] wegen "
          "der Transzendenz von e^pi, dann Subadditivität. Allgemein ist das offen. (c) Die Lücke pi ist keine eigene Eingabe: "
          "Auf dem Qubit (gerade Zustände (1 + b Gamma)/2) ist sigma_-1/2 genau für g in pi(4Z + 1) der Viertelschritt der "
          "Binäroperationen, minimal g = pi gibt rho_1 und R = diag(1, i) auf V_1; geometrisch ist sie tau = i (Marken als "
          "Quadrat, j = 1728, E6(i) = 0) mit Zeit als zweitem Torus-Zyklus: beta_angle = 2 pi, beta = 4 pi, q = e^(-2 pi), "
          "Saatgewicht e^(i pi tau) = e^(-pi). Übrig bleiben zwei benannte Prämissen, P_sym (Vorzeichensymmetrie in derselben "
          "irreduziblen Darstellung) und P_clock (Delta^(-i/2) = geometrische Vierteldrehung, äquivalent KMS bei beta = L), "
          "dazu der Rahmen der Hardy-Naht mit Ladungssuperauswahl; die Gaußprämisse ist in den Abschneidungen Folge von P_clock. "
          "K19 wird nicht hergeleitet, nur verengt (Lesart: Naht = Modulbahn der Marken auf dem Saat-Qubit).")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z87_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
