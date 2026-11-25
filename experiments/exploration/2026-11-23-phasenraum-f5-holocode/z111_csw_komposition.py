"""Z111: Komposition der Nähte = MHV-Diagramme (Cachazo-Svrcek-Witten)? Tiefe Welle T2.

Leitlinie 'Rekonstruktion statt Auswahl': Die CSW-Regeln, die BCFW-Rekursion und die Parke-Taylor-Formel sind hier
VERGLEICHSZIELE, keine Prämissen. Gefragt ist, ob native TFPT-Strukturen beim Zusammensetzen zweier Nähte die CSW-Zutaten
liefern: den Propagator 1/P^2 zwischen zwei MHV-Vertizes, die Off-shell-Fortsetzung lambda_P = P|eta] und die Klebung.
Firewall DIMENSION.UPLIFT.FIREWALL.01: Die Naht ist eine chirale 2D-Struktur; 4D gibt es nur als Lesart (Naht =
Twistorgerade L_x, Lage = Raumzeitpunkt x, z61/z85) und als Kinematik freier Einteilchenzustände. Keine wechselwirkende
3+1-Schließung wird behauptet; gelöste Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt, nichts wird befördert.
Bestand: z61 (Naht = Twistorgerade, L_x und L_y schneiden sich genau für det(x - y) = 0), z85 (Fockbasis = Penroses
Elementarzustände, Vakuumklasse [1/(z1 z2)] -> 2 kappa^2/det(x - i kappa^2 1), Tomita n -> -n - 4), z89 (I_pi = eps auf der
pi-Ebene, nativ bis auf die Modulphase; Lorentzbüschel I_pi, I_omega; Zerlegung B = m_Z o (eps c^+ (x) c^+); invariante Masse
zweier Nähte |<c(1) c(2)>|^2), z91 (bosonische Komposition = freier 4D-Fockraum, keine native Fusion gibt Wechselwirkung),
z99 (Vertexoperatoren = su(2)_1-Ströme, Delta = 1 starr, verbundene 1->2-Amplitude braucht einen Zweikörperterm), z106
(gemeinsame CAR = Hardy-Komposition h_0 (+) h_1 mit einer Fermisee-Projektion, kanonischer Tausch in der U(2) der
Transferströme). Hypothese T (eine Naht trägt einen MHV-Vertex) prüft z110 parallel; hier ist der MHV-Vertex Vergleichsziel.

Bezeichnungen: lambda_i, lambda~_i in C^2, p_i = lambda_i lambda~_i^T, <ij> = det(lambda_i, lambda_j),
[ij] = det(lambda~_i, lambda~_j), P^2 = det P, also det(p_i + p_j) = <ij>[ij] (paritätssymmetrische Konvention).
Parke-Taylor A_MHV = <ab>^4/(<12><23>...<n1>), Dreipunkt-Anti-MHV = [ab]^4/([12][23][31]); BCFW-Verschiebung [i, j>:
lambda~_i -> lambda~_i + z lambda~_j, lambda_j -> lambda_j - z lambda_i. CSW: MHV-Vertizes, Propagator 1/P^2,
lambda_P = P eta~ (eta~ ein fester Referenzspinor), jede innere Linie + an einem, - am anderen Ende.

Teil 1  Vergleichsmaschine: komplexe Kinematik (Impulserhaltung, Schouten, P^2 = Summe der Paarinvarianten); BCFW aus
        Dreipunktamplituden gegen Parke-Taylor (n = 4, 5, 6), gegen die geschlossene Anti-MHV-Formel und die geschlossene
        Split-Helizitätsformel bei sechs Punkten, Unabhängigkeit von der Verschiebung, Negativkontrolle (schlechte
        Verschiebung); CSW für A(1-,2-,3-,4+,5+), die Split-Helizität A(1-,2-,3-,4+,5+,6+) und zwei weitere NMHV-Sechspunkt-
        amplituden gegen BCFW, eta-Unabhängigkeit und Negativkontrollen (Diagramm weglassen, Helizität der inneren Linie
        tauschen, je Diagramm ein anderes eta).
Teil 2  Native Komposition, je Kandidat gegen Teil 1 gerechnet: (a) Propagator: Kreuzkovarianz der gemeinsamen CAR
        (Vakuum und KMS bei beta = 4 pi), Hüpfterm der Transferströme, Szegő-Kern einer Naht an den Inzidenzpunkten der
        anderen (Penrose-Transformation, euklidische Fouriertransformation), thermischer Nahtkern gegen den Cauchykern,
        Residuen an schneidenden Nähten; (b) Off-shell-Fortsetzung: Referenzort auf I_pi, I_omega, I_AdS, I_dS (z89-Koordinaten),
        Translationen, Modulphase, rahmenabhängige native Referenz, Referenz der anderen Chiralität; (c) Klammerzerlegung
        aus z89 als Klebung: wörtlich auf der Diagonale, mit Referenztwistor, mit einer Partnernaht; Farbklebung über die
        U(2) der Transferströme.
Teil 3  Bilanz (E/C/L/O), Urteil zu 'Komposition = CSW', Folgen für K15, T3/T5 und die Szenarioachse Wechselwirkung.
"""
import json
import time
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.sparse as sps

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(111)
BETA = 4*np.pi
TOL = 1e-10


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def r(x, n=12):
    return float(np.round(float(x), n))


def laufzeit():
    say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")


# ------------------------------------------------------------------ Spinor-Helizität (Vergleichsmaschine)
def ang(a, b):
    return a[0]*b[1] - a[1]*b[0]


def cz(*shape):
    return rng.normal(size=shape) + 1j*rng.normal(size=shape)


def kinematik(n):
    """Zufällige komplexe Spinoren, lambda~_(n-2), lambda~_(n-1) aus der Impulserhaltung gelöst."""
    lam, lt = cz(n, 2), cz(n, 2)
    S = sum(np.outer(lam[i], lt[i]) for i in range(n - 2))
    sol = np.linalg.solve(np.array([lam[n - 2], lam[n - 1]]).T, -S)
    lt[n - 2], lt[n - 1] = sol[0], sol[1]
    return lam, lt


def impuls(lam, lt, idx):
    return sum(np.outer(lam[i], lt[i]) for i in idx)


def pt(lam, h):
    """Parke-Taylor (Vergleichsziel): <ab>^4/(<12>...<n1>), a, b die beiden Beine mit h = -1."""
    n = len(h)
    m = [i for i in range(n) if h[i] < 0]
    return ang(lam[m[0]], lam[m[1]])**4/np.prod([ang(lam[i], lam[(i + 1) % n]) for i in range(n)])


def ptq(lt, h):
    """Anti-MHV in eckigen Klammern: [ab]^4/([12]...[n1]), a, b die beiden Beine mit h = +1."""
    n = len(h)
    p = [i for i in range(n) if h[i] > 0]
    return ang(lt[p[0]], lt[p[1]])**4/np.prod([ang(lt[i], lt[(i + 1) % n]) for i in range(n)])


def a3(lam, lt, h):
    """Dreipunktamplituden: MHV in <>, Anti-MHV in []; die verschwindende Klammersorte gibt 0."""
    nm = sum(1 for x in h if x < 0)
    sa = max(abs(ang(lam[i], lam[(i + 1) % 3])) for i in range(3))
    ss = max(abs(ang(lt[i], lt[(i + 1) % 3])) for i in range(3))
    if nm == 2:
        return 0.0 if sa < 1e-9*ss else pt(lam, h)
    if nm == 1:
        return 0.0 if ss < 1e-9*sa else ptq(lt, h)
    return 0.0


def nullzerlegung(Ph):
    """Ph = lambda lambda~^T für eine Rang-1-Matrix."""
    k = int(np.argmax(np.abs(Ph).max(0)))
    zeile = int(np.argmax(np.abs(Ph[:, k])))
    lP = Ph[:, k].copy()
    return lP, Ph[zeile]/lP[zeile]


def amp(lam, lt, h, schieb=None):
    """Baumamplitude aus BCFW (Rekursion bis zu den Dreipunktamplituden), Standardverschiebung (h_i, h_j) = (-, +)."""
    n = len(h)
    if n == 3:
        return a3(lam, lt, h)
    nm = sum(1 for x in h if x < 0)
    if nm < 2 or n - nm < 2:
        return 0.0
    art = schieb or "-+"
    for i in range(n):
        j = (i + 1) % n
        if (h[i], h[j]) == {"-+": (-1, 1), "--": (-1, -1), "++": (1, 1), "+-": (1, -1)}[art]:
            return bcfw(lam, lt, h, i, j)
    return None


def bcfw(lam, lt, h, i, j):
    n = len(h)
    folge = [(j + k) % n for k in range(n)]                                           # j zuerst, i zuletzt
    summe = 0.0
    for m in range(2, n - 1):
        L, R = folge[:m], folge[m:]
        P = impuls(lam, lt, L)
        adj = np.array([[P[1, 1], -P[0, 1]], [-P[1, 0], P[0, 0]]])
        z = np.linalg.det(P)/(lt[j] @ adj @ lam[i])
        lam2, lt2 = lam.copy(), lt.copy()
        lt2[i] = lt[i] + z*lt[j]
        lam2[j] = lam[j] - z*lam[i]
        lP, ltP = nullzerlegung(impuls(lam2, lt2, L))
        for hp in (1, -1):
            aL = amp(np.array([lam2[k] for k in L] + [lP]), np.array([lt2[k] for k in L] + [-ltP]), [h[k] for k in L] + [hp])
            if aL == 0:
                continue
            aR = amp(np.array([lP] + [lam2[k] for k in R]), np.array([ltP] + [lt2[k] for k in R]), [-hp] + [h[k] for k in R])
            summe += aL/np.linalg.det(P)*aR
    return summe


def bogen(n):
    """Alle Zerlegungen der zyklischen Ordnung in zwei Bögen L | R mit je >= 2 Beinen (L ohne das letzte Bein)."""
    for s in range(n - 1):
        for l in range(2, n - 1):
            if s + l <= n - 1:
                L = list(range(s, s + l))
                yield L, list(range(s + l, n)) + list(range(0, s))


def csw(lam, lt, h, eta, prop=None, fort=None, nur=None, tausch=False):
    """MHV-Diagramme für NMHV (drei negative Helizitäten): Summe V_L(lambda_P) prop(P) V_R(lambda_P).
    eta: fester Referenzspinor (oder Funktion des Diagramms), prop: Propagator (Standard 1/det P),
    fort: Fortsetzung P -> lambda_P (Standard P eta), nur: Liste der verwendeten Diagramme, tausch: innere Helizität vertauscht."""
    n = len(h)
    summe, d = 0.0, 0
    for L, R in bogen(n):
        mL, mR = sum(h[k] < 0 for k in L), sum(h[k] < 0 for k in R)
        if (mL, mR) not in ((2, 1), (1, 2)):
            continue
        hL = 1 if mL == 2 else -1
        if tausch:
            hL = -hL
        d += 1
        if nur is not None and d not in nur:
            continue
        P = impuls(lam, lt, L)
        e = eta(L, R) if callable(eta) else eta
        lP = fort(P, L, R) if fort else P @ e
        VL = pt(np.array([lam[k] for k in L] + [lP]), [h[k] for k in L] + [hL]) if sum(
            x < 0 for x in [h[k] for k in L] + [hL]) == 2 else 0.0
        VR = pt(np.array([lP] + [lam[k] for k in R]), [-hL] + [h[k] for k in R]) if sum(
            x < 0 for x in [-hL] + [h[k] for k in R]) == 2 else 0.0
        summe += VL*(prop(P) if prop else 1/np.linalg.det(P))*VR
    return summe, d


def csw_quer(lam, lt, h, zeta):
    """Paritätsgespiegelte Regeln: Anti-MHV-Vertizes in [], lambda~_P = zeta^T P (Referenz der anderen Chiralität)."""
    n = len(h)
    summe, d = 0.0, 0
    for L, R in bogen(n):
        pL, pR = sum(h[k] > 0 for k in L), sum(h[k] > 0 for k in R)
        if (pL, pR) not in ((2, 1), (1, 2)):
            continue
        hL = -1 if pL == 2 else 1
        d += 1
        P = impuls(lam, lt, L)
        ltP = zeta @ P
        VL = ptq(np.array([lt[k] for k in L] + [ltP]), [h[k] for k in L] + [hL])
        VR = ptq(np.array([ltP] + [lt[k] for k in R]), [-hL] + [h[k] for k in R])
        summe += VL/np.linalg.det(P)*VR
    return summe, d


def rel(a, b):
    return float(abs(a - b)/abs(b))


HEL = {"5: 1-2-3-4+5+": [-1, -1, -1, 1, 1], "6: 1-2-3-4+5+6+ (split)": [-1, -1, -1, 1, 1, 1],
       "6: 1-2+3-4+5-6+": [-1, 1, -1, 1, -1, 1], "6: 1-2-3+4-5+6+": [-1, -1, 1, -1, 1, 1]}
NPUNKTE = 4                                                                            # Kinematikpunkte je Amplitude
KIN = {name: [kinematik(len(h)) for _ in range(NPUNKTE)] for name, h in HEL.items()}
REF = {name: [amp(lam, lt, h) for (lam, lt) in KIN[name]] for name, h in HEL.items()}

# =================================================================== Teil 1
say("Teil 1: Vergleichsmaschine (Spinor-Helizität, BCFW, CSW)")
say("  (a) Kinematik")
kin_def = {"Impulserhaltung": 0.0, "p_i^2": 0.0, "Schouten": 0.0, "P^2 = Summe <ij>[ij]": 0.0}
for n in (4, 5, 6):
    for _ in range(5):
        lam, lt = kinematik(n)
        kin_def["Impulserhaltung"] = max(kin_def["Impulserhaltung"], np.abs(impuls(lam, lt, range(n))).max())
        kin_def["p_i^2"] = max(kin_def["p_i^2"], max(abs(np.linalg.det(np.outer(lam[i], lt[i]))) for i in range(n)))
        i, j, k, l = 0, 1, 2, 3
        kin_def["Schouten"] = max(kin_def["Schouten"], abs(ang(lam[i], lam[j])*ang(lam[k], lam[l]) + ang(lam[i], lam[k])*ang(
            lam[l], lam[j]) + ang(lam[i], lam[l])*ang(lam[j], lam[k])))
        L = [0, 1, 2]
        kin_def["P^2 = Summe <ij>[ij]"] = max(kin_def["P^2 = Summe <ij>[ij]"], abs(np.linalg.det(impuls(lam, lt, L)) - sum(
            ang(lam[a], lam[b])*ang(lt[a], lt[b]) for a in L for b in L if a < b)))
ZAHLEN["kinematik_defekte"] = {k: float(v) for k, v in kin_def.items()}
check("[E] (a) Zufällige komplexe Kinematik (n = 4, 5, 6, je fünf Punkte, lambda~_(n-1), lambda~_n aus der Impulserhaltung "
      f"gelöst): Impulserhaltung {kin_def['Impulserhaltung']:.0e}, p_i^2 = det p_i {kin_def['p_i^2']:.0e}, Schouten "
      f"{kin_def['Schouten']:.0e}, und P^2 = det P_L = Summe_(i<j in L) <ij>[ij] ({kin_def['P^2 = Summe <ij>[ij]']:.0e}): Jede "
      "Kanalinvariante ist die Summe der Zwei-Naht-Invarianten (P(i) + P(j))^2 = <ij>[ij] aus z89",
      all(v < 1e-12 for v in kin_def.values()))

say("  (b) BCFW aus Dreipunktamplituden")
mhv_def, mhv_n = 0.0, 0
for n in (4, 5, 6):
    for a_ in range(n):
        for b_ in range(a_ + 1, n):
            h = [1]*n
            h[a_] = h[b_] = -1
            lam, lt = kinematik(n)
            mhv_def = max(mhv_def, rel(amp(lam, lt, h), pt(lam, h)))
            mhv_n += 1
lam5, lt5 = KIN["5: 1-2-3-4+5+"][0]
anti5 = max(rel(REF["5: 1-2-3-4+5+"][k], ptq(KIN["5: 1-2-3-4+5+"][k][1], HEL["5: 1-2-3-4+5+"])) for k in range(NPUNKTE))
schieb_def = {}
for name, h in HEL.items():
    d_ = 0.0
    for (lam, lt), A in zip(KIN[name], REF[name]):
        for art in ("--", "++"):
            v = amp(lam, lt, h, art)
            if v is not None:
                d_ = max(d_, rel(v, A))
    schieb_def[name] = d_


def split6(lam, lt):
    """Geschlossene Formel der Split-Helizität (Vergleichsziel): 1/<5|3+4|2] (<1|2+3|4]^3/([23][34]<56><61>s_234)
    + <3|4+5|6]^3/([61][12]<34><45>s_345)), Beine 1..6 = Indizes 0..5."""
    def sw(a_, Q, b_):
        return sum(ang(lam[a_], lam[k])*ang(lt[k], lt[b_]) for k in Q)
    s234, s345 = np.linalg.det(impuls(lam, lt, [1, 2, 3])), np.linalg.det(impuls(lam, lt, [2, 3, 4]))
    t1 = sw(0, [1, 2], 3)**3/(ang(lt[1], lt[2])*ang(lt[2], lt[3])*ang(lam[4], lam[5])*ang(lam[5], lam[0])*s234)
    t2 = sw(2, [3, 4], 5)**3/(ang(lt[5], lt[0])*ang(lt[0], lt[1])*ang(lam[2], lam[3])*ang(lam[3], lam[4])*s345)
    return (t1 + t2)/sw(4, [2, 3], 1)


nm6 = "6: 1-2-3-4+5+6+ (split)"
split_def = max(rel(split6(*KIN[nm6][k]), REF[nm6][k]) for k in range(NPUNKTE))
schlecht = [abs(amp(lam, lt, HEL[nm6], "+-")/A - 1) for (lam, lt), A in zip(KIN[nm6], REF[nm6])]
ZAHLEN["bcfw"] = {"MHV gegen Parke-Taylor (max rel., Anzahl)": [r(mhv_def, 15), mhv_n],
                  "A(1-2-3-4+5+) gegen [45]^4/([12][23][34][45][51])": r(anti5, 15),
                  "Verschiebungen [--> und [++> gegen [-+>": {k: r(v, 15) for k, v in schieb_def.items()},
                  "Split-Helizität gegen geschlossene Formel": r(split_def, 15),
                  "schlechte Verschiebung [+-> (|Verhältnis - 1|)": [r(x, 6) for x in schlecht],
                  "Beispielwerte A(1-2-3-4+5+), A_split": [str(np.round(REF['5: 1-2-3-4+5+'][0], 10)),
                                                         str(np.round(REF[nm6][0], 10))]}
check(f"[E] (b) BCFW aus den Dreipunktamplituden (Propagator 1/det P, keine weitere Konvention frei) gibt für alle {mhv_n} "
      f"MHV-Helizitätslagen n = 4, 5, 6 exakt Parke-Taylor (max rel. {mhv_def:.0e}) und für A(1-,2-,3-,4+,5+) exakt die "
      f"Anti-MHV-Formel [45]^4/([12][23][34][45][51]) (max rel. {anti5:.0e}); die geschlossene Split-Helizitätsformel bei sechs "
      f"Punkten stimmt auf {split_def:.0e}", mhv_def < TOL and anti5 < TOL and split_def < TOL)
check("[E] (b) Unabhängigkeit von der Verschiebung: [-,-> und [+,+> gegen [-,+> für alle vier NMHV-Amplituden "
      f"(max rel. {max(schieb_def.values()):.0e}). Negativkontrolle: Die Verschiebung [+,-> hat einen Randterm bei z -> "
      f"unendlich und verfehlt die Split-Amplitude um |Verhältnis - 1| = {[f'{x:.2f}' for x in schlecht]}: der Test trennt",
      max(schieb_def.values()) < TOL and min(schlecht) > 1e-3)
laufzeit()

say("  (c) CSW-Regeln (MHV-Diagramme)")
csw_def, eta_spread, ndiag = {}, {}, {}
for name, h in HEL.items():
    d_, werte = 0.0, []
    for (lam, lt), A in zip(KIN[name], REF[name]):
        for _ in range(5):
            c, ndiag[name] = csw(lam, lt, h, cz(2))
            d_ = max(d_, rel(c, A))
    csw_def[name] = d_
lam, lt = KIN[nm6][0]
eta0 = cz(2)
c0 = csw(lam, lt, HEL[nm6], eta0)[0]
skal = max(rel(csw(lam, lt, HEL[nm6], t*eta0)[0], c0) for t in (np.exp(0.7j), 3.1, -0.2 + 1.3j))
diagramme = [csw(lam, lt, HEL[nm6], eta0, nur=[d])[0] for d in range(1, ndiag[nm6] + 1)]
eta1 = cz(2)
diag_eta = [abs(csw(lam, lt, HEL[nm6], eta1, nur=[d])[0]/diagramme[d - 1] - 1) for d in range(1, ndiag[nm6] + 1)]
ohne = rel(c0 - diagramme[0], REF[nm6][0])
tausch = rel(csw(lam, lt, HEL[nm6], eta0, tausch=True)[0], REF[nm6][0])
vorz = rel(csw(lam, lt, HEL[nm6], eta0, prop=lambda P: -1/np.linalg.det(P))[0], REF[nm6][0])
ZAHLEN["csw"] = {"Diagramme je Amplitude": ndiag, "CSW gegen BCFW (max rel., 4 Punkte x 5 eta)": {k: r(v, 15) for k, v in csw_def.items()},
                 "eta -> t eta (t = e^0.7i, 3.1, -0.2+1.3i)": r(skal, 15),
                 "Einzeldiagramme hängen von eta ab (|D(eta1)/D(eta0) - 1|)": [r(x, 6) for x in diag_eta],
                 "Negativkontrollen (rel. Abweichung): ohne Diagramm 1, innere Helizität getauscht, Propagator -1/P^2":
                     [r(ohne, 6), r(tausch, 6), r(vorz, 6)]}
check(f"[E] (c) CSW = BCFW: A(1-,2-,3-,4+,5+) aus {ndiag['5: 1-2-3-4+5+']} MHV-Diagrammen, die Split-Helizität aus "
      f"{ndiag[nm6]}, 1-2+3-4+5-6+ aus {ndiag['6: 1-2+3-4+5-6+']} und 1-2-3+4-5+6+ aus {ndiag['6: 1-2-3+4-5+6+']} (max rel. "
      f"{max(csw_def.values()):.0e} über je vier Kinematikpunkte und fünf zufällige eta); eta -> t eta mit Phase oder Betrag "
      f"ändert nichts ({skal:.0e}), denn jede innere Linie trägt in eta den Grad -2h + 2h = 0",
      max(csw_def.values()) < TOL and skal < TOL)
check("[E] (c) Die eta-Unabhängigkeit ist eine Eigenschaft der Summe: Einzeldiagramme der Split-Amplitude ändern sich mit eta "
      f"um {[f'{x:.2f}' for x in diag_eta]}; Negativkontrollen verfehlen die Amplitude (ein Diagramm weg: {ohne:.2f}, "
      f"innere Helizität vertauscht: {tausch:.2f}, Propagator -1/P^2: {vorz:.2f})",
      min(diag_eta) > 1e-3 and min(ohne, tausch, vorz) > 1e-2)
befund("Teil 1: Die Vergleichsmaschine steht. BCFW aus den Dreipunktamplituden reproduziert Parke-Taylor, die Anti-MHV-Formel und "
       "die geschlossene Split-Helizitätsformel, unabhängig von der Verschiebung; die CSW-Regeln (MHV-Vertizes, 1/P^2, "
       "lambda_P = P eta) geben alle vier NMHV-Amplituden exakt und unabhängig von eta, während Einzeldiagramme von eta "
       "abhängen. Konvention: P^2 = det P, paritätssymmetrisch, ohne freie Vorzeichen", True)
laufzeit()

# =================================================================== Teil 2
say("\nTeil 2: Native Komposition zweier Nähte gegen die CSW-Zutaten")
say("  (a) Propagator zwischen zwei MHV-Vertizes")


def cdag_matrix(j, nb):
    """Jordan-Wigner-Erzeuger der Mode j auf nb Moden (wie z99)."""
    dim = 1 << nb
    s = np.arange(dim)
    src = s[((s >> j) & 1) == 0]
    tief = src & ((1 << j) - 1)
    par = np.zeros(len(src), dtype=np.int64)
    for b in range(j):
        par ^= (tief >> b) & 1
    return sps.csr_matrix(((1 - 2*par).astype(float), (src | (1 << j), src)), shape=(dim, dim))


XS = [Fr(1, 4) + k for k in (-2, -1, 0, 1)]                                          # Quarter a, Moden -7/4 .. 5/4
MODEN = [(i, x) for i in (0, 1) for x in XS]                                           # zwei Nähte in einer CAR (z106)
NB = len(MODEN)
CD = [cdag_matrix(j, NB).toarray() for j in range(NB)]
C_ = [m.T for m in CD]
fF = lambda x: 0.5*(1 - np.tanh(BETA*np.asarray(x, dtype=float)/2))                    # c(x) = 1/(1 + e^{4 pi x})


def einteilchen(t):
    h = np.zeros((NB, NB))
    for a_, (i, x) in enumerate(MODEN):
        h[a_, a_] = float(x)
        for b_, (i2, x2) in enumerate(MODEN):
            if i2 != i and x2 == x:
                h[a_, b_] = t                                                          # t (E^+_0 + E^-_0), z106-Hüpfterm
    return h


def kovarianz(t, beta):
    """<c^+_a c_b> im Gibbszustand von dGamma(h) (beta = inf: Grundzustand), im Fockraum gerechnet."""
    h = einteilchen(t)
    H = sum(h[a_, b_]*CD[a_] @ C_[b_] for a_ in range(NB) for b_ in range(NB) if h[a_, b_] != 0)
    w, V = np.linalg.eigh(H)
    if np.isinf(beta):
        rho = np.outer(V[:, 0], V[:, 0])
    else:
        p = np.exp(-beta*(w - w.min()))
        rho = (V*p) @ V.T/p.sum()
    return np.array([[np.trace(rho @ CD[a_] @ C_[b_]) for b_ in range(NB)] for a_ in range(NB)]), h, H


kreuz = {}
for t in (0.0, 0.15, 0.4):
    for beta in (BETA, np.inf):
        G, h, H = kovarianz(t, beta)
        ew, ev = np.linalg.eigh(h)
        f1 = (ev*(1/(1 + np.exp(beta*ew)) if not np.isinf(beta) else (ew < 0).astype(float))) @ ev.T
        X = G[:4, 4:]
        kreuz[(t, "beta = 4 pi" if not np.isinf(beta) else "Vakuum")] = {
            "max |Kreuzkern|": r(np.abs(X).max()), "max |abseits der Diagonale|": r(np.abs(X - np.diag(np.diag(X))).max(), 15),
            "Diagonale": [r(v, 10) for v in np.diag(X)], "Fock gegen f(h)": r(np.abs(G - f1.T).max(), 15)}
formel = [r(0.5*(fF(float(x) + 0.15) - fF(float(x) - 0.15)), 10) for x in XS]
Hn = sum(float(x)*CD[a_] @ C_[a_] for a_, (i, x) in enumerate(MODEN))
Hh = einteilchen(1.0) - einteilchen(0.0)
Hhop = sum(Hh[a_, b_]*CD[a_] @ C_[b_] for a_ in range(NB) for b_ in range(NB) if Hh[a_, b_] != 0)
komm = float(np.abs(Hn @ Hhop - Hhop @ Hn).max())
ZAHLEN["kreuzkern_gemeinsame_car"] = {f"t = {k[0]}, {k[1]}": v for k, v in kreuz.items()}
ZAHLEN["kreuzkern_formel_t015_beta4pi"] = formel
ZAHLEN["[E^nat, Hüpfterm]"] = komm
k0 = [kreuz[(0.0, b)]["max |Kreuzkern|"] for b in ("beta = 4 pi", "Vakuum")]
check("[E negativ] (a) Kreuzkern der gemeinsamen CAR (zwei Nähte, je vier Moden des Quarters a, Fockraum 256): Im nativen "
      f"Zustand ist <psi*_(0,x) psi_(1,y)> = 0 für alle x, y, im KMS-Zustand bei beta = 4 pi und im Vakuum ({k0}). Die "
      "Hardy-Komposition hat eine blockdiagonale Fermisee-Projektion P_0 (+) P_1: Zwischen zwei Nähten gibt es keinen "
      "Szegő-Kern, die Komposition ist frei (Produkt), und eine CSW-Klebung mit diesem Kern gibt A = 0 statt A(1-,2-,3-,4+,5+)",
      max(k0) == 0.0 and all(v["Fock gegen f(h)"] < 1e-12 for v in kreuz.values()))
d15 = kreuz[(0.15, "beta = 4 pi")]
check("[E negativ] (a) Hüpfterm t (E^+_0 + E^-_0) aus der U(2) der Transferströme (z106): Der Kreuzkern ist exakt "
      f"diagonal in den Moden (abseits {d15['max |abseits der Diagonale|']:.0e}) und gleich (c(x + t) - c(x - t))/2, bei t = 0.15 "
      f"und beta = 4 pi {d15['Diagonale']} (Formel {formel}); im Vakuum ist er für |t| < 1/4 exakt 0 "
      f"({kreuz[(0.15, 'Vakuum')]['max |Kreuzkern|']}: die Quarterlücke schützt den Fermisee) und erst für t = 0.4 ungleich "
      f"null ({kreuz[(0.4, 'Vakuum')]['Diagonale']}). Er vertauscht mit E^nat ([., .] = {komm:.0e}): Er setzt Mode x der einen "
      "Naht mit derselben Mode der anderen gleich, ohne Stufe und ohne Lage. Ein Kern, der von der Relativlage der Nähte nicht "
      "abhängt, hat seinen Impulsträger bei P = 0 und trägt zu keinem Kanal mit P^2 != 0 bei",
      d15["max |abseits der Diagonale|"] < 1e-13 and np.allclose(d15["Diagonale"], formel, atol=1e-10)
      and kreuz[(0.15, "Vakuum")]["max |Kreuzkern|"] < 1e-13 and kreuz[(0.4, "Vakuum")]["max |Kreuzkern|"] > 0.4 and komm < 1e-12)
laufzeit()

# --- Szegő-Kern einer Naht an den Inzidenzpunkten der anderen (Penrose-Transformation)
EPS2 = np.array([[0, 1], [-1, 0]])


def inzidenz(x, y, a, b):
    """Ebenen A = (-y^T a, a), B = (-y^T b, b) enthalten L_y; auf L_x ist A.Z = a^T (x - y) lambda. Schnittpunkte s = l1/l0."""
    al, be = (x - y).T @ a, (x - y).T @ b
    lA, lB = EPS2 @ al, EPS2 @ be                                                      # al . lA = 0, be . lB = 0
    return al, be, lA, lB


def penrose(x, y, a, b, M=64):
    """(1/2 pi i) Kontur um die Nullstelle von A.Z auf L_x, lambda = (1, s): Trapezregel auf einem Kreis."""
    al, be, lA, lB = inzidenz(x, y, a, b)
    sA, sB = -al[0]/al[1], -be[0]/be[1]
    rad = 0.4*abs(sA - sB)
    th = 2*np.pi*np.arange(M)/M
    s = sA + rad*np.exp(1j*th)
    f = 1/((al[0] + al[1]*s)*(be[0] + be[1]*s))
    return np.mean(f*rad*np.exp(1j*th)), lA, lB, sA, sB


def kern_quarter(s1, s2, alpha, beta=np.inf, K=200):
    """Nahtkern sum_x (1 - c(x)) s1^(-x-1/2) s2^(x-1/2), x in Z + alpha (Vakuum: x > 0)."""
    xs = np.arange(-K, K) + alpha
    w = (xs > 0).astype(float) if np.isinf(beta) else 1 - fF(xs)
    return np.sum(w*np.exp((-xs - 0.5)*np.log(s1 + 0j) + (xs - 0.5)*np.log(s2 + 0j)))


y0, a0, b0 = cz(2, 2), cz(2), cz(2)
pen, cau, qa, qb, qab, therm = [], [], [], [], [], []
zweig = []
for _ in range(16):
    x = cz(2, 2)
    phi, lA, lB, sA, sB = penrose(x, y0, a0, b0)
    dxy = np.linalg.det(x - y0)
    pen.append(phi*dxy*ang(a0, b0))
    cau.append(dxy/ang(lA, lB))
    if 0.2 < abs(sB/sA) < 0.75 or 0.2 < abs(sA/sB) < 0.75:
        s1, s2 = (sA, sB) if abs(sB) < abs(sA) else (sB, sA)
        c_ = 1/(s1 - s2)
        qa.append(kern_quarter(s1, s2, 0.25)/c_)
        zweig.append(abs(qa[-1] - np.exp(0.25*(np.log(s1 + 0j) - np.log(s2 + 0j)))))
        qb.append(kern_quarter(s1, s2, 0.75)/c_)
        qab.append(qa[-1]*qb[-1])
        therm.append(kern_quarter(s1, s2, 0.25, BETA)/kern_quarter(s1, s2, 0.25) - 1)
spread = lambda v: float(np.abs(np.array(v) - v[0]).max()/abs(v[0])) if len(v) > 1 else float("nan")
ZAHLEN["inzidenz_szegoe"] = {"Penrose * det(x - y) <ab> (soll -1)": [str(np.round(v, 12)) for v in pen[:3]],
                             "det(x - y)/<l_A l_B> (soll 1/<ab>)": r(spread(cau), 15),
                             "Quarter a: Kern/Cauchy (Streuung)": r(spread(qa), 6), "Quarter b": r(spread(qb), 6),
                             "Produkt a*b gegen (s1/s2)^0": r(max(abs(v - 1) for v in qab), 12),
                             "Quarter a gegen (s1/s2)^(1/4)": r(max(zweig), 12),
                             "thermisch/Vakuum - 1 (Quarter a, beta = 4 pi)": [r(abs(v), 6) for v in therm],
                             "c(1/4) = 1/(1 + e^pi)": r(fF(0.25), 8), "Punkte mit konvergenter Modenreihe": len(qa)}
check("[E] (a) Szegő-Kern über die Inzidenz: Für zwei Ebenen A, B durch L_y ist die Penrose-Transformation von "
      "1/((A.Z)(B.Z)) auf L_x (Konturintegral, 64 Punkte) gleich -1/(<ab> det(x - y)) (16 zufällige komplexe x, Abweichung "
      f"{max(abs(v + 1) for v in pen):.0e}), und das ist der homogene Cauchykern 1/<l_A l_B> der Naht L_x zwischen ihren "
      f"Schnittpunkten mit A und B (det(x - y)/<l_A l_B> konstant, Streuung {spread(cau):.0e}). Der 4D-Propagator zwischen den "
      "Lagen zweier Nähte ist also der Szegő-Kern der einen an den Inzidenzpunkten der anderen; die Normierung 1/<ab> hängt "
      "an der Wahl der Ebenen im Büschel durch L_y", max(abs(v + 1) for v in pen) < 1e-10 and spread(cau) < 1e-12)
check("[E negativ] (a) Der native Nahtkern ist nicht der Cauchykern: Die Moden x in Z + 1/4 geben im Vakuum Cauchy mal "
      f"(s1/s2)^(1/4) (Abweichung {max(zweig):.0e}), dessen Verhältnis zu 1/det(x - y) über die Lagen streut ({spread(qa):.2f}; Quarter b mit "
      f"(s1/s2)^(-1/4): {spread(qb):.2f}, {len(qa)} Lagen mit konvergenter Modenreihe). Der Twist hängt an den Marken 0 und "
      f"unendlich der Naht; nur das Produkt der Quarter a und b ist twistfrei (|a*b/Cauchy^2 - 1| = "
      f"{max(abs(v - 1) for v in qab):.0e}), und das ist ein Bilinear (Gewicht 2), kein Propagator. Bei beta = 4 pi weicht der "
      f"Quarterkern zusätzlich um {[f'{abs(v):.3f}' for v in therm[:3]]} vom Vakuum ab (Leitterm c(1/4) = 1/(1 + e^pi) = "
      f"{fF(0.25):.4f}): Ein Propagator als Funktion von (x - y)^2 entsteht nur im Vakuumgrenzfall ohne Quartertwist",
      len(qa) >= 3 and max(zweig) < 1e-10 and spread(qa) > 1e-2 and spread(qb) > 1e-2 and max(abs(v - 1) for v in qab) < 1e-9
      and min(abs(v) for v in therm) > 1e-3)

mp.mp.dps = 20
hankel = {}
for Pm in (0.5, 1.0, 2.0):
    integral = mp.quadosc(lambda u: mp.besselj(1, u), [0, mp.inf], zeros=lambda k: mp.besseljzero(1, k))
    hankel[Pm] = float(4*mp.pi**2/Pm*integral/Pm*Pm**2)                                # FT[1/|x|^2](P) * P^2
ZAHLEN["euklidische_FT_1_durch_x2_mal_P2"] = {str(k): r(v, 12) for k, v in hankel.items()}
check("[E] (a) Euklidische Fouriertransformation (quaternionische Realform, z61/z84): int d^4x e^(-iPx)/|x|^2 = (4 pi^2/|P|) "
      f"int_0^inf J_1(|P| r) dr = 4 pi^2/|P|^2; P^2 mal Transformierte = {[f'{v:.10f}' for v in hankel.values()]} für |P| = "
      f"0.5, 1, 2 (4 pi^2 = {4*np.pi**2:.10f}). Der Kern 1/(x - y)^2 der Inzidenzrechnung ist im Impulsraum genau 1/P^2",
      all(abs(v - 4*np.pi**2) < 1e-8 for v in hankel.values()))
laufzeit()

# --- Residuen an schneidenden Nähten (z61: L_x und L_y schneiden sich genau für (x - y)^2 = 0)
lam, lt = KIN[nm6][1]
hh = HEL[nm6]
i_, j_ = 1, 4                                                                          # lambda~_2 verschoben, lambda_5 verschoben
L234 = [1, 2, 3]


def verschoben(w):
    l2, t2 = lam.copy(), lt.copy()
    t2[i_] = lt[i_] + w*lt[j_]
    l2[j_] = lam[j_] - w*lam[i_]
    return l2, t2


P0 = impuls(lam, lt, L234)
adj = np.array([[P0[1, 1], -P0[0, 1]], [-P0[1, 0], P0[0, 0]]])
wstern = -np.linalg.det(P0)/(lt[j_] @ adj @ lam[i_])                                     # det P_234(w) linear in w
l2, t2 = verschoben(wstern)
lP, ltP = nullzerlegung(impuls(l2, t2, L234))
AL = pt(np.array([l2[k] for k in L234] + [lP]), [hh[k] for k in L234] + [1])
AR = pt(np.array([lP] + [l2[k] for k in (4, 5, 0)]), [-1] + [hh[k] for k in (4, 5, 0)])
res = {}
for eps in (1e-4, 1e-6):
    l3, t3 = verschoben(wstern + eps)
    werte = [csw(l3, t3, hh, cz(2))[0]*np.linalg.det(impuls(l3, t3, L234)) for _ in range(3)]
    res[eps] = max(rel(v, AL*AR) for v in werte)
ZAHLEN["residuum_P234"] = {f"eps = {k}": r(v, 10) for k, v in res.items()}
check("[E/L] (a) Schneidende Nähte = Residuen: Nahe dem Pol P_234^2 -> 0 (Verschiebung lambda~_2 + w lambda~_5, "
      f"lambda_5 - w lambda_2) gibt P^2 A_CSW mit drei zufälligen eta das Produkt A(2-,3-,4+,-P+) A(P-,5+,6+,1-) der "
      f"On-shell-MHV-Amplituden (rel. {res[1e-4]:.0e} bei eps = 1e-4, {res[1e-6]:.0e} bei 1e-6, linear in eps). Am Pol ist "
      "der innere Impuls null, also schneiden sich die Geraden der beiden Vertexlagen (z61): Die native Inzidenzklebung liefert "
      "genau die Residuen, die Off-shell-Fortsetzung zwischen den Polen nicht", res[1e-6] < 1e-4 and res[1e-6] < res[1e-4]/10)
befund("Teil 2 (a): Die gemeinsame CAR hat keinen Kern zwischen Nähten (blockdiagonal; der Hüpfterm ist modendiagonal, "
       "lagenblind und im Vakuum durch die Quarterlücke unterdrückt). Die Form 1/P^2 ist dagegen native Kinematik: Sie ist der "
       "Szegő-Kern der einen Naht an den Inzidenzpunkten der anderen (Penrose-Transformation, euklidische Fouriertransformation), "
       "aber nur im Vakuumgrenzfall ohne Quartertwist; die nativen Quarterkerne tragen den Twist (s1/s2)^(+-1/4) und bei "
       "beta = 4 pi die Besetzung c(1/4). Schneidende Nähte geben genau die Residuen. Was fehlt, ist der Träger: ein Feld, das "
       "zwischen den Lagen zweier Nähte propagiert, mit einer festen Normierung", True)
laufzeit()

say("  (b) Off-shell-Fortsetzung lambda_P = P|eta]: Wer liefert eta?")
KAPPA = np.sqrt(np.tanh(np.pi))                                                        # kappa^2 = tanh(pi m), m = 1 (z85)
U_PI = np.array([[1, 0, -1, 0], [0, 1, 0, -1]], dtype=complex)                       # z89: span{u} = {pi = 0}
V_OM = np.array([[1, 0, 1, 0], [0, 1, 0, 1]], dtype=complex)                         # z89: span{v} = {omega = 0}
pi_von = lambda Z: np.array([Z[0] + Z[2], Z[1] + Z[3]])/(np.sqrt(2)*KAPPA)             # z85: d +- e^+ = sqrt2 (kappa pi, omega/kappa)
om_von = lambda Z: KAPPA*np.array([Z[0] - Z[2], Z[1] - Z[3]])/np.sqrt(2)


def zurueck(pi_, om_):
    p, o = np.sqrt(2)*KAPPA*pi_, np.sqrt(2)*om_/KAPPA
    return np.array([(p[0] + o[0])/2, (p[1] + o[1])/2, (p[0] - o[0])/2, (p[1] - o[1])/2])


def keil(a_, b_):
    return np.outer(a_, b_) - np.outer(b_, a_)


I_PI, I_OM = keil(*U_PI), keil(*V_OM)
BUESCHEL = {"I_pi (flach)": I_PI, "I_omega (Ursprung)": I_OM, "I_AdS = I_pi - I_omega": I_PI - I_OM, "I_dS = I_pi + I_omega": I_PI + I_OM}
pf = lambda B: B[0, 1]*B[2, 3] - B[0, 2]*B[1, 3] + B[0, 3]*B[1, 2]
TRIPEL = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]


def kern_dim(B):
    """Z mit Z ^ B = 0 (Punkte auf der Geraden des Bivektors): (Z ^ B)_abc = Z_a B_bc - Z_b B_ac + Z_c B_ab."""
    M = np.zeros((4, 4), dtype=complex)
    for row, (a_, b_, c_) in enumerate(TRIPEL):
        M[row, a_] += B[b_, c_]
        M[row, b_] -= B[a_, c_]
        M[row, c_] += B[a_, b_]
    return 4 - np.linalg.matrix_rank(M, tol=1e-10)


ort = {k: {"Pfaffsche": str(np.round(pf(B), 12)), "dim{Z: Z ^ I = 0}": int(kern_dim(B))} for k, B in BUESCHEL.items()}
koord = max(np.abs(pi_von(u)).max() for u in U_PI) + max(np.abs(om_von(v)).max() for v in V_OM)
roundtrip = max(np.abs(zurueck(pi_von(Z), om_von(Z)) - Z).max() for Z in cz(5, 4))
verschiebung = {"I_pi": 0.0, "I_omega": 0.0}
for _ in range(20):
    a_ = cz(2, 2)
    for k, Ub in (("I_pi", U_PI), ("I_omega", V_OM)):
        Z = cz(2) @ Ub
        Zt = zurueck(pi_von(Z), om_von(Z) + a_ @ pi_von(Z))                              # Translation: omega -> omega + a pi
        verschiebung[k] = max(verschiebung[k], float(np.abs(Zt - Z).max()/np.abs(Z).max()))
ZAHLEN["referenzort"] = {"Büschel": ort, "Koordinaten pi(u) = 0, omega(v) = 0": r(koord, 15), "Rücktransformation": r(roundtrip, 15),
                         "max relative Verschiebung unter Translationen": {k: r(v, 12) for k, v in verschiebung.items()}}
check("[E] (b) Referenzort im Büschel aus z89 (Koordinaten von z85, pi(u_j) = 0, omega(v_j) = 0, Fehler "
      f"{koord:.0e}): Punkte Z mit Z ^ I = 0 gibt es genau auf den einfachen Bivektoren: I_pi {ort['I_pi (flach)']}, I_omega "
      f"{ort['I_omega (Ursprung)']}, I_AdS {ort['I_AdS = I_pi - I_omega']}, I_dS {ort['I_dS = I_pi + I_omega']}. Die native "
      "E-invariante Wahl I_AdS (z89/z98) trägt keinen Referenztwistor. Unter den Translationen omega -> omega + a pi bleibt "
      f"jeder Punkt von I_pi fest (Verschiebung {verschiebung['I_pi']:.0e}), Punkte von I_omega wandern "
      f"({verschiebung['I_omega']:.2f}): Nur I_pi gibt eine translationsinvariante Eichung, also eine Impulsraumregel",
      kern_dim(I_PI) == 2 and kern_dim(I_OM) == 2 and kern_dim(I_PI - I_OM) == 0 and kern_dim(I_PI + I_OM) == 0
      and abs(pf(I_PI)) < 1e-14 and abs(pf(I_PI - I_OM)) > 1 and koord < 1e-14 and roundtrip < 1e-13
      and verschiebung["I_pi"] < 1e-13 and verschiebung["I_omega"] > 0.1)

eta_I = {}
for name in ("5: 1-2-3-4+5+", nm6):
    d_ = 0.0
    for (lam, lt), A in zip(KIN[name], REF[name]):
        for _ in range(4):
            W = cz(4)
            Zst = I_PI @ W                                                             # Punkt auf I_pi: Z* = I_pi . W
            for th in (0.0, 0.9, 2.3):
                et = om_von(np.exp(2j*th)*Zst)                                          # Modulphase: I -> e^(2 i theta) I
                d_ = max(d_, rel(csw(lam, lt, HEL[name], EPS2 @ et)[0], A))
    eta_I[name] = d_
auf_I = max(np.abs(pi_von(I_PI @ cz(4))).max() for _ in range(10))
ZAHLEN["eta_aus_I_pi"] = {"pi(Z*) für Z* = I_pi W": r(auf_I, 15), "CSW gegen BCFW": {k: r(v, 15) for k, v in eta_I.items()}}
check("[E/C] (b) eta aus I_pi: Z* = I_pi W liegt für jeden dualen Twistor W auf der Geraden im Unendlichen (pi(Z*) = "
      f"{auf_I:.0e}); eta~ = omega(Z*) mit der Modulphase e^(2 i theta) (theta = 0, 0.9, 2.3, z98: Ladung 2) gibt in den "
      f"CSW-Regeln A(1-,2-,3-,4+,5+) und die Split-Amplitude exakt (max rel. {max(eta_I.values()):.0e}, je vier "
      "Kinematikpunkte und vier W). Gebraucht wird nur die Gerade I_pi; welcher Punkt und welche Phase, ist Eichung",
      max(eta_I.values()) < TOL and auf_I < 1e-14)


def rahmen(P, L, R):
    """Native Zeit T = 1 (Herm(2), z61) als Rahmen: lambda_P = dominanter Eigenvektor von P, je Propagator."""
    w, V = np.linalg.eig(P)
    k = int(np.argmax(np.abs(w)))
    return V[:, k]*np.sqrt(w[k])


rahmen_def = {name: max(rel(csw(lam, lt, HEL[name], None, fort=rahmen)[0], A) for (lam, lt), A in zip(KIN[name], REF[name]))
              for name in ("5: 1-2-3-4+5+", nm6)}
falsch_chiral = {name: max(rel(csw(lam, lt, HEL[name], None, fort=lambda P, L, R, z=cz(2): z @ P)[0], A)
                           for (lam, lt), A in zip(KIN[name], REF[name])) for name in ("5: 1-2-3-4+5+", nm6)}
quer = {}
for name in (nm6, "6: 1-2+3-4+5-6+", "6: 1-2-3+4-5+6+"):
    d_ = 0.0
    for (lam, lt), A in zip(KIN[name], REF[name]):
        for _ in range(3):
            v, nq = csw_quer(lam, lt, HEL[name], cz(2))
            d_ = max(d_, rel(v, A))
    quer[name] = (d_, nq)
ZAHLEN["fortsetzung_kandidaten"] = {"Rahmen T = 1 je Propagator (rel.)": {k: r(v, 6) for k, v in rahmen_def.items()},
                                    "MHV-Vertizes mit zeta^T P (falsche Chiralität)": {k: r(v, 6) for k, v in falsch_chiral.items()},
                                    "Anti-MHV-Vertizes mit lambda~_P = zeta^T P": {k: [r(v[0], 15), v[1]] for k, v in quer.items()}}
check("[E negativ] (b) Rahmenabhängige native Referenz: lambda_P = dominanter Eigenvektor von P relativ zur nativen Zeit T = 1 "
      "(gleichbedeutend mit einem eta~_P je Propagator) verfehlt die Amplituden um "
      f"{[f'{v:.2f}' for v in rahmen_def.values()]}; ebenso MHV-Vertizes mit einem Spinor der anderen Chiralität, "
      f"zeta^T P ({[f'{v:.2f}' for v in falsch_chiral.values()]}). Die Referenz muss global und gestrichen sein, ein Punkt auf I_pi",
      min(rahmen_def.values()) > 1e-2 and min(falsch_chiral.values()) > 1e-2)
check("[E] (b) Gespiegelte Regeln: Anti-MHV-Vertizes mit lambda~_P = zeta^T P (zeta beliebig, dualer Punkt der Geraden im "
      "Unendlichen) geben die drei Sechspunkt-NMHV-Amplituden exakt "
      f"({ {k: f'{v[0]:.0e} aus {v[1]} Diagrammen' for k, v in quer.items()} }); A(1-,2-,3-,4+,5+) ist dort ein einziger "
      "Anti-MHV-Vertex (Teil 1b). I_pi trägt beide chiralen Entwicklungen, über seine Gerade und ihre duale",
      max(v[0] for v in quer.values()) < TOL)
laufzeit()

say("  (c) Klammerzerlegung B = m_Z o (eps c^+ (x) c^+) aus z89 als Klebung")
diag_eps = 0.0
for name in ("5: 1-2-3-4+5+", nm6):
    for (lam, lt) in KIN[name]:
        n = len(HEL[name])
        for m in range(2, n - 1):
            L = list(range(m))
            P = impuls(lam, lt, L)
            Ph = P - (np.linalg.det(P)/(lt[n - 1] @ np.array([[P[1, 1], -P[0, 1]], [-P[1, 0], P[0, 0]]]) @ lam[0]))*np.outer(lam[0], lt[n - 1])
            lP, ltP = nullzerlegung(Ph)
            diag_eps = max(diag_eps, abs(ang(lP, lP)), abs(ang(ltP, -ltP)), abs(np.linalg.det(Ph)))
ZAHLEN["klammer_diagonal"] = r(diag_eps, 15)
check("[E negativ] (c) Wörtliche Klebung zweier Vertizes durch eps auf den inneren Beinen: Auf der Diagonale (dasselbe "
      "On-shell-Bein P auf beiden Seiten, BCFW-Kanäle) sind <lambda_P lambda_-P> und [lambda~_P lambda~_-P] exakt null "
      f"(max {diag_eps:.0e}): die z89-Aussage 'Identifizieren der Impulse ist der kollineare Rand'. Die so geklebte Amplitude "
      "ist 0", diag_eps < 1e-12)

eps_Z, eps_P = 0.0, {}
for name in ("5: 1-2-3-4+5+", nm6):
    for (lam, lt), A in zip(KIN[name], REF[name]):
        W = cz(4)
        et = om_von(I_PI @ W)

        def fort_B(P, L, R, lam=lam, lt=lt, et=et):
            return sum(lam[k]*ang(lt[k], et) for k in L)                               # m_Z(eps(c^+(k), Z*) c(k))

        Lt = [0, 1, 2]
        eps_Z = max(eps_Z, float(np.abs(fort_B(impuls(lam, lt, Lt), Lt, None) - impuls(lam, lt, Lt) @ (EPS2 @ et)).max()),
                    rel(csw(lam, lt, HEL[name], None, fort=fort_B)[0], A))
    eps_P[name] = max(rel(csw(lam, lt, HEL[name], lambda L, R, lt=lt: EPS2 @ lt[R[0]])[0], A)
                      for (lam, lt), A in zip(KIN[name], REF[name]))
ZAHLEN["klammer_als_fortsetzung"] = {"mit Z* auf I_pi (Identität und Amplitude)": r(eps_Z, 15),
                                     "mit der ersten Naht jenseits des Schnitts (rel.)": {k: r(v, 6) for k, v in eps_P.items()}}
check("[E/C] (c) Dieselbe Struktur mit dem Referenztwistor statt der Partnernaht: lambda_P = sum_(k in L) c(k) eps(c^+(k), Z*), "
      "also die eps-Kontraktion jeder Naht der einen Seite mit Z* auf I_pi und danach die holomorphe Fusion, ist exakt P|eta] "
      f"und gibt beide Amplituden (max {eps_Z:.0e}). Mit der ersten Naht jenseits des Schnitts als zweitem Argument (je Diagramm "
      f"eine andere Referenz) fällt die Summe um {[f'{v:.2f}' for v in eps_P.values()]} daneben. Die z89-Zerlegung ist die "
      "CSW-Fortsetzung, wenn ein Faktor der Referenztwistor ist; als Klebung zweier Nähte ist sie null oder falsch",
      eps_Z < TOL and min(eps_P.values()) > 1e-2)

# --- Farbklebung über die U(2) der Transferströme (z99/z106), aus dem Fockmodell abgelesen
vak = np.zeros(1 << NB)
vak[int(sum(1 << j for j, (i, x) in enumerate(MODEN) if x < 0))] = 1.0
ein = [CD[MODEN.index((i, Fr(1, 4)))] @ vak for i in (0, 1)]                            # psi*_(i, 1/4) Omega
Ep0 = sum(CD[MODEN.index((0, x))] @ C_[MODEN.index((1, x))] for x in XS)
Q0 = sum(CD[MODEN.index((i, x))] @ C_[MODEN.index((i, x))] for i in (0, 1) for x in XS)
ops = {"E+": Ep0, "E-": Ep0.T, "H": sum((CD[MODEN.index((0, x))] @ C_[MODEN.index((0, x))]
                                        - CD[MODEN.index((1, x))] @ C_[MODEN.index((1, x))]) for x in XS), "Q": Q0}
mat = {k: np.array([[ein[a_] @ (X @ ein[b_]) for b_ in (0, 1)] for a_ in (0, 1)]) for k, X in ops.items()}
mat["Q"] = mat["Q"] - (vak @ (Q0 @ vak))*np.eye(2)
Tb = [(mat["E+"] + mat["E-"])/np.sqrt(2), 1j*(mat["E+"] - mat["E-"])/np.sqrt(2), mat["H"]/np.sqrt(2), mat["Q"]/np.sqrt(2)]
gram = np.array([[np.trace(A_ @ B_) for B_ in Tb] for A_ in Tb])
voll = np.einsum("bij,bkl->ijkl", np.array(Tb), np.array(Tb))
soll = np.einsum("il,jk->ijkl", np.eye(2), np.eye(2))
Xr, Yr = cz(2, 2), cz(2, 2)
kleb_u2 = abs(sum(np.trace(Xr @ T) * np.trace(T @ Yr) for T in Tb) - np.trace(Xr @ Yr))
kleb_su2 = sum(np.trace(Xr @ T) * np.trace(T @ Yr) for T in Tb[:3]) - np.trace(Xr @ Yr)
ZAHLEN["farbklebung"] = {"Gram Tr(T^a T^b) - delta": r(np.abs(gram - np.eye(4)).max(), 15),
                         "Vollständigkeit u(2)": r(np.abs(voll - soll).max(), 15), "Spurklebung u(2)": r(kleb_u2, 15),
                         "su(2) allein: Abweichung + Tr X Tr Y/2": r(abs(kleb_su2 + np.trace(Xr)*np.trace(Yr)/2), 15)}
check("[E/L] (c) Farbhälfte der Klebung: Die Nullmoden E^+_0, E^-_0, H_0 und die Gesamtladung Q auf den Einteilchenzuständen "
      "psi*_(i,1/4) Omega (Fockmodell oben) geben eine Orthonormalbasis T^b von u(2) (Gram-Defekt "
      f"{np.abs(gram - np.eye(4)).max():.0e}) mit sum_b T^b_ij T^b_kl = delta_il delta_jk ({np.abs(voll - soll).max():.0e}); damit "
      f"klebt sum_b Tr(X T^b) Tr(T^b Y) = Tr(XY) zwei Spuren zu einer ({kleb_u2:.0e}), wie es die Farbordnung der MHV-Diagramme "
      "verlangt. Ohne die Gesamtladung (nur su(2)) fehlt -Tr X Tr Y/2 (Defekt "
      f"{abs(kleb_su2 + np.trace(Xr)*np.trace(Yr)/2):.0e}). Lesart: Farbe = Nahtmarke",
      np.abs(gram - np.eye(4)).max() < 1e-12 and np.abs(voll - soll).max() < 1e-12 and kleb_u2 < 1e-12
      and abs(kleb_su2 + np.trace(Xr)*np.trace(Yr)/2) < 1e-12 and abs(kleb_su2) > 1e-3)
befund("Teil 2 (b, c): Die Referenz der Off-shell-Fortsetzung ist nativ vorhanden: Jeder Punkt der flachen Geraden I_pi (z89) "
       "gibt eta, Phase und Punkt sind Eichung, die Gegenwahl I_omega ist nicht translationsinvariant und I_AdS hat keinen Punkt. "
       "Die Klammerzerlegung aus z89 ist genau die CSW-Fortsetzung, wenn ihr zweiter Faktor der Referenztwistor ist; als Klebung "
       "zweier Nähte auf der Diagonale verschwindet sie. Die Farbklebung ist die U(2) der Transferströme. Nicht geliefert wird "
       "der Propagator selbst: Zwischen zwei Nähten gibt es nativ keinen Austausch", True)
laufzeit()

# =================================================================== Teil 3
say("\nTeil 3: Bilanz für T2 'Komposition der Nähte = MHV-Diagramme'")
STATUS = {
    "Vergleichsmaschine": ("E", "BCFW aus Dreipunktamplituden = Parke-Taylor (31 MHV-Lagen, n <= 6), Anti-MHV- und Split-Formel, "
                                "verschiebungsunabhängig; CSW = BCFW für vier NMHV-Amplituden (4, 6, 9, 8 Diagramme), eta-unabhängig, "
                                "Einzeldiagramme eta-abhängig; Negativkontrollen trennen"),
    "Kreuzkern der gemeinsamen CAR": ("E (negativ)", "blockdiagonale Fermisee-Projektion: <psi*_0 psi_1> = 0 im Vakuum und bei "
                                                     "beta = 4 pi; die Komposition ist ein Produkt, CSW-Klebung damit A = 0"),
    "Hüpfterm der Transferströme": ("E (negativ)", "Kreuzkern (c(x + t) - c(x - t))/2, modendiagonal, vertauscht mit E^nat, "
                                                   "lagenblind; im Vakuum für |t| < 1/4 exakt 0 (Quarterlücke)"),
    "Propagatorform": ("E/C", "Szegő-(Cauchy-)Kern von L_x an den Schnittpunkten mit zwei Ebenen durch L_y = -1/(<ab> det(x - y)) "
                              "= Penrose-Transformation von 1/((A.Z)(B.Z)); euklidisch 4 pi^2/P^2; Normierung (Büschel, 4 pi^2) "
                              "nicht nativ, durch Faktorisierung festgelegt"),
    "native Nahtkerne": ("E (negativ)", "Quarter a, b: Cauchy mal (s1/s2)^(+-1/4) (Twist an den Marken), bei beta = 4 pi zusätzlich "
                                        "c(1/4) = 0.041; nur Vakuum ohne Twist gibt eine Funktion von (x - y)^2"),
    "Residuen": ("E/L", "P^2 A -> A_L A_R am Pol, für jedes eta; schneidende Nähte (z61) = On-shell-Klebung"),
    "Referenz eta": ("E/C", "Z* = I_pi W auf der flachen Geraden aus z89, eta~ = omega(Z*); Punkt und Modulphase sind Eichung; "
                            "I_omega nicht translationsinvariant, I_AdS und I_dS ohne Punkt (Pfaffsche -+4)"),
    "falsche Referenzen": ("E (negativ)", "rahmenabhängig (native Zeit T = 1 je Propagator) und falsche Chiralität scheitern"),
    "gespiegelte Regeln": ("E", "Anti-MHV-Vertizes mit lambda~_P = zeta^T P exakt (drei Sechspunktamplituden)"),
    "Klammerzerlegung z89": ("E/C", "m_Z o (eps c^+ (x) Z*) = CSW-Fortsetzung; wörtlich zwischen zwei Nähten auf der Diagonale null, "
                                    "mit der Partnernaht als Referenz falsch"),
    "Farbklebung": ("E/L", "U(2) der Transferströme mit Gesamtladung: Vollständigkeit, Spurklebung; su(2) allein nicht"),
    "Austauschfeld": ("O", "fehlt nativ: freies Twistorfeld a (Homogenität 0), b (Homogenität -4) mit Kinetik tr(b dbar a), "
                           "Propagator dbar^(-1) in CSW-Eichung zwischen den Nahtgeraden, Kopplung J.a an die Nahtströme, Summe über "
                           "die Lagen; Kandidat: a <-> b ist Tomitas n -> -n - 4 (z85), die bilineare Paarung <J psi|phi>"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("T2 ist bedingt beantwortet: Die Komposition der Nähte ist nicht nativ CSW, aber bis auf eine genau benannte Zutat. "
          "Die Vergleichsmaschine steht: BCFW aus Dreipunktamplituden gibt Parke-Taylor (31 MHV-Lagen bis n = 6), die "
          "Anti-MHV-Formel und die geschlossene Split-Helizitätsformel, unabhängig von der Verschiebung, und die CSW-Regeln geben "
          "A(1-,2-,3-,4+,5+) (4 Diagramme) und drei Sechspunkt-NMHV-Amplituden (6, 9, 8 Diagramme) exakt und unabhängig von eta, "
          "während einzelne Diagramme von eta abhängen. Nativ vorhanden sind vier Zutaten. Erstens die Referenz der "
          "Off-shell-Fortsetzung: Jeder Punkt Z* = I_pi W der flachen Geraden aus z89 gibt eta~ = omega(Z*), und Punkt und "
          "Modulphase sind Eichung, sodass die offene Phase aus z98 nicht stört; I_omega ist nicht translationsinvariant, die native "
          "E-invariante Wahl I_AdS hat keinen Punkt (Pfaffsche -4), rahmenabhängige Referenzen und Spinoren der falschen Chiralität "
          "scheitern, die gespiegelten Regeln mit Anti-MHV-Vertizes funktionieren mit dem dualen Punkt. Zweitens die "
          "z89-Zerlegung B = m_Z o (eps c^+ (x) c^+): Sie ist genau die CSW-Fortsetzung lambda_P = sum c(k) eps(c^+(k), Z*), wenn "
          "ihr zweiter Faktor der Referenztwistor ist; als Klebung zweier Nähte auf der Diagonale ist sie null (kollinear), mit der "
          "Partnernaht als Referenz falsch. Drittens die Form des Propagators als Kinematik: Der Cauchykern der Naht L_x an ihren "
          "Schnittpunkten mit zwei Ebenen durch L_y ist -1/(<ab> det(x - y)), die Penrose-Transformation von 1/((A.Z)(B.Z)), "
          "euklidisch 4 pi^2/P^2; schneidende Nähte (z61) geben genau die Residuen P^2 A -> A_L A_R. Viertens die Farbklebung über "
          "die U(2) der Transferströme einschließlich der Gesamtladung. Nicht nativ ist der Träger des Propagators. Die gemeinsame "
          "CAR hat zwischen Nähten keinen Kern (blockdiagonale Fermisee-Projektion, exakt 0 im Vakuum und bei beta = 4 pi); der "
          "Hüpfterm der Transferströme ist modendiagonal, lagenblind und im Vakuum durch die Quarterlücke unterdrückt; die nativen "
          "Quarterkerne tragen den Twist (s1/s2)^(+-1/4) und bei beta = 4 pi die Besetzung c(1/4) = 0.041, sodass der Inzidenzkern "
          "nur im Vakuumgrenzfall ohne Twist eine Funktion von (x - y)^2 ist; die Normierung legt erst die Faktorisierung fest. Die "
          "fehlende Zutat ist ein freies Austauschfeld auf dem Twistorraum zwischen den Nahtgeraden: das Paar a (Homogenität 0), "
          "b (Homogenität -4) mit der Kinetik tr(b dbar a) und dem Propagator dbar^(-1) in CSW-Eichung (Referenztwistor auf I_pi), "
          "linear an die Nahtströme gekoppelt (J.a, eine Kopplung g) und über die Nahtlagen summiert. Unter Hypothese T (z110) "
          "liefern die Nähte dann die Vertexhälfte der Twistorwirkung (log det auf den Geraden), das Austauschfeld die gaußsche "
          "Kinetikhälfte. Folgen: K15 ist präzisiert (Yang-Mills-Bäume = Nahtvertizes unter T + Austauschfeld + I_pi; die "
          "Kinetikhälfte fehlt). Für T3/T5 ist der gesuchte Zweikörperterm nicht lokal in einem Nahtpaar (Doppeltransfer, z99), "
          "sondern der Austausch J^a Delta J^a zwischen Nähten an verschiedenen Lagen, der beim Ausintegrieren eines gaußschen "
          "Feldes entsteht; die Nichtgaußizität käme ganz aus den Nähten. Die Szenarioachse Wechselwirkung bekommt einen bedingten "
          "Wert 'CSW-Komposition' mit den Prämissen T, P_BF (Austauschfeld) und I_pi, die Kopplung bleibt frei. Kandidat für die "
          "native Kinetik: a <-> b ist genau Tomitas n -> -n - 4 (z85), also die bilineare Paarung <J psi|phi> der Sektoren "
          "D = -+2. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\n   (Laufzeit {time.time() - T0:.1f} s)")
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z111_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
