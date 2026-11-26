"""Z176: Schritt 4 der Kernstrategie aus z173. Ursprung: K19 (Naht mit Defekt im Zustand), K20 (Elternoperator), Nahtpaar.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Bausteine bleiben unberührt, nichts wird befördert.
Statuskürzel: E exakt, N numerisch, C bedingt, L Lesart, O offen.
Bestand: z168 (Zustand e^(-2 pi K), K = L_0 - 1/3, <K> = 0 bei tau = i), z169/z170 (Uhrsektoren von E8 = D5 x A3,
Diskriminante Z4 x Z4 mit q = (5x^2 + 3y^2)/8; 27 Spiegelvervollständigungen; Grundniveau des Defekts +1/24; Vorschlag
M13), z130/z174 (Standardform Naht + Tomita-Partner, Delta^(-i/2) = Paaruhr; Z4-Verklebung anomal auf einem Paar, frei über
zwei: P_paar), z113 (Quartermodell).

Teil A  M13: alle Randbedingungen des Sektorenobjekts (16 einseitige Uhr-Einschübe in Raum und Zeit, 27
        Spiegelvervollständigungen) gegen <T> = 0 und Positivität; Grundniveaus exakt; Eichsumme der Paaruhr.
Teil B  Zwei Blätter als Nahtpaar: Blatttausch auf den Uhrsektoren, Vier-Marken-Funktion der E8-Naht, Anomalie der
        einseitigen Uhr im Sektorenobjekt, Eichsumme über ein Paar von Nähten.
Teil C  K20: Was (Delta, J) der Standardform rekonstruieren (endliches Modell von z174) und der einseitige modulare
        Operator des Sektorenobjekts als Elternoperator.
"""
import itertools
import json
import time
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import expm, null_space

HERE = Path(__file__).parent
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


mp.mp.dps = 40
I, PI = mp.mpc(0, 1), mp.pi
HALB = mp.mpf(1)/2
R_ = lambda x: mp.mpf(x.numerator)/x.denominator if isinstance(x, Fr) else mp.mpf(x)
ex = lambda r: mp.exp(2*I*PI*R_(r))
ns = lambda x, n=12: mp.nstr(x, n)


def theta(a, b, tau):
    """theta[a;b](tau) = sum_n exp(i pi tau (n+a)^2 + 2 pi i (n+a) b)."""
    a, b = R_(a), R_(b)
    nmax = int(mp.sqrt((mp.mp.dps*2.31 + 12)/(PI*tau.imag))) + 3
    return mp.fsum(mp.exp(I*PI*tau*(n + a)**2 + 2*I*PI*(n + a)*b) for n in range(-nmax, nmax + 1))


def eta(tau):
    kmax = int(mp.sqrt((mp.mp.dps*2.31 + 12)/(3*PI*tau.imag))) + 3
    return mp.exp(2*I*PI*tau/24)*mp.fsum((-1)**k*mp.exp(I*PI*tau*k*(3*k - 1)) for k in range(-kmax, kmax + 1))


def bausteine(tau):
    return theta(HALB, 0, tau), theta(0, 0, tau), theta(0, HALB, tau), eta(tau)       # theta_2, theta_3, theta_4, eta


def dn(n, t2, t3, t4, e):
    """Charaktere von D_n Stufe 1 nach Klasse 0 (o), 1 (s), 2 (v), 3 (c)."""
    o, v, s = (t3**n + t4**n)/(2*e**n), (t3**n - t4**n)/(2*e**n), t2**n/(2*e**n)
    return {0: o, 1: s, 2: v, 3: s}


def chars(tau):
    """Die 16 Uhrsektoren (x, y) von D5 x A3 (A3 = D3): x Klasse von D5, y Klasse von A3."""
    b = bausteine(tau)
    d5, d3 = dn(5, *b), dn(3, *b)
    return {(x, y): d5[x]*d3[y] for x in range(4) for y in range(4)}


SEK = [(x, y) for x in range(4) for y in range(4)]
HD5 = {0: Fr(0), 1: Fr(5, 8), 2: Fr(1, 2), 3: Fr(5, 8)}
HA3 = {0: Fr(0), 1: Fr(3, 8), 2: Fr(1, 2), 3: Fr(3, 8)}
C_NAHT = Fr(1, 3)                                                          # c/24 für c = 8
qf = lambda a: Fr(5*a[0]**2 + 3*a[1]**2, 8)                                 # quadratische Form der Diskriminante (z170)
bl = lambda a, b: Fr(5*a[0]*b[0] + 3*a[1]*b[1], 4)                          # ihre Bilinearform
neg = lambda a: tuple((-t) % 4 for t in a)
plus = lambda g, h: tuple((i + j) % 4 for i, j in zip(g, h))
TP = [mp.mpc('0.137', '1.093'), mp.mpc('-0.211', '0.947')]


def erzeuge(gens, n):
    S = {tuple([0]*n)}
    neu = True
    while neu:
        neu = False
        for s in list(S):
            for g in gens:
                x = plus(s, g)
                if x not in S:
                    S.add(x)
                    neu = True
    return frozenset(S)


# =================================================================== Teil A
say("Teil A: M13, Randbedingungen des Sektorenobjekts gegen <T> = 0 und Positivität")

# --- A1 Sektorenobjekt
a1 = []
for tau in TP:
    c, c1, cs = chars(tau), chars(tau + 1), chars(-1/tau)
    t2, t3, t4, e = bausteine(tau)
    E4 = (t2**8 + t3**8 + t4**8)/2
    a1.append(max(abs(sum(c[(k, k)] for k in range(4)) - E4/e**8), abs(sum(c[(k, (-k) % 4)] for k in range(4)) - E4/e**8),
                  abs(sum(I**k*c[(k, k)] for k in range(4)) - (eta(2*tau)/eta(4*tau))**4),
                  max(abs(c1[a] - ex(qf(a) - C_NAHT)*c[a]) for a in SEK),
                  max(abs(cs[a] - sum(ex(-bl(a, b))*c[b] for b in SEK)/4) for a in SEK)))
ZAHLEN["A1_max_Fehler"] = float(max(a1))
check("A1 [E/N, 40 Stellen] Sektorenobjekt: 16 Uhrsektoren (x, y) von D5 x A3 mit h = (5x^2 + 3y^2)/8 mod 1. Die E8-Naht ist "
      "die Summe der vier Sektoren (k, k) (ebenso (k, -k)) und gleich E_4/eta^8; die Uhr R = i^x gibt Sp R q^(L_0 - 1/3) = "
      "(eta(2 tau)/eta(4 tau))^4 (z169). T wirkt mit e^(2 pi i (q - 1/3)), S mit e^(-2 pi i b(a, b))/4, b die Bilinearform "
      "von q: Beide Matrizen folgen aus der quadratischen Form allein",
      max(a1) < mp.mpf(10)**-35 and all((qf(a) - HD5[a[0]] - HA3[a[1]]) % 1 == 0 for a in SEK))

# --- A2 Grundniveaus exakt
K0 = {a: HD5[a[0]] + HA3[a[1]] - C_NAHT for a in SEK}
klasse = lambda a: (a[1] - a[0]) % 4                                        # Verdrillung relativ zur Verklebung (k, k)
grund = {m: min(K0[a] for a in SEK if klasse(a) == m) for m in range(4)}
spiegelnull = [(a, b) for a in SEK for b in SEK if K0[a] + K0[b] == 0]
ZAHLEN["A2_Grundniveaus L_0 - 1/3 je Sektor"] = {str(a): str(K0[a]) for a in SEK}
ZAHLEN["A2_Grundniveau je Verdrillung R^m"] = {str(m): str(v) for m, v in grund.items()}
check("A2 [E] Grundniveaus L_0 - 1/3 aller 16 Sektoren, exakt: unverdrillt -1/3, R- und R^3-verdrillt +1/24 (die 4 bei "
      "h = 3/8; die 16 liegt bei +7/24), R^2-verdrillt +1/6. Kein Sektor hat das Grundniveau 0, und kein Paar aus Sektor und "
      "Spiegelsektor hat L_0 + L_0-quer - 2/3 = 0 (h liegt in Z/8, 2/3 nicht). Die Energiebilanz null gibt es in keinem "
      "Sektor als Grundzustand",
      grund == {0: Fr(-1, 3), 1: Fr(1, 24), 2: Fr(1, 6), 3: Fr(1, 24)} and K0[(1, 0)] == Fr(7, 24)
      and all(v != 0 for v in K0.values()) and not spiegelnull)


# --- A3 einseitige positive Zustände
def Zab(a, b, tau):
    """Sp über den R^a-verdrillten Sektor mit Einschub R^b: sum_x i^(b x) chi_(x, x+a)."""
    c = chars(tau)
    return sum(I**(b*x)*c[(x, (x + a) % 4)] for x in range(4))


def boost(F, y):
    """Re <K> = -(1/2 pi) d/dy log |F(i y)|."""
    return -mp.diff(lambda yy: mp.log(abs(F(mp.mpc(0, yy)))), y)/(2*PI)


YG = [mp.mpf(k)/10 for k in range(4, 31)]
a3 = {a: [boost(lambda t, a=a: Zab(a, 0, t), y) for y in YG] for a in range(4)}
null0 = mp.findroot(lambda y: boost(lambda t: Zab(0, 0, t), y), mp.mpf('1.1'))
ZAHLEN["A3_min <K> je Verdrillung (y = 0.4..3.0)"] = {str(a): float(min(v)) for a, v in a3.items()}
ZAHLEN["A3_Nullstelle unverdrillt"] = float(null0)
check("A3 [E/N] Einseitige positive Zustände (Spur über den R^a-verdrillten Sektor ohne Einschub, Gewicht e^(-2 pi y K)): "
      f"<K> = 0 gibt es nur unverdrillt, genau bei y = 1 (tau = i; Nullstelle {ns(null0, 6)}). In den verdrillten Sektoren "
      "ist <K> für jede Temperatur größer als das Grundniveau (+1/24 bzw. +1/6), also nie null: Für sich genommen "
      "erfüllt der Defektsektor A* nicht",
      abs(null0 - 1) < mp.mpf(10)**-25 and all(min(a3[a]) > R_(grund[a]) for a in (1, 2, 3))
      and min(a3[0]) < 0 < max(a3[0]))
say(f"   {zeit()}")

# --- A4 alle 16 Randbedingungen
a4, nullst = {}, {}
for a in range(4):
    for b in range(4):
        F = lambda t, a=a, b=b: Zab(a, b, t)
        YV = [y + mp.mpf('0.03') for y in YG]                               # versetztes Gitter: keine Nullstelle auf einem Punkt
        werte = [boost(F, y) for y in YV]
        a4[(a, b)] = (F(mp.mpc(0, 1)), boost(F, mp.mpf(1)))
        ws = [YV[i] for i in range(len(YV) - 1) if werte[i]*werte[i + 1] < 0]
        nullst[(a, b)] = [float(mp.findroot(lambda y: boost(F, y), w + mp.mpf('0.05'))) for w in ws]
diag = {k: a4[(k, k)][0] for k in range(4)}
quer = {(a, b): a4[(a, b)][0] for (a, b) in ((1, 3), (3, 1))}
bilanz0 = sorted(ab for ab, v in a4.items() if abs(v[1]) < mp.mpf(10)**-25)
ZAHLEN["A4_Z[a;b](i)"] = {str(k): ns(v[0], 10) for k, v in a4.items()}
ZAHLEN["A4_Re<K> bei tau = i"] = {str(k): ns(v[1], 8) for k, v in a4.items()}
ZAHLEN["A4_Nullstellen von Re<K> auf der imaginären Achse"] = {str(k): [round(x, 9) for x in v] for k, v in nullst.items()}
check("A4 [N, 40 Stellen] Universum der 16 einseitigen Randbedingungen (Verdrillung R^a im Raum, Einschub R^b in der Zeit): "
      "Bei tau = i verschwindet Re <K> genau für die sechs Bedingungen a = +-b, mit den Spuren 12 (0, 0), 4 (2, 2), 2 - 2i "
      "((1, 1), (3, 3)) und 2 + 2i ((1, 3), (3, 1)). Nur (0, 0) ist ein positives Funktional; jeder Einschub R^b mit "
      "b != 0 ist eine uhrgraduierte Spur und kein Zustand. Die Bilanz null im Defektsektor (a = 1, b = +-1) gibt es "
      "also nur als graduierte Spur",
      bilanz0 == sorted({(a, b) for a in range(4) for b in range(4) if (a - b) % 4 == 0 or (a + b) % 4 == 0})
      and max(abs(diag[0] - 12), abs(diag[1] - (2 - 2*I)), abs(diag[2] - 4), abs(diag[3] - (2 - 2*I)),
              abs(quer[(1, 3)] - (2 + 2*I)), abs(quer[(3, 1)] - (2 + 2*I))) < mp.mpf(10)**-30)

# --- A5 Kompensation des Grundniveaus
E_LEITER = Fr(-1, 24)                                                       # eine chirale Bose-Leiter (z113: osz)
E_QUARTERNAHT = 2*Fr(1, 2)*Fr(1, 4)**2 + 2*E_LEITER                         # Ladungen (Q -+ 1/4)^2/2 und zwei Leitern: -1/48
E_NS, E_R = Fr(-1, 48), Fr(1, 24)                                           # reelles Fermion, NS und R
zusatz_c = {str(m): str(24*grund[m]) for m in (1, 2)}
ZAHLEN["A5"] = {"Bose-Leiter": str(E_LEITER), "chirale Hälfte einer Quarternaht (z113)": str(E_QUARTERNAHT),
                "10 NS + 6 R reelle Fermionen": str(10*E_NS + 6*E_R), "nötige Zusatzladung c je Verdrillung": zusatz_c}
check("A5 [E, Zuordnung L] Kompensation: Das Grundniveau des Defekts ist +1/24 = 10 (-1/48) + 6 (+1/24) (zehn reelle "
      "Fermionen NS, sechs R). Es hebt sich exakt gegen eine einzelne Bose-Leiter (-1/24): mit c = 9 ist h = 3/8 = c/24. "
      "Ebenso gegen die chiralen Hälften zweier Quarternähte (je -1/48, z113). Innerhalb des Sektorenobjekts (c = 8) gibt "
      "es diesen Beitrag nicht; der R^2-Sektor bräuchte c = 12",
      grund[1] + E_LEITER == 0 and grund[1] + 2*E_QUARTERNAHT == 0 and 10*E_NS + 6*E_R == grund[1]
      and E_QUARTERNAHT == Fr(-1, 48) and Fr(3, 8) == Fr(9, 24) and zusatz_c == {"1": "1", "2": "4"})
befund("A5: Eine native einzelne Leiter (oder ein natives Paar von Quarternähten) neben der E8-Naht ist im Bestand "
       "nachgewiesen", False)

# --- A6 Spiegelvervollständigungen
G4 = list(itertools.product(range(4), repeat=4))
QQ = lambda g: (qf(g[:2]) - qf(g[2:])) % 1
BB = lambda g, h: (bl(g[:2], h[:2]) - bl(g[2:], h[2:])) % 1
gefunden, rand = set(), {frozenset({(0, 0, 0, 0)})}
gesehen = set(rand)
while rand:
    neu = set()
    for Sg in rand:
        if len(Sg) == 16:
            gefunden.add(Sg)
            continue
        for g in G4:
            if g in Sg or BB(g, g) != 0 or any(BB(g, s) != 0 for s in Sg):
                continue
            S2 = erzeuge(list(Sg) + [g], 4)
            if len(S2) <= 16 and S2 not in gesehen and all(BB(p, r) == 0 for p in S2 for r in S2):
                gesehen.add(S2)
                neu.add(S2)
    rand = neu
VERV = sorted(gefunden, key=lambda L: sorted(L))
CI = chars(mp.mpc(0, 1))


def ZL(L, y):
    c = chars(mp.mpc(0, y))
    return mp.fsum(mp.re(c[g[:2]])*mp.re(c[g[2:]]) for g in L)                        # Charaktere sind auf der imaginären Achse reell


zeilen = []
for L in VERV:
    Z1 = ZL(L, mp.mpf(1))
    zeilen.append({"bosonisch": all(QQ(g) == 0 for g in L),
                   "jeder Sektor einmal": all(sum(1 for g in L if g[:2] == a) == 1 for a in SEK),
                   "Ströme links": any(g[2:] == (0, 0) and g[:2] != (0, 0) for g in L),
                   "16 lokal": any(g[:2] in ((1, 0), (3, 0)) for g in L),
                   "Z": Z1, "E": -mp.diff(lambda yy: mp.log(ZL(L, yy)), 1)/(2*PI),
                   "w": mp.fsum(mp.re(CI[g[:2]])*mp.re(CI[g[2:]]) for g in L if (g[0] + g[1]) % 2)/Z1})
W2 = mp.sqrt(2)
bos = [z for z in zeilen if z["bosonisch"]]
voll = [z for z in bos if z["jeder Sektor einmal"]]
rest = [z for z in bos if not z["jeder Sektor einmal"]]
W_DEF = (12*W2 - 13)/14
ZAHLEN["A6"] = {"Vervollständigungen": len(VERV), "bosonisch": len(bos),
                "Tabelle (bosonisch, jeder Sektor einmal, Ströme links, 16 lokal, Z(i), Gewicht uhrverdrillt)":
                    sorted([str(z["bosonisch"]), str(z["jeder Sektor einmal"]), str(z["Ströme links"]), str(z["16 lokal"]),
                            ns(z["Z"], 10), ns(z["w"], 10)] for z in zeilen),
                "max |<L_0 + L_0-quer - 2/3>| bei tau = i": float(max(abs(z["E"]) for z in zeilen))}
check("A6 [E/N] Mit Spiegelpartner: Es gibt 27 Vervollständigungen, 10 bosonisch (z170). Jede ist S-invariant, also gilt "
      "<L_0 + L_0-quer - 2/3> = 0 bei tau = i für alle 27, mit positiven Gewichten: A* allein wählt keine aus. Sechs "
      "bosonische haben Z(i) = 144 und keinen uhrverdrillten Sektor (dort sind Ströme der Erweiterung lokal). Die vier, in "
      "denen jeder Sektor einmal vorkommt, haben Z(i) = 80 + 48 sqrt(2), enthalten die 16 lokal und tragen im Zustand bei "
      f"tau = i das Gewicht (12 sqrt(2) - 13)/14 = {ns(W_DEF, 8)} in den uhrverdrillten Sektoren",
      len(VERV) == 27 and len(bos) == 10 and max(abs(z["E"]) for z in zeilen) < mp.mpf(10)**-25
      and len(rest) == 6 and all(abs(z["Z"] - 144) < mp.mpf(10)**-30 and z["w"] == 0 and not z["16 lokal"] for z in rest)
      and len(voll) == 4 and all(abs(z["Z"] - (80 + 48*W2)) < mp.mpf(10)**-30 and abs(z["w"] - W_DEF) < mp.mpf(10)**-30
                                 and z["16 lokal"] and not z["Ströme links"] for z in voll))
say(f"   {zeit()}")

# --- A7 Eichsumme der Paaruhr
L_E8 = frozenset((k, k, l, l) for k in range(4) for l in range(4))          # E8 (x) E8-quer: Naht mit Tomita-Partner
DIAG = frozenset((x, y, x, y) for x in range(4) for y in range(4))
CDIAG = frozenset((x, y, (-x) % 4, (-y) % 4) for x in range(4) for y in range(4))


def eichsumme(L, v, B):
    return erzeuge([v] + sorted(g for g in L if B(v, g) == 0), len(v))


paaruhr = {v: (QQ(v), eichsumme(L_E8, v, BB)) for v in ((1, 0, 1, 0), (1, 0, 3, 0), (1, 0, 0, 0))}
pk = {a: mp.re(CI[a])**2/(80 + 48*W2) for a in SEK}
ek = {a: -mp.diff(lambda yy, a=a: mp.log(mp.re(chars(mp.mpc(0, yy))[a])), 1)/(2*PI) for a in SEK}
kl = {m: (mp.fsum(pk[a] for a in SEK if klasse(a) == m), mp.fsum(2*pk[a]*ek[a] for a in SEK if klasse(a) == m))
      for m in range(4)}
ZAHLEN["A7_Gewicht und Energiebeitrag je Verdrillung (diagonal, tau = i)"] = \
    {str(m): [ns(v[0], 10), ns(v[1], 10)] for m, v in kl.items()}
check("A7 [E, entscheidend] Eichsumme der Paaruhr: Die Naht mit Tomita-Partner ist E8 (x) E8-quer (keine Uhrsektoren). Die "
      "Paaruhr R (x) R-quer^-1 ist anomaliefrei (q - q-quer = 0); summiert man über sie (neutrale Felder plus modulare "
      "Ergänzung), entsteht genau die diagonale Vervollständigung mit allen 16 Sektoren, je nach Spiegelbit a (x) a-quer "
      "oder a (x) Ca-quer. Die einseitige Uhr allein ist anomal (q = 5/8). Im Zustand bei tau = i hebt der Beitrag der "
      f"uhrverdrillten Sektoren ({ns(kl[1][1] + kl[3][1], 6)}) zusammen mit dem R^2-Anteil ({ns(kl[2][1], 6)}) den der "
      f"unverdrillten ({ns(kl[0][1], 6)}) auf",
      paaruhr[(1, 0, 1, 0)] == (0, DIAG) and paaruhr[(1, 0, 3, 0)] == (0, CDIAG) and paaruhr[(1, 0, 0, 0)][0] == Fr(5, 8)
      and DIAG in gefunden and CDIAG in gefunden and abs(mp.fsum(v[1] for v in kl.values())) < mp.mpf(10)**-25
      and abs(kl[1][0] + kl[3][0] - W_DEF) < mp.mpf(10)**-30 and kl[0][1] < 0 < kl[1][1])
befund("M13: Eine Randbedingung wählt den Defektsektor als Grundzustand oder als einseitigen positiven Zustand mit "
       "<T> = 0", False)
befund("K19 bedingt: Unter P_std (Paaruhr superselektiert, gelesen als Eichsumme) und A* (Basiswechsel sind Symmetrien, "
       "<T> = 0) enthält der Zustand die uhrverdrillten Sektoren zwangsläufig, mit Gewicht (12 sqrt(2) - 13)/14", True)

# =================================================================== Teil B
say("\nTeil B: Zwei Blätter als Nahtpaar")

# --- B1 Blatttausch auf den Uhrsektoren
w = lambda r: np.exp(2j*np.pi*float(r))
nrm = lambda X: float(np.abs(X).max())
SM = np.array([[w(-bl(a, b))/4 for b in SEK] for a in SEK])
TM = np.diag([w(qf(a) - C_NAHT) for a in SEK])
perm = lambda f: np.array([[1.0 if a == f(b) else 0.0 for b in SEK] for a in SEK])
CM = perm(neg)
RM = np.diag([1j**a[0] for a in SEK])
EINS = perm(lambda a: (a[0], (-a[1]) % 4))                                  # Ladungskonjugation einer Seite (Bit P_IF)
nichtentartet = all(any(bl(d, b) % 1 != 0 for b in SEK) for d in SEK if d != (0, 0))
b1 = {"S^2 - C": nrm(SM @ SM - CM), "(ST)^3 - S^2": nrm(np.linalg.matrix_power(SM @ TM, 3) - SM @ SM),
      "S^4 - 1": nrm(np.linalg.matrix_power(SM, 4) - np.eye(16)), "C R C - R^-1": nrm(CM @ RM @ CM - RM.conj()),
      "C - R^2": nrm(CM - RM @ RM), "C - einseitige Konjugation": nrm(CM - EINS)}
fix = [a for a in SEK if neg(a) == a]
ZAHLEN["B1"] = dict(b1, Fixsektoren_von_C=[str(a) for a in fix])
check("B1 [E] Blatttausch: Die Decktransformation der zweiblättrigen Sphäre ist die elliptische Involution des Torus, "
      "also S^2. Auf den 16 Uhrsektoren ist S^2 = (ST)^3 = C die volle Ladungskonjugation (x, y) -> (-x, -y), mit "
      "C R C = R^-1 und vier Fixsektoren. Der Blatttausch ist weder R^2 = (-1)^F (diagonal, keine Permutation) noch die "
      "Konjugation einer Seite (das Bit P_IF) noch der Seitentausch F aus z174 (der Original und Partner tauscht)",
      nichtentartet and max(b1["S^2 - C"], b1["(ST)^3 - S^2"], b1["S^4 - 1"], b1["C R C - R^-1"]) < 1e-12
      and b1["C - R^2"] > 0.5 and b1["C - einseitige Konjugation"] > 0.5 and len(fix) == 4)

# --- B2 Vier-Marken-Funktion
b2 = []
for tau in [mp.mpc(0, 1)] + TP:
    t2, t3, t4, e = bausteine(tau)
    lam = t2**4/t3**4
    E4 = (t2**8 + t3**8 + t4**8)/2
    wick = 1/lam + 1/(1 - lam) - 1
    b2.append((lam, wick, max(abs(E4 - t3**8*(1 - lam + lam**2)), abs(e**24 - (t2*t3*t4)**8/256),
                              abs(E4/(t2*t4)**4 - wick), abs(wick - (1 - lam + lam**2)/(lam*(1 - lam))))))
ZAHLEN["B2 (lambda, Vier-Marken-Funktion, Fehler)"] = [[ns(x[0], 10), ns(x[1], 10), float(x[2])] for x in b2]
check("B2 [E, Überlagerungsformel aus der Literatur; 40 Stellen] Eine E8-Naht auf dem Torus ist dasselbe wie ein Paar von "
      "E8-Nähten auf der Sphäre mit vier Blatttausch-Marken vom Gewicht c/16 = 1/2: j^(1/3)/(2^8 x (1 - x))^(1/3) = "
      "(1 - x + x^2)/(x (1 - x)) = 1/x + 1/(1 - x) - 1 mit x = lambda(tau). Das ist die Vierpunktfunktion eines freien "
      "Majorana-Fermions (drei Wick-Paarungen). tau = i ist die quadratische Lage x = 1/2 der Marken; der Wert dort ist 3",
      max(x[2] for x in b2) < mp.mpf(10)**-35 and abs(b2[0][0] - HALB) < mp.mpf(10)**-35
      and abs(b2[0][1] - 3) < mp.mpf(10)**-35 and Fr(8, 16) == Fr(1, 2))

# --- B3 Quartermodell: Decktransformation der Holonomie
mp.mp.dps = 50
QU = Fr(1, 4)
L_lad = lambda u, t, tau: mp.exp(-2*I*PI*R_(t)*R_(u))*theta(u, t, tau)
osz = lambda r, tau: mp.exp(2*I*PI*R_(r)/24)/eta(tau + R_(r))


def naht(ein):
    def hol(par):
        ta, tb, ra, rb = par
        return lambda tau: L_lad(-QU, ta, tau)*L_lad(QU, tb, tau)*osz(ra, tau)*osz(rb, tau)
    fl, fr = hol(ein[0]), hol(ein[1])
    return lambda tau: fl(tau)*mp.conj(fr(tau))


TP3 = [mp.mpc('0.137', '1.093'), mp.mpc('-0.211', '0.947'), mp.mpc('0.389', '1.211')]
FX = naht(((QU, -QU, -QU, -QU), (0, 0, -QU, -QU)))                          # X = R i^(-2h) (z174)
FD = naht(((-QU, QU, -QU, -QU), (0, 0, -QU, -QU)))                          # R^-1 i^(-2h): Bild unter der Decktransformation
FI = naht(((-QU, QU, QU, QU), (0, 0, QU, QU)))                              # X^-1
b3 = {"Deckbild gegen X": float(min(abs(FX(t) - FD(t))/abs(FX(t)) for t in TP3)),
      "X^-1 gegen X": float(min(abs(FX(t) - FI(t))/abs(FX(t)) for t in TP3)),
      "X^-1 gegen Paritätsbild": float(max(abs(FI(t) - mp.conj(FX(-mp.conj(t))))/abs(FX(t)) for t in TP3))}
mp.mp.dps = 40
z174 = json.loads((HERE / "z174_ergebnis.json").read_text())
c174 = lambda frag: any(k.startswith(frag) and v for k, v in z174["checks"].items())
ch174 = z174["zahlen"]["B2_charaktere"]
frei174 = {k: [N for N in ("1", "2", "3", "4") if ch174[k][N]] for k in ch174}
ZAHLEN["B3"] = b3
check("B3 [N, 50 Stellen] Quartermodell (z113, z174): Die Decktransformation kehrt die Holonomie um, lässt aber die "
      "Drehung i^(-2h) stehen. Die Spur mit R^-1 i^(-2h) ist eine andere Funktion als die mit X = R i^(-2h) (relative "
      f"Abweichung mindestens {b3['Deckbild gegen X']:.2f}); X^-1 ist das Paritätsbild von X. Die Uhr des Quartermodells "
      "ist gerade unter der Ladungskonjugation, die die Quarter tauscht. Die zwei Blätter einer Naht ergeben also eine "
      "einzelne Torusspur und nicht die diagonale Summe X_1 X_2 aus z174",
      b3["Deckbild gegen X"] > 0.1 and b3["X^-1 gegen X"] > 0.1 and b3["X^-1 gegen Paritätsbild"] < 1e-40)

# --- B4 Anomalie der einseitigen Uhr im Sektorenobjekt
HEB = [((1 + k) % 4, k) for k in range(4)]                                  # Hebungen der Uhr: (1, 0) + (k, k)
iso = {N: [vs for vs in itertools.product(HEB, repeat=N) if sum(qf(v) for v in vs) % 1 == 0] for N in range(1, 5)}
spiegel_frei = all((N*(qf(v) - qf(v))) % 1 == 0 for N in range(1, 5) for v in HEB)
ZAHLEN["B4"] = {"q der Hebungen": [str(qf(v) % 1) for v in HEB], "isotrope Hebungen je N": {str(N): len(v) for N, v in iso.items()},
                "z174 frei für N": frei174}
check("B4 [E] Anomalie im Sektorenobjekt: Die Uhr als Verschiebung hat die Hebungen (1, 0) + (k, k) mit q = 5/8, 7/8, 1/8, "
      "3/8; keine ist isotrop, 2q liegt in Z/4: eine Anomalie der Ordnung 2. Über N Nähte gibt es isotrope Hebungen genau "
      f"für gerades N ({ {N: len(v) for N, v in iso.items()} }). Die Paaruhr mit Spiegelpartner ist für jedes N frei. Das "
      "ist dasselbe Muster wie im Quartermodell von z174 (X und X' frei genau für N = 2, 4; X_p für N = 1..4), dort "
      "numerisch, hier exakt",
      all(len(iso[N]) == 0 for N in (1, 3)) and len(iso[2]) == 4 and len(iso[4]) == 64 and spiegel_frei
      and all((2*qf(v)) % Fr(1, 4) == 0 and qf(v) % Fr(1, 4) != 0 for v in HEB)
      and frei174["X = R i^(-2h)"] == ["2", "4"] and frei174["X' = R i^(2h)"] == ["2", "4"]
      and frei174["X_p = R Rq^-1 i^(-2h)"] == ["1", "2", "3", "4"] and c174("B2 [E, entscheidend] Anomalie"))

# --- B5 Eichsumme über ein Paar von Nähten
B2N = lambda g, h: (bl(g[:2], h[:2]) + bl(g[2:], h[2:])) % 1
KREUZ = frozenset((k, l, (-l) % 4, (-k) % 4) for k in range(4) for l in range(4))
paar = {vs[0] + vs[1]: eichsumme(L_E8, vs[0] + vs[1], B2N) for vs in iso[2]}
STUFEN = {("D", 0): [(Fr(0), 1), (Fr(1), 45)], ("D", 2): [(Fr(1, 2), 10)], ("D", 1): [(Fr(5, 8), 16)], ("D", 3): [(Fr(5, 8), 16)],
          ("A", 0): [(Fr(0), 1), (Fr(1), 15)], ("A", 2): [(Fr(1, 2), 6)], ("A", 1): [(Fr(3, 8), 4)], ("A", 3): [(Fr(3, 8), 4)]}


def stroeme(L):
    """Zahl der Zustände vom Gewicht 1 (Ströme) eines Gitters aus Klassen (x1, y1, x2, y2)."""
    n = 0
    for g in L:
        for kombi in itertools.product(STUFEN[("D", g[0])], STUFEN[("A", g[1])], STUFEN[("D", g[2])], STUFEN[("A", g[3])]):
            if sum(h for h, _ in kombi) == 1:
                n += int(np.prod([m for _, m in kombi]))
    return n


grund_paar = 2*Fr(3, 8) - 2*C_NAHT
ZAHLEN["B5"] = {"isotrope Hebungen": [str(v) for v in sorted(paar)], "Ströme E8 x E8": stroeme(L_E8),
                "Ströme nach Eichsumme": stroeme(KREUZ), "Grundniveau des verdrillten Paarsektors": str(grund_paar)}
check("B5 [E] Eichsumme über ein Paar von Nähten: Für jede der vier isotropen Hebungen ist die Eichsumme von E8 x E8 nach "
      "der diagonalen Uhr wieder E8 x E8 (496 Ströme), aber über Kreuz verklebt: Klassen (k, l, -l, -k), die D5 jeder Naht "
      "hängt an der A3 der anderen. Die 16 der einen Naht ist dort zusammen mit der 4-quer der anderen ein Strom vom "
      "Gewicht 5/8 + 3/8 = 1; der verdrillte Paarsektor hat das Grundniveau 3/4 - 2/3 = +1/12. Das Paar hebt die Anomalie "
      "auf, senkt aber die Energiebilanz des Defekts nicht auf null",
      len(paar) == 4 and all(L == KREUZ for L in paar.values()) and KREUZ != L_E8 and stroeme(L_E8) == 496
      and stroeme(KREUZ) == 496 and (1, 0, 0, 3) in KREUZ and grund_paar == Fr(1, 12))
befund("P_paar aus den zwei Blättern von P1: Das Nahtpaar aus z174 ist die Blattverdopplung", False)
befund("P_paar strukturell: Die einseitige Uhr hat in beiden Modellen eine Anomalie der Ordnung 2; wer sie summiert, "
       "braucht eine gerade Zahl von Nähten, mindestens ein Paar", True)
say(f"   {zeit()}")

# =================================================================== Teil C
say("\nTeil C: K20, ein Elternoperator für Naht, Uhr und Partner")


def modell(K, kfun):
    """Standardform von z174 (Teil A): K Moden je Seite, rho = e^(-kfun(Q))/Z, Omega = vec(rho^(1/2))."""
    d = 2**K
    qa = np.array([bin(i).count("1") for i in range(d)])
    kk = kfun(qa)
    p = np.exp(-kk)
    p = p/p.sum()
    F = np.zeros((d*d, d*d))
    for i in range(d):
        for j in range(d):
            F[j*d + i, i*d + j] = 1
    return {"d": d, "qa": qa, "p": p, "Om": np.diag(np.sqrt(p)).reshape(-1), "F": F, "QA": np.repeat(qa, d),
            "QB": np.tile(qa, d), "logD": np.array([kk[j] - kk[i] for i in range(d) for j in range(d)])}   # log Delta


KMS = lambda q: np.pi*q
KMS2 = lambda q: 2*np.pi*q**2 - np.pi*q
MOD = {"K = pi Q": modell(2, KMS), "K = 2 pi Q^2 - pi Q": modell(2, KMS2)}
rng = np.random.default_rng(176)
c1, c2, c3 = {}, {}, {}
for name, m in MOD.items():
    d, D, F, logD = m["d"], m["d"]**2, m["F"], m["logD"]
    paaruhr_f = 1j**((m["QA"] - m["QB"]) % 4)
    c1[name] = (nrm(np.exp(-0.5j*logD) - paaruhr_f), abs(float(m["Om"] @ (logD*m["Om"]))))
    # reell-lineares System: h antihermitesch, [h, log Delta] = 0, J h J = h mit J = F o Konjugation
    basis = []
    for i in range(D):
        for j in range(D):
            for c in (1, 1j):
                Ee = np.zeros((D, D), complex)
                Ee[i, j] = c
                basis.append(Ee)

    def zwang(h):
        v = np.concatenate([(h + h.conj().T).reshape(-1), (h*logD[None, :] - logD[:, None]*h).reshape(-1),
                            (F @ h.conj() @ F - h).reshape(-1)])
        return np.concatenate([v.real, v.imag])

    kern = null_space(np.array([zwang(Ee) for Ee in basis]).T, rcond=1e-10)
    vek = lambda X: np.array([t for z in X.reshape(-1) for t in (z.real, z.imag)])
    abl, kan = [], []
    for i in range(d):
        for j in range(d):
            for c in (1, 1j):
                k = np.zeros((d, d), complex)
                k[i, j] = c
                k = k - k.conj().T
                if nrm(k) == 0:
                    continue
                gen = np.kron(k, np.eye(d)) + np.kron(np.eye(d), k.conj())    # Erzeuger von Y -> u Y u^+
                abl.append(vek(gen))
                if nrm(k*m["p"][None, :] - m["p"][:, None]*k) < 1e-14:
                    kan.append(vek(gen))
    rg = lambda M: int(np.linalg.matrix_rank(M, tol=1e-9))
    abl, kan = np.array(abl).T, np.array(kan).T
    schnitt = rg(abl) + rg(kern) - rg(np.concatenate([abl, kern], axis=1))
    # ein Element des Kommutanten senkrecht zu den kanonischen: bricht es den Kegel?
    Uk, sk, _ = np.linalg.svd(kan, full_matrices=False)
    Qk = Uk[:, sk > 1e-9]
    restk = kern - Qk @ (Qk.T @ kern)
    col = restk[:, int(np.argmax(np.linalg.norm(restk, axis=0)))]
    h = (col[0::2] + 1j*col[1::2]).reshape(D, D)
    h = h/nrm(h)
    Wk = expm(0.7*h)
    mind = 0.0
    for _ in range(40):
        Xr = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
        Y2 = (Wk @ (Xr @ Xr.conj().T).reshape(-1)).reshape(d, d)
        mind = min(mind, -nrm(Y2 - Y2.conj().T), float(np.linalg.eigvalsh((Y2 + Y2.conj().T)/2).min()))
    c2[name] = {"Kommutant von (Delta, J)": int(kern.shape[1]), "kanonisch": rg(kan), "Kegelableitungen": rg(abl),
                "Schnitt": schnitt, "Antihermitizität": nrm(h + h.conj().T), "Kegelbruch": mind}
    # einseitige Uhr
    RA = 1j**(m["QA"] % 4)
    eig = {}
    for i in range(D):
        eig.setdefault(round(float(logD[i]), 9), set()).add(complex(np.round(RA[i], 9)))
    rho_uhr = m["p"].astype(complex)**(-0.5j)
    phase = rho_uhr/(1j**(m["qa"] % 4))
    c3[name] = (all(len(v) == 1 for v in eig.values()), nrm(phase - phase[0]))
ZAHLEN["C1 (Delta^(-i/2) - Paaruhr, <log Delta>)"] = {k: list(v) for k, v in c1.items()}
ZAHLEN["C2"] = c2
ZAHLEN["C3 (R_A Funktion von Delta, rho^(-i/2) = Phase mal i^Q)"] = {k: [v[0], v[1]] for k, v in c3.items()}
check("C1 [E] Endliches Modell (zwei Moden je Seite, beide Zustände aus z174): Delta^(-i/2) ist die Paaruhr i^(Q_A - Q_B), "
      "J der Partner. Der modulare Erzeuger log Delta = K_B - K_A hat im Standardvektor den Erwartungswert null, für "
      "jeden Zustand: In der Partnerlesart ist <T> = 0 automatisch und wählt keinen Sektor aus",
      all(v[0] < 1e-12 and v[1] < 1e-12 for v in c1.values()))
check("C2 [E] (Delta, J) allein legen die Teilung in Naht und Partner nicht fest: Die Unitären, die mit Delta und J "
      f"vertauschen, bilden eine Gruppe der Dimension {c2['K = pi Q']['Kommutant von (Delta, J)']} (zweiter Zustand: "
      f"{c2['K = 2 pi Q^2 - pi Q']['Kommutant von (Delta, J)']}); nur 5 Richtungen davon erhalten den natürlichen Kegel "
      "(die kanonischen Implementierungen u (x) u-quer mit [u, rho] = 0). Ein Element außerhalb bricht den Kegel. Zur "
      "Rekonstruktion der Naht braucht man zusätzlich den Kegel, also die Algebra",
      c2["K = pi Q"]["Kommutant von (Delta, J)"] == 32 and c2["K = 2 pi Q^2 - pi Q"]["Kommutant von (Delta, J)"] == 24
      and all(v["kanonisch"] == 5 and v["Schnitt"] == 5 and v["Kegelableitungen"] == 15 and v["Antihermitizität"] < 1e-9
              and v["Kegelbruch"] < -1e-6 for v in c2.values()))
check("C3 [E] Die einseitige Uhr R_A = i^(Q_A) ist keine Funktion von Delta (sie nimmt auf Eigenräumen von Delta "
      "verschiedene Werte an). Mit der Algebra ist sie rho^(-i/2) bis auf eine konstante Phase, für beide Zustände; "
      "eingesteckt sind dabei der Zustand selbst (beta = 4 pi, Quarterenergie 1/4) und die Teilung in Naht und Partner",
      all((not v[0]) and v[1] < 1e-12 for v in c3.values()))

# --- C4 Sektorenobjekt: der einseitige modulare Operator
versch = lambda dd: perm(lambda a: plus(a, dd))
c4 = {"S R S^-1 - Fusion mit (1, 0)": nrm(SM @ RM @ SM.conj().T - versch((1, 0))),
      "S^-1 R S - Fusion mit (3, 0)": nrm(SM.conj().T @ RM @ SM - versch((3, 0)))}
monodromie = [d_ for d_ in SEK if nrm(np.diag([w(bl(a, d_)) for a in SEK]) - RM) < 1e-12]
S_aus_T = nrm(np.array([[w(-(qf(plus(a, b)) - qf(a) - qf(b)))/4 for b in SEK] for a in SEK]) - SM)
bild = sorted({plus((k, k), (1, 0)) for k in range(4)})
isotrop = [a for a in SEK if qf(a) % 1 == 0]
lagr = sorted({tuple(sorted(erzeuge([a], 2))) for a in isotrop if len(erzeuge([a], 2)) == 4})
ZAHLEN["C4"] = dict(c4, Monodromievektor_der_Uhr=[str(d_) for d_ in monodromie], S_aus_T=S_aus_T,
                    Bild_der_E8_Sektoren=[str(a) for a in bild], isotrope_Untergruppen_der_Ordnung_4=[str(L) for L in lagr])
check("C4 [E] Sektorenobjekt: Der einseitige modulare Operator des Zustands e^(-2 pi K) ist rho^(-i) = e^(2 pi i K) = T (bis "
      "auf eine konstante Phase). Aus T allein folgt mit der Fusion die Bilinearform und daraus S (A1), damit S^2 = C "
      "(Blatttausch) und S R S^-1 = Fusion mit dem Sektor (1, 0) der 16: Die Uhr und das Erzeugen des Defekts sind "
      "dieselbe Antwort, um die Vierteldrehung S gedreht; die E8-Sektoren (k, k) gehen dabei auf die R-verdrillten. Die "
      "Uhr ist die Monodromie mit genau einem Sektor, (1, 0). Die E8-Verklebung ist eine der zwei isotropen Untergruppen "
      "von T, (k, k) oder (k, -k)",
      max(c4.values()) < 1e-12 and monodromie == [(1, 0)] and S_aus_T < 1e-12
      and bild == [(0, 3), (1, 0), (2, 1), (3, 2)] and all(klasse(a) == 3 for a in bild) and len(lagr) == 2
      and tuple(sorted((k, k) for k in range(4))) in lagr and tuple(sorted((k, (-k) % 4) for k in range(4))) in lagr)
befund("K20: Das Paar (Delta, J) allein rekonstruiert Naht, Uhr und Partner", False)
befund("K20 bedingt: Im Sektorenobjekt ist der einseitige modulare Operator T = rho^(-i) ein gemeinsamer Elternoperator "
       "für Uhr, Defekt und Blatttausch; eingesteckt sind die 16 Sektoren mit ihrer Fusion (dieselbe Bedingung wie bei "
       "K19), der Sektor der 16 und die Teilung in Naht und Partner", True)
befund("Komplettlösung", False)

urteil = (
    "Schritt 4 ist gerechnet. (1) M13 ist in der angekündigten Form negativ: Kein Sektor hat das Grundniveau null, und "
    "kein einseitiger positiver Zustand im Defektsektor erfüllt <T> = 0; mit Uhr-Einschub verschwindet die Bilanz genau "
    "für die Randbedingungen a = +-b, aber das sind graduierte Spuren, keine Zustände. Das Grundniveau +1/24 hebt "
    "sich exakt gegen eine einzelne Bose-Leiter (c = 9, h = c/24); eine native Leiter ist nicht nachgewiesen. (2) Positiv, "
    "bedingt: Die Frage 'warum ist der Defekt im Grundzustand' ist falsch gestellt, weil der Zustand von A* kein "
    "Grundzustand ist, sondern der S-symmetrische Zustand bei tau = i. Jede der 27 Spiegelvervollständigungen erfüllt dort "
    "<T> = 0. Welche gilt, entscheidet die Paaruhr: Summiert man über sie (P_std: die Paaruhr ist superselektiert), wird "
    "aus E8 (x) E8-quer genau die diagonale Vervollständigung mit allen 16 Sektoren, und der Zustand trägt das Gewicht "
    "(12 sqrt(2) - 13)/14 = 28,4 % in den uhrverdrillten Sektoren. K19 ist damit auf P_std und A* zurückgeführt, mit der "
    "Lesart 'Superselektion ist Eichsumme'; geschlossen ist es nicht. (3) Zwei Blätter: Der Blatttausch ist S^2 = C, die "
    "volle Ladungskonjugation, nicht R^2 und nicht der Seitentausch. Eine E8-Naht auf dem Torus ist ein Paar von E8-Nähten "
    "auf der Sphäre mit vier Marken, die freie Majorana-Fermionen vom Gewicht 1/2 sind; tau = i ist die quadratische Lage. "
    "Das Nahtpaar aus z174 ist aber nicht diese Blattverdopplung: Die zwei Blätter ergeben eine einzelne Torusspur. "
    "Strukturell bleibt: Die einseitige Uhr hat in beiden Modellen eine Anomalie der Ordnung 2; wer sie summiert, braucht "
    "ein Paar. P_paar schrumpft auf 'die Verklebung ist einseitig und wird summiert', hergeleitet ist es nicht. (4) K20: "
    "(Delta, J) allein reichen nicht, die Teilung in Naht und Partner braucht den Kegel. Im Sektorenobjekt ist der "
    "einseitige modulare Operator T = rho^(-i) ein Elternoperator für Uhr, Defekt und Blatttausch, unter derselben "
    "Bedingung wie K19. Keine Komplettlösung, nichts befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z176_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil},
                                                    indent=1, ensure_ascii=False, default=str))
(HERE / "ausgabe_z176.txt").write_text("\n".join(x for x in OUT if "Laufzeit bisher" not in x) + "\n")
say(f"{zeit()}")
