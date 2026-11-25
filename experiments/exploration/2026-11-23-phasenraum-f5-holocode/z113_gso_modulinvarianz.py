"""Z113: Tiefe Welle T4 'Folgt die Uhr-verdrehte Hebung P_GSO aus der Modulinvarianz der Zustandssumme über N Nähte?'

Leitlinie 'Rekonstruktion statt Auswahl': Eine Projektion zählt nur dann als nativ, wenn die Modulinvarianz der
Zustandssumme sie erzwingt, also wenn die Spur mit Einschub in derselben SL(2,Z)-Bahn liegt wie die Spur ohne Einschub und
die Bahnsumme den relativen Koeffizienten festlegt (so entsteht im Vergleichsziel Stringtheorie die GSO-Projektion).
Firewalls: DIMENSION.UPLIFT.FIREWALL.01 (Helizität ist Kinematik der radialen Einpaarschicht), GRAV.NONCIRCULAR.01. Gelöste
Bausteine (16 je Generation, E8 aus der Spinstruktursumme, Flavour) werden benutzt und nicht neu geöffnet; nichts wird
befördert. Bestand: z97 (G = U_2pi (x) R^2, (-1)^(2h) = R^2, Prämisse P_GSO; Z_SU = i^(2h-2); J-kovariant nur
G_tot = U_2pi (-1)^(F + F-quer); kein Element von Modulgruppe mal Ladungsdrehungen ist G), z106 (P_CAR gibt R^2-Statistik),
z83 (Paarschalen H_m mit L_0 = m in beiden Quartersektoren, Uhr R = e^(2 pi i L_0)), z92 (E^nat, Spin x_Z2 Z8 an den Marken),
z94 (|theta_3/eta|^2 ist S-invariant), v1/z10/z38 (16 Majorana-Weyl, GSO-Summe, E8 mit flachem (P,P)-Sektor).

Modell (exakt): Das Original trägt zwei komplexe chirale Fermionen, Quarter a (Moden Z + 1/4) und b (Z + 3/4), also den
Spektralfluss u = -1/4, +1/4 der NS-Fermionen: Gewicht (Q + u)^2/2 = E_lad(alpha, Q) + 1/32. Der Partner (Tomita-Bild)
ist antiholomorph mit denselben Gewichten. Einschübe als exakte Spuren: R = i^(Q_a - Q_b), R^2 = (-1)^F (Charakteristik der
Thetafunktion), U_2pi = (-1)^(L_0^osc - L-quer_0^osc) auf den Oszillatoren. Auf der radialen Einpaarschicht
H_m (x) H-quer_mquer ist das (-1)^(m - mquer) = (-1)^(2h) aus z97; auf Mehrpaarzuständen ist es multiplikativ. Der Einschub
ist auf dem Oszillatorfaktor die halbe Verschiebung tau -> tau + 1/2. Bahnen werden numerisch exakt bestimmt (mpmath, eigene
Theta- und Eta-Reihen): Zwei Bahnelemente F(g tau) gelten als gleich, wenn ihr Quotient an drei Testpunkten dieselbe
Einheitswurzel ist. Daraus folgen die Monomialdarstellung von S und T, die Charaktere (Anomalie), die erzwungenen relativen
Koeffizienten und die Positivität der vollständigen Zustandssumme (spinaufgelöste Fouriermoden).

Teil 1  Bausteine und Sektoren: Theta mit Charakteristik und eta (S, T, Identität für eta(tau + 1/2)); Ladungssektor als
        Spektralfluss, Uhr = 2D-Drehung (T fügt R ein, T^2 fügt R^2 ein); die vier Spinstrukturen; U_2pi auf der
        Einpaarschicht (Weyl-Zählung) und als halbe Drehung der Oszillatoren; ihre Bahn (lange Nähte, Moden in Z/4).
Teil 2  Welche Projektionen sind modulinvariant? Kalibrierung 16 -> E8. Chirales Original über N Nähte (N = 0 mod 4, die
        c = 8-Naht = (E8)_1, Uhr dort eingefroren). Nichtchirale Naht: Bahnzugehörigkeit, Charaktere und erzwungene
        Koeffizienten für U (x) 1, G = U R^2, G_tot, U R R-quer, einseitiges R^2 und die geometrische Hebung; Positivität der
        vollständigen Zustandssumme für N = 1..4 und mit nur einer radialen Leiter; die U-Summe auf der c = 8-Naht.
Teil 3  Bilanz, Urteil für P_GSO, Folgen für die Statistikachse und für T4.
"""
import json
import time
from collections import defaultdict, deque
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
mp.mp.dps = 25
I = mp.mpc(0, 1)
PI = mp.pi
TP = [mp.mpc('0.137', '1.093'), mp.mpc('-0.211', '0.947'), mp.mpc('0.389', '1.211')]   # Testpunkte
NULL = mp.mpf(10)**-16


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
    return f"(Laufzeit bisher {time.time() - T0:.0f} s)"


# ------------------------------------------------------------------ Theta, eta, Modulgruppe
def R_(x):
    return mp.mpf(x.numerator)/x.denominator if isinstance(x, Fr) else mp.mpf(x)


def theta(a, b, tau):
    """theta[a;b](tau) = sum_n exp(i pi tau (n+a)^2 + 2 pi i (n+a) b), q = e^(2 pi i tau)."""
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


GEN = {"S": (0, -1, 1, 0), "T": (1, 1, 0, 1)}
EINS = (1, 0, 0, 1)


def wurzel(c, nmax=4608):
    """c = exp(2 pi i r) mit rationalem r; gibt r mod 1 oder None."""
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
    return c if all(abs(x/y - c) < mp.mpf(10)**-12*abs(c) for x, y in zip(v[1:], w[1:])) else None


def bahn(F, maxn=600):
    """SL(2,Z)-Bahn {F o g}: Vertreter, Werte, Permutationen und Phasen (Einheitswurzeln) von S und T."""
    reps, vals = [EINS], [werte(F, EINS)]
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
                perm[gn][k], ph[gn][k] = len(reps) - 1, Fr(0)
                assert len(reps) <= maxn
        k += 1
    assert all(x is not None for d in ph.values() for x in d.values())
    return {"F": F, "reps": reps, "vals": vals, "perm": perm, "ph": ph, "n": len(reps)}


def lokalisiere(B, H):
    """H = e^(2 pi i r) (F o g_j)? -> (j, r) oder None."""
    v = werte(H, EINS)
    for j, w in enumerate(B["vals"]):
        c = proportional(v, w)
        if c is not None:
            return j, wurzel(c)
    return None


def t_zyklen(B):
    seen, out = set(), []
    for i in range(B["n"]):
        if i not in seen:
            c, j = 0, i
            while j not in seen:
                seen.add(j)
                c += 1
                j = B["perm"]["T"][j]
            out.append(c)
    return sorted(out)


def gemeinsam(faktoren):
    """Bahn von prod_k (F_k o g)^(e_k) (e_k < 0: komplex konjugiert) über die Faktorbahnen; Phasen addiert."""
    start = tuple(s for _, s, _ in faktoren)
    idx, states, kanten = {start: 0}, [start], []
    dq = deque([start])
    while dq:
        st = dq.popleft()
        for gn in ("S", "T"):
            nst = tuple(B["perm"][gn][s] for (B, _, _), s in zip(faktoren, st))
            phase = sum((e*B["ph"][gn][s] for (B, _, e), s in zip(faktoren, st)), Fr(0)) % 1
            if nst not in idx:
                idx[nst] = len(states)
                states.append(nst)
                dq.append(nst)
            kanten.append((idx[st], gn, idx[nst], phase))
    return states, idx, kanten


def konsistenz(n, kanten, chi):
    """Invariante Kombination sum_j e^(2 pi i V_j) F_j mit F|gen = e^(2 pi i chi(gen)) F; bad = Widersprüche (Anomalie)."""
    V = [None]*n
    V[0] = Fr(0)
    adj = defaultdict(list)
    for (i, gn, j, r) in kanten:
        adj[i].append((gn, j, r))
    stack = [0]
    while stack:
        i = stack.pop()
        for (gn, j, r) in adj[i]:
            if V[j] is None:
                V[j] = (V[i] + r - chi[gn]) % 1
                stack.append(j)
    return V, {(V[i] + r - chi[gn] - V[j]) % 1 for (i, gn, j, r) in kanten} - {Fr(0)}


def charaktere(n, kanten):
    """alle Charaktere chi(T) = k/24, chi(S) = j/4, für die die Bahnsumme widerspruchsfrei ist."""
    out = []
    for k in range(24):
        for j in range(4):
            chi = {"T": Fr(k, 24), "S": Fr(j, 4)}
            V, bad = konsistenz(n, kanten, chi)
            if not bad:
                out.append((chi, V))
    return out


def chistr(ch):
    return [f"T {c['T']}, S {c['S']}" for c, _ in ch]


# ------------------------------------------------------------------ Bausteine der Naht
QU = Fr(1, 4)
HA = Fr(1, 2)


def L_lad(u, t, tau):
    """sum_Q q^((Q+u)^2/2) e^(2 pi i t Q): exakte Spur von e^(2 pi i t Q) über die Ladungen eines Quarters."""
    return mp.exp(-2*I*PI*R_(t)*R_(u))*theta(u, t, tau)


def osz(s, tau):
    """Oszillatorspur einer chiralen Bose-Leiter: s = 0: 1/eta; s = 1: mit (-1)^(L_0^osc): q^(-1/24)/prod(1 - (-q)^n)."""
    return 1/eta(tau) if s == 0 else mp.exp(I*PI/24)/eta(tau + HA)


def h_orig(ta, tb, sa=0, sb=0, geo=0):
    """holomorphe Spur des Originals: Ladungseinschub e^(2 pi i (ta Q_a + tb Q_b)), U auf den Oszillatoren von a und/oder b,
    geo = volle halbe Drehung (auch der Ladungsgitter)."""
    def f(tau):
        tl = tau + HA if geo else tau
        return L_lad(-QU, ta, tl)*L_lad(QU, tb, tl)*osz(max(sa, geo), tau)*osz(max(sb, geo), tau)
    return f


TYP = {"1": (0, 0, 0, 0), "R": (QU, -QU, 0, 0), "R2": (HA, HA, 0, 0), "R3": (-QU, QU, 0, 0),
       "U": (0, 0, 1, 1), "UR": (QU, -QU, 1, 1), "UR2": (HA, HA, 1, 1), "UR3": (-QU, QU, 1, 1)}

# =================================================================== Teil 1
say("Teil 1: Bausteine und Sektoren der Nahtzustandssumme auf dem Torus")

# --- 1.1 Transformationsgesetze (exakt, mpmath)
fehler = []
for t in TP:
    for a, b in ((QU, Fr(0)), (-QU, HA), (HA, QU), (Fr(1, 8), Fr(3, 8))):
        fehler.append(abs(theta(a, b, t + 1) - mp.exp(-I*PI*R_(a)*(R_(a) + 1))*theta(a, a + b + HA, t)))
        fehler.append(abs(theta(a, b, -1/t) - mp.sqrt(-I*t)*mp.exp(2*I*PI*R_(a)*R_(b))*theta(b, -a, t)))
        fehler.append(abs(theta(a, b + 1, t) - mp.exp(2*I*PI*R_(a))*theta(a, b, t)))
    fehler.append(abs(eta(-1/t) - mp.sqrt(-I*t)*eta(t)))
    fehler.append(abs(eta(t + 1) - mp.exp(I*PI/12)*eta(t)))
    fehler.append(abs(eta(t + HA) - mp.exp(I*PI/24)*eta(2*t)**3/(eta(t)*eta(4*t))))
    q = mp.exp(2*I*PI*t)
    prod = mp.mpc(1)
    for n in range(1, 400):
        prod *= 1 - (-q)**n
    fehler.append(abs(osz(1, t) - q**(-mp.mpf(1)/24)/prod))
f_max = float(max(fehler))
ZAHLEN["transformationsgesetze_max_fehler"] = f_max
check("[E] Bausteine (mpmath, 25 Stellen, drei Testpunkte): theta[a;b](tau + 1) = e^(-i pi a(a+1)) theta[a; a+b+1/2], "
      "theta[a;b](-1/tau) = sqrt(-i tau) e^(2 pi i ab) theta[b;-a], theta[a;b+1] = e^(2 pi i a) theta[a;b], eta(-1/tau) = "
      "sqrt(-i tau) eta, eta(tau+1) = e^(i pi/12) eta, eta(tau + 1/2) = e^(i pi/24) eta(2tau)^3/(eta(tau) eta(4tau)) und die "
      f"Oszillatorspur mit (-1)^N = q^(-1/24)/prod(1 - (-q)^n) = e^(i pi/24)/eta(tau + 1/2); größter Fehler {f_max:.0e}",
      f_max < 1e-18)

# --- 1.2 Ladungssektor = Spektralfluss, Uhr = 2D-Drehung
E_lad = lambda al, Q: Fr(Q*Q, 2) + (al - HA)*Q                      # z82/z92/z97
spektral = all(Fr((4*Q - 1)**2, 32) - Fr(1, 32) == E_lad(Fr(1, 4), Q) and
               Fr((4*Q + 1)**2, 32) - Fr(1, 32) == E_lad(Fr(3, 4), Q) for Q in range(-8, 9))
uhr = all(((E_lad(Fr(1, 4), qa) + E_lad(Fr(3, 4), qb)) - Fr((qa - qb) % 4, 4)) % 1 == 0
          for qa in range(-8, 9) for qb in range(-8, 9))
uhr2 = all((2*(E_lad(Fr(1, 4), qa) + E_lad(Fr(3, 4), qb)) - Fr((qa + qb) % 2, 2)) % 1 == 0
           for qa in range(-8, 9) for qb in range(-8, 9))
B0 = bahn(h_orig(*TYP["1"]))
BU = bahn(h_orig(*TYP["U"]))
BGEO = bahn(h_orig(0, 0, 0, 0, geo=1))
LOK = {}
for p, par in TYP.items():
    for nm, B in (("B0", B0), ("BU", BU)):
        x = lokalisiere(B, h_orig(*par))
        if x is not None:
            LOK[p] = (nm, x[0], x[1])
            break
T1 = B0["perm"]["T"][0]
T2 = B0["perm"]["T"][T1]
ZAHLEN["ladungssektor"] = {"Bahn h_1 (Original, Ladung + Oszillatoren)": B0["n"], "T-Zyklen": t_zyklen(B0),
                           "Lage der Einschübe": {k: [v[0], v[1], str(v[2])] for k, v in LOK.items()},
                           "T h_1 -> Element": T1, "T^2 h_1 -> Element": T2}
check("[E] Ladungssektor als Spektralfluss: (Q - 1/4)^2/2 - 1/32 = E_lad(1/4, Q) und (Q + 1/4)^2/2 - 1/32 = E_lad(3/4, Q) "
      "exakt (|Q| <= 8). Die Uhr ist die 2D-Drehung: e^(2 pi i L_0) = i^(Q_a - Q_b) = R und e^(4 pi i L_0) = (-1)^(Q_a + Q_b) "
      "= R^2 exakt auf allen Ladungen (z83 Teil 5), und auf den Spuren bildet T die Spur ohne Einschub auf die mit R ab "
      f"(Element {T1} = Lage von R {LOK['R'][1]}), T^2 auf die mit R^2 (Element {T2} = Lage von R^2 {LOK['R2'][1]}). Die "
      f"Bahn hat {B0['n']} Elemente (T-Zyklen {t_zyklen(B0)}): R^2 = (-1)^F ist auf dem Original eine Modultransformation, "
      "wie (-1)^F = T im NS-Sektor des Vergleichsziels",
      spektral and uhr and uhr2 and LOK["R"][:2] == ("B0", T1) and LOK["R2"][:2] == ("B0", T2) and B0["n"] == 6)

# --- 1.3 Spinstrukturen (Dirac-Paar, z94)
z_ns = lambda tau: abs(theta(0, 0, tau)/eta(tau))**2
BNS = bahn(z_ns)
s_inv = max(float(abs(z_ns(t) - z_ns(-1/t))/z_ns(t)) for t in TP)
ns_tilde = lokalisiere(BNS, lambda tau: abs(theta(0, HA, tau)/eta(tau))**2)
r_ns = lokalisiere(BNS, lambda tau: abs(theta(HA, 0, tau)/eta(tau))**2)
rr = float(max(abs(theta(HA, HA, t)) for t in TP))
ZAHLEN["spinstrukturen"] = {"S-Abweichung |theta3/eta|^2": s_inv, "Bahn": BNS["n"],
                            "NS mit (-1)^F, R ohne Einschub in der Bahn": [ns_tilde is not None, r_ns is not None],
                            "|theta_1| (R mit (-1)^F)": rr}
check(f"[E] Spinstrukturen: |theta_3/eta|^2 ist S-invariant (rel. Abweichung {s_inv:.0e}, z94); ihre Bahn hat genau "
      f"{BNS['n']} Elemente, NS-NS, NS mit (-1)^F und R-NS; R mit (-1)^F verschwindet (|theta_1| = {rr:.0e}, Nullmode). Die "
      "Spinstruktursumme ist eine Bahnsumme, ihr relatives Vorzeichen ist festgelegt",
      s_inv < 1e-20 and BNS["n"] == 3 and ns_tilde is not None and r_ns is not None and rr < 1e-20)

# --- 1.4 U_2pi: Einpaarschicht, halbe Drehung, lange Nähte
M = 12
weyl, schicht = {}, {}
for n in range(1, 6):
    weyl[n] = sum(m*mb for m in range(1, M) for mb in range(1, M) if m + mb == 2*n + 1 and m - mb == 1)
gen_fehler = 0
for m in range(1, M):
    for mb in range(1, M):
        u = (-1)**(m - mb)
        gen_fehler += abs(u - (-1)**((m + mb) % 2))              # e^(2 pi i E) = (-1)^(2h) auf der Schicht
ungerade_halb = all(((m - mb) % 2 == 1) == ((-1)**(m - mb) == -1) for m in range(1, M) for mb in range(1, M))
u_in_b0 = lokalisiere(B0, h_orig(*TYP["U"]))
s_bild = []
for t in TP:
    s_bild.append(abs(osz(1, -1/t) - mp.sqrt(2)*(-I*t)**(-mp.mpf('0.5'))*eta(t)*eta(t/4)/eta(t/2)**3))
s_bild = float(max(s_bild))
ZAHLEN["U_2pi"] = {"Weyl-Zählung h = +-1/2 bei E = n + 1/2": [weyl[n] for n in range(1, 6)],
                   "Bahn h_U": BU["n"], "T-Zyklen h_U": t_zyklen(BU), "h_U in Bahn(h_1)": u_in_b0 is not None,
                   "Bahn geometrische halbe Drehung": BGEO["n"], "T-Zyklen geometrisch": t_zyklen(BGEO),
                   "S-Bild osz_U = sqrt2 (-i tau)^-1/2 eta(tau) eta(tau/4)/eta(tau/2)^3, Fehler": s_bild}
check("[E] U_2pi auf der Einpaarschicht H_m (x) H-quer_mquer (L_0^osc = m in beiden Quartersektoren, z83): "
      "(-1)^(m - mquer) = (-1)^(2h) = e^(2 pi i E), ungerade genau für halbzahlige h; Weyl-Zählung h = +-1/2 bei "
      f"E = n + 1/2: {[weyl[n] for n in range(1, 6)]} = n(n+1) (z97). Fortgesetzt ist U_2pi = (-1)^(L_0^osc - L-quer_0^osc), "
      "multiplikativ, weil L_0^osc die Summe der Paarenergien ist",
      [weyl[n] for n in range(1, 6)] == [n*(n + 1) for n in range(1, 6)] and gen_fehler == 0 and ungerade_halb)
check("[E, entscheidend] U_2pi ist keine Modultransformation der Naht: Der Einschub ist auf den Oszillatoren tau -> tau + 1/2; "
      f"h_U liegt nicht in der Bahn von h_1 ({B0['n']} Elemente), seine Bahn hat {BU['n']} Elemente mit T-Zyklen "
      f"{t_zyklen(BU)}. Das S-Bild der Oszillatorspur mit U ist exakt sqrt2 (-i tau)^(-1/2) eta(tau) eta(tau/4)/eta(tau/2)^3 "
      f"(Fehler {s_bild:.0e}), Moden in Z/4: Die verdrehten Sektoren der 2 pi-Drehung sind lange Nähte (4-fache Überlagerung "
      f"des Kreises, Hecke-Sektoren). Die volle geometrische halbe Drehung hat eine eigene Bahn ({BGEO['n']} Elemente, "
      f"T-Zyklen {t_zyklen(BGEO)})",
      u_in_b0 is None and BU["n"] == 24 and t_zyklen(BU) == [4]*6 and s_bild < 1e-18 and BGEO["n"] == 24)
befund("Teil 1: Auf dem Original sind die Uhr R und R^2 = (-1)^F die 2D-Drehungen T und T^2, also Modultransformationen; die "
       "Spinstruktursumme ist eine Bahnsumme. Die 4D-2 pi-Drehung ist dagegen die halbe Drehung der Oszillatoren, keine "
       "Modultransformation; ihre Sektoren sind lange Nähte", True)
say(f"   {zeit()}")

# =================================================================== Teil 2
say("\nTeil 2: Welche Projektionen sind modulinvariant? Kalibrierung, N Nähte, Positivität")

# --- 2.1 Kalibrierung 16 -> E8
kal = {}
for n in (2, 4, 8, 12, 16):
    F = (lambda n: (lambda tau: (theta(0, 0, tau)/eta(tau))**n))(n)
    B = bahn(F)
    _, idx, kn = gemeinsam([(B, 0, 1)])
    ch = charaktere(B["n"], kn)
    j, r = lokalisiere(B, (lambda n: (lambda tau: (theta(0, HA, tau)/eta(tau))**n))(n))
    eintrag = []
    for chi, V in ch:
        lam = (V[idx[(j,)]] - r) % 1
        summe = max(abs(mp.fsum(mp.exp(2*I*PI*R_(V[idx[(k,)]]))*B["vals"][k][p] for k in range(B["n"])))
                    for p in range(3))
        eintrag.append({"chi": f"T {chi['T']}", "lambda (-1)^F": str(lam), "Summe": float(summe)})
    kal[n] = eintrag
ZAHLEN["kalibrierung_16_E8 (n komplexe Fermionen)"] = kal
gso_ok = [n for n, e in kal.items() if any(x["lambda (-1)^F"] == "0" and x["Summe"] > 1e-6 for x in e)]
check("[E] Kalibrierung 16 -> E8 (gelöst, nur benutzt): Für n komplexe chirale Fermionen ist die Bahn von (theta_3/eta)^n "
      "{theta_3, theta_4, theta_2}^n; eine widerspruchsfreie Summe mit dem GSO-Vorzeichen +1 für die Spur mit (-1)^F und "
      f"nicht verschwindender Summe gibt es genau für n in {gso_ok} (16 bzw. 32 Majorana-Weyl, Charakter T = -c/24); n = 4 "
      "ist widerspruchsfrei, aber die Summe verschwindet (Jacobi), n = 2 ist anomal, n = 12 nur mit umgekehrtem Vorzeichen. "
      "Dort erzwingt die Modulinvarianz die Projektion, weil (-1)^F = T in der Bahn liegt",
      gso_ok == [8, 16] and not kal[2] and kal[4] and all(x["Summe"] < 1e-15 for x in kal[4])
      and all(x["lambda (-1)^F"] == "1/2" for x in kal[12]))

# --- 2.2 chirales Original über N Nähte, die c = 8-Naht
chiral = {}
for N in range(1, 9):
    _, idx, kn = gemeinsam([(B0, 0, N)])
    chiral[N] = chistr(charaktere(B0["n"], kn))
frei = [N for N, c in chiral.items() if c]
_, idx4, kn4 = gemeinsam([(B0, 0, 4)])
(chi4, V4), = charaktere(B0["n"], kn4)
Z_uhr = lambda tau: mp.fsum(mp.exp(2*I*PI*R_(V4[idx4[(k,)]]))*B0["F"](act(B0["reps"][k], tau))**4 for k in range(B0["n"]))
Z_E8 = lambda tau: (theta(HA, 0, tau)**8 + theta(0, 0, tau)**8 + theta(0, HA, tau)**8)/(2*eta(tau)**8)
e8_quot = [Z_uhr(t)/Z_E8(t) for t in TP]
e8_fehler = float(max(abs(x - 1) for x in e8_quot))
lam_e8 = {p: str((V4[idx4[(LOK[p][1],)]] - 4*LOK[p][2]) % 1) for p in ("R", "R2", "R3")}
_, idxu4, knu4 = gemeinsam([(BU, 0, 4)])
chu4 = [(c, V) for c, V in charaktere(BU["n"], knu4) if c["T"] == chi4["T"] and c["S"] == chi4["S"]]
Vu4 = chu4[0][1]
Z_Uc = lambda tau: mp.fsum(mp.exp(2*I*PI*R_(Vu4[idxu4[(k,)]]))*BU["F"](act(BU["reps"][k], tau))**4 for k in range(BU["n"]))
u_quot = [Z_Uc(t)/Z_E8(t) for t in TP]
u_fehler = float(max(abs(x - u_quot[0]) for x in u_quot))
lam_u4 = {p: str((Vu4[idxu4[(LOK[p][1],)]] - 4*LOK[p][2]) % 1) for p in ("UR", "UR2", "UR3")}
ZAHLEN["chirales_original"] = {"Charaktere je N": chiral, "anomaliefrei für N": frei,
                               "Uhrensumme N = 4 / Z_E8": [complex(x) for x in e8_quot],
                               "relative Koeffizienten (Phase/2pi) von R, R^2, R^3 bei N = 4": lam_e8,
                               "U-Bahnsumme N = 4 / Z_E8": complex(u_quot[0]), "Abweichung": u_fehler,
                               "relative Koeffizienten UR, UR^2, UR^3 bei N = 4": lam_u4}
check("[E] Chirales Original (Ladung und Oszillatoren, c = 2 je Naht) über N Nähte: Die Bahnsumme ist genau für N = 0 mod 4 "
      f"widerspruchsfrei (N = 1..8: {frei}), also c = 8 mal k. Bei N = 4 ist die Uhrensumme exakt (E8)_1: Z_Uhr/Z_E8 = 1 "
      f"(Abweichung {e8_fehler:.0e}, Koeffizienten 1, 248, 4124, ...). Ihre Koeffizienten {lam_e8} (Phase/2 pi von R, R^2, "
      "R^3) projizieren auf R = -i: Auf der c = 8-Naht ist die Uhr eingefroren, R^2 = -1, also G = U R^2 = -U",
      frei == [4, 8] and e8_fehler < 1e-15 and lam_e8 == {"R": "1/4", "R2": "1/2", "R3": "3/4"})
check("[E] U-Summe auf der c = 8-Naht: Die Bahnsumme von h_U^4 mit dem E8-Charakter ist genau proportional zu Z_E8 "
      f"(Faktor {mp.nstr(u_quot[0].real, 8)}, Abweichung {u_fehler:.0e}; der Raum dieser Modulfunktionen ist eindimensional). "
      "Z_E8 + beta Z_U = (1 - 16 beta) Z_E8 hat für keine Projektion beta = +-1 ein einfaches Vakuum: Keine Hebung der "
      "2 pi-Drehung ist auf (E8)_1 modulverträglich, die Cartanströme (U = -1) und die Gitterströme (U = +1) liegen in "
      f"einem Multiplett. Relative Koeffizienten in der U-Bahn: {lam_u4}",
      abs(u_quot[0] + 16) < 1e-12 and u_fehler < 1e-12 and lam_u4["UR2"] == "1/2")
say(f"   {zeit()}")

# --- 2.3 nichtchirale Naht (Original (x) Partner): Bahnen, Charaktere, erzwungene Koeffizienten
KAND = {"1": ("1", "1"), "(-1)^(F+Fq)": ("R2", "R2"), "R Rq^-1": ("R", "R"), "U (x) 1": ("U", "U"),
        "G_tot = U (-1)^(F+Fq)": ("UR2", "UR2"), "U R Rq^-1": ("UR", "UR"), "G = U R^2 (einseitig)": ("UR2", "U"),
        "U R Rq": ("UR", "UR3"), "R^2 einseitig": ("R2", "1"), "U R einseitig": ("UR", "U")}
BB = {"B0": B0, "BU": BU}


def naht(pl, pr, N, lok=LOK, bb=BB):
    nl, sl, phl = lok[pl]
    nr, sr, phr = lok[pr]
    fak = [(bb[nl], sl, N), (bb[nr], sr, -N)]
    states, idx, kn = gemeinsam(fak)
    return fak, states, idx, kn, (N*phl - N*phr) % 1


tabelle = {}
for N in (1, 2, 4):
    tab = {}
    for name, (pl, pr) in KAND.items():
        fak, states, idx, kn, ph0 = naht(pl, pr, N)
        ch = charaktere(len(states), kn)
        triv = [V for c, V in ch if c["T"] == 0 and c["S"] == 0]
        enth = {}
        for b, (ql, qr) in KAND.items():
            ml, il, pl_ = LOK[ql]
            mr, ir, pr_ = LOK[qr]
            if ml == LOK[pl][0] and mr == LOK[pr][0] and (il, ir) in idx and triv:
                enth[b] = str((triv[0][idx[(il, ir)]] - (N*pl_ - N*pr_) - (triv[0][0] - ph0)) % 1)
            elif ml == LOK[pl][0] and mr == LOK[pr][0] and (il, ir) in idx:
                enth[b] = "in Bahn, anomal"
        tab[name] = {"Bahn": len(states), "Charaktere": chistr(ch), "enthält (rel. Koeffizient)": enth}
    tabelle[N] = tab
ZAHLEN["nichtchirale_naht"] = tabelle
t1 = tabelle[1]
check("[E] Nichtchirale Naht, Bahn der Spur ohne Einschub (6 Elemente, Charakter trivial): Sie enthält R R-quer^-1 und "
      f"(-1)^(F+F-quer) mit festem Koeffizienten {t1['1']['enthält (rel. Koeffizient)']} (Phase/2 pi). Die Modulinvarianz "
      "erzwingt also Niveaugleichheit R = R-quer und die Gesamtparität (-1)^(F+F-quer) = +1, J-kovariant; das ist die native "
      "GSO-Projektion der Naht. Keine Hebung der 2 pi-Drehung (U (x) 1, G, G_tot, U R R-quer) liegt in dieser Bahn, für N = 1, "
      "2, 4: Ihr Koeffizient ist frei, die Modulinvarianz erzwingt keine von ihnen",
      all(tabelle[N]["1"]["Bahn"] == 6 and tabelle[N]["1"]["Charaktere"] == ["T 0, S 0"]
          and tabelle[N]["1"]["enthält (rel. Koeffizient)"] == {"1": "0", "(-1)^(F+Fq)": "0", "R Rq^-1": "0"}
          for N in (1, 2, 4)))
check("[E] Die Hebungen ordnen sich in zwei getrennte Bahnen (je 24 Elemente, Charakter trivial und (1/2, 1/2)): "
      f"Bahn(U (x) 1) enthält {sorted(t1['U (x) 1']['enthält (rel. Koeffizient)'])} mit Koeffizient 0, also G_tot = U auf allen "
      "niveaugleichen Zuständen (dort ist F + F-quer gerade): Die J-kovariante Fassung der Regel ist die unverdrehte Hebung. "
      f"Bahn(G) enthält {sorted(t1['G = U R^2 (einseitig)']['enthält (rel. Koeffizient)'])}, getrennt von Bahn(U). Das "
      f"einseitige R^2 allein ist anomal für N = 1 ({t1['R^2 einseitig']['Charaktere']}), für N = 2 nur mit Charakter "
      f"{tabelle[2]['R^2 einseitig']['Charaktere']}, frei erst bei N = 4 ({tabelle[4]['R^2 einseitig']['Charaktere']})",
      sorted(t1["U (x) 1"]["enthält (rel. Koeffizient)"]) == ["G_tot = U (-1)^(F+Fq)", "U (x) 1", "U R Rq^-1"]
      and set(t1["U (x) 1"]["enthält (rel. Koeffizient)"].values()) == {"0"}
      and sorted(t1["G = U R^2 (einseitig)"]["enthält (rel. Koeffizient)"]) == ["G = U R^2 (einseitig)", "U R Rq"]
      and t1["G = U R^2 (einseitig)"]["Charaktere"] == ["T 0, S 0", "T 1/2, S 1/2"]
      and t1["R^2 einseitig"]["Charaktere"] == [] and tabelle[2]["R^2 einseitig"]["Charaktere"] == ["T 1/2, S 1/2"]
      and tabelle[4]["R^2 einseitig"]["Charaktere"] == ["T 0, S 0"])

# Gegenprobe: gemeinsame Bahn gegen direkte Bahn (N = 1)
direkt = {}
for name in ("1", "U (x) 1", "G = U R^2 (einseitig)"):
    pl, pr = KAND[name]
    fl, fr = h_orig(*TYP[pl]), h_orig(*TYP[pr])
    direkt[name] = bahn(lambda tau, fl=fl, fr=fr: fl(tau)*mp.conj(fr(tau)))["n"]
ZAHLEN["gegenprobe_direkte_bahn"] = direkt
check(f"[E] Gegenprobe: direkte Bahnbestimmung der nichtchiralen Spuren gibt dieselben Größen {direkt} wie die "
      "Zusammensetzung aus den holomorphen Faktorbahnen",
      direkt == {"1": 6, "U (x) 1": 24, "G = U R^2 (einseitig)": 24})
say(f"   {zeit()}")

# --- 2.4 Positivität der vollständigen Zustandssumme (spinaufgelöste Fouriermoden)
mp.mp.dps = 20


def summe(pl, pr, N, lok=LOK, bb=BB):
    """Bahnsumme mit trivialem Charakter, normiert so, dass die exakte Spur Tr(X q^L_0 qq^Lq_0) den Koeffizienten 1 hat."""
    fak, states, idx, kn, ph0 = naht(pl, pr, N, lok, bb)
    ch = [V for c, V in charaktere(len(states), kn) if c["T"] == 0 and c["S"] == 0]
    if not ch:
        return None
    V = ch[0]

    def Z(tau):
        cache, tot = {}, mp.mpc(0)
        for st, j in idx.items():
            val = mp.mpc(1)
            for (B, _, e), s in zip(fak, st):
                if (id(B), s) not in cache:
                    cache[(id(B), s)] = B["F"](act(B["reps"][s], tau))
                x = cache[(id(B), s)]
                val *= x**e if e > 0 else mp.conj(x)**(-e)
            tot += mp.exp(2*I*PI*R_(V[j]))*val
        return tot*mp.exp(2*I*PI*R_(ph0))
    return Z


YS = (mp.mpf('0.35'), mp.mpf('0.6'), mp.mpf('1.0'))
K = 24


def moden(Z):
    out = {}
    for Y in YS:
        v = [Z(mp.mpc(mp.mpf(k)/K, Y)) for k in range(K)]
        out[Y] = {s: mp.fsum(x*mp.exp(-2*I*PI*s*k/K) for k, x in enumerate(v))/K for s in range(-5, 6)}
    return out


def minimum(m0, m1, beta):
    """kleinster normierter Spinmode von Z_1 + beta Z_X (negativ: keine Hilbertraumspur)."""
    return float(min(((m0[Y][s] + beta*m1[Y][s]).real/abs(m0[Y][0]).real) for Y in YS for s in m0[Y]))


positiv = {}
for N in (1, 2, 3, 4):
    m0 = moden(summe("1", "1", N))
    zeile = {"Z_1 allein": minimum(m0, m0, 0)}
    for name, (pl, pr) in (("U (x) 1", ("U", "U")), ("G = U R^2", ("UR2", "U"))):
        m1 = moden(summe(pl, pr, N))
        zeile[name] = {"beta = +1": minimum(m0, m1, 1), "beta = -1": minimum(m0, m1, -1)}
    positiv[N] = zeile
    say(f"      N = {N}: min Z_1 {zeile['Z_1 allein']:.1e}; U: +1 {zeile['U (x) 1']['beta = +1']:.2e}, "
        f"-1 {zeile['U (x) 1']['beta = -1']:.2e}; G: +1 {zeile['G = U R^2']['beta = +1']:.2e}, "
        f"-1 {zeile['G = U R^2']['beta = -1']:.2e}  {zeit()}")
mg = moden(summe("geo", "geo", 1, {"geo": ("BG", 0, Fr(0)), "1": LOK["1"]}, {"BG": BGEO, "B0": B0}))
m01 = moden(summe("1", "1", 1))
positiv["geometrisch N = 1"] = {"beta = +1": minimum(m01, mg, 1), "beta = -1": minimum(m01, mg, -1)}
ZAHLEN["positivitaet (min normierter Spinmode, Y = 0.35, 0.6, 1.0, |s| <= 5)"] = positiv
EPS = -1e-8
ok = lambda x: x > EPS
check("[E, entscheidend] Positivität der vollständigen Zustandssumme Z_1 + beta Z_X (Untergrund: Bahnsumme ohne Einschub; "
      "beta = +1 heißt Projektion X = +1, beta = -1 X = -1; geprüft an allen Spinmoden |s| <= 5 bei Y = 0,35, 0,6, 1,0): "
      "Die unverdrehte Hebung U (x) 1 (gleich G_tot auf niveaugleichen Zuständen) ist für N = 1..4 nur mit U = +1 eine "
      "Hilbertraumspur, also ohne halbzahlige Helizitäten. Die Uhr-verdrehte Hebung G = U R^2 mit G = +1 (P_GSO) gibt bei "
      "N = 1 negative Vielfachheiten in den langen Nähten; konsistent ist dort nur G = -1 (Spin gegen Statistik vertauscht). "
      "Bei N = 2 (diagonal, G_1 G_2 = +1) ist G = +1 konsistent, bei N = 3 und 4 kein Vorzeichen. Die geometrische halbe "
      "Drehung ist bei N = 1 mit +1 konsistent",
      all(ok(positiv[N]["Z_1 allein"]) and ok(positiv[N]["U (x) 1"]["beta = +1"])
          and not ok(positiv[N]["U (x) 1"]["beta = -1"]) for N in (1, 2, 3, 4))
      and not ok(positiv[1]["G = U R^2"]["beta = +1"]) and ok(positiv[1]["G = U R^2"]["beta = -1"])
      and ok(positiv[2]["G = U R^2"]["beta = +1"]) and not ok(positiv[2]["G = U R^2"]["beta = -1"])
      and all(not ok(positiv[N]["G = U R^2"][b]) for N in (3, 4) for b in ("beta = +1", "beta = -1"))
      and ok(positiv["geometrisch N = 1"]["beta = +1"]))

# Robustheit: U nur auf einer radialen Leiter (Oszillatoren des Quarters a)
TYP1 = {"1": (0, 0, 0, 0), "U": (0, 0, 1, 0), "UR2": (HA, HA, 1, 0)}
B0_1 = bahn(h_orig(*TYP1["1"]))
BU_1 = bahn(h_orig(*TYP1["U"]))
LOK1 = {p: (nm, *lokalisiere(B, h_orig(*par))) for p, par in TYP1.items()
        for nm, B in (("B0", B0_1), ("BU", BU_1)) if lokalisiere(B, h_orig(*par)) is not None}
BB1 = {"B0": B0_1, "BU": BU_1}
leiter = {}
for N in (1, 2):
    m0 = moden(summe("1", "1", N, LOK1, BB1))
    leiter[N] = {}
    for name, (pl, pr) in (("U (x) 1", ("U", "U")), ("G = U R^2", ("UR2", "U"))):
        m1 = moden(summe(pl, pr, N, LOK1, BB1))
        leiter[N][name] = {"beta = +1": minimum(m0, m1, 1), "beta = -1": minimum(m0, m1, -1)}
ZAHLEN["positivitaet_eine_leiter"] = {"Bahnen": [B0_1["n"], BU_1["n"]], "N": leiter}
check("[E] Robustheit: Wirkt U_2pi nur auf einer radialen Leiter (Oszillatoren eines Quarters je Seite, wie die eine "
      f"Standardform in z97; Bahnen {B0_1['n']} und {BU_1['n']}), ist das Muster gleich: N = 1 nur U = +1 oder G = -1, "
      "N = 2 U = +1 oder G = +1",
      B0_1["n"] == 6 and BU_1["n"] == 24 and ok(leiter[1]["U (x) 1"]["beta = +1"])
      and not ok(leiter[1]["U (x) 1"]["beta = -1"]) and not ok(leiter[1]["G = U R^2"]["beta = +1"])
      and ok(leiter[1]["G = U R^2"]["beta = -1"]) and ok(leiter[2]["G = U R^2"]["beta = +1"])
      and not ok(leiter[2]["G = U R^2"]["beta = -1"]))
befund("Teil 2: Die Modulinvarianz erzwingt auf der Naht nur Niveaugleichheit R = R-quer und die Gesamtparität "
       "(-1)^(F+F-quer) = +1. Keine Hebung der 2 pi-Drehung liegt in der Bahn der Spur ohne Einschub. Verlangt man die Summe "
       "über die 2 pi-Drehung, dann lassen Modulinvarianz und Positivität pro Naht nur U = +1 (keine Spinoren) oder G = -1 "
       "(vertauschte Statistik) zu. P_GSO pro Naht ist ausgeschlossen und nur diagonal über zwei Nähte verträglich; auf der "
       "c = 8-Naht ist die Uhr eingefroren und keine Hebung verträglich", True)
say(f"   {zeit()}")

# =================================================================== Teil 3
say("\nTeil 3: Bilanz und Urteil für T4 (P_GSO aus Modulinvarianz?)")
STATUS = {
    "Uhr und R^2 auf dem Original": ("E", "e^(2 pi i L_0) = R, e^(4 pi i L_0) = R^2: T und T^2 sind Modultransformationen; "
                                          "Bahn 6, T-Zyklen [1, 1, 4]"),
    "U_2pi auf dem Torus": ("E", "halbe Drehung der Oszillatoren, keine Modultransformation; Bahn 24, T-Zyklen [4]^6, "
                                 "S-Bild eta(tau) eta(tau/4)/eta(tau/2)^3: lange Nähte"),
    "Kalibrierung 16 -> E8": ("E", "GSO-Vorzeichen +1 erzwungen genau für 16 und 32 Majorana-Weyl"),
    "chirales Original": ("E", "anomaliefrei genau für N = 0 mod 4; N = 4: Uhrensumme = (E8)_1, R = -i eingefroren, G = -U"),
    "U-Summe auf (E8)_1": ("E (negativ)", "U-Bahnsumme = -16 Z_E8: keine Projektion mit einfachem Vakuum"),
    "native Projektion der Naht": ("E", "Niveaugleichheit R = R-quer und (-1)^(F+F-quer) = +1, Koeffizienten erzwungen, "
                                        "J-kovariant"),
    "Hebungen der 2 pi-Drehung": ("E (negativ als Herleitung)", "keine in der Bahn der Spur ohne Einschub: nie erzwungen"),
    "G_tot (J-kovariant)": ("E", "in Bahn(U); auf niveaugleichen Zuständen gleich U: die unverdrehte Hebung"),
    "Positivität pro Naht": ("E (negativ für P_GSO)", "U = +1 oder G = -1 zulässig, G = +1 (P_GSO) negativ in den langen "
                                                       "Nähten; robust bei einer Leiter"),
    "N-Abhängigkeit": ("E", "G = +1 nur diagonal bei N = 2; N = 3, 4 kein Vorzeichen; U = +1 für N = 1..4"),
    "P_GSO gesamt": ("nein (nicht nativ)", "nicht aus Modulinvarianz mit P_CAR; pro Naht im Widerspruch zu Modulinvarianz "
                                           "und Positivität; verträglich nur als diagonale Zweinahtregel"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
ZAHLEN["urteil_kurz"] = {"P_GSO": "nicht nativ; pro Naht durch Modulinvarianz + Positivität ausgeschlossen",
                         "native Projektion": "R = R-quer und (-1)^(F+F-quer) = +1",
                         "erzwungene Hebung bei Rotationssumme": "unverdreht (U = +1) oder vertauscht (G = -1)",
                         "Chiralität": "unberührt (nichtchirale Naht vektorartig; chirales Original nur als (E8)_1 mit "
                                       "eingefrorener Uhr)"}
urteil = ("T4 (Folgt P_GSO aus der Modulinvarianz der Zustandssumme über N Nähte?): nein. Auf dem Original ist die Uhr die "
          "2D-Drehung: e^(2 pi i L_0) = R und e^(4 pi i L_0) = R^2 = (-1)^F gelten exakt, die Spuren mit R und R^2 liegen in "
          "der Bahn der Spur ohne Einschub, wie (-1)^F = T im NS-Sektor des Vergleichsziels. Die 4D-2 pi-Drehung U_2pi = "
          "(-1)^(m - mquer) ist dagegen die halbe Drehung der Oszillatoren (tau -> tau + 1/2) und keine Modultransformation. "
          "Ihre Spur hat eine eigene Bahn mit 24 Elementen; die verdrehten Sektoren sind lange Nähte mit Moden in Z/4 (S-Bild "
          "sqrt2 eta(tau) eta(tau/4)/eta(tau/2)^3). Die Kalibrierung reproduziert 16 -> E8: Dort erzwingt die Bahnsumme das "
          "GSO-Vorzeichen für 16 und 32 Majorana-Weyl. Das chirale Original ist genau für N = 0 mod 4 Nähte anomaliefrei; bei "
          "N = 4 (c = 8) ist seine Uhrensumme exakt (E8)_1, die Uhr ist dort auf R = -i eingefroren, also G = -U, und die "
          "U-Bahnsumme ist -16 Z_E8, so dass keine Hebung der 2 pi-Drehung mit einfachem Vakuum verträglich ist. Auf der "
          "nichtchiralen Naht (Original mit Partner) erzwingt die Modulinvarianz genau die Niveaugleichheit R = R-quer und die "
          "Gesamtparität (-1)^(F+F-quer) = +1 mit festen Koeffizienten. Das ist die native, J-kovariante GSO-Projektion der "
          "Naht. Keine Hebung der 2 pi-Drehung liegt in dieser Bahn, für N = 1, 2, 4; ihr Koeffizient ist frei, erzwungen wird "
          "keine. G_tot = U_2pi (-1)^(F+F-quer) liegt in der Bahn von U (x) 1 und ist auf niveaugleichen Zuständen gleich U: "
          "Die J-kovariante Fassung aus z97 ist die unverdrehte Hebung. Die einseitige Hebung G = U R^2 hat eine eigene, "
          "anomaliefreie Bahn. Verlangt man die Summe über die 2 pi-Drehung, entscheidet die Positivität der vollständigen "
          "Zustandssumme. Pro Naht (N = 1) sind nur U = +1 (keine halbzahligen Helizitäten) oder G = -1 (Spin gegen R^2-"
          "Statistik vertauscht) Hilbertraumspuren; G = +1, also P_GSO, gibt negative Vielfachheiten in den langen Nähten. "
          "Das gilt auch, wenn U nur auf einer radialen Leiter wirkt. Diagonal über zwei Nähte (G_1 G_2 = +1) ist G = +1 "
          "verträglich, bei N = 3 und 4 kein Vorzeichen; U = +1 ist für N = 1..4 verträglich. P_CAR ändert die Spuren nicht, "
          "weil alle Einschübe gerade sind. Urteil: P_GSO ist nicht nativ und folgt nicht aus Modulinvarianz mit P_CAR. Als "
          "Regel pro Naht steht P_GSO sogar im Widerspruch zu Modulinvarianz und Positivität, sobald über die 2 pi-Drehung "
          "summiert wird; verträglich ist sie nur als diagonale Zweinahtregel, die nicht die Regel pro Zustand aus z97 ist. "
          "Folge für die Statistikachse der Szenario-Suche: Der Wert 'Uhr-verdrehte Statistik P_GSO' ist durch die "
          "Modulinvarianz nicht gestützt. Wer die Torussumme verlangt, bekommt pro Naht 'U = +1' (keine Spinoren) oder "
          "'G = -1'; die native Projektion ist R = R-quer mit (-1)^(F+F-quer) = +1. Folge für T4: Die Chiralität bleibt "
          "getrennt und offen; die nichtchirale Naht ist vektorartig, das chirale Original ist nur als (E8)_1 mit "
          "eingefrorener Uhr konsistent. Nichts ist befördert; 3+1-Aussagen bleiben Kinematik der Einpaarschicht.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z113_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
