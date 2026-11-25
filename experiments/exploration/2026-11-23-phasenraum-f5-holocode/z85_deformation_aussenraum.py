"""Z85: Gravitation als Deformation des Außenraums (Strategie S4) und die Gegenhypothese S4b 'Marken = Weylkrümmung'.

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (kein Schluss 'modular => Bisognano-Wichmann'),
GRAV.PARENT.GAP.NO_GO.01 (gleichmäßig gegappter Elternoperator => kein masseloser Spin-2-Pol; T7 braucht eine kritische Folge
mit lokalen Tensor-Eichbedingungen), DIMENSION.UPLIFT.FIREWALL.01 (alles hier ist Kinematik oder lineare Twistortheorie).
Bestand: z58 (Ladung x Oszillatoren, Kumulanten), z61 (Marken = Hauptspinoren von 1 - z^4, I = -1, J = 0), z65 (googly,
Grad +2 / -6, Kerr-Schild), z75 (<J J J> = A f^abc, kein Graviton aus T), z78/z79 (Araki-Woods, Twistorraum aus den
Quarteroszillatoren mit Tomita-Partner), z80 (zusammengesetzter Spin 2 ohne Pol), z81 (R = i^(Q_a - Q_b) innen).
Statuskürzel: E exakt, C bedingt, L Lesart, O offen, N numerisch.

Teil 1  Außenraum und Wörterbuch: Twistoroperatoren Z = (d_1, d_2, e_1^+, e_2^+) aus den nativen Oszillatoren; Euler-Operator
        Upsilon = -(N_a - N_b) - 2 = -K/(4 pi m) - 2 (Homogenität = Modulgrad); Tomita J: n -> -n - 4. Die Fockbasis ist
        Penroses Elementarzustandsbasis (Čech-Repräsentanten modulo Korändern), alle 16 Bilineare von u(2,2) exakt
        verflochten; Vakuum = [1/(z1 z2)] mit Zentrum i kappa^2 im Tubus. Tafel n = 2, 0, -1, -2, -3, -4, -6 mit K-Typen;
        Kodaira-Spencer: nur Hamiltonfunktionen vom Grad 2 geben Vektorfelder auf PT.
Teil 2  Native Quellen einer Deformation: (a) gaußscher Sektor, (b) Vertexoperatoren e^{iQ phi} (Verschiebungen des
        Außenraums), (c) nichtabelsche Ströme (Frenkel-Kac, so(4)_1 der zwei Quarter), (d) Ladungssektor. Helizitätsgehalt,
        Grad-2-Klasse, Universalität; zweite Ordnung gegen die Plebanski-Penrose-Klammer; GRAV.PARENT.GAP.NO_GO.01.
Teil 3  S4b: Weylinvarianten und Stabilisatoren der H6-Quartik; Bianchi IX nach Cartan (geprüft an rund, Euler-Kreisel,
        Halphen); Isometrie R => biaxial => Typ D/O; Tods Form (Omega-System, W- = 0 exakt): J = 0 nur auf isolierten Bahnen,
        J identisch null nur konform flach; Hitchins vier Punkte als Symmetrielesart; Gegenprobe S1.
Teil 4  Bilanz und Urteil für S4 und S4b.
"""
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, least_squares

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")
    return bool(cond)


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


rng = np.random.default_rng(85)
mp.mp.dps = 30
I2 = np.eye(2)


def fock4(n):
    a = sps.diags(np.sqrt(np.arange(1, n)), 1, format="csr")
    Id = sps.identity(n, format="csr")
    ops = []
    for k in range(4):
        mats = [Id]*4
        mats[k] = a
        op = mats[0]
        for Mx in mats[1:]:
            op = sps.kron(op, Mx, format="csr")
        ops.append(op)
    return ops


# =================================================================== Teil 1
say("Teil 1: Außenraum -- Twistoroperatoren, Homogenität und Penroses Elementarzustände")
NM = 8
DIM = NM**4
d1, d2, e1, e2 = fock4(NM)
DD, EE = [d1, d2], [e1, e2]
OCC = np.array(list(itertools.product(range(NM), repeat=4)))      # (N_d1, N_d2, N_e1, N_e2) je Basisvektor
LOW = np.where(OCC.sum(1) <= 4)[0]
DVAL = OCC[:, 0] + OCC[:, 1] - OCC[:, 2] - OCC[:, 3]               # D = N_d - N_e
DOP = sps.diags(DVAL.astype(float), format="csr")
OMEGA = np.zeros(DIM)
OMEGA[0] = 1.0
ID = sps.identity(DIM, format="csr")


def dag(A):
    return A.T.conj().tocsr()


def komm(A, B):
    return (A @ B - B @ A).tocsr()


def klein(A, tol=1e-10):
    A = A.tocsc()[:, LOW]
    return (abs(A).max() if A.nnz else 0.0) < tol


qt = 0.3                                                          # Prüfwert; nativ q = e^{-2 pi m}
st = math.sqrt(1 - qt**2)
AA = [((DD[i] + qt*dag(EE[i]))/st).tocsr() for i in range(2)]     # Stromoszillatoren des Originals
BB = [((EE[i] + qt*dag(DD[i]))/st).tocsr() for i in range(2)]     # Tomita-Partner
ok_ccr = all(klein(komm(AA[i], dag(AA[j])) - (i == j)*ID) and klein(komm(BB[i], dag(BB[j])) - (i == j)*ID)
             and klein(komm(AA[i], BB[j])) and klein(komm(AA[i], dag(BB[j]))) for i in range(2) for j in range(2))
ok_tfd = all(np.allclose(AA[i] @ OMEGA, qt*(dag(BB[i]) @ OMEGA)) and np.allclose(BB[i] @ OMEGA, qt*(dag(AA[i]) @ OMEGA))
             for i in range(2))
nbar = float(OMEGA @ (dag(AA[0]) @ (AA[0] @ OMEGA)))
ZT = [d1, d2, dag(e1), dag(e2)]                                   # Z^alpha
ZB = [dag(d1), dag(d2), -e1, -e2]                                 # Z-quer_alpha
ok_tw = all(klein(komm(ZT[al], ZB[be]) - (al == be)*ID) and klein(komm(ZT[al], ZT[be])) for al in range(4) for be in range(4))
check("a_i = (d_i + q e_i^+)/s und b_i = (e_i + q d_i^+)/s erfüllen die CCR und a_i Omega = q b_i^+ Omega (Purifikation, "
      f"<a^+ a> = {nbar:.6f} = q^2/(1-q^2)); die Twistoroperatoren Z = (d_1, d_2, e_1^+, e_2^+) aus span{{a_i, b_i^+}} (z79) "
      "mit Z-quer = (d_1^+, d_2^+, -e_1, -e_2) erfüllen [Z^alpha, Z-quer_beta] = delta, [Z, Z] = 0 (Fockraum 8^4, alle Zustände "
      "bis Stufe 4)", ok_ccr and ok_tfd and ok_tw and abs(nbar - qt**2/(1 - qt**2)) < 1e-12)

UPS = -sum((ZT[al] @ ZB[al] for al in range(4)), sps.csr_matrix((DIM, DIM)))
NaNb = sum((dag(AA[i]) @ AA[i] - dag(BB[i]) @ BB[i] for i in range(2)), sps.csr_matrix((DIM, DIM)))
ok_eu = klein(UPS + DOP + 2*ID) and klein(NaNb - DOP) and all(klein(komm(DOP, ZT[al]) + ZT[al]) for al in range(4))
kms = all(abs(np.linalg.norm((qt**DVAL)*(AA[i] @ OMEGA)) - np.linalg.norm(dag(AA[i]) @ OMEGA)) < 1e-12 for i in range(2))
perm = np.array([np.ravel_multi_index((o[2], o[3], o[0], o[1]), (NM,)*4) for o in OCC])
SW = sps.csr_matrix((np.ones(DIM), (perm, np.arange(DIM))), shape=(DIM, DIM))
ok_J = all(klein(SW @ AA[i] @ SW - BB[i]) for i in range(2)) and klein(SW @ DOP @ SW + DOP) and np.allclose(SW @ OMEGA, OMEGA)
check("Homogenität = Modulgrad: Upsilon = -Z^alpha Z-quer_alpha = -(N_d - N_e) - 2 = -(N_a - N_b) - 2 als Operatoridentität, "
      "[D, Z^alpha] = -Z^alpha; mit Delta^(1/2) = q^D (KMS geprüft) ist K = -2 log q (N_a - N_b), nativ q = e^{-2 pi m}: "
      "Upsilon = -K/(4 pi m) - 2. Die Tomita-Konjugation (a <-> b) gibt J D J = -D, also n -> -n - 4 (Twistorkonjugation)",
      ok_eu and kms and ok_J)


# --- 1.3 Fockbasis = Elementarzustände (exakt, sqrt-Darstellung: Koeffizient = Vorzeichen * sqrt(Bruch))
def fock_op(kind, idx, st_):
    key, sg, Q = st_
    k = list(key)
    n = k[idx]
    if kind == "ab":                                              # Vernichter auf Mode idx
        if n == 0:
            return None
        k[idx] = n - 1
        return tuple(k), sg, Q*n
    k[idx] = n + 1                                                # Erzeuger
    return tuple(k), sg, Q*(n + 1)


def Z_fock(al, st_):                                              # Z = (d1, d2, e1^+, e2^+)
    return fock_op("ab", al, st_) if al < 2 else fock_op("auf", al, st_)


def Zb_fock(be, st_):                                             # Z-quer = (d1^+, d2^+, -e1, -e2)
    if be < 2:
        return fock_op("auf", be, st_)
    r = fock_op("ab", be, st_)
    return None if r is None else (r[0], -r[1], r[2])


def cech(st_):                                                    # |n1 n2 m1 m2> -> sqrt(n1! n2!/(m1! m2!)) z3^m1 z4^m2/(z1^(n1+1) z2^(n2+1))
    if st_ is None:
        return None
    (n1, n2, m1, m2), sg, Q = st_
    return (n1, n2, m1, m2), sg, Q*Fraction(math.factorial(n1)*math.factorial(n2), math.factorial(m1)*math.factorial(m2))


def Z_cech(al, mono):                                             # Multiplikation mit z^alpha, Korand (kein Pol) -> 0
    if mono is None:
        return None
    k, sg, Q = mono
    k = list(k)
    if al < 2:
        if k[al] == 0:
            return None
        k[al] -= 1
    else:
        k[al] += 1
    return tuple(k), sg, Q


def Zb_cech(be, mono):                                            # -d/dz^beta
    if mono is None:
        return None
    k, sg, Q = mono
    k = list(k)
    if be < 2:
        c = k[be] + 1                                             # -d z^-(n+1) = (n+1) z^-(n+2)
        k[be] += 1
        return tuple(k), sg, Q*c*c
    if k[be] == 0:
        return None
    c = k[be]
    k[be] -= 1
    return tuple(k), -sg, Q*c*c


basis = [(tuple(o), 1, Fraction(1)) for o in itertools.product(range(5), repeat=4) if sum(o) <= 4]
ok_int, n_vgl = True, 0
for st_ in basis:
    for al in range(4):
        ok_int &= cech(Z_fock(al, st_)) == Z_cech(al, cech(st_))
        ok_int &= cech(Zb_fock(al, st_)) == Zb_cech(al, cech(st_))
        n_vgl += 2
        for be in range(4):
            r = Zb_fock(be, st_)
            lhs = cech(None if r is None else Z_fock(al, r))
            rhs_ = Z_cech(al, Zb_cech(be, cech(st_)))
            ok_int &= lhs == rhs_
            n_vgl += 1
ZAHLEN["verflechtung_vergleiche"] = n_vgl
check(f"Fockbasis = Penroses Elementarzustände: |n1 n2; m1 m2> <-> sqrt(n1! n2!/(m1! m2!)) z3^m1 z4^m2/(z1^(n1+1) z2^(n2+1)) "
      f"auf U1 n U2 modulo Korändern verflicht Z^alpha (Multiplikation), Z-quer_beta (-d/dz^beta) und alle 16 Bilineare "
      f"Z^alpha Z-quer_beta von u(2,2) exakt ({n_vgl} Vergleiche auf allen {len(basis)} Zuständen bis Stufe 4); d_i Omega = 0 "
      "ist die Korandbedingung (z_i/z_i = 1 hat keinen Pol in z_i)", ok_int)


def cech_falsch(st_):                                             # Negativkontrolle: Normierung sqrt(n! m!)
    if st_ is None:
        return None
    (n1, n2, m1, m2), sg, Q = st_
    return (n1, n2, m1, m2), sg, Q*math.factorial(n1)*math.factorial(n2)*math.factorial(m1)*math.factorial(m2)


ok_kontr = all(cech_falsch(Zb_fock(be, s_)) == Zb_cech(be, cech_falsch(s_)) and
               cech_falsch(Z_fock(be, s_)) == Z_cech(be, cech_falsch(s_)) for s_ in basis for be in range(4))
check("Negativkontrolle: mit der Normierung sqrt(n! m!) statt sqrt(n!/m!) scheitert die Verflechtung (der Test trennt)",
      not ok_kontr)

# --- 1.4 Koordinaten und Vakuumklasse
kap = math.sqrt((1 - qt)/(1 + qt))
PI_ = [((AA[i] + dag(BB[i]))/math.sqrt(2)).tocsr() for i in range(2)]
OM_ = [((AA[i] - dag(BB[i]))/math.sqrt(2)).tocsr() for i in range(2)]
ok_ko = all(klein(DD[i] - (kap*PI_[i] + OM_[i]/kap)/math.sqrt(2)) and klein(dag(EE[i]) - (kap*PI_[i] - OM_[i]/kap)/math.sqrt(2))
            for i in range(2))
Xs = sp.Matrix(2, 2, sp.symbols("x00 x01 x10 x11"))
ks = sp.Symbol("kappa", positive=True)
zeta = sp.Symbol("zeta")
Mz = (ks*sp.eye(2) + sp.I/ks*Xs)/sp.sqrt(2)                       # z_i = (kappa pi_i + omega_i/kappa)/sqrt2, omega = i x pi
l1, l2 = Mz[0, 0] + Mz[0, 1]*zeta, Mz[1, 0] + Mz[1, 1]*zeta
z0 = sp.solve(l1, zeta)[0]
phi = sp.simplify(1/(Mz[0, 1]*l2.subs(zeta, z0)))                 # Residuum von 1/(z1 z2) auf L_x
konst = sp.simplify(phi*(Xs - sp.I*ks**2*sp.eye(2)).det())
kap2_nativ = (1 - mp.exp(-2*mp.pi))/(1 + mp.exp(-2*mp.pi))
ZAHLEN["vakuumklasse"] = {"phi*det(x - i kappa^2)": str(konst), "kappa^2 nativ (m=1)": mp.nstr(kap2_nativ, 12),
                          "tanh(pi)": mp.nstr(mp.tanh(mp.pi), 12)}
check(f"Koordinaten: d_i = (kappa pi_i + omega_i/kappa)/sqrt2, e_i^+ = (kappa pi_i - omega_i/kappa)/sqrt2 (pi, omega aus z79); "
      f"die Vakuumklasse [1/(z1 z2)] (Homogenität -2) hat die Penrose-Transformierte phi = {konst}/det(x - i kappa^2 I): das "
      f"skalare Elementarzustandsfeld mit Zentrum i kappa^2 I im Vorwärtstubus (z79: Re M = kappa^2 I); nativ kappa^2 = "
      f"tanh(pi m) = {mp.nstr(kap2_nativ, 8)} für m = 1", ok_ko and konst.free_symbols <= {ks} and konst != 0
      and abs(kap2_nativ - mp.tanh(mp.pi)) < mp.mpf(10)**-25)

# --- 1.5 Wörterbuch Homogenität <-> Helizität <-> Focksektor
tafel = []
ok_tab = True
for D in (-4, -2, -1, 0, 1, 2, 4):
    n_h = -D - 2
    zaehl = []
    for k in range(6):
        Ntot = abs(D) + 2*k
        c = sum(1 for o in itertools.product(range(Ntot + 1), repeat=4)
                if sum(o) == Ntot and (o[0] + o[1]) - (o[2] + o[3]) == D)
        zaehl.append(c)
        ok_tab &= c == (k + 1)*(k + 1 + abs(D))
    tafel.append({"D = N_d - N_e": D, "Homogenität n": n_h, "Helizität h": D/2, "E_min": 1 + abs(D)/2, "K-Typen k=0..5": zaehl})
rolle = {2: "Penrose-Deformation von PT: H1(PT, O(2)), Psi_ABCD (z65: Helizität -2), Kodaira-Spencer-Klasse",
         0: "Helizität -1 (ASD-Maxwell bzw. Wards Bündeldeformation), NICHT der Skalar",
         -1: "Helizität -1/2 (Weyl, eine Chiralität)",
         -2: "Skalar (konformer Skalar), Sektor des Vakuums Omega",
         -3: "Helizität +1/2 (z61: Grad -3)",
         -4: "Maxwell phi_A'B' (Helizität +1)",
         -6: "Psi~_A'B'C'D' (z61, Helizität +2); googly: Deformation von PT* = Tomita-Partner (z65)"}
for zeile in tafel:
    zeile["Rolle"] = rolle[zeile["Homogenität n"]]
    say(f"   n = {zeile['Homogenität n']:+d}  D = {zeile['D = N_d - N_e']:+d}  h = {zeile['Helizität h']:+.1f}  "
        f"E_min = {zeile['E_min']}  K-Typen {zeile['K-Typen k=0..5']}  -- {zeile['Rolle']}")
ZAHLEN["woerterbuch"] = tafel
check("Wörterbuch n = -D - 2 = -2h - 2: die Elementarzustände vom Grad n bilden genau den Fock-Sektor D = N_d - N_e = -n - 2 "
      "mit den 4D-K-Typ-Dimensionen (k+1)(k+1+|D|) und niedrigster Energie 1 + |D|/2 (k = 0..5, alle sieben Grade); n = 2 "
      "(Graviton, Kodaira-Spencer) ist D = -4, n = -6 ist D = +4 (Tomita-Bild), der Skalar ist n = -2 (D = 0), n = 0 ist "
      "Helizität -1", ok_tab)

# --- 1.6 Kodaira-Spencer: Hamiltonsches Vektorfeld V_h = eps^AB d_B h d_A auf den omega-Fasern
w0, w1, p0, p1 = sp.symbols("w0 w1 p0 p1")
KO = [w0, w1, p0, p1]


def V_h(h):
    return [sp.diff(h, w1), -sp.diff(h, w0), 0, 0]               # Komponenten auf (w0, w1, p0, p1), eps^01 = 1


def euler_gewicht(V):
    """[Upsilon, V] = c V, Upsilon = sum Z d_Z: liefert c (oder None)."""
    br = [sp.simplify(sum(KO[j]*sp.diff(V[i], KO[j]) for j in range(4)) - V[i]) for i in range(4)]
    cs = {sp.simplify(b/v) for b, v in zip(br, V) if v != 0}
    return cs.pop() if len(cs) == 1 else None


ok_ks = True
gew = {}
for n_h in (-6, -2, 0, 2, 4):
    for _ in range(2):
        ex = [int(t) for t in rng.integers(-3, 4, size=3)]
        e4 = n_h - sum(ex)
        h = w0**ex[0]*w1**ex[1]*p0**ex[2]*p1**e4*(1 + w0*p1/(w1*p0))   # Homogenität n_h, rational (Pole wie Elementarzustände)
        c = euler_gewicht(V_h(h))
        gew.setdefault(n_h, set()).add(c)
        ok_ks &= c == n_h - 2
y_ = sp.symbols("y0:4")
h_tr = w0*(y_[2]*p0 + y_[3]*p1) - w1*(y_[0]*p0 + y_[1]*p1)       # Grad 2, polynomial
Vt = V_h(h_tr)
ok_tr = all(sp.diff(Vt[i], KO[j]) == 0 for i in range(2) for j in range(2)) and Vt[2] == 0 and Vt[3] == 0
ZAHLEN["euler_gewicht_V_h"] = {str(k): [str(c) for c in v] for k, v in gew.items()}
check("Kodaira-Spencer: V_h = eps^AB (dh/d omega^B) d/d omega^A hat das Euler-Gewicht n_h - 2 (rationale h mit Polen, "
      "n_h = -6, -2, 0, 2, 4): nur Hamiltonfunktionen vom Grad 2 geben Vektorfelder auf PT; ein polynomiales h vom Grad 2 gibt "
      "omega -> omega + y pi (Translation, global, also trivial in H1)", ok_ks and ok_tr)
befund("E: Der Außenraum von z79 ist als Fockraum die Elementarzustandsdarstellung der Twistortheorie: Fockvektoren = "
       "Čech-Klassen, Vernichtung = Korand, Homogenität = -(Modulgenerator)/(4 pi m) - 2, Tomita J = Konjugation n -> -n - 4. "
       "Eine Kodaira-Spencer-Klasse erster Ordnung (Penroses selbstduales Graviton, Grad 2) ist genau ein Vektor im Sektor "
       "D = N_a - N_b = -4, also Modulfrequenz K = -16 pi m; das googly-Gegenstück (Grad -6) ist sein Tomita-Bild D = +4", True)

# --- 1.7 Casimir je Helizitätssektor (Abgleich mit S2/z83: Helizität h = (m - m-quer)/2 als Einteilchensektor der Leiter)
EOP = [[(ZT[al] @ ZB[be]).tocsr() for be in range(4)] for al in range(4)]
NOP = sum((EOP[al][al] for al in range(4)), sps.csr_matrix((DIM, DIM)))
CU = sum((EOP[al][be] @ EOP[be][al] for al in range(4) for be in range(4)), sps.csr_matrix((DIM, DIM)))
CSU = (CU - 0.25*(NOP @ NOP)).tocsc()
SAFE = np.where((OCC.max(1) <= 5) & (OCC.sum(1) <= 6))[0]          # zwei Bilineare heben jede Mode höchstens um 2
res_cas = abs((CSU - sps.diags(3*(DVAL**2/4 - 1))).tocsc()[:, SAFE]).max()
sek = sorted({int(DVAL[i]) for i in SAFE})
ZAHLEN["casimir_su22"] = {"Formel": "C2 = sum E^a_b E^b_a - (tr E)^2/4 = 3(h^2 - 1), h = D/2", "Residuum": float(res_cas),
                          "geprüfte D-Sektoren": sek, "z83 (S2)": "3(h^2 - 1)/8 in der dortigen Erzeugernormierung"}
check(f"su(2,2)-Casimir C2 = sum E^a_b E^b_a - (tr E)^2/4 = 3(h^2 - 1) mit h = D/2 auf jedem Sektor D = {sek[0]}..{sek[-1]} "
      f"(alle Zustände mit Besetzung <= 5 und Stufe <= 6, Residuum {res_cas:.1e}): jeder D-Sektor ist ein irreduzibler "
      "masseloser Einteilchensektor, Helizität +-2 eingeschlossen; das ist dieselbe Leiter wie in S2/z83 (h = (m - m-quer)/2, "
      "Casimir 3(h^2 - 1)/8 in der dortigen Normierung, niedrigster K-Typ bei E = |h| + 1)", res_cas < 1e-9 and min(sek) <= -4
      and max(sek) >= 4)
befund("E/C: Die Helizität-2-Sektoren existieren kinematisch als Einteilchensektoren (z75, z83, hier): D = -4 trägt die "
       "Kodaira-Spencer-Klassen vom Grad 2, D = +4 die googly-Klassen vom Grad -6. Gerealisiert sind sie mit mindestens vier "
       "Oszillatorquanten, wie in jeder Twistorquantisierung. Die Frage von S4 ist daher nicht, ob der Sektor existiert, "
       "sondern ob eine native Quelle ihn mit universeller Kopplung bevölkert und ob die zweite Ordnung Einstein ist", True)

# =================================================================== Teil 2
say("\nTeil 2: Native Quellen einer Deformation des Außenraums")


def d_gewichte(v):
    w = {}
    for Dv in np.unique(DVAL):
        p = float(np.sum(np.abs(v[DVAL == Dv])**2))
        if p > 1e-28:
            w[int(Dv)] = p
    return w


# --- (a) gaußscher Sektor
LIN = [d1, d2, dag(d1), dag(d2), e1, e2, dag(e1), dag(e2)]
ok_a, maxD = True, 0
for i in range(8):
    w = d_gewichte(LIN[i] @ OMEGA)
    maxD = max([maxD] + [abs(k) for k in w])
    for j in range(i, 8):
        w = d_gewichte(LIN[i] @ (LIN[j] @ OMEGA))
        maxD = max([maxD] + [abs(k) for k in w])
        ok_a &= all(abs(k) <= 2 for k in w)
dDelta = {"a": -1, "a+": 1, "b": 1, "b+": -1}
OPS8 = {("a", i): AA[i] for i in range(2)} | {("a+", i): dag(AA[i]) for i in range(2)} | \
       {("b", i): BB[i] for i in range(2)} | {("b+", i): dag(BB[i]) for i in range(2)}
ok_zp, n_nz = True, 0
for (k1, X), (k2, Y) in itertools.product(OPS8.items(), repeat=2):
    val = OMEGA @ (X @ (Y @ OMEGA))
    if abs(val) > 1e-12:
        n_nz += 1
        ok_zp &= dDelta[k1[0]] + dDelta[k2[0]] == 0
cc = [(DD[i] + dag(EE[i])).tocsr() for i in range(2)]             # c = d + e^+ = kappa (a + b^+) (z78)
SIGM = [I2, np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1., -1.])]
PMU = [0.5*sum((s_[i, j]*(dag(cc[i]) @ cc[j]) for i in range(2) for j in range(2)), sps.csr_matrix((DIM, DIM)))
       for s_ in SIGM]
P2 = (PMU[0] @ PMU[0] - sum((PMU[k] @ PMU[k] for k in range(1, 4)), sps.csr_matrix((DIM, DIM)))).tocsc()
ok_P = all(klein(komm(DOP, Pm)) for Pm in PMU) and abs(P2[:, SAFE]).max() < 1e-9
check(f"(a) gaußsch: alle 8 linearen und 36 quadratischen Oszillatormonome auf Omega haben nur D in {{-2..2}} (max |D| = "
      f"{maxD}), der Sektor |D| = 4 bleibt exakt leer; die nativen Gaußdaten sind D-neutral: alle {n_nz} nichtverschwindenden "
      "Zweipunktfunktionen <X Y> (X, Y aus a, a^+, b, b^+) haben Delta D = 0, und P_mu = (1/2) c^+ sigma_mu c, Lorentz, "
      "Dilatation und K vertauschen mit D (P^2 = 0 exakt)", ok_a and ok_zp and n_nz > 0 and ok_P)
befund("E (a): Der gaußsche Sektor allein deformiert nichts. Seine nativen Daten sind modular invariant (D = 0) und erzeugen "
       "nur u(2,2), also projektive Abbildungen von PT (Translationen, Lorentz, Dilatation): globale, triviale Klassen. Selbst "
       "die nicht nativen Quetscher erreichen in erster Ordnung nur |D| = 2 (Helizität +-1), nie den Grad 2", True)

# --- (b) Vertexoperatoren: Oszillatorfaktor von :e^{iQ phi}: an Mode m = Verschiebung mit alpha_i = Q_i e^{i m theta}/sqrt(m)
NV = 10
DIMV = NV**4
dv1, dv2, ev1, ev2 = fock4(NV)
OCCV = np.array(list(itertools.product(range(NV), repeat=4)))
AV = [((d_ + qt*dag(e_))/st).tocsr() for d_, e_ in ((dv1, ev1), (dv2, ev2))]
OMV = np.zeros(DIMV, dtype=complex)
OMV[0] = 1


def W_omega(alpha):
    G = sum((alpha[i]*dag(AV[i]) - np.conj(alpha[i])*AV[i] for i in range(2)), sps.csr_matrix((DIMV, DIMV)))
    return spla.expm_multiply(G.tocsc(), OMV), G


def koh_analytisch(alpha, q):
    s_ = math.sqrt(1 - q**2)
    be = np.array(alpha)/s_
    ga = -q*np.conj(np.array(alpha))/s_
    fac = np.exp(-(np.sum(np.abs(be)**2) + np.sum(np.abs(ga)**2))/2)
    v = np.array([be[0]**o[0]*be[1]**o[1]*ga[0]**o[2]*ga[1]**o[3]/math.sqrt(math.factorial(o[0])*math.factorial(o[1])
                                                                           * math.factorial(o[2])*math.factorial(o[3]))
                  for o in OCCV])
    return fac*v


def a_niedrig(alpha, q, D):
    s_ = math.sqrt(1 - q**2)
    n2 = float(np.sum(np.abs(np.array(alpha))**2))
    return math.exp(-(1 + q**2)*n2/(2*s_**2))*n2**2/(s_**4*math.sqrt(24))*(q**4 if D < 0 else 1.0)


ok_b1, ok_b2, fid = True, True, 1.0
for _ in range(3):
    al = (rng.normal(size=2) + 1j*rng.normal(size=2))*0.35
    v, G = W_omega(al)
    va = koh_analytisch(al, qt)
    fid = min(fid, abs(np.vdot(va, v))/np.linalg.norm(va)/np.linalg.norm(v))
    g1 = G @ OMV
    ok_b1 &= all(abs((OCCV[i, 0] + OCCV[i, 1]) - (OCCV[i, 2] + OCCV[i, 3])) == 1 for i in np.nonzero(np.abs(g1) > 1e-14)[0])
    for D, sel in ((4, (OCCV[:, 0] + OCCV[:, 1] == 4) & (OCCV[:, 2] + OCCV[:, 3] == 0)),
                   (-4, (OCCV[:, 0] + OCCV[:, 1] == 0) & (OCCV[:, 2] + OCCV[:, 3] == 4))):
        num_ = np.linalg.norm(v[sel])
        ok_b2 &= abs(num_/a_niedrig(al, qt, D) - 1) < 1e-8
check(f"(b) der Oszillatorfaktor eines Vertexoperators wirkt als W(alpha) = D_d(alpha/s) (x) D_e(-q alpha-quer/s) auf Omega "
      f"(Fockraum 10^4, q = 0,3, drei zufällige alpha, Treue 1 - {1 - fid:.1e}); in erster Ordnung (Erzeuger auf Omega) nur "
      "D = +-1 (Helizität +-1/2); die niedrigsten K-Typen von D = +4 und D = -4 haben exakt die Amplituden "
      "e^{-(1+q^2)|alpha|^2/(2 s^2)} |alpha|^4/(s^4 sqrt24) bzw. q^4 mal das: der Grad-2-Sektor wird erst in vierter Ordnung "
      "erreicht", fid > 1 - 1e-9 and ok_b1 and ok_b2)
# nativ: q = e^{-2 pi m}; Universalität über Ladungssektoren (Quarter alpha = 1/4, 3/4)
mp_q = {m: mp.exp(-2*mp.pi*m) for m in (1, 2, 3)}
ratio = {m: mp.nstr(mp_q[m]**4, 6) for m in mp_q}
ok_rat = all(abs(mp_q[m]**4/mp.exp(-8*mp.pi*m) - 1) < mp.mpf(10)**-25 for m in mp_q)
qn = float(mp_q[1])
ladungen = [(1, 0), (0, 1), (1, 1), (1, -1), (2, 0), (2, 1), (3, 0)]
univ = {}
for Qv in ladungen:
    al = np.array(Qv, dtype=float)                                # m = 1, theta = 0
    h_w = 0.5*float(al @ al)
    E_tw = sum(Qi**2/2 + (ai - 0.5)*Qi for Qi, ai in zip(Qv, (0.25, 0.75)))
    univ[str(Qv)] = {"a(+4)": a_niedrig(al, qn, 4), "a(+4)/h": a_niedrig(al, qn, 4)/h_w, "h": h_w, "E(Q)": E_tw}
quot = [u["a(+4)/h"] for u in univ.values()]
v_pol = {}
for Qv in ((1, 0), (0, 1), (1, 1), (1, -1)):
    v, _ = W_omega(0.3*np.array(Qv, dtype=complex)/np.linalg.norm(Qv))
    sel = (OCCV[:, 0] + OCCV[:, 1] == 4) & (OCCV[:, 2] + OCCV[:, 3] == 0)
    v_pol[Qv] = v[sel]/np.linalg.norm(v[sel])
ov = {"(1,0)|(0,1)": abs(np.vdot(v_pol[(1, 0)], v_pol[(0, 1)])), "(1,1)|(1,-1)": abs(np.vdot(v_pol[(1, 1)], v_pol[(1, -1)])),
      "(1,0)|(1,1)": abs(np.vdot(v_pol[(1, 0)], v_pol[(1, 1)]))}
mittel = sum(W_omega(0.3*np.exp(2j*np.pi*k/8)*np.array([1, 0]))[0][(OCCV[:, 0] + OCCV[:, 1] == 4) & (OCCV[:, 2] + OCCV[:, 3] == 0)]
             for k in range(8))/8
ZAHLEN["vertex_grad2"] = {"a(-4)/a(+4) = q^4 = e^{-8 pi m}": ratio, "Ladungssektoren m=1": univ, "Polarisation_Überlapp": ov,
                          "theta_mittel_norm": float(np.linalg.norm(mittel))}
check(f"(b) nativ ist a(D = -4)/a(D = +4) = q^4 = e^(-8 pi m) ({ratio[1]} für m = 1): Penroses Grad-2-Klasse ist gegen das "
      f"googly-Bild KMS-unterdrückt. Nicht universell: a(+4)/h schwankt über die Ladungssektoren zwischen {min(quot):.3g} und "
      f"{max(quot):.3g} (h = |Q|^2/2), die Polarisation ist der Ladungsvektor (Überlapp (1,0)|(0,1) = {ov['(1,0)|(0,1)']:.1e}, "
      f"(1,1)|(1,-1) = {ov['(1,1)|(1,-1)']:.1e}, (1,0)|(1,1) = {ov['(1,0)|(1,1)']:.4f} = 1/4), und das Mittel über die Nahtlage "
      f"theta (Phase e^(4 i m theta)) verschwindet ({np.linalg.norm(mittel):.1e})", ok_rat and max(quot)/min(quot) > 2
      and ov["(1,0)|(0,1)"] < 1e-10 and ov["(1,1)|(1,-1)"] < 1e-10 and abs(ov["(1,0)|(1,1)"] - 0.25) < 1e-8
      and np.linalg.norm(mittel) < 1e-12)
befund("E (b): Vertexoperatoren (Ladung x Oszillator) sind auf dem Außenraum Twistorverschiebungen. Sie erreichen den Grad-2-"
       "Sektor nur in vierter Ordnung, als Produkt von vier Helizität-1/2-Quanten, mit Amplitude ~ h^2 e^{-c h} statt ~ Energie, "
       "mit dem Ladungsvektor als Polarisationsspinor und mit dem Tomita-Gewicht q^4: Eine Grad-2-Klasse entsteht, aber nicht "
       "universell", True)

# --- (c) nichtabelsche Ströme: Frenkel-Kac, so(4)_1 = su(2)_L + su(2)_R auf den zwei Quarterfermionen
wurzeln = {"E_L": (1, -1), "F_L": (-1, 1), "E_R": (1, 1), "F_R": (-1, -1)}
v_w = {}
for nm_, r_ in wurzeln.items():
    v, _ = W_omega(0.25*np.array(r_, dtype=complex))
    sel = (OCCV[:, 0] + OCCV[:, 1] == 4) & (OCCV[:, 2] + OCCV[:, 3] == 0)
    v_w[nm_] = v[sel]
betr = {k: float(np.linalg.norm(v_)) for k, v_ in v_w.items()}
ov_LR = abs(np.vdot(v_w["E_L"], v_w["E_R"]))/betr["E_L"]/betr["E_R"]
ov_LL = abs(np.vdot(v_w["E_L"], v_w["F_L"]))/betr["E_L"]/betr["F_L"]
x_, y_ = np.array([0.3, 0.0]), np.array([0.0, 0.3])               # reelle, orthogonale Verschiebungen
Gx, Gy = W_omega(x_)[1], W_omega(y_)[1]
kom = spla.expm_multiply(Gx.tocsc(), spla.expm_multiply(Gy.tocsc(), OMV)) - \
    spla.expm_multiply(Gy.tocsc(), spla.expm_multiply(Gx.tocsc(), OMV))
al_c, be_c = np.array([0.2 + 0.1j, -0.15j]), np.array([0.1, 0.25 - 0.05j])
lhs = spla.expm_multiply(W_omega(al_c)[1].tocsc(), W_omega(be_c)[0])
rhs_w = np.exp(1j*np.imag(al_c @ np.conj(be_c)))*W_omega(al_c + be_c)[0]
# Plebanski-Penrose-Klammer zweier Grad-2-Klassen h_a = (a1 z3 + a2 z4)^4/(z1 z2) auf den omega-Fasern
zc = sp.symbols("z1:5")
ka_ = sp.Symbol("kappa", positive=True)
al_s, be_s = sp.symbols("al1 al2"), sp.symbols("be1 be2")
Aform = al_s[0]*zc[2] + al_s[1]*zc[3]
Bform = be_s[0]*zc[2] + be_s[1]*zc[3]
h_al, h_be = Aform**4/(zc[0]*zc[1]), Bform**4/(zc[0]*zc[1])


def d_om(f, j):                                                   # d/d omega_j bei festem pi
    return (sp.diff(f, zc[j]) - sp.diff(f, zc[j + 2]))/(sp.sqrt(2)*ka_)


klammer = sp.simplify(d_om(h_al, 0)*d_om(h_be, 1) - d_om(h_al, 1)*d_om(h_be, 0))
ab = al_s[0]*be_s[1] - al_s[1]*be_s[0]
soll = 2/ka_**2*ab*Aform**3*Bform**3*(zc[2]/(zc[0]**3*zc[1]**2) + zc[3]/(zc[0]**2*zc[1]**3) + 4/(zc[0]**2*zc[1]**2))
eul = sp.simplify(sum(z_*sp.diff(klammer, z_) for z_ in zc) - 2*klammer)
pole = all(sp.degree(sp.denom(sp.together(t)), zc[0]) >= 1 and sp.degree(sp.denom(sp.together(t)), zc[1]) >= 1
           for t in sp.Add.make_args(sp.expand(soll)))
spez = klammer.subs({al_s[0]: 1, al_s[1]: 0, be_s[0]: 0, be_s[1]: 1})
ZAHLEN["zweite_ordnung"] = {"Klammer {h_a, h_b}": str(sp.factor(soll)), "nativer Kommutator reell orthogonal":
                            float(np.linalg.norm(kom)), "Weyl-Relation Residuum": float(np.linalg.norm(lhs - rhs_w)),
                            "Wurzelbeträge D=+4": betr, "Überlapp su2_L|su2_R": float(ov_LR)}
check(f"(c) Frenkel-Kac: die Wurzelströme von so(4)_1 sind Verschiebungen mit |alpha|^2 = 2 (h = 1): gleicher Betrag des "
      f"D = +4-Anteils für alle vier Wurzeln, aber su(2)_L und su(2)_R polarisieren orthogonal (Überlapp {ov_LR:.1e}, innerhalb "
      f"su(2)_L {ov_LL:.4f}); die native Komposition ist Heisenberg: W(a)W(b) = e^(i Im a.b-quer) W(a + b) (Residuum "
      f"{np.linalg.norm(lhs - rhs_w):.1e}, Abschneidegrenze 10 Stufen), für reelle orthogonale a, b ist der Kommutator null ({np.linalg.norm(kom):.1e})",
      max(betr.values())/min(betr.values()) - 1 < 1e-9 and ov_LR < 1e-10 and abs(ov_LL - 1) < 1e-9
      and np.linalg.norm(lhs - rhs_w) < 1e-8 and np.linalg.norm(kom) < 1e-10)
check("(c) die Plebanski-Penrose-Klammer zweier Grad-2-Klassen ist {h_a, h_b} = (2/kappa^2) [a b] A^3 B^3 (z3/(z1^3 z2^2) + "
      "z4/(z1^2 z2^3) + 4/(z1^2 z2^2)) (symbolisch exakt): Grad 2, Pole in z1 und z2 (nichttriviale Klasse, Fock-Sektor D = -4 "
      "mit N_d = 2, 3 statt 0), antisymmetrisch und != 0 für a = (1,0), b = (0,1), wo der native Kommutator verschwindet: die "
      "zweite Ordnung der nativen Quellen ist nicht der Einstein-Vertex", sp.simplify(klammer - soll) == 0 and eul == 0 and pole
      and sp.simplify(spez) != 0)
befund("E/L (c): Die Nichtabelschheit von <J J J> = A f^abc (z75) sitzt im Kozyklus des Ladungsgitters (innen); außen sind die "
       "Ströme Verschiebungen, und deren Algebra ist die Heisenberg-Algebra mit zentralem Kommutator. Die native Seite liefert "
       "damit die Farbhälfte eines Doppelkopie-Schemas (f^abc), nicht die kinematische Hälfte (die Poissonklammer auf den "
       "omega-Fasern), die die Gravitation braucht", True)

# --- (d) Ladungssektor: sektoraufgelöste Oszillatordaten im fermionischen Fockraum


def sektor_aufgeloest(beta, alpha, R):
    xs = np.array([r + alpha for r in range(-R, R)])
    nm = len(xs)
    idx = np.arange(2**nm)
    bits = ((idx[:, None] >> np.arange(nm)[::-1]) & 1).astype(float)
    lc = -np.logaddexp(0, beta*xs)                                # log c(x)
    l1c = beta*xs - np.logaddexp(0, beta*xs)                      # log(1 - c(x))
    logw = bits @ lc + (1 - bits) @ l1c
    w = np.exp(logw - logw.max())
    Q = np.rint((bits - (xs < 0)).sum(1)).astype(int)
    f1 = (bits[:, 1:]*(1 - bits[:, :-1])).sum(1)                  # Diagonalteil von J_-1 J_1
    f2 = (bits[:, 2:]*(1 - bits[:, :-2])).sum(1)                  # Diagonalteil von J_-2 J_2
    out = {}
    for Qv in range(-2, 3):
        sel = Q == Qv
        Ws = w[sel].sum()
        out[Qv] = ((w[sel]*f1[sel]).sum()/Ws, (w[sel]*f2[sel]).sum()/Ws, Ws)
    return out


ok_d, sek_tab = True, {}
for beta_, R_, tol in ((4*np.pi, 7, 1e-8), (3.0, 7, 1e-5), (3.0, 8, 1e-6)):
    res = sektor_aufgeloest(beta_, 0.25, R_)
    soll1, soll2 = 1/np.expm1(beta_), 2/np.expm1(2*beta_)
    dev = max(max(abs(v[0]/soll1 - 1), abs(v[1]/soll2 - 1)) for v in res.values())
    ok_d &= dev < tol
    sek_tab[f"beta = {beta_:.4g}, R = {R_}"] = {"max rel. Abweichung über Q = -2..2": dev,
                                      "<J_-1 J_1>_Q": [res[Qv][0] for Qv in range(-2, 3)]}
    if beta_ > 10:
        pQ = np.array([res[Qv][2] for Qv in range(-2, 3)])
        pQ /= pQ.sum()
        qv_ = np.arange(-2, 3)
        mu_ = pQ @ qv_
        k2_ = pQ @ (qv_ - mu_)**2
        k4_ = pQ @ (qv_ - mu_)**4 - 3*k2_**2
        theta = np.exp(-beta_*(qv_**2/2 - 0.25*qv_))
        theta /= theta.sum()
        sek_tab["Ladung bei 4 pi"] = {"Exzess-Kurtosis": float(k4_/k2_**2), "Theta-Abweichung": float(np.abs(pQ - theta).max())}
        ok_d &= k4_/k2_**2 > 10 and np.abs(pQ - theta).max() < 1e-10
ZAHLEN["ladungssektor"] = sek_tab
# innere Operatoren (Q, Uhr R = e^{2 pi i L_0}) vertauschen mit allen Strommoden, also mit den daraus gebauten P_mu
xs_k = [r + 0.25 for r in range(-4, 4)]
nk = len(xs_k)
Zp = sps.diags([1.0, -1.0], format="csr")
sm = sps.csr_matrix(np.array([[0.0, 1.0], [0.0, 0.0]]))
psi = []
for j in range(nk):
    mats = [Zp]*j + [sm] + [sps.identity(2, format="csr")]*(nk - j - 1)
    op = mats[0]
    for Mx in mats[1:]:
        op = sps.kron(op, Mx, format="csr")
    psi.append(op)
num_k = [np.real((p_.T @ p_).diagonal()) for p_ in psi]
sea_k = [1.0 if x < 0 else 0.0 for x in xs_k]
Qk = sps.diags(sum(n_ - s_ for n_, s_ in zip(num_k, sea_k)))
L0k = sum(x*(n_ - s_) for x, n_, s_ in zip(xs_k, num_k, sea_k))
Rk = sps.diags(np.exp(2j*np.pi*L0k))


def Jk(m):
    op = sps.csr_matrix((2**nk, 2**nk))
    for j, x in enumerate(xs_k):
        if j + m in range(nk):
            op = op + psi[j].T @ psi[j + m]
    return op.tocsr()


ok_inn = all(abs(komm(Qk, Jk(m))).max() < 1e-12 and abs(komm(Rk, Jk(m))).max() < 1e-12 for m in (1, -1, 2, -2, 3))
check(f"(d) sektoraufgelöst (fermionischer Fockraum, 14 bzw. 16 Moden, Quarter alpha = 1/4): <J_-1 J_1>_Q = 1/(e^beta - 1) und <J_-2 J_2>_Q = "
      f"2/(e^(2 beta) - 1) in jedem Ladungssektor Q = -2..2 (beta = 4 pi: rel. Abweichung "
      f"{sek_tab['beta = 12.57, R = 7']['max rel. Abweichung über Q = -2..2']:.1e}; beta = 3: "
      f"{sek_tab['beta = 3, R = 7']['max rel. Abweichung über Q = -2..2']:.1e} -> "
      f"{sek_tab['beta = 3, R = 8']['max rel. Abweichung über Q = -2..2']:.1e} mit einem Modenpaar mehr, Abschneidefehler "
      f"~ e^(-beta R)), obwohl die Ladung selbst stark nichtgaußsch ist "
      f"(Exzess-Kurtosis {sek_tab['Ladung bei 4 pi']['Exzess-Kurtosis']:.1f}, Thetaverteilung); Q und die Uhr R = e^(2 pi i L_0) "
      "vertauschen mit allen Strommoden, also mit P_mu", ok_d and ok_inn)
befund("E (d): Die Kumulanten des Ladungssektors tragen keine Außendaten. Jeder rein innere Operator vertauscht mit P_mu, hat "
       "also nur Impulsübertrag null und kann keinen propagierenden masselosen Spin-2-Pol erzeugen (Hinweis des Koordinators, "
       "hier gerechnet); auf dem Außenraum wirkt er als c-Zahl im Sektor D = 0. Gravitation müsste aus einer gemischten "
       "Antwort innen/außen ((b), (c)) oder aus einem eigenen geometrischen Sektor kommen", True)

# --- GRAV.PARENT.GAP.NO_GO.01: der native Generator ist auf dem Graviton-Sektor gleichmäßig gegappt
K_grav = sorted({abs(4*np.pi*m_*D_) for m_ in range(1, 11) for D_ in (-4, 4)})
sel4 = SAFE[np.abs(DVAL[SAFE]) == 4]
p2_4 = abs(P2[:, sel4]).max()
check(f"GRAV.PARENT.GAP.NO_GO.01: im Sektor D = +-4 ist die Modulfrequenz |K| = 16 pi m >= 16 pi = {K_grav[0]:.4f} für alle "
      f"Moden m >= 1 (kein m -> 0, keine kritische Folge), die konforme Energie ist >= 3; die Masselosigkeit P^2 = 0 gilt dort "
      f"wie in jedem Sektor ({p2_4:.1e}), ist also kinematisch und kein Pol der nativen Dynamik",
      abs(K_grav[0] - 16*np.pi) < 1e-12 and p2_4 < 1e-9 and len(sel4) > 0)
befund("E/L: Im Sinn von GRAV.PARENT.GAP.NO_GO.01 ist der native Elterngenerator (Modulgenerator je Mode) auf dem Graviton-"
       "Sektor gleichmäßig gegappt; einen masselosen Spin-2-Pol der nativen Dynamik gibt es damit nicht. Die Čech-Koränder "
       "sind die lineare Eichfreiheit des Grad-2-Sektors (kinematisch), aber ohne kritische Folge; der zugelassene Ausweg "
       "(kritische Folge mit lokalen Tensor-Eichbedingungen) wird von keiner der Quellen (a)-(d) geliefert", True)

# =================================================================== Teil 3
say("\nTeil 3: S4b -- können die Marken die selbstduale Weylkrümmung einer gekrümmten Metrik sein?")


def IJ(a, b, c, d, e):                                            # Quartik a z^4 + 4b z^3 + 6c z^2 + 4d z + e (wie z61)
    return a*e - 4*b*d + 3*c*c, a*c*e + 2*b*c*d - a*d*d - b*b*e - c**3


# --- 3.1 Spin-2-Zustand <-> Quadrupol, Vierteldrehung (Entscheidungstest des Koordinators)
ex_, ey_, ez_ = np.eye(3)
pl, mi = ex_ + 1j*ey_, ex_ - 1j*ey_
ebas = {2: 0.5*np.outer(pl, pl), -2: 0.5*np.outer(mi, mi), 1: -0.5*(np.outer(pl, ez_) + np.outer(ez_, pl)),
        -1: 0.5*(np.outer(mi, ez_) + np.outer(ez_, mi)), 0: (2*np.outer(ez_, ez_) - np.outer(ex_, ex_) - np.outer(ey_, ey_))/math.sqrt(6)}
c_h6 = {2: 1/math.sqrt(2), -2: -1/math.sqrt(2)}
Q_h6 = sum(c*ebas[m] for m, c in c_h6.items())
ew_h6 = np.sort(np.linalg.eigvals(Q_h6/(1j)).real)
Rz = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)    # Vierteldrehung um die Markenachse
anti = np.allclose(Rz @ Q_h6 @ Rz.T, -Q_h6)
Bm = np.array([ebas[m].ravel() for m in (2, 1, 0, -1, -2)]).T     # 9 x 5
Lmat = np.linalg.lstsq(Bm, np.array([(Rz @ ebas[m] @ Rz.T).ravel() for m in (2, 1, 0, -1, -2)]).T, rcond=None)[0]
ew_R, ev_R = np.linalg.eig(Lmat)
inv = ev_R[:, np.abs(ew_R - 1) < 1e-9]
binom4 = {m: math.sqrt(math.comb(4, 2 + m)) for m in range(-2, 3)}


def quartik(cvec):                                                # Koeffizienten von z^4 .. z^0 (Majorana-Sterne wie z61)
    return [binom4[m]*cvec[i] for i, m in enumerate((2, 1, 0, -1, -2))]


def IJ_von(cvec):
    k = quartik(cvec)                                             # a z^4 + 4b z^3 + 6c z^2 + 4d z + e
    return IJ(k[0], k[1]/4, k[2]/6, k[3]/4, k[4])


I_h6, J_h6 = IJ_von([c_h6.get(m, 0) for m in (2, 1, 0, -1, -2)])
I_d, J_d = IJ_von(inv[:, 0]) if inv.shape[1] == 1 else (0, 0)
stab_t, stab_p = 0, 0
Qd = np.diag([1.0, -1.0, 0.0])                                    # H6 in Hauptachsen
for perm_ in itertools.permutations(range(3)):
    for sg in itertools.product((1, -1), repeat=3):
        Mo = np.zeros((3, 3))
        for i, j in enumerate(perm_):
            Mo[i, j] = sg[i]
        if abs(np.linalg.det(Mo) - 1) > 1e-9:
            continue
        stab_t += np.allclose(Mo @ Qd @ Mo.T, Qd)
        stab_p += np.allclose(Mo @ Qd @ Mo.T, Qd) or np.allclose(Mo @ Qd @ Mo.T, -Qd)
ZAHLEN["h6_quadrupol"] = {"Eigenwerte (bis Phase i)": ew_h6.tolist(), "I, J H6": [float(np.real(I_h6)), float(abs(J_h6))],
                          "R-invarianter Unterraum": int(inv.shape[1]), "I, J invariant": [complex(I_d), complex(J_d)],
                          "Stabilisator Tensor/Satz (Oktaedergruppe)": [stab_t, stab_p]}
check(f"Entscheidungstest S4b (Koordinator, gerechnet): der H6-Zustand (|2> - |-2>)/sqrt2 ist der Quadrupol diag(a, -a, 0) "
      f"(Eigenwerte {np.round(ew_h6, 6).tolist()} bis auf die Phase i; für 1 - z^4 ist I = -1, J = 0, hier normiert I = -1/2); die Vierteldrehung um die Markenachse "
      f"gibt R Q R^T = -Q; der R-invariante Unterraum des Spin-2-Raums ist eindimensional (m = 0, diag(a, a, -2a)) mit "
      f"I^3 = 27 J^2 != 0, also Typ D; in der Oktaedergruppe hat Q den Tensorstabilisator D2 ({stab_t} Drehungen) und den "
      f"Satzstabilisator D4 ({stab_p})", anti and np.allclose(ew_h6/ew_h6.max(), [-1, 0, 1], atol=1e-12) and abs(I_h6 + 0.5) < 1e-12
      and abs(J_h6) < 1e-12 and inv.shape[1] == 1 and abs(I_d**3 - 27*J_d**2) < 1e-9 and abs(J_d) > 1e-6
      and stab_t == 4 and stab_p == 8)
befund("E: S4b in der Form 'Marken = invariante selbstduale Weylkrümmung an einem Fixpunkt von R' ist ausgeschlossen. Eine "
       "Isometrie erhält W; die Vierteldrehung kehrt das Vorzeichen des H6-Quadrupols um (Eigenwert e^{-i pi m/2} = -1 auf "
       "m = +-2), ein R-invariantes W ist axialsymmetrisch (Typ D, J != 0) oder null. Die Marken haben als Satz D4, als Tensor "
       "nur D2", True)

# --- 3.2 Bianchi IX nach Cartan: e0 = N dt, ei = a_i sigma_i, d sigma_1 = sigma_2 ^ sigma_3 (zyklisch)
tt = sp.Symbol("t")
af = [sp.Function(f"a{i}")(tt) for i in range(1, 4)]
Nf = sp.Function("N")(tt)


def kruemmungsoperator(N, a):
    E = [N] + list(a)
    Dm = [sp.zeros(4, 4) for _ in range(4)]                       # de^a = (1/2) Dm[a]_bc e^b ^ e^c
    for i, j, k in ((1, 2, 3), (2, 3, 1), (3, 1, 2)):
        Dm[i][0, i] = sp.diff(E[i], tt)/(N*E[i])
        Dm[i][i, 0] = -Dm[i][0, i]
        Dm[i][j, k] = E[i]/(E[j]*E[k])
        Dm[i][k, j] = -Dm[i][j, k]
    w, unb = {}, []
    for p_ in range(4):
        for r_ in range(p_ + 1, 4):
            v = sp.symbols(f"w{p_}{r_}_0:4")
            unb += list(v)
            w[(p_, r_)], w[(r_, p_)] = list(v), [-x for x in v]
    for p_ in range(4):
        w[(p_, p_)] = [0]*4
    gl = []
    for p_ in range(4):
        M = sp.zeros(4, 4)
        for r_ in range(4):
            u = w[(p_, r_)]
            for c_ in range(4):
                for d_ in range(4):
                    M[c_, d_] -= u[c_]*(1 if d_ == r_ else 0) - u[d_]*(1 if c_ == r_ else 0)
        gl += [sp.Eq(M[c_, d_], Dm[p_][c_, d_]) for c_ in range(4) for d_ in range(c_ + 1, 4)]
    lsg = sp.solve(gl, unb, dict=True)[0]
    W_ = {k: [sp.simplify(sp.sympify(x).subs(lsg)) for x in v] for k, v in w.items()}

    def keil(u, v):
        return sp.Matrix(4, 4, lambda c_, d_: u[c_]*v[d_] - u[d_]*v[c_])

    def dform(u):
        M = sp.zeros(4, 4)
        for c_ in range(4):
            M += keil([sp.diff(u[c_], tt)/N, 0, 0, 0], [1 if x == c_ else 0 for x in range(4)]) + u[c_]*Dm[c_]
        return M
    Om = {(p_, r_): (dform(W_[(p_, r_)]) + sum((keil(W_[(p_, c_)], W_[(c_, r_)]) for c_ in range(4)), sp.zeros(4, 4)))
          .applyfunc(sp.simplify) for p_ in range(4) for r_ in range(4)}
    Rm = lambda p_, r_, c_, d_: Om[(p_, r_)][c_, d_]
    skal = sp.simplify(sum(Rm(p_, r_, p_, r_) for p_ in range(4) for r_ in range(4)))

    def Eb(k, sg):
        i, j = {1: (2, 3), 2: (3, 1), 3: (1, 2)}[k]
        M = sp.zeros(4, 4)
        M[0, k], M[k, 0], M[i, j], M[j, i] = 1, -1, sg, -sg
        return M

    def block(s1, s2):
        return sp.Matrix(3, 3, lambda k, l: sp.simplify(sum(Eb(k + 1, s1)[p_, r_]*Rm(p_, r_, c_, d_)*Eb(l + 1, s2)[c_, d_]
                                                            for p_, r_, c_, d_ in itertools.product(range(4), repeat=4))/8))
    return block(1, 1), block(-1, -1), block(1, -1), skal


A0, A1, A2 = sp.symbols("A0_1:4", positive=True), sp.symbols("A1_1:4"), sp.symbols("A2_1:4")
N0, N1 = sp.symbols("N0 N1", positive=True)


def jet(e):
    e = e.subs({sp.diff(af[i], tt, 2): A2[i] for i in range(3)})
    e = e.subs({**{sp.diff(af[i], tt): A1[i] for i in range(3)}, sp.diff(Nf, tt): N1})
    return sp.simplify(e.subs({af[i]: A0[i] for i in range(3)}).subs(Nf, N0))


Mpp, Mmm, Mpm, skal = kruemmungsoperator(Nf, af)
WP = (Mpp - skal/12*sp.eye(3)).applyfunc(jet)                     # W+ in der Basis e0^ei + ej^ek
WM = (Mmm - skal/12*sp.eye(3)).applyfunc(jet)
BJ, SJ = Mpm.applyfunc(jet), jet(skal)
diag_ok = all(WP[i, j] == 0 and WM[i, j] == 0 and BJ[i, j] == 0 for i in range(3) for j in range(3) if i != j)


def auswerten(Nexpr, aexpr, var, fluss):
    """Jets entlang eines Flusses d var_j/dx = fluss_j; x = var[0] ist die Zeit."""
    xv, Pv = var[0], var[1:]

    def Dx(f):
        return sp.diff(f, xv) + sum(sp.diff(f, Pv[j])*fluss[j] for j in range(len(Pv)))
    a1_ = [Dx(ai) for ai in aexpr]
    rep = {N0: Nexpr, N1: Dx(Nexpr)}
    for i in range(3):
        rep[A0[i]], rep[A1[i]], rep[A2[i]] = aexpr[i], a1_[i], Dx(a1_[i])
    return rep


rund = {A0[1]: A0[0], A0[2]: A0[0], A1[1]: A1[0], A1[2]: A1[0], A2[1]: A2[0], A2[2]: A2[0]}
ok_rund = all(sp.simplify(WP[i, i].subs(rund)) == 0 and sp.simplify(WM[i, i].subs(rund)) == 0 for i in range(3))
xT = sp.Symbol("x", positive=True)
Wv = sp.symbols("W1:4", positive=True)
kal = {}
for name, fl in (("Euler-Kreisel", [Wv[1]*Wv[2], Wv[2]*Wv[0], Wv[0]*Wv[1]]),
                 ("Halphen", [Wv[1]*Wv[2] - Wv[0]*(Wv[1] + Wv[2]), Wv[2]*Wv[0] - Wv[1]*(Wv[2] + Wv[0]),
                              Wv[0]*Wv[1] - Wv[2]*(Wv[0] + Wv[1])])):
    rep = auswerten(sp.sqrt(Wv[0]*Wv[1]*Wv[2]), [sp.sqrt(Wv[1]*Wv[2]/Wv[0]), sp.sqrt(Wv[2]*Wv[0]/Wv[1]),
                                                sp.sqrt(Wv[0]*Wv[1]/Wv[2])], (xT, *Wv), fl)
    kal[name] = (sp.simplify(SJ.subs(rep)) == 0 and all(sp.simplify(BJ[i, i].subs(rep)) == 0 for i in range(3))
                 and all(sp.simplify(WP[i, i].subs(rep)) == 0 for i in range(3))
                 and any(sp.simplify(WM[i, i].subs(rep)) != 0 for i in range(3)))
check(f"Bianchi-IX-Krümmung nach Cartan (allgemeine a_i(t), N(t), Jets): W+, W- und der spurfreie Ricci sind diagonal in "
      f"e0^ei +- ej^ek; Kontrollen: rund (a1 = a2 = a3) ist konform flach, Euler-Kreisel und Halphen (Gibbons-Pope, "
      f"Atiyah-Hitchin) sind Ricci-flach mit W+ = 0 und W- != 0 ({kal})", diag_ok and ok_rund and all(kal.values()))

# --- 3.3 R als Isometrie: sigma_1 -> sigma_2 -> -sigma_1 ist Isometrie nur für a1 = a2
biax = {A0[1]: A0[0], A1[1]: A1[0], A2[1]: A2[0]}
ok_biax = all(sp.simplify((W_[0, 0] - W_[1, 1]).subs(biax)) == 0 for W_ in (WP, WM))
lam_b = sp.Symbol("lambda")
Jb = sp.Matrix(sp.diag(lam_b, lam_b, -2*lam_b)).det()
check("R als Isometrie (Vierteldrehung sigma_1 -> sigma_2 -> -sigma_1) verlangt a1 = a2; dann ist für beliebige a(t), N(t) in "
      f"beiden Hälften W = diag(lambda, lambda, -2 lambda), also Typ D oder O; J = det W = {Jb} verschwindet nur für W = 0: "
      "harmonische Marken (J = 0, Typ I) und isometrisches R schließen sich in Bianchi IX aus",
      ok_biax and sp.solve(Jb, lam_b) == [0])

# --- 3.4 Tods Form der SU(2)-invarianten selbstdualen Konformstrukturen (Omega-System; Hitchin: Painlevé VI)
Pv = sp.symbols("P1:4", positive=True)
flT = [-Pv[1]*Pv[2]/(xT*(1 - xT)), -Pv[2]*Pv[0]/xT, -Pv[0]*Pv[1]/(1 - xT)]
repT = auswerten(1/sp.sqrt(xT*(1 - xT)), [1/Pv[0], sp.sqrt(1 - xT)/Pv[1], sp.sqrt(xT)/Pv[2]], (xT, *Pv), flT)
wmT = [sp.simplify(WM[i, i].subs(repT)) for i in range(3)]
lamT = [sp.factor(sp.simplify(WP[i, i].subs(repT))) for i in range(3)]
I0T = Pv[0]**2 - Pv[1]**2 - Pv[2]**2
dI0 = sp.simplify(sum(sp.diff(I0T, Pv[j])*flT[j] for j in range(3)))
pt_ = {xT: sp.Rational(1, 3), Pv[0]: sp.Rational(7, 10), Pv[1]: sp.Rational(1, 2), Pv[2]: sp.Rational(3, 10)}
lam_pt = [lv.subs(pt_) for lv in lamT]
check(f"Tods Form g0 = dx^2/(x(1-x)) + sigma_1^2/Omega_1^2 + (1-x) sigma_2^2/Omega_2^2 + x sigma_3^2/Omega_3^2 mit "
      f"Omega_1' = -Omega_2 Omega_3/(x(1-x)), Omega_2' = -Omega_3 Omega_1/x, Omega_3' = -Omega_1 Omega_2/(1-x): W- = 0 exakt "
      f"(symbolisch auf dem Fluss), I0 = Omega_1^2 - Omega_2^2 - Omega_3^2 ist Erhaltungsgröße, W+ = diag(lambda_k) mit "
      f"Summe 0; an einem triaxialen Punkt drei verschiedene Eigenwerte {[str(v) for v in lam_pt]} (Typ I)",
      all(v == 0 for v in wmT) and dI0 == 0 and sp.simplify(sum(lamT)) == 0 and len(set(lam_pt)) == 3)

# --- 3.5 J = lambda_1 lambda_2 lambda_3 identisch null? (Invarianz von {lambda_k = 0} unter dem Fluss)
numT = [sp.expand(sp.numer(sp.together(lv))) for lv in lamT]


def Xfl(f):
    g = sp.diff(f, xT) + sum(sp.diff(f, Pv[j])*flT[j] for j in range(3))
    return sp.expand(sp.numer(sp.together(g*xT*(1 - xT))))


ok_res, kleine_faktoren = True, set()
X123 = []
for k in range(3):
    X1 = Xfl(numT[k])
    X2 = Xfl(X1)
    X123.append((X1, X2, Xfl(X2)))
    g_ = sp.gcd(sp.resultant(numT[k], X1, xT), sp.resultant(numT[k], X2, xT))
    for fa, _ in sp.factor_list(g_)[1]:
        ok_res &= sp.Poly(fa, *Pv).total_degree() <= 2
        kleine_faktoren.add(str(fa))
lin_nicht_inv = True
for sg_, j_ in ((1, 1), (-1, 1), (1, 2), (-1, 2)):
    v_ = flT[0] - sg_*flT[j_]                                     # X(Omega_1 - sg Omega_j)
    lin_nicht_inv &= sp.simplify(v_.subs(Pv[0], sg_*Pv[j_])) != 0
check(f"exakt: für jedes k haben Res_x(N_k, X N_k) und Res_x(N_k, X^2 N_k) nur die gemeinsamen Faktoren "
      f"{sorted(kleine_faktoren)} (Grad <= 2): entartet (Omega_i = 0), nicht reell (Omega_2^2 + Omega_3^2) oder nicht "
      "flussinvariant (Omega_1 = +-Omega_2, +-Omega_3); eine zweidimensionale invariante Familie mit lambda_k = 0 gibt es nicht",
      ok_res and lin_nicht_inv)


def absterme_fn(e):
    """Skala sum_t |t| der Monome von e, numerisch (die Symbole P sind als positiv deklariert, Abs also erst nach Einsetzen)."""
    fs = sp.lambdify((xT, *Pv), list(sp.Add.make_args(e)), "mpmath")
    return lambda *v: mp.fsum(abs(t_) for t_ in fs(*v))


mp.mp.dps = 40
inv_typO, nicht_inv, unklar, n_C, nahe_O = 0, 0, 0, 0, 0
min_rel_nichtinv = mp.mpf(1)
for k in range(3):
    fN = sp.lambdify((xT, *Pv), [numT[k], X123[k][0], X123[k][1]], "mpmath")
    f3 = sp.lambdify((xT, *Pv), X123[k][2], "mpmath")
    a3 = absterme_fn(X123[k][2])
    fall = [sp.lambdify((xT, *Pv), numT[j], "mpmath") for j in range(3)]
    aall = [absterme_fn(numT[j]) for j in range(3)]
    fF = sp.lambdify((xT, *Pv), [numT[k], X123[k][0], X123[k][1]], "numpy")
    jF = sp.lambdify((xT, *Pv), sp.Matrix([numT[k], X123[k][0], X123[k][1]]).jacobian(sp.Matrix(Pv)), "numpy")
    for _ in range(250):
        xx0 = rng.uniform(0.05, 0.95)
        vv, info_, ier_, msg_ = fsolve(lambda v: fF(xx0, *v), rng.uniform(-2, 2, size=3),
                                       fprime=lambda v: np.array(jF(xx0, *v), dtype=float), full_output=True, xtol=1e-13)
        if ier_ != 1 or np.min(np.abs(vv)) < 1e-3 or np.abs(fF(xx0, *vv)).max() > 1e-8:
            continue
        xx = mp.mpf(xx0)
        try:
            p1, p2, p3 = mp.findroot(lambda a_, b_, c_: fN(xx, a_, b_, c_), [mp.mpf(t_) for t_ in vv], tol=mp.mpf(10)**-35)
        except (ValueError, ZeroDivisionError):
            continue
        n_C += 1
        rel3 = abs(f3(xx, p1, p2, p3))/a3(xx, p1, p2, p3)
        wgr = max(abs(fall[j](xx, p1, p2, p3))/aall[j](xx, p1, p2, p3) for j in range(3))   # Größe von W+ (relativ)
        if wgr < mp.mpf(10)**-20:                                 # konform flach (Typ O)
            inv_typO += rel3 < mp.mpf(10)**-15 and abs(p1**2 - p2**2 - p3**2 + mp.mpf(1)/4) < mp.mpf(10)**-12
            unklar += not (rel3 < mp.mpf(10)**-15)
        elif wgr > mp.mpf(10)**-3:                                # W+ deutlich != 0
            nicht_inv += rel3 > mp.mpf(10)**-9
            unklar += not (rel3 > mp.mpf(10)**-9)
            min_rel_nichtinv = min(min_rel_nichtinv, rel3)
        else:
            nahe_O += 1                                           # nahe am konform flachen Ort, entartete Konvergenz
mp.mp.dps = 30
ZAHLEN["J_null_invarianz"] = {"Punkte auf C = {N = XN = X^2N = 0}": n_C, "X^3N = 0 und W+ = 0 (I0 = -1/4)": inv_typO, "nahe Typ O": nahe_O,
                              "nicht invariant": nicht_inv, "unklar": unklar, "min |X^3 N| rel (nicht inv.)": mp.nstr(min_rel_nichtinv, 4)}
check(f"numerisch (fsolve, dann mpmath mit 40 Stellen nachpoliert): {n_C} Punkte der Kurve C = {{N_k = X N_k = X^2 N_k = 0}}; "
      f"wo W+ deutlich != 0 ist ({nicht_inv} Punkte), ist X^3 N_k != 0 (relativ >= {mp.nstr(min_rel_nichtinv, 2)}), die Kurve "
      f"also nicht invariant; X^3 N_k verschwindet nur zusammen mit ganz W+ ({inv_typO} Punkte, alle auf I0 = -1/4: konform "
      f"flach, Typ O; {nahe_O} Punkte nahe daran mit entarteter Konvergenz): J = 0 entlang einer Lösung gibt es nur konform "
      "flach [N]", n_C >= 40 and unklar == 0 and nicht_inv >= 10 and inv_typO >= 5)


def omega_fluss(xx, y):
    return [-y[1]*y[2]/(xx*(1 - xx)), -y[2]*y[0]/xx, -y[0]*y[1]/(1 - xx)]


def explodiert(xx, y):
    return 50 - np.max(np.abs(y))


explodiert.terminal = True
lamN = sp.lambdify((xT, *Pv), lamT, "numpy")
nullstellen, beschr = [], 0
for _ in range(40):
    y0 = rng.uniform(0.2, 1.2, size=3)*rng.choice([-1, 1], size=3)
    stuecke = []
    for xe in (0.97, 0.03):
        so = solve_ivp(omega_fluss, (0.5, xe), y0, rtol=1e-11, atol=1e-13, events=explodiert, dense_output=True)
        stuecke.append(so)
    if any(so.status != 0 for so in stuecke):
        continue
    beschr += 1
    xs_ = np.linspace(0.04, 0.96, 1500)
    ys_ = np.array([stuecke[0].sol(x_) if x_ >= 0.5 else stuecke[1].sol(x_) for x_ in xs_])
    L_ = np.array([lamN(x_, *y_) for x_, y_ in zip(xs_, ys_)])
    jr = L_[:, 0]*L_[:, 1]*L_[:, 2]/((L_**2).sum(1)/2)**1.5
    wechsel = np.where(np.sign(jr[1:]) != np.sign(jr[:-1]))[0]
    nullstellen.append((len(wechsel), float(np.max(np.abs(jr))), float(min(abs(jr[i + 1] - jr[i]) for i in wechsel)) if len(wechsel) else 0.0))
ZAHLEN["J_entlang_loesungen"] = {"beschränkte Lösungen": beschr, "(Nullstellen, max|J|/(I/2)^1.5, min Sprung)": nullstellen[:8]}
check(f"entlang {beschr} beschränkter triaxialer Lösungen auf (0,04; 0,96) wechselt J/(I/2)^(3/2) das Vorzeichen nur an "
      f"isolierten Stellen ({[n[0] for n in nullstellen]} Nullstellen) und erreicht dazwischen den Typ-D-Rand (max "
      f"{max(n[1] for n in nullstellen):.3f}): harmonische Marken (J = 0, tau = i) liegen auf einzelnen Bahnen, nicht überall",
      beschr >= 3 and all(n[0] >= 1 for n in nullstellen) and all(n[1] > 0.1 for n in nullstellen))
befund("E/N: In der SU(2)-symmetrischen Klasse (Bianchi IX, Tods Form) ist J = 0 mit W+ != 0 kein Dauerzustand: Die harmonische "
       "Konfiguration (tau = i) tritt auf einzelnen Bahnen auf, die Lösung läuft durch sie hindurch. J identisch null kommt nur "
       "konform flach vor (Niveau I0 = -1/4; Lesart: die runde S^4 mit der irreduziblen SO(3)-Wirkung auf Sym_0(3))", True)

# --- 3.6 Hitchins vier Punkte: Twistorraum der S^4 = P(Sym^3 C^2), Geraden = Büschel binärer Kubiken (f, sigma f)


def sig3(f):                                                      # (sigma f)(u, v) = conj(f)(-v, u), sigma^2 = -1
    g = np.zeros(4, complex)
    for k in range(4):
        g[3 - k] += (-1)**k*np.conj(f[k])
    return g


def jac33(f, g):                                                  # f_u g_v - f_v g_u, Koeffizienten von u^j v^(4-j)
    J_ = np.zeros(5, complex)
    for k in range(4):
        for l_ in range(4):
            if k >= 1 and l_ <= 2:
                J_[k - 1 + l_] += k*f[k]*(3 - l_)*g[l_]
            if l_ >= 1 and k <= 2:
                J_[k + l_ - 1] -= (3 - k)*f[k]*l_*g[l_]
    return J_


def disk3(c):
    a, b, c2, d = c[3], c[2], c[1], c[0]
    return b*b*c2*c2 - 4*a*c2**3 - 4*b**3*d - 27*a*a*d*d + 18*a*b*c2*d


def vier_punkte(f):
    g = sig3(f)
    J_ = jac33(f, g)
    zJ = np.roots(J_[::-1])
    ts = 1.3*np.exp(2j*np.pi*np.arange(9)/9)
    troots = np.roots(np.polyfit(ts, [disk3(f + t_*g) for t_ in ts], 4))
    doppel = []
    for t_ in troots:
        r_ = np.roots((f + t_*g)[::-1])
        i_, j_ = min(((i, j) for i in range(3) for j in range(i + 1, 3)), key=lambda pq: abs(r_[pq[0]] - r_[pq[1]]))
        doppel.append((r_[i_] + r_[j_])/2)
    err = max(min(abs(z_ - w_) for w_ in zJ) for z_ in doppel)
    anti_ = max(min(abs(-1/np.conj(z_) - w_) for w_ in zJ) for z_ in zJ)
    return err, anti_, J_


ok_hp = all(max(vier_punkte(rng.normal(size=4) + 1j*rng.normal(size=4))[:2]) < 1e-9 for _ in range(5))
ziel = np.array([-1, 0, 0, 0, 1], complex)/math.sqrt(2)          # u^4 - v^4: Wurzeln 1, i, -1, -i


def h6_res(p_):
    f = p_[:4] + 1j*p_[4:]
    J_ = jac33(f, sig3(f))
    ph = np.vdot(ziel, J_)/(abs(np.vdot(ziel, J_)) + 1e-300)
    r_ = J_/np.linalg.norm(J_) - ph*ziel
    return np.concatenate([r_.real, r_.imag, [np.linalg.norm(f) - 1]])


bestes = min((least_squares(h6_res, rng.normal(size=8), xtol=1e-15, ftol=1e-15, gtol=1e-15) for _ in range(20)),
             key=lambda r_: r_.cost)
f_h6 = bestes.x[:4] + 1j*bestes.x[4:]
err_h6, anti_h6, J_h6q = vier_punkte(f_h6)
wurz_h6 = np.roots(J_h6q[::-1])
ok_h6 = bestes.cost < 1e-25 and err_h6 < 1e-9 and np.allclose(np.sort(np.round(wurz_h6**4, 9).real), [1, 1, 1, 1])
ZAHLEN["hitchin_vier_punkte"] = {"Kubik f der H6-Geraden": [str(np.round(c, 8)) for c in f_h6],
                                 "Doppelwurzeln = Wurzeln der Jacobi-Quartik (Fehler)": err_h6}
check(f"Hitchins vier Punkte: im Büschel f + t sigma(f) binärer Kubiken (Twistorgerade der S^4 in P(Sym^3 C^2)) haben die vier "
      f"entarteten Glieder ihre Doppelwurzeln genau in den Wurzeln der Jacobi-Quartik J(f, sigma f) (5 Zufallsgeraden, Fehler "
      f"< 1e-9), und diese Quartik ist antipodisch reell (ein Punkt der S^4 = S(Sym^4 R)); es gibt eine Gerade mit "
      f"J(f, sigma f) ~ u^4 - v^4 (Fehler {err_h6:.1e}): die H6-Marken sind die vier Punkte, in denen die Twistorgerade über "
      "m_H6 die Diskriminantenfläche der SL(2,C)-Bahnen trifft", ok_hp and ok_h6)
befund("E (z65 + hier): Die Krümmungsquartik einer halbflachen Metrik liegt nie auf ihren eigenen Twistorgeraden: In Tods Form "
       "verschwindet W- (die Hälfte der Twistorfaser), die Marken einer Quartik von W+ sitzen auf der anderen Sphäre. Liest man "
       "die Naht als Twistorgerade, so lässt eine Weylquartik vom H6-Typ dort nach z65 (R(m_pi, m_pi) = 4 Psi~(pi^4)) nur die vier "
       "Marken als ungehinderte alpha-Ebenen übrig: Die Gerade selbst überlebt die Krümmung nicht", True)
befund("L/C: Die tragfähige Lesart von S4b ist eine Symmetrielesart, keine Krümmungslesart: In SO(3)-invarianten selbstdualen "
       "Einsteinräumen (Hitchin, Painlevé VI) trifft jede Twistorgerade den antikanonischen Divisor der entarteten Bahnen in vier "
       "Punkten; für die runde S^4 sind das genau die Wurzeln der Quartik, die den Punkt darstellt, und die H6-Marken gehören "
       "zum Punkt m_H6 auf der harmonischen Bahn. Das sind Daten der Symmetriewirkung, sie variieren mit der Bahn (Kreuzverhältnis "
       "= Bahnparameter) und sind harmonisch nur auf einer Bahn", True)
befund("E/L Gegenprobe S1: R ist auf den Außenoszillatoren trivial und innen ein Spin-Z4-Erzeuger (z79, z81). S4b verlangt R als "
       "Raumzeitsymmetrie: an Fixpunkten unmöglich (Quadrupol wechselt das Vorzeichen), als freie Isometrie erzwingt R Typ D. Die "
       "Entscheidung fällt damit für S1: Die Marken sind innen", True)

# =================================================================== Teil 4
say("\nTeil 4: Bilanz und Urteil")
STATUS = {
    "Außenraum = Elementarzustände": ("E", "Fockvektoren = Čech-Klassen, Homogenität = -K/(4 pi m) - 2, Tomita J = n -> -n - 4, "
                                      "Casimir 3(h^2 - 1) je Sektor; Kodaira-Spencer-Klassen = Sektor D = -4"),
    "S4 (a) gaußscher Sektor": ("E, negativ", "nur D = 0 (u(2,2), triviale Klassen), keine Deformation"),
    "S4 (b) Vertexoperatoren": ("E, negativ für Universalität", "Grad-2-Klassen erst in vierter Ordnung, Amplitude ~ h^2 e^{-ch}, "
                                "Polarisation = Ladungsvektor, KMS-Gewicht q^4 = e^{-8 pi m} gegen googly"),
    "S4 (c) nichtabelsche Ströme": ("E, negativ", "Verschiebungen mit |alpha|^2 = 2, su(2)_L und su(2)_R orthogonal polarisiert; "
                                    "f^abc innen; native zweite Ordnung Heisenberg, nicht die Plebanski-Penrose-Klammer"),
    "S4 (d) Ladungssektor": ("E, negativ", "Außendaten in jedem Sektor gleich, innere Operatoren vertauschen mit P_mu"),
    "GRAV.PARENT.GAP.NO_GO.01": ("E", "nativer Generator auf D = +-4 gleichmäßig gegappt (|K| >= 16 pi), keine kritische Folge"),
    "S4b Fixpunkt": ("E, ausgeschlossen", "R-invariantes W ist Typ D oder null; H6 hat J = 0"),
    "S4b Bianchi IX": ("E/N, ausgeschlossen", "isometrisches R => Typ D/O; J = 0 nur auf Bahnen, J identisch null nur konform flach"),
    "S4b Twistorgerade": ("E, ausgeschlossen", "Krümmungsquartik nicht auf der eigenen Twistorgeraden (googly, z65)"),
    "Hitchins vier Punkte": ("E/L", "H6-Marken = Schnitt der Twistorgeraden über m_H6 mit der Diskriminantenfläche; Symmetriedatum"),
    "Gravitation nichtlinear, G, Universalität": ("O", "braucht eine native Poissonstruktur auf den omega-Fasern oder einen eigenen "
                                                       "geometrischen Sektor; hier nicht vorhanden"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("Der Außenraum aus den Quarteroszillatoren mit Tomita-Partner ist als Fockraum genau die Elementarzustandsdarstellung "
          "der Twistortheorie: Fockvektoren sind Čech-Klassen, Vernichtung ist Korand, die Homogenität ist der negative "
          "Modulgenerator (n = -K/(4 pi m) - 2), die Tomita-Konjugation ist die Twistorkonjugation n -> -n - 4, und jeder "
          "Sektor D = 2h ist ein masseloser Einteilchensektor mit Casimir 3(h^2 - 1), Helizität +-2 eingeschlossen (wie S2/z83). "
          "Eine Kodaira-Spencer-Klasse erster Ordnung ist ein Vektor im Sektor D = -4. Keine native Quelle bevölkert ihn "
          "universell: Der gaußsche Sektor bleibt bei D = 0; Vertexoperatoren und nichtabelsche Ströme sind Twistorverschiebungen "
          "und erreichen den Grad-2-Sektor erst in vierter Ordnung, mit Amplitude ~ h^2 e^{-ch}, dem Ladungsvektor als "
          "Polarisation und dem KMS-Gewicht q^4 = e^{-8 pi m} gegen das googly-Bild; der Ladungssektor trägt keine Außendaten "
          "und hat nur Impulsübertrag null. Die native zweite Ordnung ist die Heisenberg-Komposition mit zentralem Kommutator, "
          "nicht die Plebanski-Penrose-Klammer {h_a, h_b} ~ [a b], die der Einstein-Vertex verlangt; der native Generator ist "
          "auf dem Graviton-Sektor gleichmäßig gegappt (GRAV.PARENT.GAP.NO_GO.01). S4 ist damit für die Quellen (a)-(d) gescopt "
          "negativ: Die Klasse koppelt nicht universell und die zweite Ordnung ist nicht Einstein. S4b ist in der getesteten Form "
          "ausgeschlossen: Ein R-invariantes Weyl-W an einem Fixpunkt ist axialsymmetrisch (Typ D) oder null, isometrisches R "
          "erzwingt in Bianchi IX Typ D/O, in Tods Form (W- = 0 exakt verifiziert) ist J = 0 nur auf einzelnen Bahnen und "
          "identisch nur konform flach, und die Krümmungsquartik liegt nie auf der eigenen Twistorgeraden. Tragfähig bleibt eine "
          "Symmetrielesart: Die H6-Marken sind die vier Punkte, in denen die Twistorgerade über dem Punkt m_H6 der S^4 die "
          "Diskriminantenfläche der SL(2,C)-Bahnen trifft (Hitchin). Die Entscheidung S1 gegen S4b fällt für S1: Die Marken sind "
          "innen.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z85_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
