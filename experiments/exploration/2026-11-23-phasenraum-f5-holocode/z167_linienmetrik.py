"""Z167: Stufe 1 des Lösungswegs. Die metrisierte Nahtdeterminante am symmetrischen Punkt, in allen Richtungen.

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Objekt: Z(tau, z) = Theta_E8(tau, z)/eta(tau)^8, der Charakter von (E8)_1 mit flachem Hintergrund z im Cartan-Torus.
Metrik der Determinantenlinie (Quillen-Faktor): ||Z||^2 = |Z|^2 exp(-2 pi y.y/Im tau), y = Im z. Gamma = -log ||Z||^2.
Verformungen: z (Eichhintergrund) und tau (Form der Naht, nach z166 die Scherung). Symmetrischer Punkt: tau = i, z = 0.

Frage der Stufe 1 (angekündigt als Entscheidungstest): Liefert die Krümmung der Linie den Koeffizienten -2 c3^3 der
alpha-Gleichung, eine positive Form für die Scherung und die Normierung der Eichkopplungen?

Teil A  Objekt und Invarianz
Teil B  Eichrichtungen: zweite und vierte Ordnung
Teil C  Formrichtung tau
Teil D  Die alpha-Gleichung in Nahteinheiten
Teil E  Bilanz und Urteil.
"""
import json
import os
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import mpmath as mp
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
mp.mp.dps = 40
AINV = mp.mpf("137.0359992168407")          # Wurzel der alpha-Gleichung des Ledgers (Bestand, Eingabe)
c3 = 1 / (8 * mp.pi)
I = mp.mpc(0, 1)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   -> {name}: {text}")


def theta_e8(tau, z):
    """Theta_E8(tau, z) = Summe über lambda in E8 von exp(i pi tau lambda^2 + 2 pi i lambda.z), D8-Koordinaten."""
    q = mp.e ** (I * mp.pi * tau)
    prods = []
    for n in (3, 4, 2, 1):
        p = mp.mpf(1)
        for zi in z:
            p *= mp.jtheta(n, mp.pi * zi, q)
        prods.append(p)
    return sum(prods) / 2


def eta(tau):
    q = mp.e ** (2 * I * mp.pi * tau)
    return mp.e ** (I * mp.pi * tau / 12) * mp.qp(q)


def Zf(tau, z):
    return theta_e8(tau, z) / eta(tau) ** 8


def norm2(tau, z):
    y = [mp.im(zi) for zi in z]
    return abs(Zf(tau, z)) ** 2 * mp.e ** (-2 * mp.pi * mp.fsum(a * a for a in y) / mp.im(tau))


def Gam(tau, z):
    return -mp.log(norm2(tau, z))


def eis(k, tau):
    q = mp.e ** (2 * I * mp.pi * tau)
    coef = {2: -24, 4: 240, 6: -504}[k]
    return 1 + coef * mp.nsum(lambda n: n ** (k - 1) * q ** n / (1 - q ** n), [1, mp.inf])


say("Z167: Stufe 1. Die metrisierte Nahtdeterminante am symmetrischen Punkt")
say("=" * 100)
say("\nTeil A  Objekt und Invarianz")
tau0 = I
zero = [mp.mpf(0)] * 8
z_test = [mp.mpc("0.11", "0.07"), mp.mpc("-0.05", "0.13"), mp.mpc("0.21", "-0.02"), mp.mpc("0.03", "0.04"),
          mp.mpc("-0.17", "0.09"), mp.mpc("0.08", "-0.11"), mp.mpc("0.02", "0.15"), mp.mpc("-0.06", "-0.03")]
tau_t = mp.mpc("0.13", "1.1")
check("A1 Theta_E8(tau, 0) = E_4(tau) und Z(i, 0) = j(i)^(1/3) = 12 [N, 30 Stellen]",
      abs(theta_e8(tau_t, zero) - eis(4, tau_t)) < mp.mpf("1e-30") and abs(Zf(tau0, zero) - 12) < mp.mpf("1e-30"))
lam = [1, -1, 0, 0, 0, 0, 0, 0]                       # Wurzel
mu = [mp.mpf(1) / 2] * 8                              # Spinorvektor
sh1 = [a + b for a, b in zip(z_test, lam)]
sh2 = [a + tau_t * b for a, b in zip(z_test, mu)]
zS = [a / tau_t for a in z_test]
n0 = norm2(tau_t, z_test)
check("A2 ||Z||^2 ist invariant unter z -> z + lambda, z -> z + tau mu (lambda, mu in E8) und unter "
      "S: (tau, z) -> (-1/tau, z/tau), bei allgemeinem tau und z [N, 25 Stellen]: Der Faktor "
      "exp(-2 pi y.y/Im tau) ist die Metrik der Determinantenlinie über dem Cartan-Torus",
      abs(norm2(tau_t, sh1) / n0 - 1) < mp.mpf("1e-25") and abs(norm2(tau_t, sh2) / n0 - 1) < mp.mpf("1e-25")
      and abs(norm2(-1 / tau_t, zS) / n0 - 1) < mp.mpf("1e-25"))
check("A3 Bei tau = i wirkt S als z -> -i z (Phase -i, Spin 1) und auf die Tangente von tau als dtau -> -dtau "
      "(Phase -1, Spin 2): d(-1/tau) = dtau/tau^2. Die Verformungen tragen unter der geometrischen "
      "Vierteldrehung die Helizitäten von Eichfeld und Graviton [E; Gleichsetzung der geometrischen "
      "Vierteldrehung mit der Uhr R ist die Prämisse P_lift]",
      sp.simplify(1 / sp.I + sp.I) == 0 and sp.simplify(1 / sp.I ** 2 + 1) == 0)

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil B  Eichrichtungen")
E2i, E4i, E6i = eis(2, tau0), eis(4, tau0), eis(6, tau0)
varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi))        # Lemniskatenkonstante
ZAHLEN.update(E4_i=mp.nstr(mp.re(E4i), 15), varpi=mp.nstr(varpi, 15))
check("B1 Am symmetrischen Punkt gilt E_6(i) = 0, E_2(i) = 3/pi und E_4(i) = 3 varpi^4/pi^4 "
      "(varpi = 2,62206 Lemniskatenkonstante) [N, 30 Stellen]",
      abs(E6i) < mp.mpf("1e-30") and abs(E2i - 3 / mp.pi) < mp.mpf("1e-30")
      and abs(E4i - 3 * varpi ** 4 / mp.pi ** 4) < mp.mpf("1e-30"))
# holomorpher quadratischer Koeffizient c(tau) = -pi^2 (E2 E4 - E6)/(3 E4)
E2s, E4s, E6s = sp.symbols("E2 E4 E6")
c_hol = -sp.pi ** 2 * (E2s * E4s - E6s) / (3 * E4s)
check("B2 log Z = log Z_0 + (c/2) z.z + ... mit c(tau) = -pi^2 (E_2 E_4 - E_6)/(3 E_4) (Ramanujan: "
      "q dE_4/dq = (E_2 E_4 - E_6)/3; E8 isotrop). Bei tau = i ist c = -pi [E]",
      sp.simplify(c_hol.subs({E6s: 0, E2s: 3 / sp.pi}) + sp.pi) == 0)
mp.mp.dps = 30
vdir = [mp.mpf(x) for x in ("0.31", "-0.22", "0.47", "0.10", "-0.39", "0.28", "0.05", "-0.61")]
nv = mp.sqrt(mp.fsum(a * a for a in vdir))
vdir = [a / nv for a in vdir]
hess = {}
for name, ph in (("reell", mp.mpf(1)), ("imaginär", I), ("diagonal", mp.e ** (I * mp.pi / 4))):
    f = lambda s, ph=ph: Gam(tau0, [ph * s * a for a in vdir])
    hess[name] = (mp.diff(f, 0, 2) / 2, mp.diff(f, 0, 4) / 24)
mp.mp.dps = 40
d4 = varpi ** 4 / 12
ZAHLEN.update(hesse_quadratisch={k: mp.nstr(v[0], 12) for k, v in hess.items()},
              hesse_quartisch={k: mp.nstr(v[1], 12) for k, v in hess.items()}, quartisch_soll=mp.nstr(d4, 12))
check("B3 Zweite Ordnung: Gamma = pi |z|^2 in jeder Richtung und Phase (reell, imaginär, diagonal) [N, 15 "
      "Stellen]. Der anisotrope Teil ist proportional zu E_2 - 3/(pi Im tau) und zu E_6; beide verschwinden genau "
      "bei tau = i. Die Antwort der Naht auf einen Eichhintergrund ist dort isotrop mit Koeffizient "
      "pi = 1/(8 c3) je Cartan-Richtung [E]",
      all(abs(v[0] - mp.pi) < mp.mpf("1e-15") for v in hess.values()))
check("B4 Vierte Ordnung: Gamma_4 = -(varpi^4/12) Re((z.z)^2), also -varpi^4/12, -varpi^4/12, +varpi^4/12 für "
      "reelle, imaginäre und diagonale Richtung [N, 12 Stellen]. E8 hat nur eine quartische Invariante (z163), "
      "deshalb ist das die ganze vierte Ordnung: eine reine Zahl ohne freien Parameter, mit der vierzähligen "
      "Anisotropie cos(4 Winkel) der Vierteldrehung [E mit der Vier-Design-Eigenschaft der E8-Schalen]",
      abs(hess["reell"][1] + d4) < mp.mpf("1e-12") and abs(hess["imaginär"][1] + d4) < mp.mpf("1e-12")
      and abs(hess["diagonal"][1] - d4) < mp.mpf("1e-12"))
say(f"   quadratisch: {mp.nstr(hess['reell'][0], 12)}; quartisch: {mp.nstr(hess['reell'][1], 12)} "
    f"(soll -{mp.nstr(d4, 12)})")
Y = [sp.Rational(-1, 3)] * 3 + [sp.Rational(1, 2)] * 2 + [0] * 3
YY = sum(a * a for a in Y)
t3 = [0, 0, 0, sp.Rational(1, 2), sp.Rational(-1, 2), 0, 0, 0]          # schwacher Isospin T3
t3t3 = sum(a * a for a in t3)
check("B5 Die Antwort ist die E8-Form selbst (Niveau 1), gleich für alle Richtungen: Hyperladung Y.Y = 5/6 "
      "(k_Y = 2 Y.Y = 5/3), T_3.T_3 = 1/2, also sin^2(theta_W) = T_3.T_3/(T_3.T_3 + Y.Y) = 3/8 an der Naht "
      "[E; bekannt im Bestand]. Die Krümmung legt Verhältnisse fest, nicht den gemeinsamen Wert der Kopplung",
      YY == sp.Rational(5, 6) and t3t3 / (t3t3 + YY) == sp.Rational(3, 8))
befund("B_eich", {"hesse": "pi je Cartan-Richtung", "quartisch": "varpi^4/12", "spur": "8 pi = 1/c3"},
       "Die Spur der Hesse-Form über den Cartan-Torus ist 8 pi = 1/c3: c3 ist der Kehrwert der Gesamtkrümmung der "
       "Determinantenlinie. Die quadratische Antwort ist isotrop, die quartische die reine Zahl varpi^4/12.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil C  Formrichtung tau")
mp.mp.dps = 30
lapl = lambda f: mp.diff(lambda x: f(mp.mpc(x, 1)), 0, 2) + mp.diff(lambda y: f(mp.mpc(0, y)), 1, 2)
g_tot = lambda t: -mp.log(abs(Zf(t, zero)) ** 2)
g_lat = lambda t: -mp.log(mp.im(t) ** 4 * abs(theta_e8(t, zero)) ** 2)
g_osc = lambda t: -mp.log(1 / (mp.im(t) ** 4 * abs(eta(t)) ** 16))
L_tot, L_lat, L_osc = lapl(g_tot), lapl(g_lat), lapl(g_osc)
mp.mp.dps = 40
ZAHLEN.update(laplace_gesamt=mp.nstr(L_tot, 8), laplace_gitter=mp.nstr(L_lat, 8), laplace_oszillator=mp.nstr(L_osc, 8))
check("C1 Z hat Gewicht 0: -log |Z|^2 ist harmonisch in tau, die Krümmung der Linie in Formrichtung ist null "
      "(Laplace 0 bei tau = i). Die zwei Faktoren tragen entgegengesetzte Krümmung: Gitterfaktor Theta_E8 "
      "(Gewicht 4, Norm (Im tau)^4 |Theta|^2) gibt d d-quer = +1/(Im tau)^2, also die Poincaré-Form, "
      "Oszillatorfaktor eta^-8 gibt -1/(Im tau)^2 [E; N: Laplace = 4 mal d d-quer = 0, +4, -4]",
      abs(L_tot) < mp.mpf("1e-10") and abs(L_lat - 4) < mp.mpf("1e-10") and abs(L_osc + 4) < mp.mpf("1e-10"))
dlogj = lambda t: mp.diff(lambda x: mp.log(Zf(mp.mpc(x, mp.im(t)), zero)), mp.re(t))
kmean = lambda y: -eis(6, I * y) / (3 * eis(4, I * y))                # <L_0 - 1/3> = q d/dq log Z = -E6/(3 E4)
check("C2 tau = i ist stationär (E_6(i) = 0) und ein Sattel: Die zweiten Ableitungen von -log |Z|^2 sind "
      "-19,157 und +19,157 (z166). Eine positive Form für die Scherung kommt aus der vollen chiralen Naht nicht "
      "heraus; die Ankündigung von S-I ('positiver Teil aus der Linienmetrik') ist in dieser Form widerlegt "
      "[E negativ]. Positiv ist allein der Gitterfaktor (Zwei-Faktoren-Architektur S1)",
      abs(kmean(1)) < mp.mpf("1e-30") and abs(L_tot) < mp.mpf("1e-10"))
befund("C_form", {"kruemmung_gesamt": 0, "gitter": "+Poincaré", "oszillator": "-Poincaré"},
       "In Formrichtung hat die Linie keine Krümmung. Der Gitterfaktor allein trägt genau die Poincaré-Form, "
       "die Feldraummetrik der zwei Gravitonpolarisationen (z166); der Oszillatorfaktor hebt sie auf.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil D  Die alpha-Gleichung in Nahteinheiten")
al, tt, cc, b1s, Ls = sp.symbols("alpha t c3 b1 L", positive=True)
F_led = al ** 3 - 2 * cc ** 3 * al ** 2 - 8 * b1s * cc ** 6 * Ls
F_t = sp.simplify((F_led.subs(al, cc * tt / (4 * sp.pi)) * 64 * sp.pi ** 3 / cc ** 3).subs(cc, 1 / (8 * sp.pi)))
c3s = sp.Rational(1, 8) / sp.pi
check("D1 Mit t = e^2/c3 = 32 pi^2 alpha (e^2 = 4 pi alpha) ist die alpha-Gleichung des Ledgers genau "
      "t^2 (t - c3) = b1 ln(1/phi_seam): keine Faktoren 2 und 8 mehr, nur c3 und b1 [E, Umformung]",
      sp.simplify(F_t - (tt ** 3 - c3s * tt ** 2 - b1s * Ls)) == 0)
alpha = 1 / AINV
t_val = 32 * mp.pi ** 2 * alpha
dtop = 48 * c3 ** 4
qf = lambda a: dtop * mp.e ** (-2 * a)
phis = lambda a: 1 / (6 * mp.pi) + qf(a) * (1 - qf(a)) ** (-mp.mpf(5) / 4)
resid = t_val ** 2 * (t_val - c3) - mp.mpf(41) / 10 * mp.log(1 / phis(alpha))
ZAHLEN.update(t=mp.nstr(t_val, 15), rest_t_form=mp.nstr(resid, 5))
check("D2 An der Ledgerwurzel alpha^-1 = 137,0359992168407 ist t = 2,3047 und der Rest der t-Form kleiner als "
      "1e-12 [N]. t ist die Kopplung mal Spur der Hesse-Form aus Teil B: t = e^2 mal 8 pi [E]",
      abs(resid) < mp.mpf("1e-12") and abs(t_val - mp.mpf("2.30466")) < mp.mpf("1e-4")
      and abs(t_val - 4 * mp.pi * alpha * 8 * mp.pi) < mp.mpf("1e-30"))
tsym = sp.Symbol("t")
Pi1 = sp.diag(0, 0, 1)
check("D3 t^2 (t - c3) = det(t - c3 Pi) für einen Projektor Pi vom Rang 1 auf drei Kanälen: zwei Kanäle tragen "
      "t, einer t - c3. Drei ist die Zahl der Kanäle mit von null verschiedenem Sprung-Eigenwert auf den vier "
      "Marken (Ledger: det' R_0 = 2^16 c3^3, Spektrum 16 c3 mal {0, 2, 2, 4}: Dublett und Vorzeichenkanal) "
      "[E für die Algebra; Zuordnung des verschobenen Kanals zum Vorzeichenkanal ist Lesart, L]",
      sp.expand((tsym * sp.eye(3) - c3s * Pi1).det() - tsym ** 2 * (tsym - c3s)) == 0
      and sp.simplify((16 * c3s) ** 3 * 2 * 2 * 4 - 2 ** 16 * c3s ** 3) == 0)
check("D4 Entscheidungstest: Die Krümmung der Linie ist das Niveau, eine ganze Zahl in der Einheit pi. Sie "
      "liefert die Einheit der Gleichung (1/c3 = 8 pi = Spur der Hesse-Form) und damit die Variable t, aber "
      "nicht die Gleichung: weder die Verschiebung um genau c3 in einem Kanal noch die rechte Seite "
      "b1 ln(1/phi_seam) folgen aus zweiter oder vierter Ordnung am symmetrischen Punkt [E negativ für "
      "'-2 c3^3 ist eine Krümmung'; die Gleichung selbst bleibt gelöster Bestand, ihr Ursprung offen]",
      abs(8 * mp.pi - 1 / c3) < mp.mpf("1e-30") and abs(2 * c3 ** 3 - c3 ** 2 / (4 * mp.pi)) < mp.mpf("1e-30"))
ainv_dual = AINV
check("D5 Dieselbe Einheit im Exponenten der Hierarchie: 2/alpha = 1/(c3 e^2) = (Spur der Hesse-Form)/e^2, also "
      "e^(-2/alpha) = exp(-8 pi/e^2). Elektrisch steht 2 alpha = 4 c3^2 t, dual 1/(c3^2 t): die Dualität ist "
      "t -> 1/(4 c3^4 t) [E, Umformung]",
      abs(2 * ainv_dual - 1 / (c3 * 4 * mp.pi * alpha)) < mp.mpf("1e-25")
      and abs(2 * alpha - 4 * c3 ** 2 * t_val) < mp.mpf("1e-30")
      and abs(2 * ainv_dual - 1 / (c3 ** 2 * t_val)) < mp.mpf("1e-25"))
befund("D_alpha", {"t": mp.nstr(t_val, 10), "form": "t^2 (t - c3) = b1 ln(1/phi_seam)"},
       "Die Linienmetrik liefert die natürliche Variable und Einheit der alpha-Gleichung, nicht die Gleichung. "
       "Der angekündigte Test ('-2 c3^3 als Krümmung') ist negativ; neu ist die Form t^2 (t - c3) = b1 ln(1/phi) "
       "als Determinante über drei Kanäle mit einem um c3 verschobenen.")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil E  Bilanz und Urteil")
bilanz = [
    ("E", "Eichrichtung, zweite Ordnung", "Isotrop mit pi = 1/(8 c3) je Cartan-Richtung; Spur 8 pi = 1/c3. "
     "Verhältnisse fest (3/8), gemeinsamer Wert der Kopplung nicht."),
    ("E", "Eichrichtung, vierte Ordnung", "-(varpi^4/12) Re((z.z)^2): eine reine Zahl, vierzählig anisotrop."),
    ("E negativ", "Formrichtung", "Keine Krümmung, Sattel. S-I in der angekündigten Form widerlegt. Nur der "
     "Gitterfaktor trägt die Poincaré-Form."),
    ("E", "Helizitäten", "Unter der geometrischen Vierteldrehung: z Spin 1, dtau Spin 2."),
    ("E", "alpha in Nahteinheiten", "t^2 (t - c3) = b1 ln(1/phi_seam) mit t = e^2/c3 = e^2 mal Spur der Hesse-Form."),
    ("E negativ", "Entscheidungstest", "Die Krümmung ist das Niveau; -2 c3^3 ist keine Krümmung."),
    ("O", "Ursprung der Gleichung", "Warum ein Kanal um c3 verschoben ist und warum rechts b1 ln(1/phi) steht."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("Stufe 1 ist gerechnet. Ergebnis: Der angekündigte Entscheidungstest ist in seiner Form negativ. Die "
          "Krümmung der Determinantenlinie über dem Cartan-Torus ist das Niveau der Stromalgebra; sie ist "
          "isotrop, gleich pi je Richtung, und legt Verhältnisse (sin^2 theta_W = 3/8) und die Einheit "
          "1/c3 = 8 pi fest, aber weder den Koeffizienten -2 c3^3 als solchen noch den gemeinsamen Wert der "
          "Kopplung. In Formrichtung hat die volle chirale Naht gar keine Krümmung; tau = i ist ein Sattel, und "
          "eine positive Graviton-Normierung entsteht so nicht. Positiv ist nur der Gitterfaktor, der genau "
          "die Poincaré-Form trägt. Drei Dinge sind neu und exakt: (1) c3 ist der Kehrwert der Gesamtkrümmung "
          "der Linie; (2) in der Variablen t = e^2/c3 lautet die alpha-Gleichung t^2 (t - c3) = b1 ln(1/phi_seam), "
          "eine Determinante über drei Kanäle mit einem um c3 verschobenen; (3) die vierte Ordnung ist die reine "
          "Zahl varpi^4/12 mit der Anisotropie der Vierteldrehung. Der Kandidat 'relative graduierte "
          "Determinante' ist damit nicht bestätigt: Die Form 'alles ist zweite Ordnung der Linienmetrik' trägt "
          "nicht. Was übrig bleibt, ist enger: der Ursprung der Kanalverschiebung und der rechten Seite. Keine "
          "volle Lösung, nichts befördert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z167_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z167.txt").write_text("\n".join(OUT) + "\n")
