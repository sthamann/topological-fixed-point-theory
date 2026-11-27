"""Z140: Maßbilanz -- Erzeugt die Theorie (Korpus plus H6) aus ihrer Topologie, Mathematik und Geometrie selbst Maße:
Planck-Länge, Zeiteinheiten, Skalenverhältnisse?

Leitlinie 'Rekonstruktion statt Auswahl' und die Anti-Numerologie-Regel des Ledgers (v354/v355): Eine Zahl durch Probieren
von Ausdrücken zu treffen ist kein Ergebnis; eine Relation zählt nur, wenn ein Mechanismus sie erzeugt. Wo hier überhaupt
Ausdrücke abgetastet werden, stehen Größe des Raums und erwartete Zufallstreffer dabei. Gelöste Blöcke (alpha^-1 =
137,03599921684 aus der Fixpunktgleichung, c_3 = 1/(8 pi), phi_0 = 1/(6 pi) + 48 c_3^4, die 16 in drei Familien, E8, die
Massenverhältnisse) werden wie gegeben benutzt und nicht neu angepasst. Nichts wird befördert; es wird keine vollständige
Theorie behauptet. Gelesen (nicht verändert): verification/v3, v7, v18, v20, v46, v60, v153, v164, v274, v725, das Ledger,
z33, z48, z68, z96, z98, z99, z104, z112, z118, z119, z121. Nicht wiederholt: z136 (Verschränkungsgleichgewicht, G) und
z137 (induzierte Gravitation, Nahtzählung). Statuskürzel: E exakt, C bedingt, L Lesart, O offen, N numerisch.

Teil 1  Inventar: alle nativen dimensionslosen Zahlen mit Maßbedeutung (Ort, Bedeutung, exakt?) und alle nativen Größen
        mit Dimension; zwei Dimensionsmatrizen: (A) physikalische Dimension modulo hbar, c, k_B (ein Torsor?), (B)
        unabhängige Skalen nach Einsetzen nur der vorhandenen Relationen (wie v725); Negativkontrollen.
Teil 2  Zeit: was macht den dimensionslosen Modulparameter zur physikalischen Zeit? beta = 4 pi mit l = R_Vak/2; Vergleich
        mit Gibbons-Hawking/Unruh; natives Energiequant und nativer Takt der Uhr R in Einheiten R_Vak/c.
Teil 3  Tabelle aller Skalenverhältnisse mit Wert, Quelle, ehrlichem Status; mit einem Anker (CODATA G) in SI umgerechnet.
Teil 4  Gibt es einen nativen Mechanismus für die eine große Zahl G/R_Vak^2 (Ziel e^(-2/alpha)/(128 pi^3))?
        (a) Boltzmann/KMS, (b) dimensionale Transmutation, (c) Instanton, (d) Zählung; Wirkung des Vorzeichenproblems.
Teil 5  Urteil in drei Teilen.
"""
import itertools
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.linalg import expm

HERE = Path(__file__).parent
REPO = HERE.parents[2]
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(140)
mp.mp.dps = 40


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


def lade(z):
    return json.loads((HERE / f"{z}_ergebnis.json").read_text())


def hat(d, teil, feld="checks"):
    return any(teil in k and v for k, v in d.get(feld, {}).items())


def ns(x, n=6):
    return mp.nstr(mp.mpf(x), n)


# ------------------------------------------------------------------ feste Eingaben (nur gelesen, nicht angepasst)
PI = mp.pi
c3 = 1/(8*PI)
dtop = 48*c3**4
phi0 = 1/(6*PI) + dtop
lamY = mp.sqrt(phi0*(1 - phi0))


def F_alpha(a, M=41, rueck=None):
    """F_U(1) wie verification/v3_em_alpha.py (make_F); rueck ersetzt die Nahtrückwirkung e^(-2 a) (nur für Kontrollen)."""
    Q = dtop*(mp.e**(-2*a) if rueck is None else rueck(a))
    phiseam = 1/(6*PI) + Q*(1 - Q)**(mp.mpf(-5)/4)
    return a**3 - 2*c3**3*a**2 - (mp.mpf(4)/5)*c3**6*M*mp.log(1/phiseam)


AINV = 1/mp.findroot(F_alpha, mp.mpf("0.0073"))
r_nat = mp.e**(-4*PI)

# SI (hbar, c, k_B, eV exakt; G gemessen) -- der eine Anker ist CODATA G
HBAR, CL, KB, EV = mp.mpf("1.054571817e-34"), mp.mpf("299792458"), mp.mpf("1.380649e-23"), mp.mpf("1.602176634e-19")
G_SI = mp.mpf("6.67430e-11")                       # CODATA 2018/2022, rel. Unsicherheit 2,2e-5
LJ, JAHR = mp.mpf("9.4607304725808e15"), mp.mpf("31557600")
GEV = EV*mp.mpf("1e9")
lPl = mp.sqrt(HBAR*G_SI/CL**3)
tPl = lPl/CL
Mbar_kg = mp.sqrt(HBAR*CL/(8*PI*G_SI))
Mbar_GeV = Mbar_kg*CL**2/GEV
MPl_GeV = Mbar_GeV*mp.sqrt(8*PI)
# externe Messwerte (keine Theorieausgaben; nur für Vergleiche und für die Umrechnung der Flavourzeile)
V_EW = mp.mpf("246.21965")                          # GeV, (sqrt2 G_F)^(-1/2)
PDG = {"e": mp.mpf("0.51099895e-3"), "mu": mp.mpf("0.1056583755"), "tau": mp.mpf("1.77686")}   # GeV
H0_PLANCK, DH0, OL_PLANCK, DOL = mp.mpf("67.36"), mp.mpf("0.54"), mp.mpf("0.6847"), mp.mpf("0.0073")  # Planck 2018
MPC = mp.mpf("3.0856775814913673e22")

say("Z140: Maßbilanz -- erzeugt die Theorie selbst Maße?")
say("\nTeil 0: feste Eingaben")
check(f"[E] alpha^-1 als Wurzel von F_U(1) (v3, Budget 41) = {ns(AINV, 15)}; Soll 137,03599921684; c_3 = 1/(8 pi), delta_top = "
      f"48 c_3^4 = {ns(dtop, 8)}, phi_0 = {ns(phi0, 10)}. Kontrolle: Budget 40 gibt {ns(1/mp.findroot(lambda a: F_alpha(a, 40), mp.mpf('0.0073')), 8)}",
      abs(AINV - mp.mpf("137.03599921684")) < mp.mpf("1e-10")
      and abs(1/mp.findroot(lambda a: F_alpha(a, 40), mp.mpf("0.0073")) - AINV) > 0.5
      and abs(dtop - 3/(256*PI**4)) < mp.mpf("1e-35"))
v3_text = (REPO / "verification" / "v3_em_alpha.py").read_text()
check("[E] Quelle gelesen: v3 enthält die Rückwirkung als 'dt * mp.e**(-2 * a)' und den Exponenten -5/4 (die hier benutzte "
      "Funktion ist dieselbe)", "Q = dt * mp.e**(-2 * a)" in v3_text and "mp.mpf(-5) / 4" in v3_text)
check(f"[E] Anker: aus CODATA G folgen l_Pl = {ns(lPl, 7)} m, t_Pl = {ns(tPl, 7)} s, Mbar = (8 pi G)^(-1/2) = {ns(Mbar_GeV, 7)} GeV "
      f"= {ns(Mbar_kg, 6)} kg (Korpus: 2,435323e18 GeV)", abs(Mbar_GeV/mp.mpf("2.435323203e18") - 1) < mp.mpf("2e-5")
      and abs(lPl/mp.mpf("1.616255e-35") - 1) < mp.mpf("1e-6"))

# =================================================================== Teil 1
say("\nTeil 1: Inventar der nativen Zahlen und Größen")
z33, z48, z68, z98, z104, z112, z118, z119, z121 = (lade(z) for z in
                                                   ("z33", "z48", "z68", "z98", "z104", "z112", "z118", "z119", "z121"))
kappa2 = mp.gamma(mp.mpf(1)/4)**8/(384*PI**4)
ddv_nat = 4*mp.e**(-PI/2)
# (Name, Wert, Ort, was es misst, exakt?)
INVENTAR = [
    ("beta = 4 pi", 4*PI, "z46/z58/z94, VERDICT 'KMS bei beta = 4 pi'", "KMS-Wert der Naht zur Uhr E = L_0", "exakt"),
    ("r = e^(-4 pi)", r_nat, "z47/z48/z118", "Boltzmann-Faktor je Stromquant; Gewicht der massiven Linie", "exakt"),
    ("Uhrordnung 4 (R^4 = 1)", mp.mpf(4), "z82/z113/z118", "R = e^(2 pi i L_0), L_0 in Z/4", "exakt"),
    ("Lücke 4 pi L_0 = pi (L_0 = 1/4)", PI, "VERDICT (Jacobi-Tripelprodukt)", "niedrigste Anregung in Moduleinheiten", "exakt"),
    ("c_3 = 1/(8 pi)", c3, "tfpt_constants.py (P1)", "Nahtnormierung; 8 pi G = 1/c_3 in v_geo-Einheiten", "exakt (Eingabe P1)"),
    ("delta_top = 48 c_3^4", dtop, "tfpt_constants.py", "topologischer Term des Seeds", "exakt"),
    ("phi_0 = 1/(6 pi) + 48 c_3^4", phi0, "tfpt_constants.py", "Seed der Flavourleiter (Schrittweite der Massen)", "exakt"),
    ("lambda_Y = sqrt(phi_0 (1 - phi_0))", lamY, "v18", "Yukawa-Schritt je Wortlänge", "exakt"),
    ("alpha^-1", AINV, "v3_em_alpha.py", "Kopplung; einzige Wurzel von F_U(1)", "Wurzel exakt bestimmt, Lesart C"),
    ("Indexspur 64 (32 je J-Summand)", mp.mpf(64), "z118", "ganzzahliges Gewicht des eps-Terms", "exakt"),
    ("|Lambda| R_Vak^2 = 12", mp.mpf(12), "z98", "Krümmung des nativen AdS in Vakuumeinheiten", "exakt"),
    ("l/R_Vak = 1/2", mp.mpf(1)/2, "z98", "AdS-Radius gegen Vakuumgröße Im x_Vak", "exakt"),
    ("beta l/R_Vak = 2 pi", 2*PI, "z140 Teil 2", "physikalische inverse Temperatur am Zentrum", "exakt"),
    ("beta_Original/beta_dS = 2", mp.mpf(2), "z104", "KMS-Wert der Uhr gegen Gibbons-Hawking der Kappe", "exakt, ohne Mechanismus"),
    ("S je Quarter = ln 2/4", mp.log(2)/4, "z68", "Entropie am KMS-Punkt (Naht ln 2/2, E8-Naht ln 12)", "exakt"),
    ("kappa_2 = Gamma(1/4)^8/(384 pi^4)", kappa2, "z68/VERDICT", "Hesse-Form der Modulfunktion", "exakt"),
    ("log 3", mp.log(3), "z116/z119", "Vielfachheitsterm der Spur (dim Isospin 1), keine Energie", "exakt"),
    ("M_f/M_g = 1/sqrt(e^(4 pi) - 1)", 1/mp.sqrt(mp.e**(4*PI) - 1), "z48/z96", "bimetrische Lesart der Gewichte", "exakt als Zahl, Lesart L"),
    ("DdV: m a = 4 e^(-pi Theta), Theta = 1/2", ddv_nat, "z121", "Massenformel des Lichtkegelgitters am nativen Punkt", "Formel L, Zahl exakt"),
    ("v = 2 tanh(pi/2)", 2*mp.tanh(PI/2), "z121", "Spinonengeschwindigkeit der Nahtkette (Plätze je Periode)", "Lesart L"),
    ("Leiterexponenten K = (4,2,0), (4,3,2), (5,3,2)", mp.mpf(25), "v46", "phi_0-Potenzen der neun Massen; Summe 25", "exakt (innerhalb der Annahmen)"),
    ("Leptonkoeffizienten (16/7, 4/3, 7/6)", mp.mpf(32)/9, "v20", "Produkt 32/9", "exakt (innerhalb der Annahmen)"),
    ("c_3^(7/2)", c3**(mp.mpf(7)/2), "v7 (GRAV.SCAL.01)", "M_scal/Mbar", "Zahl exakt, Zuordnung N"),
    ("48 c_3^2 e^(-2/alpha)", 48*c3**2*mp.e**(-2*AINV), "v7/v60 (COSMO.LAM.01)", "rho_Lambda/Mbar^4", "Zahl exakt, Zuordnung N/C"),
    ("(1 - 1/(4 pi)) phi_0", (1 - 1/(4*PI))*phi0, "v7 (COSMO.OMB.01)", "Omega_b (Anteil, kein Skalenverhältnis)", "Zahl exakt, Zuordnung N"),
    ("e^(3/4)", mp.e**mp.mpf("0.75"), "v152", "m/mu der Naht-EH-Normierung", "exakt"),
]
say(f"   {'Zahl':46s} {'Wert':>16s}  Ort; Bedeutung; exakt?")
for nm, val, ort, bed, ex in INVENTAR:
    say(f"   {nm:46s} {ns(val, 10):>16s}  {ort}; {bed}; {ex}")
ZAHLEN["inventar_dimensionslos"] = [{"Zahl": nm, "Wert": ns(val, 15), "Ort": ort, "misst": bed, "exakt": ex}
                                    for nm, val, ort, bed, ex in INVENTAR]
check(f"[E] Inventar ({len(INVENTAR)} native dimensionslose Zahlen), Abgleich mit dem Bestand: r = e^(-4 pi) = {ns(r_nat, 8)}; "
      f"M_f/M_g = {ns(1/mp.sqrt(mp.e**(4*PI) - 1), 10)} (z48: {z48['zahlen']['alpha_bimetrik']:.10f}); S je Quarter = ln 2/4 = "
      f"{ns(mp.log(2)/4, 10)} (z68: {z68['zahlen']['kms_punkt']['S_je_quarter']}); kappa_2 = {ns(kappa2, 10)} (z68: "
      f"{z68['zahlen']['kms_punkt']['kappa2']}); DdV m a = {ns(ddv_nat, 10)} (z121); l = 1/2, Lambda = -12 (z98); 64 je Index (z118)",
      abs(1/mp.sqrt(mp.e**(4*PI) - 1)/mp.mpf(z48["zahlen"]["alpha_bimetrik"]) - 1) < mp.mpf("1e-10")
      and abs(mp.log(2)/4 - mp.mpf(z68["zahlen"]["kms_punkt"]["S_je_quarter"])) < mp.mpf("1e-10")
      and abs(kappa2 - mp.mpf(z68["zahlen"]["kms_punkt"]["kappa2"])) < mp.mpf("1e-10")
      and abs(ddv_nat - mp.mpf(z121["zahlen"]["kontinuum"]["DdV-Masse m a = 4 e^(-pi Theta) bei Theta = 1/2 (Lesart)"])) < mp.mpf("1e-12")
      and z98["zahlen"]["radius"]["l (Zustandseinheiten, Im x_Vak = 1)"] == "1/2"
      and z98["zahlen"]["radius"]["Lambda AdS"] == "-12"
      and hat(z118, "Die relative Spur liefert nur ganze Zahlen"))
lam_ = sp.Symbol("lambda", positive=True)
Rv_, l_, Lam_, bet_, eps_, rho_, Mb_, Gi_, al_ = sp.symbols("R_Vak l Lambda beta_phys epsilon rho_Lambda Mbar Ginv alpha", positive=True)
MDIM = {Rv_: -1, l_: -1, Lam_: 2, bet_: -1, eps_: 1, rho_: 4, Mb_: 1, Gi_: 2, al_: 0}     # Massendimension (hbar = c = k_B = 1)


def skaliere(expr, dim=MDIM):
    return expr.subs({v: lam_**(-d)*v for v, d in dim.items()}, simultaneous=True)


RELATIONEN = [   # (Name, linke Seite = 0, Homogenitätsgrad)
    ("l = R_Vak/2 (z98)", l_ - Rv_/2, -1), ("Lambda = -12/R_Vak^2 (z98)", Lam_ + 12/Rv_**2, 2),
    ("beta_phys = 2 pi R_Vak (Teil 2)", bet_ - 2*sp.pi*Rv_, -1), ("Stromquant = 2/R_Vak (Teil 2)", eps_ - 2/Rv_, 1),
    ("1/G = 8 pi Mbar^2 (v274)", Gi_ - 8*sp.pi*Mb_**2, 2),
    ("rho_L = 48 c_3^2 e^(-2/alpha) Mbar^4 (v7)", rho_ - sp.Rational(3, 4)/sp.pi**2*sp.exp(-2/al_)*Mb_**4, 4),
    ("P_Lambda: R_Vak^2/G = 128 pi^3 e^(2/alpha) (z96/z98)", Rv_**2*Gi_ - 128*sp.pi**3*sp.exp(2/al_), 0),
]
hom = {n: sp.simplify(skaliere(e) - lam_**(-d)*e) == 0 for n, e, d in RELATIONEN}
mut = Rv_ - sp.Rational(5, 2)                                          # 'R_Vak = reine Zahl' würde lambda festlegen
mut_hom = [sp.simplify(skaliere(mut) - lam_**(-d)*mut) == 0 for d in (-1, 0, 1)]
falsch = sp.simplify(skaliere(Gi_ - 8*sp.pi*Mb_**2, {**MDIM, Gi_: 1}) - lam_**(-2)*(Gi_ - 8*sp.pi*Mb_**2)) == 0
check(f"[E] No-Unit in Maschinenform für H6 und Korpus: Alle {len(RELATIONEN)} Relationen zwischen Größen mit Dimension sind "
      f"unter L -> lambda L homogen (Grade {[d for _, _, d in RELATIONEN]}); keine legt lambda fest, auch P_Lambda nicht (Grad 0: "
      "sie fixiert nur die reine Zahl G/R_Vak^2). Die Zahlen des Inventars haben Massendimension 0. Negativkontrollen: 'R_Vak = 5/2' (Länge "
      f"gleich reine Zahl) ist in keinem Grad homogen ({mut_hom}), und mit falsch typisiertem 1/G (Dimension 1) fällt "
      f"die Newton-Relation durch ({falsch})",
      all(hom.values()) and len(hom) == 7 and not any(mut_hom) and not falsch
      and sp.simplify(skaliere(Rv_**2*Gi_) - Rv_**2*Gi_) == 0 and sp.simplify(skaliere(Gi_/Rv_**2) - Gi_/Rv_**2) != 0)

# --- Matrix A: physikalische Dimension (L, T, M, Theta, I) modulo der exakt festgelegten Konstanten hbar, c, k_B
DIM = [
    ("R_Vak (Vakuumgröße Im x_Vak)", [1, 0, 0, 0, 0]), ("l = R_Vak/2 (AdS-Radius)", [1, 0, 0, 0, 0]),
    ("Lambda = -12/R_Vak^2", [-2, 0, 0, 0, 0]), ("Takt der Uhr R", [0, 1, 0, 0, 0]),
    ("T_phys der Naht", [0, 0, 0, 1, 0]), ("Energiequant", [2, -2, 1, 0, 0]),
    ("Nahtlänge L, Kragenbreite L/2, Gitterabstand a", [1, 0, 0, 0, 0]),
    ("G (Koeffizient der Einsteinwirkung)", [3, -2, -1, 0, 0]), ("Mbar = v_geo", [0, 0, 1, 0, 0]),
    ("rho_Lambda (Energiedichte)", [-1, -2, 1, 0, 0]), ("M_scal, f_a, m_f, v_EW, M_R (Massen)", [0, 0, 1, 0, 0]),
    ("Lambda_QCD, T_reh (Energie, Temperatur)", [2, -2, 1, 0, 0]),
]
KONST = [("hbar", [2, -1, 1, 0, 0]), ("c", [1, -1, 0, 0, 0]), ("k_B", [2, -2, 1, -1, 0])]
MA = sp.Matrix([d for _, d in DIM] + [d for _, d in KONST])
MK = sp.Matrix([d for _, d in KONST])
rang_A = MA.rank() - MK.rank()
MA_laenge = sp.Matrix([d for _, d in DIM] + [[1, 0, 0, 0, 0]] + [d for _, d in KONST])
MA_ladung = sp.Matrix([d for _, d in DIM] + [[0, 1, 0, 0, 1]] + [d for _, d in KONST])
rang_A_ohne_c = sp.Matrix([d for _, d in DIM] + [KONST[0][1], KONST[2][1]]).rank() - 2
ZAHLEN["matrix_A"] = {"Zeilen": [n for n, _ in DIM], "Rang modulo (hbar, c, k_B)": int(rang_A),
                      "mit zusätzlicher fiktiver Länge": int(MA_laenge.rank() - MK.rank()),
                      "mit neuer Basisdimension (Ladung A s)": int(MA_ladung.rank() - MK.rank()),
                      "ohne c als Umrechnung": int(rang_A_ohne_c)}
check(f"[E] Matrix A (physikalische Dimension, {len(DIM)} Zeilen über (L, T, M, Theta, I)): Rang modulo der exakt gesetzten "
      f"hbar, c, k_B = {rang_A}: genau ein Einheiten-Torsor R_+. Eine zusätzliche fiktive Länge ändert daran nichts (Rang "
      f"{MA_laenge.rank() - MK.rank()}): in A ist jede Länge dieselbe Einheit. Negativkontrollen: eine neue Basisdimension "
      f"(Ladung) hebt den Rang auf {MA_ladung.rank() - MK.rank()}, und ohne c als Umrechnung sind Länge und Zeit getrennt (Rang "
      f"{rang_A_ohne_c})", rang_A == 1 and MA_laenge.rank() - MK.rank() == 1 and MA_ladung.rank() - MK.rank() == 2
      and rang_A_ohne_c == 2)

# --- Matrix B: unabhängige Skalen. Spalten = Längen, die keine vorhandene Relation aufeinander zurückführt:
#     (R_Vak, l_G = sqrt(hbar G/c^3), l_EW = hbar/(v_EW c), l_fiktiv). Zeilen = Exponenten.
SKALEN = ["R_Vak", "l_G", "l_EW", "l_fiktiv"]
H6_NATIV = {"l": [1, 0, 0, 0], "Lambda": [-2, 0, 0, 0], "Takt": [1, 0, 0, 0], "beta_phys": [1, 0, 0, 0],
            "Energiequant": [-1, 0, 0, 0], "Markenschritt": [1, 0, 0, 0]}
H6_G = {"G": [0, 2, 0, 0]}
KORPUS = {"Mbar": [0, -1, 0, 0], "M_scal": [0, -1, 0, 0], "f_a": [0, -1, 0, 0], "rho_Lambda": [0, -4, 0, 0]}
FLAVOUR = {"m_f (neun Massen)": [0, 0, -1, 0]}


def rang(*bloecke, ersetze=None):
    rows = [list(v) for b in bloecke for v in b.values()]
    if ersetze:                                              # Relation: Spalte j ist Zahl mal Spalte i
        for (j, i) in ersetze:
            for row in rows:
                row[i] += row[j]
                row[j] = 0
    return int(sp.Matrix(rows).rank())


rB = {
    "H6 nativ (ohne G)": rang(H6_NATIV),
    "H6 mit G (MM-Koeffizient frei)": rang(H6_NATIV, H6_G),
    "H6 mit G unter P_Lambda (COSMO.LAM.01) und P_cap": rang(H6_NATIV, H6_G, ersetze=[(1, 0)]),
    "Korpus, zertifizierter Block": rang(KORPUS),
    "Korpus mit Flavour-Absolutwert": rang(KORPUS, FLAVOUR),
    "alles, nur vorhandene Relationen": rang(H6_NATIV, H6_G, KORPUS, FLAVOUR),
    "alles unter P_Lambda": rang(H6_NATIV, H6_G, KORPUS, FLAVOUR, ersetze=[(1, 0)]),
    "alles unter P_Lambda und abgeleitetem Mbar/v_EW": rang(H6_NATIV, H6_G, KORPUS, FLAVOUR, ersetze=[(1, 0), (2, 0)]),
    "Negativkontrolle: H6 nativ plus fiktive unabhängige Länge": rang(H6_NATIV, {"fiktiv": [0, 0, 0, 1]}),
    "Negativkontrolle: alles plus fiktive Länge": rang(H6_NATIV, H6_G, KORPUS, FLAVOUR, {"fiktiv": [0, 0, 0, 1]}),
}
ZAHLEN["matrix_B"] = {"Spalten": SKALEN, "Ränge": rB}
for k, v in rB.items():
    say(f"      Rang {v}: {k}")
check("[E] Matrix B (unabhängige Skalen nach Einsetzen nur der vorhandenen Relationen): H6 nativ hat Rang "
      f"{rB['H6 nativ (ohne G)']} (alles ist Zahl mal R_Vak); mit G Rang {rB['H6 mit G (MM-Koeffizient frei)']}, also genau "
      "eine freie reine Zahl G/R_Vak^2; der zertifizierte Korpusblock hat Rang "
      f"{rB['Korpus, zertifizierter Block']} (wie v725), mit dem Flavour-Absolutwert {rB['Korpus mit Flavour-Absolutwert']} "
      f"(freie Zahl Mbar/v_EW); alles zusammen {rB['alles, nur vorhandene Relationen']}, unter P_Lambda "
      f"{rB['alles unter P_Lambda']}, erst mit abgeleitetem Mbar/v_EW 1. Negativkontrolle: eine fiktive unabhängige Länge hebt "
      f"den Rang ({rB['H6 nativ (ohne G)']} -> {rB['Negativkontrolle: H6 nativ plus fiktive unabhängige Länge']}, "
      f"{rB['alles, nur vorhandene Relationen']} -> {rB['Negativkontrolle: alles plus fiktive Länge']})",
      rB["H6 nativ (ohne G)"] == 1 and rB["H6 mit G (MM-Koeffizient frei)"] == 2
      and rB["H6 mit G unter P_Lambda (COSMO.LAM.01) und P_cap"] == 1 and rB["Korpus, zertifizierter Block"] == 1
      and rB["Korpus mit Flavour-Absolutwert"] == 2 and rB["alles, nur vorhandene Relationen"] == 3
      and rB["alles unter P_Lambda"] == 2 and rB["alles unter P_Lambda und abgeleitetem Mbar/v_EW"] == 1
      and rB["Negativkontrolle: H6 nativ plus fiktive unabhängige Länge"] == 2
      and rB["Negativkontrolle: alles plus fiktive Länge"] == 4)
befund("E: Es gibt keine zweite unabhängige Einheit: Die physikalische Dimension aller nativen Größen ist ein einziger Torsor "
       "(Matrix A, Rang 1). Was fehlt, sind keine Einheiten, sondern reine Zahlen: G/R_Vak^2 (Gravitationssektor von H6) und "
       "Mbar/v_EW (Flavour-Absolutwert), dazu die im Ledger als externe Eingaben geführten M_R, T_reh und alpha_s(M_Z) "
       "(Matrix B, Rang 3 statt 1)", rang_A == 1 and rB["alles, nur vorhandene Relationen"] == 3)

# =================================================================== Teil 2
say("\nTeil 2: Zeit")
t, x, y, z = sp.symbols("t x y z", real=True)
CO = [t, x, y, z]
r2 = x**2 + y**2 + z**2
ETA = [-1, 1, 1, 1]


def ricci_skalar(Om, eta):
    """Ricci-Skalar von g = eta/Om^2, direkt aus den Christoffel-Symbolen."""
    g = sp.diag(*[e/Om**2 for e in eta])
    gi = sp.diag(*[Om**2/e for e in eta])
    Gam = [[[sp.cancel(sum(gi[a, d]*(sp.diff(g[d, b], CO[c]) + sp.diff(g[d, c], CO[b]) - sp.diff(g[b, c], CO[d]))
                           for d in range(4))/2) for c in range(4)] for b in range(4)] for a in range(4)]

    def ric(b, c):
        return sum(sp.diff(Gam[a][b][c], CO[a]) - sp.diff(Gam[a][b][a], CO[c])
                   + sum(Gam[a][a][d]*Gam[d][b][c] - Gam[a][c][d]*Gam[d][b][a] for d in range(4)) for a in range(4))
    return sp.cancel(sum(gi[b, b]*ric(b, b) for b in range(4)))


def lie(xi, Om, eta):
    g = sp.diag(*[e/Om**2 for e in eta])
    return sp.Matrix(4, 4, lambda a, b: sp.cancel(sum(xi[c]*sp.diff(g[a, b], CO[c]) + g[c, b]*sp.diff(xi[c], CO[a])
                                                      + g[a, c]*sp.diff(xi[c], CO[b]) for c in range(4))))


OM_ADS, OM_DS = 1 + t**2 - r2, 1 - t**2 + r2                         # z98: Omega = 1 +- det x, Einheit Im x_Vak = 1
R_ads, R_ds, R_kugel = ricci_skalar(OM_ADS, ETA), ricci_skalar(OM_DS, ETA), ricci_skalar((1 + t**2 + r2)/6, [1, 1, 1, 1])
l_nat = sp.sqrt(sp.Rational(-12)/R_ads)
check(f"[E] natives AdS (z98, Einheit R_Vak = Im x_Vak = 1): g = eta/Omega^2 mit Omega = 1 + det x hat Ricci-Skalar {R_ads} = "
      f"-12/l^2, also l = {l_nat} R_Vak und Lambda = -3/l^2 = {-3/l_nat**2}/R_Vak^2; die dS-Wahl Omega = 1 - det x hat {R_ds}. "
      f"Kontrolle der Routine: runde S^4 mit Radius 3 gibt {R_kugel} = 12/9",
      R_ads == -48 and R_ds == 48 and l_nat == sp.Rational(1, 2) and R_kugel == sp.Rational(4, 3))
XI_E = [(1 + t**2 + r2)/2, t*x, t*y, t*z]                            # E = (P_0 + K_0)/2, z104
XI_Q = [(1 - t**2 - r2)/2, -t*x, -t*y, -t*z]                         # Q^- = (P_0 - K_0)/2, z104
Z4 = sp.zeros(4, 4)
kill = {"E in AdS": lie(XI_E, OM_ADS, ETA) == Z4, "Q^- in dS": lie(XI_Q, OM_DS, ETA) == Z4,
        "E in dS": lie(XI_E, OM_DS, ETA) == Z4, "Q^- in AdS": lie(XI_Q, OM_ADS, ETA) == Z4}
norm_E = sp.cancel(sum(ETA[a]*XI_E[a]**2 for a in range(4))/OM_ADS**2)
norm_Q = sp.cancel(sum(ETA[a]*XI_Q[a]**2 for a in range(4))/OM_DS**2)
nE_achse, nQ_achse = sp.cancel(norm_E.subs({x: 0, y: 0, z: 0})), sp.cancel(norm_Q.subs({x: 0, y: 0, z: 0}))
check(f"[E] Die Uhr E = (P_0 + K_0)/2 ist Killingfeld genau des nativen AdS, die Kappe Q^- = (P_0 - K_0)/2 genau des dS "
      f"(Lie-Ableitung der Metrik = 0: {kill}). Auf der Achse x = 0 ist g(xi_E, xi_E) = {nE_achse} und g(xi_Q, xi_Q) = "
      f"{nQ_achse}, konstant längs der ganzen Bahn: Eigenzeit = l mal Flussparameter, exakt, mit l = 1/2",
      kill["E in AdS"] and kill["Q^- in dS"] and not kill["E in dS"] and not kill["Q^- in AdS"]
      and nE_achse == sp.Rational(-1, 4) and nQ_achse == sp.Rational(-1, 4))
th, s_ = sp.symbols("theta s", real=True)
ode_E = sp.simplify(sp.diff(sp.tan(th/2), th) - (1 + sp.tan(th/2)**2)/2) == 0
ode_Q = sp.simplify(sp.diff(sp.tanh(s_/2), s_) - (1 - sp.tanh(s_/2)**2)/2) == 0
per = {"tan(theta/2): reelle Periode 2 pi": abs(mp.tan((mp.mpf("0.7") + 2*PI)/2) - mp.tan(mp.mpf("0.35"))),
       "tanh(s/2): imaginäre Periode 2 pi": abs(mp.tanh((mp.mpf("0.7") + 2j*PI)/2) - mp.tanh(mp.mpf("0.35"))),
       "tan(theta/2) bei theta + 4 pi i (keine Periode)": abs(mp.tan((mp.mpf("0.7") + 4j*PI)/2) - mp.tan(mp.mpf("0.35"))),
       "tanh(s/2) bei s + 4 pi i (zweite Periode)": abs(mp.tanh((mp.mpf("0.7") + 4j*PI)/2) - mp.tanh(mp.mpf("0.35")))}
ZAHLEN["zeit_perioden"] = {k: ns(v, 4) for k, v in per.items()}
check("[E] Bahnen auf der Achse: Uhr t = tan(theta/2) (reelle Periode 2 pi, keine imaginäre: bei theta + 4 pi i Abweichung "
      f"{ns(per['tan(theta/2) bei theta + 4 pi i (keine Periode)'], 3)}), Kappe t = tanh(s/2) (imaginäre Periode 2 pi, Abweichung "
      f"{ns(per['tanh(s/2): imaginäre Periode 2 pi'], 2)}). Die 2 pi der Kappe ist eine geometrische Periode (Gibbons-Hawking), "
      "die 4 pi der Uhr ist ein Datum des Zustands (KMS), keine Periode des Flusses",
      ode_E and ode_Q and per["tan(theta/2): reelle Periode 2 pi"] < mp.mpf("1e-30")
      and per["tanh(s/2): imaginäre Periode 2 pi"] < mp.mpf("1e-30")
      and per["tan(theta/2) bei theta + 4 pi i (keine Periode)"] > mp.mpf("0.5"))

# dieselben Erzeuger in der Vektordarstellung von so(4,2) (Konstruktion wie z96)
Gm = np.diag([-1., -1, 1, 1, 1, 1])
E6 = np.eye(6)


def Mg(A, B):
    return np.outer(E6[A], Gm @ E6[B]) - np.outer(E6[B], Gm @ E6[A])


P0v, K0v = Mg(1, 5) + Mg(1, 0), Mg(1, 5) - Mg(1, 0)
kand = {s: (P0v + s*K0v)/2 for s in (+1, -1)}
klasse = {}
for s in (+1, -1):
    ev = np.linalg.eigvals(kand[s])
    klasse[s] = ("elliptisch" if np.max(np.abs(ev.real)) < 1e-12 and np.max(np.abs(ev.imag)) > 0.5 else
                 "hyperbolisch" if np.max(np.abs(ev.imag)) < 1e-12 and np.max(np.abs(ev.real)) > 0.5 else "?")
s_E = [s for s in (+1, -1) if klasse[s] == "elliptisch"]
s_Q = [s for s in (+1, -1) if klasse[s] == "hyperbolisch"]
Ev, Qv = kand[s_E[0]], kand[s_Q[0]]
per_v = {"|exp(2 pi E) - 1|": float(np.abs(expm(2*np.pi*Ev) - E6).max()),
         "|exp(2 pi i Q^-) - 1|": float(np.abs(expm(2j*np.pi*Qv) - E6).max()),
         "|exp(4 pi i E) - 1|": float(np.linalg.norm(expm(4j*np.pi*Ev) - E6, 2)),
         "|exp(pi E) - 1|": float(np.abs(expm(np.pi*Ev) - E6).max())}
check(f"[E] In so(4,2) ist genau eine der beiden Kombinationen (P_0 +- K_0)/2 elliptisch (E, Eigenwerte +-i) und die andere "
      f"hyperbolisch (Q^-, Eigenwerte +-1): {sorted(klasse.values())}; exp(2 pi E) = 1 und exp(2 pi i Q^-) = 1 "
      f"({per_v['|exp(2 pi E) - 1|']:.0e}, {per_v['|exp(2 pi i Q^-) - 1|']:.0e}), aber exp(4 pi i E) ist nicht 1 "
      f"(Operatornorm {per_v['|exp(4 pi i E) - 1|']:.6g} = e^(4 pi) - 1; z104: 286750) und exp(pi E) auch nicht ({per_v['|exp(pi E) - 1|']:.0f}). Es gilt "
      "K_cap = 2 pi Q^- = 2 pi (P_0 - E) = 2 pi P_0 - K_Uhr/2 mit K_Uhr = 4 pi E",
      len(s_E) == 1 and len(s_Q) == 1 and per_v["|exp(2 pi E) - 1|"] < 1e-12 and per_v["|exp(2 pi i Q^-) - 1|"] < 1e-12
      and per_v["|exp(4 pi i E) - 1|"] > 1e3 and per_v["|exp(pi E) - 1|"] > 1
      and np.allclose(2*np.pi*Qv, 2*np.pi*(P0v - Ev))
      and abs(per_v["|exp(4 pi i E) - 1|"]/z104["zahlen"]["temperatur"]["Perioden"]["exp(4 pi i V_E) - 1 (beta_Original imaginär)"] - 1) < 1e-6)

Rv = sp.Symbol("R_Vak", positive=True)
l_s = Rv/2
beta_phys = 4*sp.pi*l_s                                             # rho = e^(-4 pi E), H = E/l
beta_GH_l = 2*sp.pi*l_s                                             # Gibbons-Hawking zum Horizontradius l (Kappe, P_cap)
beta_GH_R = 2*sp.pi*Rv                                              # dieselbe Formel zum Radius R_Vak
check(f"[E] Umrechnung: Der Zustand ist rho = e^(-4 pi E); die Eigenzeit am Zentrum ist tau = l theta, also H = E/l und "
      f"beta_phys = 4 pi l = {beta_phys} (hbar = c = k_B = 1), T = 1/(2 pi R_Vak). Das ist formal die Gibbons-Hawking-/"
      f"Unruh-Formel zum Radius R_Vak (beta l = 2 pi R_Vak gilt exakt: {sp.simplify(beta_phys - beta_GH_R) == 0}); der native "
      f"Horizont der Kappe hat aber den Radius l = R_Vak/2 und beta_GH = 2 pi l = {beta_GH_l}. Verhältnis beta_phys/beta_GH = "
      f"{sp.simplify(beta_phys/beta_GH_l)} (z104: {z104['zahlen']['temperatur']['beta_Original/beta_dS (Uhr)']})",
      sp.simplify(beta_phys - beta_GH_R) == 0 and sp.simplify(beta_phys/beta_GH_l) == 2
      and z104["zahlen"]["temperatur"]["beta_Original/beta_dS (Uhr)"] == "2"
      and z104["zahlen"]["temperatur"]["beta_dS Eigenzeit"] == "pi")
tolman = sp.cancel(sp.sqrt(-norm_E.subs({t: 0, y: 0, z: 0}))**2)
check(f"[E] Die Angabe gilt am Zentrum: abseits ist -g(xi_E, xi_E) = {sp.factor(tolman)} >= 1/4 (bei t = 0), die lokale "
      "inverse Temperatur also 4 pi sqrt(-g(xi, xi)) >= 2 pi R_Vak (Tolman); das Zentrum ist der heißeste Punkt, und zum Rand "
      "r -> 1 geht die lokale Temperatur gegen null",
      sp.simplify(tolman - (1 + x**2)**2/(4*(1 - x**2)**2)) == 0 and tolman.subs(x, 0) == sp.Rational(1, 4)
      and tolman.subs(x, sp.Rational(1, 2)) > sp.Rational(1, 4))
befund("E/L: beta l = 2 pi R_Vak ist eine exakte Folge zweier nativer Daten (beta = 4 pi aus dem Zustand, l = R_Vak/2 aus "
       "der Kinematik), aber keine erzwungene Gibbons-Hawking-Beziehung: Die Uhr E hat keine imaginäre Periode, und der native "
       "Horizont (Kappe, nur unter P_cap) hat Radius l, nicht R_Vak. Zwischen beiden steht der Faktor 2 ohne Mechanismus "
       "(S^1 x S^3 gegen S^4, z104). Die Übereinstimmung mit 'Horizont vom Radius R_Vak' ist ein Zufall der Normierung",
       True)

ZEITEN = [   # (Name, Faktor vor R_Vak/c, Herkunft)
    ("Markenschritt (theta = pi/2, Marken z^4 = 1)", sp.pi/4, "z104: Marken bei theta = 0, pi/2, pi, 3 pi/2"),
    ("Takt der Uhr R = e^(2 pi i L_0) (theta = 2 pi)", sp.pi, "z113/z118"),
    ("inverse Temperatur hbar/(k_B T) (beta = 4 pi)", 2*sp.pi, "z46/z58, Teil 2"),
    ("Rückkehr R^4 = 1 (theta = 8 pi)", 4*sp.pi, "L_0 in Z/4"),
]
ENERGIEN = [("Stromquant (Delta L_0 = 1)", sp.Integer(2), "Leiter 4 pi je Quant"),
            ("Lücke (L_0 = 1/4, Twistzustand)", sp.Rational(1, 2), "Jacobi-Tripelprodukt")]
check("[E] Native Zeiten sind reine Zahlen mal R_Vak/c: Markenschritt pi/4, Takt der Uhr R pi, inverse Temperatur 2 pi, "
      "Rückkehr R^4 = 1 nach 4 pi; native Energien sind Zahlen mal hbar c/R_Vak: Stromquant 2, Lücke 1/2. Konsistenz: "
      "beta mal Stromquant = 4 pi, beta mal Lücke = pi, Rückkehrzeit = 2 beta, und e^(2 pi i/4) hat genau die Ordnung 4 "
      "(Quellen: Marken z104, R = e^(2 pi i L_0) z118, L_0 = 1/4 mit Lücke pi in VERDICT.md)",
      [f for _, f, _ in ZEITEN] == [sp.pi/4, sp.pi, 2*sp.pi, 4*sp.pi]
      and sp.simplify(ZEITEN[2][1]*ENERGIEN[0][1] - 4*sp.pi) == 0 and sp.simplify(ZEITEN[2][1]*ENERGIEN[1][1] - sp.pi) == 0
      and sp.simplify(ZEITEN[3][1] - 2*ZEITEN[2][1]) == 0
      and [k for k in range(1, 9) if sp.exp(2*sp.pi*sp.I*sp.Rational(k, 4)) == 1] == [4, 8]
      and "pi/2" in z104["zahlen"]["marken"]["Marken z^4 = 1"] and hat(z118, "R = e^(2 pi i L_0)")
      and "L₀ = 1/4 (Lücke π)" in (HERE / "VERDICT.md").read_text())
R_lam = 2*mp.sqrt(32*PI**3*mp.e**(2*AINV))*lPl                       # R_Vak unter P_Lambda (COSMO.LAM.01), P_cap, MM
R_o1 = [mp.mpf(z112["zahlen"]["prämissen_G"]["Spanne R_Vak/l_Pl"][0])*lPl, mp.mpf(z112["zahlen"]["prämissen_G"]["Spanne R_Vak/l_Pl"][1])*lPl]
ZAHLEN["zeit"] = {"Faktoren vor R_Vak/c": {n: str(f) for n, f, _ in ZEITEN}, "Energien in hbar c/R_Vak": {n: str(f) for n, f, _ in ENERGIEN},
                  "unter P_Lambda (R_Vak = 2 l_dS)": {"R_Vak [m]": ns(R_lam, 5), "R_Vak/c [s]": ns(R_lam/CL, 5),
                                                      "Takt [Mrd. Jahre]": ns(PI*R_lam/CL/JAHR/1e9, 5),
                                                      "inverse Temperatur [Mrd. Jahre]": ns(2*PI*R_lam/CL/JAHR/1e9, 5),
                                                      "T [K]": ns(HBAR*CL/(2*PI*R_lam*KB), 4),
                                                      "Stromquant [eV]": ns(2*HBAR*CL/R_lam/EV, 4)},
                  "unter O(1)-Prämissen (z112, R_Vak = 0,40 bis 10,03 l_Pl)": {
                      "R_Vak [m]": [ns(v, 4) for v in R_o1], "Takt [s]": [ns(PI*v/CL, 4) for v in R_o1],
                      "Stromquant [GeV]": [ns(2*HBAR*CL/v/GEV, 4) for v in R_o1[::-1]]}}
befund("E/O: Die Theorie erzeugt keine Zeiteinheit, sondern exakte Zeitverhältnisse: jede native Zeit ist (pi/4, pi, 2 pi "
       "oder 4 pi) mal R_Vak/c. Welche SI-Zeit das ist, hängt an der einen offenen Zahl G/R_Vak^2: unter P_Lambda ist der Takt "
       f"{ns(PI*R_lam/CL/JAHR/1e9, 4)} Mrd. Jahre, unter den O(1)-Prämissen von z112 {ns(PI*R_o1[0]/CL, 2)} bis "
       f"{ns(PI*R_o1[1]/CL, 2)} s; dazwischen liegen 61 Größenordnungen", True)

# =================================================================== Teil 3
say("\nTeil 3: Tabelle der Skalenverhältnisse")
rho_ratio = 48*c3**2*mp.e**(-2*AINV)
Mscal_r = c3**(mp.mpf(7)/2)
fa_r = (c3/4)**(mp.mpf(7)/2)
l2G = 32*PI**3*mp.e**(2*AINV)
G_R2 = mp.e**(-2*AINV)/(128*PI**3)
c_lep, K_lep = {"e": mp.mpf(16)/7, "mu": mp.mpf(4)/3, "tau": mp.mpf(7)/6}, {"e": 5, "mu": 3, "tau": 2}
m_lep = {f: V_EW*PI/mp.sqrt(2)*c_lep[f]*phi0**K_lep[f] for f in ("e", "mu", "tau")}
TAB18 = {"u": (7, "0.440"), "d": (7, "0.936"), "s": (5, "0.943"), "c": (3, "0.646"), "b": (2, "0.480"), "t": (0, "0.986")}
m_q = {f: V_EW/mp.sqrt(2)*lamY**TAB18[f][0]*mp.mpf(TAB18[f][1]) for f in ("u", "d", "s", "c", "b", "t")}
H_EW = mp.log(Mbar_GeV/(V_EW/mp.sqrt(2)))
check(f"[E] gelöster Block, Leptonen (v20): m_mu/m_tau = (c_mu/c_tau) phi_0 = (8/7) phi_0 = {ns(c_lep['mu']/c_lep['tau']*phi0, 8)} "
      f"(PDG {ns(PDG['mu']/PDG['tau'], 6)}), m_e/m_mu = (12/7) phi_0^2 = {ns(c_lep['e']/c_lep['mu']*phi0**2, 8)} (PDG "
      f"{ns(PDG['e']/PDG['mu'], 6)}); Produkt der c = 32/9. Der Vergleich mit Polmassen braucht die Transfervorschrift [C]",
      abs(c_lep["mu"]/c_lep["tau"] - mp.mpf(8)/7) < mp.mpf("1e-35") and abs(c_lep["e"]/c_lep["mu"] - mp.mpf(12)/7) < mp.mpf("1e-35")
      and abs(c_lep["e"]*c_lep["mu"]*c_lep["tau"] - mp.mpf(32)/9) < mp.mpf("1e-35")
      and abs(c_lep["mu"]/c_lep["tau"]*phi0/(PDG["mu"]/PDG["tau"]) - 1) < mp.mpf("0.03")
      and abs(c_lep["e"]/c_lep["mu"]*phi0**2/(PDG["e"]/PDG["mu"]) - 1) < mp.mpf("0.03"))
check(f"[N] Quarks (v18, Lambda-Tabelle dreistellig): m_u/m_d = {ns(m_q['u']/m_q['d'], 6)} (55/117 = {ns(mp.mpf(55)/117, 6)}), "
      f"m_c/m_s = {ns(m_q['c']/m_q['s'], 5)}, m_t/m_b = {ns(m_q['t']/m_q['b'], 5)} (v18: 0,470085, 13,61, 40,80); exakte "
      "Quarkkoeffizienten und die absolute Quarkamplitude sind offen (TFPT-STAND §5)",
      abs(m_q["u"]/m_q["d"] - mp.mpf("0.470085")) < mp.mpf("1e-5") and abs(m_q["c"]/m_q["s"]/mp.mpf("13.61") - 1) < mp.mpf("2e-3")
      and abs(m_q["t"]/m_q["b"]/mp.mpf("40.80") - 1) < mp.mpf("2e-3"))
check(f"[O] Die Massen stehen in Einheiten von v_EW, nicht von Mbar: H_EW = ln(Mbar/(v_EW/sqrt2)) = {ns(H_EW, 6)} (v725: "
      f"37,1776) ist nicht abgeleitet; als Leiterstufe wäre es {ns(H_EW/mp.log(1/c3), 4)} Stufen von c_3 oder "
      f"{ns(H_EW/mp.log(1/phi0), 4)} von phi_0, nicht ganzzahlig (v725: 11,53 und 12,67). Die Teilchenmassen über Mbar sind "
      "damit je eine bekannte Zahl mal der offenen Zahl v_EW/Mbar",
      abs(H_EW - mp.mpf("37.1776")) < mp.mpf("2e-3") and abs(H_EW/mp.log(1/c3) - mp.mpf("11.53")) < mp.mpf("0.01")
      and abs(H_EW/mp.log(1/phi0) - mp.mpf("12.67")) < mp.mpf("0.01")
      and min(abs(H_EW/mp.log(1/c3) - mp.nint(H_EW/mp.log(1/c3))), abs(H_EW/mp.log(1/phi0) - mp.nint(H_EW/mp.log(1/phi0)))) > mp.mpf("0.3"))
rho_qrt_meV = rho_ratio**mp.mpf("0.25")*Mbar_GeV*mp.mpf("1e12")
H0_SI = H0_PLANCK*1000/MPC
rho_obs = OL_PLANCK*3*H0_SI**2/(8*PI*G_SI)*CL**2                       # J/m^3
rho_obs_qrt_meV = (rho_obs*(HBAR*CL)**3)**mp.mpf("0.25")/EV*1000
drho = mp.sqrt((2*DH0/H0_PLANCK)**2 + (DOL/OL_PLANCK)**2)               # ohne Korrelation, grobe 1-sigma-Breite
abw = rho_qrt_meV/rho_obs_qrt_meV - 1
Mbar_cosmo = mp.mpf("2.24e-12")/rho_ratio**mp.mpf("0.25")
ZAHLEN["lambda"] = {"rho_L/Mbar^4": ns(rho_ratio, 7), "rho_L^(1/4) [meV]": ns(rho_qrt_meV, 7),
                    "beobachtet (Planck 2018, Omega_L und H_0 hier eingetippt) [meV]": ns(rho_obs_qrt_meV, 6),
                    "1 sigma der vierten Wurzel [%]": ns(100*drho/4, 3), "Abweichung [%]": ns(100*abw, 3),
                    "Abweichung [sigma]": ns(abw/(drho/4), 3), "mit H_0 = 73,04 statt 67,36 [meV]": ns(rho_obs_qrt_meV*mp.sqrt(mp.mpf("73.04")/H0_PLANCK), 5),
                    "Mbar aus rho_L (2,24 meV) [GeV]": ns(Mbar_cosmo, 5), "gegen Mbar aus G [%]": ns(100*(Mbar_cosmo/Mbar_GeV - 1), 3),
                    "l^2/G": ns(l2G, 6), "l_dS/l_Pl": ns(mp.sqrt(l2G), 6), "G/R_Vak^2": ns(G_R2, 6), "log10 G/R_Vak^2": ns(mp.log10(G_R2), 7)}
check(f"[N] COSMO.LAM.01 nachgerechnet: rho_L/Mbar^4 = 48 c_3^2 e^(-2/alpha) = {ns(rho_ratio, 6)} (v7: 7,12533e-121), "
      f"rho_L^(1/4) = {ns(rho_qrt_meV, 6)} meV; aus Planck-Werten (Omega_L = 0,6847, H_0 = 67,36) folgt {ns(rho_obs_qrt_meV, 5)} meV "
      f"mit etwa {ns(100*drho/4, 2)} % Breite: Abweichung {ns(100*abw, 2)} % ({ns(abs(abw)/(drho/4), 2)} sigma). v274: Mbar aus "
      f"2,24 meV = {ns(Mbar_cosmo, 5)} GeV, {ns(100*(Mbar_cosmo/Mbar_GeV - 1), 2)} % über dem Wert aus G. Die 0,11 % liegen "
      "innerhalb der Messbreite; mit H_0 = 73,04 wären es rund 4 %",
      abs(rho_ratio/mp.mpf("7.12533e-121") - 1) < mp.mpf("1e-5") and abs(rho_qrt_meV - mp.mpf("2.23747")) < mp.mpf("1e-4")
      and abs(rho_obs_qrt_meV - mp.mpf("2.24")) < mp.mpf("0.006") and abs(abw) < drho/4
      and abs(Mbar_cosmo/Mbar_GeV - 1 - mp.mpf("0.0011")) < mp.mpf("2e-4"))
check(f"[E] COSMO.LAM.01 in H6-Variablen (z96/z98): l_dS^2/G = 32 pi^3 e^(2/alpha) = {ns(l2G, 5)}, l_dS = {ns(mp.sqrt(l2G), 5)} "
      f"l_Pl; mit R_Vak = 2 l: G/R_Vak^2 = e^(-2/alpha)/(128 pi^3) = 10^{ns(mp.log10(G_R2), 6)} (z112: "
      f"{z112['zahlen']['prämissen_G']['Ziel log10 G/R_Vak^2']}); |rho| = 12/(8 pi G R_Vak^2) gibt dasselbe rho_L/Mbar^4 "
      "(Betrag; das native Vorzeichen ist negativ)",
      abs(mp.log10(G_R2) - mp.mpf(z112["zahlen"]["prämissen_G"]["Ziel log10 G/R_Vak^2"])) < mp.mpf("1e-3")
      and abs(96*PI*G_R2/rho_ratio - 1) < mp.mpf("1e-30") and abs(mp.sqrt(l2G)/mp.mpf("1.03e61") - 1) < mp.mpf("0.01"))

# Anti-Numerologie: Wie viel trägt der Vorfaktor? Abgetastete Familie n c_3^k, n = 1..64, k = 0..6
A_obs = rho_obs/GEV**4*(HBAR*CL)**3/Mbar_GeV**4*mp.e**(2*AINV)          # beobachteter Vorfaktor vor e^(-2/alpha)
fam = sorted(((n*c3**k, n, k) for n in range(1, 65) for k in range(0, 7)), key=lambda v: (v[0], v[1], v[2]))
tr1 = [(n, k) for a, n, k in fam if abs(a/A_obs - 1) < drho]
tr2 = [(n, k) for a, n, k in fam if abs(a/A_obs - 1) < 2*drho]
dichte = sum(1 for a, _, _ in fam if abs(mp.log(a/A_obs)) < mp.log(10))/(2*mp.log(10))
erw1, erw2 = dichte*2*drho, dichte*4*drho
k_eff = (mp.log(Mbar_GeV**4/(rho_obs/GEV**4*(HBAR*CL)**3)))/AINV
ZAHLEN["anti_numerologie_vorfaktor"] = {"Familie": "n c_3^k, n = 1..64, k = 0..6", "Größe": len(fam),
                                        "beobachteter Vorfaktor A": ns(A_obs, 5), "48 c_3^2": ns(48*c3**2, 5),
                                        "Treffer in 1 sigma (n, k)": tr1, "Treffer in 2 sigma (n, k)": tr2,
                                        "lokale Dichte je ln-Einheit": ns(dichte, 4),
                                        "erwartete Zufallstreffer 1 sigma, 2 sigma": [ns(erw1, 3), ns(erw2, 3)],
                                        "ln(Mbar^4/rho_obs)/alpha^-1": ns(k_eff, 6)}
check(f"[N] Anti-Numerologie am Bestand: Der beobachtete Vorfaktor vor e^(-2/alpha) ist A = {ns(A_obs, 4)}; 48 c_3^2 = "
      f"{ns(48*c3**2, 4)} liegt im 1-sigma-Fenster (+-{ns(100*drho, 2)} % in rho). In der Familie n c_3^k (n <= 64, k <= 6, "
      f"{len(fam)} Ausdrücke) liegen {len(tr1)} im 1-sigma- und {len(tr2)} im 2-sigma-Fenster, (n, k) = {tr2}; bei der lokalen Dichte "
      f"{ns(dichte, 3)} je ln-Einheit sind {ns(erw1, 2)} bzw. {ns(erw2, 2)} Zufallstreffer zu erwarten. Die Messung legt "
      f"den Vorfaktor nur auf etwa 2 % fest (n = 48 +- 1); er allein ist nicht signifikant. Tragend ist der Exponent: "
      f"ln(Mbar^4/rho_obs)/alpha^-1 = {ns(k_eff, 5)}",
      (48, 2) in tr1 and len(fam) == 448 and erw2 > mp.mpf("0.3") and abs(k_eff - 2) < mp.mpf("0.03"))

# Tabelle: (Verhältnis, Formel, Wert, Quelle, Status)
VERH = [
    ("alpha^-1", "Wurzel von F_U(1)", AINV, "v3", "M: Fixpunktgleichung; Lesart als Elektromagnetismus C"),
    ("l/R_Vak", "1/2", mp.mpf(1)/2, "z98", "M (Kinematik des Kandidaten H6)"),
    ("|Lambda| R_Vak^2", "12", mp.mpf(12), "z98", "M (H6); Vorzeichen nativ negativ (AdS)"),
    ("c beta_phys/R_Vak", "2 pi", 2*PI, "z140", "M (H6): beta = 4 pi und l = R_Vak/2"),
    ("c Takt/R_Vak", "pi", PI, "z140", "M (H6)"),
    ("Energiequant R_Vak/(hbar c)", "2", mp.mpf(2), "z140", "M (H6)"),
    ("m_mu/m_tau", "(8/7) phi_0", c_lep["mu"]/c_lep["tau"]*phi0, "v20", "M innerhalb der Annahmen; Polvergleich C (+2,2 %)"),
    ("m_e/m_mu", "(12/7) phi_0^2", c_lep["e"]/c_lep["mu"]*phi0**2, "v20", "M innerhalb der Annahmen; Polvergleich C (+0,2 %)"),
    ("m_u/m_d", "55/117", mp.mpf(55)/117, "v18/v49", "M (Quellverhältnis); Vergleich C"),
    ("m_c/m_s", "lambda_Y^-2 L_c/L_s", m_q["c"]/m_q["s"], "v18", "N (Tabelle dreistellig)"),
    ("m_t/m_b", "lambda_Y^-2 L_t/L_b", m_q["t"]/m_q["b"], "v18", "N (Tabelle dreistellig)"),
    ("m_tau/v_EW", "(pi/sqrt2)(7/6) phi_0^2", m_lep["tau"]/V_EW, "v20/v46", "M innerhalb der Annahmen (Quellmasse)"),
    ("m_t/v_EW", "0,986/sqrt2", m_q["t"]/V_EW, "v18", "N"),
    ("v_EW/Mbar (also jede Masse über Mbar)", "nicht abgeleitet", V_EW/Mbar_GeV, "v725, Ledger ANCHOR.VGEO.01", "O (Wert ist Messung)"),
    ("M_Pl/v_EW", "nicht abgeleitet", MPl_GeV/V_EW, "Ledger", "O (Wert ist Messung)"),
    ("M_scal/Mbar", "c_3^(7/2)", Mscal_r, "v7 GRAV.SCAL.01", "N; ungemessen (Vorhersage)"),
    ("f_a/Mbar", "(c_3/4)^(7/2)", fa_r, "v25", "Vermutung; ungemessen"),
    ("rho_Lambda/Mbar^4", "48 c_3^2 e^(-2/alpha)", rho_ratio, "v7/v60 COSMO.LAM.01", "N, bedingt auf Lambda-Zweig; kein Mechanismus für e^(-2/alpha)"),
    ("R_Vak/l_Pl", "sqrt(128 pi^3) e^(1/alpha)", 1/mp.sqrt(G_R2), "z96/z98", "C: P_Lambda (COSMO.LAM.01), P_cap (Vorzeichen), MM-Wirkung; nativ O"),
    ("R_Vak/l_Pl (Alternative)", "O(1)-Prämissen P_L", mp.mpf(z112["zahlen"]["prämissen_G"]["Spanne R_Vak/l_Pl"][1]), "z112", "C: 0,40 bis 10,03 (Obergrenze gezeigt)"),
    ("Energiequant/Mbar", "2 sqrt(8 pi) l_Pl/R_Vak", 2*mp.sqrt(8*PI)*mp.sqrt(G_R2), "z140", "C wie R_Vak/l_Pl"),
    ("Lambda_QCD/v_EW", "b_3 = -7, alpha_s(M_Z) gemessen", mp.mpf("0.5")/V_EW, "v164 QCD.LAMBDA.01", "O: braucht alpha_s(M_Z) als Eingabe (mu_conf ~ 0,5 GeV)"),
    ("M_R/M_scal", "Kandidat 3", mp.mpf(3), "FLAV.NUSCALE.05", "Kandidat unter y_nu = y_t; O"),
    ("T_reh/M_scal", "importierte Reheating-Physik", mp.mpf("9.55e9")/(Mscal_r*Mbar_GeV), "COSMO.NSTAR.01", "P (importiert); O"),
    ("M_f/M_g", "1/sqrt(e^(4 pi) - 1)", 1/mp.sqrt(mp.e**(4*PI) - 1), "z48", "L (bimetrische Lesart)"),
    ("m a (Nahtkette)", "4 e^(-pi/2)", ddv_nat, "z121", "L; O(1), keine Hierarchie"),
]
say(f"   {'Verhältnis':40s} {'Formel':30s} {'Wert':>14s}  Quelle; Status")
for nm, fo, val, qu, st in VERH:
    say(f"   {nm:40s} {fo:30s} {ns(val, 7):>14s}  {qu}; {st}")
ZAHLEN["verhaeltnisse"] = [{"Verhältnis": nm, "Formel": fo, "Wert": ns(val, 10), "Quelle": qu, "Status": st}
                           for nm, fo, val, qu, st in VERH]
n_M = sum(1 for v in VERH if v[4].startswith("M"))
n_N = sum(1 for v in VERH if v[4].startswith("N"))
n_C = sum(1 for v in VERH if v[4].startswith(("C", "L", "Vermutung", "Kandidat", "P")))
n_O = sum(1 for v in VERH if v[4].startswith("O"))
check(f"[E] Tabelle vollständig typisiert: {len(VERH)} Zeilen, davon {n_M} durch einen Mechanismus (innerhalb der jeweiligen "
      f"Annahmen), {n_N} numerische Treffer, {n_C} bedingt/Lesart/Kandidat, {n_O} offen; M_scal = {ns(Mscal_r*Mbar_GeV, 5)} GeV "
      "(v7: 3,06e13)", n_M + n_N + n_C + n_O == len(VERH) and n_M == 10 and n_O == 3
      and abs(Mscal_r*Mbar_GeV/mp.mpf("3.06e13") - 1) < mp.mpf("5e-3"))


def si(E_GeV):
    """Masse [kg], Compton-Länge [m], Zeit [s] zu einer Energie in GeV."""
    m = E_GeV*GEV/CL**2
    return m, HBAR/(m*CL), HBAR/(m*CL**2)


l_dS = mp.sqrt(l2G)*lPl
SI = []


def zeile(name, wert, einheit, grundlage):
    SI.append({"Größe": name, "Wert": ns(wert, 5), "Einheit": einheit, "Grundlage": grundlage})
    say(f"   {name:46s} {ns(wert, 5):>12s} {einheit:8s} {grundlage}")


say("   SI-Tabelle mit dem einen Anker CODATA G:")
zeile("Planck-Länge l_Pl", lPl, "m", "Definition des Ankers, keine Vorhersage")
zeile("Planck-Zeit t_Pl", tPl, "s", "Definition des Ankers, keine Vorhersage")
zeile("Mbar", Mbar_kg, "kg", "Definition des Ankers, keine Vorhersage")
zeile("M_scal (Masse)", si(Mscal_r*Mbar_GeV)[0], "kg", "N: c_3^(7/2), ungemessen")
zeile("M_scal (Compton-Länge)", si(Mscal_r*Mbar_GeV)[1], "m", "N")
zeile("M_scal (Zeit)", si(Mscal_r*Mbar_GeV)[2], "s", "N")
zeile("f_a (Masse)", si(fa_r*Mbar_GeV)[0], "kg", "Vermutung")
zeile("rho_Lambda (Massendichte)", rho_ratio*Mbar_GeV**4*GEV**4/(HBAR*CL)**3/CL**2, "kg/m^3", "NUR numerischer Treffer (COSMO.LAM.01)")
zeile("Länge zu rho_L^(1/4)", HBAR*CL/(rho_qrt_meV*EV/1000), "m", "NUR numerischer Treffer")
zeile("l_dS = sqrt(3/Lambda)", l_dS, "m", "NUR numerischer Treffer (dieselbe Aussage)")
zeile("l_dS in Mrd. Lichtjahren", l_dS/LJ/1e9, "Gly", "NUR numerischer Treffer")
zeile("H6: R_Vak = 2 l_dS", R_lam, "m", "bedingt: P_Lambda, P_cap, MM-Wirkung")
zeile("H6: Takt der Uhr R = pi R_Vak/c", PI*R_lam/CL, "s", "bedingt wie R_Vak")
zeile("H6: inverse Temperatur 2 pi R_Vak/c", 2*PI*R_lam/CL, "s", "bedingt wie R_Vak")
zeile("H6: T = hbar c/(2 pi k_B R_Vak)", HBAR*CL/(2*PI*R_lam*KB), "K", "bedingt wie R_Vak")
zeile("H6: Stromquant 2 hbar c/R_Vak", 2*HBAR*CL/R_lam/EV, "eV", "bedingt wie R_Vak")
zeile("H6 Alternative: R_Vak (z112, Obergrenze)", R_o1[1], "m", "bedingt: O(1)-Prämissen P_L")
zeile("H6 Alternative: Takt (z112, Obergrenze)", PI*R_o1[1]/CL, "s", "bedingt: O(1)-Prämissen P_L")
for f in ("e", "mu", "tau"):
    zeile(f"m_{f} (Quellmasse)", si(m_lep[f])[0], "kg", "braucht ZWEITEN Messwert v_EW; Verhältnisse M")
for f in ("u", "d", "s", "c", "b", "t"):
    zeile(f"m_{f} (Quellmasse)", si(m_q[f])[0], "kg", "braucht v_EW; N (Tabelle dreistellig)")
zeile("Lambda_QCD ~ mu_conf", si(mp.mpf("0.5"))[0], "kg", "braucht alpha_s(M_Z) und v_EW-Skala (v164)")
ZAHLEN["si_tabelle"] = SI
ZAHLEN["massen_GeV"] = {**{f: ns(m_lep[f], 6) for f in ("e", "mu", "tau")}, **{f: ns(m_q[f], 6) for f in ("u", "d", "s", "c", "b", "t")}}
check(f"[E/N] SI-Umrechnung mit einem Anker: l_dS = {ns(l_dS, 5)} m = {ns(l_dS/LJ/1e9, 4)} Mrd. Lichtjahre (z96: 17,6), M_scal = "
      f"{ns(si(Mscal_r*Mbar_GeV)[0], 4)} kg, Länge zu rho_L^(1/4) = {ns(1e6*HBAR*CL/(rho_qrt_meV*EV/1000), 4)} Mikrometer. "
      f"Die neun Teilchenmassen lassen sich mit G allein nicht umrechnen (v_EW/Mbar offen); mit dem Messwert v_EW sind die "
      f"Quellmassen m_tau = {ns(m_lep['tau'], 5)}, m_mu = {ns(m_lep['mu'], 5)}, m_e = {ns(m_lep['e']*1000, 5)}e-3 GeV "
      f"({ns(100*(m_lep['tau']/PDG['tau'] - 1), 2)}, {ns(100*(m_lep['mu']/PDG['mu'] - 1), 2)}, {ns(100*(m_lep['e']/PDG['e'] - 1), 2)} % "
      "über den Polmassen)",
      abs(l_dS/LJ/1e9 - mp.mpf("17.6")) < mp.mpf("0.05") and all(abs(m_lep[f]/PDG[f] - 1) < mp.mpf("0.05") for f in ("e", "mu", "tau"))
      and len(SI) == 28)

# =================================================================== Teil 4
say("\nTeil 4: Gibt es einen nativen Mechanismus für die eine große Zahl?")
ziel = mp.e**(-2*AINV)
say(f"   Ziel: G/R_Vak^2 = e^(-2/alpha)/(128 pi^3); e^(-2/alpha) = 10^{ns(mp.log10(ziel), 8)}")
# --- (a) Boltzmann/KMS
E_stern, E_kappe = AINV/(2*PI), AINV/PI
nachbar = [mp.floor(4*E_stern)/4, mp.ceil(4*E_stern)/4]
lage = (E_stern - nachbar[0])/(nachbar[1] - nachbar[0])
gew = [r_nat**e for e in nachbar]
ZAHLEN["mechanismus_a"] = {"E* = alpha^-1/(2 pi) (beta = 4 pi)": ns(E_stern, 12), "alpha^-1/pi (beta = 2 pi)": ns(E_kappe, 12),
                           "nächste Niveaus in Z/4": [ns(e, 6) for e in nachbar], "Lage im Niveauabstand": ns(lage, 5),
                           "log10 r^E der Nachbarn": [ns(mp.log10(g), 7) for g in gew], "log10 Ziel": ns(mp.log10(ziel), 7),
                           "Ziel/Nachbar": [ns(ziel/g, 4) for g in gew]}
check(f"[E] (a) Als Boltzmann-Gewicht: e^(-2/alpha) = r^E* mit r = e^(-4 pi) und E* = alpha^-1/(2 pi) = {ns(E_stern, 9)}. "
      f"Die nativen Niveaus liegen in Z/4 (Stromquant 1, Lücke 1/4); E* liegt bei {ns(lage, 3)} des Niveauabstands zwischen "
      f"{ns(nachbar[0], 5)} und {ns(nachbar[1], 5)}, die nächsten nativen Gewichte sind 10^{ns(mp.log10(gew[0]), 6)} und "
      f"10^{ns(mp.log10(gew[1]), 6)} (Faktor {ns(gew[0]/ziel, 3)} und {ns(ziel/gew[1], 3)} daneben). Mit der Kappe (beta = 2 pi) "
      f"wäre E = alpha^-1/pi = {ns(E_kappe, 7)}, ebenfalls kein Niveau",
      abs(r_nat**E_stern/ziel - 1) < mp.mpf("1e-30") and mp.mpf("0.05") < lage < mp.mpf("0.95")
      and gew[0]/ziel > 2 and ziel/gew[1] > 2 and abs(4*E_kappe - mp.nint(4*E_kappe)) > mp.mpf("0.05"))
a_s = sp.Symbol("a", positive=True)
c3_s = 1/(8*sp.pi)
Q_s = 48*c3_s**4*sp.exp(-2*a_s)
F_s = a_s**3 - 2*c3_s**3*a_s**2 - sp.Rational(4, 5)*c3_s**6*41*sp.log(1/(1/(6*sp.pi) + Q_s*(1 - Q_s)**sp.Rational(-5, 4)))
exp_args = sorted(str(e.args[0]) for e in F_s.atoms(sp.exp))
lam_s = 48*c3_s**2*sp.exp(-2/a_s)
dual = sp.simplify((8*sp.pi)**2*Q_s.subs(a_s, 1/a_s) - lam_s) == 0
ainv_ohne = 1/mp.findroot(lambda a: F_alpha(a, rueck=lambda a: mp.mpf(1)), mp.mpf("0.0073"))
ZAHLEN["mechanismus_a"]["Exponentialargumente in F_U(1)"] = exp_args
ZAHLEN["mechanismus_a"]["alpha^-1 ohne die Rückwirkung e^(-2 alpha)"] = ns(ainv_ohne, 12)
check(f"[E] (a) In der Fixpunktgleichung kommt 2/alpha nicht vor: F_U(1) enthält genau die Exponentialargumente {exp_args} "
      "(die Nahtrückwirkung Q(alpha) = delta_top e^(-2 alpha), bei alpha = 1/137 nahe 1; ohne sie wäre alpha^-1 = "
      f"{ns(ainv_ohne, 9)}). 2 alpha^-1 tritt dort weder als Wirkung noch als beta E auf. Exakte Identität: rho_L G^2 = "
      f"delta_top e^(-2/alpha) = Q(1/alpha), dieselbe Funktion am inversen Argument ({dual}). Kontrolle: die Lambda-Formel "
      f"selbst hat das Argument {sorted(str(e.args[0]) for e in lam_s.atoms(sp.exp))}",
      exp_args == ["-2*a"] and dual and sorted(str(e.args[0]) for e in lam_s.atoms(sp.exp)) == ["-2/a"]
      and abs(ainv_ohne - AINV) > mp.mpf("1e-4"))
befund("(a) Boltzmann/KMS: Mechanismus vorhanden (r = e^(-4 pi) ist nativ), aber er gibt andere Zahlen: r^E mit E in Z/4, "
       "nächstgelegen 10^-118,70 und 10^-120,07 statt 10^-119,03, und nichts zeichnet ein Niveau bei 22 aus. "
       "alpha^-1/(2 pi) ist keine native Energie. Die Identität rho_L G^2 = Q(1/alpha) ist eine Schreibweise des "
       "Ledgereintrags, kein Mechanismus: in den gelesenen Quellen wertet kein Operator Q am inversen Argument aus", True)

# --- (b) dimensionale Transmutation
def gd(u):
    return 2*mp.atan(mp.e**u) - PI/2


def kern(lam):
    f = lambda k: mp.e**(-k/2)/(1 + mp.e**(-k))*mp.cos(k*lam)
    return mp.quad(f, [0, 10, 40, 160])/PI                             # (1/2 pi) int_-oo^oo ... dk


kern_abw = max(abs(kern(mp.mpf(l)) - 1/(2*mp.cosh(PI*mp.mpf(l)))) for l in ("0", "0.3", "1.1"))
kern_kontr = abs(kern(mp.mpf("0.3")) - 1/(2*mp.cosh(mp.mpf("0.3"))))


def E_ddv(lam, Th):
    return gd(PI*(lam + Th)) - gd(PI*(lam - Th))


def P_ddv(lam, Th):
    return gd(PI*(lam + Th)) + gd(PI*(lam - Th))


def v_ddv(Th, lam=mp.mpf(6)):
    return -mp.diff(lambda u: E_ddv(u, Th), lam)/mp.diff(lambda u: P_ddv(u, Th), lam)


v_z121 = mp.mpf(z121["zahlen"]["kontinuum"]["v(tau) = 2 tanh(pi tau/2), tau = 1"])
masse_ex = {Th: PI - E_ddv(mp.mpf(0), mp.mpf(Th)) for Th in ("0.5", "3")}
masse_as = {Th: 4*mp.e**(-PI*mp.mpf(Th)) for Th in ("0.5", "3")}
Th_noetig = (AINV + mp.log(4))/PI
ZAHLEN["mechanismus_b"] = {"Kernabweichung (Fourier gegen 1/(2 cosh pi lambda))": ns(kern_abw, 3),
                           "v(Theta = 1/2) je Platzpaar": ns(2*v_ddv(mp.mpf("0.5")), 10), "z121": ns(v_z121, 10),
                           "pi - E(0) exakt (Theta = 1/2, 3)": [ns(masse_ex[k], 8) for k in ("0.5", "3")],
                           "4 e^(-pi Theta)": [ns(masse_as[k], 8) for k in ("0.5", "3")],
                           "nötiges Theta für m a = e^(-1/alpha)": ns(Th_noetig, 8), "alpha^-1/pi": ns(AINV/PI, 8),
                           "natives Theta": "1/2", "nötige Kopplung g = 1/Theta (b = h-dual = 2)": ns(1/Th_noetig, 6),
                           "nötige Kopplung bei b = 4 (z99: sum f f = 4)": ns(1/(2*Th_noetig), 6), "native g": "2",
                           "pi alpha": ns(PI/AINV, 6)}
check(f"[E/L] (b) Lichtkegelgitter (z121: R(u) = u + i P, Inhomogenitäten -+Theta): Der gekleidete Kern ist 1/(2 cosh pi lambda) "
      f"(Fourier-Quadratur, Abweichung {ns(kern_abw, 2)}; Kontrolle mit cosh(lambda) statt cosh(pi lambda): {ns(kern_kontr, 2)}), "
      "also E a = gd(pi(lambda + Theta)) - gd(pi(lambda - Theta)). Daraus v = tanh(pi Theta) je Platz; bei Theta = 1/2 je "
      f"Platzpaar {ns(2*v_ddv(mp.mpf('0.5')), 8)} = z121 ({ns(v_z121, 8)}), Kontrolle Theta = 1: {ns(2*v_ddv(mp.mpf(1)), 6)}. "
      f"Massenskala pi - E(0) = {ns(masse_ex['0.5'], 5)} bei Theta = 1/2 (Asymptotik 4 e^(-pi Theta) = {ns(masse_as['0.5'], 5)}), "
      f"bei Theta = 3 {ns(masse_ex['3'], 6)} gegen {ns(masse_as['3'], 6)}",
      kern_abw < mp.mpf("1e-12") and kern_kontr > mp.mpf("0.05") and abs(2*v_ddv(mp.mpf("0.5")) - v_z121) < mp.mpf("1e-9")
      and abs(2*v_ddv(mp.mpf(1)) - v_z121) > mp.mpf("0.1") and abs(masse_ex["3"]/masse_as["3"] - 1) < mp.mpf("1e-6")
      and abs(masse_ex["0.5"]/masse_as["0.5"] - 1) > mp.mpf("0.01") and abs(masse_as["0.5"] - ddv_nat) < mp.mpf("1e-30"))
b_nativ = {"b_1 = 41/10": mp.mpf(41)/10, "|b_2| = 19/6": mp.mpf(19)/6, "|b_3| = 7": mp.mpf(7), "2 h-dual su(2) = 4 (z99)": mp.mpf(4),
           "h-dual su(2) = 2": mp.mpf(2)}
exp_b = {k: 2*PI*AINV/v for k, v in b_nativ.items()}
ZAHLEN["mechanismus_b"]["Einschleifenexponent 2 pi/(b alpha) für native b"] = {k: ns(v, 7) for k, v in exp_b.items()}
ZAHLEN["mechanismus_b"]["nötig"] = {"b für e^(-2/alpha)": "pi", "b für e^(-1/alpha)": "2 pi", "2/alpha": ns(2*AINV, 7)}
naechst = min(exp_b.items(), key=lambda kv: (abs(kv[1] - 2*AINV), kv[0]))
check(f"[E] (b) Native Zahl: m a = 4 e^(-pi/2) = {ns(ddv_nat, 6)}, also O(1) und keine Hierarchie (z121: der Vakuumsektor ist "
      f"bei Theta = 1/2 kritisch, kein Parameter schickt m a gegen null). Für m a = e^(-1/alpha) wäre Theta = (alpha^-1 + ln 4)/pi "
      f"= {ns(Th_noetig, 6)} nötig (nativ 1/2), als Kopplung g = 1/Theta = {ns(1/Th_noetig, 4)} nahe pi alpha = {ns(PI/AINV, 4)} "
      f"(nativ 2). Mit der Kopplung alpha selbst verlangt exp(-2 pi/(b alpha)) = e^(-2/alpha) den Koeffizienten b = pi (für e^(-1/alpha): b = 2 pi); die "
      f"nativen b sind rational ({len(b_nativ)} geprüft). Am nächsten liegt {naechst[0]} mit Exponent {ns(naechst[1], 6)} gegen "
      f"{ns(2*AINV, 6)} (19/6 liegt 0,8 % neben pi), also e^(-2/alpha) um den Faktor {ns(mp.e**(2*AINV - naechst[1]), 3)} verfehlt",
      mp.mpf("0.5") < ddv_nat < 1 and Th_noetig > 40 and abs(1/Th_noetig/(PI/AINV) - 1) < mp.mpf("0.02")
      and naechst[0] == "|b_2| = 19/6" and abs(mp.e**(2*AINV - naechst[1]) - 1) > mp.mpf("0.5")
      and all(abs(v - 2*AINV) > 1 for v in exp_b.values()))
befund("(b) Dimensionale Transmutation: Mechanismus vorhanden (die native Kette hat die Destri-de-Vega-Formel und eine "
       "laufende Strom-Strom-Kopplung), aber er gibt eine andere Zahl: m a = 0,83 am nativen Punkt Theta = 1/2. Für e^(-1/alpha) "
       "bräuchte es Theta = 44,06 bzw. g nahe pi alpha; weder Theta noch g sind nativ mit alpha verbunden (z99: die Kopplung ist "
       "gewählt; z121: Theta = 1/2 kommt aus P_viertel). b = pi ist für rationale native Koeffizienten ausgeschlossen", True)

# --- (c) Instanton
gew16 = [w for w in itertools.product((sp.Rational(1, 2), -sp.Rational(1, 2)), repeat=5) if sum(1 for v in w if v < 0) % 2 == 0]
t3_16 = [(w[0] + w[1])/2 for w in gew16]                                # Wurzel-SU(2) zur Wurzel e_1 + e_2
gew10 = [tuple(sgn*(1 if i == j else 0) for i in range(5)) for j in range(5) for sgn in (1, -1)]
t3_10 = [sp.Rational(w[0] + w[1], 2) for w in gew10]
nm16, nm10 = 2*sum(v**2 for v in t3_16), 2*sum(v**2 for v in t3_10)     # Nullmoden = 2 T(R), T = Summe T_3^2
S_inst_em = 2*PI*AINV
unterdr = {k: mp.e**(-2*PI*mp.mpf(k)) for k in sorted(z33["zahlen"]["unterdrueckung"], key=float)}
ok_z33 = all(abs(mp.log(unterdr[k]) - mp.log(mp.mpf(z33["zahlen"]["unterdrueckung"][k]))) < mp.mpf("1e-6") for k in unterdr)
kopp = {"induziert stark (z33)": mp.mpf(4), "z33, Wert 23": mp.mpf(23), "(3/8) alpha^-1 (sin^2 theta_W = 3/8, ohne Laufen)": mp.mpf(3)/8*AINV,
        "nackte SO(10) (z25)": mp.mpf(60), "alpha^-1 selbst": AINV}
ZAHLEN["mechanismus_c"] = {"48 c_3^2 = 3/(4 pi^2) = (8 pi)^2 delta_top": True, "Nullmoden 16": str(nm16), "Nullmoden 10": str(nm10),
                           "drei Familien": str(3*nm16), "z33": z33["zahlen"]["nullmoden"],
                           "Instantonwirkung 2 pi/alpha_g, nötig alpha_g^-1 = alpha^-1/pi": ns(AINV/PI, 8),
                           "log10 e^(-2 pi alpha_g^-1) für vorhandene Kopplungen": {k: ns(-2*PI*v/mp.log(10), 6) for k, v in kopp.items()},
                           "log10 Ziel": ns(mp.log10(ziel), 6), "2/alpha in Einheiten der Instantonwirkung 2 pi/alpha": ns(1/PI, 8)}
check("[E] (c) Vorfaktor: 3/(4 pi^2) = 48 c_3^2 = (8 pi)^2 delta_top exakt (sympy), und 48 = 3 x 16 ist die Zahl der "
      f"Weyl-Zustände. Die Nullmodenzahl eines Einheitsinstantons in einer Wurzel-SU(2) von Spin(10) ist aber 2 T(R): {nm16} "
      f"je 16 (vier Dubletts, acht Singuletts), {nm10} für die 10, {3*nm16} für drei Familien (z33: {z33['zahlen']['nullmoden']['drei Familien']} "
      f"plus {z33['zahlen']['nullmoden']['Eltern (psi, chi je (10,1)+(1,6))']} Eltern). 48 ist keine Nullmodenzahl",
      sp.simplify(sp.Rational(3, 4)/sp.pi**2 - 48*c3_s**2) == 0 and sp.simplify((8*sp.pi)**2*48*c3_s**4 - 48*c3_s**2) == 0
      and len(gew16) == 16 and nm16 == 4 and nm10 == 2 and 3*nm16 == z33["zahlen"]["nullmoden"]["drei Familien"]
      and z33["zahlen"]["nullmoden"]["Spinor 16"] == 4 and 3*nm16 != 48)
check(f"[E] (c) Wirkung: Ein Einheitsinstanton hat S = 8 pi^2/g^2 = 2 pi/alpha_g (z33: e^(-2 pi 60) = {ns(unterdr['60.0'], 3)}, "
      f"nachgerechnet: {ok_z33}). 2/alpha = 8 pi/e^2 = 1/(c_3 e^2) ist das 1/pi-fache der Instantonwirkung zur Kopplung alpha, "
      f"kein ganzzahliges oder 1/N-Vielfaches. Nötig wäre alpha_g^-1 = alpha^-1/pi = {ns(AINV/PI, 7)}; die im Bestand vorhandenen "
      f"Kopplungen ({len(kopp)}) geben log10 = "
      f"{[ns(-2*PI*v/mp.log(10), 5) for v in kopp.values()]} gegen {ns(mp.log10(ziel), 6)}",
      ok_z33 and abs(2*AINV/S_inst_em - 1/PI) < mp.mpf("1e-30") and abs(2*AINV - 1/(c3*4*PI/AINV)) < mp.mpf("1e-30")
      and all(abs(-2*PI*v/mp.log(10) - mp.log10(ziel)) > 15 for v in kopp.values())
      and all(abs(1/PI - mp.mpf(1)/N) > mp.mpf("0.01") for N in range(1, 13)))
befund("(c) Instanton/Tunneln: kein solcher Mechanismus. Der Vorfaktor 48 c_3^2 ist eine exakte Identität, aber 48 zählt "
       "Weyl-Zustände, nicht Nullmoden (12 bzw. 16). Kein Instanton des Bestands hat die Wirkung 2/alpha = 1/(c_3 e^2): die "
       "Instantonwirkung zur Kopplung alpha ist 2 pi/alpha (10^-374), und alpha^-1/pi = 43,62 ist keine vorhandene Kopplung", True)

# --- (d) Zählung
S_dS = PI*l2G
S_einheit = {"ln 2/2 (Quarterquelle je Naht)": mp.log(2)/2, "ln 3 (Vielfachheit Isospin 1)": mp.log(3), "ln 12 (E8-Naht)": mp.log(12)}
N_noetig = {k: S_dS/v for k, v in S_einheit.items()}
d_werte = [2, 3, 4, 12, 16]
n_plaetze = {d: 2*AINV/mp.log(d) for d in d_werte}
abst = {d: abs(n_plaetze[d] - mp.nint(n_plaetze[d])) for d in d_werte}
ZAHLEN["mechanismus_d"] = {"S_dS = pi l^2/G = 32 pi^4 e^(2/alpha)": ns(S_dS, 6), "log10": ns(mp.log10(S_dS), 7),
                           "nötige Nahtzahl je Entropieeinheit": {k: ns(v, 5) for k, v in N_noetig.items()},
                           "n mit d^n = e^(2/alpha)": {str(d): ns(n_plaetze[d], 7) for d in d_werte},
                           "Abstand zur nächsten ganzen Zahl": {str(d): ns(abst[d], 3) for d in d_werte},
                           "Abtastgröße": len(d_werte), "erwartete Zufallstreffer bei Fenster 0,05": ns(len(d_werte)*mp.mpf("0.1"), 2)}
check(f"[E] (d) Zählung: S_dS = pi l^2/G = 32 pi^4 e^(2/alpha) = 10^{ns(mp.log10(S_dS), 6)} (z112: "
      f"{z112['zahlen']['prämissen_G']['de-Sitter-Lesart pi l^2/G']['Ziel log10 pi l^2/G']}); als Nahtzahl wären das "
      f"{ns(N_noetig['ln 2/2 (Quarterquelle je Naht)'], 3)} Nähte zu ln 2/2. Als Zustandszahl d^n = e^(2/alpha) für d in {d_werte} "
      f"(Abtastgröße {len(d_werte)}): n = {[ns(n_plaetze[d], 6) for d in d_werte]}, keine ganze Zahl innerhalb 0,05 "
      f"(erwartete Zufallstreffer {ns(len(d_werte)*mp.mpf('0.1'), 2)}). Die Nahtzahl N ist nativ frei (z119: K_H(S) für jedes N)",
      abs(mp.log10(S_dS) - mp.mpf(z112["zahlen"]["prämissen_G"]["de-Sitter-Lesart pi l^2/G"]["Ziel log10 pi l^2/G"])) < mp.mpf("1e-3")
      and all(abst[d] > mp.mpf("0.05") for d in d_werte) and hat(z119, "N Nähte mit je einer Anregung"))
befund("(d) Zählung: kein solcher Mechanismus. l^2/G proportional N ist eine Buchhaltung, solange nichts N festlegt; die "
       "Theorie fixiert weder die Nahtzahl noch die Fläche je Naht (z68, z90, z119; die Zählung selbst rechnet z137)", True)

# --- Vorzeichenproblem
lam_ab = lambda a, b: -12*sp.Rational(a - b, a + b)                    # z104: Lambda R_Vak^2 für Q_0 = a E + b Q^-
check(f"[E] Vorzeichen: Lambda R_Vak^2 = -12 (a - b)/(a + b) für Q_0 = a E + b Q^- (z104) gibt {lam_ab(1, 0)} für die Uhr E "
      f"(nativ, AdS), {lam_ab(1, 1)} parabolisch und {lam_ab(0, 1)} für die Kappe (dS, Prämisse P_cap). COSMO.LAM.01 verlangt "
      "rho_L > 0. Die native Kombination wäre rho_L/Mbar^4 = -48 c_3^2 e^(-2/alpha)",
      lam_ab(1, 0) == -12 and lam_ab(1, 1) == 0 and lam_ab(0, 1) == 12 and hat(z104, "Lambda = -12"))
VORZ = {
    "(a) Boltzmann/KMS": "r = e^(-4 pi) gehört zur elliptischen Uhr E (AdS, diskrete positive Niveaus). Unter P_cap ist die Zeit "
                         "der Kappenfluss mit beta = 2 pi und kontinuierlichem, symmetrischem Spektrum (z104): es gibt dann keine "
                         "Niveaus, an denen ein Boltzmann-Gewicht r^E abgelesen werden könnte. Der Mechanismus und das beobachtete "
                         "Vorzeichen schließen einander aus",
    "(b) Transmutation": "legt nur einen Betrag m a fest; das Vorzeichen von Lambda bleibt die Zeitklasse (P_cap). Unberührt, aber "
                         "auch nicht geholfen",
    "(c) Instanton": "ein verdünntes Instantongas senkt die Vakuumenergie (-2 K e^(-S) cos theta): das natürliche Vorzeichen "
                     "ist negativ, passend zum nativen AdS und gegen die Beobachtung",
    "(d) Zählung": "S_dS = pi l^2/G zählt einen kosmologischen Horizont, den es nur für Lambda > 0 gibt; im nativen AdS fehlt "
                   "die Fläche, auf der gezählt würde. Die Zählung setzt P_cap voraus",
}
ZAHLEN["vorzeichen"] = VORZ
for k, v in VORZ.items():
    say(f"      {k}: {v}")
befund("Vorzeichen: Kein Mechanismus (a)-(d) liefert zugleich die große Zahl und Lambda > 0. (a) und P_cap schließen "
       "einander aus, (c) gibt natürlich das AdS-Vorzeichen, (d) setzt P_cap voraus, (b) ist neutral", True)
ERG4 = {"(a) Boltzmann/KMS": "vorhanden, andere Zahl (10^-118,70 oder 10^-120,07 statt 10^-119,03; kein Niveau ausgezeichnet)",
        "(b) Transmutation": "vorhanden, andere Zahl (m a = 0,83 statt e^(-1/alpha) = 10^-59,5)",
        "(c) Instanton": "kein solcher Mechanismus (Wirkung 2 pi/alpha statt 2/alpha; 12 Nullmoden statt 48)",
        "(d) Zählung": "kein solcher Mechanismus (N frei)"}
ZAHLEN["mechanismen"] = ERG4

# =================================================================== Teil 5
say("\nTeil 5: Bilanz und Urteil")
bilanz = [
    "[E] Ein Torsor: Die physikalische Dimension aller nativen Größen hat Rang 1 modulo hbar, c, k_B; eine absolute Einheit "
    "kann die Theorie prinzipiell nicht ausgeben (No-Unit, v153).",
    "[E] Offene reine Zahlen statt offener Einheiten: G/R_Vak^2 (H6), Mbar/v_EW (Flavour-Absolutwert), dazu alpha_s(M_Z), M_R, "
    "T_reh als externe Eingaben (Matrix B: Rang 3).",
    "[E] Zeit: tau = l theta, beta_phys = 2 pi R_Vak/c, Takt der Uhr R = pi R_Vak/c, Stromquant 2 hbar c/R_Vak; exakte "
    "Verhältnisse, keine Einheit. beta l = 2 pi R_Vak ist exakt, als Gibbons-Hawking-Beziehung aber Zufall der Normierung "
    "(Faktor 2 zum nativen Horizontradius l).",
    "[E/C] Durch einen Mechanismus (innerhalb der jeweiligen Annahmen): alpha^-1, die Massenverhältnisse der Leiter, "
    "l/R_Vak = 1/2, |Lambda| R_Vak^2 = 12, die Zeit- und Energieverhältnisse der Naht.",
    "[N/C] Nur numerisch oder bedingt: rho_L/Mbar^4 = 48 c_3^2 e^(-2/alpha) (0,09 % neben Planck, innerhalb der Messbreite; "
    "Vorfaktor allein nicht signifikant), R_Vak/l_Pl = 2,06e61 unter P_Lambda und P_cap, M_scal/Mbar = c_3^(7/2), "
    "Quarkverhältnisse aus dreistelliger Tabelle, M_f/M_g als Lesart.",
    "[O] Die große Zahl hat keinen nativen Mechanismus: (a) und (b) geben andere Zahlen, (c) und (d) gibt es nicht; das "
    "Vorzeichen (nativ AdS) verschärft (a), (c), (d).",
]
for zl in bilanz:
    say("   " + zl)
urteil = ("Die Theorie erzeugt Maßverhältnisse, aber kein Maß. (1) Als abgeleitete Ergebnisse (durch einen Mechanismus, "
          "innerhalb der jeweils benannten Annahmen) liefert sie reine Zahlen: alpha^-1 = 137,03599921684 als Wurzel der "
          "Fixpunktgleichung; die Massenverhältnisse der phi_0-Leiter (m_mu/m_tau = (8/7) phi_0, m_e/m_mu = (12/7) phi_0^2, "
          "m_u/m_d = 55/117) und jede geladene Masse in Einheiten von v_EW; im Kandidaten H6 die Kinematik l = R_Vak/2, "
          "|Lambda| R_Vak^2 = 12 und die Zeit: Eigenzeit = l mal Uhrwinkel, inverse Temperatur 2 pi R_Vak/c, Takt der Uhr R "
          "gleich pi R_Vak/c, Rückkehr R^4 = 1 nach 4 pi R_Vak/c, Stromquant 2 hbar c/R_Vak, Lücke hbar c/(2 R_Vak). Das sind "
          "exakte Zeit- und Energieverhältnisse, keine Zeiteinheit. beta l = 2 pi R_Vak gilt exakt, ist aber keine erzwungene "
          "Gibbons-Hawking-Beziehung: der native Horizont hätte den Radius l, der Faktor 2 hat keinen Mechanismus (z104). "
          "(2) Nur als numerische Treffer oder unter benannten Prämissen: rho_Lambda/Mbar^4 = 48 c_3^2 e^(-2/alpha) = "
          "7,125e-121 (2,2375 meV gegen 2,2395 meV aus Planck-Werten, 0,09 %, innerhalb der Messbreite von etwa 0,5 %; Ledger: "
          "Numerical, bedingt auf den Lambda-Zweig) und damit l_dS = 1,03e61 l_Pl = 17,6 Mrd. Lichtjahre; im Kandidaten H6 "
          "R_Vak = 2,06e61 l_Pl, Takt 110,4 Mrd. Jahre, T = 1,10e-30 K unter den Prämissen P_Lambda (COSMO.LAM.01 als "
          "Eingabe), P_cap (dS-Vorzeichen) und MacDowell-Mansouri-Wirkung, alternativ R_Vak = 0,4 bis 10 l_Pl unter den "
          "O(1)-Prämissen von z112; M_scal/Mbar = c_3^(7/2) (ungemessen); Quarkverhältnisse m_c/m_s, m_t/m_b aus dreistelliger "
          "Tabelle; M_f/M_g als Lesart. In der Familie n c_3^k (448 Ausdrücke) liegen vier im 2-sigma-Fenster des Vorfaktors "
          "(n = 47 bis 50 bei k = 2; zufällig erwartet 1,3): die Messung legt ihn nur auf etwa 2 % fest, er allein trägt "
          "nicht, tragend ist der Exponent 2/alpha. "
          "(3) Prinzipiell nicht erzeugbar ist die eine Einheit: Alle nativen Größen mit Dimension liegen auf einem einzigen "
          "Torsor (Rang 1 modulo hbar, c, k_B), eine zweite unabhängige Einheit gibt es nicht. Offen sind reine Zahlen: "
          "G/R_Vak^2 im Gravitationssektor von H6, Mbar/v_EW (H_EW = 37,18; die neun Massen lassen sich mit G allein nicht in "
          "Kilogramm umrechnen), dazu alpha_s(M_Z), M_R und T_reh als externe Eingaben. Für die große Zahl e^(-2/alpha) hat "
          "H6 keinen Mechanismus: (a) Das native Boltzmann-Gewicht r = e^(-4 pi) verlangt E = alpha^-1/(2 pi) = 21,81, kein "
          "Niveau in Z/4 (nächste Gewichte 10^-118,70 und 10^-120,07), und in der Fixpunktgleichung steht e^(-2 alpha), nicht "
          "e^(-2/alpha); rho_L G^2 = Q(1/alpha) ist eine Schreibweise, kein Mechanismus. (b) Die native Kette hat dimensionale "
          "Transmutation, gibt aber m a = 4 e^(-pi/2) = 0,83; nötig wäre Theta = 44,06 statt 1/2 bzw. ein Koeffizient b = pi. "
          "(c) Kein Instanton hat die Wirkung 2/alpha (vorhanden ist 2 pi/alpha), und 48 zählt Weyl-Zustände, nicht Nullmoden "
          "(12). (d) Die Nahtzahl ist frei. Das Vorzeichenproblem verschärft (a), (c), (d): die Niveaus gehören zur AdS-Uhr, "
          "Instantonen senken die Vakuumenergie, und eine Horizontzählung braucht P_cap. Entscheiden würde eine einzige "
          "Rechnung: der native Wert des MacDowell-Mansouri-Koeffizienten b = R_Vak^2/(256 pi G), also die Normierung des "
          "Einsteinterms der relativen Uhrspur für das volle, wechselwirkende Nahtsystem. Heute liefert sie nur ganze Zahlen "
          "(64 je Index, z118); COSMO.LAM.01 verlangt |b| = (pi^2/2) e^(2/alpha). Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z140_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil, "bilanz": bilanz}, indent=1,
    ensure_ascii=False, default=str) + "\n")
