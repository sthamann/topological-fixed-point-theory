"""Z142: Die Planck-Wurzel rho_* des Archivmanuskripts tfpt-42 -- ist p_G = G/R_Vak^2 aus dem dort gestellten
Funktional W(rho) mit nativen H6-Daten berechenbar?

Auftrag (Strategie S4 der Gravitationsdurchsicht R1, Punkt 10): tfpt-42 ('Canonical spectral Einstein functional and
geometric Hodge closure', 'Absolute spectral Planck closure') ist die einzige Stelle des Korpus, an der die reine Zahl
p_G als berechenbar gestellt wird: M_Pl^2 = chi_0^2/(2 pi^2), chi_0 = lambda_Sigma sqrt(rho_*), rho_* die eindeutige
stationäre Wurzel von W(rho) = -6 E^_D rho + Gamma^(4)(rho) + rho^2 L^_Sigma(rho) unter angenommener strenger
Konvexität, also G_N lambda_Sigma^2 = pi/(4 rho_*). rho_* wurde nie berechnet. Dieses Skript rekonstruiert das
Funktional aus den Quellzeilen, prüft die Algebra, setzt für die nur benannten Zutaten die nativen H6-Objekte ein
(als Identifikation gekennzeichnet) und rechnet rho_*, p_G und Lambda G mit Sensitivitäten.

Leitplanken: H6 und die gelösten Blöcke (alpha^-1 = 137.03599921684, c_3 = 1/(8 pi), phi_0, die 16 in drei Familien,
E8) sind fest; nichts wird befördert; keine vollständige Theorie wird behauptet. Anti-Numerologie-Regel (Ledger
v354/v355): eine passende Zahl ohne Mechanismus ist kein Ergebnis. Andere z-Skripte, ihr JSON, Ledger und Manuskripte
werden nur gelesen. Statuskürzel: E exakt (bzw. exakt gelesen), C bedingt (Identifikation), L Lesart, O offen.

Bestand (nur gelesen): z107 (Tr e^(-K_beta) = q^(1/24) theta_3/eta bei beta = 4 pi), z98/z104 (Lambda = -12/R_Vak^2,
AdS, l = R_Vak/2), z112 (R_Vak = 0,40 bis 10,03 Planck-Längen unter P_L), z127 (Lambda G = -3/(1024 pi kappa)),
z90 (xi = 1/6 induziert 1/G = 0; Flächenkoeffizient ist Regulatorfunktional), z27 (Gitterkopplung nicht universell),
z41 mit Ledgernotiz v1037 (Spin-2-Geist der Spektralwirkung bei x* = 6,768), v68 (a_2 = -d R/(192 pi^2)),
v259 (f_2/f_0 = 1 für das KMS-Gewicht), v164 (QCD.LAMBDA.01), z99 (laufende Strom-Strom-Kopplung).

Teil 1  Rekonstruktion: Quellzeilen jeder Zutat (spezifiziert oder nur benannt); Algebra a_2 -> M_Pl^2 =
        chi_0^2/(2 pi^2) -> G chi_0^2 = pi/4 (48 Komponenten, f_2 = 1) mit Negativkontrollen; a_4 des Dirac-Operators
        (R^2-Koeffizient, Weyl-Koeffizient, Anschluss an z41); Vorzeichen und Normierung von W.
Teil 2  rho_*: spezifizierte Zutaten gerechnet, benannte nativ identifiziert [C]; stationäre Punkte, Konvexität,
        Sensitivitäten; p_G = G/R_Vak^2 und Lambda G; Vergleich mit z127, z112, Flächenzustand und Ledgerziel.
Teil 3  Hierarchieformel v_geo = Mbar g_car beta_rad^2 exp[-(alpha^-1 + delta_ph)/5]: Herkunft, delta_ph, Zahl,
        dimensionale Transmutation, Anti-Numerologie.
Teil 4  Bilanz und Urteil (Folge für die Prämisse P_G).
"""
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.optimize import brentq

HERE = Path(__file__).parent
ROOT = HERE.parents[2]
MS = ROOT / "documents" / "archive--paper-latex--tfpt-42" / "sections"
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(142)                                   # nur für die Ableitungsprobe (feste Saat)
mp.mp.dps = 40
PI = mp.pi

# feste Blöcke (nur gelesen, nichts wird befördert)
AINV = mp.mpf("137.03599921684")
ALPHA = 1/AINV
C3 = 1/(8*PI)
PHIBASE = 1/(6*PI)
DTOP = 48*C3**4                                                    # 3/(256 pi^4)
PHI0 = PHIBASE + DTOP
B1 = mp.mpf(41)/10
GCAR = 5
MBAR = mp.mpf("2.435323203e18")                                    # GeV, reduzierte Planck-Masse (tfpt_constants, gelesen)
V_EW = mp.mpf("246.21965")                                         # GeV, (sqrt2 G_F)^(-1/2)
H_EW_LEDGER = mp.mpf("37.1776")                                    # Ledger L. 238/670 (gelesen)


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


def fz(x, n=6):
    return f"{float(x):.{n}g}"


def lies_json(name):
    try:
        return json.loads((HERE / name).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def ms_datei(prefix):
    c = sorted(MS.glob(prefix + "-*.tex"))
    return c[0] if len(c) == 1 else None


_MS_CACHE = {}


def ms_zeilen(prefix):
    if prefix not in _MS_CACHE:
        f = ms_datei(prefix)
        _MS_CACHE[prefix] = f.read_text(encoding="utf-8").split("\n") if f else []
    return _MS_CACHE[prefix]


def ms_hat(prefix, nr, token):
    z = ms_zeilen(prefix)
    return 1 <= nr <= len(z) and token in z[nr - 1]


# =================================================================== Teil 1
say("Teil 1: Rekonstruktion des Funktionals aus tfpt-42 (documents/archive--paper-latex--tfpt-42/sections)")

# --- 1.1 Quellzeilen: (Zutat, Datei-Präfix, Zeile, Marke in der Zeile, Status)
QUELLEN = [
    ("E_rel = -(1/4) R 1_48 (48 Komponenten)", "013", 26, r"\mathbf 1_{48}", "spezifiziert"),
    ("a_4-Zeile der Residuentafel: Skalierung chi^0", "013", 45, r"$\chigeo^0$", "spezifiziert"),
    ("E_D: -R/(24 pi^2)", "013", 108, r"\frac{1}{24\pi^2}R", "spezifiziert"),
    ("E_D: -Phi^+Phi/(48 pi^2)", "013", 110, r"\frac{1}{48\pi^2}\Phi^\dagger\Phi", "spezifiziert"),
    ("E_D: Nahtanteil E_Sigma", "013", 113, r"\mathcal E_\Sigma", "nur benannt"),
    ("Gamma_D^(4): R^2/(6 M^2)", "013", 128, r"\frac{1}{6M^2}R^2", "spezifiziert (M nicht an chi gebunden)"),
    ("Gamma_D^(4): Nahtanteil A_Sigma", "013", 131, r"\mathcal A_\Sigma", "nur benannt"),
    ("Gamma_grav := -6 chi^2 E_D + Gamma_D^(4)", "013", 141, r"-6\chigeo^2\mathcal E_D+\Gamma_D^{(4)}", "spezifiziert"),
    ("derselbe R^2-Term als Mbar^2 R^2/(12 M^2)", "013", 151, r"\frac{\Mpl^2}{12M^2}R^2", "spezifiziert (widerspricht Zeile 128 um chi^2)"),
    ("Mbar^2 = chi_0^2/(2 pi^2)", "013", 161, r"\Mpl^2=\frac{\chigeozero^2}{2\pi^2}", "spezifiziert"),
    ("G_N chi_0^2 = pi/4", "013", 195, r"G_N\chigeozero^2=\frac{\pi}{4}", "spezifiziert"),
    ("E_D := Res_(s=1) Z_rel", "011", 376, r"\operatorname*{Res}_{s=1}", "spezifiziert (D_geo, B_rel abstrakt)"),
    ("Gamma_D^(4) := FP_(s=0) Z_rel", "011", 380, r"\operatorname{FP}_{s=0}", "spezifiziert (D_geo, B_rel abstrakt)"),
    ("Wirkung trägt -Lambda_(Sigma,ren)", "015", 32, r"-\Lambda_{\Seam,\mathrm{ren}}", "spezifiziert"),
    ("Behauptung: positiv auf TT-Moden", "015", 117, "positive on transverse traceless modes", "behauptet"),
    ("lambda_Sigma := lambda_1^+(|B_Sigma|)", "015", 158, r"\lambda_\Sigma:=\lambda_1^+(|B_\Sigma|)", "definiert (Spektrum von B_Sigma nicht angegeben)"),
    ("rho := chi^2/lambda_Sigma^2", "015", 182, r"\rho&:=\frac{\chigeo^2}{\lambda_\Sigma^2}", "spezifiziert"),
    ("Auswertung bei g = g_vac, Phi = Phi_vac", "015", 190, r"g=g_{\mathrm{vac}}", "nur benannt (g_vac nicht angegeben)"),
    ("E^_D := lambda_Sigma^-2 E_D", "015", 194, r"\widehat{\mathcal E}_D:=\lambda_\Sigma^{-2}\mathcal E_D", "spezifiziert"),
    ("U^_Sigma(rho) := U_Sigma(lambda_Sigma sqrt(rho))", "015", 196, r"U_\Sigma(\lambda_\Sigma\sqrt{\rho})", "nur benannt (chi-Abhängigkeit nie definiert)"),
    ("L^_Sigma := -log det_adm(1 - U^_Sigma)", "015", 199, r"-\log\det\nolimits_{\mathrm{adm}}(1-\widehat U_\Sigma(\rho))", "spezifiziert bei gegebenem U"),
    ("Gamma^(4)(rho) := lambda^-4 Gamma_D^(4)(lambda sqrt(rho))", "015", 201, r"\lambda_\Sigma^{-4}\Gamma_D^{(4)}", "spezifiziert"),
    ("W: -6 E^_D rho", "015", 206, r"-6\widehat{\mathcal E}_D\rho", "spezifiziert"),
    ("W: + rho^2 L^_Sigma(rho)", "015", 208, r"+\rho^2\widehat L_\Sigma(\rho)", "spezifiziert"),
    ("Annahme d^2 Gamma^(4) > 0", "015", 216, r"\partial_\rho^2\widehat\Gamma_D^{(4)}(\rho)>0", "angenommen"),
    ("Annahme d^2 [rho^2 L] >= 0", "015", 218, r"\rho^2\widehat L_\Sigma(\rho)\bigr]\ge 0", "angenommen"),
    ("Annahme dV(0+) < 0", "015", 272, r"(0^+)<0", "angenommen"),
    ("Annahme strenge Konvexität", "015", 274, r">0\ \text{for all }\rho>0", "angenommen"),
    ("G_N lambda_Sigma^2 = pi/(4 rho_*)", "015", 341, r"\frac{\pi}{4\rho_\star}", "spezifiziert"),
    ("Mbar^-2 = 8 pi G_N", "015", 364, r"\Mpl^{-2}=8\pi G_N", "spezifiziert"),
    ("Vakuumzweig: Schwarzschild-de Sitter", "015", 367, "Schwarzschild de Sitter", "spezifiziert (dS)"),
    ("U_Sigma := P_adm e^(-l_Sigma |B_Sigma|) P_adm", "029", 14, r"\Padm e^{-\ell_\Seam|B_\Seam|}\Padm", "spezifiziert bis auf l_Sigma"),
    ("S_Sigma = alpha^2 (alpha - 2 c_3^3)/(8 c_3^6 b_1)", "029", 44, r"{8\cthree^6\,b_1}", "spezifiziert"),
    ("Lambda_IR = M_Pl^4 (unreduziert) [...]", "029", 81, r"\Mplunred^4", "spezifiziert"),
    ("... [-log det_adm(1 - U_Sigma)]", "029", 82, r"-\log\det\nolimits_{\mathrm{adm}}(1-U_\Sigma)", "spezifiziert"),
    ("phi_seam(alpha) = phi_base + q (1 - q)^(-5/4)", "007", 46, r"\phibase+q(\alpha)(1-q(\alpha))^{-5/4}", "spezifiziert"),
    ("IR-Fortsetzung: S_IR = 2/alpha (Vermutung)", "040", 332, r"\frac{2}{\alpha_\star}", "Vermutung"),
    ("bedingte Zahl 1.128e-123", "040", 368, "1.1280419085743317e-123", "bedingt"),
    ("Lemma zum Zweiblattgewicht w_H", "013", 206, r"\begin{lemma}[Two-sheet seam-even zero-mode normalization]", "behauptet"),
    ("Exponent -(alpha^-1(0) + delta_ph)/5", "013", 213, r"-\frac{\alpha^{-1}(0)+\delta_{\mathrm{ph}}}{5}", "behauptet"),
    ("mu_Phi^2 = chi_0^2 w_H^2/(8 pi^2)", "013", 217, r"\frac{\chigeozero^2}{8\pi^2}w_H^2", "behauptet"),
    ("nur Beweisskizze", "013", 231, r"\begin{proof}[Proof sketch]", "Beweisskizze"),
    ("beta_rad = phi_0/(4 pi)", "036", 66, r"\beta_{\mathrm{rad}}=\frac{\phi_{\mathrm{ret}}^0}{4\pi}", "spezifiziert"),
    ("z_* = (794 - 7 sqrt(9961))/2187", "017", 112, r"\frac{794-7\sqrt{9961}}{2187}", "spezifiziert"),
    ("delta_ph = 0.593277280133224", "017", 117, "0.593277280133224", "spezifiziert"),
    ("delta_ph in (1/3, 2/3)", "017", 37, r"\left(\frac13,\frac23\right)", "spezifiziert"),
]
fehl = [f"{p}:{n} {name}" for (name, p, n, tok, st) in QUELLEN if not ms_hat(p, n, tok)]
ZAHLEN["quellen"] = [{"zutat": name, "datei": (ms_datei(p).name if ms_datei(p) else p), "zeile": n, "status": st}
                     for (name, p, n, tok, st) in QUELLEN]
ZAHLEN["quellen_fehlend"] = fehl
for (name, p, n, tok, st) in QUELLEN:
    say(f"      {p}:{n:<4d} {name}  [{st}]")
check(f"[E] Quellenbeleg: alle {len(QUELLEN)} zitierten Formeln stehen in den angegebenen Zeilen von tfpt-42 "
      f"(Abschnitte 007, 011, 013, 015, 017, 029, 036, 040); fehlend: {fehl}", not fehl)

# l_Sigma und die chi-Abhängigkeit von U_Sigma: nur benannt
alle_tex = sorted(MS.glob("*.tex"))
ell_orte = []
uchi_orte = []
exp_orte = []
for f in alle_tex:
    for i, z in enumerate(f.read_text(encoding="utf-8").split("\n"), start=1):
        if r"\ell_\Seam" in z or r"\ell_\Sigma" in z:
            ell_orte.append(f"{f.name[:3]}:{i}")
        if r"U_\Sigma(\lambda_\Sigma" in z:
            uchi_orte.append(f"{f.name[:3]}:{i}")
        if r"\alpha^{-1}(0)+\delta_{\mathrm{ph}}" in z:
            exp_orte.append(f"{f.name[:3]}:{i}")
ZAHLEN["nur_benannt"] = {"l_Sigma kommt vor in": ell_orte, "U_Sigma(chi) kommt vor in": uchi_orte,
                         "Exponent (alpha^-1(0) + delta_ph) kommt vor in": exp_orte}
check(f"[E] nur benannt: die Nahtlänge l_Sigma kommt im ganzen Manuskript nur in der Definition von U_Sigma vor "
      f"({ell_orte}), eine chi-abhängige Naht U_Sigma(chi) nur in der einen Zeile {uchi_orte}; weder l_Sigma noch die "
      f"chi-Abhängigkeit von U_Sigma ist irgendwo definiert. Ebenso bleiben E_Sigma, A_Sigma, g_vac und das Spektrum von "
      f"B_Sigma ohne Angabe", ell_orte == ["029:14", "029:57", "029:73"] and uchi_orte == ["015:196"])

# --- 1.2 Algebra a_2 -> Mbar^2 -> G chi^2
Rs, Ns, chis, f2s, Gs = sp.symbols("R N chi f2 G", positive=True)
E_dirac = -Rs/4                                                    # Lichnerowicz: D^2 = nabla^* nabla + R/4
a2_dichte = sp.Rational(1, 16)/sp.pi**2*(Ns*Rs/6 + Ns*E_dirac)     # (4 pi)^-2 tr(E + R/6)
a2_koeff = sp.simplify(a2_dichte/Rs)
N_ED = sp.solve(sp.Eq(a2_koeff, -sp.Rational(1, 24)/sp.pi**2), Ns)[0]
a2_skalar_konform = sp.simplify(sp.Rational(1, 16)/sp.pi**2*(Rs/6 - Rs/6))   # E = -xi R mit xi = 1/6


def G_chi2(N, vorfaktor=1, f2=1):
    """G chi^2 aus (Mbar^2/2) R = -vorfaktor f2 chi^2 a_2(N) und Mbar^-2 = 8 pi G."""
    mbar2 = sp.simplify(-2*vorfaktor*f2*a2_koeff.subs(Ns, N))      # in Einheiten chi^2
    return sp.simplify(1/(8*sp.pi*mbar2)), mbar2


G48, M48 = G_chi2(48)
G_6x8, M_6x8 = G_chi2(8, vorfaktor=6)
check(f"[E] a_2 des Dirac-Typs (sympy): (4 pi)^-2 tr(E + R/6) mit E = -R/4 gibt {a2_koeff} R je N Komponenten "
      f"(wie v68); der Manuskriptkoeffizient -1/(24 pi^2) (013:108) entspricht N = {N_ED}, der Vorfaktor 6 in "
      f"-6 chi^2 E_D (013:141) macht daraus 6 x 8 = 48. Ein konformer Skalar (xi = 1/6) trägt {a2_skalar_konform} bei "
      f"(z90: nur Fermionen zählen)", a2_koeff == -Ns/(192*sp.pi**2) and N_ED == 8 and a2_skalar_konform == 0)
check(f"[E] Algebra (sympy): (Mbar^2/2) R = -chi^2 a_2(48) gibt Mbar^2 = {M48} chi^2 und mit Mbar^-2 = 8 pi G "
      f"G chi_0^2 = {G48}; die Manuskriptform 6 x 8 gibt dasselbe ({G_6x8}). f_2 kommt im Manuskript nicht vor: das "
      f"Residuum bei s = 1 ist a_2 selbst, das entspricht f_2 = 1",
      M48 == 1/(2*sp.pi**2) and G48 == sp.pi/4 and G_6x8 == sp.pi/4)
kontrolle_N = {N: G_chi2(N)[0] for N in (16, 32, 96)}
G_wortlich = G_chi2(48, vorfaktor=6)[0]                            # 1_48 wörtlich im Residuum und zusätzlich der Faktor 6
G_gauss = G_chi2(48, f2=sp.sqrt(sp.pi)/2)[0]                       # Gauß-Profil e^(-u^2): f_2 = sqrt(pi)/2 (v259)
ZAHLEN["algebra"] = {"a_2 je N": str(a2_koeff), "N im Manuskript-E_D": int(N_ED), "Mbar^2/chi^2 (N = 48)": str(M48),
                     "G chi^2 (N = 48)": str(G48), "G chi^2 bei N = 16, 32, 96": {str(k): str(v) for k, v in kontrolle_N.items()},
                     "G chi^2 mit 1_48 und Vorfaktor 6": str(G_wortlich), "G chi^2 mit Gauß-Profil": str(G_gauss)}
check(f"[E] Negativkontrolle: falsche Komponentenzahl ändert die Konstante, N = 16, 32, 96 geben G chi^2 = "
      f"{[str(kontrolle_N[k]) for k in (16, 32, 96)]}; 1_48 wörtlich im Residuum mit dem Vorfaktor 6 gäbe {G_wortlich}; "
      f"ein Gauß-Profil (f_2 = sqrt(pi)/2) gäbe {G_gauss}. pi/4 gilt genau für N x Vorfaktor = 48 und f_2 = 1. Die 16 in "
      f"drei Familien sind 48 Weyl-Fermionen, als Spinorspur gezählt 96: dann pi/8 (Faktor 2 offen, [C])",
      kontrolle_N[16] == 3*sp.pi/4 and kontrolle_N[32] == 3*sp.pi/8 and kontrolle_N[96] == sp.pi/8
      and G_wortlich == sp.pi/24 and G_gauss != sp.pi/4 and all(v != sp.pi/4 for v in kontrolle_N.values()))

rho_s, lam_s, wH = sp.symbols("rho lambda w_H", positive=True)
chi0 = lam_s*sp.sqrt(rho_s)
Mbar_s = sp.sqrt(chi0**2/(2*sp.pi**2))
G_s = 1/(8*sp.pi*Mbar_s**2)
mP = 1/sp.sqrt(G_s)
ident = [sp.simplify(Mbar_s/lam_s - sp.sqrt(rho_s)/(sp.sqrt(2)*sp.pi)), sp.simplify(G_s*lam_s**2 - sp.pi/(4*rho_s)),
         sp.simplify(mP/lam_s - 2*sp.sqrt(rho_s)/sp.sqrt(sp.pi)), sp.simplify(lam_s/mP - sp.sqrt(sp.pi)/(2*sp.sqrt(rho_s)))]
ZH, lamPhi = sp.Rational(1, 16)/sp.pi**2, sp.Rational(1, 16)/sp.pi**2
v2 = ZH*(chi0**2*wH**2/(8*sp.pi**2))/lamPhi
rang_P2 = sp.Rational(1, 48)*6*16                                  # 6 x 1/(48 pi^2) = rang/(16 pi^2)
check(f"[E] Folgerungen (sympy): Mbar/lambda = sqrt(rho)/(sqrt2 pi), G lambda^2 = pi/(4 rho), m_P/lambda = "
      f"2 sqrt(rho)/sqrt(pi), lambda l_P = sqrt(pi)/(2 sqrt(rho)) (015:327-351); v_geo^2 = Z_H mu^2/lambda_Phi = "
      f"Mbar^2 w_H^2/4 (013:232-238); der Higgs-Term 6/(48 pi^2) entspricht Rang P_2 = {rang_P2}",
      all(x == 0 for x in ident) and sp.simplify(v2 - Mbar_s**2*wH**2/4) == 0 and rang_P2 == 2)

# --- 1.3 a_4 des Dirac-Operators: R^2- und Weyl-Koeffizient, Anschluss an z41
R2, Ric2, Riem2, C2, E4, e_s = sp.symbols("R2 Ric2 Riem2 C2 E4 e")
lsg = sp.solve([sp.Eq(C2, Riem2 - 2*Ric2 + R2/3), sp.Eq(E4, Riem2 - 4*Ric2 + R2)], [Ric2, Riem2], dict=True)[0]
# Gilkey: 360 (4 pi)^2 a_4 = tr(60 R E + 180 E^2 + 30 Omega^2) + N (5 R^2 - 2 Ric^2 + 2 Riem^2), E = -e R, tr Omega^2 = -N Riem^2/8
a4_360 = sp.expand((Ns*(5*R2 - 2*Ric2 + 2*Riem2) - 60*e_s*Ns*R2 + 180*e_s**2*Ns*R2 - sp.Rational(30, 8)*Ns*Riem2).subs(lsg))
kR2 = sp.factor(a4_360.coeff(R2))
kC2 = a4_360.coeff(C2)
kE4 = a4_360.coeff(E4)
kR2_dirac = kR2.subs(e_s, sp.Rational(1, 4))
c_C = kC2/(360*16*sp.pi**2)                                        # Koeffizient von C^2 in a_4 (t^0)
c_R = a2_koeff                                                     # Koeffizient von R in a_2 (t^-1)
verh = sp.simplify(c_C/c_R)                                        # C^2 gegen R, in Einheiten der Eigenzeit
lin_z41 = -2*verh                                                  # TT: R -> (1/4) h Box h, C^2 -> (1/2) h Box^2 h (Standard, zitiert)
z41 = lies_json("z41_ergebnis.json")
z41_sp = (z41 or {}).get("zahlen", {}).get("spektralwirkung_spin2", {})
x_trunk = sp.Rational(1)/(2*verh)
ZAHLEN["a4"] = {"R^2-Koeffizient (360 (4 pi)^2 a_4/N) als Funktion von e (E = -e R)": str(kR2), "bei e = 1/4": str(kR2_dirac),
                "C^2-Koeffizient": str(kC2), "E_4-Koeffizient": str(kE4), "C^2 in a_4": str(c_C),
                "Verhältnis C^2/R": str(verh), "linearer TT-Koeffizient": str(lin_z41),
                "Nullstelle der a_2 + a_4-Trunkierung": str(x_trunk), "z41 x*": z41_sp.get("x_stern"), "z41 Reihe": z41_sp.get("reihe")}
check(f"[E] a_4 des Dirac-Operators (Gilkey, sympy, Basis C^2, E_4, R^2): 360 (4 pi)^2 a_4 = N ({kC2/Ns} C^2 + {kE4/Ns} E_4 "
      f"+ {kR2/Ns} R^2); bei E = -R/4 (e = 1/4) ist der R^2-Koeffizient {kR2_dirac} (Negativkontrolle e = 0: "
      f"{kR2.subs(e_s, 0)}). Der im Manuskript behaltene Term R^2/(6 M^2) (013:128) entsteht also nicht aus dem "
      f"Residuum von D^2 allein; was entsteht, ist der Weyl-Term {c_C} C^2 mit negativem Vorzeichen, den das Manuskript "
      f"unter '...' führt",
      kR2_dirac == 0 and kR2.subs(e_s, 0) != 0 and kC2 == -sp.Rational(9, 2)*Ns and kE4 == sp.Rational(11, 4)*Ns)
check(f"[E] Anschluss an z41 (gelesen): das Verhältnis Weyl zu Einstein ist {verh} (in Einheiten der Eigenzeit), der "
      f"lineare Koeffizient des TT-Formfaktors damit {lin_z41}, genau der erste Term der z41-Reihe "
      f"'{z41_sp.get('reihe')}'. Schon die a_2 + a_4-Trunkierung hat eine euklidische Nullstelle bei x = {x_trunk}, die "
      f"volle KMS-Resummation bei x* = {fz(z41_sp.get('x_stern', float('nan')), 5)} (z41, Residuum "
      f"{fz(z41_sp.get('residuum', float('nan')), 5)}): 'positiv auf TT-Moden' (015:117) gilt nicht",
      verh == sp.Rational(3, 20) and lin_z41 == -sp.Rational(3, 10) and "3x/10" in str(z41_sp.get("reihe", ""))
      and abs(float(z41_sp.get("x_stern", 0)) - 6.768124) < 1e-5 and x_trunk == sp.Rational(10, 3))

# --- 1.4 Vorzeichen und Normierung von W
a_s, b_s, r_s, k_s, V4 = sp.symbols("a b r k V4", positive=True)
chi_s = sp.Symbol("chi", positive=True)
S_wirk = V4*r_s**4*(a_s*chi_s**2*12/r_s**2 - b_s*chi_s**4)        # maximal symmetrisch, R = 12/r^2 (dS), Wirkung a chi^2 R - b chi^4
gleich = sp.simplify(chi_s*sp.diff(S_wirk, chi_s) - r_s*sp.diff(S_wirk, r_s))
Rsym = sp.Symbol("Rk", real=True)
R_kons = sp.solve(sp.diff(a_s*chi_s**2*Rsym - b_s*chi_s**4, chi_s), Rsym)[0]        # konsistent: V = -(a chi^2 R - b chi^4)
R_wort = sp.solve(sp.diff(a_s*chi_s**2*Rsym + b_s*chi_s**4, chi_s), Rsym)[0]        # wörtlich: +6 E_D-Term und +rho^2 L
kappa_unred = sp.simplify((8*sp.pi)**2*(1/(2*sp.pi**2))**2)        # M_Pl^4 = (8 pi)^2 Mbar^4, Mbar^2 = chi^2/(2 pi^2)
kappa_red = sp.simplify((1/(2*sp.pi**2))**2)
ZAHLEN["vorzeichen_normierung"] = {"R der chi-Stationarität (konsistent)": str(R_kons), "R (wörtliches W)": str(R_wort),
                                   "kappa wie geschrieben": "1", "kappa aus Lambda_IR = M_Pl^4 L (029:81)": str(kappa_unred),
                                   "kappa aus rho_Lambda = Mbar^4 L (Ledgerform)": str(kappa_red)}
check(f"[E] Vorzeichen (sympy; gelesen Lambda_(Sigma,ren) = Lambda_IR > 0): Für die Wirkung a chi^2 R - b chi^4 (015:27-45: "
      f"Mbar^2 R/2 - Lambda) ist die "
      f"chi-Stationarität bei festem g dieselbe Gleichung wie die Spur der Einsteingleichung (chi d_chi S = r d_r S, "
      f"Rest {gleich}), mit R = {R_kons} > 0 (dS). Im gedruckten W = -6 E^_D rho + ... + rho^2 L^ (015:205-208) haben "
      f"Einstein- und Nahtterm dasselbe relative Vorzeichen (in der Wirkung entgegengesetzte, in jeder Konvention); "
      f"seine Wurzel liegt bei R = {R_wort}, dem Negativen davon. W ist so, wie es gedruckt ist, nicht das Potential der "
      f"eigenen Wirkung",
      gleich == 0 and sp.simplify(R_kons - 2*b_s*chi_s**2/a_s) == 0 and sp.simplify(R_wort + R_kons) == 0)
check(f"[E] Normierung (sympy): Lambda_IR = M_Pl^4 [-log det(1 - U)] mit unreduziertem M_Pl (029:79-82) ist in "
      f"Randeinheiten kappa rho^2 L^ mit kappa = {kappa_unred}, mit reduziertem Mbar^4 kappa = {kappa_red}; gedruckt ist "
      f"rho^2 L^ (kappa = 1). Drei Normierungen, Spanne {fz(kappa_unred/kappa_red, 5)}",
      kappa_unred == 16/sp.pi**2 and kappa_red == 1/(4*sp.pi**4) and kappa_unred != 1)
befund("E: Vollständig spezifiziert sind die lokalen Koeffizienten von E_D und Gamma_D^(4), die Definitionen von rho, "
       "E^_D, L^_Sigma, W, die Planck-Algebra und der Spektralradius von U_Sigma (029:34-44). Nur benannt sind E_Sigma, "
       "A_Sigma, g_vac, l_Sigma, die chi-Abhängigkeit von U_Sigma, das Spektrum von B_Sigma und die chi-Bindung von M. "
       "Dazu drei innere Unstimmigkeiten: Vorzeichen von W, Normierung des Nahtterms (kappa), R^2-Term ohne Quelle im "
       "Residuum von D^2", True)

# =================================================================== Teil 2
say("\nTeil 2: rho_* mit spezifizierten und nativ identifizierten Zutaten")

# --- 2.1 Zutaten
q_seam = DTOP*mp.e**(-2*ALPHA)
phi_seam = PHIBASE + q_seam*(1 - q_seam)**(-mp.mpf(5)/4)
S_ms = mp.log(1/phi_seam)
S_ms_rhs = ALPHA**2*(ALPHA - 2*C3**3)/(8*C3**6*B1)
L_ms = -mp.log(1 - phi_seam)                                       # Rang eins, Manuskript
x_bed = DTOP*mp.e**(-2*AINV)
L_bedingt = -mp.log1p(-x_bed)                                      # bedingte Schicht 040:358-381
L_bed_ms = mp.mpf("1.1280419085743317e-123")                       # gedruckt in 040:368
abw_bed = abs(L_bedingt/L_bed_ms - 1)
ZAHLEN["manuskript_naht"] = {"phi_seam(alpha)": mp.nstr(phi_seam, 15), "S_Sigma = log(1/phi_seam)": mp.nstr(S_ms, 15),
                             "rechte Seite 029:44": mp.nstr(S_ms_rhs, 15), "L (Rang eins)": mp.nstr(L_ms, 12),
                             "bedingte Schicht (040:368)": mp.nstr(L_bedingt, 12),
                             "log10 Abstand Satzschicht zu bedingter Schicht": mp.nstr(mp.log10(phi_seam/L_bedingt), 6)}
check(f"[E] Spezifiziert im Manuskript: der Spektralradius der Naht ist rho(U_Sigma) = e^(-S_Sigma) = phi_seam(alpha) = "
      f"{mp.nstr(phi_seam, 10)} (029:29-45; die Identität S_Sigma = alpha^2 (alpha - 2 c_3^3)/(8 c_3^6 b_1) ist die "
      f"alpha-Fixpunktgleichung, relativer Rest {mp.nstr(abs(S_ms - S_ms_rhs)/S_ms, 3)} beim festen alpha), also "
      f"L^_Sigma >= tr U >= {mp.nstr(phi_seam, 6)} nach der eigenen Schranke 029:96-106, im Rang-eins-Fall "
      f"L^ = {mp.nstr(L_ms, 8)}", abs(S_ms - S_ms_rhs)/S_ms < mp.mpf("1e-9") and mp.mpf("0.05") < phi_seam < mp.mpf("0.06"))
check(f"[E] Die bedingte Schicht braucht dagegen rho(U_Sigma)|_IR = delta_top e^(-2/alpha) und gibt "
      f"{mp.nstr(L_bedingt, 12)} (040:368 nachgerechnet, relative Abweichung {mp.nstr(abw_bed, 2)} aus den letzten "
      f"Stellen von alpha): {mp.nstr(mp.log10(phi_seam/L_bedingt), 5)} Größenordnungen "
      f"unter dem Spektralradius der Satzschicht. Dieselbe Größe trägt im Manuskript zwei Werte; der kleine ist als "
      f"Vermutung geführt (040:324-340)",
      abw_bed < mp.mpf("1e-9") and mp.log10(phi_seam/L_bedingt) > 121)

QN = mp.e**(-2*PI)                                                 # q = e^(-beta/2) bei beta = 4 pi (z107)
prod_nat = mp.nprod(lambda n: (1 + QN**(n - mp.mpf(1)/2))**2, [1, mp.inf])
theta_eta = mp.jtheta(3, 0, mp.e**(-PI))/mp.nprod(lambda n: 1 - QN**n, [1, mp.inf])   # q^(1/24) theta_3/eta
z107 = lies_json("z107_ergebnis.json")
z107_wert = mp.mpf(str((z107 or {}).get("zahlen", {}).get("Tr_e^-K_bei_beta_4pi", {}).get("Produkt", "0")))
L_F = mp.log(prod_nat)                                             # fermionisch: log Tr e^(-K_beta)
L_B = -2*mp.nsum(lambda n: mp.log(1 - QN**(n - mp.mpf(1)/2)), [1, mp.inf])           # wörtlich -log det(1 - U)
L_FC = L_F + 2*PI/24                                               # mit Casimir: log(theta_3/eta)
L_VARIANTEN = {"nativ fermionisch log(q^(1/24) theta_3/eta)": L_F, "nativ, wörtlich -log det(1 - U)": L_B,
               "nativ mit Casimir log(theta_3/eta)": L_FC, "Manuskript Rang eins -log(1 - phi_seam)": L_ms}
ZAHLEN["naht_nativ"] = {"Tr e^(-K) bei beta = 4 pi": mp.nstr(prod_nat, 25), "z107": mp.nstr(z107_wert, 25),
                        "L-Varianten": {k: mp.nstr(v, 12) for k, v in L_VARIANTEN.items()},
                        "nativer Spektralradius e^(-pi)": mp.nstr(mp.e**(-PI), 10), "native Dämpfung S_Sigma": "pi",
                        "nötige Dämpfung 2/alpha": mp.nstr(2*AINV, 12), "Verhältnis (2/alpha)/pi": mp.nstr(2*AINV/PI, 10)}
check(f"[C] Identifikation I3 (native Naht): U_Sigma <-> Einteilchengewichte q^(n-1/2), je zweifach, q = e^(-2 pi) "
      f"(z107, beta = 4 pi). Tr e^(-K_beta) = prod (1 + q^(n-1/2))^2 = q^(1/24) theta_3/eta = {mp.nstr(prod_nat, 20)} "
      f"(Abweichung Produkt gegen Thetaquotient {mp.nstr(abs(prod_nat - theta_eta), 3)}, gegen z107 "
      f"{mp.nstr(abs(prod_nat - z107_wert), 3)}); L^ = {mp.nstr(L_F, 10)} (fermionisch), {mp.nstr(L_B, 10)} (wörtliche "
      f"bosonische Determinante), {mp.nstr(L_FC, 10)} (mit Casimir). Alle Lesarten sind positiv und von Ordnung 0,1 wie der "
      f"Manuskriptwert {mp.nstr(L_ms, 6)}",
      abs(prod_nat - theta_eta) < mp.mpf("1e-35") and abs(prod_nat - z107_wert) < mp.mpf("1e-23")
      and all(mp.mpf("0.05") < v < mp.mpf("0.4") for v in L_VARIANTEN.values()))
verh_daempf = 2*AINV/PI
check(f"[E] Native Dämpfung: -log rho(U_nat) = pi (unterstes Gewicht q^(1/2)); das Ledgerziel braucht 2/alpha = "
      f"{mp.nstr(2*AINV, 9)}, das {mp.nstr(verh_daempf, 8)}-Fache, kein ganzzahliges Vielfaches (Abstand zur nächsten "
      f"ganzen Zahl {mp.nstr(abs(verh_daempf - mp.nint(verh_daempf)), 4)}). Keine native Potenz von e^(-pi) gibt e^(-2/alpha)",
      abs(-mp.log(QN**(mp.mpf(1)/2)) - PI) < mp.mpf("1e-35") and abs(verh_daempf - mp.nint(verh_daempf)) > mp.mpf("0.1"))

wH_num = 2*GCAR*(PHI0/(4*PI))**2*mp.e**(-(AINV + ((794 - 7*mp.sqrt(9961))/2187)**(mp.mpf(1)/6))/5)
higgs_rho2 = wH_num**4/(16*PI**2)                                  # |V_min| = mu^4/(4 lambda) = rho^2 w_H^4/(16 pi^2)
g4_konst = sp.diff(sp.Symbol("A0") + 0*rho_s, rho_s, 2)
g4_higgs = sp.diff(-sp.Symbol("h", positive=True)*rho_s**2, rho_s, 2)
ZAHLEN["gamma4"] = {"w_H": mp.nstr(wH_num, 8), "Higgs-Anteil |rho^2-Koeffizient|": mp.nstr(higgs_rho2, 6),
                    "d^2 (Konstante)": str(g4_konst), "d^2 (Higgs-Vakuumwert, Potentialvorzeichen)": str(g4_higgs)}
check(f"[E] Gamma^(4)(rho): nach der eigenen Residuentafel skaliert der a_4-Sektor mit chi^0 (013:45), auf dem "
      f"homogenen Vakuum ist er also rho-unabhängig, d^2/drho^2 = {g4_konst}; der einzige rho-abhängige Rest, der "
      f"Higgs-Vakuumwert rho^2 w_H^4/(16 pi^2), ist im Potentialvorzeichen konkav und mit {mp.nstr(higgs_rho2, 3)} um rund "
      f"{mp.nstr(mp.log10(L_F/higgs_rho2), 3)} Größenordnungen unter L^. Die Annahme d^2 Gamma^(4) > 0 (015:216) wird von "
      f"der spezifizierten Form nicht erfüllt; die Konvexität von W kann nur vom Nahtterm kommen",
      g4_konst == 0 and g4_higgs.is_negative and higgs_rho2 < mp.mpf("1e-60"))


# --- 2.2 Modell und Stationaritätssuche
def turm(t, eta):
    """Naht-Turm: Energien n - 1/2, je zweifach; eta = -1 fermionisch (log det(1 + U)), +1 bosonisch (-log det(1 - U))."""
    nmax = int(min(4e5, 45.0/t + 3))
    e = np.arange(1, nmax + 1) - 0.5
    y = np.exp(-t*e)
    den = 1 - eta*y
    return (-eta*2*float(np.sum(np.log1p(-eta*y))), -2*float(np.sum(e*y/den)), 2*float(np.sum(e*e*y/den**2)))


def rang_eins(t, eta=1):
    y = math.exp(-t)
    return -math.log1p(-y), -y/(1 - y), y/(1 - y)**2


T0_NAT = 2*math.pi
L_MODELLE = {"F": (lambda t: turm(t, -1), T0_NAT), "B": (lambda t: turm(t, +1), T0_NAT),
             "R1": (rang_eins, float(S_ms))}


def g_abl(rho, lesart, modell, kappa):
    """g(rho) = kappa rho^2 L(t(rho)) mit erster und zweiter Ableitung. Lesart A: t fest; B: t = t0/sqrt(rho)
    (l_Sigma ~ 1/chi); C: t = t0 sqrt(rho) (l_Sigma ~ chi/lambda^2)."""
    Lfun, t0 = L_MODELLE[modell]
    if lesart == "A":
        t, t1, t2 = t0, 0.0, 0.0
    elif lesart == "B":
        t, t1, t2 = t0*rho**-0.5, -0.5*t0*rho**-1.5, 0.75*t0*rho**-2.5
    else:
        t, t1, t2 = t0*rho**0.5, 0.5*t0*rho**-0.5, -0.25*t0*rho**-1.5
    L, Lt, Ltt = Lfun(t)
    return (kappa*rho**2*L, kappa*(2*rho*L + rho**2*Lt*t1),
            kappa*(2*L + 4*rho*Lt*t1 + rho**2*(Ltt*t1**2 + Lt*t2)))


def a_koeff(c, N=48):
    """|6 E^_D| = N |R^|/(192 pi^2) mit R^ = R_vac/lambda^2 = 48/c^2 (l = R_Vak/2, lambda_Sigma = c/R_Vak)."""
    return N*48.0/(192*math.pi**2*c**2)


GITTER = np.logspace(-4, 4, 1601)


def stationaer(fp, gitter=GITTER):
    """Alle Vorzeichenwechsel von fp auf dem Gitter, mit brentq verfeinert."""
    w = np.array([fp(x) for x in gitter])
    wurzeln = []
    for i in range(len(gitter) - 1):
        if w[i] == 0.0:
            wurzeln.append(float(gitter[i]))
        elif w[i]*w[i + 1] < 0:
            wurzeln.append(float(brentq(fp, gitter[i], gitter[i + 1], xtol=1e-15, rtol=1e-13)))
    return wurzeln


# Proben der Maschinerie
test_w = stationaer(lambda x: (x - 1.0)*(x - 3.0)*(x - 10.0))
test_0 = stationaer(lambda x: x*x + 1.0)
fd_fehler = 0.0
for lesart in ("A", "B", "C"):
    for modell in ("F", "B", "R1"):
        for rho in np.exp(rng.uniform(math.log(0.05), math.log(20.0), 3)):
            h = 1e-4*rho
            gm, gp, g0 = g_abl(rho - h, lesart, modell, 1.0), g_abl(rho + h, lesart, modell, 1.0), g_abl(rho, lesart, modell, 1.0)
            d1 = (gp[0] - gm[0])/(2*h)
            d2 = (gp[1] - gm[1])/(2*h)
            fd_fehler = max(fd_fehler, abs(d1 - g0[1])/(abs(g0[1]) + 1e-12), abs(d2 - g0[2])/(abs(g0[2]) + 1e-12))
LF_float = turm(T0_NAT, -1)[0]
check(f"[E] Maschinerie: die Suche findet genau die Nullstellen einer Probe ({[fz(x, 6) for x in test_w]}) und keine, wo "
      f"keine ist ({test_0}); analytische Ableitungen von g = kappa rho^2 L(t(rho)) stimmen mit Differenzenquotienten "
      f"(größter relativer Fehler {fd_fehler:.1e} über drei Lesarten, drei Nahtmodelle); der Gleitkomma-Turm gibt L^ = "
      f"{LF_float:.12f} wie mpmath",
      len(test_w) == 3 and max(abs(test_w[0] - 1), abs(test_w[1] - 3), abs(test_w[2] - 10)) < 1e-10 and test_0 == []
      and fd_fehler < 1e-5 and abs(LF_float - float(L_F)) < 1e-13)

# --- 2.3 Lesart A (U_Sigma ist Randdatum nach 029:14, hängt nicht von chi ab)
C_KAND = {"1": 1.0, "2": 2.0, "2 pi": 2*math.pi, "4 pi": 4*math.pi}
KAPPA = {"wie gedruckt (1)": 1.0, "M_Pl^4 unreduziert (16/pi^2)": 16/math.pi**2, "Mbar^4 reduziert (1/(4 pi^4))": 1/(4*math.pi**4)}
A_tab = {}
A_ok = True
for cn, c in C_KAND.items():
    a = a_koeff(c)
    w = stationaer(lambda x: -a + g_abl(x, "A", "F", 1.0)[1])
    geschl = a/(2*float(L_F))
    kruemm = g_abl(1.0, "A", "F", 1.0)[2]
    A_ok = A_ok and len(w) == 1 and abs(w[0]/geschl - 1) < 1e-10 and kruemm > 0
    A_tab[cn] = {"rho_*": w[0] if w else None, "geschlossen 6/(pi^2 c^2 L)": geschl, "G lambda^2": math.pi/(4*geschl),
                 "p_G": math.pi/(4*geschl)/c**2}
pG_A = [A_tab[cn]["p_G"] for cn in C_KAND]
pG_A_geschl = float(PI**3*L_F/24)
ZAHLEN["lesart_A"] = {"rho_* je c (N = 48, kappa = 1, L fermionisch)": A_tab, "p_G = pi^3 L/24": pG_A_geschl,
                      "W'' = 2 kappa L": 2*float(L_F)}
check(f"[C] Lesart A (U_Sigma hängt nicht von chi ab, wie es seine Definition 029:14 als Randdatum sagt; I1: lambda_Sigma "
      f"= c/R_Vak, I2: g_vac = natives AdS mit R = -48/R_Vak^2, E_Sigma = 0, wörtliches Vorzeichen): W = -a rho + kappa "
      f"rho^2 L^ ist streng konvex (W'' = 2 kappa L^ = {2*float(L_F):.6f}), genau eine Wurzel rho_* = a/(2 kappa L^) = "
      f"6/(pi^2 c^2 kappa L^): {', '.join(f'c = {cn}: {A_tab[cn]["rho_*"]:.6g}' for cn in C_KAND)} (Suche = "
      f"geschlossene Form)", A_ok)
check(f"[E] In Lesart A fällt lambda_Sigma heraus: p_G = G/R_Vak^2 = (pi/(4 rho_*))/c^2 = pi^3 kappa L^/24 = "
      f"{pG_A_geschl:.8f} für alle vier c (Streuung {max(pG_A) - min(pG_A):.1e}). Die Stationarität ist die Spur der "
      f"Vakuum-Einsteingleichung, |Lambda| G = (pi^3/2) kappa L^; die 'Randeinheit' trägt nichts bei, und rho_* selbst "
      f"ist nur die Umrechnung c -> rho_* ~ c^-2. Negativkontrolle: G lambda_Sigma^2 allein ist nicht c-unabhängig "
      f"({', '.join(f'{A_tab[cn]["G lambda^2"]:.4g}' for cn in C_KAND)})",
      max(pG_A) - min(pG_A) < 1e-12 and abs(pG_A[0] - pG_A_geschl) < 1e-12
      and max(A_tab[cn]["G lambda^2"] for cn in C_KAND)/min(A_tab[cn]["G lambda^2"] for cn in C_KAND) > 100)

# Vorzeichentafel: s = Vorzeichen des linearen Terms
vz = {}
for lesung, vak, s in (("wörtlich", "AdS", -1), ("wörtlich", "dS", +1), ("konsistent", "dS", -1), ("konsistent", "AdS", +1)):
    w = stationaer(lambda x, s=s: s*a_koeff(1.0) + g_abl(x, "A", "F", 1.0)[1])
    vz[f"{lesung}, {vak}"] = len(w)
ZAHLEN["vorzeichentafel (Zahl positiver Wurzeln)"] = vz
check(f"[E] Vorzeichentafel (L^ > 0 in jeder Lesart): positive Wurzeln {vz}. Mit dem Vakuum des Manuskripts selbst "
      f"(Schwarzschild-de Sitter, 015:367) hat das gedruckte W keine positive Wurzel, weil dV(0+) = -6 E^_D > 0 ist "
      f"(Annahme 015:272 verletzt, solange E_Sigma nicht eingreift); mit dem nativen AdS hat es eine, widerspricht dann "
      f"aber der Spur der eigenen Einsteingleichung; das vorzeichenkonsistente Potential hat seine Wurzel nur bei dS",
      vz == {"wörtlich, AdS": 1, "wörtlich, dS": 0, "konsistent, dS": 1, "konsistent, AdS": 0})

# --- 2.4 Lesarten B und C (chi-abhängige Naht), je c und kappa
BC = {}
B_konvex, B_eindeutig, C_nichtkonvex, C_zahlen = True, True, True, []
for lesart in ("B", "C"):
    g2 = np.array([g_abl(x, lesart, "F", 1.0)[2] for x in GITTER])
    g1 = np.array([g_abl(x, lesart, "F", 1.0)[1] for x in GITTER])
    BC[lesart] = {"min d2g": float(g2.min()), "max d1g": float(g1.max()), "d1g am rechten Rand": float(g1[-1]), "je kappa, c": {}}
    if lesart == "B":
        B_konvex = B_konvex and g2.min() > 0
    else:
        C_nichtkonvex = C_nichtkonvex and g2.min() < 0 and g1[-1] < 0
    for kn, kap in KAPPA.items():
        for cn, c in C_KAND.items():
            a = a_koeff(c)
            w = stationaer(lambda x: -a + g_abl(x, lesart, "F", kap)[1])
            art = ["Minimum" if g_abl(x, lesart, "F", kap)[2] > 0 else "Maximum" for x in w]
            BC[lesart]["je kappa, c"][f"{kn}; c = {cn}"] = {"Wurzeln": w, "Art": art,
                                                           "p_G (Minimum)": [math.pi/(4*x)/c**2 for x, t in zip(w, art) if t == "Minimum"]}
            if lesart == "B":
                B_eindeutig = B_eindeutig and len(w) == 1 and art == ["Minimum"]
            else:
                C_zahlen.append(len(w))
                C_nichtkonvex = C_nichtkonvex and len(w) in (0, 2) and (len(w) == 0 or art == ["Minimum", "Maximum"])
pG_B = [v["p_G (Minimum)"][0] for v in BC["B"]["je kappa, c"].values()]
pG_C = [p for v in BC["C"]["je kappa, c"].values() for p in v["p_G (Minimum)"]]
ZAHLEN["lesarten_B_C"] = BC
B_gedruckt = {cn: BC["B"]["je kappa, c"][f"wie gedruckt (1); c = {cn}"] for cn in C_KAND}
check(f"[C/L] Lesart B (l_Sigma ~ 1/chi, t = 2 pi/sqrt(rho); natives Gewicht bei rho = 1): rho^2 L^(rho) ist konvex "
      f"(min g'' = {BC['B']['min d2g']:.3e} > 0 auf rho = 1e-4..1e4), genau eine Wurzel je c und kappa; bei kappa = 1: "
      f"{', '.join(f'c = {cn}: rho_* = {B_gedruckt[cn]["Wurzeln"][0]:.5g}, p_G = {B_gedruckt[cn]["p_G (Minimum)"][0]:.4g}' for cn in C_KAND)}. "
      f"Hier hängt p_G von c ab; über alle c und kappa {min(pG_B):.3g} bis {max(pG_B):.3g}", B_konvex and B_eindeutig)
check(f"[C/L] Lesart C (l_Sigma ~ chi, t = 2 pi sqrt(rho)): rho^2 L^(rho) ist nicht konvex (min g'' = "
      f"{BC['C']['min d2g']:.3e} < 0), g' hat ein Maximum {BC['C']['max d1g']:.4g} und wird negativ; W fällt für rho -> "
      f"unendlich wie -a rho, hat also kein globales Minimum. Je nach c und kappa gibt es keine oder zwei stationäre "
      f"Punkte (Minimum, dann Maximum): Anzahlen {C_zahlen}. Die im Manuskript angenommene Konvexität ist eine "
      f"Eigenschaft der nicht definierten chi-Abhängigkeit von U_Sigma, nicht des Funktionals", C_nichtkonvex and 0 in C_zahlen and 2 in C_zahlen)

# --- 2.5 p_G und Lambda G, Sensitivitäten, Vergleich
p_ziel = mp.e**(-2*AINV)/(128*PI**3)
lg_ziel = float(mp.log10(p_ziel))
TAB = {}
for ln_, Lv in L_VARIANTEN.items():
    for kn, kap in KAPPA.items():
        for N in (48, 96):
            pG = 96*PI**3*mp.mpf(kap)*Lv/N**2
            if kn.startswith("Mbar^4") or kn.startswith("M_Pl^4"):
                pG = pG*(mp.mpf(N)/48)**2                          # kappa ~ Mbar^4 ~ N^2: N fällt heraus
            TAB[f"{ln_}; {kn}; N = {N}"] = {"p_G": float(pG), "Lambda G (AdS, Lambda = -12/R_Vak^2)": float(-12*pG),
                                            "R_Vak/l_Pl": float(1/mp.sqrt(pG)), "log10 p_G": float(mp.log10(pG)),
                                            "Größenordnungen über dem Ziel": float(mp.log10(pG)) - lg_ziel,
                                            "implizites kappa_MM (z127)": float(1/(4096*PI*pG))}
pG_A_alle = [v["p_G"] for v in TAB.values()]
alle_pG = pG_A_alle + pG_B + pG_C
haupt = TAB["nativ fermionisch log(q^(1/24) theta_3/eta); wie gedruckt (1); N = 48"]
haupt_u = TAB["nativ fermionisch log(q^(1/24) theta_3/eta); M_Pl^4 unreduziert (16/pi^2); N = 48"]
haupt_r = TAB["nativ fermionisch log(q^(1/24) theta_3/eta); Mbar^4 reduziert (1/(4 pi^4)); N = 48"]
sens = {"d ln rho_*/d ln L^": -1, "d ln rho_*/d ln c": -2, "d ln rho_*/d ln N": 1, "d ln rho_*/d ln kappa": -1,
        "d ln p_G/d ln L^": 1, "d ln p_G/d ln c": 0, "d ln p_G/d ln N (kappa fest)": -2, "d ln p_G/d ln kappa": 1,
        "d rho_*/d E^_Sigma = 3/(kappa L^) (kappa = 1)": float(3/L_F)}
ZAHLEN["p_G"] = {"Tafel Lesart A": TAB, "Sensitivitäten Lesart A": sens, "Ziel p_G": mp.nstr(p_ziel, 10), "log10 Ziel": lg_ziel,
                 "Spanne p_G über alle Lesarten": [min(alle_pG), max(alle_pG)], "Spanne p_G Lesart A": [min(pG_A_alle), max(pG_A_alle)],
                 "Spanne p_G Lesart B": [min(pG_B), max(pG_B)], "Spanne p_G Lesart C (Nebenminima)": [min(pG_C), max(pG_C)]}
say(f"      Lesart A, nativ fermionisch, N = 48: p_G = {haupt['p_G']:.5f} (kappa = 1), {haupt_u['p_G']:.5f} (16/pi^2), "
    f"{haupt_r['p_G']:.4e} (1/(4 pi^4)); Lambda G = {haupt['Lambda G (AdS, Lambda = -12/R_Vak^2)']:.4f}, "
    f"{haupt_u['Lambda G (AdS, Lambda = -12/R_Vak^2)']:.4f}, {haupt_r['Lambda G (AdS, Lambda = -12/R_Vak^2)']:.4e}; "
    f"R_Vak = {haupt['R_Vak/l_Pl']:.3f}, {haupt_u['R_Vak/l_Pl']:.3f}, {haupt_r['R_Vak/l_Pl']:.2f} Planck-Längen")
check(f"[C] p_G über alle Identifikationen (vier Nahtlesarten, drei Normierungen, N = 48 und 96, Lesarten A, B, C, "
      f"c = 1, 2, 2 pi, 4 pi): {min(alle_pG):.3g} bis {max(alle_pG):.3g}, also R_Vak = {1/math.sqrt(max(alle_pG)):.2f} bis "
      f"{1/math.sqrt(min(alle_pG)):.0f} Planck-Längen (das obere Ende von p_G kommt aus den Nebenminima der Lesart C). Jede "
      f"Fassung liegt mindestens "
      f"{min(math.log10(p) for p in alle_pG) - lg_ziel:.1f} Größenordnungen über dem Ledgerziel 10^{lg_ziel:.2f}: ein "
      f"planckgroßes Vakuum, in keiner Lesart die Hierarchie",
      min(math.log10(p) for p in alle_pG) - lg_ziel > 115 and max(abs(math.log10(p)) for p in alle_pG) < 4)

# Ziel als Forderung an L^: exakt die Vermutung des Manuskripts bzw. COSMO.LAM.01
pt_s = sp.exp(-2*sp.Symbol("ainv", positive=True))/(128*sp.pi**3)
L_noetig = {kn: sp.simplify(24*pt_s/(sp.pi**3*ks)) for kn, ks in
            (("wie gedruckt (1)", sp.Integer(1)), ("M_Pl^4 unreduziert (16/pi^2)", kappa_unred), ("Mbar^4 reduziert (1/(4 pi^4))", kappa_red))}
dtop_s = sp.Rational(3, 256)/sp.pi**4
ex = sp.exp(-2*sp.Symbol("ainv", positive=True))
rhoL_ledger = 3/(4*PI**2)*mp.e**(-2*AINV)
ZAHLEN["ziel_als_L"] = {k: str(v) for k, v in L_noetig.items()}
ZAHLEN["ziel_als_L"]["COSMO.LAM.01 rho_Lambda/Mbar^4"] = mp.nstr(rhoL_ledger, 10)
check(f"[E] Das Ledgerziel p_G = e^(-2/alpha)/(128 pi^3) ist in Lesart A identisch mit einer Forderung an die Naht: "
      f"mit der Normierung des Manuskripts (M_Pl^4) L^_Sigma = {L_noetig['M_Pl^4 unreduziert (16/pi^2)']} = delta_top "
      f"e^(-2/alpha), wörtlich die Vermutung 040:324-340; mit reduziertem Mbar^4 L^ = "
      f"{L_noetig['Mbar^4 reduziert (1/(4 pi^4))']} = COSMO.LAM.01 ({mp.nstr(rhoL_ledger, 7)}, Status 'Numerical'). Die "
      f"Wurzelgleichung ist COSMO.LAM.01 in anderen Variablen und fügt nichts hinzu. Negativkontrolle: mit kappa = 1 "
      f"wäre es {L_noetig['wie gedruckt (1)']}, keine der beiden Formen",
      sp.simplify(L_noetig["M_Pl^4 unreduziert (16/pi^2)"] - dtop_s*ex) == 0
      and sp.simplify(L_noetig["Mbar^4 reduziert (1/(4 pi^4))"] - 3*ex/(4*sp.pi**2)) == 0
      and sp.simplify(L_noetig["wie gedruckt (1)"] - dtop_s*ex) != 0
      and abs(rhoL_ledger/mp.mpf("7.12533e-121") - 1) < mp.mpf("1e-5"))

# Vergleich mit dem Bestand
z112 = lies_json("z112_ergebnis.json")
z127 = lies_json("z127_ergebnis.json")
sp112 = (z112 or {}).get("zahlen", {}).get("prämissen_G", {}).get("Spanne log10 G/R_Vak^2", [float("nan")]*2)
mm127 = (z127 or {}).get("zahlen", {}).get("mm_normierung", {})
nat127 = mm127.get("natives kappa: [Lambda G, G/R_Vak^2]", {})
pG_index = nat127.get("Index (1)", [float("nan")]*2)[1]
pG_2c3 = nat127.get("2 c_3", [float("nan")]*2)[1]
pG_flaeche = float(1/(72*mp.zeta(3)*C3))                           # G = l^2/(18 zeta(3) c_3), l = R_Vak/2
naechste = min(math.exp(abs(math.log(h["p_G"]/x))) for h in (haupt, haupt_u, haupt_r) for x in (pG_flaeche, pG_index, pG_2c3))
in112 = {k: bool(sp112[0] <= v["log10 p_G"] <= sp112[1]) for k, v in TAB.items()}
ZAHLEN["vergleich"] = {"z112 Spanne log10 p_G": sp112, "z127 p_G (Index, kappa = 1)": pG_index, "z127 p_G (2 c_3)": pG_2c3,
                       "z127 Formel": mm127.get("G/R_Vak^2"), "Flächenzustand p_G = 1/(72 zeta(3) c_3)": pG_flaeche,
                       "nächste Übereinstimmung (Faktor)": naechste,
                       "in der z112-Spanne": in112,
                       "implizites kappa_MM (nativ fermionisch, N = 48)": {"kappa = 1": haupt["implizites kappa_MM (z127)"],
                                                                         "16/pi^2": haupt_u["implizites kappa_MM (z127)"],
                                                                         "1/(4 pi^4)": haupt_r["implizites kappa_MM (z127)"]}}
check(f"[E/C] Vergleich: z112 (gelesen) log10 p_G in {sp112}; z127 (gelesen) p_G = 1/(4096 pi kappa) = "
      f"{pG_index:.4e} (Index), {pG_2c3:.4e} (2 c_3); Flächenzustand G = l^2/(18 zeta(3) c_3) mit l = R_Vak/2: p_G = "
      f"{pG_flaeche:.4f}. Diese Route gibt {haupt['p_G']:.4f} bzw. {haupt_u['p_G']:.4f} (in der z112-Spanne) und "
      f"{haupt_r['p_G']:.3e} (reduzierte Normierung, darunter, zwischen den beiden z127-Werten; implizites kappa_MM = "
      f"{haupt_r['implizites kappa_MM (z127)']:.3f}). Alle Routen liegen innerhalb von rund vier Größenordnungen "
      f"beieinander und rund 119 bis 123 über dem Ziel; die nächste Übereinstimmung dieser Route mit einer anderen ist "
      f"ein Faktor {naechste:.2f}, daraus folgt nichts",
      mm127.get("G/R_Vak^2") == "1/(4096*pi*kappa)" and abs(pG_index - 7.7712e-5) < 1e-8 and abs(sp112[0] + 2.0023) < 1e-9
      and in112["nativ fermionisch log(q^(1/24) theta_3/eta); wie gedruckt (1); N = 48"]
      and in112["nativ fermionisch log(q^(1/24) theta_3/eta); M_Pl^4 unreduziert (16/pi^2); N = 48"]
      and not in112["nativ fermionisch log(q^(1/24) theta_3/eta); Mbar^4 reduziert (1/(4 pi^4)); N = 48"]
      and pG_index < haupt_r["p_G"] < pG_2c3 and abs(pG_flaeche - 0.2904) < 1e-4
      and naechste > 1.5)
befund("C: rho_* ist nur nach Identifikationen berechenbar. In Lesart A (U_Sigma als Randdatum) ist rho_* = 6/(pi^2 c^2 kappa L^) "
       "eindeutig und W konvex, aber die Gleichung ist die Spur der Vakuum-Einsteingleichung: lambda_Sigma fällt heraus, "
       "und p_G = pi^3 kappa L^/24 ist die Vakuumenergie der Naht in Planck-Einheiten. Nativ ist L^ von Ordnung 0,1, das "
       "Vakuum planckgroß", True)
befund("E: Die Konvexitätsannahme des Manuskripts ist nicht allgemein wahr: der Anteil Gamma^(4) ist nicht streng konvex "
       "(null bzw. winzig konkav), der Nahtterm ist konvex in den Lesarten A und B und nicht konvex in Lesart C", True)

# =================================================================== Teil 3
say("\nTeil 3: Hierarchieformel v_geo = Mbar g_car beta_rad^2 exp[-(alpha^-1 + delta_ph)/5]")
zs = sp.Symbol("z")
P_pol = (zs - 1)*(zs - sp.Rational(64, 729))*(zs - sp.Rational(1, 729))
z_stern = (794 - 7*sp.sqrt(9961))/2187
dP = sp.simplify(sp.diff(P_pol, zs).subs(zs, z_stern))
delta_ph = mp.mpf(sp.N(z_stern, 45))**(mp.mpf(1)/6)
kontrollform = 1594323*delta_ph**12 - 1157652*delta_ph**6 + 47449
check(f"[E] delta_ph des Manuskripts (017:41-118): P'(z_*) = {dP} für z_* = (794 - 7 sqrt(9961))/2187 (kleinere "
      f"kritische Stelle der Kubik), delta_ph = z_*^(1/6) = {mp.nstr(delta_ph, 16)}, im zulässigen Intervall (1/3, 2/3); "
      f"Kontrollform zwölften Grades {mp.nstr(kontrollform, 3)}. delta_ph ist die Transportpolstelle des C_6-Sechsecks, "
      f"eine algebraische Zahl ohne Bezug zur Skala",
      dP == 0 and abs(delta_ph - mp.mpf("0.593277280133224")) < mp.mpf("1e-15") and mp.mpf(1)/3 < delta_ph < mp.mpf(2)/3
      and abs(kontrollform) < mp.mpf("1e-28"))

beta_rad = PHI0/(4*PI)
sockel = -mp.log(GCAR) - 2*mp.log(beta_rad) + AINV/5
H_ms = sockel + delta_ph/5
H_246 = mp.log(MBAR/V_EW)
H_174 = mp.log(MBAR/(V_EW/mp.sqrt(2)))
d_noetig_ledger = 5*(H_EW_LEDGER - sockel)
d_noetig_246 = 5*(H_246 - sockel)
v_geo = MBAR*mp.e**(-H_ms)
ledger_text = (ROOT / "verification" / "status_ledger.csv").read_text(encoding="utf-8")
ZAHLEN["hierarchie"] = {"beta_rad": mp.nstr(beta_rad, 12), "Sockel -ln 5 - 2 ln beta_rad + alpha^-1/5": mp.nstr(sockel, 10),
                        "ln(Mbar/v_geo) mit delta_ph des Manuskripts": mp.nstr(H_ms, 10), "v_geo [GeV]": mp.nstr(v_geo, 7),
                        "ln(Mbar/246.22)": mp.nstr(H_246, 10), "ln(Mbar/174.10)": mp.nstr(H_174, 10), "H_EW Ledger": "37.1776",
                        "delta_ph nötig für H_EW (Ledger)": mp.nstr(d_noetig_ledger, 6),
                        "delta_ph nötig für v = 246.22": mp.nstr(d_noetig_246, 6),
                        "Abweichung zu H_EW": mp.nstr(H_ms - H_EW_LEDGER, 5), "Abweichung zu ln(Mbar/246.22)": mp.nstr(H_ms - H_246, 5)}
check(f"[E] Zahl: mit dem delta_ph des Manuskripts ist ln(Mbar/v_geo) = {mp.nstr(sockel, 8)} + delta_ph/5 = "
      f"{mp.nstr(H_ms, 8)}, v_geo = {mp.nstr(v_geo, 6)} GeV. Der Ledgerwert H_EW = 37,1776 ist ln(Mbar/(v/sqrt2)) mit "
      f"v/sqrt2 = 174,1 GeV ({mp.nstr(H_174, 8)}); dafür bräuchte die Formel delta_ph = {mp.nstr(d_noetig_ledger, 5)}, "
      f"außerhalb von (1/3, 2/3): der Wert 2,25 wird von nichts erzeugt, er ist die Rückrechnung aus H_EW. Gegen "
      f"H_EW fehlen {mp.nstr(H_EW_LEDGER - H_ms, 4)} im Logarithmus (Faktor {mp.nstr(mp.e**(H_EW_LEDGER - H_ms), 4)}); "
      f"gegen v = 246,22 GeV (ln = {mp.nstr(H_246, 8)}) liegt die Formel {mp.nstr((v_geo/V_EW - 1)*100, 3)} % daneben "
      f"(nötig wäre delta_ph = {mp.nstr(d_noetig_246, 4)})",
      abs(H_174 - H_EW_LEDGER) < mp.mpf("2e-3") and "H_EW = ln(Mbar/v_EW) = 37.1776" in ledger_text
      and abs(d_noetig_ledger - mp.mpf("2.2467")) < mp.mpf("1e-3") and d_noetig_ledger > mp.mpf(2)/3
      and abs(H_ms - mp.mpf("36.8469144")) < mp.mpf("1e-6") and abs(v_geo/V_EW - 1) < mp.mpf("0.02")
      and abs(H_EW_LEDGER - H_ms) > mp.mpf("0.3"))

skizze = ms_zeilen("013")[231] if len(ms_zeilen("013")) > 231 else ""
exp_dateien = sorted({o[:3] for o in exp_orte})
check(f"[E] Herkunft des Exponenten: (alpha^-1(0) + delta_ph)/5 steht nur in den Abschnitten {exp_dateien} ({exp_orte}); "
      f"022, 036 und 037 wiederholen ihn und verweisen auf das Lemma (022:394). Eingeführt wird er im Lemma 013:206-229 als "
      f"Zweiblattgewicht w_H, dessen Beweisskizze (013:231-239) nur sagt, die "
      f"Zweiblattnormierung 'doubles the unnormalized zero-mode weight and therefore produces the stated w_H', und dann "
      f"v_geo = Mbar w_H/2 ausrechnet. Der Exponent, die 5 und der Vorfaktor g_car beta_rad^2 sind behauptet, nicht "
      f"hergeleitet; das Manuskript nennt v_geo zudem ausdrücklich nicht die physikalische elektroschwache Skala "
      f"(013:265-277, Faktor Z_EW offen)",
      exp_dateien == ["013", "022", "036", "037"] and ms_hat("022", 394, "lem:two-sheet-zero-mode-normalization")
      and "therefore produces the stated $w_H$" in skizze
      and ms_hat("013", 270, "not identified with") and ms_hat("013", 231, "Proof sketch"))

# dimensionale Transmutation: Skalenverhältnis exp(-(2 pi/b) alpha^-1) aus d alpha^-1/d ln mu = -b/(2 pi)
b_noetig = 2*sp.pi*5
B_KAND = {"U(1)_Y, b_1 = 41/10": sp.Rational(41, 10), "SU(2), |b_2| = 19/6": sp.Rational(19, 6), "SU(3), |b_3| = 7 (QCD.LAMBDA.01, v164)": sp.Integer(7),
          "QED, alle SM-Fermionen: 32/3": sp.Rational(32, 3), "z99 Strom-Strom: 2 h-dual = 4": sp.Integer(4), "g_car = 5": sp.Integer(5)}
trans = {k: float(2*sp.pi/b) for k, b in B_KAND.items()}
abst = min(abs(x - 0.2) for x in trans.values())
lam_z99 = float(5*sp.pi/AINV)                                      # exp(-pi/lambda) = exp(-alpha^-1/5) verlangte lambda = 5 pi alpha
ZAHLEN["transmutation"] = {"nötiges b für exp(-alpha^-1/5)": str(b_noetig), "Exponentfaktor 2 pi/b je Kandidat": trans,
                           "kleinster Abstand zu 1/5": abst, "z99: nötiges lambda = 5 pi alpha": lam_z99}
check(f"[E] Mechanismusprüfung (dimensionale Transmutation): ein Einschleifenlauf d alpha^-1/d ln mu = -b/(2 pi) gibt "
      f"Skalenverhältnisse exp(-(2 pi/b) alpha^-1); für exp(-alpha^-1/5) bräuchte es b = {b_noetig} = "
      f"{float(b_noetig):.4f}, irrational, während Einschleifenkoeffizienten rational sind. Kandidaten 2 pi/b: "
      f"{ {k: round(v, 4) for k, v in trans.items()} }; kleinster Abstand zu 1/5 ist {abst:.3f}. QCD.LAMBDA.01 hat die Form "
      f"exp(-2 pi/(7 alpha_s)); die native laufende Strom-Strom-Kopplung aus z99 (Koeffizient 2 h-dual = 4, Betrag und "
      f"Vorzeichen von lambda nicht nativ) bräuchte lambda = 5 pi alpha = {lam_z99:.4f} als gesetzten Wert. Die 5 ist "
      f"g_car, kein Betakoeffizient. Negativkontrolle: b = 10 pi gibt genau 1/5",
      abst > 0.3 and sp.simplify(2*sp.pi/b_noetig - sp.Rational(1, 5)) == 0 and not b_noetig.is_rational
      and all(b.is_rational for b in B_KAND.values()))

# Anti-Numerologie: wie viele Formeln derselben Bauart treffen ebenso gut?
tol = abs(H_ms - H_246)
treffer_246, treffer_HEW, enthalten = [], [], False
for n in range(1, 9):
    for a_ in range(-2, 3):
        for b_ in range(0, 5):
            for dn, dv in (("0", mp.mpf(0)), ("delta_ph", delta_ph)):
                val = -a_*mp.log(GCAR) - b_*mp.log(beta_rad) + (AINV + dv)/n
                marke = f"g_car^{a_} beta_rad^{b_} exp(-(alpha^-1 + {dn})/{n})"
                if (n, a_, b_, dn) == (5, 1, 2, "delta_ph"):
                    enthalten = abs(val - H_ms) < mp.mpf("1e-30")
                if abs(val - H_246) <= tol:
                    treffer_246.append(marke)
                if abs(val - H_EW_LEDGER) <= tol:
                    treffer_HEW.append(marke)
ZAHLEN["anti_numerologie"] = {"Raster": "n = 1..8, g_car-Potenz -2..2, beta_rad-Potenz 0..4, delta in {0, delta_ph}: 400 Formeln",
                              "Toleranz |Delta ln|": mp.nstr(tol, 4), "Treffer bei ln(Mbar/246.22)": treffer_246,
                              "Treffer bei H_EW (Ledger)": treffer_HEW}
check(f"[L] Anti-Numerologie: in einem Raster von 400 Formeln derselben Bauart (n = 1..8, Potenzen von g_car und "
      f"beta_rad, mit und ohne delta_ph) treffen ln(Mbar/246,22) innerhalb der Abweichung der Manuskriptformel "
      f"(|Delta ln| = {mp.nstr(tol, 3)}) {len(treffer_246)} Formel(n): {treffer_246}; den Ledgerwert H_EW trifft in "
      f"derselben Toleranz {len(treffer_HEW)}. "
      + ("Eine weitere, von niemandem gemeinte Formel trifft also ebenso gut: die Nähe zu 246 GeV zeichnet die "
         "Manuskriptformel nicht aus. " if len(treffer_246) > 1 else "Im Raster trifft nur die Manuskriptformel. ")
      + "Ohne Herleitung des Exponenten, mit einer Bezugsgröße, die das Manuskript selbst nicht als v_phys liest, und "
        "mit 1,6 % Abstand ist sie nach der Regel kein Ergebnis; H_EW bleibt offen",
      enthalten and len(treffer_246) >= 1 and "g_car^1 beta_rad^2 exp(-(alpha^-1 + delta_ph)/5)" in treffer_246
      and len(treffer_HEW) == 0)
befund("E/O: Die Hierarchieformel ist im Manuskript behauptet (Lemma mit Beweisskizze), nicht hergeleitet. delta_ph ist "
       "dort 0,5933 (Transportpol); der Wert 2,25 ist eine Rückrechnung aus H_EW und liegt außerhalb des zulässigen "
       "Intervalls. Ein nativer laufender Kopplungsmechanismus mit Koeffizient 10 pi existiert nicht", True)

# =================================================================== Teil 4
say("\nTeil 4: Bilanz und Urteil")
IDENT = [
    "I1 [C] lambda_Sigma = c/R_Vak mit c in {1, 2, 2 pi, 4 pi} (plausibel, nicht belegt); in Lesart A ohne Einfluss auf p_G",
    "I2 [C] g_vac = natives AdS (l = R_Vak/2, R = -48/R_Vak^2, z98/z104), je Volumen; E_Sigma = 0, A_Sigma = 0",
    "I3 [C] U_Sigma = native Nahtgewichte q^(n-1/2), q = e^(-2 pi), L^ = log(q^(1/24) theta_3/eta) (z107); Varianten: "
    "bosonische Determinante, mit Casimir, Manuskriptwert phi_seam",
    "I4 [L] chi-Abhängigkeit von U_Sigma: keine (A, Randdatum nach 029:14), l_Sigma ~ 1/chi (B), l_Sigma ~ chi (C)",
    "I5 [L] Vorzeichen von W wie gedruckt (sonst mit AdS keine Wurzel) und Normierung kappa in {1, 16/pi^2, 1/(4 pi^4)}",
    "I6 [C] N = 48 Komponenten wie im Manuskript (Spinorspur der 48 Weyl-Fermionen wäre 96), f_2 = 1 (Residuum)",
]
bilanz = [
    "[E] Algebra: a_2 = -N R/(192 pi^2), N x Vorfaktor = 48, f_2 = 1 geben Mbar^2 = chi_0^2/(2 pi^2), G chi_0^2 = pi/4; "
    "N = 16, 32, 96 geben 3 pi/4, 3 pi/8, pi/8",
    "[E] spezifiziert: lokale Koeffizienten, Definitionen, Spektralradius phi_seam = 0,0532; nur benannt: E_Sigma, "
    "A_Sigma, g_vac, l_Sigma, U_Sigma(chi), Spektrum von B_Sigma",
    "[E] innere Unstimmigkeiten: Vorzeichen von W, Normierung des Nahtterms (Faktor bis 4 pi^4), R^2-Term ohne Quelle in "
    "a_4(D^2), Weyl-Term mit Geistvorzeichen unter '...' (z41: 1 - 3x/10)",
    f"[C] Lesart A: rho_* = 6/(pi^2 c^2 kappa L^), konvex, eindeutig; p_G = pi^3 kappa L^/24 = {haupt['p_G']:.4f} "
    f"(kappa = 1), {haupt_u['p_G']:.4f} (16/pi^2), {haupt_r['p_G']:.3e} (1/(4 pi^4)); lambda_Sigma fällt heraus",
    f"[C/L] Lesart B: konvex, eindeutig, p_G = {min(pG_B):.3g} bis {max(pG_B):.3g}; Lesart C: nicht konvex, keine oder "
    f"zwei stationäre Punkte, kein globales Minimum",
    "[E] Ziel = Forderung L^ = delta_top e^(-2/alpha) (Vermutung des Manuskripts) = COSMO.LAM.01 in anderen Variablen",
    f"[E] Vorzeichen: gedrucktes W + dS keine Wurzel; gedrucktes W + AdS Wurzel gegen die eigene Einsteinspur; "
    f"konsistentes Potential nur dS",
    f"[E/O] Hierarchie: behauptet; delta_ph = 0,5933 gibt ln(Mbar/v_geo) = {mp.nstr(H_ms, 7)} (v_geo = "
    f"{mp.nstr(v_geo, 5)} GeV), H_EW = 37,1776 verlangte 2,25; kein Laufmechanismus mit b = 10 pi",
    "[O] nicht geprüft: z27 (Universalität der Gitterkopplung der 48 Komponenten), E_Sigma und A_Sigma, das Spektrum "
    "von B_Sigma, eine renormierte Fassung Gamma^ren_TFPT",
]
for z in IDENT + bilanz:
    say("   " + z)
ZAHLEN["identifikationen"] = IDENT
ZAHLEN["bilanz"] = bilanz
urteil = (
    "Z142 rechnet die Planck-Wurzel aus tfpt-42 unter benannten Identifikationen aus (im Bestand lag nach der "
    "Durchsicht R1 keine Rechnung vor); nichts ist befördert, eine vollständige Theorie wird nicht behauptet. Antwort auf die Frage, ob p_G = G/R_Vak^2 aus diesem Funktional mit nativen Eingaben berechenbar ist: "
    "(b) mit Anteilen von (c). (1) Die Algebra stimmt: a_2 = -N R/(192 pi^2) mit N x Vorfaktor = 48 und f_2 = 1 gibt "
    "Mbar^2 = chi_0^2/(2 pi^2) und G chi_0^2 = pi/4; andere Komponentenzahlen geben andere Konstanten (N = 96: pi/8). "
    "Das ist die Definition der Skala chi, kein Wert von G. (2) Das Funktional W ist im Manuskript nicht vollständig "
    "definiert: E_Sigma, A_Sigma, g_vac, die Nahtlänge l_Sigma und vor allem die chi-Abhängigkeit von U_Sigma sind nur "
    "benannt; dazu ist das gedruckte Vorzeichen nicht das Potential der eigenen Wirkung, die Normierung des Nahtterms "
    "ist dreifach (kappa = 1, 16/pi^2, 1/(4 pi^4)), und der R^2-Term hat im Residuum von D^2 keine Quelle (dort steht "
    "der Weyl-Term mit dem Geistvorzeichen aus z41). (3) Mit den Identifikationen I1 bis I6 ist rho_* berechenbar. In "
    "Lesart A (U_Sigma als Randdatum nach 029:14, ohne chi-Abhängigkeit) ist W streng konvex mit genau einer Wurzel rho_* = "
    "6/(pi^2 c^2 kappa L^); für die native Naht L^ = log(q^(1/24) theta_3/eta) = 0,08477 und kappa = 1 ist rho_* = "
    f"{A_tab['1']['rho_*']:.4g}/c^2. Die Gleichung ist aber die Spur der Vakuum-Einsteingleichung: lambda_Sigma fällt "
    f"heraus, p_G = pi^3 kappa L^/24 = {haupt['p_G']:.4f} (kappa = 1), {haupt_u['p_G']:.4f} (Manuskriptnormierung), "
    f"{haupt_r['p_G']:.2e} (reduzierte Normierung), Lambda G = -12 p_G = {haupt['Lambda G (AdS, Lambda = -12/R_Vak^2)']:.3f} "
    f"bis {haupt_r['Lambda G (AdS, Lambda = -12/R_Vak^2)']:.2e}. Sensitivität: p_G ist linear in L^ und kappa, "
    "quadratisch invers in N, unabhängig von c; rho_* geht mit c^-2 und linear mit jedem E_Sigma (Steigung 3/(kappa L^) "
    f"= {float(3/L_F):.0f}). Über alle Lesarten liegt p_G zwischen {min(alle_pG):.1e} und {max(alle_pG):.2g}. Die "
    "Konvexität, die das Manuskript annimmt, gilt für Gamma^(4) nicht (null bzw. winzig konkav) und für den Nahtterm "
    "nur in den Lesarten A und B; in Lesart C gibt es keine oder zwei stationäre Punkte und kein globales Minimum. "
    "(4) Die Route bestätigt das native Ergebnis: Das Vakuum ist planckgroß (R_Vak = "
    f"{1/math.sqrt(max(pG_A_alle)):.1f} bis {1/math.sqrt(min(pG_A_alle)):.0f} Planck-Längen in Lesart A), mindestens {min(math.log10(p) for p in alle_pG) - lg_ziel:.0f} Größenordnungen über dem Ledgerziel "
    "e^(-2/alpha)/(128 pi^3). Das Ziel ist in dieser Sprache exakt die Forderung L^_Sigma = delta_top e^(-2/alpha), die "
    "das Manuskript als Vermutung führt und die seinem eigenen Spektralradius phi_seam = 0,0532 um 122 Größenordnungen "
    "widerspricht; in reduzierter Normierung ist es COSMO.LAM.01. Die native Dämpfung ist pi, nötig wäre 2/alpha. "
    "Vorzeichen: Mit dem nativen AdS gibt nur das gedruckte, mit der eigenen Einsteinspur unverträgliche W eine Wurzel; "
    "das konsistente Potential verlangt dS, also weiter P_cap. Falsifikatoren: z41 trifft den Gravitationszweig (kein "
    "positives TT-Spektrum), nicht die Wurzel selbst; z90 ist verträglich (nur Fermionen, E = -R/4) und sagt, dass "
    "G chi^2 = pi/4 die Skala nur benennt; z27 (Universalität der Kopplung) wurde nicht geprüft. (5) Hierarchieformel: "
    "behauptet, nicht hergeleitet. Das delta_ph des Manuskripts ist 0,5933 und gibt ln(Mbar/v_geo) = "
    f"{mp.nstr(H_ms, 7)} (v_geo = {mp.nstr(v_geo, 4)} GeV, 1,6 % unter 246,22 GeV, das Manuskript liest v_geo nicht als "
    "v_phys); der Ledgerwert H_EW = 37,1776 bräuchte delta_ph = 2,25 außerhalb des zulässigen Intervalls. Ein "
    "Laufmechanismus exp(-alpha^-1/5) bräuchte den irrationalen Koeffizienten 10 pi; keiner der nativen oder "
    f"Standardkoeffizienten gibt ihn, und im Raster von 400 Formeln derselben Bauart treffen {len(treffer_246)} die 246 GeV "
    "ebenso gut. Nach der Anti-Numerologie-Regel ist das kein Ergebnis, H_EW bleibt offen. "
    "Folge für P_G: P_G bleibt Prämisse. Die Route liefert keinen nativen Wert jenseits von O(1) und keine zweite "
    "Länge; sie zeigt, dass p_G und Lambda G dieselbe eine Zahl sind (die Vakuumenergie der Naht in Planck-Einheiten, "
    "gleichwertig zu kappa in z127), und dass die Hierarchie in diesem Funktional nur als gesetzter Nahtwert vorkommt.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z142_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n", encoding="utf-8")
(HERE / "ausgabe_z142.txt").write_text("\n".join(OUT) + "\n", encoding="utf-8")
