"""Z147: Legt die Planck-Wurzel rho_* aus tfpt-42 (G_N lambda_Sigma^2 = pi/(4 rho_*)) die freie Zahl P_G fest?

Leitlinie 'Rekonstruktion statt Auswahl'; H6 und die gelösten Blöcke (Flavour, Massen, alpha, 16, E8) bleiben unberührt,
nichts wird befördert. Statuskürzel: E exakt, C bedingt (Prämisse genannt), L Lesart, N numerisch, O offen.

Frage: Der Bestand enthält G_N lambda_Sigma^2 = pi/(4 rho_*) (tfpt-42; in der Recherche R1 von H6h50 als ungerechnet
geführt). Liefert H6 native Werte für rho_* und lambda_Sigma, sodass G/R_Vak^2 ohne freie Wahl folgt?

Bestand: z136 (Überbestimmungstabelle für G/R_Vak^2, s_0 = ln 2/2 je Naht, P_S, P_eta), z127 (G/R_Vak^2 =
1/(4096 pi kappa)), z137 (l^2/G = (4 c/pi) L N), z112 (Doppelkopie), z96 (l^2/G = 32 pi^3 e^(2/alpha) aus COSMO.LAM.01),
z90 (1/G_ind = 0 beim konformen Leiterfeld), z2 (Verschränkungsdichte je Gitterfläche). z142 ist ein ungeprüfter Entwurf
zum selben Thema; hier ist alles neu gerechnet, nichts daraus übernommen.

Teil 1  Bestand: Definition von rho_* und lambda_Sigma, Herleitung, Annahmen, Nichtzirkularität
  1.1  Quellenzeilen in tfpt-42 (Abschnitte 003, 007, 013, 015, 029, 036, 040) und in der Gesamtsynthese.
  1.2  Kein Zahlenwert von rho_* im Bestand; die chi-Abhängigkeit von U_Sigma und l_Sigma sind nirgends definiert.
  1.3  Algebra: a_2 = -N R/(192 pi^2), N x Vorfaktor = 48 gibt Mbar^2 = chi_0^2/(2 pi^2), G chi_0^2 = pi/4.
  1.4  Profil: das Residuum entspricht dem Wärmekern in D^2; ein Boltzmann-Profil in |D| verdoppelt den Koeffizienten.
  1.5  Nichtzirkularität: weder G = 1 noch eine Flächenregel; aber G wird nur gegen die Skala chi_0 getauscht.
Teil 2  Native Kandidaten
  2.1  Native Nahtdaten neu gerechnet (KMS bei beta = 4 pi, q = e^(-2 pi)): ln Z, s_0, E8-Naht.
  2.2  Identität: rho_* (lambda_Sigma R_Vak)^2 = (chi_0 R_Vak)^2 = pi R_Vak^2/(4 G) = S_dS. Nur das Produkt zählt.
  2.3  Familie M: rho_* als Verhältnis zweier nativer Massen. Familie S: rho_* als native reine Zahl.
  2.4  Definitionstreue Wurzel: Stationarität = Spur der Einsteingleichung, lambda_Sigma fällt heraus.
  2.5  Empfindlichkeit: eine chi-abhängige Naht (Lesart B) ändert den Wert und macht ihn von lambda_Sigma abhängig.
Teil 3  Überbestimmungstabelle aus z136, ergänzt; exakte Übereinstimmungen und ihr Gewicht
Teil 4  Ledger-Ziel COSMO.LAM.01: nötiges rho_*, native Bedeutung
Teil 5  Bilanz und Urteil.
"""
import json
import os
import re
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # deterministische Rundung (einfädig)

import mpmath as mp
import sympy as sp

HERE = Path(__file__).parent
REPO = HERE.parents[2]
T42 = REPO / "documents" / "archive--paper-latex--tfpt-42" / "sections"
SYN = (REPO / "documents" / "tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28" / "sections")
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
mp.mp.dps = 50
PI = mp.pi


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val, text):
    BEFUNDE[name] = val
    say(f"   [Befund: {val}] {text}")


def f6(x, n=6):
    """Zahl mit n geltenden Stellen."""
    return mp.nstr(mp.mpf(x), n)


def lg(x):
    """Zehnerlogarithmus mit zwei Nachkommastellen."""
    return f"{float(mp.log10(mp.mpf(x))):+.2f}"


def N(expr):
    """sympy-Ausdruck als mpmath-Zahl (50 Stellen)."""
    return mp.mpf(sp.N(expr, 50))


def sec_lines(num):
    files = sorted(T42.glob(f"{num}-*.tex"))
    assert len(files) == 1, num
    return files[0].read_text().splitlines()


# ------------------------------------------------------------------------------------------------------------------
say("Teil 1: Bestand -- Definition, Herleitung, Annahmen (tfpt-42, documents/archive--paper-latex--tfpt-42/sections)")

ZITATE = [
    ("003", 20, r"D_+=\gamma_n(\partial_n+B_\Sigma);", "B_Sigma = tangentialer Randoperator im Kragen (Masse)"),
    ("003", 251, r"\tr f\!\left(\frac{D_{\mathrm{rel}}}{\chiseed}\right)", "chi_seed als Spektralskala (Masse)"),
    ("007", 89, r"\deltatop=\frac{3}{256\pi^4},", "delta_top = 3/(256 pi^4)"),
    ("007", 132, r"\phiseam(\alpha)=\frac{1}{6\pi}+\frac{3e^{-2\alpha}}{256\pi^4}", "phi_seam(alpha)"),
    ("011", 368, r"-|D_{\mathrm{ref}}|^{-2s}", "relative Zetafunktion von |D|^(-2s)"),
    ("011", 376, r"\operatorname*{Res}_{s=1}\mathcal Z_{\mathrm{rel}}(s;\sigma),", "E_D := Res_(s=1) Z_rel"),
    ("013", 26, r"-\frac14 R\,\mathbf 1_{48}", "E_rel = -(1/4) R 1_48"),
    ("013", 45, r"\chigeo^0", "a_4-Sektor skaliert mit chi^0"),
    ("013", 108, r"\frac{1}{24\pi^2}R", "E_D enthält -R/(24 pi^2)"),
    ("013", 141, r"-6\chigeo^2\mathcal E_D+\Gamma_D^{(4)}", "Gamma_grav := -6 chi^2 E_D + Gamma_D^(4)"),
    ("013", 161, r"\Mpl^2=\frac{\chigeozero^2}{2\pi^2},", "Mbar^2 = chi_0^2/(2 pi^2)"),
    ("013", 195, r"G_N\chigeozero^2=\frac{\pi}{4}.", "G_N chi_0^2 = pi/4"),
    ("013", 197, r"a rewriting of the Einstein-branch closure, not a second", "ausdrücklich nur eine Umschreibung"),
    ("013", 201, r"G_N\lambda_\Sigma^2=\frac{\pi}{4\rho_\star},", "G_N lambda_Sigma^2 = pi/(4 rho_*)"),
    ("015", 32, r"-\Lambda_{\Seam,\mathrm{ren}}", "Wirkung trägt -Lambda_(Sigma,ren)"),
    ("015", 158, r"\lambda_\Sigma:=\lambda_1^+(|B_\Sigma|)=\chiseed^{-1}", "lambda_Sigma := lambda_1^+(|B_Sigma|)"),
    ("015", 182, r"\rho&:=\frac{\chigeo^2}{\lambda_\Sigma^2},", "rho := chi^2/lambda_Sigma^2"),
    ("015", 196, r"\widehat U_\Sigma(\rho):=U_\Sigma(\lambda_\Sigma\sqrt{\rho}),", "U_Sigma(chi), nur hier"),
    ("015", 199, r"-\log\det\nolimits_{\mathrm{adm}}(1-\widehat U_\Sigma(\rho)),", "L^_Sigma := -log det(1 - U^)"),
    ("015", 206, r"-6\widehat{\mathcal E}_D\rho", "W: -6 E^_D rho"),
    ("015", 208, r"+\rho^2\widehat L_\Sigma(\rho).", "W: + rho^2 L^_Sigma"),
    ("015", 272, r"\partial_\rho\widehat V_{\mathrm{Pl}}^{\mathrm{ren}}(0^+)<0,", "Annahme V'(0+) < 0"),
    ("015", 274, r"\partial_\rho^2\widehat V_{\mathrm{Pl}}^{\mathrm{ren}}(\rho)>0", "Annahme strenge Konvexität"),
    ("015", 276, r"\lim_{\rho\to\infty}", "Annahme V' -> unendlich"),
    ("015", 341, r"\frac{\pi}{4\rho_\star},", "G_N lambda_Sigma^2 = pi/(4 rho_*) (Korollar)"),
    ("015", 364, r"$\Mpl^{-2}=8\pi G_N$", "Standardbeziehung Mbar^-2 = 8 pi G_N"),
    ("015", 367, r"[Schwarzschild de Sitter vacuum branch]", "Vakuumzweig des Manuskripts: dS"),
    ("029", 14, r"U_\Sigma:=\Padm e^{-\ell_\Seam|B_\Seam|}\Padm", "U_Sigma := P e^(-l_Sigma |B_Sigma|) P"),
    ("029", 81, r"\Mplunred^4", "Lambda_IR = M_Pl^4 (unreduziert) x ..."),
    ("029", 82, r"\Bigl[-\log\det\nolimits_{\mathrm{adm}}(1-U_\Sigma)\Bigr]", "... [-log det(1 - U_Sigma)]"),
    ("036", 204, r"\frac{\pi}{4\rho_\star},", "Metrologie-Korollar wiederholt die Relation"),
    ("040", 332, r"\frac{2}{\alpha_\star},", "IR-Fortsetzung S_IR = 2/alpha (Vermutung)"),
    ("040", 336, r"\deltatop e^{-2/\alpha_\star},", "rho_IR = delta_top e^(-2/alpha) (Vermutung)"),
    ("040", 368, r"\num{1.1280419085743317e-123},", "bedingter Zahlenwert Lambda_IR/M_Pl^4"),
]
fehl = []
for num, zeile, nadel, was in ZITATE:
    ls = sec_lines(num)
    ok = zeile <= len(ls) and nadel in ls[zeile - 1]
    if not ok:
        fehl.append(f"{num}:{zeile}")
    say(f"      {num}:{zeile:<4} {was}")
syn = sorted(SYN.glob("004-*.tex"))[0].read_text().splitlines()
syn_ok = "Eine Konvexitätsannahme allein berechnet" in syn[12]
check(f"[E] Quellenbeleg: alle {len(ZITATE)} zitierten Formeln stehen in den angegebenen Zeilen von tfpt-42 "
      f"(fehlend: {fehl}); Gesamtsynthese 004:13 sagt selbst 'Eine Konvexitätsannahme allein berechnet rho_* noch "
      f"nicht' ({syn_ok})", not fehl and syn_ok)

# 1.2 Kein Zahlenwert; nicht definierte Stücke
pat_zahl = re.compile(r"\\rho_\\star\s*(=|\\approx|\\simeq)\s*(\\num\{|[0-9])")
treffer_zahl, dateien_rho = [], 0
for p in sorted((REPO / "documents").rglob("*.tex")):
    txt = p.read_text(errors="replace")
    if "\\rho_\\star" not in txt:
        continue
    dateien_rho += 1
    if pat_zahl.search(txt):
        treffer_zahl.append(str(p.relative_to(REPO)))
u_chi, l_naht = [], []
for p in sorted(T42.glob("*.tex")):
    for i, z in enumerate(p.read_text().splitlines(), 1):
        if re.search(r"U_\\Sigma\(\\lambda_\\Sigma|U_\\Sigma\(\\chi", z):
            u_chi.append(f"{p.name[:3]}:{i}")
        if "\\ell_\\Seam" in z:
            l_naht.append(f"{p.name[:3]}:{i}")
check(f"[E] Kein Zahlenwert: In den {dateien_rho} Quelldateien unter documents/, die rho_* enthalten, steht nirgends "
      f"'rho_* = Zahl' (Treffer: {treffer_zahl}). Eine chi-abhängige Naht U_Sigma(chi) kommt in tfpt-42 nur in "
      f"{u_chi} vor, die Nahtlänge l_Sigma nur in {l_naht} (Definition von U_Sigma und zwei Wiederholungen); beide "
      "sind nirgends bestimmt. Das Funktional, dessen Wurzel rho_* ist, ist im Bestand nicht ausgewertet",
      not treffer_zahl and u_chi == ["015:196"] and l_naht == ["029:14", "029:57", "029:73"] and dateien_rho >= 10)
ZAHLEN["bestand"] = {"Dateien mit rho_*": dateien_rho, "Zahlenwert rho_*": treffer_zahl, "U_Sigma(chi)": u_chi,
                     "l_Sigma": l_naht}

# 1.3 Algebra
Rs, Ns, chi, lam, rho, Gs = sp.symbols("R N chi lambda rho G", positive=True)
a2 = sp.simplify(Ns * (-Rs / 4 + Rs / 6) / (16 * sp.pi**2))          # (4 pi)^-2 tr(E + R/6), E = -R/4
a2_8 = a2.subs(Ns, 8)
Mbar2 = sp.solve(sp.Eq(sp.Symbol("M") / 2 * Rs, -6 * chi**2 * a2_8), sp.Symbol("M"))[0]
Gchi2 = sp.simplify(1 / (8 * sp.pi * Mbar2) * chi**2)
Mbar2_N = sp.simplify(-2 * chi**2 * a2 / Rs)                          # N Komponenten insgesamt, Vorfaktor 1
G_N_allg = sp.simplify(1 / (8 * sp.pi * Mbar2_N))
werte_N = {n: sp.simplify(G_N_allg.subs(Ns, n) * chi**2) for n in (16, 32, 48, 96)}
check(f"[E] Algebra (sympy): a_2 = (4 pi)^-2 tr(E + R/6) mit E = -R/4 gibt {a2} je N Komponenten; der "
      f"Manuskriptkoeffizient -R/(24 pi^2) (013:108) ist N = 8 ({sp.simplify(a2_8 + Rs / (24 * sp.pi**2)) == 0}), "
      f"mit dem Vorfaktor 6 (013:141) also 6 x 8 = 48 Komponenten. Daraus Mbar^2 = {Mbar2} und mit Mbar^-2 = 8 pi G "
      f"G chi_0^2 = {Gchi2}. Andere Komponentenzahlen: N = 16, 32, 96 geben "
      f"{[str(werte_N[n]) for n in (16, 32, 96)]}; N = 96 wäre die Zählung der 48 Weyl-Fermionen mit ihrer "
      "Spinorspur (je 2)",
      sp.simplify(a2 + Ns * Rs / (192 * sp.pi**2)) == 0 and sp.simplify(Mbar2 - chi**2 / (2 * sp.pi**2)) == 0
      and Gchi2 == sp.pi / 4 and werte_N[48] == sp.pi / 4 and werte_N[96] == sp.pi / 8)
Mb2 = chi**2 / (2 * sp.pi**2)
rel = {
    "Mbar^2/lambda^2 = rho/(2 pi^2)": sp.simplify((Mb2 / lam**2).subs(chi, lam * sp.sqrt(rho)) - rho / (2 * sp.pi**2)),
    "G lambda^2 = pi/(4 rho)": sp.simplify((lam**2 / (8 * sp.pi * Mb2)).subs(chi, lam * sp.sqrt(rho))
                                           - sp.pi / (4 * rho)),
    "m_P/lambda = 2 sqrt(rho)/sqrt(pi)": sp.simplify((sp.sqrt(8 * sp.pi * Mb2) / lam).subs(chi, lam * sp.sqrt(rho))
                                                     - 2 * sp.sqrt(rho) / sp.sqrt(sp.pi)),
}
check("[E] Die Relation ist eine Umschreibung: Mit rho := chi_0^2/lambda_Sigma^2 folgen Mbar^2/lambda^2 = "
      "rho/(2 pi^2), G lambda^2 = pi/(4 rho) und m_P/lambda = 2 sqrt(rho)/sqrt(pi) identisch aus G chi_0^2 = pi/4 "
      f"(Reste {[str(v) for v in rel.values()]}). Sie enthält genau eine Aussage, G chi_0^2 = pi/4; rho_* ist der "
      "Name für (chi_0/lambda_Sigma)^2. Das Manuskript sagt das selbst (013:197: 'a rewriting ..., not a second "
      "independent prediction')", all(v == 0 for v in rel.values()))

# 1.4 Profil
p_s = sp.Symbol("p", positive=True)
faktor = sp.gamma(p_s + sp.Rational(1, 2)) * 4**(p_s + sp.Rational(1, 2)) / (2 * sp.sqrt(sp.pi)) / sp.gamma(p_s)
fak1 = sp.simplify(faktor.subs(p_s, 1))                               # Koeffizient von t^-2 je a_2
tt, la = mp.mpf("0.7"), mp.mpf("2.3")
sub_num = mp.quad(lambda s: tt / (2 * mp.sqrt(PI)) * s**(-1.5) * mp.exp(-tt**2 / (4 * s)) * mp.exp(-s * la),
                  [0, 0.05, 0.5, 5, mp.inf])
sub_res = abs(sub_num - mp.exp(-tt * mp.sqrt(la)))
pot_num = mp.quad(lambda s: tt / (2 * mp.sqrt(PI)) * s**(-2.5) * mp.exp(-tt**2 / (4 * s)), [0, 0.05, 0.5, 5, mp.inf])
pot_res = abs(pot_num - 2 / tt**2)
f2_heat = sp.integrate(sp.exp(-sp.Symbol("u", positive=True)**2) * sp.Symbol("u", positive=True),
                       (sp.Symbol("u", positive=True), 0, sp.oo))
f2_boltz = sp.integrate(sp.exp(-sp.Symbol("u", positive=True)) * sp.Symbol("u", positive=True),
                        (sp.Symbol("u", positive=True), 0, sp.oo))
check("[E] Profilabhängigkeit: Das profilfreie Residuum Res_(s=1) der Zetafunktion von D^2 (011:368-376) ist a_2 "
      "selbst und "
      f"entspricht dem Wärmekernprofil e^(-D^2/chi^2) (2 f_2 = {2 * f2_heat}). Ein Boltzmann-Profil e^(-|D|/chi) hat "
      f"2 f_2 = {2 * f2_boltz}: Über die Subordination e^(-t sqrt(x)) = int K_t(s) e^(-s x) ds (Rest {f6(sub_res, 2)}) "
      f"wird aus s^-1 der Term {fak1} t^-2 (numerisch Rest {f6(pot_res, 2)}). Dann ist G chi_0^2 = pi/8 statt pi/4. "
      "Die Konstante pi/4 gilt genau für N x Vorfaktor = 48 und dieses Profil",
      2 * f2_heat == 1 and 2 * f2_boltz == 2 and fak1 == 2 and sub_res < mp.mpf("1e-20") and pot_res < mp.mpf("1e-20"))

# 1.5 Nichtzirkularität
muster = re.compile(r"G_N\s*=\s*1[^0-9]|Bekenstein|A/\(4|area law|Bisognano|modular flow", re.I)
nc = []
for num in ("013", "015", "036"):
    for i, z in enumerate(sec_lines(num), 1):
        if muster.search(z):
            nc.append(f"{num}:{i}")
dim_konflikt = ("\\chiseed^{-1}" in sec_lines("015")[157]) and ("{\\chiseed}" in sec_lines("003")[250])
check("[E] Nichtzirkularität (GRAV.NONCIRCULAR.01): Die Herleitung in 013/015/036 setzt weder G = 1 noch benutzt sie "
      f"eine Flächenregel, Bisognano-Wichmann oder einen Modulfluss (Treffer: {nc}). Sie ist eine induzierte "
      "Einsteinkopplung nach Sakharov: 48 Komponenten mit E = -R/4, profilfreies Residuum, dazu die Definition "
      "Mbar^-2 = 8 pi G (015:364). Im Sinn der Brandmauer ist sie also nicht zirkulär. Sie setzt aber einen "
      "4D-Dirac-Operator mit Wärmekernentwicklung voraus (Import nach DIMENSION.UPLIFT.FIREWALL.01) und tauscht G "
      "nur gegen die Skala chi_0", not nc)
check("[E] Unstimmigkeit im Bestand: 015:158 setzt lambda_Sigma = chi_seed^-1, während 003:251 chi_seed als "
      "Spektralskala (Masse) in tr f(D/chi_seed) führt und 015:182 lambda_Sigma als Masse braucht (rho dimensionslos). "
      "Hier wird lambda_Sigma als Masse gelesen, lambda_1^+(|B_Sigma|) des Randoperators aus D_+ = "
      "gamma_n (d_n + B_Sigma) (003:20)", dim_konflikt)
befund("Bestand", "ja",
       "E: rho_* ist in tfpt-42 als eindeutige stationäre Wurzel eines randnormierten Potentials definiert, unter drei "
       "Annahmen (015:270-277), die dort nicht geprüft werden. Ausgewertet ist das Potential nirgends; die Relation "
       "G lambda_Sigma^2 = pi/(4 rho_*) ist G chi_0^2 = pi/4 in anderen Variablen. Annahmen: 48 Komponenten, "
       "Residuenprofil, Existenz der stationären Skala chi_0, ein Randoperator B_Sigma mit Spektrum. Weder G = 1 noch "
       "eine Flächenregel geht ein")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil 2: Native Kandidaten")
q = mp.exp(-2 * PI)
NMAX = 60
lnZ_f = 2 * mp.fsum(mp.log(1 + q**(n - mp.mpf(1) / 2)) for n in range(1, NMAX))
lnZ_b = -2 * mp.fsum(mp.log(1 - q**(n - mp.mpf(1) / 2)) for n in range(1, NMAX))


def h_ferm(eps):
    return mp.log(1 + mp.exp(-eps)) + eps / (mp.exp(eps) + 1)


s_quarter = mp.fsum(h_ferm(4 * PI * abs(mp.mpf(k) + mp.mpf(1) / 4)) for k in range(-NMAX, NMAX))
s_quarter34 = mp.fsum(h_ferm(4 * PI * abs(mp.mpf(k) + mp.mpf(3) / 4)) for k in range(-NMAX, NMAX))
s_naht = s_quarter + s_quarter34
bE = s_naht - lnZ_f
spek14 = sorted(abs(mp.mpf(k) + mp.mpf(1) / 4) for k in range(-6, 6))[:4]


def lnZ_voll(y):
    qq = mp.exp(-2 * PI * y)
    return PI * y / 12 + 2 * mp.fsum(mp.log(1 + qq**(n - mp.mpf(1) / 2)) for n in range(1, NMAX))


def lnZ_e8(y):
    qq = mp.exp(-2 * PI * y)
    e4 = 1 + 240 * mp.fsum(n**3 * qq**n / (1 - qq**n) for n in range(1, NMAX))
    eta8 = 8 * (-2 * PI * y / 24 + mp.fsum(mp.log(1 - qq**n) for n in range(1, NMAX)))
    return mp.log(e4) - eta8


theta_gamma = 2 * PI / (mp.gamma(mp.mpf(1) / 4) * mp.gamma(mp.mpf(3) / 4))      # theta_3(i)/eta(i)
d_voll = mp.diff(lnZ_voll, 1)
d_e8 = mp.diff(lnZ_e8, 1)
s_e8 = lnZ_e8(1)
res_nativ = {
    "ln Z_ferm - (ln 2/2 - pi/12)": abs(lnZ_f - (mp.log(2) / 2 - PI / 12)),
    "-ln det(1 - U) - (pi/12 - ln 2/4)": abs(lnZ_b - (PI / 12 - mp.log(2) / 4)),
    "s_Quarter - ln 2/4": abs(s_quarter - mp.log(2) / 4),
    "s_Naht - ln 2/2": abs(s_naht - mp.log(2) / 2),
    "beta<E> - pi/12": abs(bE - PI / 12),
    "theta_3(i)/eta(i) - sqrt 2": abs(theta_gamma - mp.sqrt(2)),
    "ln Z_voll(i) - ln 2/2": abs(lnZ_voll(1) - mp.log(2) / 2),
    "S_E8 - ln 12": abs(s_e8 - mp.log(12)),
}
check("[E] Native Nahtdaten, neu gerechnet (KMS bei beta = 4 pi, Einteilchenenergien 4 pi |x|, x in Z + 1/4 bzw. "
      f"Z + 3/4; tiefste |x| = {[str(sp.nsimplify(float(v))) for v in spek14]}, also Gewichte q^(n-1/2), q = e^(-2 pi), "
      f"je Naht zweifach): ln Z_ferm = ln 2/2 - pi/12 = {f6(lnZ_f, 10)}, wörtliche bosonische Determinante "
      f"-ln det(1 - U) = pi/12 - ln 2/4 = {f6(lnZ_b, 10)}, beta<E> = pi/12, s_0 = ln 2/4 je Quarter und ln 2/2 = "
      f"{f6(s_naht, 10)} je Naht (wie z136). Geschlossen: mit Casimir ist Z = theta_3(i)/eta(i) = "
      "2 pi/(Gamma(1/4) Gamma(3/4)) = sqrt 2 und d ln Z/dy = 0 am selbstdualen Punkt, also s_0 = ln Z_voll = ln 2/2. "
      f"E8-Naht: Z = E_4/eta^8 = 12, S = ln 12. Größter Rest {f6(max(res_nativ.values()), 2)}; Ableitungen "
      f"{f6(abs(d_voll), 2)}, {f6(abs(d_e8), 2)}",
      max(res_nativ.values()) < mp.mpf("1e-40") and abs(d_voll) < mp.mpf("1e-30") and abs(d_e8) < mp.mpf("1e-30"))
ZAHLEN["nativ"] = {"ln Z_ferm": f6(lnZ_f, 20), "-ln det(1 - U) bosonisch": f6(lnZ_b, 20),
                   "s_0 Quarter": f6(s_quarter, 20), "s_0 Naht": f6(s_naht, 20), "S E8-Naht": f6(s_e8, 20),
                   "beta<E>": f6(bE, 20), "Reste": {k: f6(v, 3) for k, v in sorted(res_nativ.items())}}

# 2.2 Identität rho_* c^2 = S_dS
Rv, cs, ks, Ls = sp.symbols("R_Vak c kappa L", positive=True)
G_von_chi = sp.pi / (4 * chi**2)
S_dS = sp.simplify(sp.pi * (Rv / 2)**2 / G_von_chi)                    # 4 pi l^2/(4 G), l = R_Vak/2
S_flach = sp.simplify(4 * sp.pi * Rv**2 / (4 * G_von_chi))
rho_c2 = (chi / lam)**2 * (lam * Rv)**2
S_dS_N = sp.simplify(sp.pi * (Rv / 2)**2 / G_N_allg)
check("[E] Identität: rho_* (lambda_Sigma R_Vak)^2 = (chi_0 R_Vak)^2 = pi R_Vak^2/(4 G) = S_dS, die "
      "Bekenstein-Hawking-Entropie des dS-Horizonts vom Radius l = R_Vak/2 (Fläche pi R_Vak^2, z136); im flachen "
      f"Rahmen (Fläche 4 pi R_Vak^2) ist A/(4G) = 4 rho_* c^2 (Reste {sp.simplify(S_dS - rho_c2)}, "
      f"{sp.simplify(S_flach - 4 * rho_c2)}). In G/R_Vak^2 = pi/(4 rho_* c^2) gehen rho_* und lambda_Sigma nur als "
      "Produkt ein: eine Zahl, nicht zwei. Die Gleichheit der Faktoren pi/4 hängt an der Normierung 48 und an "
      f"l = R_Vak/2: allgemein S_dS = {sp.simplify(S_dS_N / (chi * Rv)**2)} (chi_0 R_Vak)^2, bei N = 96 das Doppelte",
      sp.simplify(S_dS - rho_c2) == 0 and sp.simplify(S_flach - 4 * rho_c2) == 0
      and sp.simplify(S_dS_N.subs(Ns, 48) - (chi * Rv)**2) == 0)

# 2.3 Kandidaten
C_KAND = [
    ("c1", "Viertelmode alpha = 1/4, tiefster Eigenwert des Nahtgenerators (P_B)", sp.Rational(1, 4)),
    ("c2", "NS-Mode (L_0^NS = 2 L_0)", sp.Rational(1, 2)),
    ("c3", "1/R_Vak (Leiterstufe Delta = 1)", sp.Integer(1)),
    ("c4", "MM-Länge l = R_Vak/2", sp.Integer(2)),
    ("c5", "KMS-Länge beta R_Vak = 4 pi R_Vak", 1 / (4 * sp.pi)),
    ("c6", "Nahtumfang 2 pi R_Vak", 1 / (2 * sp.pi)),
    ("c7", "Markenabstand pi R_Vak/2", 2 / sp.pi),
]
kappa_z2 = json.loads((HERE / "z2_ergebnis.json").read_text())["kappa_naive_dirac"]
R_KAND = [
    ("r1", "s_0 Quarter = ln 2/4", sp.log(2) / 4, "E"),
    ("r2", "s_0 Naht = ln 2/2", sp.log(2) / 2, "E"),
    ("r3", "S E8-Naht = ln 12", sp.log(12), "E"),
    ("r4", "ln Z_ferm = ln 2/2 - pi/12", sp.log(2) / 2 - sp.pi / 12, "E"),
    ("r5", "-ln det(1 - U) = pi/12 - ln 2/4", sp.pi / 12 - sp.log(2) / 4, "E"),
    ("r6", "kappa je Gitterfläche (z2) [N]", sp.Float(round(kappa_z2, 6), 12), "N"),
    ("r7", "beta = 4 pi", 4 * sp.pi, "E"),
    ("r8", "-ln q = 2 pi", 2 * sp.pi, "E"),
]
ell_lam = sp.simplify(4 * sp.pi * C_KAND[0][2])                       # l_Sigma lambda_Sigma mit l_Sigma = beta R_Vak
check("[C] lambda_Sigma nativ (Prämisse P_B: B_Sigma ist der Einteilchenerzeuger der Naht, l_Sigma = beta R_Vak): "
      f"lambda_1^+ = 1/(4 R_Vak), l_Sigma lambda_Sigma = {ell_lam}, tiefstes Gewicht e^(-pi) = q^(1/2) "
      f"({f6(abs(mp.exp(-PI) - mp.sqrt(q)), 2)}). Das ist der einzige Kandidat, der der Definition lambda_1^+(|B_Sigma|) "
      "folgt; die sechs anderen (NS-Mode, 1/R_Vak, MM-Länge, KMS-Länge, Umfang, Markenabstand) sind Lesarten",
      ell_lam == sp.pi and abs(mp.exp(-PI) - mp.sqrt(q)) < mp.mpf("1e-45"))

say("   Familie M: rho_* = (chi_0/lambda_Sigma)^2 als Verhältnis zweier nativer Massen; lambda_Sigma fällt heraus,")
say("   G/R_Vak^2 = pi/(4 (chi_0 R_Vak)^2):")
fam_M = []
for key, name, x in C_KAND:
    val = sp.simplify(sp.pi / (4 * x**2))
    fam_M.append((f"M/{key}", f"chi_0 R_Vak = {x}: {name}", val))
    say(f"      chi_0 R_Vak = {str(x):<9} {name:<70} G/R_Vak^2 = {str(val):<10} = {f6(N(val))}")
vals_M = [N(v) for _, _, v in fam_M]
check(f"[L] Familie M (7 Werte, jede Wahl von chi_0 eine Prämisse P_chi): G/R_Vak^2 von {f6(min(vals_M))} bis "
      f"{f6(max(vals_M))}, also R_Vak/l_P = {f6(1 / mp.sqrt(max(vals_M)), 3)} bis {f6(1 / mp.sqrt(min(vals_M)), 3)}. "
      "Jede native Masse ist eine Zahl mal 1/R_Vak (konforme Kinematik, eine Länge, z136), deshalb ist jedes so "
      "gebildete rho_* von der Ordnung 1 und das Vakuum planckgroß. Ein natives chi_0 mit eigenem Rang gibt es nicht: "
      "Das Spektrum der Leiter ist unbeschränkt, das konforme Leiterfeld induziert 1/G = 0 (z90)",
      len(set(str(v) for _, _, v in fam_M)) == 7 and max(vals_M) < 200 and min(vals_M) > mp.mpf("0.1"))

say("   Familie S: rho_* = native reine Zahl r, lambda_Sigma = c/R_Vak; G/R_Vak^2 = pi/(4 r c^2):")
kopf = "      " + " " * 36 + " | ".join(f"c = {str(c):<9}" for _, _, c in C_KAND)
say(kopf)
fam_S = []
for rk, rname, r, rst in R_KAND:
    zeile = []
    for ck, cname, c in C_KAND:
        val = sp.simplify(sp.pi / (4 * r * c**2))
        fam_S.append((f"S/{rk}/{ck}", f"rho_* = {rname}; lambda_Sigma R_Vak = {c}", val, r * c**2, rst))
        zeile.append(f"{f6(N(val), 5):<13}")
    say(f"      {rname:<36}" + " | ".join(zeile))
vals_S = [N(v) for _, _, v, _, _ in fam_S]
produkte = sorted({f6(N(p_), 30) for _, _, _, p_, _ in fam_S})
check(f"[L] Familie S ({len(fam_S)} Paarungen, jede eine Prämisse P_rho und eine Prämisse P_lambda): G/R_Vak^2 von "
      f"{f6(min(vals_S))} bis {f6(max(vals_S))} (10^{lg(min(vals_S))} bis 10^{lg(max(vals_S))}; R_Vak/l_P = "
      f"{f6(1 / mp.sqrt(max(vals_S)), 3)} bis {f6(1 / mp.sqrt(min(vals_S)), 3)}); {len(produkte)} "
      "verschiedene Produkte r c^2, keine zwei Paarungen fallen zusammen. Die Definition zeichnet nur die Spalte "
      f"c = 1/4 aus (P_B); dort {f6(N(fam_S[0][2]))} (Quarter), {f6(N(fam_S[7][2]))} (Naht), {f6(N(fam_S[14][2]))} "
      "(E8-Naht). Eine Regel, die eine Zeile auswählt, gibt es nicht",
      len(fam_S) == 56 and len(produkte) == 56 and min(vals_S) > mp.mpf("1e-3") and max(vals_S) < mp.mpf("1e4"))
eta_s, a_s, k2_s = sp.symbols("eta a kappa_2", positive=True)
X_dichte = sp.simplify(eta_s * sp.pi * Rv**2)                          # S_dS = eta mal Fläche pi R_Vak^2
pG_dichte = sp.simplify((sp.pi / (4 * X_dichte)).subs(eta_s, k2_s / a_s**2))
a_lP = 2 * mp.sqrt(mp.mpf(kappa_z2))
check("[E/L] Dichte-Lesart: Eine Entropie je Fläche eta gibt rho_* c^2 = S_dS = pi eta R_Vak^2, also G = 1/(4 eta), "
      "die Bekenstein-Dichte selbst. Für z2 (eta = kappa/a^2 je Gitterfläche) ist G/R_Vak^2 = "
      f"{pG_dichte}, mit kappa = {f6(kappa_z2, 6)} also {f6(1 / (4 * mp.mpf(kappa_z2)), 6)} "
      f"(a/R_Vak)^2 und a/l_P = 2 sqrt(kappa) = {f6(a_lP, 5)} (z2: 1,77). Der Gitterabstand a ist eine dritte, nicht "
      "native Länge (Prämisse D aus z2 und a/R_Vak). Der Tabellenwert 0,9985 bei c = 1 ist pi/(4 kappa), die schon in "
      f"z2 verworfene Nähe kappa/(pi/4) = {f6(mp.mpf(kappa_z2) / (PI / 4), 6)}",
      sp.simplify(pG_dichte - a_s**2 / (4 * k2_s * Rv**2)) == 0 and abs(a_lP - mp.mpf("1.7738")) < mp.mpf("1e-3")
      and abs(mp.mpf(kappa_z2) / (PI / 4) - 1) > mp.mpf("1e-3"))
ZAHLEN["familie_M"] = {k: {"Lesart": n_, "G/R_Vak^2": str(v), "Wert": f6(N(v), 12)} for k, n_, v in fam_M}
ZAHLEN["familie_S"] = {k: {"Lesart": n_, "G/R_Vak^2": str(v), "Wert": f6(N(v), 12), "Status": st}
                       for k, n_, v, _, st in fam_S}

# 2.4 Definitionstreue Wurzel
sg = sp.Symbol("sigma")                                                # Vorzeichen der Skalarkrümmung
Rhat = sg * 48 / cs**2                                                 # R/lambda^2 mit |R| = 48/R_Vak^2
W_druck = rho * Rhat / (4 * sp.pi**2) + ks * rho**2 * Ls               # -6 E^_D rho + kappa rho^2 L^
V_kons = -rho * Rhat / (4 * sp.pi**2) + ks * rho**2 * Ls               # -(Mbar^2 R/2 - Lambda_Sigma) je Volumen
w_druck = sp.solve(sp.diff(W_druck, rho), rho)[0]
w_kons = sp.solve(sp.diff(V_kons, rho), rho)[0]
rho_stern_A = sp.simplify(w_druck.subs(sg, -1))
pG_A = sp.simplify(sp.pi / (4 * rho_stern_A * cs**2))
R_einstein = sp.simplify(4 * ks * chi**4 * Ls / (chi**2 / (2 * sp.pi**2)))      # R = 4 Lambda_Sigma/Mbar^2
R_stat = sp.solve(sp.diff(chi**2 * Rs / (4 * sp.pi**2) - ks * chi**4 * Ls, chi), Rs)[0]
check("[E] Die Stationarität in chi bei festem g ist die Spur der Vakuum-Einsteingleichung: Für die Wirkungsdichte "
      "Mbar^2 R/2 - Lambda_Sigma (015:27-45) mit Mbar^2 = chi^2/(2 pi^2) und Lambda_Sigma = kappa chi^4 L^ gibt "
      f"d/dchi = 0 R = {R_stat}, die Spur R = 4 Lambda_Sigma/Mbar^2 dasselbe ({sp.simplify(R_stat - R_einstein)}). "
      "Die Wurzelgleichung fügt der Einsteingleichung nichts hinzu",
      sp.simplify(R_stat - R_einstein) == 0)
check(f"[C] Lesart A (U_Sigma hängt nicht von chi ab; |R| = 48/R_Vak^2 aus Lambda = -12/R_Vak^2, z98/z118; "
      f"E_Sigma = A_Sigma = 0, Higgs-Anteil vernachlässigt; a_4-Sektor chi^0): rho_* = {rho_stern_A}, also "
      f"rho_* c^2 = S_dS = {sp.simplify(rho_stern_A * cs**2)} und G/R_Vak^2 = {pG_A}, |Lambda| G = "
      f"{sp.simplify(12 * pG_A)}. lambda_Sigma fällt heraus (d/dc = {sp.diff(pG_A, cs)}); rho_* allein ist nur die "
      "Umrechnung c -> rho_* ~ c^-2. Die drei Annahmen des Satzes (015:270-277) gelten hier: W'' = 2 kappa L^ > 0, "
      "W'(0+) < 0 genau für R < 0, W' -> unendlich",
      sp.simplify(rho_stern_A - 6 / (sp.pi**2 * cs**2 * ks * Ls)) == 0
      and sp.simplify(pG_A - sp.pi**3 * ks * Ls / 24) == 0 and sp.diff(pG_A, cs) == 0
      and sp.simplify(sp.diff(W_druck, rho, 2) - 2 * ks * Ls) == 0
      and sp.simplify(sp.diff(W_druck, rho).subs(rho, 0) - sg * 12 / (sp.pi**2 * cs**2)) == 0
      and sp.limit(sp.diff(W_druck, rho).subs(sg, -1), rho, sp.oo) == sp.oo)
check("[E] Vorzeichen: Im gedruckten W (015:205-208) haben Einstein- und Nahtterm dasselbe relative Vorzeichen, in "
      f"der Wirkung (015:27-45) entgegengesetzte. Wurzeln: gedruckt rho = {w_druck}, konsistent rho = {w_kons} "
      "(sigma = Vorzeichen von R). Das gedruckte W hat eine positive Wurzel nur für AdS (sigma = -1), das konsistente "
      "nur für dS (sigma = +1); der Betrag ist gleich. Mit dem eigenen Vakuum des Manuskripts (Schwarzschild-de "
      "Sitter, 015:367) hat das gedruckte W keine Wurzel; mit dem nativen AdS hat es eine, die der Spur der eigenen "
      "Einsteingleichung im Vorzeichen widerspricht",
      sp.simplify(w_druck + w_kons) == 0 and w_druck.subs({sg: -1, cs: 1, ks: 1, Ls: 1}) > 0
      and w_kons.subs({sg: 1, cs: 1, ks: 1, Ls: 1}) > 0 and w_druck.subs({sg: 1, cs: 1, ks: 1, Ls: 1}) < 0)

ainv = mp.mpf("137.0359992168")                                        # EM.FP.01 (gelöster Block, nur benutzt)
alpha = 1 / ainv
delta_top = 3 / (256 * PI**4)
qa = delta_top * mp.exp(-2 * alpha)
phi_seam = 1 / (6 * PI) + qa * (1 - qa)**(-mp.mpf(5) / 4)
L_ms = -mp.log(1 - phi_seam)
L_KAND = [
    ("L1", "ln Z_ferm = ln 2/2 - pi/12 (native Gewichte, fermionisch)", sp.log(2) / 2 - sp.pi / 12),
    ("L2", "-ln det(1 - U) = pi/12 - ln 2/4 (wörtliche Determinante)", sp.pi / 12 - sp.log(2) / 4),
    ("L3", "ln Z_voll = s_0 = ln 2/2 (mit Casimir)", sp.log(2) / 2),
    ("L4", "Manuskript, Satzschicht: -ln(1 - phi_seam(alpha)), Rang eins", sp.Float(mp.nstr(L_ms, 15), 15)),
]
K_KAND = [
    ("k1", "kappa = 1 (gedruckt, rho^2 L^)", sp.Integer(1)),
    ("k2", "kappa = 16/pi^2 (Lambda_IR = M_Pl^4 L^, unreduziert, 029:79-82)", 16 / sp.pi**2),
    ("k3", "kappa = 1/(4 pi^4) (reduziert, Mbar^4 L^)", 1 / (4 * sp.pi**4)),
]
kap_unred = sp.simplify((8 * sp.pi)**2 * (1 / (2 * sp.pi**2))**2)      # M_Pl^4 = (8 pi)^2 Mbar^4, Mbar^2 = chi^2/2pi^2
say("   Wurzelroute, Lesart A: G/R_Vak^2 = pi^3 kappa L^/24 (lambda_Sigma beliebig):")
fam_W = []
for lk, lname, Lv in L_KAND:
    for kk, kname, kv in K_KAND:
        val = sp.simplify(sp.pi**3 * kv * Lv / 24)
        fam_W.append((f"W/{lk}/{kk}", f"{lname}; {kname}", val))
        say(f"      {lname:<62} {kname.split(' (')[0]:<18} G/R_Vak^2 = {f6(N(val))}   S_dS = {f6(PI / (4 * N(val)))}")
vals_W = [N(v) for _, _, v in fam_W]
check(f"[C] Wurzelroute mit nativen Nahtgewichten (Prämissenbündel P_rho*: P_48, P_U, P_kappa, P_vac): 12 Werte, "
      f"G/R_Vak^2 von {f6(min(vals_W))} bis {f6(max(vals_W))} (R_Vak/l_P = {f6(1 / mp.sqrt(max(vals_W)), 3)} bis "
      f"{f6(1 / mp.sqrt(min(vals_W)), 3)}). Die Normierung des Nahtterms steht im Manuskript dreifach (kappa = 1, "
      f"{kap_unred}, 1/(4 pi^4); Spanne {f6(N(kap_unred * 4 * sp.pi**4))}). Auch die Satzschicht des Manuskripts selbst "
      f"(Spektralradius phi_seam(alpha) = {f6(phi_seam, 10)}, L^ >= {f6(L_ms, 8)}) gibt ein planckgroßes Vakuum",
      kap_unred == 16 / sp.pi**2 and min(vals_W) > mp.mpf("1e-4") and max(vals_W) < 1 and len(fam_W) == 12)
ZAHLEN["wurzelroute"] = {k: {"Lesart": n_, "G/R_Vak^2": str(v), "Wert": f6(N(v), 12)} for k, n_, v in fam_W}
ZAHLEN["wurzelroute"]["Formel"] = {"rho_* (Lesart A)": str(rho_stern_A), "G/R_Vak^2": str(pG_A),
                                   "phi_seam(alpha)": f6(phi_seam, 15), "L^ Manuskript": f6(L_ms, 15)}

# Wurzelroute gegen Jacobson: ganze Nahtzahl?
s0 = mp.log(2) / 2
n_tab = {}
for lk, lname, Lv in L_KAND[:3]:
    for kk, kname, kv in K_KAND:
        n2 = 6 / (PI**2 * N(kv) * N(Lv) * s0)                          # N^2 aus pi^3 kappa N L/24 = pi/(4 N s_0)
        n_tab[f"{lk}/{kk}"] = mp.sqrt(n2)
abst = {k: abs(v - mp.nint(v)) for k, v in n_tab.items()}
check("[E negativ, gescopt] Wurzelroute gegen Jacobson (z136, dS-Rahmen): Beide zugleich verlangen S_dS L^ = "
      "6/(pi^2 kappa). Mit additiven Nähten (P_ord: L^ = N L_1, S_dS = N s_0) wäre N = "
      f"{ {k: f6(v, 5) for k, v in sorted(n_tab.items())} }: in keiner der neun Fassungen eine ganze Zahl "
      f"(kleinster Abstand {f6(min(abst.values()), 3)}). Die beiden bedingten Routen widersprechen sich für jede "
      "ganze Nahtzahl; G/R_Vak^2 wächst in der einen mit dem Nahtinhalt und fällt in der anderen",
      min(abst.values()) > mp.mpf("0.05"))
ZAHLEN["nahtzahl_wurzel_gegen_jacobson"] = {k: f6(v, 10) for k, v in sorted(n_tab.items())}

# 2.5 Empfindlichkeit: Lesart B
mp.mp.dps = 25


def L_B(r):
    t = PI / mp.sqrt(r)
    nmax = int(40 * mp.sqrt(r) / float(PI)) + 12
    return 2 * mp.fsum(mp.log(1 + mp.exp(-t * (2 * n - 1))) for n in range(1, nmax))


def dL_B(r):
    t = PI / mp.sqrt(r)
    nmax = int(40 * mp.sqrt(r) / float(PI)) + 12
    return 2 * mp.fsum((2 * n - 1) * PI / (2 * r**mp.mpf("1.5")) / (mp.exp(t * (2 * n - 1)) + 1)
                       for n in range(1, nmax))


def g_B(r):
    return 2 * r * L_B(r) + r**2 * dL_B(r)


gitter = [mp.mpf(10)**(mp.mpf(k) / 4) for k in range(-12, 13)]
gw = [g_B(r) for r in gitter]
monoton = all(gw[i + 1] > gw[i] for i in range(len(gw) - 1))
les_B = {}
for key, name, c in C_KAND:
    a = 12 / (PI**2 * N(c)**2)
    lo, hi = mp.mpf("1e-3"), mp.mpf("1e3")
    for _ in range(90):
        mid = mp.sqrt(lo * hi)
        if g_B(mid) < a:
            lo = mid
        else:
            hi = mid
    r_st = mp.sqrt(lo * hi)
    les_B[key] = (r_st, PI / (4 * r_st * N(c)**2))
mp.mp.dps = 50
pB = [v[1] for v in les_B.values()]
check("[L] Empfindlichkeit (Lesart B: Nahtlänge in Einheiten der Skala, l_Sigma ~ 1/chi, natives Gewicht bei "
      f"rho = 1, kappa = 1, fermionisch): d/drho [rho^2 L^] ist monoton ({monoton}), je c genau eine Wurzel; "
      f"rho_* = { {k: f6(v[0], 6) for k, v in sorted(les_B.items())} }, G/R_Vak^2 = "
      f"{ {k: f6(v[1], 5) for k, v in sorted(les_B.items())} }. Hier hängt der Wert von lambda_Sigma ab (Faktor "
      f"{f6(max(pB) / min(pB), 4)} über die sieben c) und weicht von Lesart A ({f6(vals_W[0])}) ab. Die im Manuskript "
      "nicht definierte chi-Abhängigkeit von U_Sigma ist also eine wertbestimmende Wahl (P_U)",
      monoton and max(pB) / min(pB) > 5 and all(mp.mpf("1e-3") < v < 10 for v in pB))
ZAHLEN["lesart_B"] = {k: {"rho_*": f6(v[0], 10), "G/R_Vak^2": f6(v[1], 10)} for k, v in sorted(les_B.items())}
befund("Native Werte für rho_* und lambda_Sigma", "nein",
       "E/C: H6 liefert kein Paar ohne Wahl. (1) In G/R_Vak^2 zählt nur das Produkt rho_* c^2 = (chi_0 R_Vak)^2 = "
       "S_dS; die Relation ist die Definition dieser einen Zahl. (2) lambda_Sigma hat einen definitionstreuen "
       "Kandidaten (1/(4 R_Vak), P_B), rho_* keinen: Es bräuchte eine zweite native Skala chi_0, und die konforme "
       "Kinematik hat nur R_Vak. (3) Die Wurzel selbst ist die Spur der Einsteingleichung mit der freien Energie der "
       "Naht als Vakuumenergie; lambda_Sigma fällt heraus, der Wert hängt an drei Lesewahlen (Nahtgewicht, kappa, "
       "chi-Abhängigkeit) und am Vorzeichen")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil 3: Überbestimmungstabelle (z136), ergänzt")
z136 = json.loads((HERE / "z136_ergebnis.json").read_text())
z127 = json.loads((HERE / "z127_ergebnis.json").read_text())
z137 = json.loads((HERE / "z137_ergebnis.json").read_text())
alt = []
for row in z136["zahlen"]["tabelle"]:
    if row["G/R_Vak^2"] is not None:
        alt.append((row["Quelle"], row["Lesart"], mp.mpf(row["G/R_Vak^2"]), row["Status"]))
for name, w in sorted(z127["zahlen"]["mm_normierung"]["natives kappa: [Lambda G, G/R_Vak^2]"].items()):
    alt.append(("z127", f"1/(4096 pi kappa), kappa = {name}", mp.mpf(w[1]), "C"))
# exakte Formen der Zeilen, die hier als Vergleich dienen (gegen die Tabellenwerte geprüft)
EXAKT_ALT = {
    "(a2) reines Gewicht w = 1 (Index) auf 64 eps, lambda = 4/R_Vak^2": 1 / (4096 * sp.pi),
    "(a2) reines Gewicht w = 2 c_3 auf 64 eps, lambda = 4/R_Vak^2": sp.Rational(1, 1024),
    "(a2) reines Gewicht w = c_3 (Majorana-Weyl)": sp.Rational(1, 512),
    "l = R_Vak/2 (MM-Länge), asymptotisch": sp.pi / (9 * sp.zeta(3)),
}
for Lname, Lw in (("R_Vak/2", sp.Rational(1, 2)), ("R_Vak", sp.Integer(1)), ("R_Vak/tanh(pi)", 1 / sp.tanh(sp.pi))):
    for gname, g2 in (("1", sp.Integer(1)), ("2", sp.Integer(2)), ("2 pi", 2 * sp.pi), ("4 pi", 4 * sp.pi)):
        EXAKT_ALT[f"L = kappa/(2 g^2) = {Lname}, g^2 = {gname}"] = g2**2 * Lw**2 / (8 * sp.pi)
for sname, sw in (("Quarter (ein Majorana-NS-Spektrum, alpha = 1/4)", sp.log(2) / 4),
                  ("Naht (beide Quarter)", sp.log(2) / 2),
                  ("Träger 5 (g_car = 5 Nähte)", 5 * sp.log(2) / 2),
                  ("16 Quarter ohne GSO (= c = 8 Nähte)", 4 * sp.log(2)),
                  ("E8-Naht (16 Quarter mit GSO)", sp.log(12)),
                  ("drei Familien: 3 E8-Nähte", 3 * sp.log(12)),
                  ("drei Familien: 48 Quarter ohne GSO", 12 * sp.log(2))):
    EXAKT_ALT[f"s_0 = {sname}; A = 4 pi R_Vak^2 (flach/ESU)"] = sp.pi / sw
    EXAKT_ALT[f"s_0 = {sname}; A = pi R_Vak^2 (dS)"] = sp.pi / (4 * sw)
rek = {}
for quelle, lesart, wert, st in alt:
    if lesart in EXAKT_ALT:
        rek[lesart] = abs(N(EXAKT_ALT[lesart]) / wert - 1)
check(f"[E] z136-Tabelle gelesen: {len(alt)} Zeilen mit Zahlenwert (darunter die drei Gewichte aus z127). "
      f"{len(rek)} Zeilen sind hier aus ihrer geschlossenen Form nachgerechnet (z118/z127: 1/(4096 pi w); "
      "Flächenzustand pi/(9 zeta(3)); Jacobson pi/s_0 und pi/(4 s_0); z112: g^4 L^2/(8 pi)); größte relative "
      f"Abweichung {f6(max(rek.values()), 2)}", len(rek) == len(EXAKT_ALT) and max(rek.values()) < mp.mpf("1e-9"))

neu = ([(k, n_, N(v), "L", v) for k, n_, v in fam_M] + [(k, n_, N(v), "L", v) for k, n_, v, _, _ in fam_S]
       + [(k, n_, N(v), "C", v) for k, n_, v in fam_W])
treffer, naechst = [], None
for nk, nname, nv, nst, nex in neu:
    for quelle, lesart, av, ast in alt:
        d = abs(nv / av - 1)
        if d < mp.mpf("1e-9"):
            ex = lesart in EXAKT_ALT and sp.simplify(nex - EXAKT_ALT[lesart]) == 0
            treffer.append((nk, nname, quelle, lesart, bool(ex)))
        elif naechst is None or d < naechst[0]:
            naechst = (d, nk, nname, quelle, lesart)
say(f"   Neue Zeilen: {len(neu)} (Familie M 7, Familie S 56, Wurzelroute 12); Vergleich mit {len(alt)} alten Zeilen "
    f"({len(neu) * len(alt)} Paare)")
for nk, nname, quelle, lesart, ex in treffer:
    say(f"      gleich: {nk} [{nname}]  =  {quelle} [{lesart}]  (symbolisch: {ex})")
# dieselbe Aussage: neue Zeile 'rho_* = native Entropie' gegen alte Zeile 'Jacobson mit nativer Entropie'
taut = [t for t in treffer if t[0].split("/")[0] == "S" and t[0].split("/")[1] in ("r1", "r2", "r3")
        and t[2].startswith("(e) Jacobson")]
zaehlend = [t for t in treffer if t not in taut]
check(f"[E] Exakte Übereinstimmungen neuer mit alten Zeilen: {len(treffer)}, alle symbolisch bestätigt "
      f"({all(t[4] for t in treffer)}). Alle {len(taut)} sind dieselbe Aussage in zwei Schreibweisen: 'rho_* ist eine "
      "native Entropie' gegen 'Jacobson mit nativer Entropie'. Weil rho_* c^2 = S_dS ist (2.2), ist rho_* = s_0 mit "
      "c = 1 Jacobson im dS-Rahmen und mit c = 1/2 Jacobson im flachen Rahmen; andere c geben ganzzahlige Vielfache "
      "derselben Entropie. Beide Seiten tragen dieselbe Prämisse P_S; die Übereinstimmung prüft nichts. Zählende "
      f"Übereinstimmungen (beide Seiten ohne Wahlfreiheit festgelegt): {len(zaehlend)}",
      len(treffer) >= 6 and len(taut) == len(treffer) and all(t[4] for t in treffer) and not zaehlend)
neu_exakt = {k: ex for k, _, _, _, ex in neu}
neu_wert = {k: v for k, _, v, _, _ in neu}
alt_wert = {lesart: av for _, lesart, av, _ in alt}
form, verh = "--", alt_wert[naechst[4]] / neu_wert[naechst[1]]          # Verhältnis alt/neu
if naechst[4] in EXAKT_ALT:
    form = str(sp.simplify(EXAKT_ALT[naechst[4]] / neu_exakt[naechst[1]]))
innerhalb5 = sum(1 for nk, nname, nv, nst, nex in neu for quelle, lesart, av, ast in alt
                 if mp.mpf("1e-9") <= abs(nv / av - 1) < mp.mpf("0.05"))
check(f"[E] Nächstes nicht identisches Paar: {naechst[1]} [{naechst[2]}] gegen {naechst[3]} [{naechst[4]}], Verhältnis "
      f"{form} = {f6(verh, 8)}, nicht 1. Es vergleicht die Wurzelroute mit einer Naht (kappa = 1, mit "
      f"Casimir) mit Jacobson für fünf Nähte; als Nahtzahl gelesen 24/(pi^2 ln^2 2) = "
      f"{f6(24 / (PI**2 * mp.log(2)**2), 6)} statt 5. Von {len(neu) * len(alt)} Paaren liegen {innerhalb5} innerhalb "
      "5 %, wie bei so vielen Zahlen der Ordnung 1 zu erwarten; Zahlennähe ist kein Befund (z86, z136)",
      naechst[0] > mp.mpf("1e-6") and naechst[1] == "W/L3/k1"
      and abs(N(sp.sympify(form)) - verh) < mp.mpf("1e-9") and abs(verh - 1) > mp.mpf("0.01"))

# Wörterbuch der Zeilen
kap_s, L_s, N_s, c1_s, s_s, Nn_s = sp.symbols("kappa_MM L N c_1 s_0 N_Naht", positive=True)
X = sp.Symbol("X", positive=True)                                       # X = rho_* c^2 = S_dS
woerter = {
    "z127: G/R_Vak^2 = 1/(4096 pi kappa)": sp.solve(sp.Eq(sp.pi / (4 * X), 1 / (4096 * sp.pi * kap_s)), X)[0],
    "z137: l^2/G = (4 c_1/pi) L N, l = R_Vak/2": sp.solve(sp.Eq(sp.pi / (4 * X), sp.pi / (16 * c1_s * L_s * N_s)),
                                                          X)[0],
    "z136 Jacobson (dS): G/R_Vak^2 = pi/(4 N s_0)": sp.solve(sp.Eq(sp.pi / (4 * X), sp.pi / (4 * Nn_s * s_s)), X)[0],
    "Wurzelroute: G/R_Vak^2 = pi^3 kappa L^/24": sp.solve(sp.Eq(sp.pi / (4 * X), sp.pi**3 * ks * Ls / 24), X)[0],
}
check("[E] Wörterbuch: Jede Zeile der Tabelle ist dieselbe freie Zahl X = rho_* c^2 = S_dS = pi/(4 N_G) in einer "
      f"anderen Parametrisierung: { {k: str(v) for k, v in woerter.items()} }. Die rho_*-Relation fügt eine fünfte "
      "Schreibweise hinzu, keine Gleichung zwischen zweien",
      woerter["z127: G/R_Vak^2 = 1/(4096 pi kappa)"] == 1024 * sp.pi**2 * kap_s
      and woerter["z137: l^2/G = (4 c_1/pi) L N, l = R_Vak/2"] == 4 * c1_s * L_s * N_s
      and woerter["z136 Jacobson (dS): G/R_Vak^2 = pi/(4 N s_0)"] == Nn_s * s_s
      and sp.simplify(woerter["Wurzelroute: G/R_Vak^2 = pi^3 kappa L^/24"] - 6 / (sp.pi**2 * ks * Ls)) == 0)
ZAHLEN["tabelle_neu"] = {"Zeilen neu": len(neu), "Zeilen alt": len(alt),
                         "exakte Übereinstimmungen": [list(map(str, t)) for t in treffer],
                         "zählend": len(zaehlend),
                         "nächstes Paar": [naechst[1], naechst[2], naechst[3], naechst[4], f6(verh, 10), form],
                         "Paare innerhalb 5 %": innerhalb5,
                         "Wörterbuch X = rho_* c^2": {k: str(v) for k, v in woerter.items()}}
befund("Überbestimmung", "nein",
       f"E: Die rho_*-Zeilen stimmen mit keiner vorhandenen Zeile in einer Weise überein, die zählt. Die {len(treffer)} "
       "exakten Gleichheiten sind Jacobson-Zeilen selbst (rho_* c^2 = S_dS, Prämisse P_S auf beiden Seiten). Die "
       "Wurzelroute gibt neue, von allen anderen verschiedene Werte, wieder von der Ordnung 0,001 bis 1")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil 4: Ledger-Ziel COSMO.LAM.01")
A = sp.Symbol("A", positive=True)                                      # A = alpha^-1
rhoL = 3 / (4 * sp.pi**2) * sp.exp(-2 * A)                             # rho_Lambda/Mbar^4
GLam = sp.simplify(rhoL / (8 * sp.pi))                                 # G Lambda = rho_L/(8 pi Mbar^4)
pG_ziel = sp.simplify(GLam / 12)                                       # Lambda = 12/R_Vak^2 (dS, l = R_Vak/2)
X_ziel = sp.simplify(sp.pi / (4 * pG_ziel))
l2G = sp.simplify(sp.Rational(1, 4) / pG_ziel)                         # l^2/G = R_Vak^2/(4 G)
check(f"[E] Ziel: rho_Lambda/Mbar^4 = (3/4 pi^2) e^(-2/alpha) gibt G Lambda = {GLam}, mit Lambda = 12/R_Vak^2 "
      f"G/R_Vak^2 = {pG_ziel} und l^2/G = {l2G} (z96). Nötig wäre rho_* c^2 = S_dS = {X_ziel}. Dieselbe Zahl in den "
      "anderen Zeilen: kappa_MM = X/(1024 pi^2) = (pi^2/32) e^(2/alpha) (z127), N L c_1 = X/4 = 8 pi^4 e^(2/alpha) "
      "(z137), N s_0 = X (z136)",
      sp.simplify(pG_ziel - sp.exp(-2 * A) / (128 * sp.pi**3)) == 0
      and sp.simplify(l2G - 32 * sp.pi**3 * sp.exp(2 * A)) == 0
      and sp.simplify(X_ziel - 32 * sp.pi**4 * sp.exp(2 * A)) == 0
      and sp.simplify(X_ziel / (1024 * sp.pi**2) - sp.pi**2 * sp.exp(2 * A) / 32) == 0
      and sp.simplify(X_ziel / 4 - 8 * sp.pi**4 * sp.exp(2 * A)) == 0)
Xz = 32 * PI**4 * mp.exp(2 * ainv)
pz = mp.exp(-2 * ainv) / (128 * PI**3)
rho_noetig = {key: Xz / N(c)**2 for key, _, c in C_KAND}
chi0R = mp.sqrt(Xz)
z136_ziel = [w for qn, ls_, w, st in alt if qn.startswith("(b')")][0]
z137_lg = z137["zahlen"]["vergleichsziel"]["log10 S_dS (Ziel)"]
check(f"[N] Zahlen (alpha^-1 = 137,0359992168, EM.FP.01): G/R_Vak^2 = {f6(pz)} (z136: {f6(z136_ziel)}), "
      f"S_dS = rho_* c^2 = {f6(Xz)} = 10^{lg(Xz)} (z137: 10^{z137_lg:.2f}). Nötiges rho_* je lambda_Sigma: "
      f"{ {k: '10^' + lg(v) for k, v in sorted(rho_noetig.items())} }; definitionstreu (c = 1/4) rho_* = "
      f"512 pi^4 e^(2/alpha) = {f6(rho_noetig['c1'])}. Als Skala: chi_0 R_Vak = 4 sqrt(2) pi^2 e^(1/alpha) = "
      f"{f6(chi0R)}, die Hierarchie e^(1/alpha) in der Länge (z112)",
      abs(pz / z136_ziel - 1) < mp.mpf("1e-5") and abs(float(mp.log10(Xz)) - z137_lg) < 1e-5
      and abs(rho_noetig["c1"] / (512 * PI**4 * mp.exp(2 * ainv)) - 1) < mp.mpf("1e-40")
      and abs(chi0R / (4 * mp.sqrt(2) * PI**2 * mp.exp(ainv)) - 1) < mp.mpf("1e-40"))
# Wurzelroute: nötiges L^
L_ziel = {kk: sp.simplify(6 / (sp.pi**2 * kv * X_ziel)) for kk, _, kv in K_KAND}
check(f"[E] In der Wurzelroute ist das Ziel eine Forderung an die Naht: L^ = 6/(pi^2 kappa S_dS) = "
      f"{ {k: str(v) for k, v in sorted(L_ziel.items())} }. Mit der Normierung des Manuskripts (k2, M_Pl^4) ist das "
      "wörtlich delta_top e^(-2/alpha), die Vermutung 040:324-340; mit reduzierter Normierung (k3) ist es "
      "(3/4 pi^2) e^(-2/alpha), also COSMO.LAM.01 selbst. Die Wurzelgleichung ist das Ledger-Ziel in anderen "
      "Variablen und sagt es nicht vorher",
      sp.simplify(L_ziel["k2"] - 3 * sp.exp(-2 * A) / (256 * sp.pi**4)) == 0
      and sp.simplify(L_ziel["k3"] - rhoL) == 0
      and sp.simplify(L_ziel["k1"] - 3 * sp.exp(-2 * A) / (16 * sp.pi**6)) == 0)
L_verm = -mp.log1p(-delta_top * mp.exp(-2 * ainv))
gedruckt = mp.mpf("1.1280419085743317e-123")
luecke = mp.log10(phi_seam / (delta_top * mp.exp(-2 * ainv)))
daempf = 2 * ainv - mp.log(delta_top)
check(f"[E] Das Manuskript führt den Spektralradius von U_Sigma zweifach: Satzschicht phi_seam(alpha) = "
      f"{f6(phi_seam, 8)} (029:29-45), bedingte Schicht delta_top e^(-2/alpha) = {f6(delta_top * mp.exp(-2 * ainv), 8)} "
      f"(040:324-340; -ln(1 - .) = {f6(L_verm, 12)}, gedruckt 1.1280419086e-123). Dazwischen liegen "
      f"{f6(luecke, 5)} Größenordnungen. Nativ ist die Dämpfung des tiefsten Gewichts -ln q^(1/2) = pi; das Ziel "
      f"braucht 2/alpha - ln delta_top = {f6(daempf, 8)}, das {f6(daempf / PI, 6)}-Fache",
      abs(L_verm / gedruckt - 1) < mp.mpf("1e-9") and 121 < luecke < 122.5)
# native Bedeutung der Zielzahl
teiler = {"pi (tiefstes Gewicht e^-pi)": PI, "2 pi (-ln q)": 2 * PI, "4 pi (beta)": 4 * PI, "ln 2/2 (s_0 Naht)": s0,
          "ln 2": mp.log(2), "ln 12 (E8-Naht)": mp.log(12)}
vielf = {k: 2 * ainv / v for k, v in teiler.items()}
vielf_ln = {k: mp.log(Xz) / v for k, v in teiler.items()}
abst_v = {k: abs(v - mp.nint(v)) for k, v in vielf.items()}
abst_ln = {k: abs(v - mp.nint(v)) for k, v in vielf_ln.items()}
naehte = {"Quarter": Xz / (mp.log(2) / 4), "Naht": Xz / s0, "E8-Naht": Xz / mp.log(12)}
k_min = min(sorted(abst_ln), key=lambda k: abst_ln[k])
d_min = min(min(abst_v.values()), min(abst_ln.values()))
zufall = 1 - (1 - 2 * d_min)**(len(abst_v) + len(abst_ln))
pot2 = Xz / mp.mpf(2)**int(mp.nint(vielf_ln["ln 2"]))
check("[E negativ, gescopt] Native Bedeutung von S_dS = 32 pi^4 e^(2/alpha): Die nativen Boltzmann-Exponenten sind "
      "Vielfache von pi, die nativen Entropien Vielfache von ln 2 und ln 12. 2/alpha in diesen Einheiten: "
      f"{ {k: f6(v, 7) for k, v in sorted(vielf.items())} }; ln S_dS in diesen Einheiten: "
      f"{ {k: f6(v, 7) for k, v in sorted(vielf_ln.items())} }. Kein Wert ist ganz. Der kleinste Abstand ist "
      f"{f6(d_min, 3)} bei ln S_dS/{k_min.split(' (')[0]}: S_dS = {f6(pot2, 6)} x 2^{int(mp.nint(vielf_ln['ln 2']))}. "
      f"Bei {len(abst_v) + len(abst_ln)} geprüften Verhältnissen ist ein so kleiner Abstand mit Wahrscheinlichkeit "
      f"{f6(zufall, 2)} zufällig; die Abweichung von {f6(abs(pot2 - 1) * 100, 2)} % liegt innerhalb der Streuung des "
      "Bestands für rho_Lambda (v60 gegen v206: 3,8 %, aus R1 gelesen), und eine Entropie als Zweierpotenz hat keine "
      "native Lesart "
      "(eine Zustandszahl wäre e^S). Das ist Zahlennähe, kein Befund. Als Zählung wären es "
      f"{ {k: '10^' + lg(v) for k, v in sorted(naehte.items())} } Einheiten je Horizont (z136). Geprüft sind nur "
      "diese sechs Einheiten; e^(2/alpha) kommt aus dem gelösten alpha-Block und hat in den KMS-Daten der Naht kein "
      "Gegenstück", d_min > mp.mpf("1e-6") and abs(pot2 - 1) > mp.mpf("1e-4") and zufall > mp.mpf("0.05"))
check("[E] Vorzeichen des Ziels: rho_Lambda > 0 ist dS. Nativ ist AdS (z96, z127, z137). In der Wurzelroute hat das "
      "gedruckte W bei dS keine Wurzel, das vorzeichenkonsistente keine bei AdS (Teil 2.4). Ein natives rho_* und das "
      "Ziel verlangen also verschiedene Fassungen desselben Funktionals",
      w_druck.subs({sg: 1, cs: 1, ks: 1, Ls: 1}) < 0 and w_kons.subs({sg: -1, cs: 1, ks: 1, Ls: 1}) < 0)
abstand_ziel = {"Familie M": (mp.log10(min(vals_M) / pz), mp.log10(max(vals_M) / pz)),
                "Familie S": (mp.log10(min(vals_S) / pz), mp.log10(max(vals_S) / pz)),
                "Wurzelroute": (mp.log10(min(vals_W) / pz), mp.log10(max(vals_W) / pz))}
check("[N] Abstand zum Ziel in Größenordnungen: "
      f"{ {k: [f6(v[0], 5), f6(v[1], 5)] for k, v in sorted(abstand_ziel.items())} }. Jede native Fassung liegt "
      "mindestens 118 Größenordnungen darüber",
      min(v[0] for v in abstand_ziel.values()) > 118)
ZAHLEN["ziel"] = {"G/R_Vak^2": f6(pz, 12), "S_dS = rho_* c^2": f6(Xz, 12), "log10 S_dS": lg(Xz),
                  "rho_* je c": {k: f6(v, 12) for k, v in sorted(rho_noetig.items())},
                  "chi_0 R_Vak": f6(chi0R, 12), "L^ nötig": {k: str(v) for k, v in sorted(L_ziel.items())},
                  "2/alpha in nativen Einheiten": {k: f6(v, 12) for k, v in sorted(vielf.items())},
                  "ln S_dS in nativen Einheiten": {k: f6(v, 12) for k, v in sorted(vielf_ln.items())},
                  "Lücke Satzschicht/Vermutung (log10)": f6(luecke, 8),
                  "Abstand zum Ziel (log10)": {k: [f6(v[0], 6), f6(v[1], 6)] for k, v in sorted(abstand_ziel.items())}}
befund("Ziel", "nein",
       "E: Das Ziel verlangt rho_* c^2 = 32 pi^4 e^(2/alpha) = 10^122,52, die dS-Entropie des Ziels. In der "
       "Wurzelroute ist das die Vermutung des Manuskripts (delta_top e^(-2/alpha)) bzw. COSMO.LAM.01 selbst. Eine "
       "native Bedeutung hat die Zahl nicht: kein ganzes Vielfaches nativer Exponenten oder Entropien, und das "
       "Vorzeichen (dS) ist dem nativen entgegengesetzt")

# ------------------------------------------------------------------------------------------------------------------
say("\nTeil 5: Bilanz und Urteil")
bilanz = [
    ("E", "Bestand",
     "rho_* ist in tfpt-42 die stationäre Wurzel eines randnormierten Potentials (015:243-288), lambda_Sigma der "
     "erste positive Eigenwert des Randoperators (015:158). Ein Zahlenwert steht nirgends. G lambda_Sigma^2 = "
     "pi/(4 rho_*) ist G chi_0^2 = pi/4 in anderen Variablen (013:192-204)."),
    ("E", "Annahmen",
     "48 Komponenten mit E = -R/4 und Residuenprofil (N = 96 oder Boltzmann-Profil: pi/8), Existenz einer "
     "stationären Skala chi_0 (drei ungeprüfte Annahmen), Mbar^-2 = 8 pi G. Weder G = 1 noch eine Flächenregel; "
     "nicht zirkulär im Sinn von GRAV.NONCIRCULAR.01, aber ein 4D-Wärmekern ist Import."),
    ("E", "Identität",
     "rho_* (lambda_Sigma R_Vak)^2 = (chi_0 R_Vak)^2 = pi R_Vak^2/(4 G) = S_dS. Nur das Produkt zählt; P_G ist "
     "gleichwertig zu dieser einen Zahl."),
    ("C/L", "Native Kandidaten",
     f"lambda_Sigma = 1/(4 R_Vak) ist definitionstreu (P_B). Für rho_* gibt es kein natives chi_0. 7 Massenverhältnisse "
     f"geben G/R_Vak^2 = {f6(min(vals_M), 3)} bis {f6(max(vals_M), 3)}, 56 Paarungen mit reinen Zahlen "
     f"{f6(min(vals_S), 3)} bis {f6(max(vals_S), 3)}."),
    ("C", "Wurzelroute (P_rho*)",
     f"Die Stationarität ist die Spur der Einsteingleichung: G/R_Vak^2 = pi^3 kappa L^/24, lambda_Sigma fällt heraus. "
     f"Mit nativen Nahtgewichten {f6(min(vals_W), 3)} bis {f6(max(vals_W), 3)}; fermionisch mit kappa = 1: "
     f"pi^3 (ln 2/2 - pi/12)/24 = {f6(vals_W[0])}."),
    ("E", "Unstimmigkeiten im Funktional",
     "Vorzeichen (gedrucktes W hat nur bei AdS eine Wurzel, die Wirkung nur bei dS), dreifache Normierung des "
     "Nahtterms, nicht definierte chi-Abhängigkeit von U_Sigma (Lesart B ändert den Wert und macht ihn von "
     "lambda_Sigma abhängig), lambda_Sigma = chi_seed^-1 gegen chi_seed als Masse."),
    ("E", "Überbestimmung",
     f"{len(treffer)} exakte Gleichheiten, alle identisch mit Jacobson-Zeilen (P_S auf beiden Seiten); keine zählt. "
     "Wurzelroute und Jacobson sind für keine ganze Nahtzahl zugleich erfüllt."),
    ("E", "Ziel",
     "rho_* c^2 = 32 pi^4 e^(2/alpha) = 10^122,52. In der Wurzelroute gleich der IR-Vermutung des Manuskripts bzw. "
     "COSMO.LAM.01. Keine native Bedeutung; Vorzeichen entgegengesetzt."),
    ("O", "nicht geprüft",
     "E_Sigma und A_Sigma (im Manuskript nur benannt), eine renormierte Fassung Gamma^ren, die Universalität der 48 "
     "(z27), weitere chi-Abhängigkeiten von U_Sigma außer Lesart A und B."),
]
for st, was, text in bilanz:
    say(f"   [{st}] {was}: {text}")
urteil = ("rho_* legt P_G nicht fest; es gibt nur eine weitere Zeile, und sie ist eine Umschreibung. Nichts ist "
          "befördert. (1) Bestand: tfpt-42 definiert rho_* als stationäre Wurzel eines randnormierten Potentials und "
          "lambda_Sigma als ersten positiven Eigenwert des Randoperators. Ausgewertet ist das nirgends. Die Relation "
          "G lambda_Sigma^2 = pi/(4 rho_*) ist G chi_0^2 = pi/4, die induzierte Einsteinkopplung von 48 Komponenten "
          "mit Residuenprofil; das Manuskript nennt sie selbst eine Umschreibung. Sie setzt weder G = 1 noch eine "
          "Flächenregel voraus, tauscht G aber nur gegen die Skala chi_0. (2) In G/R_Vak^2 = pi/(4 rho_* c^2) zählt "
          "nur das Produkt rho_* c^2 = (chi_0 R_Vak)^2, und das ist exakt die dS-Entropie S_dS = pi R_Vak^2/(4 G). "
          "Die Frage nach rho_* ist also dieselbe wie die nach P_eta aus z136. (3) Nativ: lambda_Sigma hat einen "
          "definitionstreuen Kandidaten, 1/(4 R_Vak) (Prämisse P_B). Für rho_* fehlt die zweite Skala; jede native "
          "Masse ist eine Zahl mal 1/R_Vak. Alle 63 Kandidatenpaarungen geben ein planckgroßes Vakuum. (4) Die "
          "Wurzel selbst ist die Spur der Vakuum-Einsteingleichung mit der freien Energie der Naht als "
          "Vakuumenergie: G/R_Vak^2 = pi^3 kappa L^/24, lambda_Sigma fällt heraus. Mit nativen Gewichten sind das "
          f"{f6(min(vals_W), 3)} bis {f6(max(vals_W), 3)}, je nach Nahtgewicht und Normierung. Das ist bedingt fest "
          "nur unter dem Bündel P_rho* (P_48, P_U, P_kappa, P_vac), vier Prämissen statt einer freien Zahl, dazu ein "
          "Vorzeichenkonflikt im gedruckten Funktional. (5) Überbestimmung: Keine Übereinstimmung zählt. Die "
          f"{len(treffer)} exakten Gleichheiten sind Jacobson-Zeilen in anderer Schreibweise. Wurzelroute und Jacobson "
          "schließen sich für jede ganze Nahtzahl aus. (6) Das Ledger-Ziel verlangt rho_* c^2 = 32 pi^4 e^(2/alpha); in der "
          "Wurzelroute ist das die IR-Vermutung des Manuskripts bzw. COSMO.LAM.01 selbst, ohne native Bedeutung und "
          "mit dem falschen Vorzeichen. K13 bleibt L, P_G bleibt eine freie Zahl, die neun Prämissen sind "
          "unverändert.")
say("\nUrteil\n   " + urteil)
print(f"\n   (Laufzeit {time.time() - T0:.1f} s, nicht Teil der Ausgabe)", flush=True)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z147_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "bilanz": [list(b) for b in bilanz], "urteil": urteil},
    indent=1, ensure_ascii=False, default=str) + "\n")
(HERE / "ausgabe_z147.txt").write_text("\n".join(OUT) + "\n")
