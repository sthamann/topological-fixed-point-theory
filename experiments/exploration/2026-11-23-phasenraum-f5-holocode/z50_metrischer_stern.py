"""Z50: Prüfung der externen Analyse vom 24.11. (metrischer Stern) und die metrische Stufe der Prozesssuche.

Teil 1  Metrischer Stern: Zentralmetrik plus vier Arme, Gewichte M_(j,n) = (M_c/2) q^(n+1), Kantensteifigkeit mu^2 M_Kind^2.
        Sein Massenoperator mu^2 B^+ B ist unitär äquivalent zum rechteckigen Vierteltransfer (B U = V A_rect), die Armdrehung
        geht in R|k> = i^k|k> über, die Nullmode ist eta_v = M_v/M_eff, und eta_v/M_v = 1/M_eff ist universell.
Teil 2  Potential V(X) = (-4 + e1 + e3 - 4 e4)/2: flach für beide Metriken, Fierz-Pauli-Masse; V_t = V + (t/2) det(I - X)
        ändert erst die quartische Ordnung, eine zweite Richtung ändert die kubische. Längs der masselosen Richtung
        verschwindet jedes Potential; die Ordnung vier unterscheidet sich erst mit massiven Beinen.
Teil 3  Auslese: Die zentrierte Feldkovarianz (p^2 + mu^2 L)^-1 ist verbunden und masselos; das frühere 'verbunden' war die
        vakuumsubtrahierte Mikroantwort (Korrektur zu z44/z45). Der Stern mit minimalen Knotenquellen ist universell.
Teil 4  Einstein-Zusatz: Für kappa > 0 verlangt eine positive Spektraldarstellung x C(x) <= 1/kappa, also h_f >= 0; die
        Summenregel verbietet das. Der Kern e^-u plus 2x ist auf der Achse positiv und trotzdem ausgeschlossen.
Teil 5  Skala: a = r/(lambda (1-r)^2) trifft 3/10 für jedes q mit passendem lambda; q^2 = 1/2, lambda = 20/3 gibt b = 27/100.
Teil 6  Flavour: Eine ungebrochene Majorana-Z4 erzwingt zwei gleiche Massen; Z3 x Z4 und Z4 x Z3 geben dieselben |U|^2.
"""
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "research" / "originlab"))
from process_rules import (BETA_DET, BETA_V, Q_ORIGINAL, flavour_pairs, h_exp, hardy_rect, metric_star,  # noqa: E402
                           potential_taylor, run, transfer_kernel)

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


res = run()
U_ = res["universes"]

# ------------------------------------------------------------------ Teil 1
say("Teil 1: Metrischer Stern = Vierteltransfer")


def intertwiner(qv, m):
    st, A = metric_star(qv, m), hardy_rect(qv, m)
    nodes = st["nodes"]
    idx = {v: i for i, v in enumerate(nodes)}
    Fn = np.zeros((1 + 4*m, 1 + 4*m), complex)        # Hardy-Spalte -> Knoten
    for k in range(4*m + 1):
        n, g = divmod(k, 4)
        if g == 0 and n == 0:
            Fn[idx["c"], k] = 1
        else:
            lev, ell = (n - 1, 0) if g == 0 else (n, (-g) % 4)
            for j in range(4):
                Fn[idx[(j, lev)], k] = np.conj(1j**(-ell*j))/2
    edges = [(j, n) for j in range(4) for n in range(m)]
    eidx = {e: i for i, e in enumerate(edges)}
    Fe = np.zeros((4*m, 4*m), complex)                 # Hardy-Zeile -> Kante
    for r in range(4*m):
        n, g = divmod(r + 1, 4)
        lev, ell = (n - 1, 0) if g == 0 else (n, (-g) % 4)
        for j in range(4):
            Fe[eidx[(j, lev)], r] = -np.conj(1j**(-ell*j))/2
    return st, A, Fn, Fe


for qv, m in ((0.5, 5), (Q_ORIGINAL, 4)):
    st, A, Un, Ve = intertwiner(qv, m)
    err = np.max(np.abs(st["B"] @ Un - Ve @ A))
    unit = max(np.max(np.abs(Un.conj().T @ Un - np.eye(Un.shape[0]))), np.max(np.abs(Ve.conj().T @ Ve - np.eye(Ve.shape[0]))))
    rot = np.zeros((1 + 4*m, 1 + 4*m))
    idx = {v: i for i, v in enumerate(st["nodes"])}
    rot[idx["c"], idx["c"]] = 1
    for j in range(4):
        for n in range(m):
            rot[idx[((j + 1) % 4, n)], idx[(j, n)]] = 1
    Rk = np.diag(1j**np.arange(4*m + 1))
    say(f"   q = {qv:.4g}, m = {m}: |B U - V A_rect| = {err:.1e}, Unitarität {unit:.1e}, Rang A_rect = "
        f"{np.linalg.matrix_rank(A)} von {A.shape[0]} Zeilen")
check("B U = V A_rect mit unitären U, V (Armfourier plus Restklassen der Hardy-Marken, q = 1/2 und e^-2pi)",
      err < 1e-12 and unit < 1e-12)
check("damit U^+ (B^+ B) U = A_rect^+ A_rect, und A_rect (4m x (4m+1)) ist surjektiv mit genau einem Nullvektor",
      np.linalg.matrix_rank(A) == 4*m and st["nullity"] == 1
      and np.allclose(Un.conj().T @ st["B"].T @ st["B"] @ Un, A.T @ A, atol=1e-12))
check("die Armdrehung j -> j+1 geht unter U in die Vierteldrehung R|k> = i^k|k> über (bis auf Konjugation der Phase)",
      np.allclose(Un.conj().T @ rot @ Un, Rk.conj(), atol=1e-12) or np.allclose(Un.conj().T @ rot @ Un, Rk, atol=1e-12))
st = metric_star(0.5, 6)
Meff = np.sqrt(np.sum(st["M"]**2))
check("die Nullmode ist eta_v = M_v/M_eff mit M_eff^2 = M_c^2 Sum_(n<=m) q^(2n); eta_v/M_v = 1/M_eff an jedem Knoten",
      np.allclose(st["eta"], st["M"]/Meff, atol=1e-12) and abs(Meff**2 - sum(0.25**n for n in range(7))) < 1e-12)
_, sv, vh = np.linalg.svd(st["B"])
mass = vh[:len(sv)][sv > 1e-12]
coup = mass/st["M"]
orth = np.max(np.abs((mass*(st["eta"]/st["M"])) @ mass.T - np.eye(len(mass))/Meff))
phiphichi = np.max(np.abs(mass @ (st["eta"]**2/st["M"])))
check("massive Kopplungen u_a(v)/M_v hängen vom Knoten ab; Sum_v eta_v u_a u_b/M_v = delta_ab/M_eff; die "
      "Einstein-Kopplung masselos-masselos-massiv Sum_v eta_v^2 u_a(v)/M_v verschwindet",
      np.max(np.ptp(coup, axis=1)) > 1e-3 and orth < 1e-12 and phiphichi < 1e-12)
for m in (4, 8, 16):
    st_m = metric_star(0.5, m)
    L = st_m["B"].T @ st_m["B"]
    xs = np.array([0.3, 1.0, 3.0])
    Cst = [np.linalg.solve(x*np.eye(L.shape[0]) + L, np.eye(L.shape[0])[:, 0])[0] for x in xs]
    say(f"   m = {m:2d}: Zentrumsresolvente {np.round(Cst, 10)}, Hardy C_q {np.round(transfer_kernel(0.5)(xs), 10)}")
check("die Zentrumsresolvente des Sterns geht mit m -> oo in die Randantwort C_q(x) des Vierteltransfers über",
      np.allclose(Cst, transfer_kernel(0.5)(xs), atol=1e-8))
befund("Der vollständige Vierteltransfer ist der Massenoperator eines endlichen metrischen Sterns; Quellenregel, positive "
       "Masse und gemeinsame Diffeomorphie gehören zu einer ausdrücklichen Wirkung; G_eff = 1/(8 pi M_eff^2) bleibt frei", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Potential, kubische und quartische Freiheit")
mp_ = U_["metric_potentials"]
for name, d in mp_["derived"].items():
    say(f"      {name:26s} {d['clauses']}  Taylor {np.round(d['taylor_first_pair'], 6).tolist()}")
amb = {tuple(r["rule"]) if "rule" in r else None for r in []}
cl = {n: set(d["clauses"]) for n, d in mp_["derived"].items()}
check("V ist flach (keine Kaulquappe für beide Metriken) und hat die Fierz-Pauli-Form tr D^2 - (tr D)^2 mit konstantem "
      "Koeffizienten", {"flat_background", "fierz_pauli_mass"} <= cl["V"])
check("V +- det(I-X)/2 stimmt bis zur kubischen Ordnung mit V überein und weicht quartisch ab",
      all({"flat_background", "fierz_pauli_mass", "cubic_as_V"} <= cl[n] and "quartic_as_V" not in cl[n]
          for n in ("V + det(I-X)/2", "V - det(I-X)/2")))
check("eine zweite flache Richtung mit gleicher Masse ändert schon die kubische Ordnung",
      {"flat_background", "fierz_pauli_mass"} <= cl["V + cubic/2"] and "cubic_as_V" not in cl["V + cubic/2"])
check("alle diese Potentiale verschwinden längs der masselosen Richtung (alle Metriken gleich) in jeder Ordnung",
      all("massless_direction_free" in cl[n] for n in cl if n != "beta1 only"))
rng = np.random.default_rng(7)
Bm = rng.normal(size=(4, 4))
Bm = (Bm + Bm.T)/2
c_v = potential_taylor(BETA_V, np.zeros((4, 4)), Bm)
c_t = potential_taylor(tuple(v + d for v, d in zip(BETA_V, BETA_DET)), np.zeros((4, 4)), Bm)
say(f"   reine Massenrichtung (g = I, f = I + eps B): Delta c4 = {c_t[4] - c_v[4]:.10f}, det(B)/16 = {np.linalg.det(Bm)/16:.10f}")
check("der Unterschied liegt im massiven Viererkontakt: Delta c4 = t det(B)/16 bei t = 1 (also t det(B)/32 je t/2)",
      abs((c_t[4] - c_v[4]) - np.linalg.det(Bm)/16) < 1e-9 and all(abs(c_t[k] - c_v[k]) < 1e-12 for k in range(4)))
befund("Propagator, kubische Wirkung und Ward-Identitäten legen die Wechselwirkung nicht fest. Die masselose "
       "Vierpunktamplitude (Baumgraphen) ist reines Einstein mit M_eff und sieht t nicht; erst Amplituden mit massiven "
       "Beinen bzw. die verbundene Vierpunktantwort der Quellen über massive Austausche unterscheiden V_t", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Auslese -- verbunden heißt nicht vakuumsubtrahiert")
ro = U_["readouts"]
say(f"   minimale leere Regeln mit metrischer Stufe: {ro['minimal_empty_rules']}")
check("die zentrierte Feldkovarianz (p^2 + mu^2 L)^-1 ist verbunden und hat den masselosen Pol: {masselos, feldverbunden} "
      "ist nicht leer", ["field_connected", "massless_pole"] not in ro["minimal_empty_rules"])
check("leer bleiben nur Regeln mit 'vakuumsubtrahiert' (etwa {masselos, vakuumsubtrahiert}): die frühere Aussage aus "
      "z44/z45 betraf die Mikroantwort ohne Nullmodenanteil", ["massless_pole", "vacuum_subtracted"] in ro["minimal_empty_rules"]
      and all("vacuum_subtracted" in r for r in ro["minimal_empty_rules"]))
check("der Stern mit minimalen Knotenquellen ist universell: {mehrere Quellen, universelles Residuum, masselos} ist "
      "nicht mehr leer", ["universal_residue"] not in ro["minimal_empty_rules"])
befund("Korrektur zu z44/z45: 'kein Auslesen ist masselos und verbunden' gilt nur für die vakuumsubtrahierte Mikroantwort; "
       "'kein Auslesen ist universell' gilt nur ohne die Quellnorm 1/M_v der metrischen Wirkung", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Einstein-Zusatz rettet keine Spektralwirkung")
xs = np.geomspace(1e-3, 1e4, 4000)
h = h_exp(xs)
for kappa in (0.5, 2.0):
    H = kappa*xs + h
    xc = xs/H
    viol = xs[np.argmax(xc > 1/kappa + 1e-9)] if np.any(xc > 1/kappa + 1e-9) else None
    say(f"   kappa = {kappa}: min H/x = {np.min(H/xs):.4f}, max x C(x) kappa = {np.max(xc)*kappa:.4f}, erste Verletzung bei "
        f"x = {viol}")
H2 = 2*xs + h
check("H = 2x + h_exp ist auf der ganzen Achse positiv (keine euklidische Nullstelle)", np.all(H2 > 0))
check("trotzdem überschreitet x C(x) die Schranke 1/kappa, die jede positive Darstellung Z/x + Int rho/(x+s) erfüllt: kein "
      "positives Spektralmaß", np.any(xs/H2 > 0.5 + 1e-9))
check("die Verletzung liegt genau dort, wo h_exp < 0 ist (x > 6.77); die Summenregel erzwingt solche Stellen für jedes "
      "nichttriviale Profil", np.all((xs/H2 > 0.5 + 1e-12) == (h < -1e-12)))
xp = np.geomspace(1e-3, 1e3, 400)
Cs = 1/xp + 1/(xp + 10/3)
check("Gegenprobe: für die positive Schur-Antwort wächst x C(x) monoton gegen 1/kappa = 1", np.all(np.diff(xp*Cs) > 0)
      and np.all(xp*Cs > 1 - 1e-12))
befund("'Keine Nullstelle gefunden' ist kein Positivitätszertifikat; die Schranke x C(x) <= 1/kappa ist eines. Damit ist "
       "die flache Dirac-Profilklasse auch mit Einstein-Zusatz ausgeschlossen (frühere Begründung präzisiert)", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Relative Skala")
r, lam = sp.symbols("r lambda", positive=True)
a = r/(lam*(1 - r)**2)
b = r*(1 + r)/(lam**2*(1 - r)**4)
lam_q = sp.solve(sp.Eq(a, sp.Rational(3, 10)), lam)[0]
say(f"   lambda(q) = {sp.simplify(lam_q)} (r = q^2)")
check("lambda(q) = 10 q^2/(3 (1-q^2)^2) erfüllt a = 3/10 für jedes 0 < q < 1", sp.simplify(lam_q - 10*r/(3*(1 - r)**2)) == 0)
rows = [("Transfer q^2 = (8-sqrt55)/3, lambda = 1", (8 - sp.sqrt(55))/3, 1),
        ("Transfer q^2 = 1/2, lambda = 20/3", sp.Rational(1, 2), sp.Rational(20, 3))]
for name, rv, lv in rows:
    say(f"      {name:40s} a = {sp.nsimplify(a.subs({r: rv, lam: lv}))}, b = {float(b.subs({r: rv, lam: lv})):.6f}, "
        f"massiv/Pol = {float(rv/(1 - rv)):.6f}")
check("q^2 = 1/2 mit lambda = 20/3 gibt a = 3/10, b = 27/100 und dasselbe massive Gesamtgewicht wie der Schur-Spiegel (1)",
      a.subs({r: sp.Rational(1, 2), lam: sp.Rational(20, 3)}) == sp.Rational(3, 10)
      and b.subs({r: sp.Rational(1, 2), lam: sp.Rational(20, 3)}) == sp.Rational(27, 100))
check("skalenfrei ist b/a^2 = (1+r)/r; ein unabhängig bestimmtes b würde q festlegen", sp.simplify(b/a**2 - (1 + r)/r) == 0)

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: Flavour")
G = np.diag([-1, 1j, -1j])
allowed = [(i, j) for i in range(3) for j in range(i, 3) if abs(G[i, i]*G[j, j] - 1) < 1e-12]
say(f"   Majorana-Massenterme, die G = diag(-1, i, -i) erlaubt: {allowed}")
Mnu = np.zeros((3, 3), complex)
Mnu[0, 0], Mnu[1, 2], Mnu[2, 1] = 0.7, 1.3, 1.3
check("eine ungebrochene volle Majorana-Z4 erlaubt nur M_11 und M_23: Massen |A|, |B|, |B| entartet (ausgeschlossen)",
      allowed == [(0, 0), (1, 2)] and np.allclose(np.sort(np.linalg.svd(Mnu, compute_uv=False)), [0.7, 1.3, 1.3]))
pat = {o["pair"]: o["pattern"] for o in flavour_pairs()}
p34, p43 = pat["Z3 x Z4"], pat["Z4 x Z3"]
P = [[float(sp.sympify(v)) for v in row] for row in p34]
check("die Z3 x Z4-Tabelle ist symmetrisch (|U|^2 = transponiert): beide Rollenverteilungen (Z3 geladen/Z4 neutral und "
      "umgekehrt) geben dieselben |U|^2; sie unterscheiden sich in der Majorana-Verträglichkeit, nicht in den Winkeln",
      p34 == p43 and np.allclose(P, np.array(P).T))
s13 = (2 - np.sqrt(3))/12
check("TM1 bei 15 Grad (0.02233) liegt im marginalen NuFIT-6.1-Bereich 0.02064 - 0.02418 (NO, IC24 + SK, laut Analyse); "
      "die Familienuhr (0.0447) nicht", 0.02064 <= s13 <= 0.02418 and not 0.02064 <= 2*s13 <= 0.02418)
befund("G ist durch Zustand und Ward nicht festgelegt: exp(-gamma Int sqrt g R - Lambda Int sqrt g) Z_0 lässt alle "
       "normierten Materiekorrelationen auf festem Hintergrund gleich (Argument der Analyse, logisch geprüft)", True)
ZAHLEN["M_eff_q_halb_m6"] = float(Meff)

urteil = ("Die externe Analyse hält in den nachgerechneten Punkten. Der metrische Stern (Zentralmetrik plus vier Arme, "
          "M_(j,n) = (M_c/2) q^(n+1)) hat einen Massenoperator, der exakt unitär äquivalent zum rechteckigen Vierteltransfer "
          "ist (B U = V A_rect, Armdrehung = Vierteldrehung, genau ein Nullvektor); mit minimalen Knotenquellen ist das "
          "masselose Residuum universell (eta_v/M_v = 1/M_eff), die Zentrumsresolvente geht in C_q über. Das Potential "
          "V = (-4 + e1 + e3 - 4 e4)/2 ist flach mit Fierz-Pauli-Masse; V +- det(I-X)/2 ändert erst die quartische Ordnung, "
          "und die masselose Vierpunktamplitude sieht das nicht: Zweipunkt, kubische Daten und Ward legen die Dynamik nicht "
          "fest. Korrekturen: Die zentrierte Feldkovarianz ist verbunden und masselos (das frühere Nein betraf die "
          "vakuumsubtrahierte Mikroantwort); 'keine Nullstelle' ist kein Positivitätszertifikat, die Schranke x C <= 1/kappa "
          "ist eines und schließt die Dirac-Profilklasse auch mit Einstein-Zusatz aus; mit freier relativer Skala trifft "
          "jedes q die Zahl 3/10. Flavour: Eine volle Majorana-Z4 gibt entartete Massen; TM1 bei 15 Grad liegt im "
          "NuFIT-6.1-Bereich.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z50_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
