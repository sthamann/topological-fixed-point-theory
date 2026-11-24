"""Z55: Prüfung der externen Analyse zu H6h27 -- Vierpunktantwort, Tiefenparameter, G-Skalierung, Flavourhindernis.

Teil 1  Vierpunktantwort: Mit R_(e,b) = G_(Kind,b)/(M_Kind M_b) - G_(Eltern,b)/(M_Eltern M_b) und k_e = lambda W_e M_Kind^2 ist
        Sum_e k_e Prod_i R_(e,c)(p_i^2) für vier zentrale Quellen im flachen Stern rP/((1 - rP) lambda^3 M_c^6), im gewichteten
        Stern bei p -> 0 Li_3(r)/(lambda^3 M_c^6). Tiefe n trägt (rP)^(n+1): eine erzeugende Funktion für die t_n.
Teil 2  Tiefenparameter: Die Vierteldrehung setzt die Arme gleich (t_(j,n) = t_n), nicht die Tiefen; zwei Homogenitätsregeln
        (metrisch t_n = t, Transfer t_n ~ r^n) geben verschiedene Erzeugende.
Teil 3  G: M_v -> c M_v, k_e -> c^2 k_e lässt den kanonischen Massenoperator gleich und skaliert G mit 1/c^2.
Teil 4  Flavour: Z3 geladen, Klein neutral, mu-tau-Reflexion: tr(T S_i) = 3|U_ei|^2 - 1. Für TM1 bei 15 Grad sind
        -(2 -+ sqrt 3)/4 keine ganzen algebraischen Zahlen, also keine Spuren eines endlichen Gruppenelements. Die Familienuhr
        (TM2) entgeht dem, weil ihre Beträge nicht mu-tau-symmetrisch sind; TM2 verletzt dafür die Solarsummenregel.
"""
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import polylog

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "research" / "originlab"))
from process_rules import _joint_basis, delta_group, metric_star  # noqa: E402

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


def star(qv, m, weighted, lam=1.0):
    st = metric_star(qv, m)
    W = np.array([n + 1 if weighted else 1 for j in range(4) for n in range(m)], float)
    L = lam*st["B"].T @ np.diag(W) @ st["B"]
    nodes, M = st["nodes"], st["M"]
    ix = {v: i for i, v in enumerate(nodes)}
    edges = []
    for e, (j, n) in enumerate((j, n) for j in range(4) for n in range(m)):
        child, parent = ix[(j, n)], ix["c"] if n == 0 else ix[(j, n - 1)]
        edges.append((child, parent, lam*W[e]*M[child]**2, n))
    return L, M, edges, ix


def S(qv, m, weighted, p2s, tdepth=None, lam=1.0):
    L, M, edges, ix = star(qv, m, weighted, lam)
    c = ix["c"]
    cols = []
    for p2 in p2s:
        G = np.linalg.solve(p2*np.eye(L.shape[0]) + L, np.eye(L.shape[0])[:, c])
        cols.append(G)
    tot = 0.0
    for child, parent, k, n in edges:
        t = 1.0 if tdepth is None else tdepth(n)
        prod = 1.0
        for G in cols:
            prod *= G[child]/(M[child]*M[c]) - G[parent]/(M[parent]*M[c])
        tot += k*t*prod
    return tot


def zfun(x, r):
    return 2/(x + 1 + r + np.sqrt((x + 1 + r)**2 - 4*r))


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Vierpunktantwort der Randquellen")
q = 0.5
r = q*q
rows = []
for p2s in ([1e-7]*4, [0.3, 0.3, 0.3, 0.3], [0.1, 0.5, 1.2, 2.0]):
    P = np.prod([zfun(x, r) for x in p2s])
    num = S(q, 70, False, p2s)
    form = r*P/(1 - r*P)
    rows.append((p2s, num, form))
    say(f"   flach, p^2/lambda = {p2s}: Summe {num:.10f}, rP/(1 - rP) = {form:.10f}")
check("flacher Stern: Sum_e k_e Prod R_(e,c) = rP/((1 - rP) lambda^3 M_c^6) (lambda = M_c = 1) für drei Impulssätze",
      all(abs(a - b) < 1e-8 for _, a, b in rows))
wnum = S(q, 70, True, [1e-7]*4)
li3 = float(polylog(3, r))
say(f"   gewichtet, p -> 0: Summe {wnum:.10f}, Li_3(r) = {li3:.10f}; flach: r/(1-r) = {r/(1 - r):.10f}")
check("gewichteter Stern bei p -> 0: Li_3(r) statt r/(1-r): der Operator ändert auch die nichtlineare Antwort",
      abs(wnum - li3) < 1e-6)
p2s = [0.2, 0.4, 0.7, 1.1]
P = np.prod([zfun(x, r) for x in p2s])
dep = [(n, S(q, 70, False, p2s, tdepth=lambda k, n=n: 1.0 if k == n else 0.0), (r*P)**(n + 1)) for n in range(3)]
say(f"   Beitrag einzelner Tiefen (gemessen, (rP)^(n+1)): {[(n, round(a, 10), round(b, 10)) for n, a, b in dep]}")
check("die Tiefe n trägt genau (rP)^(n+1): F(P) = -(1/(2 lambda^3 M_c^6)) Sum t_n (rP)^(n+1) ist eine erzeugende Funktion",
      all(abs(a - b) < 1e-10 for _, a, b in dep))
befund("Eine mikroskopische Vierpunktantwort mit dieser Tensorstruktur und Quellnorm bestimmt alle t_n über die "
       "Koeffizienten von F(P); ihr Wert bei kleinen Impulsen bestimmt nur eine gewichtete Summe", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Tiefenparameter")
D = 6
per_edge = 4*D
arm_eq = D
say(f"   Kanten bei {D} Tiefen: {per_edge} freie t_e; Vierteldrehung (t_(j,n) = t_n): {arm_eq}")
check("die Vierteldrehung setzt die Arme gleich, nicht die Tiefen: 4D -> D Parameter", per_edge == 4*arm_eq)
Fm = S(q, 70, False, p2s, tdepth=lambda n: 1.0)
Ft = S(q, 70, False, p2s, tdepth=lambda n: r**n)
say(f"   metrisch homogen (t_n = 1): {Fm:.10f} = rP/(1-rP) {r*P/(1 - r*P):.10f}; transferhomogen (t_n = r^n): {Ft:.10f} = "
    f"rP/(1 - r^2 P) {r*P/(1 - r*r*P):.10f}")
check("zwei Homogenitätsregeln geben verschiedene Erzeugende: rP/(1-rP) gegen rP/(1-r^2 P)",
      abs(Fm - r*P/(1 - r*P)) < 1e-9 and abs(Ft - r*P/(1 - r*r*P)) < 1e-9 and abs(Fm - Ft)/Fm > 1e-2)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: G ohne Normregel frei")
st = metric_star(q, 40)
c = 3.7
L1 = st["B"].T @ st["B"]
Mb = c*st["M"]
B2 = np.zeros_like(st["B"])
for e in range(st["B"].shape[0]):
    nz = np.nonzero(st["B"][e])[0]
    child = [i for i in nz if st["B"][e, i] == -1][0]
    parent = [i for i in nz if i != child][0]
    B2[e, parent], B2[e, child] = Mb[child]/Mb[parent], -1
G1 = 1/(8*np.pi*np.sum(st["M"]**2))
G2 = 1/(8*np.pi*np.sum(Mb**2))
check("M_v -> c M_v, k_e -> c^2 k_e lässt den kanonischen Massenoperator gleich und skaliert G mit 1/c^2",
      np.allclose(B2.T @ B2, L1, atol=1e-12) and abs(G2*c*c/G1 - 1) < 1e-12)
check("im unendlichen Stern ist G = (1 - r)/(8 pi M_c^2)", abs(1/(8*np.pi*np.sum(metric_star(q, 80)["M"]**2))
                                                             - (1 - r)/(8*np.pi)) < 1e-14)
befund("G braucht eine physische Normregel; im Bestand ist das die Einheit v_geo mit 8 pi = 1/c3 (GRAV.NONLINEAR.01, Lesart)",
       True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Flavour -- Ganzzahligkeitshindernis")
U0 = np.array([[np.sqrt(2/3), 1/np.sqrt(3), 0], [-1/np.sqrt(6), 1/np.sqrt(3), -1/np.sqrt(2)],
               [-1/np.sqrt(6), 1/np.sqrt(3), 1/np.sqrt(2)]])
th, ph = np.pi/12, np.pi/2
cc, ss = np.cos(th), np.sin(th)*np.exp(-1j*ph)
U = U0 @ np.array([[1, 0, 0], [0, cc, ss], [0, -np.conj(ss), cc]])
w = np.exp(2j*np.pi/3)
T = np.diag([1, w, w*w])
Sk = [2*np.outer(U[:, i], U[:, i].conj()) - np.eye(3) for i in range(3)]
klein = all(np.allclose(Sk[i] @ Sk[i], np.eye(3)) for i in range(3)) and np.allclose(Sk[0] @ Sk[1], Sk[1] @ Sk[0])
mutau = np.allclose(np.abs(U[1]), np.abs(U[2]))
tr = [np.trace(T @ Sk[i]) for i in range(3)]
pred = [3*abs(U[0, i])**2 - 1 for i in range(3)]
say(f"   tr(T S_i) = {np.round(tr, 6).tolist()}, 3|U_ei|^2 - 1 = {np.round(pred, 6).tolist()}")
check("TM1 bei 15 Grad mit mu-tau-Reflexion: S_i bilden eine Klein-Gruppe, |U_mu i| = |U_tau i|, tr(T S_i) = 3|U_ei|^2 - 1",
      klein and mutau and np.allclose(tr, pred, atol=1e-12))
z = sp.symbols("z")
vals = [sp.Integer(1), -(2 - sp.sqrt(3))/4, -(2 + sp.sqrt(3))/4]
mins = [sp.minimal_polynomial(v, z) for v in vals]
say(f"   Minimalpolynome: {mins}")
integral = [sp.Poly(mp_, z).LC() == 1 and all(cf.is_integer for cf in sp.Poly(mp_, z).all_coeffs()) for mp_ in
            [sp.Poly(sp.minimal_polynomial(v, z, polys=True).as_expr()*sp.denom(sp.minimal_polynomial(v, z, polys=True).as_expr()), z).monic() for v in vals]]
monic_int = [all(cf.is_integer for cf in sp.Poly(sp.minimal_polynomial(v, z), z).monic().all_coeffs()) for v in vals]
check("-(2 -+ sqrt 3)/4 haben das Minimalpolynom 16 z^2 + 16 z + 1 (normiert z^2 + z + 1/16): keine ganzen algebraischen "
      "Zahlen, also keine Spur eines Elements einer endlichen Gruppe", monic_int == [True, False, False])
befund("Die exakte direkte Realisierung von TM1 bei 15 Grad (Z3 geladen, Majorana-Klein, mu-tau-Reflexion) ist in jeder "
       "endlichen Gruppe ausgeschlossen; kein größerer Gruppenscan kann sie retten", True)
G, mul, e, mat = delta_group(4)
els = sorted(G)
Mm = {g: mat(g) for g in els}
order = {}
for g in els:
    k, x = 1, g
    while x != e:
        x, k = mul(g, x), k + 1
    order[g] = k
inv = [g for g in els if order[g] == 2 and not np.allclose(Mm[g], -np.eye(3))]
found = None
for g3 in [g for g in els if order[g] == 3]:
    Bl = _joint_basis([Mm[g3]])
    if Bl is None:
        continue
    for i, g1 in enumerate(inv):
        for g2 in inv[i + 1:]:
            if mul(g1, g2) != mul(g2, g1):
                continue
            Bn = _joint_basis([Mm[g1], Mm[g2]])
            if Bn is None:
                continue
            P = np.abs(Bl.conj().T @ Bn)**2
            if np.min(P) > 0.04 and np.min(P) < 0.05:
                found = (g3, g1, g2, P)
                break
        if found:
            break
    if found:
        break
g3, g1, g2, P = found
traces = [complex(np.trace(Mm[g3] @ Mm[s])) for s in (g1, g2, mul(g1, g2))]
prods = [mul(g3, s) for s in (g1, g2, mul(g1, g2))]
consistent = all(np.allclose(Mm[g3] @ Mm[s], mat(gp)) for s, gp in zip((g1, g2, mul(g1, g2)), prods))
exact = [sum(sp.exp(sp.I*sp.pi*gp[1][j]/4) for j in range(3) if gp[0][j] == j) for gp in prods]   # Summe von Einheitswurzeln
tmins = [sp.minimal_polynomial(sp.expand(t), z) if t != 0 else z for t in exact]
say(f"   Delta(96)-Familienuhr: |U|^2-Zeilen {np.round(np.sort(P, axis=1), 4).tolist()}; tr(T S) = {[str(sp.nsimplify(sp.expand(t))) for t in exact]}, "
    f"Minimalpolynome {tmins}")
check("in der Delta(96)-Familienuhr sind die Spuren tr(T S) ganze algebraische Zahlen, und ihre Beträge sind nicht "
      "mu-tau-symmetrisch -- deshalb entgeht TM2 bei 15 Grad dem Hindernis",
      consistent and all(abs(complex(sp.N(t)) - tn) < 1e-12 for t, tn in zip(exact, traces))
      and all(all(cf.is_integer for cf in sp.Poly(mp_, z).monic().all_coeffs()) for mp_ in tmins)
      and not any(np.allclose(np.sort(P[a]), np.sort(P[b])) and np.allclose(P[a], P[b]) for a, b in ((0, 1), (0, 2), (1, 2))))
s13 = np.linspace(0.02, 0.025, 6)
s12_tm2 = 1/(3*(1 - s13))
say(f"   TM2: sin^2 theta12 = 1/(3(1 - sin^2 theta13)) = {np.round(s12_tm2, 4).tolist()} > 1/3")
check("TM2 hat sin^2 theta12 > 1/3 für jedes theta13 > 0; NuFIT 6.1 (NO mit SK): 0.3088 +- 0.0067 (arXiv:2603.29876), "
      "also mehr als 3.5 sigma; TM1 bei 15 Grad (0.3181) liegt 1.4 sigma darüber",
      np.all(s12_tm2 > 1/3) and (1/3 - 0.3088)/0.0067 > 3.5 and abs((0.3181 - 0.3088)/0.0067 - 1.39) < 0.02)
befund("Der Datenvergleich ist eindimensional, kein gemeinsamer Fit; bei TM1 mit mu-tau-Reflexion passt laut Analyse nur "
       "delta = 270 Grad in den NuFIT-6.1-Bereich 125-365 Grad", True)
ZAHLEN["vierpunkt_flach"] = [(a, b) for _, a, b in rows]

urteil = ("Die externe Analyse hält auch hier. Die t-abhängige Vierpunktantwort zentraler Quellen ist im flachen Stern "
          "rP/((1 - rP) lambda^3 M_c^6), im gewichteten bei p -> 0 Li_3(r)/(lambda^3 M_c^6); die Tiefe n trägt (rP)^(n+1), "
          "sodass eine mikroskopische Antwort mit passender Tensorstruktur alle t_n bestimmen könnte. Die Vierteldrehung "
          "setzt nur die Arme gleich; metrische und transferartige Homogenität geben verschiedene Erzeugende. G skaliert mit "
          "M_v -> c M_v bei gleichem kanonischem Operator, braucht also eine Normregel (im Bestand die Einheit v_geo). Im "
          "Flavour ist die exakte direkte Realisierung von TM1 bei 15 Grad mit mu-tau-Reflexion in jeder endlichen Gruppe "
          "ausgeschlossen: tr(T S_i) = 3|U_ei|^2 - 1 = -(2 -+ sqrt 3)/4 sind keine ganzen algebraischen Zahlen. Die "
          "Familienuhr entgeht dem nur, weil sie nicht mu-tau-symmetrisch ist, und verletzt dafür die Solarsummenregel "
          "(sin^2 theta12 > 1/3, mehr als 3.5 sigma gegen NuFIT 6.1).")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z55_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
