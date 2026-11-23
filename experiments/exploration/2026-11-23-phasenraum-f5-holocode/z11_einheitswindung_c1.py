"""Z11: Ist TFPTs Einheitswindung die erste Chernklasse des Bündels der Paarungsmasse Phi?

Frage (Folgeaufgabe 1 aus H6h8): Ein nichttriviales Bündel des Massenfelds erzwingt eine Nullstellenfläche
[Sigma] = PD(c1). Kann die Einheitswindung von P1 genau dieses c1 sein und so die Naht erzwingen?

A. Was TFPT unter Einheitswindung versteht (Zitate mit Datei und Zeile, im Skript nachgeprüft).
B. Die lokale Einheit im Nahtmodell von Z7: Windung von det M um die Naht = Index = Modenzahl.
C. Grad- und Paritätsbuchhaltung unter G (exakt): Wie viele verdrillte Durchstoßpunkte verlangt eine
   geschlossene Querfläche bei gegebenem c1?
D. Gitter-Ginzburg-Landau auf dem Torus mit Hintergrundbündel vom Grad d (Peierls-Phasen): erzwungene
   Nullstellen, robust aus Zufallsstarts, gezählt über Plakettenwindungen.
E. Phi mit freiem Blatt (Vorzeichen als Eichfreiheit, Schnitt unsichtbar): gewöhnlicher Wirbel gegen zwei
   verdrillte Defekte.
F. Dasselbe auf der Kugel (Würfeloberfläche), Lage von zwei und vier erzwungenen Defekten.
G. Urteil.
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
LINES = []
CHECKS = {}
BEFUNDE = {}
DATA = {}


def say(s=""):
    print(s)
    LINES.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    return "ja" if cond else "NEIN"


# ----------------------------------------------------------------------------------------------------------
say("A. Was TFPT unter Einheitswindung versteht\n")
QUOTES = [
    ("docs/research/TFPT-STAND.md", 26, "Randkern mit Einheitswindung",
     "P1 im Stand: orientierter, reflexionspositiver Randkern mit Einheitswindung."),
    ("documents/tfpt_1_architecture_e8/sections/001-source.tex", 72, "unit winding $[u_\\Sigma]=1$",
     "Architekturpapier: Randkern mit Einheitswindung [u_Sigma] = 1, primitives Postulat."),
    ("documents/archive--paper-latex--tfpt-42/sections/006-subsection-character-trace-identities-and-the-abelian-in.tex",
     253, "[u_\\Sigma]=1\\in K^1(S^1)\\cong \\mathbb Z.",
     "Original (tfpt-42): [u_Sigma] = 1 in K^1(S^1) = Z, also eine Schleifenklasse auf dem Nahtumlauf."),
    ("documents/archive--paper-latex--tfpt-42/sections/006-subsection-character-trace-identities-and-the-abelian-in.tex",
     280, "spectral flow $SF(U_\\Seam)=1$",
     "Dort: ein positiver Nahtumlauf ändert den Spektralfluss um eins; Delta = 2 pi, c3 = Delta/(16 pi^2)."),
    ("documents/archive--paper-latex--tfpt-42/sections/012-subsection-bosonic-relative-index-and-higgs-selection.tex",
     12, "compactified oriented normal sphere $S^2=D_N\\cup_{S^1}D_S$",
     "Original: L2 = det E2 über der kompaktifizierten orientierten Normalensphäre der Naht."),
    ("documents/archive--paper-latex--tfpt-42/sections/012-subsection-bosonic-relative-index-and-higgs-selection.tex",
     27, "c_1(L_2)=\\deg g_2=1",
     "Dort: c1(L2) = Grad der Verklebungsfunktion g2 = e^{i phi} = 1, c1(L3) = 0."),
    ("documents/archive--paper-latex--tfpt-42/sections/012-subsection-bosonic-relative-index-and-higgs-selection.tex",
     50, "half-angle representation",
     "Dort: der Spinlift wirkt im Zweierblock mit halbem Winkel, erst die Determinante windet einmal."),
    ("documents/archive--paper-latex--tfpt-42/sections/029-section-frw-reduction-seam-transfer-and-the-scalaron-bra.tex",
     205, "N_{\\mathrm{DW}}=n",
     "Original: Domänenwandzahl N_DW = n = Nahtwindung, also N_DW = 1."),
    ("documents/archive--paper-latex--tfpt-42/sections/011-subsection-family-monodromy-and-occupancy.tex",
     80, "N_{\\mathrm{corner}}\\cdot \\frac{\\pi}{2}=2\\pi",
     "Original: jede Ecke trägt eine Vierteldrehung des Spinlifts, vier Ecken ergeben eine volle Windung."),
    ("docs/research/fundamental-plan/TFPT-MARKIERTER-STABILISATOR-20261112.txt", 39, "pi2(C/K) ≅ pi1(K) ≅ Z.",
     "Plan 12.11.: markierte Einheitswindung = Endliniengrad in pi2(C/K) = pi1(K) = Z."),
    ("docs/research/fundamental-plan/TFPT-MARKIERTER-STABILISATOR-20261112.txt", 40,
     "Der Endliniengrad ist dieselbe ganze Zahl als erste Chernzahl",
     "Dort: der Endliniengrad ist die erste Chernzahl des Linienbündels über einer inneren Wurzelsphäre."),
    ("docs/research/fundamental-plan/TFPT-MARKIERTER-STABILISATOR-20261112.txt", 90,
     "Die Topologie bestimmt die kleinste nichttriviale Klasse",
     "Dort: welche Klasse physisch realisiert wird, muss eine Regel noch begründen."),
    ("experiments/theory-contracts/landkarte-hypothesen-20261119/check_h03_fluss.py", 42,
     "Fluss 1 auf dem Torus ist kein Rückzug von der Kugel",
     "Atlas 19.11.: Fluss 1 auf dem Markentorus ist kein Rückzug von der Kugel (Rückzüge haben geraden Grad)."),
    ("verification/status_ledger.csv", 135, "8pi=|mu4|*2pi seam winding",
     "Ledger SEAM.HORIZON.01: 8 pi = |mu4| mal 2 pi Nahtwindung."),
]
found = []
for path, line, needle, text in QUOTES:
    rows = (ROOT / path).read_text(encoding="utf-8").splitlines()
    ok = needle in rows[line - 1]
    found.append(ok)
    say(f"   {path}:{line}")
    say(f"      {text}")
check("Zitate an den angegebenen Zeilen vorhanden", all(found))
say("\n   Drei Orte für die Einheit. (i) Lokal: Klasse der Schleife auf dem Nahtumlauf, über Verklebung c1 = 1")
say("   auf der kompaktifizierten Normalensphäre (Determinantenlinie, Halbwinkel im Zweierblock).")
say("   (ii) Auf dem Nahtquotienten: vier Ecken mit je einer Vierteldrehung. (iii) Intern: Wurzelsphäre in C/K.")
say("   Keiner der drei Orte ist eine geschlossene Querfläche im Raum, die die Naht von außen umschließt.")

# ----------------------------------------------------------------------------------------------------------
say("\nB. Lokale Einheit im Nahtmodell von Z7\n")
_src = (HERE / "z7_ursprung_naht.py").read_text(encoding="utf-8")
_ns = {"__file__": str(HERE / "z7_ursprung_naht.py")}
exec(_src[:_src.index("L = 14")], _ns)
build = _ns["build"]
z7 = json.loads((HERE / "z7_ergebnis.json").read_text(encoding="utf-8"))


def det_winding(n_plus, n_minus, flux, copies, L=14, R=5):
    """Windung von det M (Produkt der Majoranamassen aller Arten, Wilson-Glied an der physischen Ecke null)."""
    H, (rad, ns) = build(L, 0.0, n_plus, n_minus, flux, copies=copies)
    dim = L * L * ns * 2
    D = H[:dim, dim:]
    idx = lambda x, y, a, s: ((x * L + y) * ns + a) * 2 + s
    c = (L - 1) / 2.0
    lo, hi = int(np.floor(c - R)), int(np.ceil(c + R))
    ring = [(x, lo) for x in range(lo, hi)] + [(hi, y) for y in range(lo, hi)] + \
           [(x, hi) for x in range(hi, lo, -1)] + [(lo, y) for y in range(hi, lo, -1)]
    dets = [np.prod([D[idx(x, y, a, 0), idx(x, y, a, 1)] for a in range(ns)]) for x, y in ring]
    ph = np.angle(np.array(dets + dets[:1]))
    return int(round(np.sum(np.angle(np.exp(1j * np.diff(ph)))) / (2 * np.pi)))


cases_b = [("verdrillter Defekt (J, Phase windet um pi), ein Paar", 1, 0, 0.25, 1),
           ("gewöhnlicher Wirbel (ohne J, Windung 1), ein Paar", 1, 1, 0.0, 1),
           ("nur J, keine Massenwindung", 0, 0, 0.25, 1),
           ("verdrillter Defekt, zwei Paare", 1, 0, 0.25, 2)]
match = []
for name, npl, nmi, fl, cp in cases_b:
    w = det_winding(npl, nmi, fl, cp)
    modes = z7[name]["net_core"]
    match.append(w == modes)
    say(f"   {name}: Windung von det M = {w}, chirale Moden aus Z7 = {modes}")
check("det-Windung um die Naht = Modenzahl aus Z7 (vier Fälle)", all(match))
say("   Die Klasse der Schleife M(phi) in K^1(S^1) ist die Windung von det M. Der verdrillte Defekt hat genau")
say("   die Einheit 1, der gewöhnliche Wirbel 2. Lokal wählt die Einheitswindung also den verdrillten Defekt.")
say("   Das ist die Lesart von [u_Sigma] = 1 und SF = 1 im Nahtmodell; c1 = 1 auf der Normalensphäre ist")
say("   dieselbe Zahl über die Verklebung. Sie setzt die Naht voraus und erzwingt sie nicht.")
say("   Dazu passt das Original: N_DW = n (Domänenwandzahl = Nahtwindung, also 1). Z9 fand für den verdrillten")
say("   Defekt auf dem Gitter genau eine Wand (Periode pi der Massenphase, ein Minimum je Periode).")
BEFUNDE["lokal: Einheit = det-Windung 1 = verdrillter Defekt (gewöhnlicher Wirbel hat 2)"] = True

# ----------------------------------------------------------------------------------------------------------
say("\nC. Grad- und Paritätsbuchhaltung unter G (exakt)\n")
say("   Ladungen: U(1)_V (gemeinsame Phase der Fermionen) und J (Z4). eta+ = (1, +1), eta- = (1, -1),")
say("   Phi = (-2, 2), Wilson eta+ eta- = (2, 0), det M ~ Phi^2 = (-4, 0). Ein verdrillter Durchstoßpunkt")
say("   trägt die J-Holonomie i auf eta+ und -i auf eta-; Phi windet dort um pi, det M einmal.")
say("   Auf einer geschlossenen Querfläche S mit k gleich orientierten Durchstoßpunkten und U(1)_V-Fluss f")
say("   (Fermioneinheiten) sind die Fermionlinien echt genau dann, wenn f + k/4 und f - k/4 ganz sind.\n")


def admissible(k, f, with_minus=True):
    ok = (f + Fr(k, 4)).denominator == 1
    if with_minus:
        ok = ok and (f - Fr(k, 4)).denominator == 1
    return ok


fs = [Fr(j, 8) for j in range(-24, 25)]
adm_k = sorted({k for k in range(0, 9) for f in fs if admissible(k, f)})
adm_k_plus_only = sorted({k for k in range(0, 9) for f in fs if admissible(k, f, with_minus=False)})
say(f"   zulässige k (mit beiden J-Ladungen): {adm_k}")
say(f"   zulässige k, wenn es nur eta+ gäbe:   {adm_k_plus_only}")
check("unter G ist die Zahl verdrillter Durchstoßpunkte gerade", adm_k == [0, 2, 4, 6, 8])
BEFUNDE["ohne Partner eta- wäre ein einzelner Durchstoßpunkt zulässig (G liefert aber beide J-Ladungen)"] = \
    1 in adm_k_plus_only
# Lesarten der Einheit als globale Zahl auf S. Phi-Windung gesamt = -2 f = k/2 (nur verdrillte Punkte).
readings = {}
for name, cond in (("R1: Grad der Determinantenlinie det M = 1", lambda k, f: k == 1),
                   ("R2: c1 des Phi-Bündels = 1", lambda k, f: Fr(k, 2) == 1),
                   ("R3: Vierteldrehungen summieren zu 2 pi (vier Ecken)", lambda k, f: Fr(k, 4) == 1),
                   ("R4: c1 der Paarungsphase in Fermioneinheiten = 1", lambda k, f: -f == 1)):
    ks = sorted({k for k in range(1, 9) for f in fs if admissible(k, f) and cond(k, f) and -2 * f == Fr(k, 2)})
    readings[name] = ks
    say(f"   {name}: verdrillte Durchstoßpunkte {ks if ks else 'keine zulässige Lösung'}")
check("R1 (eine einzelne Naht als Determinantengrad 1) ist unzulässig", readings[
    "R1: Grad der Determinantenlinie det M = 1"] == [])
check("R2 verlangt genau zwei, R3 und R4 genau vier Durchstoßpunkte",
      readings["R2: c1 des Phi-Bündels = 1"] == [2]
      and readings["R3: Vierteldrehungen summieren zu 2 pi (vier Ecken)"] == [4]
      and readings["R4: c1 der Paarungsphase in Fermioneinheiten = 1"] == [4])
say("   Grad von det M = 2 c1(Phi) ist immer gerade. Eine einzelne verdrillte Naht durch eine geschlossene")
say("   Querfläche ist unter G topologisch verboten. Das passt zum Atlas: Rückzüge haben geraden Grad.")
# Allgemein: M ist eine symmetrische Bilinearform auf dem Paarbündel E, det M liegt in (det E*)^2, Grad gerade.
# Nichtabelsche Holonomien helfen nicht: die Symmetrie von Phi psi chi auf dem Paar ist D4, J ist kein Kommutator.
sx_ = np.array([[0, 1], [1, 0]])
grp = []
for kk in range(16):
    th = kk * np.pi / 8
    rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    for g_ in (rot, rot @ np.diag([1, -1])):
        form = g_.T @ sx_ @ g_
        if np.allclose(form, sx_) or np.allclose(form, -sx_):
            grp.append(np.round(g_).astype(int))
comms = {tuple(int(c) for c in (a_ @ b_ @ np.linalg.inv(a_) @ np.linalg.inv(b_)).round().ravel())
         for a_ in grp for b_ in grp}
Jm = (0, -1, 1, 0)
say(f"   Symmetrie von psi chi (bis aufs Vorzeichen) in O(2): {len(grp)} Elemente (D4); Kommutatoren:"
    f" {sorted(comms)}; J darunter: {Jm in comms}")
check("D4-Holonomien umgehen die Parität nicht (J ist kein Kommutator)", len(grp) == 8 and Jm not in comms
      and comms == {(1, 0, 0, 1), (-1, 0, 0, -1)})
say("   Für jedes echte Paarbündel E gilt: M liegt in Sym^2 E*, det M in (det E*)^2, also ist die Summe der")
say("   det-Windungen -2 c1(E) gerade. Auch nichtabelsche Holonomien auf dem Paar ändern das nicht.")

# exhaustive Z2 parity on small closed lattices (Blatt-Monodromie = Z2-Eichfeld auf Links)
def torus(L):
    idx = lambda x, y: (x % L) * L + (y % L)
    links, pos, plaq = [], [], []
    for x in range(L):
        for y in range(L):
            links.append((idx(x, y), idx(x + 1, y)))
            links.append((idx(x, y), idx(x, y + 1)))
    lid = lambda x, y, d: 2 * idx(x, y) + d
    for x in range(L):
        for y in range(L):
            plaq.append(([lid(x, y, 0), lid(x + 1, y, 1), lid(x, y + 1, 0), lid(x, y, 1)], [1, 1, -1, -1]))
            pos.append((x + 0.5, y + 0.5))
    return L * L, np.array(links), plaq, np.array(pos, float)


def cube(N):
    pts = [p for p in itertools.product(range(N + 1), repeat=3) if any(c in (0, N) for c in p)]
    index = {p: i for i, p in enumerate(pts)}
    links, lindex = [], {}
    for p in pts:
        for ax in range(3):
            qn = list(p)
            qn[ax] += 1
            qn = tuple(qn)
            if qn in index:
                lindex[(index[p], index[qn])] = len(links)
                links.append((index[p], index[qn]))
    plaq, pos = [], []
    E = np.eye(3, dtype=int)
    for k in range(3):
        u, w = (k + 1) % 3, (k + 2) % 3
        for side in (0, N):
            uu, ww = (u, w) if side == N else (w, u)
            for a in range(N):
                for b in range(N):
                    c0 = np.zeros(3, int)
                    c0[k] = side
                    c0[uu], c0[ww] = a, b
                    corners = [c0, c0 + E[uu], c0 + E[uu] + E[ww], c0 + E[ww]]
                    ls, ss = [], []
                    for s in range(4):
                        pa, pb = index[tuple(corners[s])], index[tuple(corners[(s + 1) % 4])]
                        if (pa, pb) in lindex:
                            ls.append(lindex[(pa, pb)]); ss.append(1)
                        else:
                            ls.append(lindex[(pb, pa)]); ss.append(-1)
                    plaq.append((ls, ss))
                    ctr = c0 + 0.5 * (E[uu] + E[ww]) - N / 2.0
                    pos.append(ctr / np.linalg.norm(ctr))
    return len(pts), np.array(links), plaq, np.array(pos)


def z2_parities(nl, plaq, configs):
    out = np.ones((configs.shape[0], len(plaq)), int)
    for p, (ls, _) in enumerate(plaq):
        for l in ls:
            out[:, p] *= configs[:, l]
    return np.sum(out < 0, axis=1)


n3, l3, p3, _ = torus(3)
cfg = 1 - 2 * ((np.arange(2 ** len(l3))[:, None] >> np.arange(len(l3))) & 1)
par_t = z2_parities(len(l3), p3, cfg)
nc, lc, pc, _ = cube(1)
cfgc = 1 - 2 * ((np.arange(2 ** len(lc))[:, None] >> np.arange(len(lc))) & 1)
par_c = z2_parities(len(lc), pc, cfgc)
say(f"\n   Alle {len(cfg)} Z2-Blattfelder auf dem 3x3-Torus und alle {len(cfgc)} auf der Würfelkugel:")
say(f"   Zahl der Plaketten mit Monodromie -1 ist immer gerade: Torus {bool(np.all(par_t % 2 == 0))},"
    f" Kugel {bool(np.all(par_c % 2 == 0))} (vorkommend: {sorted(set(par_t.tolist()))})")
check("Blattmonodromie -1 tritt auf geschlossenen Flächen nur paarweise auf (erschöpfend)",
      np.all(par_t % 2 == 0) and np.all(par_c % 2 == 0) and par_t.max() > 0)

# ----------------------------------------------------------------------------------------------------------
say("\nD. Ginzburg-Landau auf dem Gitter-Torus mit Bündel vom Grad d\n")


def incidence(nl, plaq):
    B = np.zeros((len(plaq), nl))
    for p, (ls, ss) in enumerate(plaq):
        for l, s in zip(ls, ss):
            B[p, l] += s
    return B


def background(B, d):
    """Peierls-Phasen A mit gleichförmigem Fluss 2 pi d / P je Plakette; Diracfaden in Plakette 0."""
    P = B.shape[0]
    F = np.full(P, 2 * np.pi * d / P)
    rhs = F.copy()
    rhs[0] -= 2 * np.pi * d
    A = np.linalg.lstsq(B, rhs, rcond=None)[0]
    assert np.allclose(B @ A, rhs, atol=1e-8)
    return A, F


def energy(x, n, links, U, lam, v, nem=False, eps=1e-3):
    phi = x[:n] + 1j * x[n:]
    i, j = links[:, 0], links[:, 1]
    r = np.abs(phi) ** 2
    E = lam * np.sum((r - v * v) ** 2)
    g = 2 * lam * (r - v * v) * phi
    if not nem:
        diff = phi[j] - U * phi[i]
        E += np.sum(np.abs(diff) ** 2)
        a, b = diff, -np.conj(U) * diff
    else:                                   # Vorzeichen auf jedem Link frei: min_sigma |phi_j - sigma U phi_i|^2
        R = np.real(np.conj(phi[j]) * U * phi[i])
        s = np.sqrt(R * R + eps * eps)
        E += np.sum(r[i] + r[j] - 2 * s)
        t = R / s
        a, b = phi[j] - t * U * phi[i], phi[i] - t * np.conj(U) * phi[j]
    g = g + np.bincount(j, weights=a.real, minlength=n) + 1j * np.bincount(j, weights=a.imag, minlength=n)
    g = g + np.bincount(i, weights=b.real, minlength=n) + 1j * np.bincount(i, weights=b.imag, minlength=n)
    return E, 2 * np.concatenate([g.real, g.imag])


def relax(x, n, links, U, lam, v, nem=False, rounds=3):
    for _ in range(rounds):
        res = minimize(energy, x, args=(n, links, U, lam, v, nem), jac=True, method="L-BFGS-B",
                       options={"maxiter": 20000, "maxfun": 40000, "gtol": 1e-10, "ftol": 1e-16})
        x = res.x
    return res.x, float(res.fun)


def windings(phi, links, plaq, A, F, q, nem=False):
    """Plakettenwindungen von Phi (Ladung q) und Psi = Phi^2 (Ladung 2q); Z2-Monodromie des Blatts."""
    i, j = links[:, 0], links[:, 1]
    U = np.exp(1j * q * A)
    z = phi[j] * np.conj(U * phi[i])
    sig = np.sign(z.real) if nem else np.ones(len(z))
    w = np.angle(sig * z)
    wpsi = np.angle(phi[j] ** 2 * np.conj(U ** 2 * phi[i] ** 2))
    nphi, npsi, z2 = np.zeros(len(plaq)), np.zeros(len(plaq)), np.zeros(len(plaq), int)
    for p, (ls, ss) in enumerate(plaq):
        nphi[p] = (sum(s * w[l] for l, s in zip(ls, ss)) + q * F[p]) / (2 * np.pi)
        npsi[p] = (sum(s * wpsi[l] for l, s in zip(ls, ss)) + 2 * q * F[p]) / (2 * np.pi)
        z2[p] = int(np.prod([sig[l] for l in ls]) < 0)
    dev = np.abs(npsi - np.round(npsi)).max()
    return np.round(nphi * 2) / 2, np.round(npsi).astype(int), z2, dev


def corner_min(phi, plaq, links, p):
    sites = {s for l in plaq[p][0] for s in links[l]}
    return min(abs(phi[s]) for s in sites)


LT = 24
n, links, plaq, pos = torus(LT)
B = incidence(len(links), plaq)
lam, v = 0.25, 1.0
runs = {}
say(f"   Torus {LT} x {LT}, Energie sum |Phi_j - e^(i q A) Phi_i|^2 + lambda (|Phi|^2 - v^2)^2, lambda = {lam}.")
say("   Fünf Zufallsstarts je Fall, dazu ein Start mit konstantem Phi = v.\n")
cases_d = [(1, 0, lam), (1, 1, lam), (1, 1, 2.0), (2, 1, lam), (4, 1, lam), (1, -2, lam)]
for q, d, lm in cases_d:
    A, F = background(B, d)
    U = np.exp(1j * q * A)
    rec = []
    starts = [("zufall", s) for s in range(5)] + [("konstant", None)]
    for kind, seed in starts:
        if kind == "zufall":
            x0 = np.random.default_rng(seed).normal(size=2 * n) * 0.7
        else:                               # konstantes Phi = v mit Rauschen 1e-3 (bricht die Symmetrie)
            x0 = np.concatenate([np.full(n, v), np.zeros(n)]) + 1e-3 * np.random.default_rng(99).normal(size=2 * n)
        phi0 = x0[:n] + 1j * x0[n:]
        nphi0, _, _, _ = windings(phi0, links, plaq, A, F, q)
        x, E = relax(x0, n, links, U, lm, v)
        phi = x[:n] + 1j * x[n:]
        nphi, npsi, z2, dev = windings(phi, links, plaq, A, F, q)
        vort = np.nonzero(nphi)[0]
        rec.append({"start": kind, "E": E, "vor": int(np.sum(np.abs(nphi0))), "netto": float(nphi.sum()),
                    "windungen": sorted(nphi[vort].astype(int).tolist()), "zahl": int(len(vort)),
                    "kern_min": float(max([corner_min(phi, plaq, links, p) for p in vort], default=np.nan)),
                    "median": float(np.median(np.abs(phi))), "dev": float(dev)})
    runs[(q, d, lm)] = rec
    Es = [r["E"] for r in rec]
    outs = sorted({tuple(r["windungen"]) for r in rec})
    say(f"   q = {q}, d = {d:+d}, lambda = {lm}: Windungen am Minimum {[list(o) for o in outs]} (netto"
        f" {sorted({r['netto'] for r in rec})}, q d = {q * d:+d}); vorher {[r['vor'] for r in rec[:5]]}"
        f" Wirbelplaketten; E = {min(Es):.4f} bis {max(Es):.4f}")
    if len(rec[0]["windungen"]):
        say(f"      |Phi| am Kern höchstens {max(r['kern_min'] for r in rec):.3f} v, im Volumen (Median)"
            f" {min(r['median'] for r in rec):.3f} v")
DATA["torus_gl"] = {f"q={q},d={d},lambda={lm}": rec for (q, d, lm), rec in runs.items()}
ok_forced = all(all(r["netto"] == q * d and r["zahl"] == abs(q * d) and all(w == np.sign(q * d) for w in r["windungen"])
                    for r in rec) for (q, d, lm), rec in runs.items() if d != 0)
ok_zero = all(r["zahl"] == 0 for r in runs[(1, 0, lam)]) and min(r["vor"] for r in runs[(1, 0, lam)][:5]) > 20
ok_core = all(r["kern_min"] < 0.5 * r["median"] for (q, d, lm), rec in runs.items() if d != 0 and lm == lam
              for r in rec)
ok_energy = min(r["E"] for r in runs[(1, 1, lam)]) > max(r["E"] for r in runs[(1, 0, lam)]) + 1.0
ok_int = all(r["dev"] < 1e-8 for rec in runs.values() for r in rec)
say(f"\n   d = 0: alle Wirbel aus den Zufallsstarts verschwinden: {check('d = 0: Nullstellen verschwinden', ok_zero)}")
say(f"   d != 0: genau |q d| Elementarwirbel gleichen Vorzeichens in jedem Lauf, auch vom konstanten Start und bei"
    f" lambda = 2: {check('d != 0: genau |q d| erzwungene Elementarwirbel in jedem Lauf', ok_forced)}")
say(f"   Die Wirbel sind echte Nullstellen (|Phi| an den Plakettenecken unter halbem Volumenwert, lambda = {lam};"
    f" bei lambda = 2 ist der Kern kleiner als eine Plakette): {check('erzwungene Wirbel sind Nullstellen des Betrags', ok_core)}")
say(f"   Sie kosten Energie und bleiben trotzdem: {check('erzwungene Wirbel kosten Energie', ok_energy)};"
    f" Plakettenwindungen ganzzahlig auf 1e-8: {check('Peierls-Phasen konsistent: Plakettenwindungen ganzzahlig', ok_int)}")
say("   Lesart für H6: Ist Psi = Phi^2 (det M) das Feld mit Bündelgrad 2 c1(Phi), dann ist jeder Elementarwirbel")
say("   von Psi ein verdrillter Defekt. q d = 2 ist R2, q d = 4 ist R3/R4. q d = 1 wäre die einzelne Naht (R1),")
say("   die als Bündel von Psi allein existiert, aber nach C mit Phi und den Fermionen unter G unvereinbar ist.")

# Lage der vier erzwungenen Wirbel (q d = 4) auf dem quadratischen Torus
A4, F4 = background(B, 1)
best4 = min(range(5), key=lambda s: runs[(4, 1, lam)][s]["E"])
x0 = np.random.default_rng(best4).normal(size=2 * n) * 0.7
x4, _ = relax(x0, n, links, np.exp(4j * A4), lam, v)
phi4 = x4[:n] + 1j * x4[n:]
nphi4, _, _, _ = windings(phi4, links, plaq, A4, F4, 4)
P4 = pos[np.nonzero(nphi4)[0]]
devs = []
for a_, b_ in itertools.combinations(range(len(P4)), 2):
    dlt = (P4[b_] - P4[a_]) % LT
    devs.append(max(min(abs(c), abs(c - LT / 2), abs(c - LT)) for c in dlt))
half_dev = max(devs) if devs else np.nan
say(f"\n   Vier erzwungene Wirbel auf dem quadratischen Torus (Plakettenmitten) bei {P4.tolist()}; größte Abweichung der")
say(f"   Abstände vom Halbperiodengitter {{0, L/2}}^2: {half_dev:.1f} Gitterabstände.")
say("   Sie liegen auf einem Halbperiodengitter." if half_dev <= 1.0 else
    "   Sie liegen nicht auf dem Halbperiodengitter, sondern versetzt (rhombisch).")
BEFUNDE["vier erzwungene Wirbel sitzen auf einem Halbperiodengitter des Torus"] = bool(half_dev <= 1.0)

# ----------------------------------------------------------------------------------------------------------
say("\nE. Phi mit freiem Blatt: gewöhnlicher Wirbel gegen zwei verdrillte Defekte (c1(Phi) = 1)\n")
say("   Unter G ist das Vorzeichen von Phi das Blatt und die Lage des Schnitts eine Eichwahl. Auf dem Gitter")
say("   heißt das: Auf jedem Link ist das Vorzeichen frei, die Energie ist min_sigma |Phi_j - sigma U Phi_i|^2.")
A1, F1 = background(B, 1)
U1 = np.exp(1j * A1)
ords = []
for seed in range(3):
    x, E = relax(np.random.default_rng(seed).normal(size=2 * n) * 0.7, n, links, U1, lam, v)
    ords.append((E, x))
E_ord, x_ord = min(ords, key=lambda t: t[0])
_, npsi_ord, _, _ = windings(x_ord[:n] + 1j * x_ord[n:], links, plaq, A1, F1, 1)
E_ord_nem = float(energy(x_ord, n, links, U1, lam, v, nem=True)[0])
# Start mit zwei verdrillten Defekten: Wurzel des Psi-Minimums (Ladung 2, d = 1); das Vorzeichen ist frei
A2, F2 = background(B, 1)
xp, _ = relax(np.random.default_rng(7).normal(size=2 * n) * 0.7, n, links, np.exp(2j * A2), lam, v * v)
psi = xp[:n] + 1j * xp[n:]
phi_t = np.sqrt(psi.astype(complex))
x_t, E_tw = relax(np.concatenate([phi_t.real, phi_t.imag]), n, links, U1, lam, v, nem=True)
ph_t = x_t[:n] + 1j * x_t[n:]
_, npsi_t, z2_t, _ = windings(ph_t, links, plaq, A1, F1, 1, nem=True)
# gewöhnlicher Wirbel mit kleinem Rauschen unter freiem Blatt relaxiert
rng = np.random.default_rng(11)
x_o2, E_o2 = relax(x_ord + 0.05 * rng.normal(size=2 * n), n, links, U1, lam, v, nem=True)
ph_o2 = x_o2[:n] + 1j * x_o2[n:]
_, npsi_o2, z2_o2, _ = windings(ph_o2, links, plaq, A1, F1, 1, nem=True)
say(f"   gewöhnlicher Wirbel (Blatt fest): E = {E_ord:.4f}; Psi-Windungen {sorted(npsi_ord[npsi_ord != 0].tolist())};"
    f" mit freiem Blatt bewertet E = {E_ord_nem:.4f}")
say(f"   derselbe Wirbel mit freiem Blatt nachrelaxiert: E = {E_o2:.4f}; Psi-Windungen"
    f" {sorted(npsi_o2[npsi_o2 != 0].tolist())}, Blattmonodromien {int(z2_o2.sum())}")
say(f"   zwei verdrillte Defekte (Blatt frei): E = {E_tw:.4f}; Psi-Windungen {sorted(npsi_t[npsi_t != 0].tolist())},"
    f" Blattmonodromien {int(z2_t.sum())}")
pair_ok = sorted(npsi_t[npsi_t != 0].tolist()) == [1, 1] and int(z2_t.sum()) == 2 and \
    all(z2_t[p] == 1 for p in np.nonzero(npsi_t)[0])
check("c1(Phi) = 1 mit freiem Blatt: genau zwei verdrillte Defekte, je Psi-Windung 1", pair_ok)
check("zwei verdrillte Defekte haben weniger Energie als der gewöhnliche Wirbel",
      E_tw < min(E_ord, E_ord_nem, E_o2) - 0.5)
BEFUNDE["gewöhnlicher Wirbel zerfällt unter freiem Blatt von selbst in zwei verdrillte"] = \
    sorted(npsi_o2[npsi_o2 != 0].tolist()) == [1, 1]
DATA["blatt_frei"] = {"E_gewoehnlich": E_ord, "E_gewoehnlich_blatt_frei_bewertet": E_ord_nem,
                      "E_gewoehnlich_nachrelaxiert": E_o2, "E_zwei_verdrillt": E_tw,
                      "psi_windungen_verdrillt": sorted(npsi_t[npsi_t != 0].tolist()),
                      "psi_windungen_nachrelaxiert": sorted(npsi_o2[npsi_o2 != 0].tolist())}
say("   Mit c1(Phi) = 1 erzwingt das Bündel also zwei verdrillte Nähte, nicht eine; der gewöhnliche Wirbel")
say("   (Einheit 2, zwei Moden je Paar) ist teurer. Jede der beiden Nähte trägt lokal die Einheit 1.")

# ----------------------------------------------------------------------------------------------------------
say("\nF. Kugel (Würfeloberfläche, N = 10)\n")
NS = 10
ns_, links_s, plaq_s, pos_s = cube(NS)
Bs = incidence(len(links_s), plaq_s)
sph = {}
for q, d in ((1, 0), (1, 1), (2, 1), (4, 1)):
    As, Fs = background(Bs, d)
    Us = np.exp(1j * q * As)
    rec = []
    for seed in range(3):
        x, E = relax(np.random.default_rng(100 + seed).normal(size=2 * ns_) * 0.7, ns_, links_s, Us, lam, v)
        phi = x[:ns_] + 1j * x[ns_:]
        nphi, _, _, dev = windings(phi, links_s, plaq_s, As, Fs, q)
        rec.append((E, sorted(nphi[nphi != 0].astype(int).tolist()), pos_s[np.nonzero(nphi)[0]]))
    sph[(q, d)] = rec
    say(f"   q = {q}, d = {d}: Windungen {[r[1] for r in rec]}, E = {min(r[0] for r in rec):.4f}")
check("Kugel: genau |q d| Elementarwirbel, d = 0 ohne Wirbel",
      all(all(len(r[1]) == abs(q * d) and all(w == 1 for w in r[1]) for r in rec) for (q, d), rec in sph.items()))
E2, _, P2 = min(sph[(2, 1)], key=lambda t: t[0])
ang = float(np.degrees(np.arccos(np.clip(P2[0] @ P2[1], -1, 1))))
say(f"   Zwei erzwungene Defekte (q d = 2) stehen im Winkel {ang:.1f} Grad: Zur Naht im Nordpol der")
say("   Normalensphäre gehört ein zweiter verdrillter Punkt, bei Gegenlage im Unendlichen (D_S).")
BEFUNDE["Kugel: zwei erzwungene Defekte liegen gegenüber (Winkel > 150 Grad)"] = ang > 150
E4, _, P4s = min(sph[(4, 1)], key=lambda t: t[0])
pole = np.array([0.3, -0.5, 0.81])
pole /= np.linalg.norm(pole)
ax1 = np.cross(pole, [1.0, 0, 0])
ax1 /= np.linalg.norm(ax1)
ax2 = np.cross(pole, ax1)
zs = [(p @ ax1 + 1j * (p @ ax2)) / (1 - p @ pole) for p in P4s]
lamb = (zs[0] - zs[2]) * (zs[1] - zs[3]) / ((zs[0] - zs[3]) * (zs[1] - zs[2]))
jinv = 256 * (lamb ** 2 - lamb + 1) ** 3 / (lamb ** 2 * (lamb - 1) ** 2)
say(f"   Vier erzwungene Defekte (q d = 4): j-Invariante der Lage im Betrag {abs(jinv):.1f}")
say("   (Quadrat, harmonisch, mu4: 1728; Tetraeder: 0).")
BEFUNDE["Kugel: vier erzwungene Defekte liegen quadratisch (j nahe 1728)"] = bool(abs(jinv - 1728) < 200)
BEFUNDE["Kugel: vier erzwungene Defekte liegen tetraedrisch (j nahe 0)"] = bool(abs(jinv) < 200)
DATA["kugel"] = {"winkel_zwei": ang, "j_vier": [float(jinv.real), float(jinv.imag)],
                 "positionen_vier": P4s.tolist()}

# ----------------------------------------------------------------------------------------------------------
say("\nG. Urteil\n")
say("   1. Lokal ist die Identifikation sauber: [u_Sigma] = 1 in K^1(S^1) ist die Windung von det M um die")
say("      Naht, gleich dem Index und der Modenzahl je Paar. Sie wählt den verdrillten Defekt (1) und")
say("      schließt den gewöhnlichen Wirbel (2) aus. c1 = 1 auf der Normalensphäre ist dieselbe Zahl über die")
say("      Verklebung. Das ist eine Lesart, die die Naht voraussetzt; erzwungen wird nichts.")
say("   2. Erzwingen kann nur ein Bündel auf einer geschlossenen Querfläche. Das verlangt zweierlei, was TFPT")
say("      nicht enthält: einen nicht zusammenziehbaren 2-Zyklus quer zur Naht (kompakter Raum) und eine exakte,")
say("      geeichte Phase der Paarung. Auf dem H6-Gitter bricht das Wilson-Glied diese Phase explizit (Z9),")
say("      dann gibt es kein Bündel mit c1, nur das Z2 des Blatts.")
say("   3. Auch wenn man beides zugesteht: Unter G ist der Grad von det M auf jeder geschlossenen Fläche gerade.")
say("      Eine einzelne verdrillte Naht (R1) ist verboten. c1(Phi) = 1 erzwingt zwei verdrillte Nähte,")
say("      die Viereckenbilanz (Vierteldrehungen zu 2 pi) und c1 = 1 in Fermioneinheiten erzwingen vier.")
say("   4. Die Satzzahl (c = 8 gegen 24) bleibt davon unberührt: Jede erzwungene Naht trägt die Moden aller")
say("      windenden Sätze. Die Einheitswindung legt die Zahl der Nähte fest, nicht die Zahl der Sätze.")
say(f"   5. Lage: Zwei erzwungene Defekte stehen auf der Kugel im Winkel {ang:.0f} Grad. Vier liegen auf der Kugel bei")
say(f"      j = {abs(jinv):.1f} im Betrag (Quadrat 1728, Tetraeder 0), auf dem Torus "
    f"{'auf' if half_dev <= 1.0 else 'nicht auf'} den Halbperioden.")
if abs(jinv - 1728) > 200:
    say("      Die quadratische Lage der Marken folgt also nicht aus der Abstoßung erzwungener Defekte;")
    say("      in H6 kommt sie aus A* (E6 = 0 auf dem Nahttorus).")
say("   Urteil: Lesart lokal; gescheitert als Erzwingung genau einer Naht. Ein Ausweg verlangt eine Änderung")
say("   von G: ein Paar ohne Partner entgegengesetzter J-Ladung (dann fiele das J-invariante Wilson-Glied)")
say("   oder ein Paarungsfeld, dessen det M kein Quadrat ist. Beides ist hier nicht gerechnet.")

say("\nZusammenfassung")
for k_, c_ in CHECKS.items():
    say(f"   {'bestanden' if c_ else 'NICHT bestanden'}: {k_}")
for k_, c_ in BEFUNDE.items():
    say(f"   Befund ({'ja' if c_ else 'nein'}): {k_}")
(HERE / "ausgabe_z11.txt").write_text("\n".join(LINES) + "\n", encoding="utf-8")
(HERE / "z11_ergebnis.json").write_text(json.dumps({
    "checks": CHECKS, "befunde": BEFUNDE, "zulaessige_k": adm_k, "zulaessige_k_nur_eta_plus": adm_k_plus_only,
    "lesarten": readings, "torus_halbperioden_abweichung": float(half_dev), "daten": DATA,
    "urteil": "lokal Lesart (Einheit = det-Windung = Index, wählt den verdrillten Defekt); als Erzwingung genau "
              "einer Naht gescheitert (Paritätsschranke unter G, fehlende geschlossene Querfläche und explizite "
              "Brechung der Paarungsphase); c1(Phi) = 1 erzwingt zwei, Viereckenbilanz vier verdrillte Nähte"},
    ensure_ascii=False, indent=1, default=float), encoding="utf-8")
assert all(CHECKS.values())
