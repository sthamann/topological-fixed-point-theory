"""Z17: Sektorwahl. Gibt es eine Regel, die den primitiven Sektor (Ladung bzw. Windung 1) und genau einen
windenden 16er-Block als Grundzustandseigenschaft wählt, statt beides vorauszusetzen?

Offene Verbindung 1 aus H6h10. Geprüft werden sieben Kandidaten:
 A. Quellen (K5): Was TFPT selbst über die Einheit sagt; Hopf-Hebung und Linsenraum X = S^3/<D> (v957).
 B. K1  Verzweigte Überlagerung: Ist das Ordnungsfeld das Datum w^2 = q der Doppelüberlagerung?
 C. K2  Spin und Statistik: Erzwingt ein spinorielles Phi eine ungerade Ladung von q = Phi^2?
 D. K3  Energie und Entropie der Hopf-Sektoren k = 0..3 (Gitter-Ginzburg-Landau, Hopf-reduziert).
 E. K3' Nullmittel: Ordnungsfeld ohne Hintergrundwert; tiefste zulässige Mode l = 1 auf S^3.
 F. K4  Kopiedeterminante: Smith-Normalform, Generizität, J-Symmetrie, Familiencharaktere.
 G. K6  Gemeinsame GSO-Summe (A*) mit n windenden Blöcken: D_{8n}^+ statt E8^n.
 H. Urteil.
Alles mit .venv/bin/python; Laufzeit etwa eine Minute.
"""
import itertools
import json
import textwrap
from math import comb
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
LINES, CHECKS, BEFUNDE, DATA = [], {}, {}, {}


def say(s=""):
    print(s)
    LINES.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    return "ja" if cond else "NEIN"


def find_line(rel, needle):
    p = ROOT / rel
    if not p.exists():
        return None
    for i, line in enumerate(p.read_text(errors="replace").splitlines(), 1):
        if needle in line:
            return i
    return None


# ==========================================================================================================
say("A. Quellen (K5): Was TFPT über die Einheit sagt\n")
T42 = "documents/archive--paper-latex--tfpt-42/sections/"
QUOTES = [
    ("documents/tfpt_1_architecture_e8/sections/001-source.tex", "primitive postulate",
     "P1 im Architekturpapier: Einheitswindung [u_Sigma] = 1, Status 'primitive postulate'."),
    (T42 + "006-subsection-character-trace-identities-and-the-abelian-in.tex", "[u_\\Sigma]=1\\in K^1(S^1)",
     "tfpt-42: [u_Sigma] = 1 in K^1(S^1) = Z."),
    (T42 + "012-subsection-bosonic-relative-index-and-higgs-selection.tex", "c_1(L_2)=\\deg g_2=1",
     "tfpt-42: c1(L2) = deg g2 = 1 auf der kompaktifizierten Normalensphäre."),
    (T42 + "029-section-frw-reduction-seam-transfer-and-the-scalaron-bra.tex", "N_{\\mathrm{fam}}=3n,",
     "tfpt-42, Satz 'Boundary winding controls family': N_fam = 3n, N_DW = n."),
    (T42 + "029-section-frw-reduction-seam-transfer-and-the-scalaron-bra.tex", "canonical rigid branch has $n=1$",
     "Dort: n = 1 ist der kanonische Zweig, gesetzt über den primitiven Generator."),
    ("documents/tfpt_research_contracts/sections/014-source.tex", "the seam sphere \\emph{is} the Hopf base",
     "Forschungsverträge (v957): die Nahtsphäre ist die Hopf-Basis von S^3."),
    ("verification/status_ledger.csv", "four marks lift to four Hopf fibre circles",
     "Ledger CCC.SEAM.KINEMATICS.01 (4): vier Marken = vier Hopf-Fasern, in zwei Blattbahnen, verschlungen."),
    ("verification/status_ledger.csv", "pole strings carry the full mu_4 monodromy",
     "Ledger CCC.SEAM.KINEMATICS.01 (6): Polfäden tragen chi(D) = i, Nahtfäden chi(D^2) = -1 [C]."),
    ("experiments/tfpt-discovery/ccc_uplift_4d_probe.py", "TYPED FORK (honest, not decided)",
     "Probe v957: offene Gabel, Deck auf der Raumzeit (X = S^3/Z4) oder intern (Raum S^3)."),
    ("docs/research/fundamental-plan/TFPT-MARKIERTER-STABILISATOR-20261112.txt", "physische Regel muss noch begr",
     "Plan 12.11.: Welche Klasse physisch realisiert wird, muss eine Regel noch begründen."),
]
found = []
for rel, needle, comment in QUOTES:
    ln = find_line(rel, needle)
    found.append(ln is not None)
    say(f"   {rel}:{ln}")
    say(f"      {comment}")
check("A: alle zitierten Stellen im Repo gefunden", all(found))
say("\n   Befund: TFPT setzt die Einheit (n = 1) als primitives Postulat. Der Satz N_fam = 3n verknüpft sie")
say("   mit der Familienzahl. Die Hopf-Hebung (v957) legt die vier Marken als vier Hopf-Fasern auf S^3.")

say("\n   Rechnung: Deck D = diag(i, -i) und Uhr rho = diag(i, 1) aus der Probe v957 auf Monomen z1^a z2^b.")
say("   Hopf-Ladung k = a + b, Deck-Charakter chi(D) = i^(a - b).")
inv_by_k = {}
for k in range(1, 5):
    mons = [(a, k - a) for a in range(k, -1, -1)]
    chars = {m: sp.I ** (m[0] - m[1]) for m in mons}
    inv = [m for m in mons if sp.simplify(chars[m] - 1) == 0]
    inv_by_k[k] = inv
    say(f"   k = {k}: Charaktere {[str(sp.simplify(chars[m])) for m in mons]}; D-invariant: "
        f"{['z1^%d z2^%d' % m for m in inv] or 'keine'}")
check("A: kein D-invariantes Feld mit k = 1 oder 3; kleinstes invariantes k = 2 (z1 z2)",
      inv_by_k[1] == [] and inv_by_k[3] == [] and inv_by_k[2] == [(1, 1)])
chi1 = sorted(str(sp.I ** (a - (1 - a))) for a in (0, 1))
say(f"   Ein Feld mit k = 1 trägt chi(D) in {chi1}: das ist die Holonomie eines Spin-Z4-Fermions, nicht eines")
say("   J-neutralen Bosons. z1 verschwindet auf der Polfaser z = 0; dort sitzt nach Ledger (6) chi(D) = i.")
check("A: Ladung-1-Felder tragen chi(D) = +-i", chi1 == ["-I", "I"])

# clock-invariant quartic, marks, D-orbits, linking
z = sp.symbols("z")
clock_inv = [m for m in inv_by_k[4] if sp.simplify(sp.I ** m[0] - 1) == 0]
s_marks = sp.Poly(z ** 4 - 1, z)                     # s = z1^4 - z2^4 in z = z1/z2
roots = [complex(r) for r in sp.nroots(s_marks.as_expr())]
orbits = {frozenset((round(r.real, 9) + 1j * round(r.imag, 9), round(-r.real, 9) + 1j * round(-r.imag, 9))) for r in roots}
say(f"   k = 4: uhrinvariant {['z1^%d z2^%d' % m for m in clock_inv]}; s = z1^4 - z2^4 hat einfache Nullfasern über")
say(f"   z in mu4; D bildet z -> -z ab, also {len(orbits)} Bahnen: im Linsenraum X zwei Nahtfäden.")
rng = np.random.default_rng(17)
n_str = []
for _ in range(50):
    a_, b_, c_ = rng.normal(size=3) + 1j * rng.normal(size=3)
    us = np.roots([a_, b_, c_])                       # a u^2 + b u + c, u = z^2
    zs = np.concatenate([np.sqrt(us), -np.sqrt(us)])
    orb = {tuple(sorted([np.round(zz, 8), np.round(-zz, 8)], key=lambda w: (w.real, w.imag))) for zz in zs}
    n_str.append(len(orb))
check("A: generisches D-invariantes q mit k = 4 gibt im Linsenraum genau zwei Nahtfäden",
      len(orbits) == 2 and set(n_str) == {2})


def fibre(m, t):
    v = np.array([m, 1.0], complex) / np.sqrt(abs(m) ** 2 + 1)
    return np.exp(1j * t)[:, None] * v[None, :]


tt = np.linspace(0, 2 * np.pi, 4001)
fi = fibre(1j, tt)
vals = fi[:, 0] - 1.0 * fi[:, 1]                     # vanishes on the fibre over z = 1
lk = np.sum(np.diff(np.unwrap(np.angle(vals)))) / (2 * np.pi)
say(f"   Verschlingung der Markenfasern über 1 und i (Windung von z1 - z2 längs der Faser über i): {lk:+.4f}")
check("A: Markenfasern paarweise einfach verschlungen", abs(abs(lk) - 1) < 1e-9)
say("   Folgerung (exakt, unter G): q = Phi^2 = det M ist J-neutral. Wirkt das Deck auf den Raum (X = S^3/<D>),")
say("   muss q D-invariant sein. Dann ist k gerade, und Nahtfäden kommen in X paarweise (wie Z11).")
say("   Eine einzelne Naht (k = 1) ist nur möglich, wenn das Deck intern ist und der Raum S^3 bleibt.")
BEFUNDE["Deck auf der Raumzeit: J-neutrales q hat gerades k, einzelne Naht in X unmöglich"] = True

say("\n   Familien je Naht (Lesart C'' je Naht): Index von K^(1/2) x O(D) ist H^0(K(D)) = 3 je Markensphäre.")
P_ = sp.symbols("p0:4")
t = sp.symbols("t")
ok_idx = []
for deg in (2, 3):
    Pz = sum(P_[i] * z ** i for i in range(deg + 1))
    expr = sp.simplify(-Pz.subs(z, 1 / t) / t ** 2 / ((1 / t) ** 4 - 1))
    regular = sp.limit(expr.subs({P_[i]: 1 for i in range(deg + 1)}), t, 0)
    ok_idx.append(regular.is_finite)
check("A: P(z)dz/(z^4-1) im Unendlichen regulär genau für deg P <= 2 (h^0 = 3)", ok_idx == [True, False])
Nnu, dNnu = 2.9840, 0.0082
pulls = {k: (3 * k - Nnu) / dNnu for k in (1, 2, 3)}
say(f"   Index additiv über k Nähte: N_fam = 3k, wie TFPTs Satz N_fam = 3n. Mit N_nu = {Nnu} +- {dNnu}")
say(f"   (LEP, externe Zahl, nicht im Repo geprüft): Abstand {', '.join(f'k = {k}: {p:+.0f} sigma' for k, p in pulls.items())}.")
check("A: nur k = 1 verträglich mit N_nu (externe Zahl, Lesart Familien je Naht)",
      abs(pulls[1]) < 5 and all(abs(pulls[k]) > 100 for k in (2, 3)))
DATA["N_nu_pulls"] = {str(k): float(p) for k, p in pulls.items()}

# ==========================================================================================================
say("\nB. K1: Das Ordnungsfeld als Datum einer verzweigten Doppelüberlagerung w^2 = q\n")
say("   q sei ein Schnitt von O(k) auf der Normalensphäre (Hopf-Ladung k). w = sqrt(q) = Phi lebt auf der")
say("   Doppelüberlagerung mit r = k + (k mod 2) Verzweigungspunkten, Geschlecht g = r/2 - 1.")
gen = {}
for k in range(1, 9):
    r = k + (k % 2)
    gg = sp.solve(sp.Eq(2 * sp.Symbol("g") - 2, 2 * (-2) + r), sp.Symbol("g"))[0]
    gen[k] = int(gg)
say(f"   Geschlecht je k: {gen}")
torus_k = [k for k, g_ in gen.items() if g_ == 1]
check("B: Torus (nirgends verschwindende 1-Form) nur für k in {3, 4}", torus_k == [3, 4])
w = sp.symbols("w")
f4 = sp.series(-1 / sp.sqrt(1 - t ** 4), t, 0, 3).removeO()   # dz/w at infinity for w^2 = z^4 - 1
say(f"   k = 4, w^2 = z^4 - 1: dz/w bei z = 1/t ist {f4} dt + ..., regulär und ohne Nullstelle.")
k1form = sp.simplify(sp.diff(w ** 2, w) / w)                # z = w^2: dz/w = 2 dw
inf_order = sp.simplify((k1form * sp.diff(1 / t, t)).subs(w, 1 / t))
say(f"   k = 1, w^2 = z: dz/w = {k1form} dw; bei w = 1/t: {inf_order} dt, Pol zweiter Ordnung.")
check("B: k = 1 hat keine holomorphe 1-Form, k = 4 eine ohne Nullstelle",
      sp.limit(inf_order * t ** 2, t, 0) != 0 and f4.subs(t, 0) == -1)
a4, b4 = -1, 0
j33 = sp.Rational(1728) * 4 * a4 ** 3 / (4 * a4 ** 3 + 27 * b4 ** 2)
say(f"   k = 3 mit Punkt im Unendlichen, z -> -z-symmetrisch: w^2 = z^3 - z, j = {j33}. Auch das ist tau = i.")
check("B: j(z^3 - z) = 1728", j33 == 1728)
say("   Vierteldrehung (z, w) -> (iz, w) auf w^2 = z^4 - 1: dz/w -> i dz/w; Quadrat -1 auf der 1-Form.")
say("   Die Bedingung J^2 = -1 auf der Form ist für k = 4 automatisch und für k = 1 leer (keine Form).")
say("   Ergebnis: Als Normalendatum wählt die Torusbedingung k = 3 oder 4, nicht 1. Mit mu4 ist k = 4 die")
say("   Hopf-Hebung der vier Marken aus v957 (vier verdrillte Fasern). Gescheitert für Primitivität.")
BEFUNDE["K1: Torusbedingung am Normalendatum wählt k in {3,4}, nicht 1"] = True

# ==========================================================================================================
say("\nC. K2: Spin und Statistik\n")
gam = [np.array(x, complex) for x in (
    [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]],
    [[0, 0, 0, -1j], [0, 0, -1j, 0], [0, 1j, 0, 0], [1j, 0, 0, 0]],
    [[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]],
    [[0, 0, -1j, 0], [0, 0, 0, 1j], [1j, 0, 0, 0], [0, -1j, 0, 0]])]
g5 = gam[0] @ gam[1] @ gam[2] @ gam[3]
I4 = np.eye(4)
cliff = all(np.allclose(gam[a] @ gam[b] + gam[b] @ gam[a], 2 * (a == b) * I4) for a in range(4) for b in range(4))


def spin_rot(th01, th23):
    return expm(th01 / 2 * gam[0] @ gam[1]) @ expm(th23 / 2 * gam[2] @ gam[3])


S_hopf = spin_rot(2 * np.pi, 2 * np.pi)
S_one = spin_rot(2 * np.pi, 0)
say(f"   Clifford-Relationen: {cliff}. Hopf-Drehung um 2 pi (zwei Ebenen) auf Spinoren: "
    f"{'+1' if np.allclose(S_hopf, I4) else '?'}; Drehung in einer Ebene: {'-1' if np.allclose(S_one, -I4) else '?'}.")
check("C: Hopf-Drehung 2 pi hebt zu +1, Einzelebene zu -1 = (-1)^F", cliff and np.allclose(S_hopf, I4) and np.allclose(S_one, -I4))
S_D = spin_rot(np.pi / 2, -np.pi / 2)
P_plus, P_minus = (I4 + g5) / 2, (I4 - g5) / 2
blocks = []
for Pp in (P_plus, P_minus):
    U_, s_, _ = np.linalg.svd(Pp)
    B_ = U_[:, :2]
    blocks.append(B_.conj().T @ S_D @ B_)
evs = [np.round(np.linalg.eigvals(b), 6) for b in blocks]
dets = [np.linalg.det(b) for b in blocks]
say(f"   Deck D gehoben: Eigenwerte je Chiralität {[sorted(e.tolist(), key=lambda x: (x.real, x.imag)) for e in evs]},")
say(f"   Determinanten {np.round(dets, 6).tolist()}. Der Lorentz-skalare Bilinear eps_ab eta^a eta^b trägt det = 1.")
check("C: geometrisches Deck wirkt auf skalare Bilineare trivial", all(abs(d - 1) < 1e-9 for d in dets))
say("   Damit: Jede geometrische Drehung wirkt auf Phi (Boson, J^2-neutral) trivial. Sie gibt q keine")
say("   ungerade Ladung. Ungerade Ladung heißt J-Holonomie längs der Faser, also eine verschlungene Naht.")
mono = []
for k in range(1, 5):
    cs = rng.normal(size=k) + 1j * rng.normal(size=k)
    base = np.array([0.3 + 0.8j, 1.0]) / np.linalg.norm([0.3 + 0.8j, 1.0])
    pts = np.exp(1j * tt)[:, None] * base[None, :]
    qv = np.prod([pts[:, 0] - c * pts[:, 1] for c in cs], axis=0)
    ph = np.unwrap(np.angle(qv))
    sq_end = np.exp(0.5j * (ph[-1] - ph[0]))
    mono.append(int(np.round(sq_end.real)))
say(f"   Monodromie von sqrt(q) längs einer generischen Hopf-Faser für k = 1..4: {mono} = (-1)^k.")
check("C: sqrt(q) spinoriell längs der Faser genau für ungerades k", mono == [(-1) ** k for k in range(1, 5)])
say("   'Phi spinoriell' ist also gleichwertig zu 'ungerade viele verschlungene Nähte'. Das ist zirkulär.")
say("   Spin^c-Buchhaltung unter G (Fermionen Ladung 1/4 relativ zu q): Nähte durch einen Zyklus C sind")
tab = {cc: (2 * cc) % 4 for cc in range(4)}
say(f"   kongruent 2 (C.C) mod 4: {tab}. Nie ungerade. Nur ohne J (q = Phi, Ladung 1/2) wäre k = C.C mod 2.")
check("C: Spin^c unter G erlaubt nie eine ungerade Nahtzahl", all(v_ % 2 == 0 for v_ in tab.values()))
BEFUNDE["K2: Spin-Statistik gibt keine ungerade Ladung, Spinorialität ist zirkulär"] = True

# ==========================================================================================================
say("\nD. K3: Energie und Entropie der Hopf-Sektoren\n")
say("   Ladung-k-Felder auf S^3 (Radius 1) sind Schnitte von O(k) auf der Basis S^2 (Radius 1/2) mit")
say("   Monopolzusammenhang; dazu kommt die Faserenergie k^2 |q|^2. Gitter: Würfeloberfläche, Peierls-Phasen.")
say("   E = 2 pi [ sum |q_j - U q_i|^2 + a^2 sum (k^2 |q|^2 + lam/2 (|q|^2 - v^2)^2) ], lam = 1.")


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
    plaq = []
    E3 = np.eye(3, dtype=int)
    for kk in range(3):
        u, ww = (kk + 1) % 3, (kk + 2) % 3
        for side in (0, N):
            uu, w2 = (u, ww) if side == N else (ww, u)
            for a in range(N):
                for b in range(N):
                    c0 = np.zeros(3, int)
                    c0[kk] = side
                    c0[uu], c0[w2] = a, b
                    corners = [c0, c0 + E3[uu], c0 + E3[uu] + E3[w2], c0 + E3[w2]]
                    ls, ss = [], []
                    for s_ in range(4):
                        pa, pb = index[tuple(corners[s_])], index[tuple(corners[(s_ + 1) % 4])]
                        if (pa, pb) in lindex:
                            ls.append(lindex[(pa, pb)]); ss.append(1)
                        else:
                            ls.append(lindex[(pb, pa)]); ss.append(-1)
                    plaq.append((ls, ss))
    return len(pts), np.array(links), plaq


from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsqr


class Lattice:
    def __init__(self, N):
        self.N = N
        self.n, self.LK, self.PQ = cube(N)
        rows, cols, vals = [], [], []
        for p, (ls, ss) in enumerate(self.PQ):
            for l_, s_ in zip(ls, ss):
                rows.append(p); cols.append(l_); vals.append(s_)
        self.B = csr_matrix((vals, (rows, cols)), shape=(len(self.PQ), len(self.LK)))
        self.a2 = np.pi / self.n

    def background(self, k):
        P = len(self.PQ)
        F = np.full(P, 2 * np.pi * k / P)
        rhs = F.copy()
        rhs[0] -= 2 * np.pi * k
        A = lsqr(self.B, rhs, atol=1e-14, btol=1e-14, iter_lim=200000)[0]
        assert np.allclose(self.B @ A, rhs, atol=1e-8)
        return np.exp(1j * A), F

    def energy(self, x, U, k, lam_, v):
        n, LK, a2 = self.n, self.LK, self.a2
        q = x[:n] + 1j * x[n:]
        i, j = LK[:, 0], LK[:, 1]
        diff = q[j] - U * q[i]
        r = np.abs(q) ** 2
        E = np.sum(np.abs(diff) ** 2) + a2 * np.sum(k * k * r + 0.5 * lam_ * (r - v * v) ** 2)
        g = a2 * (k * k * q + lam_ * (r - v * v) * q)
        g = g + np.bincount(j, weights=diff.real, minlength=n) + 1j * np.bincount(j, weights=diff.imag, minlength=n)
        b = -np.conj(U) * diff
        g = g + np.bincount(i, weights=b.real, minlength=n) + 1j * np.bincount(i, weights=b.imag, minlength=n)
        return 2 * np.pi * E, 4 * np.pi * np.concatenate([g.real, g.imag])

    def windings(self, q, U, F):
        i, j = self.LK[:, 0], self.LK[:, 1]
        w_ = np.angle(q[j] * np.conj(U * q[i]))
        out = np.array([(sum(s_ * w_[l_] for l_, s_ in zip(ls, ss)) + F[p]) / (2 * np.pi)
                        for p, (ls, ss) in enumerate(self.PQ)])
        return np.round(out).astype(int)

    def sectors(self, v, trials):
        Es, ok = [], True
        for k in range(4):
            U, F = self.background(k)
            best = None
            for _ in range(trials):
                res = minimize(self.energy, rng.normal(size=2 * self.n) * v, args=(U, k, 1.0, v), jac=True,
                               method="L-BFGS-B", options={"maxiter": 30000, "gtol": 1e-9, "ftol": 1e-15})
                if best is None or res.fun < best.fun:
                    best = res
            q = best.x[:self.n] + 1j * best.x[self.n:]
            wn = self.windings(q, U, F)
            ok &= sorted(wn[wn != 0].tolist()) == [1] * k
            Es.append(float(best.fun))
        return Es, ok


L12 = Lattice(12)
ev_ok = []
for k in range(4):
    U, F = L12.background(k)
    H = np.zeros((L12.n, L12.n), complex)
    for l_, (a, b) in enumerate(L12.LK):
        H[a, a] += 1; H[b, b] += 1; H[b, a] -= U[l_]; H[a, b] -= np.conj(U[l_])
    lo = np.linalg.eigvalsh(H)[0] / L12.a2 + k * k
    ev_ok.append(abs(lo - k * (k + 2)) < 0.05 * max(1, k * (k + 2)))
    say(f"   k = {k}: tiefster Eigenwert des Gitteroperators {lo:7.3f}, Kontinuum k(k+2) = {k * (k + 2)}")
check("D: Gitter trifft die Ladung-k-Spektren k(k+2) auf 5 %", all(ev_ok))

EK, seams_ok = {}, True
L24 = Lattice(24)
for v, lat, tr in ((4.0, L12, 3), (8.0, L12, 3), (12.0, L12, 3), (16.0, L24, 2)):
    Es, ok = lat.sectors(v, tr)
    seams_ok &= ok
    EK[v] = Es
    d1, d2 = Es[1] - Es[0], Es[2] - Es[1]
    say(f"   v R = {v:4.1f} (Würfel N = {lat.N}): E(k) = {[round(e_, 1) for e_ in Es]};"
        f" E(1)-E(0) = {d1:.1f}, E(2)-E(1) = {d2:.1f}")
DATA["E_k_hopf"] = {str(v): e_ for v, e_ in EK.items()}
check("D: im Sektor k genau k einfache Nähte (Windung +1 je Naht)", seams_ok)
check("D: E(k) steigt streng mit k, also k = 1 minimal unter k != 0",
      all(all(np.diff(e_) > 0) for e_ in EK.values()))
widths = {v: (EK[v][2] - 2 * EK[v][1] + EK[v][0]) / (EK[v][1] - EK[v][0]) for v in EK}
say("   Mit konstanter Entropie mu je Naht wählt F(k) = beta E(k) - k mu genau k = 1 nur im Fenster")
say("   E(1)-E(0) < mu/beta < E(2)-E(1). Relative Fensterbreite (E(2)-2E(1)+E(0))/(E(1)-E(0)):")
say(f"   {', '.join(f'v R = {v:.0f}: {w_:+.4f}' for v, w_ in widths.items())}.")
say("   Negativ heißt leer (E konkav). Wo es existiert, ist es unter 1 % breit: Feinabstimmung von mu.")
check("D: Fenster für genau eine Naht leer oder relativ unter 1 %",
      all(w_ < 0.01 for w_ in widths.values()) and any(w_ <= 0 for w_ in widths.values()))
BEFUNDE["K3: Fenster für genau eine Naht bei konstanter Nahtentropie leer oder unter 1 % breit"] = True
DATA["entropy_window_width"] = {str(v): float(w_) for v, w_ in widths.items()}
ratio = {v: 2 * np.pi * EK[v][1] / np.log(12) for v in EK}
say(f"   Größenordnung bei beta = 2 pi (tau = i): beta E(1) / log 12 = {[round(r_) for r_ in ratio.values()]}.")
say("   Die Nahtentropie log chi_E8(i) = log 12 ist gegen die Volumenkosten verschwindend; k = 0 gewinnt.")
BEFUNDE["K3: k = 0 bleibt Grundzustand, log 12 je Naht reicht nicht"] = True
check("D: k = 0 hat die kleinste freie Energie auch mit log 12 je Naht",
      all(2 * np.pi * EK[v][k] - k * np.log(12) > 0 for v in EK for k in (1, 2, 3)))

# ==========================================================================================================
say("\nE. K3': Nullmittel. Ordnungsfeld ohne Hintergrundwert, tiefste zulässige Mode l = 1\n")
say("   Lesart: A* auf das Ordnungsfeld übertragen, also kein konstanter Anteil, int q dV = 0. Dann ist die")
say("   tiefste Mode l = 1 (Eigenwert 3), und jede komplexe l = 1-Mode ist q = a.x mit a in C^4.")
circ_ok, wind_vals = True, []
for _ in range(200):
    a = rng.normal(size=4) + 1j * rng.normal(size=4)
    Mreal = np.vstack([a.real, a.imag])
    ns = np.linalg.svd(Mreal)[2][2:]                   # orthonormal basis of the null plane
    tc = np.linspace(0, 2 * np.pi, 200)[:, None]
    pts = np.cos(tc) * ns[0] + np.sin(tc) * ns[1]
    if np.max(np.abs(pts @ a)) > 1e-12:
        circ_ok = False
    nrm_basis = np.linalg.svd(Mreal)[2][:2]
    p0 = ns[0]
    eps = 1e-3
    loop = p0[None, :] * np.sqrt(1 - eps ** 2) + eps * (np.cos(tt)[:, None] * nrm_basis[0] + np.sin(tt)[:, None] * nrm_basis[1])
    wv = np.sum(np.diff(np.unwrap(np.angle(loop @ a)))) / (2 * np.pi)
    wind_vals.append(int(np.round(wv)))
say(f"   200 zufällige a: Nullmenge = Großkreis (Kern von Re a, Im a): {circ_ok}; Windungen {sorted(set(wind_vals))}.")
check("E: jede generische l = 1-Mode hat genau eine Naht mit Windung +-1", circ_ok and set(wind_vals) <= {1, -1})
X = rng.normal(size=(400000, 4))
X /= np.linalg.norm(X, axis=1)[:, None]
m4, m22 = np.mean(X[:, 0] ** 4), np.mean(X[:, 0] ** 2 * X[:, 1] ** 2)
say(f"   Momente auf S^3 (Monte Carlo): <x^4> = {m4:.4f} (1/8), <x1^2 x2^2> = {m22:.4f} (1/24).")
check("E: Momente 1/8 und 1/24", abs(m4 - 1 / 8) < 2e-3 and abs(m22 - 1 / 24) < 1e-3)


def quartic_ratio(a):
    Pm = np.outer(a.real, a.real) + np.outer(a.imag, a.imag)
    p = np.sort(np.linalg.eigvalsh(Pm))[-2:]
    return (2 * (p[0] ** 2 + p[1] ** 2) + 4 / 3 * p[0] * p[1]) / (p[0] + p[1]) ** 2


rs = [quartic_ratio(rng.normal(size=4) + 1j * rng.normal(size=4)) for _ in range(2000)]
r_circ = quartic_ratio(np.array([1, 1j, 0, 0]))
r_real = quartic_ratio(np.array([1, 0, 0, 0], complex))
mc_circ = np.mean(np.abs(X[:, 0] + 1j * X[:, 1]) ** 4) / np.mean(np.abs(X[:, 0] + 1j * X[:, 1]) ** 2) ** 2
say(f"   <|q|^4>/<|q|^2>^2: kreisförmig (x0 + i x1) {r_circ:.4f} (MC {mc_circ:.4f}), reell (Wand) {r_real:.4f},")
say(f"   Minimum über 2000 zufällige a: {min(rs):.4f}. Nahe der Schwelle lam v^2 = 3 wählt das Quartische die")
say("   kreisförmige Mode: genau ein Großkreis, Windung 1, Lage frei (Modulraum der orientierten Ebenen).")
check("E: Quartik minimal genau bei der kreisförmigen Mode (eine Naht), Wand teurer",
      abs(r_circ - 4 / 3) < 1e-12 and r_real == 2 and min(rs) >= 4 / 3 - 1e-12 and abs(mc_circ - 4 / 3) < 0.01)

say("\n   Tief geordnet: Naht-Ansatz q = f(eta) e^(i xi1) gegen nullstellenfreie Textur q = g(chi) mit int q = 0.")
lam = 1.0


def seam_energy(v, M1=300, M2=400):
    xi = 1 / (v * np.sqrt(lam))
    rho = np.unique(np.concatenate([np.linspace(0, min(12 * xi, 0.5), M1), np.geomspace(min(12 * xi, 0.5), np.pi / 2, M2)]))
    h = np.diff(rho)
    mid = 0.5 * (rho[1:] + rho[:-1])
    wmid = np.sin(mid) * np.cos(mid)
    dual = np.zeros_like(rho)
    dual[:-1] += h / 2
    dual[1:] += h / 2
    wnod = np.sin(rho) * np.cos(rho)
    cent = np.zeros_like(rho)
    cent[1:] = np.cos(rho[1:]) / np.sin(rho[1:])

    def E(f_in):
        f = np.concatenate([[0.0], f_in])
        df = np.diff(f)
        e = np.sum(df ** 2 / h * wmid) + np.sum(dual * (cent * f ** 2 + wnod * 0.5 * lam * (f ** 2 - v * v) ** 2))
        g = np.zeros_like(f)
        g[1:] += 2 * df / h * wmid
        g[:-1] -= 2 * df / h * wmid
        g += dual * (2 * cent * f + wnod * 2 * lam * (f ** 2 - v * v) * f)
        return 4 * np.pi ** 2 * e, 4 * np.pi ** 2 * g[1:]
    r = minimize(E, v * np.tanh(rho[1:] / xi), jac=True, method="L-BFGS-B",
                 options={"maxiter": 50000, "gtol": 1e-9, "ftol": 1e-14})
    return r.fun


def texture_energy(v, M=800):
    chi = np.linspace(0, np.pi, M)
    h = chi[1] - chi[0]
    s_ = np.sin(chi) ** 2
    wmid = np.sin(0.5 * (chi[1:] + chi[:-1])) ** 2
    dual = np.full(M, h)
    dual[0] = dual[-1] = h / 2
    cw = dual * s_
    nrm = cw @ s_
    proj = lambda hh: hh - s_ * (cw @ hh) / nrm

    def E(x):
        g = proj(x[:M] + 1j * x[M:])
        dg = np.diff(g)
        r = np.abs(g) ** 2
        e = 4 * np.pi * (np.sum(np.abs(dg) ** 2 / h * wmid) + np.sum(cw * 0.5 * lam * (r - v * v) ** 2))
        gr = np.zeros(M, complex)
        gr[1:] += 2 * dg / h * wmid
        gr[:-1] -= 2 * dg / h * wmid
        gr = 4 * np.pi * (gr + cw * 2 * lam * (r - v * v) * g)
        gr = gr - cw * (s_ @ gr) / nrm
        return e, np.concatenate([gr.real, gr.imag])
    g0 = v * np.exp(4j * (chi - np.pi / 2))
    r = minimize(E, np.concatenate([g0.real, g0.imag]), jac=True, method="L-BFGS-B",
                 options={"maxiter": 50000, "gtol": 1e-9, "ftol": 1e-14})
    g = proj(r.x[:M] + 1j * r.x[M:])
    return r.fun, abs(cw @ g) / (v * nrm), np.abs(g).min() / v


VS = [2.0, 4.0, 8.0, 16.0, 32.0, 64.0]
EA, EB = [], []
for v in VS:
    ea = seam_energy(v)
    eb, res_, mn = texture_energy(v)
    EA.append(ea / v ** 2)
    EB.append(eb / v ** 2)
    say(f"   v R = {v:5.1f}: Naht E/v^2 = {ea / v ** 2:7.2f}; Textur E/v^2 = {eb / v ** 2:7.2f} (Mittel {res_:.0e}, min|g|/v = {mn:.3f})")
slope = (EA[-1] - EA[-3]) / np.log(VS[-1] / VS[-3])
vc = VS[-1] * np.exp((EB[-1] - EA[-1]) / slope)
say(f"   Naht wächst wie {slope:.1f} log(vR) (4 pi^2 = {4 * np.pi ** 2:.1f}), die Textur sättigt bei {EB[-1]:.0f}.")
say(f"   Übergang bei v R = {vc:.0f}, also R/xi ~ {vc:.0f} (lam = 1). Darüber ist eine nullstellenfreie Textur billiger.")
check("E: Naht unterhalb, Textur oberhalb des Übergangs (Ansatzschranken)",
      all(a_ < b_ for a_, b_ in zip(EA, EB)) and abs(slope - 4 * np.pi ** 2) < 1.0 and 64 < vc < 1000)
DATA["zero_mean"] = {"v": VS, "seam_E_over_v2": EA, "texture_E_over_v2": EB, "crossover_vR": float(vc)}
BEFUNDE["K3': Nullmittel erzwingt die Naht nur für R/xi unter etwa %d" % round(vc)] = True

# ==========================================================================================================
say("\nF. K4: Kopiedeterminante\n")
U1 = sp.Matrix([[1, z, 0], [0, 1, z ** 2], [0, 0, 1]])
U2 = sp.Matrix([[1, 0, 0], [2 * z, 1, 0], [1 - z, z, 1]])
V1 = sp.Matrix([[1, 0, 3], [0, 1, 0], [0, 1 - z, 1]])
C = sp.expand(U1 * U2 * sp.diag(1, 1, z) * V1)
detC = sp.factor(C.det())
rk0 = C.subs(z, 0).rank()
say(f"   C(z) = U diag(1,1,z) V mit unimodularen U, V: det C = {detC}, Rang bei z = 0: {rk0}, Kern {3 - rk0}.")
check("F: primitive Determinante gibt genau einen Rangverlust an genau einer Stelle", detC in (z, -z) and rk0 == 2)
ok_ord, eq_cnt = True, 0
for r_ in (1, 2, 3):
    for _ in range(10):
        if r_ < 3:
            A0 = sp.Matrix(rng.integers(-3, 4, size=(3, 3 - r_)).tolist()) * sp.Matrix(rng.integers(-3, 4, size=(3 - r_, 3)).tolist())
        else:
            A0 = sp.zeros(3, 3)
        A1 = sp.Matrix(rng.integers(-3, 4, size=(3, 3)).tolist())
        A2 = sp.Matrix(rng.integers(-3, 4, size=(3, 3)).tolist())
        d_ = sp.Poly(sp.expand((A0 + A1 * z + A2 * z ** 2).det()), z)
        if d_.is_zero:
            continue
        kerdim = 3 - A0.rank()
        ordv = min(m[0] for m in d_.monoms())
        ok_ord &= ordv >= kerdim
        eq_cnt += ordv == kerdim
say(f"   Ungleichung dim ker C(z0) <= ord_z0 det C an 30 Zufallsfällen: {ok_ord} (Gleichheit in {eq_cnt}).")
check("F: dim ker <= Nullstellenordnung der Determinante", ok_ord)
say("   Drei gleiche Kopien C = z 1_3: det = z^3, Kern 3 bei z = 0, alle drei winden (c = 24).")

say("\n   Ohne J (Konstanten erlaubt): C(z) = Y z + m, Y und m komplex symmetrisch.")
cor_ok, sep_min = True, np.inf
for _ in range(500):
    Y = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); Y = Y + Y.T
    m = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); m = m + m.T
    zs = np.linalg.eigvals(-np.linalg.solve(Y, m))
    for z0 in zs:
        sv = np.linalg.svd(Y * z0 + m, compute_uv=False)
        if not (sv[-1] < 1e-8 * sv[0] and sv[-2] > 1e-6 * sv[0]):
            cor_ok = False
    sep_min = min(sep_min, min(abs(zs[i] - zs[j]) for i in range(3) for j in range(i + 1, 3)))
mm = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); mm = mm + mm.T
spl = []
for e_ in (1e-2, 1e-4):
    zs = np.linalg.eigvals(-e_ * mm)
    spl.append(min(abs(zs[i] - zs[j]) for i in range(3) for j in range(i + 1, 3)) / e_)
say(f"   500 Zufallsfälle: drei getrennte einfache Nullstellen, jede mit Korang 1: {cor_ok} (kleinster Abstand {sep_min:.1e}).")
say(f"   Störung von z 1_3 durch eps m: Aufspaltung / eps = {spl[0]:.4f} und {spl[1]:.4f} (linear in eps).")
check("F: ohne J generisch Korang 1 an getrennten Stellen; Dreifachwindung spaltet linear auf",
      cor_ok and abs(spl[0] - spl[1]) < 1e-6 * max(spl))


def var_dim(fun, nparam, sym=False):
    p0 = rng.normal(size=nparam) + 1j * rng.normal(size=nparam)
    base = fun(p0).ravel()
    J = []
    for i in range(nparam):
        dp = np.zeros(nparam, complex); dp[i] = 1e-7
        J.append((fun(p0 + dp).ravel() - base) / 1e-7)
    J = np.array(J).T
    if sym:
        iu = np.triu_indices(3)
        J = np.array([fun_col.reshape(3, 3)[iu] for fun_col in J.T]).T
    return np.linalg.matrix_rank(J, tol=1e-4)


d_gen1 = var_dim(lambda p: p[:6].reshape(3, 2) @ p[6:].reshape(2, 3), 12)
d_gen2 = var_dim(lambda p: np.outer(p[:3], p[3:]), 6)
d_sym1 = var_dim(lambda p: p[:6].reshape(3, 2) @ np.diag(p[6:8]) @ p[:6].reshape(3, 2).T, 8, sym=True)
d_sym2 = var_dim(lambda p: np.outer(p, p), 3, sym=True)
say(f"   Dimension der Korang-r-Mengen (komplex): allgemein r=1: {d_gen1} von 9, r=2: {d_gen2};"
    f" symmetrisch r=1: {d_sym1} von 6, r=2: {d_sym2}.")
say(f"   Reelle Kodimension von Korang 2: allgemein {2 * (9 - d_gen2)}, symmetrisch {2 * (6 - d_sym2)}; beides > 4.")
check("F: Kodimensionen 1, 4 (allgemein) und 1, 3 (symmetrisch)",
      (d_gen1, d_gen2, d_sym1, d_sym2) == (8, 5, 5, 3))
say("   Ohne Symmetrie trägt also in 4D jede Naht generisch genau einen windenden Block.")

say("\n   Mit J (Prämisse G): J wirkt auf eta+ mit -i, auf eta- mit +i.")
Jc = {"eta+ eta+": (-1j) ** 2, "eta+ eta-": (-1j) * (1j), "eta- eta-": (1j) ** 2}
say(f"   J-Ladungen der Bilineare: {', '.join(f'{k_}: {v_.real:+.0f}' for k_, v_ in Jc.items())}.")
say("   Nackte Majoranamassen im windenden Kanal sind verboten; dort ist C = Y Phi ohne Konstante.")
check("F: J verbietet nackte Masse im windenden Kanal, erlaubt die Diracmasse eta+ eta-",
      abs(Jc["eta+ eta+"] + 1) < 1e-12 and abs(Jc["eta+ eta-"] - 1) < 1e-12)
ranks = []
for _ in range(500):
    Y = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); Y = Y + Y.T
    ranks.append(int(np.linalg.matrix_rank(Y)))
say(f"   det(Y Phi) = det Y Phi^3: alle gekoppelten Kopien winden an derselben Stelle. Rang generischer Y: {set(ranks)}.")
say(f"   Rang 1 (genau ein windender Block) hat komplexe Kodimension {6 - d_sym2} in Sym^2 C^3: Feinabstimmung.")
check("F: unter J ist der generische Rang 3, Rang 1 kodimensional", set(ranks) == {3} and 6 - d_sym2 == 3)
BEFUNDE["K4: unter J winden alle an Phi gekoppelten Kopien gemeinsam; genau einer verlangt Rang Y = 1"] = True

say("\n   Kopien mit Familiencharakteren (-1, i, -i) unter der Vierteldrehung, Phi mit Charakter chi_Phi:")
fam = [-1, 1j, -1j]
rk_by = {}
for cname, cphi in (('1', 1), ('i', 1j), ('-1', -1), ('-i', -1j)):
    allowed = [(a, b) for a in range(3) for b in range(a, 3) if abs(cphi * fam[a] * fam[b] - 1) < 1e-12]
    Yf = np.zeros((3, 3), complex)
    for (a, b) in allowed:
        val = rng.normal() + 1j * rng.normal()
        Yf[a, b] = Yf[b, a] = val
    rk_by[cname] = (int(np.linalg.matrix_rank(Yf)), allowed)
    say(f"   chi_Phi = {cname}: erlaubte Einträge {[(a + 1, b + 1) for a, b in allowed]}, generischer Rang {rk_by[cname][0]}")
check("F: mit Familiencharakteren ist der Rang nie 1", all(r_[0] >= 2 for r_ in rk_by.values()))
say("   Tragen die schweren Kopien die Familiencharaktere, winden stets mindestens zwei Blöcke zusammen.")
say("   Genau ein windender Block verlangt also eine schwere Kopie ohne Familienladung (TFPTs ein E8).")
BEFUNDE["K4: familiengeladene schwere Kopien können nicht genau einen Block winden"] = True

# ==========================================================================================================
say("\nG. K6: Gemeinsame GSO-Summe (A*) für n windende Blöcke\n")
say("   Unter A* ist die GSO-Summe die Summe über die Spinstrukturen des einen Nahttorus. Sie ist allen")
say("   windenden Fermionen gemeinsam. N = 16n Majorana-Weyl geben Z = (th3^8n + th4^8n + th2^8n)/(2 eta^8n).")
PMAX = 24                                                   # powers of s = q^(1/8)


def smul(a, b):
    c = [0] * (PMAX + 1)
    for i, ai in enumerate(a):
        if ai:
            for j in range(PMAX + 1 - i):
                if b[j]:
                    c[i + j] += ai * b[j]
    return c


def spow(a, e):
    r = [1] + [0] * PMAX
    for _ in range(e):
        r = smul(r, a)
    return r


th3, th4, th2 = [0] * (PMAX + 1), [0] * (PMAX + 1), [0] * (PMAX + 1)
for n_ in range(-4, 5):
    if 4 * n_ * n_ <= PMAX:
        th3[4 * n_ * n_] += 1
        th4[4 * n_ * n_] += (-1) ** n_
    if (2 * n_ + 1) ** 2 <= PMAX:
        th2[(2 * n_ + 1) ** 2] += 1
etainv = [1] + [0] * PMAX
for n_ in range(1, PMAX // 8 + 1):
    geo = [1 if m_ % (8 * n_) == 0 else 0 for m_ in range(PMAX + 1)]
    etainv = smul(etainv, geo)
cur, hR, Zi = {}, {}, {}
mp.mp.dps = 30
nome = mp.e ** (-mp.pi)
T2, T3, T4 = (mp.jtheta(k_, 0, nome) for k_ in (2, 3, 4))
ETA = mp.gamma(0.25) / (2 * mp.pi ** 0.75)
gso_int = True
for nb in (1, 2, 3):
    N = 16 * nb
    t2p = spow(th2, N // 2)
    Z_ = [(x + y + w_) // 2 for x, y, w_ in zip(spow(th3, N // 2), spow(th4, N // 2), t2p)]
    Z_ = smul(Z_, spow(etainv, N // 2))
    gso_int &= all(Z_[i] == 0 for i in range(PMAX + 1) if i % 8)
    cur[nb] = int(Z_[8])
    hR[nb] = min(i for i in range(PMAX + 1) if t2p[i]) / 8
    Zi[nb] = float((T3 ** (N // 2) + T4 ** (N // 2) + T2 ** (N // 2)) / (2 * ETA ** (N // 2)))
    say(f"   n = {nb}: Ströme vom Gewicht 1: {cur[nb]} (so({N}): {comb(N, 2)}, E8^{nb}: {248 * nb}); h_R = {hR[nb]:.0f};"
        f" Z(i) = {Zi[nb]:.4f}")
check("G: GSO-Summe hat nur ganze q-Potenzen für N = 16, 32, 48", gso_int)
check("G: Ströme 248, 496, 1128; Spinorströme nur für n = 1",
      cur == {1: 248, 2: 496, 3: 1128} and all(cur[n_] - comb(16 * n_, 2) == 128 * (n_ == 1) for n_ in cur))
check("G: Twistgewicht h_R = n", hR == {1: 1.0, 2: 2.0, 3: 3.0})
check("G: Z(i) = 12, 144, 2112 (= j(i) + 384)", abs(Zi[1] - 12) < 1e-9 and abs(Zi[2] - 144) < 1e-8 and abs(Zi[3] - 2112) < 1e-6)
roots24 = cur[3] - 24
check("G: D24+ hat 1104 = 2*24*23 Wurzeln, keine aus der Spinorklasse", roots24 == 2 * 24 * 23)
say(f"   n = 3 mit gemeinsamer GSO ist D24+ (Niemeier, Wurzelsystem D24, {roots24} Wurzeln), nicht E8^3.")
say("   Unterwurzelsysteme von D24 sind vom Typ A und D (Borel-de Siebenthal, zitiert): kein E8.")
say("   E8^3 verlangt drei unabhängige GSO-Projektionen, also zwei zusätzliche geeichte Z2. A* liefert nur eine.")
say("   Die sichtbare Sonde SO(2) hat in so(48)_1 Niveau 1, c = 24: dieselben Zahlen wie Z14s verdeckter Zweig.")
say("   Alpha (EM1) schließt n = 3 so nicht aus; der Ausschluss kommt allein aus der Forderung E8 an der Naht.")
say(f"   Entropie: log Z(i) = {np.log(Zi[1]):.3f}, {np.log(Zi[2]):.3f}, {np.log(Zi[3]):.3f}. Sie bevorzugt mehr Blöcke.")
BEFUNDE["K6: Nahtentropie bevorzugt mehr windende Blöcke (2112 > 12)"] = True
BEFUNDE["K6: n = 3 mit gemeinsamer GSO ist D24+, besteht alpha (Niveau 1), enthält aber kein E8"] = True
DATA["gso"] = {"currents": cur, "h_R": hR, "Z_i": Zi}

# ==========================================================================================================
say("\nH. Urteil\n")
URTEIL = {
    "K5_quellen": ("Lesart", "TFPT setzt n = 1 als primitives Postulat; N_fam = 3n verknüpft n mit der Familienzahl. "
                   "Mit C'' je Naht ist N_fam = 3k, also wählt N_nu = 3 (externe Zahl) k = 1 empirisch, nicht als "
                   "Grundzustand. In TFPTs eigener Hopf-Hebung (v957) sind die vier Marken vier Fasern; wirkt das Deck "
                   "auf den Raum, hat das J-neutrale q gerades k, und eine einzelne Naht ist unmöglich."),
    "K1_verzweigung": ("gescheitert", "Die Torusbedingung am Normalendatum w^2 = q wählt k = 3 oder 4, nicht 1; "
                       "J^2 = -1 auf der 1-Form ist bei k = 4 automatisch und bei k = 1 leer."),
    "K2_spin_statistik": ("gescheitert", "Geometrische Drehungen heben zu (-1)^F oder 1 und wirken auf Phi trivial. "
                          "Spinorielles Phi ist gleichwertig zu ungerader Verschlingung (zirkulär). Spin^c unter G gibt "
                          "nur gerade Nahtzahlen."),
    "K3_energie_entropie": ("gescheitert", "E(k) steigt streng: unter k != 0 ist k = 1 minimal (hergeleitet im Modell). "
                            "k = 0 bleibt aber der Grundzustand; log 12 je Naht ist um Größenordnungen zu klein. Das Fenster "
                            "für genau eine Naht ist für v R <= 12 leer (E konkav) und bei v R = 16 nur 0,1 % breit."),
    "K3s_nullmittel": ("Lesart", "Ohne Hintergrundwert des Ordnungsfelds (int q = 0) ist die tiefste Mode l = 1; jede "
                       "generische l = 1-Mode hat genau eine Naht mit Windung 1, und das Quartische wählt sie nahe der "
                       "Schwelle (hergeleitet). Gilt nur für R/xi unter etwa %d; darüber ist eine nullstellenfreie Textur "
                       "billiger. Das Nullmittel selbst ist eine übertragene Lesart von A*, keine Folgerung." % round(vc)),
    "K4_kopiedeterminante": ("Lesart", "Primitive Determinante gibt genau einen Rangverlust (hergeleitet, Smith). Ohne J "
                             "trägt jede Naht generisch genau einen Block. Unter G verbietet J aber Konstanten im windenden "
                             "Kanal; dann winden alle gekoppelten Kopien gemeinsam, und genau einer verlangt Rang Y = 1. "
                             "Mit Familiencharakteren ist Rang 1 unmöglich."),
    "K6_gemeinsame_GSO": ("Lesart", "Unter A* bilden n windende Blöcke D_{8n}^+; E8 nur für n = 1. Genau ein Block ist "
                          "damit gleichwertig zur Forderung E8 an der Naht (H6 wählte N = 16 als kleinstes N). Keine "
                          "Grundzustandsauswahl: Die Entropie bevorzugt mehr Blöcke."),
    "gesamt": ("offen", "Keine Regel wählt den primitiven Sektor und genau einen windenden Block als Grundzustand. "
               "Erreicht ist eine Reduktion: (a) k = 1 ist minimal unter k != 0 und folgt aus einem Nullmittel des "
               "Ordnungsfelds nahe der Schwelle, oder empirisch aus N_fam = 3k; im Linsenraum X ist es ausgeschlossen. "
               "(b) Genau ein Block ist unter A* gleichwertig zu E8 und verlangt eine schwere Kopie ohne Familienladung."),
}
for k_, (v_, txt) in URTEIL.items():
    say(f"   {k_}: {v_}.")
    for chunk in textwrap.wrap(txt, 100):
        say(f"      {chunk}")

say("\nZusammenfassung")
for n_, c_ in CHECKS.items():
    say(f"   {'bestanden' if c_ else 'NICHT bestanden'}: {n_}")
say(f"   Befunde (negativ oder begrenzend): {len(BEFUNDE)}")
(HERE / "ausgabe_z17.txt").write_text("\n".join(LINES) + "\n")
(HERE / "z17_ergebnis.json").write_text(json.dumps({
    "checks": CHECKS, "befunde": BEFUNDE,
    "urteil": {k_: {"stufe": v_, "text": txt} for k_, (v_, txt) in URTEIL.items()},
    "data": DATA}, ensure_ascii=False, indent=1, default=str))
assert all(CHECKS.values())
