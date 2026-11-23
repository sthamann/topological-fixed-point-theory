"""Z14: Satzzahl n = 1 gegen n = 3 über das Niveau eines gemeinsamen Stroms.

A. Gitter (Modell aus Z7, Sektortechnik aus Z9). Jeder windende Satz wird auf zwei Komponenten seines
   reellen 16-Pletts reduziert; eine SO(2) dreht die beiden. Ein Probestrom koppelt an
   Q_probe = sum_s w_s Q_s. Gemessen wird im Sektor Q_probe = +1 der Sprung der Kernbesetzung bei k_z = 0
   (gepumpte Ladung je Flussquant), die Zahl der chiralen Kernmoden und ihre Geschwindigkeit. Varianten:
   n = 1, 2, 3 ohne Kopplung, n = 2, 3 mit einer Kopplung zwischen den Sätzen, die nur die gemeinsame
   SO(2) erhält (die einzelnen Satzladungen sind gebrochen), und ein Zweig mit verborgenen Sätzen
   (Probe nur an Satz 1).
B. E8-Einbettung: Dynkin-Index der SO(2) in E8, Weyl-Dimension von V(2 theta), Sugawara und Coset.
C. TFPT-Beziehungen bei Niveau 3 beziehungsweise c = 24: c3 je Cartanrichtung, chi_E8(i) = 12 und j(i),
   Holomorphie-Eindeutigkeit, c = g_car + N_fam, das alpha^3-Glied als U(1)-Niveau, der Seedterm 48 c3^4,
   R102, T-Phase.
D. Urteil.
"""
import itertools
import json
import time
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.linalg import block_diag, eigh

T0 = time.time()
HERE = Path(__file__).resolve().parent
LINES = []
CHECKS = {}
BEFUNDE = {}
DATA = {}


def P(s=""):
    print(s)
    LINES.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    return "ja" if cond else "NEIN"


def befund(name, val, text):
    BEFUNDE[name] = {"wert": val, "text": text}


# --------------------------------------------------------------------------------------------------
P("A. Gitter: Niveau eines gleichförmig koppelnden Probestroms\n")
_src = (HERE / "z7_ursprung_naht.py").read_text()
_ns = {"__file__": str(HERE / "z7_ursprung_naht.py")}
exec(_src[:_src.index("L = 14")], _ns)
build, MASS = _ns["build"], _ns["MASS"]

L, R = 12, 3.5
nsite = L * L
NP = nsite * 4 * 2                     # Teilchenblock eines Satzes (2 Komponenten x (eta+, eta-) x Spin)
qs = np.kron(np.array([[0, -1j], [1j, 0]]), np.eye(2))
q1 = np.kron(np.eye(nsite), np.kron(qs, np.eye(2)))
QB = np.block([[q1, np.zeros_like(q1)], [np.zeros_like(q1), -q1.T]])
wq, UQ = np.linalg.eigh(QB)
UP, UM = UQ[:, wq > 0], UQ[:, wq < 0]
_, (rad, nsp) = build(L, 0.0, 1, 0, 0.25, copies=2)
CORE1 = np.tile(np.repeat(rad <= R, nsp * 2), 2).astype(float)
OUTER1 = np.tile(np.repeat(rad >= L / 2 - 2.0, nsp * 2), 2).astype(float)
rng = np.random.default_rng(20261123)
TSITE = rng.uniform(0.15, 0.35, nsite)
ISY = np.array([[0, 1], [-1, 0]], complex)


def set_H(kz, m):
    MASS.clear()
    MASS["m"] = m
    H, _ = build(L, kz, 1, 0, 0.25, copies=2)
    MASS.clear()
    return H


def mixing(n, pairs, t=1.0, mu=0.15):
    """Kopplung zwischen Sätzen: Hüpfen t(r) psi_s^dag psi_s' und Paarung mu eta-_s eta-_s', beides
    diagonal im Komponentenindex, also nur unter der gemeinsamen SO(2) invariant."""
    dim = 2 * NP
    V = np.zeros((n * dim, n * dim), complex)
    T = np.kron(np.diag(TSITE * t), np.eye(8))
    sel = np.kron(np.eye(2), np.diag([0.0, 1.0]))          # (Komponente, eta+/eta-): nur eta-
    Dm = mu * np.kron(np.eye(nsite), np.kron(sel, ISY))
    for (s, r) in pairs:
        for (a, b) in ((s, r), (r, s)):
            pa, pb = slice(a * dim, a * dim + NP), slice(b * dim, b * dim + NP)
            ha, hb = slice(a * dim + NP, (a + 1) * dim), slice(b * dim + NP, (b + 1) * dim)
            V[pa, pb] += T
            V[ha, hb] += -T.conj()
            V[pa, hb] += Dm
            V[ha, pb] += Dm.conj().T
    return V


def total_H(kz, masses, V):
    return block_diag(*[set_H(kz, m) for m in masses]) + V


def measure(label, masses, pairs, weights, t=1.0, mu=0.15, dk=1e-3):
    n = len(masses)
    V = mixing(n, pairs, t, mu) if pairs else np.zeros((n * 2 * NP, n * 2 * NP), complex)
    Qs = [block_diag(*[QB if j == s else np.zeros_like(QB) for j in range(n)]) for s in range(n)]
    Qp = sum(w * Q for w, Q in zip(weights, Qs))
    H03 = total_H(0.3, masses, V)
    comm_p = float(np.max(np.abs(H03 @ Qp - Qp @ H03)))
    comm_1 = float(np.max(np.abs(H03 @ Qs[0] - Qs[0] @ H03)))
    blocks = [UP if w == 1 else (UM if w == -1 else np.zeros((2 * NP, 0))) for w in weights]
    Up = block_diag(*blocks)
    assert np.allclose(Qp @ Up, Up)
    core, outer = np.tile(CORE1, n), np.tile(OUTER1, n)
    Pc = Up.conj().T @ (core[:, None] * Up)
    Po = Up.conj().T @ (outer[:, None] * Up)
    dH = block_diag(*[(set_H(dk, m) - set_H(-dk, m)) / (2 * dk) for m in masses])
    dHs = Up.conj().T @ dH @ Up
    res, low = {}, {}
    for kz in (-0.1, -0.05, 0.05, 0.1):
        E, W = eigh(Up.conj().T @ total_H(kz, masses, V) @ Up, driver="evr")
        fc = np.real((W.conj() * (Pc @ W)).sum(0))
        fo = np.real((W.conj() * (Po @ W)).sum(0))
        occ = E < 0
        res[kz] = (fc[occ].sum(), fo[occ].sum(), int(occ.sum()))
        if kz in (-0.05, 0.05):
            sel = (np.abs(E) < 0.15) & (fc > 0.5)              # am Kern lokalisierte Eigenzustände nahe E = 0
            vel = np.real((W[:, sel].conj() * (dHs @ W[:, sel])).sum(0))   # Hellmann-Feynman: dE/dk_z
            low[kz] = (E[sel], fc[sel], vel, W[:, sel], int(np.sum(np.abs(E) < 0.15)))
    jc = (2 * res[0.05][0] - res[0.1][0]) - (2 * res[-0.05][0] - res[-0.1][0])
    jo = (2 * res[0.05][1] - res[0.1][1]) - (2 * res[-0.05][1] - res[-0.1][1])
    Ep, cw, vel, Wc, nlow = low[0.05]
    # Spektralfluss: Kernzustände, die zwischen k_z = -0,05 und +0,05 durch E = 0 laufen (Vorzeichen = Richtung)
    ncross = int(np.sum(low[0.05][0] > 0) - np.sum(low[-0.05][0] > 0))
    nz = sum(1 for bl in blocks if bl.shape[1])
    setw = [[float(np.sum(np.abs(Wc[j * NP:(j + 1) * NP, i]) ** 2)) for j in range(nz)] for i in range(Wc.shape[1])]
    out = {"n": n, "weights": list(weights), "masses": list(masses), "pairs": pairs, "t_scale": t, "mu": mu,
           "comm_probe": comm_p, "comm_Q1": comm_1, "sector_dim": int(Up.shape[1]), "core_jump": float(jc),
           "edge_jump": float(jo), "core_weights": [float(x) for x in cw], "occupied": [res[kk][2] for kk in res],
           "velocities": [round(float(v), 4) for v in vel], "n_low": nlow,
           "net_chiral": int(np.sum(vel > 0.1) - np.sum(vel < -0.1)), "spectral_flow": ncross,
           "jump_per_weight": float(abs(jc) / np.sum(cw)) if len(cw) else 0.0}
    P(f"   {label}:")
    P(f"      [H, Q_probe] = {comm_p:.1e}, [H, Q_1] = {comm_1:.2f}; Sektor Q_probe = +1 hat Dimension {Up.shape[1]}.")
    P(f"      Spektralfluss am Kern zwischen k_z = -0,05 und +0,05: {ncross}. Geschwindigkeiten dE/dk_z "
      f"{np.round(vel, 3)}; Kerngewichte {np.round(cw, 3)}; Zustände mit |E| < 0,15: {nlow}.")
    P(f"      Sprung der Kernbesetzung bei k_z = 0: {jc:+.4f}, Außenrand {jo:+.4f}; Sprung / Summe der Kerngewichte "
      f"= {out['jump_per_weight']:.3f}.")
    if pairs:
        P(f"      Satzgewichte der Kernmoden: {[[round(x, 2) for x in ww] for ww in setw]}")
    out["set_weights"] = setw
    return out


runs = {}
runs["n1"] = measure("n = 1", [0.6], [], [1])
runs["n2"] = measure("n = 2, ungekoppelt, Massen 0,6 und 0,5", [0.6, 0.5], [], [1, 1])
runs["n3"] = measure("n = 3, ungekoppelt, Massen 0,6, 0,5, 0,7", [0.6, 0.5, 0.7], [], [1, 1, 1])
runs["n2mix"] = measure("n = 2, gekoppelt (t = 0,15..0,35, mu = 0,15)", [0.6, 0.5], [(0, 1)], [1, 1])
runs["n2mix_mehr"] = measure("n = 2, stärker gekoppelt (t = 0,19..0,44, mu = 0,25)", [0.6, 0.5], [(0, 1)], [1, 1],
                             t=1.25, mu=0.25)
runs["n2mix_grenze"] = measure("n = 2, an der Lückengrenze (t = 0,3..0,7 > m = 0,5, mu = 0,3)", [0.6, 0.5], [(0, 1)],
                               [1, 1], t=2.0, mu=0.3)
runs["n3mix"] = measure("n = 3, alle Paare gekoppelt", [0.6, 0.5, 0.7], [(0, 1), (1, 2), (0, 2)], [1, 1, 1])
P("\n   Verborgener Zweig: drei Sätze winden, die Probe koppelt nur an Satz 1 (Sätze 2 und 3 gekoppelt).")
runs["hidden_probe"] = measure("Probe w = (1, 0, 0)", [0.6, 0.5, 0.7], [(1, 2)], [1, 0, 0])
runs["hidden_all"] = measure("dieselbe Konfiguration, gemeinsame SO(2) (zählt alle Moden)", [0.6, 0.5, 0.7], [(1, 2)],
                             [1, 1, 1])

k = {key: runs[key]["spectral_flow"] for key in runs}
jr = {key: abs(runs[key]["core_jump"]) / abs(runs["n1"]["core_jump"]) for key in runs}
P("\n   Niveau der Probe = Summe q^2 über den Spektralfluss (hier q = 1), daneben das rohe Sprungverhältnis:")
for key in runs:
    P(f"      {key:13s}: Niveau {k[key]}, |Sprung| / |Sprung bei n = 1| = {jr[key]:.3f}, Sprung je Kerngewicht "
      f"{runs[key]['jump_per_weight']:.3f}")
DATA["lattice"] = runs
gapped = [x for x in runs if x != "n2mix_grenze"]
check("[H, Q_probe] = 0 in allen Läufen", all(r["comm_probe"] < 1e-10 for r in runs.values()))
check("Kopplung bricht die einzelnen Satzladungen", runs["n3mix"]["comm_Q1"] > 0.1 and runs["n2mix"]["comm_Q1"] > 0.1)
check("Niveau = n für n = 1, 2, 3, ungekoppelt und gekoppelt",
      [k[x] for x in ("n1", "n2", "n3", "n2mix", "n2mix_mehr", "n2mix_grenze", "n3mix")] == [1, 2, 3, 2, 2, 2, 3])
check("Sprung der Kernbesetzung = Summe der Kerngewichte der kreuzenden Moden (gelückte Läufe)",
      all(0.92 < runs[x]["jump_per_weight"] < 1.02 for x in gapped))
check("alle Kernmoden laufen in dieselbe Richtung", all(runs[x]["net_chiral"] == len(runs[x]["velocities"]) == k[x]
                                                         for x in runs))
check("ungekoppelt v nahe 1 (0,9 bis 1), gekoppelt 0,1 < v < 0,9",
      all(0.9 < v < 1.01 for x in ("n1", "n2", "n3") for v in runs[x]["velocities"])
      and all(0.1 < v < 0.9 for v in runs["n3mix"]["velocities"]))
check("Rand nimmt gegenläufig auf", all(abs(runs[x]["core_jump"] + runs[x]["edge_jump"]) < 0.1 * max(1, runs[x]["n"])
                                        for x in gapped))
check("Kernmoden über die Sätze verteilt (gekoppelt)", any(max(w) < 0.6 for w in runs["n3mix"]["set_weights"]))
check("Besetzungszahl über k_z konstant", all(len(set(r["occupied"])) == 1 for r in runs.values()))
check("verborgener Zweig: Probe sieht Niveau 1, gemeinsame SO(2) sieht 3 Moden",
      k["hidden_probe"] == 1 and k["hidden_all"] == 3)
befund("Lueckengrenze", {"core_weights": runs["n2mix_grenze"]["core_weights"],
                         "velocities": runs["n2mix_grenze"]["velocities"],
                         "jump_per_weight": runs["n2mix_grenze"]["jump_per_weight"]},
       "Bei t > m nähert sich die Volumenlücke dem Schluss: Kerngewicht und Geschwindigkeit fallen, der Sprung ist "
       "nicht mehr am Kern lokalisiert. Der Spektralfluss bleibt 2.")
P("\n   Ungekoppelt ist v = 0,93 bis 1,00 bei k_z = 0,05 (Endlichkeit, am kleinsten für m = 0,5).")
P("   Die Kopplung zwischen Sätzen wirkt für die Satzkombinationen wie ein chemisches Potential: Die")
P("   Geschwindigkeit sinkt, Zahl und Richtung der Kernmoden bleiben. Das Niveau ist eine Eigenschaft der")
P("   Kopplung, nicht der Zerlegung in Sätze: Auch wenn die Kopplung die einzelnen Satzladungen bricht und die")
P("   Kernmoden über alle Sätze verteilt, laufen je Flussquant genau n Moden durch E = 0, und der Kernsprung")
P("   ist ihre volle Ladung. Koppelt die Probe nur an einen Satz, sieht sie Niveau 1, während drei Sätze winden.")

# --------------------------------------------------------------------------------------------------
P("\nB. Vom SO(2)-Probestrom zum E8-Niveau\n")
roots = []
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product((1, -1), repeat=2):
        v = [Fr(0)] * 8
        v[i], v[j] = Fr(si), Fr(sj)
        roots.append(tuple(v))
for signs in itertools.product((1, -1), repeat=8):
    if signs.count(-1) % 2 == 0:
        roots.append(tuple(Fr(s, 2) for s in signs))
dot = lambda a, b: sum(x * y for x, y in zip(a, b))
assert len(roots) == 240 and all(dot(r, r) == 2 for r in roots)
idx_adj = sum(r[0] ** 2 for r in roots)                 # Tr_adj(H^2) für H = e_1 (Plan 1-2)
idx_vec = 2                                             # Gewichte +-e_1 des SO(16)-Vektors
P(f"   H = e_1 dreht zwei Komponenten des 16-Pletts, |H|^2 = 1. Tr_adj(H^2) über die 240 Wurzeln = {idx_adj}")
P(f"   = 2 h^v |H|^2 mit h^v = 30. Im Vektor 16: Tr(H^2) = {idx_vec} = 2 T(16), T(16) = 1.")
P("   Also hat die SO(2) Einbettungsindex 1: Niveau der SO(2) = Niveau von E8. Gemessen: 1 je Satz, n bei")
P("   gleichförmiger Kopplung. Ein E8-Strom, der in jedem Satz gleich wirkt, hat Niveau n.")
check("E8: Tr_adj(H^2) = 2 h^v = 60", idx_adj == 60)
w = (1, Fr(1, 10), Fr(1, 100), Fr(1, 1000), Fr(1, 10 ** 4), Fr(1, 10 ** 5), Fr(1, 10 ** 6), Fr(1, 10 ** 7))
pos = [r for r in roots if dot(r, w) > 0]
rho = tuple(sum(r[i] for r in pos) / 2 for i in range(8))
theta = max(pos, key=lambda r: dot(r, rho))


def weyl_dim(lam):
    num, den = Fr(1), Fr(1)
    for a in pos:
        num *= dot(tuple(x + y for x, y in zip(lam, rho)), a)
        den *= dot(rho, a)
    return num / den


dV1, dV2 = weyl_dim(theta), weyl_dim(tuple(2 * x for x in theta))
w2_generic = 248 * 249 // 2 + 248
P(f"   Weyl-Dimension: V(theta) = {dV1}, V(2 theta) = {dV2}, Höhe von theta = {dot(theta, rho)}.")
P(f"   Gewicht 2 im Vakuummodul von E8 auf Niveau k >= 2: Sym^2(248) + 248 = {w2_generic}.")
P(f"   Nur bei k = 1 fällt V(2 theta) heraus (Z: CELEST.WP5B.01, Faktor 2(1-k)): {w2_generic} - {dV2} = "
  f"{w2_generic - int(dV2)}.")
check("dim V(theta) = 248, dim V(2 theta) = 27000", dV1 == 248 and dV2 == 27000)


def c_sug(kk):
    return Fr(248 * kk, kk + 30)


c3sug, coset = c_sug(3), 24 - c_sug(3)
P(f"   Sugawara: c(E8_1) = {c_sug(1)}, c(E8_3) = {c3sug} = {float(c3sug):.4f}. Bei n = 3 bleibt der Coset")
P(f"   E8_1^3 / E8_3 mit c = {coset} = 1/2 + 21/22 (Ising und M(11,12); E8_1^2/E8_2 hat {16 - c_sug(2)}).")
mm = lambda m: 1 - Fr(6, m * (m + 1))                  # Minimalmodell M(m, m+1)
check("Coset c = 16/11 = c(M(3,4)) + c(M(11,12))", coset == Fr(16, 11) and 16 - c_sug(2) == mm(3)
      and 8 + c_sug(2) - c_sug(3) == mm(11) and mm(3) + mm(11) == coset)
P("   Bei n = 3 mit gleichförmiger Kopplung ist die Naht also keine reine Stromalgebra des sichtbaren E8:")
P("   T = T_Sugawara(E8_3) + T_Coset. Bei n = 1 gilt T = T_Sugawara(E8_1) mit c = 8.")

# --------------------------------------------------------------------------------------------------
P("\nC. TFPT-Beziehungen bei Niveau 3 beziehungsweise c = 24\n")
mp.mp.dps = 30
c3 = 1 / (8 * mp.pi)
qi = mp.e ** (-2 * mp.pi)
E2i = 1 - 24 * mp.nsum(lambda nn: mp.mpf(sum(d for d in range(1, int(nn) + 1) if int(nn) % d == 0)) * qi ** nn, [1, 40])
P(f"   C1. E2(i) = {mp.nstr(E2i, 20)}, 3/pi = {mp.nstr(3 / mp.pi, 20)}; E2(i)/24 = c3.")
check("E2(i) = 3/pi", abs(E2i - 3 / mp.pi) < mp.mpf("1e-25"))


def theta_E8(y):
    x = mp.e ** (-mp.pi * y)
    return (mp.jtheta(2, 0, x) ** 8 + mp.jtheta(3, 0, x) ** 8 + mp.jtheta(4, 0, x) ** 8) / 2


def theta_D8(y):
    x = mp.e ** (-mp.pi * y)
    return (mp.jtheta(3, 0, x) ** 8 + mp.jtheta(4, 0, x) ** 8) / 2


def mean_energy(f):
    return -mp.diff(lambda y: mp.log(f(y)), 1) / (2 * mp.pi)


eE8 = mean_energy(theta_E8)
eE83 = mean_energy(lambda y: theta_E8(y) ** 3)
eD8 = mean_energy(theta_D8)
eS3 = mean_energy(lambda y: theta_E8(3 * y))
P(f"   Mittlere Gitterladungsenergie bei tau = i je Richtung: E8 {mp.nstr(eE8 / 8, 15)}, E8^3 {mp.nstr(eE83 / 24, 15)},")
P(f"   c3 = {mp.nstr(c3, 15)}. Gesamt: E8 {mp.nstr(eE8, 10)} = 1/pi, E8^3 {mp.nstr(eE83, 10)} = 3/pi.")
P(f"   Gegenproben: D8 {mp.nstr(eD8 / 8, 8)} je Richtung; das Cartangitter sqrt(3) E8 des Niveaus 3 für sich "
  f"{mp.nstr(eS3 / 8, 6)}.")
# isotrope zweite Momente aus Aufzählung, dann die Diagonalrichtung in E8^3
def shell(vals, maxn=12):
    arr = np.zeros((1, 0))
    for _ in range(8):
        arr = np.array([np.append(a, v) for a in arr for v in vals if (a ** 2).sum() + v * v <= maxn + 1e-9])
    return arr


ints = shell(range(-3, 4))
ints = ints[ints.sum(1) % 2 == 0]
halfs = shell((-2.5, -1.5, -0.5, 0.5, 1.5, 2.5))
halfs = halfs[np.round(halfs.sum(1)) % 2 == 0]
vecs = np.vstack([ints, halfs])
wts = np.exp(-np.pi * (vecs ** 2).sum(1))
Z = wts.sum()
M2 = (vecs.T * wts) @ vecs / Z
m1 = (vecs.T * wts).sum(1) / Z
e = np.ones(8) / np.sqrt(8)
diag_dir = (3 * e @ M2 @ e + 6 * (e @ m1) ** 2) / 6
P(f"   Aufzählung von {len(vecs)} E8-Vektoren (Norm <= 12): zweites Moment = {M2[0, 0]:.12f} x Eins "
  f"(Abweichung {np.max(np.abs(M2 - M2[0, 0] * np.eye(8))):.1e}), 1/(4 pi) = {1 / (4 * np.pi):.12f}.")
P(f"   Diagonale Cartanrichtung des gemeinsamen E8 in E8^3: Energie {diag_dir:.12f} = c3.")
check("c3 je Cartanrichtung für E8 und E8^3", abs(eE8 / 8 - c3) < mp.mpf("1e-20") and abs(eE83 / 24 - c3) < mp.mpf("1e-20"))
check("Diagonalrichtung in E8^3 trägt c3", abs(diag_dir - 1 / (8 * np.pi)) < 1e-12)
check("Gegenproben D8 und sqrt(3) E8 weichen ab", abs(eD8 / 8 - c3) > 0.01 and abs(eS3 / 8 - c3) > 0.01)
befund("c3_trennt_n_nicht", True,
       "c3 = E2(i)/24 je Cartanrichtung gilt für jedes unimodulare Gitter am selbstdualen Punkt, also für E8 "
       "und E8^3, auch entlang der Diagonalen. Es unterscheidet n = 1 und n = 3 nicht; nur die Summe 1/pi "
       "gegen 3/pi ändert sich, und sie geht in keine TFPT-Zahl ein. Ein für sich stehendes Niveau-3-Gitter "
       "sqrt(3) E8 trägt c3 nicht.")

P("\n   C2. Charaktere.")
NQ = 5
sig3 = [0] + [sum(d ** 3 for d in range(1, nn + 1) if nn % d == 0) for nn in range(1, NQ + 1)]
E4s = sp.Poly(1 + 240 * sum(sig3[nn] * sp.Symbol("q") ** nn for nn in range(1, NQ + 1)), sp.Symbol("q"))
qq = sp.Symbol("q")
prod8 = sp.series(sp.prod([(1 - qq ** nn) ** (-8) for nn in range(1, NQ + 1)]), qq, 0, NQ).removeO()
chi = sp.series(E4s.as_expr() * prod8, qq, 0, NQ).removeO()
chic = [int(sp.Poly(chi, qq).coeff_monomial(qq ** nn)) for nn in range(4)]
chi3 = sp.series(chi ** 3, qq, 0, 3).removeO()
chi3c = [int(sp.Poly(chi3, qq).coeff_monomial(qq ** nn)) for nn in range(3)]
eta_i = mp.gamma(mp.mpf(1) / 4) / (2 * mp.pi ** (mp.mpf(3) / 4))
E4i = 1 + 240 * mp.nsum(lambda nn: mp.mpf(sum(d ** 3 for d in range(1, int(nn) + 1) if int(nn) % d == 0)) * qi ** nn,
                        [1, 40])
chi_i = E4i / eta_i ** 8
j_i = E4i ** 3 / eta_i ** 24
P(f"   chi_E8 = E4/eta^8: Koeffizienten {chic}; chi_E8^3 = j: {chi3c} (744 = 3 x 248, 196884 = 3 x 4124 + 3 x 248^2).")
P(f"   chi_E8(i) = {mp.nstr(chi_i, 15)}, j(i) = chi^3 = {mp.nstr(j_i, 15)}.")
P("   Bei n = 3 ist die Zustandssumme der Naht j(tau) selbst; ihr Wert 1728 ist TFPTs Kennzeichen des")
P("   quadratischen Moduls (QGEO.PILLOW.01). Bei n = 1 ist es die Kubikwurzel 12. Keine TFPT-Zahl hängt")
P("   an der Wahl; in Z9 verschiebt sich nur die Schwelle für unterdrückte Nähte von log 12 = 2,48 auf")
P("   log 1728 = 7,45.")
P(f"   Die Kill-Probe SEAM.EQUIV.KILLTEST.01 nennt 248 Ströme vom Gewicht 1. Bei n = 3 sind es 744; das")
P(f"   sichtbare E8 auf Niveau 3 hat im Vakuummodul bei Gewicht 2 {w2_generic} statt {chic[2]} Zustände.")
check("chi = 1, 248, 4124, 34752", chic == [1, 248, 4124, 34752])
check("chi^3 = j: 1, 744, 196884", chi3c == [1, 744, 196884] and 196884 == 3 * 4124 + 3 * 248 ** 2)
check("chi(i) = 12, j(i) = 1728", abs(chi_i - 12) < mp.mpf("1e-20") and abs(j_i - 1728) < mp.mpf("1e-17"))
check("4124 = Gewicht 2 nur bei Niveau 1", chic[2] == w2_generic - int(dV2))
befund("chi_i_trennt_n_nicht", True,
       "chi_E8(i) = 12 und j(i) = 1728 sind beide TFPT-Zahlen; die Satzzahl wählt nur, welche die Zustandssumme "
       "der Naht ist. Kein Vorhersagewert hängt daran.")

P("\n   C3. Holomorphe Eindeutigkeit (SEAM.EQUIV.UNIQ.01).")
roots_D = lambda nn: 2 * nn * (nn - 1)
v1 = {"E8^3 (GSO je Satz)": 24 + 3 * 240, "D16+ x E8 (GSO über 32 + 16)": 24 + roots_D(16) + 240,
      "D24+ (GSO über alle 48)": 24 + roots_D(24)}
P("   Bei c = 8 gibt es genau eine holomorphe Theorie, E8_1. Aus 48 Majorana-Weyl Fermionen entstehen je nach")
P("   GSO-Projektion mindestens drei verschiedene holomorphe Theorien mit c = 24:")
for kk, vv in v1.items():
    P(f"      {kk:30s}: dim V_1 = {vv}")
P("   (Schellekens: 71 Einträge bei c = 24.) Bei n = 3 braucht die Wahl E8^3 eine eigene Regel: eine GSO-Summe je")
P("   Satz. Bei n = 1 legt die Holomorphie allein E8 fest.")
check("c = 24: verschiedene freie Fermionkonstruktionen", len(set(v1.values())) == 2 and v1["D24+ (GSO über alle 48)"] == 1128
      and v1["E8^3 (GSO je Satz)"] == 744)

P("\n   C4. c = 8 = g_car + N_fam (GRAV.ENTROPY.EQUILIBRIUM.01, SEAM.EQUIV.EDGE.ENTANGLEMENT.01).")
g_car, N_fam = 5, 3
P(f"   n = 1: c = 8 = g_car + N_fam. n = 3: c = 24 = N_fam (g_car + N_fam) = {N_fam * (g_car + N_fam)}.")
P("   Der Einsteinkoeffizient 1/c3 = 8 pi und die Dichte 1/4 = 1/|mu4| hängen nicht an c. c tritt nur als")
P("   Lesart der Verschränkungsnormierung auf; v_geo ist offen.")
befund("c_gleich_gcar_plus_Nfam_ist_Lesart", True,
       "Bei n = 3 wird aus c = g_car + N_fam die Beziehung c = N_fam (g_car + N_fam); keine Zahl der "
       "Gravitationsherleitung ändert sich. Die Beziehung bevorzugt n = 1 nur als Lesart.")

P("\n   C5. Das alpha^3-Glied als U(1)-Niveau (EM.WARD.01, ALPHA.QUILLEN.PROGRESS.04, ALPHA.QUILLEN.INFLOW.01).")
P("   F(alpha) = k0 alpha^3 - 2 c3^3 alpha^2 - 8 b1 c3^6 log(1/phi_seam(alpha)), wie in verification/v3_em_alpha.py.")
P("   TFPT liest k0 = 1 als Niveau beziehungsweise Chernzahl |C| = 1 des Kragens; 2 = |Z2|; 8 b1 = (4/5) 41 aus")
P("   dem Materieinhalt.")
CODATA, SIG = mp.mpf("137.035999177"), mp.mpf("0.000000021")


def alpha_inv(k0=1, kc=1, kt=1, dtf=48, b1=mp.mpf(41) / 10):
    pb, dt = 1 / (6 * mp.pi), dtf * c3 ** 4

    def ps(a):
        Q = dt * mp.e ** (-2 * a)
        return pb + Q * (1 - Q) ** (mp.mpf(-5) / 4)

    F = lambda a: k0 * a ** 3 - kc * 2 * c3 ** 3 * a ** 2 - kt * 8 * b1 * c3 ** 6 * mp.log(1 / ps(a))
    return 1 / mp.findroot(F, mp.mpf("0.0073") / mp.mpf(k0) ** (mp.mpf(1) / 3))


variants = [
    ("Referenz (n = 1, k0 = 1)", {}),
    ("n = 3 gleichförmig: k0 = 3", {"k0": 3}),
    ("k0 = 3 und Calderónglied x 3", {"k0": 3, "kc": 3}),
    ("k0 = 3, b1 = (3/(5k)) 41/6 = 41/30", {"k0": 3, "b1": mp.mpf(41) / 30}),
    ("alle drei Glieder x 3 (nur Umnormierung)", {"k0": 3, "kc": 3, "kt": 3}),
    ("verborgene Sätze: sichtbares k0 = 1", {}),
]
alpha_tab = {}
for name, kw in variants:
    ai = alpha_inv(**kw)
    zz = (ai - CODATA) / SIG
    alpha_tab[name] = {"alpha_inv": float(ai), "sigma": float(zz)}
    P(f"      {name:42s}: alpha^-1 = {mp.nstr(ai, 13):>16s}, Zug {mp.nstr(zz, 4):>10s} sigma")
check("Referenz reproduziert alpha^-1 = 137,0359992168", abs(alpha_tab["Referenz (n = 1, k0 = 1)"]["alpha_inv"]
                                                         - 137.0359992168) < 1e-9)
check("Niveau 3 im alpha^3-Glied zerstört alpha (> 10^6 sigma)", abs(alpha_tab["n = 3 gleichförmig: k0 = 3"]["sigma"]) > 1e6
      and abs(alpha_tab["k0 = 3 und Calderónglied x 3"]["sigma"]) > 1e6)
P("   alpha übersteht Niveau 3 nur, wenn alle drei Glieder mit der Satzzahl wachsen. Das widerspricht ihrer")
P("   Herkunft: 2 = |Z2| (zwei Blätter) und 41 (drei Generationen und Higgs) zählen keine Nahtsätze.")
befund("alpha_ueberlebt_bei_gleichmaessiger_Umnormierung", True,
       "Werden alle drei Glieder der Wardidentität mit 3 multipliziert, bleibt alpha^-1 unverändert. Nur die "
       "Typisierung der Glieder (EM.WARD.01) schließt das aus; die EM1-Lesart k0 = Niveau ist [C], nicht bewiesen.")

P("\n   C6. Seedterm 48 c3^4 in phi0 = 1/(6 pi) + 48 c3^4.")
seed_tab = {}
for name, dtf in (("48 (Kumulant je normierter Richtung, niveaublind)", 48), ("16 (16 Nahtmoden, n = 1 gezählt)", 16),
                  ("48 = 3 x 16 (Nahtmoden bei n = 3)", 48), ("144 = 3 x 48 (Summe über drei Sätze)", 144),
                  ("3888 = 48 x 3^4 (Diagonalladung unnormiert)", 48 * 81)):
    ai = alpha_inv(dtf=dtf)
    zz = (ai - CODATA) / SIG
    seed_tab[name] = {"alpha_inv": float(ai), "sigma": float(zz)}
    P(f"      {name:46s}: alpha^-1 = {mp.nstr(ai, 13):>16s}, Zug {mp.nstr(zz, 4):>10s} sigma")
check("nur der Faktor 48 trägt alpha", all(abs(v["sigma"]) > 1e5 for kk, v in seed_tab.items() if not kk.startswith("48")))
P("   Die Zahl 48 ist starr. Ihre Zerlegung legt die Satzzahl nicht fest: TFPT liest 48 = 2^3 3! (Kumulant")
P("   von X^2, X ~ N(0, c3)) oder 48 = 16 x 3 (Generation mal Familien, radical-global-hypotheses check_three,")
P("   W1). Als Zählung von Nahtmoden passte sie zu n = 3 (nu = 48), nicht zu n = 1 (16). Nach Z10 sind die")
P("   16 der Naht aber nicht die 16 einer Generation; die Zählung als Nahtmoden ist nicht belegt.")
befund("seed_48_entscheidet_n_nicht", True,
       "48 c3^4 ist starr (16, 144, 3888 scheitern um > 10^5 sigma). Als Zählung der Nahtmoden würde 48 = nu(n=3) "
       "n = 3 stützen; TFPTs eigene Lesarten (Kumulant, 16 x N_fam Materie) sind satzzahlblind.")

P("\n   C7. R102: D5 + A3 mit Z4-Glue gibt E8.")
cart = lambda name: {"A3": sp.Matrix([[2, -1, 0], [-1, 2, -1], [0, -1, 2]]),
                     "D5": sp.Matrix([[2, -1, 0, 0, 0], [-1, 2, -1, 0, 0], [0, -1, 2, -1, -1], [0, 0, -1, 2, 0],
                                      [0, 0, -1, 0, 2]]),
                     "E8": sp.Matrix([[2, -1, 0, 0, 0, 0, 0, 0], [-1, 2, -1, 0, 0, 0, 0, 0], [0, -1, 2, -1, 0, 0, 0, -1],
                                      [0, 0, -1, 2, -1, 0, 0, 0], [0, 0, 0, -1, 2, -1, 0, 0], [0, 0, 0, 0, -1, 2, -1, 0],
                                      [0, 0, 0, 0, 0, -1, 2, 0], [0, 0, -1, 0, 0, 0, 0, 2]])}[name]
dA3, dD5, dE8 = cart("A3").det(), cart("D5").det(), cart("E8").det()
idx = sp.sqrt(dA3 * dD5 / dE8)
d3E8 = (3 * cart("E8")).det()
P(f"   det A3 = {dA3}, det D5 = {dD5}, det E8 = {dE8}: Index von D5 + A3 in E8 = {idx}. Ein Satz (10,1)+(1,6) ist")
P("   genau D5 + D3 = D5 + A3 in E8 ⊃ D8: ein Träger (Rang 5 = g_car) und ein Familiengitter (Rang 3 = N_fam).")
P(f"   Drei Sätze brauchen Rang 24 = 3 (5 + 3), also drei Träger. Das Cartangitter eines Niveau-3-E8 ist sqrt(3) E8")
P(f"   mit Wurzelnorm 6 und det = {d3E8} = 3^8; R102 setzt Wurzelnorm 2 voraus und gilt dort nicht.")
check("R102-Index 4 und det(sqrt3 E8) = 3^8", idx == 4 and dE8 == 1 and d3E8 == 3 ** 8)

P("\n   C8. Eine Phase, zwei Antworten (ALPHA.QUILLEN.INFLOW.01): c_- = 16 x 1/2 x |C|, also c_-/|C| = 8.")
ratio = {"n = 1": Fr(8, 1), "n = 3 gleichförmig": Fr(24, 3), "n = 3 verborgen": Fr(24, 1)}
for kk, vv in ratio.items():
    P(f"      {kk:20s}: c_- / k_U(1) = {vv}")
lat_ratio_hidden = k["hidden_all"] / k["hidden_probe"]
P(f"   Im Gitter (2 von 16 Komponenten je Satz, c_red = Zahl der komplexen Kernmoden): gleichförmig ist c_red/k = 1,")
P(f"   weil dieselbe Ladung alle Moden zählt; verborgen {lat_ratio_hidden:.0f}. Mal 8 für volle Sätze: 8 gegen 24.")
P("   Das Verhältnis c_-/k trennt n = 1 vom verborgenen Zweig; k allein trennt n = 1 vom gleichförmigen.")

P("\n   C9. T-Phase der Naht allein: e^(-2 pi i c/24).")
for nn in (1, 3):
    ph = mp.exp(-2j * mp.pi * 8 * nn / 24)
    P(f"      n = {nn}: {mp.nstr(ph, 6)}")
P("   n = 1 braucht Zufluss (gerechnet für SO(2) in Z9, gravitativ offen); n = 3 ist allein modular invariant.")

# --------------------------------------------------------------------------------------------------
P("\nD. Urteil\n")
P("   1. Gitter: Ein Strom, der in jedem windenden Satz gleich wirkt, hat auf der Naht Niveau n. Gemessen für")
P(f"      n = 1, 2, 3 (Spektralfluss {k['n1']}, {k['n2']}, {k['n3']}; gekoppelt {k['n2mix']}, {k['n3mix']}; rohe "
  f"Sprungverhältnisse {jr['n2']:.2f}, {jr['n3']:.2f}, {jr['n3mix']:.2f}). Mit Einbettungsindex 1 gilt das für E8.")
P("   2. Gleichförmiges n = 3 (drei identische Sätze, an die dieselben Eichfelder koppeln) hat sichtbares Niveau 3.")
P(f"      TFPT liest das alpha^3-Glied als dieses Niveau. Mit k0 = 3 ist alpha^-1 = "
  f"{alpha_tab['n = 3 gleichförmig: k0 = 3']['alpha_inv']:.4f}. Ebenso verlangen")
P("      CELEST.WP5B.01 (nur k = 1 tilgt V(2 theta), 4124), SEAM.EQUIV.EDGE.VIRASORO.01 (T = T_Sugawara, c = 8),")
P("      SEAM.EQUIV.UNIQ.01 (Holomorphie legt E8 nur bei c = 8 fest) und R102 (Wurzelnorm 2) Niveau 1.")
P("      Dieser Zweig scheitert unter der EM1-Lesart [C].")
P("   3. n = 3 mit einem sichtbaren und zwei verborgenen Sätzen hat sichtbares Niveau 1 und c = 24. alpha, phi0,")
P("      c3 und chi(i) überstehen das. Dagegen sprechen nur Lesarten: c = g_car + N_fam, c_-/|C| = 8 (eine Phase),")
P("      ein Träger und ein Familiengitter in R102, Minimalität.")
P("   4. TFPT legt also das sichtbare Niveau auf 1 fest, nicht c. Entscheiden würde die gravitative Antwort der")
P("      Naht: c_- = 8 oder 24 bei k = 1, also c_-/k = 8 oder 24. Rechnung: der gravitative Zufluss aus der")
P("      Determinante der schweren 4D Fermionen (Koeffizient c/24 von p1), oder im Gitter die thermische")
P("      Hallantwort der Kernmoden, die auch ungeladene Moden zählt.")
verdict = {
    "gleichfoermig_n3": "gescheitert unter der EM1-Lesart (k0 = Niveau, [C]); alpha^-1 = 198,24",
    "verborgen_n3": "offen; alle Zahlen bestehen, nur Lesarten (c = g_car + N_fam, c_-/|C| = 8, R102-Zählung) dagegen",
    "n1": "Lesart; mit allen geprüften TFPT-Beziehungen verträglich, braucht gravitativen Zufluss (nicht gerechnet)",
    "gesamt": "Lesart: TFPT fixiert das sichtbare E8-Niveau 1 (über alpha unter EM1), nicht die Satzzahl; "
              "entscheidend ist c_-/k der Naht (8 gegen 24)",
    "entscheidende_rechnung": "gravitativer Zufluss der schweren 4D Fermionen beziehungsweise thermische Hallantwort "
                              "der Naht: c_- = 8 (n = 1) oder 24 (n = 3) bei sichtbarem Niveau 1",
}
befund("verborgener_Zweig_besteht_alle_Zahlen", True, verdict["verborgen_n3"])
befund("gravitativer_Zufluss_nicht_gerechnet", True, "c_- der Naht ist im Gitter nur gezählt, nicht als Antwort "
       "(Energiestrom, p1-Koeffizient) aus dem 4D Volumen gerechnet.")
befund("Gitter_reduziert", True, "Je Satz 2 von 16 Komponenten; die volle 16 ist über Additivität und den "
       "Einbettungsindex 1 argumentiert, nicht gerechnet. Die Probe ist eine SO(2), kein nichtabelscher E8-Strom.")

P("\nZusammenfassung")
for nm, c in CHECKS.items():
    P(f"   {'bestanden' if c else 'NICHT bestanden'}: {nm}")
P("Befunde (nicht als Prüfung gezählt)")
for nm, b in BEFUNDE.items():
    P(f"   {nm}: {b['wert']}")
P(f"\nLaufzeit {time.time() - T0:.0f} s")

(HERE / "ausgabe_z14.txt").write_text("\n".join(LINES) + "\n")
(HERE / "z14_ergebnis.json").write_text(json.dumps({
    "checks": CHECKS, "befunde": BEFUNDE, "urteil": verdict,
    "niveau": k, "sprungverhaeltnis": jr, "gitter": DATA["lattice"],
    "alpha": alpha_tab, "seed": seed_tab,
    "c3_je_richtung": {"E8": float(eE8 / 8), "E8^3": float(eE83 / 24), "diag_in_E8^3": float(diag_dir),
                       "D8": float(eD8 / 8), "sqrt3_E8": float(eS3 / 8), "c3": float(c3)},
    "charaktere": {"chi": chic, "chi3": chi3c, "chi_i": float(chi_i), "j_i": float(j_i), "gewicht2_niveau_ge2": w2_generic,
                   "dim_V_2theta": int(dV2)},
    "sugawara": {"E8_1": str(c_sug(1)), "E8_3": str(c3sug), "coset_n3": str(coset)},
    "holomorph_c24_dimV1": v1,
}, ensure_ascii=False, indent=1, default=str))
assert all(CHECKS.values())
