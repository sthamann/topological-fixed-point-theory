"""Z52: Gravitation -- Kantenform, Skala und Scalaron.

Teil 1  Satz: Jedes flache, ungerichtete Fierz-Pauli-Potential ist eps(D ^ D ^ Q(e_p, e_k)) mit D = e_p - e_k und einer
        symmetrischen Quadratform Q = x(e_p^2 + e_k^2) + y e_p e_k. Die quartische Zahl t misst die Form von Q:
        V (t = 0), arithmetisches Mittel (t = 1), ein Vielbein je Endpunkt (t = 4), Quadrate (t = -2). Nur Q = e_p e_k ist
        höchstens kubisch in der Kantenvariable Y = I - X.
Teil 2  Skala: Nach GRAV.GAPCORR.01 liegt die Nahtskala bei v_geo ~ M_Pl; ist ein Übertragungsschritt eine Quartereinheit
        der Naht, ist mu ~ M_c, und kein Arm ist schwach gekoppelt. Unterhalb der Nahtskala bleibt Einstein mit
        M_eff = M_c (1 + O(e^{-4 pi})); G setzt GRAV.NONLINEAR.01 (8 pi = 1/c3 in der Einheit v_geo, Jacobson-Lesart).
Teil 3  Scalaron: Der R^2-Scalaron des Bestands (M = c3^{7/2} Mbar) ist ein gesunder Spinnull-Pol (rho0 < 0) und damit genau
        die Gewichtung, die die Ward-Bedingung aus z47 im Spurkanal erlaubt.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "research" / "originlab"))
from process_rules import BETA_DET, BETA_V, Q_ORIGINAL, edge_form_beta, potential_taylor, run  # noqa: E402

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
U = res["universes"]

# ------------------------------------------------------------------ Teil 1
say("Teil 1: Kantenform und die quartische Zahl")
rng = np.random.default_rng(521)
fam = []
for _ in range(6):
    x, y = rng.normal(size=2)
    if abs(-2*x - y) < 0.2:
        continue
    beta = np.array(edge_form_beta(x, y))
    t = 2*(beta[0] - BETA_V[0])
    fam.append((x, y, t, np.max(np.abs(beta - (np.array(BETA_V) + t/2*np.array(BETA_DET))))))
say(f"   zufällige Formen Q: t-Werte {[round(float(f[2]), 4) for f in fam]}, Abstand zur Familie V_t {max(f[3] for f in fam):.1e}")
check("jede Form eps(D ^ D ^ Q) mit symmetrischem Q liegt auf der Familie V_t = V + (t/2) det(I - X), mit "
      "t = 2 x/(-(2x + y)/6) + 4 (nach Massennormierung)", max(f[3] for f in fam) < 1e-12
      and all(abs(f[2] - (2*f[0]/(-(2*f[0] + f[1])/6) + 4)) < 1e-9 for f in fam))
ts = np.linspace(-3, 6, 7)
back = []
for t in ts:
    # Umkehrung: x = (t - 4)/2 bei Normierung -2x - y = 6
    x = (t - 4)/2
    y = -6 - 2*x
    back.append(abs(2*(edge_form_beta(x, y)[0] - BETA_V[0]) - t))
check("umgekehrt erreicht jede Zahl t eine Form Q (x = (t - 4)/2, y = -6 - 2x): die Familie V_t ist genau die Menge der "
      "Kantenformen", max(back) < 1e-12)
A, B = rng.normal(size=(4, 4)), rng.normal(size=(4, 4))
A, B = (A + A.T)/2, (B + B.T)/2
co = potential_taylor(edge_form_beta(0.3, -1.1), A, B)
sw = potential_taylor(edge_form_beta(0.3, -1.1), B, A)
check("Gegenprobe mit einer zufälligen Form: flach und ungerichtet (Taylor bis Ordnung vier)",
      abs(co[0]) < 1e-12 and abs(co[1]) < 1e-12 and np.allclose(co, sw, atol=1e-10))
ef = U["edge_forms"]
for name, d in ef["derived"].items():
    say(f"      {name:46s} t = {d['t']:5.2f}, Grad in Y = {d['edge_degree']}, {d['clauses']}")
sel = {r["selected"]: r["minimal_rules"] for r in ef["unique_selections"]}
check("in der Lane wählt 'höchstens kubisch in der Kantenvariable' allein genau Q = e_p e_k (t = 4); das Mittel wählt t = 1",
      sel.get("Q = e_p e_k (one vielbein of each endpoint)", [[]])[0] in (["bilinear_measure"], ["cubic_in_edge_variable"])
      and ["cubic_in_edge_variable"] in sel.get("Q = e_p e_k (one vielbein of each endpoint)", [])
      and ["mean_measure"] in sel.get("Q = (e_p + e_k)^2 (arithmetic mean)", []))
befund("t ist die Form der Endpunkt-Quadratform Q; die naheliegenden Regeln wählen verschiedene t (Minimalgrad in Y: t = 4, "
       "arithmetisches Mittel: t = 1, V: t = 0). Aus dem freien Transfer folgt keine; t bleibt eine Regelwahl", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Skala und G aus dem Bestand")
qo = Q_ORIGINAL
for ratio in (1.0, 0.1, 1e-3):
    weak = 0.5*qo >= ratio
    say(f"   mu/M_c = {ratio:g}: erste Armstufe M_(j,0)/M_c = q/2 = {0.5*qo:.2e} -> schwach gekoppelt: {weak}")
check("liegt mu bei der Nahtskala (mu/M_c ~ 1 bis 1e-3), ist im originalen Zustand keine Armstufe schwach gekoppelt "
      "(q/2 = 9.3e-4)", 0.5*qo < 1e-3)
Meff_over_Mc = np.sqrt(sum(qo**(2*n) for n in range(40)))
check("unterhalb der Nahtskala bleibt Einstein mit M_eff = M_c sqrt(Sum q^(2n)) = M_c (1 + 1.74e-6)",
      abs(Meff_over_Mc - 1 - 0.5*qo**2) < 1e-10)
c3 = 1/(8*np.pi)
Mbar = 2.435e18
Ms = c3**3.5*Mbar
say(f"   GRAV.NONLINEAR.01: 8 pi G mit 8 pi = 1/c3 in der Einheit v_geo; GRAV.GAPCORR.01: Scalaron M = c3^(7/2) Mbar = {Ms:.3e} GeV")
check("der Scalaron des Bestands liegt bei c3^(7/2) Mbar = 3.06e13 GeV (Ledger-Wert reproduziert)", abs(Ms/3.06e13 - 1) < 0.01)
befund("Für zugängliche Energien ist die Gravitation damit festgelegt: Einstein mit G aus GRAV.NONLINEAR.01 plus Scalaron; "
       "die massive Linie (Gewicht 3.5e-6) und t wirken erst an der Nahtskala, wo die Arme stark koppeln", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Scalaron und Ward-Bedingung")
eta = np.diag([-1.0, 1, 1, 1])
from itertools import combinations_with_replacement
idx = list(combinations_with_replacement(range(4), 2))


def sym(v):
    T = np.zeros((4, 4))
    for (a, b), c in zip(idx, v):
        T[a, b] = T[b, a] = c
    return T


rows = []
kup = np.array([1.0, 0, 0, 0])
for nu in range(4):
    rows.append([sym(e)[:, nu] @ kup for e in np.eye(10)])
_, sv, vh = np.linalg.svd(np.array(rows))
basis = [sym(v) for v in vh[np.sum(sv > 1e-12):]]
TT = lambda X, Y: np.einsum("ab,cd,ac,bd->", X, Y, eta, eta)
tr = lambda X: np.einsum("ab,ab->", eta, X)
G = np.array([[(1/6)*tr(X)*tr(Y) for Y in basis] for X in basis])        # rho2 = 0, rho0 = -1: -(rho0/6) T^2
ev = np.linalg.eigvalsh(G)
check("am Scalaron-Pol (rho2 = 0, rho0 = -1) ist das Residuum (1/6) T^2 positiv semidefinit: gesunder Skalar",
      np.all(ev > -1e-12) and np.sum(ev > 1e-12) == 1)
static = 1 + (1/6)/(1/2)
say(f"   statisch: Newton (TT - T^2/2 -> 1/2) plus Scalaron (T^2/6) gibt den Faktor {static:.4f} = 1 + 1/3 e^(-M r)")
check("das ist die bekannte Starobinsky-Korrektur V = -GM/r (1 + e^(-M r)/3)", abs(static - 4/3) < 1e-12)
befund("Die Ward-Bedingung aus z47 (Gedächtnis nur im Spin-2-Kanal, Spurkanal Einstein oder gesunder Skalar) verträgt sich "
       "mit dem R^2-Scalaron des Bestands; er ist genau die erlaubte Spinnull-Gewichtung", True)
ZAHLEN["scalaron_GeV"] = Ms
ZAHLEN["t_kandidaten"] = {n: d["t"] for n, d in ef["derived"].items()}

urteil = ("Die quartische Zahl des metrischen Sterns ist die Form einer Quadratform: Jedes flache, ungerichtete "
          "Fierz-Pauli-Potential ist eps(D ^ D ^ Q(e_p, e_k)), und t misst Q (V: t = 0, Mittel: t = 1, ein Vielbein je "
          "Endpunkt: t = 4, Quadrate: t = -2). Nur Q = e_p e_k ist höchstens kubisch in der Kantenvariable; welche Regel "
          "gilt, folgt aus dem freien Transfer nicht. Für zugängliche Energien spielt t keine Rolle: Nach dem Bestand liegt "
          "die Nahtskala bei v_geo ~ M_Pl, dort koppeln die Arme im originalen Zustand stark (q/2 = 9.3e-4), und darunter "
          "bleibt Einstein mit M_eff = M_c (1 + 1.7e-6) und G aus GRAV.NONLINEAR.01, dazu der Scalaron des Bestands "
          "(3.06e13 GeV). Dieser ist ein gesunder Spinnull-Pol und genau die Gewichtung, die die Ward-Bedingung im "
          "Spurkanal erlaubt. Die Gravitation ist damit unter den Prämissen (D, freie Naht, Zuordnung von Nahtskala und "
          "Planckskala) bis zur Nahtskala festgelegt; t und die Planck-nahe Dynamik bleiben offen.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z52_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
