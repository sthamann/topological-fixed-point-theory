"""Z47: Ward-Struktur für das Spin-2-Gedächtnis der Viertelübertragung.

Teil 1  Index-Schutz: Der masselose Pol ist der Nullvektor eines Fredholmoperators mit Index +1. Kleine Störungen des
        Symbols und endlichrangige Störungen lassen ihn bestehen (normierbarer Nullvektor), ändern aber sein Gewicht Z:
        Masselosigkeit ist geschützt, G nicht.
Teil 2  Universalität: Weinbergs Softbedingung erzwingt gleiche Kopplungen aller Arten. Die Leiterquellen |4j> koppeln
        mit sqrt(Z) q^j; universell ist nur eine einzige Quelle, der Randzustand |0> (Materie auf der Naht).
Teil 3  Ward plus Positivität: Mit erhaltener Quelle ist die Austauschamplitude T[C2 P2 - C0 P0/2]T. Ein massiver Pol
        mit gleichem Formfaktor in beiden Kanälen trägt einen Skalargeist (Gewicht -rho/2 auf der Spur); gesund ist das
        Gedächtnis nur im Spin-2-Kanal, der Spurkanal bleibt reiner Einsteinpol mit demselben Residuum.
Teil 4  Kovariante Form: Das linearisierte Weyl-Quadrat ist proportional zu k^4 h P2 h; ein Term W F(box) W trägt das
        Gedächtnis diffeomorphismeninvariant nur im Spin-2-Kanal.
"""
import json
from itertools import combinations_with_replacement
from pathlib import Path

import numpy as np

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
HERE = Path(__file__).parent
rng = np.random.default_rng(47)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Index-Schutz des masselosen Pols")
N = 240


def section(g, shift=1, K=None):
    """(N - shift) x N section of S^{+shift} T_g (exact: every image component is kept)."""
    T = np.zeros((N + shift, N), complex)
    for j in range(N):
        for k, gk in enumerate(g):
            if j + k < N + shift:
                T[j + k, j] += gk
    M = T[shift:shift + N - 1, :] if shift else T[:N - 1, :]
    if K is not None:
        M = M.copy()
        M[:K.shape[0], :K.shape[1]] += K
    return M


def kernel(M):
    _, sv, vh = np.linalg.svd(M)
    v = vh[-1].conj()
    tail = float(np.sum(np.abs(v[N//2:])**2))
    return v, tail, float(sv[-1])


def base(qv):
    g = np.zeros(5, complex)
    g[0], g[4] = 1, -qv
    return g


rows = []
for qv in (0.5, float(np.exp(-2*np.pi))):
    v, tail, gap = kernel(section(base(qv)))
    Z0 = abs(v[0])**2
    for trial in range(12):
        g = base(qv)
        eps = np.zeros(7, complex)
        eps[1:] = 0.04*(rng.normal(size=6) + 1j*rng.normal(size=6))
        g2 = np.concatenate([g, np.zeros(2)]) + eps
        K = 0.05*(rng.normal(size=(6, 6)) + 1j*rng.normal(size=(6, 6)))/np.sqrt(6)
        M = section(g2, K=K)
        vv, tt, gp = kernel(M)
        L = M.conj().T @ M
        x = 1e-9
        Cx = np.linalg.solve(x*np.eye(N) + L, np.eye(N)[:, 0])[0].real
        roots = np.roots(g2[::-1])                         # Nullstellen von g2(z): alle außerhalb der Scheibe -> Windung -1
        rows.append((qv, tt, gp, abs(vv[0])**2, x*Cx, float(np.min(np.abs(roots)))))
    say(f"   q = {qv:.4g}: ungestört Z = {Z0:.8f} (1 - q^2 = {1 - qv**2:.8f}), Randgewicht des Nullvektors {tail:.1e}, "
        f"kleinster Singulärwert {gap:.3f}")
zs = {qv: [r[3] for r in rows if r[0] == qv] for qv in (0.5, float(np.exp(-2*np.pi)))}
say(f"   24 Störungen (Symbol bis z^6 und Rang-6-Block): größtes Randgewicht {max(r[1] for r in rows):.1e}, "
    f"kleinster Abstand {min(r[2] for r in rows):.3f}, kleinster Nullstellenbetrag des Symbols {min(r[5] for r in rows):.3f}")
say(f"   Polgewicht Z unter Störung: q = 1/2 von {min(zs[0.5]):.4f} bis {max(zs[0.5]):.4f}; "
    f"q = e^-2pi von {min(zs[float(np.exp(-2*np.pi))]):.4f} bis {max(zs[float(np.exp(-2*np.pi))]):.4f}")
check("ungestört ist Z = 1 - q^2 (Polgewicht des Randzustands im Nullvektor)",
      all(abs(abs(kernel(section(base(qv)))[0][0])**2 - (1 - qv**2)) < 1e-12 for qv in (0.5, np.exp(-2*np.pi))))
check("alle gestörten Symbole haben ihre Nullstellen außerhalb der Einheitsscheibe (Windung -1, Index +1 auch mit dem "
      "kompakten Block)", all(r[5] > 1 for r in rows))
check("unter allen 24 Störungen bleibt ein normierbarer Nullvektor (Randgewicht < 1e-6, Gegenprobe ~ 1) mit Abstand "
      "> 0.2 zum übrigen Spektrum", all(r[1] < 1e-6 and r[2] > 0.2 for r in rows))
check("x C(x) -> |<0|eta'>|^2: der masselose Pol der Randantwort bleibt, mit dem gestörten Gewicht",
      all(abs(r[4] - r[3]) < 1e-6 for r in rows))
check("das Gewicht Z ändert sich unter der Störung (Masselosigkeit geschützt, Kopplung nicht)",
      max(zs[0.5]) - min(zs[0.5]) > 1e-2)
v0, tail0, _ = kernel(section(base(0.5), shift=0))
check("Gegenprobe Windung 0 (ohne S^+): der 'Nullvektor' sitzt am Abschneiderand (Randgewicht ~ 1), kein Pol",
      tail0 > 0.5)
befund("Der masselose Spin-2-Pol ist topologisch geschützt (Fredholmindex +1), sein Residuum und damit G nicht", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Universalität (Weinbergs Softbedingung)")
rows2 = []
for _ in range(30):
    sp_ = rng.integers(0, 3, size=4)
    p = rng.normal(size=(3, 4))
    p = np.vstack([p, -p.sum(axis=0)])                    # alle einlaufend, Impulserhaltung
    A = np.zeros((4, 3))
    for i in range(4):
        A[:, sp_[i]] += p[i]
    rows2.append(A)
A = np.vstack(rows2)
_, sv, vh = np.linalg.svd(A)
null = vh[-1]/vh[-1][0]
say(f"   Bedingung Sum_i kappa_(Art i) p_i = 0 für 30 Prozesse mit drei Arten: Singulärwerte {np.round(sv, 4)}")
check("die Softbedingung hat genau eine Lösung, kappa proportional zu (1, 1, 1): universelle Kopplung",
      sv[-1] < 1e-12 and sv[-2] > 1e-3 and np.allclose(null, 1, atol=1e-10))
qo = np.exp(-2*np.pi)
kap = [np.sqrt(1 - qo**2)*qo**j for j in range(3)]
say(f"   Kopplungen der Leiterquellen |4j> an den Pol (q = e^-2pi): {[f'{k:.3e}' for k in kap]}")
check("die Leiterquellen koppeln nicht universell (kappa_1/kappa_0 = e^-2pi = 0.0019): Materie auf tieferen Leiterstufen "
      "verletzte das Äquivalenzprinzip in Ordnung eins", abs(kap[1]/kap[0] - qo) < 1e-15)
befund("Universell ist nur eine einzige erhaltene Quelle: alle Materie koppelt über den Randzustand |0> (die Naht); "
       "unter Störungen (Teil 1) renormiert sich dann nur G, universell", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Ward plus Positivität -- welcher Kanal darf das Gedächtnis tragen?")
eta = np.diag([-1.0, 1, 1, 1])
idx = list(combinations_with_replacement(range(4), 2))


def sym(vec):
    T = np.zeros((4, 4))
    for (a, b), c in zip(idx, vec):
        T[a, b] = T[b, a] = c
    return T


def conserved_basis(kup):
    rows = []
    for nu in range(4):
        rows.append([sym(e)[:, nu] @ kup for e in np.eye(10)])
    _, s, vh = np.linalg.svd(np.array(rows))
    return [sym(v) for v in vh[np.sum(s > 1e-12):]]


def form(basis, c_tt, c_tr):
    TT = lambda A, B: np.einsum("ab,cd,ac,bd->", A, B, eta, eta)
    tr = lambda A: np.einsum("ab,ab->", eta, A)
    G = np.array([[c_tt*TT(A, B) + c_tr*tr(A)*tr(B) for B in basis] for A in basis])
    ev = np.linalg.eigvalsh((G + G.T)/2)
    return (int(np.sum(ev > 1e-10)), int(np.sum(ev < -1e-10)), int(np.sum(abs(ev) <= 1e-10))), ev


massive = conserved_basis(np.array([1.0, 0, 0, 0]))
massless = conserved_basis(np.array([1.0, 0, 0, 1]))
gr0, _ = form(massless, 1, -0.5)
fp, _ = form(massive, 1, -1/3)
same, ev_same = form(massive, 1, -0.5)
say(f"   Trägheit (+, -, 0): masselos TT - T^2/2 {gr0}; massiv TT - T^2/3 {fp}; massiv TT - T^2/2 {same}")
check("masselos: TT - T^2/2 ist positiv vom Rang 2 (zwei Helizitäten), der Rest Eichrichtungen", gr0 == (2, 0, 4))
check("massiv: TT - T^2/3 (Fierz-Pauli) ist positiv vom Rang 5", fp == (5, 0, 1))
check("massiv mit gleichem Formfaktor in Spin 2 und Spin 0 (TT - T^2/2): ein negativer Eigenwert -1/2 auf der Spur "
      "(Skalargeist)", same == (5, 1, 0) and abs(min(ev_same) + 0.5) < 1e-12)
ok = []
for rho0 in (-0.3, 0.0, 0.2):
    inert, _ = form(massive, 1.0, -1/3 - rho0/6)
    ok.append((rho0, inert))
say(f"   massiver Pol mit Spin-2-Gewicht 1 und Spin-0-Gewicht rho0: {ok}")
check("Residuum rho2 (TT - T^2/3) - (rho0/6) T^2 ist genau dann positiv, wenn rho0 <= 0", [o[1][1] for o in ok] == [0, 0, 1])
dz = []
for z0 in (0.8, 1.0, 1.2):
    inert, _ = form(massless, 1.0, -1/3 - z0/6)
    dz.append((z0, inert))
say(f"   masseloser Pol mit Z2 = 1 und Spur-Residuum Z0: {dz}")
check("Z0 = Z2 gibt den reinen Graviton (Rang 2); Z0 < Z2 einen zusätzlichen masselosen Skalar, Z0 > Z2 einen Geist",
      dz[1][1] == (2, 0, 4) and dz[0][1][0] == 3 and dz[2][1][1] == 1)
q2 = np.exp(-4*np.pi)
say(f"   im originalen Zustand (z46): Bandgewicht {q2:.3e}; mit gleichem Formfaktor in beiden Kanälen wäre das ein "
    f"Skalargeist der Stärke {q2/2:.2e}")
befund("Ward plus Positivität legen die Kanalstruktur fest: das Übertragungsgedächtnis nur im Spin-2-Kanal, der "
       "Spurkanal reiner Einsteinpol mit demselben Residuum 1 - q^2; gleiche Formfaktoren in beiden Kanälen sind "
       "ausgeschlossen, sobald massive Pole auftreten", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Kovariante Form -- Weyl-Quadrat trägt nur Spin 2")
d = np.eye(4)


def lin_weyl(h, k):
    R = -0.5*(np.einsum("n,r,ms->mnrs", k, k, h) + np.einsum("m,s,nr->mnrs", k, k, h)
              - np.einsum("m,r,ns->mnrs", k, k, h) - np.einsum("n,s,mr->mnrs", k, k, h))
    Ric = np.einsum("mnms->ns", R)
    Rs = np.trace(Ric)
    g = d
    C = (R - 0.5*(np.einsum("mr,ns->mnrs", g, Ric) - np.einsum("ms,nr->mnrs", g, Ric)
                  - np.einsum("nr,ms->mnrs", g, Ric) + np.einsum("ns,mr->mnrs", g, Ric))
         + Rs/6*(np.einsum("mr,ns->mnrs", g, g) - np.einsum("ms,nr->mnrs", g, g)))
    return C, Rs


def p2(h, k):
    th = d - np.outer(k, k)/(k @ k)
    hp = th @ h @ th
    return np.sum(hp*hp) - np.trace(hp)**2/3


ratios, rs = [], []
for _ in range(8):
    k = rng.normal(size=4)
    h = rng.normal(size=(4, 4))
    h = h + h.T
    C, Rs = lin_weyl(h, k)
    ratios.append(np.sum(C*C)/((k @ k)**2*p2(h, k)))
    th = d - np.outer(k, k)/(k @ k)
    rs.append(Rs**2/((k @ k)**2*np.trace(th @ h)**2))
k = rng.normal(size=4)
xi = rng.normal(size=4)
Cg, _ = lin_weyl(np.outer(k, xi) + np.outer(xi, k), k)
Ct, _ = lin_weyl(np.eye(4), k)
say(f"   W.W/(k^4 h P2 h) über 8 Zufallsproben: {np.round(ratios, 12)}; R^2/(k^4 (theta.h)^2): {np.round(rs, 12)}")
check("linearisiertes W.W = (1/2) k^4 h P2 h (euklidisch, konstant über Proben); Eichmoden und reine Spur geben null",
      np.allclose(ratios, 0.5, atol=1e-10) and np.max(abs(Cg)) < 1e-12 and np.max(abs(Ct)) < 1e-12)
check("R^2 = k^4 (theta.h)^2 trägt nur den Spinnull-Kanal", np.allclose(rs, 1, atol=1e-10))
th = d - np.outer(k, k)/(k @ k)
P2op = lambda X: 0.5*(th @ X @ th + th @ X.T @ th) - th*np.trace(th @ X)/3
P0op = lambda X: th*np.trace(th @ X)/3
hg = np.outer(k, xi) + np.outer(xi, k)
check("P2 und P0 vernichten reine Eichmoden h = k xi + xi k: jeder Formfaktor F2 P2 + F0 P0 ist linear eichinvariant",
      np.max(abs(P2op(hg))) < 1e-12 and np.max(abs(P0op(hg))) < 1e-12)
befund("Kovariante Realisierung: S = (1/16 pi G) Int sqrt g [R + W F(box) W] mit F aus dem Übertragungsgedächtnis trägt "
       "die Spinzwei-Erinnerung und lässt den Spurkanal Einstein; die nichtlineare Vervollständigung der massiven Linie "
       "(bimetrische Wechselwirkungen) ist nicht festgelegt", True)

urteil = ("Die Ward-Struktur ist dreiteilig. Erstens ist der masselose Pol topologisch geschützt: Er ist der Nullvektor "
          "eines Fredholmoperators mit Index +1 und überlebt Symbol- und endlichrangige Störungen; sein Gewicht Z und "
          "damit G ändern sich. Zweitens verlangt Weinbergs Softbedingung gleiche Kopplungen; die Leiterquellen koppeln "
          "mit sqrt(Z) q^j, universell ist nur eine einzige Quelle, der Randzustand (Materie auf der Naht). Drittens "
          "erlauben Ward-Identität und Positivität das Gedächtnis nur im Spin-2-Kanal: ein massiver Pol mit gleichem "
          "Formfaktor in Spin 2 und Spin 0 trägt einen Skalargeist (-1/2 auf der Spur), und der Spurkanal muss reiner "
          "Einsteinpol mit demselben Residuum sein. Kovariant leistet das ein Term W F(box) W, denn das linearisierte "
          "Weyl-Quadrat ist (1/2) k^4 h P2 h. Offen bleiben die nichtlineare Vervollständigung der massiven Linie und "
          "warum alle Materie auf dem Randzustand sitzt.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
ZAHLEN["z_stoerung_q_halb"] = [min(zs[0.5]), max(zs[0.5])]
(HERE / "z47_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
