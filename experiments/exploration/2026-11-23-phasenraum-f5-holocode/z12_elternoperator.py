"""Z12: Ein gemeinsamer Elternoperator für Naht, E8 und Familien (Folgeauftrag 2 der externen Kritik).

Konstruktion. Produktgeometrie (Gitterquerschnitt L x L mit verdrilltem Defekt) x (Impuls k_z) x
(Markensphäre CP^1 mit den vier Marken D = mu4). Ein einziger BdG-Operator

    H_P(k_z) = [[ A(k_z) (x) Gamma + 1 (x) D_int ,  B(k_z) (x) C_int          ],
                [ (B (x) C_int)^dag              , -A(-k_z)^T (x) Gamma - 1 (x) D_int ]]

mit A, B aus dem Z7-Defektmodell (Weyl-Kinetik, windende Paarungsmasse, Viertelfluss, Wilson-Glied),
D_int dem Spin-Dirac auf der Sphäre, gekoppelt an L = K^{1/2} (x) O(D) (Monopolharmonische, Grad
deg L = -1 + deg D = 3), Gamma der inneren Chiralität und C_int einer symmetrischen, mit D_int und
Gamma vertauschenden Paarungsform. Das ist die Form eines 6D-Weyl-Felds auf R^3 x S^2 (A (x) Gamma +
1 (x) D_int antikommutiert blockweise). Die Vierteldrehung z -> i z ist R_int = lambda_D exp(i pi J_z/2);
lambda_D ist der Lift, der den kanonischen Markenschnitt s_D = u^4 - v^4 festhält.

Geprüft in einem Lauf: (i) chirale Majorana-Weyl Moden am Defekt, (ii) Twistgewicht und GSO des
entstehenden Nahtsektors, (iii) ker D_+ = 3, ker D_- = 0 mit Charakteren {-1, i, -i}; dazu die
Vierteldrehung auf Familien und Nahtmoden, J^2 = -1, und der Einwand aus der Gesamtsynthese
("Die vier ursprünglichen Marken wählen eine echte lokale Zweiteilung").
"""
import itertools
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import sympy as sp
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).parent
_src = (HERE / "z7_ursprung_naht.py").read_text()
_ns = {"__file__": str(HERE / "z7_ursprung_naht.py")}
exec(_src[:_src.index("L = 14")], _ns)
build = _ns["build"]

LINES, CHECKS, BEFUNDE, DATA = [], {}, {}, {}


def say(s=""):
    print(s)
    LINES.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    return "ja" if cond else "NEIN"


def befund(name, val):
    BEFUNDE[name] = bool(val)
    return val


# ---------------------------------------------------------------------------------------------------
# 1. Innerer Operator: Spin-Dirac auf der Markensphäre, gekoppelt an K^{1/2} (x) O(D)
# ---------------------------------------------------------------------------------------------------
def internal(deg, lmax, Rint=2.0):
    """Monopolharmonische Basis. deg = Grad von K^{1/2} (x) L_twist = 2q. Komponenten u (Spingewicht
    -q-1/2, Gamma = -1) und v (Spingewicht -q+1/2, Gamma = +1). D_int = [[0, A], [A^T, 0]] mit dem
    Absteiger eth-bar: v_{lm} -> u_{lm}, Koeffizient sqrt((l+s_v)(l-s_v+1))/R. Reell symmetrisch."""
    q = sp.Rational(deg, 2)
    s_u, s_v = -q - sp.Rational(1, 2), -q + sp.Rational(1, 2)
    states = []                                     # (komponente, l, m)
    for comp, s in (("v", s_v), ("u", s_u)):
        l = abs(s)
        while l <= lmax:
            for k in range(int(2 * l) + 1):
                states.append((comp, l, -l + k))
            l += 1
    n = len(states)
    pos = {st: i for i, st in enumerate(states)}
    D = np.zeros((n, n))
    for (comp, l, m), i in pos.items():
        if comp == "v" and ("u", l, m) in pos:
            c = float(sp.sqrt((l + s_v) * (l - s_v + 1))) / Rint
            j = pos[("u", l, m)]
            D[i, j] = D[j, i] = c
    Gam = np.diag([1.0 if st[0] == "v" else -1.0 for st in states])
    C = np.zeros((n, n))
    for (comp, l, m), i in pos.items():
        C[i, pos[(comp, l, -m)]] = 1.0             # symmetrische Paarung m <-> -m, mit D und Gamma vertauschend
    Jz = np.array([float(st[2]) for st in states])
    return dict(states=states, D=D, Gam=Gam, C=C, Jz=Jz, n=n, q=q)


def mark_lift(marks, inf_mult=0, trials=5, seed=1):
    """SU(2)-Bild von z -> i z: (u, v) -> (e^{i pi/4} u, e^{-i pi/4} v). Eigenwert des homogenen
    Markenschnitts S_D. lambda_D = 1/Eigenwert hält s_D fest. None, wenn S_D kein Eigenvektor ist."""
    rng = np.random.default_rng(seed)
    g = np.exp(1j * np.pi / 4)

    def S(u, v):
        out = v ** inf_mult
        for p in marks:
            out = out * (u - p * v)
        return out

    rats = []
    for _ in range(trials):
        u, v = rng.normal(size=2) + 1j * rng.normal(size=2)
        rats.append(S(g * u, v / g) / S(u, v))
    rats = np.array(rats)
    if np.max(np.abs(rats - rats[0])) > 1e-9:
        return None
    return 1.0 / rats[0]


def chars(lam, I):
    zero = [i for i, st in enumerate(I["states"]) if st[0] == "v" and abs(I["D"][i]).sum() < 1e-12]
    return [lam * np.exp(1j * np.pi * I["Jz"][i] / 2) for i in zero], zero


def fmt(zs):
    out = []
    for z in sorted(zs, key=lambda x: (round(np.angle(x), 6))):
        z = complex(np.round(z, 9))
        lab = {1: "1", -1: "-1", 1j: "i", -1j: "-i"}.get(z, f"{z:.3f}")
        out.append(lab)
    return "{" + ", ".join(out) + "}"


def same_set(a, b):
    a = sorted(np.round(np.array(a, complex), 9), key=lambda x: (x.real, x.imag))
    b = sorted(np.round(np.array(b, complex), 9), key=lambda x: (x.real, x.imag))
    return len(a) == len(b) and np.allclose(a, b)


say("1. Innerer Operator auf der Markensphäre\n")
degD = 4
degL = -1 + degD
say(f"   Grad: deg K^(1/2) + deg O(D) = -1 + {degD} = {degL}. Monopolladung q = {degL}/2.")
I3 = internal(degL, lmax=2)
ev = np.linalg.eigvalsh(I3["D"])
nz = int(np.sum(np.abs(ev) < 1e-10))
kerp = sum(1 for i, st in enumerate(I3["states"]) if st[0] == "v" and abs(I3["D"][i]).sum() < 1e-12)
kerm = sum(1 for i, st in enumerate(I3["states"]) if st[0] == "u" and abs(I3["D"][i]).sum() < 1e-12)
gap = np.min(np.abs(ev[np.abs(ev) > 1e-10]))
say(f"   D_int ({I3['n']} Zustände bis l = 2): {nz} Nullmoden, ker D+ = {kerp}, ker D- = {kerm}, nächstes |lambda| = {gap:.3f}")
say(f"   {{Gamma, D_int}} = 0: {np.allclose(I3['Gam'] @ I3['D'] + I3['D'] @ I3['Gam'], 0)}; [C, D] = 0: "
    f"{np.allclose(I3['C'] @ I3['D'], I3['D'] @ I3['C'])}")
check("Index D_int = 3, ker D- = 0", nz == 3 and kerp == 3 and kerm == 0)

say("\n   Vierteldrehung z -> i z. Der Lift ist durch den Markenschnitt festgelegt (s_D bleibt fest).")
mu4 = [1, 1j, -1, -1j]
lam_mu4 = mark_lift(mu4)
ch_mu4, zero_idx = chars(lam_mu4, I3)
say(f"   D = mu4: SU(2)-Eigenwert von s_D = {complex(np.round(1 / lam_mu4, 9))}, Lift lambda_D = "
    f"{complex(np.round(lam_mu4, 9))}, Familiencharaktere {fmt(ch_mu4)}")
z = sp.symbols("z")
MK = sp.Matrix([[sp.Poly(sp.expand(sp.I * b.subs(z, sp.I * z)), z).coeff_monomial(c) for b in (sp.Integer(1), z, z ** 2)]
                for c in (sp.Integer(1), z, z ** 2)])
ev_sym = [complex(e) for e in MK.eigenvals(multiple=True)]
say(f"   Gegenrechnung wie Z10 (omega = P dz/(z^4 - 1) -> i P(i z)): {fmt(ev_sym)}")
check("Familien {-1, i, -i} aus dem Markenlift", same_set(ch_mu4, [-1, 1j, -1j]) and same_set(ch_mu4, ev_sym))

say("\n   Gegenproben zur Lage der Marken (Grad 4, invariant unter z -> i z):")
scan = {}
for name, mk, inf in (("freie Bahn mu4", mu4, 0),
                      ("freie Bahn e^(i pi/4) mu4", list(np.exp(1j * np.pi / 4) * np.array(mu4)), 0),
                      ("freie Bahn 2 mu4", list(2 * np.array(mu4)), 0),
                      ("2[0] + 2[oo]", [0, 0], 2), ("3[0] + [oo]", [0, 0, 0], 1), ("4[0]", [0, 0, 0, 0], 0),
                      ("nicht invariant {1, 2, -1, -2}", [1, 2, -1, -2], 0)):
    lam = mark_lift(mk, inf)
    if lam is None:
        scan[name] = None
        say(f"   {name:32s}: keine Liftung, z -> i z ist keine Symmetrie")
        continue
    c, _ = chars(lam, I3)
    scan[name] = fmt(c)
    distinct = len(set(np.round(mk, 9))) + (1 if inf else 0) == 4
    say(f"   {name:32s}: lambda_D = {complex(np.round(lam, 9))}, Charaktere {fmt(c)}"
        + ("  (vier verschiedene Marken)" if distinct else ""))
ok_scan = (scan["freie Bahn e^(i pi/4) mu4"] == scan["freie Bahn mu4"] == scan["freie Bahn 2 mu4"]
           and scan["2[0] + 2[oo]"] != scan["freie Bahn mu4"] and scan["3[0] + [oo]"] != scan["freie Bahn mu4"]
           and scan["nicht invariant {1, 2, -1, -2}"] is None)
check("Charaktere hängen an der Markenlage (freie Bahn nötig)", ok_scan)
say("   Vier verschiedene invariante Marken bilden zwangsläufig eine freie Bahn; dann ist der Eigenwert von")
say("   s_D immer -1 und die Charaktere sind {-1, i, -i}. Nicht reduzierte Divisoren an den Fixpunkten geben")
say("   andere Werte (4[0] fällt zufällig zusammen). Hier gehen die Marken wirklich in den Operator ein.")
lam_flip = -lam_mu4
say(f"   Eingebaut: Spinor und Kopplungsbündel benutzen dieselbe Wurzel K^(1/2). Mit entgegengesetztem")
say(f"   Spinlift der Kopplungswurzel wären es {fmt(chars(lam_flip, I3)[0])}.")
DATA["family_chars"] = fmt(ch_mu4)
DATA["mark_scan"] = scan


# ---------------------------------------------------------------------------------------------------
# 2. Der Elternoperator auf Querschnitt x k_z x Markensphäre
# ---------------------------------------------------------------------------------------------------
_cache = {}


def perp_blocks(L, kz, copies, npl=1, nmi=0, fl=0.25):
    key = (L, round(kz, 12), copies, npl, nmi, fl)
    if key not in _cache:
        H, (rad, ns) = build(L, kz, npl, nmi, fl, copies=copies)
        d = H.shape[0] // 2
        _cache[key] = (sps.csr_matrix(H[:d, :d]), sps.csr_matrix(H[:d, d:]), rad, ns, d)
    return _cache[key]


def parent(L, kz, I, copies=1, pair=True):
    A, B, rad, ns, d = perp_blocks(L, kz, copies)
    Am = perp_blocks(L, -kz, copies)[0]
    Id = sps.identity(d, format="csr")
    G, Dn, C = sps.csr_matrix(I["Gam"]), sps.csr_matrix(I["D"]), sps.csr_matrix(I["C"])
    TL = sps.kron(A, G) + sps.kron(Id, Dn)
    BR = -sps.kron(Am.T, G) - sps.kron(Id, Dn.T)
    TR = sps.kron(B, C) if pair else sps.csr_matrix((TL.shape[0], TL.shape[0]))
    return sps.bmat([[TL, TR], [TR.conj().T, BR]], format="csc"), rad, ns, d


def quarter(I, lam, d):
    R = np.diag(lam * np.exp(1j * np.pi * I["Jz"] / 2))
    Id = sps.identity(d, format="csr")
    return sps.block_diag([sps.kron(Id, sps.csr_matrix(R)), sps.kron(Id, sps.csr_matrix(R.conj()))], format="csr")


def seam(L, I, copies=1, lam=1.0, Rc=3.5, window=0.15, K=40, dk=1e-3, sigma=0.0137):
    H0, rad, ns, d = parent(L, 0.0, I, copies)
    n_int = I["n"]
    vals, vecs = eigsh(H0, k=K, sigma=sigma, which="LM")
    reach = np.max(np.abs(vals - sigma))
    Q, _ = np.linalg.qr(vecs)                        # ARPACK-Vektoren in Clustern nachorthonormieren
    e, w = np.linalg.eigh(Q.conj().T @ (H0 @ Q))     # Rayleigh-Ritz im gefundenen Unterraum
    S = (Q @ w)[:, np.abs(e) < window]
    core_perp = np.repeat(rad <= Rc, ns * 2)
    mask = np.tile(np.repeat(core_perp, n_int), 2).astype(float)
    cv, cvec = np.linalg.eigh(S.conj().T @ (mask[:, None] * S))
    Vc, Ve = S @ cvec[:, cv > 0.7], S @ cvec[:, cv < 0.3]
    dH = (parent(L, dk, I, copies)[0] - parent(L, -dk, I, copies)[0]) / (2 * dk)
    vc = np.linalg.eigvalsh(Vc.conj().T @ (dH @ Vc)) if Vc.shape[1] else np.array([])
    ve = np.linalg.eigvalsh(Ve.conj().T @ (dH @ Ve)) if Ve.shape[1] else np.array([])
    zero_int = np.array([abs(I["D"][i]).sum() < 1e-12 for i in range(n_int)])
    m0 = np.tile(np.tile(zero_int, d), 2).astype(float)
    w0 = np.linalg.eigvalsh(Vc.conj().T @ (m0[:, None] * Vc)) if Vc.shape[1] else np.array([])
    RP = quarter(I, lam, d)
    comm = sps.linalg.norm(RP @ H0 - H0 @ RP)
    if Vc.shape[1]:
        M = Vc.conj().T @ (RP @ Vc)
        inv = np.linalg.norm(RP @ Vc - Vc @ M)
        rch = np.linalg.eigvals(M)
    else:
        inv, rch = 0.0, np.array([])
    net_c = int(np.sum(vc > 0.2) - np.sum(vc < -0.2))
    net_e = int(np.sum(ve > 0.2) - np.sum(ve < -0.2))
    return dict(net_core=net_c, net_edge=net_e, vc=vc, w0=w0, comm=comm, inv=inv, rch=rch, reach=reach,
                window_ok=reach > window + sigma, dim=H0.shape[0], nlow=S.shape[1])


say("\n2. Der Elternoperator: Querschnitt mit verdrilltem Defekt x k_z x Markensphäre\n")
Lc = 12
say(f"   Querschnitt {Lc} x {Lc}, Kernradius 3,5; Defekt wie Z7 (n+ = 1, n- = 0, Viertelfluss), m = 0,6, r_w = 1.")
say("   Ein Operator H_P auf dem ganzen Produkt; gezählt wird bei k_z = 0 mit Kernprojektor und <dH/dk_z>.")
H_test, _, _, dtest = parent(Lc, 0.3, I3)
Hm_test = parent(Lc, -0.3, I3)[0]
nn = H_test.shape[0] // 2
tx = sps.bmat([[None, sps.identity(nn)], [sps.identity(nn), None]], format="csr")
herm = sps.linalg.norm(H_test - H_test.conj().T)
phs = sps.linalg.norm(tx @ H_test.conj() @ tx + Hm_test)
say(f"   Hermitesch: {herm:.1e}; Teilchen-Loch-Symmetrie tau_x H(k)^* tau_x = -H(-k): {phs:.1e}")
check("H_P hermitesch und teilchen-loch-symmetrisch", herm < 1e-10 and phs < 1e-10)

runs = {}
for label, deg, lmax, cp, lam in (("Grad 3 (Familienbündel), ein Paar", 3, 2, 1, lam_mu4),
                                  ("Grad 1, ein Paar", 1, 1, 1, 1.0),
                                  ("Grad 1, zwei Paare", 1, 1, 2, 1.0),
                                  ("Grad 0 (ohne Monopol), ein Paar", 0, 1.5, 1, 1.0),
                                  ("Grad -1, ein Paar", -1, 1, 1, 1.0)):
    I = I3 if deg == 3 else internal(deg, lmax)
    ind = int(round(sum(I["Gam"][i, i] for i, st in enumerate(I["states"]) if abs(I["D"][i]).sum() < 1e-12)))
    r = seam(Lc, I, copies=cp, lam=lam)
    r["ind"], r["copies"] = ind, cp
    runs[label] = r
    say(f"   {label:36s}: ind D_int = {ind:+d}, Dimension {r['dim']}, netto Kern {r['net_core']:+d}, Rand "
        f"{r['net_edge']:+d}, v_Kern = {np.round(r['vc'], 3)}")
    say(f"   {'':36s}  Gewicht der Kernmoden im inneren Nullraum {np.round(r['w0'], 4)}; [R_P, H_P] = {r['comm']:.1e}")
mult = all(r["net_core"] == r["copies"] * r["ind"] and r["net_edge"] == -r["copies"] * r["ind"] for r in runs.values())
check("Naht im Produkt: netto Kernmoden = Paare x ind(D_int)", mult and all(r["window_ok"] for r in runs.values()))
check("Kernmoden mit Lichtgeschwindigkeit", all(np.all(np.abs(np.abs(r["vc"]) - 1) < 0.03) for r in runs.values()))
check("KK-Niveaus tragen keine Nahtmoden", all(np.all(r["w0"] > 0.99) for r in runs.values() if len(r["w0"])))
check("Vierteldrehung ist Symmetrie von H_P", all(r["comm"] < 1e-9 for r in runs.values()))
r3 = runs["Grad 3 (Familienbündel), ein Paar"]
say(f"\n   Vierteldrehung R_P auf den drei Nahtmoden (Grad 3): {fmt(r3['rch'])}, Invarianzfehler {r3['inv']:.1e}")
check("R_P auf den Nahtmoden = Familiencharaktere", same_set(r3["rch"], [-1, 1j, -1j]) and r3["inv"] < 1e-6)
say(f"   R_P^2 auf den Nahtmoden: {fmt(r3['rch'] ** 2)}. Fermionparität wäre -1 auf allen drei.")
befund("R_P^2 = (-1)^F auf den Nahtmoden", np.allclose(r3["rch"] ** 2, -1))
say("   Die Nahtmoden eines Multipletts mit Familienbündel sind Naht (x) Familie: Jede Familie trägt ihren")
say("   eigenen Satz chiraler Moden, und die Vierteldrehung wirkt auf ihnen mit den Familiencharakteren.")
DATA["seam_runs"] = {k: {kk: (vv if not isinstance(vv, np.ndarray) else [str(np.round(x, 4)) for x in vv])
                         for kk, vv in v.items() if kk not in ("comm", "inv", "reach")} for k, v in runs.items()}


# ---------------------------------------------------------------------------------------------------
# 3. J als Monodromie der Nahtfermionen, und ob eine einzige Vierteldrehung beides sein kann
# ---------------------------------------------------------------------------------------------------
say("\n3. J aus dem Operator und die Frage nach einer einzigen Vierteldrehung\n")
A0, _, rad0, ns0, d0 = perp_blocks(Lc, 0.0, 1)
A0 = A0.toarray()
sxm, sym_ = np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]])
ix = lambda x, y, a, s: ((x * Lc + y) * ns0 + a) * 2 + s


def loop_phase(a, x0, x1):
    tot = 0.0
    path = [((x, x0), (x + 1, x0), sxm, +1) for x in range(x0, x1)] + \
           [((x1, y), (x1, y + 1), sym_, +1) for y in range(x0, x1)] + \
           [((x, x1), (x + 1, x1), sxm, -1) for x in range(x0, x1)] + \
           [((x0, y), (x0, y + 1), sym_, -1) for y in range(x0, x1)]
    for (r, r2, sig, sgn) in path:
        ph = A0[ix(*r, a, 0), ix(*r2, a, 1)] / (sig[0, 1] / 2j)
        tot += sgn * np.angle(ph)
    return np.exp(1j * tot)


Jm = {}
for a, nm in ((0, "eta+"), (1, "eta-")):
    w1, w2 = loop_phase(a, 3, 8), loop_phase(a, 1, 10)
    Jm[nm] = w2
    say(f"   Monodromie um den Defekt für {nm}: Schleife 5x5 {complex(np.round(w1, 9))}, 9x9 {complex(np.round(w2, 9))}")
Jd = np.diag([Jm["eta+"], Jm["eta-"]])
say(f"   J = diag{tuple(complex(np.round(x, 9)) for x in np.diag(Jd))}, J^2 = {np.round(np.diag(Jd @ Jd).real, 9)}"
    " (eingebaut über den Viertelfluss).")
check("J^2 = -1 als Monodromie, unabhängig von der Schleife",
      np.allclose(Jd @ Jd, -np.eye(2)) and abs(loop_phase(0, 3, 8) - loop_phase(0, 1, 10)) < 1e-9)
rk = np.array([-1, 1j, -1j])
us = np.exp(1j * np.linspace(0, 2 * np.pi, 3601))
best = min(np.max(np.abs((u * rk) ** 2 + 1)) for u in us)
say(f"   Jede Produktsymmetrie U = U_perp (x) U_int wirkt auf Naht (x) Familie k mit u r_k. Verlangt man")
say(f"   U^2 = -1 auf allen Nahtmoden und U = {{-1, i, -i}} auf den Familien: kleinster Fehler über alle u {best:.3f}.")
befund("eine Vierteldrehung mit {-1,i,-i} auf Familien und J^2=-1 auf allen Nahtmoden", best < 1e-6)
detc = np.prod(rk)
say(f"   Grund: r_k^2 = {fmt(rk ** 2)}. Die Charaktere enthalten -1, und (-1)^2 = +1. Ein Spin-Z4 mit J^2 = (-1)^F")
say(f"   hat nur die Eigenwerte +-i. det = {complex(detc)}: kein Element von SU(3); (-1)^F mal die Charaktere gibt")
say(f"   {fmt(-rk)} mit det {complex(np.prod(-rk))}, also diag(1, -i, i) in SU(3).")
befund("Familiencharaktere als Element von SU(3) (det 1)", abs(detc - 1) < 1e-12)

# ---------------------------------------------------------------------------------------------------
# 4. Volumen: Wer bleibt leicht? Die Paarung, die die Naht macht, lückt die Felder, die sie trägt.
# ---------------------------------------------------------------------------------------------------
say("\n4. Volumen ohne Defekt: Lücke der Familienmoden mit und ohne Paarung\n")


def bulk(k, I, pair=True, m=0.6, rw=1.0):
    hs = lambda kk: np.sin(kk[0]) * sxm + np.sin(kk[1]) * sym_ + np.sin(kk[2]) * np.diag([1, -1]).astype(complex)
    W = rw * (3 - np.cos(k[0]) - np.cos(k[1]) - np.cos(k[2]))
    hp, hm = np.kron(np.eye(2), hs(k)), np.kron(np.eye(2), hs(-np.array(k)))
    Dl = np.kron(np.array([[m, W], [W, m]]), 1j * sym_) if pair else np.zeros((4, 4))
    TL = np.kron(hp, I["Gam"]) + np.kron(np.eye(4), I["D"])
    BR = -np.kron(hm.T, I["Gam"]) - np.kron(np.eye(4), I["D"])
    TR = np.kron(Dl, I["C"])
    return np.linalg.eigvalsh(np.block([[TL, TR], [TR.conj().T, BR]]))


e_pair = bulk([0, 0, 0], I3, True)
e_free = bulk([0, 0, 0], I3, False)
e_free_k = bulk([0.02, 0, 0], I3, False)
gap_pair = np.min(np.abs(e_pair))
nzero = int(np.sum(np.abs(e_free) < 1e-9))
lin = np.sort(np.abs(e_free_k))[:nzero]
say(f"   Mit Paarung (Grad 3): kleinste |E| bei k = 0 ist {gap_pair:.3f} (Masse m = 0,6). Die drei Familien sind schwer.")
say(f"   Ohne Paarung: {nzero} Nullzustände bei k = 0 (3 Familien x 2 Arten x 2 Spin x BdG 2); bei k = 0,02 |E| = "
    f"{np.round(np.unique(np.round(lin, 6)), 4)}")
check("Paarung lückt die Familien, ungepaarte Materie hat 3 Weyl-Kegel",
      abs(gap_pair - 0.6) < 1e-6 and nzero == 24 and np.allclose(lin, np.sin(0.02), atol=1e-6))
befund("Familien bleiben leicht, wenn dasselbe Multiplett die Naht trägt", gap_pair < 1e-6)
say("   Die windende Paarungsmasse ist eine Volumenmasse. Ein Multiplett, das die Naht trägt, ist im Volumen")
say("   schwer. Leichte chirale Familien brauchen ein zweites Multiplett ohne diese Paarung.")

# ---------------------------------------------------------------------------------------------------
# 4b. Zwei Multipletts: Naht am halben Markenbündel, Materie am Markenbündel
# ---------------------------------------------------------------------------------------------------
say("\n4b. Zweite Variante: Nahtmultiplett am halben Markenbündel M mit M^2 = O(D)\n")
say("   Nahtfelder koppeln an K^(1/2) (x) M (Grad -1 + 2 = 1), Materie an K^(1/2) (x) O(D) (Grad 3). Die")
say("   Nullmode der Naht ist ein Schnitt von K (x) M = (K^2(D))^(1/2): die Wurzel des quadratischen")
say("   Differentials q = dz^2/(z^4 - 1) der flachen zweiblättrigen Sphäre, also omega_T = dz/sqrt(z^4 - 1).")
rng = np.random.default_rng(7)
zs = rng.normal(size=6) * 0.4 + 1j * rng.normal(size=6) * 0.4
rat_q, rat_w = [], []
for z0 in zs:
    y0 = np.sqrt(z0 ** 4 - 1 + 0j)                 # Punkt (z0, y0) der Doppelüberlagerung y^2 = z^4 - 1
    assert abs(y0 ** 2 - ((1j * z0) ** 4 - 1)) < 1e-12   # (i z0, y0) liegt wieder auf der Kurve
    rat_w.append((1j / y0) / (1 / y0))             # Lift (z, y) -> (i z, y): g^* omega_T / omega_T
    rat_q.append(((1j) ** 2 / ((1j * z0) ** 4 - 1)) / (1 / (z0 ** 4 - 1)))
r_s = complex(np.mean(rat_w))
say(f"   Gemessen: g^* q / q = {complex(np.round(np.mean(rat_q), 9))} (q wechselt das Vorzeichen: das Blatt);")
say(f"   Lift (z, y) -> (i z, y): g^* omega_T / omega_T = {complex(np.round(r_s, 9))}; der andere Lift (i z, -y) gibt "
    f"{complex(np.round(-r_s, 9))}.")
say(f"   Gegenrechnung über den Markenlift: lambda_M^2 = lambda_D = {complex(np.round(lam_mu4, 9))}, also lambda_M = +-i.")
check("Nahtnullmode hat Charakter +-i, q wechselt das Vorzeichen",
      np.allclose(rat_q, -1) and np.allclose(rat_w, r_s) and abs(abs(r_s.imag) - 1) < 1e-9
      and np.isclose(r_s ** 2, lam_mu4))
I1 = internal(1, 1)
rs_run = seam(Lc, I1, copies=1, lam=r_s)
H1, _, _, d1 = parent(Lc, 0.0, I1)
RP1 = quarter(I1, r_s, d1)
H1flip = sps.bmat([[H1[:H1.shape[0] // 2, :H1.shape[0] // 2], -H1[:H1.shape[0] // 2, H1.shape[0] // 2:]],
                   [-H1[H1.shape[0] // 2:, :H1.shape[0] // 2], H1[H1.shape[0] // 2:, H1.shape[0] // 2:]]], format="csr")
flip_err = sps.linalg.norm(RP1 @ H1 @ RP1.conj().T - H1flip)
say(f"   Im Operator (Grad 1, ein Paar): netto Kern {rs_run['net_core']:+d}; R_P H_P R_P^dag = H_P mit Paarung -B: "
    f"Fehler {flip_err:.1e}")
z1 = np.array([abs(I1["D"][i]).sum() < 1e-12 for i in range(I1["n"])])
sec = np.tile(np.tile(z1, d1), 2)
R2 = (RP1 @ RP1).diagonal()[sec]
Ug = np.concatenate([np.full(d1 * I1["n"], 1j), np.full(d1 * I1["n"], -1j)])[sec]
say(f"   Auf dem Nullmodensektor der Naht ist R_P die globale Phase exp(i pi N/2): Abweichung "
    f"{np.max(np.abs(RP1.diagonal()[sec] - Ug)):.1e}; R_P^2 = {np.unique(np.round(R2, 9))} = (-1)^F.")
say(f"   Eine Majoranamode des Defekts ist daher kein Eigenvektor (Überlapp {np.round(np.abs(rs_run['rch']), 6)}):")
say("   R_P bildet sie auf die Mode des umgekehrten Blatts ab. Das ist dieselbe Struktur wie bei J.")
check("Zwei Multipletts: eine Vierteldrehung, +-i auf der Naht (R^2 = (-1)^F), {-1, i, -i} auf Familien",
      rs_run["net_core"] == 1 and flip_err < 1e-9 and np.allclose(R2, -1)
      and np.allclose(RP1.diagonal()[sec], Ug) and same_set(ch_mu4, [-1, 1j, -1j]))
say("   Die Vierteldrehung ist dann nur zusammen mit dem Blattwechsel Phi -> -Phi eine Symmetrie, genau wie J.")
say("   Unterschied zu J: R_int multipliziert eta+ und eta- beide mit i und kippt damit auch das Wilson-Glied")
say("   W eta+ eta-; J = diag(-i, i) lässt es fest. Gleich sind beide auf der windenden Paarmasse.")
Bp = perp_blocks(Lc, 0.0, 1)[1].toarray()
spc = (np.arange(Bp.shape[0]) // 2) % ns0
jv, rv = np.diag(Jd)[spc], np.full(Bp.shape[0], r_s)
off = spc[:, None] != spc[None, :]
Wil = np.where(off, Bp, 0)
Win = np.where(~off, Bp, 0)
dJ = (np.linalg.norm(jv[:, None] * Wil * jv[None, :] - Wil), np.linalg.norm(jv[:, None] * Win * jv[None, :] + Win))
dR = (np.linalg.norm(rv[:, None] * Wil * rv[None, :] - Wil), np.linalg.norm(rv[:, None] * Win * rv[None, :] + Win))
say(f"   Gerechnet am Paarungsblock B: J hält das Wilson-Glied fest (Fehler {dJ[0]:.1e}) und kippt die Paarmassen")
say(f"   ({dJ[1]:.1e}); R_int kippt die Paarmassen ({dR[1]:.1e}), hält das Wilson-Glied aber nicht fest ({dR[0]:.2f}).")
befund("R_int = J auch auf dem Wilson-Glied", dR[0] < 1e-9 and dJ[0] < 1e-9)
befund("R_int ist Spin-Z4 (R^2 = (-1)^F) auch auf den Familien", np.allclose(np.array(ch_mu4) ** 2, -1))
say("   Auf den Familien gilt R^2 = {1, -1, -1}. Eine Familie hat gerade Z4-Ladung; die Spin-Z4-Zählung aus Z8")
say("   (alle Fermionen mit Ladung 1 mod 4) gilt für R_int also nicht.")
DATA["seam_half_bundle_char"] = str(np.round(r_s, 9))


# ---------------------------------------------------------------------------------------------------
# 5. E8 oder nicht: GSO-Summe des entstehenden Nahtsektors
# ---------------------------------------------------------------------------------------------------
say("\n5. Twistgewicht und GSO des Nahtsektors (Reihen in t = q^(1/2))\n")
T = 10


def pmul(a, b):
    out = [0] * (T + 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b[:T + 1 - i]):
                out[i + j] += x * y
    return out


def ppow(a, n):
    out = [1] + [0] * T
    for _ in range(n):
        out = pmul(out, a)
    return out


def gso(N):
    nsp, nsm, rr = [1] + [0] * T, [1] + [0] * T, [1] + [0] * T
    for n in range(1, T + 1):
        if 2 * n - 1 <= T:
            f = [0] * (T + 1); f[0] = 1; f[2 * n - 1] = 1; nsp = pmul(nsp, ppow(f, N))
            f[2 * n - 1] = -1; nsm = pmul(nsm, ppow(f, N))
        if 2 * n <= T:
            f = [0] * (T + 1); f[0] = 1; f[2 * n] = 1; rr = pmul(rr, ppow(f, N))
    shift = N // 8                                   # q^(N/16) = t^(N/8)
    R = [0] * (T + 1)
    for i in range(T + 1 - shift):
        R[i + shift] = 2 ** (N // 2) * rr[i]
    tot = [(a + b) // 2 + c // 2 for a, b, c in zip(nsp, nsm, R)]
    return tot, sp.Rational(N, 16)


def e8_char(nmax=4):
    qq = sp.symbols("qq")
    E4 = 1 + 240 * sum(sp.divisor_sigma(n, 3) * qq ** n for n in range(1, nmax + 1))
    eta8 = sp.prod([(1 - qq ** n) ** 8 for n in range(1, nmax + 1)])
    s = sp.series(E4 / eta8, qq, 0, nmax).removeO()
    return [int(s.coeff(qq, n)) for n in range(nmax)]


chi = e8_char()
g16, h16 = gso(16)
g48, h48 = gso(48)
say(f"   16 MW (Grad 1, 16 Paare): Gewicht 1, 2, 3: {g16[2]}, {g16[4]}, {g16[6]}; h_R = {h16}")
say(f"   E4/eta^8:                Gewicht 1, 2, 3: {chi[1]}, {chi[2]}, {chi[3]}")
check("GSO von 16 MW = E8-Charakter bis Gewicht 3", [g16[2], g16[4], g16[6]] == chi[1:4] and h16 == 1)
say(f"   48 MW (Grad 3, 16 Paare, eine gemeinsame Spinstruktur der Naht): Gewicht 1: {g48[2]}, h_R = {h48}")
say(f"   Vergleich: E8^3 hätte 3 x 248 = 744 Ströme. Gemeinsame GSO über 48 Fermionen gibt so(48), kein E8.")
check("48 MW: 1128 Ströme, h_R = 3", g48[2] == 1128 and h48 == 3)
befund("dasselbe Multiplett (Grad 3) gibt E8 auf der Naht", g48[2] == 248)
befund("dasselbe Multiplett (Grad 3) gibt E8^3 auf der Naht", g48[2] == 744)
say("   Die Kopien im Z7-Modell sind entkoppelt; die Zahl der Kernmoden ist daher exakt additiv (gerechnet für")
say("   1 und 2 Paare). 16 Paare mit Grad 1 geben 16 Moden und E8; 16 Paare mit Grad 3 geben 48 Moden.")

# ---------------------------------------------------------------------------------------------------
# 6. Lokale innere Paarung: Grad des Massenbündels auf der Markensphäre
# ---------------------------------------------------------------------------------------------------
say("\n6. Lokale Paarung zweier Nullmoden auf der Markensphäre\n")
mb = {}
for dl in (3, 1):
    e = dl - 1                                       # Nullmoden: Schnitte von E = K^(1/2) (x) L, deg E = deg L - 1
    degPhi = (+2) + (-2) - 2 * e                     # Phi in Kbar (x) K (x) E^-2 (glatt: deg Kbar = +2, deg K = -2)
    mb[dl] = degPhi
    say(f"   deg L = {dl}: psi psi' liegt in E^2 (Grad {2 * e}), Phi muss in Kbar (x) K (x) E^-2 liegen: Grad {degPhi}")
DATA["mass_bundle_degree"] = mb                      # reine Gradarithmetik, daher keine Prüfung
befund("lokale Paarung der Familienmoden ohne Nullstellen auf der Markensphäre", mb[3] == 0)
say("   Für das Familienbündel muss die lokale Paarungsmasse vier Nullstellen haben (Grad -4 = -deg D). Der")
say("   kanonische Schnitt ist conj(s_D): Er verschwindet genau an den vier Marken. Die Marken wären dann")
say("   Wirbel der inneren Masse; sie liegen aber in der inneren Sphäre und füllen die ganze Raumzeit. Eine")
say("   Naht in der Raumzeit (Kodimension 2 dort) folgt daraus nicht. Für Grad 1 gibt es keine Pflichtnullstellen.")

# ---------------------------------------------------------------------------------------------------
# 7. Einwand der Gesamtsynthese: Vierteldrehung, Pin+, 248 Ströme
# ---------------------------------------------------------------------------------------------------
say("\n7. Der Einwand \"Die vier ursprünglichen Marken wählen eine echte lokale Zweiteilung\"\n")
rho_V1 = 1j ** 1
say(f"   Geometrisch auf der Naht: rho = exp(i pi L0/2) hat auf V1 (248 Ströme) den Eigenwert {rho_V1}, rho^2 = {rho_V1 ** 2}.")
say("   Intern auf der Naht: 16 MW = 8 komplexe Fermionen c_j, J c_j = i c_j (Ladung 1). Abzählung der")
say("   Gewicht-1-Zustände im Fockraum (NS-Bilineare und GSO-behaltene Ramond-Grundzustände):")
ph_ns = []
for a, b in itertools.combinations(range(16), 2):    # Moden 0..7 = c_j (Phase i), 8..15 = c_j^dag (Phase -i)
    ph_ns.append((1j if a < 8 else -1j) * (1j if b < 8 else -1j))
spec = {}
for phi0 in (1, -1):
    ph_r = [phi0 * 1j ** sum(occ) for occ in itertools.product((0, 1), repeat=8) if sum(occ) % 2 == 0]
    keys = [(int(round(v.real)), int(round(v.imag))) for v in np.array(ph_ns + ph_r)]
    spec[phi0] = {k: keys.count(k) for k in sorted(set(keys))}
    lab = {(1, 0): "+1", (-1, 0): "-1", (0, 1): "i", (0, -1): "-i"}
    say(f"   Vakuumphase {phi0:+d} der Ramond-Grundzustände: J auf V1 = "
        f"{ {lab[k]: v for k, v in spec[phi0].items()} } (Summe {len(keys)})")
dims_fix = sorted(spec[p].get((1, 0), 0) for p in (1, -1))
say(f"   J ist auf den Strömen eine Involution, Fixalgebra der Dimension {dims_fix}: so(16) (120) oder e7 + su(2)")
say("   (136), die beiden Involutionsklassen von E8. Der Eigenwert i kommt nie vor; rho hat auf allen 248 den Wert i.")
check("J auf den 248 Strömen: nur +-1 (Involution von E8), rho: i",
      all(set(s) <= {(1, 0), (-1, 0)} for s in spec.values()) and dims_fix == [120, 136] and rho_V1 == 1j)
DATA["J_on_V1"] = {str(p): {str(k): v for k, v in s.items()} for p, s in spec.items()}
say("   Innere Sphäre: Im Nahtmultiplett am halben Markenbündel wirkt R_int auf allen 16 Fermionen mit derselben")
say("   Phase +-i (Abschnitt 4b), auf den Strömen also mit +-1, auf V1 ebenfalls ohne den Wert i.")
say(f"   Antiunitäre Konjugation, die mit rho vertauscht, bildet V1 (Eigenwert i) in Gewichte 4m+3 ab;")
say(f"   dim V3 = {chi[3]}. Die 248 Ströme sind unter ihr nicht abgeschlossen: bestätigt.")
befund("eine einzige Vierteldrehung wirkt auf Naht-Strömen, Nahtfermionen und Familien gleich",
       any((0, 1) in s for s in spec.values()) and same_set(np.diag(Jd), ch_mu4))
say("   Der Elternoperator macht den Einwand konkret. Er hat drei verschiedene Z4: die geometrische Drehung der")
say("   Naht (rho), die Monodromie der Nahtfermionen (J) und die Drehung der Markensphäre (R_int). Mit einem")
say("   Multiplett leistet keine Produktsymmetrie {-1, i, -i} und J^2 = -1 zugleich (Abschnitt 3). Mit zwei")
say("   Multipletts tut es R_int (Abschnitt 4b), ist aber auf den Familien kein Spin-Z4 und auf den Strömen")
say("   nicht rho. Eine Pin+-Identifikation ist also nicht automatisch; ein Intertwiner bleibt Zusatzannahme.")

# ---------------------------------------------------------------------------------------------------
# 8. Urteil
# ---------------------------------------------------------------------------------------------------
say("\n8. Was hergeleitet, was eingebaut ist, und das Urteil\n")
say("   Eingebaut: Viertelfluss und windende Masse des Defekts (G, wie Z7); Produktform des Operators;")
say("   gleiche Wurzel K^(1/2) für Spinor und Kopplung; die symmetrische innere Paarungsform C_int.")
say("   Hergeleitet im Lauf: Grad 3 aus vier Marken; ker D+ = 3, ker D- = 0; Charaktere {-1, i, -i} aus dem")
say("   Markenlift; Nahtmoden = Paare x Index, mit Vorzeichen und Lichtgeschwindigkeit; KK-Niveaus tragen nichts;")
say("   R_P auf Nahtmoden = Familiencharaktere; GSO-Zahlen; Grad des Massenbündels.")
say("   Obstruktion in Produktform:")
say("   (a) Trägt das Familienmultiplett die Naht, gibt es 3 x 16 = 48 Moden, so(48) und kein E8.")
say("   (b) Die Paarung, die die Naht macht, macht dieselben Familien schwer (Lücke m).")
say("   (c) Mit einem Multiplett leistet keine Produktsymmetrie {-1, i, -i} auf Familien und J^2 = -1 auf allen")
say("       Nahtmoden.")
say("   Ausweg mit zwei Multipletts (Abschnitt 4b): Naht am halben Markenbündel K^(1/2) (x) M, M^2 = O(D)")
say("   (Index 1, Nullmode omega_T = sqrt(q), Charakter +-i), Materie am Markenbündel (Index 3, ungepaart).")
say("   Dann gelten (i), (ii) mit 16 Paaren und (iii) zugleich, und eine Vierteldrehung wirkt als +-i auf der")
say("   Naht und als {-1, i, -i} auf den Familien. Eingesetzt bleiben: zwei Multipletts als direkte Summe,")
say("   welches an M und welches an O(D) koppelt, und dass nur die Naht gepaart ist. R_int ist auf den")
say("   Familien kein Spin-Z4 und weicht auf dem Wilson-Glied von J ab.")
say("   Urteil: gescheitert für einen Produkt-Elternoperator mit einem Multiplett (Obstruktionen a bis c).")
say("   Lesart für die Variante mit zwei Multipletts am Markenbündel und an seiner Wurzel: alle drei")
say("   Eigenschaften in einem Lauf, aber als direkte Summe mit eingesetzter Zuordnung. Offen: ein Elternoperator")
say("   ohne diese Zuordnung, etwa mit Familien als Zuständen im Spinorteil 128 des Nahtkerns (H6h7).")
DATA["verdict"] = ("gescheitert (Produktform, ein Multiplett); Lesart (zwei Multipletts an O(D) und an "
                   "O(D)^(1/2), direkte Summe); offen (Nicht-Produkt)")

say("\nZusammenfassung der Prüfungen")
for n, c in CHECKS.items():
    say(f"   {'bestanden' if c else 'NICHT bestanden'}: {n}")
say("Befunde (negativ = ehrlich negativ, nicht erzwungen)")
for n, c in BEFUNDE.items():
    say(f"   {'ja' if c else 'nein'}: {n}")
(HERE / "ausgabe_z12.txt").write_text("\n".join(LINES) + "\n")
(HERE / "z12_ergebnis.json").write_text(json.dumps({"checks": CHECKS, "befunde": BEFUNDE, "data": DATA},
                                                   ensure_ascii=False, indent=1, default=str))
assert all(CHECKS.values())
