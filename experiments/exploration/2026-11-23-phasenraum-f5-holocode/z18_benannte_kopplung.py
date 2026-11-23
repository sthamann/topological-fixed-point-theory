"""Z18: Ein Operator mit benannten Kopplungen (offene Verbindung 2 aus H6h10).

Ansatz. Statt zweier Multipletts auf der Markensphäre (Z12) ein einziges Multiplett auf dem Nahttorus
E: w^2 = z^4 - 1 mit geeichtem Blattwechsel sigma: w -> -w. Die Nullmoden des Dirac-Operators auf E,
gekoppelt so, dass K_E^(1/2) (x) L = K_E(mR) (R = Verzweigungsdivisor = die vier Marken auf E), sind
meromorphe 1-Formen mit Polen höchstens der Ordnung m an R. Für m = 1 ist das genau das zurückgezogene
Familienbündel pi^* K(D). Geprüft wird:
 1. Automorphismen von E: sigma, Markendrehung rho (z -> iz), iota (z -> 1/z, w -> i w/z^2) mit iota^2 = sigma.
 2. Nullmodenräume K_E(mR), m = 0..3: Basis, Polordnungen, Riemann-Roch, Deckteilung, Charaktere,
    Lefschetz-Gegenrechnung aus den Fixpunktdaten.
 3. Eine Struktur statt zwei Bündel: pi^* K(D) = K_E (x) O(R), O(R)^2 = pi^* O(D); Träger und T.
 4. J als Spinlift des Blattwechsels; Kopplung von J-Ladung und Decksektor durch die Eichbedingung;
    Vierteldrehungen und Verflechter.
 5. Paarung: welche Massen der Elternoperator erzwingt (exakte Identität q = f2^2 - f0^2).
 6. E8 auf der Naht und Gegenprobe am Z7-Defekt.
 7. Gitter: tiefstes Landauband auf dem quadratischen Torus mit Deckprojektion (Index je Sektor).
"""
import json
from collections import deque
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

HERE = Path(__file__).parent
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


def lab(x, tol=1e-9):
    for v, s in ((1, "1"), (-1, "-1"), (1j, "i"), (-1j, "-i")):
        if abs(x - v) < tol:
            return s
    return f"{complex(np.round(x, 4))}"


def fmt(vals):
    return "{" + ", ".join(lab(v) for v in vals) + "}"


def same_multiset(a, b, tol=1e-8):
    a, b = list(np.array(a, complex)), list(np.array(b, complex))
    if len(a) != len(b):
        return False
    for x in a:
        j = next((k for k, y in enumerate(b) if abs(x - y) < tol), None)
        if j is None:
            return False
        b.pop(j)
    return True


# ---------------------------------------------------------------------------------------------------
# 1. Der Nahttorus und seine Automorphismen
# ---------------------------------------------------------------------------------------------------
say("1. Nahttorus E: w^2 = z^4 - 1 und seine Automorphismen\n")
z, w, t = sp.symbols("z w t")
curve = w ** 2 - (z ** 4 - 1)
AUT = {  # (Z(z, w), W(z, w), dZ/dz)
    "sigma": (z, -w, sp.Integer(1)),
    "rho": (sp.I * z, w, sp.I),
    "iota": (1 / z, sp.I * w / z ** 2, -1 / z ** 2),
}


def on_curve(Z, W):
    return sp.simplify((W ** 2 - (Z ** 4 - 1)).subs(w ** 2, z ** 4 - 1).subs(w, sp.sqrt(z ** 4 - 1))) == 0


def compose(a, b):          # a nach b: (z, w) -> a(b(z, w))
    Za, Wa, _ = AUT[a]
    Zb, Wb, _ = AUT[b]
    return (sp.simplify(Za.subs({z: Zb, w: Wb}, simultaneous=True)),
            sp.simplify(Wa.subs({z: Zb, w: Wb}, simultaneous=True)))


ok_curve = all(on_curve(Z, W) for Z, W, _ in AUT.values())
ii = compose("iota", "iota")
rs, sr = compose("rho", "sigma"), compose("sigma", "rho")
rr = compose("rho", "rho")
say(f"   sigma, rho, iota bilden E auf sich ab: {ok_curve}")
say(f"   iota^2 = {ii}, also iota^2 = sigma; rho sigma = sigma rho: {rs == sr}; rho^2 = {rr}")
check("iota^2 = sigma (Blattwechsel ist Quadrat einer Vierteldrehung)", ok_curve and ii == (z, -w) and rs == sr)
omega = 1 / w                                     # omega = dz / w, Koeffizient von dz
pull = {}
for nm, (Z, W, dZ) in AUT.items():
    pull[nm] = sp.simplify(omega.subs({z: Z, w: W}, simultaneous=True) * dZ / omega)
say(f"   omega = dz/w = sqrt(q), q = dz^2/(z^4 - 1): sigma* omega = {pull['sigma']} omega, rho* omega = "
    f"{pull['rho']} omega, iota* omega = {pull['iota']} omega")
check("omega -> -omega (Blatt), i omega (rho und iota)",
      pull["sigma"] == -1 and pull["rho"] == sp.I and pull["iota"] == sp.I)

# Fixpunkte. Endliche Punkte direkt, die zwei Punkte über z = oo über w/z^2 -> +-1.
marks = [1, sp.I, -1, -sp.I]


def fixed_points(nm):
    Z, W, _ = AUT[nm]
    eqs = [e for e in (sp.numer(sp.together(Z - z)), sp.numer(sp.together(W - w)), curve) if e != 0]
    pts = []
    for sol in sp.solve(eqs, [z, w], dict=True):
        p = (sp.simplify(sol[z]), sp.simplify(sol[w]))
        if p not in pts:
            pts.append(p)
    for sgn in (1, -1):                               # Punkt oo_sgn: w ~ sgn z^2
        lim_z = sp.limit(1 / Z.subs(w, sgn * z ** 2 * sp.sqrt(1 - z ** -4)), z, sp.oo)
        if lim_z == 0:
            ratio = sp.limit((W / Z ** 2).subs(w, sgn * z ** 2 * sp.sqrt(1 - z ** -4)), z, sp.oo)
            if ratio == sgn:
                pts.append(("oo", sgn))
    return pts


FIX = {nm: fixed_points(nm) for nm in AUT}
AUT["tau"] = (-z, -w, sp.Integer(-1))                # tau = sigma rho^2
FIX["tau"] = fixed_points("tau")
for nm in ("sigma", "rho", "iota", "tau"):
    say(f"   Fixpunkte von {nm:5s}: {len(FIX[nm])}  {FIX[nm]}")
sig_fix_marks = sorted(str(p[0]) for p in FIX["sigma"]) == sorted(str(m) for m in marks) and all(p[1] == 0 for p in FIX["sigma"])
check("Fixpunktzahlen sigma 4 (die Marken), rho 2, iota 2, tau = sigma rho^2 keine",
      sig_fix_marks and len(FIX["rho"]) == 2 and len(FIX["iota"]) == 2 and len(FIX["tau"]) == 0)
say("   tau = sigma rho^2 ist fixpunktfrei, also eine Verschiebung um einen Halbperiodenpunkt.")
say("   Auf den Marken: rho dreht zyklisch (1 -> i -> -1 -> -i); iota hält 1 und -1 fest und vertauscht i, -i.")


# Ableitung in den Fixpunkten über lokale Parameter (numerisch, mpmath)
mp.mp.dps = 40


def local_point(p, tt):
    """lokaler Parameter tt an p: Marke (t = w), Punkt über z = 0 (t = z), oo_sgn (t = 1/z)."""
    z0, w0 = p
    if z0 == "oo":
        zz = 1 / tt
        return zz, w0 * zz ** 2 * mp.sqrt(1 - tt ** 4)
    z0c, w0c = complex(z0), complex(w0)
    if abs(w0c) < 1e-30:                              # Marke
        return z0c * (1 + tt ** 2) ** mp.mpf(0.25), tt
    zz = z0c + tt
    return zz, w0c * mp.sqrt((zz ** 4 - 1) / (z0c ** 4 - 1))


def local_param(p, zz, ww):
    z0, w0 = p
    if z0 == "oo":
        return 1 / zz
    return ww if abs(complex(w0)) < 1e-30 else zz - complex(z0)


def derivative_at(nm, p):
    Z, W, _ = AUT[nm]
    fZ, fW = sp.lambdify((z, w), Z, "mpmath"), sp.lambdify((z, w), W, "mpmath")
    eps = mp.mpf("1e-15")
    zz, ww = local_point(p, eps)
    return complex(local_param(p, fZ(zz, ww), fW(zz, ww)) / eps)


DER = {nm: [derivative_at(nm, p) for p in FIX[nm]] for nm in ("sigma", "rho", "iota")}
for nm, d in DER.items():
    say(f"   Drehung in den Fixpunkten von {nm}: {fmt(d)}")
check("sigma dreht in allen vier Marken um pi (d sigma = -1); rho, iota um pi/2 (d = i)",
      np.allclose(DER["sigma"], -1, atol=1e-10) and np.allclose(DER["rho"], 1j, atol=1e-10)
      and np.allclose(DER["iota"], 1j, atol=1e-10))
say("   Folgerung. Ein Spinlift von sigma wirkt in jedem Fixpunkt als Wurzel aus d sigma = -1, also als +-i.")
say("   Damit gilt sigma~^2 = -1 auf jedem Spinor, für jede Spinstruktur: Der Spinlift des Blattwechsels hat")
say("   Ordnung 4 auf Fermionen und Ordnung 2 auf Bosonen. Das ist die Struktur von J (J^2 = (-1)^F), hier")
say("   erzwungen durch die vier Fixpunkte, nicht durch einen eingesetzten Viertelfluss.")
DATA["fixpunkte"] = {k: [str(p) for p in v] for k, v in FIX.items()}


# ---------------------------------------------------------------------------------------------------
# 2. Nullmodenräume K_E(mR): meromorphe 1-Formen g(z, w) dz mit Polen höchstens der Ordnung m an R
# ---------------------------------------------------------------------------------------------------
say("\n2. Nullmoden auf E: H^0(K_E(mR)), Deckteilung und Charaktere\n")
say("   Ein Dirac-Operator auf E mit K_E^(1/2) (x) L = K_E(mR) hat ker D+ = H^0(K_E(mR)) und ker D- = 0")
say("   (Grad 4m > 0, Geschlecht 1). Ansatz für die Basis: gerade Formen P(z) dz / w^(2j), ungerade")
say("   Formen Q(z) dz / w^(2j'+1). Geprüft werden Polordnungen, Unabhängigkeit und Riemann-Roch.")
rng = np.random.default_rng(18)
SAMPLES = []
for _ in range(24):
    z0 = complex(rng.normal() + 1j * rng.normal())
    w0 = np.sqrt(z0 ** 4 - 1 + 0j) * (1 if rng.random() < 0.5 else -1)
    SAMPLES.append((z0, w0))
SPECIAL = [(m_, 0) for m_ in marks] + [("oo", 1), ("oo", -1)]


def basis_K(m):
    j, jp = (m + 1) // 2, m // 2
    even = [z ** k / w ** (2 * j) for k in range(4 * j - 1)] if j > 0 else []
    odd = [z ** k / w ** (2 * jp + 1) for k in range(4 * jp + 1)]
    return even, odd


def order_form(g, p, is_form=True):
    """Ordnung von g dz (bzw. der Funktion g) im lokalen Parameter am Punkt p."""
    f = sp.lambdify((z, w), g, "mpmath")
    vals = []
    for tt in (mp.mpf("1e-9"), mp.mpf("1e-11")):
        zz, ww = local_point(p, tt)
        if p[0] == "oo":
            dz = -1 / tt ** 2
        else:
            z0c = complex(p[0])
            dz = z0c * tt / (2 * (1 + tt ** 2) ** mp.mpf(0.75))
        vals.append(f(zz, ww) * (dz if is_form else 1))
    return int(mp.nint(mp.log(abs(vals[0] / vals[1])) / mp.log(100)))


def evalmat(gs, pts=None):
    pts = SAMPLES if pts is None else pts
    return np.array([[complex(sp.lambdify((z, w), g, "numpy")(*p)) for p in pts] for g in gs])


def pull_matrix(gs, nm, is_form=True, weight=None):
    """Matrix von nm^* auf span(gs): nm^* g_j = sum_i M_ij g_i (Ausgleich über Stichproben).
    weight = Zahl der dz-Faktoren (0 Funktion, 1 Form, 2 quadratisches Differential)."""
    Z, W, dZ = AUT[nm]
    wt = (1 if is_form else 0) if weight is None else weight
    B = evalmat(gs).T                                  # Stichproben x Basis
    img = evalmat([g.subs({z: Z, w: W}, simultaneous=True) * dZ ** wt for g in gs]).T
    M, *_ = np.linalg.lstsq(B, img, rcond=None)
    res = np.linalg.norm(B @ M - img) / max(1.0, np.linalg.norm(img))
    Mr = np.round(M.real, 8) + 1j * np.round(M.imag, 8)
    return Mr, res


TABLE = {}
for m in range(4):
    ev, od = basis_K(m)
    gs = ev + od
    allowed = {**{(mk, 0): -m for mk in marks}, ("oo", 1): 0, ("oo", -1): 0}
    ords = {str(p): min(order_form(g, p) for g in gs) for p in SPECIAL}
    pole_ok = all(min(order_form(g, p) for g in gs) >= allowed[p] for p in SPECIAL)
    rank = np.linalg.matrix_rank(evalmat(gs), tol=1e-8)
    rr = 4 * m if m > 0 else 1                          # Riemann-Roch auf dem Torus
    mats = {nm: pull_matrix(gs, nm) for nm in ("sigma", "rho", "iota")}
    res = max(r for _, r in mats.values())
    S = mats["sigma"][0]
    n_even = int(round(np.real(np.trace(S) + len(gs)) / 2))
    n_odd = len(gs) - n_even
    ne = len(ev)
    R = mats["rho"][0]
    rho_even = np.diag(R)[:ne] if ne else np.array([])
    rho_odd = np.diag(R)[ne:]
    diag_ok = np.allclose(R, np.diag(np.diag(R))) and np.allclose(S, np.diag(np.diag(S)))
    iota_ev = np.linalg.eigvals(mats["iota"][0])
    TABLE[m] = dict(dim=len(gs), rank=rank, rr=rr, even=n_even, odd=n_odd, pole_ok=pole_ok, res=res,
                    rho_even=rho_even, rho_odd=rho_odd, iota=iota_ev, diag=diag_ok, gs=gs, ne=ne, mats=mats)
    say(f"   m = {m}: Basis {len(gs)}, Rang {rank}, Riemann-Roch {rr}, Polordnungen zulässig {pole_ok}; "
        f"gerade {n_even}, ungerade {n_odd}")
    say(f"          rho auf gerade {fmt(rho_even)}, auf ungerade {fmt(rho_odd)}; Fehler der Darstellung {res:.1e}")
ok_tab = all(T_["rank"] == T_["dim"] == T_["rr"] and T_["pole_ok"] and T_["res"] < 1e-8 and T_["diag"]
             for T_ in TABLE.values())
check("Basen von H^0(K_E(mR)), m = 0..3: Pole zulässig, unabhängig, vollständig (Riemann-Roch)", ok_tab)
T1 = TABLE[1]
fam_ok = (T1["even"], T1["odd"]) == (3, 1) and same_multiset(T1["rho_even"], [1j, -1, -1j]) \
    and same_multiset(T1["rho_odd"], [1j]) and T1["gs"][3] == 1 / w
say(f"\n   m = 1 (pi^* K(D)): gerade {T1['even']} mit rho-Charakteren {fmt(T1['rho_even'])}, ungerade {T1['odd']}: "
    f"die Form dz/w = sqrt(q) mit {fmt(T1['rho_odd'])}.")
check("m = 1: drei gerade Nullmoden {-1, i, -i} (Familien) und genau eine ungerade, sqrt(q) mit Charakter i", fam_ok)
say("   Die geraden Formen P(z) dz/(z^4 - 1) sind genau H^0(CP^1, K(D)); die ungerade ist die Nahtnullmode aus Z12.")

# Lefschetz aus Fixpunktdaten (unabhängig von den Basen): L(g) = sum_p a_p^(1 - m_p) / (1 - a_p)
say("\n   Gegenrechnung mit der holomorphen Lefschetzformel aus den Fixpunktdaten allein:")
lef_ok = True
for m in (1, 2, 3):
    row = []
    for nm in ("sigma", "rho", "iota"):
        tot = 0
        for p, a in zip(FIX[nm], DER[nm]):
            mp_ = m if (p[0] != "oo" and abs(complex(p[1])) < 1e-12) else 0
            tot += a ** (1 - mp_) / (1 - a)
        tr = np.trace(TABLE[m]["mats"][nm][0])
        lef_ok &= abs(tot - tr) < 1e-8
        row.append(f"{nm}: {complex(np.round(tot, 6))} gegen Spur {complex(np.round(tr, 6))}")
    say(f"   m = {m}: " + "; ".join(row))
check("Lefschetz aus Fixpunkten = Spur auf den Nullmoden (sigma, rho, iota; m = 1, 2, 3)", lef_ok)
say("   Mit gerade + ungerade = 4m und gerade - ungerade = 2(-1)^(m+1) (vier Marken, je +-1/2) folgt")
say("   ungerade = 2m + (-1)^m und gerade = 2m - (-1)^m für alle m >= 1.")
counts = {m: (2 * m - (-1) ** m, 2 * m + (-1) ** m) for m in range(1, 40)}
uniq = [m for m, (e_, o_) in counts.items() if (e_, o_) == (3, 1)]
check("Formel stimmt mit m = 1..3 überein, und nur m = 1 gibt 3 gerade + 1 ungerade",
      all(counts[m] == (TABLE[m]["even"], TABLE[m]["odd"]) for m in (1, 2, 3)) and uniq == [1]
      and (TABLE[0]["even"], TABLE[0]["odd"]) == (0, 1))
say(f"   m = 0 gibt (0, 1): nur die Nahtmode, keine Familien. Drei Familien und genau eine Nahtmode verlangen m = 1.")

# Vergleich: Funktionen statt Formen (gleiches Bündel O(R), andere Wirkung)
LR = [sp.Integer(1)] + [z ** k / w for k in range(3)]
SR = pull_matrix(LR, "sigma", is_form=False)[0]
RR_ = pull_matrix(LR, "rho", is_form=False)[0]
say(f"\n   Gegenprobe Funktionen in O(R) = M_E: sigma {fmt(np.diag(SR))}, rho {fmt(np.diag(RR_))}.")
say("   Als Funktionen sind die Familien die ungeraden Schnitte z^k/w (Charaktere {1, i, -1}) und die Naht die")
say("   gerade Konstante. Gerade/ungerade hängt an der Linearisierung; invariant ist nur: Familien und Naht")
say("   liegen in entgegengesetzten Sektoren, 3 gegen 1. Die Familiencharaktere {-1, i, -i} verlangen Formen.")
befund("Familiencharaktere {-1,i,-i} schon auf Funktionen (ohne K_E)", same_multiset(np.diag(RR_)[1:], [1j, -1, -1j]))
DATA["tabelle_KmR"] = {m: {k: (fmt(v) if isinstance(v, np.ndarray) else v) for k, v in T_.items()
                           if k in ("dim", "even", "odd", "rho_even", "rho_odd", "iota")} for m, T_ in TABLE.items()}


# ---------------------------------------------------------------------------------------------------
# 3. Eine Struktur statt zwei Bündel
# ---------------------------------------------------------------------------------------------------
say("\n3. Eine Struktur: pi^* K(D) = K_E (x) O(R) mit O(R)^2 = pi^* O(D)\n")
gfun = [sp.Integer(1)] + [z ** k / w for k in range(3)]          # H^0(O(R)) = L(R)
prod_om = [sp.simplify(g * omega) for g in gfun]
same = set(map(sp.simplify, prod_om)) == set(map(sp.simplify, T1["gs"]))
say(f"   H^0(O(R)) = span{{1, z^k/w}}; mal omega ergibt {prod_om}: genau die Basis von H^0(K_E(R)): {same}")
L2R_even = [z ** k / (z ** 4 - 1) for k in range(5)]            # Träger V = H^0(CP^1, O(D))
L2R_odd = [z ** k / w for k in range(3)]
ords2 = [min(order_form(g, p, is_form=False) for g in L2R_even + L2R_odd) for p in SPECIAL]
rank2 = np.linalg.matrix_rank(evalmat(L2R_even + L2R_odd), tol=1e-8)
say(f"   L(2R) = H^0(O(2R)) = H^0(pi^* O(D)): Rang {rank2} = Grad 8; gerade 5 = Träger V, ungerade 3; "
    f"Polordnungen an Marken/oo {ords2}")
prods_odd = [sp.simplify(gfun[1 + k] * gfun[1 + l]) for k in range(3) for l in range(k, 3)]
rk_sym = np.linalg.matrix_rank(evalmat(prods_odd), tol=1e-8)
rk_with = np.linalg.matrix_rank(evalmat(prods_odd + L2R_even), tol=1e-8)
conic = sp.simplify(gfun[1] * gfun[3] - gfun[2] ** 2)
say(f"   Produkte der drei ungeraden Schnitte z^k/w: 6 Stück, Rang {rk_sym}, spannen den Träger auf (Rang mit V: "
    f"{rk_with}); Kern = Kegelschnitt g0 g2 - g1^2 = {conic}.")
check("Ein Bündel: K_E(R) = O(R) (x) omega; Träger = Sym^2(Familienschnitte von O(R)) modulo Kegelschnitt",
      same and rank2 == 8 and min(ords2) >= -2 and rk_sym == 5 and rk_with == 5 and conic == 0)
say("   Lesart: Das Nahtbündel M aus Z12 ist auf E das echte Bündel O(R) (O(R)^2 = pi^* O(D)), das")
say("   Materiebündel ist O(R) (x) omega. Die Nahtmode ist 1 (x) omega = sqrt(q), die Familien sind")
say("   (z^k/w) (x) omega. Der Träger V ist das Quadrat der Familienschnitte geteilt durch q: V q = Sym^2 F / Kegel.")

# T als Deckblock einer Abbildung auf E
fam_basis = [z ** k / (z ** 4 - 1) for k in range(3)]


def P0(f):
    return sp.simplify(sum(f.subs(z, sp.I ** j * z) for j in range(4)) / 4)


def T_map(f):
    return sp.simplify((f - P0(f)) / z)                            # Koeffizient von dz


Tmat = np.zeros((3, 5), complex)
for j, f in enumerate(L2R_even):
    img = sp.expand(sp.cancel(T_map(f) * (z ** 4 - 1)))
    for k in range(3):
        Tmat[k, j] = complex(sp.Poly(img, z).coeff_monomial(z ** k)) if img != 0 else 0
theta_comm = all(sp.simplify((f / z).subs(w, -w) - f / z) == 0 for f in L2R_even) and \
    all(sp.simplify((g / z).subs(w, -w) + g / z) == 0 for g in L2R_odd)
Rc = np.diag([1j ** k for k in range(5)])
Rf = np.diag([1j * 1j ** k for k in range(3)])
equiv = np.allclose(Tmat @ Rc, Rf @ Tmat)
chars_V, chars_F = [1j ** k for k in range(5)], [1j ** (k + 1) for k in range(3)]
dim_equiv = sum(1 for a in chars_V for b in chars_F if abs(a - b) < 1e-12)
say(f"   T(f) = (1 - P0)(f) dz/z: Matrix {np.round(Tmat.real, 3).tolist()}, Rang {np.linalg.matrix_rank(Tmat)}, "
    f"äquivariant {equiv}.")
say(f"   Auf E ist Theta: f -> f dz/z mit sigma vertauschbar ({theta_comm}); T ist sein gerader Block, gefolgt")
say(f"   von der Projektion weg von den rho-invarianten Moden. Äquivariante Abbildungen V -> F: {dim_equiv} Parameter.")
check("T = gerader Block von Theta = (.) dz/z, surjektiv, äquivariant",
      np.linalg.matrix_rank(Tmat) == 3 and equiv and theta_comm)
befund("T ist durch Deckstruktur und Vierteldrehung eindeutig festgelegt", dim_equiv == 1)
say("   Die Deckstruktur legt T nicht fest: Theta ist ein Rückzug von CP^1 (dz/z ist deckgerade), und P0 kommt")
say("   von der Markendrehung. T bleibt die äquivariante Wahl mit drei gleichen Koeffizienten.")

# ---------------------------------------------------------------------------------------------------
# 4. J als Spinlift des Blattwechsels; Eichbedingung koppelt J-Ladung an den Decksektor
# ---------------------------------------------------------------------------------------------------
say("\n4. J als Spinlift des Blattwechsels und die Eichbedingung\n")
say("   Fermion = Spinor (x) Eichbündel. Spinstruktur (P,P): K_E^(1/2) trivial mit sqrt(omega), sigma~ sqrt(omega)")
say("   = eps i sqrt(omega). Eichbündel O(R) mit der ehrlichen Wirkung sigma^* (Ordnung 2, bosonisch).")
S_fun = pull_matrix(gfun, "sigma", is_form=False)[0]
J_geo = {eps: eps * 1j * S_fun for eps in (1, -1)}
for eps, Jg in J_geo.items():
    say(f"   eps = {eps:+d}: sigma^ auf (Naht, f0, f1, f2) = {fmt(np.diag(Jg))}, sigma^2 = {fmt(np.diag(Jg @ Jg))}")
check("sigma^ = sigma~ (x) sigma^*: Ordnung 4, Naht und Familien mit entgegengesetzter Ladung (+-i, 1 gegen 3)",
      all(np.allclose(Jg @ Jg, -np.eye(4)) and same_multiset(np.diag(Jg), [eps * 1j, -eps * 1j, -eps * 1j, -eps * 1j])
          for eps, Jg in J_geo.items()))

say("\n   Eichung des Blattwechsels auf Fermionen braucht eine innere Phase J_int mit J_int^2 = -1, damit")
say("   sigma^ (x) J_int^(-1) die Ordnung 2 hat. Für die Paare (psi, chi) aus Z7 ist das J mit xi+- = +-i.")
say("   Physikalische Zustände: Eigenwert 1 von sigma^ (x) J_int^(-1).")
Jint = np.diag([1j, -1j])                                         # (xi+, xi-)
kept = {}
for eps, Jg in J_geo.items():
    G = np.kron(np.linalg.inv(Jint), Jg)
    vals, vecs = np.linalg.eig(G)
    sel = np.abs(vals - 1) < 1e-9
    K = vecs[:, sel]
    w_xi = np.sum(np.abs(K.reshape(2, 4, -1)) ** 2, axis=1)       # Gewicht in xi+ / xi-
    seam_xi = [int(np.argmax(w_xi[:, c])) for c in range(K.shape[1])]
    mode_w = np.sum(np.abs(K.reshape(2, 4, -1)) ** 2, axis=0)
    seam_mode = [int(np.argmax(mode_w[:, c])) for c in range(K.shape[1])]
    kept[eps] = dict(n=K.shape[1], xi=seam_xi, mode=seam_mode)
    nseam = sum(1 for a, b in zip(seam_xi, seam_mode) if b == 0)
    say(f"   eps = {eps:+d}: je Paar {K.shape[1]} Moden; Naht sqrt(q) in xi{'+' if seam_xi[seam_mode.index(0)] == 0 else '-'}, "
        f"Familien in xi{'+' if seam_xi[seam_mode.index(1)] == 0 else '-'}; 16 Paare: {16 * nseam} Naht + {16 * (K.shape[1] - nseam)} Familien")
lock_ok = all(v["n"] == 4 for v in kept.values()) and \
    kept[1]["xi"][kept[1]["mode"].index(0)] == 0 and kept[-1]["xi"][kept[-1]["mode"].index(0)] == 1
check("Eichbedingung: je Paar 1 Nahtmode + 3 Familien, Naht und Familien in entgegengesetzten J-Sektoren", lock_ok)
say("   Mit 16 Paaren: 16 Nahtmoden und 48 Familienkomponenten, aus einem Multiplett und einem Bündel.")
say("   Welcher J-Sektor die Naht trägt, legt die Orientierung eps fest. Z7 lässt xi+ winden; E8 verlangt dann,")
say("   dass xi+ die Naht trägt (sonst winden 48 Moden, Abschnitt 6): eps = +1 ist durch E8 festgelegt.")

say("\n   Vierteldrehungen und J:")
Z1 = T1["gs"]
Ms = {nm: pull_matrix(Z1, nm)[0] for nm in ("sigma", "rho", "iota", "tau")}
r2 = Ms["rho"] @ Ms["rho"]
i2 = Ms["iota"] @ Ms["iota"]
say(f"   Basis (f0, f1, f2, Naht). rho^*^2 = {fmt(np.diag(r2))}, sigma^* = {fmt(np.diag(Ms['sigma']))}, "
    f"tau^* = {fmt(np.diag(Ms['tau']))}")
say(f"   iota^*^2 = sigma^*: {np.allclose(i2, Ms['sigma'])}; rho^*^2 = sigma^* tau^*: {np.allclose(r2, Ms['sigma'] @ Ms['tau'])}")
check("iota^2 = sigma und rho^2 = sigma tau auf den Nullmoden", np.allclose(i2, Ms["sigma"]) and np.allclose(r2, Ms["sigma"] @ Ms["tau"]))
befund("Markendrehung rho^2 = Blattwechsel auf allen Nullmoden", np.allclose(r2, Ms["sigma"]))
say(f"   Auf der Nahtmode gilt rho^2 = sigma (tau = +1). Auf den Familien f0, f2 (Charaktere i, -i) ist tau = -1.")
say("   Verflechter: J = sigma~ = iota~^2 exakt (iota: z -> 1/z, w -> i w/z^2, Drehung um die Marken +-1).")
say("   Die Markendrehung rho erfüllt rho~^2 = J tau~; tau ist die Verschiebung um eine Halbperiode und wirkt")
say("   auf zwei Familien mit -1. iota ist keine Markendrehung (auf den Marken eine Spiegelung) und gibt auf")
say(f"   den Familien {fmt(np.linalg.eigvals(Ms['iota'][:3, :3]))} statt {{-1, i, -i}}.")

# Gruppe aller hier betrachteten Automorphismen (erzeugt von sigma, rho, iota) und ihre Wurzeln aus sigma
def comp_maps(A, B):                                   # A nach B
    return (sp.simplify(A[0].subs({z: B[0], w: B[1]}, simultaneous=True)),
            sp.simplify(A[1].subs({z: B[0], w: B[1]}, simultaneous=True)))


gens = [(AUT[n][0], AUT[n][1]) for n in ("sigma", "rho", "iota")]
group, frontier = [(z, w)], [(z, w)]
while frontier:
    new_ = []
    for g in frontier:
        for h in gens:
            gh = comp_maps(h, g)
            if not any(sp.simplify(gh[0] - e[0]) == 0 and sp.simplify(gh[1] - e[1]) == 0 for e in group):
                group.append(gh)
                new_.append(gh)
    frontier = new_
roots = [g for g in group if comp_maps(g, g) == (z, -w)]


def mark_perm(g):
    return tuple(marks.index(sp.nsimplify(sp.simplify(g[0].subs(z, mk)))) for mk in marks)


cyc = [g for g in roots if len(set(mark_perm(g))) == 4 and mark_perm(comp_maps(g, g)) == (0, 1, 2, 3)
       and all(mark_perm(g)[k] != k for k in range(4))]
perms = sorted(set(mark_perm(g) for g in roots))
say(f"   Gruppe aus sigma, rho, iota: Ordnung {len(group)}. Wurzeln aus sigma (g^2 = sigma): {len(roots)}; ihre")
say(f"   Wirkung auf die Marken {perms}. Jede solche Wurzel wird auf CP^1 eine Involution und kann die vier")
say("   Marken nicht zyklisch drehen.")
check("Gruppe der Ordnung 16 mit Wurzeln aus sigma", len(group) == 16 and len(roots) > 0)
befund("eine Vierteldrehung mit Quadrat J, die zugleich die Marken zyklisch dreht", len(cyc) > 0)
DATA["gruppe"] = dict(ordnung=len(group), wurzeln_aus_sigma=len(roots), markenpermutationen=[str(p_) for p_ in perms])
DATA["J_sektoren"] = {str(k): v for k, v in kept.items()}


# ---------------------------------------------------------------------------------------------------
# 5. Paarung: was ein Multiplett mit einem Paarungsprofil erzwingt
# ---------------------------------------------------------------------------------------------------
say("\n5. Paarung der Nullmoden: Nahtkanal q und Familienkanäle\n")
say("   Eine 6D-Majoranamasse koppelt zwei Nullmoden über ihr Produkt (ein quadratisches Differential) an ein")
say("   inneres Profil s. 4D-Masse M_ab = <s, psi_a psi_b> im L^2-Produkt auf E. Produkte der Familien:")
say("   f_k f_l = z^(k+l) dz^2/(z^4 - 1)^2 = p_(k+l); Nahtkanal sqrt(q)^2 = q = dz^2/(z^4 - 1).")
ident = sp.simplify(z ** 4 / (z ** 4 - 1) ** 2 - 1 / (z ** 4 - 1) ** 2 - 1 / (z ** 4 - 1))
say(f"   Exakt: f2^2 - f0^2 - q = {ident}. Der Nahtkanal liegt im Aufspann der Familienkanäle.")
Pbasis = [z ** k_ / w ** 4 for k_ in range(5)]
PM = {nm: pull_matrix(Pbasis, nm, weight=2)[0] for nm in ("sigma", "rho", "iota")}
say(f"   Auf P = span(p_0..p_4): rho^* = {fmt(np.diag(PM['rho']))}, sigma^* = {fmt(np.diag(PM['sigma']))}, "
    f"iota^*: p_m -> p_(4-m) {np.allclose(PM['iota'], np.eye(5)[::-1])}")
qv = np.array([-1, 0, 0, 0, 1], complex)
cov = np.vstack([PM["rho"] + np.eye(5), PM["iota"] + np.eye(5), PM["sigma"] - np.eye(5)])
_, sv, Vh = np.linalg.svd(cov)
null = Vh[np.sum(sv > 1e-9):].conj().T
say(f"   Profile mit den Charakteren von q (rho -1, iota -1, sigma +1): Dimension {null.shape[1]}, "
    f"parallel zu q: {np.linalg.matrix_rank(np.column_stack([null, qv]), tol=1e-9) == 1}")
check("q = f2^2 - f0^2 exakt; das kovariante Paarungsprofil ist eindeutig s ~ q",
      ident == 0 and null.shape[1] == 1 and np.linalg.matrix_rank(np.column_stack([null, qv]), tol=1e-9) == 1)

# alle G-invarianten hermiteschen Produkte auf P
nP = 5
av = sp.symbols(f"a0:{nP * nP}", real=True)
bv = sp.symbols(f"b0:{nP * nP}", real=True)
Hs = sp.Matrix(nP, nP, lambda i, j: av[i * nP + j] + sp.I * bv[i * nP + j])
eqs = list(Hs - Hs.H)
for nm in ("rho", "iota"):
    Rm = sp.Matrix(nP, nP, lambda i, j: sp.nsimplify(complex(PM[nm][i, j]).real) + sp.I * sp.nsimplify(complex(PM[nm][i, j]).imag))
    eqs += list(Rm.H * Hs * Rm - Hs)
real_eqs = []
for e in eqs:
    e = sp.expand(e)
    real_eqs += [sp.re(e), sp.im(e)]
sol = sp.solve([e for e in real_eqs if e != 0], list(av) + list(bv), dict=True)[0]
Hinv = sp.simplify(Hs.subs(sol))
free = sorted(Hinv.free_symbols, key=str)
qS = sp.Matrix([-1, 0, 0, 0, 1])
ip = lambda x, y: sp.simplify((x.H * Hinv * y)[0])
e_ = lambda k_: sp.Matrix([1 if i == k_ else 0 for i in range(nP)])
qq = ip(qS, qS)
fam_ch = {(k_, l_): ip(qS, e_(k_ + l_)) for k_ in range(3) for l_ in range(k_, 3)}
say(f"   G-invariante hermitesche Produkte auf P: {len(free)} freie Parameter {free}.")
say(f"   <q, q> = {qq}; <q, f_k f_l>: " + ", ".join(f"f{k_}f{l_}: {v}" for (k_, l_), v in fam_ch.items()))
exact_ok = (sp.simplify(fam_ch[(2, 2)] - qq / 2) == 0 and sp.simplify(fam_ch[(0, 0)] + qq / 2) == 0
            and all(fam_ch[kl] == 0 for kl in fam_ch if kl not in ((0, 0), (2, 2))))
check("Für jedes G-invariante Produkt: <q,f2^2> = -<q,f0^2> = <q,q>/2, alle anderen Familienkanäle 0", exact_ok)

# 4D-Massen mit der Paarung psi chi aus G = (i/2)(xi+ xi+ - xi- xi-): Naht in xi+, Familien in xi-
rngH = np.random.default_rng(5)
while True:
    vals = {s_: rngH.normal() for s_ in free}
    Hn = np.array(Hinv.subs(vals).evalf(), dtype=complex)
    if np.all(np.linalg.eigvalsh(Hn) > 0):
        break
qn = float(sp.re(qq.subs(vals)))
Mfam = np.array([[-complex(sp.N(ip(qS, e_(k_ + l_)).subs(vals))) for l_ in range(3)] for k_ in range(3)])
svf = np.sort(np.linalg.svd(Mfam, compute_uv=False) / qn)
say(f"   Mit G (psi chi, beide xi gleich stark) und s = q: Nahtmasse 1 (in Einheiten m <q,q>), Familienmassen "
    f"{np.round(svf, 6)}.")
say("   Leicht bleibt nur die Familie f1 (Charakter -1); f0 und f2 (Charaktere i, -i) erhalten genau m/2.")
befund("drei leichte Familien mit der Paarung psi chi aus G (Z7) im einen Multiplett", np.sum(svf < 1e-9) == 3)
check("Paarung psi chi: Familienmassen exakt {0, 1/2, 1/2} der Nahtmasse", np.allclose(svf, [0, 0.5, 0.5], atol=1e-9))
say("   Mit der J-projizierten Paarung Phi xi+ xi+ (nur der Nahtsektor) bleiben alle drei Familien leicht.")

# Auswahlregeln: kann irgendeine Ladung von Phi die Naht paaren und alle Familien verbieten?
say("\n   Auswahlregeln für alle Ladungen von Phi unter (sigma, rho, iota) auf den vier Nullmoden:")
Rz = {nm: Ms[nm] for nm in ("sigma", "rho", "iota")}
sym_basis = []
for i_ in range(4):
    for j_ in range(i_, 4):
        E = np.zeros((4, 4))
        E[i_, j_] = E[j_, i_] = 1
        sym_basis.append(E)
roots4 = [1, 1j, -1, -1j]
sep, seam_ok = [], []
for cs in roots4[:2] + roots4[2:]:
    for cr in roots4:
        for ci in roots4:
            rows = []
            for nm, chi in (("sigma", cs), ("rho", cr), ("iota", ci)):
                Rg = Rz[nm]
                rows.append(np.column_stack([(chi * Rg.T @ Eb @ Rg - Eb).ravel() for Eb in sym_basis]))
            A_ = np.vstack(rows)
            _, sv_, Vh_ = np.linalg.svd(A_)
            ns = Vh_[np.sum(sv_ > 1e-9):].conj().T
            if ns.shape[1] == 0:
                continue
            Ms_space = [sum(c * Eb for c, Eb in zip(ns[:, c_], sym_basis)) for c_ in range(ns.shape[1])]
            seam_allowed = any(abs(M_[3, 3]) > 1e-9 for M_ in Ms_space)
            fam_allowed = any(np.abs(M_[:3, :3]).max() > 1e-9 for M_ in Ms_space)
            if seam_allowed:
                seam_ok.append((lab(cs), lab(cr), lab(ci)))
                if not fam_allowed:
                    sep.append((lab(cs), lab(cr), lab(ci)))
say(f"   Ladungen, die die Nahtmode paaren dürfen: {seam_ok}; davon ohne jede Familienmasse: {sep}")
befund("eine Symmetrieladung von Phi paart die Naht und verbietet alle Familienmassen", len(sep) > 0)
check("Nahtpaarung ist für genau eine Ladung erlaubt, und dort sind Familienmassen erlaubt",
      len(seam_ok) == 1 and len(sep) == 0)
say("   Die Geometrie von (E, Marken) trennt Naht- und Familienpaarung also nicht. Getrennt wird erst durch")
say("   die Kopplung an den J-Projektor P+ = (1 - iJ)/2: Phi koppelt nur an xi+, das ist nach Abschnitt 4 der")
say("   Nahtsektor. Das ist eine benannte Kopplung, aber eine Eingabe; G (Z7) setzt psi chi, nicht P+.")
DATA["paarung"] = dict(familienmassen_psichi=[float(x) for x in svf], erlaubte_ladungen=seam_ok)


# ---------------------------------------------------------------------------------------------------
# 6. E8 auf der Naht und Gegenprobe am Gitterdefekt
# ---------------------------------------------------------------------------------------------------
say("\n6. E8-Bedingung für den Nahtsektor\n")
TT = 8


def pmul(a, b):
    out = [0] * (TT + 1)
    for i_, x in enumerate(a):
        if x:
            for j_, y in enumerate(b[:TT + 1 - i_]):
                out[i_ + j_] += x * y
    return out


def gso(N):
    nsp, nsm, rr_ = [1] + [0] * TT, [1] + [0] * TT, [1] + [0] * TT
    for n in range(1, TT + 1):
        if 2 * n - 1 <= TT:
            for _ in range(N):
                f = [0] * (TT + 1); f[0] = 1; f[2 * n - 1] = 1; nsp = pmul(nsp, f)
                f = [0] * (TT + 1); f[0] = 1; f[2 * n - 1] = -1; nsm = pmul(nsm, f)
        if 2 * n <= TT:
            for _ in range(N):
                f = [0] * (TT + 1); f[0] = 1; f[2 * n] = 1; rr_ = pmul(rr_, f)
    shift = N // 8 if N % 8 == 0 else None
    R_ = [0] * (TT + 1)
    if shift is not None:
        for i_ in range(TT + 1 - shift):
            R_[i_ + shift] = 2 ** (N // 2) * rr_[i_]
    return [(a + b) // 2 + c // 2 for a, b, c in zip(nsp, nsm, R_)], sp.Rational(N, 16)


qq_ = sp.symbols("qq")
E4 = 1 + 240 * sum(sp.divisor_sigma(n, 3) * qq_ ** n for n in range(1, 5))
eta8 = sp.prod([(1 - qq_ ** n) ** 8 for n in range(1, 5)])
chi_e8 = sp.series(E4 / eta8, qq_, 0, 4).removeO()
chi_e8 = [int(chi_e8.coeff(qq_, n)) for n in range(4)]
seam_counts = {"m = 1, eps = +1 (xi+ trägt sqrt(q))": 16 * TABLE[1]["odd"],
               "m = 2 (Träger x sqrt(q))": 16 * TABLE[2]["odd"],
               "m = 1, eps = -1 (xi+ trägt die Familien)": 16 * TABLE[1]["even"]}
res6 = {}
for lbl, N in seam_counts.items():
    g_, h_ = gso(N)
    res6[lbl] = (N, g_[2], h_)
    say(f"   {lbl:42s}: {N} Majorana-Weyl, Gewicht-1-Zustände {g_[2]}, Twistgewicht {h_}")
g16, _ = gso(16)
say(f"   E4/eta^8: {chi_e8[1]}, {chi_e8[2]}, {chi_e8[3]}; 16 MW: {g16[2]}, {g16[4]}, {g16[6]}")
check("Nahtsektor aus einem Multiplett (m = 1, eps = +1): 16 MW, Charakter = E8 bis Gewicht 3, Twistgewicht 1",
      seam_counts["m = 1, eps = +1 (xi+ trägt sqrt(q))"] == 16 and [g16[2], g16[4], g16[6]] == chi_e8[1:4]
      and gso(16)[1] == 1)
check("Gegenproben: m = 2 gibt 80 MW (3160 Ströme), falsche Orientierung 48 MW (1128)",
      res6["m = 2 (Träger x sqrt(q))"][:2] == (80, 3160) and res6["m = 1, eps = -1 (xi+ trägt die Familien)"][:2] == (48, 1128))

say("\n   Gegenprobe am Gitterdefekt aus Z7 (verdrillter Defekt, Viertelfluss): Kernmoden je Paar.")
_src = (HERE / "z7_ursprung_naht.py").read_text()
_ns = {"__file__": str(HERE / "z7_ursprung_naht.py")}
exec(_src[:_src.index("L = 14")], _ns)
cb = {}
for cp in (1, 2):
    nlow, cval, vc, ve = _ns["core_branches"](12, 1, 0, 0.25, copies=cp)
    cb[cp] = int(np.sum(vc > 0.2) - np.sum(vc < -0.2))
    say(f"   {cp} Paar(e): netto {cb[cp]} chirale Kernmoden, Geschwindigkeiten {np.round(vc, 3)}")
check("Z7-Defekt: eine chirale Mode je windendem Paar (additiv)", cb == {1: 1, 2: 2})
say("   Mit einem inneren Index 1 im Nahtsektor (sqrt(q)) gibt jedes der 16 Paare genau eine Mode (vgl. Z12,")
say("   Nahtmoden = Paare x Index). Die 48 Familienkomponenten in xi- winden nicht und tragen keine Nahtmode.")
DATA["e8"] = {k: [str(x) for x in v] for k, v in res6.items()}


# ---------------------------------------------------------------------------------------------------
# 7. Gitter: tiefstes Landauband auf dem quadratischen Torus mit Deckprojektion
# ---------------------------------------------------------------------------------------------------
say("\n7. Gitterprobe: quadratischer Torus (tau = i), N Flussquanten, Blattwechsel u -> -u\n")
say("   Hofstadter-Gitter L x L mit Landau-Eichung; das tiefste Band ist exakt N-fach entartet (Bild der")
say("   holomorphen Schnitte, Index N). Die Randphasen (0 oder pi) wählen eine der vier Wurzeln bzw.")
say("   Spinstrukturen. Die Eichtransformation zu u -> -u und u -> iu wird über einen Spannbaum bestimmt.")


def hofstadter(L, N, tx=0.0, ty=0.0):
    phi = N / L ** 2
    n = L * L
    idx = lambda x, y: (x % L) * L + (y % L)
    H = np.zeros((n, n), complex)
    for x in range(L):
        for y in range(L):
            H[idx(x, y + 1), idx(x, y)] += -np.exp(2j * np.pi * phi * x) * (np.exp(1j * ty) if y == L - 1 else 1)
            H[idx(x + 1, y), idx(x, y)] += -((np.exp(-2j * np.pi * N * y / L) * np.exp(1j * tx)) if x == L - 1 else 1)
    return H + H.conj().T


def lattice_symmetry(H, L, perm):
    n = L * L
    Pm = np.zeros((n, n))
    for x in range(L):
        for y in range(L):
            x2, y2 = perm(x, y)
            Pm[(x2 % L) * L + (y2 % L), x * L + y] = 1
    Hp = Pm @ H @ Pm.T
    G = np.zeros(n, complex)
    G[0] = 1
    seen = np.zeros(n, bool)
    seen[0] = True
    dq = deque([0])
    while dq:
        i_ = dq.popleft()
        for j_ in np.nonzero(np.abs(H[:, i_]) > 1e-12)[0]:
            if not seen[j_]:
                G[j_] = H[j_, i_] / (Hp[j_, i_] * np.conj(G[i_]))
                seen[j_] = True
                dq.append(j_)
    U = np.diag(G) @ Pm
    return U, np.linalg.norm(U @ H @ U.conj().T - H)


Lg = 24
lat = {}
for N in (4, 8):
    for tw in ((0, 0), (np.pi, 0), (0, np.pi), (np.pi, np.pi)):
        H = hofstadter(Lg, N, *tw)
        E, V = np.linalg.eigh(H)
        low = V[:, :N]
        U, res = lattice_symmetry(H, Lg, lambda x, y: (-x, -y))
        ev = np.linalg.eigvals(low.conj().T @ U.conj().T @ low)
        npl, nmi = int(np.sum(np.abs(ev - 1) < 1e-8)), int(np.sum(np.abs(ev + 1) < 1e-8))
        lat[(N, tw)] = (npl, nmi, E[N] - E[N - 1], E[N - 1] - E[0], res)
        say(f"   N = {N}, Randphasen ({tw[0] / np.pi:.0f}, {tw[1] / np.pi:.0f}) pi: Entartung {E[N - 1] - E[0]:.1e}, Lücke "
            f"{E[N] - E[N - 1]:.3f}; Blattwechsel: {npl} gerade, {nmi} ungerade (Spur {npl - nmi:+d}); Symmetriefehler {res:.1e}")
ok_lat = all(v[3] < 1e-10 and v[2] > 0.05 and v[4] < 1e-10 for v in lat.values())
p4 = sorted(lat[(4, (0, 0))][:2])
p8 = sorted(lat[(8, (0, 0))][:2])
mixed = all(lat[(N, tw)][:2] == (N // 2, N // 2) for N in (4, 8) for tw in ((np.pi, 0), (0, np.pi), (np.pi, np.pi)))
check("Gitter: Teilung 3+1 (N = 4) und 5+3 (N = 8) nur für gleiche Fixpunktphasen, sonst N/2 + N/2",
      ok_lat and p4 == sorted([TABLE[1]["even"], TABLE[1]["odd"]]) and p8 == sorted([TABLE[2]["even"], TABLE[2]["odd"]])
      and mixed)
say("   Die Teilung 3 + 1 verlangt, dass der Spinlift an allen vier Marken dieselbe Phase hat. Lesart: das ist die")
say("   Spinstruktur (P,P) (sqrt(omega) global); die drei geraden Spinstrukturen geben 2 + 2.")
rot_match = {}
for N, m in ((4, 1), (8, 2)):
    H = hofstadter(Lg, N)
    E, V = np.linalg.eigh(H)
    low = V[:, :N]
    U, res = lattice_symmetry(H, Lg, lambda x, y: (-y, x))
    evr = np.linalg.eigvals(low.conj().T @ U.conj().T @ low)       # Rückzug = Inverse der Verschiebung
    alg = TABLE[m]["iota"]
    rot_match[N] = any(same_multiset(evr, c * alg) for c in (1, 1j, -1, -1j))
    say(f"   N = {N}: Drehung u -> iu um die Marke im Ursprung: Gitter {fmt(sorted(evr, key=np.angle))}, algebraisch "
        f"iota^* auf K_E({m}R) {fmt(sorted(alg, key=np.angle))}; gleich bis auf eine globale Phase: {rot_match[N]}")
check("Gitterdrehung um eine Marke = iota auf den Nullmoden (bis auf globale Phase)", all(rot_match.values()))
DATA["gitter"] = {f"N={k[0]},tw={k[1][0] / np.pi:.0f}{k[1][1] / np.pi:.0f}": list(v[:2]) for k, v in lat.items()}


# ---------------------------------------------------------------------------------------------------
# 8. Was hergeleitet ist, was Eingabe bleibt, Urteil
# ---------------------------------------------------------------------------------------------------
say("\n8. Bilanz und Urteil\n")
say("   Konstruktion. Ein Multiplett: 16 Paare (psi, chi) wie in Z7, als 6D-Weylfelder auf (Querschnitt mit")
say("   Defekt) x E. Innerer Dirac auf E mit K_E^(1/2) (x) L = pi^* K(D) = K_E(R), Spinstruktur (P,P).")
say("   Der Blattwechsel ist geeicht: physikalisch ist nur, was unter sigma~ (x) J^(-1) invariant ist.")
say("   Hergeleitet in diesem Lauf:")
say("   - J^2 = -1 auf Fermionen folgt aus den vier Fixpunkten des Blattwechsels (Spinlift), ohne Viertelfluss.")
say("   - Ein Bündel statt zwei: pi^* K(D) = O(R) (x) omega, O(R)^2 = pi^* O(D). Nullmoden 3 deckgerade")
say("     (Familien, {-1, i, -i}) + 1 deckungerade (sqrt(q), Charakter i). Nur m = 1 gibt 3 + 1 (Lefschetz,")
say("     Riemann-Roch, Gitter). Die Teilung 3 + 1 verlangt gleiche Fixpunktphasen (Lesart: (P,P)).")
say("   - Die Eichbedingung koppelt J-Ladung und Decksektor: xi+ trägt nur sqrt(q), xi- nur die Familien.")
say("     16 Paare geben 16 Nahtmoden und 48 Familienkomponenten. Die Zuordnung aus Z12 (welches Multiplett")
say("     an M, welches an O(D)) ist durch den Eigenwert von J ersetzt. E8 legt die Orientierung fest.")
say("   - Träger V = Sym^2(Familienschnitte von O(R)) modulo Kegelschnitt; T ist der deckgerade Block von")
say("     f -> f dz/z (Lesart; T ist dadurch nicht eindeutig, drei Parameter).")
say("   - Verflechter: J = sigma~ = iota~^2 exakt. Die Markendrehung erfüllt rho~^2 = J tau~, tau ist eine")
say("     Halbperiodenverschiebung mit -1 auf zwei Familien. Keine Wurzel aus sigma dreht die Marken zyklisch.")
say("   Gescheitert bzw. Eingabe:")
say("   - Paarung: q = f2^2 - f0^2 exakt, das kovariante Profil ist eindeutig s ~ q. Mit der Paarung psi chi")
say("     aus G bekommen die Familien f0, f2 genau die halbe Nahtmasse; leicht bleibt nur eine Familie. Keine")
say("     Symmetrieladung von Phi trennt Naht- und Familienpaarung.")
say("   - Drei leichte Familien verlangen die benannte Kopplung Phi xi+ xi+ = Phi psi^T C P+ psi mit dem")
say("     J-Projektor P+ = (1 - iJ)/2. Sie ist eine Eingabe, die G (psi chi) widerspricht.")
say("   - Eichdarstellungen sind nicht gerechnet: Familien brauchen eine chirale 16 von SO(10), die Nahteltern")
say("     eine reelle 16 von SO(16). Ein Multiplett aus Paaren gibt xi+ und xi- dieselbe Darstellung.")
say("   - Weiterhin eingesetzt: 16 Paare, die Paarstruktur mit J (G), die Spinstruktur (P,P), die Wahl m = 1")
say("     (durch die Ziele 3 + 1 ausgewählt, nicht aus einem Prinzip).")
URTEIL = ("Lesart: Zählung (16 Naht + 48 Familien aus einem Multiplett), Charaktere, J^2 = -1, E8 und die "
          "Zuordnung Naht/Familien über den J-Eigenwert sind hergeleitet. Die Paarung bleibt eine benannte "
          "Eingabe (J-Projektor P+); mit der Paarung psi chi aus G gescheitert (eine leichte Familie, zwei mit m/2).")
say(f"   Urteil: {URTEIL}")

say("\nZusammenfassung der Prüfungen")
for n_, c_ in CHECKS.items():
    say(f"   {'bestanden' if c_ else 'NICHT bestanden'}: {n_}")
say("Befunde (negativ = ehrlich negativ, nicht erzwungen)")
for n_, c_ in BEFUNDE.items():
    say(f"   {'ja' if c_ else 'nein'}: {n_}")
(HERE / "ausgabe_z18.txt").write_text("\n".join(LINES) + "\n")
(HERE / "z18_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "urteil": URTEIL,
     "hergeleitet": ["J^2 = -1 aus den Fixpunkten des Blattwechsels",
                     "pi^* K(D) = O(R) (x) omega: 3 gerade Familien {-1,i,-i} + 1 ungerade sqrt(q); nur m = 1",
                     "Eichbedingung: xi+ <-> sqrt(q), xi- <-> Familien; 16 + 48 aus einem Multiplett",
                     "E8 (16 MW, Twistgewicht 1) und Orientierung eps = +1",
                     "J = sigma~ = iota~^2; rho~^2 = J tau~",
                     "q = f2^2 - f0^2; Familienmassen {0, 1/2, 1/2} bei Paarung psi chi"],
     "eingaben": ["Paarung nur im J-Sektor +i (P+ = (1 - iJ)/2), widerspricht psi chi aus G",
                  "16 Paare und Paarstruktur mit J (G)", "Spinstruktur (P,P)", "m = 1 durch die Ziele gewählt",
                  "Eichdarstellungen von Familien und Nahteltern (nicht gerechnet)"],
     "data": DATA}, ensure_ascii=False, indent=1, default=str))
assert all(CHECKS.values())
