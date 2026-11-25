"""Z89: Gravitation über die Komposition (Folgeauftrag F1): Infinity-Twistor und Einstein-Vertex.

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (die Minkowskizeit ist der parabolische Erzeuger,
nicht der Modulfluss), GRAV.PARENT.GAP.NO_GO.01 (ein gleichmäßig gegappter Elternoperator schließt einen masselosen Spin-2-Pol
aus; T7 braucht eine kritische Folge mit lokalen Tensor-Eichbedingungen), DIMENSION.UPLIFT.FIREWALL.01 (alles ist Kinematik bzw.
Twistorklammer auf Elementarzuständen; keine wechselwirkende 3+1-Schließung). Gelöste Blöcke (Flavour, Massen, alpha, 16, E8)
bleiben unberührt; nichts wird befördert. Bestand: z78 (Araki-Woods, P_mu = c^+ sigma_mu c/2 mit c = d + e^+ = kappa (a + b^+)),
z79 (T = span{a_i, b_i^+}, tau = (Ad J)^+, pi- und omega-Ebene, Kommutant gl(2,C)), z83 (radiale su(2,2), Moden = Stufen von E,
parabolische Minkowskienergie), z85 (Fockbasis = Elementarzustände, Grad 2 = Sektor D = -4, native zweite Ordnung Heisenberg
statt Plebanski-Penrose, D = +-4 gleichmäßig gegappt), z88 (viele Teilchen brauchen viele Nähte). Geprüft werden die zehnte
externe Analyse (Auswahlregel D(-4) x D(-4) -> D(-4), Energieschranke 6 > 3, Bedarf an I^(alpha beta)) und die Kandidaten der
elften (Zwei-Naht-Bivektor Z_i ^ Z_j, diagonale E8-Komposition mit Mischterm).
Statuskürzel: E exakt, C bedingt, L Lesart, O offen, N numerisch.

Teil 1  Auswahlregel und Energieschranke: Modulladung der Kompositionen (Heisenberg mit hermiteschem Kozyklus, Normalprodukt,
        Produkt der Elementarzustände, Klammer aus z85); tiefste Gewichte und kein su(2,2)-Verflechter H x H -> H auf
        Abschneidungen, Kontrollen Yang-Mills (D = -2) und Skalar; Derivationsalgebra der Klammer in u(2,2), invariante
        Bivektoren, Eulergewicht.
Teil 2  Native Kandidaten für I und ihre Klammer auf dem Grad-2-Sektor D = -4 im Fockwörterbuch: (A) eps auf der pi-Ebene,
        (B) eps auf der omega-Ebene, (C) Uhr, Ladung und Vertexoperatoren, (D) zwei Nähte relativ über der kommutativen
        c-Algebra und die Zerlegung der Klammer in zwei Nähte, (E) Zwei-Naht-Bivektor Z_i ^ Z_j aus Zuständen zweier Nähte
        (elfte Analyse), (F) diagonale E8-Komposition mit dem Mischterm Summe J_1^a J_2^a.
Teil 3  Skala, Büschel und Universalität: Lorentzbüschel lambda I_pi + mu I_omega (flach, AdS, dS), der radiale Kanal ist AdS;
        Normierungen (Zustand, Gitter, kappa^2 = tanh(pi m)); Ladungssektoren, Vertexoperatoren, Helizitäten; kritische Folge
        der Minkowskizeit gegen E und den Modulgenerator; was G wäre.
Teil 4  Bilanz und Urteil für F1.
"""
import itertools
import json
import math
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
rng = np.random.default_rng(89)
mp.mp.dps = 30
I2 = np.eye(2)
SIGMA = [np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]),
         np.diag([1., -1.]).astype(complex)]
KAP2_NATIV = {m: float(mp.tanh(mp.pi*m)) for m in (1, 2, 3)}             # kappa^2 = (1 - q)/(1 + q), q = e^{-2 pi m}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")
    return bool(cond)


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def laufzeit():
    say(f"   (Laufzeit bisher {time.time() - T0:.1f} s)")


def nullraum(M, tol=1e-9):
    M = np.atleast_2d(M)
    _, s, vh = np.linalg.svd(M)
    rank = int(np.sum(s > tol*max(1.0, s[0]))) if s.size else 0
    return vh[rank:]


# ------------------------------------------------------------------ Twistorfunktionen und Fockwörterbuch (z85)
Z = sp.symbols("z1:5")
KAP = sp.Symbol("kappa", positive=True)
A_ = sp.symbols("a1 a2")
B_ = sp.symbols("b1 b2")
U_PI = [sp.Matrix([1, 0, -1, 0]), sp.Matrix([0, 1, 0, -1])]      # [c_j^+, f] = -(u_j . grad) f; span{u} = {pi = 0}
V_OM = [sp.Matrix([1, 0, 1, 0]), sp.Matrix([0, 1, 0, 1])]       # [omega''_j^+, f] = -(v_j . grad) f; span{v} = {omega = 0}


def bivektor(r1, r2):
    return r1*r2.T - r2*r1.T


I_PI, I_OM = bivektor(*U_PI), bivektor(*V_OM)


def grad(f):
    return sp.Matrix([sp.diff(f, z) for z in Z])


def klammer(I, f, g):
    """B_I(f, g) = I^(alpha beta) d_alpha f d_beta g (Poissonklammer mit konstantem Bivektor I)."""
    return sp.expand((grad(f).T*I*grad(g))[0, 0])


def norm_fak(N):
    n1, n2, m1, m2 = N
    return sp.sqrt(sp.Rational(math.factorial(n1)*math.factorial(n2), math.factorial(m1)*math.factorial(m2)))


def elementar(N):
    return norm_fak(N)*Z[2]**N[2]*Z[3]**N[3]/(Z[0]**(N[0] + 1)*Z[1]**(N[1] + 1))


def fock_von(expr):
    """Laurentfunktion -> Fockvektor {(n1, n2, m1, m2): c} über |n1 n2; m1 m2> <-> sqrt(n1! n2!/(m1! m2!)) z3^m1 z4^m2 /
    (z1^(n1+1) z2^(n2+1)) (z85); Monome ohne Pol in z1 oder z2 sind Koränder und fallen weg."""
    out = {}
    for t in sp.Add.make_args(sp.expand(expr)):
        if t == 0:
            continue
        c, ex = t, []
        for z in Z:
            c, e = c.as_coeff_exponent(z)
            ex.append(int(e))
        if ex[0] >= 0 or ex[1] >= 0:
            continue
        N = (-ex[0] - 1, -ex[1] - 1, ex[2], ex[3])
        if N[2] < 0 or N[3] < 0:
            raise ValueError(f"kein Elementarzustand: {t}")
        out[N] = sp.expand(out.get(N, 0) + c/norm_fak(N))
    return {k: v for k, v in out.items() if v != 0}


def funktion_von(vec):
    return sp.expand(sum((c*elementar(N) for N, c in vec.items()), sp.Integer(0)))


def D_von(N):
    return N[0] + N[1] - N[2] - N[3]


def E_von(N):
    return sp.Rational(sum(N), 2) + 1


def c_plus(j, vec):
    """c_j^+ = d_j^+ + e_j auf Fockvektoren (exakt)."""
    out = {}
    for N, c in vec.items():
        n = list(N)
        M1 = n.copy()
        M1[j] += 1
        out[tuple(M1)] = out.get(tuple(M1), 0) + c*sp.sqrt(n[j] + 1)
        if n[j + 2] > 0:
            M2 = n.copy()
            M2[j + 2] -= 1
            out[tuple(M2)] = out.get(tuple(M2), 0) + c*sp.sqrt(n[j + 2])
    return out


def normalprodukt(vx, vy):
    """Produkt der Erzeugermonome: |N> (x) |M> -> prod sqrt(C(N+M, N)) |N + M> (Zustandsform der Heisenberg-Komposition)."""
    out = {}
    for N, c in vx.items():
        for M, d in vy.items():
            S = tuple(a + b for a, b in zip(N, M))
            out[S] = out.get(S, 0) + c*d*sp.prod([sp.sqrt(sp.binomial(S[k], N[k])) for k in range(4)])
    return out


# ------------------------------------------------------------------ Fockraum einer Naht (vier Oszillatoren, wie z85)
def fock4(n):
    a = sps.diags(np.sqrt(np.arange(1, n)), 1, format="csr")
    Id = sps.identity(n, format="csr")
    ops = []
    for k in range(4):
        mats = [Id]*4
        mats[k] = a
        op = mats[0]
        for Mx in mats[1:]:
            op = sps.kron(op, Mx, format="csr")
        ops.append(op)
    return ops


def dag(A):
    return A.T.conj().tocsr()


def komm(A, B):
    return (A @ B - B @ A).tocsr()


NV = 8
DIMV = NV**4
d1, d2, e1, e2 = fock4(NV)
DD, EE = [d1, d2], [e1, e2]
OCC = np.array(list(itertools.product(range(NV), repeat=4)))
DVAL = OCC[:, 0] + OCC[:, 1] - OCC[:, 2] - OCC[:, 3]
LOW = np.where(OCC.sum(1) <= NV - 3)[0]
OMV = np.zeros(DIMV, complex)
OMV[0] = 1


def klein(A, tol=1e-10):
    A = sps.csc_matrix(A)[:, LOW]
    return (abs(A).max() if A.nnz else 0.0) < tol


qt = 0.3                                                          # Prüfwert; nativ q = e^{-2 pi m}
st = math.sqrt(1 - qt**2)
kt = math.sqrt((1 - qt)/(1 + qt))
AA = [((DD[i] + qt*dag(EE[i]))/st).tocsr() for i in range(2)]     # Stromoszillatoren des Originals
BB = [((EE[i] + qt*dag(DD[i]))/st).tocsr() for i in range(2)]     # Tomita-Partner
ZT = [d1, d2, dag(e1), dag(e2)]                                   # Z^alpha
CC = [(DD[i] + dag(EE[i])).tocsr() for i in range(2)]             # c = d + e^+ = kappa (a + b^+)
WO = [(DD[i] - dag(EE[i])).tocsr() for i in range(2)]             # omega'' = d - e^+ = (a - b^+)/kappa


def U_theta(th):
    return sps.diags(np.exp(1j*th*DVAL)).tocsr()


def gen(alpha):
    return sum((alpha[i]*dag(AA[i]) - np.conj(alpha[i])*AA[i] for i in range(2)), sps.csr_matrix((DIMV, DIMV))).tocsc()


def W_auf(alpha, v):
    return spla.expm_multiply(gen(alpha), v)


def erwartung(v, X):
    return complex(np.vdot(v, X @ v))


# ------------------------------------------------------------------ u(2,2) auf dem Twistorraum (z-Koordinaten)
G4 = np.diag([1., 1., -1., -1.])
_rows = []
for _k in range(32):
    _v = np.zeros(32)
    _v[_k] = 1
    _M = (_v[:16] + 1j*_v[16:]).reshape(4, 4)
    _C = _M.conj().T @ G4 + G4 @ _M
    _rows.append(np.concatenate([_C.real.ravel(), _C.imag.ravel()]))
U22 = [(v[:16] + 1j*v[16:]).reshape(4, 4) for v in nullraum(np.array(_rows).T)]


def in_u22(M):
    return np.abs(M.conj().T @ G4 + G4 @ M).max() < 1e-12


def teilalgebra(bed, basis=None):
    basis = U22 if basis is None else basis
    rows = []
    for M in basis:
        parts = []
        for f in bed:
            C = f(M)
            parts += [C.real.ravel(), C.imag.ravel()]
        rows.append(np.concatenate(parts))
    ns = nullraum(np.array(rows).T)
    return [sum(c*M for c, M in zip(v, basis)) for v in ns]


def derivationen(I):
    return teilalgebra([lambda M: M @ I + I @ M.T])


def im_spann(M, basis):
    if not basis:
        return np.abs(M).max() < 1e-12
    A = np.array([np.concatenate([B.real.ravel(), B.imag.ravel()]) for B in basis]).T
    b = np.concatenate([M.real.ravel(), M.imag.ravel()])
    x = np.linalg.lstsq(A, b, rcond=None)[0]
    return np.linalg.norm(A @ x - b) < 1e-9*max(1.0, np.linalg.norm(b))


ANTI = []
for _i, _j in itertools.combinations(range(4), 2):
    for _ph in (1, 1j):
        _X = np.zeros((4, 4), complex)
        _X[_i, _j], _X[_j, _i] = _ph, -_ph
        ANTI.append(_X)


def inv_bivektoren(gens):
    rows = []
    for X in ANTI:
        parts = []
        for M in gens:
            C = M @ X + X @ M.T
            parts += [C.real.ravel(), C.imag.ravel()]
        rows.append(np.concatenate(parts))
    ns = nullraum(np.array(rows).T)
    return [sum(c*X for c, X in zip(v, ANTI)) for v in ns]


def signatur(basis):
    Gm = np.array([[np.trace(X @ Y).real for Y in basis] for X in basis])
    ew = np.linalg.eigvalsh((Gm + Gm.T)/2)
    sk = max(1.0, np.abs(ew).max())
    return int(np.sum(ew > 1e-9*sk)), int(np.sum(ew < -1e-9*sk)), int(np.sum(np.abs(ew) <= 1e-9*sk))


def pfaff(I):
    return I[0, 1]*I[2, 3] - I[0, 2]*I[1, 3] + I[0, 3]*I[1, 2]


uvec = [np.array([1, 0, -1, 0.]), np.array([0, 1, 0, -1.])]
vvec = [np.array([1, 0, 1, 0.]), np.array([0, 1, 0, 1.])]
Inp = np.array(I_PI.tolist(), dtype=complex)
Ion = np.array(I_OM.tolist(), dtype=complex)
S4 = np.block([[np.zeros((2, 2)), I2], [I2, np.zeros((2, 2))]]).astype(complex)   # tau = (Ad J)^+ auf Z: d <-> e^+
NP = [-0.5*sum(s[i, j]*np.outer(uvec[i], vvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]  # [P_mu, Z] = N Z
NK = [-0.5*sum(s[i, j]*np.outer(vvec[i], uvec[j]) for i in range(2) for j in range(2)) for s in SIGMA]  # omega''-Translationen
NE = np.diag([-0.5, -0.5, 0.5, 0.5]).astype(complex)                 # [E, Z] = NE Z, E = (N_d + N_e)/2 + 1
MP = [1j*N for N in NP]
MK = [1j*N for N in NK]
MEt, MDz, MDil = 1j*NE, -1j*np.eye(4, dtype=complex), S4.copy()      # konforme Energie, Zentrum (Modulfluss), Dilatation

# =================================================================== Teil 1
say("Teil 1: Auswahlregel und Energieschranke (zehnte Analyse, nachgerechnet)")

# --- 1.1 Modulladung: Heisenberg-Kozyklus (hermitesch) gegen eps-Form (komplex)
th = 0.7
U = U_theta(th)
ok_ph = all(klein(U @ X @ dag(U) - np.exp(-1j*th)*X, 1e-12) for X in ZT + AA + [dag(b) for b in BB] + CC + WO)
al = (rng.normal(size=2) + 1j*rng.normal(size=2))*0.15
be = (rng.normal(size=2) + 1j*rng.normal(size=2))*0.15
WWv = W_auf(al, W_auf(be, OMV))
res_kov = float(np.linalg.norm(U @ WWv - W_auf(np.exp(1j*th)*al, W_auf(np.exp(1j*th)*be, OMV))))
koz = float(np.imag(al @ np.conj(be)))
koz_rot = float(np.imag((np.exp(1j*th)*al) @ np.conj(np.exp(1j*th)*be)))
res_weyl = float(np.linalg.norm(WWv - np.exp(1j*koz)*W_auf(al + be, OMV)))
va, vb = W_auf(al, OMV), W_auf(be, OMV)
dca = np.array([erwartung(va, CC[j]) for j in range(2)])
dcb = np.array([erwartung(vb, CC[j]) for j in range(2)])


def eps2(x, y):
    return x[0]*y[1] - x[1]*y[0]


e0 = eps2(dca, dcb)
ok_dc = np.allclose(dca, kt*al, atol=1e-6) and np.allclose(dcb, kt*be, atol=1e-6)
ok_eps = (abs(eps2(kt*np.exp(1j*th)*al, kt*np.exp(1j*th)*be) - np.exp(2j*th)*eps2(kt*al, kt*be)) < 1e-14
          and abs(e0 - eps2(kt*al, kt*be)) < 1e-6 and abs(e0) > 1e-4)
ZAHLEN["modulladung"] = {"Residuum Kovarianz W(a)W(b)": res_kov, "Kozyklus vor/nach Drehung": [koz, koz_rot],
                         "Weyl-Residuum": res_weyl, "eps(delta c_a, delta c_b)": str(np.round(e0, 10))}
check(f"Modulladung: U_theta = e^(i theta D) (D = N_a - N_b) dreht T = span{{a_i, b_i^+}} = span{{Z^alpha}} samt c und omega'' "
      f"um die gemeinsame Phase e^(-i theta) (exakt auf 8^4). Die native zweite Ordnung ist kovariant: U W(a) W(b) Omega = "
      f"W(e^(i theta) a) W(e^(i theta) b) Omega (Residuum {res_kov:.1e}); ihr Kozyklus Im a.b-quer ist die hermitesche "
      f"CCR-Form und phaseninvariant ({koz:.6f} -> {koz_rot:.6f}; Weyl-Relation {res_weyl:.1e}). Die komplexe eps-Form auf "
      f"den pi-Komponenten derselben Verschiebungen (delta c = <c> = kappa alpha, numerisch) trägt die Ladung 2: "
      f"eps -> e^(2 i theta) eps", ok_ph and res_kov < 1e-10 and abs(koz - koz_rot) < 1e-14 and res_weyl < 1e-5
      and ok_dc and ok_eps)

# --- 1.2 D-Verschiebung der Bilinearabbildungen auf Elementarzuständen
PROBEN = {-4: [(0, 0, 4, 0), (0, 0, 1, 3), (1, 0, 3, 2), (0, 1, 5, 0)], -2: [(0, 0, 2, 0), (0, 0, 1, 1), (1, 0, 1, 2)],
          0: [(0, 0, 0, 0), (1, 0, 0, 1), (0, 1, 1, 0)], 2: [(2, 0, 0, 0), (1, 1, 0, 0), (2, 1, 1, 0)],
          4: [(4, 0, 0, 0), (2, 2, 0, 0)]}
assert all(D_von(N) == D for D, L in PROBEN.items() for N in L)
ABB = {
    "Normalprodukt (Zustandsform von Heisenberg und Zwei-Naht-Tensorprodukt)": lambda N, M: normalprodukt({N: 1}, {M: 1}),
    "Produkt der Elementarzustände (Komposition des Yang-Mills-Vertex b[a,a])":
        lambda N, M: fock_von(elementar(N)*elementar(M)),
    "Klammer B_I mit I = eps auf der pi-Ebene (z85)": lambda N, M: fock_von(klammer(I_PI, elementar(N), elementar(M))),
    "Klammer B_I mit I_AdS = I_pi - I_omega (Teil 3)": lambda N, M: fock_von(klammer(I_PI - I_OM, elementar(N), elementar(M))),
}
SEKTORPAARE = [(-4, -4), (-2, -2), (0, 0), (-4, 0), (-4, 4), (-2, 2), (-4, -2)]
verschiebung = {}
for name, f in ABB.items():
    dsh, esh, n_nz = set(), set(), 0
    for s1, s2 in SEKTORPAARE:
        for N in PROBEN[s1]:
            for M in PROBEN[s2]:
                for K, c in f(N, M).items():
                    if c != 0:
                        dsh.add(D_von(K) - s1 - s2)
                        esh.add(E_von(K) - E_von(N) - E_von(M))
                        n_nz += 1
    verschiebung[name] = {"D-Verschiebung": sorted(dsh), "E-Verschiebung": sorted(esh), "Komponenten": n_nz}
ZAHLEN["verschiebungen"] = {k: {kk: [str(x) for x in vv] if isinstance(vv, list) else vv for kk, vv in v.items()}
                            for k, v in verschiebung.items()}
vl = list(verschiebung.values())
for k, v in verschiebung.items():
    say(f"   {k}: D-Verschiebung {v['D-Verschiebung']}, E-Verschiebung {[str(x) for x in v['E-Verschiebung']]}, "
        f"{v['Komponenten']} Komponenten")
check("D-Verschiebung auf Elementarzuständen (exakt, sieben Sektorpaare): Normalprodukt 0, Produkt der Elementarzustände "
      "(die Komposition, die der kubische Yang-Mills-Vertex b[a,a] braucht) +2, Klammer +4. Modulkovariant (U_theta-äquivariant) "
      "ist eine Bilinearabbildung genau bei Verschiebung 0. Bestätigt: Keine native Komposition bildet D(-4) x D(-4) -> D(-4) "
      "(nötig +4). Ebenso wenig aber D(-2) x D(-2) -> D(-2) (nötig +2, Yang-Mills): Die nötige Verschiebung ist "
      "2 (h_3 - h_1 - h_2), und masselose Dreipunktvertizes haben stets eine Helizitätssumme != 0; die Auswahlregel trennt "
      "Gravitation also nicht von Yang-Mills, gravitationsspezifisch ist nur die Differenz +2",
      vl[0]["D-Verschiebung"] == [0] and vl[1]["D-Verschiebung"] == [2] and vl[2]["D-Verschiebung"] == [4]
      and vl[3]["D-Verschiebung"] == [4] and all(v["Komponenten"] > 0 for v in vl))


# --- 1.3 Energieschranke: tiefste Gewichte (Kern der Senker K_ij = d_i e_j) je Sektor
def senker_kerne(D, Emax):
    kern, dims = {}, {}
    n, prev = max(D, 0), None
    while n - D/2 + 1 <= Emax + 1e-9:
        mN = n - D
        basis = [(n1, n - n1, m1, mN - m1) for n1 in range(n + 1) for m1 in range(mN + 1)]
        if prev is None:
            k = len(basis)
        else:
            pidx = {s: i for i, s in enumerate(prev)}
            M = np.zeros((4*len(prev), len(basis)))
            for c, s in enumerate(basis):
                for i in range(2):
                    for j in range(2):
                        if s[i] > 0 and s[2 + j] > 0:
                            t = list(s)
                            t[i] -= 1
                            t[2 + j] -= 1
                            M[(2*i + j)*len(prev) + pidx[tuple(t)], c] = math.sqrt(s[i]*s[2 + j])
            k = len(basis) - np.linalg.matrix_rank(M)
        E = n - D/2 + 1
        kern[E], dims[E] = int(k), len(basis)
        prev = basis
        n += 1
    return kern, dims


ener = {}
ok_en = True
for D in (-4, -2, 0):
    kern, dims = senker_kerne(D, 9)
    Emin = 1 + abs(D)/2
    tw = {E: k for E, k in kern.items() if k}
    ener[D] = {"tiefste Gewichte {E: Anzahl}": tw, "Dimensionen je E": dims, "E_min": Emin,
               "E_min von H x H": 2*Emin, "Tiefstgewichte von H x H bei 2 E_min": (abs(D) + 1)**2,
               "Kern der Senker von H bei 2 E_min": kern.get(2*Emin)}
    ok_en &= tw == {Emin: abs(D) + 1} and kern.get(2*Emin) == 0 and kern.get(2*Emin + 1) == 0
ZAHLEN["energieschranke"] = {str(k): v for k, v in ener.items()}
say(f"   tiefste Gewichte je Sektor (E <= 9): D = -4: {ener[-4]['tiefste Gewichte {E: Anzahl}']}, D = -2: "
    f"{ener[-2]['tiefste Gewichte {E: Anzahl}']}, D = 0: {ener[0]['tiefste Gewichte {E: Anzahl}']}")
check("Energieschranke (exakt auf Abschneidungen E <= 9): Im Sektor D = -4 gibt es tiefste Gewichte nur bei E = 3 (fünf, "
      "K-Typ (0, 2)); H x H beginnt bei E = 6 mit 25 Tiefstgewichten, und der Kern der Senker von H ist bei E = 6 und 7 null. "
      "Ein su(2,2)-Verflechter T: H x H -> H bildet Tiefstgewichte auf Tiefstgewichte gleicher Energie ab, also T = 0 Stufe "
      "für Stufe: Hom(D(-4) x D(-4), D(-4)) = 0, wie behauptet. Dieselbe Rechnung gibt Hom = 0 für Yang-Mills (D = -2: E = 2, "
      "Produkt ab 4) und den Skalar (E = 1, Produkt ab 2): Die Schranke schließt jede konform äquivariante Fusion "
      "positiver Einteilchendarstellungen aus, nicht nur die gravitative", ok_en)

# --- 1.4 Derivationsalgebren in u(2,2), invariante Bivektoren, Eulergewicht
Pop = [0.5*sum(s[i, j]*(dag(CC[i]) @ CC[j]) for i in range(2) for j in range(2)) for s in SIGMA]
Kop = [0.5*sum(s[i, j]*(dag(WO[i]) @ WO[j]) for i in range(2) for j in range(2)) for s in SIGMA]
Eop = 0.5*sum((dag(X) @ X for X in DD + EE), sps.csr_matrix((DIMV, DIMV))) + sps.identity(DIMV)
ok_N = all(klein(komm(Pop[mu], ZT[a_]) - sum(NP[mu][a_, b_]*ZT[b_] for b_ in range(4))) and
           klein(komm(Kop[mu], ZT[a_]) - sum(NK[mu][a_, b_]*ZT[b_] for b_ in range(4))) for mu in range(4) for a_ in range(4))
ok_N &= all(klein(komm(Eop, ZT[a_]) - NE[a_, a_]*ZT[a_]) for a_ in range(4))
ok_N &= all(in_u22(M) for M in MP + MK + [MEt, MDz, MDil])
DER_PI, DER_OM = derivationen(Inp), derivationen(Ion)
KOMM_TAU = teilalgebra([lambda M: M @ S4 - S4 @ M])
LOR = teilalgebra([lambda M: M @ S4 - S4 @ M, lambda M: M @ Inp + Inp @ M.T])
SCHNITT = teilalgebra([lambda M: M @ Inp + Inp @ M.T, lambda M: M @ Ion + Ion @ M.T])
inv_u22 = inv_bivektoren(U22)
inv_poin = inv_bivektoren(DER_PI)
inv_lor = inv_bivektoren(LOR)
ok_der = (len(U22) == 16 and len(DER_PI) == 10 and len(DER_OM) == 10 and len(KOMM_TAU) == 8 and len(LOR) == 6
          and len(SCHNITT) == 6 and all(im_spann(M, DER_PI) for M in MP) and all(im_spann(M, DER_OM) for M in MK)
          and not any(im_spann(M, DER_PI) for M in MK + [MEt, MDz, MDil])
          and not any(im_spann(M, DER_OM) for M in MP) and len(inv_u22) == 0 and len(inv_poin) == 2
          and all(im_spann(X, [Inp, 1j*Inp]) for X in inv_poin) and len(inv_lor) == 4
          and all(im_spann(X, [Inp, 1j*Inp, Ion, 1j*Ion]) for X in inv_lor))
sig_poin = signatur(DER_PI)
ZAHLEN["derivationen"] = {"dim u(2,2)": len(U22), "Der(I_pi)": len(DER_PI), "Der(I_omega)": len(DER_OM),
                          "Kommutant tau": len(KOMM_TAU), "Lorentz = Kommutant tau n Der(I_pi)": len(LOR),
                          "Der(I_pi) n Der(I_omega)": len(SCHNITT), "u(2,2)-invariante Bivektoren": len(inv_u22),
                          "Poincaré-invariante Bivektoren (reell)": len(inv_poin), "Lorentz-invariante (reell)": len(inv_lor),
                          "Spurform-Signatur Der(I_pi) (+, -, 0)": sig_poin}
check(f"Derivationsalgebren (u(2,2) = {len(U22)} reelle Dimensionen, G = diag(1,1,-1,-1)): [P_mu, Z] = N_mu Z für die "
      f"nativen P_mu = (1/2) c^+ sigma_mu c, ebenso für die omega''-Translationen und E (Fock 8^4). Die Klammer mit I_pi "
      f"(eps auf der pi-Ebene) hat die Derivationen Poincaré: Dimension {len(DER_PI)}, alle vier P_mu, Lorentz = Kommutant "
      f"von tau (Dimension {len(KOMM_TAU)} = gl(2,C)) geschnitten mit Der ({len(LOR)}), Spurform {sig_poin} (entartet: "
      f"Translationen nilpotent); nicht E, nicht die Dilatation, nicht das Zentrum (Modulfluss), keine K_mu. I_omega hat die "
      f"gespiegelte Poincaréalgebra (K_mu statt P_mu), Schnitt {len(SCHNITT)} = Lorentz. Kein Bivektor ist u(2,2)-invariant "
      f"({len(inv_u22)}); Poincaré-invariant ist genau C I_pi ({len(inv_poin)} reell), Lorentz-invariant genau "
      f"C I_pi + C I_omega ({len(inv_lor)} reell)", ok_N and ok_der)
fs = [elementar((0, 0, 4, 0)), elementar((1, 0, 1, 2)), elementar((0, 1, 0, 1))]
eu = lambda f: sp.expand(sum(z*sp.diff(f, z) for z in Z))            # Euler-Feld X_1 = Z . grad
ok_eu = all(sp.expand(eu(klammer(I_PI, f, g)) - klammer(I_PI, eu(f), g) - klammer(I_PI, f, eu(g))
                      + 2*klammer(I_PI, f, g)) == 0 for f, g in itertools.combinations(fs, 2))
Mr = sp.Matrix(4, 4, [sp.Rational(int(x), 7) for x in rng.integers(-6, 7, size=16)])
XM = lambda f: sp.expand(sum(Mr[a_, b_]*Z[b_]*sp.diff(f, Z[a_]) for a_ in range(4) for b_ in range(4)))
ok_lb = all(sp.expand(XM(f*g) - XM(f)*g - f*XM(g)) == 0 for f, g in itertools.combinations(fs, 2))
check("Eulergewicht: X_1 B_I(f, g) - B_I(X_1 f, g) - B_I(f, X_1 g) = -2 B_I(f, g) (X_1 = Z.grad, Zentrum geometrisch), "
      "während das Produkt der Elementarzustände für jedes lineare Vektorfeld (hier eine zufällige rationale Matrix) die "
      "Leibnizregel erfüllt, also die volle u(2,2) (16) als Derivationen hat: Der Bivektor I bricht su(2,2) auf Poincaré und "
      "trägt gegenüber dem Yang-Mills-Produkt genau die Zusatzladung (Gewicht -2, D-Verschiebung +2)", ok_eu and ok_lb)
befund("E: Die zehnte Analyse stimmt in beiden Aussagen (Auswahlregel, Energieschranke 6 > 3), beide sind aber nicht "
       "gravitationsspezifisch: Sie schließen jede kubische Kopplung masseloser Felder als modulkovariante bzw. konform "
       "äquivariante Abbildung H x H -> H aus, Yang-Mills eingeschlossen. Gravitationsspezifisch ist allein der Bivektor I: "
       "Derivationen 10 (Poincaré) statt 16, Zusatzladung +2 gegenüber dem konformen Yang-Mills-Produkt", True)
laufzeit()

# =================================================================== Teil 2
say("\nTeil 2: Native Kandidaten für I und ihre Klammer auf dem Grad-2-Sektor D = -4")
ha = (A_[0]*Z[2] + A_[1]*Z[3])**4/(Z[0]*Z[1])
hb = (B_[0]*Z[2] + B_[1]*Z[3])**4/(Z[0]*Z[1])
Af, Bf = A_[0]*Z[2] + A_[1]*Z[3], B_[0]*Z[2] + B_[1]*Z[3]
ab = A_[0]*B_[1] - A_[1]*B_[0]
FORM_PI = Z[2]/(Z[0]**3*Z[1]**2) + Z[3]/(Z[0]**2*Z[1]**3) + 4/(Z[0]**2*Z[1]**2)


def d_om(f, j):                                                   # z85: d/d omega_j bei festem pi (CCR-Normierung)
    return (sp.diff(f, Z[j]) - sp.diff(f, Z[j + 2]))/(sp.sqrt(2)*KAP)


kl_z85 = sp.expand(d_om(ha, 0)*d_om(hb, 1) - d_om(ha, 1)*d_om(hb, 0))
br_pi = klammer(I_PI, ha, hb)
br_om = klammer(I_OM, ha, hb)
SUB = {A_[0]: 1, A_[1]: 0, B_[0]: 0, B_[1]: 1}
v_pi, v_om = fock_von(br_pi.subs(SUB)), fock_von(br_om.subs(SUB))

# --- (A) eps auf der pi-Ebene
ok_A1 = (sp.expand(kl_z85 - 2/KAP**2*ab*Af**3*Bf**3*FORM_PI) == 0 and sp.expand(br_pi - 2*KAP**2*kl_z85) == 0
         and sp.expand(br_pi - 4*ab*Af**3*Bf**3*FORM_PI) == 0)
ha1, hb1 = ha.subs(SUB), hb.subs(SUB)
norm2_pi = sum(abs(complex(c))**2 for c in v_pi.values())
norm_in = float(sp.sqrt(sum(c**2 for c in fock_von(ha1).values())))
ZAHLEN["klammer_pi"] = {"B_I_pi(h_a, h_b)": str(sp.factor(br_pi)), "= 2 kappa^2 mal z85": True,
                        "Fock (a = (1,0), b = (0,1))": {str(k): str(v) for k, v in v_pi.items()},
                        "D": sorted({D_von(k) for k in v_pi}), "E": sorted({str(E_von(k)) for k in v_pi}),
                        "Norm für normierte Eingänge": math.sqrt(norm2_pi)/norm_in**2}
check(f"(A) eps auf der pi-Ebene: B_I_pi(h_a, h_b) = 2 kappa^2 {{h_a, h_b}}_z85 = 4 [a b] A^3 B^3 (z3/(z1^3 z2^2) + "
      f"z4/(z1^2 z2^3) + 4/(z1^2 z2^2)) (symbolisch exakt; die z85-Form (2/kappa^2) [a b] ... ist reproduziert): genau die "
      f"Plebanski-Penrose-Klammer, in Fockkoordinaten ohne kappa. Fockbild für a = (1,0), b = (0,1): "
      f"{ {str(k): str(v) for k, v in v_pi.items()} }, alles D = -4 bei E = 5 und 6 (nicht E-äquivariant, E liegt nicht in "
      f"Poincaré); Norm für normierte Eingänge {math.sqrt(norm2_pi)/norm_in**2:.6f} = sqrt(20)",
      ok_A1 and {D_von(k) for k in v_pi} == {-4} and {E_von(k) for k in v_pi} == {5, 6}
      and abs(math.sqrt(norm2_pi)/norm_in**2 - math.sqrt(20)) < 1e-12)
# Nativität: tau-Eigenraum, Zustandsnormierung, Phase
SW = sps.csr_matrix((np.ones(DIMV), ([np.ravel_multi_index((o[2], o[3], o[0], o[1]), (NV,)*4) for o in OCC],
                                     np.arange(DIMV))), shape=(DIMV, DIMV))
tau = lambda X: dag(SW @ X @ SW)                                  # tau(A) = (J A J)^+, J = Vertauschung d <-> e (reell)
ok_tau = all(klein(tau(CC[j]) - CC[j]) and klein(tau(WO[j]) + WO[j]) for j in range(2))
PIo = [((AA[i] + dag(BB[i]))/np.sqrt(2)).tocsr() for i in range(2)]
Gc = np.array([[np.vdot(CC[i] @ OMV, CC[j] @ OMV) for j in range(2)] for i in range(2)])
Gw = np.array([[np.vdot(WO[i] @ OMV, WO[j] @ OMV) for j in range(2)] for i in range(2)])
Gpi = np.array([[np.vdot(PIo[i] @ OMV, PIo[j] @ OMV) for j in range(2)] for i in range(2)])
ok_cpi = all(klein(CC[i] - np.sqrt(2)*kt*PIo[i]) for i in range(2))
ok_ccomm = all(klein(komm(CC[i], dag(CC[j]))) and klein(komm(CC[i], CC[j])) and klein(komm(WO[i], dag(WO[j])))
               for i in range(2) for j in range(2))
P0_id = klein(Pop[0] - Eop - 0.5*sum((dag(DD[i]) @ dag(EE[i]) + EE[i] @ DD[i] for i in range(2)),
                                     sps.csr_matrix((DIMV, DIMV))))
K0_id = klein(Kop[0] - Eop + 0.5*sum((dag(DD[i]) @ dag(EE[i]) + EE[i] @ DD[i] for i in range(2)),
                                     sps.csr_matrix((DIMV, DIMV))))
ths = 2*np.pi*np.arange(8)/8
mittel_phase = abs(np.mean(np.exp(-2j*ths)))
ZAHLEN["zustandsnormierung"] = {"<c^+ c>": np.round(Gc.real, 12).tolist(), "<omega''^+ omega''>": np.round(Gw.real, 12).tolist(),
                                "<pi^+ pi> (CCR-normiert)": np.round(Gpi.real, 12).tolist(), "1/(2 kappa^2)": 1/(2*kt**2),
                                "|I_Zustand|/|I_CCR| = 2 kappa^2 (nativ m = 1)": 2*KAP2_NATIV[1]}
check(f"(A) Nativität: die pi-Ebene ist der Eigenraum tau = +1 (c = d + e^+ = sqrt2 kappa pi, tau(c) = c, tau(omega'') = "
      f"-omega'', Fock 8^4), c und c^+ vertauschen (kommutative c-Algebra). Der KMS-Zustand normiert sie: <c_i^+ c_j> = "
      f"delta_ij, <pi^+ pi> = 1/(2 kappa^2) (q = 0,3), also |eps(c_1, c_2)| = 1 und I_Zustand = I_pi = 2 kappa^2 I_CCR "
      f"(nativ 2 tanh(pi) = {2*KAP2_NATIV[1]:.6f}). Die Zeitkomponente der nativen Translationen hat genau die parabolische "
      f"Form von z83: (1/2) c^+ c = E + (P^0 + K^0)/2, hier mit P^0 = Summe d_i^+ e_i^+ (Operatoridentität), die "
      f"omega-Ebene gibt das Inversionsbild (1/2) omega''^+ omega'' = E - (P^0 + K^0)/2. Die Phase von eps ist die "
      f"Modulphase: eps -> e^(-2 i theta) eps, Mittel über eine Periode {mittel_phase:.1e}, kein modulinvariantes Element "
      f"!= 0", ok_tau and ok_cpi and ok_ccomm and np.allclose(Gc, I2, atol=1e-12) and np.allclose(Gw, I2, atol=1e-12)
      and np.allclose(Gpi, I2/(2*kt**2), atol=1e-10) and P0_id and K0_id and mittel_phase < 1e-12)
befund("C: Die Linie [I] ist nativ: Sie ist die pi-Ebene von tau, der Fixpunkt im Unendlichen der nativen P_mu (z78) und "
       "der parabolischen Minkowskizeit (z83). Der Betrag von I ist durch den KMS-Zustand festgelegt (c ist zustandsnormiert), "
       "die Phase nicht: Sie ist die Modulphase (Ladung 2). Die induzierte Klammer ist exakt die Plebanski-Penrose-Klammer "
       "aus z85 (Faktor 2 kappa^2). Nativ ist I damit bis auf die Modulphase, nicht als Element", True)

# --- (B) eps auf der omega-Ebene
FORM_OM = 4/(Z[0]**2*Z[1]**2) - Z[2]/(Z[0]**3*Z[1]**2) - Z[3]/(Z[0]**2*Z[1]**3)
ok_B1 = sp.expand(br_om - 4*ab*Af**3*Bf**3*FORM_OM) == 0
keys = sorted(set(v_pi) | set(v_om))
rang_po = np.linalg.matrix_rank(np.array([[complex(v_pi.get(k, 0)), complex(v_om.get(k, 0))] for k in keys]))
ZAHLEN["klammer_omega"] = {"B_I_omega(h_a, h_b)": str(sp.factor(br_om)), "Fock": {str(k): str(v) for k, v in v_om.items()},
                           "Rang (pi, omega)": int(rang_po)}
check(f"(B) eps auf der omega-Ebene: B_I_omega(h_a, h_b) = 4 [a b] A^3 B^3 (4/(z1^2 z2^2) - z3/(z1^3 z2^2) - z4/(z1^2 z2^3)) "
      f"(exakt): gleicher E = 5-Teil, entgegengesetzter E = 6-Teil, als Fockvektor von der pi-Klammer unabhängig (Rang "
      f"{rang_po}); ihre Derivationen enthalten die K_mu, nicht die P_mu (Teil 1): I_omega ist der Infinity-Twistor des "
      f"invertierten Raums, also der Ursprungs-Twistor der nativen Minkowskizeit, nicht deren Infinity-Twistor",
      ok_B1 and rang_po == 2 and {D_von(k) for k in v_om} == {-4})

# --- (C) Uhr, Ladung und Vertexoperatoren
xs_k = [r + 0.25 for r in range(-4, 4)]
nk = len(xs_k)
Zp = sps.diags([1.0, -1.0], format="csr")
sm = sps.csr_matrix(np.array([[0.0, 1.0], [0.0, 0.0]]))
psi = []
for j in range(nk):
    mats = [Zp]*j + [sm] + [sps.identity(2, format="csr")]*(nk - j - 1)
    op = mats[0]
    for Mx in mats[1:]:
        op = sps.kron(op, Mx, format="csr")
    psi.append(op)
num_k = [np.real((p_.T @ p_).diagonal()) for p_ in psi]
sea_k = [1.0 if x < 0 else 0.0 for x in xs_k]
Qk = sps.diags(sum(n_ - s_ for n_, s_ in zip(num_k, sea_k)))
L0k = sum(x*(n_ - s_) for x, n_, s_ in zip(xs_k, num_k, sea_k))
Rk = sps.diags(np.exp(2j*np.pi*L0k))


def Jk(m):
    op = sps.csr_matrix((2**nk, 2**nk))
    for j in range(nk):
        if j + m in range(nk):
            op = op + psi[j].T @ psi[j + m]
    return op.tocsr()


ok_inn = all(abs(komm(Qk, Jk(m))).max() < 1e-12 and abs(komm(Rk, Jk(m))).max() < 1e-12 for m in (1, -1, 2, -2, 3))
# Vertexoperatoren: Verschiebungstwistor A_Q = (alpha, -q alpha)/s, alpha = Q e^{i m theta}/sqrt(m)
ok_AQ, x_ok, gram_ok = True, True, True
for Qv in ((1, 0), (0, 1), (1, 1), (1, -1)):
    alv = 0.2*np.array(Qv, complex)*np.exp(0.4j)
    v = W_auf(alv, OMV)
    Zerw = np.array([erwartung(v, ZT[k]) for k in range(4)])
    ok_AQ &= np.allclose(Zerw, np.concatenate([alv, -qt*alv])/st, atol=1e-6)
    pi_ = (Zerw[:2] + Zerw[2:])/(np.sqrt(2)*kt)
    om_ = kt*(Zerw[:2] - Zerw[2:])/np.sqrt(2)
    x_ok &= np.allclose(om_, pi_, atol=1e-6)                   # omega = i x pi mit x = -i 1
vert_tab = {}
for m in (1, 2, 3):
    k2 = KAP2_NATIV[m]
    for Q1, Q2 in (((1, 0), (0, 1)), ((1, 1), (1, -1)), ((1, 0), (1, 1)), ((1, 0), (2, 0))):
        det = Q1[0]*Q2[1] - Q1[1]*Q2[0]
        vert_tab[f"m={m} {Q1}|{Q2}"] = {"eps_pi(delta c) / e^(2 i m theta)": k2*det/m, "det": det}
AQ1 = np.concatenate([np.array([1, 0]), -qt*np.array([1, 0])])/st
AQ2 = np.concatenate([np.array([1, 1]), -qt*np.array([1, 1])])/st
I_V = np.outer(AQ1, AQ2) - np.outer(AQ2, AQ1)
DER_V = derivationen(I_V.astype(complex))
plane = np.array([AQ1, AQ2]).T
gram = plane.T @ G4 @ plane
sig_V = signatur(DER_V)
ZAHLEN["vertexoperatoren"] = {"Tabelle (nativ kappa^2 = tanh(pi m))": vert_tab, "Der(A_Q1 ^ A_Q2)": len(DER_V),
                              "Gram-Matrix der Ebene (CCR-Form)": np.round(gram, 10).tolist(), "Spurform-Signatur": sig_V,
                              "Punkt der Linie": "x = -i 1 (alle Q, theta)"}
check(f"(C) Uhr und Ladung: Q und R = e^(2 pi i L_0) vertauschen mit allen Strommoden J_m (fermionisches Modell, 8 Moden), "
      f"also mit allen Oszillatoren; ihr Stabilisator ist ganz u(2,2), und der einzige u(2,2)-invariante Bivektor ist 0: "
      f"Innere Daten liefern kein I. Vertexoperatoren verschieben Z um A_Q = (alpha, -q alpha)/s (numerisch <Z>), mit pi = "
      f"omega = alpha/sqrt2: Alle A_Q liegen auf der Linie des komplexen Punkts x = -i 1. Zwei Ladungen geben eps_pi(delta c) "
      f"= tanh(pi m) det(Q1, Q2) e^(2 i m theta)/m (nicht universell über m, Phase der Nahtlage, Mittel 0); der Bivektor "
      f"A_Q1 ^ A_Q2 ist einfach, aber seine Ebene ist definit (Gram {np.round(np.linalg.eigvalsh(gram), 4).tolist()}), "
      f"Derivationen {len(DER_V)} statt 10, Spurform {sig_V} (kompakt): kein Infinity-Twistor eines reellen Minkowskiraums",
      ok_inn and ok_AQ and x_ok and len(DER_V) == 7 and np.all(np.linalg.eigvalsh(gram) > 0) and sig_V[0] == 0
      and abs(pfaff(I_V)) < 1e-12 and not any(im_spann(M, DER_V) for M in MP))

# --- (D) zwei Nähte relativ über der kommutativen c-Algebra
g_ = sp.symbols("g1 g2")
gq_ = sp.symbols("gq1 gq2")
psi0 = sp.exp(-(g_[0]*gq_[0] + g_[1]*gq_[1])/2)
GAUSS_INV = sp.exp((g_[0]*gq_[0] + g_[1]*gq_[1])/2)


def opc(kind, i, F):
    if kind == "d":
        return sp.diff(F, gq_[i]) + g_[i]*F/2
    if kind == "e":
        return sp.diff(F, g_[i]) + gq_[i]*F/2
    if kind == "dp":
        return -sp.diff(F, g_[i]) + gq_[i]*F/2
    return -sp.diff(F, gq_[i]) + g_[i]*F/2                           # e^+


def ist_null(expr):
    return sp.expand(sp.powsimp(sp.expand(expr*GAUSS_INV))) == 0


Ft = (1 + g_[0]**2*gq_[1] + 3*g_[1]*gq_[0]**2 - 2*g_[0]*g_[1]*gq_[0])*psi0


def kc(k1, i1, k2, i2, F):
    return opc(k1, i1, opc(k2, i2, F)) - opc(k2, i2, opc(k1, i1, F))


ok_crep = (ist_null(kc("d", 0, "dp", 0, Ft) - Ft) and ist_null(kc("e", 1, "ep", 1, Ft) - Ft) and ist_null(kc("d", 0, "e", 0, Ft))
           and ist_null(kc("d", 0, "ep", 0, Ft)) and ist_null(kc("d", 0, "dp", 1, Ft)) and ist_null(opc("d", 0, psi0))
           and ist_null(opc("e", 1, psi0)))
ok_crep &= all(ist_null(opc("d", i, Ft) + opc("ep", i, Ft) - g_[i]*Ft) and ist_null(opc("dp", i, Ft) + opc("e", i, Ft) - gq_[i]*Ft)
               for i in range(2))
Fg = psi0
for _ in range(4):
    Fg = A_[0]*opc("ep", 0, Fg) + A_[1]*opc("ep", 1, Fg)
ok_grav = ist_null(Fg - (A_[0]*g_[0] + A_[1]*g_[1])**4*psi0)
tt = sp.Symbol("t", real=True)
dreh = {g_[0]: sp.exp(-sp.I*tt)*g_[0], g_[1]: sp.exp(-sp.I*tt)*g_[1], gq_[0]: sp.exp(sp.I*tt)*gq_[0],
        gq_[1]: sp.exp(sp.I*tt)*gq_[1]}
grav_a = (A_[0]*g_[0] + A_[1]*g_[1])**4*psi0
grav_b = (B_[0]*g_[0] + B_[1]*g_[1])**4*psi0
fus = grav_a*grav_b
ok_fus = (ist_null(sp.expand(grav_a.subs(dreh, simultaneous=True) - sp.exp(-4*sp.I*tt)*grav_a))
          and sp.expand(sp.powsimp(sp.expand((fus.subs(dreh, simultaneous=True) - sp.exp(-8*sp.I*tt)*fus)*GAUSS_INV**2)))
          == 0)
g1s, g2s = sp.symbols("u1 u2"), sp.symbols("w1 w2")
g1q, g2q = sp.symbols("uq1 uq2"), sp.symbols("wq1 wq2")
eps12 = g1s[0]*g2s[1] - g1s[1]*g2s[0]
diag_null = sp.expand(eps12.subs({g2s[0]: g1s[0], g2s[1]: g1s[1]})) == 0
Ssy = [sp.Matrix(s.tolist()).applyfunc(sp.nsimplify) for s in SIGMA]
Pm = lambda gg, ggq, s: sp.Rational(1, 2)*sum(ggq[i]*s[i, j]*gg[j] for i in range(2) for j in range(2))
Ptot = [Pm(g1s, g1q, s) + Pm(g2s, g2q, s) for s in Ssy]
M2 = sp.expand(Ptot[0]**2 - Ptot[1]**2 - Ptot[2]**2 - Ptot[3]**2)
ok_m2 = sp.expand(M2 - eps12*(g1q[0]*g2q[1] - g1q[1]*g2q[0])) == 0
check("(D) c-Darstellung: d = d/d gq + g/2, e = d/d g + gq/2 (und Adjungierte) erfüllen die CCR, c = d + e^+ ist "
      "Multiplikation mit g, c^+ mit gq, psi0 = e^(-|g|^2/2) ist das Vakuum, (a.e^+)^4 psi0 = (a.g)^4 psi0: Gravitonen "
      "tiefster Energie sind holomorphe Quartiken im Impulsspinor (Modulphase e^(-4 i t), D = -4). Die Fusion über der "
      "c-Algebra (punktweises Produkt) gibt (a.g)^4 (b.g)^4 psi0^2 mit e^(-8 i t), also D = -8: additiv, modulkovariant. Die "
      "Kreuzkontraktion <c(1) c(2)> = eps^jk g_j(1) g_k(2) verschwindet auf der Diagonale g(1) = g(2) identisch, und die "
      "invariante Masse zweier Nähte ist (P(1) + P(2))^2 = |<c(1) c(2)>|^2 (Fierz-Lagrange, exakt): Das Identifizieren der "
      "Impulse ist der kollineare Rand M^2 = 0, auf dem I verschwindet", ok_crep and ok_grav and ok_fus and diag_null and ok_m2)
# Zerlegung der Klammer in zwei Nähte
basis5 = [(0, 0, m, 4 - m) for m in range(5)]
basis12 = [(1, 0, m, 5 - m) for m in range(6)] + [(0, 1, m, 5 - m) for m in range(6)]
paare = [(N, M) for N in basis5 for M in basis5] + [(N, M) for N in basis5[:2] for M in basis12]
ok_zer = True
for N, M in paare:
    F = [funktion_von(c_plus(j, {N: 1})) for j in range(2)]
    Gf = [funktion_von(c_plus(j, {M: 1})) for j in range(2)]
    zwei = fock_von(F[0]*Gf[1] - F[1]*Gf[0])
    direkt = fock_von(klammer(I_PI, elementar(N), elementar(M)))
    ok_zer &= all(sp.simplify(zwei.get(k, 0) - direkt.get(k, 0)) == 0 for k in set(zwei) | set(direkt))
ok_cd = all(sp.expand(funktion_von(c_plus(j, {N: 1})) + sum(U_PI[j][a_]*sp.diff(elementar(N), Z[a_]) for a_ in range(4))) == 0
            for N in basis5 + basis12 for j in range(2))
check(f"(D) Zerlegung: Auf allen {len(paare)} Paaren (Gravitonen E = 3 und E = 4) ist B_I_pi(f, g) = m_Z(eps^jk (c_j^+ f) (x) "
      f"(c_k^+ g)) exakt: c_j^+ wirkt auf Elementarzuständen als -(u_j . grad) (omega-Ableitung), die Kreuzkontraktion zweier "
      f"getrennter Nähte mit eps (D +2) und danach das Produkt der Elementarzustände m_Z (Fusion über der holomorphen "
      f"Twistoralgebra A_T = C[a_i, b_i^+], D +2 wie bei Yang-Mills) geben die Klammer. Die Impulse der Nähte müssen bis "
      f"nach der eps-Kontraktion getrennt bleiben; die Fusion über A_c identifiziert sie vorher und löscht I", ok_zer and ok_cd)
befund("E/C: Kandidat (D) in der gestellten Form (relatives Tensorprodukt über der c-Algebra) trägt kein I: Es ist additiv "
       "(D(-4) x D(-4) -> D(-8)) und setzt die Kreuzkontraktion <c(1) c(2)> = 0 (kollinear, M^2 = 0). Die Klammer braucht "
       "stattdessen zwei getrennte Nähte, ihre eps-Kontraktion (nativ bis auf die Modulphase) und danach die Fusion über der "
       "holomorphen Twistoralgebra (wie Yang-Mills)", True)

# --- (E) Zwei-Naht-Bivektor Z_i ^ Z_j aus Zuständen zweier Nähte (elfte Analyse)
NV2 = 3
dd1, dd2, ee1, ee2 = fock4(NV2)
Z1 = [dd1.toarray(), dd2.toarray(), dag(ee1).toarray(), dag(ee2).toarray()]
OCC2 = np.array(list(itertools.product(range(NV2), repeat=4)))
DV2 = OCC2[:, 0] + OCC2[:, 1] - OCC2[:, 2] - OCC2[:, 3]
dim1 = NV2**4


def biv_erw(C):
    """B^(alpha beta) = <Psi| Z(1)^alpha Z(2)^beta - (alpha <-> beta) |Psi>, Psi = Summe C_ij |i>|j>."""
    B = np.array([[np.trace(C.conj().T @ Z1[a_] @ C @ Z1[b_].T) for b_ in range(4)] for a_ in range(4)])
    return B - B.T


C_vak = np.zeros((dim1, dim1), complex)
C_vak[0, 0] = 1
jidx = [np.ravel_multi_index((o[2], o[3], o[0], o[1]), (NV2,)*4) for o in OCC2]
C_tfd = np.zeros((dim1, dim1), complex)
for i_, o in enumerate(OCC2):
    C_tfd[i_, jidx[i_]] = qt**(o.sum()/2)
C_tfd /= np.linalg.norm(C_tfd)
mask0 = (DV2[:, None] + DV2[None, :]) == 0
mask2 = (DV2[:, None] + DV2[None, :]) == -2
C_inv = (rng.normal(size=(dim1, dim1)) + 1j*rng.normal(size=(dim1, dim1)))*mask0
C_inv /= np.linalg.norm(C_inv)
C_br = (rng.normal(size=(dim1, dim1)) + 1j*rng.normal(size=(dim1, dim1)))*(mask0 | mask2)
C_br /= np.linalg.norm(C_br)
B_vak, B_tfd, B_inv, B_br = (np.abs(biv_erw(C)).max() for C in (C_vak, C_tfd, C_inv, C_br))
Ub2 = np.diag(np.exp(1j*0.9*DV2))
inv_tfd = np.abs(Ub2 @ C_tfd @ Ub2.T - C_tfd).max()
# kohärente Zwei-Naht-Zustände: Vertexverschiebung gegen Verschiebung in der c-Algebra
beta1, beta2 = np.array([0.12, 0.05j]), np.array([-0.03, 0.1 + 0.02j])


def c_koh(beta):
    Xg = sum((beta[j]*dag(CC[j]) - np.conj(beta[j])*CC[j] for j in range(2)), sps.csr_matrix((DIMV, DIMV))).tocsc()
    return spla.expm_multiply(Xg, OMV)


Zb = [np.array([erwartung(c_koh(bb), ZT[k]) for k in range(4)]) for bb in (beta1, beta2)]
ok_ck = all(np.allclose(Zb[s], np.concatenate([bb, -bb]), atol=1e-6) for s, bb in enumerate((beta1, beta2)))
I_beta = np.outer(Zb[0], Zb[1]) - np.outer(Zb[1], Zb[0])
det_b = eps2(beta1, beta2)
ok_ib = np.allclose(I_beta, det_b*Inp, atol=1e-6) and abs(pfaff(I_beta)) < 1e-10
plane_b = np.array([np.concatenate([beta1, -beta1]), np.concatenate([beta2, -beta2])]).T
null_b = np.abs(plane_b.conj().T @ G4 @ plane_b).max() < 1e-14
qs = sp.Rational(3, 10)
ss = sp.sqrt(1 - qs**2)
I_Vs = bivektor(sp.Matrix([1, 0, -qs, 0])/ss, sp.Matrix([1, 1, -qs, -qs])/ss)    # A_Q1 ^ A_Q2, Q1 = (1,0), Q2 = (1,1)
v_v = fock_von(klammer(I_Vs, ha1, hb1))
keys_v = sorted(set(v_pi) | set(v_v))
rang_pv = np.linalg.matrix_rank(np.array([[complex(v_pi.get(k, 0)), complex(v_v.get(k, 0))] for k in keys_v]))
ZAHLEN["zwei_naht_bivektor"] = {"|B| Omega x Omega": float(B_vak), "|B| Purifikation zweier Saaten": float(B_tfd),
                                "|B| zufällig modulinvariant": float(B_inv), "|B| modulbrechend": float(B_br),
                                "det(beta1, beta2)": str(np.round(det_b, 6)), "Klammer": "det(beta1, beta2) 2 kappa^2 {.,.}_z85"}
check(f"(E) Zwei-Naht-Bivektor <Z(1) ^ Z(2)> (zwei Nähte, je 3^4): {B_vak:.1e} in Omega x Omega, {B_tfd:.1e} in der "
      f"Purifikation zweier Saaten (Zwei-Naht-TFD, modulinvariant: {inv_tfd:.1e}), {B_inv:.1e} in einem zufälligen "
      f"modulinvarianten Zustand, aber {B_br:.2f} in einem Zustand mit D_ges in {{0, -2}}: Ein Bivektor wird nur von "
      f"einem modulbrechenden Zwei-Naht-Zustand gewählt (Ladung 2). Kohärente Nähte: Vertexverschiebungen geben die "
      f"definite Ebene aus (C), und ihre Klammer B_(A_Q1 ^ A_Q2)(h_a, h_b) ist von der Plebanski-Penrose-Klammer linear "
      f"unabhängig (Rang {rang_pv}); Verschiebungen aus der c-Algebra, exp(beta.c^+ - h.c.), geben <Z> = (beta, -beta) "
      f"(pi = 0, numerisch) und <Z(1)> ^ <Z(2)> = det(beta1, beta2) I_pi exakt: einfach, null, Derivationen Poincaré, "
      f"Klammer = det(beta1, beta2) 2 kappa^2 {{.,.}}_z85", B_vak < 1e-14 and B_tfd < 1e-14 and B_inv < 1e-14
      and B_br > 1e-3 and inv_tfd < 1e-14 and ok_ck and ok_ib and null_b and rang_pv == 2
      and {D_von(k) for k in v_v} == {-4})
befund("C/E: Die elfte Analyse hat recht, dass ein kollektiver Zustand I wählen muss, und präzisiert: Kein modulinvarianter "
       "Zwei-Naht-Zustand (Vakuum, Purifikation zweier Saaten) wählt einen Bivektor. Ein kohärenter Zustand wählt einen "
       "einfachen Bivektor; Minkowski-tauglich (null, Poincaré) ist er nur für Verschiebungen aus der c-Algebra, also "
       "Impulsphasen exp(beta.c^+ - h.c.) zweier Nähte mit det(beta1, beta2) != 0; dann ist die Klammer genau die "
       "Plebanski-Penrose-Klammer mal det(beta1, beta2). Gebrochen werden die Modul-U(1), die Dilatation und die K_mu; die "
       "Amplitude det(beta1, beta2) setzt die Kopplung und ist nicht nativ", True)

# --- (F) diagonale E8-Komposition: Mischterm Summe J_1^a J_2^a
D2s = sps.diags(DV2.astype(float)).tocsr()
Id1 = sps.identity(dim1, format="csr")
Dtot = (sps.kron(D2s, Id1) + sps.kron(Id1, D2s)).tocsr()
s2 = math.sqrt(1 - qt**2)
a2 = [((dd1 + qt*dag(ee1))/s2).tocsr(), ((dd2 + qt*dag(ee2))/s2).tocsr()]   # a_i = (d_i + q e_i^+)/s auf 3^4
Xmix = sum((sps.kron(a2[i], dag(a2[i])) + sps.kron(dag(a2[i]), a2[i]) for i in range(2)),
           sps.csr_matrix((dim1**2, dim1**2))).tocsr()
Xeps = (sps.kron(a2[0], a2[1]) - sps.kron(a2[1], a2[0])).tocsr()
mix0 = float(abs(Dtot @ Xmix - Xmix @ Dtot).max())
epsl = float(abs(Dtot @ Xeps - Xeps @ Dtot + 2*Xeps).max())
c_coset = 8 + 8 - 2*248/(2 + 30)
ZAHLEN["e8_diagonal"] = {"c(E8_1 x E8_1 / E8_2)": c_coset, "[D_ges, Mischterm]": float(mix0),
                         "[D_ges, eps-Term] + 2 eps-Term": float(epsl)}
check(f"(F) diagonale E8-Komposition: T_rel = (T_1 + T_2 - Summe J_1^a J_2^a)/32 ist der Coset (E8)_1 x (E8)_1/(E8)_2 mit "
      f"c = {c_coset} (Ising). Der Mischterm paart die Nähte über die symmetrische Killingform delta^ab; auf den Oszillatoren "
      f"(Cartananteil, J_(1,1) J_(2,-1) + h.c. ~ a_1 a_2^+ + a_1^+ a_2) vertauscht er mit D_ges = D_1 + D_2 ({mix0:.1e}), "
      f"seine Nullmode ist modulneutral. Die eps-Paarung eps^ij a_1i a_2j, die I bräuchte, hat die Ladung -2 ({epsl:.1e}) und "
      f"ist nicht E8-invariant (eine invariante symmetrische Form hat keinen antisymmetrischen Teil): Der Mischterm kann die "
      f"Graduierung der Klammer nicht tragen", abs(c_coset - 0.5) < 1e-12 and mix0 < 1e-12 and epsl < 1e-12)
laufzeit()

# =================================================================== Teil 3
say("\nTeil 3: Skala, Büschel und Universalität")
# --- 3.1 Lorentzbüschel lambda I_pi + mu I_omega: flach, AdS, dS
I_ADS, I_DS = Inp - Ion, Inp + Ion
DER_ADS, DER_DS = derivationen(I_ADS), derivationen(I_DS)
inv_lorE = inv_bivektoren(LOR + [MEt])
inv_lorP = inv_bivektoren(LOR + [MP[0]])
sig_ads, sig_ds = signatur(DER_ADS), signatur(DER_DS)
br_ads, br_ds = sp.expand(br_pi - br_om), sp.expand(br_pi + br_om)
E_ads = sorted({str(E_von(k)) for k in fock_von(br_ads.subs(SUB))})
E_ds = sorted({str(E_von(k)) for k in fock_von(br_ds.subs(SUB))})
ok_bue = (len(inv_lorE) == 2 and all(im_spann(X, [I_ADS, 1j*I_ADS]) for X in inv_lorE) and len(inv_lorP) == 2
          and all(im_spann(X, [Inp, 1j*Inp]) for X in inv_lorP) and len(DER_ADS) == 10 and len(DER_DS) == 10
          and im_spann(MEt, DER_ADS) and not im_spann(MEt, DER_DS) and not im_spann(MP[0], DER_ADS)
          and sig_ads == (6, 4, 0) and sig_ds == (4, 6, 0) and abs(pfaff(I_ADS)) > 0.5 and abs(pfaff(I_DS)) > 0.5
          and abs(pfaff(Inp)) < 1e-14 and abs(pfaff(Ion)) < 1e-14
          and sp.expand(br_ads - 8*ab*Af**3*Bf**3*(Z[2]/(Z[0]**3*Z[1]**2) + Z[3]/(Z[0]**2*Z[1]**3))) == 0
          and E_ads == ["6"] and E_ds == ["5"] and verschiebung[list(ABB)[3]]["E-Verschiebung"] == [0])
ZAHLEN["buesche"] = {"Lorentz + E invariant": "C (I_pi - I_omega) (AdS)", "Lorentz + P_0 invariant": "C I_pi (flach)",
                     "Signatur Der(I_AdS)": sig_ads, "Signatur Der(I_dS)": sig_ds, "Pfaffsche": {"I_pi": 0, "I_omega": 0,
                     "I_AdS": complex(pfaff(I_ADS)).real, "I_dS": complex(pfaff(I_DS)).real},
                     "B_AdS(h_a, h_b)": str(sp.factor(br_ads)), "E von B_AdS": E_ads, "E von B_dS": E_ds}
check(f"Lorentzbüschel: Die Lorentz-invarianten Bivektoren sind lambda I_pi + mu I_omega. Zusätzlich P_0-invariant ist "
      f"genau C I_pi (einfach, flach); zusätzlich E-invariant ist genau C (I_pi - I_omega) (nicht einfach, Pfaffsche "
      f"{complex(pfaff(I_ADS)).real:.0f}; Derivationen 10 mit Spurform {sig_ads} = so(3,2) mit kompaktem E: AdS); I_pi + "
      f"I_omega gibt {sig_ds} = so(4,1) (dS). Die AdS-Klammer ist B_AdS(h_a, h_b) = 8 [a b] A^3 B^3 (z3/(z1^3 z2^2) + "
      f"z4/(z1^2 z2^3)), rein E = {E_ads} (E-additiv auf allen Probenpaaren), die dS-Klammer rein E = {E_ds}; die flache "
      f"Klammer ist ihr Mittel B_pi = (B_AdS + B_dS)/2", ok_bue)
w_, p_ = sp.symbols("w0 w1"), sp.symbols("p0 p1")
zsub = {Z[0]: (KAP*p_[0] + w_[0]/KAP)/sp.sqrt(2), Z[1]: (KAP*p_[1] + w_[1]/KAP)/sp.sqrt(2),
        Z[2]: (KAP*p_[0] - w_[0]/KAP)/sp.sqrt(2), Z[3]: (KAP*p_[1] - w_[1]/KAP)/sp.sqrt(2)}
ftest = Z[2]**3*Z[3]/(Z[0]**2*Z[1]) + Z[3]**2/(Z[0]*Z[1]**2)
ok_kk = all(sp.simplify(sp.diff(ftest.subs(zsub), w_[j]) - (sum(U_PI[j][a_]*sp.diff(ftest, Z[a_]) for a_ in range(4))/(
    sp.sqrt(2)*KAP)).subs(zsub)) == 0 and sp.simplify(sp.diff(ftest.subs(zsub), p_[j]) - (KAP*sum(
        V_OM[j][a_]*sp.diff(ftest, Z[a_]) for a_ in range(4))/sp.sqrt(2)).subs(zsub)) == 0 for j in range(2))
Lam_hat = sp.simplify((2/KAP**2)/(2*KAP**2))
ok_lam = ok_kk and sp.simplify(Lam_hat*(KAP**2)**2 - 1) == 0
ZAHLEN["ads_skala"] = {"I_AdS in CCR-Ableitungen": "2 kappa^2 eps d_omega d_omega - (2/kappa^2) eps d_pi d_pi",
                       "Lambda-hat (Konvention I = eps_omega + Lambda eps_pi)": str(Lam_hat),
                       "Im x des Vakuums (z85)": "kappa^2", "|Lambda-hat| (Im x)^2": 1}
check("Skala des radialen Kanals: d/d omega_j = (u_j . grad)/(sqrt2 kappa), d/d pi_j = kappa (v_j . grad)/sqrt2 (symbolisch, "
      "Koordinaten aus z85), also I_AdS = 2 kappa^2 eps d_omega d_omega - (2/kappa^2) eps d_pi d_pi: Lambda-hat = "
      "-kappa^(-4) in CCR-Koordinaten, und mit dem Vakuum bei Im x = kappa^2 (z85) gilt |Lambda-hat| (Im x)^2 = 1 exakt. "
      "Eine E-graduierte (radiale) Komposition sieht nur den E-erhaltenden Teil der flachen Klammer, B_AdS/2: Der radiale "
      "Kanal gibt höchstens AdS mit einem Krümmungsradius von der Größe des KMS-Vakuums, nicht die flache Klammer", ok_lam)
befund("C/L: Der radiale Kanal (Abschnitt 11, Rechnung 1) ist damit entschieden: Mit E (Moden = Stufen) und Lorentz ist der "
       "Infinity-Twistor eindeutig I_pi - I_omega (AdS, |Lambda| R_Vakuum^2 = 1). Die flache Plebanski-Penrose-Klammer "
       "verlangt die parabolische Zeit P_0 statt E; ihr E = 5-Anteil ist auf der radialen Leiter unsichtbar", True)

# --- 3.2 Normierungen und Universalität


def sektor_aufgeloest(beta, alpha, R):
    xs = np.array([r + alpha for r in range(-R, R)])
    nm = len(xs)
    idx = np.arange(2**nm)
    bits = ((idx[:, None] >> np.arange(nm)[::-1]) & 1).astype(float)
    lc = -np.logaddexp(0, beta*xs)
    l1c = beta*xs - np.logaddexp(0, beta*xs)
    logw = bits @ lc + (1 - bits) @ l1c
    w = np.exp(logw - logw.max())
    Q = np.rint((bits - (xs < 0)).sum(1)).astype(int)
    f1 = (bits[:, 1:]*(1 - bits[:, :-1])).sum(1)
    out = {}
    for Qv in range(-2, 3):
        sel = Q == Qv
        out[Qv] = float((w[sel]*f1[sel]).sum()/w[sel].sum())
    return out


sek = sektor_aufgeloest(4*np.pi, 0.25, 7)
soll1 = 1/np.expm1(4*np.pi)
dev_sek = max(abs(v/soll1 - 1) for v in sek.values())
Asym = sp.symbols("A1:5")
fsh = lambda f: f.subs({Z[k]: Z[k] + Asym[k] for k in range(4)}, simultaneous=True)
f1_, f2_ = elementar((0, 0, 3, 1)), elementar((1, 0, 2, 3))
ok_versch = sp.cancel(sp.together(klammer(I_PI, fsh(f1_), fsh(f2_)) - fsh(klammer(I_PI, f1_, f2_)))) == 0
hel_ok = True
for D, L in PROBEN.items():
    for N in L:
        for K in fock_von(klammer(I_PI, ha1, elementar(N))):
            hel_ok &= D_von(K) == D
h2s = [elementar((0, 0, 2, 2)), elementar((1, 0, 3, 2)), (Z[2] + 2*Z[3])**4/(Z[0]*Z[1])]
jac = sp.expand(klammer(I_PI, h2s[0], klammer(I_PI, h2s[1], h2s[2])) + klammer(I_PI, h2s[1], klammer(I_PI, h2s[2], h2s[0]))
                + klammer(I_PI, h2s[2], klammer(I_PI, h2s[0], h2s[1])))
lb2 = sp.expand(klammer(I_PI, ha1, f1_*f2_) - klammer(I_PI, ha1, f1_)*f2_ - f1_*klammer(I_PI, ha1, f2_))
ZAHLEN["universalitaet"] = {"<J_-1 J_1>_Q bei 4 pi (Q = -2..2)": sek, "max. rel. Abweichung": dev_sek,
                            "Normierung Zustand/Gitter": {f"m={m}": KAP2_NATIV[m]/m for m in (1, 2, 3)}}
check(f"Universalität bei gegebenem I: (i) Ladungssektoren: <J_-1 J_1>_Q = 1/(e^(4 pi) - 1) in jedem Sektor Q = -2..2 (rel. "
      f"{dev_sek:.1e}), also derselbe Oszillatorzustand und dieselbe Zustandsnormierung von I; (ii) Vertexoperatoren: "
      f"B_I vertauscht mit jeder Twistorverschiebung Z -> Z + A (symbolisch); (iii) Helizitäten: V_h = B_I(h, .) mit Grad 2 "
      f"erhält jeden Sektor D = -4..4, erfüllt Leibniz mit dem Produkt und Jacobi (exakt): minimale Kopplung derselben "
      f"Klammer an alle Felder. Die Gitternormierung (eps(Q1, Q2) = det) weicht davon um tanh(pi m)/m ab "
      f"({ {m: round(KAP2_NATIV[m]/m, 6) for m in (1, 2, 3)} }): nicht über die Moden universell",
      dev_sek < 1e-8 and ok_versch and hel_ok and jac == 0 and lb2 == 0)


# --- 3.3 kritische Folge: Minkowskizeit gegen E und Modulgenerator
def P0_sektor(D, Emax):
    basis = []
    n = max(D, 0)
    while n - D/2 + 1 <= Emax + 1e-9:
        mN = n - D
        basis += [(n1, n - n1, m1, mN - m1) for n1 in range(n + 1) for m1 in range(mN + 1)]
        n += 1
    idx = {s: i for i, s in enumerate(basis)}
    r, c, v = [], [], []
    for s, i in idx.items():
        r.append(i)
        c.append(i)
        v.append(sum(s)/2 + 1)
        for k in range(2):
            t = list(s)
            t[k] += 1
            t[k + 2] += 1
            j = idx.get(tuple(t))
            if j is not None:
                amp = 0.5*math.sqrt((s[k] + 1)*(s[k + 2] + 1))
                r += [j, i]
                c += [i, j]
                v += [amp, amp]
    return sps.csr_matrix((v, (r, c)), shape=(len(basis),)*2)


krit = {}
for Emax in (10, 20, 40):
    P0m = P0_sektor(-4, Emax)
    lam = float(spla.eigsh(P0m.tocsc(), k=1, sigma=-0.01, which="LM", return_eigenvectors=False)[0])
    krit[Emax] = {"dim": P0m.shape[0], "lambda_min": lam, "lambda_min * E_max": lam*Emax}
lams = [krit[E]["lambda_min"] for E in (10, 20, 40)]
prod_ = [krit[E]["lambda_min * E_max"] for E in (10, 20, 40)]
grenze = float(mp.besseljzero(5, 1))**2/8                         # tiefster K-Typ E_0 = 3: Leiter D_3, Bessel J_(2 E_0 - 1)
ZAHLEN["kritische_folge"] = {"P_0 auf D = -4": {str(k): v for k, v in krit.items()}, "j_(5,1)^2/8": grenze,
                             "E auf D = -4": "Spektrum ab 3", "Modulgenerator auf D = +-4": "|K| >= 16 pi m"}
check(f"Kritische Folge (GRAV.PARENT.GAP.NO_GO.01): Die Minkowskizeit P_0 = (1/2) c^+ c ist auf dem Gravitonsektor D = -4 "
      f"lückenlos: kleinste Eigenwerte {[round(x, 5) for x in lams]} bei E_max = 10, 20, 40 (Dimensionen "
      f"{[krit[E]['dim'] for E in (10, 20, 40)]}), lambda_min E_max = {[round(x, 4) for x in prod_]} steigt gegen "
      f"j_(5,1)^2/8 = {grenze:.4f} (tiefster K-Typ E_0 = 3, Bessel J_(2 E_0 - 1); für den Skalar gibt dieselbe Regel "
      f"j_(1,1)^2/8 = 1,8353 aus z83), also lambda_min ~ 9,62/E_max -> 0. Die konforme Energie E ist dort ab 3 gegappt, der "
      f"Modulgenerator ab 16 pi. Die kritische Folge gehört zu genau der Zeit, deren Infinity-Twistor I_pi ist; die "
      f"AdS-Wahl des radialen Kanals (E) hat sie nicht",
      lams[0] > lams[1] > lams[2] > 0 and prod_[0] < prod_[1] < prod_[2] < grenze and abs(prod_[2]/grenze - 1) < 0.03)
befund("C/O: Mit I gibt es Universalität (eine Klammer für alle Ladungssektoren, Vertexoperatoren und Helizitäten), eine "
       "lückenlose Zeit P_0 auf dem Gravitonsektor und die Eichalgebra der Grad-2-Hamiltonfunktionen (Jacobi). Eine "
       "Gravitationskonstante folgt nicht: Der Koeffizient g der Klammer (dbar h + g {h, h}) ist frei; I legt nur Poincaré und "
       "die Längeneinheit des KMS-Vakuums fest. In Nahteinheiten wäre kappa_G = g mal Vakuumlänge mit unbestimmtem g; im "
       "Zwei-Naht-Kondensat (E) wäre g = det(beta1, beta2), ebenfalls frei", True)
laufzeit()

# =================================================================== Teil 4
say("\nTeil 4: Bilanz und Urteil")
STATUS = {
    "Auswahlregel (zehnte Analyse)": ("E, bestätigt, nicht spezifisch", "native Kompositionen D-additiv; Klammer +4, "
                                      "Yang-Mills-Produkt +2; schließt jeden Vertex mit Summe h != 0 aus"),
    "Energieschranke (zehnte Analyse)": ("E, bestätigt, nicht spezifisch", "Hom(D(-4) x D(-4), D(-4)) = 0; ebenso D(-2), D(0)"),
    "Zusatzstruktur I": ("E", "Derivationen der Klammer = Poincaré (10), Eulergewicht -2, kein u(2,2)-invarianter Bivektor"),
    "(A) eps auf der pi-Ebene": ("C", "Linie nativ (tau, P_mu, parabolische Zeit), Betrag aus dem KMS-Zustand, Phase = "
                                 "Modulphase; Klammer = 2 kappa^2 z85 exakt"),
    "(B) eps auf der omega-Ebene": ("E, negativ als flaches I", "Ursprungs-Twistor (K_mu), Klammer verschieden"),
    "(C) Uhr, Ladung, Vertexoperatoren": ("E, negativ", "innen trivial auf T; Vertex-Bivektoren definit (Punkt x = -i 1), "
                                          "kompakter Stabilisator, Normierung tanh(pi m)/m"),
    "(D) relativ über der c-Algebra": ("E, negativ", "D-additiv, Kreuzkontraktion auf der Diagonale 0 (kollinear, M^2 = 0)"),
    "Zerlegung der Klammer": ("E", "B_pi = m_Z o (eps c^+ (x) c^+): zwei getrennte Nähte, dann holomorphe Fusion"),
    "(E) Zwei-Naht-Bivektor": ("E/C", "modulinvariante Zustände wählen keinen; Impulsphasen-Kondensat wählt I_pi mit "
                               "g = det(beta1, beta2)"),
    "(F) diagonale E8": ("E, negativ", "Mischterm modulneutral (Killingform symmetrisch), Coset c = 1/2"),
    "radialer Kanal": ("C/L", "E + Lorentz erzwingen I_AdS = I_pi - I_omega, |Lambda| R_Vakuum^2 = 1; flach braucht P_0"),
    "Universalität bei gegebenem I": ("E/C", "Ladungssektoren, Vertexoperatoren, Helizitäten"),
    "kritische Folge": ("E/C", "P_0 lückenlos auf D = -4 (lambda_min E_max -> j_(5,1)^2/8 = 9,62), E und Modulgenerator "
                               "gegappt"),
    "G, Kopplung g": ("O", "nicht nativ; I gibt nur Poincaré und die Vakuumlänge"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
urteil = ("F1 ist bedingt beantwortet. Die zehnte Analyse stimmt in beiden Rechnungen: Native Kompositionen (Heisenberg mit "
          "hermiteschem Kozyklus, Zwei-Naht-Tensorprodukt, Fusion über der c-Algebra) sind modulkovariant und D-additiv, die "
          "Plebanski-Penrose-Klammer hat die D-Verschiebung +4, und Hom_su(2,2)(D(-4) x D(-4), D(-4)) = 0, weil die tiefsten "
          "Gewichte von H x H bei E = 6 liegen und H nur bei E = 3 eines hat. Beide Aussagen sind aber nicht "
          "gravitationsspezifisch: Sie schließen ebenso die Yang-Mills-Komposition D(-2) x D(-2) -> D(-2) (+2) und jede "
          "konform äquivariante Fusion positiver Einteilchenräume aus. Spezifisch ist allein I: Die Klammer hat als Derivationen "
          "genau Poincaré (10 statt 16), Eulergewicht -2, und kein Bivektor ist u(2,2)-invariant. Ein nativer Infinity-Twistor "
          "existiert bedingt: Die Linie [I] ist die pi-Ebene von tau, der Fixpunkt im Unendlichen der nativen Translationen "
          "P_mu = c^+ sigma c/2 und der parabolischen Minkowskizeit ((1/2) c^+ c = E + (P^0 + K^0)/2 exakt); der KMS-Zustand "
          "normiert den Betrag (<c^+ c> = 1, I_Zustand = 2 kappa^2 I_CCR), die Phase ist die Modulphase mit Ladung 2 und wird "
          "von keiner nativen Größe festgelegt. Mit diesem I ist die Einstein-Klammer exakt reproduziert, B_I(h_a, h_b) = "
          "2 kappa^2 {h_a, h_b}_z85 = 4 [a b] A^3 B^3 (...), aber keine native Komposition erzeugt sie: Die Klammer zerfällt "
          "exakt in die eps-Kreuzkontraktion zweier getrennter Nähte (D +2) und die Fusion über der holomorphen "
          "Twistoralgebra (D +2, wie Yang-Mills). Das relative Tensorprodukt über der c-Algebra ist dagegen additiv und "
          "löscht I, denn die invariante Masse zweier Nähte ist |<c(1) c(2)>|^2 und das Identifizieren der Impulse ist der "
          "kollineare Rand. Die omega-Ebene gibt den Ursprungs-Twistor, Uhr und Ladung wirken auf T trivial, "
          "Vertexoperatoren geben definite Bivektoren (Punkt x = -i 1, kompakter Stabilisator), und der Mischterm der "
          "diagonalen E8 ist modulneutral (symmetrische Killingform, Coset c = 1/2). Zur elften Analyse: Kein "
          "modulinvarianter Zwei-Naht-Zustand (Vakuum, Purifikation zweier Saaten) wählt einen Bivektor; kohärente "
          "Vertex-Zustände wählen einen einfachen, aber definiten Bivektor mit einer von Plebanski-Penrose unabhängigen "
          "Klammer; erst ein Impulsphasen-Kondensat exp(beta.c^+ - h.c.) zweier Nähte wählt genau det(beta1, beta2) I_pi und "
          "damit die Plebanski-Penrose-Klammer mit Kopplung det(beta1, beta2), also Symmetriebrechung der Modul-U(1) mit "
          "freier Amplitude. Neu ist die Entscheidung des radialen Kanals: Mit E (Moden = Stufen) und Lorentz ist der "
          "Infinity-Twistor eindeutig I_pi - I_omega (so(3,2), AdS) mit |Lambda| R_Vakuum^2 = 1; die flache Klammer braucht "
          "die parabolische Zeit P_0, deren Spektrum auf dem Gravitonsektor lückenlos ist (lambda_min E_max -> j_(5,1)^2/8, "
          "kritische Folge), während E ab 3 und der Modulgenerator ab 16 pi gegappt sind. Universalität gilt bei gegebenem "
          "I für Ladungssektoren, "
          "Vertexoperatoren und Helizitäten. Eine Gravitationskonstante folgt nicht; offen bleiben die Kopplung g (eine "
          "Kompositionsregel mit festem Koeffizienten), die Festlegung der Modulphase und die Ward-Identitäten einer "
          "wechselwirkenden Theorie. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z89_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
