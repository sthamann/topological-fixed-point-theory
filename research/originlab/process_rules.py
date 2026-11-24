"""Process rule search: which structural rule singles out the TFPT process?

Same engine as the genesis lane (`genesis_rules.search`): each universe is a finite class of candidate process
pieces, every clause is computed from the candidate, every subset of the clause pool is a rule, and the smallest
solutions decide EMPTY, UNIQUE or AMBIGUOUS. The universes are the open slots of the H6 exploration
(experiments/exploration/2026-11-23-phasenraum-f5-holocode, z39-z42):

  matrix_kernels  positive kernels on 2x2 matrices: det(I - W^+ Z)^-nu, optionally plus eps det Z conj det W.
  spin2_kernels   inverse spin-2 propagators H(x): Dirac spectral action with several cutoffs, the ledger's entire
                  ansatz, local truncations, Einstein, and Schur complements of positive two-channel pencils.
  seam_layers     stacks of n Chern layers, real (Majorana) or complex, with the second sheet half equally or
                  reversely oriented.
  flavour_pairs   pairs of residual symmetries (abelian subgroups of S4 with three distinct characters on the
                  family triplet); their eigenbases give a mixing pattern.
  transfer_processes  Hardy transfer functions a(z) = z^-d (1 - q z^k) and their compressions S^{+d}(1 - q S^k) (z44).
  readouts        source maps of the primitive quarter transfer: which pole, connected or not, which residues (z44).

The spin-2 slot also holds the Marchenko-Pastur memory of the quarter transfer for several q. `run(transfer=False)`
reproduces the lane as of H6h23 (z43): no transfer universes, no transfer memory, nine composed principles.

Only one clause carries a number: `dirac_low_energy` asks for the first correction a = 3/10 of the Dirac heat kernel
(a structural property of the Dirac operator with the e^-u cutoff, z41). Record structures are named afterwards in
`RECORD`. A unique minimal solution is a statement about the universe and its clause pool; it derives no physical
map and promotes nothing.
"""
from fractions import Fraction as Fr
from itertools import combinations, permutations
from math import gamma
import argparse
import json

import numpy as np
import sympy
from numpy.polynomial.legendre import leggauss
from scipy.special import roots_genlaguerre

from genesis_rules import search

SCOPE = ("exhaustive finite enumeration of candidate process pieces with computed clauses; a unique minimal solution "
         "names the shortest structural rule inside the universe, not a derivation of physics")


# ---------------------------------------------------------------- A: matrix kernels
def two_row_partitions(top):
    return [(a, b) for n in range(top+1) for b in range(n//2+1) for a in [n-b] if a >= b]


def schur_coefficient(nu, lam):
    """Coefficient of s_lambda(X) in det(I - X)^-nu: prod over cells (nu + j - i)/hook (2 x 2, two rows)."""
    a, b = lam
    rows = [a, b]
    cols = [sum(1 for r in rows if r > j) for j in range(a)]
    out = Fr(1)
    for i, r in enumerate(rows):
        for j in range(r):
            hook = (r - j - 1) + (cols[j] - i - 1) + 1
            out *= (nu + j - i) / Fr(hook)
    return out


def cayley_vanishes(nu, eps):
    z = sympy.symbols("z11 z12 z21 z22")
    Z = sympy.Matrix(2, 2, z)
    W = sympy.Matrix([[sympy.Rational(1, 3), sympy.Rational(1, 5)], [sympy.Rational(-1, 7), sympy.Rational(1, 4)]])
    K = (sympy.eye(2) - W.T * Z).det() ** (-sympy.Rational(nu.numerator, nu.denominator)) + eps * Z.det() * W.det()
    box = sympy.diff(K, z[0], z[3]) - sympy.diff(K, z[1], z[2])
    val = box.subs({z[0]: sympy.Rational(1, 9), z[1]: sympy.Rational(1, 11), z[2]: sympy.Rational(-1, 13), z[3]: sympy.Rational(1, 6)})
    return abs(complex(sympy.N(val, 30))) < 1e-20


def moebius_covariant(nu, eps, seed=5):
    """K(gZ, gW) = det(CZ+D)^nu conj(det(CW+D))^nu K(Z, W) for g in U(2,2) near the identity (principal branch)."""
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
    d = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
    b = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
    X = 0.15*np.block([[a - a.conj().T, b], [b.conj().T, d - d.conj().T]])
    w, V = np.linalg.eig(X)
    g = V @ np.diag(np.exp(w)) @ np.linalg.inv(V)
    A, B, C, D = g[:2, :2], g[:2, 2:], g[2:, :2], g[2:, 2:]
    act = lambda M: (A @ M + B) @ np.linalg.inv(C @ M + D)
    nuf = float(nu)
    K = lambda Z, W: np.linalg.det(np.eye(2) - W.conj().T @ Z) ** (-nuf) + eps*np.linalg.det(Z)*np.conj(np.linalg.det(W))
    worst = 0.0
    for _ in range(4):
        Z = 0.25*(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))
        W = 0.25*(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))
        lhs = K(act(Z), act(W))
        rhs = np.linalg.det(C @ Z + D) ** nuf * np.conj(np.linalg.det(C @ W + D)) ** nuf * K(Z, W)
        worst = max(worst, abs(lhs - rhs) / abs(rhs))
    return worst < 1e-9


def matrix_kernels(top=8):
    out = []
    for nu in (Fr(0), Fr(1, 2), Fr(1), Fr(3, 2), Fr(2), Fr(3)):
        for eps in (0, 1):
            coeff = {lam: schur_coefficient(nu, lam) + (eps if lam == (1, 1) else 0) for lam in two_row_partitions(top)}
            clauses = set()
            if all(c >= 0 for c in coeff.values()):
                clauses.add("positive")
            if all(coeff[(k, 0)] == 1 for k in range(top+1)):
                clauses.add("hardy_slice")
            if cayley_vanishes(nu, eps):
                clauses.add("wave_equation")
            if moebius_covariant(nu, eps):
                clauses.add("conformally_covariant")
            out.append({"nu": nu, "eps": eps, "clauses": clauses,
                        "negative_channels": sorted(str(l) for l, c in coeff.items() if c < 0)[:4]})
    return out


# ---------------------------------------------------------------- B: spin-2 kernels
_XG, _WG = leggauss(200)


def _bv_F(x):
    """F(x) = Int_0^1 exp(-x xi(1-xi)) d xi."""
    xi = (_XG + 1)/2
    return np.sum(_WG/2*np.exp(-np.multiply.outer(np.asarray(x, float), xi*(1-xi))), axis=-1)


def h_exp(x):
    F = _bv_F(x)
    return 3*x*F + 12*(F - 1)                    # x a(x) of z41/z42, H ~ x - 3x^2/10


def h_laplace(k):
    """Cutoff (1+u)^-k = Int s^{k-1} e^-s e^{-su}/Gamma(k): H_f(x) = Int mu(s) H_exp(s x)/s^2 ds, normalised to H ~ x."""
    s, w = roots_genlaguerre(80, k - 3)
    w = w/gamma(k)
    norm = float(np.sum(w*s))                  # Int mu/s ds  (H_exp(y) ~ y)
    return lambda x: np.array([np.sum(w*h_exp(s*xx)) for xx in np.atleast_1d(x)])/norm


def h_step(x):
    """Cutoff 1_[0,1]: F_1 = (1-v)_+, F_2 = (1-v)_+^2/2, integrated exactly over the region x xi(1-xi) <= 1."""
    out = []
    for xx in np.atleast_1d(x):
        pieces = [(0.0, 1.0)] if xx <= 4 else [(0.0, (1 - np.sqrt(1 - 4/xx))/2), ((1 + np.sqrt(1 - 4/xx))/2, 1.0)]
        i1 = i2 = 0.0
        for lo, hi in pieces:
            xi = lo + (hi - lo)*(_XG + 1)/2
            v = 1 - xx*xi*(1 - xi)
            i1 += np.sum((hi - lo)/2*_WG*v)
            i2 += np.sum((hi - lo)/2*_WG*v*v/2)
        out.append(3*xx*i1 + 12*(i2 - 0.5))
    return np.array(out)


Q_DIRAC = float(np.sqrt((8 - np.sqrt(55))/3))      # q^2/(1-q^2)^2 = 3/10


def transfer_kernel(qv):
    """Boundary response of A_q = S^+(1 - q S^4): C_q(x) = [sqrt((x+a)(x+b)) - x + 1 - q^2]/(2x), a, b = (1 -+ q)^2 (z44)."""
    a_, b_ = (1 - qv)**2, (1 + qv)**2
    return lambda x: (np.sqrt((np.asarray(x, float) + a_)*(np.asarray(x, float) + b_)) - np.asarray(x, float) + 1 - qv**2)/(2*np.asarray(x, float))


def spin2_kernels(transfer=True):
    rat = lambda num_den: (lambda x: num_den(np.asarray(x, float)))
    cands = [
        {"name": "Einstein H = x", "H": rat(lambda x: x), "poles": 1, "schur": None},
        {"name": "ledger ansatz H = x e^x", "H": rat(lambda x: x*np.exp(np.minimum(x, 700))), "poles": 1, "schur": None},
        {"name": "Stelle local H = x(1+x)", "H": rat(lambda x: x*(1 + x)), "poles": 2, "schur": None},
        {"name": "Dirac local H = x(1-3x/10)", "H": rat(lambda x: x*(1 - 0.3*x)), "poles": 2, "schur": None},
        {"name": "spectral action e^-u", "H": h_exp, "poles": 99, "schur": None},
        {"name": "spectral action step", "H": h_step, "poles": 99, "schur": None},
        {"name": "spectral action (1+u)^-3", "H": h_laplace(3), "poles": 99, "schur": None},
    ]
    for w in (Fr(1, 4), Fr(1, 2), Fr(3, 4)):
        for m in (Fr(5, 6), Fr(5, 3), Fr(5, 2), Fr(10, 3)):
            cands.append({"name": f"Schur w={w} m={m}", "H": (lambda ww, mm: rat(lambda x: x - ww*x*x/(x + mm)))(float(w), float(m)),
                          "poles": 2, "schur": (w, m)})
    cands.append({"name": "Schur two channels w=1/4,1/4 m=5/3,10/3",
                  "H": rat(lambda x: x - 0.25*x*x/(x + 5/3) - 0.25*x*x/(x + 10/3)), "poles": 3, "schur": "two"})
    for label, qv in (("1/4", 0.25), ("1/2", 0.5), ("3/4", 0.75), ("q_dirac", Q_DIRAC), ("e^-pi", float(np.exp(-np.pi))))[:5*transfer]:
        Cq = transfer_kernel(qv)
        cands.append({"name": f"transfer memory q={label}", "H": (lambda C, Z: lambda x: Z/C(x))(Cq, 1 - qv**2), "poles": 50,
                      "schur": None, "transfer": True})
    xs = np.geomspace(1e-3, 1e4, 900)
    out = []
    for c in cands:
        H = np.asarray(c["H"](xs), float)
        clauses = set()
        if np.all(H > 0):
            clauses.add("no_euclidean_zero")
            R = 1/H - 1/xs
            dR = np.diff(R)/np.diff(xs)
            if np.all(R > -1e-12) and np.all(dR < 1e-12) and np.all(np.diff(dR)/np.diff(xs[:-1]) > -1e-9):
                clauses.add("kallen_lehmann")
        h1, h2 = (float(v) for v in np.asarray(c["H"](np.array([1e-3, 2e-3])), float))
        # H/x = 1 - a x + b x^2: two points give a
        r1, r2 = h1/1e-3, h2/2e-3
        b = ((r2 - 1)/2e-3 - (r1 - 1)/1e-3)/(2e-3 - 1e-3)
        a = -((r1 - 1)/1e-3 - b*1e-3)
        if abs(a - 0.3) < 1e-4:
            clauses.add("dirac_low_energy")
        if isinstance(c["schur"], tuple) and c["schur"][0]/(1 - c["schur"][0]) == 1:
            clauses.add("mirror_weight")
        if c["poles"] <= 3 and "e^x" not in c["name"]:
            clauses.add("rational")
        if c.get("transfer"):
            clauses.add("transfer_memory")
        zero = None
        if np.any(H <= 0):
            k = int(np.argmax(H <= 0))
            zero = float(xs[k])
        out.append({"name": c["name"], "poles": c["poles"], "clauses": clauses, "a": round(a, 6), "first_zero_near": zero})
    return out


# ---------------------------------------------------------------- B2: transfer functions and readouts
Q_VALUES = (("1/4", 0.25), ("1/2", 0.5), ("3/4", 0.75), ("e^-pi", float(np.exp(-np.pi))), ("e^-2pi", float(np.exp(-2*np.pi))))


def _compression(d, k, qv, n=400):
    S = np.diag(np.ones(n - 1), -1)
    return np.linalg.matrix_power(S.T, d) @ (np.eye(n) - qv*np.linalg.matrix_power(S, k))


def transfer_processes():
    """Transfer functions a(z) = z^-d (1 - q z^k) on the Hardy space; compression A = S^{+d}(1 - q S^k).

    The square truncation is exact for the kernel: 1 - q S^k is unit lower triangular, so ker A = (1 - q S^k)^-1 ker S^{+d}
    has dimension d with no truncation artefact. The cokernel vanishes on the Hardy space (S^+ onto, 1 - q S^k invertible
    for q < 1), so the index is the kernel dimension.
    """
    out = []
    zs = np.exp(1j*np.array([0.3, 1.1, 2.0]))*0.9
    for d in (1, 2):
        for k in (2, 3, 4, 8):
            for qi, (label, qv) in enumerate(Q_VALUES):
                a = lambda z: z**(-d)*(1 - qv*z**k)
                ratios = [a(1j*zz)/a(zz) for zz in zs]
                clauses = set()
                if max(abs(r - ratios[0]) for r in ratios) < 1e-12:
                    clauses.add("quarter_covariant")
                A = _compression(d, k, qv, n=120)
                _, sv, vh = np.linalg.svd(A)
                kernel_dim = int(np.sum(sv < 1e-10))
                if kernel_dim == 1:
                    clauses.add("index_one")
                grades = {n % 4 for vec in vh[len(sv) - kernel_dim:] for n in np.nonzero(np.abs(vec) > 1e-8)[0]}
                if kernel_dim and len(grades) == 1:
                    clauses.add("kernel_one_grade")
                if label in ("1/4", "1/2", "3/4"):
                    clauses.add("rational_pole_weight")
                if label.startswith("e^"):
                    clauses.add("tau_i_nome")
                out.append({"d": d, "k": k, "q": label, "qi": qi, "clauses": clauses, "kernel_dim": kernel_dim,
                            "kernel_grades": sorted(int(g) for g in grades)})
    return out


def readouts():
    """Readouts of the primitive quarter transfer (d = 1, k = 4): which source map gives which pole."""
    out = []
    for label, qv in (("1/2", 0.5), ("e^-pi", float(np.exp(-np.pi)))):
        Z = 1 - qv**2
        C = transfer_kernel(qv)
        kinds = {
            "boundary source |0>": {"pole": Z, "multi": False, "connected": False, "relativistic": False, "residues": [Z]},
            "connected vacuum response": {"pole": 0.0, "multi": False, "connected": True, "relativistic": False, "residues": [0.0]},
            "4D field, one source": {"pole": Z, "multi": False, "connected": False, "relativistic": True, "residues": [Z]},
            "4D field, sources |4j>": {"pole": Z, "multi": True, "connected": False, "relativistic": True,
                                       "residues": [Z*qv**(2*j) for j in range(4)]},
            "4D field connected, sources |4j>": {"pole": 0.0, "multi": True, "connected": True, "relativistic": True,
                                                 "residues": [0.0]*4},
        }
        for kind, r in kinds.items():
            clauses = {"kl_positive"}                       # every readout is a resolvent of the positive L_q (z44)
            if r["pole"] > 1e-12:
                clauses.add("massless_pole")
            if r["connected"]:
                clauses.add("connected")
            if r["relativistic"]:
                clauses.add("relativistic")
            if r["multi"]:
                clauses.add("several_sources")
                if max(r["residues"]) - min(r["residues"]) < 1e-12 and r["residues"][0] > 0:
                    clauses.add("universal_residue")
            out.append({"q": label, "kind": kind, "clauses": clauses, "size": len(r["residues"])})
    return out


# ---------------------------------------------------------------- C: seam layers
def seam_layers(top=32):
    out = []
    for n in range(1, top+1):
        for kind in ("real", "complex"):
            for orient in ("same", "reversed"):
                c = Fr(n, 2) if kind == "real" else Fr(n)
                majorana = n if kind == "real" else 2*n
                clauses = set()
                if c == 8:
                    clauses.add("c_eight")
                if majorana % 16 == 0:
                    clauses.add("anyon_free")               # Kitaev's sixteenfold way
                if kind == "real" and n >= 10:
                    clauses.add("spin10_level_one")         # ten real modes in the vector: level 1; complex modes: even level
                if orient == "reversed":
                    clauses.add("sheet_neutral")            # c - c = 0 across the two sheet halves
                if kind == "real":
                    clauses.add("u1_free")
                out.append({"n": n, "kind": kind, "orient": orient, "c": str(c), "clauses": clauses})
    return out


# ---------------------------------------------------------------- D: flavour pairs
def _standard_rep():
    basis = np.array([[1, -1, 0, 0], [1, 1, -2, 0], [1, 1, 1, -3]], float)
    basis = (basis.T/np.linalg.norm(basis, axis=1)).T
    def rep(p):
        Pm = np.zeros((4, 4))
        for i, j in enumerate(p):
            Pm[j, i] = 1
        return basis @ Pm @ basis.T
    return rep


def _subgroups():
    els = list(permutations(range(4)))
    mul = lambda p, q: tuple(p[q[i]] for i in range(4))
    found = set()
    for gens in combinations(els, 2):
        S = {tuple(range(4))}
        frontier = list(S)
        while frontier:
            new = []
            for x in frontier:
                for g in gens:
                    y = mul(g, x)
                    if y not in S:
                        S.add(y)
                        new.append(y)
            frontier = new
        found.add(frozenset(S))
    for g in els:
        S, x = {tuple(range(4))}, g
        while x not in S:
            S.add(x)
            x = mul(g, x)
        found.add(frozenset(S))
    return found


def _kind(G):
    order = len(G)
    mul = lambda p, q: tuple(p[q[i]] for i in range(4))
    abelian = all(mul(a, b) == mul(b, a) for a in G for b in G)
    def elorder(g):
        k, x = 1, g
        while x != tuple(range(4)):
            x, k = mul(g, x), k + 1
        return k
    if not abelian:
        return None
    if order == 3:
        return "Z3"
    if order == 4:
        if any(elorder(g) == 4 for g in G):
            return "Z4"
        fixed = all(sum(1 for i in range(4) if g[i] == i) in (0, 4) for g in G)
        return "V4" if fixed else "V4'"
    return None


def flavour_pairs():
    rep = _standard_rep()
    rng = np.random.default_rng(3)
    groups = []
    for G in _subgroups():
        kind = _kind(G)
        if kind is None:
            continue
        mats = [rep(g) for g in G]
        gen = sum((rng.normal() + 1j*rng.normal())*M for M in mats)
        ev, vec = np.linalg.eig(gen)
        if min(abs(ev[i] - ev[j]) for i in range(3) for j in range(i+1, 3)) < 1e-6:
            continue                                # not three distinct characters
        groups.append({"kind": kind, "G": G, "basis": vec/np.linalg.norm(vec, axis=0)})
    seen, out = {}, []
    for g1 in groups:
        for g2 in groups:
            U = g1["basis"].conj().T @ g2["basis"]
            P = np.round(np.abs(U)**2, 9)
            canon = min(tuple(tuple(P[r][list(cp)]) for r in rp) for rp in permutations(range(3)) for cp in permutations(range(3)))
            J = abs((U[0, 0]*U[1, 1]*np.conj(U[0, 1])*np.conj(U[1, 0])).imag)
            key = (g1["kind"], g2["kind"], canon, g1["G"] == g2["G"])
            if key in seen:
                continue
            seen[key] = True
            flat = [v for row in canon for v in row]
            clauses = set()
            if any(0 < v < 1 for v in flat):
                clauses.add("mixes")
            if all(v > 1e-9 for v in flat):
                clauses.add("no_zero_entry")
            if any(all(abs(canon[r][c] - 1/3) < 1e-9 for r in range(3)) for c in range(3)):
                clauses.add("trimaximal_column")
            if J > 1e-9:
                clauses.add("cp_odd")
            if g1["G"] == g2["G"]:
                clauses.add("same_group")
            frac = tuple(tuple(str(Fr(v).limit_denominator(48)) for v in row) for row in canon)
            out.append({"pair": f"{g1['kind']} x {g2['kind']}", "index": 24//len(g1["G"]) + 24//len(g2["G"]),
                        "pattern": frac, "jarlskog": round(float(J), 9), "clauses": clauses})
    return out


PRINCIPLES = {
    "hardy": "the kernel restricts to the Hardy kernel on every rank-one slice",
    "covariance": "the kernel is covariant under the conformal group U(2,2)",
    "positivity": "kernel positive, spin-2 kernel Kaellen-Lehmann positive, seam free of anyons (unitarity everywhere)",
    "dirac_spin2": "the spin-2 kernel keeps the Dirac heat-kernel correction a = 3/10",
    "mirror_sheet": "the second sheet half is a mirror copy: seam orientation reversed AND hidden spin-2 channel of equal weight",
    "real_seam": "the seam carries Spin(10) at level one (real Majorana layers)",
    "e8_seam": "the seam has chiral central charge 8",
    "family_mixing": "two residual family symmetries mix all three families (no zero entry)",
    "family_clock": "one residual family symmetry is the quarter-turn clock Z4 (characters -1, i, -i)",
    "orientation_reversal": "the second sheet half reverses the chiral orientation (8 - 8 = 0)",
    "one_transfer": "the spin-2 kernel is the memory of the same quarter-turn Hardy transfer that gives the bundle and the index",
}
TRANSFER_PRINCIPLES = ("orientation_reversal", "one_transfer")      # added with the transfer universes (z45)


def composed_process(transfer=True):
    """Whole processes: one kernel, one spin-2 kernel, one seam with c in {4, 5, 8, 16}, one flavour pair."""
    ks = matrix_kernels()
    s2 = spin2_kernels(transfer)
    sl = [o for o in seam_layers() if o["n"] in (8, 10, 16, 32)]
    fp = [f for f in flavour_pairs() if f["pair"].split(" x ")[0] <= f["pair"].split(" x ")[1]]   # unordered: (A, B) ~ (B, A)^T
    out = []
    for k in ks:
        for t in s2:
            for z in sl:
                for f in fp:
                    c = set()
                    if "hardy_slice" in k["clauses"]:
                        c.add("hardy")
                    if "conformally_covariant" in k["clauses"]:
                        c.add("covariance")
                    if "positive" in k["clauses"] and "kallen_lehmann" in t["clauses"] and "anyon_free" in z["clauses"]:
                        c.add("positivity")
                    if "dirac_low_energy" in t["clauses"]:
                        c.add("dirac_spin2")
                    if "sheet_neutral" in z["clauses"] and "mirror_weight" in t["clauses"]:
                        c.add("mirror_sheet")
                    if "spin10_level_one" in z["clauses"]:
                        c.add("real_seam")
                    if "c_eight" in z["clauses"]:
                        c.add("e8_seam")
                    if {"mixes", "no_zero_entry"} <= f["clauses"]:
                        c.add("family_mixing")
                    if "Z4" in f["pair"]:
                        c.add("family_clock")
                    if transfer and "sheet_neutral" in z["clauses"]:
                        c.add("orientation_reversal")
                    if "transfer_memory" in t["clauses"]:
                        c.add("one_transfer")
                    out.append({"clauses": c, "size": float(k["nu"]) + t["poles"] + z["n"] + f["index"],
                                "label": f"kernel nu={k['nu']} eps={k['eps']} | spin2 {t['name']} | seam n={z['n']} {z['kind']} "
                                         f"{z['orient']} | flavour {f['pair']} {f['pattern']}"})
    return out


RECORD = {
    "matrix_kernels": {"nu=1 eps=0": "the Hardy determinant kernel det(I - W^+ Z)^-1 (z40/z42: massless, Wallach edge, oscillator pairs)"},
    "transfer_processes": {"d=1 k=4 q=1/4": "the quarter-turn transfer a_q(z) = z^-1(1 - q z^4) of the eleventh analysis (z44)"},
    "readouts": {},
    "spin2_kernels": {"transfer memory q=q_dirac": "Marchenko-Pastur memory of the quarter transfer with q^2/(1-q^2)^2 = 3/10 (z45)",
                      "Schur w=1/2 m=5/3": "positive memory kernel C = 1/x + 1/(x + 10/3) of H6h23/z42 (keeps a = 3/10)",
                      "Einstein H = x": "pure Einstein graviton (Jacobson reading, premise D)"},
    "seam_layers": {"n=16 real reversed": "SEAM.S3 (p+ip x 16, c = 8) with an orientation-reversing second sheet half (z41/z42)",
                    "n=8 complex reversed": "the complex flux seam G' of z39 (U(1) charge, no Spin(10)_1)"},
    "flavour_pairs": {"(('0', '1/3', '2/3'), ('1/2', '1/3', '1/6'), ('1/2', '1/3', '1/6'))":
                      "tri-bimaximal mixing (S4 leading order of the TFPT flavour channels, K18)"},
}


def run(root=None, transfer=True):
    """transfer=False reproduces the lane as of H6h23 (z43): no transfer universes, no transfer memory, nine principles."""
    out = {"status": "PROCESS_RULES_SEARCHED", "physical_solution": False, "scope": SCOPE, "universes": {}}
    mk = matrix_kernels()
    res = search(mk, ["positive", "hardy_slice", "wave_equation", "conformally_covariant"],
                 lambda o: float(o["nu"]), lambda o: f"nu={o['nu']} eps={o['eps']}")
    res["derived"] = {"negative_channels": {f"nu={o['nu']} eps={o['eps']}": o["negative_channels"] for o in mk if "positive" not in o["clauses"]}}
    res["bound"] = "nu in {0, 1/2, 1, 3/2, 2, 3}, eps in {0, 1}; Schur channels up to degree 8"
    out["universes"]["matrix_kernels"] = res

    sk = spin2_kernels(transfer)
    res = search(sk, ["no_euclidean_zero", "kallen_lehmann", "dirac_low_energy", "mirror_weight", "rational"]
                 + ["transfer_memory"]*transfer,
                 lambda o: o["poles"], lambda o: o["name"])
    res["derived"] = {o["name"]: {"a": o["a"], "first_zero_near": o["first_zero_near"], "clauses": sorted(o["clauses"])} for o in sk}
    res["bound"] = ("Dirac spectral action with e^-u, step and (1+u)^-3; x e^x; Stelle and Dirac local truncations; Einstein; "
                    "Schur H = x - w x^2/(x + m) on a 3 x 4 grid and one two-channel kernel; x in [1e-3, 1e4]")
    out["universes"]["spin2_kernels"] = res

    if transfer:
        tp = transfer_processes()
        res = search(tp, ["quarter_covariant", "index_one", "kernel_one_grade", "rational_pole_weight", "tau_i_nome"],
                     lambda o: (o["d"], o["k"]), lambda o: f"d={o['d']} k={o['k']} q={o['q']}")
        res["bound"] = "a(z) = z^-d (1 - q z^k), d in {1, 2}, k in {2, 3, 4, 8}, q in {1/4, 1/2, 3/4, e^-pi, e^-2pi}"
        out["universes"]["transfer_processes"] = res

        ro = readouts()
        res = search(ro, ["massless_pole", "kl_positive", "connected", "relativistic", "several_sources", "universal_residue"],
                     lambda o: o["size"], lambda o: f"q={o['q']} {o['kind']}")
        pool = ["massless_pole", "kl_positive", "connected", "relativistic", "several_sources", "universal_residue"]
        empty = [set(r) for k in range(1, len(pool) + 1) for r in combinations(pool, k)
                 if not any(set(r) <= o["clauses"] for o in ro)]
        res["minimal_empty_rules"] = sorted(sorted(r) for r in empty if not any(e < r for e in empty))
        res["bound"] = "readouts of the primitive quarter transfer at q = 1/2 and e^-pi"
        out["universes"]["readouts"] = res

    sl = seam_layers()
    res = search(sl, ["c_eight", "anyon_free", "spin10_level_one", "sheet_neutral", "u1_free"],
                 lambda o: o["n"], lambda o: f"n={o['n']} {o['kind']} {o['orient']}")
    res["bound"] = "1 to 32 layers, real or complex, second sheet half equally or reversely oriented"
    out["universes"]["seam_layers"] = res

    fp = flavour_pairs()
    res = search(fp, ["mixes", "no_zero_entry", "trimaximal_column", "cp_odd", "same_group"],
                 lambda o: o["index"], lambda o: str(o["pattern"]))
    res["derived"] = {"patterns": sorted({(o["pair"], str(o["pattern"]), o["jarlskog"]) for o in fp if "mixes" in o["clauses"]})}
    res["bound"] = "ordered pairs of abelian subgroups of S4 with three distinct characters on the standard triplet"
    out["universes"]["flavour_pairs"] = res

    principles = {k: v for k, v in PRINCIPLES.items() if transfer or k not in TRANSFER_PRINCIPLES}
    cp = composed_process(transfer)
    res = search(cp, list(principles), lambda o: o["size"], lambda o: o["label"])
    longest = max(len(r["minimal_rules"][0]) for r in res["unique_selections"]) if res["unique_selections"] else 0
    full = [r for r in res["unique_selections"] if len(r["minimal_rules"][0]) == longest]
    out["composed_process"] = {"principles": principles, "objects": res["objects"], "rules": res["rules"],
                               "verdicts": res["verdicts"], "longest_minimal_rule": longest,
                               "selected_by_the_longest_rules": full,
                               "unique_selections": res["unique_selections"][:40],
                               "bound": "all kernels, all spin-2 kernels, seams with 8, 10, 16, 32 layers, unordered flavour pairs"}
    met = []
    for key, res in out["universes"].items():
        for sel in res["unique_selections"]:
            for pattern, meaning in RECORD[key].items():
                if pattern == sel["selected"]:
                    met.append({"universe": key, "selected": sel["selected"], "shortest_rule": sel["minimal_rules"][0],
                                "record": meaning})
    out["selections_that_meet_the_record"] = met
    out["rules_generated"] = sum(r["rules"] for r in out["universes"].values())
    out["objects_enumerated"] = sum(r["objects"] for r in out["universes"].values())
    return out


if __name__ == "__main__":
    argparse.ArgumentParser().parse_args()
    print(json.dumps(run(), indent=1, default=str, ensure_ascii=False))
