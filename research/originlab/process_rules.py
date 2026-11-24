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
uses the historical H6h23 candidate switch with current certificates: no transfer universes, no transfer memory, nine composed principles.

The original finite comparison also uses thermal and flavour inputs. In particular, `dirac_low_energy` asks for the first correction a = 3/10 of the Dirac heat kernel
(a structural property of the Dirac operator with the e^-u cutoff, z41). Record structures are named afterwards in
`RECORD`. A unique minimum is conditional on the candidate universe, clause pool and cost. The source_dynamics audit adds
exact operator and interaction checks, including positive parent alternatives. No physical map is promoted.
"""
from fractions import Fraction as Fr
from functools import lru_cache
from itertools import combinations, permutations
from math import gamma, gcd
import argparse
import json

import numpy as np
import sympy
from numpy.polynomial.legendre import leggauss
from scipy.special import expit, roots_genlaguerre

from genesis_rules import search
from process_matrix_certificates import matrix_certificate
from process_spin2_certificates import candidate_certificate, rational_certificate
from process_transfer_certificates import transfer_certificate
from process_metric_potential_certificates import metric_potential_certificate

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
            cert = matrix_certificate(nu, eps)
            clauses = set(cert["clauses"])
            out.append({"nu": nu, "eps": eps, "clauses": clauses,
                        "negative_channels": sorted(str(l) for l, c in coeff.items() if c < 0)[:4],
                        "certificate": cert})
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


def original_current_norm(alpha=0.25, b=4*np.pi, cut=40):
    """omega(J_-1 J_1) in the original labelled quarter Gaussian state: c_r = (1 + e^{b(r + alpha)})^-1, b = 4 pi
    (ELLIPTIC-CURRENT-SOURCE, native-arc-charge-source-interface-20261026, eq. 1 and 6), by the direct Fermi sum
    Sum_r c_{r+1}(1 - c_r). The Hardy state has J_1 Omega = 0, i.e. norm 0."""
    r = np.arange(-cut, cut, dtype=float)
    c = lambda rr: expit(-b*(rr + alpha))
    return float(np.sum(c(r + 1)*(1 - c(r))))


def ladder_occupation(qv, n=200):
    """Mean ladder level <n> of the normalised kernel vector of S^+(1 - q S^4) (computed from the compression)."""
    if qv == 0:
        return 0.0
    A = _compression(1, 4, qv, n=n)
    vec = np.linalg.svd(A)[2][-1]
    w = np.abs(vec)**2/np.sum(np.abs(vec)**2)
    return float(np.sum(np.arange(n)/4*w))


Q_ORIGINAL = float(np.exp(-2*np.pi))                # ladder occupation q^2/(1-q^2) = omega(J_-1 J_1) = 1/(e^{4 pi} - 1)


def spin2_kernels(transfer=True, state=True):
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
    memories = (("1/4", 0.25), ("1/2", 0.5), ("3/4", 0.75), ("q_dirac", Q_DIRAC), ("e^-pi", float(np.exp(-np.pi))))
    if state:
        memories += (("e^-2pi", Q_ORIGINAL), ("1/sqrt2", float(np.sqrt(0.5))))
    for label, qv in memories[:len(memories)*transfer]:
        Cq = transfer_kernel(qv)
        # 1/H - 1/x = 2 q^2 / (Z [sqrt((x+a)(x+b)) + x + Z]): the closed form avoids cancellation for small q (z46)
        R = (lambda qq, Z: lambda x: 2*qq**2/(Z*(np.sqrt((np.asarray(x, float) + (1 - qq)**2)*(np.asarray(x, float) + (1 + qq)**2))
                                                + np.asarray(x, float) + Z)))(qv, 1 - qv**2)
        cands.append({"name": f"transfer memory q={label}", "H": (lambda C, Z: lambda x: Z/C(x))(Cq, 1 - qv**2), "poles": 50,
                      "schur": None, "transfer": True, "q": qv, "R": R})
    xs = np.geomspace(1e-3, 1e4, 900)
    out = []
    for c in cands:
        H = np.asarray(c["H"](xs), float)
        kinds = {"Einstein H = x": "einstein", "ledger ansatz H = x e^x": "ledger_entire",
                 "Stelle local H = x(1+x)": "stelle", "Dirac local H = x(1-3x/10)": "dirac_local",
                 "spectral action e^-u": "exp", "spectral action step": "step",
                 "spectral action (1+u)^-3": "laplace3"}
        if c.get("transfer"):
            label = c["name"].split("q=", 1)[1]
            exact_q = sympy.sqrt(sympy.Rational(1, 2)) if label == "1/sqrt2" else label
            cert = transfer_certificate(exact_q)
        elif isinstance(c["schur"], tuple):
            w, m = c["schur"]
            cert = rational_certificate((w,), (m,))
        elif c["schur"] == "two":
            cert = rational_certificate((Fr(1, 4), Fr(1, 4)), (Fr(5, 3), Fr(10, 3)))
        else:
            cert = candidate_certificate(kinds[c["name"]])
        clauses = set(cert["clauses"])
        if c.get("transfer"):
            if state and label == "e^-2pi":
                clauses.add("original_state")
        zero = None
        if np.any(H <= 0):
            k = int(np.argmax(H <= 0))
            zero = float(xs[k])
        out.append({"name": c["name"], "search_cost": c["poles"],
                    "poles": 1 if c.get("transfer") else (None if c["poles"] == 99 else c["poles"]),
                    "continuum": bool(c.get("transfer")), "clauses": clauses,
                    "a": float(cert["a"]), "a_exact": str(cert["a"]),
                    "b_exact": str(cert["b"]), "c_exact": str(cert["c"]),
                    "certificate": _json_exact(cert), "first_zero_near": zero,
                    "sample_scope": "first nonpositive sampled value, not a certified root",
                    "cost_scope": "50 and 99 are declared search penalties, not numbers of poles"})
    return out


def _json_exact(value):
    if isinstance(value, (Fr, sympy.Basic)):
        return str(value)
    if isinstance(value, dict):
        return {key: _json_exact(item) for key, item in value.items()}
    if isinstance(value, set):
        return sorted(value)
    if isinstance(value, (tuple, list)):
        return [_json_exact(item) for item in value]
    return value


# ---------------------------------------------------------------- B2: transfer functions and readouts
Q_VALUES = (("1/4", 0.25), ("1/2", 0.5), ("3/4", 0.75), ("e^-pi", float(np.exp(-np.pi))), ("e^-2pi", float(np.exp(-2*np.pi))))


def _compression(d, k, qv, n=400):
    S = np.diag(np.ones(n - 1), -1)
    return np.linalg.matrix_power(S.T, d) @ (np.eye(n) - qv*np.linalg.matrix_power(S, k))


def transfer_processes():
    """Exact identities for A=S*^d(1-qS^k), 0<q<1; no sampled-rank inference.

    ker(A)=(1-qS^k)^(-1) span(e_0,...,e_(d-1)); the cokernel is zero.
    a(iz)/a(z) is constant exactly when k is divisible by four. Kernel
    grades are the union of j+m*k modulo four for 0<=j<d.
    """
    out = []
    for d in (1, 2):
        for k in (2, 3, 4, 8):
            grades = sorted({(j+m*k) % 4 for j in range(d) for m in range(4//gcd(k, 4))})
            for qi, (label, q) in enumerate(Q_VALUES):
                clauses = set()
                if k % 4 == 0:
                    clauses.add("quarter_covariant")
                if d == 1:
                    clauses.add("index_one")
                if len(grades) == 1:
                    clauses.add("kernel_one_grade")
                if label in ("1/4", "1/2", "3/4"):
                    clauses.add("rational_pole_weight")
                if label.startswith("e^"):
                    clauses.add("tau_i_nome")
                out.append({"d": d, "k": k, "q": label, "qi": qi, "clauses": clauses,
                            "kernel_dim": d, "kernel_grades": grades,
                            "certificate": {"method": "Hardy inverse and exact cyclic grading",
                                            "kernel_basis": "sum_(m>=0) q^m e_(j+m*k), 0<=j<d",
                                            "hardy_cokernel_dimension": 0,
                                            "q_selected": False}})
    return out


def readouts(metric=False):
    """Declared source maps, with Hamilton vacuum subtraction separate from field centering.

    Centering a Hamilton-vacuum source removes its zero-state atom. Centering
    a Gaussian spacetime field leaves its connected two-point covariance.
    Equal residues obtained by separately normalizing sources are inputs,
    not a gravitational Ward identity or a native source-map derivation.
    """
    out = []
    for label, q in (("1/2", sympy.Rational(1, 2)), ("e^-pi", sympy.exp(-sympy.pi))):
        Z = 1-q*q
        raw = [Z*q**(2*j) for j in range(4)]
        kinds = (
            ("boundary source |0>", Z, [Z], False, False, False, "boundary_resolvent", False),
            ("connected vacuum response", 0, [0], True, False, False, "Hamilton_vacuum", False),
            ("4D field, one source", Z, [Z], False, True, False, "Gaussian_field", False),
            ("4D field, sources |4j>", Z, raw, False, True, True, "Gaussian_field", False),
            ("4D field connected, sources |4j>", Z, raw, True, True, True, "Gaussian_field", False),
            ("4D field connected, pole-normalized sources", Z, [Z]*4, True, True, True, "Gaussian_field", True),
        )
        for kind, pole, residues, connected, relativistic, multiple, interpretation, normalized in kinds:
            clauses = {"kl_positive"}
            if pole != 0:
                clauses.add("massless_pole")
            if connected:
                clauses.add("connected")
            if metric and interpretation == "Hamilton_vacuum":
                clauses.discard("connected")
                clauses.add("vacuum_subtracted")
            if metric and connected and interpretation == "Gaussian_field":
                clauses.add("field_connected")
            if relativistic:
                clauses.add("relativistic")
            if multiple:
                clauses.add("several_sources")
                if pole != 0 and all(sympy.simplify(r-residues[0]) == 0 for r in residues):
                    clauses.add("universal_residue")
            out.append({"q": label, "kind": kind, "clauses": clauses, "size": len(residues),
                        "interpretation": interpretation, "pole_exact": str(pole),
                        "residues_exact": [str(sympy.simplify(r)) for r in residues],
                        "source_normalization_is_input": normalized,
                        "physical_source_map_derived": False, "gravitational_ward_identity": False,
                        "scope": "A declared positive scalar readout; relativistic embedding and source conventions are additional inputs."})
        if metric:
            # B M = 0 and the connected tree has one kernel direction exactly.
            # This is a consequence of the supplied metric action and its
            # minimal source coupling, not a source map derived from TFPT.
            arm_length = 6
            meff2 = sympy.simplify(sum(q**(2*n) for n in range(arm_length+1)))
            out.append({
                "q": label, "kind": "metric star, minimal node sources",
                "size": 1+4*arm_length,
                "clauses": {"kl_positive", "relativistic", "connected", "field_connected",
                            "several_sources", "massless_pole", "universal_residue"},
                "interpretation": "metric_field", "pole_exact": str(1/meff2),
                "residues_exact": [str(1/meff2)]*(1+4*arm_length),
                "source_normalization_is_input": False, "source_action_is_input": True,
                "physical_source_map_derived": False, "gravitational_ward_identity": False,
                "certificate": {"nullity": 1, "identity": "B M = 0; eta_v/M_v = 1/M_eff",
                                "M_eff_squared_at_Mc_one": str(meff2)},
                "scope": "Universal massless coupling follows from the assumed minimally coupled metric tree; the microscopic stress tensor dictionary remains open.",
            })
    return out


# ---------------------------------------------------------------- B3: metric star and potentials (z50)
def metric_star(qv, m, Mc=1.0):
    """Central metric plus four arms of m metrics; Einstein weights M_c and M_(j,n) = (M_c/2) q^(n+1); edge stiffness
    mu^2 M_child^2. Canonical mass operator mu^2 B^+ B with (Bh)_(j,0) = (q/2) h_c - h_(j,0), (Bh)_(j,n) = q h_(j,n-1) - h_(j,n)."""
    nodes = ["c"] + [(j, n) for j in range(4) for n in range(m)]
    idx = {v: i for i, v in enumerate(nodes)}
    B = np.zeros((4*m, 1 + 4*m))
    for e, (j, n) in enumerate((j, n) for j in range(4) for n in range(m)):
        B[e, idx["c"] if n == 0 else idx[(j, n - 1)]] = qv/2 if n == 0 else qv
        B[e, idx[(j, n)]] = -1
    M = np.array([Mc] + [Mc/2*qv**(n + 1) for j in range(4) for n in range(m)])
    _, sv, vh = np.linalg.svd(B)
    # the kernel vector is M itself (B M = 0 exactly); the SVD vector loses relative precision on deep nodes for small q
    assert np.max(np.abs(B @ M)) < 1e-14*np.max(np.abs(M))
    eta = M/np.linalg.norm(M)
    return {"B": B, "M": M, "nodes": nodes, "eta": eta, "nullity": 1, "sv": sv,
            "nullity_certificate": "The child columns form an invertible triangular matrix; B has full row rank.",
            "svd_nullity_diagnostic": int(B.shape[1] - np.sum(sv > 1e-12)),
            "svd_kernel_overlap": float(abs(vh[-1] @ eta))}


def hardy_rect(qv, m):
    """Finite rectangular quarter transfer: rows 0..4m-1, columns 0..4m, (A h)_k = h_(k+1) - q h_(k-3)."""
    A = np.zeros((4*m, 4*m + 1))
    for k in range(4*m):
        A[k, k + 1] = 1
        if k >= 3:
            A[k, k - 3] = -qv
    return A


BETA_V = (-2.0, 0.5, 0.0, 0.5, -2.0)          # V(X) = (-4 + e1 + e3 - 4 e4)/2: flat for both metrics, unit Fierz-Pauli mass
BETA_DET = (1.0, -1.0, 1.0, -1.0, 1.0)        # det(I - X) = Sum (-1)^n e_n(X): starts at fourth order
BETA_CUB = (-2.0, 1.0, 0.0, -1.0, 2.0)        # keeps flatness and mass, changes the cubic vertex


def potential_taylor(beta, A, B, order=4, dps=40):
    """Taylor coefficients in eps of sqrt(det g) Sum_n beta_n e_n(sqrt(g^-1 f)), g = I + eps A, f = I + eps B (Euclidean)."""
    import mpmath as mp
    mp.mp.dps = dps
    Am, Bm = mp.matrix(A.tolist()), mp.matrix(B.tolist())

    def F(eps):
        g = mp.eye(4) + eps*Am
        f = mp.eye(4) + eps*Bm
        w, Vg = mp.eigsy(g)
        gih = Vg*mp.diag([1/mp.sqrt(x) for x in w])*Vg.T
        lam, _ = mp.eigsy(gih*f*gih)
        x = [mp.sqrt(v) for v in lam]
        e = [mp.mpf(1), sum(x), sum(x[i]*x[j] for i in range(4) for j in range(i + 1, 4)),
             sum(x[i]*x[j]*x[k] for i in range(4) for j in range(i + 1, 4) for k in range(j + 1, 4)), x[0]*x[1]*x[2]*x[3]]
        return mp.sqrt(mp.fprod(w))*sum(mp.mpf(b)*en for b, en in zip(beta, e))
    return [float(c) for c in mp.taylor(F, 0, order)]


@lru_cache(maxsize=1)
def metric_potentials():
    """Exact HR coefficient certificates; finite Taylor samples are diagnostics."""
    rng = np.random.default_rng(50)
    pairs = []
    for _ in range(2):
        A, B = rng.normal(size=(4, 4)), rng.normal(size=(4, 4))
        pairs.append(((A + A.T)/2, (B + B.T)/2))
    cands = (("V", BETA_V), ("V + det(I-X)/2", tuple(v + d/2 for v, d in zip(BETA_V, BETA_DET))),
             ("V - det(I-X)/2", tuple(v - d/2 for v, d in zip(BETA_V, BETA_DET))),
             ("V + cubic/2", tuple(v + c/2 for v, c in zip(BETA_V, BETA_CUB))),
             ("V + cubic/2 + det(I-X)/2", tuple(v + c/2 + d/2 for v, c, d in zip(BETA_V, BETA_CUB, BETA_DET))),
             ("beta1 only", (0.0, 1.0, 0.0, 0.0, 0.0)))
    out = []
    for name, beta in cands:
        co = [potential_taylor(beta, A, B) for A, B in pairs]
        # The historical literal halves have a declared exact rational meaning.
        exact_beta = tuple(sympy.Rational(str(b)) for b in beta)
        certificate = metric_potential_certificate(exact_beta)
        out.append({"name": name, "beta": beta, "clauses": set(certificate["clauses"]),
                    "taylor": co, "certificate": certificate,
                    "sample_scope": "Two perturbation pairs are diagnostics; they decide no clause."})
    return out


def q_rules():
    """Exact candidate thermal weights; a shared occupation is not a full KMS proof.

    The original_state label means the candidate matches the specified canonical
    thermal purification. Its physical current and time hypotheses are recorded
    in source_dynamics, rather than inferred from one mean occupation.
    """
    cands = (("0 (Hardy state)", sympy.Integer(0)), ("e^-2pi", sympy.exp(-2*sympy.pi)),
             ("e^-4pi", sympy.exp(-4*sympy.pi)), ("e^-pi", sympy.exp(-sympy.pi)),
             ("1/sqrt2", sympy.sqrt(2)/2), ("q_dirac", sympy.sqrt((8-sympy.sqrt(55))/3)),
             ("1/2", sympy.Rational(1, 2)))
    out = []
    for label, q in cands:
        r = sympy.simplify(q*q)
        occ, residue = sympy.simplify(r/(1-r)), sympy.simplify(1-r)
        clauses = set()
        if label == "e^-2pi":
            clauses.add("original_state")
        if q == 0:
            clauses.add("hardy_state")
        else:
            clauses.add("memory")
        if sympy.simplify(r-sympy.Rational(1, 2)) == 0:
            clauses.add("equal_weight")
        if bool(r < sympy.Rational(1, 100)):
            clauses.add("einstein_dominant")
        if label == "q_dirac":
            clauses.add("dirac_same_scale")
        if q.func == sympy.exp:
            clauses.add("modular_nome")
        out.append({"q": label, "value": float(q), "clauses": clauses,
                    "occupation": float(occ), "residue": float(residue),
                    "occupation_exact": str(occ), "residue_exact": str(residue),
                    "scope": "canonical-purification parameter certificate under the stated current KMS dictionary"})
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
                    # An even unimodular D_(majorana/2)^+ extension exists.
                    # SO(16)_1 itself still has four sectors: no condensation is performed here.
                    clauses.add("bosonic_extension_available")
                if kind == "real" and n >= 10:
                    clauses.add("spin10_level_one")         # sufficient ten-real-mode current-algebra realization
                if orient == "reversed":
                    clauses.add("sheet_neutral")            # c - c = 0 across the two sheet halves
                if kind == "real":
                    clauses.add("u1_free")
                out.append({"n": n, "kind": kind, "orient": orient, "c": str(c), "clauses": clauses,
                            "condensation_constructed": False,
                            "spin10_scope": "existence of a ten-real-mode current subalgebra, not a bulk matter identification; complex alternatives depend on the retained charge structure",
                            "orientation_scope": "net chiral-central-charge cancellation only; no protection of coexisting opposite edges is proved"})
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


def _exact(v):
    """Exact label in Q(sqrt 2, sqrt 3); the first run used limit_denominator(48), which printed (2 - sqrt 3)/6 as 2/45."""
    e = sympy.nsimplify(float(v), [sympy.sqrt(2), sympy.sqrt(3)], tolerance=1e-8)
    small = (len(sympy.Add.make_args(e)) <= 2
             and all(r.q <= 12 and abs(r.p) <= 12 for r in e.atoms(sympy.Rational)))
    if abs(float(e) - float(v)) > 1e-8 or not small:          # no spurious large-denominator fits (z49)
        return repr(round(float(v), 9))
    return str(e)


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
            frac = tuple(tuple(_exact(v) for v in row) for row in canon)
            out.append({"pair": f"{g1['kind']} x {g2['kind']}", "index": 24//len(g1["G"]) + 24//len(g2["G"]),
                        "pattern": frac, "jarlskog": round(float(J), 9), "clauses": clauses})
    return out


# ---------------------------------------------------------------- D2: larger family groups (z49)
DATA_RANGES = ((0.27, 0.35), (0.020, 0.024), (0.43, 0.60))     # approximate 3-sigma, NuFIT 6.0 normal ordering, rounded (z43)


def delta_group(n):
    """Delta(6 n^2) in its faithful 3-dim representation: monomial matrices (perm, phase exponents of e^{i pi/n})."""
    N2 = 2*n

    def mul(g, h):
        pg, xg = g
        ph, xh = h
        return (tuple(pg[ph[j]] for j in range(3)), tuple((xh[j] + xg[ph[j]]) % N2 for j in range(3)))
    gens = [((1, 2, 0), (0, 0, 0)), ((2, 1, 0), (n, n, n)), ((0, 1, 2), (2, N2 - 2, 0)), ((0, 1, 2), (0, 2, N2 - 2))]
    e = ((0, 1, 2), (0, 0, 0))
    G, frontier = {e}, [e]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = mul(g, x)
                if y not in G:
                    G.add(y)
                    new.append(y)
        frontier = new
    mat = lambda g: np.array([[np.exp(1j*np.pi*g[1][j]/n) if g[0][j] == i else 0 for j in range(3)] for i in range(3)])
    return G, mul, e, mat


def _joint_basis(mats, seed=1):
    rng = np.random.default_rng(seed)
    gen = sum((rng.normal() + 1j*rng.normal())*M for M in mats)
    ev, vec = np.linalg.eig(gen)
    if min(abs(ev[i] - ev[j]) for i in range(3) for j in range(i+1, 3)) < 1e-6:
        return None
    return vec/np.linalg.norm(vec, axis=0)


def _angles_fit(P):
    """Some row/column assignment puts sin^2 theta12, theta13, theta23 into DATA_RANGES."""
    (r12, r13, r23) = DATA_RANGES
    for rp in permutations(range(3)):
        for cp_ in permutations(range(3)):
            M = [[P[r][c] for c in cp_] for r in rp]
            s13 = M[0][2]
            if s13 >= 1:
                continue
            s12, s23 = M[0][1]/(1 - s13), M[1][2]/(1 - s13)
            if r12[0] <= s12 <= r12[1] and r13[0] <= s13 <= r13[1] and r23[0] <= s23 <= r23[1]:
                return True
    return False


def flavour_groups(ns=(2, 3, 4, 6, 8)):
    """Residual symmetries in Delta(6 n^2): Klein group for Majorana neutrinos, cyclic group with three distinct
    eigenvalues for charged leptons; |U|^2 of the two eigenbases."""
    out = []
    clock = None
    for n in ns:
        G, mul, e, mat = delta_group(n)
        els = sorted(G)
        M = {g: mat(g) for g in els}
        order = {}
        for g in els:
            k, x = 1, g
            while x != e:
                x, k = mul(g, x), k + 1
            order[g] = k
        inv = [g for g in els if order[g] == 2 and not np.allclose(M[g], -np.eye(3))]
        klein, seen = [], set()
        for g1, g2 in combinations(inv, 2):
            if mul(g1, g2) == mul(g2, g1):
                key = frozenset([g1, g2, mul(g1, g2)])
                B = _joint_basis([M[g1], M[g2]])
                if key not in seen and B is not None:
                    seen.add(key)
                    klein.append(B)
        cyc = [(order[g], B) for g in els if order[g] >= 3 for B in [_joint_basis([M[g]])] if B is not None]
        pats = {}
        for Bn in klein:
            for o, Bl in cyc:
                P = np.abs(Bl.conj().T @ Bn)**2
                canon = min(tuple(tuple(round(float(P[r][c]), 9) for c in cp_) for r in rp)
                            for rp in permutations(range(3)) for cp_ in permutations(range(3)))
                pats.setdefault(canon, set()).add(o)
        for canon, orders in sorted(pats.items()):
            flat = [v for row in canon for v in row]
            cols = [tuple(canon[r][c] for r in range(3)) for c in range(3)]
            clauses = set()
            if all(v > 1e-9 for v in flat):
                clauses.add("mixes")
            if any(sorted(col) == sorted((2/3, 1/6, 1/6)) or np.allclose(sorted(col), sorted((2/3, 1/6, 1/6)), atol=1e-8)
                   for col in cols):
                clauses.add("tm1")
            if any(np.allclose(col, 1/3, atol=1e-8) for col in cols):
                clauses.add("tm2")
            if any(DATA_RANGES[1][0] <= v <= DATA_RANGES[1][1] for v in flat):
                clauses.add("theta13_window")
            if _angles_fit(canon):
                clauses.add("fits_data")
            if n % 4 == 0:
                clauses.add("quarter_phases")
            if any(abs(v - (2 - np.sqrt(3))/6) < 1e-8 for v in flat) and "tm2" in clauses:
                clauses.add("clock_pattern")                 # Z3 x Z4 of the S4 lane = TM2 at 15 degrees (z46)
            out.append({"n": n, "group": f"Delta({6*n*n})", "pattern": tuple(tuple(_exact(v) for v in row) for row in canon),
                        "charged_lepton_orders": sorted(orders), "clauses": clauses})
    return out


PRINCIPLES = {
    "hardy": "the kernel restricts to the Hardy kernel on every rank-one slice",
    "covariance": "the kernel is covariant under the conformal group U(2,2)",
    "positivity": "positive matrix kernel, positive scalar TT spectral measure, and availability of a bosonic seam extension; condensation and full quantum unitarity are not certified",
    "dirac_spin2": "the spin-2 kernel keeps the Dirac heat-kernel correction a = 3/10 (in units of the e^-u cutoff: a scale "
                   "convention, z46)",
    "mirror_sheet": "reversed seam orientation AND equality of total massless and massive scalar spectral weights; no microscopic mirror map is inferred",
    "real_seam": "the seam carries Spin(10) at level one (real Majorana layers)",
    "e8_seam": "the seam has chiral central charge 8",
    "family_mixing": "two residual family symmetries mix all three families (no zero entry)",
    "family_clock": "one residual family symmetry is the quarter-turn clock Z4 (characters -1, i, -i)",
    "orientation_reversal": "the second sheet half reverses the chiral orientation (8 - 8 = 0)",
    "one_transfer": "membership in the declared flat Hardy response family; identification as the physical spin-2 process is conditional",
    "original_state": "q matches the canonical purification of the specified lowest current KMS mode; current algebra, time and conjugate-sheet dictionary are hypotheses",
}
TRANSFER_PRINCIPLES = ("orientation_reversal", "one_transfer")      # added with the transfer universes (z45)
STATE_PRINCIPLES = ("original_state",)                               # added with the state rule for q (z46)


def allowed_principles(transfer=True, state=True):
    return {k: v for k, v in PRINCIPLES.items()
            if (transfer or k not in TRANSFER_PRINCIPLES) and (state and transfer or k not in STATE_PRINCIPLES)}


def composed_process(transfer=True, state=True):
    """Whole processes: one kernel, one spin-2 kernel, one seam with c in {4, 5, 8, 16}, one flavour pair."""
    ks = matrix_kernels()
    s2 = spin2_kernels(transfer, state)
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
                    if "positive" in k["clauses"] and "kallen_lehmann" in t["clauses"] and "bosonic_extension_available" in z["clauses"]:
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
                    if "original_state" in t["clauses"]:
                        c.add("original_state")
                    out.append({"clauses": c, "size": float(k["nu"]) + t["search_cost"] + z["n"] + f["index"],
                                "label": f"kernel nu={k['nu']} eps={k['eps']} | spin2 {t['name']} | seam n={z['n']} {z['kind']} "
                                         f"{z['orient']} | flavour {f['pair']} {f['pattern']}"})
    return out


RECORD = {
    "matrix_kernels": {"nu=1 eps=0": "the Hardy determinant kernel det(I - W^+ Z)^-1 (z40/z42: massless, Wallach edge, oscillator pairs)"},
    "transfer_processes": {"d=1 k=4 q=1/4": "the quarter-turn transfer a_q(z) = z^-1(1 - q z^4) of the eleventh analysis (z44)"},
    "readouts": {},
    "q_rules": {"e^-2pi": "transfer state = purification of the lowest current mode of the original quarter Gaussian state (z46)"},
    "spin2_kernels": {"transfer memory q=q_dirac": "Marchenko-Pastur memory of the quarter transfer with q^2/(1-q^2)^2 = 3/10 (z45)",
                      "transfer memory q=e^-2pi": "transfer memory in the original quarter Gaussian state, hidden weight e^-4pi (z46)",
                      "Schur w=1/2 m=5/3": "positive memory kernel C = 1/x + 1/(x + 10/3) of H6h23/z42 (keeps a = 3/10)",
                      "Einstein H = x": "pure Einstein graviton (Jacobson reading, premise D)"},
    "seam_layers": {"n=16 real reversed": "SEAM.S3 (p+ip x 16, c = 8) with an orientation-reversing second sheet half (z41/z42)",
                    "n=8 complex reversed": "the complex flux seam G' of z39 (U(1) charge, no Spin(10)_1)"},
    "flavour_groups": {},
    "metric_potentials": {"V": "Hassan-Rosen potential (-4 + e1 + e3 - 4 e4)/2 of the metric star (z50)"},
    "flavour_pairs": {"(('0', '1/3', '2/3'), ('1/2', '1/3', '1/6'), ('1/2', '1/3', '1/6'))":
                      "tri-bimaximal mixing (S4 leading order of the TFPT flavour channels, K18)"},
}


def run(root=None, transfer=True, state=True, metric=True, source_dynamics=True):
    """Historical candidate switches with current certificates; archived files stay unchanged."""
    metric = metric and transfer and state
    out = {"status": "PROCESS_RULES_SEARCHED", "physical_solution": False, "scope": SCOPE, "universes": {}}
    mk = matrix_kernels()
    res = search(mk, ["positive", "hardy_slice", "wave_equation", "conformally_covariant"],
                 lambda o: float(o["nu"]), lambda o: f"nu={o['nu']} eps={o['eps']}")
    res["derived"] = {"negative_channels": {f"nu={o['nu']} eps={o['eps']}": o["negative_channels"] for o in mk if "positive" not in o["clauses"]}}
    res["bound"] = "nu in {0, 1/2, 1, 3/2, 2, 3}, eps in {0, 1}; all-degree analytic Schur and covariance certificates; diagnostic channels up to degree 8"
    out["universes"]["matrix_kernels"] = res

    sk = spin2_kernels(transfer, state)
    res = search(sk, ["no_euclidean_zero", "kallen_lehmann", "dirac_low_energy", "mirror_weight", "rational"]
                 + ["transfer_memory"]*transfer + ["original_state"]*(transfer and state),
                 lambda o: o["search_cost"], lambda o: o["name"])
    res["derived"] = {o["name"]: {
        **{key: o[key] for key in ("a", "a_exact", "b_exact", "c_exact", "search_cost", "poles",
                                  "continuum", "certificate", "first_zero_near", "sample_scope", "cost_scope")},
        "clauses": sorted(o["clauses"])} for o in sk}
    res["bound"] = ("analytic scalar certificates for Dirac profiles e^-u, step and (1+u)^-3; x e^x; local truncations; "
                    "Einstein; Schur H on the historical finite parameter grid; exact transfer measures. "
                    "A finite x mesh only locates diagnostic sign changes; no tensor or quantum gravity claim.")
    out["universes"]["spin2_kernels"] = res

    if transfer:
        tp = transfer_processes()
        res = search(tp, ["quarter_covariant", "index_one", "kernel_one_grade", "rational_pole_weight", "tau_i_nome"],
                     lambda o: (o["d"], o["k"]), lambda o: f"d={o['d']} k={o['k']} q={o['q']}")
        res["bound"] = "a(z) = z^-d (1 - q z^k), d in {1, 2}, k in {2, 3, 4, 8}, q in {1/4, 1/2, 3/4, e^-pi, e^-2pi}"
        out["universes"]["transfer_processes"] = res

        ro = readouts(metric)
        pool = (["massless_pole", "kl_positive", "vacuum_subtracted", "field_connected", "relativistic", "several_sources",
                 "universal_residue"] if metric else
                ["massless_pole", "kl_positive", "connected", "relativistic", "several_sources", "universal_residue"])
        res = search(ro, pool, lambda o: o["size"], lambda o: f"q={o['q']} {o['kind']}")
        empty = [set(r) for k in range(1, len(pool) + 1) for r in combinations(pool, k)
                 if not any(set(r) <= o["clauses"] for o in ro)]
        res["minimal_empty_rules"] = sorted(sorted(r) for r in empty if not any(e < r for e in empty))
        res["bound"] = "declared microscopic and 4D field sources at q = 1/2 and e^-pi; field centering does not remove a massless particle pole"
        out["universes"]["readouts"] = res

    if metric:
        mp_ = metric_potentials()
        pool = ["flat_background", "fierz_pauli_mass", "massless_direction_free", "undirected_edge", "cubic_as_V", "quartic_as_V"]
        res = search(mp_, pool, lambda o: 0, lambda o: o["name"])
        res["tie_flat_mass_undirected"] = sorted(o["name"] for o in mp_
                                                 if {"flat_background", "fierz_pauli_mass", "undirected_edge"} <= o["clauses"])
        res["derived"] = {o["name"]: {"beta": o["beta"], "clauses": sorted(o["clauses"]),
                                      "certificate": o["certificate"],
                                      "taylor_first_pair": o["taylor"][0]} for o in mp_}
        res["bound"] = "six HR candidates with exact tensor-coefficient certificates in the full five-parameter HR class; two Taylor samples are diagnostic"
        out["universes"]["metric_potentials"] = res
        diag = {}
        for qv in (0.5, Q_ORIGINAL):
            st, A = metric_star(qv, 6), hardy_rect(qv, 6)
            sa = np.linalg.svd(A, compute_uv=False)
            diag[str(qv)] = {"singular_values_agree": bool(np.allclose(np.sort(st["sv"]), np.sort(sa), atol=1e-12)),
                             "nullity": st["nullity"],
                             "eta_over_M_spread": float(np.ptp(st["eta"]/st["M"])/np.mean(st["eta"]/st["M"]))}
        out["metric_star"] = diag

    if transfer and state:
        qr = q_rules()
        pool = ["original_state", "hardy_state", "memory", "equal_weight", "einstein_dominant", "dirac_same_scale",
                "modular_nome"]
        res = search(qr, pool, lambda o: 0, lambda o: o["q"])
        empty = [set(r) for k in range(1, len(pool) + 1) for r in combinations(pool, k)
                 if not any(set(r) <= o["clauses"] for o in qr)]
        res["minimal_empty_rules"] = sorted(sorted(r) for r in empty if not any(e < r for e in empty))
        res["derived"] = {o["q"]: {"occupation": o["occupation"], "massless_residue": o["residue"],
                                   "clauses": sorted(o["clauses"])} for o in qr}
        res["original_current_norm"] = original_current_norm()
        res["bound"] = "seven candidate states of the quarter transfer; no size preference (ties stay AMBIGUOUS)"
        out["universes"]["q_rules"] = res

    if transfer and state:
        fg = flavour_groups()
        pool = ["mixes", "tm1", "tm2", "theta13_window", "fits_data", "quarter_phases", "clock_pattern"]
        res = search(fg, pool, lambda o: o["n"], lambda o: f"{o['group']} {o['pattern']}")
        empty = [set(r) for k in range(1, len(pool) + 1) for r in combinations(pool, k)
                 if not any(set(r) <= o["clauses"] for o in fg)]
        res["minimal_empty_rules"] = sorted(sorted(r) for r in empty if not any(e < r for e in empty))
        res["derived"] = {"patterns": [(o["group"], str(o["pattern"]), sorted(o["clauses"])) for o in fg]}
        res["bound"] = ("Delta(6 n^2), n in {2, 3, 4, 6, 8}: Klein residual group for Majorana neutrinos, cyclic residual "
                        "group with three distinct eigenvalues for charged leptons; data ranges rounded (z43)")
        out["universes"]["flavour_groups"] = res

    sl = seam_layers()
    res = search(sl, ["c_eight", "bosonic_extension_available", "spin10_level_one", "sheet_neutral", "u1_free"],
                 lambda o: o["n"], lambda o: f"n={o['n']} {o['kind']} {o['orient']}")
    res["bound"] = "1 to 32 layers, real or complex, second sheet half equally or reversely oriented"
    out["universes"]["seam_layers"] = res

    fp = flavour_pairs()
    res = search(fp, ["mixes", "no_zero_entry", "trimaximal_column", "cp_odd", "same_group"],
                 lambda o: o["index"], lambda o: str(o["pattern"]))
    res["derived"] = {"patterns": sorted({(o["pair"], str(o["pattern"]), o["jarlskog"]) for o in fp if "mixes" in o["clauses"]})}
    res["bound"] = "ordered pairs of abelian subgroups of S4 with three distinct characters on the standard triplet"
    out["universes"]["flavour_pairs"] = res

    principles = allowed_principles(transfer, state)
    cp = composed_process(transfer, state)
    res = search(cp, list(principles), lambda o: o["size"], lambda o: o["label"])
    longest = max(len(r["minimal_rules"][0]) for r in res["unique_selections"]) if res["unique_selections"] else 0
    full = [r for r in res["unique_selections"] if len(r["minimal_rules"][0]) == longest]
    for r in full:                              # which objects satisfy the rule before the size order decides (z46)
        ties = sorted((o for o in cp if set(r["minimal_rules"][0]) <= o["clauses"]), key=lambda o: o["size"])
        r["tie_set_without_size"] = [o["label"] for o in ties]
    out["composed_process"] = {"principles": principles, "objects": res["objects"], "rules": res["rules"],
                               "verdicts": res["verdicts"], "longest_minimal_rule": longest,
                               "longest_minimal_rule_definition": "maximum of the shortest selecting rule length for each selected label; an exploration statistic, not a physical selection principle",
                               "selected_by_the_longest_rules": full,
                               "unique_selections": res["unique_selections"][:40],
                               "bound": "all kernels, all spin-2 kernels, seams with 8, 10, 16, 32 layers, unordered flavour pairs"}
    if transfer and state:
        # the same search without the one principle that carries a scale convention (z46)
        pool = [k for k in principles if k != "dirac_spin2"]
        res = search(cp, pool, lambda o: o["size"], lambda o: o["label"])
        top = max(len(r["minimal_rules"][0]) for r in res["unique_selections"])
        best = [r for r in res["unique_selections"] if len(r["minimal_rules"][0]) == top]
        for r in best:
            r["tie_sets_without_size"] = {"+".join(m): [o["label"] for o in cp if set(m) <= o["clauses"]]
                                          for m in r["minimal_rules"]}
        out["composed_process"]["without_scale_convention"] = {
            "principles": pool, "verdicts": res["verdicts"], "longest_minimal_rule": top,
            "selected_by_the_longest_rules": best}
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
    if source_dynamics and transfer and state:
        from process_source_dynamics import source_dynamics_audit
        out["source_dynamics"] = source_dynamics_audit()
        out["global_process_unique"] = False
        out["historical_composed_search_scope"] = (
            "The composed table keeps the named H6h26 finite universe plus the H6h27 metric stage. The parent-dynamics audit "
            "includes further positive operators with the same thermofield state. "
            "A minimum-size selection in the old table is not a unique physical process.")
    return out


if __name__ == "__main__":
    argparse.ArgumentParser().parse_args()
    print(json.dumps(run(), indent=1, default=str, ensure_ascii=False))
