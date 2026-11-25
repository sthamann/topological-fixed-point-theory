"""Reconstruction search for the open gaps (A) and (B) of the bird's-eye view (H6 exploration H6h34).

Guideline "Rekonstruktion statt Auswahl" (24 November 2026): no new selection principle. Every clause is a computed
property of its candidate: a splitting type, an orbit under the Moebius symmetries of the four marks, an exact KMS
identity, the signature of a real form, a charge count. The engine is the one of the process rules
(`genesis_rules.search`) with one change of use: every candidate has size 0. A rule is UNIQUE only if exactly one label
satisfies it; every other non-empty rule is a tie and is reported with its survivors. Gaps of
docs/research/VOGELPERSPEKTIVE.md:

  normal_bundles  (B) bundles on the seam CP^1 (variable z = S, marks {1, i, -1, -i}) obtainable from TFPT data:
                  O(1) (x) r_*O for the branched double cover E and the (Z2)^3 cover C (characters = even subsets of the
                  marks, z34), the eigen-halves of the latter under every involution of its deck group, duals, the Euler
                  bundle O(1) (x) H^0(O(1)), Hecke modifications O(-1) -> O(-1)(D) of the anti-invariant summand of the
                  double-cover bundle at mark subsets and at subsets of the poles {0, oo} (the fixed points of the clock
                  z -> iz, the only exceptional orbit of size 2), and the Tomita-partner candidate O(1) (x) Vbar.
  doublet_states  (A) 2x2 states built from the original covariance c(x) = 1/(1 + e^{4 pi x}), x in Z + 1/4 or Z + 3/4
                  (ELLIPTIC-CURRENT-SOURCE eq. 1), Hardy mark k = 4x. Enumerated by m = 1..4, the number of Hardy marks
                  spanned by the unit winding (z = S^m). The seam variable of z44/z48/z61 is z = S (m = 1); z = S^4 is
                  the quotient variable w = z^4 (m = 4). Candidates: N = 1 corners and one-particle restrictions of mode
                  pairs (x, x + m/4), and corners {psi0, psi1} of the Fock state of one quarter sector at marks M0 and
                  M0 + m; controls: trace state and Hardy vacuum.
  spin2_family    minimal: quartics of the pencil fixed by the Klein group V4 of the marks through the H6 marks 1 - z^4
                  (z61 Teil 3), with the invariants I, J and the spin-1 matrix on H^0(O(2)). No flavour reading and no
                  data comparison: flavour is settled in TFPT and is not searched here.
  doublet_identifications, position_spaces  (H6h40, z79) the quarter rotation R of the marks (character i^k on the
                  Hardy mark k) on the seam doublet, the quarter pair, the quarter oscillators and the eta doublet of the
                  spin bridge (exact transport on E: w^2 = 1 - z^4), and on the four-dimensional position spaces built
                  from them; clauses: projective order, spinor, R-equivalence, Lorentz or Euclidean compatibility.
  scenarios       (H6h42, S7) whole-theory scenarios: every combination of spinor index, role of the marks, time,
                  modes, gravity and origin, with clauses computed from the doublet/position universes and from recorded
                  exploration results (EVIDENCE: result file, named check, meaning); pending evidence gives no clause.

Exactness: marks and Moebius maps in Gaussian rationals; doublet states are rational functions of u = e^pi (u is
transcendental by Gelfond's theorem, so a rational identity in u holds at u = e^pi exactly when it holds identically);
quartics in Q(i, sqrt 3). Floats only decide signs of a nondegenerate 4x4 form whose determinant is checked exactly.
`run_reconstruction()` is not wired into `process_rules.run()`; the historical outputs of z43-z62 stay unchanged.
"""
from fractions import Fraction as Fr
from functools import lru_cache
from itertools import combinations, permutations, product
from math import gcd
import json

import numpy as np
import sympy as sp

from genesis_rules import powerset, search

SCOPE = ("finite exhaustive enumeration of reconstruction candidates for the open gaps (A) and (B) with computed clauses "
         "and no size order; a unique minimal rule names what the universe and its clause pool fix, not physics")
U = sp.Symbol("u", positive=True)                 # u = e^pi
Z = sp.Symbol("z")


# ---------------------------------------------------------------- Gaussian rationals and Moebius maps of the marks
def _g(re, im=0):
    return (Fr(re), Fr(im))


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def _neg(a):
    return (-a[0], -a[1])


def _mul(a, b):
    return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])


def _conj(a):
    return None if a is None else (a[0], -a[1])


def _div(a, b):
    n = b[0]*b[0] + b[1]*b[1]
    return _mul(a, (b[0]/n, -b[1]/n))


ZERO, ONE = _g(0), _g(1)
MARKS = (_g(1), _g(0, 1), _g(-1), _g(0, -1))
MARK_NAMES = ("1", "i", "-1", "-i")
POINTS = MARKS + (ZERO, None)                     # the four marks, then the poles 0 and oo (None stands for infinity)
POINT_NAMES = MARK_NAMES + ("0", "oo")


def _apply(m, z):
    """Moebius map m = (a, b, c, d): z -> (a z + b)/(c z + d), with None for infinity."""
    a, b, c, d = m
    if z is None:
        return None if c == ZERO else _div(a, c)
    den = _add(_mul(c, z), d)
    return None if den == ZERO else _div(_add(_mul(a, z), b), den)


def _compose(m, n):
    a, b, c, d = m
    e, f, g, h = n
    return (_add(_mul(a, e), _mul(b, g)), _add(_mul(a, f), _mul(b, h)),
            _add(_mul(c, e), _mul(d, g)), _add(_mul(c, f), _mul(d, h)))


def _normal(m):
    lead = next(x for x in m if x != ZERO)
    return tuple(_div(x, lead) for x in m)


def _three_point(z1, z2, z3):
    """The Moebius map sending z1, z2, z3 to 0, infinity, 1."""
    return (_add(z3, _neg(z2)), _mul(_neg(z1), _add(z3, _neg(z2))), _add(z3, _neg(z1)), _mul(_neg(z2), _add(z3, _neg(z1))))


def _sym(x):
    return sp.Rational(x[0].numerator, x[0].denominator) + sp.I*sp.Rational(x[1].numerator, x[1].denominator)


def _perm_mul(p, q):
    return tuple(p[q[k]] for k in range(len(q)))


@lru_cache(maxsize=None)
def mark_symmetries():
    """All Moebius maps that permute the four H6 marks {1, i, -1, -i}. A Moebius map is fixed by three images, so the 24
    assignments of the first three marks are tried and kept when the fourth mark lands on the fourth image. `perm` is
    the induced permutation of POINTS (marks, then 0 and oo)."""
    src = _three_point(*MARKS[:3])
    out = []
    for img in permutations(range(4)):
        a, b, c, d = _three_point(*(MARKS[k] for k in img[:3]))
        m = _normal(_compose((d, _neg(b), _neg(c), a), src))          # T_img^-1 o T_src
        if _apply(m, MARKS[3]) == MARKS[img[3]]:
            a_, b_, c_, d_ = (_sym(x) for x in m)
            out.append({"matrix": m, "perm": tuple(POINTS.index(_apply(m, p)) for p in POINTS),
                        "name": "z -> " + sp.sstr(sp.cancel((a_*Z + b_)/(c_*Z + d_)))})
    return tuple(sorted(out, key=lambda e: e["perm"]))


def _order(perm):
    k, p, one = 1, perm, tuple(range(len(perm)))
    while p != one:
        p, k = _perm_mul(perm, p), k + 1
    return k


def symmetry_structure():
    """Order, element orders and the Klein subgroup acting regularly (without fixed marks) on the marks."""
    group = mark_symmetries()
    perms = [e["perm"] for e in group]
    orders = sorted(_order(p) for p in perms)
    abelian = all(_perm_mul(p, q) == _perm_mul(q, p) for p in perms for q in perms)
    klein = [e for e in group if _order(e["perm"]) == 1
             or (_order(e["perm"]) == 2 and all(e["perm"][k] != k for k in range(4)))]
    closed = all(_perm_mul(p["perm"], q["perm"]) in [k["perm"] for k in klein] for p in klein for q in klein)
    name = "D4" if (len(perms), orders.count(4), abelian) == (8, 2, False) else f"order {len(perms)}"
    return {"group": name, "order": len(perms), "element_orders": orders, "abelian": abelian,
            "elements": [e["name"] for e in group], "regular_klein_group": [e["name"] for e in klein],
            "regular_klein_is_group": closed and len(klein) == 4}


def _same_point(p, q):
    if p is sp.oo or q is sp.oo:
        return p is q
    return sp.simplify(p - q) == 0


def exceptional_orbits():
    """Points of the seam with a nontrivial stabiliser in the mark symmetry group (fixed points of its non-identity
    elements, exact), grouped into orbits."""
    mats = [tuple(_sym(x) for x in e["matrix"]) for e in mark_symmetries()]

    def mob(m, z):
        a, b, c, d = m
        if z is sp.oo:
            return sp.oo if c == 0 else sp.radsimp(a/c)
        den = sp.expand(c*z + d)
        return sp.oo if sp.simplify(den) == 0 else sp.radsimp(sp.expand_complex((a*z + b)/den))
    pts = []
    for a, b, c, d in mats:
        if (a, b, c, d) == (1, 0, 0, 1):
            continue
        fixed = ([sp.oo] + ([sp.radsimp(b/(d - a))] if d != a else [])) if c == 0 else \
            [sp.radsimp(sp.expand_complex(r)) for r in sp.solve(c*Z**2 + (d - a)*Z - b, Z)]
        for p in fixed:
            if not any(_same_point(p, q) for q in pts):
                pts.append(p)
    orbits = []
    for p in pts:
        if any(_same_point(p, q) for o in orbits for q in o):
            continue
        orb = []
        for m in mats:
            q = mob(m, p)
            if not any(_same_point(q, r) for r in orb):
                orb.append(q)
        orbits.append(orb)
    return sorted((sorted(sp.sstr(p) for p in o) for o in orbits), key=lambda o: (len(o), o))


def _orbit(subset):
    return sorted({tuple(sorted(e["perm"][k] for k in subset)) for e in mark_symmetries()})


def _class_orbit(subset):
    """Orbit of a deck involution (mark subset S modulo complement) under the mark symmetries."""
    def canon(s):
        s = tuple(sorted(s))
        t = tuple(sorted(set(range(4)) - set(s)))
        return min(s, t, key=lambda x: (len(x), x))
    return sorted({canon([e["perm"][k] for k in subset]) for e in mark_symmetries()})


# ---------------------------------------------------------------- (B) normal bundles on the seam
CHARACTERS = tuple(frozenset(t) for k in (0, 2, 4) for t in combinations(range(4), k))   # even subsets T of the marks
DOUBLE_COVER = (frozenset(), frozenset(range(4)))                                       # E -> P1, branched at all marks
HECKE_DIVISORS = tuple(D for k in range(1, 5) for D in combinations(range(4), k)) + ((4,), (5,), (4, 5))


def _h0(degrees):
    return sum(max(a + 1, 0) for a in degrees)


def _h1(degrees):
    return sum(max(-a - 1, 0) for a in degrees)


def pushforward(chars, twist=1):
    """Splitting type of O(twist) (x) r_*O for the abelian (Z2)-cover with character group `chars`: the character with
    branch set T contributes L_T^-1 = O(-|T|/2) (L_T^2 = O(T))."""
    return sorted((twist - len(t)//2 for t in chars), reverse=True)


def cover_genus(chars):
    """Riemann-Hurwitz for the cover with character group `chars` (monodromy of order 2 at every ramified mark)."""
    n = len(chars)
    ramified = sum(1 for p in range(4) if any(p in t for t in chars))
    return Fr(-2*n + ramified*n//2 + 2, 2)


def _seam_involutions():
    """Antiholomorphic involutions z -> g(conj z) of the seam with g a mark symmetry, their lift J on V (sigma[v] = [J vbar])
    and type: J Jbar = +1 (real, fixed circle) or -1 (quaternionic, no fixed point)."""
    out = []
    for e in mark_symmetries():
        a, b, c, d = e["matrix"]
        sq = _normal(_compose(e["matrix"], tuple(_conj(x) for x in e["matrix"])))
        if sq != (ONE, ZERO, ZERO, ONE):
            continue
        J = sp.Matrix([[_sym(d), _sym(c)], [_sym(b), _sym(a)]])          # [1 : z] -> [c z + d : a z + b]
        lam = sp.expand((J*J.conjugate())[0, 0])
        assert sp.expand(J*J.conjugate() - lam*sp.eye(2)) == sp.zeros(2) and lam in (1, -1)
        image = lambda p, m=e["matrix"]: POINTS.index(_apply(m, _conj(POINTS[p])))
        out.append({"name": "z -> " + e["name"].split("-> ")[1].replace("z", "conj(z)"), "J": J, "type": 1 if lam > 0 else -1,
                    "perm": tuple(image(p) for p in range(len(POINTS)))})
    return out


EPS = sp.Matrix([[0, -1], [1, 0]])


@lru_cache(maxsize=None)
def _signature(kind, J_key=None, K_key=None):
    """Signature of the resultant form det M on the real form {M : R(M) = M} of H^0(O(1) (x) C^2) = C^{2x2}
    (columns: the two fibre components as linear forms). Seam structure: R(M) = Jbar^T conj(M) K^T. Hermitian pairing of a
    Vbar fibre: R(M) = conj(M)^T, or its epsilon-twisted form eps conj(M)^T eps^T."""
    if kind == "seam":
        J, K = sp.Matrix(2, 2, list(J_key)), sp.Matrix(2, 2, list(K_key))
        R = lambda M: (J.conjugate().T*M.conjugate()*K.T).applyfunc(sp.expand)
    elif kind == "hermitian":
        R = lambda M: M.conjugate().T
    else:
        R = lambda M: (EPS*M.conjugate().T*EPS.T).applyfunc(sp.expand)
    basis = []
    for pos in ((0, 0), (0, 1), (1, 0), (1, 1)):
        for f in (1, sp.I):
            B = sp.zeros(2)
            B[pos] = f
            basis.append(B)
    if any(R(R(B)) != B for B in basis):
        return None                                                  # not a real structure
    chosen, rows = [], []
    for B in basis:
        F = (B + R(B)).applyfunc(sp.expand)
        trial = rows + [[c for x in F for c in (sp.re(x), sp.im(x))]]
        if sp.Matrix(trial).rank() == len(trial):
            chosen, rows = chosen + [F], trial
    assert len(chosen) == 4
    det = lambda M: sp.expand(M.det())
    G = [[sp.expand((det(A + B) - det(A) - det(B))/2) for B in chosen] for A in chosen]
    ref = sp.conjugate(next(x for row in G for x in row if x != 0))
    G = sp.Matrix(G).applyfunc(lambda x: sp.expand(x*ref))             # all entries lie on one line through 0
    assert all(sp.im(x) == 0 for x in G) and G.det() != 0
    ev = np.linalg.eigvalsh(np.array(G.applyfunc(float).tolist(), dtype=float))
    return (int(np.sum(ev > 0)), int(np.sum(ev < 0)))


def reality_signatures(cand):
    """Real forms of H^0(N) for a twistor-line candidate N = O(1) (x) C^2 and the signature of the null cone on each.
    Seam structures: the mark-preserving antiholomorphic involutions that also preserve the modification divisor, lifted
    with a fibre structure K of the same type (K = 1 for real, K = eps for quaternionic; for the Euler bundle K is induced
    on H^0(O(1)) itself). A Vbar fibre also carries the Hermitian pairing V* (x) Vbar. A Hermitian form on V (a faithful
    doublet state) would identify H^0(O(1)) = V* with Vbar and turn the Euler bundle into the Tomita candidate; that
    identification is gap (A) and is not assumed here."""
    out = []
    D = set(cand["divisor"])
    for s in _seam_involutions():
        if {s["perm"][p] for p in D} != D:
            continue
        if cand["fibre"] == "H0(O(1))":
            K = s["J"].conjugate().T                                 # sigma on sections: s -> Jbar^T conj(s)
        elif cand["fibre"] == "Vbar":
            K = s["J"].conjugate()                                   # on Vbar, w = conj(v): w -> Jbar conj(w)
        else:
            K = sp.eye(2) if s["type"] == 1 else EPS
        sig = _signature("seam", tuple(s["J"]), tuple(K))
        if sig is not None:
            out.append({"structure": s["name"], "type": "real" if s["type"] == 1 else "quaternionic", "signature": list(sig)})
    if cand["fibre"] == "Vbar":
        for kind, name in (("hermitian", "Hermitian pairing M -> M^+"), ("hermitian_eps", "Hermitian pairing eps M^+ eps^T")):
            out.append({"structure": name, "type": "hermitian", "signature": list(_signature(kind))})
    return out


def normal_bundle_candidates():
    """Bundles on the seam obtainable from TFPT data, with computed clauses:
      twistor_line    splitting type O(1) + O(1) (normal bundle of a twistor line);
      mark_symmetric  the modification divisor (empty if none) is invariant under every Moebius symmetry of the marks;
      choice_free     the whole construction datum (divisor, deck involution) is alone in its orbit under those symmetries;
      lorentz_cone    twistor line, and one of its reality structures makes the null cone of H^0(N) Lorentzian (1,3)."""
    cands = []

    def add(name, degrees, kind, datum=(), fibre=None, divisor=()):
        cands.append({"name": name, "degrees": list(degrees), "kind": kind, "datum": tuple(datum), "fibre": fibre,
                      "divisor": tuple(divisor), "cover_derived": kind in ("cover", "half", "hecke", "dual")})

    add("O(1) x r_*O_E (double cover)", pushforward(DOUBLE_COVER), "cover")
    add("O(1) x r_*O_C ((Z2)^3 cover)", pushforward(CHARACTERS), "cover")
    halves = []
    for k in (1, 2):
        for S in combinations(range(4), k):
            if k == 2 and 0 not in S:
                continue                                             # a pair and its complement are one involution
            kind = "point involution at mark " + MARK_NAMES[S[0]] if k == 1 else \
                "half period {" + ",".join(MARK_NAMES[p] for p in S) + "}|{" + \
                ",".join(MARK_NAMES[p] for p in range(4) if p not in S) + "}"
            for sign, label in ((0, "E+"), (1, "E-")):
                chars = [t for t in CHARACTERS if len(t & set(S)) % 2 == sign]
                halves.append((f"{label} under the {kind}", pushforward(chars), S))
    for name, deg, S in halves:
        add(name, deg, "half", datum=S)
    for c in list(cands):
        add("dual of " + c["name"], sorted((-a for a in c["degrees"]), reverse=True), "dual", datum=c["datum"])
    add("O(1) x H^0(O(1)) (Euler bundle)", [1 + 0]*_h0([1]), "doublet", fibre="H0(O(1))")
    add("O(1) x Vbar (Tomita partner)", [1, 1], "doublet", fibre="Vbar")
    add("dual of O(1) x H^0(O(1))", [-1, -1], "doublet")
    add("dual of O(1) x Vbar", [-1, -1], "doublet")
    for D in HECKE_DIVISORS:
        deg = sorted([1, -1 + len(D)], reverse=True)                   # summand O(-1) of O(1) x r_*O_E replaced by O(-1)(D)
        add("Hecke O(-1) -> O(-1)(D), D = {" + ",".join(POINT_NAMES[p] for p in D) + "}", deg, "hecke", datum=D, divisor=D)
    for c in cands:
        orbit = _class_orbit(c["datum"]) if c["kind"] in ("half", "dual") and c["datum"] else _orbit(c["datum"])
        c["orbit_size"] = len(orbit)
        c["rank"], c["h0"], c["h1"] = len(c["degrees"]), _h0(c["degrees"]), _h1(c["degrees"])
        clauses = set()
        if c["degrees"] == [1, 1]:
            clauses.add("twistor_line")
        if len(_orbit(c["divisor"])) == 1:
            clauses.add("mark_symmetric")
        if c["orbit_size"] == 1:
            clauses.add("choice_free")
        c["reality"] = reality_signatures(c) if "twistor_line" in clauses else []
        if any(sorted(r["signature"]) == [1, 3] for r in c["reality"]):
            clauses.add("lorentz_cone")
        c["clauses"] = clauses
    return cands


# ---------------------------------------------------------------- (A) doublet states from the original covariance
QUARTERS = (Fr(1, 4), Fr(3, 4))
PAIR_BASES = (-2, -1, 0, 1)        # base points x = alpha + r of mode pairs
M_VALUES = (1, 2, 3, 4)            # Hardy marks spanned by the unit winding, z = S^m (m = 1: seam variable z = S)
FOCK_BASES = (0, 1, 2, 3, 4)       # base marks M0 of Fock corners
MARK_MAX = 8


def cov(x):
    """Original quarter covariance c(x) = (1 + e^{4 pi x})^-1 with u = e^pi; the Hardy mark k = 4x is odd."""
    k = 4*Fr(x)
    assert k.denominator == 1
    return 1/(1 + U**int(k))


def cov_hardy(x):
    """Hardy (Szegő) vacuum: the zero-temperature limit, every mode below the Fermi level x = 0 filled."""
    return sp.Integer(1) if x < 0 else sp.Integer(0)


def _q(x):
    x = Fr(x)
    return sp.Rational(x.numerator, x.denominator)


def original_symmetries(window=range(-3, 4)):
    """Which label maps leave the original covariance invariant (exact identities in u)."""
    xs = [r + a for a in QUARTERS for r in window]
    same = lambda f: all(sp.cancel(f(x)) == 0 for x in xs)

    def rotation_commutes(x):                      # quarter rotation z -> iz: character i^k on the Hardy mark k = 4x
        C = sp.diag(cov(x), cov(x + 1), cov(x + Fr(1, 2)))
        R = sp.diag(*(sp.I**int(4*y) for y in (x, x + 1, x + Fr(1, 2))))
        return (R*C*R.conjugate().T - C).applyfunc(sp.simplify) == sp.zeros(3)
    return {"mirror x -> -x with particle-hole": same(lambda x: cov(-x) - (1 - cov(x))),
            "reflection x -> -x": same(lambda x: cov(-x) - cov(x)),
            "shift x -> x + 1": same(lambda x: cov(x + 1) - cov(x)),
            "sector exchange x -> x + 1/2": same(lambda x: cov(x + Fr(1, 2)) - cov(x)),
            "quarter rotation z -> iz": all(rotation_commutes(x) for x in xs)}


def _state(kind, pair, c):
    """2x2 density matrix in the basis ordered by the mode label (the order of the generator 4 pi x = pi k)."""
    x1, x2 = pair
    if kind == "covariance":                                      # one-particle restriction, normalised by <N>
        w = (c(x1), c(x2))
    else:                                                          # N = 1 sector of the two-mode state
        w = (c(x1)*(1 - c(x2)), (1 - c(x1))*c(x2))
    tot = sp.cancel(w[0] + w[1])
    if tot == 0:
        return None
    return sp.diag(sp.cancel(w[0]/tot), sp.cancel(w[1]/tot))


def _positive(expr):
    """Exact sign at u = e^pi > 1. Numerator and denominator are nonzero there (u is transcendental); if neither has a root
    in (1, oo) (Sturm count), the sign on (1, oo) is the sign at u = 2, computed in rationals."""
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    if num == 0:
        return False
    for p in (num, den):
        if sp.Poly(p, U).count_roots(1, sp.oo) != 0:
            raise ValueError(f"sign of {expr} on (1, oo) is not constant")
    return bool(sp.sign(num.subs(U, 2))*sp.sign(den.subs(U, 2)) > 0)


def _mode_ratio(x):
    """Weight of flipping mode x against the Fermi sea, from the covariance: particle at x > 0, hole at x < 0."""
    return sp.cancel(cov(x)/(1 - cov(x)) if x > 0 else (1 - cov(x))/cov(x))


def fock_patterns(alpha, mark_max=MARK_MAX):
    """Occupation patterns of one quarter sector relative to its Fermi sea Omega up to the mark M = 4 Delta L0 = sum 4|x|
    of the flipped modes; Q = particles - holes, parity = number of flips mod 2."""
    modes = [r + alpha for r in range(-3, 3) if 4*abs(r + alpha) <= mark_max]
    out = []
    for k in range(len(modes) + 1):
        for flips in combinations(modes, k):
            M = sum(int(4*abs(x)) for x in flips)
            if M <= mark_max:
                out.append({"flips": tuple(flips), "M": M, "Q": sum(1 if x > 0 else -1 for x in flips), "parity": k % 2})
    return sorted(out, key=lambda s: (s["M"], s["Q"], s["flips"]))


def _partitions(n):
    table = [1] + [0]*n
    for part in range(1, n + 1):
        for s in range(part, n + 1):
            table[s] += table[s - part]
    return table[n]


def fock_multiplicities(alpha, mark_max=MARK_MAX):
    """Level multiplicities counted from the patterns and from bosonization: M = 4(N + F(Q)), F(Q) = Q^2/2 +
    (alpha - 1/2) Q (ELLIPTIC-CURRENT-SOURCE eq. 2), p(N) oscillator states."""
    pats = fock_patterns(alpha, mark_max)
    counted = [sum(1 for s in pats if s["M"] == M) for M in range(mark_max + 1)]
    formula = []
    for M in range(mark_max + 1):
        tot = 0
        for Q in range(-4, 5):
            N = Fr(M, 4) - (Fr(Q*Q, 2) + (alpha - Fr(1, 2))*Q)
            if N >= 0 and N.denominator == 1:
                tot += _partitions(int(N))
        formula.append(tot)
    return counted, formula


def _weight(flips):
    return sp.cancel(sp.Mul(*[_mode_ratio(x) for x in flips]))


def _corner(w0, w1):
    tot = sp.cancel(w0 + w1)
    return None if tot == 0 else sp.diag(sp.cancel(w0/tot), sp.cancel(w1/tot))


def _pattern_name(s):
    if not s["flips"]:
        return "Omega (M=0, Q=0)"
    return " + ".join(("particle " if x > 0 else "hole ") + str(x) for x in s["flips"]) + f" (M={s['M']}, Q={s['Q']})"


def _r_order(k0, k1):
    """Projective order of the quarter rotation R = diag(i^k0, i^k1) on P(V) (character i^k on the Hardy mark k)."""
    return next(n for n in range(1, 5) if sp.I**(n*(k1 - k0)) == 1)


def _doublet_families():
    """(construction, m, items); an item carries its state, the state the same construction gives on the mirror image,
    the Hardy marks (k0, k1) of its two basis vectors, and its charge bookkeeping."""
    fams = []
    for m in M_VALUES:
        pairs = [(r + a, r + a + Fr(m, 4)) for a in QUARTERS for r in PAIR_BASES]
        pairs = [p for p in pairs if (4*p[1]).denominator == 1 and int(4*p[1]) % 2 == 1]   # both labels are modes
        for kind, name, c in (("covariance", "one-particle restriction", cov), ("n1", "N=1 corner of a mode pair", cov),
                              ("n1", "N=1 corner of a mode pair, Hardy vacuum (control)", cov_hardy)):
            items = []
            for pair in pairs:
                rho = _state(kind, pair, c)
                if rho is None:
                    continue
                image = tuple(sorted(-x for x in pair))
                items.append({"base": [str(x) for x in pair], "rho": rho, "image_rho": _state(kind, image, c),
                              "marks": (int(4*pair[0]), int(4*pair[1])), "mirror_closed": image == pair,
                              "neutral": (pair[1] - pair[0]).denominator == 1})
            if items:
                fams.append((f"{name}, m = {m}", m, items))
    pats = {a: fock_patterns(a) for a in QUARTERS}
    for m in M_VALUES:
        items = []
        for a in QUARTERS:
            level = {}
            for s in pats[a]:
                level.setdefault(s["M"], []).append(s)
            for M0 in FOCK_BASES:
                for s0 in level.get(M0, []):
                    for s1 in level.get(M0 + m, []):
                        mirror = [tuple(sorted(-x for x in s["flips"])) for s in (s0, s1)]
                        items.append({"base": [f"alpha = {a}", _pattern_name(s0), _pattern_name(s1)],
                                      "rho": _corner(_weight(s0["flips"]), _weight(s1["flips"])),
                                      "image_rho": _corner(_weight(mirror[0]), _weight(mirror[1])),
                                      "marks": (M0, M0 + m), "mirror_closed": False,
                                      "neutral": s0["Q"] == s1["Q"] and s0["parity"] == s1["parity"]})
        fams.append((f"Fock corner, m = {m}", m, items))
    return fams


def _ratio_label(rho):
    if rho[0, 0] == 0:
        return "diag(0, 1)"
    r = sp.cancel(rho[1, 1]/rho[0, 0])
    return "diag(1, " + sp.sstr(r.subs(U, sp.exp(sp.pi))) + ")/trace"


def doublet_state_candidates():
    """Doublet states with computed clauses (Hardy marks k0 < k1 of the two basis vectors, m = k1 - k0):
      faithful                both eigenvalues > 0 (exact);
      kms                     rho = e^{-h}/Tr for the restricted modular generator h = pi diag(k0, k1) (4 pi x = pi k):
                              rho11/rho00 = u^-(k1 - k0);
      symmetric               the same construction on the mirror image (x -> -x with particle-hole conjugation, which
                              maps a quarter sector to the other) gives the same matrix, and the quarter rotation
                              R = diag(i^k0, i^k1) leaves it invariant;
      base_point_independent  every admissible base point of the construction (mode pair, sector, base mark, level state)
                              gives this matrix;
      mirror_closed           the mirror maps the doublet to itself (a mode and its own mirror partner);
      marks_separated         R has projective order 4 on P(V): the four marks are four points of the doublet's sphere;
      neutral                 the matrix units of End(V) keep the charge of each quarter sector and the fermion parity
                              (they lie among the original's current responses)."""
    out = []
    for name, m, items in _doublet_families():
        same = all((it["rho"] - items[0]["rho"]).applyfunc(sp.cancel) == sp.zeros(2) for it in items)
        for it in items:
            rho, (k0, k1) = it["rho"], it["marks"]
            R = sp.diag(sp.I**k0, sp.I**k1)
            clauses = set()
            if _positive(rho[0, 0]) and _positive(rho[1, 1]):
                clauses.add("faithful")
            if sp.cancel(rho[1, 1] - U**(k0 - k1)*rho[0, 0]) == 0:
                clauses.add("kms")
            if it["image_rho"] is not None and (it["image_rho"] - rho).applyfunc(sp.cancel) == sp.zeros(2) \
                    and (R*rho*R.conjugate().T - rho).applyfunc(sp.simplify) == sp.zeros(2):
                clauses.add("symmetric")
            if same:
                clauses.add("base_point_independent")
            if it["mirror_closed"]:
                clauses.add("mirror_closed")
            if _r_order(k0, k1) == 4:
                clauses.add("marks_separated")
            if it["neutral"]:
                clauses.add("neutral")
            out.append({"construction": name, "m": m, "base": it["base"], "label": _ratio_label(rho), "clauses": clauses,
                        "r_order": _r_order(k0, k1),
                        "rho": [sp.sstr(rho[0, 0].subs(U, sp.exp(sp.pi))), sp.sstr(rho[1, 1].subs(U, sp.exp(sp.pi)))],
                        "rho_numeric": [float(rho[0, 0].subs(U, sp.exp(sp.pi))), float(rho[1, 1].subs(U, sp.exp(sp.pi)))]})
    half = sp.Rational(1, 2)
    out.append({"construction": "trace state (control)", "m": None, "base": [], "label": "diag(1, 1)/trace", "r_order": None,
                "clauses": {"faithful", "symmetric", "base_point_independent"} if _positive(half) else set(),
                "rho": ["1/2", "1/2"], "rho_numeric": [0.5, 0.5]})
    return out


# ---------------------------------------------------------------- spin-2 pencil (minimal, no flavour reading)
def _coeffs(psi):
    p = sp.expand(psi)
    return [sp.expand(p.coeff(Z, k)) for k in range(5)]


def np_components(psi):
    """Psi(z) = Psi0 + 4 Psi1 z + 6 Psi2 z^2 + 4 Psi3 z^3 + Psi4 z^4 (pi = (1, z), Newman-Penrose labels as in z61)."""
    c = _coeffs(psi)
    return (c[0], c[1]/4, c[2]/6, c[3]/4, c[4])


def invariants(psi):
    P0, P1, P2, P3, P4 = np_components(psi)
    I_ = sp.expand(P0*P4 - 4*P1*P3 + 3*P2**2)
    J_ = sp.expand(sp.Matrix([[P4, P3, P2], [P3, P2, P1], [P2, P1, P0]]).det())
    return I_, J_, sp.expand(I_**3 - 27*J_**2)


def spin1_matrix(psi):
    """Psi_AB^CD on symmetric 2-spinors phi_CD (basis phi00, phi01, phi11; phi^CD = eps^CE eps^DF phi_EF), i.e. on
    H^0(O(2)). Trace 0, tr M^2 = 2 I, det M = 2 J."""
    P0, P1, P2, P3, P4 = np_components(psi)
    return sp.Matrix([[P2, -2*P1, P0], [P3, -2*P2, P1], [P4, -2*P3, P2]])


def _pullback(m, psi):
    a, b, c, d = (_sym(x) for x in m)
    return sp.expand(sp.cancel((c*Z + d)**4*psi.subs(Z, (a*Z + b)/(c*Z + d))))


def _proportional(p, q):
    return sp.Matrix([_coeffs(p), _coeffs(q)]).rank() <= 1


def v4_pencil():
    """Joint eigenspace of the regular Klein group V4 of the marks on H^0(O(4)) that contains the H6 quartic 1 - z^4."""
    klein = [e for e in mark_symmetries() if _order(e["perm"]) == 2 and all(e["perm"][k] != k for k in range(4))]
    h6 = 1 - Z**4
    rows = []
    for e in klein:
        lam = sp.cancel(_pullback(e["matrix"], h6)/h6)
        assert lam.is_number
        A = sp.Matrix([_coeffs(_pullback(e["matrix"], Z**k)) for k in range(5)]).T
        rows.append(A - lam*sp.eye(5))
    null = sp.Matrix.vstack(*rows).nullspace()
    span = [sp.expand(sum(v[k]*Z**k for k in range(5))) for v in null]
    assert len(span) == 2 and sp.Matrix([_coeffs(p) for p in span] + [_coeffs(1 - Z**4), _coeffs(Z**2)]).rank() == 2
    return span


def _petrov(psi):
    poly = sp.Poly(psi, Z)
    mult = sorted(list(sp.roots(poly).values()) + ([4 - poly.degree()] if poly.degree() < 4 else []), reverse=True)
    assert sum(mult) == 4
    return {(1, 1, 1, 1): "I", (2, 1, 1): "II", (2, 2): "D", (3, 1): "III", (4,): "N"}[tuple(mult)]


def spin2_family_candidates():
    """Quartics Psi_t = 1 - z^4 + t z^2 of the V4 pencil and Psi_inf = z^2: every parameter where I, J or the discriminant
    vanishes along the pencil, plus t = +-1 as generic controls. Computed clauses: j_zero, i_zero, zero_eigenvalue (of the
    spin-1 matrix), nondegenerate (four distinct roots), d4_fixed (fixed up to a scalar by every mark symmetry)."""
    t = sp.Symbol("t")
    span = v4_pencil()
    psi_t = 1 - Z**4 + t*Z**2
    special = set()
    for f in invariants(psi_t):
        special |= set(sp.solve(f, t))
    params = sorted(special | {sp.Integer(1), sp.Integer(-1)},
                    key=lambda v: (float(abs(complex(v))), float(complex(v).real), float(complex(v).imag)))
    members = [(sp.sstr(v), sp.expand(psi_t.subs(t, v))) for v in params] + [("oo", Z**2)]
    out = []
    for label, psi in members:
        I_, J_, D_ = invariants(psi)
        M = spin1_matrix(psi)
        ev = M.eigenvals()
        evs = sorted(ev, key=lambda v: (float(sp.re(v)), float(sp.im(v))))
        clauses = set()
        if J_ == 0:
            clauses.add("j_zero")
        if I_ == 0:
            clauses.add("i_zero")
        if any(sp.simplify(v) == 0 for v in ev):
            clauses.add("zero_eigenvalue")
        if sp.simplify(D_) != 0:
            clauses.add("nondegenerate")
        if all(_proportional(psi, _pullback(e["matrix"], psi)) for e in mark_symmetries()):
            clauses.add("d4_fixed")
        out.append({"t": label, "quartic": sp.sstr(psi), "label": "Psi = " + sp.sstr(psi), "clauses": clauses,
                    "I": sp.sstr(I_), "J": sp.sstr(J_), "discriminant": sp.sstr(sp.factor(D_)),
                    "eigenvalues": [sp.sstr(sp.simplify(v)) for v in evs for _ in range(ev[v])],
                    "kernel_phi00_phi01_phi11": [[sp.sstr(sp.simplify(x)) for x in v] for v in M.nullspace()],
                    "petrov": _petrov(psi),
                    "trace_identities": bool(sp.expand(M.trace()) == 0 and sp.expand((M*M).trace() - 2*I_) == 0
                                             and sp.expand(M.det() - 2*J_) == 0)})
    return out, [sp.sstr(p) for p in span]


# ---------------------------------------------------------------- search over each universe
# ---------------------------------------------------------------- doublet identifications (H6 exploration H6h40, z79)
F_ETA = (sp.Integer(0), 1/(Z**2 - 1))          # F = w/(z^2 - 1) on E: w^2 = 1 - z^4, div F = R_i + R_-i - R_1 - R_-1


def _hmul(f, g):
    """Product of functions a + b w on E (w^2 = 1 - z^4), stored as pairs (a, b) of rational functions of z."""
    return (sp.cancel(f[0]*g[0] + f[1]*g[1]*(1 - Z**4)), sp.cancel(f[0]*g[1] + f[1]*g[0]))


def _hpull(f, sign):
    """f o g^-1 for the quarter rotation g: (z, w) -> (i z, sign w); sign = 0 is the deck (z, w) -> (z, -w)."""
    if sign == 0:
        return (f[0], -f[1])
    return (sp.cancel(f[0].subs(Z, -sp.I*Z)), sp.cancel(sign*f[1].subs(Z, -sp.I*Z)))


def _hcoords(f, d):
    """Coordinates of f in the basis 1, .., z^d, F, .., z^d F of L(D); None if f is not in L(D)."""
    a, b = sp.cancel(f[0]), sp.cancel(f[1]*(Z**2 - 1))
    if not (a.is_polynomial(Z) and b.is_polynomial(Z)) or sp.degree(a, Z) > d or sp.degree(b, Z) > d:
        return None
    pa, pb = sp.Poly(a, Z), sp.Poly(b, Z)
    return [pa.coeff_monomial(Z**k) for k in range(d + 1)] + [pb.coeff_monomial(Z**k) for k in range(d + 1)]


def _transport(d, sign):
    """Quarter rotation (sign = +-1, lift w -> sign w) or deck (sign = 0) on L(D): d = 0 gives the eta doublet
    H^0(E, eta (x) pi^*O(1)) = L(R_1 + R_-1), d = 1 the sections H^0(N) = L(R_1 + R_-1 + oo_+ + oo_-) of the spin
    bridge N = pi_*(eta (x) pi^*O(2)) (z70). A section s goes to F (s o g^-1), since g(D) - D = div F."""
    basis = [(Z**k, sp.Integer(0)) for k in range(d + 1)] + [(sp.Integer(0), Z**k/(Z**2 - 1)) for k in range(d + 1)]
    cols = []
    for s in basis:
        img = _hpull(s, 0) if sign == 0 else _hmul(F_ETA, _hpull(s, sign))
        c = _hcoords(img, d)
        if c is None:
            raise ValueError("the transport leaves L(D)")
        cols.append(c)
    return sp.Matrix(cols).T


def _kron(A, B):
    return sp.Matrix(A.rows*B.rows, A.cols*B.cols, lambda i, j: A[i // B.rows, j // B.cols]*B[i % B.rows, j % B.cols])


def _proj_order(M):
    P = sp.eye(M.rows)
    for k in range(1, 17):
        P = (P*M).applyfunc(sp.expand)
        if (P - P[0, 0]*sp.eye(M.rows)).is_zero_matrix:
            return k
    return None


def _turns(M):
    """Eigenvalues of M as exponents k of e^(2 pi i k/16); M^16 = I is checked exactly, floats only read off k."""
    if not ((M**16).applyfunc(sp.expand) - sp.eye(M.rows)).is_zero_matrix:
        raise ValueError("not of order dividing 16")
    ks = np.angle(np.linalg.eigvals(np.array(M.evalf(), dtype=complex)))*16/(2*np.pi)
    if np.abs(ks - np.round(ks)).max() > 1e-9:
        raise ValueError("eigenvalue is not a 16th root of unity")
    return sorted(int(round(k)) % 16 for k in ks)


def _same_up_to_phase(s, t):
    return any(sorted((k + p) % 16 for k in s) == sorted(t) for p in range(16))


def _lorentz_compatible(t):
    """Spectrum (up to phase) of g (x) conj(g) for g = diag(l, 1/l): {1, 1, l^2, l^-2}."""
    return any(_same_up_to_phase(t, [0, 0, (2*l) % 16, (-2*l) % 16]) for l in range(16))


def _euclidean_compatible(t):
    """Spectrum (up to phase) of g (x) h for g = diag(l, 1/l), h = diag(m, 1/m) in SU(2) x SU(2)."""
    return any(_same_up_to_phase(t, [(a*l + b*m) % 16 for a in (1, -1) for b in (1, -1)])
               for l in range(16) for m in range(16))


SEAM_DOUBLET = "seam doublet V_1 (Fock corner {Omega, psi*_1/4 Omega}, rho_1)"
QUARTER_PAIR = "quarter pair (one-particle modes x = +1/4, -1/4)"
QUARTER_OSC = "quarter oscillators (a_1, a_2), every mode m"
ETA_SPINOR = "eta doublet H^0(E, eta (x) pi^*O(1)), lift w -> -w"
SPIN_BRIDGE = "spin bridge H^0(N) = L(R_1 + R_-1 + oo_+ + oo_-)"


def doublet_identification_candidates():
    """Doublets on which the quarter rotation R of the marks acts (character i^k on the Hardy mark k), and the
    four-dimensional position spaces built from them. Returns (doublets, position spaces)."""
    hardy = lambda k0, k1: sp.diag(sp.I**k0, sp.I**k1)
    eta = {"-": _transport(0, -1), "+": _transport(0, 1), "deck": _transport(0, 0)}
    vstar = sp.diag(1, -sp.I)                    # H^0(O(1)) = span{1, z}: s -> s o g^-1, z -> -i z
    doublets = [(SEAM_DOUBLET, hardy(0, 1)), (QUARTER_PAIR, hardy(1, -1)), (QUARTER_OSC, hardy(4, 4)),
                (ETA_SPINOR, eta["-"]), ("eta doublet, lift w -> w", eta["+"]), ("eta doublet, deck", eta["deck"])]
    seam, pair = _turns(hardy(0, 1)), _turns(hardy(1, -1))
    ds = []
    for name, M in doublets:
        t = _turns(M)
        cl = {f"projective_order_{_proj_order(M)}"}
        if sp.expand(M.det()) == 1 and (M*M + sp.eye(2)).applyfunc(sp.expand).is_zero_matrix:
            cl.add("spinor")
        if _same_up_to_phase(t, seam):
            cl.add("r_equivalent_to_seam")
        if _same_up_to_phase(t, pair):
            cl.add("r_equivalent_to_quarter_pair")
        ds.append({"label": name, "matrix": str(M.tolist()), "turns": t, "clauses": cl})
    N = _transport(1, -1)
    spaces = [("Tomita reading T': Hom(V, Vbar^*) = V^* (x) Vbar^*", _kron(vstar, vstar.conjugate()), False),
              (SPIN_BRIDGE, N, True),
              ("relative states of z78: Herm of the quarter oscillators", _kron(hardy(4, 4), hardy(4, 4).conjugate()), True),
              ("Herm of the quarter pair", _kron(hardy(1, -1), hardy(1, -1).conjugate()), True)]
    ps = []
    for name, M, native in spaces:
        t = _turns(M)
        cl = {c for c, ok in (("native_r", native), ("lorentz_compatible", _lorentz_compatible(t)),
                              ("euclidean_compatible", _euclidean_compatible(t)),
                              ("marked", any((a - b) % 16 in (4, 12) for a in t for b in t))) if ok}
        ps.append({"label": name, "matrix": str(M.tolist()), "turns": t, "clauses": cl,
                   "factorises_as_vstar_x_eta": _same_up_to_phase(t, _turns(_kron(vstar, eta["-"])))})
    return ds, ps


# ---------------------------------------------------------------- whole-theory scenarios (H6 exploration H6h42, S7)
from pathlib import Path as _Path

EXPLORATION = _Path(__file__).resolve().parents[2] / "experiments/exploration/2026-11-23-phasenraum-f5-holocode"
SCENARIO_AXES = {
    "spinor_index": ["quarter oscillators", "step type", "eta doublet", "seam doublet V_1"],
    "marks": ["internal spin-Z4", "Euclidean Weyl curvature", "celestial sphere"],
    "time": ["modular flow", "parabolic Wiesbrock", "OS reconstruction", "radial conformal"],
    "modes": ["species", "energy levels"],
    "gravity": ["composite spin 2", "Virasoro T(z)", "twistor deformation: native sources",
                "twistor deformation: radial channel", "entanglement: original|partner cut", "entanglement: spatial cut"],
    "origin": ["seed rho_1 + spin-Z4", "two binary observations"],
}


def _check_named(d, fragment):
    hits = [v for k, v in d["checks"].items() if fragment in k]
    if len(hits) != 1:
        raise KeyError(f"check '{fragment}' found {len(hits)} times")
    return bool(hits[0])


# (axis, value, status, result file, meaning, predicate on the result json). A missing file leaves the value pending.
EVIDENCE = [
    ("time", "modular flow", "excluded", "z70_ergebnis.json", "spec(K) in pi Z: periodic modular flow, no time translation",
     lambda d: _check_named(d, "spec(K) in pi Z gibt exp(2 i K) = 1")),
    ("time", "parabolic Wiesbrock", "supported", "z72_ergebnis.json", "parabolic Wiesbrock difference in the vacuum sector, "
     "spectrum [0, oo) (one-dimensional)", lambda d: all(d["checks"].values())),
    ("gravity", "Virasoro T(z)", "excluded", "z75_ergebnis.json", "T(z) maps to Delta = 2 and (j, j), no spin-2 level",
     lambda d: all(d["checks"].values())),
    ("gravity", "composite spin 2", "excluded", "z80_ergebnis.json", "composite spin 2 of a free scalar: continuum, no pole",
     lambda d: _check_named(d, "kein isolierter masseloser Pol")),
    ("gravity", "entanglement: original|partner cut", "excluded", "z86_ergebnis.json", "the cut original|partner is a "
     "thermofield cut of a type I factor with periodic modular flow: no fixed S/area, no horizon without BW (scoped)",
     lambda d: _check_named(d, "Typfrage: spec K in pi N")),
    ("marks", "internal spin-Z4", "conditional", "z82_ergebnis.json", "the clock is the Z4 of the 16 on the carrier modulo "
     "G_SM (spinor lift z_X, mod-16 balance); in the original it reaches only the R^2 = -1 '2', not the '3'",
     lambda d: _check_named(d, "mod-16-Bilanz der Klasse von R") and _check_named(d, "Spinorhebung R~ = det(R)^(-1/2)")),
    ("spinor_index", "quarter oscillators", "excluded", "z82_ergebnis.json", "as spinor index of matter: vertex operators "
     "carry no definite outer 2 pi sign under the quarter-oscillator Lorentz group of z79 (kill criterion 2 of S1, scoped)",
     lambda d: _check_named(d, "Gegenprobe Außenbild (z79)")),
    ("origin", "seed rho_1 + spin-Z4", "conditional", "z87_ergebnis.json", "seed family rho_m = rho_1^m/Tr and Fermi statistics "
     "forced; only the original passes the five criteria; origin reduced to two named premises (P_sym, P_clock); the Gaussian "
     "premise follows from P_clock in all computed truncations",
     lambda d: _check_named(d, "nur das Original erfüllt alle fünf Kriterien") and _check_named(d, "S6-Kriterien")),
    ("modes", "energy levels", "conditional", "z83_ergebnis.json", "su(2,2) closes exactly with the CCR step amplitudes, "
     "the Tomita square roots of the native pair-conserving raiser L^nat; nu = 1 ladder; modes are energy levels E = L_0 = m "
     "(one identification: doubled diagonal step = native single step)",
     lambda d: _check_named(d, "[E] (iii) CCR-Amplituden") and _check_named(d, "[E] K-Typen: auf der Schale m")),
    ("time", "radial conformal", "conditional", "z83_ergebnis.json", "Minkowski energy E + (P0 + K0)/2 of the radial su(2,2) "
     "is parabolic with absolutely continuous spectrum [0, oo) (same identification)",
     lambda d: _check_named(d, "[E] Minkowskienergie P^M_0")),
    ("spinor_index", "step type", "conditional", "z83_ergebnis.json", "spinor index of the radial su(2,2): the two native "
     "step types (particle up, hole down); the clock R is 1 on both (R-neutral spacetime); matter spin open",
     lambda d: _check_named(d, "[E] Uhr: R = e^{2 pi i L_0}")),
    ("gravity", "twistor deformation: native sources", "excluded", "z85_ergebnis.json", "Gaussian data, vertex operators, "
     "non-abelian currents and the charge sector give no universal degree-2 class; the native second order is Heisenberg, "
     "not the Plebanski-Penrose bracket, hence not the Einstein vertex; the D = +-4 sector is uniformly gapped (scoped)",
     lambda d: _check_named(d, "(c) die Plebanski-Penrose-Klammer") and _check_named(d, "GRAV.PARENT.GAP.NO_GO.01")),
    ("marks", "Euclidean Weyl curvature", "excluded", "z85_ergebnis.json", "the H6 state is the quadrupole diag(a, -a, 0) and "
     "changes sign under the quarter turn; an R-invariant Weyl quadrupole is type D with J != 0 (Bianchi IX, Tod form)",
     lambda d: _check_named(d, "Entscheidungstest S4b")),
    ("time", "OS reconstruction", "conditional", "z84_ergebnis.json", "native reflection-positive reflections exist; for the "
     "thermal original the OS reconstruction returns the modular flow; a positive continuous time arises only in the vacuum "
     "sector with a boundary point t* (as z72)",
     lambda d: _check_named(d, "thermische OS-Rekonstruktion") and _check_named(d, "Halbraum-Transfer x_0 -> x_0 + t")),
    ("marks", "celestial sphere", "excluded", "z84_ergebnis.json", "no factor-swapping map normalises <R> (Satz B) and R has no "
     "fixed point on U(2): no mark-invariant Minkowski chart, T' holds only pointwise",
     lambda d: _check_named(d, "Satz B: eine faktortauschende") and _check_named(d, "auf D(V, eta) ist R = (g_V, g_eta)")),
    ("origin", "two binary observations", "conditional", "z80_ergebnis.json", "sign symmetry of the law forces Cl_2 = M_2, "
     "tau = i, p = 1/2 if implemented in the same irreducible representation",
     lambda d: _check_named(d, "irreduzible Darstellungen von <A, B")),
]


def _evidence(root):
    out = {}
    for axis, value, status, fname, meaning, pred in EVIDENCE:
        path = root / fname
        if not path.exists():
            continue
        ok = pred(json.loads(path.read_text()))
        out[(axis, value)] = {"status": status if ok else "contradicted", "source": fname, "meaning": meaning}
    return out


def scenario_candidates(root=None):
    """All combinations of the six axes with clauses computed from OriginLab's own universes (R orders, Lorentz or
    Euclidean position spaces) and from recorded exploration results (EVIDENCE). Pending evidence gives no clause."""
    root = EXPLORATION if root is None else _Path(root)
    ev = _evidence(root)
    hardy = lambda k0, k1: sp.diag(sp.I**k0, sp.I**k1)
    doublet = {"quarter oscillators": hardy(4, 4), "step type": hardy(4, 4),     # J_-1 = S^4: integer L_0 steps
               "eta doublet": _transport(0, -1), "seam doublet V_1": hardy(0, 1)}
    r_order = {k: _proj_order(M) for k, M in doublet.items()}
    vstar = sp.diag(1, -sp.I)
    space = {"quarter oscillators": (_kron(doublet["quarter oscillators"], doublet["quarter oscillators"].conjugate()), True),
             "step type": (_kron(doublet["step type"], doublet["step type"].conjugate()), True),
             "eta doublet": (_transport(1, -1), True),
             "seam doublet V_1": (_kron(vstar, vstar.conjugate()), False)}             # T': R on Vbar^* imposed
    native_lorentz = {k: native and _lorentz_compatible(_turns(M)) for k, (M, native) in space.items()}
    status = lambda axis, value: ev.get((axis, value), {}).get("status", "pending")
    objs = []
    for combo in product(*SCENARIO_AXES.values()):
        sc = dict(zip(SCENARIO_AXES, combo))
        st = {a: status(a, v) for a, v in sc.items()}
        coherent = not (
            (sc["marks"] == "celestial sphere" and sc["spinor_index"] != "seam doublet V_1")
            or (sc["marks"] == "internal spin-Z4" and r_order[sc["spinor_index"]] != 1)
            or (sc["marks"] == "Euclidean Weyl curvature" and (sc["spinor_index"] != "eta doublet"
                                                                or sc["time"] != "OS reconstruction"))
            or (sc["modes"] == "energy levels" and sc["spinor_index"] != "step type"))
        lorentzian = native_lorentz[sc["spinor_index"]] or (sc["time"] == "OS reconstruction"
                                                            and st["time"] == "supported")
        cl = {c for c, ok in (
            ("admissible", "excluded" not in st.values() and "contradicted" not in st.values()),
            ("coherent", coherent), ("lorentzian_native", lorentzian),
            ("positive_time", st["time"] in ("supported", "conditional")),
            ("one_field", sc["modes"] == "energy levels" and st["modes"] in ("supported", "conditional")),
            ("graviton", st["gravity"] == "supported"),
            ("origin_reduced", st["origin"] in ("supported", "conditional")),
            ("marks_placed", st["marks"] in ("supported", "conditional"))) if ok}
        objs.append({"label": " | ".join(combo), "scenario": sc, "status": st, "clauses": cl})
    return objs, {"evidence": {f"{a}: {v}": e for (a, v), e in sorted(ev.items())}, "r_order": r_order,
                  "native_lorentz": native_lorentz}


def _rules_report(objects, pool, label):
    res = search(objects, pool, lambda o: 0, label)
    hits = {frozenset(r): sorted({label(o) for o in objects if set(r) <= o["clauses"]}) for r in powerset(pool)}
    nonempty = [r for r, h in hits.items() if h]
    empty = [r for r, h in hits.items() if not h]
    order = lambda r: (len(r), sorted(r))
    res["maximal_rules"] = [{"rule": sorted(r), "survivors": hits[r]}
                            for r in sorted(nonempty, key=order) if not any(r < s for s in nonempty)]
    res["minimal_empty_rules"] = [sorted(r) for r in sorted(empty, key=order) if not any(e < r for e in empty)]
    ext = {c: tuple(i for i, o in enumerate(objects) if c in o["clauses"]) for c in pool}
    res["equivalent_clauses"] = [[a, b] for a, b in combinations(pool, 2) if ext[a] == ext[b]]
    res["selected"] = {s["selected"]: s["minimal_rules"] for s in res["unique_selections"]}
    return res, hits


def _question(hits, rule):
    h = hits[frozenset(rule)]
    return {"rule": sorted(rule), "survivors": h, "verdict": "EMPTY" if not h else "UNIQUE" if len(h) == 1 else "AMBIGUOUS"}


BUNDLE_POOL = ["twistor_line", "mark_symmetric", "choice_free", "lorentz_cone"]
DOUBLET_POOL = ["faithful", "kms", "symmetric", "base_point_independent", "mirror_closed", "marks_separated", "neutral"]
SPIN2_POOL = ["j_zero", "i_zero", "zero_eigenvalue", "nondegenerate", "d4_fixed"]
DOUBLET_QUESTION = ["faithful", "kms", "symmetric", "base_point_independent"]
IDENT_POOL = ["spinor", "r_equivalent_to_seam", "r_equivalent_to_quarter_pair", "projective_order_1",
              "projective_order_2", "projective_order_4"]
POSITION_POOL = ["native_r", "lorentz_compatible", "euclidean_compatible", "marked"]
SCENARIO_POOL = ["admissible", "coherent", "lorentzian_native", "positive_time", "one_field", "graviton",
                 "origin_reduced", "marks_placed"]


def run_reconstruction(root=None):
    out = {"status": "RECONSTRUCTION_SEARCHED", "physical_solution": False, "scope": SCOPE, "universes": {},
           "symmetries": symmetry_structure() | {"exceptional_orbits": exceptional_orbits()}}

    nb = normal_bundle_candidates()
    res, hits = _rules_report(nb, BUNDLE_POOL, lambda o: o["name"])
    res["questions"] = {"twistor_line+choice_free": _question(hits, ["twistor_line", "choice_free"]),
                        "twistor_line+mark_symmetric": _question(hits, ["twistor_line", "mark_symmetric"]),
                        "twistor_line": _question(hits, ["twistor_line"])}
    res["derived"] = {o["name"]: {"splitting": o["degrees"], "h0": o["h0"], "h1": o["h1"], "orbit_size": o["orbit_size"],
                                  "cover_derived": o["cover_derived"], "clauses": sorted(o["clauses"]),
                                  "reality": o["reality"]} for o in nb}
    res["cover_genus"] = {"double cover": str(cover_genus(DOUBLE_COVER)), "(Z2)^3 cover": str(cover_genus(CHARACTERS)),
                          "h1(r_*O_E)": _h1(pushforward(DOUBLE_COVER, 0)), "h1(r_*O_C)": _h1(pushforward(CHARACTERS, 0))}
    res["bound"] = ("covers branched at the four marks, halves under all 7 deck involutions, duals, Euler and Tomita "
                    "doublet bundles, up-modifications of the O(-1) summand at the 15 non-empty mark subsets and at the "
                    "non-empty subsets of the poles {0, oo}")
    out["universes"]["normal_bundles"] = res
    twistor = [o for o in nb if "twistor_line" in o["clauses"]]

    ds = doublet_state_candidates()
    res, hits_d = _rules_report(ds, DOUBLET_POOL, lambda o: o["label"])
    res["questions"] = {"+".join(DOUBLET_QUESTION): _question(hits_d, DOUBLET_QUESTION),
                        "kms+neutral": _question(hits_d, ["kms", "neutral"]),
                        "kms+marks_separated": _question(hits_d, ["kms", "marks_separated"]),
                        "kms+mirror_closed": _question(hits_d, ["kms", "mirror_closed"])}
    rho_m = {}
    for m in M_VALUES:
        corners = [o for o in ds if o["construction"] == f"Fock corner, m = {m}"]
        labels = sorted({o["label"] for o in corners})
        rho_m[str(m)] = {"state": labels[0] if len(labels) == 1 else labels, "corners": len(corners),
                         "faithful_kms_base_point_independent": all({"faithful", "kms", "base_point_independent"}
                                                                    <= o["clauses"] for o in corners),
                         "projective_order_of_R": sorted({o["r_order"] for o in corners}),
                         "neutral": all("neutral" in o["clauses"] for o in corners),
                         "same_kms_state_from": sorted({o["construction"] for o in ds if o["label"] in labels
                                                        and "kms" in o["clauses"]})}
    res["derived"] = {"original_symmetries": original_symmetries(), "rho_m": rho_m,
                      "lowest_quarter_pair": [str(min((r + a for r in range(-3, 4)), key=abs)) for a in QUARTERS],
                      "one_particle_pairs_exist_for_m": sorted({o["m"] for o in ds if o["construction"].startswith("one-")}),
                      "fock_multiplicities": {str(a): dict(zip(("patterns", "bosonization"), fock_multiplicities(a)))
                                              for a in QUARTERS},
                      "candidates": [{"construction": o["construction"], "m": o["m"], "base": o["base"], "rho": o["rho"],
                                      "rho_numeric": o["rho_numeric"], "clauses": sorted(o["clauses"])} for o in ds]}
    res["bound"] = ("m = 1..4; mode pairs (x, x + m/4) with x = alpha + r, r in {-2,-1,0,1}, alpha in {1/4, 3/4} (one-particle "
                    "restriction and N = 1 corner); Fock corners of each quarter sector at base marks M0 = 0..4 with every "
                    "state of the two levels; controls: trace state, Hardy vacuum")
    out["universes"]["doublet_states"] = res

    s2, span = spin2_family_candidates()
    res, hits_s = _rules_report(s2, SPIN2_POOL, lambda o: o["label"])
    res["derived"] = {o["t"]: {k: o[k] for k in ("quartic", "I", "J", "discriminant", "eigenvalues", "kernel_phi00_phi01_phi11",
                                                  "petrov", "trace_identities")} | {"clauses": sorted(o["clauses"])} for o in s2}
    res["pencil_basis"] = span
    res["bound"] = "V4 pencil 1 - z^4 + t z^2 at all zeros of I, J, I^3 - 27 J^2 along the pencil, t = +-1, and z^2"
    out["universes"]["spin2_family"] = res

    di, ps = doublet_identification_candidates()
    res, hits_i = _rules_report(di, IDENT_POOL, lambda o: o["label"])
    res["questions"] = {"r_equivalent_to_seam": _question(hits_i, ["r_equivalent_to_seam"]),
                        "spinor+r_equivalent_to_quarter_pair": _question(hits_i, ["spinor", "r_equivalent_to_quarter_pair"]),
                        "r_equivalent_to_quarter_pair": _question(hits_i, ["r_equivalent_to_quarter_pair"])}
    res["derived"] = {o["label"]: {"turns_of_16": o["turns"], "matrix": o["matrix"], "clauses": sorted(o["clauses"])}
                      for o in di}
    res["bound"] = ("R on the seam doublet (Fock corner), the one-particle quarter pair, the current oscillators of both "
                    "quarters, and the eta doublet of the spin bridge (both lifts of z -> iz and the deck), exact on E")
    out["universes"]["doublet_identifications"] = res
    res, hits_p = _rules_report(ps, POSITION_POOL, lambda o: o["label"])
    res["questions"] = {"native_r+lorentz_compatible+marked": _question(hits_p, ["native_r", "lorentz_compatible", "marked"]),
                        "native_r+euclidean_compatible+marked": _question(hits_p, ["native_r", "euclidean_compatible",
                                                                                     "marked"])}
    res["derived"] = {o["label"]: {"turns_of_16": o["turns"], "clauses": sorted(o["clauses"]),
                                   "factorises_as_vstar_x_eta": o["factorises_as_vstar_x_eta"]} for o in ps}
    res["bound"] = ("Tomita reading T' (R induced on Vbar^*), the spin-bridge sections H^0(N) (native R), Herm of the quarter "
                    "oscillators (z78) and of the quarter pair; Lorentz: g (x) conj g, Euclidean: SU(2) x SU(2)")
    out["universes"]["position_spaces"] = res
    sc, sc_info = scenario_candidates(root)
    res, hits_c = _rules_report(sc, SCENARIO_POOL, lambda o: o["label"])
    base = ["admissible", "coherent", "lorentzian_native", "positive_time"]
    res["questions"] = {"complete": _question(hits_c, SCENARIO_POOL), "+".join(base): _question(hits_c, base)}
    best = max(len(o["clauses"]) for o in sc)
    res["derived"] = sc_info | {"most_clauses": best,
                                "best_scenarios": sorted(o["label"] for o in sc if len(o["clauses"]) == best)[:40],
                                "values_in_base_survivors": {a: sorted({o["scenario"][a] for o in sc if set(base) <= o["clauses"]})
                                                             for a in SCENARIO_AXES}}
    res["bound"] = ("product of six axes (spinor index, role of the marks, time, modes, gravity, origin) = "
                    f"{len(sc)} scenarios; clauses from OriginLab's R orders and position spaces and from the recorded "
                    "results named in EVIDENCE; pending evidence gives no clause")
    out["universes"]["scenarios"] = res
    h6 = next(o for o in s2 if o["t"] == "0")
    dq = out["universes"]["doublet_states"]["questions"]["+".join(DOUBLET_QUESTION)]
    bq = out["universes"]["normal_bundles"]["questions"]["twistor_line+choice_free"]
    cover_free = sorted(o["name"] for o in twistor if o["cover_derived"] and "choice_free" in o["clauses"])
    out["answers"] = {
        "some_bundle_is_twistor_line_and_choice_free": bool(bq["survivors"]),
        "cover_derived_choice_free_twistor_lines": cover_free,
        "choice_free_twistor_line_from_a_mark_divisor": any(o["kind"] == "hecke" and max(o["divisor"]) < 4
                                                            and "choice_free" in o["clauses"] for o in twistor),
        "lorentzian_real_form_only_with_conjugate_doublet": all(
            o["fibre"] == "Vbar" for o in nb if "lorentz_cone" in o["clauses"]),
        "rho_m_faithful_kms_base_point_independent_for_every_m": all(v["faithful_kms_base_point_independent"]
                                                                     for v in rho_m.values()),
        "doublet_state_fixed_by_faithful_kms_symmetric_base_point_independent": dq["verdict"] == "UNIQUE",
        "doublet_states_left_by_faithful_kms_symmetric_base_point_independent": len(dq["survivors"]),
        "neutral_doublet_with_separated_marks_exists": bool(hits_d[frozenset(["neutral", "marks_separated"])]),
        "h6_has_j_zero_and_a_vanishing_spin1_eigenvalue": {"j_zero", "zero_eigenvalue"} <= h6["clauses"],
        "quarter_pair_is_the_eta_doublet_not_the_seam_doublet": (
            hits_i[frozenset(["r_equivalent_to_seam"])] == [SEAM_DOUBLET]
            and sorted(hits_i[frozenset(["spinor", "r_equivalent_to_quarter_pair"])]) == sorted([QUARTER_PAIR, ETA_SPINOR])),
        "no_native_marked_lorentzian_position_space": not hits_p[frozenset(["native_r", "lorentz_compatible", "marked"])],
        "spin_bridge_positions_euclidean_not_lorentzian": hits_p[frozenset(["native_r", "euclidean_compatible",
                                                                             "marked"])] == [SPIN_BRIDGE]}
    out["summary"] = {
        "normal_bundles": {
            "twistor_line_and_choice_free": bq["survivors"],
            "cover_derived_twistor_lines": {o["name"]: o["orbit_size"] for o in twistor if o["cover_derived"]},
            "lorentz_cone": sorted(o["name"] for o in nb if "lorentz_cone" in o["clauses"]),
            "seam_reality_signatures": [list(s) for s in sorted({tuple(sorted(r["signature"])) for o in twistor
                                                                 for r in o["reality"] if r["type"] != "hermitian"})],
            "selected": out["universes"]["normal_bundles"]["selected"]},
        "doublet_states": {"rho_m": {m: v["state"] for m, v in rho_m.items()},
                           "faithful_kms_symmetric_base_point_independent": dq["survivors"],
                           "selected": out["universes"]["doublet_states"]["selected"],
                           "minimal_empty_rules": out["universes"]["doublet_states"]["minimal_empty_rules"]},
        "spin2_family": {"h6": {"I": h6["I"], "J": h6["J"], "eigenvalues": h6["eigenvalues"],
                                "zero_mode": h6["kernel_phi00_phi01_phi11"]}}}
    out["rules_generated"] = sum(r["rules"] for r in out["universes"].values())
    out["objects_enumerated"] = sum(r["objects"] for r in out["universes"].values())
    return out


run = run_reconstruction


if __name__ == "__main__":
    print(json.dumps(run_reconstruction(), indent=1, sort_keys=True, ensure_ascii=False, allow_nan=False))
