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
                  Since H6h47 also the statistics of matter (lift of the 2 pi rotation, z92, z97) and the interaction
                  (z91, z99), with clauses spin_statistics, chiral_matter and interacting.

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

from genesis_rules import powerset

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
                "twistor deformation: radial channel", "entanglement: original|partner cut", "entanglement: spatial cut",
                "P0-graded two-seam composition", "MacDowell-Mansouri: modular Cartan connection",
                "double copy of seam amplitudes (P_diag)", "moved algebra at fixed state"],
    "origin": ["seed rho_1 + spin-Z4", "two binary observations"],
    "statistics": ["untwisted lift", "clock-twisted lift (P_GSO)", "E8 glue with A3 = conformal (P_A3)",
                   "Z4 glue on the spinor module (P_glue)", "Z4 glue summed over a pair of seams (P_glue, P_paar)"],
    "interaction": ["fusion over even or linear algebras", "vertex operators = level-1 currents", "chosen two-body term",
                    "twistor action from seam currents (T)", "transfer exchange (twirl)", "local coherent transfer gate"],
}
# the statistics rules of z97 live on charge lattice (x) radial standard form (step type = spinor index, modes = levels)
RADIAL_STATISTICS = ("clock-twisted lift (P_GSO)", "E8 glue with A3 = conformal (P_A3)",
                     "Z4 glue on the spinor module (P_glue)", "Z4 glue summed over a pair of seams (P_glue, P_paar)")
CHIRAL_STATISTICS = {"untwisted lift": False, "clock-twisted lift (P_GSO)": False,       # z97: index 0 per R class
                     "E8 glue with A3 = conformal (P_A3)": True,                         # z97: chiral Z4, 16 one helicity
                     "Z4 glue on the spinor module (P_glue)": True,                      # z117: one chiral 16 with the cut
                     "Z4 glue summed over a pair of seams (P_glue, P_paar)": True}                         # the same glue (z117)
# z113 (H6h49): summing the 2 pi rotation U_2pi over the seam torus, modular invariance plus positivity per seam keeps
# only U = +1 (no spinors); the per-seam rule G = U R^2 = +1 has negative long-seam multiplicities, and no lift is
# compatible with (E8)_1, where the clock is frozen. z174 (H6h50) ran the same test for the Z4 glue X = R i^(-2h) and its
# mirror: on one pair seam + partner the glue is anomalous (no consistent character), so the per-seam value stays without
# the clause; summed diagonally over two pairs it is modular and positive exactly for the eigenvalue -1 (R = Z_SU^-1 on
# all states, normal statistics) and for one eigenvalue with exchanged statistics, which is the pair-summed value
MODULAR_STATISTICS = {"untwisted lift": True, "clock-twisted lift (P_GSO)": False,
                      "E8 glue with A3 = conformal (P_A3)": False, "Z4 glue on the spinor module (P_glue)": False,
                      "Z4 glue summed over a pair of seams (P_glue, P_paar)": True}


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
    ("gravity", "twistor deformation: radial channel", "excluded", "z89_ergebnis.json", "the E-invariant (radial) bivector in "
     "the Lorentz pencil gives only AdS with curvature radius equal to the vacuum size; no flat Einstein bracket, no G (scoped)",
     lambda d: _check_named(d, "Lorentzbüschel") and _check_named(d, "Skala des radialen Kanals")),
    ("gravity", "entanglement: spatial cut", "excluded", "z90_ergebnis.json", "native weights do not regulate the UV, the area "
     "coefficient depends on the localisation, and at the forced xi = 1/6 the induced 1/G vanishes (scoped)",
     lambda d: _check_named(d, "beim nativen xi = 1/6")),
    ("gravity", "P0-graded two-seam composition", "conditional", "z89_ergebnis.json", "the infinity twistor I_pi is native up "
     "to its modular phase and reproduces the Plebanski-Penrose bracket exactly; P_0 is gapless (critical sequence); missing: "
     "a native composition rule with fixed coefficient g and a fixed phase of I",
     lambda d: _check_named(d, "Kritische Folge (GRAV.PARENT.GAP.NO_GO.01)") and _check_named(d, "Universalität bei gegebenem I")),
    ("origin", "seed rho_1 + spin-Z4", "conditional", "z94_ergebnis.json", "seed family and Fermi statistics forced (z87); "
     "rigidity carried by non-degeneracy (S(2) proven, three and four pairs exhaustive in windows); P_clock equivalent to "
     "time neutrality of the law P_tn; two independent premises P_sym and P_tn remain",
     lambda d: _check_named(d, "F6-Kriterien") and _check_named(d, "Unabhängigkeit: P_sym ohne Zeitneutralität")),
    ("origin", "two binary observations", "conditional", "z80_ergebnis.json", "sign symmetry of the law forces Cl_2 = M_2, "
     "tau = i, p = 1/2 if implemented in the same irreducible representation",
     lambda d: _check_named(d, "irreduzible Darstellungen von <A, B")),
    # H6h47: follow-up round G1-G3 and the twelfth external analysis (z96 checks it independently)
    ("gravity", "P0-graded two-seam composition", "conditional", "z98_ergebnis.json", "the phase of I is a gauge direction "
     "(all native data are modular-U(1) invariant, Theta_1 = tau fixes nothing); the second-order OS form is positive for "
     "every g and phase, so positivity fixes reality, not |g|; g stays free",
     lambda d: _check_named(d, "Modul-U(1) gegen die nativen Daten") and _check_named(d, "Reflexionspositivität zweiter Ordnung")),
    ("gravity", "MacDowell-Mansouri: modular Cartan connection", "conditional", "z98_ergebnis.json", "the (A)dS algebra of "
     "Q = (P + sigma K)/2 is native (z89's pencil, sigma = +1, so(3,2) with Q_0 = E, radius = half the vacuum size); the "
     "family of Tomita conjugations is AdS4 with the Levi-Civita connection and flat MacDowell-Mansouri curvature F = 0 "
     "(rigid vacuum); the two sheets give Sigma_+, Sigma_- and a flat tetrad; missing: a non-flat deformation with fixed b, "
     "the Plebanski sector (fixed only up to Z4) and the dS sign",
     lambda d: _check_named(d, "MacDowell-Mansouri-Algebra nativ") and _check_named(d, "tau als Punktspiegelung von AdS4")
     and _check_named(d, "Plebanski aus zwei Blättern") and _check_named(d, "Vorzeichen und Skala")),
    ("statistics", "untwisted lift", "excluded", "z92_ergebnis.json", "with the untwisted radial 2 pi rotation the "
     "half-integer helicity sectors are R-bosons; no native Lorentz action whose -1 is R^2 carries the vertex operators",
     lambda d: _check_named(d, "auf der ganzen Leiter ist L_0 = m") and _check_named(d, "Vertexoperatoren (274 Paare")),
    ("statistics", "clock-twisted lift (P_GSO)", "conditional", "z97_ergebnis.json", "G = U_2pi R^2 is central and closes "
     "under su(2,2), fusion and the dressed vertices (charge step x step type, the step type is the spinor index); not "
     "native (Tomita J only compatible, no modular element equals G, eta leaves two lifts); vector-like (index 0 per R "
     "class); across seams it needs a graded composition",
     lambda d: _check_named(d, "[E] Definition: auf Ladungsgitter") and _check_named(d, "[E] gekleidete Vertizes")
     and _check_named(d, "[E] Fusion: G ist multiplikativ") and _check_named(d, "[negativ, E] Modulgruppe bei beta = 4 pi")),
    ("statistics", "E8 glue with A3 = conformal (P_A3)", "excluded", "z97_ergebnis.json", "reading A3 as the complexified "
     "radial su(2,2) gives the rule and a chiral Z4 rule (the 16 with one helicity), but forbids the neutral radial scalar "
     "ladder of z83 and reinterprets the solved family lattice A3 (scoped)",
     lambda d: _check_named(d, "[C, mit Konflikt] Liest man A3")),
    ("interaction", "fusion over even or linear algebras", "excluded", "z91_ergebnis.json", "no native fusion over even or "
     "linear shared algebras gives a connected 4D three-point structure; cubic vertices between scalar ladders vanish",
     lambda d: _check_named(d, "(ii-d) Kubische Vertizes")),
    ("interaction", "vertex operators = level-1 currents", "excluded", "z99_ergebnis.json", "a literally shared envelope is "
     "non-local or empty; the only local coupling of two seams (lattice Z^2) transfers the su(2)_1 current, so vertex and "
     "current routes coincide; connected three-points only as form factors or sl(2,R) fusion, no covariant scalar 4D "
     "vertex; Delta = 1 rigid, the marginal radius deformation lacks the long-multiplet descendants (scoped)",
     lambda d: _check_named(d, "(a) Wörtlich geteilte Hülle") and _check_named(d, "(b) Frenkel-Kac")
     and _check_named(d, "(c) Delta = 1 bleibt") and _check_named(d, "(c-ii) Die verschobene Kreuzleiter")),
    ("interaction", "chosen two-body term", "conditional", "z99_ergebnis.json", "every native bilinear gives only "
     "disconnected 1->2 parts; the smallest connected 1->2 amplitude is the four-fermion double transfer "
     "e^{2i(phi_1 - phi_2)} (h = 4) with a chosen coupling; whether it supplies the long-multiplet descendants is open",
     lambda d: _check_named(d, "(e) Jeder Bilinearoperator")),
    # H6h48/H6h49: origin on P_sym (z107), deep wave T1-T6 and U3 (z110-z115, z119), eighteenth analysis (z120)
    ("origin", "seed rho_1 + spin-Z4", "conditional", "z107_ergebnis.json", "P_tn follows from the seam premise A* without "
     "circularity (law torus = seam torus, the only stress-free point is beta = 4 pi); the origin rests on P_sym; S(3) holds "
     "iff one pair splits off, rigidity from n = 3 on is arithmetic",
     lambda d: _check_named(d, "N8-Kriterien") and _check_named(d, "Nichtzirkularität (Prämissen)")),
    ("gravity", "double copy of seam amplitudes (P_diag)", "conditional", "z112_ergebnis.json", "one native infinity "
     "twistor gives both brackets; Hodges = KLT for n = 3..6; original x partner gives only the dilaton-axion part, the "
     "sheet-diagonal form needs P_diag; G reduces to one missing length L = kappa/(2 g^2), the Lambda-G relation to a second",
     lambda d: _check_named(d, "Hodges = KLT") and _check_named(d, "Prämissenabtastung")),
    ("gravity", "MacDowell-Mansouri: modular Cartan connection", "conditional", "z118_ergebnis.json", "the (A)dS algebra "
     "is native (z98); the clock-weighted index trace over the joint spinor gives the MacDowell-Mansouri term 64 eps: the "
     "admissible weights are span{Gamma_2, R}, the weight choice is the P_vol pair choice (premise P_R: the clock acts on the "
     "internal factor only) and one pairing bit fixes the sign; Lambda = -12/R_vac^2 (AdS); the trace gives integers only, "
     "G follows only with the weight 2 c_3 put in",
     lambda d: _check_named(d, "MacDowell-Mansouri in der Spur") and _check_named(d, "Normierung (sympy)")),
    ("statistics", "Z4 glue on the spinor module (P_glue)", "conditional", "z117_ergebnis.json", "the native Tomita "
     "J_H6 is the Majorana-Weyl structure of Lorentz x internal Cl(3,11) on the index cut; gluing I_T = I_F, total "
     "chirality +1 and the z97 rule R = Z_SU^-1 restricted to the spinor module are one condition (one bit P_glue, symmetry "
     "(SU(2,2) x Spin(10))/Z4_diag); with the cut one chiral 16 remains; whether the neutral scalar ladder stays outside "
     "the spinor module is open",
     lambda d: _check_named(d, "Verklebung = Gesamtchiralität = Z4-Regel") and _check_named(d, "Schnitt und Helizität")),
    ("statistics", "Z4 glue summed over a pair of seams (P_glue, P_paar)", "conditional", "z174_ergebnis.json",
     "the Z4 glue of z117 summed over the seam torus by the method of z113: anomalous on one pair seam + partner, free on "
     "two; there the projections on X = -1 (R = Z_SU^-1 on all states: vectors clock-neutral, scalars at R = -1, Weyl "
     "fermions at R = +-i, normal statistics) and on one eigenvalue with exchanged statistics are Hilbert space traces, the "
     "sector of the scalar vacuum is not; over four pairs as a whole only exchanged statistics survives",
     lambda d: _check_named(d, "Anomalie: Die einseitige Verklebung")
     and _check_named(d, "Positivität der verklebten Zustandssumme") and _check_named(d, "Was der Eigenwert -1 bedeutet")),
    ("gravity", "moved algebra at fixed state", "excluded", "z115_ergebnis.json", "the family J(M_g, Omega) is an orbit of "
     "the state's stabiliser (compact dual U(2)), flat in its own model; the quadratic response vanishes on the family "
     "directions, condensates spread over the orbit with a free amplitude (scoped)",
     lambda d: _check_named(d, "Lage: Jedes J(M_g, Omega)") and _check_named(d, "Quadratische Antwort im Vakuum")),
    ("statistics", "clock-twisted lift (P_GSO)", "conditional", "z113_ergebnis.json", "G = U_2pi R^2 closes (z97) but is "
     "not modular-native: modular invariance forces only R = Rbar and (-1)^(F + Fbar) = +1, U_2pi is no modular "
     "transformation (its twisted sectors are long seams); if U_2pi is summed, per-seam G = +1 violates positivity "
     "(only a diagonal two-seam rule survives)",
     lambda d: _check_named(d, "U_2pi ist keine Modultransformation") and
     _check_named(d, "Positivität der vollständigen Zustandssumme")),
    ("interaction", "twistor action from seam currents (T)", "conditional", "z111_ergebnis.json", "under T (seam = "
     "twistor line) the su(2)_1 seam currents give the Parke-Taylor denominator with coefficient fixed by the level "
     "(z110); the CSW reference spinor is native from I_pi; missing: a free exchange field on twistor space (kinetic "
     "term, P_BF), the numerator <ij>^4 (P_N4) and a suppression of the multi-trace terms",
     lambda d: _check_named(d, "CSW = BCFW") and _check_named(d, "eta aus I_pi")),
    ("interaction", "transfer exchange (twirl)", "excluded", "z119_ergebnis.json", "averaging over the transfer group "
     "(P_twirl) gives a number-conserving modular weight log(dim W/f) (log 3 = log dim isospin 1), no vertex, gamma = 0; "
     "the global average gives no consistent family of local states (z120); as an interaction mechanism excluded (scoped)",
     lambda d: _check_named(d, "(b-ii) Vergleich z105/z108") and _check_named(d, "(c-ii) Keine verbundene 1->2-Amplitude")),
    ("interaction", "local coherent transfer gate", "conditional", "z120_ergebnis.json", "U = P_T - i P_S (the XXX "
     "R-matrix at u = 1) on one excitation per seam; leak-freeness forces the quarter turn (= P_sym), the phase is the root "
     "of the clock, equivalently the Z4 mutual phase i^(Q_A Q_B); with the free energy the gate is trivial, so it needs "
     "the spin-blind charge energy as native seam response; continuum and cubic vertex open",
     lambda d: _check_named(d, "Transferdrehung B_theta") and _check_named(d, "Auf dem Sektor N = 2")),
    # H6h50: wave V (z121, z123, z124)
    ("interaction", "local coherent transfer gate", "conditional", "z121_ergebnis.json", "natively trivial: E^nat is "
     "one-body in all four readings of the seam modes; under P_null (composition sees only the charge lattice) the gate is "
     "the square root of the seam exchange, U^2 = -Gamma(pi); only theta = pi/2 is forced, every phase is leak-free and "
     "gives the XXX family R(u), u = 1 needs P_viertel; the brickwork is the integrable Destri-de Vega light-cone lattice "
     "(su(2)_1 reading); three seams give a T-odd f-channel c s_A.(s_B x s_C), no production, no phi^3",
     lambda d: _check_named(d, "Vier Lesarten der zwei Moden je Naht") and _check_named(d, "Inhomogene Transfermatrix")
     and _check_named(d, "Drei Nähte A, B, C")),
    ("interaction", "twistor action from seam currents (T)", "conditional", "z124_ergebnis.json", "the kinetic half is "
     "native in form, uniqueness and normalisation modulus: the Tomita pairing with J E J = -(tau E tau)^+ gives the CSW "
     "propagator and CSW = BCFW (z123), leaving P_X; the numerator <ij>^4 has no native source (no native degree twist, the "
     "Pati-Salam quartet lacks homogeneity +1), and the multi-trace sector is conformal supergravity of equal strength "
     "(z124): pure Yang-Mills is not natively reachable on this route",
     lambda d: _check_named(d, "Ladungsanomalie im gaußschen") and _check_named(d, "Der native Doppelspuranteil bei n = 4")),
    # H6h50: waves W-Z and G (z127-z137) and the nineteenth analysis (z144); the last entry per value is the current one
    ("statistics", "Z4 glue on the spinor module (P_glue)", "conditional", "z129_ergebnis.json", "the index route to the "
     "bit is closed: the leaf swap is a spatial mirror, the bit P_IF is P-odd, C-odd and CP-even while the seam index is "
     "invariant; the MacDowell-Mansouri sign is the same bit. In the lattice reading (z144) the bit is the choice between "
     "the two isomorphic gluings D5 + A3 -> E8, exchanged by the charge conjugation of one side",
     lambda d: all(d["checks"].values())),
    ("gravity", "MacDowell-Mansouri: modular Cartan connection", "conditional", "z127_ergebnis.json", "premises of the "
     "clock trace: P_R is mostly a consequence (the finite modular weights on S+ are +-R and +-i Gamma_2, the zero mode is "
     "the KMS ladder), left are the reading choices P_mod and P_idx; the pairing bit is not fixed by the KMS direction "
     "(the native half-turn keeps all KMS data and flips it); Lambda G = -3/(1024 pi kappa)",
     lambda d: all(d["checks"].values())),
    ("gravity", "MacDowell-Mansouri: modular Cartan connection", "conditional", "z136_ergebnis.json", "G/R_vac^2 stays a "
     "free pure number: clock trace (z127), entanglement equilibrium (z136: K_cap = 2 pi Q^- native, but the seam K is "
     "elliptic and the first law holds only globally), induced gravity (z137: sign necessarily AdS, G Lambda = -3 pi/(4 c L), "
     "no graviton pole) and the double copy (z112) give rows that do not agree, all of Planck size; the ledger value needs "
     "9.6e122 seams per horizon. State profiles carry no curvature (z128, one stabiliser orbit); the sign of the action is "
     "the chirality bit P_IF (z129)",
     lambda d: all(d["checks"].values())),
    ("interaction", "local coherent transfer gate", "conditional", "z133_ergebnis.json", "P_null does not follow from the "
     "two-factor split of z82; in the phase family it is equivalent to P_sl2 (the gate commutes strictly with the native "
     "sl(2), b in 4Z); on Q_A = Q_B the gate is 1 (x) sqrt(SWAP) on the oscillators; U^2 = exchange only at even total "
     "charge; P_viertel and P_ord stay", lambda d: all(d["checks"].values())),
    ("interaction", "local coherent transfer gate", "conditional", "z135_ergebnis.json", "loops of the gate: gauge "
     "covariance and the fundamental Wilson functional of the 2D quarter-rotation network are native, sin^2 theta_W = 3/8 "
     "at the matching scale; the absolute coupling (three readings over a factor 16384) and link dynamics are not",
     lambda d: all(d["checks"].values())),
    ("interaction", "local coherent transfer gate", "conditional", "z144_ergebnis.json", "nineteenth analysis: the mutual "
     "phase of the gate at phi = pi/2 is the Z4 bicharacter i^(Q_A Q_B) = B_D5 on the clock values; a composition phase is "
     "a function of the clock values exactly for phi in (pi/2)Z without oscillator part (P_uhr, containing P_sl2), and from "
     "transfer spin 2 on only these phases keep the pair sector Q_A = Q_B leak-free and inside the XXX family (z121's "
     "'every phase is leak-free' holds for one excitation per seam); k = 0 trivial, k = 2 the one-body exchange, k = +-1 "
     "interacting. The constant braid relation holds only for one-body operators (seam clock and beam splitter = Ising "
     "braid matrices, finite Clifford-type group), the interacting gate violates it with defect sqrt(3/8); Yang-Baxter "
     "with spectral parameter (z121) needs level matching mod 4 for the full gates",
     lambda d: _check_named(d, "Funktion der Uhrwerte (Q mod 4)") and _check_named(d, "Leckfreiheit des Paarsektors")
     and _check_named(d, "verletzt die Zopfrelation; im ausgeglichenen Sektor")),
]

# results that feed a clause instead of an axis value: (clause, result file, meaning, predicate)
CLAUSE_EVIDENCE = [
    ("modular_statistics", "z113_ergebnis.json", "statistics rule compatible with modular invariance and positivity "
     "per seam when the 2 pi rotation is summed (MODULAR_STATISTICS)",
     lambda d: _check_named(d, "Positivität der vollständigen Zustandssumme")),
    ("twistor_line_index", "z114_ergebnis.json", "chirality as index on the twistor lines: among 256 line rules only "
     "P_RJ and its mirror are J-covariant; under T, P_GSO, K2 and P_RJ the net index is 3 per 16",
     lambda d: _check_named(d, "Alle 256 Linienregeln") and _check_named(d, "Nullmodenzählung je 16 unter P_RJ")),
]


# H6h50: the premises a scenario still rests on. Each value names what its reading assumes (a supported value only if it
# needs one); a conditional value without entry keeps its own name, so nothing drops out silently. One bit under two
# names counts once (P_IF = P_glue, z117, z122). A recorded implication replaces a premise when its result holds.
VALUE_PREMISES = {
    ("spinor_index", "step type"): ("I_step",),
    ("marks", "internal spin-Z4"): ("P_lift",),
    ("time", "parabolic Wiesbrock"): ("P_vac",),
    ("time", "OS reconstruction"): ("P_vac",),
    ("time", "radial conformal"): ("I_step",),
    ("modes", "energy levels"): ("I_step",),
    ("gravity", "P0-graded two-seam composition"): ("P_G", "P_curv"),
    ("gravity", "MacDowell-Mansouri: modular Cartan connection"): ("P_R", "P_pair", "P_G", "P_curv"),
    ("gravity", "double copy of seam amplitudes (P_diag)"): ("T", "P_diag", "P_G"),
    ("origin", "seed rho_1 + spin-Z4"): ("A*", "P_sym"),
    ("origin", "two binary observations"): ("A*", "P_sym"),
    ("statistics", "clock-twisted lift (P_GSO)"): ("P_CAR", "P_GSO"),
    ("statistics", "Z4 glue on the spinor module (P_glue)"): ("P_IF",),
    ("statistics", "Z4 glue summed over a pair of seams (P_glue, P_paar)"): ("P_IF", "P_paar"),                    # z174
    ("interaction", "chosen two-body term"): ("P_2body",),
    ("interaction", "twistor action from seam currents (T)"): ("T", "P_N4", "P_norm", "P_X"),          # z123, z124
    ("interaction", "local coherent transfer gate"): ("P_null", "P_viertel", "P_ord"),               # z121
}
CLAUSE_PREMISES = {"twistor_line_index": ("T", "K2", "P_RJ")}
# (condition on the scenario, premise replaced, replacement, result file, meaning, predicate on the result json)
IMPLICATIONS = [
    (lambda sc: True, "P_RJ", ("P_IF",), "z122_ergebnis.json", "the native pairing P_IF forces the line rule P_RJ with "
     "fixed sign (256 -> 16 -> 3 -> 1 rules)",
     lambda d: _check_named(d, "[E, bedingt] Native Paarung als Bedingung") and _check_named(d, "[E unter P_IF] Vorzeichen")),
    (lambda sc: sc["interaction"] == "local coherent transfer gate", "P_sym", ("P_univ",), "z125_ergebnis.json",
     "a nontrivial leak-free composition under seam-parity superselection forces the quarter turn, i.e. P_sym",
     lambda d: _check_named(d, "Superauswahl durch Z4-Gradierungen") and _check_named(d, "[C] Auftragskriterien")),
    (lambda sc: sc["interaction"] == "local coherent transfer gate", "P_null", ("P_uhr",), "z144_ergebnis.json",
     "a composition phase that is a function of the clock values has no oscillator part (chi = 0 mod 2 pi), which is "
     "P_sl2 = P_null in the phase family of z133",
     lambda d: _check_named(d, "Funktion der Uhrwerte (Q mod 4)")),
    (lambda sc: sc["interaction"] == "local coherent transfer gate", "P_viertel", ("P_uhr",), "z144_ergebnis.json",
     "a composition phase that is a function of the clock values has phi in (pi/2)Z; k = 0 is trivial and k = 2 the "
     "one-body exchange, so an interacting gate has phi = +-pi/2 (equivalently: strict leak-freedom P_Q from transfer "
     "spin 2 on)",
     lambda d: _check_named(d, "Funktion der Uhrwerte (Q mod 4)") and _check_named(d, "Leckfreiheit des Paarsektors")
     and _check_named(d, "k gerade: das Gatter normalisiert")),
    (lambda sc: sc["interaction"] == "local coherent transfer gate", "P_univ", ("P_std",), "z130_ergebnis.json",
     "with the Tomita partner as second seam the conservation of seam parity is level matching of the doubled seam, "
     "native from modular invariance (z113, z130)",
     lambda d: _check_named(d, "[E] Physikalischer Raum der verdoppelten Naht")
     and _check_named(d, "[E] Besetzung: Der KMS-Vektor")),
    (lambda sc: sc["interaction"] == "local coherent transfer gate", "P_uhr", ("P_std",), "z174_ergebnis.json",
     "the gate maps the standard vector of seam + partner into the clock orbit of a native standard vector (bosonic or "
     "fermionic) exactly for phi in (pi/2)Z, and with oscillators for chi in pi Z, where chi = pi is the 2 pi rotation "
     "that fixes the vector; the quarter phase is the Klein twist of the graded standard form",
     lambda d: _check_named(d, "Uhrbahn der Standardvektoren") and _check_named(d, "Mit Oszillatoranteil")
     and _check_named(d, "Klein-Identität")),
    (lambda sc: True, "P_lift", ("P_zentral",), "z175_ergebnis.json",
     "among the mu_4 elements of the torus exactly 1, (-1)^F and the two lifts c, c^-1 commute with so(10) + so(6) of the "
     "two deck halves; order 4 plus commuting gives the clock with Tr R = 0 and R^2 = (-1)^F, the lift branch is a "
     "convention and 'modulo G_SM' is no longer a separate choice",
     lambda d: _check_named(d, "A5 [E] Zentralisatorsatz") and _check_named(d, "A9 [E] Was P_lift")
     and _check_named(d, "A4 [E] Die gewählte Klasse liefert")),
    (lambda sc: sc["gravity"] == "MacDowell-Mansouri: modular Cartan connection", "P_R", ("P_mod", "P_idx"),
     "z127_ergebnis.json", "the clock as trace weight follows from the finite modular weights and the KMS ladder up to "
     "the reading choices P_mod and P_idx", lambda d: all(d["checks"].values())),
    (lambda sc: sc["gravity"] == "MacDowell-Mansouri: modular Cartan connection", "P_pair", ("P_IF",),
     "z129_ergebnis.json", "the sign of the MacDowell-Mansouri action is odd under the leaf swap: the same bit as "
     "chirality", lambda d: all(d["checks"].values())),
]
PREMISE_MEANING = {
    "A*": "seam premise: no background geometry, base changes are symmetries and <T> = 0 (H6h4); gives P_tn (z107)",
    "P_sym": "the quarter turn A -> B -> -A of the two binary operations is inner, tau = i (z80, z94)",
    "P_univ": "every physical composition conserves the fermion parity of each seam, the pair clock is superselected (z125); "
              "its native form is level matching, native when the second seam is the Tomita partner (P_TP, z130)",
    "P_std": "the composition is a morphism of the graded standard form of one seam: its second seam is the Tomita "
             "partner, and it maps the standard vector into the clock orbit of a native standard vector (z130, z174); "
             "gives P_univ and P_uhr; not derived from A*, and the gate is then a holomorphic-antiholomorphic coupling; "
             "its pair-clock part is under A* equivalent to 'no one-sided current mixes the deck halves' and to 'the 16 is "
             "a local field of the pair'; A*, time neutrality and the covering geometry do not give it: six counter-models "
             "in two spectra (z183)",
    "P_paar": "the statistics rule (the Z4 glue) is one-sided and summed diagonally over separately closed seams, each with "
              "its partner; then anomaly, positivity and normal statistics leave exactly two pairs among N = 1..12 (z174, "
              "z177; the same anomaly of order 2 in the sector object, z176); the pair is neither the sheet doubling of P1 "
              "(z176) nor a cyclic gluing of two copies (z177); the evenness of N is the associator class 2 in Z4 of the "
              "one-sided clock and holds for every N, 'exactly two' is numerical for N <= 24 (z178)",
    "I_step": "identification doubled diagonal step = native single step: spinor index = step type, modes = levels (z83)",
    "P_lift": "the clock is the Z4 of the 16 on the carrier through a spinor lift branch (z82)",
    "P_zentral": "the clock has order 4 and commutes with the base changes of the sections of both deck halves (so(10) + "
                 "so(6) of the 16 seam parents, family-number U(1) included); gives R = c^+-1, Tr R = 0, R^2 = (-1)^F and "
                 "the statement of P_lift (z175); not derived from A*, and Tr R = 0 does not follow from <T> = 0; "
                 "equivalent to P_deck, 'the clock is the deck parity of the sections extended to the seam': that "
                 "extension exists in exactly two mutually inverse forms, both of order 4 (z179)",
    "P_vac": "physical time is read in the vacuum sector; the thermal original has no positive continuous time (z72, z84)",
    "P_CAR": "shared CAR of the seams (z106)",
    "P_GSO": "clock-twisted lift G = U_2pi R^2 = +1 (z97); not modular-native, per seam incompatible with a summed U_2pi "
             "(z113)",
    "P_RJ": "J-covariant line rule: R = +-i couple through opposite Tomita lines (z114)",
    "P_IF": "native pairing Ad Z_SU = Ad R on fermions (z122), the same bit as the spinor-module glue P_glue (z117) and as "
            "the sign of the MacDowell-Mansouri action (z129); P- and C-odd, not fixed by the index (z129); in the lattice "
            "reading the choice between two isomorphic gluings D5 + A3 -> E8 (z144)",
    "T": "hypothesis T: the seam is a twistor line with current algebra (z110-z112, z114)",
    "K2": "product coupling of the family zero modes with the twistor line (z114)",
    "P_N4": "four Grassmann-odd O(1) directions for the numerator <ij>^4 (z110); no native degree twist or quartet (z124)",
    "P_norm": "4D states normalised by the current two-point function (z110, z123)",
    "P_X": "a free field with the native Tomita kinetic form is exchanged between seams (z123; replaces P_BF of z111)",
    "P_null": "the composition sees only the charge lattice (inner energy Q^2/2, not the oscillators); natively the gate "
              "is trivial (z121); in the phase family equivalent to P_sl2, strict commutation with the native sl(2) (z133)",
    "P_uhr": "the composition phase is a function of the seam clock values (it brings no phase independent of the clock); "
             "gives P_null and, for an interacting gate, P_viertel; equivalent from transfer spin 2 on to strict "
             "leak-freedom P_Q: seam charges are exchanged at most, never redistributed (z144)",
    "P_mod": "reading choice of the clock trace: which modular weight is traced (z127)",
    "P_idx": "reading choice of the clock trace: the index convention of the zero-mode plane (z127)",
    "P_viertel": "the gate phase is phi = pi/2 (u = 1, the self-dual point of the XXX family); only theta is forced (z121)",
    "P_ord": "distinguishable, ordered seams in the chain (z121)",
    "P_2body": "a chosen two-body term with chosen coupling (z99)",
    "P_G": "the value of the gravitational coupling: g (z98), the length L (z112) or the weight 2 c_3 (z118); a free pure "
           "number G/R_vac^2 that no native route fixes (z127, z136, z137), minimal form P_eta (a finite entropy per area)",
    "P_curv": "a non-orbit family carrying curvature with fixed coefficient (z100, z115); state profiles do not give one "
              "(z128)",
    "P_R": "the clock is the trace weight on the internal factor only (z118)",
    "P_pair": "pairing bit: the light-cone sector of the zero mode fixes the sign of the action (z118)",
    "P_diag": "sheet-diagonal double copy (z112)",
}


def _implications(root):
    out = []
    for cond, old, new, fname, meaning, pred in IMPLICATIONS:
        path = root / fname
        ok = path.exists() and pred(json.loads(path.read_text()))
        out.append({"replaces": old, "by": list(new), "source": fname, "meaning": meaning, "holds": ok, "_cond": cond})
    return out


# premises that set a free continuous number (K14: no free parameters) weigh before any discrete premise
FREE_PARAMETER_PREMISES = frozenset({"P_2body", "P_G"})


def _premise_closure(objects, rule):
    """Premises of the scenarios that satisfy a rule: first the fewest free parameters, then the fewest premises any of
    them needs, the sets that achieve it, how many survivors need each premise, and the clauses no survivor reaches."""
    surv = [o for o in objects if set(rule) <= o["clauses"]]
    out = {"rule": list(rule), "survivors": len(surv)}
    if not surv:
        return out
    key = lambda o: (len(o["premises"] & FREE_PARAMETER_PREMISES), len(o["premises"]))
    least = min(key(o) for o in surv)
    count = {}
    for o in surv:
        for p in o["premises"]:
            count[p] = count.get(p, 0) + 1
    return out | {"least_free_parameters": least[0], "least_premises": least[1],
                  "minimal_premise_sets": sorted({tuple(sorted(o["premises"])) for o in surv if key(o) == least}),
                  "best": sorted(o["label"] for o in surv if key(o) == least),
                  "premise_frequency": dict(sorted(count.items())),
                  "meaning": {p: PREMISE_MEANING.get(p, "unnamed (conditional value without entry)") for p in sorted(count)},
                  "clauses_no_survivor_has": [c for c in SCENARIO_POOL if not any(c in o["clauses"] for o in surv)]}


def _clause_evidence(root):
    out = {}
    for clause, fname, meaning, pred in CLAUSE_EVIDENCE:
        path = root / fname
        if path.exists():
            out[clause] = {"status": "conditional" if pred(json.loads(path.read_text())) else "contradicted",
                           "source": fname, "meaning": meaning}
    return out


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
    """All combinations of the eight axes with clauses computed from OriginLab's own universes (R orders, Lorentz or
    Euclidean position spaces) and from recorded exploration results (EVIDENCE, CLAUSE_EVIDENCE). Pending evidence gives
    no clause."""
    root = EXPLORATION if root is None else _Path(root)
    ev = _evidence(root)
    cev = _clause_evidence(root)
    holds = lambda clause: cev.get(clause, {}).get("status") == "conditional"
    impl = _implications(root)
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
            or (sc["modes"] == "energy levels" and sc["spinor_index"] != "step type")
            or (sc["statistics"] in RADIAL_STATISTICS and (sc["spinor_index"] != "step type"
                                                          or sc["modes"] != "energy levels")))
        lorentzian = native_lorentz[sc["spinor_index"]] or (sc["time"] == "OS reconstruction"
                                                            and st["time"] == "supported")
        cl = {c for c, ok in (
            ("admissible", "excluded" not in st.values() and "contradicted" not in st.values()),
            ("coherent", coherent), ("lorentzian_native", lorentzian),
            ("positive_time", st["time"] in ("supported", "conditional")),
            ("one_field", sc["modes"] == "energy levels" and st["modes"] in ("supported", "conditional")),
            ("graviton", st["gravity"] == "supported"),
            ("origin_reduced", st["origin"] in ("supported", "conditional")),
            ("marks_placed", st["marks"] in ("supported", "conditional")),
            ("spin_statistics", st["statistics"] in ("supported", "conditional")),
            # z114 (H6h49): under P_GSO the twistor-line index with P_RJ gives net index 3 per 16 (premises T, K2, P_RJ)
            ("chiral_matter", CHIRAL_STATISTICS[sc["statistics"]] or (sc["statistics"] == "clock-twisted lift (P_GSO)"
                                                                       and holds("twistor_line_index"))),
            ("interacting", st["interaction"] in ("supported", "conditional")),
            ("modular_statistics", holds("modular_statistics") and MODULAR_STATISTICS[sc["statistics"]])) if ok}
        prem = set()
        for a, v in sc.items():
            if (a, v) in VALUE_PREMISES:
                prem.update(VALUE_PREMISES[(a, v)])
            elif st[a] == "conditional":
                prem.add(f"{a}: {v}")
        if sc["statistics"] == "clock-twisted lift (P_GSO)" and "chiral_matter" in cl:
            prem.update(CLAUSE_PREMISES["twistor_line_index"])
        for im in impl:
            if im["holds"] and im["replaces"] in prem and im["_cond"](sc):
                prem = (prem - {im["replaces"]}) | set(im["by"])
        objs.append({"label": " | ".join(combo), "scenario": sc, "status": st, "clauses": cl, "premises": frozenset(prem)})
    return objs, {"evidence": {f"{a}: {v}": e for (a, v), e in sorted(ev.items())}, "clause_evidence": cev,
                  "implications": [{k: v for k, v in im.items() if k != "_cond"} for im in impl],
                  "r_order": r_order, "native_lorentz": native_lorentz}


# H6h50 (z173): the complete solution as the target of the same machine. A scenario is complete when it has every clause
# of SCENARIO_POOL, rests on the axiom alone and sets no free number. The distance of a scenario is the triple (missing
# clauses, free numbers, premises beyond the axiom). A convention is a bit whose two values give isomorphic scenarios;
# it is reported, not counted.
AXIOM = "A*"
CONVENTIONS = [
    ("P_IF", "z144_ergebnis.json", "the two gluings D5 + A3 -> E8 are exchanged by the charge conjugation of one side and "
     "both give E8: the bit says which of two mirror-isomorphic worlds is called left",
     lambda d: _check_named(d, "tauschen die beiden Lagrangeschen") and _check_named(d, "beide Verklebungen geben 240 Wurzeln")),
]
# what would close each open item, and what is recorded about it (question: composition, identification, gravity,
# statistics); recorded candidates are not implications until a result shows them
FRONTIER = {
    "graviton": ("gravity", "a spin-2 response derived from the seam with fixed coefficient (a gravity value with status "
                 "supported); closes together with P_curv and P_G; recorded: the glued single object has no clock-neutral "
                 "spin 2, h = 2 with R = 1 exists only as a state of two objects (z177); products of massless "
                 "representations of su(2,2) contain no massless representation, so no pole arises as a collective "
                 "excitation, and the elementary h = +-2 state carries no bilinear stress tensor (z182)"),
    "modular_statistics": ("statistics", "a statistics rule that is chiral and modular; recorded candidate: the 16 as "
                           "clock-defect endpoint of a left-right completion glued by the Tomita conjugation (z170, z130), "
                           "not yet a scenario value"),
    "I_step": ("identification", "the doubled diagonal step is the native single step as a theorem (z83); recorded: "
               "half-sided modular inclusions fix the unit step (the mark spacing) only in the vacuum of the base, the "
               "standard form of the state of A* is type I and has none (z180)"),
    "P_lift": ("identification", "the seam clock is the Z4 of the 16 as a theorem; recorded candidate: trace zero pins "
               "the clock among the classes of E8 (z170), which trades P_lift for the requirement Tr R = 0"),
    "P_zentral": ("identification", "centrality of the clock as a theorem: that 'base changes are symmetries' of A* also "
                  "means the section bases of the two deck halves (z175); recorded: Tr R = 0 follows from centrality and "
                  "from the mark count n = 4, not from <T> = 0; restated as P_deck: what is left is the identification of "
                  "the clock of the composition with the deck parity of the sections (z179)"),
    "P_univ": ("composition", "seam parity conserved by every composition as a theorem; native when the second seam is "
               "the Tomita partner (z130); one partner is the only real Lorentzian net of the catalogue (z172)"),
    "P_uhr": ("composition", "the composition phase is a function of the clock values as a theorem (z144)"),
    "P_ord": ("composition", "ordered distinguishable seams as a theorem (z121); recorded: S reverses the orientation of "
              "the circle through the four marks and fixes the state of A*, so the order needs a bit that breaks S (z180), "
              "and pentagon, hexagon and Yang-Baxter for three and four seams require no order (z181)"),
    "P_std": ("composition", "the composition as a morphism of the graded standard form derived from A*; recorded: no "
              "member of the gate family that fixes the standard vector interacts, and the only non-product symmetry of "
              "the standard form is the one-body side exchange (z174), so the interaction needs this premise; the "
              "standard form of a pair of seams admits the interacting gate as a symmetry but selects neither angle nor "
              "phase (z177); the pair-clock part is the statement that every continuous symmetry of the pair preserves the "
              "split 10 | 6, and the covering of the four marks gives the joint spin structure, not the gauge sum (z183)"),
    "P_paar": ("statistics", "why the glue is one-sided and summed over separately closed seams; recorded: the number two "
               "then follows (N = 2 alone among 1..12 has normal statistics, z177), the sheets of P1 and the cyclic gluing "
               "are excluded as its source (z176, z177); evenness holds for every N, the number two is numerical for "
               "N <= 24 (z178)"),
    "P_curv": ("gravity", "a family carrying curvature with fixed coefficient derived from the seam (z100, z115, z128)"),
    "P_G": ("gravity", "the pure number G/R_vac^2 derived; recorded: one relation eps^2 G/a^2 = 4 c3^3 exp(-2/alpha) from "
            "the native infinity twistor and the ledger rung (z164, z172), the second relation is missing; G Lambda = "
            "12 X is already fixed, the free number is G/a^2 alone; under the area rule G/a^2 = 1/(4 s_1) with the "
            "entropy per seam s_1, which the spectrum gives in three physically different values (z182)"),
}


# rows of the H6 bench that the scenario pool does not carry as clauses
BENCH_OPEN = {
    "K19": ("origin", "what forces the seam with its defect in the state of the original process; the defect sector is "
            "modularly required but not energetically preferred (z170); no one-sided boundary condition selects it, and "
            "under P_std and A* the state contains the clock-twisted sectors if superselection of the pair clock is read "
            "as the gauge sum (z176)"),
    "K20": ("origin", "one parent operator from which seam, clock and partner follow; recorded: (Delta, J) alone do not fix "
            "the split seam | partner, and in the sector object T = rho^(-i) gives clock, defect and sheet exchange under "
            "the same reading (z176)"),
}


# H6h51 (z178): bench rows that follow from inventory and from premises the best scenarios already carry. Such a row
# is no independent open item any more: it closes when its premises close. It is not closed, and it is listed with the
# premises it rests on. Inventory used: the GNS rule of the process structure (status E) and A*.
BENCH_DERIVED = [
    ("K19", ("P_std",), (("z178_ergebnis.json",
     "superselection of the pair clock (P_std) and the GNS rule give the Gauss constraint U = 1 on the state space as a "
     "theorem of the standard form; S and T invariance from A* then leave exactly one completion over the clock-neutral "
     "algebra, the diagonal one with all 16 clock sectors, with fixed weights (12 sqrt(2) - 13)/14 in the sectors of odd "
     "clock flux. The reading 'superselection = gauge sum' of z176 is replaced by this derivation",
     lambda d: _check_named(d, "A1 [E] Die uhrneutrale Algebra") and _check_named(d, "A3 [E] Modulgleichung")
     and _check_named(d, "B1 [E, drei Modelle] Standardform") and _check_named(d, "B2 [E] Gegenmodelle")
     and _check_named(d, "F1 [E] Kette")),)),
    ("K20", ("P_std", "P_zentral"), (("z178_ergebnis.json",
     "the parent is the state of A* on the torus: its one-sided modulator gives T = rho^(-i), with the fusion that is now "
     "derived (K19) this gives S, the sheet exchange S^2 = C and the clock as the monodromy with one sector (z176 C4); "
     "the split seam | partner is the pair of opposite arcs of the four marks, two parallel cycles at the distance of "
     "the half period (z180). Which sector carries the monodromy is the identification of the clock, so the row rests on "
     "P_zentral as well",
     lambda d: _check_named(d, "F1 [E] Kette")),
     ("z176_ergebnis.json", "", lambda d: _check_named(d, "C4 [E] Sektorenobjekt")),
     ("z180_ergebnis.json", "", lambda d: _check_named(d, "B1 [N] Zweipunktfunktion") and _check_named(d, "B2 [E] Modulare Flüsse")))),
]


# H6h50 and H6h51 (z175 to z182): candidates for one object behind several open items, each with a computed test.
# Verdicts: "excluded" (the test fails, with the scope in the meaning), "holds" (the test passes and the content of the
# item shrinks or is restated, the count stays), "reading" (the items close under a named reading that is not derived),
# "derived" (the item follows from inventory and from premises that are already counted, see BENCH_DERIVED; it stops
# being an independent open item and is not closed).
UNIFICATION_CANDIDATES = [
    ("pair of seams = the two sheets of P1", ("P_paar",), "z176_ergebnis.json", "excluded",
     "the sheet exchange is the full charge conjugation S^2 = C, not R^2, not the one-sided conjugation and not the side "
     "exchange; the two sheets give one torus trace",
     lambda d: _check_named(d, "B1 [E] Blatttausch") and _check_named(d, "B3 [N")),
    ("pair of seams = two cyclically glued copies (opening angle 4 pi, arXiv:2611.34924)", ("P_paar",),
     "z177_ergebnis.json", "excluded",
     "the cyclic gluing is consistent without clock and with G = U R^2, but has no consistent character with the Z4 glue",
     lambda d: _check_named(d, "Zyklische Verklebung zweier Kopien")),
    ("number of seams from anomaly, positivity and normal statistics", ("P_paar",), "z177_ergebnis.json", "holds",
     "N = 2 alone among N = 1..12; what remains of P_paar is that the glue is one-sided and summed over separately "
     "closed seams",
     lambda d: _check_named(d, "Nahtzahl:")),
    ("standard form of the pair of seams carries the composition", ("P_std", "P_ord"), "z177_ergebnis.json", "excluded",
     "the pair admits the interacting gate as a symmetry of its standard form in the free quarter state but selects "
     "neither angle nor phase and does not conserve the parity of each seam; the partner reading selects the phase and "
     "carries no interaction (z174)",
     lambda d: _check_named(d, "Kanonische Implementierung") and _check_named(d, "Wechselwirkung: Die Mitglieder")),
    ("superselection of the pair clock = gauge sum", ("K19", "K20"), "z176_ergebnis.json", "reading",
     "under P_std and A* the state then contains the clock-twisted sectors (diagonal completion with all 16 sectors) and "
     "T = rho^(-i) is the parent of clock, defect and sheet exchange; the reading, the sector of the 16 and the split "
     "seam | partner are put in",
     lambda d: _check_named(d, "A7 [E, entscheidend] Eichsumme der Paaruhr") and _check_named(d, "C4 [E] Sektorenobjekt")
     and _check_named(d, "C2 [E] (Delta, J) allein")),
    ("defect as ground state of a boundary condition (M13)", ("K19",), "z176_ergebnis.json", "excluded",
     "none of the 16 one-sided boundary conditions of the sector object has ground level zero or a positive state with "
     "<T> = 0 in a twisted sector",
     lambda d: _check_named(d, "A4 [") and _check_named(d, "A3 [E/N] Einseitige positive")),
    ("trace zero of the clock from <T> = 0", ("P_lift",), "z175_ergebnis.json", "excluded",
     "sector by sector the condition fails for all 138 classes of E8, pairwise it holds for all, and no class has a "
     "vanishing ground-level balance",
     lambda d: _check_named(d, "B3 [N]") and _check_named(d, "B4 [N]") and _check_named(d, "B6 [E]")),
    ("clock from centrality (P_zentral)", ("P_lift",), "z175_ergebnis.json", "holds",
     "order 4 and commuting with the base changes of both deck halves gives the clock, Tr R = 0 and R^2 = (-1)^F; the "
     "premise is restated, the count stays",
     lambda d: _check_named(d, "A5 [E] Zentralisatorsatz") and _check_named(d, "A6 [E] Markenzahl")),
    ("step identification from the clock sectors", ("I_step",), "z175_ergebnis.json", "excluded",
     "the clock is the identity on the doubled ladder, so no condition in clock values sees the shell factor, and a "
     "single spinor index lies only in representations on which the clock is +-i",
     lambda d: _check_named(d, "C2 [E, negativ]") and _check_named(d, "C3 [E, negativ]")),
    ("corner of the null cone (arXiv:2611.34924) fixes the gravity items", ("P_G", "P_curv", "graviton"),
     "z177_ergebnis.json", "excluded",
     "G and Lambda are inputs there and the discrete area spectrum is an existence statement; the paper gives the "
     "dictionary seam = corner (the clock is the boost by pi, 4 pi is the spinor period) and the test separate against "
     "cyclic gluing",
     lambda d: _check_named(d, "Abgleich 2 pi gegen 4 pi") and _check_named(d, "Wörterbuch im endlichen Modell")),
    ("graviton as an excitation of one glued object", ("graviton",), "z177_ergebnis.json", "excluded",
     "the spin table R = i^(2-2h) has clock-neutral states only at odd h; h = 2 with R = 1 exists only as a state of two "
     "objects, which is the double-copy route",
     lambda d: _check_named(d, "Erweiterte Spintabelle")),
    ("gauge completion of the pair clock from superselection, the GNS rule and A*", ("K19",), "z178_ergebnis.json",
     "derived",
     "over the clock-neutral algebra lie five of the 27 completions (three bosonic, two spin forms); the GNS space of the "
     "state over the neutral observables is the U-invariant subspace, and S and T leave one completion without charged "
     "field. K19 rests on P_std; the counter-models are E8 (x) E8-bar with a selection rule (fails the GNS rule) and a "
     "spin form (fails T)",
     lambda d: _check_named(d, "A1 [E] Die uhrneutrale Algebra") and _check_named(d, "A3 [E] Modulgleichung")
     and _check_named(d, "B1 [E, drei Modelle] Standardform") and _check_named(d, "B2 [E] Gegenmodelle")),
    ("one parent: the state of A* on the torus, its modulator and the opposite arcs of the marks", ("K20",),
     "z180_ergebnis.json", "derived",
     "the reading is no longer needed for clock, defect and sheet exchange from T = rho^(-i) (z178, z176 C4), and the "
     "split seam | partner is the pair of opposite arcs of the four marks; the sector that carries the monodromy is the "
     "identification of the clock. K20 rests on P_std and P_zentral",
     lambda d: _check_named(d, "B1 [N] Zweipunktfunktion") and _check_named(d, "B2 [E] Modulare Flüsse")),
    ("evenness of the number of seams from the associator class of the clock", ("P_paar",), "z178_ergebnis.json", "holds",
     "the anomaly of the one-sided clock is the class 2 in Z4 of an explicit abelian 3-cocycle and is additive, so N is "
     "even for every N; 'exactly two' is numerical for N <= 24, and the one-sided summed glue stays the premise",
     lambda d: _check_named(d, "C1 [E] Expliziter abelscher 3-Kozykel") and _check_named(d, "E1 [N, 80 Stellen]")),
    ("clock = deck parity of the sections extended to the seam (P_deck)", ("P_zentral",), "z179_ergebnis.json", "holds",
     "the deck parity (+1 on ten, -1 on six parents) has exactly two extensions to E8, mutually inverse and of order 4; "
     "they commute with the base changes of both deck halves and are the only elements of order 4 that do. P_zentral is "
     "equivalent to the identification of the clock with the deck parity; the count stays",
     lambda d: _check_named(d, "A2 [E] Eindeutigkeit der Hebung") and _check_named(d, "B1 [E] Fortsetzung auf E8")
     and _check_named(d, "C1 [E] Natürlichkeit") and _check_named(d, "D1 [E] Umkehrung")),
    ("spacetime and order from the relative modularity of the mark arcs", ("P_ord", "I_step", "P_std"),
     "z180_ergebnis.json", "excluded",
     "in the state of A* the four arcs have only two modular flows, the commuting translations of the torus, and no "
     "half-sided modular inclusion; S reverses the orientation of the marks and fixes the state, so no time direction "
     "follows. Translation, positivity, adjacency and the unit step come from the same inclusions only in the vacuum of "
     "the base",
     lambda d: _check_named(d, "B2 [E] Modulare Flüsse") and _check_named(d, "C2 [N] Zustand von A*")
     and _check_named(d, "D1 [E] Zeitrichtung")),
    ("interaction from the gluing of three and four seams", ("P_ord", "P_std"), "z181_ergebnis.json", "excluded",
     "associator and braiding of the clock are classified by H^3_ab = Z8 (Z4 for the parity) without a continuous "
     "parameter; the gate family is associative and Yang-Baxter for every phase on integer charges, and under P_uhr the "
     "phase lies on a grid of eight values that still contains zero. Topology fixes discrete phases, neither the point "
     "of the grid nor an order nor a coupling",
     lambda d: _check_named(d, "A1 [E] Abelsche 3-Kozykel") and _check_named(d, "B1 [E] Gatterfamilie")
     and _check_named(d, "C1 [E, N = 3 und 4")),
    ("graviton as a collective excitation of several seams", ("graviton", "P_curv"), "z182_ergebnis.json", "excluded",
     "the product of two massless representations of su(2,2) starts one unit above every massless representation, so "
     "two-seam states are conserved currents and the spin-2 channel is a continuum; the elementary h = +-2 state has no "
     "bilinear stress tensor (Weinberg-Witten as a statement on characters), so universal coupling needs a gauge "
     "principle",
     lambda d: _check_named(d, "A1 [E, Potenzreihen") and _check_named(d, "A2 [E] Kein masseloser Pol")
     and _check_named(d, "B1 [E] Weinberg-Witten")),
    ("second gravity relation from the entropy of the state of A*", ("P_G",), "z182_ergebnis.json", "holds",
     "the spectrum of one seam induces no Einstein term at fixed tau, and its entropy is exactly log 12 because <T> = 0; "
     "under the area rule G/a^2 = 1/(4 s_1), with three physically different values of s_1. The free real number becomes "
     "a discrete ambiguity under a named premise; P_G stays",
     lambda d: _check_named(d, "C1 [N, 40 Stellen]") and _check_named(d, "C2 [E/N] Entropie")
     and _check_named(d, "D2 [C, Mehrdeutigkeit]")),
    # H6h52 (z183): the attack on P_std. Its pair-clock part is restated as a locality statement and stays a premise.
    ("pair-clock part of P_std as locality: no one-sided current mixes the deck halves", ("P_std",),
     "z183_ergebnis.json", "holds",
     "the parents see the deck parity as monodromy exactly around the eight sectors with h in odd eighths; a completion "
     "contains such a mark field exactly if the mixing current (10, 6) is local on neither side (all 27), and under A* "
     "these are the four graphs with one spectrum, the gauge sums of z178. Superselection of the pair clock with the "
     "GNS rule, a mark field, no mixing current and 'the 16 is a local field of the pair' are equivalent under A*; the "
     "standard form is not needed. The carrier of the response power 5/4 (z170) is local only there. P_std stays",
     lambda d: _check_named(d, "A2 [E] Lemma") and _check_named(d, "A3 [E] Satz")
     and _check_named(d, "A4 [E] Gleichwertig") and _check_named(d, "C3 [E] Gelöster alpha-Block")),
    ("P_std from A*, time neutrality or the covering of the four marks", ("P_std",), "z183_ergebnis.json", "excluded",
     "under A* there are ten completions in three spectra; the six without gauge sum (four of type E8 (x) E8-bar, two of "
     "type so(16)) satisfy A* with the deck parity as a symmetry. Time neutrality removes two of the ten and no "
     "spectrum; orientation reversal of the section bases removes E8 (x) E8-bar and leaves so(16). 'No distinguished "
     "cycles' admits the relative spin structure zero on its own, and every cycle of the seam torus encircles two marks, "
     "so odd sections are single-valued there: the covering gives the joint sum, not the gauge sum",
     lambda d: _check_named(d, "B1 [E] Die zehn") and _check_named(d, "B3 [E] Orientierungswechsel")
     and _check_named(d, "B5 [E] Zeitneutralität") and _check_named(d, "C1 [E/N] A allein")
     and _check_named(d, "C2 [N] Überlagerungsgeometrie")),
]


# H6h51: ambiguities recorded by exploration results: sets of physically different models that satisfy the same
# conditions. An ambiguity is reported and never resolved by choice.
RECORDED_AMBIGUITIES = [
    ("value of G/a^2 under the area rule", ("P_G",), "z182_ergebnis.json", 3,
     "entropy per cell: one seam (log 12), a gauged pair (half of log(80 + 48 sqrt(2))) or one seam minus the "
     "topological term log 4",
     lambda d: _check_named(d, "D2 [C, Mehrdeutigkeit]")),
    ("gate phase under P_uhr", ("P_std",), "z181_ergebnis.json", 8,
     "pentagon and hexagon on the clock values leave the eight classes of H^3_ab(Z4), among them zero (no gate) and "
     "P_viertel",
     lambda d: _check_named(d, "B1 [E] Gatterfamilie")),
    ("orientation of the circle through the four marks", ("P_ord",), "z180_ergebnis.json", 2,
     "the two orientations are exchanged by S, which fixes the state of A*; a time direction or an order of the arcs is "
     "a bit that breaks S",
     lambda d: _check_named(d, "D1 [E] Zeitrichtung")),
    ("natural Z4 of the geometry with the fixed algebra of H6", ("P_zentral",), "z180_ergebnis.json", 1,
     "of the quarter turn of the marks, the half turn of the seam and the deck parity of the sections only the last has "
     "the fixed algebra so(10) + so(6)",
     lambda d: _check_named(d, "D2 [E] Drei geometrische Kandidaten")),
    ("spectrum of the pair of seam and partner under A*", ("P_std",), "z183_ergebnis.json", 3,
     "the gauge sum with all 16 clock sectors (Z(i) = 80 + 48 sqrt(2)), the so(16) form and E8 (x) E8-bar (both 144 at "
     "tau = i); the last two have the same partition function on the whole imaginary axis and differ in spin only",
     lambda d: _check_named(d, "B1 [E] Die zehn") and _check_named(d, "B2 [N, 40 Stellen]")),
]


def counter_models(root=None):
    """Automatic counter-models. Part one enumerates OriginLab's own universe of the 27 left-right completions of the 16
    clock sectors and counts, for nested sets of conditions, the physically different models that satisfy all of them.
    Part two lists the ambiguities recorded by exploration results. More than one model is an ambiguity."""
    import sys as _sys
    here = str(_Path(__file__).resolve().parent)
    if here not in _sys.path:
        _sys.path.insert(0, here)
    import clock_sector_search as cs
    root = EXPLORATION if root is None else _Path(root)
    comp = [frozenset(x) for x in cs.isotropic_subgroups() if len(x) == 16]
    glue = [(k, k) for k in range(4)]
    untwisted = frozenset((a, b) for a in glue for b in glue)
    neutral = frozenset((a, a) for a in glue)
    charge = lambda g: (g[0][0] - g[1][0]) % 4
    flux = lambda g: (g[0][0] - g[0][1]) % 4

    def signature(L):
        return ("bosonic" if all(cs.Q_lr(g) == 0 for g in L) else "spin", len(L & untwisted),
                tuple(sorted({charge(g) for g in L & untwisted})),
                tuple(sorted({(charge(g), str(cs.Q_lr(g))) for g in L if flux(g) == 1})))

    # z183: conditions without a mirror bit. The parents see the deck parity around a sector (x, y) exactly if x + y is
    # odd; a reflection of one parent acts as x -> -x (or y -> -y) on both sides; the spectrum is the multiset of weights.
    odd = lambda a: (a[0] + a[1]) % 2 == 1
    image = lambda L, fx, fy: frozenset((((fx * g[0][0]) % 4, (fy * g[0][1]) % 4), ((fx * g[1][0]) % 4, (fy * g[1][1]) % 4))
                                        for g in L)
    spectrum = lambda L: tuple(sorted((str(cs.q_form(g[0])), str(cs.q_form(g[1]))) for g in L))
    cond = {"S": lambda L: True,
            "T": lambda L: all(cs.Q_lr(g) == 0 for g in L),
            "superselection": lambda L: all(charge(g) == 0 for g in L & untwisted),
            "GNS": lambda L: L & untwisted == neutral,
            "mark": lambda L: any(odd(g[0]) and odd(g[1]) for g in L),
            "reflection": lambda L: image(L, 3, 1) == L and image(L, 1, 3) == L,
            "exchange": lambda L: frozenset((g[1], g[0]) for g in L) == L}
    levels = [("base changes S alone (fixed spin structure)", ["S"], "signature"),
              ("A*: S and T", ["S", "T"], "signature"),
              ("A* and no local field with pair-clock charge", ["S", "T", "superselection"], "signature"),
              ("A*, superselection and the GNS rule", ["S", "T", "superselection", "GNS"], "signature"),
              ("S, superselection and the GNS rule, without T", ["S", "superselection", "GNS"], "signature"),
              ("A* and time neutrality (exchange of seam and partner)", ["S", "T", "exchange"], "spectrum"),
              ("A* and orientation reversal of the section bases", ["S", "T", "reflection"], "spectrum"),
              ("A* and a local field with the deck parity as monodromy of the parents", ["S", "T", "mark"], "spectrum")]
    reports = []
    for name, cs_, by in levels:
        models = [L for L in comp if all(cond[c](L) for c in cs_)]
        sigs = sorted({signature(L) for L in models})
        spectra = len({spectrum(L) for L in models})
        n = len(sigs) if by == "signature" else spectra
        reports.append({"conditions": name, "models": len(models), "physically_different": len(sigs), "spectra": spectra,
                        "counted_by": by,
                        "signatures (statistics, untwisted sectors, charges, (charge, spin) at clock flux 1)":
                            [[a, b, list(c), [list(x) for x in d]] for a, b, c, d in sigs],
                        "verdict": "unique" if n == 1 else ("empty" if not n else "ambiguous")})
    recorded = []
    for name, items, fname, n, meaning, pred in RECORDED_AMBIGUITIES:
        path = root / fname
        ok = path.exists() and pred(json.loads(path.read_text()))
        recorded.append({"question": name, "items": list(items), "source": fname, "models": n if ok else None,
                         "verdict": "pending" if not ok else ("unique" if n == 1 else "ambiguous"), "meaning": meaning})
    return {"completions_of_the_clock_sectors": reports, "recorded": recorded,
            "ambiguous": [r["question"] for r in recorded if r["verdict"] == "ambiguous"]}


def unification_candidates(root=None):
    """Verdicts of the recorded unification candidates and the open items that would remain under the readings."""
    root = EXPLORATION if root is None else _Path(root)
    rows = []
    for name, items, fname, verdict, meaning, pred in UNIFICATION_CANDIDATES:
        path = root / fname
        ok = path.exists() and pred(json.loads(path.read_text()))
        rows.append({"candidate": name, "would_close": list(items), "source": fname,
                     "verdict": verdict if ok else "pending", "meaning": meaning})
    d = complete_solution_distance(root)
    open_items = set(d["open_items"]) | set(d["bench_open"])
    read = {x for r in rows if r["verdict"] == "reading" for x in r["would_close"]}
    counts = {}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    derived = sorted({x for r in rows if r["verdict"] == "derived" for x in r["would_close"]} & set(d["bench_derived"]))
    return {"candidates": rows, "counts": counts, "closed_by_derivation": [],
            "derived_under_premises": {k: d["bench_derived"][k]["rests_on"] for k in derived},
            "open_total": d["open_total"],
            "open_under_readings": sorted(open_items - read), "open_total_under_readings": len(open_items - read)}


def complete_solution_distance(root=None):
    """Distance of the admissible coherent scenarios to the complete solution, the best scenarios and the open items."""
    root = EXPLORATION if root is None else _Path(root)
    sc, _ = scenario_candidates(root)
    conv = {}
    for name, fname, meaning, pred in CONVENTIONS:
        path = root / fname
        conv[name] = {"source": fname, "meaning": meaning, "holds": path.exists() and pred(json.loads(path.read_text()))}
    not_counted = {AXIOM} | {k for k, v in conv.items() if v["holds"]}
    rows = []
    for o in sc:
        if not {"admissible", "coherent"} <= o["clauses"]:
            continue
        missing = [c for c in SCENARIO_POOL if c not in o["clauses"]]
        extra = o["premises"] - not_counted
        rows.append(((len(missing), len(extra & FREE_PARAMETER_PREMISES), len(extra)), o["label"], missing, sorted(extra)))
    least = min(r[0] for r in rows)
    best = [r for r in rows if r[0] == least]
    items = sorted({x for r in best for x in r[2] + r[3]})
    derived = {}
    for row, prem, sources in BENCH_DERIVED:
        paths = [(root / fname, pred) for fname, _, pred in sources]
        if set(prem) <= set(items) and all(p_.exists() and pred(json.loads(p_.read_text())) for p_, pred in paths):
            derived[row] = {"rests_on": list(prem), "inventory": [AXIOM, "GNS rule"],
                            "sources": [fname for fname, _, _ in sources], "meaning": sources[0][1]}
    bench = {k: v for k, v in BENCH_OPEN.items() if k not in derived}
    by_question = {"origin": sorted(bench)}
    for x in items:
        by_question.setdefault(FRONTIER.get(x, ("unnamed", ""))[0], []).append(x)
    return {"target": "all clauses, axiom alone, no free number, bench rows K19 and K20 closed", "scenarios": len(rows),
            "distance": {"missing_clauses": least[0], "free_numbers": least[1], "premises_beyond_axiom": least[2]},
            "complete": least == (0, 0, 0) and not bench, "best": [r[1] for r in best],
            "open_items": {x: {"question": FRONTIER.get(x, ("unnamed", ""))[0], "closes_with": FRONTIER.get(x, ("", ""))[1],
                               "recorded": PREMISE_MEANING.get(x, "clause of the scenario pool")} for x in items},
            "bench_open": {k: {"question": q, "closes_with": t} for k, (q, t) in bench.items()},
            "bench_derived": derived,
            "open_total": len(items) + len(bench), "by_question": by_question, "conventions": conv}


class _Survivors:
    """Sorted distinct labels of the objects that satisfy a rule, read lazily from the rule's bit mask."""

    def __init__(self, masks, labels):
        self._masks, self._labels = masks, labels

    def labels(self, m, limit=None):
        out = []
        while m and (limit is None or len(out) < limit):
            low = m & -m
            name = self._labels[low.bit_length() - 1]
            if not out or out[-1] != name:
                out.append(name)
            m ^= low
        return out

    def __getitem__(self, rule):
        return self.labels(self._masks[frozenset(rule)])


def _rules_report(objects, pool, label):
    """All rules over the clause pool, as genesis_rules.search with constant size, computed on bit masks: bit p is the
    p-th object in label order, so the set bits of a rule's mask list its survivors sorted, duplicates adjacent."""
    labels = [label(o) for o in objects]
    order = sorted(range(len(objects)), key=lambda i: labels[i])
    lab = [labels[i] for i in order]
    unique_labels = len(set(lab)) == len(lab)
    ext = {c: 0 for c in pool}
    for p, i in enumerate(order):
        for c in objects[i]["clauses"]:
            if c in ext:
                ext[c] |= 1 << p
    rules = powerset(pool)
    masks = {}
    for r in rules:
        m = (1 << len(objects)) - 1
        for c in r:
            m &= ext[c]
        masks[r] = m
    hits = _Survivors(masks, lab)
    rows = []
    for r in rules:
        m = masks[r]
        if not m:
            rows.append({"rule": sorted(r), "verdict": "EMPTY", "selected": [], "size": None})
            continue
        count = m.bit_count() if unique_labels else len(hits.labels(m))
        rows.append({"rule": sorted(r), "verdict": "UNIQUE" if count == 1 else "AMBIGUOUS",
                     "selected": hits.labels(m, 12), "selected_count": count, "size": 0})
    by_rule = {frozenset(row["rule"]): row for row in rows}
    selections = {}
    for row in rows:
        if row["verdict"] == "UNIQUE":
            selections.setdefault(row["selected"][0], []).append(frozenset(row["rule"]))
    summary = []
    for chosen, makers in selections.items():
        minimal = [m for m in makers if not any(o < m for o in makers)]
        bearing = {}
        for m in minimal:
            for clause in m:
                other = by_rule[m - {clause}]
                bearing[f"{'+'.join(sorted(m))} without {clause}"] = (other["verdict"], other["selected"][:4], other["size"])
        summary.append({"selected": chosen, "shortest_rule_length": min(len(m) for m in makers),
                        "minimal_rules": [sorted(m) for m in minimal], "ablation": bearing})
    res = {"objects": len(objects), "clauses": list(pool), "rules": len(rules),
           "verdicts": {v: sum(row["verdict"] == v for row in rows) for v in ("EMPTY", "UNIQUE", "AMBIGUOUS")},
           "unique_selections": sorted(summary, key=lambda s: (s["shortest_rule_length"], str(s["selected"])))}
    order_key = lambda r: (len(r), sorted(r))
    # monotone in the rule: maximal non-empty = every one-clause extension empty, minimal empty = every one-clause
    # reduction non-empty
    res["maximal_rules"] = [{"rule": sorted(r), "survivors": hits[r]} for r in sorted(rules, key=order_key)
                            if masks[r] and not any(masks[r | {c}] for c in pool if c not in r)]
    res["minimal_empty_rules"] = [sorted(r) for r in sorted(rules, key=order_key)
                                  if not masks[r] and all(masks[r - {c}] for c in r)]
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
                 "origin_reduced", "marks_placed", "spin_statistics", "chiral_matter", "interacting", "modular_statistics"]


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
    matter = base + ["one_field", "spin_statistics", "interacting"]
    res["questions"] = {"complete": _question(hits_c, SCENARIO_POOL), "+".join(base): _question(hits_c, base),
                        "+".join(matter): _question(hits_c, matter),
                        "+".join(matter + ["chiral_matter"]): _question(hits_c, matter + ["chiral_matter"]),
                        "admissible+chiral_matter": _question(hits_c, ["admissible", "chiral_matter"]),
                        "admissible+modular_statistics": _question(hits_c, ["admissible", "modular_statistics"]),
                        "admissible+graviton": _question(hits_c, ["admissible", "graviton"])}
    best = max(len(o["clauses"]) for o in sc)
    res["premise_closure"] = _premise_closure(sc, matter + ["chiral_matter"])
    res["derived"] = sc_info | {"most_clauses": best,
                                "best_scenarios": sorted(o["label"] for o in sc if len(o["clauses"]) == best)[:40],
                                "values_in_base_survivors": {a: sorted({o["scenario"][a] for o in sc if set(base) <= o["clauses"]})
                                                             for a in SCENARIO_AXES}}
    res["bound"] = ("product of eight axes (spinor index, role of the marks, time, modes, gravity, origin, statistics, "
                    "interaction) = "
                    f"{len(sc)} scenarios; clauses from OriginLab's R orders and position spaces and from the recorded "
                    "results named in EVIDENCE; pending evidence gives no clause")
    out["universes"]["scenarios"] = res
    out["complete_solution"] = complete_solution_distance(root)
    out["unification"] = unification_candidates(root)
    out["counter_models"] = counter_models(root)
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
