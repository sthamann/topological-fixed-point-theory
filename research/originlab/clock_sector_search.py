"""Clock-sector search (H6h50, z170): finite universes behind the four clock sectors of (E8)_1.

Guideline "Rekonstruktion statt Auswahl": every clause is a computed property of its candidate; a rule is UNIQUE only
if exactly one label satisfies it (engine: `genesis_rules.search` with constant size). Nothing here is a physical
solution; the search states which native properties pin an object and which freedom is left.

  clocks       every conjugacy class of elements of order 2 to MAX_ORDER of E8 (Kac coordinates s, coprime, with
               s_0 + sum a_i s_i = N on the affine diagram), with the multiplicities of R on the 248, Tr R (exact, as a
               cyclotomic integer), the fixed root system, the lowest weight v.v/2 of the R-twisted sector and the
               anomaly of the cyclic group on the chiral factor.
  completions  every self-orthogonal subgroup H of order 16 of A (+) A-bar, A = Z4 x Z4 the discriminant group of the
               fixed lattice D5 + D3 with q(x, y) = (5 x^2 + 3 y^2)/8 mod 1 (z144). H is a left-right completion
               (sum over (a, a') in H of V_a (x) conj V_a'); h - h' is in Z/2 on H. Clauses: bosonic or spin,
               factorised or graph of an isometry, mirror symmetric, clock sectors present, scalar dimensions.
  responses    exact power laws (1 - q)^(-p) offered by the objects above for the seam response q (1 - q)^(-5/4) of the
               ledger (EM.WARD.01): coherent kernels of the twisted ground states, scalar two-point functions of the
               completions, the positive quartic measure of analysis 22. Clauses compare the power and the insertion.
"""
from fractions import Fraction as Fr
from itertools import product

import sympy as sp

from genesis_rules import search

H = Fr(1, 2)
# Bourbaki simple roots of E8 and the marks of the highest root
SIMPLE = [tuple(Fr(x, 2) for x in (1, -1, -1, -1, -1, -1, -1, 1)),
          (1, 1, 0, 0, 0, 0, 0, 0), (-1, 1, 0, 0, 0, 0, 0, 0), (0, -1, 1, 0, 0, 0, 0, 0), (0, 0, -1, 1, 0, 0, 0, 0),
          (0, 0, 0, -1, 1, 0, 0, 0), (0, 0, 0, 0, -1, 1, 0, 0), (0, 0, 0, 0, 0, -1, 1, 0)]
SIMPLE = [tuple(Fr(x) for x in r) for r in SIMPLE]
MARKS = (2, 3, 4, 6, 5, 4, 3, 2)
MAX_ORDER = 8
CLOCK_POOL = ["order_4", "order_at_most_4", "trace_zero", "square_is_fermion_parity", "semisimple_fixed",
              "anomaly_free", "fixed_has_D5"]
COMPLETION_POOL = ["bosonic", "spin", "factorised", "graph", "mirror_symmetric", "clock_sectors_present",
                   "scalar_16", "scalar_4", "holomorphic_fermions", "charge_conjugate", "diagonal"]
RESPONSE_POOL = ["power_5/4", "positive", "native", "insertion_matches", "no_extra_structure"]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def e8_roots():
    out = []
    for i in range(8):
        for j in range(i + 1, 8):
            for a in (1, -1):
                for b in (1, -1):
                    v = [Fr(0)] * 8
                    v[i], v[j] = Fr(a), Fr(b)
                    out.append(tuple(v))
    for sg in product((1, -1), repeat=8):
        if sum(1 for x in sg if x < 0) % 2 == 0:
            out.append(tuple(Fr(x, 2) for x in sg))
    return out


ROOTS = e8_roots()


def coweights():
    """Fundamental coweights w_i with w_i . alpha_j = delta_ij (E8 is simply laced and self-dual)."""
    M = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in r] for r in SIMPLE])
    W = M.inv().T
    return [tuple(Fr(int(x.p), int(x.q)) for x in W.row(i)) for i in range(8)]


def _components(roots):
    """Irreducible components of a root subsystem as sorted names, from (number of roots, rank)."""
    left, names = set(roots), []
    while left:
        comp, todo = set(), [next(iter(left))]
        while todo:
            r = todo.pop()
            if r in comp:
                continue
            comp.add(r)
            todo.extend(s for s in left if s not in comp and dot(r, s) != 0)
        left -= comp
        rank = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in r] for r in comp]).rank()
        n = len(comp)
        name = {(72, 6): "E6", (126, 7): "E7", (240, 8): "E8"}.get((n, rank))
        if name is None:
            name = f"A{rank}" if n == rank * (rank + 1) else f"D{rank}" if n == 2 * rank * (rank - 1) else f"?{n}/{rank}"
        names.append(name)
    return sorted(names)


def _trace_zero(cls):
    """Tr R = sum n_k zeta^k = 0 exactly: the polynomial sum n_k x^k is divisible by the cyclotomic polynomial."""
    x = sp.Symbol("x")
    n = len(cls)
    return sp.rem(sum(c * x ** k for k, c in enumerate(cls)), sp.cyclotomic_poly(n, x), x) == 0


def _kac(N, i=0, left=None):
    """Non-negative s_1..s_8 with sum a_i s_i <= N."""
    left = N if left is None else left
    if i == 8:
        yield ()
        return
    for k in range(left // MARKS[i] + 1):
        for tail in _kac(N, i + 1, left - k * MARKS[i]):
            yield (k,) + tail


def clock_candidates(max_order=MAX_ORDER):
    """All classes of order 2..max_order in E8. The shift vector v = (1/N) sum s_i w_i lies in the closed alcove, so
    v.v/2 is the lowest conformal weight of the R-twisted sector of the lattice theory; the cyclic group is free of
    anomaly on the chiral factor exactly if N v.v/2 is an integer."""
    from math import gcd
    W = coweights()
    out = []
    for N in range(2, max_order + 1):
        for s in _kac(N):
            rest = N - sum(a * b for a, b in zip(MARKS, s))
            g = rest
            for x in s:
                g = gcd(g, x)
            if g != 1:
                continue                                                     # order smaller than N
            v = tuple(sum(Fr(s[i], N) * W[i][k] for i in range(8)) for k in range(8))
            cls = [8] + [0] * (N - 1)
            for r in ROOTS:
                cls[int(N * dot(v, r)) % N] += 1
            fixed = [r for r in ROOTS if dot(v, r).denominator == 1]
            comps = _components(fixed)
            rank = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in r] for r in fixed]).rank() if fixed else 0
            h = dot(v, v) / 2
            cl = {c for c, ok in (
                ("order_4", N == 4), ("order_at_most_4", N <= 4), ("trace_zero", _trace_zero(cls)),
                # for N = 4: R^2 has the multiplicities (120, 128) of the fermion parity (D8 involution)
                ("square_is_fermion_parity", N == 4 and (cls[0] + cls[2], cls[1] + cls[3]) == (120, 128)),
                ("semisimple_fixed", rank == 8), ("anomaly_free", (N * h).denominator == 1),
                ("fixed_has_D5", "D5" in comps)) if ok}
            name = ("+".join(comps) if comps else "torus") + ("" if rank == 8 else f"+U(1)^{8 - rank}")
            out.append({"label": f"Z{N}: {name} {(rest,) + s}", "order": N, "fixed": name, "kac": (rest,) + s,
                        "multiplicities": cls, "twist_ground": h, "anomaly_type": int(N * N * h) % N, "clauses": cl})
    return out


# ---------------------------------------------------------------------------------------------- completions
A = [(x, y) for x in range(4) for y in range(4)]
q_form = lambda a: Fr(5 * a[0] ** 2 + 3 * a[1] ** 2, 8) % 1
b_form = lambda a, c: Fr(5 * a[0] * c[0] + 3 * a[1] * c[1], 4) % 1
G = [(a, c) for a in A for c in A]
Q_lr = lambda g: (q_form(g[0]) - q_form(g[1])) % 1
B_lr = lambda g, k: (b_form(g[0], k[0]) - b_form(g[1], k[1])) % 1
g_add = lambda g, k: (((g[0][0] + k[0][0]) % 4, (g[0][1] + k[0][1]) % 4), ((g[1][0] + k[1][0]) % 4, (g[1][1] + k[1][1]) % 4))
E8_GLUE = [frozenset((k % 4, (s * k) % 4) for k in range(4)) for s in (1, 3)]      # the two Lagrangians of z144


def _span(sub, g):
    out = set(sub)
    cur = g
    while cur not in out:
        out |= {g_add(h, cur) for h in sub}
        cur = g_add(cur, g)
    return frozenset(out)


def isotropic_subgroups():
    """All subgroups of A (+) A-bar on which the bilinear form vanishes."""
    zero = ((0, 0), (0, 0))
    self_orth = [g for g in G if B_lr(g, g) == 0]
    seen, todo = {frozenset([zero])}, [frozenset([zero])]
    while todo:
        sub = todo.pop()
        for g in self_orth:
            if g in sub or any(B_lr(g, h) != 0 for h in sub):
                continue
            new = _span(sub, g)
            if new not in seen:
                seen.add(new)
                todo.append(new)
    return seen


def completion_candidates():
    """Self-orthogonal subgroups of order 16: the left-right completions of the sixteen sectors."""
    out = []
    for sub in sorted(isotropic_subgroups(), key=lambda s: sorted(s)):
        if len(sub) != 16:
            continue
        left = frozenset(g[0] for g in sub if g[1] == (0, 0))
        right = frozenset(g[1] for g in sub if g[0] == (0, 0))
        graph = len(left) == 1 and len(right) == 1
        fact = len(left) == 4 and len(right) == 4
        spins = {Q_lr(g) for g in sub}
        mp_ = dict(sub) if graph else None
        iso = None
        if graph:
            iso = (mp_[(1, 0)], mp_[(0, 1)])
        kind = ("graph " + str(iso)) if graph else ("factorised" if fact else f"mixed {len(left)}x{len(right)}")
        cl = {c for c, ok in (
            ("bosonic", spins == {Fr(0)}), ("spin", Fr(1, 2) in spins), ("factorised", fact), ("graph", graph),
            ("mirror_symmetric", all((g[1], g[0]) in sub for g in sub)),
            ("clock_sectors_present", any(not any(g[0] in L for L in E8_GLUE) for g in sub)),
            ("scalar_16", any(g in sub for g in (((1, 0), (1, 0)), ((1, 0), (3, 0))))),
            ("scalar_4", any(g in sub for g in (((0, 1), (0, 1)), ((0, 1), (0, 3))))),
            ("holomorphic_fermions", any(g[1] == (0, 0) and Q_lr(g) == H for g in sub)),
            ("charge_conjugate", graph and iso == ((3, 0), (0, 3))), ("diagonal", graph and iso == ((1, 0), (0, 1)))) if ok}
        label = f"{len(out):02d} " + kind + " | " + ("bosonic" if spins == {Fr(0)} else "spin") + " | left " \
            + str(sorted(left)) + " right " + str(sorted(right))
        out.append({"label": label, "kind": kind, "isometry": iso, "left": sorted(left), "right": sorted(right),
                    "clauses": cl})
    return out


# ---------------------------------------------------------------------------------------------- responses
def response_candidates():
    """Exact power laws (1 - q)^(-p) against the seam response q (1 - q)^(-5/4). A lowest-weight state of weight h has
    the coherent kernel sum_n |L_-1^n h|^2 q^n/(n!)^2 = (1 - q)^(-2h) (norms n! (2h)_n, positive for h > 0)."""
    q = sp.Symbol("q")
    target = q * (1 - q) ** sp.Rational(-5, 4)
    rows = [
        # name, power p, insertion F as a function of the kernel K = (1-q)^(-p), native, extra structure
        ("coherent kernel of the 16 (R sector, h = 5/8), F = q K", sp.Rational(5, 4), lambda K: q * K, True, False),
        ("coherent kernel of the 4 (R sector, h = 3/8), F = q K", sp.Rational(3, 4), lambda K: q * K, True, False),
        ("coherent kernel of the vector (R^2 sector, h = 1/2), F = q K", sp.Integer(1), lambda K: q * K, True, False),
        ("coherent kernel of weight h = 1/8, F = 4 q dK/dq", sp.Rational(1, 4), lambda K: 4 * q * sp.diff(K, q), False, True),
        ("coherent kernel of the 16, F = q d log K/dq", sp.Rational(5, 4), lambda K: q * sp.diff(sp.log(K), q), True, False),
        ("quartic measure on five real directions (analysis 22), F = q K", sp.Rational(5, 4), lambda K: q * K, False, True),
        ("quartic measure, F = q d log K/dq", sp.Rational(5, 4), lambda K: q * sp.diff(sp.log(K), q), False, True),
        ("clock character without vacuum factor", None, None, True, False),
    ]
    out = []
    for name, p, ins, native, extra in rows:
        if p is None:
            cl = {"native", "no_extra_structure"}                           # starts at q^2: not a power law in 1 - q
            out.append({"label": name, "power": None, "clauses": cl})
            continue
        K = (1 - q) ** (-p)
        norms_positive = all(sp.rf(p, n) > 0 for n in range(8))
        match = sp.simplify(ins(K) - target) == 0
        cl = {c for c, ok in (("power_5/4", p == sp.Rational(5, 4)), ("positive", norms_positive), ("native", native),
                              ("insertion_matches", match), ("no_extra_structure", not extra)) if ok}
        out.append({"label": name, "power": str(p), "clauses": cl})
    return out


def run_clock_sector_search():
    clocks = clock_candidates()
    comps = completion_candidates()
    resp = response_candidates()
    uni = {}
    for name, objs, pool in (("clocks", clocks, CLOCK_POOL), ("completions", comps, COMPLETION_POOL),
                             ("responses", resp, RESPONSE_POOL)):
        res = search(objs, pool, lambda o: 0, lambda o: o["label"])
        uni[name] = {"candidates": len(objs), "labels": sorted({o["label"] for o in objs}),
                     "clauses": {o["label"]: sorted(o["clauses"]) for o in objs}, "search": res}
    h6 = [o for o in clocks if o["order"] == 4 and o["fixed"] == "A3+D5"]
    uni["clocks"]["derived"] = {o["label"]: {"multiplicities": o["multiplicities"], "twist_ground": str(o["twist_ground"]),
                                              "anomaly_type": o["anomaly_type"]} for o in clocks}
    uni["clocks"]["trace_zero_classes"] = [o["label"] for o in clocks if "trace_zero" in o["clauses"]]
    uni["clocks"]["classes_per_order"] = {N: sum(o["order"] == N for o in clocks) for N in range(2, MAX_ORDER + 1)}
    uni["clocks"]["h6_clock"] = h6[0]["label"] if len(h6) == 1 else None
    both = [o["label"] for o in comps if {"factorised", "scalar_16"} <= o["clauses"]]
    uni["completions"]["e8_currents_and_local_16"] = both
    uni["completions"]["counts"] = {c: sum(c in o["clauses"] for o in comps) for c in COMPLETION_POOL}
    return {"status": "CLOCK_SECTORS_SEARCHED", "physical_solution": False, "universes": uni}
