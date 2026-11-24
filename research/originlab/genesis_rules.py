"""Genesis rule search: which simple rule, in which ontology, singles out one object?

Each ontology is a finite universe of small mathematical structures, enumerated exhaustively up to a
stated bound, with a size and a pool of clauses (predicates computed from the structure). A rule is a
subset of the clause pool. The search generates every rule, keeps the solutions of minimal size and
classifies the rule as EMPTY, UNIQUE or AMBIGUOUS. For every unique selection it reports the smallest
rules that make it and which clauses are load bearing (dropping one changes the selection).

Nothing here uses a TFPT number as input. Structures of the TFPT record are named only afterwards, in
`RECORD`, to state which selections meet the record. A selection that meets the record shows that a
short rule in that universe has this object as its smallest solution. It does not show that nature
obeys the rule, and it derives no physical map. No ledger status, no Atlas entry, no promotion.
"""
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations, product
from math import comb, factorial, gcd, lcm
import argparse
import json

import sympy

SCOPE = ("exhaustive finite enumeration inside stated bounds; a unique minimal solution is a statement about "
         "the universe and its clause pool, not a derivation of physics")


def powerset(pool):
    return [frozenset(c) for k in range(len(pool)+1) for c in combinations(pool, k)]


def search(objects, pool, size, label):
    """All rules over the clause pool. objects: dicts with 'clauses' (set), size(obj) sortable, label(obj) hashable."""
    rules = []
    for rule in powerset(pool):
        hits = [o for o in objects if rule <= o["clauses"]]
        if not hits:
            rules.append({"rule": sorted(rule), "verdict": "EMPTY", "selected": [], "size": None,
                          "matching_objects": 0, "matching_labels": 0, "minimal_objects": 0})
            continue
        least = min(size(o) for o in hits)
        minimal_hits = [o for o in hits if size(o) == least]
        chosen = sorted({label(o) for o in minimal_hits})
        rules.append({"rule": sorted(rule), "verdict": "UNIQUE" if len(chosen) == 1 else "AMBIGUOUS",
                      "selected": chosen[:12], "selected_count": len(chosen), "size": least,
                      "matching_objects": len(hits), "matching_labels": len({label(o) for o in hits}),
                      "minimal_objects": len(minimal_hits)})
    by_rule = {frozenset(r["rule"]): r for r in rules}
    selections = {}
    for row in rules:
        if row["verdict"] != "UNIQUE":
            continue
        selections.setdefault(row["selected"][0], []).append(frozenset(row["rule"]))
    summary = []
    for chosen, makers in selections.items():
        short = min(len(m) for m in makers)
        minimal = sorted((m for m in makers if not any(o < m for o in makers)),
                         key=lambda m: (len(m), tuple(sorted(m))))
        bearing = {}
        for m in minimal:
            for clause in sorted(m):
                other = by_rule[m-{clause}]
                bearing[f"{'+'.join(sorted(m))} without {clause}"] = (
                    other["verdict"], other["selected"][:4], other["size"])
        summary.append({"selected": chosen, "shortest_rule_length": short,
                        "minimal_rules": [sorted(m) for m in minimal], "ablation": bearing,
                        "selection_support": [{"rule": sorted(m),
                            "matching_objects": by_rule[m]["matching_objects"],
                            "matching_labels": by_rule[m]["matching_labels"],
                            "minimal_objects": by_rule[m]["minimal_objects"],
                            "minimal_labels": by_rule[m]["selected_count"],
                            "requires_size_minimization": by_rule[m]["matching_labels"] > 1}
                            for m in minimal]})
    counts = {v: sum(r["verdict"] == v for r in rules) for v in ("EMPTY", "UNIQUE", "AMBIGUOUS")}
    return {"objects": len(objects), "clauses": list(pool), "rules": len(rules), "verdicts": counts,
            "selection_semantics": "minimum size, then distinct labels; not all satisfying objects",
            "unique_selections": sorted(summary, key=lambda r: (r["shortest_rule_length"], str(r["selected"])))}


# ---------------------------------------------------------------- A: binary codes (distinctions that test each other)
def subspaces(n, k):
    """Every k-dimensional subspace of F_2^n as a tuple of row bitmasks in reduced echelon form."""
    for pivots in combinations(range(n), k):
        free = [[c for c in range(p+1, n) if c not in pivots] for p in pivots]
        for bits in product(*[range(1 << len(f)) for f in free]):
            rows = []
            for p, f, b in zip(pivots, free, bits):
                row = 1 << p
                for j, c in enumerate(f):
                    if b >> j & 1:
                        row |= 1 << c
                rows.append(row)
            yield tuple(rows)


def weights(rows, n):
    dist = [0]*(n+1)
    words = [0]
    for row in rows:
        words += [w ^ row for w in words]
    for w in words:
        dist[bin(w).count("1")] += 1
    return tuple(dist)


def codes(bound=8):
    out = []
    for n in range(1, bound+1):
        for k in range(1, n+1):
            for rows in subspaces(n, k):
                w = weights(rows, n)
                if w[1] or any(not any(r >> c & 1 for r in rows) for c in range(n)):
                    continue                      # a position no word uses, or a word of weight one, is no test
                orth = all(bin(a & b).count("1") % 2 == 0 for a in rows for b in rows)
                clauses = set()
                if orth:
                    clauses.add("self_orthogonal")
                if orth and 2*k == n:
                    clauses.add("self_dual")
                if all(c == 0 for i, c in enumerate(w) if i % 2):
                    clauses.add("even")
                if all(c == 0 for i, c in enumerate(w) if i % 4):
                    clauses.add("doubly_even")
                if min(i for i, c in enumerate(w) if i and c) >= 3:
                    clauses.add("corrects_one_error")
                out.append({"n": n, "k": k, "w": w, "clauses": clauses, "rows": rows})
    return out


def construction_a(code):
    """Lattice (1/sqrt 2)(C + 2 Z^n): unimodular iff C self dual, even iff doubly even; vectors of norm 2 counted."""
    n, w = code["n"], code["w"]
    return {"rank": n, "unimodular": "self_dual" in code["clauses"], "even": "doubly_even" in code["clauses"],
            "norm_two_vectors": 2*n + 16*w[4] if n >= 4 else 2*n}


def one_orbit(rows_list, n):
    """Do the given codes form one orbit under coordinate permutations? (breadth first over transpositions)"""
    def canon(rows):
        words = [0]
        for row in rows:
            words += [w ^ row for w in words]
        return frozenset(words)
    todo, seen, target = [canon(rows_list[0])], {canon(rows_list[0])}, {canon(r) for r in rows_list}
    while todo:
        c = todo.pop()
        for a, b in combinations(range(n), 2):
            d = frozenset(w if (w >> a & 1) == (w >> b & 1) else w ^ (1 << a | 1 << b) for w in c)
            if d not in seen:
                seen.add(d)
                todo.append(d)
    return seen == target


# ---------------------------------------------------------------- B: two-part cuts of a self-dual even lattice
def disc(t, n):
    """Discriminant group of a root lattice: list of elements, quadratic form mod 2, negation."""
    if t == "A":
        return [(k,) for k in range(n+1)], lambda x: Fraction(x[0]*(n+1-x[0]), n+1) % 2, lambda x: ((-x[0]) % (n+1),)
    if t == "D" and n % 2:
        q = {0: Fraction(0), 1: Fraction(n, 4) % 2, 2: Fraction(1), 3: Fraction(n, 4) % 2}
        return [(k,) for k in range(4)], lambda x: q[x[0]], lambda x: ((-x[0]) % 4,)
    if t == "D":
        q = {(0, 0): Fraction(0), (1, 0): Fraction(n, 4) % 2, (0, 1): Fraction(n, 4) % 2, (1, 1): Fraction(1)}
        return list(q), lambda x: q[x], lambda x: x
    if t == "E":
        m = {6: 3, 7: 2, 8: 1}[n]
        v = {6: Fraction(4, 3), 7: Fraction(3, 2), 8: Fraction(0)}[n]
        return [(k,) for k in range(m)], lambda x: (v*x[0]*x[0]) % 2, lambda x: ((-x[0]) % m,)


def add(t, n, x, y):
    if t == "A":
        return ((x[0]+y[0]) % (n+1),)
    if t == "D" and n % 2:
        return ((x[0]+y[0]) % 4,)
    if t == "D":
        return ((x[0]+y[0]) % 2, (x[1]+y[1]) % 2)
    return ((x[0]+y[0]) % {6: 3, 7: 2, 8: 1}[n],)


IRREDUCIBLE = [("A", n) for n in range(1, 9)] + [("D", n) for n in range(4, 9)] + [("E", n) for n in (6, 7, 8)]


def root_lattices(rank):
    """All multisets of irreducible root lattices of the given total rank."""
    def rec(left, start):
        if left == 0:
            yield ()
            return
        for i in range(start, len(IRREDUCIBLE)):
            t, n = IRREDUCIBLE[i]
            if n <= left:
                for rest in rec(left-n, i):
                    yield ((t, n),)+rest
    return list(rec(rank, 0))


def group(parts):
    data = [disc(t, n) for t, n in parts]
    elements = list(product(*[d[0] for d in data]))
    q = lambda x: sum((d[1](c) for d, c in zip(data, x)), Fraction(0)) % 2
    plus = lambda x, y: tuple(add(t, n, a, b) for (t, n), a, b in zip(parts, x, y))
    neg = lambda x: tuple(d[2](c) for d, c in zip(data, x))
    return elements, q, plus, neg


def anti_isometric(left, right):
    """Is there a group isomorphism disc(left) -> disc(right) with q_right(f(x)) = -q_left(x)? (generators searched)"""
    a, qa, pa, _ = group(left)
    b, qb, pb, _ = group(right)
    if len(a) != len(b):
        return False
    zero_a, zero_b = a[0], b[0]

    def span(gens, plus, zero):
        seen, todo = {zero}, [zero]
        while todo:
            x = todo.pop()
            for g in gens:
                y = plus(x, g)
                if y not in seen:
                    seen.add(y)
                    todo.append(y)
        return seen
    gens = []
    while len(span(gens, pa, zero_a)) < len(a):
        gens.append(next(x for x in a if x not in span(gens, pa, zero_a)))
    if not gens:
        return True
    order = lambda x, plus, zero: next(k for k in range(1, 400) if _times(x, k, plus, zero) == zero)
    choices = [[y for y in b if order(y, pb, zero_b) == order(g, pa, zero_a) and (qb(y)+qa(g)) % 2 == 0] for g in gens]
    for images in product(*choices):
        table, todo, good = {zero_a: zero_b}, [zero_a], True
        while todo and good:
            x = todo.pop()
            for g, y in zip(gens, images):
                s, t = pa(x, g), pb(table[x], y)
                if s in table:
                    good = good and table[s] == t
                else:
                    table[s] = t
                    todo.append(s)
        if good and len(set(table.values())) == len(a) and all((qb(table[x])+qa(x)) % 2 == 0 for x in a):
            return True
    return False


def _times(x, k, plus, zero):
    y = zero
    for _ in range(k):
        y = plus(y, x)
    return y


def name(parts):
    return " + ".join(f"{t}{n}" for t, n in sorted(parts, key=lambda p: (-p[1], p[0]))) or "0"


def cuts(rank=8):
    """Primitive complementary pairs (K, K') of root lattices inside an even unimodular lattice of the rank."""
    out = []
    for r in range(1, rank//2+1):
        for left in root_lattices(r):
            for right in root_lattices(rank-r):
                if r == rank-r and left > right:
                    continue
                if not anti_isometric(left, right):
                    continue
                elements, q, plus, neg = group(left)
                zero = elements[0]
                cyclic = any(len({_times(x, k, plus, zero) for k in range(len(elements))}) == len(elements) for x in elements)
                clauses = set()
                if len(left) == 1 and len(right) == 1:
                    clauses.add("both_irreducible")
                if cyclic:
                    clauses.add("cyclic_glue")
                if any(neg(x) != x for x in elements):
                    clauses.add("glue_has_conjugate_pair")
                if any(neg(x) == x and x != zero for x in elements):
                    clauses.add("glue_has_neutral_class")
                out.append({"pair": f"{name(left)} | {name(right)}", "glue": len(elements), "clauses": clauses})
    return out


# ---------------------------------------------------------------- C: one vector and one clock generate a lattice
def bareiss(m):
    """Exact determinant and leading minors of an integer matrix."""
    m = [row[:] for row in m]
    n, prev, minors = len(m), 1, []
    for k in range(n):
        if m[k][k] == 0:
            return 0, minors+[0]
        minors.append(m[k][k])
        for i in range(k+1, n):
            for j in range(k+1, n):
                m[i][j] = (m[i][j]*m[k][k]-m[i][k]*m[k][j])//prev
        prev = m[k][k]
    return m[-1][-1], minors


def clock_sets(n):
    """Sets of distinct m with sum phi(m) = n: the characteristic polynomials of finite order with a cyclic vector."""
    cand = [m for m in range(1, 61) if sympy.totient(m) <= n]
    def rec(left, start):
        if left == 0:
            yield ()
            return
        for i in range(start, len(cand)):
            p = int(sympy.totient(cand[i]))
            if p <= left:
                for rest in rec(left-p, i+1):
                    yield (cand[i],)+rest
    return list(rec(n, 0))


def clock_lattices(rank, norms=(1, 2, 3, 4)):
    """Positive integer Gram matrices <C^i v, C^j v> for an isometry C of finite order with cyclic vector v."""
    x = sympy.symbols("x")
    s = sympy.symbols(f"s0:{rank}")
    out = []
    for ms in clock_sets(rank):
        f = sympy.Poly(sympy.prod([sympy.cyclotomic_poly(m, x) for m in ms]), x)
        a = [int(c) for c in reversed(f.all_coeffs())][:rank]
        comp = sympy.zeros(rank, rank)
        for i in range(rank-1):
            comp[i+1, i] = 1
        for i in range(rank):
            comp[i, rank-1] = -a[i]
        gram = sympy.Matrix(rank, rank, lambda i, j: s[abs(i-j)])
        equations = [e for e in (comp.T*gram*comp-gram) if e != 0]
        system = sympy.Matrix([[e.coeff(v) for v in s] for e in equations]) if equations else sympy.zeros(1, rank)
        basis = system.nullspace()
        if not basis:
            continue
        b = sympy.Matrix.hstack(*basis)
        _, pivots = b.T.rref()
        coords = b*(b.extract(list(pivots), list(range(b.shape[1]))).inv())      # free coordinates = pivot entries
        coords = [[Fraction(int(c.p), int(c.q)) for c in coords.row(i)] for i in range(rank)]
        if 0 not in pivots:
            continue
        for s0 in norms:
            ranges = [[s0] if p == 0 else range(-s0+1, s0) for p in pivots]
            for free in product(*ranges):
                seq = [sum(c*v for c, v in zip(row, free)) for row in coords]
                if any(v.denominator != 1 or abs(v) > s0 for v in seq):
                    continue
                seq = [int(v) for v in seq]
                g = [[seq[abs(i-j)] for j in range(rank)] for i in range(rank)]
                det, minors = bareiss(g)
                if det <= 0 or any(m <= 0 for m in minors):
                    continue
                clauses = {"generated_by_one_vector_and_one_clock"}
                if det == 1:
                    clauses.add("self_dual")
                if s0 % 2 == 0:
                    clauses.add("even")
                if len(ms) == 1:
                    clauses.add("clock_without_subclock")
                if 1 not in ms:
                    clauses.add("clock_fixes_nothing")
                out.append({"rank": rank, "det": det, "norm": s0, "clock": ms, "order": lcm(*ms), "seq": seq,
                            "clauses": clauses})
    return out


def norm_two_vectors(g):
    """Count integer vectors x with x^T G x = 2 for a positive Gram matrix (Fincke and Pohst, exact final test)."""
    n = len(g)
    m = sympy.Matrix(g)
    l, d = m.LDLdecomposition()
    l = [[float(l[i, j]) for j in range(n)] for i in range(n)]
    d = [float(d[i, i]) for i in range(n)]
    count = 0
    x = [0]*n

    def rec(i, left):
        nonlocal count
        if i < 0:
            if any(x) and sum(x[a]*g[a][b]*x[b] for a in range(n) for b in range(n)) == 2:
                count += 1
            return
        centre = -sum(l[j][i]*x[j] for j in range(i+1, n))
        radius = (max(left, 0)/d[i])**.5
        lo, hi = int(-(-(centre-radius-1e-9)//1)), int((centre+radius+1e-9)//1)
        for v in range(lo, hi+1):
            x[i] = v
            rec(i-1, left-d[i]*(v-centre)**2)
        x[i] = 0
    rec(n-1, 2+1e-9)
    return count


# ---------------------------------------------------------------- D: a surface folded onto a sphere with cone points
def folds(max_points=6, max_order=12):
    out = []
    for k in range(3, max_points+1):
        for orders in combinations_with_replacement(range(2, max_order+1), k):
            chi = 2-sum(Fraction(m-1, m) for m in orders)
            clauses = set()
            if chi == 0:
                clauses.add("flat")
            if chi > 0:
                clauses.add("closed_positive_curvature")
            if k >= 4:
                clauses.add("has_a_shape_parameter")
            if len(set(orders)) == 1:
                clauses.add("all_points_alike")
            if set(orders) == {2}:
                clauses.add("fold_is_a_mirror")
            if all(gcd(a, b) == 1 for a, b in combinations(orders, 2)):
                clauses.add("turn_orders_pairwise_coprime")
            out.append({"orders": orders, "clauses": clauses})
    return out


def double_cover(n):
    """n marks on the sphere (n even), double cover branched there, line bundle pulled back from the mark divisor."""
    genus = (n-2)//2
    even, odd = n+1, n//2+1
    return {"marks": n, "genus": genus, "modes": even+odd, "even": even, "odd": odd, "loops_of_punctured_sphere": n-1,
            "cover_is_flat": genus == 1, "modes_are_twice_marks": even+odd == 2*n, "odd_modes_are_the_loops": odd == n-1,
            "difference_is_half_the_marks": even-odd == n//2}


# ---------------------------------------------------------------- E: a relabelling symmetry with a second description
def class_size(n, k):
    return factorial(n)//(2**k*factorial(k)*factorial(n-2*k))


def symmetric(bound=40):
    out = []
    for n in range(2, bound+1):
        twins = [k for k in range(2, n//2+1) if class_size(n, k) == comb(n, 2)]
        out.append({"n": n, "clauses": {"swap_has_a_double"} if twins else set(), "twins": twins})
    return out


def six():
    """The second description of the permutations of six things: duads, synthemes, totals; cycle types exchanged."""
    duads = list(combinations(range(6), 2))
    synthemes = [s for s in combinations(duads, 3) if len({p for d in s for p in d}) == 6]
    totals = [t for t in combinations(synthemes, 5) if len({d for s in t for d in s}) == 15]
    index = {t: i for i, t in enumerate(totals)}

    def image(perm):
        move = lambda t: tuple(sorted(tuple(sorted(tuple(sorted((perm[a], perm[b]))) for a, b in s)) for s in t))
        return tuple(index[move(t)] for t in totals)

    def cycle_type(perm):
        seen, out = set(), []
        for i in range(len(perm)):
            if i not in seen:
                j, c = i, 0
                while j not in seen:
                    seen.add(j)
                    j, c = perm[j], c+1
                out.append(c)
        return tuple(sorted(out, reverse=True))
    table = {}
    for perm in permutations(range(6)):
        table.setdefault(cycle_type(perm), set()).add(cycle_type(image(perm)))
    images = {image(p) for p in permutations(range(6))}
    compose = lambda p, q: tuple(p[i] for i in q)
    a, b = (1, 0, 2, 3, 4, 5), (1, 2, 3, 4, 5, 0)
    homomorphism = image(compose(a, b)) == compose(image(a), image(b))
    line_of = {d: [s for s in synthemes if d in s] for d in duads}
    quadrangle = all(sum(1 for e in s if e != d and any(e in t and d in t for t in synthemes)) == 1
                     for d in duads for s in synthemes if d not in s)
    return {"duads": len(duads), "synthemes": len(synthemes), "totals": len(totals),
            "faithful": len(images) == 720, "homomorphism_checked_on_generators": homomorphism,
            "cycle_types_exchanged": {str(k): sorted(map(str, v)) for k, v in sorted(table.items()) if {k} != v},
            "cycle_types_kept": [str(k) for k, v in sorted(table.items()) if {k} == v],
            "lines_through_a_point": sorted({len(v) for v in line_of.values()}),
            "generalised_quadrangle_15_points_15_lines": quadrangle}


# ---------------------------------------------------------------- F: signs that no probe of low order can tell from their mirror
def sign_codes(bound=7):
    """Odd functions c on {+1,-1}^n (a set of sign words with its mirror set counted negatively).

    Blind to order m means: sum over words of c(s) s_i1 ... s_im = 0 for all indices. The odd functions are spanned
    by the characters of odd degree, so the solutions blind to orders 1..m are spanned by the characters of odd
    degree above m. Each such character is one object; its support is the set of words with product +1 on its
    positions. The kernel dimension is checked against exact linear algebra for n <= 5.
    """
    out = []
    for n in range(1, bound+1):
        for degree in range(1, n+1, 2):
            for positions in combinations(range(n), degree):
                clauses = {"not_its_own_mirror"}
                if degree > 1:
                    clauses.add("blind_to_first_order")
                if degree > 3:
                    clauses.add("blind_to_third_order")
                if degree > 5:
                    clauses.add("blind_to_fifth_order")
                if degree == n:
                    clauses.add("uses_every_sign")
                out.append({"n": n, "degree": degree, "positions": positions, "states": 2**(n-1), "clauses": clauses})
    return out


def sign_kernel(n, orders=(1, 3)):
    """Exact dimension of the space of odd functions on the n-cube blind to the given odd orders."""
    words = [w for w in product((1, -1), repeat=n) if w[0] == 1]          # one word of each mirror pair
    rows = []
    for m in orders:
        for idx in combinations(range(n), m):
            rows.append([sympy.prod([w[i] for i in idx]) for w in words])
    return len(words)-sympy.Matrix(rows).rank() if rows else len(words)


def spinor_check():
    """The words with an even number of minus signs among five: weights of a half spinor of D5, charges under 3 + 2."""
    words = [w for w in product((1, -1), repeat=5) if sympy.prod(w) == 1]
    y = lambda w: Fraction(-1, 3)*sum(Fraction(w[i], 2) for i in range(3))+Fraction(1, 2)*sum(Fraction(w[i], 2) for i in (3, 4))
    table = {}
    for w in words:
        table[str(y(w))] = table.get(str(y(w)), 0)+1
    cubic = sum(y(w)**3 for w in words)
    return {"states": len(words), "first_moments": [sum(w[i] for w in words) for i in range(5)],
            "third_moments_all_zero": all(sum(w[i]*w[j]*w[k] for w in words) == 0 for i, j, k in combinations(range(5), 3)),
            "charge_3_plus_2_multiplicities": dict(sorted(table.items())), "sum_of_charges": str(sum(y(w) for w in words)),
            "sum_of_cubed_charges": str(cubic), "norm_of_a_word_over_two": str(Fraction(5, 4))}


# ---------------------------------------------------------------- G: the two rules composed: blind signs, then self-duality
def completion(max_rank=3, max_norm=8):
    """Start from the five signs. Which smallest positive even lattice completes them to a self-dual even lattice?

    The sixteen words/2 span, by differences, a lattice D; the words lie in one class of D*/D. A partner K must have a
    discriminant form anti-isometric to that of D. All positive even Gram matrices of rank <= max_rank with diagonal
    <= max_norm and reduced off-diagonal entries are enumerated; classes are told apart by rank, determinant and the
    number of vectors of norm two.
    """
    words = [w for w in product((1, -1), repeat=5) if sympy.prod(w) == 1]
    diffs = sorted({tuple((a-b)//2 for a, b in zip(u, v)) for u in words for v in words})
    basis = sympy.Matrix([d for d in diffs if any(d)]).T.columnspace()
    roots = [d for d in diffs if sum(x*x for x in d) == 2]
    hnf = sympy.Matrix(roots)
    from sympy.matrices.normalforms import smith_normal_form
    index = abs(sympy.prod([d for d in smith_normal_form(hnf, domain=sympy.ZZ).diagonal() if d != 0]))   # det of lattice basis
    det_d = int(index)**2
    norm = Fraction(sum(x*x for x in words[0]), 4)
    classes = sorted({Fraction(sum(Fraction(a, 2)*Fraction(b, 2) for a, b in zip(words[0], w))) % 1 for w in words}, key=float)
    gauss = sum(sympy.exp(sympy.I*sympy.pi*q) for q in (0, Fraction(5, 4), 1, Fraction(5, 4)))/2
    phase = sympy.nsimplify(sympy.arg(sympy.simplify(gauss))/(2*sympy.pi))
    target = (2-norm) % 2
    found = {}
    for rank in range(1, max_rank+1):
        diag = [range(2, max_norm+1, 2)]*rank
        for d in product(*diag):
            if list(d) != sorted(d):
                continue
            pairs = list(combinations(range(rank), 2))
            for off in product(*[range(-(min(d[i], d[j])//2), min(d[i], d[j])//2+1) for i, j in pairs]):
                g = [[d[i] if i == j else 0 for j in range(rank)] for i in range(rank)]
                for (i, j), v in zip(pairs, off):
                    g[i][j] = g[j][i] = v
                det, minors = bareiss(g)
                if det != det_d or any(m <= 0 for m in minors):
                    continue
                inv = sympy.Matrix(g).inv()
                if all((2*x).is_integer for x in inv):
                    continue                                   # exponent two: not a cyclic group of order four
                if not any((Fraction(int(inv[i, i].p), int(inv[i, i].q))-target) % 2 == 0 for i in range(rank)
                           if not (2*inv[i, i]).is_integer or not all((2*x).is_integer for x in inv.col(i))):
                    continue
                key = (rank, det, norm_two_vectors(g))
                found.setdefault(key, g)
    short = lambda g: min((Fraction(int(v.p), int(v.q)) for v in (sympy.Matrix(g).inv()).diagonal()), default=None)
    return {"lattice_spanned_by_the_words": {"rank": 5, "determinant": det_d, "vectors_of_norm_two": len(roots),
                                             "norm_of_a_word": str(norm), "words_in_one_class": classes == [classes[0]] or len(classes) <= 2},
            "gauss_sum_phase_in_turns": str(phase), "partner_rank_mod_8": int((-phase*8) % 8),
            "required_norm_of_the_partner_class_mod_2": str(target),
            "partners_by_rank": {r: [{"determinant": k[1], "vectors_of_norm_two": k[2], "gram": g,
                                      "smallest_dual_diagonal": str(short(g))}
                                     for k, g in sorted(found.items()) if k[0] == r] for r in range(1, max_rank+1)},
            "bound": f"all positive even Gram matrices of rank <= {max_rank}, diagonal <= {max_norm}, reduced off-diagonal entries"}


def words_and_flips(bound=24):
    """n signs: a word is the vector of its signs over two, a flip changes two signs. When are they the same kind?

    Flips have norm 2 for every n; words have norm n/4. Words and flips together form a lattice with integer products
    and even norms exactly when n is divisible by 8 (products of two even words are (n - 2k)/4 with k even).
    Inside the n = 8 solution, the split into the smallest blind chiral part and the rest is read off.
    """
    rows = []
    for n in range(1, bound+1):
        products = {Fraction(n-2*k, 4) for k in range(0, n+1, 2)}
        rows.append({"n": n, "word_norm": str(Fraction(n, 4)), "same_norm_as_a_flip": Fraction(n, 4) == 2,
                     "integral_and_even": all(p.denominator == 1 for p in products) and Fraction(n, 4) % 2 == 0})
    rest = [w for w in product((1, -1), repeat=3) if sympy.prod(w) == 1]
    table = {(a, b): tuple(x*y for x, y in zip(a, b)) for a in rest for b in rest}
    return {"n_with_words_of_flip_norm": [r["n"] for r in rows if r["same_norm_as_a_flip"]],
            "n_integral_and_even": [r["n"] for r in rows if r["integral_and_even"]],
            "eight_signs": {"even_words": 2**7, "flips": 4*comb(8, 2), "together": 2**7+4*comb(8, 2),
                            "blind_chiral_part": 5, "rest": 3},
            "rest_of_three_signs": {"even_words": len(rest), "closed_under_product": set(table.values()) == set(rest),
                                    "every_word_squares_to_one": all(table[(a, a)] == (1, 1, 1) for a in rest),
                                    "one_identity_and_three_others": [list(w) for w in rest]}}


def completed_chirality():
    """Red-team test: does any mirror-odd label survive on the completed object?

    The 128 glue vectors of the completion are the even words of eight signs. A label is an odd function on them
    (64 mirror pairs). Blindness is tested with probes of order one and three (a) in all eight directions, (b) only
    in the five directions of the blind part, (c) in the five and in one, two or three of the partner directions.
    The dimension of the space of surviving labels is the exact corank of the probe matrix.
    """
    words = [w for w in product((1, -1), repeat=8) if sympy.prod(w) == 1 and w[0] == 1]
    def survivors(directions):
        rows = [[w[i] for w in words] for i in directions]
        rows += [[w[i]*w[j]*w[k] for w in words] for i, j, k in combinations(directions, 3)]
        return len(words)-sympy.Matrix(rows).rank()
    carrier = tuple(range(5))
    every = [w for w in product((1, -1), repeat=8) if sympy.prod(w) == 1]
    four = {w: w[5]*w[6]*w[7] for w in every}                                  # 16 against its mirror, four family words
    three = {w: 0 if w[5] == w[6] == w[7] else four[w] for w in every}         # fourth family word left out
    def blind(label, directions):
        return (all(sum(label[w]*w[i] for w in every) == 0 for i in directions) and
                all(sum(label[w]*w[i]*w[j]*w[k] for w in every) == 0 for i, j, k in combinations(directions, 3)))
    families = {name: {"words_marked_plus": sum(v == 1 for v in label.values()),
                       "blind_with_k_partner_directions_probed": {k: blind(label, carrier+tuple(range(5, 5+k))) for k in range(4)}}
                for name, label in (("four_families", four), ("three_families", three))}
    return {"mirror_pairs": len(words), "family_labels": families,
            "labels_blind_in_all_eight_directions": survivors(tuple(range(8))),
            "labels_blind_in_the_five_directions_only": survivors(carrier),
            "labels_blind_in_five_plus_k_partner_directions": {k: survivors(carrier+tuple(range(5, 5+k))) for k in (1, 2, 3)}}


def current_probes():
    """Which root currents of the completed object see a family label?

    The 240 norm-two vectors split into 40 of the blind part, 12 of the partner, 60 mixed and the 128 glue vectors.
    A current of root r probes a label f on the glue vectors through its charge <r, w>: order one, order three, and
    mixed order three for several currents. It also moves w to w + r. The largest set of root currents for which a
    label stays blind is found by closing the roots with vanishing first order under the full cubic form.
    """
    from fractions import Fraction
    every = [w for w in product((1, -1), repeat=8) if sympy.prod(w) == 1]
    four = {w: w[5]*w[6]*w[7] for w in every}
    three = {w: 0 if w[5] == w[6] == w[7] else four[w] for w in every}
    roots = []
    for i, j in combinations(range(8), 2):
        for a, b in product((1, -1), repeat=2):
            r = [0]*8; r[i], r[j] = a, b; roots.append(tuple(r))
    def kind(r):
        s = [k for k in range(8) if r[k]]
        return "blind_part" if s[1] < 5 else "partner" if s[0] >= 5 else "mixed"
    def charge(r, w):
        return Fraction(sum(a*b for a, b in zip(r, w)), 2)
    def cubic(f, a, b, c):
        return sum(f[w]*charge(a, w)*charge(b, w)*charge(c, w) for w in every)
    def blind_span(f, rs):
        basis = []
        for r in rs:
            if sympy.Matrix(basis+[r]).rank() > len(basis):
                basis.append(r)
        return (all(sum(f[w]*charge(r, w) for w in every) == 0 for r in basis) and
                all(cubic(f, a, b, c) == 0 for i, a in enumerate(basis) for j, b in enumerate(basis[i:], i) for c in basis[j:]))
    out = {}
    for name, f in (("four_families", four), ("three_families", three)):
        row = {}
        for k in ("blind_part", "partner", "mixed"):
            rs = [r for r in roots if kind(r) == k]
            first = [sum(f[w]*charge(r, w) for w in every) for r in rs]
            third = [cubic(f, r, r, r) for r in rs]
            moves = set()
            for r in rs:
                keep = flip = leave = 0
                for w in every:
                    if charge(r, w) == -1 and f[w]:
                        v = tuple(-x if y else x for x, y in zip(w, r))
                        keep += f[v] == f[w]; flip += f[v] == -f[w]; leave += f[v] == 0
                moves.add((keep, flip, leave))
            row[k] = {"currents": len(rs), "seen_at_first_order": sum(x != 0 for x in first),
                      "seen_at_third_order": sum(x != 0 for x in third),
                      "first_order_values": sorted({int(x) for x in first}),
                      "marked_words_kept_flipped_removed": sorted(moves)}
        candidates = {"blind_part_only": [r for r in roots if kind(r) == "blind_part"],
                      "blind_part_and_partner": [r for r in roots if kind(r) != "mixed"],
                      "all_roots_blind_at_first_order": [r for r in roots if sum(f[w]*charge(r, w) for w in every) == 0],
                      "blind_part_and_mixed_into_two_partner_directions": [r for r in roots if not r[7]],
                      "all_240_minus_glue": roots}
        row["sets_of_currents"] = {k: {"currents": len(v), "rank": sympy.Matrix(v).rank(), "label_blind": blind_span(f, v)}
                                   for k, v in candidates.items()}
        mixed = [r for r in roots if kind(r) == "mixed" and r[0] == 1 and r[5]+r[6]+r[7] == 1]
        row["mixed_cubic_of_three_currents_into_three_partner_directions"] = int(8*cubic(f, *mixed))
        out[name] = row
    return out


# ---------------------------------------------------------------- record structures met by a selection (named afterwards)
RECORD = {
    "codes": {"n=8 k=4 weights (1, 0, 0, 0, 14, 0, 0, 0, 1)":
              "extended Hamming code; Construction A gives the E8 lattice (Atlas: aspect:website-finite-e8-code-quantum)"},
    "cuts": {"D5 | A3": "the D5 + A3 chart with Z4 glue (Atlas: structure:e8-glue, claim:BOOT.CLASS.01)",
             "E6 | A2": "E6 + A2 cut, recorded as one of the other cuts of E8", "A4 | A4": "A4 + A4 cut (two five-clocks)",
             "E7 | A1": "E7 + A1 cut", "D4 | D4": "D4 + D4 cut"},
    "clocks": {"rank 8 det 1 even": "E8 with a clock that has a cyclic vector (Atlas: Coxeter Krylov basis, "
               "finding:coxeter-seventh-power-has-complete-thirty-e8-lattice-lifts; SEAM.GALOIS.01)"},
    "folds": {"(2, 2, 2, 2)": "sphere with four marks, cover a torus (Atlas: map:shared-fourmark-origin, GESAMT2-20261120)",
              "(2, 3, 5)": "icosahedral fold; orders 2, 3, 5 with product 30 = Coxeter number of E8 (clock network)"},
    "signs": {"n=5 degree 5": "five signs, sixteen words: the half spinor of the carrier D5 (Atlas: map:carrier-d5)"},
    "symmetric": {"n=6": "six things; pairs of five (Petersen), type 3 + 2 of order six (SECHS20261119), doily"},
}


def run(root=None, code_bound=8):
    out = {"status": "GENESIS_RULES_SEARCHED", "physical_solution": False, "scope": SCOPE, "universes": {}}

    cs = codes(code_bound)
    res = search(cs, ["self_orthogonal", "self_dual", "even", "doubly_even", "corrects_one_error"],
                 lambda o: o["n"], lambda o: f"n={o['n']} k={o['k']} weights {o['w']}")
    hamming = [c for c in cs if c["n"] == 8 and {"self_dual", "doubly_even"} <= c["clauses"]]
    res["derived"] = {"self_dual_doubly_even_codes_of_length_8": len(hamming),
                      "one_orbit_under_relabelling": one_orbit([c["rows"] for c in hamming], 8) if hamming else None,
                      "construction_a": construction_a(hamming[0]) if hamming else None,
                      "self_dual_codes_by_length": {n: sum(1 for c in cs if c["n"] == n and "self_dual" in c["clauses"])
                                                    for n in range(2, code_bound+1, 2)}}
    res["bound"] = f"all binary linear codes of length at most {code_bound} without unused positions or weight-one words"
    out["universes"]["codes"] = res

    ct = cuts(8)
    res = search(ct, ["both_irreducible", "cyclic_glue", "glue_has_conjugate_pair", "glue_has_neutral_class"],
                 lambda o: o["glue"], lambda o: o["pair"])
    res["derived"] = {"all_primitive_pairs": sorted(f"{c['pair']} (glue {c['glue']})" for c in ct),
                      "ranks_below_8_with_any_pair": [r for r in range(2, 8) if cuts(r)]}
    res["bound"] = "all pairs of root lattices (reducible allowed) of total rank 8 with anti-isometric discriminant forms"
    out["universes"]["cuts"] = res

    cl = [row for rank in range(1, 9) for row in clock_lattices(rank)]
    res = search(cl, ["self_dual", "even", "clock_without_subclock", "clock_fixes_nothing"],
                 lambda o: o["rank"], lambda o: f"rank {o['rank']} det {o['det']} {'even' if o['norm'] % 2 == 0 else 'odd'}")
    e8 = [c for c in cl if c["rank"] == 8 and c["det"] == 1 and c["norm"] % 2 == 0]
    by_clock = {}
    for c in e8:
        by_clock.setdefault(c["clock"], []).append(c)
    sample = min(e8, key=lambda c: (c["norm"], c["order"])) if e8 else None
    res["derived"] = {
        "even_self_dual_clock_lattices_of_rank_8": len(e8),
        "their_clocks": [{"cyclotomic_factors": list(k), "order": lcm(*k), "generator_norms": sorted({c["norm"] for c in v}),
                          "sequences": len(v)} for k, v in sorted(by_clock.items(), key=lambda kv: -lcm(*kv[0]))],
        "clock_orders": sorted({c["order"] for c in e8}),
        "clocks_without_subclock": sorted(k[0] for k in by_clock if len(k) == 1),
        "sample": sample and {"clock": list(sample["clock"]), "sequence": sample["seq"],
                              "norm_two_vectors": norm_two_vectors([[sample["seq"][abs(i-j)] for j in range(8)] for i in range(8)])},
        "even_self_dual_below_rank_8": sum(1 for c in cl if c["rank"] < 8 and c["det"] == 1 and c["norm"] % 2 == 0)}
    res["bound"] = "ranks 1 to 8, every finite-order isometry with a cyclic vector of norm 1 to 4"
    out["universes"]["clocks"] = res

    fs = folds()
    res = search(fs, ["flat", "closed_positive_curvature", "has_a_shape_parameter", "all_points_alike", "fold_is_a_mirror",
                      "turn_orders_pairwise_coprime"],
                 lambda o: (len(o["orders"]), sum(o["orders"])), lambda o: str(o["orders"]))
    res["derived"] = {"flat_folds": [str(f["orders"]) for f in fs if "flat" in f["clauses"]],
                      "closed_folds_without_the_family_2_2_n": [str(f["orders"]) for f in fs
                                                                if "closed_positive_curvature" in f["clauses"]
                                                                and f["orders"][:2] != (2, 2)],
                      "double_covers": [double_cover(n) for n in range(2, 13, 2)]}
    res["bound"] = "spheres with 3 to 6 cone points of order 2 to 12"
    out["universes"]["folds"] = res

    sg = sign_codes()
    res = search(sg, ["not_its_own_mirror", "blind_to_first_order", "blind_to_third_order", "blind_to_fifth_order",
                      "uses_every_sign"], lambda o: o["n"], lambda o: f"n={o['n']} degree {o['degree']}")
    res["derived"] = {"kernel_dimension_blind_to_orders_1_and_3": {n: sign_kernel(n) for n in range(1, 6)},
                      "kernel_dimension_blind_to_order_1": {n: sign_kernel(n, (1,)) for n in range(1, 6)},
                      "five_signs": spinor_check(),
                      "states_by_blindness": {"order 1": 4, "orders 1, 3": 16, "orders 1, 3, 5": 64}}
    res["bound"] = "odd functions on the cube of 1 to 7 signs; one object per character of odd degree (they span all solutions)"
    out["universes"]["signs"] = res

    sy = symmetric()
    res = search(sy, ["swap_has_a_double"], lambda o: o["n"], lambda o: f"n={o['n']}")
    res["derived"] = {"all_n_with_a_double_up_to_40": [o["n"] for o in sy if o["clauses"]], "six": six()}
    res["bound"] = "permutations of 2 to 40 things; a second description needs a class of involutions as large as the swaps"
    out["universes"]["symmetric"] = res

    out["composed_rule"] = {"rule": ["signs: not_its_own_mirror + blind_to_third_order", "then: self_dual + even"],
                            "result": completion()}
    out["red_team_test"] = {"question": "is the completed object still chiral and blind?", "result": completed_chirality()}
    out["current_test"] = {"question": "which root currents of the completed object see a family label?",
                           "result": current_probes()}
    out["one_rule_variant"] = {"rule": "words and flips of n signs are the same kind of vector", "result": words_and_flips()}
    met = []
    for key, res in out["universes"].items():
        for sel in res["unique_selections"]:
            for pattern, meaning in RECORD[key].items():
                if pattern in sel["selected"]:
                    met.append({"universe": key, "selected": sel["selected"], "shortest_rule": sel["minimal_rules"][0],
                                "record": meaning})
    out["selections_that_meet_the_record"] = met
    out["rules_generated"] = sum(res["rules"] for res in out["universes"].values())
    out["objects_enumerated"] = sum(res["objects"] for res in out["universes"].values())
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--code-bound", type=int, default=8)
    parser.add_argument("--foreign", action="store_true",
                        help="universes without TFPT vocabulary, scored against measured targets (foreign_universes.py)")
    args = parser.parse_args()
    if args.foreign:
        from foreign_universes import run as foreign_run
        print(json.dumps(foreign_run(), indent=1, default=str, ensure_ascii=False))
    else:
        print(json.dumps(run(code_bound=args.code_bound), indent=1, default=str, ensure_ascii=False))
