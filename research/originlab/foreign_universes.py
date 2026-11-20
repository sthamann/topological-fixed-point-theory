"""Foreign universes for the genesis rule search, scored against measured targets with a random baseline.

The universes here use no TFPT vocabulary: small connected graphs, elementary cellular automata and regular
polytopes. Each object carries clauses (for the rule search in `genesis_rules.search`) and integer invariants.
Every object's invariants are combined by a fixed small expression grammar and compared with the measured
dimensionless targets of docs/research/EXPLORATION.md. The same grammar is applied to pseudo-objects whose
invariants are drawn at random from the universe's own pooled invariants. A target counts only where real objects
beat that baseline. A hit is a lead for a human, not a derivation. No TFPT number enters; no ledger, no Atlas.
"""
from itertools import combinations, permutations, product
from math import comb, isqrt, pi, sqrt
import random

import numpy

from genesis_rules import search

SCOPE = ("exhaustive enumeration inside stated bounds; number hits come from a fixed expression grammar and are "
         "compared with a shuffled-invariant baseline; a hit is a lead, not physics")

# Measured dimensionless targets (see docs/research/EXPLORATION.md for sources).
TARGETS = {
    "alpha_inverse": 137.035999177,
    "proton_electron_mass_ratio": 1836.15267,
    "muon_electron_mass_ratio": 206.768283,
    "tau_muon_mass_ratio": 16.8170,
    "koide_q": 0.666661,
    "sin2_theta_w": 0.23129,
    "cabibbo_vus": 0.2243,
    "pmns_sin2_theta12": 0.307,
    "pmns_sin2_theta13": 0.0220,
    "alpha_s_mz": 0.1180,
    "omega_lambda": 0.685,
}
TOLERANCE = 1e-3


# ---------------------------------------------------------------- expression grammar
def expressions(values):
    """Yield (value, text) for a fixed grammar over named positive integer invariants."""
    items = [(k, float(v)) for k, v in values.items() if v and v > 0]
    for k, v in items:
        yield v, k
        yield v*pi, f"{k}*pi"
        yield v/pi, f"{k}/pi"
        yield sqrt(v), f"sqrt({k})"
    for (a, x), (b, y) in permutations(items, 2):
        yield x/y, f"{a}/{b}"
        yield x*y, f"{a}*{b}"
        yield x+y, f"{a}+{b}"
        if x > y:
            yield x-y, f"{a}-{b}"
        yield x/y*pi, f"{a}/{b}*pi"
    for (a, x), (b, y), (c, z) in permutations(items, 3):
        yield x*y/z, f"{a}*{b}/{c}"
        yield (x+y)/z, f"({a}+{b})/{c}"
        if a < b:
            yield x/(y*z), f"{a}/({b}*{c})"


def hits(values, tolerance=TOLERANCE):
    """For each target the best expression within tolerance, else None."""
    best = {}
    for value, text in expressions(values):
        for name, target in TARGETS.items():
            err = abs(value-target)/target
            if err <= tolerance and (name not in best or err < best[name][1]):
                best[name] = (text, err)
    return best


def score(objects, label, samples=10, seed=20261120):
    """Per-target count of objects with a hit, versus pseudo-objects drawn from the pooled invariants."""
    real = {name: [] for name in TARGETS}
    for o in objects:
        for name, (text, err) in hits(o["invariants"]).items():
            real[name].append({"object": label(o), "expression": text, "relative_error": err})
    pool = [v for o in objects for v in o["invariants"].values() if v and v > 0]
    rng = random.Random(seed)
    base = {name: 0 for name in TARGETS}
    for _ in range(samples):
        for o in objects:
            keys = list(o["invariants"])
            fake = {k: rng.choice(pool) for k in keys}
            for name in hits(fake):
                base[name] += 1
    rows = {}
    for name in TARGETS:
        n_real, n_base = len(real[name]), base[name]/samples
        rows[name] = {"objects_hit": n_real, "baseline_hits": n_base,
                      "beats_baseline": n_real >= 3 and n_real > 3*max(n_base, 0.5),
                      "examples": sorted(real[name], key=lambda r: r["relative_error"])[:5]}
    return rows


# ---------------------------------------------------------------- universe 1: small connected graphs
def _degree_canonical(n, edges):
    adj = [[0]*n for _ in range(n)]
    for a, b in edges:
        adj[a][b] = adj[b][a] = 1
    deg = [sum(r) for r in adj]
    order = sorted(range(n), key=lambda v: deg[v])
    blocks, start = [], 0
    for i in range(1, n+1):
        if i == n or deg[order[i]] != deg[order[start]]:
            blocks.append(order[start:i])
            start = i
    best, autos = None, 0
    for choice in product(*[list(permutations(b)) for b in blocks]):
        perm = [v for block in choice for v in block]
        key = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
        if best is None or key > best:
            best = key
    for choice in product(*[list(permutations(b)) for b in blocks]):
        perm = [v for block in choice for v in block]
        mapping = dict(zip(order, perm))
        if all(adj[mapping[a]][mapping[b]] == adj[a][b] for a in range(n) for b in range(a+1, n)):
            autos += 1
    return best, adj, autos


def _connected(adj):
    n, seen, stack = len(adj), {0}, [0]
    while stack:
        v = stack.pop()
        for w in range(n):
            if adj[v][w] and w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == n


def graphs(max_n=6):
    """All connected simple graphs on 2..max_n vertices up to isomorphism."""
    out = []
    for n in range(2, max_n+1):
        pairs = list(combinations(range(n), 2))
        seen = set()
        for mask in range(1, 1 << len(pairs)):
            edges = [p for i, p in enumerate(pairs) if mask >> i & 1]
            if len(edges) < n-1:
                continue
            key, adj, autos = _degree_canonical(n, edges)
            if key in seen or not _connected(adj):
                continue
            seen.add(key)
            out.append(_graph_object(n, edges, adj, autos, key))
    return out


def _graph_object(n, edges, adj, autos, key):
    a = numpy.array(adj, dtype=float)
    eig = numpy.linalg.eigvalsh(a)
    lap = numpy.diag(a.sum(1))-a
    trees = round(numpy.linalg.det(lap[1:, 1:]))
    deg = [sum(r) for r in adj]
    triangles = sum(adj[i][j]*adj[j][k]*adj[i][k] for i, j, k in combinations(range(n), 3))
    dist = numpy.full((n, n), n)
    for s in range(n):
        dist[s][s], frontier, d = 0, [s], 0
        while frontier:
            d += 1
            nxt = [w for v in frontier for w in range(n) if adj[v][w] and dist[s][w] > d]
            for w in nxt:
                dist[s][w] = d
            frontier = list(set(nxt))
    two_colour = triangles == 0 and _bipartite(adj)
    complement = [1-adj[i][j] for i in range(n) for j in range(i+1, n)]
    clauses = set()
    if len(set(deg)) == 1:
        clauses.add("regular")
    if two_colour:
        clauses.add("bipartite")
    if triangles == 0:
        clauses.add("triangle_free")
    if all(d % 2 == 0 for d in deg):
        clauses.add("eulerian")
    if len(edges) == n-1:
        clauses.add("tree")
    if all(abs(x-round(x)) < 1e-9 for x in eig):
        clauses.add("integral_spectrum")
    if 2*len(edges) == n*(n-1)//2:
        comp_key = _degree_canonical(n, [p for p, c in zip(combinations(range(n), 2), complement) if c])[0]
        if comp_key == key:
            clauses.add("self_complementary")
    if autos >= n and len(set(deg)) == 1 and _vertex_transitive(adj):
        clauses.add("vertex_transitive")
    return {"n": n, "edges": len(edges), "edge_list": edges, "clauses": clauses,
            "invariants": {"vertices": n, "edges": len(edges), "automorphisms": autos, "spanning_trees": trees,
                           "triangles": triangles, "diameter": int(dist.max()), "max_degree": max(deg)},
            "spectral_radius": float(eig[-1])}


def _bipartite(adj):
    n, colour = len(adj), {0: 0}
    stack = [0]
    while stack:
        v = stack.pop()
        for w in range(n):
            if adj[v][w]:
                if w not in colour:
                    colour[w] = 1-colour[v]
                    stack.append(w)
                elif colour[w] == colour[v]:
                    return False
    return True


def _vertex_transitive(adj):
    n = len(adj)
    reach = {0}
    for perm in permutations(range(n)):
        if perm[0] in reach:
            continue
        if all(adj[perm[a]][perm[b]] == adj[a][b] for a in range(n) for b in range(a+1, n)):
            reach.add(perm[0])
            if len(reach) == n:
                return True
    return len(reach) == n


# ---------------------------------------------------------------- universe 2: elementary cellular automata
def _step(rule, state, n):
    return sum(((rule >> (((state >> ((i+1) % n)) & 1) << 2 | ((state >> i) & 1) << 1 | ((state >> ((i-1) % n)) & 1))) & 1) << i
               for i in range(n))


def _orbit_data(rule, n):
    image = [_step(rule, s, n) for s in range(1 << n)]
    gardens = (1 << n)-len(set(image))
    on_cycle, cycles, longest = set(), 0, 0
    for s in range(1 << n):
        path, x = {}, s
        while x not in path:
            path[x] = len(path)
            x = image[x]
        if x not in on_cycle:
            y, length = x, 0
            while True:
                on_cycle.add(y)
                y, length = image[y], length+1
                if y == x:
                    break
            cycles += 1
            longest = max(longest, length)
    return gardens, len(on_cycle), cycles, longest


def automata(ring=8):
    """The 88 inequivalent elementary cellular automata (mirror and complement identified)."""
    def mirror(r):
        return sum(((r >> (((k & 1) << 2) | (k & 2) | (k >> 2))) & 1) << k for k in range(8))

    def complement(r):
        return sum((1-((r >> (7-k)) & 1)) << k for k in range(8))

    classes = sorted({min(r, mirror(r), complement(r), mirror(complement(r))) for r in range(256)})
    out = []
    for r in classes:
        gardens, cyc_states, cycles, longest = _orbit_data(r, ring)
        table = [(r >> k) & 1 for k in range(8)]
        clauses = set()
        if gardens == 0 and all(_orbit_data(r, m)[0] == 0 for m in (5, 6, 7)):
            clauses.add("reversible_on_rings")
        if all(table[k] == (sum(c*((k >> i) & 1) for i, c in enumerate((table[1], table[2], table[4]))) % 2)
               for k in range(8)) and table[0] == 0:
            clauses.add("additive")
        if mirror(r) == r:
            clauses.add("mirror_symmetric")
        if complement(r) == r:
            clauses.add("complement_symmetric")
        if all(bin(_step(r, s, ring)).count("1") == bin(s).count("1") for s in range(1 << ring)):
            clauses.add("conserves_number")
        if table[0] == 0 and table[7] == 1:
            clauses.add("keeps_uniform_states")
        out.append({"rule": r, "ones": sum(table), "clauses": clauses,
                    "invariants": {"garden_states": gardens, "cycle_states": cyc_states, "cycles": cycles,
                                   "longest_cycle": longest, "table_ones": sum(table)}})
    return out


# ---------------------------------------------------------------- universe 3: regular polytopes
def polytopes(max_dim=8):
    """Regular convex polytopes of dimension 2..max_dim with f-vectors and symmetry orders."""
    out = []

    def add(name, dim, f, group, clauses):
        out.append({"name": name, "dim": dim, "clauses": set(clauses),
                    "invariants": {"dimension": dim, "vertices": f[0], "edges": f[1], "facets": f[-1],
                                   "symmetry_order": group, "faces_total": sum(f)}})

    for p in range(3, 13):
        add(f"{p}-gon", 2, [p, p], 2*p, {"self_dual"} | ({"family_in_every_dimension"} if p in (3, 4) else set()))
    platonic = [("tetrahedron", [4, 6, 4], 24, {"self_dual", "simplicial", "simple", "family_in_every_dimension"}),
                ("cube", [8, 12, 6], 48, {"simple", "family_in_every_dimension"}),
                ("octahedron", [6, 12, 8], 48, {"simplicial", "family_in_every_dimension"}),
                ("dodecahedron", [20, 30, 12], 120, {"simple", "exceptional"}),
                ("icosahedron", [12, 30, 20], 120, {"simplicial", "exceptional"})]
    for name, f, g, cl in platonic:
        add(name, 3, f, g, cl)
    four = [("5-cell", [5, 10, 10, 5], 120, {"self_dual", "simplicial", "simple", "family_in_every_dimension"}),
            ("tesseract", [16, 32, 24, 8], 384, {"simple", "family_in_every_dimension"}),
            ("16-cell", [8, 24, 32, 16], 384, {"simplicial", "family_in_every_dimension"}),
            ("24-cell", [24, 96, 96, 24], 1152, {"self_dual", "exceptional"}),
            ("120-cell", [600, 1200, 720, 120], 14400, {"simple", "exceptional"}),
            ("600-cell", [120, 720, 1200, 600], 14400, {"simplicial", "exceptional"})]
    for name, f, g, cl in four:
        add(name, 4, f, g, cl)
    for n in range(5, max_dim+1):
        fact = 1
        for k in range(2, n+2):
            fact *= k
        add(f"{n}-simplex", n, [comb(n+1, k+1) for k in range(n)], fact,
            {"self_dual", "simplicial", "simple", "family_in_every_dimension"})
        hyper = (2**n)*(fact//(n+1))
        add(f"{n}-cube", n, [2**(n-k)*comb(n, k) for k in range(n)], hyper, {"simple", "family_in_every_dimension"})
        add(f"{n}-orthoplex", n, [2**(k+1)*comb(n, k+1) for k in range(n)], hyper,
            {"simplicial", "family_in_every_dimension"})
    for o in out:
        if o["dim"] >= 3 and "exceptional" not in o["clauses"] and "family_in_every_dimension" not in o["clauses"]:
            o["clauses"].add("exceptional")
    return out


# ---------------------------------------------------------------- run
def run(root=None, graph_bound=6, ring=8):
    out = {"status": "FOREIGN_UNIVERSES_SEARCHED", "physical_solution": False, "scope": SCOPE,
           "targets": TARGETS, "tolerance": TOLERANCE, "universes": {}}
    gs = graphs(graph_bound)
    res = search(gs, ["regular", "bipartite", "triangle_free", "eulerian", "tree", "integral_spectrum",
                      "self_complementary", "vertex_transitive"],
                 lambda o: (o["n"], o["edges"]), lambda o: f"n={o['n']} {o['edge_list']}")
    res["targets"] = score(gs, lambda o: f"n={o['n']} {o['edge_list']}")
    res["bound"] = f"all connected simple graphs on 2 to {graph_bound} vertices up to isomorphism"
    out["universes"]["graphs"] = res

    ca = automata(ring)
    res = search(ca, ["reversible_on_rings", "additive", "mirror_symmetric", "complement_symmetric",
                      "conserves_number", "keeps_uniform_states"],
                 lambda o: (o["ones"], o["rule"]), lambda o: f"rule {o['rule']}")
    res["targets"] = score(ca, lambda o: f"rule {o['rule']}")
    res["bound"] = f"88 elementary cellular automata, dynamics on rings of {ring} cells"
    out["universes"]["automata"] = res

    po = polytopes()
    res = search(po, ["self_dual", "simplicial", "simple", "family_in_every_dimension", "exceptional"],
                 lambda o: (o["dim"], o["invariants"]["vertices"]), lambda o: o["name"])
    res["targets"] = score(po, lambda o: o["name"])
    res["bound"] = "regular convex polytopes of dimension 2 to 8 (polygons up to 12 sides)"
    out["universes"]["polytopes"] = res

    leads = []
    for key, res in out["universes"].items():
        for name, row in res["targets"].items():
            if row["beats_baseline"]:
                leads.append({"universe": key, "target": name, "objects_hit": row["objects_hit"],
                              "baseline_hits": row["baseline_hits"], "best": row["examples"][:1]})
    out["leads_beating_baseline"] = leads
    out["objects_enumerated"] = sum(res["objects"] for res in out["universes"].values())
    return out
