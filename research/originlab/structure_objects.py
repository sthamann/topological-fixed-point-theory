"""Generated mathematical objects with computed features for the structure scenario search.

Every function gen_*() yields objects {id, fam, name, feats}; feats is a set of
(kind, value, role) with kind 'int', 'frac' or 'type'. An object may name
other objects as its parts when it contains them by construction. No generator reads the
Atlas or any target. A new family is one more gen_*() function in this file.
"""
import itertools, math
from collections import Counter
from fractions import Fraction as Fr
from functools import reduce


def Z(v, w): return ("int", int(v), w)
def B(v, w):
    v = Fr(v); return ("int", int(v), w) if v.denominator == 1 else ("frac", f"{v.numerator}/{v.denominator}", w)
def T(v, w): return ("type", str(v), w)
def obj(id_, fam, name, feats): return {"id": id_, "fam": fam, "name": name, "feats": {f for f in feats if f is not None}}

# ---------- ADE Bausteine ----------
ROOTS = lambda t, n: {"A": n * (n + 1), "D": 2 * n * (n - 1), "E": {6: 72, 7: 126, 8: 240}.get(n, 0)}[t]
DET = lambda t, n: {"A": n + 1, "D": 4, "E": {6: 3, 7: 2, 8: 1}.get(n, 0)}[t]
COX = lambda t, n: {"A": n + 1, "D": 2 * n - 2, "E": {6: 12, 7: 18, 8: 30}.get(n, 0)}[t]
def graph(t, n, ext=False):
    if t == "A":
        e = [(i, i + 1) for i in range(1, n)]
        if ext: e += [(0, 1)] + ([(0, n)] if n > 1 else [])
    elif t == "D":
        e = [(i, i + 1) for i in range(1, n - 2)] + [(n - 2, n - 1), (n - 2, n)]
        if ext: e += [(0, 2)]
    else:
        chain = [1, 3, 4, 5, 6, 7, 8][: n - 1]; e = list(zip(chain, chain[1:])) + [(2, 4)]
        if ext: e += [(0, {6: 2, 7: 1, 8: 8}[n])]
    return set(range(0 if ext else 1, n + 1)), e
def classify(nodes, edges):
    out = []; nodes = set(nodes); adj = {v: set() for v in nodes}
    for a, b in edges:
        if a in nodes and b in nodes: adj[a].add(b); adj[b].add(a)
    seen = set()
    for s in nodes:
        if s in seen: continue
        comp = {s}; stack = [s]
        while stack:
            for w in adj[stack.pop()]:
                if w not in comp: comp.add(w); stack.append(w)
        seen |= comp; m = len(comp); br = [v for v in comp if len(adj[v]) == 3]
        if not br: out.append(("A", m)); continue
        arms = []
        for w in adj[br[0]]:
            k, prev, cur = 1, br[0], w
            while True:
                nxt = [x for x in adj[cur] if x != prev]
                if not nxt: break
                prev, cur = cur, nxt[0]; k += 1
            arms.append(k)
        arms.sort(); out.append(("D", m) if arms[:2] == [1, 1] else ("E", m))
    return tuple(sorted(out))
def subsystems(parent):
    start = (parent,); seen = {start}; todo = [start]
    while todo:
        cur = todo.pop()
        for i, (t, n) in enumerate(cur):
            rest = cur[:i] + cur[i + 1:]
            for ext in (False, True):
                nodes, edges = graph(t, n, ext)
                for d in nodes:
                    new = tuple(sorted(rest + classify(nodes - {d}, edges)))
                    if new and new not in seen: seen.add(new); todo.append(new)
    return seen
def tname(c): return " + ".join(f"{t}{n}" for t, n in sorted(c, key=lambda x: (-x[1], x[0])))
def _e8_classes():
    """Sizes of the glue classes of E8 roots for one explicit embedding of every full-rank subsystem."""
    import numpy as np
    simple = [[1, -1, 0, 0, 0, 0, 0, 0], [0, 1, -1, 0, 0, 0, 0, 0], [0, 0, 1, -1, 0, 0, 0, 0], [0, 0, 0, 1, -1, 0, 0, 0],
              [0, 0, 0, 0, 1, -1, 0, 0], [0, 0, 0, 0, 0, 1, -1, 0], [0, 0, 0, 0, 0, 1, 1, 0], [-.5]*8]
    roots = [tuple(v) for v in itertools.product((-1, 0, 1), repeat=8) if sum(x*x for x in v) == 2]
    roots += [tuple(s/2 for s in v) for v in itertools.product((-1, 1), repeat=8) if sum(x < 0 for x in v) % 2 == 0]
    R = np.array(roots)
    def kind(basis):
        B = np.array(basis); G = np.rint(B@B.T).astype(int)
        return classify(range(len(basis)), [(i, j) for i in range(len(basis)) for j in range(i) if G[i, j] == -1])
    def components(basis):
        B = np.array(basis); G = np.rint(B@B.T).astype(int); left, out = set(range(len(basis))), []
        while left:
            comp, stack = set(), [left.pop()]
            while stack:
                i = stack.pop(); comp.add(i)
                for j in list(left):
                    if G[i, j]: left.discard(j); stack.append(j)
            out.append(sorted(comp))
        return out
    def highest(basis, comp):
        B = np.array([basis[i] for i in comp]); coeff = R@B.T@np.linalg.inv(B@B.T)
        inside = (np.all(np.abs(coeff@B-R) < 1e-9, axis=1) & np.all(coeff > -1e-9, axis=1)
                  & np.all(np.abs(coeff-np.rint(coeff)) < 1e-9, axis=1))   # a root of this component, not of a glue class
        return tuple(R[np.argmax(np.where(inside, coeff.sum(1), -1))])
    start = tuple(map(tuple, simple)); seen = {kind(start): start}; todo = [start]
    while todo:
        basis = todo.pop()
        for comp in components(basis):
            extended = list(basis)+[tuple(-x for x in highest(basis, comp))]
            for drop in comp+[len(basis)]:
                new = tuple(v for i, v in enumerate(extended) if i != drop)
                if kind(new) not in seen: seen[kind(new)] = new; todo.append(new)
    sizes = {}
    for name, basis in seen.items():
        B = np.array(basis); coeff = R@B.T@np.linalg.inv(B@B.T)
        cls = Counter(tuple(np.rint((c-np.floor(c+1e-9))*240).astype(int) % 240) for c in coeff)
        sizes[name] = sorted(n for key, n in cls.items() if any(key))
    return sizes


def gen_teilgitter():
    classes = _e8_classes()
    for parent in (("E", 8), ("E", 7), ("E", 6)):
        prank = parent[1]
        for c in subsystems(parent):
            rank = sum(n for _, n in c)
            if rank < prank - 1: continue
            roots = sum(ROOTS(t, n) for t, n in c); det = reduce(lambda a, b: a * b, (DET(t, n) for t, n in c), 1)
            f = [T(f"{t}{n}", "component") for t, n in c] + [T(f"E{prank}", "parent lattice"), Z(roots, "roots"), Z(ROOTS(*parent) - roots, "roots outside"), Z(rank, "rank"), Z(roots + rank, "dimension"), Z(len(c), "number of components")]
            for t, n in c: f += [Z(n, f"rank {t}{n}"), Z(ROOTS(t, n), f"roots {t}{n}"), Z(ROOTS(t, n) + n, f"dimension {t}{n}"), Z(COX(t, n), f"Coxeter number {t}{n}"), Z(DET(t, n), f"determinant {t}{n}")]
            for t, n in c:
                if t == "A": f += [B(Fr(k * (n + 1 - k), n + 1), f"glue norm A{n} class {k}") for k in range(1, n // 2 + 2) if k <= n]
                if t == "D": f += [B(Fr(n, 4), f"glue norm D{n} spinor"), B(Fr(n, 8), f"spin field weight D{n}")]
                if t == "A" and n % 2 == 1: f.append(B(Fr(n, 8), f"weight A{n} middle class/2"))
            if rank == prank:
                idx = math.isqrt(det // DET(*parent))
                if idx * idx * DET(*parent) == det:
                    f.append(Z(idx, "index in parent lattice"))
                    if len(c) == 2 and idx > 1 and {DET(t, n) for t, n in c} == {idx}: f.append(T(f"Z{idx}", "cyclic glue of two components"))
            else: f.append(Z(1, "free direction"))
            if prank == 8 and rank == 8:
                f += [Z(n, "roots in one glue class") for n in classes.get(c, [])]
            yield obj(f"lattice:{tname(c)}<E{prank}", "Sublattice", f"{tname(c)} in E{prank}", f)

# ---------- Spiegelungsgruppen (Grade nach Shephard und Todd) ----------
DEG = {"G4": [4, 6], "G5": [6, 12], "G6": [4, 12], "G7": [12, 12], "G8": [8, 12], "G9": [8, 24], "G10": [12, 24], "G11": [24, 24], "G12": [6, 8], "G13": [8, 12], "G14": [6, 24], "G15": [12, 24],
       "G16": [20, 30], "G17": [20, 60], "G18": [30, 60], "G19": [60, 60], "G20": [12, 30], "G21": [12, 60], "G22": [12, 20], "G23 (H3)": [2, 6, 10], "G24": [4, 6, 14], "G25": [6, 9, 12], "G26": [6, 12, 18],
       "G27": [6, 12, 30], "G28 (F4)": [2, 6, 8, 12], "G29": [4, 8, 12, 20], "G30 (H4)": [2, 12, 20, 30], "G31": [8, 12, 20, 24], "G32": [12, 18, 24, 30], "G33": [4, 6, 10, 12, 18], "G34": [6, 12, 18, 24, 30, 42],
       "G35 (E6)": [2, 5, 6, 8, 9, 12], "G36 (E7)": [2, 6, 8, 10, 12, 14, 18], "G37 (E8)": [2, 8, 12, 14, 18, 20, 24, 30]}
def gen_spiegelungsgruppen():
    table = dict(DEG)
    for n in range(1, 9): table[f"W(A{n})"] = list(range(2, n + 2))
    for n in range(2, 9): table[f"W(B{n})"] = list(range(2, 2 * n + 1, 2))
    for n in range(4, 9): table[f"W(D{n})"] = sorted(list(range(2, 2 * n - 1, 2)) + [n])
    for m in (5, 6, 8, 10, 12): table[f"I2({m})"] = [2, m]
    for name, d in table.items():
        order = reduce(lambda a, b: a * b, d, 1)
        f = [Z(order, "order"), Z(len(d), "rank"), Z(sum(x - 1 for x in d), "number of reflections"), Z(reduce(math.gcd, d), "order of the centre"), Z(max(d), "highest degree")] + [Z(x, "degree") for x in d]
        yield obj(f"reflection:{name}", "Reflection group", f"Reflection group {name}", f)

# ---------- endliche Drehgruppen, McKay ----------
def gen_mckay():
    data = {"binary tetrahedral (E6)": (24, [1, 1, 1, 2, 2, 2, 3], (2, 3, 3)), "binary octahedral (E7)": (48, [1, 1, 2, 2, 2, 3, 3, 4], (2, 3, 4)), "binary icosahedral (E8)": (120, [1, 2, 2, 3, 3, 4, 4, 5, 6], (2, 3, 5))}
    for n in range(2, 13): data[f"cyclic Z{n} (A{n-1})"] = (n, [1] * n, (n, n))
    for n in range(2, 9): data[f"binary dihedral Dic{n} (D{n+2})"] = (4 * n, [1, 1, 1, 1] + [2] * (n - 1), (2, 2, n))
    for name, (order, marks, sig) in data.items():
        f = [Z(order, "order"), Z(len(marks), "number of representations"), Z(sum(marks), "Coxeter number (sum of marks)"), Z(max(marks), "largest mark")] + [Z(m, "mark") for m in set(marks) if m > 1] + [Z(s, "rotation order") for s in sig]
        if order % 2 == 0 and "cyclic" not in name: f.append(Z(order // 2, "order of the rotation group"))
        yield obj(f"mckay:{name}", "Rotation group", name, f)

# ---------- Linsenräume ----------
def gen_linsen():
    for p in range(2, 17):
        for k in range(0, p // 2 + 1):
            eta = Fr(p * p - 1, 6 * p) - Fr(k * (p - k), p)
            f = [Z(p, "order of the cut"), B(abs(eta), "boundary phase eta") if eta else Z(0, "boundary phase vanishes"), B(Fr(k * (p - k), p), "glue norm of the winding") if k else None, B(Fr(p * p - 1, 6 * p), "eta without winding")]
            if k: f.append(Z(k, "winding"))
            yield obj(f"lens:S3/Z{p},k={k}", "Lens space", f"Lens space S³/Z{p} with winding {k}", f)

# ---------- Modulkurven ----------
def _pf(n): return {p for p in range(2, n + 1) if n % p == 0 and all(p % q for q in range(2, int(p**.5) + 1))}
def _phi(n): return sum(1 for a in range(1, n + 1) if math.gcd(a, n) == 1)
def gen_modulkurven():
    for N in range(1, 26):
        mu = Fr(N) * reduce(lambda a, b: a * b, (1 + Fr(1, p) for p in _pf(N)), 1); cusps = sum(_phi(math.gcd(d, N // d)) for d in range(1, N + 1) if N % d == 0)
        def leg(a, p): return 0 if a % p == 0 else (1 if pow(a % p, (p - 1) // 2, p) == 1 else -1)
        nu2 = 0 if N % 4 == 0 else reduce(lambda a, b: a * b, (1 + (leg(-1, p) if p > 2 else 0) for p in _pf(N)), 1)
        nu3 = 0 if N % 9 == 0 else reduce(lambda a, b: a * b, (1 + (leg(-3, p) if p > 3 else (0 if p == 3 else -1)) for p in _pf(N)), 1)
        g = 1 + mu / 12 - Fr(nu2, 4) - Fr(nu3, 3) - Fr(cusps, 2)
        yield obj(f"modular:X0({N})", "Modular curve", f"Modular curve X₀({N})", [Z(N, "level"), Z(mu, "Index"), Z(cusps, "cusps"), Z(g, "genus"), Z(nu2, "points of order 2"), Z(nu3, "points of order 3")])
        if N >= 3:
            mu = Fr(N**3, 2) * reduce(lambda a, b: a * b, (1 - Fr(1, p * p) for p in _pf(N)), 1); cusps = mu / N; g = 1 + mu / 12 - cusps / 2
            yield obj(f"modular:X({N})", "Modular curve", f"Modular curve X({N})", [Z(N, "level"), Z(mu, "index = order of the symmetry group"), Z(cusps, "cusps"), Z(g, "genus"), Z(N, "width of each cusp")])
        if N >= 5:
            mu = Fr(N * N, 2) * reduce(lambda a, b: a * b, (1 - Fr(1, p * p) for p in _pf(N)), 1); cusps = Fr(sum(_phi(d) * _phi(N // d) for d in range(1, N + 1) if N % d == 0), 2); g = 1 + mu / 12 - cusps / 2
            yield obj(f"modular:X1({N})", "Modular curve", f"Modular curve X₁({N})", [Z(N, "level"), Z(mu, "Index"), Z(cusps, "cusps"), Z(g, "genus"), T(f"Z{N}", "torsion of the curve family")])

# ---------- elliptische Flächen mit endlich vielen Schnitten (Miranda und Persson) ----------
FIB = {"I": lambda n: ("A", n - 1) if n > 1 else None, "I*": lambda n: ("D", n + 4), "II*": lambda n: ("E", 8), "III*": lambda n: ("E", 7), "IV*": lambda n: ("E", 6), "III": lambda n: ("A", 1), "IV": lambda n: ("A", 2), "II": lambda n: None}
EUL = {"I": lambda n: n, "I*": lambda n: n + 6, "II*": lambda n: 10, "III*": lambda n: 9, "IV*": lambda n: 8, "III": lambda n: 3, "IV": lambda n: 4, "II": lambda n: 2}
SURF = {"X22": ([("II", 0), ("II*", 0)], [1]), "X33": ([("III", 0), ("III*", 0)], [2]), "X44": ([("IV", 0), ("IV*", 0)], [3]), "X11": ([("I*", 0), ("I*", 0)], [2, 2]), "X211": ([("II*", 0), ("I", 1), ("I", 1)], [1]),
        "X411": ([("I*", 4), ("I", 1), ("I", 1)], [2]), "X9111": ([("I", 9), ("I", 1), ("I", 1), ("I", 1)], [3]), "X8211": ([("I", 8), ("I", 2), ("I", 1), ("I", 1)], [4]), "X6321": ([("I", 6), ("I", 3), ("I", 2), ("I", 1)], [6]),
        "X5511": ([("I", 5), ("I", 5), ("I", 1), ("I", 1)], [5]), "X4422": ([("I", 4), ("I", 4), ("I", 2), ("I", 2)], [4, 2]), "X3333": ([("I", 3)] * 4, [3, 3]), "X321": ([("III*", 0), ("I", 2), ("I", 1)], [2]),
        "X431": ([("IV*", 0), ("I", 3), ("I", 1)], [3]), "X222": ([("I*", 2), ("I", 2), ("I", 2)], [2, 2]), "X141": ([("I*", 1), ("I", 4), ("I", 1)], [4])}
def gen_flaechen():
    for name, (fibs, tors) in SURF.items():
        lat = [FIB[t](n) for t, n in fibs]; lat = [x for x in lat if x]; tor = reduce(lambda a, b: a * b, tors, 1)
        f = [Z(len(fibs), "exceptional fibres"), Z(tor, "order of the torsion"), Z(sum(n for _, n in lat), "rank of fibre lattices"), Z(sum(EUL[t](n) for t, n in fibs), "Euler number")] + [T(f"{t}{n}", "fibre lattice") for t, n in lat]
        f += [Z(n, "width of a fibre") for t, n in fibs if t == "I" and n > 1]
        if len(tors) == 1 and tor > 1: f.append(T(f"Z{tor}", "cyclic torsion"))
        if sum(n for _, n in lat) == 8: f.append(T("E8", "fibre lattices fill E8"))
        yield obj(f"surface:{name}", "Elliptic surface", f"elliptic surface {name} (" + ", ".join(t.replace("I", f"I{n}", 1) if t in ("I", "I*") else t for t, n in fibs) + ")", f)

# ---------- Codes ----------
def gen_codes():
    for m in range(2, 8):
        for r in range(0, m):
            n = 2**m; k = sum(math.comb(m, i) for i in range(r + 1)); d = 2**(m - r)
            f = [Z(n, "length"), Z(k, "dimension"), Z(d, "distance"), B(Fr(k, n), "rate"), Z(m, "number of variables")]
            if 2 * r == m - 1: f += [T("self-dual", "code"), Z(m, "self-dual at m")]
            yield obj(f"code:RM({r},{m})", "code", f"Reed-Muller code RM({r},{m}) [{n},{k},{d}]", f)
    for name, n, k, d in (("Hamming [7,4,3]", 7, 4, 3), ("extended Hamming [8,4,4]", 8, 4, 4), ("Golay [23,12,7]", 23, 12, 7), ("extended Golay [24,12,8]", 24, 12, 8), ("Hexacode [6,3,4] over F4", 6, 3, 4), ("ternary Golay [12,6,6]", 12, 6, 6)):
        yield obj(f"code:{name}", "code", name, [Z(n, "length"), Z(k, "dimension"), Z(d, "distance"), B(Fr(k, n), "rate")] + ([T("self-dual", "code")] if 2 * k == n else []))
    for name, n, k, d in (("[[5,1,3]]", 5, 1, 3), ("[[7,1,3]] Steane", 7, 1, 3), ("[[9,1,3]] Shor", 9, 1, 3), ("[[4,2,2]]", 4, 2, 2), ("[[8,3,3]]", 8, 3, 3), ("[[15,1,3]]", 15, 1, 3), ("[[6,4,2]]", 6, 4, 2), ("[[16,10,3]]", 16, 10, 3), ("[[8,0,4]]", 8, 0, 4)):
        yield obj(f"qcode:{name}", "Quantum code", f"Quantum code {name}", [Z(n, "qubits"), Z(k, "logical qubits") if k else None, Z(d, "distance"), Z(n - k, "stabilizers"), Z(2**k, "code dimension") if k else None, B(Fr(k, n), "rate") if k else None])

# ---------- Qubits, Stabilisatorzustände, Clifford ----------
def gen_qudits():
    for d in (2, 3, 4, 5):
        for n in range(1, 6 if d == 2 else 3):
            D = d**n; f = [Z(D, "dimension"), Z(D * D - 1, "Pauli directions"), Z(n, "number of systems"), Z(d, "levels per system")]
            if d in (2, 3, 5):
                st = D * reduce(lambda a, b: a * b, (d**k + 1 for k in range(1, n + 1)), 1); f.append(Z(st, "stabilizer states"))
                cl = d**(n * n + 2 * n) * reduce(lambda a, b: a * b, (d**(2 * j) - 1 for j in range(1, n + 1)), 1); f.append(Z(cl, "order of the Clifford group modulo phases"))
                f.append(Z(D + 1, "unbiased bases")); f.append(Z(D * (D + 1), "rays in unbiased bases"))
            for kk in (2, 3, 4): f.append(Z(math.comb(D + kk - 1, kk), f"symmetric part of {kk} copies"))
            yield obj(f"qudit:{n}x{d}", "qubits", f"{n} system(s) with {d} levels", f)

# ---------- Anyonen und Randphasen ----------
def gen_anyonen():
    for nu in range(1, 17):
        yield obj(f"kitaev:{nu}", "Kitaev level", f"Kitaev level {nu} of sixteen", [Z(nu, "level"), B(Fr(nu, 2), "central charge"), B(Fr(nu, 16), "vortex weight"), Z(4 if nu % 2 == 0 else 3, "number of sectors"), T("Z4" if nu % 4 == 2 else ("Z2xZ2" if nu % 4 == 0 else "Ising"), "fusion")])
    for k in range(1, 13):
        f = [Z(k, "level"), B(Fr(3 * k, k + 2), "central charge"), Z(k + 1, "number of sectors"), Z(k + 2, "denominator")] + [B(Fr(j * (j + 2), 4 * (k + 2)), "weight") for j in range(1, k + 1)]
        yield obj(f"su2:{k}", "Anyons", f"SU(2) at level {k}", f)
    for n in range(2, 13):
        for k in range(1, 2 * n, 2 if n % 2 == 0 else 1):
            if n % 2 == 1 and k % 2: continue
            hs = {Fr(k * a * a, 2 * n) % 1 for a in range(1, n)}
            yield obj(f"zn:Z{n},k={k}", "Glue group", f"Glue group Z{n} with form {k}a²/{2*n}", [Z(n, "order"), T(f"Z{n}", "Glue group")] + [B(h, "weight") for h in hs if h])

# ---------- flache und runde Kegelflächen ----------
def gen_kegel():
    for sig in itertools.combinations_with_replacement(range(2, 9), 3):
        chi = 2 - sum(1 - Fr(1, m) for m in sig)
        if chi > 0: yield obj(f"cone:{sig}", "Cone surface", f"round sphere with cone points {sig}", [Z(2 / chi, "order of the rotation group"), Z(4 / chi, "order of the double group")] + [Z(m, "cone order") for m in sig])
        if chi == 0: yield obj(f"cone:{sig}", "Cone surface", f"flat sphere with cone points {sig}", [Z(3, "cone points"), Z(max(sig), "highest cone order"), T("flat", "geometry")] + [Z(m, "cone order") for m in sig])
    yield obj("cone:(2,2,2,2)", "Cone surface", "flat sphere with four half-turn points", [Z(4, "cone points"), Z(2, "cone order"), T("flat", "geometry"), Z(1, "free position parameter"), T("Z2", "rotation")])

def gen_polytope():
    data = {"Simplex 5-cell": ([5, 10, 10, 5], 120), "tesseract": ([16, 32, 24, 8], 384), "16-cell": ([8, 24, 32, 16], 384), "24-cell": ([24, 96, 96, 24], 1152), "120-cell": ([600, 1200, 720, 120], 14400), "600-cell": ([120, 720, 1200, 600], 14400),
            "icosahedron": ([12, 30, 20], 120), "dodecahedron": ([20, 30, 12], 120), "octahedron": ([6, 12, 8], 48), "cube": ([8, 12, 6], 48), "tetrahedron": ([4, 6, 4], 24),
            "Gosset 4_21": ([240, 6720, 60480], 696729600), "Gosset 3_21": ([56, 756, 4032], 2903040), "Gosset 2_21": ([27, 216, 720], 51840), "demicube 5": ([16, 80, 160], 1920), "cross polytope 5": ([10, 40, 80, 80, 32], 3840)}
    for name, (f, sym) in data.items():
        yield obj(f"polytope:{name}", "Polytope", f"Polytope {name}", [Z(f[0], "vertices"), Z(f[1], "edges"), Z(sym, "order of the symmetry group"), Z(2 * f[1] // f[0], "neighbours per vertex"), Z(len(f), "dimension")] + [Z(x, "face count") for x in f[2:]])

def gen_darstellungen():
    reps = {"A4 = SU(5)": [5, 10, 24, 15, 40, 45, 50], "D5 = SO(10)": [10, 16, 45, 54, 120, 126, 144, 210], "E6": [27, 78, 351, 650], "E7": [56, 133, 912, 1539], "E8": [248, 3875, 30380], "A3 = SU(4)": [4, 6, 10, 15, 20, 36, 45],
            "A2 = SU(3)": [3, 6, 8, 10, 15, 27], "D4 = SO(8)": [8, 28, 35, 56], "A1 = SU(2)": [2, 3, 4, 5], "G2": [7, 14, 27, 64, 77], "F4": [26, 52, 273, 324], "A5 = SU(6)": [6, 15, 20, 35, 21], "B4 = SO(9)": [9, 16, 36, 84, 126], "C3 = Sp(6)": [6, 14, 21]}
    for name, d in reps.items():
        yield obj(f"rep:{name}", "Representations", f"small representations of {name}", [Z(x, "dimension of a representation") for x in d] + [T(name.split(" ")[0], "algebra")])

def gen_modulformen():
    for k in range(2, 26, 2):
        dim = k // 12 + (0 if k % 12 == 2 else 1)
        yield obj(f"modform:k={k}", "Modular forms", f"modular forms of weight {k} for the full group", [Z(k, "weight"), Z(dim, "dimension"), B(Fr(k, 12), "degree of the bundle"), Z(dim - 1, "cusp forms") if dim > 1 else None])
    for m in range(3, 10): yield obj(f"cft:M({m+1},{m})", "Boundary theory", f"minimal model M({m+1},{m})", [B(1 - Fr(6, m * (m + 1)), "central charge"), Z(m * (m - 1) // 2, "number of sectors"), Z(m, "level")])
    for k in range(2, 13): yield obj(f"cft:PF{k}", "Boundary theory", f"parafermions Z{k}", [B(Fr(2 * (k - 1), k + 2), "central charge"), Z(k * (k + 1) // 2, "number of sectors"), Z(k, "order"), T(f"Z{k}", "symmetry")])
    for N in range(3, 6):
        for k in range(1, 7): yield obj(f"cft:SU({N})_{k}", "Boundary theory", f"SU({N}) at level {k}", [B(Fr(k * (N * N - 1), k + N), "central charge"), Z(math.comb(N + k - 1, N - 1), "number of sectors"), Z(k, "level"), Z(k + N, "denominator"), T(f"A{N-1}", "algebra")])
    for name, t, c, n in (("SO(10) level 1", "D5", 5, 4), ("E6 level 1", "E6", 6, 3), ("E7 level 1", "E7", 7, 2), ("E8 level 1", "E8", 8, 1)):
        yield obj(f"cft:{name}", "Boundary theory", f"Boundary theory {name}", [Z(c, "central charge"), Z(n, "number of sectors"), T(t, "algebra")])


# ---------- fermion modes, graphs on subsets, repeated steps, matter content ----------
def gen_fermions():
    for n in range(1, 9):
        f = [Z(n, "modes"), Z(2*n, "Majorana operators"), Z(2**n, "Fock states"), Z(4**n, "operator dimension"),
             Z(n*(2*n-1), "dimension of so(2n)"), B(Fr(n, 4), "spinor glue norm"), B(Fr(n, 8), "spin field weight"),
             Z(n*n, "operators on the modes"), B(Fr(1, n), "share of one place") if n > 1 else None, B(Fr(n-1, n), "share of all other places") if n > 2 else None]
        if n > 1: f += [Z(2**(n-1), "even Fock states"), Z(2*n*(n-1), "roots of D_n"), Z(math.comb(n, 2), "mode pairs")]
        if n == 8: f += [Z(240, "roots of D8 with the even spinor"), Z(248, "so(16) with the even spinor"), T("E8", "even unimodular closure")]
        row = obj(f"fermion:{n}", "Fermion modes", f"{n} fermion mode(s)", f)
        yield {**row, "parts": [f"moves:{n}x2"]} if 2 <= n <= 6 else row   # the Fock space of n modes is the cube of n bits
    for a in range(2, 7):
        for b in range(1, a+1):
            if a+b > 8: continue
            n = a+b
            f = [Z(n, "modes"), Z(a, "modes of the first block"), Z(b, "modes of the second block"), Z(4*a*b, "mixed roots (vector, vector)"),
                 Z(2**(a-1)*2**(b-1)*2, "even Fock states"), Z(2**(a-1), "half spinor of the first block"), Z(2**(b-1), "half spinor of the second block"),
                 Z(a*(2*a-1), "dimension of the first so"), Z(b*(2*b-1), "dimension of the second so"), Z(2*a, "vector of the first block"), Z(2*b, "vector of the second block"),
                 B(Fr(a, 4), "spinor glue norm of the first block"), B(Fr(b, 4), "spinor glue norm of the second block"), T(f"D{a}", "first block")]
            if b == 3: f.append(T("A3", "second block, D3 = A3"))
            if n == 8: f += [Z(240, "roots of D8 with the even spinor"), Z(248, "so(16) with the even spinor"), T("E8", "even unimodular closure")]
            block = lambda m: {1: [], 2: [("A", 1), ("A", 1)], 3: [("A", 3)]}.get(m, [("D", m)])
            parts = [f"fermion:{a}", f"fermion:{b}", f"fermion:{n}"] + [f"graph:K({m},2)" for m in (a, b, n) if m >= 5]
            if n == 8: parts.append(f"lattice:{tname(tuple(sorted(block(a)+block(b))))}<E8")
            yield {**obj(f"fermion:{a}+{b}", "Fermion modes", f"{n} fermion modes split {a} + {b}", f), "parts": parts}


def gen_graphs():
    for n in range(4, 10):
        for k in (2, 3):
            if n < 2*k+1: continue
            v, d = math.comb(n, k), math.comb(n-k, k)
            eig = [((-1)**i*math.comb(n-k-i, k-i), math.comb(n, i)-(math.comb(n, i-1) if i else 0)) for i in range(k+1)]
            f = [Z(v, "vertices"), Z(d, "neighbours per vertex"), Z(v*d//2, "edges"), Z(n, "points"), Z(math.factorial(n), "order of the symmetry group"), Z(v*v, "operator dimension")]
            for lam, mult in eig[1:]:
                f += [Z(abs(lam), "eigenvalue"), Z(mult, "multiplicity"), B(Fr(abs(lam), d), "normalised eigenvalue"), B(Fr(lam*lam, d*d), "squared normalised eigenvalue")]
            if k == 2:   # operators on pairs, and pairs of pairs, split under the symmetry of the n places
                f += [Z(n*n-1, "operators on pairs: adjoint part"), Z(v*v-n*n, "operators on pairs: remaining part"),
                      Z(math.comb(n, 4), "two pairs: four distinct places"), Z(3*math.comb(n+1, 4), "two pairs: mixed symmetry"),
                      Z(n*n*(n*n-1)//12, "two pairs: window symmetry")]
            yield obj(f"graph:K({n},{k})", "Subset graph", f"disjointness graph of {k}-subsets of {n} points" + (" (Petersen)" if (n, k) == (5, 2) else ""), f)
    for n in range(2, 9):
        f = [Z(2**n, "vertices"), Z(n, "neighbours per vertex"), Z(n*2**(n-1), "edges"), Z(2**n*math.factorial(n), "order of the symmetry group")] + [B(Fr(abs(n-2*i), n), "normalised eigenvalue") for i in range(1, n) if n != 2*i] + [Z(math.comb(n, i), "multiplicity") for i in range(1, n)]
        yield obj(f"graph:Q{n}", "Subset graph", f"cube graph of dimension {n}", f)


def gen_moves():
    """d places with q values each; a step draws one of the d*q moves 'set place i to value v' uniformly. One hand is d*q steps."""
    for d in range(2, 7):
        for q in range(2, 6):
            moves = d*q
            f = [Z(d, "places"), Z(q, "values per place"), Z(moves, "moves"), Z(q**d, "states"), Z(d+1, "levels with the wall"),
                 B(Fr(1, q), "probability to stay"), B(Fr(1, moves), "probability per neighbour"), B(Fr(q-1, q)-Fr(1, moves), "loss from a pair"),
                 Z(d**moves, "histories of one hand")]
            for n in range(1, d):
                f += [B(1-Fr(n, d), "rate of a level"), B((1-Fr(n, d))**moves, "rate of a level after one hand"),
                      Z((d-n)**moves, "histories of one hand that avoid a level"), Z(math.comb(d, n)*(q-1)**n, "multiplicity of a level")]
            if q == 2: f += [Z(2*d, "Majorana operators"), Z(d*(2*d-1), "dimension of so(2d)"), Z(2, "slope of the linear clock")]
            yield obj(f"moves:{d}x{q}", "Move rule", f"{d} places with {q} values, one of {moves} moves drawn per step", f)


def gen_steps():
    for b in range(2, 7):
        for k in range(1, 9):
            f = [Z(b**k, "number of histories"), Z(k, "steps"), Z(b, "outcomes per step")] + [B(Fr(a, b)**k, "weight of a repeated outcome") for a in range(1, b) if math.gcd(a, b) == 1 or k > 1]
            yield obj(f"steps:{b}^{k}", "Repeated steps", f"{k} independent step(s) with {b} outcomes", f)


def gen_matter():
    for g in range(1, 6):
        for h in range(0, 4):
            b1, b2, b3 = Fr(4*g, 3)+Fr(h, 10), Fr(22, 3)-Fr(4*g, 3)-Fr(h, 6), 11-Fr(4*g, 3)
            f = [Z(g, "families"), Z(16*g, "Weyl fermions with a right-handed neutrino"), Z(15*g, "Weyl fermions"), Z(12, "gauge fields"), Z(32*g, "fermion states with antiparticles"),
                 B(abs(b1), "hypercharge running coefficient"), B(abs(b2), "weak running coefficient"), B(abs(b3), "colour running coefficient"), B(Fr(5, 3), "hypercharge normalisation"),
                 B(Fr(3, 8), "weak angle at unification"), B(Fr(10*g, 9)*3, "family trace of hypercharge squared"), Z(h, "Higgs doublets") if h else None]
            yield obj(f"matter:{g}g{h}h", "Matter content", f"{g} family(ies) and {h} Higgs doublet(s)", f)


def gen_nested_splits():
    """A block of a modes split a1 + a2 with its traceless charge, beside b further modes read as b copies."""
    for a1 in range(1, 6):
        for a2 in range(1, a1+1):
            a = a1+a2
            if not 3 <= a <= 6: continue
            y = [Fr(-1, a1)]*a1+[Fr(1, a2)]*a2          # traceless, one unit across the split
            square = sum(x*x for x in y)                # trace of the squared charge on the a places
            spinor = Fr(2**(a-1))*square/4              # the same trace on the even Fock space
            f = [Z(a, "modes"), Z(a1, "first part"), Z(a2, "second part"), Z(a1*a2, "mixed places"), Z(2*a1*a2, "mixed roots"),
                 B(square, "trace of the squared charge"), B(Fr(a1, a2), "ratio of the parts"), B(square*Fr(a1, a2), "charge trace times ratio"),
                 B(spinor, "charge trace on the even Fock space"), B(2*square, "charge normalisation"),
                 B(Fr(a1, a), "share of the first part"), B(Fr(a2, a), "share of the second part"),
                 Z(a1*a1, "operators on the first part"), Z(a2*a2, "operators on the second part"),
                 Z(a1*a1+a2*a2, "operators that keep the split"), Z(a1*a1+a2*a2-1, "traceless operators that keep the split")]
            charges = Counter(abs(sum(y[i] for i in S)) for r in range(0, a+1, 2) for S in itertools.combinations(range(a), r))
            f += [B(q, "charge of an even Fock state") for q in charges if q] + [Z(m, "even Fock states of one charge") for m in charges.values()]
            f += [B(abs(p-q), "charge difference of two places") for p in set(y) for q in set(y) if p != q]
            yield {**obj(f"nested:{a1}+{a2}", "Nested split", f"{a} modes split {a1} + {a2} with its charge", f), "parts": [f"fermion:{a1}+{a2}"]}
            for b in range(1, 9-a):
                g = [Z(b, "copies"), Z(2**(a-1)*b, "even Fock states in all copies"), Z(2**a*b, "Fock states in all copies")]
                for h in (0, 1, 2):   # the scalar sits on the second part with charge 1/a2
                    g.append(B((Fr(2, 3)*b*spinor+Fr(h, 3*a2))/(2*square), f"running coefficient of the charge with {h} scalar(s) on the second part"))
                yield {**obj(f"nested:({a1}+{a2})+{b}", "Nested split", f"{a+b} modes split ({a1} + {a2}) + {b} with its charge", f+g),
                       "parts": [f"fermion:{a}+{b}" if a >= b else f"fermion:{b}+{a}", f"fermion:{a1}+{a2}"]}
