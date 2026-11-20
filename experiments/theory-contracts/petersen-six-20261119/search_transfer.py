"""Suche: Gibt es unter allen Operatoren auf den zehn Paaren, die die Teilung 3 + 2 respektieren, einen, der beim
Zusammenfassen zu drei Klassen den Originaltransfer (v486) ergibt? Alle Dreiteilungen, beide Lesarten des Verlusts.
Aufruf: python3 -B search_transfer.py"""
import itertools, json
import numpy as np
from scipy.optimize import linprog

F, W = ("f1", "f2", "f3"), ("w1", "w2")
places = F+W
pairs = [frozenset(p) for p in itertools.combinations(places, 2)]
idx = {p: i for i, p in enumerate(pairs)}
H = [dict(zip(places, f+w)) for f in itertools.permutations(F) for w in itertools.permutations(W)]
act = lambda g, p: frozenset(g[x] for x in p)

# Bahnen von H auf geordneten Paaren von Paaren = Basis aller Operatoren, die mit H vertauschen
orbit = {}
for i, p in enumerate(pairs):
    for j, q in enumerate(pairs):
        if (i, j) not in orbit:
            for g in H: orbit[(idx[act(g, p)], idx[act(g, q)])] = len(set(orbit.values())) if (i, j) not in orbit else orbit[(i, j)]
n = len(set(orbit.values()))
basis = np.zeros((n, 10, 10))
for (i, j), k in orbit.items(): basis[k, i, j] = 1
adjacent = np.array([[float(not p & q) for q in pairs] for p in pairs])
local = [k for k in range(n) if np.all(basis[k]*(1-adjacent-np.eye(10)) == 0)]
symmetric_pairs = [(k, next(m for m in range(n) if np.array_equal(basis[m], basis[k].T))) for k in range(n)]

LOSSY = np.array([[1, 0, 0], [0, 1/2, 1/6], [0, 1/6, 1/2]])
ABSORB = np.array([[1, 0, 0], [1/3, 1/2, 1/6], [1/3, 1/6, 1/2]])

def system(labels, target, extra=()):
    rows, rhs = [], []
    for i in range(10):
        for c in range(3):
            rows.append([basis[k, i][labels == c].sum() for k in range(n)]); rhs.append(target[labels[i], c])
    for row in extra: rows.append(row); rhs.append(0)
    return np.array(rows), np.array(rhs)

def solve(labels, target, extra=()):
    Aeq, b = system(labels, target, extra)
    x, res, rank, _ = np.linalg.lstsq(Aeq, b, rcond=None)
    if np.abs(Aeq@x-b).max() > 1e-9: return None
    lp = linprog(np.zeros(n), A_eq=Aeq, b_eq=b, bounds=[(0, None)]*n, method="highs")
    return {"free": int(n-rank), "nonnegative": bool(lp.status == 0), "x": lp.x if lp.status == 0 else x}

sym_rows = [[(1 if m == k else 0)-(1 if m == k2 else 0) for m in range(n)] for k, k2 in symmetric_pairs if k < k2]
loc_rows = [[1 if m == k else 0 for m in range(n)] for k in range(n) if k not in local]

relation = [0 if np.trace(basis[k]) else 1 if (basis[k]*adjacent).sum() else 2 for k in range(n)]
graph_rows = [[(1 if m == k else 0)-(1 if m == k2 else 0) for m in range(n)] for k in range(n) for k2 in range(k+1, n) if relation[k] == relation[k2]]

def stable(labels):
    classes = {frozenset(p for p, l in zip(pairs, labels) if l == c) for c in range(3)}
    return all({frozenset(act(g, p) for p in c) for c in classes} == classes for g in H)

def canonical(labels):
    best = None
    for g in H:
        img = [0]*10
        for i, p in enumerate(pairs): img[idx[act(g, p)]] = labels[i]
        for swap in ((0, 1, 2), (0, 2, 1)):
            key = tuple(swap[x] for x in img)
            if best is None or key < best: best = key
    return best

def describe(labels):
    return [sorted("".join(sorted(x[0]+x[1] for x in p)) for p, l in zip(pairs, labels) if l == c) for c in range(3)]

out = {"operators": n, "local_operators": len(local), "symmetric_operators": n-len(sym_rows), "targets": {}}
for name, target in (("verlust", LOSSY), ("aufgesogen", ABSORB)):
    seen, hits = set(), []
    for labels in itertools.product(range(3), repeat=10):
        if len(set(labels)) < 3: continue
        key = canonical(labels)
        if key in seen: continue
        seen.add(key); lab = np.array(labels)
        r = solve(lab, target)
        if not r: continue
        hit = {"classes": describe(labels), "stable": stable(labels), "free": r["free"], "nonnegative": r["nonnegative"]}
        for tag, extra in (("symmetric", sym_rows), ("local", loc_rows), ("local_symmetric", sym_rows+loc_rows), ("graph_only", graph_rows)):
            s = solve(lab, target, extra)
            hit[tag] = None if not s else {"free": s["free"], "nonnegative": s["nonnegative"]}
        hits.append(hit)
    out["targets"][name] = {"partitions_up_to_symmetry": len(seen), "feasible": len(hits), "hits": hits}
    print(name, "Teilungen bis auf Symmetrie:", len(seen), "lösbar:", len(hits))
    ok = lambda h, tag: bool(h[tag] and h[tag]["nonnegative"])
    counts = {"nonnegative": sum(h["nonnegative"] for h in hits), "stable": sum(h["stable"] and h["nonnegative"] for h in hits),
              **{tag: sum(ok(h, tag) for h in hits) for tag in ("symmetric", "local", "local_symmetric", "graph_only")},
              "stable_symmetric": sum(h["stable"] and ok(h, "symmetric") for h in hits)}
    out["targets"][name]["counts"] = counts; print("  ", counts)
    for h in hits:
        if h["stable"] and h["nonnegative"] or ok(h, "local_symmetric") or ok(h, "graph_only"): print("  ", h)
json.dump(out, open("results_transfer_search.json", "w"), indent=1, ensure_ascii=False)
print("Operatoren:", n, "lokal:", len(local), "symmetrisch:", n-len(sym_rows))
