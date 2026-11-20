"""H4: admissible histories. Local constraint = ternary relation R on alphabet A, |A|=n, such that
any two legs determine the third (no leg is 'time')  <=> Latin square / quasigroup z = x*y.
RULE = (S) leg symmetry: R invariant under all permutations of its 3 legs
     + (X) crossing: gluing two constraints along a hidden leg gives a 4-leg relation independent of
           the pairing:  (a*b)*c = (a*c)*b.
Exhaustive over ALL Latin squares n<=5; n=6..8 by backtracking with the rule as pruning.
Controls: S only, X only, associativity only."""
import itertools, time, sys

def latin(n, prune=None, tlimit=1e9):
    t = [[-1]*n for _ in range(n)]
    rowu = [set() for _ in range(n)]; colu = [set() for _ in range(n)]
    cells = [(x, y) for x in range(n) for y in range(n)]; t0 = time.time()
    def rec(i):
        if time.time() - t0 > tlimit: raise TimeoutError
        if i == len(cells):
            yield [r[:] for r in t]; return
        x, y = cells[i]
        for v in range(n):
            if v in rowu[x] or v in colu[y]: continue
            t[x][y] = v; rowu[x].add(v); colu[y].add(v)
            if prune is None or prune(t, n): yield from rec(i+1)
            rowu[x].discard(v); colu[y].discard(v)
        t[x][y] = -1
    yield from rec(0)

def S(t, n):   # totally symmetric: commutative and x*(x*y)=y
    return all(t[x][y] == t[y][x] and t[x][t[x][y]] == y for x in range(n) for y in range(n))
def X(t, n):
    R = range(n)
    return all(t[t[a][b]][c] == t[t[a][c]][b] for a in R for b in R for c in R)
def ASSOC(t, n):
    R = range(n)
    return all(t[t[a][b]][c] == t[a][t[b][c]] for a in R for b in R for c in R)
def partial_rule(t, n):
    R = range(n)
    for x in R:
        for y in R:
            v = t[x][y]
            if v < 0: continue
            if t[y][x] >= 0 and t[y][x] != v: return False
            if t[x][v] >= 0 and t[x][v] != y: return False
            for c in R:
                w = t[x][c]
                if w < 0: continue
                p, q = t[v][c], t[w][y]
                if p >= 0 and q >= 0 and p != q: return False
    return True
def partial_S(t, n):
    for x in range(n):
        for y in range(n):
            v = t[x][y]
            if v < 0: continue
            if t[y][x] >= 0 and t[y][x] != v: return False
            if t[x][v] >= 0 and t[x][v] != y: return False
    return True
def canon(t, n):
    best = None
    for p in itertools.permutations(range(n)):
        inv = [0]*n
        for i, pi in enumerate(p): inv[pi] = i
        c = tuple(tuple(p[t[inv[a]][inv[b]]] for b in range(n)) for a in range(n))
        if best is None or c < best: best = c
    return best

print("== ALL Latin squares, n<=5 ==")
for n in range(1, 6):
    c = dict(all=0, S=0, X=0, ASSOC=0, RULE=0); iso = {k: set() for k in ('S', 'X', 'RULE')}
    for t in latin(n):
        c['all'] += 1
        s, x = S(t, n), X(t, n)
        if s: c['S'] += 1; iso['S'].add(canon(t, n))
        if x: c['X'] += 1; iso['X'].add(canon(t, n))
        if ASSOC(t, n): c['ASSOC'] += 1
        if s and x: c['RULE'] += 1; iso['RULE'].add(canon(t, n))
    print(f"n={n} labelled {c}  iso classes: S={len(iso['S'])} X={len(iso['X'])} RULE={len(iso['RULE'])}")
    for s in sorted(iso['RULE']): print("    RULE table:", s)
    sys.stdout.flush()
print("== n=6..9: backtracking with pruning (RULE; control S only) ==")
for n in range(6, 10):
    out = []
    for name, pr, tl in (('RULE', partial_rule, 25), ('S', partial_S, 12)):
        k = 0; comp = True; first = None
        try:
            for t in latin(n, pr, tl):
                k += 1
                if first is None: first = t
        except TimeoutError: comp = False
        out.append(f"{name}: labelled={k}{'' if comp else '+ (TIMEOUT, lower bound)'}")
        if name == 'RULE' and first is not None:
            # check: is it Z2^k with x*y = x xor y after relabelling? test: exists e with e*x=x, x*x=e, assoc
            e = [a for a in range(n) if all(first[a][x] == x for x in range(n))]
            out.append(f"  [first RULE table: identity={e}, assoc={ASSOC(first, n)}, exponent2={bool(e) and all(first[x][x]==e[0] for x in range(n))}]")
    print(f"n={n} " + " ".join(out)); sys.stdout.flush()
