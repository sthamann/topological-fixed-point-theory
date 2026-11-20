"""H1: magma (X,*) that is its own interpreter:  L: x -> (y -> x*y) is a BIJECTION X -> End(X,*).
Exhaustive over all magmas of order n<=4 (via left-distributivity backtracking, which is implied
by the rule), and over all magmas with a left identity for n=5,6 (WLOG under the rule: id in End).
Controls: (c1) left-distributive only; (c2) LD + faithful (L injective); (c3) Aut-version
(L bijection X->Aut(X)); (c4) RIGHT-handed mismatch: L bijection onto End of the opposite magma
is the same up to relabel, so instead we use 'anti-interpreter': L bijects onto the ANTI-endomorphisms."""
import itertools, sys, time
import numpy as np

def ld_magmas(n, row0_id, tlimit):
    t = [[-1]*n for _ in range(n)]
    cells = [(x,y) for x in range(n) for y in range(n)]
    if row0_id:
        for y in range(n): t[0][y] = y
        cells = [c for c in cells if c[0] != 0]
    t0 = time.time(); R = range(n)
    def ok():
        for x in R:
            tx = t[x]
            for y in R:
                xy = tx[y]
                if xy < 0: continue
                for z in R:
                    yz = t[y][z]; xz = tx[z]
                    if yz < 0 or xz < 0: continue
                    a = tx[yz]; b = t[xy][xz]
                    if a >= 0 and b >= 0 and a != b: return False
        return True
    def rec(i):
        if time.time() - t0 > tlimit: raise TimeoutError
        if i == len(cells):
            yield tuple(tuple(r) for r in t); return
        x, y = cells[i]
        for v in R:
            t[x][y] = v
            if ok(): yield from rec(i+1)
        t[x][y] = -1
    yield from rec(0)

def maps(n):
    return np.array(list(itertools.product(range(n), repeat=n)), dtype=np.int8)

def count_homs(T, F, anti=False):
    n = T.shape[0]; good = np.ones(len(F), bool)
    for x in range(n):
        for y in range(n):
            lhs = F[:, T[x, y]]
            rhs = T[F[:, y], F[:, x]] if anti else T[F[:, x], F[:, y]]
            good &= (lhs == rhs)
    return good

def canon(T, perms):
    n = len(T); best = None
    for p in perms:
        inv = [0]*n
        for i, pi in enumerate(p): inv[pi] = i
        c = tuple(tuple(p[T[inv[a]][inv[b]]] for b in range(n)) for a in range(n))
        if best is None or c < best: best = c
    return best

def run(n, row0_id, tlimit):
    F = maps(n); isperm = np.array([len(set(f)) == n for f in F.tolist()])
    perms = list(itertools.permutations(range(n)))
    cnt = dict(LD=0, LD_faithful=0, RULE_End=0, RULE_Aut=0, CTRL_antiEnd=0)
    sol_end, sol_aut = set(), set(); complete = True
    try:
        for t in ld_magmas(n, row0_id, tlimit):
            cnt['LD'] += 1
            if len(set(t)) < n: continue
            cnt['LD_faithful'] += 1
            T = np.array(t)
            g = count_homs(T, F)
            if g.sum() == n:
                cnt['RULE_End'] += 1; sol_end.add(canon(t, perms))
            if all(len(set(r)) == n for r in t) and (g & isperm).sum() == n:
                cnt['RULE_Aut'] += 1; sol_aut.add(canon(t, perms))
            ga = count_homs(T, F, anti=True)
            rows = set(t)
            if ga.sum() == n and all(tuple(f) in rows for f in F[ga].tolist()):
                cnt['CTRL_antiEnd'] += 1
    except TimeoutError:
        complete = False
    print(f"n={n} row0=id:{row0_id} complete:{complete} total_magmas_in_class={'n^(n^2)='+str(n**(n*n)) if not row0_id else 'n^(n(n-1))='+str(n**(n*(n-1)))}")
    print("   labelled counts:", cnt)
    print(f"   RULE_End iso classes: {len(sol_end)}   RULE_Aut iso classes: {len(sol_aut)}")
    for s in sorted(sol_end): print("   End-solution table:", s)
    for s in sorted(sol_aut)[:6]: print("   Aut-solution table:", s)
    sys.stdout.flush()
    return sol_end

if __name__ == "__main__":
    print("== full class: ALL magmas of order n (LD pruning is implied by the rule) ==")
    for n in (1, 2, 3, 4): run(n, False, 100)
    print("== class: all magmas with 0 a left identity (WLOG under the rule) ==")
    for n, tl in ((4, 60), (5, 60), (6, 45)): run(n, True, tl)
