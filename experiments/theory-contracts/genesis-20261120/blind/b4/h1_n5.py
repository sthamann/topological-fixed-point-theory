"""H1, n=5 completion attempt: magmas with left identity 0 (WLOG), LD backtracking plus
duplicate-row pruning (faithfulness is required by the rule). Time limit 170 s."""
import itertools, time, sys
import numpy as np
from h1_self_interpreter import maps, count_homs, canon
n = 5; R = range(n); TL = 170
t = [[-1]*n for _ in R]
for y in R: t[0][y] = y
cells = [(x, y) for x in range(1, n) for y in R]
def ok(x0, y0):
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
    if y0 == n-1 and any(t[x0] == t[k] for k in range(x0)): return False
    return True
F = maps(n); perms = list(itertools.permutations(R)); t0 = time.time()
cnt = 0; sols = set(); lab = 0
def rec(i):
    global cnt, lab
    if time.time() - t0 > TL: raise TimeoutError
    if i == len(cells):
        cnt += 1
        T = np.array(t)
        if count_homs(T, F).sum() == n:
            lab += 1; sols.add(canon(tuple(tuple(r) for r in t), perms))
        return
    x, y = cells[i]
    for v in R:
        t[x][y] = v
        if ok(x, y): rec(i+1)
    t[x][y] = -1
complete = True
try: rec(0)
except TimeoutError: complete = False
print(f"n=5 left-identity class, complete={complete}, faithful LD magmas visited={cnt}, RULE labelled={lab}, iso classes={len(sols)}, last row1 prefix reached={t[1]}, time={time.time()-t0:.0f}s")
for s in sorted(sols): print("  table:", s)
