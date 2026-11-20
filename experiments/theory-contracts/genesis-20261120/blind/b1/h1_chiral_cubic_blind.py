"""H1: exhaustive search for the smallest irreducible complex (non-self-dual) representation
of a simple compact Lie group whose cubic moment Tr rho(X)^3 vanishes identically.
Class: ALL simple types admitting non-self-dual irreps: A_n (n>=2), D_n (n odd >=5), E6.
(B,C,G2,F4,E7,E8,A1,D_even have -w0 = 1: every irrep self-dual -> excluded by theorem.)
Controls: (a) drop cubic condition, (b) drop irreducibility (abelian multisets, = hypothesis H4)."""
import sys, itertools, random
from fractions import Fraction
from math import gcd
sys.setrecursionlimit(10000)
DMAX = 1000

def cartan(typ, n):
    C = [[2 if i == j else 0 for j in range(n)] for i in range(n)]
    def link(i, j): C[i][j] = C[j][i] = -1
    if typ == 'A':
        for i in range(n-1): link(i, i+1)
    elif typ == 'D':
        for i in range(n-2): link(i, i+1)
        link(n-3, n-1)
    elif typ == 'E6':  # Bourbaki: 1-3-4-5-6, 2-4
        for a, b in [(0,2),(2,3),(3,4),(4,5),(1,3)]: link(a, b)
    return C

def pos_roots(C):
    n = len(C); roots = {tuple(1 if i == j else 0 for j in range(n)) for i in range(n)}
    frontier = set(roots)
    while frontier:
        new = set()
        for r in frontier:
            for i in range(n):
                # simply laced: r+alpha_i is a root iff <r,alpha_i^vee> = -1
                if sum(r[j]*C[j][i] for j in range(n)) == -1:
                    s = list(r); s[i] += 1; s = tuple(s)
                    if s not in roots: new.add(s)
        roots |= new; frontier = new
    return sorted(roots)

def weyl_dim(lab, roots):
    num = den = 1
    for r in roots:
        num *= sum(c*(l+1) for c, l in zip(r, lab)); den *= sum(r)
    return num // den

def enum_irreps(typ, n, dmax):
    C = cartan(typ, n); R = pos_roots(C); out = []
    def rec(i, lab):
        if i == n:
            out.append((tuple(lab), weyl_dim(lab, R))); return
        k = 0
        while True:
            lab[i] = k
            if weyl_dim(lab, R) > dmax: break
            rec(i+1, lab); k += 1
        lab[i] = 0
    rec(0, [0]*n); return out, len(R)

def dual(typ, n, lab):
    if typ == 'A': return tuple(reversed(lab))
    if typ == 'D':
        return lab if n % 2 == 0 else lab[:-2] + (lab[-1], lab[-2])
    if typ == 'E6': return (lab[5], lab[1], lab[4], lab[3], lab[2], lab[0])

def dimGL(part):
    m = len(part); num = den = 1
    for i in range(m):
        for j in range(i+1, m):
            num *= part[i]-part[j]+j-i; den *= j-i
    return num // den

def cubic_index_A(lab):
    """sum over weights of (w,H)^3, H=diag(1,..,1,-(N-1)), via GL(N)->GL(N-1) branching. Exact."""
    N = len(lab)+1
    lam = [sum(lab[i:]) for i in range(N-1)] + [0]
    tot = sum(lam); A = 0; d = 0
    for mu in itertools.product(*[range(lam[i+1], lam[i]+1) for i in range(N-1)]):
        dm = dimGL(mu); k = tot - sum(mu)
        A += dm*(tot - N*k)**3; d += dm
    Afund = (N-1) + (-(N-1))**3
    return Fraction(A, Afund), d

def weyl_orbit(C, lab):
    n = len(C); seen = {tuple(lab)}; fr = [tuple(lab)]
    while fr:
        nf = []
        for w in fr:
            for i in range(n):
                if w[i] != 0:
                    v = tuple(w[j] - w[i]*C[i][j] for j in range(n))
                    if v not in seen: seen.add(v); nf.append(v)
        fr = nf
    return seen

print("=== validation of cubic index (known values: SU3 6->7, 10->27, 8->0; SU5 10->1, 15->9) ===")
for lab in [(2,0),(3,0),(1,1),(0,1,0,0),(2,0,0,0),(1,0,0,0)]:
    a, d = cubic_index_A(lab); print(" A%d %s dim %d  A=%s" % (len(lab), lab, d, a))

print("\n=== exhaustive enumeration, dim <= %d ===" % DMAX)
hits = []; chiral_all = []; nrep = 0
# A_n : for N=n+1 > 46, Lambda^2 has dim N(N-1)/2 > 1000, the only non-self-dual irreps with dim<=1000
# are the fundamental and its dual (cubic index +-1 != 0). So DFS for N<=46 is exhaustive.
for n in range(2, 46):
    reps, _ = enum_irreps('A', n, DMAX)
    for lab, d in reps:
        nrep += 1
        if dual('A', n, lab) == lab: continue
        a, d2 = cubic_index_A(lab); assert d2 == d
        chiral_all.append((d, 'A%d=SU(%d)' % (n, n+1), lab, a))
        if a == 0: hits.append((d, 'A%d=SU(%d)' % (n, n+1), lab))
nA = len(chiral_all)
print("A_n (2<=n<=45): irreps scanned %d, non-self-dual %d, of these cubic-blind: %d" % (nrep, nA, len(hits)))
print("  smallest |cubic index| among non-self-dual A_n irreps:", min(abs(x[3]) for x in chiral_all))
for typ, n in [('D',5),('D',7),('D',9),('D',11),('E6',6)]:
    reps, nr = enum_irreps(typ, n, DMAX); c = 0
    for lab, d in reps:
        if dual(typ, n, lab) == lab: continue
        c += 1
        name = 'E6' if typ == 'E6' else 'D%d=Spin(%d)' % (n, 2*n)
        chiral_all.append((d, name, lab, Fraction(0)))
        hits.append((d, name, lab))   # cubic vanishes: no degree-3 invariant (degrees D_n:2,4,..,2n-2,n ; E6:2,5,6,8,9,12)
    print("%s rank %d: pos roots %d, irreps <=%d: %d, non-self-dual: %d" % (typ, n, nr, DMAX, len(reps), c))
hits.sort()
print("\nALL cubic-blind non-self-dual irreps with dim <= %d (dual pairs listed separately):" % DMAX)
for h in hits: print("  dim %4d  %-14s %s" % h)
print("\nCONTROL (a): drop cubic condition -> smallest non-self-dual irreps:")
for h in sorted(chiral_all)[:6]: print("  dim %3d %-12s %s cubic=%s" % h)

print("\n=== direct weight check of the theorem for minuscule cases (no appeal to invariant degrees) ===")
random.seed(1)
for typ, n, lab, name in [('D',5,(0,0,0,0,1),'Spin(10) 16'),('E6',6,(1,0,0,0,0,0),'E6 27'),('D',7,(0,0,0,0,0,0,1),'Spin(14) 64'),('A',4,(0,1,0,0),'SU(5) 10 [control]')]:
    C = cartan(typ, n); W = weyl_orbit(C, lab)
    Cm = __import__('sympy').Matrix(C).inv()
    res = {}
    for k in (1,3,5,7):
        vals = []
        for _ in range(3):
            h = [random.randint(-9, 9) for _ in range(n)]   # coweight coordinates
            # pairing of weight (Dynkin labels w) with H = sum h_i alpha_i^vee  is sum w_i h_i
            vals.append(sum(sum(wi*hi for wi, hi in zip(w, h))**k for w in W))
        res[k] = vals
    print(" %-20s #weights=%d  sum(w.H)^k for 3 random H: k=1 %s  k=3 %s  k=5 %s  k=7 %s" % (name, len(W), res[1], res[3], res[5], res[7]))

print("\n=== CONTROL (b) = hypothesis H4: drop irreducibility. Multisets of nonzero integers, |q|<=10,")
print("    p1 = p3 = 0, no pair (q,-q); counted up to overall sign, gcd=1 ===")
Q = 10; vals = [q for q in range(-Q, Q+1) if q]
for k in range(3, 8):
    sols = set()
    for m in itertools.combinations_with_replacement(vals, k):
        if sum(m) or sum(q**3 for q in m): continue
        s = set(m)
        if any(-q in s for q in s): continue
        g = 0
        for q in m: g = gcd(g, abs(q))
        if g != 1: continue
        sols.add(min(m, tuple(sorted(-q for q in m))))
    ex = sorted(sols, key=lambda t: max(abs(x) for x in t))[:3]
    print("  size %d: %d solutions; first: %s" % (k, len(sols), ex))
print("  diagnostic comparison only (target, not input): multiset {1^6,-4^3,2^3,-3^2,6^1}")
sm = [1]*6 + [-4]*3 + [2]*3 + [-3]*2 + [6]
print("    p1=%d p3=%d ; proper nonempty sub-multisets with p1=p3=0:" % (sum(sm), sum(q**3 for q in sm)), end=' ')
subs = set()
for r in range(1, 15):
    for c in itertools.combinations(range(15), r):
        m = tuple(sorted(sm[i] for i in c))
        if sum(m) == 0 and sum(q**3 for q in m) == 0: subs.add(m)
print(len(subs), sorted(subs)[:4])
