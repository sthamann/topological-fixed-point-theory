"""COMMONDRESS20261119 (hypothesis 1 second part, hypothesis C).
Setting: E8MIDDLEBLOCKS20261116. One positive middle operator W for all 60 mixed
currents misses the colour and weak targets by at least 7/528; two blocks hit them.
Hypothesis: one common positive intermediate operator G, read through different
E8 couplings, produces both blocks:
    W_type(t) = delta + sum over roots beta in K_type of G(t + beta),
with K_colour, K_weak unions of gauge classes of even D8 roots (unit couplings from
the root operators, nothing fitted), G >= 0 a function of the occupation class (k, l, f),
delta >= 0 one common direct term. Tested by linear programming over the declared menu."""
from common import *
import numpy as np
from fractions import Fraction as Fr
from itertools import combinations
from scipy.optimize import linprog

c = Checks("common_middle")
act = {r: root_action(*r) for r in EVEN_ROOTS}
a = (Fr(-1, 6),) * 3 + (Fr(1, 4),) * 2
h2 = {s: float(sum(a[i] * (1 if (s >> i) & 1 else -1) for i in range(5)))**2 for s in STATES}
mu = float(sum(x * x for x in a))
target = [float((sum(x * x for x in a) - 2 * a[i]**2) / 4) for i in range(5)]
label = lambda s: min(s >> 5, 7 - (s >> 5))
classes = sorted({KLF(s) for s in STATES}); CI = {k: n for n, k in enumerate(classes)}
region = lambda m: "c" if m < 3 else ("w" if m < 5 else "f")
rclass = lambda r: (region(r[0]) + region(r[1]), "same" if r[2] == r[3] else "opp")
RC = sorted({rclass(r) for r in EVEN_ROOTS})
c.ok("12 gauge classes of even roots, 24 occupation classes", len(RC) == 12 and len(classes) == 24,
     {k: sum(rclass(r) == k for r in EVEN_ROOTS) for k in RC})

# M[kappa][t, class] = number of roots beta in class kappa with t + beta defined and of that occupation class
M = {k: np.zeros((128, 24)) for k in RC}
for r in EVEN_ROOTS:
    for s, (t, _) in act[r].items(): M[rclass(r)][INDEX[s], CI[KLF(t)]] += 1
ident = np.zeros((128, 24))
for s in STATES: ident[INDEX[s], CI[KLF(s)]] = 1


def rows(matrix, typ, extra_delta=True, sels=(None, 0, 1, 2, 3)):
    """Equations ZY = 1/4 and ZUY = target for the currents of one type; unknowns (delta, G)."""
    A, b = [], []
    for r in MIXED:
        if (r[0] < 3) != (typ == "c"): continue
        for om in sels:
            pairs = [(s, t) for s, (t, _) in act[r].items() if label(s) != om]
            n = 128 if om is None else 96
            y = sum(matrix[INDEX[t]] for _, t in pairs) / n
            uy = sum(h2[s] * matrix[INDEX[t]] for s, t in pairs) / n
            A += [np.r_[len(pairs) / n if extra_delta else 0, y], np.r_[sum(h2[s] for s, _ in pairs) / n if extra_delta else 0, uy]]
            b += [0.25, target[r[0]]]
    A, b = np.array(A), np.array(b)
    kind = np.arange(len(b)) % 2                          # 0: ZY equation, 1: ZUY equation
    _, keep = np.unique(np.round(np.c_[A, b, kind], 12), axis=0, return_index=True)
    return A[keep], b[keep], kind[keep]


def gap(blocks):
    """As recorded: ZY = 1/4 exactly, minimise the largest |ZUY - target| over x >= 0."""
    A = np.vstack([blk[0] for blk in blocks]); b = np.concatenate([blk[1] for blk in blocks])
    kind = np.concatenate([blk[2] for blk in blocks]) == 1
    n = A.shape[1]
    Au, bu = A[kind], b[kind]
    Aub = np.vstack([np.c_[Au, -np.ones(len(bu))], np.c_[-Au, -np.ones(len(bu))]])
    res = linprog(np.r_[np.zeros(n), 1], A_ub=Aub, b_ub=np.r_[bu, -bu], A_eq=np.c_[A[~kind], np.zeros((~kind).sum())],
                  b_eq=b[~kind], bounds=[(0, None)] * (n + 1), method="highs")
    return res.fun


def misfit(blocks):
    """min over x >= 0 of the largest absolute residual of the stacked equations."""
    A = np.vstack([blk[0] for blk in blocks]); b = np.concatenate([blk[1] for blk in blocks])
    n = A.shape[1]
    obj = np.r_[np.zeros(n), 1]
    Aub = np.vstack([np.c_[A, -np.ones(len(b))], np.c_[-A, -np.ones(len(b))]])
    res = linprog(obj, A_ub=Aub, b_ub=np.r_[b, -b], bounds=[(0, None)] * (n + 1), method="highs")
    return res.fun if res.status == 0 else None, res

# --- controls: reproduce the recorded statements in this reconstruction ------
free_c, free_w = rows(ident, "c", False), rows(ident, "w", False)
t_two = max(misfit([free_c])[0], misfit([free_w])[0])
c.ok("control: two separate positive blocks hit all targets (as recorded)", t_two < 1e-9, t_two)
t_full = gap([rows(ident, "c", False, (None,)), rows(ident, "w", False, (None,))])
c.ok("control: one common positive block on the full trace misses by exactly the recorded 7/528", abs(t_full - 7 / 528) < 1e-9,
     f"{t_full:.12f} vs 7/528 = {7/528:.12f}")
t_one = gap([free_c, free_w])
c.ok("with all five outer selections imposed at once the common block misses by at least as much (lower bound 7/528 respected)", t_one > 7 / 528 - 1e-9,
     f"{t_one:.12f} ~ {Fr(t_one).limit_denominator(5000)}")

# --- the hypothesis: one G, two coupling sets ---------------------------------
menu = [()] + [k for n in (1, 2, 3, 4) for k in combinations(range(12), n)]
Msum = lambda K: sum((M[RC[i]] for i in K), np.zeros((128, 24)))
blocks = {K: {typ: rows(Msum(K), typ) for typ in "cw"} for K in menu}
single = {typ: {K: misfit([blocks[K][typ]])[0] for K in menu} for typ in "cw"}
ok_c = [K for K in menu if single["c"][K] < 1e-9]
ok_w = [K for K in menu if single["w"][K] < 1e-9]
print(f"  menu: {len(menu)} coupling sets (unions of at most four root classes); alone feasible: colour {len(ok_c)}, weak {len(ok_w)}")
same = [misfit([blocks[K]["c"], blocks[K]["w"]])[0] for K in menu]
c.ok("negative control: the same coupling set for both types never works", min(same) > 1e-6, f"smallest misfit {min(same):.6f}")
best, feasible, count = (1.0, None), [], 0
for Kc in ok_c:
    for Kw in ok_w:
        if Kc == Kw: continue
        count += 1
        t = misfit([blocks[Kc]["c"], blocks[Kw]["w"]])[0]
        if t < best[0]: best = (t, (Kc, Kw))
        if t < 1e-9: feasible.append((Kc, Kw))
names = lambda K: [f"{RC[i][0]}:{RC[i][1]}" for i in K]
print(f"  pairs tested: {count}; feasible: {len(feasible)}; smallest misfit {best[0]:.3e} at {names(best[1][0])} / {names(best[1][1])}")
smallest = sorted(feasible, key=lambda p: (len(p[0]) + len(p[1]), p))[:6]
for Kc, Kw in smallest: print("   feasible:", names(Kc), "/", names(Kw))
verdict = "FEASIBLE_IN_MENU" if feasible else "EXCLUDED_IN_MENU"
c.ok("verdict recorded for: one common positive G read through different root couplings", True,
     f"{verdict}: {len(feasible)} of {count} ordered pairs")

# how tight is a feasible solution? count free directions for the simplest pair, and check delta = 0 variant
extra = {}
if feasible:
    Kc, Kw = smallest[0]
    A = np.vstack([blocks[Kc]["c"][0], blocks[Kw]["w"][0]])
    extra["simplest_pair"] = [names(Kc), names(Kw)]
    extra["free_directions_of_simplest_pair"] = int(A.shape[1] - np.linalg.matrix_rank(A))
    nodelta = [p for p in feasible if misfit([(blocks[p[0]]["c"][0][:, 1:], blocks[p[0]]["c"][1]),
                                              (blocks[p[1]]["w"][0][:, 1:], blocks[p[1]]["w"][1])])[0] < 1e-9]
    extra["feasible_without_direct_term"] = len(nodelta)
    gauge_only = [p for p in feasible if all(RC[i][0] in ("cc", "ww", "cw") for K in p for i in K)]
    extra["feasible_with_carrier_roots_only"] = len(gauge_only)
    print("  ", extra)
    c.ok("a feasible pair leaves G undetermined (free directions remain): existence, not selection",
         extra["free_directions_of_simplest_pair"] > 0, extra["free_directions_of_simplest_pair"])
c.done({"verdict": verdict, "pairs_tested": count, "feasible_pairs": len(feasible), "smallest_misfit": best[0],
        "common_block_gap": t_one, **extra})
