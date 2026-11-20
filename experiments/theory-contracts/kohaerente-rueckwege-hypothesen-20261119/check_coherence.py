"""NEUTRALCOH20261119 (hypothesis 1, first part of hypothesis C).
Object: the neutral coherence of finding d8-neutral-coherence-needs-four-root-word-20261118
in the D8 half-spinor CAR chart (five carrier modes, three family modes).
Questions: what is this coherence, which words see it, can it carry the
difference between the colour and the weak response, and is its phase fixed?"""
from common import *
import numpy as np
from fractions import Fraction as Fr
from collections import deque

c = Checks("coherence")
N = len(STATES)
act = {r: root_action(*r) for r in EVEN_ROOTS}
c.ok("chart: 128 even states, 112 even D8 roots, 60 mixed (carrier, family) roots, each mixed root on 32 states",
     N == 128 and len(EVEN_ROOTS) == 112 and len(MIXED) == 60 and all(len(act[r]) == 32 for r in MIXED))

# --- 1. what the coherence is -------------------------------------------------
Y = lambda s: Fr(-occ(s, (0, 1, 2)), 3) + Fr(occ(s, (3, 4)), 2)
# irreducible sectors of number-preserving SU(3)c x SU(2)w x SU(3)F are the occupation classes (k, l, f)
rep = lambda s: (min(KLF(s)[0], 3) if KLF(s)[0] in (0, 3) else KLF(s)[0], KLF(s)[1] % 2, Y(s))
sm_type = lambda s: ("1" if KLF(s)[0] in (0, 3) else ("3" if KLF(s)[0] == 1 else "3b"),
                     "1" if KLF(s)[1] != 1 else "2", Y(s))
fam_type = lambda s: "1" if KLF(s)[2] in (0, 3) else ("3" if KLF(s)[2] == 1 else "3b")
sectors = {}
for s in STATES: sectors.setdefault(KLF(s), []).append(s)
types = {}
for klf, members in sectors.items():
    types.setdefault((sm_type(members[0]), fam_type(members[0])), []).append(klf)
double = {t: v for t, v in types.items() if len(v) > 1}
c.ok("24 occupation sectors; exactly one type occurs twice: the total singlet (1,1)_0 x family singlet, on the states 0 and 255",
     len(sectors) == 24 and list(double.values()) == [[(0, 0, 0), (3, 2, 3)]]
     and sectors[(0, 0, 0)] == [0] and sectors[(3, 2, 3)] == [255], double)


def commutant_dim(offdiag_blocks, cartan):
    """dim of operators on the 128 commuting with number-preserving generators (numeric rank of an integer system)."""
    gens = [(i, j, 1, 0) for blk in offdiag_blocks for i in blk for j in blk if i != j]
    mats = []
    for i, j, a, b in gens:
        M = {}
        for s in STATES:
            f = letter(s, j, b)
            if f is None: continue
            g = letter(f[0], i, a)
            if g is not None: M[s] = (g[0], f[1] * g[1])
        mats.append(M)
    w = {s: tuple(h(s) for h in cartan) for s in STATES}
    unknown = {(s, t): n for n, (s, t) in enumerate((s, t) for s in STATES for t in STATES if w[s] == w[t])}
    rows = {}
    for g, M in enumerate(mats):
        for (s, t), n in unknown.items():              # contribution of X[s,t] to [E, X]
            if s in M:                                 # (E X)[M(s), t] += sign * X[s,t]
                u, sg = M[s]; rows.setdefault((g, u, t), {})[n] = rows.setdefault((g, u, t), {}).get(n, 0) + sg
        inv = {}
        for s, (u, sg) in M.items(): inv[u] = (s, sg)
        for (s, t), n in unknown.items():              # (X E)[s, v] with E[t, v]: v -> t
            if t in inv:
                v, sg = inv[t]; rows.setdefault((g, s, v), {})[n] = rows.setdefault((g, s, v), {}).get(n, 0) - sg
    Mx = np.zeros((len(rows), len(unknown)))
    for r, row in enumerate(rows.values()):
        for n, val in row.items(): Mx[r, n] = val
    return len(unknown) - np.linalg.matrix_rank(Mx)


d = lambda a, b: (lambda s: ((s >> a) & 1) - ((s >> b) & 1))
num = lambda modes: (lambda s: occ(s, modes))
sm_cartan = [d(0, 1), d(1, 2), d(3, 4), Y]
dim_sm = commutant_dim([(0, 1, 2), (3, 4)], sm_cartan)
dim_smf = commutant_dim([(0, 1, 2), (3, 4), (5, 6, 7)], sm_cartan + [d(5, 6), d(6, 7)])
dim_all = commutant_dim([(0, 1, 2), (3, 4), (5, 6, 7)], sm_cartan + [d(5, 6), d(6, 7), num((0, 1, 2, 3, 4))])
c.ok("commutant of SM x family SU(3) has dimension 26 = 22 + 2^2: the only neutral two-block is span{|0>, |255>}",
     dim_smf == 26, {"SM only": dim_sm, "SM x SU(3)_F": dim_smf, "plus carrier number": dim_all})
c.ok("adding the carrier number (the extra U(1) of D5, B-L type) splits that block: commutant 24, no neutral coherence left",
     dim_all == 24)
c.ok("the coherence joins the carrier-empty SU(5) singlet (nu^c) of the omitted family with its conjugate: carrier number changes by 5, family number by 3",
     KLF(0) == (0, 0, 0) and KLF(255) == (3, 2, 3) and Y(0) == Y(255) == 0)

# --- 2. which words reach it --------------------------------------------------
def distance(roots, start=0, goal=255):
    seen, queue = {start: 0}, deque([start])
    while queue:
        s = queue.popleft()
        for r in roots:
            if s in act[r] and act[r][s][0] not in seen:
                seen[act[r][s][0]] = seen[s] + 1; queue.append(act[r][s][0])
    return seen.get(goal)
c.ok("shortest root word from |0> to |255>: four roots with all even roots, five with the 60 mixed roots alone",
     (distance(EVEN_ROOTS), distance(MIXED)) == (4, 5), (distance(EVEN_ROOTS), distance(MIXED)))
word = [(0, 5, 1, 1), (1, 6, 1, 1), (3, 7, 1, 1), (2, 4, 1, 1)]
s, sign = 0, 1
for r in word:
    s, sg = act[r][s]; sign *= sg
c.ok("the word (0,5),(1,6),(3,7),(2,4) maps |0> to +-|255>: three mixed roots and one carrier root", s == 255, sign)

# --- 3. the colour and weak marks do not see it -------------------------------
a = (Fr(-1, 6),) * 3 + (Fr(1, 4),) * 2
h = lambda s: sum(a[i] * (1 if (s >> i) & 1 else -1) for i in range(5))
mu = sum(x * x for x in a)
target = lambda i: (mu - 2 * a[i]**2) / 4
label = lambda s: min(s >> 5, 7 - (s >> 5))
selections = [None, 0, 1, 2, 3]


def marks(W):
    """ZY and ZUY of every mixed current and outer selection for a Hermitian middle operator W (128 x 128)."""
    out = []
    for r in MIXED:
        for om in selections:
            pairs = [(s, t) for s, (t, _) in act[r].items() if label(s) != om]
            n = 128 if om is None else 96
            out.append((sum(W[INDEX[t], INDEX[t]] for _, t in pairs) / n,
                        sum(float(h(s))**2 * W[INDEX[t], INDEX[t]] for s, t in pairs) / n))
    return np.array(out)
rng = np.random.default_rng(20261119)
W0 = np.diag(rng.uniform(0.5, 1.5, N)).astype(complex)
W1 = W0.copy(); W1[INDEX[0], INDEX[255]] = 0.3 + 0.2j; W1[INDEX[255], INDEX[0]] = 0.3 - 0.2j
# explicit operator form of one mark, to confirm that only the diagonal of W enters
def mark_operator(W, r):
    X = np.zeros((N, N))
    for s_, (t, sg) in act[r].items(): X[INDEX[t], INDEX[s_]] = sg
    U = np.diag([float(h(s_))**2 for s_ in STATES])
    return np.trace(U @ X.T @ W @ X) / N
c.ok("operator form Tr(U X* W X) agrees with the diagonal formula on all 60 currents",
     all(abs(mark_operator(W1, r) - marks(W1)[5 * n, 1]) < 1e-12 for n, r in enumerate(MIXED)))
c.ok("all 600 marks ZY, ZUY (60 currents, five outer selections) are unchanged by the neutral coherence",
     np.abs(marks(W0) - marks(W1)).max() < 1e-13 and np.linalg.eigvalsh(W1).min() > 0)

# --- 4. is its phase fixed? ---------------------------------------------------
def car(mode):
    M = np.zeros((256, 256))
    for s_ in range(256):
        for create in (0, 1):
            f = letter(s_, mode, create)
            if f: M[f[0], s_] += f[1]
    return M
PH = np.eye(256)
for m in range(8): PH = PH @ car(m)                 # particle-hole map prod (a_m + a_m^+)
x_, y_ = PH[255, 0], PH[0, 255]
c.ok("occupation complement maps |0> <-> |255> with signs", abs(x_) == 1 and abs(y_) == 1, (x_, y_))
# linear invariance PH W PH^-1 = W forces W[0,255] = x*y * conj(W[0,255])
c.ok("complement pairing leaves the coherence c real (x*y = +1) or imaginary (x*y = -1): a sign, not a phase, stays free",
     x_ * y_ in (1, -1), "x*y = %+d -> c is %s" % (x_ * y_, "real" if x_ * y_ == 1 else "imaginary"))
blockmin = lambda cc: np.linalg.eigvalsh(np.array([[1.0, cc], [np.conj(cc), 0.25]])).min()
c.ok("positivity bounds only its size: with diagonal entries 1 and 1/4 the block is positive exactly for |c| <= 1/2",
     blockmin(0.5) > -1e-12 and blockmin(0.5j) > -1e-12 and blockmin(0.51) < 0)

# --- 5. a neutral 2 x 2 block as a Pauli triple (link to hypothesis 6) --------
c.ok("on span{|0>, |255>} sigma_x, sigma_y are the two quadratures of the coherence and sigma_z the population difference: SM and family-SU(3) neutral; carrier number removes exactly sigma_x and sigma_y",
     dim_smf - dim_all == 2)
c.done({"commutant_dims": {"SM": int(dim_sm), "SM_x_SU3F": int(dim_smf), "plus_carrier_number": int(dim_all)},
        "word_lengths": {"all_even_roots": 4, "mixed_only": 5}, "phase": "sign free, x*y=%d" % (x_ * y_)})
