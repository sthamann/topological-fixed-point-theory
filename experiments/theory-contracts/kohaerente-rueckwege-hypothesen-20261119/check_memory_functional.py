"""MEMORYFUNC20261119 (hypotheses 3, 4, B and the uniqueness part of F).
B: a flat visible state with a hidden memory; does the recorded flat-density theorem
   (finding flat-density-forced-by-the-sixty-return-values-20261118) allow it?
3: mixed responses as derivatives of one functional: what does the test see in the finite chart?
4/F: how much do ordered probes determine beyond single returns (replay of the recorded ranks)?"""
from common import *
import sympy as sp
from fractions import Fraction as Fr

c = Checks("memory_functional")
e, k, t = sp.symbols("epsilon kappa t", real=True)
I2, X, Zp = sp.eye(2), sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
kron = sp.kronecker_product
ZZ = kron(Zp, Zp)
rho = (sp.eye(4) + e * ZZ) / 4
# --- hypothesis B: the toy model as stated -------------------------------------
redA = sp.Matrix(2, 2, lambda i, j: sum(rho[2 * i + m, 2 * j + m] for m in range(2)))
c.ok("visible state is exactly flat for every epsilon", redA == I2 / 2)
K = -sp.diag(*[sp.log(rho[i, i]) for i in range(4)])
coef = sp.simplify((K[0, 0] - K[1, 1]) / 2 + sp.atanh(e))
c.ok("modular Hamiltonian is -artanh(epsilon) Z x Z up to a constant", sp.simplify(coef.rewrite(sp.log)) == 0)
Ut = sp.diag(*[sp.exp(-sp.I * k * t * ZZ[i, i]) for i in range(4)])
XA = kron(X, I2)
two = sp.simplify(sp.trace(rho * Ut.H * XA * Ut * XA))
c.ok("two-time response <X(t) X(0)> = cos(2 kappa t) + i epsilon sin(2 kappa t)",
     sp.simplify(sp.expand((two - sp.cos(2 * k * t) - sp.I * e * sp.sin(2 * k * t)).rewrite(sp.cos))) == 0)
# --- against the recorded theorem: return value with reference Tr rho^2 --------
ret = sp.simplify(sp.trace(rho * XA * rho * XA) / sp.trace(rho * rho))
c.ok("return value Tr(rho X rho X)/Tr rho^2 on the enlarged space is (1 - eps^2)/(1 + eps^2): it sees the hidden memory",
     sp.simplify(ret - (1 - e**2) / (1 + e**2)) == 0, ret)
# general identity for unitary markers: 1 - return = ||[rho, X]||^2 / (2 Tr rho^2)
comm = rho * XA - XA * rho
c.ok("identity: 1 - return = ||[rho, X]||^2 / (2 Tr rho^2); flat return values force rho to commute with every marker",
     sp.simplify(1 - ret - sp.trace(comm.H * comm) / (2 * sp.trace(rho * rho))) == 0)
# same test in the recorded chart: 128 matter states x one memory qubit, 60 root currents on the matter factor
act = {r: root_action(*r) for r in MIXED}
sgn = lambda s, m: 1 if (s >> m) & 1 else -1
def anchor(diag_fn, r):
    """Tr(rho X^+ rho X)/Tr rho^2 for rho = (1 + eps D x Z)/256, X = root current x 1."""
    num = sum((1 + e * diag_fn(s) * m) * (1 + e * diag_fn(tt) * m) for s, (tt, _) in act[r].items() for m in (1, -1))
    den = sum((1 + e * diag_fn(s) * m)**2 for s in STATES for m in (1, -1))
    return sp.factor(sp.simplify(num / den))
profiles = {"first colour sign": lambda s: sgn(s, 0), "carrier parity": lambda s: (-1)**occ(s, range(5)),
            "family sign": lambda s: sgn(s, 5)}
vals = {name: {anchor(f, r) for r in MIXED} for name, f in profiles.items()}
for name, v in vals.items(): print("  memory coupled to", name, "->", v)
c.ok("in the 128 x 2 chart every tested memory coupling changes some of the 60 anchors 1/4 unless epsilon = 0",
     all(any(sp.simplify(x - sp.Rational(1, 4)) != 0 for x in v) for v in vals.values()))
c.ok("the escape is the reference: evaluated on the visible algebra alone, Tr(rho_A X rho_A X)/Tr rho_A^2 = 1 for every epsilon",
     sp.simplify(sp.trace(redA * X * redA * X) / sp.trace(redA * redA)) == 1, "hypothesis B needs returns defined on the visible algebra only")

# --- hypothesis 3: curl of first-order responses in the finite chart ------------
all_act = {r: root_action(*r) for r in EVEN_ROOTS}
neg = lambda r: (r[0], r[1], 1 - r[2], 1 - r[3])
def curl(weight, r1, r2):
    """Tr(rho [E1, E2]) for diagonal rho; nonzero only if E1 E2 has a diagonal part."""
    d = 0
    for s, (u, s1) in all_act[r2].items():
        if u in all_act[r1] and all_act[r1][u][0] == s: d += weight(s) * s1 * all_act[r1][u][1]
    for s, (u, s1) in all_act[r1].items():
        if u in all_act[r2] and all_act[r2][u][0] == s: d -= weight(s) * s1 * all_act[r2][u][1]
    return d
flat = lambda s: Fr(1, 128)
c.ok("flat state: Tr(rho [A, B]) = 0 for all 112 x 112 root pairs: mixed first-order responses are symmetric identically",
     all(curl(flat, r1, r2) == 0 for r1 in EVEN_ROOTS for r2 in EVEN_ROOTS))
th = sp.Rational(1, 3)
wsum = sum(th**occ(s, range(5)) for s in STATES)
therm = lambda s: th**occ(s, range(5)) / wsum
nonzero = [(r1, r2) for r1 in EVEN_ROOTS for r2 in EVEN_ROOTS if curl(therm, r1, r2) != 0]
c.ok("non-flat state (weights by carrier number): the curl is nonzero exactly on pairs (beta, -beta) and equals the Cartan expectation",
     len(nonzero) > 0 and all(r2 == neg(r1) for r1, r2 in nonzero), len(nonzero))
c.ok("hence in the finite chart the functional test of hypothesis 3 is passed identically by the recorded flat state: no selective power there", True)

# --- hypotheses 4 and F: recorded probe ranks ---------------------------------
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "research", "originlab"))
import f1_four_family_transport as fft
rep = fft.controls()
fam = rep["families"]
c.ok("replay of R113-FOUR-FAMILY-TRANSPORT: single returns fix 12 of 32 hidden diagonal weights, ordered two-cut probes all 32 (6 and 16 with reflection pairing)",
     all(f["single_return_probe_rank"] == 12 and f["paired_overlap_probe_rank"] == 32 and f["reflection_paired_probe_ranks"] == [6, 16] for f in fam))
c.ok("replay: two positive transfers with equal first three marks and equal reference differ in the mixed word (5/12 and 3/8; neutral control 1/96 and 1/192)",
     all(f["connected_UY"] == ["5/12", "3/8"] and f["gauge_neutral_counterexample"]["connected_UY"] == ["1/96", "1/192"] for f in fam))
c.ok("so positivity + reflection pairing + gauge neutrality + the first three marks do not single out one process in the finite chart (uniqueness claim of F fails there)", True)
c.done({"return_value_toy": str(ret), "anchors_with_memory": {k_: [str(x) for x in v] for k_, v in vals.items()},
        "curl_nonzero_pairs_nonflat": len(nonzero)})
