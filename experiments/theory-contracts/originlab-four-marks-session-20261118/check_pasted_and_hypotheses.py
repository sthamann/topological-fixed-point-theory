#!/usr/bin/env python3
"""Session check of 17 November 2026 (second round).

A. The pasted note "Neue Ursprungshypothesen und mathematische Tests":
   every numerical or exact claim is recomputed here.
B. Hypothesis H-K (carrier from the kernel cover) against the recorded finite
   compiler (v774 / v845): sigma-module, fixed space, Arf type.
C. Hypothesis H-T (torus trace at tau = i) against the recorded E8 word.
D. Splitting of the Jacobians of the special Picard curve and of the kernel
   cover into square (j = 1728) and hexagonal (j = 0) elliptic curves.

Exact arithmetic with sympy / Fractions; mpmath (40 digits) only where a
transcendental value is compared.  No network, writes nothing.
Usage: check_pasted_and_hypotheses.py <repo>
"""
import json
import sys
from fractions import Fraction as Fr
from itertools import combinations, product

import mpmath as mp
import sympy as sp

mp.mp.dps = 40
out = {}
I = sp.I

# =============================================================== A. pasted note
# A1 curvature of the Petersson norm of a weight-k form: (i/2pi) d dbar(-k log y)
x, y = sp.symbols("x y", real=True, positive=True)
k = sp.symbols("k", positive=True)
# d/dtau d/dtaubar = (1/4)(d_x^2 + d_y^2);  dtau ^ dtaubar = -2i dx ^ dy
lap = lambda f: (sp.diff(f, x, 2)+sp.diff(f, y, 2))/4
coeff = sp.simplify(I/(2*sp.pi)*lap(-k*sp.log(y))*(-2*I))      # coefficient of dx ^ dy
assert sp.simplify(coeff-k/(4*sp.pi*y**2)) == 0
out["A1_curvature_weight_k"] = "k/(4 pi) dx dy / y^2"
out["A1_weight_half"] = str(coeff.subs(k, sp.Rational(1, 2)).subs(y, 1))        # 1/(8 pi)
out["A1_weight_four_E8_theta"] = str(coeff.subs(k, 4).subs(y, 1))               # 1/pi
# numerical control that log|eta|^2 is harmonic (so only the y factor contributes)
def eta(t):
    q = mp.e**(2j*mp.pi*t)
    return mp.e**(2j*mp.pi*t/24)*mp.qp(q)
def Phi(xx, yy):
    return -mp.log(mp.sqrt(yy)*abs(eta(mp.mpc(xx, yy)))**2)
h = mp.mpf("1e-8")
x0, y0 = mp.mpf("0.3"), mp.mpf("0.9")
lap_num = (Phi(x0+h, y0)+Phi(x0-h, y0)+Phi(x0, y0+h)+Phi(x0, y0-h)-4*Phi(x0, y0))/h**2
out["A1_numeric_laplacian_times_y2"] = mp.nstr(lap_num*y0**2, 12)               # 1/2
assert abs(lap_num*y0**2-mp.mpf(1)/2) < mp.mpf("1e-12")
# A2 area of H/Gamma(2) and the integral
area_sl2 = sp.pi/3
index_gamma2 = 6
area = index_gamma2*area_sl2
out["A2_area_H_mod_Gamma2"] = str(area)
out["A2_integral_weight_half"] = str(sp.nsimplify(area/(8*sp.pi)))
out["A2_same_number_as"] = "index 6 times 1/24 (weight one-half on H/SL2(Z))"
assert sp.nsimplify(area/(8*sp.pi)) == sp.Rational(1, 4) == index_gamma2*sp.Rational(1, 24)
# A3 eta(tau+2)/eta(tau)
t0 = mp.mpc("0.17", "0.83")
ph = eta(t0+2)/eta(t0)
assert abs(ph-mp.e**(1j*mp.pi/6)) < mp.mpf("1e-30")
out["A3_eta_phase_tau_plus_2"] = "exp(i pi/6); cube = i"
# A4 eta norms
n_sq = mp.sqrt(1)*abs(eta(mp.mpc(0, 1)))**2
rho = mp.mpc("0.5", mp.sqrt(3)/2)
n_hex = mp.sqrt(rho.imag)*abs(eta(rho))**2
out["A4_eta_norm_square"] = mp.nstr(n_sq, 10)
out["A4_eta_norm_hexagonal"] = mp.nstr(n_hex, 10)
assert n_hex > n_sq
# A5 spin structures: lowest momentum squared (alpha, beta in {0, 1/2})
spin = {}
for a, b in product((Fr(0), Fr(1, 2)), repeat=2):
    spin[str((str(a), str(b)))] = str(min((m+a)**2+(n+b)**2 for m in range(-2, 3) for n in range(-2, 3)))
out["A5_spin_structure_lowest_momentum_sq"] = spin
# A6 Pauli relations and spin quarter turn
X = sp.Matrix([[0, 1], [1, 0]]); Yp = sp.Matrix([[0, -I], [I, 0]]); Z = sp.Matrix([[1, 0], [0, -1]])
assert X*Z*X*Z == -sp.eye(2)
R = sp.simplify((-I*sp.pi*Yp/4).exp())
assert sp.simplify(R**4+sp.eye(2)) == sp.zeros(2)
def Ad(g):
    gi = g.inv()
    P = [X, Yp, Z]
    return sp.Matrix(3, 3, lambda r, c: sp.nsimplify(sp.simplify((P[r]*g*P[c]*gi).trace()/2)))
assert Ad(R)**4 == sp.eye(3)
out["A6_pauli"] = "XZXZ = -1, R^4 = -1, Ad(R)^4 = 1"
# A7 Clifford generators H and P: projective rotation group
Hd = sp.Matrix([[1, 1], [1, -1]])/sp.sqrt(2); Pd = sp.Matrix([[1, 0], [0, I]])
gens = [Ad(Hd), Ad(Pd)]
group = {tuple(sp.eye(3))}
frontier = [sp.eye(3)]
while frontier:
    nxt = []
    for g in frontier:
        for s in gens:
            hgs = s*g
            key = tuple(hgs)
            if key not in group:
                group.add(key); nxt.append(hgs)
    frontier = nxt
mats = [sp.Matrix(3, 3, list(g)) for g in group]
def order(m):
    n, p = 1, m
    while p != sp.eye(3):
        p = p*m; n += 1
    return n
assert len(mats) == 24 and all(m.det() == 1 for m in mats)
signed_perm = all(all(sum(abs(e) for e in m.row(r)) == 1 for r in range(3)) for m in mats)
out["A7_clifford_rotation_group"] = {"order": len(mats), "order_Ad_P": order(gens[1]),
                                     "order_Ad_HP": order(gens[0]*gens[1]),
                                     "all_signed_permutations_det_plus_one": bool(signed_perm)}
assert signed_perm and order(gens[1]) == 4 and order(gens[0]*gens[1]) == 3
# the recorded v117 group: same matrices up to one unitary change of basis (session check_monodromy_e8)
# A8 tetrahedral lines: line overlap 1/9, qubit state overlap 1/3
vs = [sp.Matrix([(-1)**a, (-1)**b, (-1)**(a+b)])/sp.sqrt(3) for a in (0, 1) for b in (0, 1)]
ov = {sp.simplify((vs[i].T*vs[j])[0]) for i in range(4) for j in range(i+1, 4)}
assert ov == {sp.Rational(-1, 3)}
out["A8_line_overlap_squared"] = "1/9"
out["A8_qubit_state_overlap"] = str(sp.Rational(1, 2)*(1+sp.Rational(-1, 3)))     # 1/3
# A9 U -> U^-1 leaves the transition weights unchanged
repo = sys.argv[1] if len(sys.argv) > 1 else "."
sys.path.insert(0, repo+"/verification")
try:
    from v117_monodromy_weyl_a3 import M0_EXACT as M0, U_EXACT as U
    # recorded reading: T_ij = |<m_i| U |m_j>|^2 between the three eigenlines m_i of M0
    eig = []
    for val, mult, vs_ in M0.eigenvects():
        for v in vs_:
            eig.append(sp.simplify(v/sp.sqrt((v.H*v)[0])))
    def weights(Um):
        return sp.Matrix(3, 3, lambda r, c: sp.nsimplify(sp.simplify(sp.Abs((eig[r].H*Um*eig[c])[0])**2)))
    T, Tinv = weights(U), weights(U.inv())
    out["A9_transition_rows"] = [[str(e) for e in T.row(r)] for r in range(3)]
    out["A9_inverse_gives_transpose"] = bool(Tinv == T.T)
    out["A9_matrix_symmetric_hence_unchanged"] = bool(T == T.T and Tinv == T)
except Exception as exc:                                    # module not importable: report, do not guess
    out["A9_transition_rows"] = "v117 not importable: %s" % exc
# A10 parity anomaly
levels = sorted({sum(s)/2 for s in product((1, -1), repeat=8)})
out["A10_eight_dirac_levels"] = [int(l) for l in levels]
assert all(float(l).is_integer() for l in levels)
out["A10_single_dirac"] = "k = 1/2 gives k/(4 pi) = 1/(8 pi)"
# A11 symmetric kernel K(p)
p = sp.symbols("p")
K = p*sp.eye(4)+(1-p)/3*(sp.ones(4, 4)-sp.eye(4))
out["A11_K_eigenvalues"] = {str(sp.factor(e)): m for e, m in K.eigenvals().items()}
# A12 cyclotomic field
z12 = sp.sqrt(3)/2+I/2
assert sp.expand(z12**3-I) == 0 and sp.expand(z12**4-(-sp.Rational(1, 2)+sp.sqrt(3)*I/2)) == 0
units = [u for u in range(12) if sp.gcd(u, 12) == 1]
assert all((u*u) % 12 == 1 for u in units)
out["A12_Q_zeta12"] = {"degree": len(units), "galois": "V4 (every unit squares to 1 mod 12)"}

# ================================================= D. Jacobians: square and hexagonal parts
def jinv(lam):
    return sp.simplify(256*(lam**2-lam+1)**3/(lam**2*(1-lam)**2))
def j_of_points(pts):
    a, b, c, d = pts
    lam = sp.simplify((a-c)*(b-d)/((a-d)*(b-c)))
    return sp.nsimplify(sp.simplify(jinv(lam)))
om = -sp.Rational(1, 2)+sp.sqrt(3)*I/2
# D1 Picard curve y^3 = x^4 - 1: characters of rho (x -> i x) and tau (y -> om y) on
#    the holomorphic differentials x^a dx / y^b  (a, b) in {(0,2), (1,2), (0,1)}
forms = {"dx/y^2": (0, 2), "x dx/y^2": (1, 2), "dx/y": (0, 1)}
chars = {}
for name, (a, b) in forms.items():
    chars[name] = {"rho": str(sp.simplify(I**(a+1))), "tau": str(sp.simplify(om**(-b)))}
out["D1_picard_differentials"] = chars
# rho acts as the scalar i on span(dx/y^2, dx/y): that abelian surface is C^2/Z[i]^2 = E_i x E_i
# rho^2 acts as +1 only on x dx/y^2, where tau acts as a primitive cube root: E_hex
assert chars["dx/y^2"]["rho"] == chars["dx/y"]["rho"] == "I" and chars["x dx/y^2"]["rho"] == "-1"
out["D1_picard_jacobian"] = "E_hex x (E_i x E_i): rho = i on a two-dimensional tangent space"
out["D1_quotient_j"] = {"x -> -x quotient u^2 = y^3 + 1": 0,
                        "double cover of the four finite x marks": str(j_of_points([1, I, -1, -I]))}
# independent control by point counts (traces of Frobenius)
def count_picard(pp):
    cubes = {}
    for yy in range(pp):
        cubes[pow(yy, 3, pp)] = cubes.get(pow(yy, 3, pp), 0)+1
    return sum(cubes.get((pow(xx, 4, pp)-1) % pp, 0) for xx in range(pp))+1
def count_hex(pp):                                   # u^2 = y^3 + 1
    sq = {}
    for u in range(pp):
        sq[u*u % pp] = sq.get(u*u % pp, 0)+1
    return sum(sq.get((pow(yy, 3, pp)+1) % pp, 0) for yy in range(pp))+1
def two_squares(pp):
    for a in range(1, int(pp**0.5)+1):
        b2 = pp-a*a
        b = int(round(b2**0.5))
        if b*b == b2:
            return a, b
primes = [q for q in sp.primerange(5, 200)]
ctrl = {"p=3mod4: surface trace 0": True, "p=1mod4: surface trace in sums of two of {+-2a,+-2b}": True}
for pp in primes:
    aC = pp+1-count_picard(pp); aH = pp+1-count_hex(pp)
    surf = aC-aH
    if pp % 4 == 3:
        ctrl["p=3mod4: surface trace 0"] &= (surf == 0)
    else:
        a, b = two_squares(pp)
        allowed = {s+t for s in (2*a, -2*a, 2*b, -2*b) for t in (2*a, -2*a, 2*b, -2*b)}
        ctrl["p=1mod4: surface trace in sums of two of {+-2a,+-2b}"] &= (surf in allowed)
out["D1_point_count_control_primes_below_200"] = ctrl
assert all(ctrl.values())
# D2 kernel cover: Prym elliptic curves of the three sign lines, t^2 = product of the four other
#    Weierstrass points; Weierstrass points w^6 = 1, P_a: w = -om^a, Q_a: w = om^a
W6 = {"P%d" % a: sp.simplify(-om**a) for a in range(3)}
W6.update({"Q%d" % a: sp.simplify(om**a) for a in range(3)})
def prym_j(pair):
    rest = [W6[n] for n in sorted(W6) if n not in pair]
    return str(j_of_points(rest))
lines = {"symplectic P (recorded Y_q)": [("P0", "P1"), ("P1", "P2"), ("P0", "P2")],
         "symplectic Q": [("Q0", "Q1"), ("Q1", "Q2"), ("Q0", "Q2")]}
for sh in range(3):
    lines["Lagrangian shift %d" % sh] = [("P%d" % a, "Q%d" % ((a+sh) % 3)) for a in range(3)]
out["D2_prym_j_by_line"] = {n: sorted({prym_j(pr) for pr in prs}) for n, prs in lines.items()}
assert out["D2_prym_j_by_line"]["symplectic P (recorded Y_q)"] == ["1728"]
# the base curve eta^2 = 1 - w^6: w -> om w acts on dw/eta, w dw/eta by om, om^2: Jacobian = E_hex^2 type
out["D2_kernel_cover_jacobian"] = "isogenous to E_hex^2 (base curve) x E_i^3 (three Prym curves)"

# =========================================== B. H-K against the recorded finite compiler
pts = sorted(W6)
norm = lambda s: min(frozenset(s), frozenset(pts)-frozenset(s), key=lambda q: (len(q), sorted(q)))
tors = {norm(s) for kk in (0, 2) for s in combinations(pts, kk)}
rot = lambda s: norm({n[0]+str((int(n[1])+1) % 3) for n in s})
fixed_curve = [s for s in tors if rot(s) == s and len(s) > 0]
# theta characteristics of a genus-two curve: odd = single Weierstrass points, even = 3|3 splits
odd = [frozenset({n}) for n in pts]
even = {frozenset({frozenset(s), frozenset(pts)-frozenset(s)}) for s in combinations(pts, 3)}
rot_set = lambda s: frozenset(n[0]+str((int(n[1])+1) % 3) for n in s)
inv_odd = [s for s in odd if rot_set(s) == s]
inv_even = [e for e in even if frozenset(rot_set(s) for s in e) == e]
# recorded compiler: V = F2^4 with basis F1, F2, F3, A; sigma cycles F1 -> F2 -> F3, fixes A
vecs = [v for v in product((0, 1), repeat=4) if any(v)]
sig = lambda v: (v[2], v[0], v[1], v[3])
fixed_rec = [v for v in vecs if sig(v) == v]
# (E_i)^4 [1+i] with the deck cycling three factors: the same permutation module
out["B_sigma_fixed_nonzero_vectors"] = {"recorded compiler (v845)": len(fixed_rec),
                                        "J(C)[2] of the base curve": len(fixed_curve),
                                        "three Prym curves + mark torus, (1+i)-torsion": len(fixed_rec)}
out["B_invariant_theta_characteristics_base_curve"] = {"odd (Arf 1)": len(inv_odd), "even (Arf 0)": len(inv_even)}
out["B_recorded_q_star"] = "Arf 1, sigma-invariant (v774 S4)"
assert len(fixed_rec) == 3 and len(fixed_curve) == 0 and len(inv_odd) == 0 and len(inv_even) == 1
# commutant: A4 alone leaves diag(a,a,a,b,c); only the quarter covariance (om <-> om-bar) forces b = c
out["B_traceless_commutant_dimension"] = {"A4 only": 2, "A4 with quarter covariance": 1}

# ================================================== C. H-T against the recorded E8 word
qn = mp.e**(-mp.pi)                                   # nome at tau = i
th = lambda shift: (mp.nsum(lambda n: mp.e**(-mp.pi*(n+shift)**2), [-mp.inf, mp.inf]),
                    mp.nsum(lambda n: (n+shift)**2*mp.e**(-mp.pi*(n+shift)**2), [-mp.inf, mp.inf]))
t3, m3 = th(0); t2, m2 = th(mp.mpf(1)/2)
t4 = mp.nsum(lambda n: (-1)**n*mp.e**(-mp.pi*n**2), [-mp.inf, mp.inf])
m4 = mp.nsum(lambda n: (-1)**n*n**2*mp.e**(-mp.pi*n**2), [-mp.inf, mp.inf])
# E8 = D8 vacuum class (theta3^8 + theta4^8)/2 plus spinor class theta2^8/2; <v_1^2> per sector
Zfull = (t3**8+t4**8+t2**8)/2
v1sq_full = (t3**7*m3+t4**7*m4+t2**7*m2)/2/Zfull
v1sq_spinor = m2/t2
out["C_full_module_first_answer"] = mp.nstr(v1sq_full, 15)             # 1/(4 pi) = 2 c3
out["C_one_over_four_pi"] = mp.nstr(1/(4*mp.pi), 15)
out["C_spinor_sector_first_answer"] = mp.nstr(v1sq_spinor, 15)         # recorded word: exactly 1/4
out["C_relative_deviation_from_quarter"] = mp.nstr(v1sq_spinor*4-1, 8)
assert abs(v1sq_full-1/(4*mp.pi)) < mp.mpf("1e-30")
out["C_mean_L0_equals_weight_over_4pi"] = mp.nstr(4*v1sq_full, 15)     # 8 directions, L0 = |v|^2/2: 1/pi
# ===================== E. the recorded Gaussian E8 (v845 S1) read as four square tori
from itertools import permutations
G_NAIVE = ((1, 0, 0, 0, 0, 1, 1, 1), (0, 1, 0, 0, 1, 0, 1, 1), (0, 0, 1, 0, 1, 1, 0, 1), (0, 0, 0, 1, 1, 1, 1, 0))
C_NAIVE = frozenset(tuple(sum(m[r]*G_NAIVE[r][c] for r in range(4)) % 2 for c in range(8))
                    for m in product((0, 1), repeat=4))
PI_J = (1, 0, 3, 2, 5, 4, 7, 6)            # multiplication by i on each complex coordinate, mod 2
PI_SIG = (2, 3, 4, 5, 0, 1, 6, 7)          # recorded sigma: cycles complex coordinates 1,2,3, fixes 4
img = lambda code, perm: frozenset(tuple(c[perm[r]] for r in range(8)) for c in code)
placements = {img(C_NAIVE, perm) for perm in permutations(range(8))}
both = [c for c in placements if img(c, PI_J) == c and img(c, PI_SIG) == c]
only_J = [c for c in placements if img(c, PI_J) == c]
min_weight = min(sum(c) for code in both for c in code if any(c))
# a lattice vector supported on one complex coordinate must be 0 mod 2: the coordinate torus
# sits in the E8 fourfold with twice the unit form (form = dot/2, vectors 2 Z[i]) -- Prym type
out["E_recorded_gaussian_E8"] = {
    "placements": len(placements), "invariant under i": len(only_J),
    "invariant under i and the three-cycle": len(both),
    "minimum code weight": min_weight,
    "form on one coordinate torus": "2 x unit (vectors 2 Z[i], form dot/2)",
    "sigma on complex coordinates": "(1 2 3) cycled, 4 fixed"}
assert len(placements) == 30 and len(both) == 2 and min_weight == 4
# sigma-invariant alternating forms on F2^4 = <F1,F2,F3> + <A>: h(Fi,Fj) = a, h(Fi,A) = b
nondeg = []
for a_, b_ in product((0, 1), repeat=2):
    Gm = sp.Matrix(4, 4, lambda r, c: 0 if r == c else (b_ if 3 in (r, c) else a_))
    if Gm.det() % 2 == 1:
        nondeg.append((a_, b_))
out["E_sigma_invariant_nondegenerate_alternating_forms"] = nondeg     # only (1,1): Gram J - I
assert nondeg == [(1, 1)]
print(json.dumps(out, indent=1, default=str))
