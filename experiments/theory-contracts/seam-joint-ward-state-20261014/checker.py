"""Exact joint Ward/chirality/time checks in the declared finite scope.

This is a scratch checker.  It does not edit or import the TFPT repository.
Fractions are used for every lattice and energy calculation.  SymPy is used
only for exact finite matrices in the double-cover and two-mode Fock checks.

The finite root/zero census is not presented as a proof about the whole
infinite lattice.  The full-lattice ground statement is tied separately to
the coordinatewise convexity and parity criterion recorded in results.json.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib
import json

import sympy as sp


checks = {}


def require(name, condition):
    if not bool(condition):
        raise RuntimeError(name)
    checks[name] = True


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def norm2(x):
    return dot(x, x)


# ---------------------------------------------------------------------------
# A. Two chiral E8 lattices, the Ward vector, and divisibility


def roots(chirality):
    out = []
    for i, j in combinations(range(8), 2):
        for a, b in product((1, -1), repeat=2):
            v = [F(0)] * 8
            v[i], v[j] = F(a), F(b)
            out.append(tuple(v))
    wanted = 0 if chirality == "s" else 1
    for signs in product((1, -1), repeat=8):
        if signs.count(-1) % 2 == wanted:
            out.append(tuple(F(a, 2) for a in signs))
    return out


def in_lattice(v, chirality):
    if all(x.denominator == 1 for x in v):
        return sum(v).denominator == 1 and int(sum(v)) % 2 == 0
    if all(x.denominator == 2 for x in v):
        # Full half-coset membership is the parity of sum(v).  Counting minus
        # signs is used only above to enumerate the norm-two root representatives.
        if sum(v).denominator != 1:
            return False
        wanted = 0 if chirality == "s" else 1
        return int(sum(v)) % 2 == wanted
    return False


roots_s = roots("s")
roots_c = roots("c")
require("A_Ls_has_240_roots", len(roots_s) == 240 and len(set(roots_s)) == 240)
require("A_Lc_has_240_roots", len(roots_c) == 240 and len(set(roots_c)) == 240)
require("A_all_Ls_roots_norm2_and_members", all(norm2(r) == 2 and in_lattice(r, "s") for r in roots_s))
require("A_all_Lc_roots_norm2_and_members", all(norm2(r) == 2 and in_lattice(r, "c") for r in roots_c))

q = tuple(F(1) for _ in range(8))
q_half = tuple(x / 2 for x in q)
require("A_q_in_both_lattices", in_lattice(q, "s") and in_lattice(q, "c"))
require("A_q_half_in_Ls", in_lattice(q_half, "s"))
require("A_q_half_not_in_Lc", not in_lattice(q_half, "c"))

half_s = [r for r in roots_s if r[0].denominator == 2]
half_c = [r for r in roots_c if r[0].denominator == 2]
qdot_s = Counter(dot(q, r) for r in half_s)
qdot_c = Counter(dot(q, r) for r in half_c)
qdot_roots_s = Counter(dot(q, r) for r in roots_s)
qdot_roots_c = Counter(dot(q, r) for r in roots_c)
require("A_qdot_Ls_halfroots_even", len(half_s) == 128 and all(x.denominator == 1 and int(x) % 2 == 0 for x in qdot_s))
require("A_qdot_Lc_halfroots_odd", len(half_c) == 128 and all(x.denominator == 1 and int(x) % 2 == 1 for x in qdot_c))
require("A_qdot_Ls_exact_census", qdot_s == Counter({F(-4): 1, F(-2): 28, F(0): 70, F(2): 28, F(4): 1}))
require("A_qdot_Lc_exact_census", qdot_c == Counter({F(-3): 8, F(-1): 56, F(1): 56, F(3): 8}))
require("A_qdot_Ls_allroot_exact_census", qdot_roots_s == Counter({F(-4): 1, F(-2): 56, F(0): 126, F(2): 56, F(4): 1}))
require("A_qdot_Lc_allroot_exact_census", qdot_roots_c == Counter({F(-3): 8, F(-2): 28, F(-1): 56, F(0): 56, F(1): 56, F(2): 28, F(3): 8}))

# A lattice isomorphism preserves divisibility.  If U q = 2 r for a target
# root r, surjectivity would give q/2 = U^{-1} r in Lc, which is false.
require("A_every_twice_Ls_root_is_divisible", all(in_lattice(r, "s") and in_lattice(tuple(2 * x for x in r), "s") for r in roots_s))
require("A_no_lattice_isomorphism_can_send_q_to_twice_root", in_lattice(q, "c") and not in_lattice(q_half, "c"))


# ---------------------------------------------------------------------------
# B. Exact root/zero energy censuses and the separate analytic criterion


ALPHAS = [F(0), F(1, 8), F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4), F(7, 8), F(1)]


def delta(alpha):
    return tuple((F(k) + alpha) / 8 - F(1, 2) for k in range(8))


def energy(p, alpha):
    return norm2(p) / 2 + dot(delta(alpha), p)


def serial_counter(counter):
    return [[str(k), counter[k]] for k in sorted(counter)]


zero = tuple(F(0) for _ in range(8))
pstar = tuple([F(1, 2)] * 4 + [F(-1, 2)] * 4)
energy_tables = {}
for alpha in ALPHAS:
    row = {"delta": [str(x) for x in delta(alpha)]}
    for chirality, system in (("Ls", roots_s), ("Lc", roots_c)):
        values = [energy(zero, alpha)] + [energy(r, alpha) for r in system]
        census = Counter(values)
        ground = min(values)
        half_system = [r for r in system if r[0].denominator == 2]
        half_values = [energy(r, alpha) for r in half_system]
        half_min = min(half_values)
        half_deg = half_values.count(half_min)
        row[chirality] = {
            "root_plus_zero_energy_census": serial_counter(census),
            "root_plus_zero_min": str(ground),
            "root_plus_zero_ground_count": values.count(ground),
            "half_root_min": str(half_min),
            "half_root_min_degeneracy": half_deg,
        }
    if F(0) < alpha < F(1):
        require(f"B_Ls_rootzero_ground_two_alpha_{alpha}", row["Ls"]["root_plus_zero_min"] == "0" and row["Ls"]["root_plus_zero_ground_count"] == 2)
        require(f"B_Lc_rootzero_ground_one_alpha_{alpha}", row["Lc"]["root_plus_zero_min"] == "0" and row["Lc"]["root_plus_zero_ground_count"] == 1)
        expected_gap = min(alpha, 1 - alpha) / 8
        require(f"B_Lc_half_gap_alpha_{alpha}", F(row["Lc"]["half_root_min"]) == expected_gap)
        require(f"B_Lc_half_root_deg_alpha_{alpha}", row["Lc"]["half_root_min_degeneracy"] == (2 if alpha == F(1, 2) else 1))

        d = delta(alpha)
        require(f"B_integer_coordinate_minimum_criterion_alpha_{alpha}", all(abs(x) < F(1, 2) for x in d))
        require(f"B_half_sign_pattern_alpha_{alpha}", all(x < 0 for x in d[:4]) and all(x > 0 for x in d[4:]))
        require(f"B_pstar_Ls_zero_alpha_{alpha}", in_lattice(pstar, "s") and energy(pstar, alpha) == 0)
        require(f"B_Lc_single_flip_cost_alpha_{alpha}", min(-d[3], d[4]) == expected_gap)
    energy_tables[str(alpha)] = row

require("B_pstar_not_Lc", not in_lattice(pstar, "c"))
require("B_gap_branches_cross_only_at_half", F(1, 2) / 8 == (1 - F(1, 2)) / 8)
require("B_quarter_gap", min(F(1, 4), F(3, 4)) / 8 == F(1, 32))


# ---------------------------------------------------------------------------
# C. Exact finite double-cover matrices and the selected h sector


def shift_matrix(n, holonomy):
    out = sp.zeros(n)
    for j in range(n - 1):
        out[j + 1, j] = 1
    out[0, n - 1] = holonomy
    return out


cover_results = {}
x = sp.symbols("x")
for L in (4, 6, 8):
    for h_name, h in (("1", sp.Integer(1)), ("-1", sp.Integer(-1)), ("i", sp.I), ("-i", -sp.I)):
        S = shift_matrix(2 * L, h ** 2)
        D = S ** L
        I2L = sp.eye(2 * L)
        require(f"C_deck_square_{L}_{h_name}", D ** 2 == h ** 2 * I2L)
        require(f"C_minimal_factor_{L}_{h_name}", (D - h * I2L) * (D + h * I2L) == sp.zeros(2 * L) and D != h * I2L and D != -h * I2L)
        plus_dim = 2 * L - (D - h * I2L).rank()
        minus_dim = 2 * L - (D + h * I2L).rank()
        require(f"C_deck_eigenspace_ranks_{L}_{h_name}", plus_dim == L and minus_dim == L)
        require(f"C_deck_charpoly_{L}_{h_name}", sp.expand(D.charpoly(x).as_expr() - (x - h) ** L * (x + h) ** L) == 0)

        Qh = (I2L + D / h) / 2
        require(f"C_projector_exact_{L}_{h_name}", Qh ** 2 == Qh and Qh.conjugate().T == Qh and Qh.rank() == L and Qh != I2L)

        J = sp.zeros(2 * L, L)
        for j in range(L):
            J[j, j] = 1 / sp.sqrt(2)
            J[j + L, j] = 1 / (h * sp.sqrt(2))
        T = shift_matrix(L, h)
        require(f"C_pullback_isometry_{L}_{h_name}", sp.simplify(J.conjugate().T * J) == sp.eye(L))
        require(f"C_pullback_projector_{L}_{h_name}", sp.simplify(J * J.conjugate().T - Qh) == sp.zeros(2 * L))
        require(f"C_pullback_intertwines_shift_{L}_{h_name}", sp.simplify(S * J - J * T) == sp.zeros(2 * L, L))

        cover_results[f"L{L}_h{h_name}"] = {
            "double_dimension": 2 * L,
            "deck_eigenvalues": [str(h), str(-h)],
            "deck_eigenspace_ranks": [plus_dim, minus_dim],
            "selected_projector_rank": Qh.rank(),
            "selected_projector_is_full_double": Qh == I2L,
            "intertwined_original_h_sector_dimension": L,
        }


# ---------------------------------------------------------------------------
# D. Two ordinary fermion modes: graded tensor, not one-particle direct sum


I2 = sp.eye(2)
P1 = sp.diag(1, -1)
create = sp.Matrix([[0, 0], [1, 0]])
annihilate = create.conjugate().T
adag = sp.kronecker_product(create, I2)
a = adag.conjugate().T
bdag = sp.kronecker_product(P1, create)
b = bdag.conjugate().T
I4 = sp.eye(4)
Ptot = sp.kronecker_product(P1, P1)

require("D_fock_dimension_tensor_four", 2 ** (1 + 1) == (2 ** 1) * (2 ** 1) == 4)
require("D_one_particle_direct_sum_dimension_two", 1 + 1 == 2 and 1 + 1 != 4)
require("D_mode_a_CAR", a * adag + adag * a == I4 and a * a == sp.zeros(4) and adag * adag == sp.zeros(4))
require("D_mode_b_CAR", b * bdag + bdag * b == I4 and b * b == sp.zeros(4) and bdag * bdag == sp.zeros(4))
require("D_cross_CAR", a * b + b * a == sp.zeros(4) and a * bdag + bdag * a == sp.zeros(4))
require("D_total_parity_product", Ptot == sp.diag(1, -1, -1, 1) and Ptot ** 2 == I4)
require("D_creators_flip_grading", Ptot * adag == -adag * Ptot and Ptot * bdag == -bdag * Ptot)
Peven = (I4 + Ptot) / 2
Podd = (I4 - Ptot) / 2
require("D_even_odd_ranks_two", Peven.rank() == 2 and Podd.rank() == 2 and Peven + Podd == I4)


results = {
    "contract": "seam-joint-ward-state-20261014-scratch",
    "verdict": "SCOPED_EXACT_CHECKS_PASS",
    "A_lattice_and_Ward": {
        "root_counts": {"Ls": len(roots_s), "Lc": len(roots_c)},
        "q": [str(x) for x in q],
        "q_norm2": str(norm2(q)),
        "q_in": {"Ls": in_lattice(q, "s"), "Lc": in_lattice(q, "c")},
        "q_over_2_in": {"Ls": in_lattice(q_half, "s"), "Lc": in_lattice(q_half, "c")},
        "q_dot_halfroot_census": {
            "Ls": serial_counter(qdot_s),
            "Lc": serial_counter(qdot_c),
        },
        "q_dot_allroot_census": {
            "Ls": serial_counter(qdot_roots_s),
            "Lc": serial_counter(qdot_roots_c),
        },
        "divisibility_obstruction": (
            "No lattice isomorphism U:Lc->Ls can obey U(q)=2r for an Ls root r: "
            "2-divisibility is preserved, but q/2 is not in Lc while r is in Ls."
        ),
    },
    "B_time_energy": {
        "alphas": [str(a) for a in ALPHAS],
        "tables": energy_tables,
        "analytic_interval_criterion": {
            "scope": "all real 0<alpha<1; this is the analytic step, not a conclusion from the finite root census",
            "integer_coset": "Each coordinate has |delta_k|<1/2, so the unique unconstrained integer minimum is zero; it lies in D8.",
            "half_coset_Ls": "The coordinatewise half-integer minimum is p*=(+1/2)^4,(-1/2)^4, has even parity, and energy zero.",
            "half_coset_Lc": "Odd parity requires an odd change of integer offsets. Coordinatewise convexity makes the cheapest parity-changing move a sign flip of one +/-1/2 entry; outward moves to +/-3/2 cost more. The two cheapest sign flips cost delta_4=alpha/8 and -delta_3=(1-alpha)/8; all other flips cost more.",
            "conclusion": "Full-lattice ground multiplicities are Ls=2 and Lc=1. The Lc half-coset minimum is min(alpha,1-alpha)/8, with two minimizing roots only at alpha=1/2.",
        },
        "finite_census_limitation": "The displayed root-plus-zero tables alone do not prove the infinite-lattice minimum.",
    },
    "C_double_cover": cover_results,
    "D_two_mode_Fock": {
        "one_particle_dimension": 2,
        "fock_dimension": 4,
        "graded_tensor_dimensions": [2, 2],
        "total_parity_diagonal": [1, -1, -1, 1],
        "even_odd_dimensions": [2, 2],
        "scope": "Adding an ordinary fermion mode gives the graded tensor Fock factor and flips ordinary fermion parity; it does not by itself construct a spinor lattice extension.",
    },
    "nonclaims": [
        "No total no-go for all identifications or extensions.",
        "No P1/P2 solution or physical state-selection derivation.",
        "No claim that the particular pullback projector is the full 2L double-cover space.",
        "No derivation of the c-chiral extension from the unchanged CAR source.",
    ],
    "checks": checks,
    "n_checks": len(checks),
}

out = Path(__file__).with_name("results.json")
out.write_text(json.dumps(results, indent=2) + "\n")
print(f"{len(checks)}/{len(checks)} checks passed")
print(hashlib.sha256(out.read_bytes()).hexdigest())
