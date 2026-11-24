"""Exact four-dimensional HR potential certificates, without random probes.

The fields g,f, square-root branch near g=f=eta, Einstein coefficients and
remaining action are fixed.  These are off-shell potential-vertex statements,
not equivalence modulo field redefinitions or a quantum positivity theorem.
For symbolic parameters clauses mean identities, not exclusions of special
parameter values.  Existing numerical perturbation samples may be retained by
the caller as diagnostics but do not decide any clause here.
"""
from collections.abc import Mapping
import sympy as sp


BETA_V = (sp.Integer(-2), sp.Rational(1, 2), sp.Integer(0), sp.Rational(1, 2), sp.Integer(-2))
BETA_DET = tuple(map(sp.Integer, (1, -1, 1, -1, 1)))
BETA_CUB = tuple(map(sp.Integer, (-2, 1, 0, -1, 2)))
FLAT_ROWS = sp.Matrix([[1, 3, 3, 1, 0], [0, 1, 3, 3, 1]])
MASS_ROW = sp.Matrix([[0, 1, 2, 1, 0]])
PALINDROME_ROWS = sp.Matrix([[1, 0, 0, 0, -1], [0, 1, 0, -1, 0]])
# e_n(I+E)=sum_(r=0)^n binomial(4-r,n-r)*e_r(E).
SHIFT_MATRIX = sp.Matrix(5, 5, lambda r, n: sp.binomial(4-r, n-r) if n >= r else 0)
CLAUSE_NAMES = ("flat_background", "fierz_pauli_mass", "massless_direction_free",
                "undirected_edge", "cubic_as_V", "quartic_as_V")


def _exact_beta(beta):
    try:
        values = tuple(beta)
    except TypeError as exc:
        raise ValueError("beta must contain exactly five exact coefficients") from exc
    if len(values) != 5 or any(isinstance(x, (bool, float)) for x in values):
        raise ValueError("beta must contain exactly five exact coefficients, without floats")
    values = tuple(sp.sympify(x) for x in values)
    if any(x.has(sp.Float) or x in (sp.oo, -sp.oo, sp.zoo, sp.nan) for x in values):
        raise ValueError("beta coefficients must be finite exact expressions")
    if any(x.is_real is False for x in values):
        raise ValueError("physical HR coefficients must be real")
    return sp.Matrix(values)


def _zero(values):
    return all(sp.simplify(x) == 0 for x in values)


def metric_vertex_probe_matrix(order):
    """Deterministic coefficient rows from full homogeneous metric vertices.

    Row 0: g=f=diag(1+eps,1,1,1), which isolates c0.
    Row r>=1: g=I,f=I+eps*diag(d1,...,d4), monomial
    d1**(order-r+1)*d2*...*dr, which isolates c_r.
    Here c=SHIFT_MATRIX*beta.  The half-binomial factors never vanish.
    """
    if isinstance(order, bool) or not isinstance(order, int) or order < 0:
        raise ValueError("order must be a nonnegative integer")
    rows = [sp.binomial(sp.Rational(1, 2), order)*SHIFT_MATRIX[0, :]]
    for r in range(1, min(order, 4)+1):
        factor = sp.binomial(sp.Rational(1, 2), order-r+1)*sp.Rational(1, 2)**(r-1)
        rows.append(factor*SHIFT_MATRIX[r, :])
    return sp.Matrix.vstack(*rows)


def metric_potential_candidates():
    def shifted(cubic=0, determinant=0):
        return tuple(v+cubic*c+determinant*d for v, c, d in zip(BETA_V, BETA_CUB, BETA_DET))
    half = sp.Rational(1, 2)
    return (("V", BETA_V), ("V + det(I-X)/2", shifted(determinant=half)),
            ("V - det(I-X)/2", shifted(determinant=-half)),
            ("V + cubic/2", shifted(cubic=half)),
            ("V + cubic/2 + det(I-X)/2", shifted(cubic=half, determinant=half)),
            ("beta1 only", tuple(map(sp.Integer, (0, 1, 0, 0, 0)))))


def metric_potential_certificate(beta):
    b = _exact_beta(beta)
    delta = b-sp.Matrix(BETA_V)
    flat = FLAT_ROWS*b
    mass = sp.simplify((MASS_ROW*b)[0])
    shifted = SHIFT_MATRIX*b
    cubic_residual = metric_vertex_probe_matrix(3)*delta
    quartic_residual = metric_vertex_probe_matrix(4)*delta
    flat_ok, mass_one = _zero(flat), sp.simplify(mass-1) == 0
    alpha, tau = sp.simplify(delta[1]+delta[2]), sp.simplify(delta[2])
    decomposition = delta-alpha*sp.Matrix(BETA_CUB)-tau*sp.Matrix(BETA_DET)
    same = _zero(delta)
    truths = {
        "flat_background": flat_ok,
        "fierz_pauli_mass": flat_ok and mass_one,
        "massless_direction_free": sp.simplify(shifted[0]) == 0,
        "undirected_edge": _zero(PALINDROME_ROWS*b),
        "cubic_as_V": _zero(cubic_residual),
        "quartic_as_V": _zero(quartic_residual),
    }
    return {
        "beta": [str(x) for x in b], "clauses": sorted(k for k, v in truths.items() if v),
        "clause_truths": truths,
        "flat_residuals": [str(sp.simplify(x)) for x in flat],
        "fierz_pauli_combination": str(mass),
        "fierz_pauli_combination_is_physical_under_flatness": flat_ok,
        "common_metric_coefficient": str(sp.simplify(shifted[0])),
        "coefficients_in_e_r_of_X_minus_I": [str(sp.simplify(x)) for x in shifted],
        "cubic_probe_residuals": [str(sp.simplify(x)) for x in cubic_residual],
        "quartic_probe_residuals": [str(sp.simplify(x)) for x in quartic_residual],
        "flat_unit_mass_deformation": {
            "applicable": flat_ok and mass_one, "alpha_cubic": str(alpha), "tau_determinant": str(tau),
            "exact_decomposition": _zero(decomposition),
            "cubic_match_condition": "alpha=0", "quartic_match_conditions": "alpha=0 and tau=0",
        },
        "all_orders_same_potential_as_V": same,
        "quartic_injectivity_check": truths["quartic_as_V"] == same,
        "undirected_scope": "Exchange of the pair potential density; unequal Planck masses need not give an exchange symmetry of the full theory.",
        "certificate_method": "exact linear identities and deterministic tensor-monomial probes, no numerical samples",
    }


def metric_potential_class_proof():
    flat_mass = FLAT_ROWS.col_join(MASS_ROW)
    undirected = flat_mass.col_join(PALINDROME_ROWS)
    cubic = metric_vertex_probe_matrix(3)
    quartic = metric_vertex_probe_matrix(4)
    directions = sp.Matrix.hstack(sp.Matrix(BETA_CUB), sp.Matrix(BETA_DET))
    z = sp.Symbol("z", real=True)
    cut = sp.expand(sum(BETA_DET[n]*sp.binomial(4, n)*z**n for n in range(5)))
    checks = {
        "two_independent_flat_conditions": FLAT_ROWS.rank() == 2,
        "flat_and_unit_mass_leave_two_directions": flat_mass.rank() == 3 and directions.rank() == 2 and flat_mass*directions == sp.zeros(3, 2),
        "palindrome_leaves_one_direction": undirected.rank() == 4 and undirected*sp.Matrix(BETA_DET) == sp.zeros(5, 1),
        "shift_matrix_is_invertible": SHIFT_MATRIX.det() == 1,
        "common_coefficient_is_sum_of_flat_rows": SHIFT_MATRIX[0, :] == FLAT_ROWS[0, :]+FLAT_ROWS[1, :],
        "cubic_vertex_has_exact_determinant_kernel": cubic.rank() == 4 and cubic*sp.Matrix(BETA_DET) == sp.zeros(4, 1),
        "quartic_vertex_is_injective_on_all_five_coefficients": quartic.rank() == 5,
        "determinant_starts_at_fourth_order": SHIFT_MATRIX*sp.Matrix(BETA_DET) == sp.Matrix([0, 0, 0, 0, 1]),
        "cubic_direction_is_e3_plus_2e4": SHIFT_MATRIX*sp.Matrix(BETA_CUB) == sp.Matrix([0, 0, 0, 1, 2]),
        "tree_edge_cut_witness": sp.expand(cut-(1-z)**4) == 0,
    }
    return {
        "checks": checks, "all_checks_pass": all(checks.values()),
        "ranks": {"flat": 2, "flat_mass": 3, "flat_mass_palindrome": 4, "full_cubic": 4, "full_quartic": 5},
        "general_family": "beta=V+alpha*(-2,1,0,-1,2)+tau*(1,-1,1,-1,1)",
        "shift_expansion": "V_beta(I+E)=sum_r c_r e_r(E), c_r=sum_(n>=r) binomial(4-r,n-r)*beta_n",
        "cubic_rule": "Full homogeneous cubic matching iff delta_beta=tau*DET; in the flat/unit-mass class edge palindromie is equivalent.",
        "quartic_rule": "Full homogeneous quartic matching iff beta=V, on the entire five-coefficient HR class in fixed metric coordinates.",
        "all_orders_rule": "beta=V gives identical potential densities and hence every vertex. For each homogeneous order k>=4 the displayed probe rows have rank five.",
        "quartic_probe_factors_c0_to_c4": ["-5/128", "-5/128", "1/32", "-1/32", "1/16"],
        "potential_density_quartic_difference": "tau*det(B-A)/16 for g=I+eps*A, f=I+eps*B",
        "canonical_action_quartic_difference": "+tau*mu^2*det((B_transfer*h)_edge)/M_child^2 in Euclidean S_E=+k*V and g=I+2h/M; Lorentzian action -k*V has the opposite sign for mixed-index determinants (det eta=-1 converts to the covariant-tensor convention)",
        "parameter_dictionary": "tau is the coefficient of det(I-X); the earlier V+(t/2)*det(I-X) uses tau=t/2",
        "scope": "Four-dimensional HR pair potentials, fixed fields and fixed remaining action; not general higher-derivative actions, field-redefinition equivalence, or quantum positivity.",
    }


def star_quartic_orbits(arm_length):
    """Each depth is a distinct C4 orbit; depth homogeneity is a new input."""
    if isinstance(arm_length, bool) or not isinstance(arm_length, int) or arm_length < 1:
        raise ValueError("arm_length must be a positive integer")
    m = arm_length
    constraints = sp.zeros(4*m)
    for n in range(m):
        for j in range(4):
            constraints[4*n+j, 4*n+j] = 1
            constraints[4*n+j, 4*n+(j+1) % 4] = -1
    return {
        "arm_length": m, "edges": 4*m,
        "independent_tau_before_mu4": 4*m,
        "mu4_constraint_rank": constraints.rank(),
        "independent_tau_after_mu4": 4*m-constraints.rank(),
        "edge_orbits": [[[j, n] for j in range(4)] for n in range(m)],
        "mu4_rule": "tau_(j,n)=tau_n; no relation between different n follows",
        "homogeneous_tau_across_depth_is_extra": True,
        "independence_proof": "Cut a chosen edge of the tree. Set all metrics on one component to g and all on the other to f=z^2*g. Every other determinant deformation vanishes at X=I; the chosen edge alone gives k_e*tau_e*sqrt(|g|)*(1-z)^4. Thus edge coefficients are independent.",
    }


def metric_potential_audit(candidates=None):
    """JSON-compatible default audit; no import of process_rules is needed.

    candidates may be a mapping name->exact beta, pairs (name,beta), or
    objects with name and beta keys.  Floats must be explicitly converted
    by the caller according to their declared rational meaning.
    """
    if candidates is None:
        candidates = metric_potential_candidates()
    elif isinstance(candidates, Mapping):
        candidates = candidates.items()
    rows = []
    for item in candidates:
        name, beta = (item["name"], item["beta"]) if isinstance(item, Mapping) else item
        rows.append({"name": str(name), **metric_potential_certificate(beta)})
    proof = metric_potential_class_proof()
    return {
        "status": "EXACT_FIVE_COEFFICIENT_HR_POTENTIAL_CLASSIFICATION",
        "reference": [str(x) for x in BETA_V], "candidates": rows,
        "class_proof": proof, "star_quartic_freedom": star_quartic_orbits(3),
        "all_checks_pass": proof["all_checks_pass"] and all(r["quartic_injectivity_check"] for r in rows),
        "samples_decide_clauses": False,
        "existing_samples_role": "Retain separately as numerical diagnostics; exact certificates decide clauses.",
        "physical_solution": False,
    }
