"""H6h26/H6h27 certificates and parent selection in the existing lane.

There is no new search engine. The explicit microscopic quadratic rule is
passed to genesis_rules.search, and its hypothesis is reported next to the
selection. A selected response is not promoted to a physical process.
"""
from functools import lru_cache

import sympy as sp

from genesis_rules import search
from process_boundary_fourpoint import boundary_fourpoint_audit
from process_depth_control import depth_control_audit
from process_flavour_obstructions import flavour_obstruction_audit
from process_gaussian_parent import gaussian_parent_audit
from process_metric_potential_certificates import metric_potential_audit
from process_source_bimetric import source_bimetric_audit
from process_spectral_no_go import spectral_no_go_audit
from process_weighted_star import weighted_star_audit


def parent_dynamics_universe():
    """Positive polynomial edge weights with the same exact geometric state.

    A current-quadratic pair operator has hopping proportional to n+1.
    Polynomial identities decide that condition for all n, not a mesh.
    Positivity follows from the displayed nonnegative coefficients in n.
    """
    n = sp.Symbol("n", nonnegative=True, integer=True)
    weights = (("flat Hardy", sp.Integer(1)), ("quadratic CCR parent", n+1),
               ("quadratic depth weight", (n+1)**2), ("affine alternative", 2*n+1))
    objects = []
    for label, weight in weights:
        polynomial = sp.Poly(weight, n)
        positive = all(c >= 0 for c in polynomial.all_coeffs()) and weight.subs(n, 0) > 0
        if not positive:
            raise ValueError("Uncertified polynomial edge weight")
        clauses = {"positive_response", "same_KMS_purification", "one_ground_state", "quarter_symmetry"}
        quadratic = sp.simplify(sp.diff(weight/(n+1), n)) == 0
        if quadratic:
            clauses.add("quadratic_in_the_two_current_modes")
        if polynomial.degree() == 0:
            clauses.add("bounded_mass_operator")
        objects.append({"label": label, "weight": str(weight), "clauses": clauses,
                        "quadratic_proof": "d/dn [w(n)/(n+1)]="+str(sp.simplify(sp.diff(weight/(n+1), n))),
                        "physical_solution": False})
    pool = ["positive_response", "same_KMS_purification", "one_ground_state", "quarter_symmetry",
            "quadratic_in_the_two_current_modes", "bounded_mass_operator"]
    result = search(objects, pool, lambda obj: 0, lambda obj: obj["label"])
    result["candidates"] = [{**obj, "clauses": sorted(obj["clauses"])} for obj in objects]
    result["bound"] = "four positive polynomial depth weights; no size preference"
    result["unqualified_unique_process"] = False
    result["status"] = "CONDITIONAL_PARENT_DYNAMICS_SEARCH"
    result["extra_rule_scope"] = (
        "Quadratic in exactly two identified microscopic CCR modes is an extra physical condition. "
        "The separate Gaussian certificate proves uniqueness up to scale in that full declared quadratic class; "
        "The full theorem also assumes conservation of N_L-N_R and exchange of the two modes. "
        "The finite enumeration alone does not prove it.")
    return result


@lru_cache(maxsize=1)
def source_dynamics_audit():
    gaussian = gaussian_parent_audit()
    weighted = weighted_star_audit()
    source = source_bimetric_audit()
    flavour = flavour_obstruction_audit()
    spectral = spectral_no_go_audit()
    boundary = boundary_fourpoint_audit()
    depth = depth_control_audit()
    potential = metric_potential_audit()
    checks = {
        "current_and_Gaussian_operator_checks": gaussian["all_checks_pass"],
        "all_four_weighted_metric_sectors": weighted["all_checks_pass"],
        "physical_source_and_quartic_checks": source["all_checks_pass"],
        "exact_flavour_obstructions": flavour["all_exact_checks_pass"],
        "connected_boundary_fourpoint": boundary["all_checks_pass"],
        "conditional_local_depth_control": depth["all_checks_pass"],
        "full_HR_potential_class": potential["all_checks_pass"],
        "spectral_action_sum_rule_beta_identity": spectral["analytic_identities"]["sum_rule_factor"] == "0",
    }
    return {
        "status": "H6H27_SOURCE_DYNAMICS_AUDITED",
        "base_commit": "b5cc3d4070a1d869cafc3c414cf0cb8b2d4f1401",
        "initial_base_commit": "0172a74bdea0d3f751372813598dfe3f777a5930",
        "previous_certificate_commit": "8c076733cf1013e2ec18acfdca83cb5447aaa58d",
        "checks": checks, "all_checks_pass": all(checks.values()),
        "physical_solution": False, "ledger_promotions": [],
        "gaussian": gaussian, "weighted_metric_parent": weighted,
        "source_bimetric": source, "flavour": flavour,
        "boundary_fourpoint": boundary, "depth_control": depth,
        "metric_potential_class": potential,
        "spectral_action_scope": spectral,
        "parent_selection": parent_dynamics_universe(),
        "normalization_freedom": {
            "transformation": "M_v -> c*M_v, k_e -> c^2*k_e for every c>0",
            "preserved": "canonical mass operator, spectral shapes, q, zero-mode direction and normalized kernel",
            "changed": "physical source factors 1/M_v -> (1/c)/M_v and G -> G/c^2",
            "consequence": "The structural and normalized spectral data do not fix G without an independent physical source normalization.",
        },
        "not_claimed": [
            "P1 or spacetime derived from no physical inputs",
            "a nonlinear quantum field integral for seam, chiral matter, gauge fields, Higgs and gravity",
            "a prediction of G, mass scale or quartic metric parameter",
            "the TM1 angle or CP branch selected by TFPT dynamics",
        ],
    }
