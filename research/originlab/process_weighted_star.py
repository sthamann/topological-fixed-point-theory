"""Depth weighted metric parent for the complete quarter transfer.

This is a conditional finite classical construction. Einstein kinetic terms,
pairwise Hassan--Rosen potentials and their graph are additional inputs.
Only the central symmetric sector has the two oscillator interpretation.
All three other quarter sectors remain in the operator and are certified.
"""
from functools import lru_cache

import sympy as sp

from process_star_metric import star_metric_data


def weighted_star_data(arm_length, q, M_c=1, scale=1):
    """Multiply the stiffness of an edge at depth n by n+1.

    The geometric Planck masses and the one exact null vector are unchanged.
    The scale has units of mass squared when this is a 4D mass operator;
    identifying it with a microscopic oscillator energy is a separate map.
    """
    d = star_metric_data(arm_length, q, M_c, scale)
    m = arm_length
    W = sp.diag(*(sp.Integer(n+1) for n in range(m) for _ in range(4)))
    B, A, U = (d[k] for k in (
        "canonical_incidence", "native_A_rect", "U_vertex_fourier"))
    d["depth_weights"] = W
    d["edge_fp_weights"] = W*d["edge_fp_weights"]
    d["mass_matrix"] = scale*B.T*W*B
    d["native_weighted_mass"] = scale*A.T*W*A
    d["fourier_mass"] = sp.simplify(U.conjugate().T*d["mass_matrix"]*U)
    d["blocks"] = [d["native_weighted_mass"].extract(
        list(range(r, 4*m+1, 4)), list(range(r, 4*m+1, 4))) for r in range(4)]
    return d


def radial_finite_matrix(length, q, scale=1):
    """Physical graph truncation with an exact zero mode, not compression."""
    if isinstance(length, bool) or not isinstance(length, int) or length < 1:
        raise ValueError("length must be a positive integer")
    A = sp.zeros(length, length+1)
    for n in range(length):
        A[n, n], A[n, n+1] = q, -1
    return scale*A.T*sp.diag(*range(1, length+1))*A


@lru_cache(maxsize=1)
def weighted_star_audit():
    q, M_c, scale = sp.symbols("q M_c scale", positive=True)
    d = weighted_star_data(3, q, M_c, scale)
    nv, ne = d["metrics"], d["edges"]
    U, V, W, A, B, D = (d[k] for k in (
        "U_vertex_fourier", "V_edge_fourier", "depth_weights",
        "native_A_rect", "canonical_incidence", "planck_masses"))
    radial = radial_finite_matrix(3, q, scale)
    exact_zero = lambda x, n, m: sp.simplify(x) == sp.zeros(n, m)
    reconstructed = D.inv()*d["incidence"].T*d["edge_fp_weights"]*d["incidence"]*D.inv()
    checks = {
        "depth_weight_commutes_with_arm_Fourier": exact_zero(W*V-V*W, ne, ne),
        "first_order_full_quarter_intertwiner": exact_zero(B*U-V*A, ne, nv),
        "weighted_full_mass_intertwiner": exact_zero(d["fourier_mass"]-d["native_weighted_mass"], nv, nv),
        "actual_EH_and_edge_weights_reproduce_mass": exact_zero(reconstructed-d["mass_matrix"], nv, nv),
        "symmetric_sector_is_Gaussian_pair_parent": exact_zero(d["blocks"][0]-radial, 4, 4),
        "three_other_sectors_are_Dirichlet_restrictions": all(
            exact_zero(block-radial[1:, 1:], 3, 3) for block in d["blocks"][1:]),
        "single_zero_by_positive_edge_weights_and_unit_minor": A[:, 1:].det() == 1,
        "same_geometric_null_vector": exact_zero(d["mass_matrix"]*d["zero"], nv, 1),
        "full_quarter_symmetry": exact_zero(d["rotation"]*d["mass_matrix"]-d["mass_matrix"]*d["rotation"], nv, nv),
        "all_physical_sources_have_same_massless_residue": exact_zero(
            d["physical_pole_matrix"]-sp.ones(nv)/d["M_eff_squared"], nv, nv),
        "massive_modes_keep_universal_massless_EH_overlap": exact_zero(
            d["cubic_zero_weight_matrix"]-sp.eye(nv)/sp.sqrt(d["M_eff_squared"]), nv, nv),
        "physical_endpoint_is_not_infinite_operator_compression": radial[-1, -1] == 3*scale,
    }
    return {
        "status": "CONDITIONAL_FINITE_CLASSICAL_WEIGHTED_METRIC_PARENT",
        "physical_solution": False,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "action": "sum_v M_v^2/2 integral sqrt(-g_v) R[g_v] - sum_edges k_e integral sqrt(-g_p) V_tau(sqrt(g_p^-1*g_c)) + sum_v S_matter,v[g_v]",
        "potential": "V_tau(X)=(-4+e1+e3-4*e4)/2+(tau/2)*det(I-X)",
        "planck_weights": "M_center=M_c; M_(j,n)=M_c*q^(n+1)/2",
        "edge_stiffness": "k_(j,n)=scale*(n+1)*M_(j,n)^2",
        "scale_dimension": "mass squared in the 4D mass operator; oscillator energy needs an explicit unit map",
        "all_four_sectors_retained": True,
        "central_sector": "K=scale*A_q^dagger*(N+1)*A_q",
        "other_three_sectors": "Dirichlet restriction of this Jacobi operator to radial n>=1; not additional copies of the two oscillator pair sector",
        "radial_infinite_entries": "K_nn=scale*[n+q^2*(n+1)], K_n,n+1=-scale*q*(n+1)",
        "self_adjoint_scope": "0<q<1. The Jacobi Carleman sum sum_n 1/[scale*q*(n+1)] diverges; the finitely supported symmetric operator has a unique self-adjoint closure.",
        "radial_spectrum": "Delta*n, n=0,1,..., Delta=scale*(1-q^2)",
        "radial_spectral_weights_at_center": "(1-q^2)*q^(2*n)",
        "other_sector_spectrum": "One simple eigenvalue lambda_n in each open interval (Delta*n, Delta*(n+1)); it is a zero of the radial Stieltjes resolvent on the negative real axis.",
        "interlacing_proof": "Delete the cyclic boundary vector from a compact-resolvent Jacobi operator with strictly positive weights at every simple eigenvalue. The scalar resolvent is strictly decreasing between its poles and changes from +infinity to -infinity.",
        "total_counting_bound": "For E>=0, N_full(E)<=4*floor(E/Delta)+4. N counts eigenvalues of mass squared, with multiplicity; no infinite number of modes below fixed finite E.",
        "universal_massless_coupling": "eta_v/M_v=1/M_eff for every vertex",
        "finite_newton_constant": "G_eff=1/[8*pi*M_c^2*sum_(n=0)^m q^(2*n)]",
        "infinite_newton_constant": "G_eff=(1-q^2)/(8*pi*M_c^2)",
        "G_is_predicted": False,
        "finite_vs_infinite": "Finite graphs have a physical terminal diagonal scale*m. Infinite central spectrum is an exact limit, not the exact spectrum of any finite graph.",
        "quantum_completion_proved": False,
        "remaining_obstructions": [
            "M_c, scale and tau remain free",
            "the physical metric and matter dictionary is an additional construction",
            "arbitrarily small local Planck coefficients along infinite arms require a controlled quantum cutoff or completion",
            "Kallen-Lehmann positivity here certifies the scalar TT coefficient, not the complete tensor and interacting quantum theory",
            "chirality, gauge and Higgs interactions are not derived by this metric construction",
        ],
        "provenance": {
            "previous_star_commit": "8c076733cf1013e2ec18acfdca83cb5447aaa58d",
            "new_step": "depth weight n+1 selected conditionally by the two CCR mode quadratic parent criterion",
        },
    }
