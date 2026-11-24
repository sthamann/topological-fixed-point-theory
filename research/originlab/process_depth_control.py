"""Conditional local depth check for flat and depth weighted metric stars.

At arm depth d>=1, M_d=M_c*q**d/2. The parent-to-child ratio is q/2
at the central edge and q afterwards. If k_d=mu**2*w_d*M_d**2,
w_d=1 or d, the isolated two-metric mass obeys
m_d**2=mu**2*w_d*(1+rho_d**2), M_rel=M_d/sqrt(1+rho_d**2).
The conventional local Lambda3 estimate is (m_d**2*M_rel)**(1/3).
Thus Lambda3/m_d>margin iff
q**(2*d)/(4*(mu/M_c)**2*w_d*(1+rho_d**2)**2)>margin**6.

These are local power-counting checks, not the strong coupling scale of
the full graph and not a dynamically derived end boundary. A decreasing
local Planck mass does not by itself prove a universal UV obstruction.
The reported native depth counts use rational bounds, not float probes.
Bounds on exp(-2*pi) follow from Machin's identity, alternating arctan
remainders, and the positive exponential series with a geometric tail.
"""
from fractions import Fraction
from functools import lru_cache

import sympy as sp


def _fraction(value):
    if isinstance(value, (bool, float)):
        raise TypeError("use exact rational input")
    return Fraction(value)


def _atan_bounds(x, terms=12):
    partial = sum(((-1)**k*x**(2*k+1)/(2*k+1) for k in range(terms)), Fraction())
    next_term = (-1)**terms*x**(2*terms+1)/(2*terms+1)
    return min(partial, partial+next_term), max(partial, partial+next_term)


def _exp_bounds(x, degree=60):
    term, partial = Fraction(1), Fraction(1)
    for k in range(1, degree+1):
        term *= x/k
        partial += term
    first_omitted = term*x/(degree+1)
    ratio_bound = x/(degree+2)
    if not 0 <= ratio_bound < 1:
        raise ValueError("exponential remainder requires degree+2 > x >=0")
    return partial, partial+first_omitted/(1-ratio_bound)


@lru_cache(maxsize=1)
def native_q_bounds():
    """Prove 1867442/10^9 < exp(-2*pi) < 1867443/10^9."""
    a_lo, a_hi = _atan_bounds(Fraction(1, 5))
    b_lo, b_hi = _atan_bounds(Fraction(1, 239))
    pi_lo, pi_hi = 16*a_lo-4*b_hi, 16*a_hi-4*b_lo
    exp_lo, _ = _exp_bounds(2*pi_lo)
    _, exp_hi = _exp_bounds(2*pi_hi)
    lower, upper = Fraction(1867442, 10**9), Fraction(1867443, 10**9)
    return {
        "lower": lower, "upper": upper,
        "lower_proved": lower < 1/exp_hi,
        "upper_proved": 1/exp_lo < upper,
        "proof": "Machin pi=16 atan(1/5)-4 atan(1/239); alternating bounds; degree60 exponential with positive geometric tail",
    }


def depth_certificate(q_lower, q_upper, mu_over_Mc, *, weighted=False,
                      margin=1, max_depth=200):
    """Exact interval certificate for the local criterion Lambda3/m > margin.

    The count is of complete arm layers d=1,...,D, not vertices, massive
    modes, or the zero-based source index n=d-1. An undecided interval is
    reported as such and never converted to a floating point predicate.
    """
    lo, hi, ratio, margin = map(_fraction, (q_lower, q_upper, mu_over_Mc, margin))
    if not 0 < lo <= hi < 1 or ratio <= 0 or margin <= 0:
        raise ValueError("require 0<q_lower<=q_upper<1 and positive scales")
    if isinstance(max_depth, bool) or not isinstance(max_depth, int) or max_depth < 1:
        raise ValueError("max_depth must be a positive integer")
    if not isinstance(weighted, bool):
        raise TypeError("weighted must be boolean")
    last_lower = None
    threshold = margin**6
    for d in range(1, max_depth+1):
        w = d if weighted else 1
        alpha = Fraction(1, 2) if d == 1 else Fraction(1)
        lower = lo**(2*d)/(4*ratio**2*w*(1+(alpha*hi)**2)**2)
        upper = hi**(2*d)/(4*ratio**2*w*(1+(alpha*lo)**2)**2)
        if lower > threshold:
            last_lower = lower
            continue
        status = "CERTIFIED_LOCAL_DEPTH" if upper <= threshold else "UNDECIDED_INTERVAL"
        return {
            "status": status, "weighted": weighted,
            "mu_over_Mc": str(ratio), "required_Lambda3_over_m": str(margin),
            "strict_criterion": "Lambda3/m_edge > margin",
            "certified_layers": d-1, "last_zero_based_source_index": d-2,
            "first_failed_depth": d if upper <= threshold else None,
            "first_unresolved_depth": d if upper > threshold else None,
            "last_pass_lower_squared_Mrel_over_m": str(last_lower) if last_lower is not None else None,
            "next_upper_squared_Mrel_over_m": str(upper),
            "threshold_squared_Mrel_over_m": str(threshold),
            "all_deeper_fail_after_first_failure": upper <= threshold,
            "monotonicity": "At fixed q the squared ratio falls by q^2*w_d/w_(d+1) in the interior; the first step has the additional factor [(1+q^2/4)/(1+q^2)]^2<1",
            "scope": "Local isolated-edge power counting only; no collective graph cutoff or dynamical termination is certified",
        }
    return {
        "status": "DEPTH_LIMIT_REACHED", "weighted": weighted,
        "certified_layers": max_depth, "first_failed_depth": None,
        "scope": "A longer local check is needed; this is not a proof of infinite control",
    }


def edge_depth_data(depth, q, M_c=1, mu=1, *, weighted=False):
    """Exact expressions; symbolic positive scales may be used for identities.

    The mass is the one nonzero isolated-edge FP eigenvalue, not an
    assertion that each edge is itself a global graph eigenmode.
    """
    if isinstance(depth, bool) or not isinstance(depth, int) or depth < 1:
        raise ValueError("depth is the positive arm-layer index")
    if not isinstance(weighted, bool):
        raise TypeError("weighted must be boolean")
    vals = []
    for name, value in (("q", q), ("M_c", M_c), ("mu", mu)):
        if isinstance(value, (bool, float)):
            raise TypeError(f"{name} must be exact")
        value = sp.sympify(value)
        if value.has(sp.Float) or value.is_positive is not True:
            raise ValueError(f"{name} must be provably positive and exact")
        vals.append(value)
    q, M_c, mu = vals
    if q.is_number and (1-q).is_positive is not True:
        raise ValueError("require q<1")
    w = sp.Integer(depth if weighted else 1)
    child = M_c*q**depth/2
    parent = M_c if depth == 1 else M_c*q**(depth-1)/2
    rho = sp.cancel(child/parent)
    relative = child/sp.sqrt(1+rho**2)
    stiffness = mu**2*w*child**2
    mass_squared = mu**2*w*(1+rho**2)
    return {
        "depth": depth, "source_index": depth-1, "weight": w,
        "M_child": child, "M_parent": parent, "rho": rho,
        "M_relative": relative, "edge_stiffness": stiffness,
        "isolated_edge_mass_squared": mass_squared,
        "Lambda3_cubed": sp.simplify(mass_squared*relative),
        "squared_Mrel_over_m": sp.cancel(relative**2/mass_squared),
        "quartic_relative_field_coefficient_per_tau": sp.simplify(stiffness/(2*relative**4)),
        "quartic_transfer_row_coefficient_per_tau": sp.simplify(mu**2*w/(2*child**2)),
    }


def truncated_preparation_data(depth, q, scale=1):
    """Positive two-mode finite parent with the physical terminal entry.

    Both oscillators contain levels 0,...,depth. They are finite matrices,
    so their commutator has a top-level defect and they are not a regular
    CCR representation. The chosen depth is an additional input.
    """
    if isinstance(depth, bool) or not isinstance(depth, int) or depth < 1:
        raise ValueError("depth must be a positive integer")
    exact = []
    for name, value in (("q", q), ("scale", scale)):
        if isinstance(value, (bool, float)):
            raise TypeError(f"{name} must be exact")
        value = sp.sympify(value)
        if value.has(sp.Float) or value.is_positive is not True:
            raise ValueError(f"{name} must be exact and provably positive")
        exact.append(value)
    q, scale = exact
    if q.is_number and (1-q).is_positive is not True:
        raise ValueError("require q<1")
    dim = depth+1
    a = sp.zeros(dim)
    for n in range(1, dim):
        a[n-1, n] = sp.sqrt(n)
    I = sp.eye(dim)
    left, right = sp.kronecker_product(a, I), sp.kronecker_product(I, a)
    dL, dR = left-q*right.T, right-q*left.T
    parent = scale*(dL.T*dL+dR.T*dR)/2
    pair_indices = [n*(dim+1) for n in range(dim)]
    vacuum = sp.zeros(dim**2, 1)
    for n, index in enumerate(pair_indices):
        vacuum[index] = q**n
    A = sp.zeros(depth, dim)
    for n in range(depth):
        A[n, n], A[n, n+1] = q, -1
    top = sp.zeros(dim)
    top[-1, -1] = 1
    return {
        "single_mode_lowering": a,
        "single_mode_commutator": a*a.T-a.T*a,
        "expected_commutator": I-(depth+1)*top,
        "full_parent": parent,
        "paired_parent": parent.extract(pair_indices, pair_indices),
        "physical_radial_parent": scale*A.T*sp.diag(*range(1, depth+1))*A,
        "unnormalized_tfd": vacuum,
        "pair_ground_weights": [q**(2*n)/sum(q**(2*j) for j in range(dim)) for n in range(dim)],
    }


@lru_cache(maxsize=1)
def depth_control_audit():
    """JSON-compatible exact local certificates and their explicit limits."""
    native = native_q_bounds()
    qlo, qhi = native["lower"], native["upper"]
    cases = []
    for label, lo, hi in (("exp(-2*pi)", qlo, qhi), ("1/2", Fraction(1, 2), Fraction(1, 2))):
        for ratio in (Fraction(1, 1000), Fraction(1, 10**10), Fraction(1, 10**16)):
            for weighted in (False, True):
                cert = depth_certificate(lo, hi, ratio, weighted=weighted)
                cert["q"] = label
                cases.append(cert)
    # This is a declared benchmark in GeV, not a prediction of M_c or mu.
    benchmark_ratio = Fraction(1000, 2435*10**15)
    benchmark = {}
    for weighted in (False, True):
        label = "weighted" if weighted else "flat"
        benchmark[label] = {
            str(margin): depth_certificate(qlo, qhi, benchmark_ratio, weighted=weighted, margin=margin)
            for margin in (1, 10, 100)
        }
    q, Mc, mu = sp.symbols("q M_c mu", positive=True)
    d2 = edge_depth_data(2, q, Mc, mu, weighted=True)
    d3 = edge_depth_data(3, q, Mc, mu, weighted=True)
    finite = truncated_preparation_data(2, sp.Rational(2, 3))
    checks = {
        "native_q_lower_bound_proved": native["lower_proved"],
        "native_q_upper_bound_proved": native["upper_proved"],
        "edge_mass_is_stiffness_over_relative_mass_squared": sp.simplify(
            d2["isolated_edge_mass_squared"]-d2["edge_stiffness"]/d2["M_relative"]**2) == 0,
        "weighted_quartic_depth_step_is_three_over_two_q_squared": sp.simplify(
            d3["quartic_relative_field_coefficient_per_tau"]/d2["quartic_relative_field_coefficient_per_tau"]-3/(2*q**2)) == 0,
        "all_reported_counts_have_rational_certificates": all(c["status"] == "CERTIFIED_LOCAL_DEPTH" for c in cases)
            and all(c["status"] == "CERTIFIED_LOCAL_DEPTH" for case in benchmark.values() for c in case.values()),
        "finite_D2_pair_parent_has_physical_boundary": finite["paired_parent"] == finite["physical_radial_parent"],
        "finite_D2_parent_annihilates_truncated_TFD": finite["full_parent"]*finite["unnormalized_tfd"] == sp.zeros(9, 1),
        "finite_D2_CCR_top_defect_is_explicit": finite["single_mode_commutator"] == finite["expected_commutator"],
    }
    return {
        "status": "CONDITIONAL_LOCAL_EFT_DEPTH_CERTIFICATES",
        "physical_solution": False, "checks": checks, "all_checks_pass": all(checks.values()),
        "native_q_interval": [str(qlo), str(qhi)], "native_q_proof": native["proof"],
        "historical_ratio_cases": cases,
        "benchmark_inputs_GeV": {"M_c": "2435000000000000000", "mu": "1000"},
        "benchmark_margin_dependence": benchmark,
        "flat_infinite_spectrum": "zero plus a bounded massive band; infinitely many massive modes in that band",
        "weighted_infinite_spectrum": "symmetric ladder mu^2*(1-q^2)*n and three interlacing Dirichlet ladders; finite counting below fixed finite mass squared",
        "growth_per_interior_depth": {
            "inverse_local_Planck_mass": "q^(-1)",
            "quartic_transfer_row_flat": "q^(-2)",
            "quartic_transfer_row_weighted": "[(d+1)/d]*q^(-2)",
            "local_Lambda3_flat": "q^(1/3)",
            "local_Lambda3_weighted": "[(d+1)/d]^(1/3)*q^(1/3)",
        },
        "normalized_zero_mode_tail_beyond_D_layers": "q^(2*(D+1)); for D=5 this is q^12",
        "zero_mode_tail_scope": "Spatial norm tail only. It does not bound all loop effects or interactions; eta_v/M_v=1/M_eff remains universal.",
        "finite_preparation_parent": {
            "formula": "H_D=lambda/2*[(a_D-q*b_D^dagger)^dagger*(a_D-q*b_D^dagger)+(b_D-q*a_D^dagger)^dagger*(b_D-q*a_D^dagger)]",
            "single_mode_commutator": "[a_D,a_D^dagger]=I-(D+1)*|D><D|",
            "interior_pair_diagonal": "lambda*[n+q^2*(n+1)], 0<=n<D",
            "terminal_pair_diagonal": "lambda*D",
            "pair_upper_diagonal": "-lambda*q*(n+1), 0<=n<D",
            "unique_ground_state": "sum_(n=0)^D q^n*|n,n>, normalized",
            "uniqueness_proof": "Vectorize a common null vector as X=q^N Y. The two annihilation equations become [a_D,Y]=[a_D^dagger,Y]=0. These matrices generate the full matrix algebra, so Y is scalar.",
            "thermal_marginal": "rho_D=(1-q^2)/(1-q^(2*(D+1)))*sum_(n=0)^D q^(2*n)*|n><n|",
            "selected_depth": False,
            "scope": "Exact finite positive preparation parent for each chosen D; not a regular CCR representation or the infinite Gaussian uniqueness theorem",
        },
        "not_derived": ["physical end boundary", "collective multigravity cutoff", "quantum completion", "M_c", "mu", "edge quartic parameters"],
        "sources": ["https://arxiv.org/abs/1511.02877", "https://arxiv.org/abs/1806.10607",
                    "https://pdg.lbl.gov/2025/reviews/rpp2025-rev-phys-constants.pdf"],
    }
