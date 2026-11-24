"""Conditional current KMS selection and its positive Gaussian preparation parent.

For one regular oscillator, full KMS for alpha_t(a)=exp(-it)*a fixes
rho=(1-exp(-beta))*exp(-beta*N). Matrix units give p[n+1]=exp(-beta)*p[n];
Gaussianity is a consequence on this oscillator, not an extra hypothesis.
The canonical purification has q=exp(-beta/2). A mean occupation alone
does not imply this state. At the supplied beta=4*pi, q=exp(-2*pi).

On two canonical oscillators let dL=(a-q*b.dag)/sqrt(1-q*q), and similarly
dR. Their common vacuum is sqrt(1-q*q)*sum q**n |n,n>. A positive quadratic
parent preserving N_L-N_R and sheet exchange is uniquely Omega*sum d.dag*d,
up to Omega>0. Its restriction to the paired sector is
K=scale*A.dag*(N+1)*A, (A f)[n]=q*f[n]-f[n+1], Omega=scale*(1-q*q)/2.
Thus K has eigenvalues Delta*m, Delta=scale*(1-q*q), and its boundary
measure is (1-q*q)*sum q**(2*m)*delta_(Delta*m). This gives the complete
Stieltjes resolvent and its norm, not just a sampled positivity check.

The flat Hardy parent scale*A.dag*A has the same null vector but a bounded
continuous band. The two parents are not unitarily equivalent for 0<q<1.
J_-1 itself is a weighted oscillator creator, not a unilateral shift;
only its polar phase is an isometry on a chosen cyclic oscillator sector.

All physical source identifications, a common clock, and a four-dimensional
gravity completion remain separate inputs. This module does not replace
the complete current algebra or the E8 source by its single first mode.
"""

from fractions import Fraction
from functools import lru_cache

import sympy as sp


EXACT_Q_ORIGINAL = sp.exp(-2*sp.pi)


def _exact(value, name):
    if isinstance(value, (float, bool)):
        raise TypeError(f"{name} must be exact")
    if isinstance(value, Fraction):
        result = sp.Rational(value.numerator, value.denominator)
    elif isinstance(value, int):
        result = sp.Integer(value)
    elif isinstance(value, str):
        aliases = {"e^-2pi": EXACT_Q_ORIGINAL, "e^-pi": sp.exp(-sp.pi)}
        if value in aliases:
            result = aliases[value]
        else:
            fraction = Fraction(value)
            result = sp.Rational(fraction.numerator, fraction.denominator)
    elif isinstance(value, sp.Expr):
        result = value
    else:
        raise TypeError(f"unsupported exact {name}")
    if result.has(sp.Float) or result.free_symbols or result.is_real is not True:
        raise ValueError(f"{name} must be an exact specified real number")
    return result


def _parameters(q, scale=1):
    q, scale = _exact(q, "q"), _exact(scale, "scale")
    if q.is_nonnegative is not True or (1-q).is_positive is not True:
        raise ValueError("require 0 <= q < 1, proved without numerical sampling")
    if scale.is_positive is not True:
        raise ValueError("scale must be provably positive")
    return q, scale


def gaussian_parent_certificate(q, scale=1):
    """JSON-compatible exact spectral data for scale*A.dag*(N+1)*A.

    Here scale multiplies an operator interpreted as mass squared only
    under an additional physical dictionary. The two-oscillator full
    preparation gap is Delta/2; the paired-sector gap is Delta.
    """
    q, scale = _parameters(q, scale)
    r, Z = q*q, 1-q*q
    delta = scale*Z
    ell = -sp.log(Z)
    a = ell/delta
    b = (ell**2+sp.polylog(2, r))/delta**2
    c = (ell**3+2*ell*sp.polylog(2, r)+sp.polylog(3, r))/delta**3
    serial = lambda x: str(sp.simplify(x))
    return {
        "status": "CONDITIONAL_EXACT_GAUSSIAN_PREPARATION_PARENT",
        "physical_solution": False,
        "q": serial(q), "r": serial(r), "Z": serial(Z), "scale": serial(scale),
        "pair_spacing": serial(delta), "full_two_mode_gap": serial(delta/2),
        "jacobi": {
            "diagonal": "scale*((1+r)*n+r), n>=0",
            "upper_diagonal": "-scale*q*(n+1)",
            "zero_vector": "sqrt(1-r)*q^n",
            "factorization": "K=scale*A^dagger*(N+1)*A; (A f)_n=q*f_n-f_(n+1)",
        },
        "self_adjointness": "The Bogoliubov unitary conjugates the two-mode number operator; its N_L=N_R sector defines K. Also sum_n 1/(scale*q*(n+1)) diverges for q>0 (Jacobi Carleman criterion).",
        "spectrum": "Delta*m, m=0,1,2,..., simple in the paired sector",
        "spectral_type": "pure point, unbounded spectrum",
        "boundary_is_cyclic": bool(q != 0),
        "raw_spectral_measure": "(1-r)*sum_(m>=0) r^m*delta_(Delta*m)",
        "raw_total_weight": "1", "raw_massless_residue": serial(Z),
        "raw_total_massive_weight": serial(r),
        "normalized_massless_residue": "1",
        "normalized_total_weight": serial(1/Z),
        "normalized_total_massive_weight": serial(r/Z),
        "raw_resolvent": "C_raw(x)=(1-r)*sum_(m>=0) r^m/(x+Delta*m)=(1-r)/Delta*lerchphi(r,1,x/Delta)",
        "normalized_resolvent": "C(x)=C_raw(x)/(1-r); its inverse has unit infrared slope",
        "raw_tail_bound_after_M_terms": "0<=tail<=r^M/(x+Delta*M), x>0; terms m=0,...,M-1 retained",
        "raw_operator_resolvent_norm": "1/x for x>0, since K>=0 has a zero eigenvalue",
        "raw_boundary_bound": "(1-r)/x <= C_raw(x) <= 1/x",
        "normalized_inverse_coefficients": {"a": serial(a), "b": serial(b), "c": serial(c)},
        "coefficient_convention": "H=x-a*x^2+b*x^3-c*x^4+O(x^5)",
        "invariant_b_over_a2": serial(1+sp.polylog(2, r)/ell**2) if q != 0 else None,
        "infrared_slope": "1", "ultraviolet_slope": serial(Z),
        "raw_mass_moments": ["1", serial(scale*r), serial(scale**2*r*(1+r)),
                             serial(scale**3*r*(1+4*r+r*r))],
        "flat_hardy_comparison": {
            "same_zero_vector": True,
            "same_boundary_total_and_massless_weights": True,
            "same_first_two_mass_moments": True,
            "third_mass_moment_difference": serial(scale**3*r*r),
            "flat_hardy_spectrum": "{0} union [scale*(1-q)^2,scale*(1+q)^2] for q>0",
            "same_process_operator": False,
            "gaussian_a_strictly_smaller_at_equal_scale_for_q_positive": bool(q != 0),
            "a_difference_proof": "a_Hardy-a_Gaussian=(r+(1-r)*log(1-r))/(scale*(1-r)^2)>0 for 0<r<1",
        },
        "preparation_parent_unique_under_named_symmetries": {
            "hypotheses": ["exactly two canonical oscillators", "the specified pure canonical TFD vacuum",
                           "positive quadratic parent with that unique ground state",
                           "conservation of N_L-N_R", "exchange of the two sheets"],
            "conclusion": "H=Omega*(dL^dagger*dL+dR^dagger*dR), Omega>0, up to an additive constant",
            "proof": "Vacuum annihilation excludes linear and anomalous terms in the d variables. Relative U(1) makes the positive matrix in d^dagger*G*d diagonal; sheet exchange equates its diagonal entries.",
            "remaining_scale_free": True,
            "modular_formulation": "In the canonical KMS standard form, Delta_mod=exp[-beta*(N_L-N_R)]. Requiring [H,Delta_mod]=0 and J_mod*H*J_mod=H gives the same selection in the declared two-mode quadratic class.",
            "modular_conjugation_scope": "Tomita J_mod is antiunitary with J_mod^2=+1. It is not identified with a fermionic spinor or deck operator of square -1. Such an operator, if commuting with K on this same space, would exclude a simple zero eigenvalue by Kramers degeneracy.",
            "physical_time_scope": "Respecting the modular data is an extra condition on the preparation parent; it does not identify that positive parent with the original current time generator.",
        },
        "scope": "The paired sector has an exact positive preparation parent. It is not identified with the original time generator, the flat Hardy transfer, all four Hardy grades, or an interacting graviton.",
    }


def gaussian_resolvent_partial(q, x, terms=8, scale=1):
    """Exact lower approximation and rigorous raw tail bound for x>0.

    This is an evaluator with a proved remainder, not the source of the
    spectral-positivity certificate. Output expressions are SymPy objects.
    """
    if isinstance(terms, bool) or not isinstance(terms, int) or terms < 1:
        raise ValueError("terms must be a positive integer")
    q, scale = _parameters(q, scale)
    x = _exact(x, "x")
    if x.is_positive is not True:
        raise ValueError("x must be provably positive")
    r, Z = q*q, 1-q*q
    delta = scale*Z
    partial = Z*sum(r**m/(x+delta*m) for m in range(terms))
    remainder = r**terms/(x+delta*terms)
    return {"partial": sp.simplify(partial), "tail_bound": sp.simplify(remainder)}


def gaussian_eigenpolynomial(n, spectral_level, q):
    """Boundary-normalized Meixner polynomial F_n(m)=q^n M_n(m;1,q^2).

    The spectral transform has measure (1-q^2)*q^(2*m). This helper is
    for exact algebraic checks; completeness also follows from the
    unitary two-mode Bogoliubov diagonalization.
    """
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError("n must be a nonnegative integer")
    q, _ = _parameters(q)
    if q == 0:
        raise ValueError("the boundary is not cyclic at q=0")
    m = sp.sympify(spectral_level)
    return sp.expand_func(q**n*sum(sp.binomial(n, k)*sp.binomial(m, k)*(1-1/q**2)**k
                                  for k in range(n+1))).expand()


def current_kms_certificate(beta=4*sp.pi):
    """One oscillator: full KMS fixes the density; one occupation does not."""
    beta = _exact(beta, "beta")
    if beta.is_positive is not True:
        raise ValueError("beta must be provably positive")
    r, q = sp.exp(-beta), sp.exp(-beta/2)
    Z = 1-r
    epsilon = Z*r/4
    serial = lambda x: str(sp.simplify(x))
    return {
        "status": "CONDITIONAL_UNIQUE_REGULAR_SINGLE_OSCILLATOR_KMS_STATE",
        "beta_times_mode_gap": serial(beta), "current_level": "1",
        "hypotheses": ["regular one-mode CCR representation, [J1,J_-1]=1",
                       "specified dynamics alpha_t(J1)=exp(-it)*J1",
                       "full KMS condition at the supplied beta, not just its two-point consequence",
                       "the current mode is physically identified with the proposed ladder"],
        "density": "rho_beta=(1-exp(-beta))*exp(-beta*N)",
        "density_matrix_recurrence": "p_(n+1)=exp(-beta)*p_n",
        "proof": "In the regular finite-mode representation the state is normal on the oscillator factor. Stationarity makes it diagonal in N; bounded analytic matrix units |n><n+1| give the displayed recurrence by KMS. Trace one fixes p_0.",
        "gaussianity_of_this_mode_is_a_consequence": True,
        "occupation": serial(r/Z), "canonical_tfd_q": serial(q),
        "canonical_tfd": "sqrt(1-q^2)*sum_(n>=0) q^n*|n>_L*|n>_R",
        "canonical_purification_requires_fixed_conjugate_sheet_dictionary": True,
        "gaussian_phase_family_with_same_reduced_state": "sqrt(1-q^2)*sum q^n*exp(i*n*theta)*|n,n>; a fixed natural-cone/reflection identification selects theta=0",
        "fermi_current_telescoping": {
            "source": "c_r=1/(1+exp(beta*(r+alpha)))",
            "identity": "c_(r+1)*(1-c_r)=(c_r-c_(r+1))/(exp(beta)-1)",
            "sum_over_all_integers": "1/(exp(beta)-1), independent of alpha",
            "boundary_conditions": "c_r tends to 1 at negative infinity and 0 at positive infinity",
        },
        "same_mean_counterexample": {
            "density": "rho_alt=rho_beta+epsilon*(P_0-2*P_1+P_2)",
            "epsilon": serial(epsilon),
            "first_three_probabilities": [serial(Z+epsilon), serial(Z*r-2*epsilon), serial(Z*r*r+epsilon)],
            "other_probabilities_unchanged": True,
            "positive_normalized_faithful": True,
            "occupation_unchanged": True, "thermal_two_point_relation_unchanged": True,
            "full_KMS_defect_p1_minus_r_p0": serial(-(2+r)*epsilon),
            "fourth_normal_ordered_moment_excess": serial(2*epsilon),
            "full_KMS_satisfied": False,
        },
        "scope": "This establishes neither the full E8 source state nor P1-to-time soldering, a two-sheet operator intertwiner, or a unique unrestricted process Hamiltonian.",
    }


def _finite_car_witness():
    # Four fermion modes, with the two lowest occupied. The first two
    # current steps do not touch a missing boundary of the filled sea.
    modes, dimension = 4, 16
    annihilators = []
    for j in range(modes):
        a = sp.zeros(dimension)
        for bits in range(dimension):
            if bits & (1 << j):
                sign = (-1)**((bits & ((1 << j)-1)).bit_count())
                a[bits ^ (1 << j), bits] = sign
        annihilators.append(a)
    creators = [a.T for a in annihilators]
    J = sum((creators[j+1]*annihilators[j] for j in range(modes-1)), sp.zeros(dimension))
    vacuum = sp.zeros(dimension, 1)
    vacuum[3] = 1
    norms = [int(((J**n*vacuum).T*(J**n*vacuum))[0]) for n in range(3)]
    return {
        "CAR_commutator_shifts_each_fermion_creation_mode": all(
            J*creators[j]-creators[j]*J == creators[j+1] for j in range(modes-1)),
        "current_powers_norm_squared": norms,
        "unilateral_shift_powers_norm_squared": [1, 1, 1],
        "operator_identity_refuted_by_second_power": norms[2] != 1,
    }


@lru_cache(maxsize=1)
def gaussian_parent_audit():
    """Complete JSON report; root lane integration does not alter any source snapshot."""
    q, n = sp.symbols("q n", positive=True)
    r = q*q
    sigma = sp.diag(1, 1, -1, -1)
    T = sp.Matrix([[1, 0, 0, -q], [0, 1, -q, 0],
                   [0, -q, 1, 0], [-q, 0, 0, 1]])
    recurrence = ((1+r)*n+r)*q**n-q*(n+1)*q**(n+1)-q*n*q**(n-1)
    y, t = sp.symbols("y t", positive=True)
    c0, c1 = 1/(1+y), 1/(1+t*y)
    fermi_identity = sp.cancel(c1*(1-c0)-(c0-c1)/(t-1))
    witness = _finite_car_witness()
    checks = {
        "unnormalized_bogoliubov_CCR": sp.simplify(T*sigma*T.T-(1-r)*sigma) == sp.zeros(4),
        "weighted_parent_annihilates_geometric_state": sp.simplify(recurrence) == 0,
        "fermi_occupation_telescoping_identity": fermi_identity == 0,
        "CAR_commutator_shift": witness["CAR_commutator_shifts_each_fermion_creation_mode"],
        "current_is_not_a_unilateral_shift": witness["operator_identity_refuted_by_second_power"],
    }
    return {
        "status": "CONDITIONAL_KMS_SELECTION_WITH_DISTINCT_GAUSSIAN_DYNAMICS",
        "physical_solution": False, "checks": checks, "all_checks_pass": all(checks.values()),
        "current_operator": {
            "finite_CAR_witness": witness,
            "CCR": "[J1,J_-1]=1; J_-1|n>=sqrt(n+1)|n+1>",
            "polar_shift": "V=J_-1*(N+1)^(-1/2); V|n>=|n+1>",
            "conditional_Hardy_map": "U|n>=e_(4*n+j) gives U*V*U^dagger=S^4 on one specified grade j",
            "not_proved_by_relabelling": "The CAR derivation [J_-1,psi_r^dagger]=psi_(r+1)^dagger does not identify J_-1 with S^4 as an operator on state vectors.",
            "full_current_level_dimensions_0_through_5": [int(sp.partition(k)) for k in range(6)],
            "isolated_J_minus_one_ladder_dimensions": [1]*6,
        },
        "one_mode_geometric_state": {
            "pure_Gaussian_for_nonzero_q": False,
            "witness": "A pure single-mode Gaussian must solve (u*a+v*a^dagger-alpha)|psi>=0. For coefficients q^n, the n=0,1,2 equations force u*q*(sqrt(3)+sqrt(2)-3)=0, incompatible with |u|^2-|v|^2=1 and q>0.",
            "two_mode_TFD_is_Gaussian": True,
            "scope": "An isometry |n> -> |n,n> changes the canonical observable dictionary; equal coefficients alone do not transport the process operator.",
        },
        "single_current_KMS": current_kms_certificate(),
        "original_q_gaussian_parent": gaussian_parent_certificate(EXACT_Q_ORIGINAL),
        "free_gaussian_counterfamily_before_KMS_is_imposed": {
            "parameters": "every 0<q<1 and fixed Omega>0",
            "parent": "H_q=Omega/(1-q^2)*[(1+q^2)*(N_L+N_R)+2*q^2-2*q*(a*b+a^dagger*b^dagger)]",
            "same_full_gap": "Omega", "positive_quadratic_and_sheet_symmetric": True,
            "does_not_contradict_full_KMS_selection": True,
        },
        "sheet_exchange": "The canonical TFD is exchange symmetric for every q. Unequal bimetric Planck weights do not follow as spontaneous sheet-exchange breaking from small q alone.",
        "four_grade_extension_scope": "With depth weight n+1 on a four-arm metric tree, the symmetric sector is K above. The other three sectors are its Dirichlet compression to n>=1; they are not claimed to have the same arithmetic spectrum.",
        "sources": ["https://arxiv.org/abs/1212.6049", "https://arxiv.org/abs/1811.11528"],
        "source_snapshot_scope": "The supplied H6h26 source fixes the one-particle Fermi covariance at beta=4*pi. Extending it to a regular physical current KMS algebra and identifying its conjugate sheet are stated hypotheses, not derived by a finite mode relabelling.",
    }
