"""Exact tree-level quartic response of physical metric sources.

Conventions: g_v=eta+2*h_v/M_v, source sum_v h_v:J_v/M_v, and
V_tau=V_0+(tau_e/2)*det(I-sqrt(g_parent^-1*g_child)).
The output is the change of the connected Euclidean four-source vertex,
with W_E=log integral exp(-S_E+source), at tree level.  A Lorentzian
contact amplitude instead uses the positive quartic Lagrangian below
and the Feynman i factors; a Euclidean W derivative is not an S-matrix.

Each external source must be transverse and traceless for its own momentum.
Then its tensor is unchanged by the TT propagator, and the scalar profiles
below suffice.  There is no extra 4! in fourpoint_prefactor: multiply it
by determinant_vertex(T1,T2,T3,T4), defined as the coefficient of z1*z2*z3*z4.

The helpers use exact algebra.  Their finite graphs are classical EFTs,
and a convergent single infinite tree contact is not a quantum completion.
"""
from fractions import Fraction
from itertools import permutations

import sympy as sp


def _exact(value, name):
    if isinstance(value, (bool, float)):
        raise TypeError(f"{name} must be exact, not a float or boolean")
    if isinstance(value, Fraction):
        value = sp.Rational(value.numerator, value.denominator)
    else:
        value = sp.sympify(value)
    if not isinstance(value, sp.Expr) or value.has(sp.Float):
        raise TypeError(f"{name} must be an exact scalar expression")
    return value


def _positive(value, name, allow_zero=False):
    value = _exact(value, name)
    valid = value.is_nonnegative if allow_zero else value.is_positive
    if valid is not True:
        raise ValueError(f"{name} must be provably {'nonnegative' if allow_zero else 'positive'}")
    return value


def _unit_q(q):
    q = _positive(q, "q")
    if sp.simplify(1-q).is_positive is not True:
        raise ValueError("the infinite-star certificate requires 0<q<1")
    return q


def tree_source_data(planck_masses, edges, stiffnesses):
    """An oriented rooted tree, root 0, with positive EH masses and FP weights.

    Every edge is (parent,child).  The physical source operator is phi_v=h_v/M_v.
    Its TT covariance is (x*diag(M_v^2)+incidence.T*K*incidence)^(-1).
    """
    masses = tuple(_positive(m, "Planck mass") for m in planck_masses)
    edges = tuple(tuple(edge) for edge in edges)
    stiff = tuple(_positive(k, "edge stiffness") for k in stiffnesses)
    nv, ne = len(masses), len(edges)
    if nv < 2 or ne != nv-1 or len(stiff) != ne:
        raise ValueError("a nontrivial tree requires V-1 edges and one stiffness per edge")
    parents, children = {}, [[] for _ in range(nv)]
    incidence = sp.zeros(ne, nv)
    for ei, edge in enumerate(edges):
        if len(edge) != 2 or any(isinstance(v, bool) or not isinstance(v, int) for v in edge):
            raise TypeError("edge endpoints must be integer vertex indices")
        parent, child = edge
        if not (0 <= parent < nv and 0 < child < nv) or parent == child or child in parents:
            raise ValueError("edges must give each nonroot vertex one parent")
        parents[child] = parent
        children[parent].append(child)
        incidence[ei, parent], incidence[ei, child] = 1, -1
    visited = set()
    def descendants(vertex):
        if vertex in visited:
            raise ValueError("the graph is not a rooted tree")
        visited.add(vertex)
        result = {vertex}
        for child in children[vertex]:
            result |= descendants(child)
        return result
    if descendants(0) != set(range(nv)):
        raise ValueError("the rooted graph must be connected")
    def subtree(vertex):
        return {vertex}.union(*(subtree(c) for c in children[vertex]))
    subtrees = tuple(tuple(sorted(subtree(c))) for _, c in edges)
    D, K = sp.diag(*masses), sp.diag(*stiff)
    physical_laplacian = incidence.T*K*incidence
    return {
        "vertices": nv, "edge_count": ne, "edges": edges,
        "planck_masses": masses, "stiffnesses": stiff,
        "incidence": incidence, "physical_laplacian": physical_laplacian,
        "canonical_mass": D.inv()*physical_laplacian*D.inv(),
        "mass_squared_weights": sp.Matrix([m*m for m in masses]),
        "M_eff_squared": sp.expand(sum(m*m for m in masses)),
        "subtrees": subtrees,
    }


def star_boundary_data(arm_length, q, M_c=1, mu_squared=1, depth_weights=None):
    """The finite physical four-arm star, including every arm and end edge."""
    if isinstance(arm_length, bool) or not isinstance(arm_length, int) or arm_length < 1:
        raise ValueError("arm_length must be a positive integer")
    q, M_c, mu_squared = (_positive(q, "q"), _positive(M_c, "M_c"),
                          _positive(mu_squared, "mu_squared"))
    weights = (sp.Integer(1),)*arm_length if depth_weights is None else tuple(
        _positive(w, "depth weight") for w in depth_weights)
    if len(weights) != arm_length:
        raise ValueError("one positive weight is required for every depth")
    masses = [M_c]+[M_c*q**(n+1)/2 for n in range(arm_length) for _ in range(4)]
    edges, stiffnesses = [], []
    for n in range(arm_length):
        for arm in range(4):
            child = 1+4*n+arm
            parent = 0 if n == 0 else 1+4*(n-1)+arm
            edges.append((parent, child))
            stiffnesses.append(mu_squared*weights[n]*masses[child]**2)
    data = tree_source_data(masses, edges, stiffnesses)
    data.update({"arm_length": arm_length, "q": q, "M_c": M_c,
                 "mu_squared": mu_squared, "depth_weights": weights})
    return data


def _source(data, source):
    nv = data["vertices"]
    if isinstance(source, int) and not isinstance(source, bool):
        if not 0 <= source < nv:
            raise ValueError("source vertex outside the graph")
        vector = sp.zeros(nv, 1)
        vector[source] = 1
        return vector
    values = tuple(source)
    if len(values) != nv:
        raise ValueError("a distributed physical source needs one coefficient per vertex")
    return sp.Matrix([_exact(v, "source coefficient") for v in values])


def edge_source_response(data, momentum_squared, source=0):
    """Return D_e/J=(h_child/M_child-h_parent/M_parent)/J, all edges.

    A vertex source already includes its physical canonical coupling 1/M_v.
    At zero the pole-free analytic limit is returned, never an inverse of
    the singular zero-momentum mass operator.
    """
    x = _positive(momentum_squared, "momentum_squared", allow_zero=True)
    source = _source(data, source)
    if x == 0:
        return edge_source_response_zero(data, source)
    kinetic = sp.diag(*data["mass_squared_weights"])
    phi = (x*kinetic+data["physical_laplacian"]).inv()*source
    return (-data["incidence"]*phi).applyfunc(sp.factor)


def edge_source_response_zero(data, source=0):
    """Exact tree flux identity after subtracting the universal massless pole.

    For the subtree S_e below an edge,
    k_e*D_e(0)=sum_(v in S_e) J_v - M_S_e^2*sum_v J_v/M_eff^2.
    """
    source = _source(data, source)
    total_source = sum(source)
    m2 = data["mass_squared_weights"]
    total_mass = data["M_eff_squared"]
    return sp.Matrix([
        sp.factor((sum(source[v] for v in sub)-sum(m2[v] for v in sub)*total_source/total_mass)/k)
        for sub, k in zip(data["subtrees"], data["stiffnesses"])
    ])


def physical_two_point(data, momentum_squared, vertex=0):
    x = _positive(momentum_squared, "momentum_squared")
    vector = _source(data, vertex)
    kinetic = sp.diag(*data["mass_squared_weights"])
    return sp.factor((vector.T*(x*kinetic+data["physical_laplacian"]).inv()*vector)[0])


def determinant_vertex(*tensors):
    """The coefficient of z1*z2*z3*z4 in det(sum z_i*T_i).

    This is symmetric and multilinear, equals 24*det(T) for four equal
    arguments, and has no assumptions about the chosen polarizations.
    The caller supplies physical TT tensors for a TT response.
    """
    if len(tensors) != 4:
        raise ValueError("four tensor arguments are required")
    tensors = tuple(sp.Matrix(t) for t in tensors)
    if any(t.shape != (4, 4) or any(v.has(sp.Float) for v in t) for t in tensors):
        raise ValueError("each tensor must be an exact 4 by 4 matrix")
    return sp.expand(sum(
        sp.Matrix.hstack(*(tensors[order[col]][:, col] for col in range(4))).det()
        for order in permutations(range(4))
    ))


def _edge_taus(data, taus):
    if isinstance(taus, (tuple, list, sp.MatrixBase)):
        values = tuple(_exact(t, "tau") for t in taus)
        if len(values) != data["edge_count"]:
            raise ValueError("one tau coefficient is required per edge")
        return values
    return (_exact(taus, "tau"),)*data["edge_count"]


def fourpoint_prefactor(data, momenta_squared, physical_sources=None, taus=1):
    """Change in the connected Euclidean tree four-source coefficient.

    Multiply by determinant_vertex(T1,T2,T3,T4) and the momentum-conserving
    delta function.  Each external physical source is a vertex index or
    a vector of source weights.  Default: four central sources.
    Distinct tau_e are allowed; a common tau is an additional input.
    """
    xs = tuple(momenta_squared)
    sources = (0, 0, 0, 0) if physical_sources is None else tuple(physical_sources)
    if len(xs) != 4 or len(sources) != 4:
        raise ValueError("four momenta and four physical source profiles are required")
    profiles = [edge_source_response(data, x, source) for x, source in zip(xs, sources)]
    taus = _edge_taus(data, taus)
    return sp.factor(-sum(k*t*sp.prod(prof[e] for prof in profiles)
                         for e, (k, t) in enumerate(zip(data["stiffnesses"], taus)))/2)


def flat_infinite_contact_data(q, momenta_squared, M_c=1, mu_squared=1, tau=1):
    """Closed single-contact sum for the infinite flat-weight star.

    A common tau at every depth is assumed here explicitly.  The returned
    generating-function variable also supports arbitrary bounded tau_n.
    Zero momentum means the continuous edge-profile limit, not C_raw(0).
    """
    q, M_c, mu_squared = (_unit_q(q), _positive(M_c, "M_c"),
                          _positive(mu_squared, "mu_squared"))
    tau = _exact(tau, "tau")
    ps = tuple(_positive(p, "momentum_squared", allow_zero=True) for p in momenta_squared)
    if len(ps) != 4:
        raise ValueError("four squared momenta are required")
    r = q*q
    z = tuple(sp.Integer(1) if p == 0 else
              2/(p/mu_squared+1+r+sp.sqrt((p/mu_squared+1+r)**2-4*r)) for p in ps)
    P = sp.prod(z)
    prefactor = -tau*r*P/(2*mu_squared**3*M_c**6*(1-r*P))
    return {
        "q": q, "r": r, "z": z, "P": P,
        "euclidean_prefactor": prefactor,
        "infrared_prefactor": -tau*r/(2*mu_squared**3*M_c**6*(1-r)),
        "per_depth_edge_profile": "-z(p^2/mu_squared)^(n+1)/(mu_squared*M_c^2)",
        "general_depth_rule": "F(P)=-sum_(n>=0) tau_n*(r*P)^(n+1)/(2*mu_squared^3*M_c^6)",
        "common_tau_is_additional": True,
    }


def flat_depth_contact_prefactor(q, momenta_squared, taus_by_depth, M_c=1, mu_squared=1):
    """Exact infinite-flat-parent contact when tau_n has finite support.

    This truncates the interaction coefficients, not the free mass operator.
    A finite physical star instead uses fourpoint_prefactor on its own graph.
    """
    data = flat_infinite_contact_data(q, momenta_squared, M_c, mu_squared, tau=0)
    M_c, mu_squared = _positive(M_c, "M_c"), _positive(mu_squared, "mu_squared")
    terms = tuple(_exact(t, "tau_n") for t in taus_by_depth)
    return sp.factor(-sum(t*(data["r"]*data["P"])**(n+1) for n, t in enumerate(terms))
                     /(2*mu_squared**3*M_c**6))


def _json(value):
    if isinstance(value, sp.Basic):
        return str(value)
    if isinstance(value, dict):
        return {k: _json(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(v) for v in value]
    return value


def boundary_fourpoint_audit():
    """Portable exact identities and scope, directly JSON serializable."""
    x, Mg, Mf, k, tau = sp.symbols("x M_g M_f k tau", positive=True)
    edge = tree_source_data((Mg, Mf), ((0, 1),), (k,))
    mass = k*(1/Mg**2+1/Mf**2)
    left = edge_source_response(edge, x, 0)[0]
    right = edge_source_response(edge, x, 1)[0]
    star = star_boundary_data(2, sp.Rational(1, 2), M_c=3, mu_squared=5)
    positive_x = sp.Integer(2)
    null_source = tuple(star["mass_squared_weights"])
    TA, TB = sp.diag(0, 1, 1, -2), sp.diag(1, 0, 1, -2)
    pA, pB = sp.Matrix([1, 0, 0, 0]), sp.Matrix([0, 1, 0, 0])
    response = edge_source_response(star, x, 0)
    soft = edge_source_response_zero(star, 0)
    flat_zero = flat_infinite_contact_data(sp.Rational(1, 2), (0,)*4, M_c=3, mu_squared=5, tau=tau)
    # R=z/(1-z), x=(1-z)*(1-r*z)/z proves the root and boundary equations.
    r, z = sp.symbols("r z", positive=True)
    xz = (1-z)*(1-r*z)/z
    raw = z/(1-z)
    checks = {
        "single_edge_parent_physical_source": sp.simplify(left+1/(Mg**2*(x+mass))) == 0,
        "single_edge_child_physical_source": sp.simplify(right-1/(Mf**2*(x+mass))) == 0,
        "finite_star_soft_flux_formula": all(sp.limit(response[e], x, 0) == soft[e] for e in range(len(soft))),
        "one_massless_distributed_source_kills_every_edge": edge_source_response(star, positive_x, null_source) == sp.zeros(star["edge_count"], 1),
        "one_massless_external_leg_kills_fourpoint": fourpoint_prefactor(
            star, (1, 2, 3, 4), (null_source, 0, 0, 0)) == 0,
        "explicit_sources_are_transverse": TA*pA == sp.zeros(4, 1) and TB*pB == sp.zeros(4, 1),
        "explicit_sources_are_traceless": sp.trace(TA) == sp.trace(TB) == 0,
        "noncollinear_TT_vertex_is_nonzero": determinant_vertex(TA, TB, TA, TB) == -16,
        "collinear_TT_vertex_vanishes": determinant_vertex(TA, TA, TA, TA) == 0,
        "determinant_vertex_includes_four_factorial": determinant_vertex(sp.eye(4), sp.eye(4), sp.eye(4), sp.eye(4)) == 24,
        "flat_root_recurrence_identity": sp.factor(r*z*z-(xz+1+r)*z+1) == 0,
        "flat_boundary_green_identity": sp.factor((xz+r-r*z)*raw-1) == 0,
        "flat_edge_cancels_massless_green_pole": sp.factor(raw*(z-1)+z) == 0,
        "flat_infinite_IR_prefactor": sp.simplify(flat_zero["euclidean_prefactor"]+tau/546750) == 0,
    }
    return _json({
        "status": "CONNECTED_BOUNDARY_FOURPOINT_AT_TREE_ORDER",
        "base_claims": "H6h27 main b5cc3d407; z50 metric star and z51 edge quartic",
        "all_checks_pass": all(checks.values()), "checks": checks,
        "physical_solution": False,
        "Euclidean_convention": "W_E=log integral exp(-S_E+sum_v h_v:J_v/M_v); Delta S_E^4=+sum_e k_e*tau_e*det(D_e)/2",
        "Lorentzian_convention": "signature(-+++); Delta L_L^4=+sum_e k_e*tau_e*det(D_e,cov)/2; tensor Wick rotation, propagator i factors and LSZ are separate",
        "canonical_source": "h_v:J_v/M_v, never h_v:J_v with the factor omitted",
        "edge_profile": "R_eb(x)=G_child,b(x)/(M_child*M_b)-G_parent,b(x)/(M_parent*M_b)",
        "connected_fourpoint_change": "-1/2 sum_e k_e*tau_e product_i R_e,b_i(p_i^2) V_det(T1,T2,T3,T4)",
        "determinant_vertex": "[z1*z2*z3*z4] det(sum_i z_i*T_i); equals 24*det(T) for equal tensors",
        "single_edge": {"mass_squared": mass, "left_source_profile": left, "right_source_profile": right},
        "soft_tree_identity": "R_eb(0)=[1_(b in subtree_e)-M_subtree,e^2/M_eff^2]/k_e",
        "soft_orders": {
            "unamputated_generic_TT_sources": "O(epsilon^0)+O(epsilon^2), p_i=epsilon*k_i",
            "four_inverse_boundary_propagators": "O(product_i p_i^2)=O(epsilon^8)",
            "pure_massless_eigenstates": "identically zero already if any one external eigenmode is massless",
        },
        "infinite_flat_star": {
            "z": "2/(x+1+r+sqrt((x+1+r)^2-4*r)), x=p^2/mu_squared",
            "edge_response": "-z^(n+1)/(mu_squared*M_c^2)",
            "common_tau": "-tau*r*P/[2*mu_squared^3*M_c^6*(1-r*P)], P=product_i z_i",
            "depth_dependent": "-sum_n tau_n*(r*P)^(n+1)/(2*mu_squared^3*M_c^6)",
            "assumptions": "0<q<1, nonnegative Euclidean momenta, bounded tau_n or an explicitly convergent series",
            "IR_common_tau": "-tau*r/[2*mu_squared^3*M_c^6*(1-r)]",
            "IR_linear_depth_weight": "-tau*polylog(3,r)/(2*mu_squared^3*M_c^6)",
        },
        "coupling_selection": {
            "tau_predicted": False,
            "quartic_matching_required": True,
            "mu4": "equal tau among the four arms at each fixed depth; it does not equate different depths",
            "central_source_blindness": "without mu4 symmetry, central sources see the mean tau of the four arms at each depth",
            "full_kernel_tomography": "the flat-star central fourpoint kernel is a generating function of tau_n; its analytic coefficients determine tau_n if the microscopic kernel and all lower normalizations are independently supplied",
            "IR_only": "one weighted sum of the tau_n, not each individual coupling",
            "Gaussian_state_is_not_a_zero_gravity_cumulant": (
                "A Gaussian state of one current mode does not determine the gravitational connected fourpoint function. "
                "A fully Gaussian process with a linear source map has zero higher cumulants, but then already conflicts "
                "with the nonzero Einstein threepoint and interactions at finite G. It does not merely select tau=0."),
            "match_the_difference_not_total_W4": (
                "The formulas compute W4(tau)-W4(0). The full tau-independent contacts and exchange diagrams must be "
                "included before matching an independently supplied microscopic tensor and momentum kernel."),
        },
        "scope": [
            "finite tree, common flat background, positive EH and FP coefficients",
            "conserved traceless external sources and the stated physical metric dictionary",
            "tree order only: exchange diagrams are unchanged because quadratic and cubic vertices are tau independent",
            "Euclidean connected response is not automatically a Lorentzian on-shell amplitude",
            "analytic low-momentum contact is not a new long-range massless force",
            "local counterterms and source conventions must be fixed in a microscopic matching",
            "convergence of this one infinite tree contact does not bound quantum corrections in the strongly coupled tail",
        ],
    })
