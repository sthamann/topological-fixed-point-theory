"""Conditional nonlinear metric parent for the complete quarter-transfer L.

A central metric and four identical arms form a tree, not an interaction
cycle.  At arm depth n=0,...,m-1 choose M_(j,n)=M_c*q**(n+1)/2 and give
each parent-child HR potential the Fierz-Pauli stiffness mu^2*M_child^2.
Use Einstein action M^2*R/2 and g=eta+2*h/M, so a minimally coupled source
has coupling h*T/M.  With V=(-4+e1+e3-4*e4)/2, the action -k*V gives
precisely the stated stiffness k in this convention.

An exact Fourier transform on the four arms sends the canonical mass
matrix to A_rect^* A_rect, where A_rect is the physical rectangular
truncation of S^*(1-q*S^4).  The symmetry, all four residue sectors, and
the single zero mode are matched; no sector is silently discarded.

EH terms, HR potentials, their weights, 4D spacetime, and source choices
are explicit extra inputs.  The finite result is a classical EFT
completion, not a proof of a UV-complete infinite quantum theory.
"""
import sympy as sp


def _positive_exact(value, name):
    if isinstance(value, (float, bool)):
        raise TypeError(f"{name} must be exact")
    value = sp.sympify(value)
    if value.has(sp.Float):
        raise TypeError(f"{name} must contain no floats")
    if value.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return value


def star_metric_data(arm_length, q, M_c=1, mu_squared=1):
    """Return exact matrices for 1+4*m metrics and 4*m interaction edges.

    Vertex ordering is center, then four vertices at each arm depth.
    Edge ordering is likewise four parent-child edges at each depth.
    U maps the native transfer column ordering (0,...,4*m) into vertices.
    V maps the native transfer row ordering (0,...,4*m-1) into edges.
    Thus B*U=V*A_rect, and U^* L_star U=mu_squared*A_rect^* A_rect.
    """
    if isinstance(arm_length, bool) or not isinstance(arm_length, int) or arm_length < 1:
        raise ValueError("arm_length must be a positive integer")
    q = _positive_exact(q, "q")
    M_c = _positive_exact(M_c, "M_c")
    mu_squared = _positive_exact(mu_squared, "mu_squared")
    m, edges, vertices = arm_length, 4*arm_length, 4*arm_length+1
    vertex = lambda n, j: 1+4*n+j
    plans = [M_c]+[M_c*q**(n+1)/2 for n in range(m) for _ in range(4)]
    D = sp.diag(*plans)
    incidence = sp.zeros(edges, vertices)
    B = sp.zeros(edges, vertices)
    weights = []
    for n in range(m):
        for j in range(4):
            row, child = 4*n+j, vertex(n, j)
            parent = 0 if n == 0 else vertex(n-1, j)
            incidence[row, parent], incidence[row, child] = 1, -1
            B[row, parent], B[row, child] = (q/2 if n == 0 else q), -1
            weights.append(mu_squared*plans[child]**2)
    K = sp.diag(*weights)
    canonical_mass = sp.simplify(D.inv()*incidence.T*K*incidence*D.inv())

    # A_rect is an actual map between distinct finite vertex/edge spaces.
    # It is not a square compression which is then assigned a fictitious index.
    native = sp.zeros(edges, vertices)
    for row in range(edges):
        native[row, row+1] = 1
        if row >= 3:
            native[row, row-3] = -q

    U = sp.zeros(vertices)
    U[0, 0] = 1
    V = sp.zeros(edges)
    for n in range(m):
        for r in range(4):
            input_col = 4*(n+1) if r == 0 else 4*n+r
            output_col = 4*n+(r-1) % 4
            for j in range(4):
                fourier = sp.I**(-r*j)/2
                U[vertex(n, j), input_col] = fourier
                V[4*n+j, output_col] = -fourier

    rotation = sp.zeros(vertices)
    rotation[0, 0] = 1
    for n in range(m):
        for j in range(4):
            rotation[vertex(n, (j+1) % 4), vertex(n, j)] = 1
    native_rotation = sp.diag(*(sp.I**n for n in range(vertices)))
    norm_squared = sp.expand(sum(q**(2*n) for n in range(m+1)))
    M_eff_squared = M_c**2*norm_squared
    zero = sp.Matrix(plans)/(M_c*sp.sqrt(norm_squared))
    source_couplings = sp.simplify(D.inv()*zero)
    native_zero = sp.Matrix([q**(n//4) if n % 4 == 0 else 0 for n in range(vertices)])/sp.sqrt(norm_squared)
    cubic_zero_weight = sp.simplify(sp.diag(*(zero[n]/plans[n] for n in range(vertices))))
    return {
        "arm_length": m, "metrics": vertices, "edges": edges,
        "q": q, "M_c": M_c, "mu_squared": mu_squared,
        "incidence": incidence, "edge_fp_weights": K,
        "planck_masses": D, "canonical_incidence": B,
        "mass_matrix": canonical_mass, "native_A_rect": native,
        "U_vertex_fourier": U, "V_edge_fourier": V,
        "rotation": rotation, "native_rotation": native_rotation,
        "zero": zero, "native_zero": native_zero,
        "M_eff_squared": sp.expand(M_eff_squared),
        "source_couplings": source_couplings,
        "physical_pole_matrix": sp.simplify(source_couplings*source_couplings.T),
        "boundary_zero_weight": sp.simplify(1/norm_squared),
        "cubic_zero_weight_matrix": cubic_zero_weight,
    }


def star_metric_audit():
    """Return a JSON-compatible analytic audit, with no floating-point gate."""
    q, M_c, mu2 = sp.symbols("q M_c mu2", positive=True)
    d = star_metric_data(3, q, M_c, mu2)
    nv, ne = d["metrics"], d["edges"]
    U, V, B, A = (d[name] for name in (
        "U_vertex_fourier", "V_edge_fourier", "canonical_incidence", "native_A_rect"
    ))
    zero_matrix = lambda expr, rows, cols: sp.simplify(expr) == sp.zeros(rows, cols)
    coupling = d["source_couplings"][0]
    checks = {
        "vertex_fourier_is_unitary": zero_matrix(U.conjugate().T*U-sp.eye(nv), nv, nv),
        "edge_fourier_is_unitary": zero_matrix(V.conjugate().T*V-sp.eye(ne), ne, ne),
        "full_first_order_intertwiner": zero_matrix(B*U-V*A, ne, nv),
        "full_mass_operator_intertwiner": zero_matrix(U.conjugate().T*d["mass_matrix"]*U-mu2*A.T*A, nv, nv),
        "canonical_normalization_from_EH_weights": zero_matrix(d["mass_matrix"]-mu2*B.T*B, nv, nv),
        "quarter_rotation_intertwiner": zero_matrix(d["rotation"]*U-U*d["native_rotation"], nv, nv),
        "native_surjectivity_unit_minor": A[:, 1:].det() == 1,
        "one_exact_zero_vector": zero_matrix(B*d["zero"], ne, 1),
        "zero_vector_is_normalized": sp.simplify((d["zero"].T*d["zero"])[0]-1) == 0,
        "native_zero_vector_intertwiner": zero_matrix(U*d["native_zero"]-d["zero"], nv, 1),
        "universal_physical_pole_matrix": zero_matrix(d["physical_pole_matrix"]-sp.ones(nv)/d["M_eff_squared"], nv, nv),
        "universal_zero_heavy_heavy_EH_overlap": zero_matrix(d["cubic_zero_weight_matrix"]-coupling*sp.eye(nv), nv, nv),
    }
    return {
        "status": "CONDITIONAL_FINITE_CLASSICAL_FULL_QUARTER_TRANSFER_COMPLETION",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "physical_solution": False,
        "graph": "one center and four identical arms; a tree with 1+4*m metrics and 4*m pairwise HR interactions",
        "action": "sum_v M_v^2/2 integral sqrt(-g_v) R[g_v] - sum_edges k_edge integral sqrt(-g_parent) V(sqrt(g_parent^-1*g_child)) + sum_v S_matter,v[g_v]",
        "potential": "V(X)=(-4+e1(X)+e3(X)-4*e4(X))/2",
        "canonical_field_convention": "g_v=eta+2*h_v/M_v; minimal source coupling is h_v*T_v/M_v",
        "parameter_assignment": {
            "center": "M_c",
            "arm_planck_mass": "M_(j,n)=(M_c/2)*q^(n+1), j=0,1,2,3, n=0,...,m-1",
            "edge_fp_weight": "k_edge=mu^2*M_child^2",
        },
        "first_order_operator": "B center edges=(q/2)*h_c-h_(j,0); B arm edges=q*h_(j,n-1)-h_(j,n)",
        "native_rectangular_operator": "A_rect[j,j+1]=1; A_rect[j,j-3]=-q for j>=3; shape (4*m,4*m+1)",
        "exact_intertwiner": "B*U=V*A_rect with U,V unitary four-arm Fourier transforms; U^dagger L_star U=mu^2*A_rect^dagger*A_rect",
        "symmetry": "The physical arm permutation is order four and U intertwines it with R*e_n=i^n*e_n.",
        "all_residue_sectors": "The invariant center-plus-symmetric-arm sector is L_q. The other three Fourier sectors tend to (1+q^2)I-q(T+T^dagger); none is discarded.",
        "finite_boundary_condition": "The last diagonal entry on every finite arm is 1, not 1+q^2. This is a physical rectangular truncation, not quadratic compression of the infinite operator.",
        "nullity": 1,
        "finite_index": "dim ker(A_rect)-dim coker(A_rect)=1; the row minor using columns 1,...,4*m is unit triangular",
        "M_eff_squared": "M_c^2*sum_(n=0)^m q^(2*n)",
        "universal_massless_coupling": "eta_v/M_v=1/M_eff at every metric; all physical source pole residues are 1/M_eff^2",
        "newton_constant_convention": "G_eff=1/(8*pi*M_eff^2), Einstein action M_eff^2*R/2",
        "infinite_limit": "For fixed 0<q<1, the bounded mass operators converge strongly to the complete quarter-transfer L; M_eff^2 -> M_c^2/(1-q^2).",
        "boundary_response": "The center resolvent is exactly the symmetric-arm response. Its infinite limit at momentum p is mu^(-2)*C_raw,q(p^2/mu^2); physical source exchange includes a further M_c^(-2). The three other sectors are present but have zero center-source overlap at quadratic order.",
        "connected_source_scope": "A centered 4D metric field retains the massless connected two-point pole. The unique-vacuum subtraction argument for a microscopic Hamiltonian is a different construction.",
        "tensor_scope": "C_q is the scalar transverse-traceless coefficient. For conserved 4D sources, massless and massive spin-2 exchange have trace factors 1/2 and 1/3 respectively; scalar covariance positivity alone does not prove nonlinear quantum gravity positivity.",
        "nonlinear_freedom": "V_t(X)=V(X)+(t/2)*det(I-X) preserves flatness, the quadratic mass operator and cubic vertices, while changing quartic interactions.",
        "loop_warning_identity": "For any orthonormal massive eigenvectors u_a,u_b: sum_v eta_v*conj(u_a(v))*u_b(v)/M_v=delta_ab/M_eff. Thus massive species do not decouple from the massless metric as their number grows.",
        "loop_warning_scope": "At fixed mu,q an unbounded number of massive modes occupy a bounded mass band. This motivates an explicit tensor and loop calculation; it is not a proved universal UV no-go.",
        "boundary_tree_vertex_bound": "For x>0 in the symmetric radial channel, C_0n(x)=t(x)^n*C_00(x), 0<t(x)<q, with t(0)=q as a limit. A physical individual arm has the additional Fourier factor 1/2. A fixed single k-point EH vertex with boundary legs has a depth ratio at most q^2 and its depth sum converges.",
        "first_new_inputs": [
            "4D dynamical metrics and Einstein kinetic terms",
            "tree of pairwise Hassan-Rosen potentials",
            "geometric Planck weights and edge stiffnesses",
            "minimal coupling of each matter sector to exactly one metric",
        ],
        "unfixed": ["q", "M_c", "mu", "arm length", "higher potential coefficients", "TFPT physical source and field dictionary"],
        "quantum_or_infinite_interacting_completion_proved": False,
        "sources": [
            "https://arxiv.org/abs/1109.3515",
            "https://arxiv.org/abs/1203.5783",
            "https://arxiv.org/abs/1207.3637",
            "https://arxiv.org/abs/1410.7774",
        ],
    }
