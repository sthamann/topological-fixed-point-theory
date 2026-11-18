"""Finite graded sewing, with distinct static and time observations."""
from itertools import combinations, product
from math import comb

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

from coherence import car_letter
from reconstruction import PAULI, require, scalar
from routes import partial_trace

I = sp.I
SIGNS = tuple(product((1,-1),repeat=4))
HADAMARD = sp.Matrix([[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]])


def complete_graphs(maximum=7):
    require(type(maximum) is int and 3<=maximum<=12,"choose three to twelve marks")
    rows = []
    for n in range(3,maximum+1):
        edges = tuple(combinations(range(n),2))
        cycles = sp.zeros(len(edges),(n-1)*(n-2)//2)
        for column,(a,b) in enumerate(combinations(range(1,n),2)):
            for edge,sign in (((0,a),1),((a,b),1),((0,b),-1)):
                cycles[edges.index(edge),column] = sign
        gram = cycles.T*cycles
        smith = smith_normal_form(gram,domain=sp.ZZ)
        factors = [int(abs(smith[k,k])) for k in range(gram.rows) if abs(smith[k,k])>1]
        require(gram.det()==n**(n-2) and factors==[n]*(n-2),"complete graph lattice mismatch")
        rows.append({"marks":n,"loops":gram.rows,"mark_relations":n-1,
                     "determinant":int(gram.det()),"critical_group_factors":factors})
    relation_trace = 2-1
    oriented_edge_trace = 1-1
    return {"rows":rows,"positive_rank_selection":[x["marks"] for x in rows if x["loops"]==x["mark_relations"]],
            "proof":"(n-1)(n-2)/2=n-1 for complete K_n with n>1 implies n=4",
            "rank_equality_also_on_self_dual_noncomplete_wheel":{"marks":5,"edges":8,"loops":4},
            "k4_transposition_traces":{"vertex_relations":relation_trace,"cycles":oriented_edge_trace-relation_trace},
            "equivariance_requires_orientation_twist":True,
            "k4_critical_group":[4,4],"a3_root_discriminant_group":[4],
            "scope":"complete graph rank selection and integral cycle lattices, no P1 completeness proof"}


def creation(modes,mode):
    result = sp.zeros(2**modes)
    for state in range(result.cols):
        target,sign = car_letter(state,mode,True)
        if sign:
            result[target,state] = sign
    return result


def clifford_cells(parity_signs=(1,1)):
    require(len(parity_signs)==2 and all(type(s) is int and s in (-1,1) for s in parity_signs),
            "declare one local parity sign per cell")
    # Independent of the occupation implementation in creation().
    one,x,y,z = PAULI
    tensor = sp.kronecker_product
    gamma = [tensor(one,x),tensor(one,y),tensor(x,z),tensor(y,z)]
    parity = -gamma[0]*gamma[1]*gamma[2]*gamma[3]
    a = [tensor(sp.eye(4),g) for g in gamma]
    b = [tensor(g,parity) for g in gamma]
    local_a = parity_signs[0]*tensor(sp.eye(4),parity)
    local_b = parity_signs[1]*tensor(parity,sp.eye(4))
    return a,b,local_a,local_b


def joint_projector(links,signs):
    eye = sp.eye(links[0].rows)
    return sp.prod((eye+s*m)/2 for s,m in zip(signs,links)).applyfunc(sp.expand)


def ground_spectrum(links,weights,selector):
    # Trace of each joint projector is its exact multiplicity.
    counts = {}
    for signs in SIGNS:
        multiplicity = (joint_projector(links,signs)*selector).trace()
        if multiplicity:
            energy = sum(abs(w)-w*s for w,s in zip(weights,signs))
            counts[energy] = counts.get(energy,0)+multiplicity
    require(sum(counts.values())==selector.trace(),"joint eigenspaces do not resolve the selector")
    energies = sorted(counts)
    return {"energy_multiplicities":[[str(e),int(counts[e])] for e in energies],
            "ground_multiplicity":int(counts[energies[0]]),
            "gap":str(energies[1]-energies[0]) if len(energies)>1 else None}


def temporal_rate(hamiltonian,ground,operator,projector):
    vector = (projector*operator*ground).applyfunc(sp.expand)
    norm = sp.expand((vector.H*vector).trace())
    if not norm:
        return None
    value = sp.expand((vector.H*hamiltonian*vector).trace())/norm
    require(not any((hamiltonian*vector-value*vector).applyfunc(sp.expand)),"temporal probe does not isolate one rate")
    return str(value)


def two_cells(data):
    g,mu,delta = (scalar(data.get(k,d)) for k,d in (("strength",1),("reference_weight",1),("gap",8)))
    require(all(x.is_Rational for x in (g,mu,delta)) and g>0 and mu>=0 and delta>0,
            "positive rational strength and gap and a nonnegative reference weight required")
    parity_signs = data.get("parity_signs",[1,1])
    a,b,pa,pb = clifford_cells(parity_signs)
    eye = sp.eye(16)
    require(all(op*op==eye for op in a+b) and
            all(x*y+y*x==sp.zeros(16) for x,y in combinations(a+b,2)),"crossed CAR failed")
    occupation = [creation(4,i) for i in range(4)]
    occupation_majoranas = [g for c in occupation for g in (c+c.H,I*(c-c.H))]
    require(a+b==occupation_majoranas,"independent occupation CAR does not match the tensor convention")
    links = [I*x*y for x,y in zip(a,b)]
    require(sp.prod(links)==parity_signs[0]*parity_signs[1]*pa*pb,"local parity orientation mismatch")
    rho = joint_projector(links,(1,)*4)
    require(rho.trace()==rho.rank()==1 and partial_trace(rho,[4,4],[0])==sp.eye(4)/4,
            "ground state is not the predicted entangled pure state")

    p = (eye+pa)*(eye+pb)/4
    v = g*sum(links,sp.zeros(16))
    q = eye-p
    early = p*v*p
    second = p*v*q*v*p
    require(early==sp.zeros(16) and second==p*v*v*p,"projection erased the ordered return")
    kk = sum(((I*a[i]*a[j])*(I*b[i]*b[j]) for i,j in combinations(range(4),2)),sp.zeros(16))
    require(second==g**2*p*(4*eye-2*kk)*p,"even return identity failed")
    basis = sp.Matrix.hstack(*p.columnspace())
    restricted = (basis.H*basis).inv()*basis.H*second*basis
    effective = [[str(-e/delta),int(mult)] for e,mult in sorted(restricted.eigenvals().items(),reverse=True)]
    energy_sizes = sorted({sp.sqrt(e) for e in restricted.eigenvals()})
    exact_blocks = [{"coupling_magnitude":str(e),
                     "low_energy":str((delta-sp.sqrt(delta**2+4*e**2))/2),
                     "high_energy":str((delta+sp.sqrt(delta**2+4*e**2))/2)} for e in energy_sizes]

    frame = HADAMARD/2
    ar = [sum((frame[i,j]*a[j] for j in range(4)),sp.zeros(16)) for i in range(4)]
    br = [sum((frame[i,j]*b[j] for j in range(4)),sp.zeros(16)) for i in range(4)]
    rotated = [(I*x*y).applyfunc(sp.expand) for x,y in zip(ar,br)]
    require(joint_projector(rotated,(1,)*4)==rho,"the weighted graph changed the positive ground ray")
    weights = [mu,4,4,4]
    h = sum((w*(eye-m) for w,m in zip(weights,rotated)),sp.zeros(16))
    require(h*rho==sp.zeros(16),"weighted reference spectrum is inconsistent")

    unit,x,_,z = PAULI
    tensor = sp.kronecker_product
    eye_link = sp.eye(32)
    gauge_a,gauge_b = tensor(pa,x),tensor(pb,x)
    gauge = (eye_link+gauge_a)*(eye_link+gauge_b)/4
    dressed = [tensor(m,z) for m in rotated]
    require(gauge.trace()==8 and all(m*ga==ga*m for m,ga in product(dressed,(gauge_a,gauge_b))),
            "dressed hopping failed the local gauge conditions")
    ground = joint_projector(dressed,(1,)*4)*gauge
    gauge_h = sum((w*(eye_link-m) for w,m in zip(weights,dressed)),sp.zeros(32))
    probes = {}
    if ground.trace():
        for label,op in (("single_reference",ar[0]),("reference_graph_pair",I*ar[0]*ar[1]),
                         ("graph_graph_pair",I*ar[1]*ar[2])):
            probes[label] = temporal_rate(gauge_h,ground,tensor(op,unit),gauge)
    expected_parity = parity_signs[0]*parity_signs[1]
    require(sp.prod(dressed)*gauge==expected_parity*gauge,"physical product is in the wrong sector")
    for signs in SIGNS:
        energy = -g*sum(signs)
        state = joint_projector(links,signs)
        projected = p*state
        if not any(projected):
            continue
        penalty = delta*(2*eye-pa-pb)/4
        exact = penalty-v
        low_polynomial = exact*exact-delta*exact-energy**2*eye
        require(low_polynomial*projected==sp.zeros(16) and low_polynomial*v*projected==sp.zeros(16),
                "exact return energy polynomial failed")
    return {"status":"EXACT_FINITE_SEWING","convention":{"crossed_car":True,"local_parity_signs":parity_signs},
            "ungauged_equal":ground_spectrum(links,[g]*4,eye),
            "ungauged_graph":ground_spectrum(rotated,weights,eye),
            "reduced_density":"I_4/4","entropy":"log(4)","ground_ray_for_all_positive_reference_weights":True,
            "unprojected_time_rates":{"reference":temporal_rate(h,rho,ar[0],eye),
                                      "graph":temporal_rate(h,rho,ar[1],eye),"ratio":str(mu/4)},
            "local_projection":{"rank":int(p.trace()),"direct_zero":not any(early),"ordered_return_rank":second.rank(),
                                "second_order_spectrum":effective,
                                "exact_coupled_blocks":exact_blocks,
                                "gap_convention":"H0 is zero on P and Delta on the reached complement",
                                "low_energy":"(Delta-sqrt(Delta^2+4*e^2))/2",
                                "energy_remainder_bound":"0 <= E_exact + e^2/Delta <= e^4/Delta^3"},
            "binary_link":{"physical_dimension":int(gauge.trace()),"dressed_product":expected_parity,
                           "equal":ground_spectrum(dressed,[g]*4,gauge),
                           "graph":ground_spectrum(dressed,weights,gauge),"projected_time_rates":probes,
                           "diagonal_link_observations_have_no_decay":True},
            "input_role":"declared_finite_model"}


def clock_and_shift(order):
    eye = sp.eye(order)
    phase = -1 if order==2 else I
    clock = sp.diag(*(phase**k for k in range(order)))
    shift = eye[:,[(k+1)%order for k in range(order)]]
    require(clock*shift==phase*shift*clock and shift**order==eye,"clock and shift covariance failed")
    return clock,shift,phase


def finite_gauge_link(order=4):
    require(type(order) is int and order in (2,4),"binary and fourth-root links are supported")
    # Two complex modes per cell, with A in the low occupation bits.
    creators = [creation(4,i) for i in range(4)]
    eye,link_eye = sp.eye(16),sp.eye(order)
    clock,shift,phase = clock_and_shift(order)
    number = [op*op.H for op in creators]
    charges_a = [int(number[0][j,j]+number[1][j,j]) for j in range(16)]
    charges_b = [int(number[2][j,j]+number[3][j,j]) for j in range(16)]
    tensor = sp.kronecker_product
    local_a = tensor(sp.diag(*(phase**k for k in charges_a)),clock)
    local_b = tensor(sp.diag(*(phase**k for k in charges_b)),clock.H)
    sector = [order*j+(-a)%order for j,(a,b) in enumerate(zip(charges_a,charges_b)) if (a+b)%order==0]
    selection = sp.eye(16*order)[:,sector]
    physical = selection*selection.H
    require(local_a*selection==selection and local_b*selection==selection,"wrong joint Gauss projector")
    gauge_average = sum((local_a**i*local_b**j for i,j in product(range(order),repeat=2)),
                        sp.zeros(16*order))/order**2
    require(gauge_average==physical,"Gauss average disagrees with the occupation census")
    a,b = creators[0],creators[2]
    bare = a*b.H
    dressed = tensor(bare,shift.H)
    require(all(dressed*gen==gen*dressed for gen in (local_a,local_b)),
            "dressed charge-one hop does not preserve both Gauss laws")
    require(not any(selection.H*tensor(bare,link_eye)*selection),"bare hop survived without the link")
    require(tensor(bare,link_eye)*local_a==sp.conjugate(phase)*local_a*tensor(bare,link_eye),
            "bare left endpoint has the wrong character")
    two_hops = tensor(creators[0]*creators[2].H*creators[1]*creators[3].H,shift.H**2)
    transfer = selection.H*(dressed+dressed.H)*selection
    require(transfer==transfer.H,"gauged transfer is not Hermitian")
    fixed = []
    for total in range(5):
        external = [order*j+(-a)%order for j,(a,b) in enumerate(zip(charges_a,charges_b)) if a+b==total]
        if not external:
            continue
        inclusion = sp.eye(16*order)[:,external]
        restored = inclusion.H*dressed*inclusion
        require(len(external)==comb(4,total) and local_a*inclusion==inclusion and
                local_b*inclusion==phase**total*inclusion,"incorrect boundary charge sector")
        fixed.append({"total_matter_number":total,"right_boundary_character":str(phase**total),
                      "dimension":len(external),"primitive_hop_rank":restored.rank()})
    return {"status":"EXACT_FINITE_GAUGE_LINK","group_order":order,"physical_dimension":selection.cols,
            "global_charge_condition":f"N_A+N_B=0 mod {order} on the closed link",
            "bare_hop_rank_after_gauss":0,"closed_primitive_hop_rank":(selection.H*dressed*selection).rank(),
            "closed_double_hop_rank":(selection.H*two_hops*selection).rank(),
            "closed_hamiltonian_spectrum":[[str(k),int(v)] for k,v in sorted(transfer.eigenvals().items())],
            "boundary_sectors":fixed,
            "majorana_primitive_i_character":"impossible by taking adjoints of a nonzero Hermitian operator",
            "r102_pairing_or_original_charge_derived":False}


def neutral_flux(inc,order):
    """Enumerate every neutral Abelian flux without enumerating tree edges."""
    basis = inc.nullspace()
    require(all(x.is_Integer for vector in basis for x in vector),"incidence cycle basis is not integral")
    cycles = sp.Matrix.hstack(*basis) if basis else sp.zeros(inc.cols,0)
    require(inc*cycles==sp.zeros(inc.rows,len(basis)),"cycle basis failed the Gauss condition")
    require(order**len(basis)<=4096,"flux enumeration exceeds the declared finite sector bound")
    labels = {tuple(sum(cycles[e,j]*value for j,value in enumerate(values))%order for e in range(inc.cols))
              for values in product(range(order),repeat=len(basis))}
    require(len(labels)==order**len(basis) and all(all(x%order==0 for x in inc*sp.Matrix(len(flux),1,flux))
                                                 for flux in labels),
            "neutral cycle space has an incorrect flux census")
    return len(basis),sorted(labels)


def link_sewing(order,graphs=None):
    require(type(order) is int and order in (2,4),"binary and fourth-root links are supported")
    clock,shift,phase = clock_and_shift(order)
    tensors = sp.kronecker_product
    eye = sp.eye(order)
    # Internal matter is fixed neutral. On two consistently oriented links,
    # its Gauss condition selects equal electric flux on the two segments.
    copy = sp.zeros(order**2,order)
    for q in range(order):
        copy[q*(order+1),q] = 1
    internal = tensors(clock.H,clock)
    physical = sum((internal**k for k in range(order)),sp.zeros(order**2))/order
    require(copy.H*copy==eye and copy*copy.H==physical,"coincident fluxes do not resolve the internal Gauss law")
    require(tensors(shift,shift)*copy==copy*shift,"charge-one shift does not sew")
    require(shift**(order//2)!=eye,"link collapsed to a proper quotient")
    require(tensors(copy,eye)*copy==tensors(eye,copy)*copy,"neutral link sewing is not associative")
    if graphs is None:
        graphs = [{"name":name,"vertices":n,"edges":edges} for name,n,edges in (
            ("path",4,[(0,1),(1,2),(2,3)]),("ring",4,[(0,1),(1,2),(2,3),(3,0)]),
            ("two_rings",5,[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]))]
    require(isinstance(graphs,list) and bool(graphs),"declare one or more oriented graphs")
    topology = []
    for graph in graphs:
        name,n,edges = (graph[k] for k in ("name","vertices","edges"))
        require(isinstance(name,str) and bool(name) and type(n) is int and 1<=n<=32 and
                isinstance(edges,list) and len(edges)<=48,"invalid finite oriented graph")
        require(all(isinstance(edge,(list,tuple)) and len(edge)==2 and all(type(v) is int and 0<=v<n for v in edge)
                    and edge[0]!=edge[1] for edge in edges),"invalid oriented edges")
        inc = sp.zeros(n,len(edges))
        for j,(left,right) in enumerate(edges):
            inc[left,j],inc[right,j] = 1,-1
        loops,labels = neutral_flux(inc,order)
        topology.append({"graph":name,"first_betti":loops,"neutral_flux_sectors":len(labels),
                         "edges":[list(e) for e in edges]})
    return {"order":order,"local_flux_constraint":"q_left=q_right","associative_isometry":True,
            "shift_intertwiner":True,"primitive_action_order":order,"topologies":topology,
            "scope":"neutral fixed internal matter, finite Abelian link kinematics"}


def return_forms(metric,selector,first,last,step):
    """Decompose an ordered transfer before the one outer physical selection."""
    from language import positive
    eye = sp.eye(metric.rows)
    require(selector*selector==selector and metric*selector==selector.H*metric,
            "return selector must be an orthogonal projector in the source metric")
    require(positive(metric*step) and positive(metric*(eye-step)),
            "return step must be a positive contraction in the source metric")
    q = eye-selector
    direct = selector*last*selector*step*selector*first*selector
    outside = selector*last*q*step*q*first*selector
    cross = selector*last*selector*step*q*first*selector + selector*last*q*step*selector*first*selector
    whole = selector*last*step*first*selector
    require((whole-direct-outside-cross).applyfunc(sp.expand)==sp.zeros(metric.rows),
            "ordered return decomposition failed")
    return whole,direct,outside,cross
