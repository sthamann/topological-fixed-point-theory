"""Exact tests for regional transport, superselection and directed information."""
from collections import Counter, deque
from itertools import product
from math import gcd, prod

import sympy as sp

from language import positive
from reconstruction import matrix, require, scalar


def phase(turns):
    angle = turns if isinstance(turns,sp.Rational) else scalar(turns)
    require(angle.is_Rational,"an angle must be a rational multiple of pi")
    return sp.exp(sp.I*sp.pi*(angle % 2))


def zero(expr):
    return sp.simplify(sp.expand_complex(expr))==0


def relative_transition(data):
    q = data["charge"]
    require(type(q) is int and q>0,"positive integer determinant charge required")
    a,b = [scalar(data[key]) for key in ("phase_a_pi","phase_b_pi")]
    require(a.is_Rational and b.is_Rational,"rational cap phases required")
    theta = phase(b-a)
    # Work in the exact 2x2 logical sectors {empty, full}; no truncation
    # of the state or C=|empty><full| is involved on these supports.
    psi = sp.kronecker_product(sp.Matrix([1,phase(a)]),sp.Matrix([1,phase(b)]))/2
    c = sp.Matrix([[0,1],[0,0]])
    op = sp.kronecker_product(c.H,c)
    rho = psi*psi.H
    charges = (0,q,q,2*q)
    common,local = sp.Matrix(rho),sp.Matrix(rho)
    for i,j in product(range(4),repeat=2):
        if charges[i]!=charges[j]:
            common[i,j] = 0
        if i!=j:
            local[i,j] = 0
    value = sp.simplify((op*rho).trace())
    require(zero(value-theta/4) and zero((op*common).trace()-value) and
            (op*local).trace()==0,"relative coherence replay failed")
    holonomies = data.get("wilson_angles_pi",[])
    require(isinstance(holonomies,list),"Wilson paths must be a list")
    angles = [scalar(x) for x in holonomies]
    require(all(x.is_Rational for x in angles),"rational Wilson angles required")
    n_a,n_b = sp.diag(0,0,q,q),sp.diag(0,q,0,q)
    op_charges = [q,-q]
    require(all(n*op-op*n==charge*op for n,charge in zip((n_a,n_b),op_charges)),
            "determinant endpoints do not have the stated regional charges")
    moduli = data.get("moduli",[4])
    require(isinstance(moduli,list) and all(type(n) is int and n>0 for n in moduli),
            "positive group moduli required")
    return {"status":"EXACT_RELATIONAL_MODEL","charge":q,
            "neutral_under_common_rotation":True,"invariant_under_independent_rotations_without_link":False,
            "absolute_cap_coherence_after_u1_average":"0",
            "absolute_cap_coherence_after_zn_average":{
                str(n):(str(phase(a)/2) if q % n==0 else "0")
                for n in moduli},
            "bare_expectation":str(theta/4),"common_twirl_expectation":str(theta/4),
            "independent_twirl_expectation":"0","operator_endpoint_charges":op_charges,
            "required_wilson_endpoint_charges":[-q,q],
            "dressed_expectations":[str(phase(b-a+q*x)/4) for x in angles],
            "loop_for_two_paths":str(phase(q*(angles[0]-angles[1]))) if len(angles)==2 else None,
            "physical_regional_field_bound":False}


def network(data):
    """Read oriented unitary edges; compare paths up to vertex frame changes."""
    vertices = data["vertices"]
    require(isinstance(vertices,list) and vertices and
            all(isinstance(x,str) and x for x in vertices) and len(set(vertices))==len(vertices),
            "vertices must be distinct names")
    edges = data["edges"]
    require(isinstance(edges,dict),"named directed edges required")
    parsed,adjacent = {},{v:[] for v in vertices}
    dim = data.get("dimension",1)
    require(type(dim) is int and dim>=1,"positive common fibre dimension required")
    identity = sp.eye(dim)
    for name,item in sorted(edges.items()):
        u,v = item["from"],item["to"]
        require(name and u in adjacent and v in adjacent,"edge outside declared graph")
        if "phase_pi" in item:
            require(dim==1,"phase edge requires a one-dimensional fibre")
            transport = sp.Matrix([[phase(item["phase_pi"])]])
        else:
            transport = matrix(item["matrix"])
        require(transport.shape==identity.shape and all(zero(x) for x in transport.H*transport-identity),
                "edge transport must be unitary on the declared fibre")
        parsed[name] = (u,v,transport)
        adjacent[u].append((name,v,transport))
        adjacent[v].append((name,u,transport.H))

    def encode(m):
        return [[str(x) for x in row] for row in m.tolist()]

    def along(start,end,route):
        require(start in adjacent and end in adjacent,"path endpoints must be declared")
        u,transport = start,identity
        for name,direction in route:
            require(name in parsed and type(direction) is int and direction in (-1,1),
                    "path edges require a name and orientation")
            a,b,edge = parsed[name]
            if direction==-1:
                a,b,edge = b,a,edge.H
            require(a==u,"discontinuous ordered path")
            transport,u = edge*transport,b
        require(u==end,"path must end at the declared vertex")
        return transport.applyfunc(sp.simplify)

    frames,tree,components = {},set(),[]
    for root in vertices:
        if root in frames:
            continue
        frames[root] = identity
        component = [root]
        queue = deque(component)
        while queue:
            u = queue.popleft()
            for name,v,transport in adjacent[u]:
                if v not in frames:
                    frames[v] = transport*frames[u]
                    tree.add(name)
                    queue.append(v)
                    component.append(v)
        components.append(component)

    cycles = []
    for name,(u,v,transport) in parsed.items():
        if name not in tree:
            holonomy = (frames[v].H*transport*frames[u]).applyfunc(sp.simplify)
            cycles.append({"edge":name,"holonomy":encode(holonomy),"trace":str(holonomy.trace()),
                           "identity":all(zero(x) for x in holonomy-identity)})
    paths = []
    for pair in data.get("path_pairs",[]):
        start,end = pair["from"],pair["to"]
        transports = [along(start,end,route) for route in pair["routes"]]
        require(len(transports)==2,"each comparison requires two paths")
        loop = (transports[1].H*transports[0]).applyfunc(sp.simplify)
        paths.append({"from":start,"to":end,"ordered_transports":[encode(m) for m in transports],
                      "holonomy":encode(loop),"trace":str(loop.trace()),
                      "identity":all(zero(x) for x in loop-identity)})
    flat = all(c["identity"] for c in cycles)
    faces = []
    for face in data.get("face_boundaries",[]):
        holonomy = along(face["base"],face["base"],face["boundary"])
        faces.append({"id":face["id"],"holonomy":encode(holonomy),
                      "identity":all(zero(x) for x in holonomy-identity)})
    return {"status":"VALID_DECLARED_TRANSPORT","components":components,
            "connected":len(components)==1,"flat_graph_transport":flat,"all_graph_loops_trivial":flat,
            "flat_frame_freedom_per_component":len(components) if flat else None,
            "spanning_tree_edges":sorted(tree),"tree_frames":{k:encode(v) for k,v in frames.items()},
            "fundamental_cycles":cycles,"path_comparisons":paths,"face_holonomies":faces,
            "flat_on_supplied_faces":all(f["identity"] for f in faces) if faces else None,
            "scope":"ordered graph and face transport in the supplied connection"}


def oriented_shifts(data):
    """GNVW indices for products of explicitly declared infinite-chain shifts."""
    dims,steps,reflection = (data[k] for k in ("dimensions","shifts","reflection_permutation"))
    count = len(dims)
    require(count and len(steps)==len(reflection)==count and
            all(type(n) is int and n>=2 for n in dims) and all(type(s) is int for s in steps) and
            all(type(r) is int for r in reflection) and set(reflection)==set(range(count)),
            "invalid labelled shift factors")
    require(all(reflection[reflection[i]]==i and dims[reflection[i]]==d for i,d in enumerate(dims)),
            "reflection must preserve factor types and square to the identity")
    reversed_steps = [-steps[reflection[i]] for i in range(count)]
    factors = [sp.Integer(d)**n for d,n in zip(dims,steps)]
    index = prod(factors)
    invariant = reversed_steps==steps
    require(not invariant or index==1,"reflection and declared index disagree")
    return {"status":"EXACT_SHIFT_INDEX","factor_indices":list(map(str,factors)),
            "index":str(index),"reflected_index":str(1/index),"reflected_shifts":reversed_steps,
            "actual_labelled_reflection_invariant":invariant,
            "scope":"products of qudit shifts on the infinite one-dimensional chain"}


def covariant_stationary(data):
    """Check a finite CPTP channel, U(1) or Z_n covariance and stationarity."""
    charges = data["charges"]
    modulus = data.get("modulus")
    require(isinstance(charges,list) and charges and all(type(x) is int for x in charges),
            "integer charges in a common basis required")
    require(modulus is None or type(modulus) is int and modulus>0,"invalid group modulus")
    d = len(charges)
    superop = matrix(data["superoperator"])
    require(superop.shape==(d*d,d*d),"wrong superoperator dimension")
    ops = [superop[:,k].reshape(d,d) for k in range(d*d)]
    choi = sp.BlockMatrix([ops[i*d:(i+1)*d] for i in range(d)]).as_explicit()
    trace_preserving = all(op.trace()==int(i//d==i%d) for i,op in enumerate(ops))
    require(positive(choi) and trace_preserving,"map is not CPTP")
    differences = [a-b if modulus is None else (a-b)%modulus for a in charges for b in charges]
    violations = [[i,j] for i,j in product(range(d*d),repeat=2)
                  if differences[i]!=differences[j] and superop[i,j]!=0]
    fixed = (superop-sp.eye(d*d)).nullspace()
    result = {"status":"EXACT_FINITE_CHANNEL","covariant":not violations,
              "first_covariance_violation":violations[0] if violations else None,
              "fixed_operator_dimension":len(fixed),"unique_stationary_density":None}
    if len(fixed)==1:
        rho = fixed[0].reshape(d,d)
        require(rho.trace()!=0,"CPTP map must have a stationary density")
        rho = rho/rho.trace()
        require(positive(rho),"stationary operator is not a density")
        result["unique_stationary_density"] = [[str(v) for v in row] for row in rho.tolist()]
        result["has_cross_sector_stationary_coherence"] = any(
            differences[i]!=0 and x!=0 for i,x in enumerate(rho))
        require(violations or not result["has_cross_sector_stationary_coherence"],
                "unique covariant fixed density cannot mix distinct charges")
    return result


def gaussian(data):
    stiffness = matrix(data["stiffness"])
    require(stiffness.is_symmetric() and all(v.is_Rational for v in stiffness) and
            positive(stiffness) and stiffness.det()>0,"real positive definite stiffness required")
    covariance = stiffness.inv()
    blocks = data.get("blocks",[])
    require(not blocks or sorted(i for b in blocks for i in b)==list(range(stiffness.rows)),
            "blocks must partition all field coordinates")
    mixed = [{"coordinates":[i,j],"connected":str(covariance[i,j])}
             for k,m in enumerate(blocks) for n in blocks[k+1:] for i in m for j in n]
    return {"status":"EXACT_GAUSSIAN_COMMON_SOURCE","hessian":[list(map(str,row)) for row in covariance.tolist()],
            "mixed_responses":mixed,"factorizes_in_supplied_partition":
                all(scalar(q["connected"])==0 for q in mixed) if blocks else None}


def controls():
    vertices = list("ABCD")
    edges = {str(i):{"from":u,"to":v,"phase_pi":x} for i,(u,v,x) in enumerate(
        [("A","B","1/7"),("B","C","-2/7"),("C","D","1/5"),("D","A","-2/35")])}
    flat = network({"vertices":vertices,"edges":edges})
    edges["3"]["phase_pi"] = "29/105"  # -2/35 + 1/3
    curved_connection = {"vertices":vertices,"edges":edges,
                         "path_pairs":[{"from":"A","to":"D","routes":[
                             [["0",1],["1",1],["2",1]],[["3",-1]]]}],
                         "face_boundaries":[{"id":"square","base":"A",
                                             "boundary":[[str(i),1] for i in range(4)]}]}
    curved = network(curved_connection)
    cases = [relative_transition({"charge":q,"phase_a_pi":"1/7","phase_b_pi":"-1/5",
                                  "moduli":[3,4,8],"wilson_angles_pi":["1/24",0]})
             for q in (4,8,16)]
    coherent_reset = [["1/2",0,0,"1/2"] for _ in range(4)]
    dephased_reset = [[x,0,0,x] for x in ("1/2",0,0,"1/2")]
    return {"status":"EXACT_RELATIONAL_CONTROLS","relative_transitions":cases,"flat":flat,"curved":curved,
            "charged_path_combinations":[transported_caps({
                "caps":{"charge":q,"phase_a_pi":"1/7","phase_b_pi":"-1/5"},
                "connection":curved_connection}) for q in (6,16)],
            "charge_visibility":[charge_visibility({"charges":charges,"flux_pi":"1/8"})
                                 for charges in ([16],[4,16],[1,16])],
            "inverse_flux":flux_inference({"fundamental_loops":["square"],"phase_order":16,
                "observations":[{"charge":16,"winding":[1],"residue":0}],
                "allowed_probes":[{"charge":q,"winding":[1]} for q in (1,4,16)]}),
            "shift_pairs":[oriented_shifts({"dimensions":[d,d],"shifts":[1,-1],
                                           "reflection_permutation":[1,0]}) for d in (2,3,4,5)],
            "stationary_symmetry":[
                covariant_stationary({"charges":[0,16],"superoperator":dephased_reset}),
                covariant_stationary({"charges":[0,16],"superoperator":coherent_reset}),
                covariant_stationary({"charges":[0,16],"superoperator":coherent_reset,"modulus":4}),
            ],
            "mixed_gaussian":gaussian({"stiffness":[[1,"1/3"],["1/3",1]],"blocks":[[0],[1]]}),
            "single_generation":anomaly_control(),
            "original_scope":"None of these declared transport models supplies the original regional endpoints."}


def transported_caps(data):
    """Compare two determinant transitions on one supplied Abelian connection."""
    cap = relative_transition(data["caps"])
    require(data["connection"].get("dimension",1)==1,
            "determinant reference uses an Abelian connection")
    connection = network(data["connection"])
    q = cap["charge"]
    edges = data["connection"]["edges"]
    a,b = (scalar(data["caps"][k]) for k in ("phase_a_pi","phase_b_pi"))
    regional = []
    for pair in data["connection"].get("path_pairs",[]):
        angles = path_angles(edges,pair["routes"])
        amps = [phase(b-a+q*x)/4 for x in angles]
        h = phase(q*(angles[0]-angles[1]))
        require(zero(amps[0]-h*amps[1]),"charged path ratio disagrees with the loop")
        regional.append({"from":pair["from"],"to":pair["to"],"responses":list(map(str,amps)),
                         "unit_charge_flux_pi":str((angles[0]-angles[1])%2),
                         "holonomy_seen_by_charge":str(h),
                         "distinguished_by_this_charge":not zero(h-1)})
    return {"status":"EXACT_COMMON_TRANSPORT_MODEL","relative_cap":cap,
            "connection":connection,"regional":regional,"physical_original_fields_derived":False}


def path_angles(edges,routes):
    require(all("phase_pi" in e for e in edges.values()),"supply the unit-charge phases")
    return [sum(direction*scalar(edges[name]["phase_pi"]) for name,direction in route) for route in routes]


def charge_visibility(data):
    charges = data["charges"]
    require(isinstance(charges,list) and charges and all(type(q) is int and q!=0 for q in charges),
            "nonzero integer charges in one U(1) dictionary required")
    primitive = gcd(*charges)
    period = sp.Rational(2,primitive)
    flux = scalar(data["flux_pi"])
    require(flux.is_Rational,"flux must be a rational multiple of pi")
    visible = [q for q in charges if (q*flux)%2!=0]
    require(bool(visible)==(flux%period!=0),"charge visibility disagrees with the common kernel")
    return {"status":"EXACT_U1_CHARACTER_KERNEL","charges":charges,"gcd":primitive,
            "first_invisible_positive_flux_pi":str(period),"flux_pi":str(flux%2),
            "visible_to_charges":visible,"indistinguishable_from_zero":not visible,
            "scope":"one supplied U(1) connection and its integer characters"}


def flux_inference(data):
    """Enumerate flux classes on a finite phase grid and choose a separator."""
    from z3 import Int, Or, Solver, sat, unsat

    loops = data["fundamental_loops"]
    n = data["phase_order"]
    limit = data.get("candidate_limit",256)
    require(isinstance(loops,list) and loops and len(set(loops))==len(loops) and
            all(isinstance(name,str) and name for name in loops),"distinct fundamental loop names required")
    require(type(n) is int and n>=2 and type(limit) is int and limit>=1,
            "phase order and candidate limit must be positive integers")

    def read(item):
        q,winding = item["charge"],item["winding"]
        require(type(q) is int and q!=0 and len(winding)==len(loops) and
                all(type(x) is int for x in winding),"integer charge and winding for each loop required")
        return q,winding

    observations = []
    for item in data["observations"]:
        q,winding = read(item)
        residue = item["residue"]
        require(type(residue) is int and 0<=residue<n,"phase residue outside the supplied grid")
        observations.append((q,winding,residue))
    solver = Solver()
    solver.set(timeout=5000)
    unknowns = [Int(f"loop_{i}") for i in range(len(loops))]
    solver.add(*(x>=0 for x in unknowns),*(x<n for x in unknowns))
    solver.add(*(q*sum(w*x for w,x in zip(winding,unknowns)) % n==r
                 for q,winding,r in observations))
    found = []
    status = solver.check()
    while status==sat and len(found)<limit:
        model = solver.model()
        flux = [model[x].as_long() for x in unknowns]
        require(all(q*sum(w*x for w,x in zip(winding,flux)) % n==r for q,winding,r in observations),
                "inverse flux solution failed exact replay")
        found.append(flux)
        solver.add(Or(*(x!=r for x,r in zip(unknowns,flux))))
        status = solver.check()
    found.sort()
    ranked = []
    total = len(found)*(len(found)-1)//2
    for item in data["allowed_probes"]:
        q,winding = read(item)
        cost = scalar(item.get("cost",1))
        require(cost.is_Rational and cost>0,"positive exact probe cost required")
        responses = [q*sum(w*x for w,x in zip(winding,flux)) % n for flux in found]
        gain = total-sum(k*(k-1)//2 for k in Counter(responses).values())
        if gain:
            ranked.append((sp.Rational(gain,cost),gain,item,responses))
    ranked.sort(key=lambda row:(-row[0],-row[1],row[2]["charge"],row[2]["winding"]))
    return {"status":"EXHAUSTED_FINITE_FLUX_CLASS" if status==unsat else
                     "PARTIAL_FINITE_FLUX_CLASS" if status==sat else "SOLVER_UNKNOWN",
            "fundamental_loops":loops,"phase_order":n,"flux_residue_candidates":found,
            "next_probes":[{**item,"pair_separation":gain,"candidate_responses":responses}
                           for _,gain,item,responses in ranked[:8]],
            "scope":"finite grid of declared Abelian loop phases and allowed probes"}


def anomaly_control():
    # All entries are left-handed Weyl multiplets (colour, weak, Y).
    # The colour sign distinguishes the fundamental from its conjugate.
    reps = [(3,2,"1/6",1),(3,1,"-2/3",-1),(3,1,"1/3",-1),
            (1,2,"-1/2",0),(1,1,1,0),(1,1,0,0)]
    reps = [(c,w,scalar(y),sign) for c,w,y,sign in reps]
    linear = sum(c*w*y for c,w,y,_ in reps)
    cubic = sum(c*w*y**3 for c,w,y,_ in reps)
    colour_y = sum(w*y/2 for c,w,y,_ in reps if c==3)
    weak_y = sum(c*y/2 for c,w,y,_ in reps if w==2)
    colour_cubic = sum(w*sign for c,w,y,sign in reps)
    doublets = sum(c for c,w,_,_ in reps if w==2)
    return {"dimension":sum(c*w for c,w,_,_ in reps),"tr_Y":str(linear),"tr_Y3":str(cubic),
            "SU3_squared_Y":str(colour_y),"SU2_squared_Y":str(weak_y),
            "SU3_cubic_relative":str(colour_cubic),"weak_doublets":doublets,
            "indistinguishable_copy_counts_for_these_tests":list(range(1,5))}
