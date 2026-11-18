"""Exact searches on the reported four marks, integral glues and Clifford cells."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product

import sympy as sp

from reconstruction import matrix, require, scalar
from routes import partial_trace
from selection import BIT_VECTORS, full_face_refinements, index, integer_rows, mod2, perm_matrix
from sewing_dynamics import SIGNS, clifford_cells, joint_projector

SOURCE_REPORT = "originlab-four-mark-and-source"


def parity_bit(value):
    return int(value) % 2


def rank_probes(candidates, probes, observations):
    require(isinstance(observations,dict) and set(observations)<=set(probes),
            "observations must name probes in the declared alphabet")
    accepted = {probe:scalar(value) for probe,value in observations.items()}
    alive = [c for c in candidates if all(probes[key](c)==value for key,value in accepted.items())]
    questions = []
    for name, read in probes.items():
        if name in accepted:
            continue
        groups = Counter(read(c) for c in alive)
        separated = (len(alive)**2-sum(n*n for n in groups.values()))//2
        if separated:
            questions.append({"id":name,"separated_pairs":separated,
                              "answers":{str(value):n for value,n in sorted(groups.items())}})
    questions.sort(key=lambda item:-item["separated_pairs"])
    @lru_cache(maxsize=None)
    def depth(group):
        if len(group)<=1:
            return 0
        remaining = []
        for read in probes.values():
            answers = {read(c) for c in group}
            if len(answers)>1:
                remaining.append(1+max(depth(tuple(c for c in group if read(c)==key)) for key in answers))
        return min(remaining,default=float("inf"))

    for question in questions:
        read = probes[question["id"]]
        challenge = 1+max(depth(tuple(c for c in alive if read(c)==answer)) for answer in {read(c) for c in alive})
        question["worst_case_total_probes"] = int(challenge) if challenge!=float("inf") else None
    complete = [q for q in questions if q["worst_case_total_probes"] is not None]
    probe = min(complete,key=lambda q:q["worst_case_total_probes"]) if complete else next(iter(questions),None)
    branches = {}
    if probe:
        read = probes[probe["id"]]
        branches = {answer:[str(c) for c in alive if str(read(c))==answer] for answer in probe["answers"]}
    return alive, probe, branches


def phase_and_clock(data):
    candidates = full_face_refinements()["candidates"]
    require([c["arf"] for c in candidates]==[1,1],"expected exactly two odd refinements")
    b = matrix(candidates[0]["pairing"])
    a = sp.ones(4,1)
    t = mod2(sp.eye(4)+a*a.T*b)
    # The boundary of all four oriented faces vanishes modulo two.
    boundary = sp.Matrix([[int(k not in (i,j)) for k in range(4)] for i,j in combinations(range(4),2)])
    q0,q1 = (c["values"] for c in candidates)
    require(mod2(t*t)==sp.eye(4) and mod2(t.T*b*t)==b and mod2(boundary*t)==boundary and
            all(q0[index(t*sp.Matrix(v))]==q1[i] for i,v in enumerate(BIT_VECTORS)) and
            all(t*perm_matrix(p)==perm_matrix(p)*t for p in permutations(range(4))),
            "refinement equivalence failed")
    # Even subsets of the five q=0 points are the parity vectors of A4.
    points = [sp.Matrix(v) for v in BIT_VECTORS if any(v) and q0[index(v)]==0]
    pairings = [parity_bit((u.T*b*v)[0]) for u,v in combinations(points,2)]
    require(len(points)==5 and mod2(sum(points,sp.zeros(4,1)))==sp.zeros(4,1) and all(pairings),
            "odd form lost its five-point frame")
    images = {}
    for bits in product((0,1),repeat=5):
        if sum(bits)%2==0:
            code = index(sum((bits[i]*points[i] for i in range(5)),sp.zeros(4,1)))
            require(q0[code]==sum(bits)//2%2,"A4 quadratic norm failed")
            images[bits] = code
    require(len(set(images.values()))==16 and all(
        parity_bit(sum(x*y for x,y in zip(u,v))) ==
        parity_bit((sp.Matrix(BIT_VECTORS[i]).T*b*sp.Matrix(BIT_VECTORS[j]))[0])
        for u,i in images.items() for v,j in images.items()),"A4/2A4 isometry failed")
    clocks = [p for p in permutations(range(5)) if sum(i==x for i,x in enumerate(p))==1 and
              all(p[p[p[p[i]]]]==i for i in range(5)) and any(p[p[i]]!=i for i in range(5))]
    splits = {tuple(p[p[i]] for i in range(5)) for p in clocks}
    require(len(clocks)==30 and len(splits)==15,"wrong clock or involution count")
    observed = data.get("observations",{})
    require(isinstance(observed,dict) and set(observed)<={"phase","clock"},"unknown observation group")
    phase_probes = {f"face_square_{i+1}":lambda k,i=i: candidates[k]["squares"][i] for i in range(4)}
    phases, phase_probe, phase_branches = rank_probes(range(2),phase_probes,observed.get("phase",{}))
    clock_probes = {f"image_of_point_{i+1}":lambda c,i=i: c[i]+1 for i in range(5)}
    survivors, clock_probe, clock_branches = rank_probes(clocks,clock_probes,observed.get("clock",{}))
    quarter = perm_matrix(clocks[0])
    y = (sp.eye(5)-5*quarter**2)/12
    y_labels = [weight for weight,multiplicity in y.eigenvals().items() for _ in range(multiplicity)]
    require(Counter(y_labels)=={sp.Rational(-1,3):3,sp.Rational(1,2):2},"wrong hypercharge lift")
    charge_counts = Counter(sum((y_labels[i] for i,bit in enumerate(bits) if bit),sp.S.Zero) for bits in images)
    z = sp.Symbol("z")
    mobius = (z-1)/(z+1)
    orbit = [z]
    for _ in range(3):
        orbit.append(sp.cancel(orbit[-1].subs(z,mobius)))
    basis = [sp.S.One,z,1/z,1/(z-1),1/(z+1)]
    points_for_interpolation = [2,3,4,5,6]
    values = sp.Matrix([[f.subs(z,x) for f in basis] for x in points_for_interpolation])
    change = values.inv()*sp.Matrix([[f.subs(z,x) for f in [sp.S.One,*orbit]]
                                      for x in points_for_interpolation])
    require(change.det()!=0 and all(sp.cancel(f-g)==0 for f,g in zip(
        [sp.S.One,*orbit],sp.Matrix([[*basis]])*change)) and sp.cancel(orbit[-1].subs(z,mobius)-z)==0,
        "meromorphic clock intertwiner failed")
    # A4 discriminant norm for class k is k(5-k)/5. A glue map must
    # cancel the sum modulo 2 on every class.
    norm = lambda k: sp.Rational(k*(5-k),5)
    glues = [r for r in range(1,5) if all((norm(k)+norm(r*k%5))%2==0 for k in range(5))]
    root_counts = [40]+[sp.binomial(5,k)*sp.binomial(5,2*k%5) for k in range(1,5)]
    require(glues==[2,3] and sum(root_counts)==240,"A4 glue replay failed")
    # D5 vector representatives are the ten signed coordinate vectors;
    # each spinor coset has sixteen half-integer sign patterns.
    a3_norms = [sp.Rational(k*(4-k),4) for k in range(4)]
    d5_norms = [0,sp.Rational(5,4),1,sp.Rational(5,4)]
    a3_d5 = [4*3+4*sp.binomial(5,2)]+[sp.binomial(4,k)*d5_min for k,d5_min in
                                        enumerate([0,16,10,16]) if k and a3_norms[k]+d5_norms[k]==2]
    require(len(a3_d5)==4 and sum(a3_d5)==240,"A3/D5 glue replay failed")
    return {"status":"EXACT_FINITE_MARKED_SEARCH" if phases and survivors else "REFUTED_SUPPLIED_MARKS",
            "report_file_id":SOURCE_REPORT,"physical_source_derived":False,
            "phase":{"boundary_equivalence":integer_rows(t),"remaining":phases,"next_probe":phase_probe,
                     "branches":phase_branches,"scope":"marked face squares in the declared frame"},
            "clock":{"point_vectors":[[int(x) for x in point] for point in points],"initial_clocks":len(clocks),
                     "initial_involutions":len(splits),"remaining_clocks":len(survivors),
                     "remaining_involutions":len({tuple(p[p[i]] for i in range(5)) for p in survivors}),
                     "next_probe":clock_probe,"branches":clock_branches},
            "spinor":{"even_subset_count":len(images),"hypercharge_multiplicities":{
                str(weight):count for weight,count in sorted(charge_counts.items())}},
            "meromorphic_clock":{"change_of_basis_determinant":str(change.det()),
                                 "matrix":[list(map(str,row)) for row in change.tolist()],
                                 "column_order":["1","z","M(z)","M^2(z)","M^3(z)"]},
            "integral_glue":{"A4_A4_multipliers":glues,"root_counts":list(map(int,root_counts)),
                             "A3_D5_root_counts":list(map(int,a3_d5)),
                             "assumptions":["declared integral lattices and their shortest coset representatives",
                                            "even unimodular gluing"]}}


def cell_search(data):
    a,b,pa,pb = clifford_cells()
    eye = sp.eye(16)
    reflect = sp.eye(4)-sp.ones(4)/2
    x,y = sp.symbols("x y",real=True)
    general = x*sp.eye(4)+y*sp.ones(4)
    # Commuting with all permutations makes every diagonal entry equal
    # and every off-diagonal entry equal. Solve on the two eigenspaces.
    choices = {}
    for on_sum,on_complement in product((1,-1),repeat=2):
        orthogonal = general.subs({x:on_complement,y:sp.Rational(on_sum-on_complement,4)})
        name = ("" if on_complement==1 else "-")+("I" if on_sum==on_complement else "R")
        require(orthogonal in (sp.eye(4),-sp.eye(4),reflect,-reflect),"wrong sewing parametrization")
        choices[name] = orthogonal
    mixed_ops = [[sp.I*x*y for y in b] for x in a]
    signatures = {}
    for name,orthogonal in choices.items():
        require(orthogonal.T*orthogonal==sp.eye(4) and
                all(orthogonal*perm_matrix(p)==perm_matrix(p)*orthogonal for p in permutations(range(4))),
                "sewing matrix fails orthogonal diagonal symmetry")
        links = [sum((orthogonal[i,j]*mixed_ops[i][j] for j in range(4)),sp.zeros(16)).applyfunc(sp.expand)
                 for i in range(4)]
        require(all(op==op.H and (op**2).applyfunc(sp.expand)==eye for op in links) and
                all(not any((x*y-y*x).applyfunc(sp.expand)) for x,y in combinations(links,2)),
                "sewing involutions failed")
        levels = Counter()
        for signs in SIGNS:
            levels[signs.count(-1)] += int(sp.expand(joint_projector(links,signs).trace()))
        rho = joint_projector(links,(1,1,1,1))
        mixed = sp.Matrix(4,4,lambda i,j:sp.expand((rho*mixed_ops[i][j]).trace()))
        parity = (rho*pa*pb).trace()
        require(levels=={i:int(sp.binomial(4,i)) for i in range(5)} and rho.rank()==1 and
                all(partial_trace(rho,[4,4],[i])==sp.eye(4)/4 for i in (0,1)) and
                parity==orthogonal.det() and mixed==orthogonal,"sewing witness failed")
        signatures[name] = {"mixed":mixed,"parity":parity}
    probes = {f"mixed_{i+1}_{j+1}":lambda key,i=i,j=j:signatures[key]["mixed"][i,j]
              for i,j in product(range(4),repeat=2)}
    probes["joint_parity"] = lambda key:signatures[key]["parity"]
    alive,probe,branches = rank_probes(choices,probes,data.get("observations",{}))
    return {"status":"EXACT_FINITE_SEWING_SEARCH" if alive else "REFUTED_SUPPLIED_MARKS",
            "report_file_id":SOURCE_REPORT,"physical_source_derived":False,
            "shared":{"local_reduced_state":"I_4/4","spectrum_in_units_of_g":[1,4,6,4,1],
                      "gap":"g","g_identified":False},
            "remaining":alive,"next_probe":probe,"branches":branches,
            "witnesses":{name:{"mixed":[list(map(str,row)) for row in signatures[name]["mixed"].tolist()],
                               "ground_parity":str(signatures[name]["parity"])} for name in alive},
            "scope":"orthogonal S4-equivariant bilinear two-cell sewings in a fixed marked frame"}


def run(data=None):
    data = {} if data is None else data
    require(isinstance(data,dict) and set(data)<={"marks","cells"},"unknown marked search inputs")
    marks = phase_and_clock(data.get("marks",{}))
    cells = cell_search(data.get("cells",{}))
    return {"status":"REFUTED_SUPPLIED_MARKS" if "REFUTED_SUPPLIED_MARKS" in (marks["status"],cells["status"])
            else "EXACT_MARKED_SEARCH","physical_source_derived":False,"marks":marks,"cells":cells}
