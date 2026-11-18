"""Executable necessary tests for six TFPT origin routes on one source."""
from itertools import combinations, product
from math import prod
import sympy as sp

from language import positive
from reconstruction import matrix, scalar, reconstruct, regional_check, require, support_rank

ROUTES = ("sewing","markov","moments","regional","inflow","scattering")
OBLIGATIONS = {
    "sewing":"Full physical gluing, original primitive local charge, stress and global anomaly uniqueness.",
    "markov":"Original overlapping states, compatible recoveries over a covering family, and controlled limit.",
    "moments":"Original complete word grammar, flat extension relations and uniqueness in the physical class.",
    "regional":"Common standard state, compatible modular inclusions, translation implementation and local continuum.",
    "inflow":"Original Dirac family, global anomaly/inflow and a controlled mirror gap.",
    "scattering":"Full analytic and unitary crossing problem, infrared control, uniqueness and local reconstruction.",
}


def partial_trace(rho,dims,keep):
    coordinates = list(product(*(range(n) for n in dims)))
    size = prod(dims[k] for k in keep)
    output = sp.zeros(size)

    def index(c):
        offset = 0
        for k in keep:
            offset = offset*dims[k]+c[k]
        return offset

    trace_out = set(range(len(dims)))-set(keep)
    for i,a in enumerate(coordinates):
        for j,b in enumerate(coordinates):
            if all(a[k]==b[k] for k in trace_out):
                output[index(a),index(b)] += rho[i,j]
    return output


def psd_power(operator,power):
    # Exact Hermitian spectral decomposition. The inverse uses the support.
    value = sp.zeros(operator.rows)
    for eigenvalue,_,vectors in operator.eigenvects():
        if not eigenvalue:
            continue
        basis = sp.Matrix.hstack(*vectors)
        value += eigenvalue**power*basis*(basis.H*basis).inv()*basis.H
    return value.applyfunc(sp.simplify)


def markov(data):
    dims = data["dimensions"]
    require(len(dims)==3 and all(type(d) is int and d>0 for d in dims),"A,B,C dimensions required")
    a,b,c = dims
    ab,bc = matrix(data["ab"]),matrix(data["bc"])
    require(ab.shape==(a*b,a*b) and bc.shape==(b*c,b*c),"marginal dimension mismatch")
    support_rank(ab)
    support_rank(bc)
    middle = partial_trace(ab,[a,b],[1])
    require(middle==partial_trace(bc,[b,c],[0]),"AB and BC have incompatible overlap")
    bridge = sp.kronecker_product(sp.eye(a),
        psd_power(bc,sp.Rational(1,2))*sp.kronecker_product(psd_power(middle,-sp.Rational(1,2)),sp.eye(c)))
    recovered = (bridge*sp.kronecker_product(ab,sp.eye(c))*bridge.H).applyfunc(sp.simplify)
    require(positive(recovered) and sp.simplify(recovered.trace()-1)==0,"invalid recovered state")
    overlap = partial_trace(recovered,dims,[0,1])==ab and partial_trace(recovered,dims,[1,2])==bc
    result = {"status":"PASS_SCOPED" if overlap else "REFUTED_SCOPED",
              "recovered_has_both_marginals":overlap,
              "scope":"finite Petz recovery on the declared overlap and support"}
    if "abc" in data:
        original = matrix(data["abc"])
        support_rank(original)
        require(partial_trace(original,dims,[0,1])==ab and partial_trace(original,dims,[1,2])==bc,
                "declared full state does not have supplied marginals")
        matches = original==recovered
        result["equals_supplied_full_state"] = matches
        result["status"] = "PASS_SCOPED" if overlap and matches else "REFUTED_SCOPED"
    return result


def sewing(data):
    a,b = matrix(data["first"]),matrix(data["second"])
    ab,reflected = matrix(data["composed"]),matrix(data["reflected_first"])
    require(b.cols==a.rows,"incompatible sewing interfaces")
    stages = data["invertible_layers"]
    require(all(type(x) is int for x in stages),"layer counts must be integers")
    charges = list(map(scalar,data["charge_labels"]))
    central = scalar(data["central_charge"])
    checks = {"composition":b*a==ab,"reflection":a.H==reflected}
    return {"status":"PASS_SCOPED" if all(checks.values()) else "REFUTED_SCOPED",
            "checks":checks,"same_charges_central_charges":[str(central+8*k) for k in stages],
            "charge_labels":list(map(str,charges)),
            "scope":"specified matrix gluings and invertible-layer ambiguity, no universal gluing theorem"}


def moments(data):
    gram = matrix(data["gram"])
    words = data["words"]
    require(len(words)==gram.rows and all(isinstance(w,list) and
            all(isinstance(k,str) for k in w) for w in words),"labelled moments required")
    tuples = [tuple(w) for w in words]
    require(len(set(tuples))==len(tuples) and () in tuples,"distinct words including the unit required")
    forms = {k:matrix(v) for k,v in data["forms"].items()}
    require(all(all(k in forms for k in w) for w in tuples),"word uses an undeclared operator")
    depth = data["base_depth"]
    require(type(depth) is int and 0<=depth and any(len(w)>depth for w in words),"two word levels required")
    low = [i for i,w in enumerate(words) if len(w)<=depth]
    require(positive(gram),"moment Gram is not positive")
    lower,upper = gram.extract(low,low).rank(),gram.rank()
    if lower!=upper:
        return {"status":"REFUTED_SCOPED","lower_rank":lower,"upper_rank":upper,
                "scope":"flatness at the supplied levels"}
    identity_index = tuples.index(())
    rep,(vacuum,row,ops) = reconstruct(gram,forms,identity_index,adjoints=data["adjoints"])
    pivots = rep["basis_indices"]
    metric = gram.extract(pivots,pivots)
    inv = metric.inv()
    for j,word in enumerate(tuples):
        vector = vacuum
        for label in reversed(word):
            vector = ops[label]*vector
        require(vector==inv*gram.extract(pivots,[j]),"declared left multiplication does not reproduce word columns")
    return {"status":"PASS_SCOPED","lower_rank":lower,"upper_rank":upper,
            "representation":rep,"scope":"finite closed representation of the supplied marked moments"}


def regional(data):
    metric = matrix(data["metric"])
    require(positive(metric) and metric.det()!=0,"positive definite regional metric required")
    checks = regional_check(metric,[matrix(a) for a in data["translations"]],
        {name:[matrix(x) for x in terms] for name,terms in data["regions"].items()},
        [{**t,"operator":matrix(t["operator"])} for t in data["transports"]])
    passes = (checks["commuting_generators"] and checks["positive_joint_spinor_form"]
              and checks["spatial_rank_modulo_identity_and_time"]==3 and bool(checks["transports"])
              and all(x["transports_full_algebra"] and x["commutes_with_supplied_translations"]
                      for x in checks["transports"]))
    return {"status":"PASS_SCOPED" if passes else "REFUTED_SCOPED","checks":checks}


def inflow(data):
    d = matrix(data["chiral_dirac_plus"])
    index = d.cols-d.rows
    charges = [scalar(q) for q in data["signed_u1_charges"]]
    require(all(q.is_Rational for q in charges),"rational signed U(1) charges required")
    linear = sum(charges)+scalar(data["inflow_linear"])
    cubic = sum(q**3 for q in charges)+scalar(data["inflow_cubic"])
    gap = scalar(data["mirror_gap_squared"])
    m = matrix(data["mirror_mass"])
    require(gap.is_Rational and gap>0,"positive rational mirror mass bound required")
    checks = {"finite_chiral_index":index,"kernel_plus":d.cols-d.rank(),
              "kernel_minus":d.rows-d.rank(),"u1_linear":str(linear),"u1_cubic":str(cubic),
              "finite_mirror_gap_verified":positive(m.H*m-gap*sp.eye(m.cols))}
    passes = index==data["target_index"] and linear==0 and cubic==0 and checks["finite_mirror_gap_verified"]
    return {"status":"PASS_SCOPED" if passes else "REFUTED_SCOPED","checks":checks,
            "scope":"finite chiral block index, abelian perturbative sums and supplied mirror block only"}


def scattering(data):
    s,t,u = sp.symbols("s t u")
    coeff = data["symmetric_contact_coefficients"]
    values = [scalar(v) for v in coeff]
    require(all(v.is_Rational for v in values),"real rational coefficients required")
    # Crossing symmetric contact basis. Contact deformations change no pole residues.
    polynomial = sum(v*(s**(i+2)+t**(i+2)+u**(i+2)) for i,v in enumerate(values))
    swapped = polynomial.xreplace({s:t,t:s})
    channel = matrix(data["finite_unitary"])
    couplings = [scalar(g) for g in data["soft_spin2_couplings"]]
    require(len(couplings)>=2,"at least two species needed for soft coupling check")
    checks = {"crossing_contact":sp.expand(polynomial-swapped)==0,
              "finite_channel_unitary":channel.rows==channel.cols and channel.H*channel==sp.eye(channel.rows),
              "universal_supplied_soft_coupling":len(set(couplings))==1,
              "contact_polynomial":str(polynomial),"changes_pole_residues":False}
    return {"status":"PASS_SCOPED" if all(checks[k] for k in (
            "crossing_contact","finite_channel_unitary","universal_supplied_soft_coupling")) else "REFUTED_SCOPED",
            "checks":checks,"scope":"independent necessary controls; no full S matrix unitarity or QFT reconstruction"}


CHECKERS = dict(zip(ROUTES,(sewing,markov,moments,regional,inflow,scattering)))


def examine(candidate):
    source = candidate["source_id"]
    require(isinstance(source,str) and source.strip(),"same candidate source identity required")
    require(set(candidate["routes"])<=set(ROUTES),"unknown origin route")
    results = {}
    for key,checker in CHECKERS.items():
        item = candidate["routes"].get(key)
        if item is None:
            results[key] = {"status":"SOURCE_REQUIRED","missing":OBLIGATIONS[key]}
            continue
        require(item["source_id"]==source,"cross-source evidence cannot be combined")
        results[key] = checker(item["input"])
    joint = (joint_markov_moments(candidate,results)
             if all(key in candidate["routes"] for key in ("markov","moments"))
             else {"status":"SOURCE_REQUIRED","missing":"same candidate Markov and word data"})
    combinations_checked = []
    for size in range(2,7):
        for group in combinations(ROUTES,size):
            statuses = [results[r]["status"] for r in group]
            if {"markov","moments"} <= set(group):
                statuses.append(joint["status"])
            verdict = ("REFUTED_SCOPED" if "REFUTED_SCOPED" in statuses else
                       "SOURCE_REQUIRED" if "SOURCE_REQUIRED" in statuses else "PASS_SCOPED")
            combinations_checked.append({"routes":group,"status":verdict})
    return {"status":"SIX_ROUTE_COMPARISON","source_id":source,"routes":results,
            "combinations":combinations_checked,"joint_markov_moments":joint,"physical_solution":False,
            "open_obligations":OBLIGATIONS}


def encoded(m):
    return [[([str(sp.re(x)),str(sp.im(x))] if sp.im(x) else str(x)) for x in row] for row in m.tolist()]


def examples():
    i2,x,z = sp.eye(2),sp.Matrix([[0,1],[1,0]]),sp.diag(1,-1)
    positive_markov = sp.diag(sp.Rational(1,2),0,0,0,0,0,0,sp.Rational(1,2))
    coherent = positive_markov.copy()
    coherent[0,7] = coherent[7,0] = sp.Rational(1,2)
    primitive = sp.Matrix([[0,0,0,0]])
    region_basis = [sp.eye(8),sp.diag(1,0,0,0,0,0,0,0)]
    x3 = sp.kronecker_product(x,x,x)

    for state_name,rho in (("classical_mixture",positive_markov),("coherent_ghz",coherent)):
        ab,bc = partial_trace(rho,[2,2,2],[0,1]),partial_trace(rho,[2,2,2],[1,2])
        for signs,clock,mirror_mass,contact in product((1,-1),("independent","spin"),
                                                       (1,0),(0,1)):
            label = f"finite-control:{state_name}:{signs}:{clock}:{mirror_mass}:{contact}"
            algebra = [sp.eye(8),x3,sp.eye(8)]
            from reconstruction import gram_of
            gram = gram_of(algebra,rho)
            form = sp.Matrix([[(rho*a.H*x3*b).trace() for b in algebra] for a in algebra])
            translations = [4*sp.eye(8),sp.kronecker_product(z,i2,i2),
                            sp.kronecker_product(i2,z,i2),sp.kronecker_product(i2,i2,z)]
            if clock=="spin":
                translations[2] = sp.kronecker_product(x,i2,i2)
            rows = {
                "sewing":{"first":encoded(x),"second":encoded(z),"composed":encoded(signs*z*x),
                          "reflected_first":encoded(x),"invertible_layers":[0,1,2],
                          "charge_labels":[0,1,2,3],"central_charge":8},
                "markov":{"dimensions":[2,2,2],"ab":encoded(ab),"bc":encoded(bc),"abc":encoded(rho)},
                "moments":{"words":[[],["X"],["X","X"]],"base_depth":1,
                           "gram":encoded(gram),"forms":{"X":encoded(form)},"adjoints":{"X":"X"}},
                "regional":{"metric":encoded(sp.eye(8)),"translations":[encoded(t) for t in translations],
                            "regions":{"A":[encoded(a) for a in region_basis],
                                       "B":[encoded(a) for a in region_basis]},
                            "transports":[{"from":"A","to":"B","operator":encoded(sp.eye(8))}]},
                "inflow":{"chiral_dirac_plus":encoded(primitive),"target_index":3,"signed_u1_charges":[1,-1],
                          "inflow_linear":0,"inflow_cubic":0,"mirror_mass":[[mirror_mass]],"mirror_gap_squared":1},
                "scattering":{"symmetric_contact_coefficients":[contact],"finite_unitary":encoded(i2),
                              "soft_spin2_couplings":[1,1]},
            }
            yield {"source_id":label,"routes":{k:{"source_id":label,"input":v} for k,v in rows.items()},
                   "markov_moment_binding":{"word_matrices":[encoded(w) for w in algebra],
                                             "generators":{"X":encoded(x3)}}}


def joint_markov_moments(candidate,results):
    data = candidate.get("markov_moment_binding")
    if not data:
        return {"status":"SOURCE_REQUIRED","missing":"word matrices on the reconstructed overlap state"}
    if any(results[k]["status"]!="PASS_SCOPED" for k in ("markov","moments")):
        return {"status":"REFUTED_SCOPED","reason":"at least one necessary route check failed"}
    markov_data = candidate["routes"]["markov"]["input"]
    if "abc" not in markov_data:
        return {"status":"SOURCE_REQUIRED","missing":"full state or explicitly bound Petz reconstruction"}
    from reconstruction import gram_of
    words = [matrix(w) for w in data["word_matrices"]]
    moments_data = candidate["routes"]["moments"]["input"]
    gram = matrix(moments_data["gram"])
    rho = matrix(markov_data["abc"])
    labels = list(moments_data["forms"])
    declarations = data.get("generators",{})
    require(set(declarations)==set(labels),"every moment generator must act on the Markov carrier")
    ops = {label:matrix(a) for label,a in declarations.items()}
    valid = len(words)==len(moments_data["words"]) and gram_of(words,rho)==gram
    for word,recorded in zip(moments_data["words"],words):
        actual = sp.eye(rho.rows)
        for label in word:
            actual *= ops[label]
        valid &= actual==recorded
    for label,form in moments_data["forms"].items():
        actual = sp.Matrix([[(rho*a.H*ops[label]*b).trace() for b in words] for a in words])
        valid &= actual==matrix(form)
    return {"status":"PASS_SCOPED" if valid else "REFUTED_SCOPED",
            "common_state_and_operator_products_replayed":bool(valid)}


def run():
    rows = []
    for candidate in examples():
        outcome = examine(candidate)
        rows.append(outcome)
    return {"status":"SIX_ROUTE_VARIATIONS","variation_count":len(rows),
            "subset_count_per_variation":57,"variations":rows,"physical_solution":False}
