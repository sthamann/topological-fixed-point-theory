"""Exact reconstruction of finite marked quantum processes from coherent data."""
from collections import deque
from itertools import combinations, product
import re

import sympy as sp

from language import positive

PAULI = (
    sp.eye(2), sp.Matrix([[0, 1], [1, 0]]),
    sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1),
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def scalar(value):
    if type(value) is int:
        return sp.Integer(value)
    if isinstance(value, str) and re.fullmatch(r"-?[0-9]+(?:/[1-9][0-9]*)?", value):
        return sp.Rational(value)
    if isinstance(value, list) and len(value) == 2:
        return scalar(value[0])+sp.I*scalar(value[1])
    raise ValueError("expected a Gaussian rational")


def matrix(rows):
    require(isinstance(rows,list) and rows and all(isinstance(r,list) and r for r in rows),
            "matrix must be a nonempty rectangular list")
    require(len({len(r) for r in rows})==1,"ragged matrix")
    return sp.Matrix([[scalar(v) for v in r] for r in rows])


def support_rank(rho):
    require(positive(rho) and rho.trace()==1,"state must be normalized and positive")
    return rho.rank()


def gram_of(words, rho):
    support_rank(rho)
    d = rho.rows
    require(all(w.shape == (d,d) for w in words),"words must use the declared carrier")
    return sp.Matrix([[(rho*a.H*b).trace() for b in words] for a in words])


def density_reconstruction(gram, words, rho):
    rank = support_rank(rho)
    require(gram_of(words,rho)==gram,"supplied density and word matrices do not reproduce the original Gram")
    dimension = rho.rows
    algebra_rank = sp.Matrix.hstack(*(w.reshape(dimension**2,1) for w in words)).rank()
    return {"dimension":dimension,"algebra_span_rank":algebra_rank,"density_rank":rank,
            "gram_rank":gram.rank(),"complete_matrix_span":algebra_rank==dimension**2,
            "rank_formula_applies":algebra_rank==dimension**2}


def matrix_units(d):
    return [sp.eye(d)[:,i]*sp.eye(d)[j,:] for i,j in product(range(d),repeat=2)]


def tomography(states, qubits):
    require(states and all(r.shape == (2**qubits,2**qubits) for r in states),
            "states must use one declared tensor factorization")
    for rho in states:
        support_rank(rho)
    for order in range(qubits+1):
        for support in combinations(range(qubits),order):
            for axes in product(range(1,4),repeat=order):
                letters = [0]*qubits
                for site,axis in zip(support,axes):
                    letters[site] = axis
                probe = sp.kronecker_product(*(PAULI[i] for i in letters))
                values = [(rho*probe).trace() for rho in states]
                if len(set(values))>1:
                    return {"status":"DISTINGUISHED","minimum_support":order,
                            "probe":"".join("IXYZ"[i] for i in letters),"values":list(map(str,values))}
    return {"status":"SAME_DENSITY_MATRIX","tested_complete_pauli_basis":True}


def quotient(gram, forms):
    """Use a basis of the null quotient and keep its metric explicitly."""
    require(positive(gram) and gram.rank()>0,"nonzero positive coherent Gram form required")
    _, pivots = gram.rref()
    null = gram.nullspace()
    metric = gram.extract(pivots,pivots)
    inv = metric.inv()
    operators = {}
    for label,form in forms.items():
        require(form.shape == gram.shape,"operator form has wrong word coordinates")
        require(all(not any(form*v) and not any(v.H*form) for v in null),
                f"form {label} fails to descend through the Gram null space")
        operators[label] = inv*form.extract(pivots,pivots)
    return pivots,metric,operators


def energy_moment_check(gram, first, second):
    """Test closure under K=H-E0 from the supplied moments of one K >= 0."""
    pivots, metric, operators = quotient(gram, {"energy":first,"energy_squared":second})
    require(positive(first), "first energy form must be positive Hermitian")
    require(all(not any(second*vector) for vector in first.nullspace()),
            "a zero-energy direction has a nonzero second energy moment")
    compressed = operators["energy"]
    defect = metric*(operators["energy_squared"]-compressed**2)
    require(positive(defect), "negative energy moment compression defect")
    _, hidden_pivots = defect.rref()
    missing_rank = len(hidden_pivots)
    report = {
        "status":"GENERATOR_INVARIANT" if not missing_rank else "ADDITIONAL_DIRECTIONS_REQUIRED",
        "observed_rank":metric.rows,
        "compressed_generator":[[str(v) for v in row] for row in compressed.tolist()],
        "defect_form":[[str(v) for v in row] for row in defect.tolist()],
        "minimum_additional_directions":missing_rank,
        "single_word_energy_variances":[
            str((gram[i,i]*second[i,i]-first[i,i]**2)/gram[i,i]**2) if gram[i,i] else None
            for i in range(gram.rows)
        ],
        "scope":"supplied moments of the same nonnegative self-adjoint generator and word vectors",
        "physical_generator_or_absolute_reference_derived":False,
    }
    if missing_rank:
        # A principal Gram of the defect gives a rational hidden carrier.
        # The Schur complement is zero for the first extension and positive
        # after raising the hidden generator by one in these energy units.
        visible_energy = metric*compressed
        hidden_metric = defect.extract(hidden_pivots,hidden_pivots)
        coupling = defect.extract(hidden_pivots,range(metric.rows))
        lower = coupling*visible_energy.pinv()*coupling.H
        extension_metric = sp.diag(metric,hidden_metric)
        first_form = visible_energy.row_join(coupling.H).col_join(coupling.row_join(lower))
        raised_form = first_form + sp.diag(sp.zeros(metric.rows),hidden_metric)
        first_generator = extension_metric.inv()*first_form
        raised_generator = extension_metric.inv()*raised_form
        third_difference = (extension_metric*(raised_generator**3-first_generator**3)).extract(
            range(metric.rows),range(metric.rows))
        require(third_difference==defect,"constructed energy extensions disagree with their third-moment test")
        report["finite_extension_witness"] = {
            "gram":[[str(v) for v in row] for row in extension_metric.tolist()],
            "energy_forms":[[[str(v) for v in row] for row in form.tolist()]
                            for form in (first_form,raised_form)],
            "visible_basis_indices":list(pivots),
            "third_moment_difference":report["defect_form"],
            "separating_word_basis_index":pivots[next(i for i in range(metric.rows) if defect[i,i])],
        }
    return report


def reconstruct(gram, forms, identity_index, step=None, twice=None, adjoints=None):
    require(type(identity_index) is int and 0<=identity_index<gram.rows,"identity word index required")
    require(gram[identity_index,identity_index]==1,"identity word must have norm one")
    require(all(isinstance(k,str) and k for k in forms),"operation labels must be nonempty strings")
    if step is not None:
        require(step.shape == gram.shape and positive(step) and positive(gram-step),
                "one-step form must be a positive contraction relative to Gram")
        require("@transfer" not in forms,"@transfer is a reserved label")
        forms = {**forms,"@transfer":step}
    pivots,metric,operators = quotient(gram,forms)
    if adjoints is not None:
        require(set(adjoints)==set(operators) and all(adjoints.get(v)==k for k,v in adjoints.items()),
                "adjoints must specify an involution on all declared labels")
        for name,adjoint in adjoints.items():
            require(operators[name].H*metric==metric*operators[adjoint],
                    f"{name} has incompatible declared adjoint")
    omega = metric.inv()*gram.extract(pivots,[identity_index])
    row = omega.H*metric
    rank = len(pivots)
    span = reachable(omega,operators)
    result = {"gram_rank":rank,"reachable_rank":len(span),
              "representation_complete_for_supplied_basis":len(span)==rank,
              "operator_labels":sorted(operators),"adjoints_checked":adjoints is not None,
              "physical_local_algebra_identified":False}
    if twice is not None:
        require(step is not None and twice.shape==gram.shape,"two-step form requires a matching one-step form")
        _,_,second = quotient(gram,{"two":twice})
        t = operators["@transfer"]
        defect = metric*(second["two"]-t*t)
        require(positive(defect),"negative two-step compression defect")
        result["two_step_compression"] = {"defect_rank":defect.rank(),
                                         "closes_on_supplied_basis_at_two_steps":not any(defect)}
    result["quotient_metric"] = [[str(v) for v in row] for row in metric.tolist()]
    result["basis_indices"] = list(pivots)
    if adjoints is not None and len(span)==rank:
        result["scope"] = "cyclic finite star representation on the declared marked alphabet"
    return result,(omega,row,operators)


def reachable(initial,operators):
    """Retain original independent word columns. Their span is invariant."""
    basis = []
    queue = deque([((),initial)])
    while queue:
        word,vector = queue.popleft()
        if sp.Matrix.hstack(*(v for _,v in basis),vector).rank() == len(basis):
            continue
        basis.append((word,vector))
        queue.extend(((label,*word),op*vector) for label,op in sorted(operators.items()))
    return basis


def equivalent(left,right):
    """Decide equality of ALL words on one explicitly matched alphabet."""
    u,l,a = left
    v,r,b = right
    require(set(a)==set(b),"comparison needs a common operational alphabet")
    ops = {k:sp.diag(a[k],b[k]) for k in a}
    combined = l.row_join(-r)
    basis = reachable(u.col_join(v),ops)
    for word,vec in basis:
        if (combined*vec)[0]!=0:
            return {"status":"DISTINGUISHED","word":list(word),
                    "responses":[str((l*vec[:u.rows,:])[0]),str((r*vec[u.rows:,:])[0])]}
    return {"status":"ALL_WORDS_EQUAL","invariant_reachable_rank":len(basis),
            "dimension_bound":u.rows+v.rows,"alphabet":sorted(a),
            "scope":"finite representations and the supplied common marked alphabet only"}


def compare_documents(left,right):
    require(left["interface"]==right["interface"],"marked interface and conventions differ")
    l_check,l_model = from_document(left)
    r_check,r_model = from_document(right)
    report = equivalent(l_model,r_model)
    if report["status"] == "ALL_WORDS_EQUAL":
        report["cyclic_star_equivalence"] = (
            l_check["adjoints_checked"] and r_check["adjoints_checked"]
            and l_check["representation_complete_for_supplied_basis"]
            and r_check["representation_complete_for_supplied_basis"]
            and left["adjoints"]==right["adjoints"]
        )
        if report["cyclic_star_equivalence"]:
            u,_,a = l_model
            v,_,b = r_model
            basis = reachable(u,a)
            pairs = []
            for word,x in basis:
                y = v
                for label in reversed(word):
                    y = b[label]*y
                pairs.append((x,y))
            left_basis = sp.Matrix.hstack(*(x for x,_ in pairs))
            right_basis = sp.Matrix.hstack(*(y for _,y in pairs))
            mapping = right_basis*left_basis.inv()
            left_metric = sp.Matrix([[sp.sympify(x) for x in row] for row in l_check["quotient_metric"]])
            right_metric = sp.Matrix([[sp.sympify(x) for x in row] for row in r_check["quotient_metric"]])
            require(mapping.H*right_metric*mapping==left_metric and mapping*u==v and
                    all(mapping*a[k]==b[k]*mapping for k in a),
                    "all-word candidate did not replay as a metric-unitary intertwiner")
            report["metric_unitary_intertwiner"] = [[str(x) for x in row] for row in mapping.tolist()]
    return report


def from_document(data):
    require(isinstance(data.get("interface"),dict) and bool(data["interface"]),
            "explicit marked interface required")
    gram = matrix(data["gram"])
    forms = {k:matrix(v) for k,v in data["forms"].items()}
    step = matrix(data["step"]) if "step" in data else None
    twice = matrix(data["twice"]) if "twice" in data else None
    return reconstruct(gram,forms,data["identity_index"],step,twice,data.get("adjoints"))


def controls():
    zero = sp.eye(8)[:,0]
    one = sp.eye(8)[:,7]
    ghz = [(zero+sgn*one)*(zero+sgn*one).H/2 for sgn in (1,-1)]
    diagnostic = tomography(ghz,3)
    require(diagnostic["minimum_support"]==3 and diagnostic["probe"]=="XXX",
            "the GHZ phase should require a three-body probe")
    ranks = []
    for r in range(1,7):
        rho = sp.diag(*([sp.Rational(1,r)]*r+[0]*(6-r)))
        rank = gram_of(matrix_units(6),rho).rank()
        require(rank==6*r,"unexpected state-support rank")
        ranks.append({"density_rank":r,"gram_rank":rank,"algebra_span":36})
    initial = sp.Matrix([1,0])
    local = sp.Matrix([[0,1],[1,0]])
    # A change of basis and an unreachable spectator give a larger presentation.
    change = sp.Matrix([[1,1],[0,1]])
    equal = equivalent((initial,initial.H,{"X":local}),
        ((change*initial).col_join(sp.zeros(1,1)),
         (initial.H*change.inv()).row_join(sp.zeros(1,1)),
         {"X":sp.diag(change*local*change.inv(),sp.Rational(2,3))}))
    different = equivalent((initial,initial.H,{"X":local}),
                           (initial,initial.H,{"X":local/2}))
    return {"status":"EXACT_RECONSTRUCTION_CONTROLS","ghz":diagnostic,"rank_controls":ranks,
            "all_word_comparison":equal,"distinguishing_comparison":different,
            "physical_solution":False}


def regional_check(metric, translations, regions, transports):
    """Test supplied regional algebras and transport on the same quotient."""
    from bridges import SIGMA
    d = metric.rows
    require(len(translations)==4 and all(t.shape==(d,d) and t.H*metric==metric*t for t in translations),
            "four metric-self-adjoint translation generators required")
    commutes = all(a*b==b*a for i,a in enumerate(translations) for b in translations[:i])
    cone = sp.kronecker_product(metric*translations[0],sp.eye(2))
    for t,sigma in zip(translations[1:],SIGMA):
        cone += sp.kronecker_product(metric*t,sigma)
    inv = metric.inv()
    identity = sp.eye(d)

    def vec(a):
        return a.reshape(d*d,1)

    spans = {}
    motion_rows = []
    for name,words in regions.items():
        require(words and all(w.shape==(d,d) for w in words),"regional operator basis is invalid")
        cols = sp.Matrix.hstack(*(vec(w) for w in words))
        rank = cols.rank()
        for word in [identity,*(inv*w.H*metric for w in words),*(a*b for a,b in product(words,repeat=2))]:
            require(cols.row_join(vec(word)).rank()==rank,f"region {name} fails unital star algebra closure")
        spans[name] = (cols,rank)
        annihilator = cols.T.nullspace()
        for word in words:
            projected = [vec(generator*word-word*generator) for generator in translations]
            motion_rows.extend([[(functional.T*column)[0] for column in projected]
                                for functional in annihilator])
    motion = sp.Matrix(motion_rows) if motion_rows else sp.zeros(0,4)
    checks = []
    for entry in transports:
        a,b,u = entry["from"],entry["to"],entry["operator"]
        require(u.shape==(d,d) and u.H*metric*u==metric,"regional transport must be metric unitary")
        cols,rank = spans[b]
        moves = (rank==spans[a][1] and all(cols.row_join(vec(u*w*u.inv())).rank()==rank for w in regions[a]))
        distinct = cols.row_join(spans[a][0]).rank() > spans[a][1] or rank != spans[a][1]
        intertwines = all(u*t==t*u for t in translations)
        checks.append({"from":a,"to":b,"transports_full_algebra":moves,
                       "source_and_target_algebras_distinct":distinct,
                       "commutes_with_supplied_translations":intertwines})
    return {"commuting_generators":commutes,"positive_joint_spinor_form":positive(cone),
            "spatial_rank_modulo_identity_and_time":(
                sp.Matrix.hstack(*(vec(t) for t in [identity,*translations])).rank()
                - sp.Matrix.hstack(vec(identity),vec(translations[0])).rank()
            ),
            "regional_motion":{
                "spatial_rank_modulo_time_and_internal_actions":motion.rank()-motion[:,0:1].rank(),
                "invisible_generator_combinations":[[str(v) for v in column] for column in motion.nullspace()],
                "scope":"first-order commutator motion out of the supplied regional algebra spans",
            },
            "regional_ranks":{k:rank for k,(_,rank) in spans.items()},"transports":checks,
            "translation_generator_transport_link_proved":False,
            "scope":"finite simultaneous tests on a supplied carrier, no continuum local net theorem"}
