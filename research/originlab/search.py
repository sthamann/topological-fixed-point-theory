"""Candidate enumeration, exact feedback, mutation and proposal import."""
from itertools import combinations, product
from random import Random
from collections import Counter, defaultdict
from hashlib import sha256
import json
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import sympy as sp

from language import Evaluator, Term, TypeMismatch, positive
from reconstruction import matrix as exact_matrix, require, scalar

UNARY = ("mirror","reverse","trace","fix")
BINARY = ("compose","glue","identify")
FAMILIES = {
    "composition": {"compose", "identify"},
    "boundary_gluing": {"glue","mirror","fix","reverse"},
    "automaton": {"compose","iterate","branch"},
    "spectral": {"trace","mirror","fix"},
    "encoding": {"lift","mirror","identify"},
    "selection": {"fix","glue","compose"},
}
# The families are finite search fragments. Their generators are explicitly supplied.
ATOMS = {
    "identity":sp.ImmutableMatrix([[1,0],[0,1]]),
    "swap":sp.ImmutableMatrix([[0,1],[1,0]]),
    "retain":sp.ImmutableMatrix([[1,0],[0,0]]),
    "read":sp.ImmutableMatrix([[1,0]]),
    "prepare":sp.ImmutableMatrix([[1],[1]]),
}


def enumerate_terms(atoms, size, allowed, parameters=(2,3)):
    layers = {1:[Term("atom",parameter=k) for k in sorted(atoms)]}
    yield from layers[1]
    for cost in range(2,size+1):
        layer = [Term(op,(arg,)) for op in UNARY if op in allowed for arg in layers[cost-1]]
        layer += [Term(op,(arg,),n) for op in ("branch","iterate","lift") if op in allowed
                  for arg in layers[cost-1] for n in parameters]
        layer += [Term(op,(a,b)) for op in BINARY if op in allowed
                  for i in range(1,cost-1) for a,b in product(layers[i],layers[cost-i-1])]
        layers[cost] = layer
        yield from layer


def mutations(terms, count, seed, operator_weights, operators, parameters):
    unary = tuple(op for op in (*UNARY,"branch","iterate","lift") if op in operators)
    binary = tuple(op for op in BINARY if op in operators)
    choices = unary+binary
    if not terms or not choices:
        return
    random = Random(seed)
    for _ in range(count):
        a = random.choice(terms)
        op = random.choices(choices,weights=[1+operator_weights[op] for op in choices])[0]
        args = (a,random.choice(terms)) if op in binary else (a,)
        parameter = random.choice(parameters) if op in ("branch","iterate","lift") else None
        yield Term(op,args,parameter)


def proposals(model, prompt, timeout):
    # Start the proposal provider in an empty read-only directory. Its
    # input is the bounded proposal prompt and its output is untrusted.
    with TemporaryDirectory(prefix="originlab-proposals-") as folder:
        out = Path(folder)/"proposal.json"
        command = ["codex","-a","never","exec","--ephemeral","--ignore-user-config","--skip-git-repo-check",
                   "-s","read-only","-C",folder,"-m",model,"-o",str(out),"-"]
        subprocess.run(command,input=prompt,text=True,capture_output=True,timeout=timeout,check=True)
        raw = json.loads(out.read_text())
        if not isinstance(raw,list) or len(raw)>32:
            raise ValueError("model must return a JSON array of at most 32 terms")
        terms = [Term.parse(data) for data in raw]
        if any(term.size>12 for term in terms):
            raise ValueError("a model proposal exceeds the prompt's term bound")
        return terms


def proposal_prompt(retained):
    return (
        'Return only a JSON array of at most 12 OriginDSL terms with at most 12 nodes each. '
        'Propose distinct positive contracting transfers with different ordered word responses. '
        'Use no tools and make no claims of physical derivation. Syntax: ["atom","swap"], '
        '["compose",A,B], ["glue",A,B], ["mirror",A], ["reverse",A], ["trace",A], '
        '["fix",A], ["branch",2,A], ["iterate",2,A], ["lift",2,A], ["identify",A,B]. '
        'Atom matrices: identity [[1,0],[0,1]], swap [[0,1],[1,0]], retain [[1,0],[0,0]], '
        'read [[1,0]], prepare [[1],[1]]. Existing finite frontier: '
        +json.dumps(retained[:8],separators=(",",":"))
    )


def search_profile(profile):
    if profile is None:
        return Evaluator(ATOMS,12),FAMILIES,(2,3),{}, {},()
    atoms,boundaries,seams,assumptions,dimensions = {},{},{},{},{}
    for name,row in profile["atoms"].items():
        require(isinstance(name,str) and name and isinstance(row.get("source"),str) and row["source"].strip(),
                "every custom primitive needs a label and a source or assumption")
        m = exact_matrix(row["matrix"])
        before,after = row["boundary"]
        for label,dim in ((before,m.cols),(after,m.rows)):
            require(isinstance(label,str) and label and dimensions.setdefault(label,dim)==dim,
                    "each named boundary must have one carrier dimension")
        atoms[name] = sp.ImmutableMatrix(m)
        boundaries[name] = (before,after)
        assumptions[name] = row["source"]
    for row in profile.get("seams",[]):
        before,after = row["boundary"]
        require(before in dimensions and after in dimensions and (before,after) not in seams and
                isinstance(row.get("source"),str) and row["source"].strip(),
                "each seam needs a unique declared boundary pair and a source")
        m = exact_matrix(row["matrix"])
        require(m.shape==(dimensions[after],dimensions[before]),"seam has wrong shape")
        seams[(before,after)] = sp.ImmutableMatrix(m)
        assumptions[f"seam:{before}:{after}"] = row["source"]
    require(atoms,"at least one primitive required")
    families = {k:set(v) for k,v in profile.get("families",FAMILIES).items()}
    allowed = set().union(*families.values())
    require(families and allowed <= set((*UNARY,*BINARY,"branch","iterate","lift")),
            "invalid generator families")
    if profile.get("combine_families",True) and len(families)>1:
        families["combined"] = allowed
    parameters = tuple(sorted(set(profile.get("parameters",[2,3]))))
    require(parameters and all(type(p) is int and 1<=p<=8 for p in parameters),"invalid expansion parameters")
    probes = {}
    for name,probe in profile.get("probes",{}).items():
        before = probe["boundary"]
        d = dimensions[before]
        state = exact_matrix(probe["state"])
        ops = {k:exact_matrix(v) for k,v in probe["operations"].items()}
        require(state.shape==(d,1) and (state.H*state)[0]==1 and "candidate" not in ops and
                all(m.shape==(d,d) for m in ops.values()),"probe requires normalized state and common carrier")
        costs = {"candidate":sp.S.One,**{k:scalar(v) for k,v in probe.get("costs",{}).items()}}
        require(set(costs)<=(set(ops)|{"candidate"}) and all(v.is_Rational and v>0 for v in costs.values()),
                "each test operation cost must be positive and rational")
        probes[name] = (before,state,ops,costs)
    bound = profile.get("max_dimension",12)
    require(type(bound) is int and bound>=1,"positive boundary dimension bound required")
    return Evaluator(atoms,bound,boundaries,seams),families,parameters,probes,assumptions,tuple(profile.get("premises",()))


def profile_job(spec, previous=None):
    previous = previous or {}
    return run(spec["size"],spec["seed"],spec.get("mutations",128),spec.get("depth",5),
               [Term.parse(x) for x in spec.get("suggested",[])],
               [Term.parse(x) for x in previous.get("retained",[])],
               previous.get("questions",[]),spec.get("observations",[]),
               spec.get("profile"),previous.get("excluded",[]))


def read_word(state, ops, word):
    vector = state
    for label in reversed(word):
        vector = ops[label]*vector
    return (state.H*vector)[0]


def pareto(survivors, evaluator, assumptions):
    # These objectives use only executable finite certificates and
    # assumptions carried by the actual term. No proposal supplies a score.
    costs = {}
    for i,(term,_,_) in enumerate(survivors):
        pending = [term]
        bindings = set()
        while pending:
            part = pending.pop()
            if part.op=="atom":
                bindings.add(part.parameter)
            elif part.op=="glue":
                seam = (evaluator.boundary(part.args[1])[1],evaluator.boundary(part.args[0])[0])
                if seam in evaluator.seams:
                    bindings.add(f"seam:{seam[0]}:{seam[1]}")
            pending.extend(part.args)
        inputs = sorted(assumptions.get(label,label) for label in bindings)
        costs[i] = {"size":term.size,"finite_assumptions":sorted(set(inputs)),
                    "finite_certificates":["positive","contractive"],
                    "independent_predictions":None,"original_derived":None,"free_physical_parameters":None}
    pairs = [(entry["size"],len(entry["finite_assumptions"])) for entry in costs.values()]
    unique = set(pairs)
    best = {p for p in unique if not any(q!=p and q[0]<=p[0] and q[1]<=p[1] for q in unique)}
    return costs,[i for i,p in enumerate(pairs) if p in best]


def default_probe(size):
    return sp.eye(size)[:,0],{1:sp.eye(size)[::-1,:]},{}


def observable_questions(survivors, max_depth, probes, model_limit):
    groups = defaultdict(list)
    for term,m,boundary in survivors:
        if not probes:
            groups[(m.rows,boundary,None)].append((term,m))
        for label,(before,state,_,_) in probes.items():
            if before==boundary[0] and state.rows==m.rows:
                groups[(m.rows,boundary,label)].append((term,m))
    ranked,coverage = [],[]
    for (size,boundary,label),models in groups.items():
        selected = models[:model_limit]
        if len(selected)<2:
            continue
        state,ops,costs = probes[label][1:] if label is not None else default_probe(size)
        candidate = "candidate" if label is not None else 0
        alphabet = [candidate,*sorted(ops)]
        model_ops = [{**ops,candidate:m} for _,m in selected]
        found = []
        for depth in range(1,max_depth+1):
            for word in product(alphabet,repeat=depth):
                if candidate not in word:
                    continue
                values = [read_word(state,letters,word) for letters in model_ops]
                gain = len(values)*(len(values)-1)//2-sum(n*(n-1)//2 for n in Counter(values).values())
                if not gain:
                    continue
                cost = sum(costs.get(k,1) for k in word)
                found.append((sp.Rational(gain,cost),gain,word,cost,values))
        found.sort(key=lambda x:(-x[0],-x[1],len(x[2]),x[2]))
        coverage.append({"dimension":size,"boundary":boundary,"probe":label,
                         "compared":len(selected),"available":len(models),"depth":max_depth,
                         "tested_word_count":sum(len(alphabet)**k-(len(alphabet)-1)**k for k in range(1,max_depth+1))})
        ranked.extend((score,gain,size,label,boundary,word,cost,values,selected)
                      for score,gain,word,cost,values in found[:12])
    ranked.sort(key=lambda x:(-x[0],-x[1],x[2],str(x[3]),len(x[5]),str(x[5])))
    return [dict(dimension=size,boundary=boundary,probe=label,word=list(word),pair_separation=gain,
                 declared_cost=str(cost),pair_separation_per_cost=str(score),
                 pair_separation_per_letter=str(sp.Rational(gain,len(word))),responses=list(map(str,values)),
                 models=[term.data() for term,_ in selected])
            for score,gain,size,label,boundary,word,cost,values,selected in ranked[:12]],coverage


def run(size, seed, evolution, max_depth, suggested=(), retained=(), prior_questions=(), observations=(),
        profile=None, known_counterexamples=()):
    evaluator,families,parameters,probes,assumptions,premises = search_profile(profile)
    allowed_ops = set().union(*families.values())
    model_limit = (profile or {}).get("experiment_models",16)
    require(type(model_limit) is int and model_limit>=2,"at least two experiment models required")
    for observation in observations:
        d,word,evidence = (observation.get(k) for k in ("dimension","word","evidence_id"))
        probe = observation.get("probe")
        require(type(d) is int and d>=1 and isinstance(word,list) and word and
                isinstance(evidence,str) and evidence.strip(),"observation needs dimension, word and evidence ID")
        require(probe is None or probe in probes,"observation uses an undeclared probe")
        alphabet = {"candidate",*probes[probe][2]} if probe is not None else {0,1}
        require(all((type(k) is str if probe is not None else type(k) is int) and k in alphabet for k in word)
                and (probe is None or d==probes[probe][1].rows),"observation word or carrier mismatches its probe")
        scalar(observation["value"])
    semantic = {}
    negatives = []
    tested = set()
    seen_counts = defaultdict(int)

    def process(family, term):
        key = json.dumps(term.data(),separators=(",",":"))
        if key in tested:
            return
        tested.add(key)
        seen_counts[family] += 1
        if profile is not None:
            queue = [term]
            while queue:
                part = queue.pop()
                if (part.op!="atom" and part.op not in allowed_ops or
                        type(part.parameter) is int and part.parameter not in parameters):
                    negatives.append(dict(term=term.data(),family=family,status="outside_declared_grammar"))
                    return
                queue.extend(part.args)
        try:
            matrix = evaluator.evaluate(term)
        except TypeMismatch as error:
            negatives.append(dict(term=term.data(),family=family,status="ill_typed",witness=str(error)))
            return
        boundary = evaluator.boundary(term)
        signature = (boundary,matrix.shape,tuple(matrix))
        previous = semantic.get(signature)
        if previous is None or term.size < previous[0].size:
            semantic[signature] = (term,matrix,boundary)

    for family, allowed in families.items():
        for term in enumerate_terms(evaluator.atoms,size,allowed,parameters):
            process(family,term)
    for term in suggested:
        process("model_proposal",term)
    for term in retained:
        process("previous_survivor",term)
    parents = [term for term,m,boundary in semantic.values() if boundary[0]==boundary[1] and positive(m)]
    informative_parents = []
    for question in prior_questions:
        frequencies = Counter(question["responses"])
        if not frequencies:
            continue
        best = min(frequencies.values())
        informative_parents.extend(Term.parse(data)
                                   for data,value in zip(question["models"],question["responses"])
                                   if frequencies[value] == best)
    operator_weights = Counter(term.op for term in informative_parents)
    for term in mutations(informative_parents+parents,evolution,seed,operator_weights,allowed_ops,parameters):
        process("mutation",term)

    def failed_observation(m,boundary):
        for observed in observations:
            if m.rows != observed["dimension"]:
                continue
            probe = observed.get("probe")
            if probe is not None and probes[probe][0]!=boundary[0]:
                continue
            state,ops,_ = probes[probe][1:] if probe is not None else default_probe(m.rows)
            actual = read_word(state,{**ops,("candidate" if probe is not None else 0):m},observed["word"])
            expected = scalar(observed["value"])
            if actual != expected:
                return dict(status="counterexample",observed=observed,predicted=str(actual))
        return None

    scope = sha256(json.dumps([profile,observations],sort_keys=True).encode()).hexdigest()
    exact_scope = sha256(json.dumps([profile,"positive_contraction"],sort_keys=True).encode()).hexdigest()
    archived = {item["semantic_key"]:item for item in known_counterexamples
                if item.get("test_scope") in (scope,exact_scope)}
    reused = 0

    survivors = []
    for term,m,boundary in semantic.values():
        if m.rows != m.cols or boundary[0]!=boundary[1]:
            continue
        key = sha256(json.dumps([boundary,m.shape,list(map(str,m))]).encode()).hexdigest()
        evidence = archived.get(key)
        if evidence is not None:
            reused += 1
        elif not positive(m):
            evidence = dict(status="fails_finite_positivity",exact_matrix=str(m))
        elif not positive(sp.eye(m.rows)-m):
            evidence = dict(status="fails_finite_contraction",exact_matrix=str(m))
        elif m not in (sp.zeros(m.rows),sp.eye(m.rows)):
            evidence = failed_observation(m,boundary)
            if evidence is None:
                survivors.append((term,m,boundary))
        if evidence is not None:
            test_scope = scope if evidence["status"]=="counterexample" else exact_scope
            negatives.append({**evidence,"term":term.data(),"semantic_key":key,"test_scope":test_scope})
    survivors.sort(key=lambda item:(item[0].size,str(item[0].data())))
    objectives,frontier = pareto(survivors,evaluator,assumptions)
    questions,coverage = observable_questions(survivors,max_depth,probes,model_limit)
    # The finite compression identity is a diagnostic, not a requirement
    # that every observation space must be closed under the full transfer.
    defects = []
    for term,t,_ in survivors:
        v = sp.eye(t.rows)[:,0]
        defect = (v.H*t*t*v-v.H*t*v*v.H*t*v)[0]
        if defect:
            defects.append(dict(term=term.data(),compressed_defect=str(defect)))
    return dict(status="finite_search", grammar_families={k:sorted(v) for k,v in families.items()},
                premises=["exact Gaussian-rational linear maps","finite depth and boundary dimension",
                          "declared test preparation and readout",*premises],
                atoms={k:str(v) for k,v in evaluator.atoms.items()}, ast_size_bound=size,
                primitive_and_seam_sources=assumptions,declared_parameters=list(parameters),
                feedback={"prior_questions":len(prior_questions),"rare_response_parents":len(informative_parents),
                          "reused_negative_certificates":reused},
                canonical_matrices=len(semantic), evaluated_terms=len(tested), counts=dict(seen_counts),
                survivors=[dict(term=x.data(),matrix=str(y),size=x.size,boundary=boundary,
                                objectives=objectives[i],
                                evidence={"finite_psd":True,"finite_contraction":True,
                                          "independent_predictions":0,"original_derived":0,
                                          "physical_premises_evaluated":False,
                                          "free_physical_parameters":None})
                           for i,(x,y,boundary) in enumerate(survivors)],
                finite_pareto_frontier=[survivors[i][0].data() for i in frontier],
                excluded=negatives, supplied_finite_observations=list(observations),
                questions=questions, experiment_coverage=coverage, compression_defects=defects)
