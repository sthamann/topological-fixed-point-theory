"""Compose operator fragments and search their joint ordered responses exactly."""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from math import prod

import sympy as sp

from language import Evaluator, TypeMismatch, positive
from reconstruction import equivalent, matrix, require, scalar
from search import enumerate_terms


def word_operator(operators, word, dimension):
    result = sp.eye(dimension)
    for letter in word:
        result = result*operators[letter]
    return result


def operational_model(weight, operators):
    """Represent Tr(weight * word) as a finite linear word model."""
    left,right = weight.rank_decomposition()
    rank = left.cols
    if not rank:
        return sp.zeros(1,1),sp.zeros(1,1),{name:sp.zeros(1,1) for name in operators}
    return (left.vec(), right.T.vec().T,
            {name:sp.kronecker_product(sp.eye(rank),op) for name,op in operators.items()})


class MetricEvaluator(Evaluator):
    def __init__(self, atoms, gram):
        super().__init__(atoms, gram.rows)
        self.gram = gram
        self.metric_inverse = gram.inv()

    @lru_cache(maxsize=16000)
    def evaluate(self, term):
        if term.op == "mirror":
            value = self.evaluate(term.args[0])
            return sp.ImmutableMatrix(self.metric_inverse*value.H*self.gram)
        if term.op == "fix":
            basis = (self.evaluate(term.args[0])-sp.eye(self.gram.rows)).nullspace()
            if not basis:
                return sp.ImmutableMatrix.zeros(self.gram.rows)
            columns = sp.Matrix.hstack(*basis)
            return sp.ImmutableMatrix(columns*(columns.H*self.gram*columns).inv()*columns.H*self.gram)
        return super().evaluate(term)


def compile_choices(slot, gram):
    dimension = gram.rows
    choices = list(slot.get("choices", []))
    if "grammar" in slot:
        require(len(slot["letters"]) == 1, "grammar slots generate one letter")
        grammar = slot["grammar"]
        atoms = {name:sp.ImmutableMatrix(matrix(rows)) for name,rows in grammar["atoms"].items()}
        size = grammar["size"]
        allowed = set(grammar["operations"])
        require(all(m.shape == (dimension, dimension) for m in atoms.values()) and atoms and
                type(size) is int and 1 <= size <= 6 and
                allowed <= {"compose","mirror","iterate","fix","reverse"},
                "invalid endomorphism grammar or size")
        evaluator = MetricEvaluator(atoms, gram)
        for term in enumerate_terms(atoms, size, allowed):
            try:
                value = evaluator.evaluate(term)
            except TypeMismatch:
                continue
            choices.append({"id":term.data(), "operators":{slot["letters"][0]:value}})
    require(choices and isinstance(slot["letters"], list) and slot["letters"] and
            len(set(slot["letters"])) == len(slot["letters"]), "empty or duplicate operation slot")
    distinct = {}
    for choice in choices:
        tags = choice.get("bindings", {})
        require(isinstance(tags, dict) and all(isinstance(k,str) and isinstance(v,str) for k,v in tags.items()),
                "shared bindings must be named strings")
        require(set(choice["operators"]) == set(slot["letters"]), "fragment must implement its whole slot")
        ops = {name:sp.ImmutableMatrix(value if isinstance(value,sp.MatrixBase) else matrix(value))
               for name,value in sorted(choice["operators"].items())}
        require(all(op.shape == (dimension, dimension) for op in ops.values()), "operator does not fit common carrier")
        key = (tuple(sorted(tags.items())), tuple(ops.values()))
        distinct.setdefault(key,(choice["id"],tags,ops))
    return list(distinct.values())


def word_tuple(raw, alphabet):
    require(isinstance(raw,list) and all(isinstance(k,str) and k in alphabet for k in raw),
            "unknown letter in ordered word")
    return tuple(raw)


def classify(survivors, budget):
    """Use exact invariant reachable spans, including trace memory."""
    classes, next_probe = [], None
    checks = 0
    complete = True
    for identity, ops, weight in survivors:
        if not complete:
            continue
        model = operational_model(weight,ops)
        for group in classes:
            if checks == budget:
                complete = False
                break
            checks += 1
            verdict = equivalent(group[1],model)
            if verdict["status"] == "ALL_WORDS_EQUAL":
                group[2] += 1
                break
            key = (len(verdict["word"]),verdict["word"])
            if next_probe is None or key < (len(next_probe["word"]),next_probe["word"]):
                next_probe = {**verdict,"left":group[0],"right":identity}
        else:
            classes.append([identity,model,1])
    return classes, next_probe, checks, complete


def compose(case, gram, weight, observations, boundary, bindings, counts, rejected):
    dimension = gram.rows
    require(dimension and gram.shape == weight.shape and positive(gram) and gram.det() != 0,
            "a nondegenerate positive common carrier and matching contraction required")
    slots = []
    alphabet = set()
    for slot in case["slots"]:
        letters = slot["letters"]
        require(isinstance(letters,list) and letters and all(isinstance(k,str) and k for k in letters)
                and not alphabet.intersection(letters), "each operation belongs to exactly one slot")
        alphabet.update(letters)
        slots.append(compile_choices(slot,gram))
    require(slots, "at least one slot required")
    slots.sort(key=len)
    require(set(case.get("positive_steps",[])) <= alphabet, "unknown positive step")
    adjoints = case.get("adjoints",{})
    require(set(adjoints) <= alphabet and all(adjoints.get(v)==k for k,v in adjoints.items()),
            "declared adjoints must be an involution on their letters")
    checks = []
    fit_words = set()
    for lane in ("fit","withheld"):
        rows = observations.get(lane,[])
        lane_words = set()
        for row in rows:
            word = word_tuple(row["word"],alphabet)
            require(word not in fit_words and word not in lane_words, "fit and withheld words must be distinct")
            lane_words.add(word)
            checks.append((lane,set(word),word,scalar(row["value"])))
        if lane == "fit":
            fit_words = lane_words
    identities = []
    for equation in case.get("identities",[]):
        terms = [(scalar(row["coefficient"]),word_tuple(row["word"],alphabet)) for row in equation["terms"]]
        require(terms, "empty operator identity")
        identities.append((set().union(*(set(w) for _,w in terms)), terms))
    step_names = set(case.get("positive_steps", []))
    cache = {}

    def admissible(name, op):
        key = (name,op)
        if key not in cache:
            form = gram*op
            cache[key] = positive(form) and positive(gram-form) if name in step_names else True
        return cache[key]

    counts["cartesian_bound"] += prod(map(len,slots))
    # Fit and operator checks run as soon as their last letter is assigned.
    # Withheld observations are replayed on complete survivors only.
    def visit(depth, identities_so_far, bindings, operators, applied):
        if depth == len(slots):
            counts["fit_survivors"] += 1
            for lane,_,word,value in checks:
                if lane == "withheld" and (weight*word_operator(operators,word,dimension)).trace() != value:
                    rejected["withheld_word"] += 1
                    return
            yield identities_so_far,operators,weight
            return
        for fragment,tags,ops in slots[depth]:
            if any(k in bindings and bindings[k]!=value for k,value in tags.items()):
                rejected["binding"] += 1
                continue
            if any(not admissible(k,op) for k,op in ops.items()):
                rejected["positive_step"] += 1
                continue
            merged = {**operators, **ops}
            if any(k in merged and v in merged and merged[k].H*gram != gram*merged[v] for k,v in adjoints.items()):
                rejected["adjoint"] += 1
                continue
            verified = set(applied)
            for i,(letters,terms) in enumerate(identities):
                if i not in applied and letters <= merged.keys():
                    defect = sum((c*word_operator(merged,w,dimension) for c,w in terms),sp.zeros(dimension))
                    if any(defect):
                        rejected["operator_identity"] += 1
                        break
                    verified.add(i)
            else:
                for i,(lane,letters,word,value) in enumerate(checks,len(identities)):
                    if lane == "fit" and i not in applied and letters <= merged.keys():
                        if (weight*word_operator(merged,word,dimension)).trace() != value:
                            rejected["fit_word"] += 1
                            break
                        verified.add(i)
                else:
                    yield from visit(depth+1,[*identities_so_far,fragment],{**bindings,**tags},merged,verified)

    yield from visit(0,[] if boundary is None else [{"boundary":boundary}],bindings,{},set())


def run(case, gram=None, weight=None, native_observations=None):
    """Join libraries on each declared carrier and classify their joint responses."""
    native = native_observations is not None
    if native:
        require(not ({"gram","trace_weight","boundary_choices","fit","withheld"} & set(case)),
                "original search receives its metric, contraction and observations from the native source")
        boundaries = [(None,gram,weight,{})]
    else:
        require("boundary_choices" not in case or not ({"gram","trace_weight"} & set(case)),
                "declare one boundary or a list of boundary alternatives")
        choices = case["boundary_choices"] if "boundary_choices" in case else [case]
        require(isinstance(choices,list) and choices, "nonempty boundary choices required")
        boundaries = [(row.get("id"),matrix(row["gram"]),matrix(row["trace_weight"]),row.get("bindings",{}))
                      for row in choices]
        require(len({g.rows for _,g,_,_ in boundaries})==1, "boundary alternatives must share coordinate size")
        if "boundary_choices" in case:
            require(all(isinstance(id,str) and id and isinstance(tags,dict) and
                        all(isinstance(k,str) and isinstance(v,str) for k,v in tags.items())
                        for id,_,_,tags in boundaries) and len({row[0] for row in boundaries})==len(boundaries),
                    "distinct boundary IDs and string bindings required")

    counts = dict(cartesian_bound=0,fit_survivors=0,withheld_survivors=0)
    rejected = Counter()
    digest = sha256()

    def candidates():
        for boundary,g,w,tags in boundaries:
            for identity,operators,contraction in compose(case,g,w,native_observations if native else case,
                                                           boundary,tags,counts,rejected):
                counts["withheld_survivors"] += 1
                digest.update(json.dumps(identity,separators=(",",":")).encode()+b"\n")
                yield identity,operators,contraction

    budget = case.get("equivalence_budget",1000)
    require(type(budget) is int and budget>=0, "nonnegative equivalence budget required")
    classes, next_probe, checked_pairs, complete = classify(candidates(),budget)
    status = "EXACT_COMPOSED_SEARCH" if counts["withheld_survivors"] else "NO_MATCH_IN_DECLARED_COMPOSITION_CLASS"
    return {"status":status,
            "input_role": "supplied_native_contractions" if native else "declared_finite_control",
            "physical_source_identity_proved":False,
            **counts,
            "pruned_prefixes":dict(rejected), "survivor_stream_sha256":digest.hexdigest(),
            "equivalence_complete":complete, "equivalence_checks":checked_pairs,
            "all_word_classes":[{"representative":name,"size":members} for name,_,members in classes],
            "unique_operational_class":complete and len(classes)==1,
            "next_distinguishing_word":next_probe}


def square_control():
    """Reuse the D4 sectors, composing time with a generated quarter and odd field."""
    def rows(value):
        return [[str(x) for x in row] for row in value.tolist()]
    one = sp.eye(4)
    stationary = sp.ones(4)/4
    sign = sp.Matrix([1,-1,1,-1])
    alternating = sign*sign.T/4
    doublet = one-stationary-alternating
    time = [{"id":[a,b],"operators":{"T":rows(stationary+a*alternating/8+b*doublet/8)}}
            for a,b in product(range(9),repeat=2)]
    quarter = rows(one[:, [1,2,3,0]])
    field = [{"id":sign,"operators":{"Y":rows(sp.diag(1,sign,0,-sign))}} for sign in (1,-1)]
    return {"gram":rows(one),"trace_weight":rows(sp.diag(1,0,0,0)),
            "slots":[{"letters":["T"],"choices":time},
                     {"letters":["U"],"grammar":{"atoms":{"quarter":quarter},"size":4,
                                                   "operations":["compose","mirror","iterate"]}},
                     {"letters":["Y"],"choices":field}],
            "positive_steps":["T"],"adjoints":{"Y":"Y"},
            "identities":[{"terms":[{"coefficient":1,"word":["U","T"]},
                                    {"coefficient":-1,"word":["T","U"]}]}],
            "fit":[{"word":["T"],"value":"1/2"},{"word":["U"],"value":0},
                   {"word":["U","U"],"value":0}],
            "withheld":[{"word":["U","T"],"value":"1/8"},
                        {"word":["U","Y","U","U","U"],"value":-1}]}
