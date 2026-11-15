"""Decide all binary radius-one tables and normalize typed compositions."""
from itertools import islice, product

try:
    from egglog import EGraph, Expr, StringLike, eq, rewrite, vars_
except ModuleNotFoundError as error:
    if error.name!="egglog":
        raise
    EGraph = None


if EGraph is not None:
    class Process(Expr):
        @classmethod
        def atom(cls, name: StringLike) -> "Process": ...
        def then(self, next: "Process") -> "Process": ...
        def adjoint(self) -> "Process": ...


def _normal_word(expression):
    if isinstance(expression,str):
        return ((expression,False),)
    if expression[0]=="then":
        return _normal_word(expression[1])+_normal_word(expression[2])
    return tuple((letter,not adjoint) for letter,adjoint in reversed(_normal_word(expression[1])))


def equivalent_examples():
    if EGraph is None:
        examples = [
            (("then",("then","x","y"),"z"),("then","x",("then","y","z"))),
            (("adjoint",("then","x","y")),("then",("adjoint","y"),("adjoint","x"))),
            (("adjoint",("adjoint","x")),"x"),
        ]
        if any(_normal_word(left)!=_normal_word(right) for left,right in examples):
            raise ArithmeticError("typed word normal forms disagree")
        return {"status":"EXACT_TYPED_IDENTITIES","backend":"typed word normalization",
                "equalities":[[str(x),str(y)] for x,y in examples],
                "typing":"finite composable linear maps; adjoint reverses boundary direction"}

    a,b,c = vars_("a b c",Process)
    g = EGraph()
    x,y,z = map(Process.atom,"xyz")
    examples = [
        (x.then(y).then(z), x.then(y.then(z))),
        (x.then(y).adjoint(),y.adjoint().then(x.adjoint())),
        (x.adjoint().adjoint(),x),
    ]
    g.register(*(a for pair in examples for a in pair))
    g.register(rewrite(a.then(b).then(c)).to(a.then(b.then(c))),
               rewrite(a.then(b).adjoint()).to(b.adjoint().then(a.adjoint())),
               rewrite(a.adjoint().adjoint()).to(a))
    g.run(8)
    for left,right in examples:
        g.check(eq(left).to(right))
    return {"status":"EXACT_TYPED_IDENTITIES","backend":"egglog",
            "equalities":[[str(x),str(y)] for x,y in examples],
            "typing":"finite composable linear maps; adjoint reverses boundary direction"}


def _expression(table):
    inputs = ("left","center","right")
    cubes = ["(and "+" ".join(name if bit else f"(not {name})" for name,bit in zip(inputs,bits))+")"
             for bits,output in table.items() if output]
    if len(cubes)<2:
        return cubes[0] if cubes else "false"
    return "(or "+" ".join(cubes)+")"


def synthesize(observations=None, limit=8, split_depth=0):
    # The default calibration supplies fixed uniform configurations only.
    # This admits many laws and introduces no TFPT target or spatial scale.
    if observations is None:
        observations = {(False,)*3:False,(True,)*3:True}
    inputs = list(product((False,True),repeat=3))
    if (type(limit) is not int or limit<0 or type(split_depth) is not int or split_depth<0 or
            any(bits not in inputs or value not in (False,True) for bits,value in observations.items())):
        raise ValueError("invalid binary observations, rule limit or split depth")
    unknown = [bits for bits in inputs if bits not in observations]
    count = 2**len(unknown)
    tables = []
    # Every table is in the given Boolean grammar by disjunctive normal form.
    # The finite set includes every completion of all currently measured bits.
    for values in islice(product((False,True),repeat=len(unknown)),min(limit,count)):
        table = {**observations,**dict(zip(unknown,values))}
        tables.append({key:table[key] for key in inputs})
    status = "EXHAUSTED" if len(tables)==count else "PARTIAL_ENUMERATION"
    questions = []
    for bits in inputs:
        if bits in observations:
            continue
        values = [t[bits] for t in tables]
        gain = sum(values)*(len(values)-sum(values))
        if gain:
            questions.append({"input":list(map(int,bits)),"outputs":list(map(int,values)),"pairs_separated":gain})
    questions.sort(key=lambda q:(-q["pairs_separated"],q["input"]))
    rows = []
    for table in tables:
        # Compare a closed periodic history against open zero boundary data.
        seed = (False,True,False,False,False)
        paths = {}
        for boundary in ("periodic","open_zero","reversed_glue"):
            state = seed
            history = [list(map(int,state))]
            for _ in range(5):
                left,right = ((state[-1],state[0]) if boundary=="periodic"
                              else (state[0],state[-1]) if boundary=="reversed_glue" else (False,False))
                extended = (left,*state,right)
                state = tuple(table[extended[j:j+3]] for j in range(len(state)))
                history.append(list(map(int,state)))
            paths[boundary] = history
        rows.append({"rule":_expression(table),"truth_table":list(map(int,table.values())),"histories":paths})
    branches = []
    if questions and split_depth:
        # Both answers remain live; choosing one would add an empirical input.
        probe = tuple(map(bool,questions[0]["input"]))
        for answer in (False,True):
            branch = synthesize({**observations,probe:answer},limit,split_depth-1)
            branches.append({"probe":list(map(int,probe)),"hypothetical_answer":int(answer),
                             "status":branch["status"],"rules":branch["rules"],
                             "questions":branch["questions"],"branches":branch["branches"]})
    return {"status":status,"solver":"exhaustive Boolean truth tables","eligible_count":count,
            "rules":rows,"questions":questions,"branches":branches,
            "observations":[[list(map(int,b)),int(v)] for b,v in observations.items()],
            "scope":"uniform binary radius-one local rules; discrete history index, no physical time"}
