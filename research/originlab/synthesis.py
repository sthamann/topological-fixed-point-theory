"""SyGuS learns local rules from exact observations; egglog simplifies compositions."""
from itertools import product
import cvc5
from cvc5 import Kind
from egglog import EGraph, Expr, StringLike, eq, rewrite, vars_


class Process(Expr):
    @classmethod
    def atom(cls, name: StringLike) -> "Process": ...
    def then(self, next: "Process") -> "Process": ...
    def adjoint(self) -> "Process": ...


def equivalent_examples():
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
    return {"status":"EXACT_TYPED_IDENTITIES",
            "equalities":[[str(x),str(y)] for x,y in examples],
            "typing":"finite composable linear maps; adjoint reverses boundary direction"}


def synthesis_solver(observations, exclude):
    tm = cvc5.TermManager()
    solver = cvc5.Solver(tm)
    solver.setOption("sygus","true")
    solver.setOption("tlimit-per","3000")
    solver.setLogic("ALL")
    sort = tm.getBooleanSort()
    inputs = [tm.mkVar(sort,x) for x in ("left","center","right")]
    nonterminal = tm.mkVar(sort,"B")
    grammar = solver.mkGrammar(inputs,[nonterminal])
    grammar.addRules(nonterminal,[*inputs,tm.mkBoolean(False),tm.mkBoolean(True),
                                  tm.mkTerm(Kind.NOT,nonterminal),
                                  *(tm.mkTerm(k,nonterminal,nonterminal)
                                    for k in (Kind.AND,Kind.OR,Kind.XOR))])
    rule = solver.synthFun("local",inputs,sort,grammar)

    def app(bits):
        return tm.mkTerm(Kind.APPLY_UF,rule,*(tm.mkBoolean(b) for b in bits))

    for bits,value in observations.items():
        solver.addSygusConstraint(tm.mkTerm(Kind.EQUAL,app(bits),tm.mkBoolean(value)))
    for table in exclude:
        solver.addSygusConstraint(tm.mkTerm(Kind.OR,*[
            tm.mkTerm(Kind.EQUAL,app(bits),tm.mkBoolean(not output))
            for bits,output in table.items()]))
    result = solver.checkSynth()
    if not result.hasSolution():
        return ("EXHAUSTED" if result.hasNoSolution() else "UNKNOWN"),None,None
    term = solver.getSynthSolution(rule)
    reducer = cvc5.Solver(tm)
    table = {}
    for bits in product((False,True),repeat=3):
        body = term[1].substitute(list(term[0]),[tm.mkBoolean(b) for b in bits])
        table[bits] = reducer.simplify(body).getBooleanValue()
    if any(table[b] != y for b,y in observations.items()) or table in exclude:
        raise ArithmeticError("synthesized truth table did not replay")
    return "SOLUTION",str(term),table


def synthesize(observations=None, limit=8, split_depth=0):
    # The default calibration supplies fixed uniform configurations only.
    # This admits many laws and introduces no TFPT target or spatial scale.
    observations = observations or {(False,)*3:False,(True,)*3:True}
    tables,expressions = [],[]
    for _ in range(limit):
        status,term,table = synthesis_solver(observations,tables)
        if status != "SOLUTION":
            break
        tables.append(table)
        expressions.append(term)
    else:
        status = "PARTIAL_ENUMERATION"
    questions = []
    for bits in product((False,True),repeat=3):
        if bits in observations:
            continue
        values = [t[bits] for t in tables]
        gain = sum(values)*(len(values)-sum(values))
        if gain:
            questions.append({"input":list(map(int,bits)),"outputs":list(map(int,values)),"pairs_separated":gain})
    questions.sort(key=lambda q:(-q["pairs_separated"],q["input"]))
    rows = []
    for term,table in zip(expressions,tables):
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
        rows.append({"rule":term,"truth_table":list(map(int,table.values())),"histories":paths})
    branches = []
    if questions and split_depth:
        # Both answers remain live; choosing one would add an empirical input.
        probe = tuple(map(bool,questions[0]["input"]))
        for answer in (False,True):
            branch = synthesize({**observations,probe:answer},limit,split_depth-1)
            branches.append({"probe":list(map(int,probe)),"hypothetical_answer":int(answer),
                             "status":branch["status"],"rules":branch["rules"],
                             "questions":branch["questions"],"branches":branch["branches"]})
    return {"status":status,"solver":"cvc5 SyGuS","rules":rows,"questions":questions,"branches":branches,
            "observations":[[list(map(int,b)),int(v)] for b,v in observations.items()],
            "scope":"uniform binary radius-one local rules; discrete history index, no physical time"}
