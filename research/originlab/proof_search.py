"""Search generative source rules and submit proofs to a fixed Lean goal."""

from contextlib import contextmanager
from hashlib import sha256
import json
from math import factorial, prod
import os
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory

import sympy as sp


HERE = Path(__file__).resolve().parent
FORMAL_FILES = ("OriginProofGoal.lean", "OriginProofRules.lean")
CONTROL_GOAL = "75dcdacdae9d4683a3b545f411a237a55e81ba82408a413429675e176b702a97"
DEFAULT_TOOLCHAIN = "leanprover/lean4:v4.29.1"
STANDARD_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}


class GaussianProducer:
    """Rank-one Fock source with one common level and arbitrary mode cutoff."""

    def __init__(self, level):
        if type(level) is not int or level < 1:
            raise ValueError("oscillator level must be a positive integer")
        self.level = level

    def gram(self, left, right, offset=0):
        if len(left) != len(right):
            raise ValueError("use the same ordered mode basis on both ends")
        if left != right:
            return 0
        return prod((self.level * mode)**n * factorial(n)
                    for mode, n in enumerate(left, offset+1))

    def endpoint(self, left, right, windows, insertions):
        """Apply ordered windows and words from the right endpoint to the left."""
        if len(windows) != len(insertions)+1:
            raise ValueError("one window on each side of every insertion")
        if any(type(n) is not int or n < 0 for word in (left, right) for n in word):
            raise ValueError("occupations must be nonnegative integers")
        if len(left) != len(right):
            raise ValueError("use the same ordered mode basis on both ends")
        occupied = list(right)
        weight = sp.S.One
        for index in range(len(windows)-1, -1, -1):
            energy = sum(mode*n for mode,n in enumerate(occupied, 1))
            weight *= windows[index]**energy
            if index == 0:
                break
            for mode, delta in reversed(insertions[index-1]):
                if type(mode) is not int or not 1 <= mode <= len(occupied) or delta not in (-1, 1):
                    raise ValueError("invalid oscillator operation")
                slot = mode-1
                if delta == -1:
                    if occupied[slot] == 0:
                        return sp.S.Zero
                    weight *= self.level * mode * occupied[slot]
                occupied[slot] += delta
        return weight * self.gram(tuple(left), tuple(occupied))

    def forms(self, left, right, windows, u, y):
        return {key:self.endpoint(left,right,windows,words) for key,words in
                (("Q0",((),())),("QU",(u,())),("QY",((),y)),("QUY",(u,y)))}


def audit_axioms(output):
    match = re.search(r"'originlab_certificate' depends on axioms: \[([^\]]*)\]", output)
    if match:
        return sorted(name.strip() for name in match[1].split(",") if name.strip())
    if "'originlab_certificate' does not depend on any axioms" in output:
        return []
    return None


class LeanKernel:
    """Compile trusted modules, then check an exact proposition and its proof."""

    def __init__(self, directory, toolchain, timeout):
        self.directory = directory
        self.env = {**os.environ, "LEAN_PATH": str(directory)}
        self.toolchain = toolchain
        self.timeout = timeout
        self.modules = []

    def invoke(self, executable, *args):
        return subprocess.run(["elan","run",self.toolchain,executable,*args],
                              cwd=self.directory, env=self.env, timeout=self.timeout,
                              text=True, capture_output=True, check=False)

    def prepare(self, modules):
        for file, expected in modules:
            contents = file.read_bytes()
            if sha256(contents).hexdigest() != expected:
                raise ValueError(f"frozen module digest mismatch: {file.name}")
            (self.directory/file.name).write_bytes(contents)
            run = self.invoke("lean","-o",f"{file.stem}.olean",file.name)
            if run.returncode:
                raise ValueError(f"trusted Lean module failed: {file.name}\n{run.stdout}{run.stderr}")
            self.modules.append(file.stem)

    def verify(self, module, statement, term, preamble="", replay=True):
        # These parameters come from a registered trusted adapter. Raw
        # model text is never executed by the autonomous queue.
        body = (f"import {module}\n{preamble}\n"
                f"theorem originlab_certificate : {statement} := {term}\n"
                "#print axioms originlab_certificate\n")
        source_path = self.directory/"ProofCandidate.lean"
        module_path = source_path.with_suffix(".olean")
        module_path.unlink(missing_ok=True)
        source_path.write_text(body)
        checked = self.invoke("lean", "-o", module_path.name, source_path.name)
        output = checked.stdout + checked.stderr
        axioms = audit_axioms(output)
        receipt = {"axioms":axioms, "toolchain":self.toolchain,
                   "statement":statement, "certificate":body,
                   "certificate_sha256":sha256(body.encode()).hexdigest()}
        if checked.returncode or axioms is None or not set(axioms) <= STANDARD_AXIOMS:
            return {**receipt,"status":"REJECTED","diagnostic":output.strip()}
        if not replay:
            return {**receipt,"status":"KERNEL_CHECKED"}
        replayed = self.invoke("leanchecker","-v",*self.modules,"ProofCandidate")
        if replayed.returncode:
            return {**receipt,"status":"INCONCLUSIVE","diagnostic":(replayed.stdout+replayed.stderr).strip()}
        return {**receipt,"status":"PROVED","replay":replayed.stdout.strip()}


@contextmanager
def control_kernel(toolchain, timeout=55):
    with TemporaryDirectory(prefix="originlab-proof-") as temp:
        kernel = LeanKernel(Path(temp),toolchain,timeout)
        goal, rules = (HERE/name for name in FORMAL_FILES)
        kernel.prepare(((goal,CONTROL_GOAL),(rules,sha256(rules.read_bytes()).hexdigest())))
        yield kernel


def candidate_pairs(cost):
    """A finite cost shell of positive levels and registered proof strategies."""
    for tactic_cost, strategy in enumerate(("finite_decision", "normal_order"),1):
        level = cost-tactic_cost
        if level > 0:
            yield level,strategy


def affine_forecast():
    """Generate a charged two-window word in the existing E8 radial chart."""
    from affine_operators import (LatticeFock, MARKED_GAMMA, MARKED_P1, MARKED_P2,
                                  dot, integral, lattice_frame, marked_word)

    frame, _ = lattice_frame()
    fock = LatticeFock([[dot(a,b) for b in frame] for a in frame])
    p1,p2,gamma = (integral(frame,r) for r in (MARKED_P1,MARKED_P2,MARKED_GAMMA))
    u,v = sp.symbols("u v",positive=True)

    def window(state, q):
        return {key:sp.sympify(weight)*q**(fock.pair(key[0],key[0])//2
                                        + sum(mode for mode,_ in key[1]))
                for key,weight in state.items()}

    m,n = 4,6
    right = fock.vertex(p2,-n,fock.vacuum())
    left = fock.vertex(p1,-m,fock.vacuum())
    q0 = sp.expand(fock.inner(right,window(window(right,u),v)))
    marked = sp.expand(fock.inner(left,window(fock.vertex(gamma,n-m,window(right,u)),v)))
    vacuum_word = marked_word(fock,frame,n)[m,n]
    # Our root-pairing convention uses the actual contravariant norm.
    # The prior vacuum word uses the raw opposite root instead of its adjoint.
    raw_adjoint_factor = fock.sign(p1,tuple(-i for i in p1))
    if q0 != n*(u*v)**n or marked != vacuum_word*raw_adjoint_factor*u**n*v**m:
        raise ArithmeticError("affine forecast does not replay the graded charge word")
    return {"status":"GENERATED_AFFINE_RADIAL_COMPARISON","mode_pair":[m,n],
            "empty_charged_form":str(q0),"ordered_charged_form":str(marked),
            "absolute_r113_reference":None,"native_source_identified":False}


def run(first_cost=2, shells=3, toolchain=DEFAULT_TOOLCHAIN, timeout=55):
    if type(first_cost) is not int or type(shells) is not int or first_cost < 2 or shells < 1:
        raise ValueError("positive proof slices start at cost two")
    attempts = []
    checked = None
    try:
        with control_kernel(toolchain,timeout) as kernel:
            for cost in range(first_cost,first_cost+shells):
                for level,strategy in candidate_pairs(cost):
                    term = ("by decide" if strategy=="finite_decision" else
                            f"closeGaussian {level} (by decide)")
                    receipt = kernel.verify("OriginProofRules",
                        f"Goal (gaussian {level})", term, "open OriginLabProof", replay=False)
                    attempts.append({"cost":cost,"level":level,"strategy":strategy,
                                     "status":receipt["status"]})
                    if receipt["status"] != "KERNEL_CHECKED":
                        continue
                    checked = kernel.verify("OriginProofRules",f"Goal (gaussian {level})",
                                            term,"open OriginLabProof")
                    if checked["status"] == "PROVED":
                        negative = kernel.verify("OriginProofRules","¬ FreeLevelUnique",
                                                 "refuteFreeLevel","open OriginLabProof")
                        controls = [
                            ("injected_axiom", "forgedGoal", "axiom forgedGoal : Goal (gaussian 1)"),
                            ("finite_test", "(by decide : gaussian 1 0 [(1,1)] = 1)", ""),
                            ("table_without_rules", "gaussianRules 1", ""),
                            ("assumed_goal", "fun h : Goal (gaussian 1) => h", ""),
                        ]
                        rejected = {}
                        for name, proof, extra in controls:
                            check = kernel.verify("OriginProofRules","Goal (gaussian 1)",proof,
                                                  f"open OriginLabProof\n{extra}",replay=False)
                            rejected[name] = {"status":check["status"],"axioms":check["axioms"]}
                        if any(case["status"] != "REJECTED" for case in rejected.values()):
                            raise ArithmeticError("invalid proof respected as a complete goal")
                        producer = GaussianProducer(level)
                        windows = tuple(map(sp.Rational,("1/2","2/3","3/5")))
                        occupations = (2,1)
                        held_out = producer.forms(occupations,occupations,windows,
                                                  ((1,-1),),((1,1),))
                        return {"status":"CONTROL_SOLVED","scope":"normal-ordered rank-one oscillator goal",
                                "goal_sha256":CONTROL_GOAL,"producer":{"level":level},
                                "proof":checked,"free_level_refutation":negative,
                                "attempts":attempts,"rejected_certificates":rejected,
                                "held_out_forms":{k:str(v) for k,v in held_out.items()},
                                "affine_forecast":affine_forecast(),
                                "physical_solution":False}
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"status":"PROOF_BACKEND_INCONCLUSIVE","error":str(error),
                "attempts":attempts,"physical_solution":False}
    return {"status":"OPEN_IN_SEARCH_SLICE","attempts":attempts,
            "next_cost":first_cost+shells,"last_kernel_replay":checked,
            "physical_solution":False}


IDENTIFIER = re.compile(r"[A-Za-z][A-Za-z0-9_]*(?:\.[A-Za-z][A-Za-z0-9_]*)*")


def read_contract(path):
    """A trusted, externally supplied list of pinned Lean modules and constructor names."""
    contract = json.loads(path.read_text())
    if set(contract) != {"modules","goal","producers","proofs"}:
        raise ValueError("proof contract requires modules, goal, producers and proofs")
    if not isinstance(contract["goal"],str) or not IDENTIFIER.fullmatch(contract["goal"]):
        raise ValueError("use the qualified name of the fixed goal predicate")
    modules = []
    for row in contract["modules"]:
        name = row["file"]
        digest = row["sha256"]
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*\.lean",name) or not re.fullmatch("[0-9a-f]{64}",digest):
            raise ValueError("contract Lean modules must have a simple name and a pinned digest")
        modules.append((path.parent/name,digest))
    if not modules or len({file.name for file,_ in modules}) != len(modules):
        raise ValueError("contract needs distinct Lean modules in dependency order")
    for kind in ("producers","proofs"):
        if not contract[kind]:
            raise ValueError("each term grammar needs at least one constructor")
        for name, arity in contract[kind].items():
            if not IDENTIFIER.fullmatch(name) or type(arity) is not int or not 0 <= arity <= 8:
                raise ValueError("use named constructors with declared arity zero through eight")
    return contract,modules


def terms_at_cost(constructors, cost):
    """Enumerate application trees by node cost, with numerals costing one plus value."""
    if cost < 1:
        return
    yield str(cost-1)
    for name, arity in constructors.items():
        for arguments in arguments_at_cost(constructors,arity,cost-1):
            yield name if not arguments else f"({name} {' '.join(arguments)})"


def arguments_at_cost(constructors, count, cost):
    if count == 0:
        if cost == 0:
            yield ()
        return
    for first in range(1,cost-count+2):
        for term in terms_at_cost(constructors,first):
            for rest in arguments_at_cost(constructors,count-1,cost-first):
                yield term,*rest


def search_contract(path, first_cost=2, shells=2, toolchain=DEFAULT_TOOLCHAIN, timeout=55):
    if type(first_cost) is not int or type(shells) is not int or first_cost < 2 or shells < 1:
        raise ValueError("positive proof slices start at cost two")
    contract, modules = read_contract(path)
    attempts = 0
    limits = []
    try:
        with TemporaryDirectory(prefix="originlab-proof-") as temp:
            kernel = LeanKernel(Path(temp),toolchain,timeout)
            kernel.prepare(modules)
            for cost in range(first_cost,first_cost+shells):
                for producer_cost in range(1,cost):
                    for producer in terms_at_cost(contract["producers"],producer_cost):
                        for proof in terms_at_cost(contract["proofs"],cost-producer_cost):
                            attempts += 1
                            try:
                                receipt = kernel.verify(modules[-1][0].stem,
                                    f"{contract['goal']} ({producer})",proof)
                            except subprocess.TimeoutExpired:
                                limits.append({"cost":cost,"producer":producer,"proof":proof})
                                continue
                            if receipt["status"] == "PROVED":
                                return {"status":"FORMAL_GOAL_PROVED","producer":producer,
                                        "proof":receipt,"attempts":attempts,"limited_pairs":limits,
                                        "trusted_modules":{file.name:digest for file,digest in modules},
                                        "physical_solution":False}
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"status":"PROOF_BACKEND_INCONCLUSIVE","error":str(error),
                "attempts":attempts,"physical_solution":False}
    return {"status":"OPEN_IN_SEARCH_SLICE","attempts":attempts,"limited_pairs":limits,
            "next_cost":first_cost+shells,"physical_solution":False}
