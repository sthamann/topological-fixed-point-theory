"""Search native F1 candidates and check proofs on one common model."""

from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory

from evidence_gates import REQUIREMENTS
from proof_search import DEFAULT_TOOLCHAIN, IDENTIFIER, LeanKernel, read_contract, terms_at_cost


HERE = Path(__file__).resolve().parent
NATIVE_REGISTRY = HERE/"f1_native_contract.json"
GATES = tuple(REQUIREMENTS)
SAFE_TERM = re.compile(r"[A-Za-z0-9_.() ]+")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read_f1_contract(path, expected):
    content = path.read_bytes()
    if sha256(content).hexdigest() != expected:
        raise ValueError("F1 contract digest mismatch")
    contract = json.loads(content)
    if set(contract) != {"lean_contract","lean_contract_sha256","gates"}:
        raise ValueError("F1 contract must pin a Lean contract, and name its gate predicates")
    name = contract["lean_contract"]
    if not isinstance(name,str) or not re.fullmatch(r"[A-Za-z0-9_.-]+\.json",name):
        raise ValueError("use a sibling JSON Lean contract")
    target = path.parent/name
    if digest(target) != contract["lean_contract_sha256"]:
        raise ValueError("nested Lean contract digest mismatch")
    gates = contract["gates"]
    if not isinstance(gates,dict) or set(gates) != set(GATES):
        raise ValueError("the F1 contract requires exactly G0 through G7")
    predicates = [gates[gate] for gate in GATES]
    if any(not isinstance(item,str) or not IDENTIFIER.fullmatch(item) for item in predicates):
        raise ValueError("each gate needs a qualified Lean predicate")
    base, modules = read_contract(target)
    if digest(target) != contract["lean_contract_sha256"]:
        raise ValueError("nested Lean contract changed during loading")
    if len(set(predicates+[base["goal"]])) != len(GATES)+1:
        raise ValueError("use distinct root and gate predicates")
    if any(digest(file) != sha for file,sha in modules):
        raise ValueError("a declared Lean module digest differs")
    return base, modules, {gate:gates[gate] for gate in GATES}, target


def native_registry():
    return json.loads(NATIVE_REGISTRY.read_text())


def resolve_contract(path=None, expected=None):
    if (path is None) != (expected is None):
        raise ValueError("supply the F1 contract and its expected SHA256 together")
    native = native_registry()["contract"]
    if path is None and native is not None:
        return HERE/native["file"],native["sha256"],True
    # External data may name a formal goal, but only the repository's
    # independently reviewed and pinned contract confers native scope.
    is_native = (native is not None and path is not None and
                 path.resolve()==(HERE/native["file"]).resolve() and expected==native["sha256"])
    return path,expected,is_native


def contract_fingerprints(path, expected):
    path,expected,_ = resolve_contract(path,expected)
    if path is None:
        return {"native_registry":digest(NATIVE_REGISTRY)}
    _,modules,_,target = read_f1_contract(path,expected)
    return {"native_registry":digest(NATIVE_REGISTRY),"contract":expected,
            "lean_contract":digest(target),"modules":{file.name:sha for file,sha in modules}}


def check_bundle(kernel, module, root, gates, producer, proof):
    if (not isinstance(producer,str) or not SAFE_TERM.fullmatch(producer)
            or not isinstance(proof,str) or not SAFE_TERM.fullmatch(proof)):
        raise ValueError("use closed terms from the constructor grammar")
    if not IDENTIFIER.fullmatch(root) or set(gates)!=set(GATES) or any(
            not IDENTIFIER.fullmatch(gates[gate]) for gate in GATES):
        raise ValueError("invalid pinned F1 predicates")
    statement = " ∧ ".join(f"{predicate} originlab_model"
                           for predicate in [root,*(gates[g] for g in GATES)])
    return kernel.verify(module,statement,proof,f"def originlab_model := {producer}")


def scheduled_budgets(stage):
    """Dovetail syntax cost and finite check budgets across resumed stages."""
    for budget in range(stage-1):
        yield stage-budget,budget+1


def candidates(base, cost):
    for producer_cost in range(1,cost):
        for producer in terms_at_cost(base["producers"],producer_cost):
            for proof in terms_at_cost(base["proofs"],cost-producer_cost):
                yield producer,proof


def original_session_checks(root):
    from e8_family_return import controls as family_return
    from f1_binary_clifford_origin import controls as binary
    from f1_four_family_transport import controls as four_family_transport
    from f1_native_pair_selection import controls as native_pair
    from f1_spherical_pairing_selection import controls as sphere
    from f1_reflected_cap import ASSEMBLER, controls as reflected_cap
    from r113_pair_search import search as pair_search
    from relative_history import control as connected_history_control
    from r113_radial_order import RADIAL_DATA, controls as radial_order
    from r113_reflected_radial import controls as reflected_radial
    from source import NATIVE

    checks = {"binary":binary(),"spherical":sphere(),"native_pair":native_pair(root),
              "e8_family_return":family_return(),
              "four_family_transport":four_family_transport(),
              "relative_connected_history":connected_history_control()}
    if (root/ASSEMBLER).is_file() and (root/NATIVE).is_file():
        checks["reflected_cap"] = reflected_cap(root)
        checks["reflected_pair_search"] = pair_search(root)
    if (root/RADIAL_DATA).is_file():
        checks["conditional_radial_order"] = radial_order(root)
        checks["conditional_reflected_radial"] = reflected_radial(root)
    return checks


def search_native(root, domain, depth, hypotheses, toolchain, seconds, stage):
    from coherent_origin import search as search_coherent
    from p1_origin_hypotheses import search as search_p1
    from q0_rule_search import native_comparison, search as search_q0
    from source import evaluate_original
    from source_return import affine_return

    q0 = search_q0(depth)
    source = evaluate_original(root, domain)
    q0["native_comparison"] = native_comparison(q0, source)
    anchors = {}
    if source.get("reason") != "ORIGINAL_SOURCE_FILE_MISSING":
        anchors["affine_return"] = affine_return(root)
    return {"status":"NATIVE_SEARCH_BOUND_REACHED","physical_solution":False,
            "formal_acceptance_registered":False,
            "acceptance_configuration":native_registry().get("contract_status","NOT_CONFIGURED"),
            "candidate_contract_path":"F1-CANDIDATE-CONTRACT.md",
            "native_gate_sources":native_registry()["gates"],
            "original_source":source,"original_anchors":anchors,
            "q0_rule_search":q0,"p1_origins":search_p1(**(hypotheses or {})),
            "coherent_origin_search":search_coherent(stage),
            "coherent_origin_certificate":check_coherent(toolchain,seconds),
            "ward_certificate":check_ward(toolchain,seconds),
            "latest_session_controls":original_session_checks(root)}


def search(path=None, expected=None, stage=2, toolchain=DEFAULT_TOOLCHAIN, seconds=20, *,
           root=HERE.parents[1], domain=None, q0_depth=3, p1_hypotheses=None):
    if type(stage) is not int or stage<2 or type(seconds) is not int or seconds<1:
        raise ValueError("F1 search starts at stage two with a positive check budget")
    path,expected,is_native = resolve_contract(path,expected)
    if path is None:
        return search_native(root,domain,q0_depth,p1_hypotheses,toolchain,seconds,stage)
    base,modules,gates,_ = read_f1_contract(path,expected)
    limited = []
    attempts = 0
    try:
        with TemporaryDirectory(prefix="originlab-f1-") as temp:
            kernel = LeanKernel(Path(temp),toolchain,seconds*(stage-1))
            kernel.prepare(modules)
            for cost,budget in scheduled_budgets(stage):
                kernel.timeout = budget*seconds
                for producer,proof in candidates(base,cost):
                    attempts += 1
                    try:
                        checked = check_bundle(kernel,modules[-1][0].stem,base["goal"],gates,producer,proof)
                    except subprocess.TimeoutExpired:
                        limited.append([producer,proof,budget])
                        continue
                    if checked["status"] == "INCONCLUSIVE":
                        limited.append([producer,proof,budget])
                    if checked["status"] == "PROVED":
                        return {"status":"F1_PROVED" if is_native else "F1_CONTRACT_PROVED",
                                "physical_solution":is_native,"contract_sha256":expected,
                                "modules":{file.name:sha for file,sha in modules},
                                "producer":producer,"same_model_gates":gates,
                                "proof":checked,"stage":stage,"attempts":attempts}
    except (OSError,subprocess.TimeoutExpired) as error:
        return {"status":"PROOF_BACKEND_INCONCLUSIVE","physical_solution":False,
                "stage":stage,"next_stage":stage+1,"attempts":attempts,"error":str(error)}
    return {"status":"F1_CONTRACT_OPEN","physical_solution":False,"stage":stage,
            "next_stage":stage+1,"attempts":attempts,"budget_limited":limited}


def check_ward(toolchain, seconds):
    ward = HERE/"OriginF1Ward.lean"
    try:
        with TemporaryDirectory(prefix="originlab-ward-") as temp:
            kernel = LeanKernel(Path(temp),toolchain,max(seconds,45))
            kernel.prepare(((ward,digest(ward)),))
            statement = ("(∀ (f : Nat → Rat), OriginF1Ward.Ward f → "
                         "∀ n, f n = OriginF1Ward.spectrum (f 0) n) ∧ "
                         "(∃ f g : Nat → Rat, OriginF1Ward.Ward f ∧ OriginF1Ward.Ward g ∧ f 0 ≠ g 0)")
            return kernel.verify(ward.stem,statement,
                "⟨OriginF1Ward.wardDeterminesRatios, OriginF1Ward.wardLeavesScale⟩")
    except (OSError,subprocess.TimeoutExpired) as error:
        return {"status":"PROOF_BACKEND_INCONCLUSIVE","error":str(error)}


def check_coherent(toolchain, seconds):
    model = HERE/"OriginCoherent.lean"
    statement = (
        "(∀ (f : Nat → Int), (∀ n, f n + f 1 = f (n + 1) + f 0) → "
        "∀ n, f n = f 0 + (n : Int) * (f 1 - f 0)) ∧ "
        "(∀ (p : Nat → Rat), p 0 ≠ 0 → (∀ n, p n * p 1 = p (n+1) * p 0) → "
        "∀ n, p n = p 0 * (p 1 / p 0)^n) ∧ "
        "(∀ s : OriginCoherent.PairState, "
        "OriginCoherent.swap (OriginCoherent.swap s) = s ∧ "
        "OriginCoherent.marked (OriginCoherent.marked s) = s ∧ "
        "OriginCoherent.swap (OriginCoherent.marked s) = OriginCoherent.ordered s ∧ "
        "OriginCoherent.swap (OriginCoherent.marked "
        "(OriginCoherent.swap (OriginCoherent.marked s))) = OriginCoherent.minus s)"
    )
    proof = ("⟨OriginCoherent.affineGenerator, OriginCoherent.geometricWeights, "
             "OriginCoherent.markerIdentities⟩")
    try:
        with TemporaryDirectory(prefix="originlab-coherent-proof-") as temp:
            kernel = LeanKernel(Path(temp),toolchain,max(seconds,45))
            kernel.prepare(((model,digest(model)),))
            return kernel.verify(model.stem,statement,proof)
    except (OSError,subprocess.TimeoutExpired) as error:
        return {"status":"PROOF_BACKEND_INCONCLUSIVE","error":str(error)}
