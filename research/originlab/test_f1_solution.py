import json
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory

import pytest

from evidence_gates import assess
import f1_solution
from f1_solution import GATES, check_bundle, digest, read_f1_contract, scheduled_budgets
from proof_search import DEFAULT_TOOLCHAIN, LeanKernel


def fixture(directory):
    lean = directory/"Fixture.lean"
    statements = ["import Std","namespace Fixture","def Goal (n : Nat) : Prop := n = 0",
                  "def G2 (n : Nat) : Prop := n*n = n"]
    statements += [f"def {gate} (n : Nat) : Prop := n + {i} = {i}"
                   for i,gate in enumerate(GATES) if gate!="G2"]
    all_at_zero = ["Goal 0",*(f"{gate} 0" for gate in GATES)]
    mixed = [claim if not claim.startswith("G2") else "G2 1" for claim in all_at_zero]
    for name, claims in (("complete",all_at_zero),("mixed",mixed)):
        statements.append(f"theorem {name} : {' ∧ '.join(claims)} := ⟨{', '.join(['rfl']*len(claims))}⟩")
    statements += ["theorem rootOnly : Goal 0 := rfl","end Fixture"]
    lean.write_text("\n".join(statements)+"\n")
    nested = directory/"lean-contract.json"
    nested.write_text(json.dumps({"modules":[{"file":lean.name,"sha256":digest(lean)}],
        "goal":"Fixture.Goal","producers":{"Nat.zero":0},"proofs":{"Fixture.complete":0}}))
    contract = directory/"contract.json"
    rows = {"lean_contract":nested.name,"lean_contract_sha256":digest(nested),
            "gates":{gate:f"Fixture.{gate}" for gate in GATES}}
    contract.write_text(json.dumps(rows))
    return contract,rows


def test_complete_gate_manifest_is_pinned(tmp_path):
    path,rows = fixture(tmp_path)
    pin = digest(path)
    _,_,gates,_ = read_f1_contract(path,pin)
    assert tuple(gates)==GATES
    rows["gates"].pop("G6")
    path.write_text(json.dumps(rows))
    with pytest.raises(ValueError,match="digest mismatch"):
        read_f1_contract(path,pin)
    with pytest.raises(ValueError,match="requires exactly"):
        read_f1_contract(path,digest(path))


def test_native_search_runs_without_registered_formal_acceptance(monkeypatch):
    monkeypatch.setattr(f1_solution,"native_registry",lambda: {"contract":None,"gates":{}})
    monkeypatch.setattr(f1_solution,"check_ward",lambda *_: {"status":"NOT_RUN"})
    monkeypatch.setattr(f1_solution,"original_session_checks",lambda _: {})
    result = f1_solution.search(q0_depth=1,p1_hypotheses={"hypotheses":["protected_defect"]})
    assert result["status"] == "NATIVE_SEARCH_BOUND_REACHED"
    assert not result["physical_solution"]
    assert result["q0_rule_search"]["word_count"] > 0
    assert result["q0_rule_search"]["native_comparison"]["status"] == "NO_SEPARATOR_AT_THIS_DEPTH"
    assert result["original_anchors"]["affine_return"]["adjoint_ordered_return"] == "1"
    assert set(result["p1_origins"]["hypotheses"]) == {"protected_defect"}


def test_scheduler_revisits_each_syntax_cost_with_growing_budgets():
    assert list(scheduled_budgets(2))==[(2,1)]
    assert list(scheduled_budgets(5))==[(5,1),(4,2),(3,3),(2,4)]
    assert all((cost, budget) in scheduled_budgets(cost+budget-1)
               for cost in range(2,10) for budget in range(1,10))


def test_coherent_search_joins_real_contractions_and_keeps_shells_disjoint():
    from coherent_origin import search
    stage2,stage3=search(2),search(3)
    assert stage2["candidate_count"]==2 and stage3["candidate_count"]==8
    key=lambda c: (c["q"],c["empty_step"],c["reference_step"],c["orientation"])
    assert not {key(c) for c in stage2["candidates"]} & {key(c) for c in stage3["candidates"]}
    forward,backward=stage2["candidates"]
    assert forward["ABAB"]==backward["ABAB"]=="-1"
    assert forward["connected_UY"]=="I/2" and backward["connected_UY"]=="-I/2"
    assert [c["same_history_defect"]["selected_rank"] for c in (forward,backward)]==[0,1]
    assert not stage2["physical_solution"]


def test_missing_original_source_blocks_only_its_own_route():
    report = assess({"status":"BLOCKED","missing":["original_blocked_fiber_packet"]})
    assert report["gates"]["G3"]["status"] == "BLOCKED"
    assert report["gates"]["G3"]["route"] == "original_source"
    assert report["gates"]["G4"]["status"] == "NOT_TESTED"
    assert report["gates"]["G4"]["route_dependencies"]["original_source"]["blocked_by"] == ["G3"]
    assert all(row["status"] == "NOT_TESTED" for gate,row in report["gates"].items() if gate != "G3")


def test_checked_candidate_contract_reaches_candidate_pass_without_promotion():
    complete = {"status":"F1_CONTRACT_PROVED","physical_solution":False,"contract_sha256":"c"*64,
                "same_model_gates":{gate:f"Candidate.{gate}" for gate in GATES},
                "proof":{"status":"PROVED","certificate_sha256":"d"*64}}
    report = assess({"status":"BLOCKED"},complete=complete)
    assert report["physical_status"] == "CANDIDATE_PROVED"
    assert not report["physical_solution"]
    assert all(row["status"] == "CANDIDATE_PASS" for row in report["gates"].values())
    incomplete = {**complete,"same_model_gates":{"G0":"Candidate.G0"}}
    assert assess({},complete=incomplete)["physical_status"] == "OPEN"


def test_proof_contract_cannot_promote_external_or_incomplete_receipt():
    assert not assess({},complete={"status":"F1_CONTRACT_PROVED","physical_solution":True})["physical_solution"]
    assert not assess({},complete={"status":"F1_PROVED","physical_solution":True,
                                  "same_model_gates":{"G0":"G0"}})["physical_solution"]


def test_kernel_requires_root_and_all_gates_on_same_model(tmp_path,monkeypatch):
    if shutil.which("elan") is None or DEFAULT_TOOLCHAIN not in subprocess.check_output(
            ["elan","toolchain","list"],text=True):
        pytest.skip("the pinned Lean toolchain is not installed")
    manifest,_ = fixture(tmp_path)
    base,modules,gates,_ = read_f1_contract(manifest,digest(manifest))
    results = {}
    with TemporaryDirectory() as compiled:
        kernel = LeanKernel(Path(compiled),DEFAULT_TOOLCHAIN,45)
        kernel.prepare(modules)
        for name in ("complete","rootOnly","mixed"):
            results[name] = check_bundle(kernel,modules[-1][0].stem,base["goal"],
                                         gates,"0",f"Fixture.{name}")["status"]
    assert results=={"complete":"PROVED","rootOnly":"REJECTED","mixed":"REJECTED"}
    registry = tmp_path/"native.json"
    registry.write_text(json.dumps({"contract":{"file":manifest.name,"sha256":digest(manifest)}}))
    monkeypatch.setattr(f1_solution,"HERE",tmp_path)
    monkeypatch.setattr(f1_solution,"NATIVE_REGISTRY",registry)
    accepted = f1_solution.search(seconds=45)
    assert accepted["status"] == "F1_PROVED"
    bound = assess({"status":"BLOCKED"},complete=accepted)
    assert bound["physical_solution"] and all(row["status"]=="PASS" for row in bound["gates"].values())
