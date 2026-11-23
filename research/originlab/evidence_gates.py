"""Report the evidence actually available for the plan's physical gates.

Gate statuses have separate meanings:

- PASS            a checked native proof covers the gate on the same producer.
- CANDIDATE_PASS  a checked proof of a supplied candidate contract covers the gate on
                  the same producer; physical acceptance still needs the contract to be
                  reviewed and registered in f1_native_contract.json.
- PARTIAL         part of the gate is computed in its stated scope.
- BLOCKED         a named prerequisite *within the same route* is missing.
- NOT_TESTED      no producer has evaluated the gate yet.
- NOT_APPLICABLE  the evaluated route does not address this gate.
- REFUTED         a computation contradicts the gate for the evaluated model.

The original R113 source is one route (``original_source``). Its missing values block
G3 of that route and, through it, that route's G4. They do not block G4 or any other
gate for a different origin model: a route dependency is recorded under
``route_dependencies`` and leaves the gate status itself untouched.
"""
import json
from pathlib import Path

from source import MARKS

REQUIREMENTS = json.loads(Path(__file__).with_name("gates.json").read_text())
STATUSES = ("PASS", "CANDIDATE_PASS", "PARTIAL", "BLOCKED", "NOT_TESTED", "NOT_APPLICABLE", "REFUTED")
ORIGINAL_ROUTE = "original_source"


def _same_model(complete, gates):
    return (set(complete.get("same_model_gates",{})) == set(gates)
            and complete.get("proof",{}).get("status") == "PROVED")


def assess(original, synthesis=None, complete=None):
    original = original or {}
    synthesis = synthesis or {}
    complete = complete or {}
    role = original.get("input_role")
    status = original.get("status")
    values = original.get("original_values") or dict.fromkeys(MARKS)
    obstruction = synthesis.get("scoped_exclusion")
    source_missing = status=="BLOCKED"
    absolute_missing = original.get("relative_observables_only") is True
    gates = {key:{"status":"NOT_TESTED","required":list(items)} for key,items in REQUIREMENTS.items()}
    gates["G0"]["source_receipts"] = original.get("source_refs",[])
    gates["G1"]["scoped_operator_exclusion"] = obstruction
    gates["G3"].update(route=ORIGINAL_ROUTE,
                       status="BLOCKED" if source_missing else "NOT_TESTED",
                       original_values=values,missing=original.get("missing",[]),
                       supplied_control_only=role=="arithmetic_control_only")
    if absolute_missing:
        connected = original.get("connected_UY")
        gates["G3"].update(status="PARTIAL" if connected is not None else "BLOCKED",
                           relative_connected_UY=connected)
    if source_missing or absolute_missing:
        # Only the original route reaches G4 through G3. Other origin models are not blocked.
        gates["G4"]["route_dependencies"] = {ORIGINAL_ROUTE:{"status":"BLOCKED","blocked_by":["G3"]}}
    if complete.get("status") == "F1_PROVED" and complete.get("physical_solution") and _same_model(complete,gates):
        for gate,predicate in complete["same_model_gates"].items():
            gates[gate].pop("route_dependencies",None)
            gates[gate].update(status="PASS",predicate=predicate,
                               certificate_sha256=complete["proof"]["certificate_sha256"])
        return {"physical_status":"PROVED","physical_solution":True,
                "scope":"pinned native contract checked on one producer","gates":gates}
    if complete.get("status") == "F1_CONTRACT_PROVED" and _same_model(complete,gates):
        for gate,predicate in complete["same_model_gates"].items():
            gates[gate].update(status="CANDIDATE_PASS",predicate=predicate,
                               contract_sha256=complete.get("contract_sha256"),
                               certificate_sha256=complete["proof"].get("certificate_sha256"))
        return {"physical_status":"CANDIDATE_PROVED","physical_solution":False,
                "scope":"supplied candidate contract checked on one producer; physical acceptance "
                        "requires review and registration of that contract as the native contract",
                "gates":gates}
    return {"physical_status":"OPEN","physical_solution":False,
            "scope":"evidence report for current campaigns; PASS requires separately reviewed physical derivations",
            "gates":gates}
