"""Report the evidence actually available for the plan's physical gates."""
import json
from pathlib import Path

from source import MARKS

REQUIREMENTS = json.loads(Path(__file__).with_name("gates.json").read_text())


def assess(original, synthesis=None, complete=None):
    original = original or {}
    synthesis = synthesis or {}
    role = original.get("input_role")
    status = original.get("status")
    values = original.get("original_values") or dict.fromkeys(MARKS)
    obstruction = synthesis.get("scoped_exclusion")
    source_missing = status=="BLOCKED"
    absolute_missing = original.get("relative_observables_only") is True
    gates = {key:{"status":"UNKNOWN","required":list(items)} for key,items in REQUIREMENTS.items()}
    gates["G0"]["source_receipts"] = original.get("source_refs",[])
    gates["G1"]["scoped_operator_exclusion"] = obstruction
    gates["G3"].update(status="BLOCKED" if source_missing else "UNKNOWN",
                       original_values=values,missing=original.get("missing",[]),
                       supplied_control_only=role=="arithmetic_control_only")
    if absolute_missing:
        connected = original.get("connected_UY")
        gates["G3"].update(status="PARTIAL" if connected is not None else "BLOCKED",
                           relative_connected_UY=connected)
    if source_missing or absolute_missing:
        gates["G4"]["status"] = "BLOCKED"
        gates["G4"]["blocked_by"] = ["G3"]
    if (complete and complete.get("status") == "F1_PROVED" and complete.get("physical_solution")
            and set(complete.get("same_model_gates",{})) == set(gates)
            and complete.get("proof",{}).get("status") == "PROVED"):
        for gate,predicate in complete["same_model_gates"].items():
            gates[gate].pop("blocked_by",None)
            gates[gate].update(status="PASS",predicate=predicate,
                               certificate_sha256=complete["proof"]["certificate_sha256"])
        return {"physical_status":"PROVED","physical_solution":True,
                "scope":"pinned native contract checked on one producer","gates":gates}
    return {"physical_status":"OPEN","physical_solution":False,
            "scope":"evidence report for current campaigns; PASS requires separately reviewed physical derivations",
            "gates":gates}
