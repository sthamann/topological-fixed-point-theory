"""P numeric source consumer with the independently reviewed Q math checks.

No primitive contribution producer is supplied by this registry wiring.
"""
from pathlib import Path
import importlib.util
from original_source_adapter import (
    evaluate_contribution, integration_plan, aggregate_contributions,
    evaluate_source, capability_report as whole_capability_report,
    contribution_capability_report as _contribution_capability,
    ensemble_test_registry as _p_registry,
)

CHECKS_PATH=Path(__file__).resolve().parents[1]/'checks/original_checks.py'
_spec=importlib.util.spec_from_file_location('tfpt_q_original_mathematical_checks',CHECKS_PATH)
CHECKS=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(CHECKS)

def ensemble_test_registry(domain):
    return {**_p_registry(domain),**CHECKS.ensemble_test_registry(domain)}

def contribution_capability_report(domain):
    report=_contribution_capability(domain)
    report['registered_ensemble_test_ids']=list(ensemble_test_registry(domain))
    report['primitive_producer_added_by_Q_registry']=False
    return report

def run_ensemble_checks(result,entries,domain,required_ids=None):
    registry=ensemble_test_registry(domain)
    ids=list(registry) if required_ids is None else required_ids
    return [registry[cid](result,entries,domain) if cid in registry else
            {'id':cid,'status':'MISSING_IMPLEMENTATION'} for cid in ids]
