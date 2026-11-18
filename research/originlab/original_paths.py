"""The original affine regional word and its first R113 binding obstruction."""
from pathlib import Path
import json
import subprocess
import sys

import sympy as sp

from reconstruction import require
from verify import digest

AFFINE = Path("docs/research/fundamental-plan/TFPT-E8-GELADENE-DREIREGIONENORDNUNG-20261111.txt")
Q_CHARGE = Path("docs/research/fundamental-plan/p1-original-marken-und-metrik-20261112/"
                "TFPT-P1-ORIGINAL-MARKEN-UND-METRIK-20261112.txt")
NATIVE_RUNTIME = Path("experiments/theory-contracts/original-paper-source-reconciliation-20261102/"
                      "contribution_20261106Q/runtime")


def native_readout(root,domain,contract):
    # Direct call to the original consumer. The separate discovery engine
    # uses none of its optimization or verification scripts.
    program = """
import hashlib, json, sys
from pathlib import Path
sys.path.insert(0,sys.argv[1])
from source_adapter import evaluate_contribution, evaluate_source
path = Path(sys.argv[2])
raw = path.read_bytes()
domain = json.loads(raw)
domain['domain_directory'] = str(path.parent)
assembled = {'domain_id':domain['domain_id'],'native_inputs_sha256':hashlib.sha256(raw).hexdigest(),
             'h_N_role':'not_used_for_original_endpoint_probe'}
whole = evaluate_source({},assembled,domain)
fibre = evaluate_contribution({'schema':'tfpt.original-source-unit.v1',
                              'original_configuration_id':None},assembled,domain)
print(json.dumps({'whole_status':whole['status'],'values':whole['values'],
                  'fibre_status':fibre['status'],'fibre_missing':fibre.get('missing_operations',[])}))
"""
    # The original runtime uses NumPy. The independent environment only
    # needs its exact backends, so allow a system interpreter for this call.
    commands = [sys.executable,"python3"]
    for python in dict.fromkeys(commands):
        try:
            completed = subprocess.run([python,"-B","-c",program,str(root/NATIVE_RUNTIME),str(domain)],
                                       cwd=root,text=True,capture_output=True,check=True,timeout=30)
            native = json.loads(completed.stdout)
            break
        except (OSError,ValueError,subprocess.SubprocessError) as error:
            failure = type(error).__name__
    else:
        return {**contract,"native_endpoint_readout":{"status":"UNAVAILABLE","error_type":failure}}
    values = native.get("values") or {}
    if contract["status"]=="BLOCKED" and any(value is not None for value in values.values()):
        return {**contract,"status":"ORIGINAL_NATIVE_VALUES_REQUIRE_REVIEW","native_endpoint_readout":native}
    return {**contract,"native_endpoint_readout":native}


def original_word(root):
    refs = {str(path):digest(root/path) for path in (AFFINE,Q_CHARGE)}
    require(list(refs.values())==[
        "0d210ea85d48ec6d99bb7079cef382e17a1ea413c1a5d94bcaa7ca380bfa9f73",
        "3cef7c6b5039ddf466868b87d5f162071724fa23c971ed55ab5647e152d7d32c",
    ],"review the changed original regional or charge source")
    half = sp.S.Half
    b1,b2 = (sp.eye(4)[:,j]-sp.ones(4,1)/4 for j in (0,1))
    alpha = sp.Matrix([half]*5).col_join(b1)
    beta = sp.Matrix([half,-half,-half,-half,-half]).col_join(b2)
    gamma = alpha+beta
    require(alpha.dot(beta)==-1 and [v.dot(v) for v in (alpha,beta,gamma)]==[2]*3,
            "wrong original affine root dictionary")
    a,b,c,w = sp.S.One,sp.I,sp.S.NegativeOne,-sp.I
    points = (a,b,c)
    value = sp.simplify(a*b*c/((a-b)*(b-c)*(a-c)))
    # The weight-one Ward operator on the same three insertions, in
    # the original Cartan convention in check_tfpt_cartan_process.
    ward = sum(1/(w-z)**2-sum(1/(z-x) for x in points if x!=z)/(w-z) for z in points)
    require(value==-sp.I/4 and sp.simplify(ward)==sp.Rational(3,4),
            "original E8 current and stress word did not replay")
    charge = sp.Matrix([-2,-2,-2,3,3,0,0,0])/12
    p1 = sp.ones(8,1)
    p2 = sp.Matrix([-1,1,1,1,1,-1,1,1])
    y_gamma,y_weak = charge.dot(p1-p2),sp.S.Zero
    require(y_gamma==-sp.Rational(1,3) and charge.dot(p2)==sp.Rational(1,3),
            "original D8 charge of the middle current changed")
    return {"status":"EXACT_ORIGINAL_AFFINE_WORD","source_refs":refs,
            "common_ordered_word":"J^(-alpha) J^gamma J^(-beta)",
            "original_scope":"given E8 level-one current vacuum with its original point and cocycle convention",
            "vacuum_word":str(value),"cartan_stress_factor":str(sp.simplify(ward)),
            "r113_binding":{"status":"BLOCKED","absolute_values":[],
                           "charge_of_original_middle":{"Y":str(y_gamma),"Q":"-1/3","triality":1},
                           "charge_of_r032_weak_comparison":{"Y":str(y_weak),"Q":"1"},
                           "interchangeable_by_the_supplied_charge_dictionary":y_gamma==y_weak,
                           "first_missing":"embed the original middle current and its compensating endpoints in the same weighted R113 word before the outer selection"}}
