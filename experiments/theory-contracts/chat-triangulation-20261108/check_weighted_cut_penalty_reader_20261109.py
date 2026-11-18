"""Read genuine penalty identities and original Yukawa/weak vertex types."""
from pathlib import Path
import contextlib
import hashlib
import io
import json
import runpy
import sys

import numpy as np
import sympy as sp

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0,str(repo/"verification"))
folder = repo / "experiments/theory-contracts"
r032_path = folder / "chat-triangulation-20261108/check_r113_dressed_open_20261109.py"
args = sys.argv
try:
    sys.argv = [str(r032_path), str(r032_path.parent)]
    with contextlib.redirect_stdout(io.StringIO()):
        r032 = runpy.run_path(str(r032_path))
finally:
    sys.argv = args

v = r032["native"]
left = [v.IDX[(f,"L",0,j)] for f in ("nu","e") for j in range(3)]
nu = [v.IDX[("nu","R",0,j)] for j in range(3)]
electron = [v.IDX[("e","R",0,j)] for j in range(3)]
quaternion = r032["quaternion"]
hi,hj = map(np.asarray,r032["history"][:2])
Y_e,Y_n = r032["qref"].yukawa[2:4]
U = np.eye(2)
Q = np.array([[0,1],[1,0]],dtype=complex)


def transverse(h, w):
    p = np.outer(h,h.conj())/np.vdot(h,h).real
    return p @ w @ (np.eye(2)-p) + (np.eye(2)-p) @ w @ p


def original_endpoint(h, right):
    return r032["inner_form"](h)[np.ix_(left,right)]


def mixed_vertex(h_left,h_right,weak_link,weak_field):
    w = transverse(h_right,weak_field)
    lift = weak_link / np.linalg.det(weak_link)
    e = original_endpoint(h_right,electron)
    n = original_endpoint(h_left,nu)
    full = e.conj().T @ np.kron(w @ lift,np.eye(3)) @ n
    scalar = np.vdot(h_right,w @ lift @ quaternion(h_left)[:,0])
    assert np.linalg.norm(full-scalar*Y_e.conj().T @ Y_n)<2e-13
    return full, scalar


vertex, coeff = mixed_vertex(hi,hj,U,Q)
# Independent SU2 and hypercharge endpoint changes. W carries only the
# adjoint action at the right endpoint; the right electron has det^-1.
u_i = np.exp(.29j)*quaternion(np.array([.6,.8j]))
u_j = np.exp(-.18j)*quaternion(np.array([.8j,.6]))
moved,_ = mixed_vertex(u_i@hi,u_j@hj,u_j@U@u_i.conj().T,u_j@Q@u_j.conj().T)
assert np.linalg.norm(moved-vertex/np.linalg.det(u_j))<2e-13
assert np.linalg.norm(vertex)>0

# E = sector even but not structurally admitted; G = structural below gap;
# P = retained; O = determinant-odd. Exact inverse of the four original
# linear insertions. a,b are diagnostic attenuation factors, not dynamics.
a,b = sp.symbols("a b",positive=True)
original_reads = sp.Matrix([[1,1,1,1],[1,1,1,a],[1,1,0,0],[1,b,0,0]])
read_inverse = sp.Matrix([[0,0,-b/(1-b),1/(1-b)],
                          [0,0,1/(1-b),-1/(1-b)],
                          [-a/(1-a),1/(1-a),-1,0],
                          [1/(1-a),-1/(1-a),0,0]])
assert (read_inverse*original_reads-sp.eye(4)).applyfunc(sp.simplify)==sp.zeros(4)

# Active original callers and every registered domain in this contract.
caller = folder / "original-paper-source-reconciliation-20261102/contribution_20261106Q"
sys.path.insert(0, str(caller/"runtime"))
from assemble import load_domain
from source_adapter import evaluate_contribution, evaluate_source, whole_capability_report

domains = []
for path in sorted((caller/"domain").glob("*domain*.json")) + [caller/"domain/native-diagnostic-v1.json"]:
    domain, _, _ = load_domain(path)
    capability = whole_capability_report(domain)
    domains.append({"name":path.name,"physical_source_identity":domain["physical_source_identity"],
                    "bundle":domain.get("whole_source_bundle"),"bundle_ready":capability["numeric_bundle_ready"]})
    if path.name=="native-diagnostic-v1.json":
        assembled = {"domain_id":domain["domain_id"],"h_N_role":"not_used_for_endpoint_reader",
                     "native_inputs_sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
        base = evaluate_source({},assembled,domain)
        fiber = evaluate_contribution({"schema":"tfpt.original-source-unit.v1",
                                       "original_configuration_id":None},assembled,domain)
        assert all(x is None for x in base["values"].values())
        assert base["missing_operations"]==["evaluate_original_full_blocked_endpoint_source"]
        assert fiber["missing_operations"]==["original_blocked_fiber_packet"]

print(json.dumps({
    "scope":"given native R032 Yukawa vertex and fixed weak insertion before physical P_adm, exact QFT002 penalty identities, actual source call readiness",
    "kinematic_weak_vertex_uses_native_first_link_and_formal_sigma_x_probe":True,
    "unweighted_endpoint_complex_coefficient":[float(coeff.real),float(coeff.imag)],
    "unweighted_endpoint_block_frobenius_square":float(np.vdot(vertex,vertex).real),
    "fixed_input_reference_chart":"R032 original free S rule, not an identified physical denominator",
    "both_diagnostic_radii_kinematic_vertex_unchanged_under_p_only_variation":[-.001,.001],
    "weak_vertex_reference_compression_is_not_computed":True,
    "observations":["I","exp(-t_theta KTheta^2)","Pstr","Pstr exp(-t_gap Kg^2)"],
    "sector_basis":["P_adm","Pstr-P_adm","PTheta-Pstr","I-PTheta"],
    "exact_reader_inverse":[[str(x) for x in row] for row in read_inverse.tolist()],
    "available_domains":domains,
    "actual_source_status":base["status"],"missing_whole":base["missing_operations"],
    "missing_fibre":fiber["missing_operations"],"weighted_original_metric_pairs":0
},indent=2))
