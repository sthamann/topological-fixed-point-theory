"""Read the original two-face collar and the installed whole-source gate."""
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp

repo = Path(sys.argv[1]).resolve()
contracts = repo / "experiments/theory-contracts"
caller = contracts / "original-paper-source-reconciliation-20261102/contribution_20261106Q"
sys.path.insert(0, str(caller / "runtime"))

from assemble import load_domain
from source_adapter import evaluate_contribution, evaluate_source

domain_file = caller / "domain/native-diagnostic-v1.json"
domain, _, _ = load_domain(domain_file)
domain_sha = hashlib.sha256(domain_file.read_bytes()).hexdigest()
assembled = {"domain_id": domain["domain_id"], "native_inputs_sha256": domain_sha,
             "h_N_role": "not_used_for_whole_source_gate"}
requests = {}
for leg in ("ramond", "higgs"):
    # candidate_id distinguishes calls. The installed interface does not
    # accept translated halfwords as executable source inputs.
    whole = evaluate_source({"candidate_id": "calderon-shift-" + leg}, assembled, domain)
    fiber = evaluate_contribution({"schema": "tfpt.original-source-unit.v1",
                                   "original_configuration_id": None}, assembled, domain)
    assert whole["status"] == "original_full_source_operation_not_supplied"
    assert whole["missing_operations"] == ["evaluate_original_full_blocked_endpoint_source"]
    assert all(x is None for x in whole["values"].values())
    assert fiber["missing_operations"] == ["original_blocked_fiber_packet"]
    assert fiber["original_measure_mass"] is None
    requests[leg] = {"candidate_id": whole["candidate_id"],
                     "whole_source_gate": whole["status"],
                     "missing_whole_operation": whole["missing_operations"][0],
                     "missing_fiber_operation": fiber["missing_operations"][0],
                     "registered_endpoint_values": whole["values"],
                     "physical_mixed_value": None}

# Original constant product collar at the actual old quarter boundary.
# Its signed h-pair and normal Clifford action are from the original
# energy-source proof, up to choice of eigenvector phases.
E = sp.Rational(1, 4)
ell = 2 * sp.pi
h = sp.diag(E, -E)
N = sp.Matrix([[0, -1], [1, 0]])
J_faces = sp.diag(N, -N)
assert h*N+N*h == sp.zeros(2)
q = sp.symbols("q", positive=True)
T = sp.diag(q,1/q)
Id = sp.eye(2)
graph = Id.col_join(T)
G = graph*(graph.H*graph).inv()*graph.H
A = (Id+T**2).inv()
green = A.row_join(T*A).col_join((T*A).row_join(T**2*A))
assert (G-green).applyfunc(sp.simplify) == sp.zeros(4)
assert (G*J_faces*G).applyfunc(sp.simplify) == sp.zeros(4)

# Test each placement of the same positive shift. Only the geometric
# collar h is used here; it is not assigned to the physical OS clock.
d = sp.Rational(1,20)
shift = sp.diag(*[sp.exp(-d*x) for x in (E,-E,E,-E)])
placed_left = shift*G
placed_right = G*shift
assert (placed_left-placed_right).applyfunc(sp.simplify) == sp.zeros(4)
q_value = sp.exp(-ell*E)
cross_before = G[0,2].subs(q,q_value)
cross_after = placed_left[0,2].subs(q,q_value)
assert sp.simplify(G[0,2]-q/(1+q*q)) == 0

# Independent use of the actual original real-energy Calderon derivative.
# This is not a window translation. It checks that a nonzero change of
# range gives no first-order direct insertion between retained traces.
a = (1-q*q)/(E*(1+q*q))
b = 2*q/(1+q*q)
Q_prime = (a*a*sp.kronecker_product(Id,N*h)+
           a*b*sp.kronecker_product(-N,N))/2
assert Q_prime == Q_prime.H
assert (Q_prime**2-a*a*sp.eye(4)/4).applyfunc(sp.simplify) == sp.zeros(4)
assert (G*Q_prime*G).applyfunc(sp.simplify) == sp.zeros(4)

print(json.dumps({
    "original_domain": domain["domain_id"],
    "domain_sha256": domain_sha,
    "whole_source_bundle": domain["whole_source_bundle"],
    "distinct_requests_at_installed_gate": requests,
    "evaluated_physical_shifted_words": 0,
    "original_geometric_quarter_pair_only": {
        "h_eigenvalues": [str(E), str(-E)], "collar_half_width": str(ell),
        "shift_in_geometric_collar_parameter": str(d),
        "cross_face_before_shift": float(cross_before),
        "cross_face_shift_on_left": float(cross_after),
        "cross_face_shift_on_right": float(placed_right[0,2].subs(q,q_value)),
        "entire_pair_left_minus_right_exact_zero": True,
        "bare_oriented_green_form_on_cauchy_range_exact_zero": True,
        "real_energy_range_derivative_operator_norm": float((a/2).subs(q,q_value)),
        "real_energy_derivative_between_retained_traces_exact_zero": True,
    },
}, indent=2))
