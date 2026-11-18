"""Check the actual native intertwiners against the cyclotomic trace metric."""
from pathlib import Path
import json
import sympy as sp

root = Path(__file__).resolve().parent
data = json.loads((root / "operator_results.json").read_text())["native_F20"]
# Tr_Q(zeta5)/Q(conj(zeta5^i) zeta5^j), in the original power basis.
# The diagonal trace is 4; every off-diagonal trace is -1.
trace_metric = 5 * sp.eye(4) - sp.ones(4)
checks = []
for key in ("v419_C5", "v419_G4"):
    operator = sp.Matrix(data[key])
    ok = sp.simplify(operator.H * trace_metric * operator - trace_metric) == sp.zeros(4)
    if not ok:
        raise ArithmeticError(f"The original operator {key} fails trace-metric invariance")
for row in data["rows"]:
    intertwiner = sp.Matrix(row["cyclotomic_intertwiner"])
    pullback = sp.simplify(intertwiner.H * intertwiner)
    scale = sp.simplify(pullback[0, 0] / 4)
    ok = scale > 0 and sp.simplify(pullback - scale * trace_metric) == sp.zeros(4)
    if not ok:
        raise ArithmeticError(f"Intertwiner {row['sylow_index']} fails metric compatibility")
    checks.append({"sylow_index": row["sylow_index"], "positive_metric_scale": str(scale), "exact": True})
result = {
    "status": "PASS_EXACT_TRACE_METRIC_INTERTWINERS",
    "trace_metric": [[str(x) for x in trace_metric.row(i)] for i in range(4)],
    "rows": checks,
    "conclusion": "A/sqrt(scale) is an isometry from the original cyclotomic trace metric to the native Euclidean metric, for each actual intertwiner.",
    "scope": "This fixes a representation metric, not a physical field or action normalization.",
}
(root / "galois_metric_results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(json.dumps(result, indent=2))
