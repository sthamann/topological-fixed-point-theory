"""Isolate the longitudinal phase in the already specified R032 comparison."""
import contextlib
import io
import json
import runpy
import sys
from pathlib import Path

import numpy as np

root = Path(sys.argv[1]).resolve()
base = root / "experiments/theory-contracts/chat-triangulation-20261108"
check = base / "check_normal_weak_curvature_bridge_20261109.py"
argv = sys.argv
try:
    sys.argv = [str(check), str(root)]
    with contextlib.redirect_stdout(io.StringIO()):
        old = runpy.run_path(str(check))
finally:
    sys.argv = argv

history = old["history"]
projectors = old["projectors"]
identity = old["identity"]
half, ref = old["half"], old["existing"]["direct_whole"]
real_action = old["source"]["real_action"]
mark = old["marks"][0]
g, u = np.linalg.eigh(mark)
basis_plus, basis_minus = u[:,g>0], u[:,g<0]
eye2 = np.eye(2)
changes = []
for pi, pj, hi, hj in zip(projectors[:2], projectors[1:3], history[:2], history[1:3]):
    overlap = np.sqrt(np.trace(pi @ pj).real)
    ray = (pj @ pi + (eye2-pj) @ (eye2-pi))/overlap
    angle = -np.angle(np.vdot(hj, hi))
    changes.append((ray, pj, angle))


def links(fractions):
    out = []
    for (ray, pj, angle), fraction in zip(changes, fractions):
        a = np.exp(1j*fraction*angle)
        u = (a*pj+a.conjugate()*(eye2-pj)) @ ray
        assert abs(np.linalg.det(u)-1)<1e-13
        out.append(real_action(u))
    return out + [identity] + [u.conj().T for u in reversed(out)]


def number_traces(b):
    assert np.linalg.norm(mark @ b-b @ mark)<1e-11
    def coefficients(block):
        values = np.linalg.eigvalsh(block)
        assert min(values)>0
        coeff = np.array([1.])
        for t in values:
            coeff = np.convolve(coeff, [1.,t])
        return coeff

    plus = coefficients(basis_plus.conj().T @ b @ basis_plus)
    minus = coefficients(basis_minus.conj().T @ b @ basis_minus)
    even = minus.copy()
    even[1::2] = 0
    return np.convolve(plus,even), np.convolve(plus,minus)


rows = {}
for m in (-.001,.001):
    reference, ref_sea = ref("S",m)
    ref_selected, ref_full = number_traces(reference)
    vals = {}
    phases = [("ray",(0,0)),("first_only",(1,0)),("second_only",(0,1)),
              ("midpoint",(.5,.5)),("vector",(1,1))]
    for name, fractions in phases:
        trial = links(fractions)
        k, sea = half(m,trial)
        k_s, sea_s = half(m,trial,"S")
        assert np.linalg.norm(k_s.conj().T @ k_s-reference)<1e-11
        assert abs(2*sea_s-ref_sea)<1e-12
        b = k.conj().T @ k
        selected, full = number_traces(b)
        weight = np.exp(2*sea-ref_sea)
        sel_value = weight*selected.sum()/ref_selected.sum()
        assert abs(sel_value/old["response"](m,trial)[0][0,0]-1)<1e-11
        vals[name] = {
            "relative_selected_empty":float(sel_value),
            "relative_full_unselected_empty":float(weight*full.sum()/ref_full.sum()),
            "one_particle_sector_trace_before_reference":float(selected[1]),
            "dominant_exterior_degree":int(np.argmax(selected)),
            "sea_factor_original_over_reference":float(weight),
        }
    first, second, ray, vector = [vals[key]["relative_selected_empty"]
        for key in ("first_only","second_only","ray","vector")]
    rows[str(m)] = {
        "phase_readouts": vals,
        "selected_vector_over_ray":vector/ray,
        "unselected_vector_over_ray":vals["vector"]["relative_full_unselected_empty"]
                                       / vals["ray"]["relative_full_unselected_empty"],
        "two_segment_nonlinearity_log":float(np.log(vector*ray/(first*second))),
    }

print(json.dumps({"scope":"previous two changed R032 transports; native local H and time, same S reference and same candidate projection",
                  "longitudinal_phase_corrections_rad":[float(x[2]) for x in changes],
                  "rows":rows},indent=2))
