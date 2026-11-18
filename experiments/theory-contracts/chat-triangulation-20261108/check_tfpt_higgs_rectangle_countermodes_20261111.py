"""Original R032 fixed history and exact old Fold packet controls."""
from contextlib import redirect_stdout
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np
from scipy.linalg import expm

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
sys.path.insert(0, str(repo / "verification"))
with redirect_stdout(io.StringIO()):
    argv = sys.argv
    try:
        sys.argv = [argv[0], str(base)]
        native = runpy.run_path(str(base / "check_r113_dressed_open_20261109.py"))
    finally:
        sys.argv = argv
    fold = runpy.run_path(str(base / "check_tfpt_cover_selector_compiler_20261110.py"))
    affine = runpy.run_path(str(base / "check_tfpt_family_middle_ordered_cut_20261110.py"))

gamma, mass, higgs, gauge, multiply = (
    native[k] for k in ("grading", "majorana", "inner_form", "real_action", "ordered")
)
sp, q, Q, vertex, phase, signed = (
    fold[k] for k in ("sp", "q", "Q", "vertex", "phase", "signed")
)
inner, combine = (affine[k] for k in ("inner", "combine"))
identity = np.eye(len(gamma))
rep = native["native"]
omega = np.exp(2j * np.pi / 3)
color = np.array([0 if rep.P[i % rep.NP][0] in ("nu", "e")
                  else (1 if i < rep.NP else -1) for i in range(len(gamma))])
center = np.diag(omega**color)
window, width = .2, .55
durations = [window / 2, window, window, window, window / 2]
total_time = sum(durations)
real, imaginary = width, 1j * width


def transfer(momentum, field, duration):
    h = momentum * gamma - 1j * gamma @ (mass + higgs(np.array([0, field])))
    np.testing.assert_allclose(h, h.conj().T, atol=1e-14)
    e, v = np.linalg.eigh(h)
    return (v * np.exp(-duration * e)) @ v.conj().T, float(duration * e[:48].sum())


def word(blocks, links=None):
    if links is None:
        links = [identity] * (len(blocks) - 1)
    return multiply(blocks, links, [np.zeros_like(identity)] * len(blocks))[0]


def reflected_center_trace(effect):
    # Exterior/Fock trace with the original finite color-center average.
    traces = []
    for k in range(3):
        sign, logabs = np.linalg.slogdet(identity + np.diag(omega**(k * color)) @ effect)
        traces.append(sign * np.exp(logabs))
    value = sum(traces) / 3
    assert value.real > 0 and abs(value.imag) < 1e-12 * value.real
    return float(value.real)


rectangle = {}
for p in (.25, .75):
    fields = [0, real, real + imaginary, imaginary, 0]
    data = [transfer(p, h, d) for h, d in zip(fields, durations, strict=True)]
    blocks, sea = zip(*data, strict=True)
    independent = identity.copy()
    for field, d in zip(fields, durations, strict=True):
        generator = p * gamma - 1j * gamma @ (mass + higgs(np.array([0, field])))
        independent = expm(-d * generator) @ independent
    fwd, rev = word(blocks), word(blocks[::-1])
    np.testing.assert_allclose(fwd, independent, atol=1e-13)
    np.testing.assert_allclose(fwd.conj().T, rev, atol=1e-13)
    np.testing.assert_allclose(fwd @ center, center @ fwd, atol=1e-13)
    plus, minus = fwd.conj().T @ fwd, rev.conj().T @ rev
    contrast = plus - minus
    # Use the same zero-Higgs reflected reference, including the same
    # Majorana and quarter data, and keep the specified sea convention.
    empty, empty_sea = transfer(p, 0, total_time)
    denominator = reflected_center_trace(empty.conj().T @ empty)
    returns = [np.exp(2 * (sum(sea) - empty_sea)) * reflected_center_trace(k) / denominator
               for k in (plus, minus)]
    np.testing.assert_allclose(*returns, rtol=2e-13)
    # Passively change the weak frame at interior cuts, with transported links.
    frames = [gauge(np.array([[np.cos(x), np.sin(x)], [-np.sin(x), np.cos(x)]]))
              for x in (0, .17, -.23, .31, 0)]
    dressed = [u @ a @ u.conj().T for u, a in zip(frames, blocks, strict=True)]
    links = [v @ u.conj().T for u, v in zip(frames[:-1], frames[1:], strict=True)]
    np.testing.assert_allclose(word(dressed, links), fwd, atol=1e-13)
    eigenvalues = np.linalg.eigvalsh(contrast)
    assert eigenvalues[0] < -1e-4 and eigenvalues[-1] > 1e-4
    rectangle[str(p)] = {
        "reflected_fixed_history_return_both_orders": returns,
        "open_order_difference_frobenius": float(np.linalg.norm(fwd-rev)),
        "effect_contrast_min_max": [float(eigenvalues[0]), float(eigenvalues[-1])],
        "covariance_residual": float(np.linalg.norm(word(dressed, links)-fwd)),
    }


def mode(momentum, step, opposite_root_sign):
    # The two terms have distinct outgoing charges for a single input.
    terms = {}
    for charge, factor in ((q, 1), (-q, opposite_root_sign)):
        a = Q.from_sympy(factor * phase(charge, momentum))
        terms[tuple(momentum + charge)] = {
            w: a*c for w, c in vertex({(): Q.one}, charge, momentum, step).items()
        }
    return terms


def norm_form(upper, lower):
    return sp.Matrix(2, 2, lambda i, j: sp.simplify(sum(
        inner((upper, lower)[i][p], (upper, lower)[j][p])
        for p in (upper, lower)[i].keys() & (upper, lower)[j].keys()
    )))


def scale(terms, multiplier):
    a = Q.from_sympy(multiplier)
    return {p: {w: a*c for w, c in x.items()} for p, x in terms.items()}


phase_q = phase(q, -q)
assert phase_q in (1, -1) and phase_q == phase(-q, q)
packet_data = {}
preps = [sp.Matrix(x) for x in ((1, 0), (0, 1), (1, 1), (1, sp.I))]
def four_reads(form):
    return [str(sp.simplify((x.H @ form @ x)[0])) for x in preps]


for alpha, rows in signed["controls"].items():
    packet_data[alpha] = {}
    for sign in (-1, 1):
        up, low = [r for r in rows if r["sign"] == sign]
        charges = [sp.Matrix([sp.sympify(x) for x in r["charge_old"]]) for r in (up, low)]
        top, bottom = [q + x for x in charges]
        assert top-bottom == 2*q
        phases = [phase(-q, top), phase(q, bottom)]
        # Verify the actual old packets isolate the outer ground components
        # with the same phase convention as the existing four-probe study.
        a, b = [sp.Matrix([[phase(c, q)*sp.sympify(z) for z in r[
            "exact_ground_amplitudes_up_to_old_cocycle"]]]) for c, r in zip(charges, (up, low))]
        i = 2 if sign == -1 else 1
        basis = sp.eye(5)
        upper = (basis[:, i] - b[i]/b[0]*basis[:, 0]) / (a[i]*phases[0])
        lower = basis[:, 0] / (b[0]*phases[1])
        coeffs = (sp.Matrix.vstack(a, b) @ sp.Matrix.hstack(upper, lower)).applyfunc(sp.simplify)
        assert coeffs == sp.diag(*(1/x for x in phases))
        matrices = {}
        for step in (1, -1):
            columns = [scale(mode(p, step, phase_q), coeffs[k, k])
                       for k, p in enumerate((top, bottom))]
            matrices[step] = norm_form(*columns)
        assert matrices[1] - matrices[-1] == 2 * sp.eye(2)
        assert matrices[1].rank() == 2 and matrices[-1].rank() == 1
        assert matrices[-1] @ preps[2] == sp.zeros(2, 1)
        packet_data[alpha][str(sign)] = {
            "raised_form": [[str(x) for x in row] for row in matrices[1].tolist()],
            "lowered_form": [[str(x) for x in row] for row in matrices[-1].tolist()],
            "raised_four_packets": four_reads(matrices[1]),
            "lowered_four_packets": four_reads(matrices[-1]),
            "difference_four_packets": four_reads(matrices[1]-matrices[-1]),
        }

print(json.dumps({
    "scope": "R032 finite fixed histories and exact old Fold outer ground channel only",
    "chosen_window_width": window,
    "chosen_higgs_square_side": width,
    "same_zero_Higgs_reference_at_reflected_time": 2*total_time,
    "rectangle": rectangle,
    "old_current_relative_root_phase": str(phase_q),
    "old_countermodes_before_time_and_physical_selection": packet_data,
}, indent=2))
