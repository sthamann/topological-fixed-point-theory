"""Open marked covariance susceptibility of the actual fixed R032 half strip."""

import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import block_diag, expm, expm_frechet


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with contextlib.redirect_stdout(io.StringIO()):
    native_cut = runpy.run_path(str(base / "check_original_interference_cut_20261109.py"))

old, source, native = (native_cut[key] for key in ("old", "src", "native"))
identity, gamma, swap = (native_cut[key] for key in ("identity", "gamma", "swap"))
history = np.array(old["history"][:3])
widths = old["durations"][:3]
e, n, w = (native_cut["yukawa"](history[2], "e"),
           native_cut["yukawa"](history[0], "nu"),
           native_cut["mark"](native_cut["weak"]))
first = identity[:, native.IDX[("e", "R", 0, 0)]]
second = identity[:, native.IDX[("nu", "R", 0, 0)]]
tri = np.triu_indices(len(gamma), 1)
radial_terms = [-1j * gamma @ source["inner_form"](h) for h in history]


def exterior_pair(a, b):
    return (a[:, None] * b[None, :] - b[:, None] * a[None, :])[tri]


def operators(radius, shifts, with_tangents=False, momentum_sign=1):
    transfer = identity.copy()
    tangents = []
    for i, (width, shift) in enumerate(zip(widths, shifts)):
        generator, _, momentum = old["existing"]["generator"](i, "A", radius)
        generator = generator + (momentum_sign - 1) * momentum * gamma
        direction = -width * radial_terms[i]
        argument = -width * generator + shift * direction
        if with_tangents:
            heat, derivative = expm_frechet(argument, direction)
            tangents = [heat @ tangent for tangent in tangents] + [derivative @ transfer]
        else:
            heat = expm(argument)
        transfer = heat @ transfer

    # The original arm factors produce E K and W K N in the checker.
    electron = (1 + shifts[2]) * e @ transfer
    neutrino = (1 + shifts[0]) * w @ transfer @ n
    if not with_tangents:
        return electron, neutrino
    d_e = [e @ tangent for tangent in tangents]
    d_n = [w @ tangent @ n for tangent in tangents]
    d_e[2] += e @ transfer
    d_n[0] += w @ transfer @ n
    return transfer, electron, neutrino, d_e, d_n, tangents


def source_frame_test(transfer, d_e, d_n, tangents):
    # Same input states, bases passively transported, and independent frames
    # on both ends. For the two words the final frame is the same.
    l = native_cut["l"]
    r = native_cut["r"]
    er = native_cut["yukawa"](native_cut["g_r"] @ history[2], "e")
    nl = native_cut["yukawa"](native_cut["g_l"] @ history[0], "nu")
    wr = native_cut["mark"](native_cut["g_r"] @ native_cut["weak"] @ native_cut["g_r"].conj().T)
    reframed = r @ transfer @ l.conj().T
    left = [er @ r @ tangent @ l.conj().T for tangent in tangents]
    right = [wr @ r @ tangent @ l.conj().T @ nl for tangent in tangents]
    left[2] += er @ reframed
    right[0] += wr @ reframed @ nl
    return max(np.max(np.abs(a - r @ b @ l.conj().T))
               for a, b in zip(left + right, d_e + d_n))


def source_test(radius):
    transfer, electron, neutrino, d_e, d_n, tangents = operators(radius, np.zeros(3), True)
    k, _ = old["half"](radius, [identity] * len(old["links"]))
    assert np.max(np.abs(k - transfer)) < 1e-13
    r, sea = old["existing"]["direct_whole"]("S", radius)
    ref = native_cut["exact_from_chains"]([], [(r, False)], native_cut["center"], sea).real
    assert ref > 0
    de = [op @ first for op in d_e]
    dn = [op @ second for op in d_n]
    l0, r0 = electron @ first, neutrino @ second
    assert min(np.linalg.norm(l0), np.linalg.norm(r0)) > 1e-8
    assert max(np.linalg.norm(c @ x - x) for c in native_cut["center"]
               for x in [l0, r0, *de, *dn]) < 1e-13
    frames = source_frame_test(transfer, d_e, d_n, tangents)
    assert frames < 1e-13

    diagonal = [exterior_pair(de[i], dn[i]) for i in range(3)]
    off_diagonal = [(exterior_pair(de[i], dn[j]) + exterior_pair(de[j], dn[i])) / np.sqrt(2.)
                    for i in range(3) for j in range(i + 1, 3)]
    coefficients = np.column_stack(diagonal + off_diagonal)
    real_matrix = np.vstack([coefficients.real, coefficients.imag])
    singular_values = np.linalg.svd(real_matrix, compute_uv=False)
    assert singular_values[-1] > 1e-3
    target = np.array([1., 1., 1., 0., 0., 0.])
    radial_trace_probe = np.linalg.lstsq(real_matrix.T, target, rcond=None)[0]
    trace_residual = np.max(np.abs(real_matrix.T @ radial_trace_probe - target))
    assert trace_residual < 1e-12

    # The unchanged partner uses negative momentum and the same Higgs history.
    # For the A-even combination of the two marked inputs the ++ and --
    # output sectors are orthogonal, so conjugate terms cannot cancel as vectors.
    _, _, _, minus_e, minus_n, _ = operators(radius, np.zeros(3), True, -1)
    for plus, minus in zip(d_e + d_n, minus_e + minus_n):
        assert np.max(np.abs(swap @ plus.conj() @ swap - minus)) < 1e-13
    minus_de = [op @ (swap @ first.conj()) for op in minus_e]
    minus_dn = [op @ (swap @ second.conj()) for op in minus_n]
    minus_diag = [exterior_pair(minus_de[i], minus_dn[i]) for i in range(3)]
    minus_off = [(exterior_pair(minus_de[i], minus_dn[j]) + exterior_pair(minus_de[j], minus_dn[i])) / np.sqrt(2.)
                 for i in range(3) for j in range(i + 1, 3)]
    partner = np.column_stack(minus_diag + minus_off)
    paired = np.vstack([coefficients, partner]) / np.sqrt(2.)
    paired_singular = np.linalg.svd(np.vstack([paired.real, paired.imag]), compute_uv=False)
    assert np.max(np.abs(paired_singular - singular_values)) < 1e-13

    # An antithetic shift checks the coefficient through a two-copy identity.
    # Its test points carry no physical probability or bosonic normalization.
    directions = [*np.eye(3), np.ones(3), np.array([1., -1., .5])]
    errors = []
    for direction in directions:
        predicted = exterior_pair(sum(v * x for v, x in zip(direction, de)),
                                  sum(v * x for v, x in zip(direction, dn)))
        step = 1e-3
        e_plus, n_plus = operators(radius, step * direction)
        e_minus, n_minus = operators(radius, -step * direction)
        x1, x2 = e_plus @ first, e_minus @ first
        y1, y2 = n_plus @ second, n_minus @ second
        direct = ((exterior_pair(x1, y1) + exterior_pair(x2, y2)) / 2 -
                  exterior_pair((x1 + x2) / 2, (y1 + y2) / 2)) / step**2
        stable = exterior_pair((x1 - x2) / (2 * step), (y1 - y2) / (2 * step))
        rounding_bound = 32 * np.finfo(float).eps * np.linalg.norm(l0) * np.linalg.norm(r0) / step**2
        assert np.linalg.norm(direct - stable) < rounding_bound
        errors.append(float(np.linalg.norm(stable - predicted) / np.linalg.norm(predicted)))
    assert max(errors) < 1e-4
    return {
        "same_original_color_reference": ref,
        "two_yukawa_arms_with_existing_weak_mark_norms": [float(np.linalg.norm(l0)), float(np.linalg.norm(r0))],
        "covariance_map_real_singular_values": singular_values.tolist(),
        "covariance_map_real_rank": int(np.linalg.matrix_rank(real_matrix)),
        "covariance_map_on_antiunitary_even_pair_minimum_singular_value": float(paired_singular[-1]),
        "radial_trace_linear_probe_residual": float(trace_residual),
        "radial_trace_linear_probe_norm": float(np.linalg.norm(radial_trace_probe)),
        "per_segment_diagonal_covariance_coefficients": list(map(float, map(np.linalg.norm, diagonal))),
        "largest_relative_fd_error": max(errors),
        "passive_common_boundary_frame_error": float(frames),
        "neutral_under_the_actual_color_center": True,
        "under_full_original_P_adm": None,
    }


def coherence_bound():
    right = [native.IDX[("nu", "R", 0, i)] for i in range(2)]
    t = identity[:, native.IDX[("nu", "L", 0, 0)]]
    target_after = source["real_action"](expm(.2j * native_cut["weak"])) @ t
    overlap = np.vdot(t, target_after)
    orthogonal = np.linalg.norm(target_after - overlap * t)
    matrix = np.array([[1., 0.], [overlap, orthogonal]]) / np.sqrt(2.)
    lower, upper = np.linalg.eigvalsh(matrix @ matrix.conj().T)
    # Singular values of the antisymmetric 4-mode matrix occur twice.
    skew = np.block([[np.zeros((2, 2)), matrix], [-matrix.T, np.zeros((2, 2))]])
    values = np.linalg.svd(skew, compute_uv=False)**2
    assert np.max(np.abs(values - [upper, upper, lower, lower])) < 1e-13
    assert np.allclose(source["real_action"](expm(.2j * native_cut["weak"])) @ identity[:, right],
                       identity[:, right])
    return {"max_fidelity_for_positive_mixture_of_fixed_decomposable_outputs": float(upper),
            "target_infidelity_floor_in_that_class": float(lower)}


print(json.dumps({
    "scope": "Fixed native R032 96-component half strip, same original open arm order, lepton color singlets; local radial susceptibility, not a physical R113 average.",
    "half_time": sum(widths),
    "rows": {str(radius): source_test(radius) for radius in (-.001, .001)},
    "control_target_incoherent_mixture_bound": coherence_bound(),
}, indent=2))
