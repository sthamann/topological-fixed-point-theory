"""Exact finite coherent kernels for even CAR operators.

Conditional algebra on the existing collar path. This reconstructs finite
operator matrices and sews their kernels; it does not select a P1/P2 source or
claim a physical field theory. Fock bit ordering is the source_sewing ordering.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re

import sympy as sp

from joint_source import add, clean, exponential, multiply, scale, variable
from source_sewing import source_lift


def _check_operator(matrix):
    t = sp.Matrix(matrix)
    if t.rows != t.cols or t.rows not in (2, 4, 8, 16, 32):
        raise ValueError("exact kernel reconstruction accepts 1..5 CAR modes")
    n = t.rows.bit_length() - 1
    for i in range(t.rows):
        for j in range(t.cols):
            if ((i.bit_count() - j.bit_count()) & 1) and t[i, j] != 0:
                raise ValueError("kernel chart accepts parity-even operators only")
    return t, n


def _sign(number):
    return (-1) ** (number * (number - 1) // 2)


def kernel(matrix, bar_start=0, field_start=None):
    """K_T(bar z,z)=sum_IJ (-1)^C(|I|,2) bar z_I T_IJ z_J."""
    t, n = _check_operator(matrix)
    field_start = bar_start + n if field_start is None else field_start
    out = {}
    for i in range(1 << n):
        for j in range(1 << n):
            coefficient = _sign(i.bit_count()) * t[i, j]
            if coefficient == 0:
                continue
            monomial = {0: sp.Integer(1)}
            for mode in range(n):
                if i >> mode & 1:
                    monomial = multiply(monomial, variable(bar_start + mode))
            for mode in range(n):
                if j >> mode & 1:
                    monomial = multiply(monomial, variable(field_start + mode))
            out = add(out, scale(monomial, coefficient))
    return clean(out)


def matrix_from_kernel(poly, modes, bar_start=0, field_start=None):
    """Invert the coefficient map for a polynomial in the standard kernel chart."""
    if not 1 <= modes <= 5:
        raise ValueError("exact kernel reconstruction accepts 1..5 CAR modes")
    field_start = bar_start + modes if field_start is None else field_start
    result = sp.zeros(1 << modes)
    for i in range(1 << modes):
        for j in range(1 << modes):
            mask = sum(1 << (bar_start + k) for k in range(modes) if i >> k & 1)
            mask |= sum(1 << (field_start + k) for k in range(modes) if j >> k & 1)
            result[i, j] = sp.simplify(_sign(i.bit_count()) * poly.get(mask, 0))
    return result


def action_split(matrix, bar_start=0, field_start=None):
    """Return (lambda, S_nil) with lambda*exp(-S_nil)=K exactly."""
    k = kernel(matrix, bar_start, field_start)
    scalar = sp.simplify(k.get(0, 0))
    if scalar == 0:
        raise ValueError("nilpotent action chart requires a nonzero scalar body")
    n = sp.Matrix(matrix).rows.bit_length() - 1
    normalized = scale(k, 1 / scalar)
    nil = dict(normalized)
    nil[0] = sp.simplify(nil.get(0, 0) - 1)
    nil = clean(nil)
    # Each nonzero term has positive Grassmann degree; 2n is a safe finite bound.
    log = {}
    power = nil
    for order in range(1, 2*n + 1):
        if not power:
            break
        coefficient = sp.Rational((-1) ** (order + 1), order)
        log = add(log, scale(power, coefficient))
        power = multiply(power, nil)
    return scalar, scale(log, -1)


def _berezin_integral(poly, modes, internal_start):
    internal_mask = ((1 << (2 * modes)) - 1) << internal_start
    orientation = (-1) ** (modes * (modes + 1) // 2)
    return clean({mask ^ internal_mask: orientation * coefficient
                  for mask, coefficient in poly.items()
                  if mask & internal_mask == internal_mask})


def sew(left, right):
    """Berezin sew K_left(bar x,w) and K_right(bar w,y), returning K_(left right)."""
    a, n = _check_operator(left)
    b, m = _check_operator(right)
    if n != m:
        raise ValueError("sewn kernels must have the same number of modes")
    # Variable blocks: bar-x, bar-w, w, y. The chosen Berezin orientation makes
    # integral exp(-bar-w w)=1, matching joint_source.berezin_source.
    left_kernel = kernel(a, bar_start=0, field_start=2*n)
    right_kernel = kernel(b, bar_start=n, field_start=3*n)
    gaussian = scale(sum_polynomials(
        multiply(variable(n + j), variable(2*n + j)) for j in range(n)), -1)
    integrand = exponential(gaussian, 4*n)
    integrand = multiply(multiply(left_kernel, integrand), right_kernel)
    integrated = _berezin_integral(integrand, n, n)
    # Collapse the external y field block from its temporary 3n offset back
    # to the canonical kernel's n offset.
    result = {}
    for mask, coefficient in integrated.items():
        bars = mask & ((1 << n) - 1)
        fields = (mask >> (3*n)) & ((1 << n) - 1)
        result[bars | (fields << n)] = coefficient
    return clean(result)


def sum_polynomials(polys):
    result = {}
    for poly in polys:
        result = add(result, poly)
    return result


def coherent_trace(matrix, insertion=None):
    """Tr(T W) via integral exp(-bar z z) K_(T W)(-bar z,z)."""
    t, n = _check_operator(matrix)
    w = sp.eye(t.rows) if insertion is None else sp.Matrix(insertion)
    w, nw = _check_operator(w)
    if n != nw:
        raise ValueError("trace insertion must have the same Fock dimension")
    tw = t * w
    k = kernel(tw, bar_start=0, field_start=n)
    substituted = {mask: coefficient * (-1) ** (mask & ((1 << n) - 1)).bit_count()
                   for mask, coefficient in k.items()}
    gaussian = scale(sum_polynomials(
        multiply(variable(j), variable(n + j)) for j in range(n)), -1)
    integral = multiply(exponential(gaussian, 2*n), substituted)
    return sp.simplify(_berezin_integral(integral, n, 0).get(0, 0))


def quarter_operator(modes, charges=None):
    """Gamma(diag(i**q_j)) in the existing source_sewing Fock ordering.

    The default is the producer's original i^N action (all q_j=+1). Supply
    explicit charges for other controls, such as a particle-hole quarter turn.
    """
    if not 1 <= modes <= 5:
        raise ValueError("quarter operator supports 1..5 CAR modes")
    charges = (1,) * modes if charges is None else tuple(charges)
    if len(charges) != modes or any(type(q) is not int for q in charges):
        raise ValueError("quarter charges must be one integer per mode")
    one_particle = sp.diag(*(sp.I**q for q in charges))
    return source_lift(one_particle)


def controls():
    """Small exact controls requested for this conditional reconstruction."""
    identity = sp.eye(4)
    assert kernel(identity) == exponential(
        sum_polynomials(multiply(variable(j), variable(2+j)) for j in range(2)), 4)
    generic = sp.Matrix([[2, 0, 0, 3], [0, 5, 7, 0],
                         [0, 11, 13, 0], [17, 0, 0, 19]])
    assert matrix_from_kernel(kernel(generic), 2) == generic
    a = sp.eye(4)
    b = sp.eye(4)
    a[0, 3] = 1
    b[3, 0] = 1
    assert a*b != b*a
    assert sew(a, b) == kernel(a*b)
    r = sp.diag(sp.Rational(1, 4), sp.Rational(1, 64))
    gaussian_fock = source_lift(r)
    body, action = action_split(gaussian_fock)
    assert body == 1 and all(mask.bit_count() <= 2 for mask in action)
    interacting = sp.diag(1, sp.Rational(1, 4), sp.Rational(1, 64), sp.Rational(1, 65536))
    ibody, iaction = action_split(interacting)
    quartic_masks = [mask for mask in iaction if mask.bit_count() == 4]
    assert ibody == 1 and quartic_masks
    assert iaction[quartic_masks[0]] == -sp.Rational(255, 65536)
    assert coherent_trace(interacting) == sp.Rational(82945, 65536)
    quarter = quarter_operator(2, charges=(1, -1))
    assert coherent_trace(interacting, quarter) == sp.Rational(65537, 65536) + sp.I*sp.Rational(15, 64)
    ell = 8 * sp.log(2)
    h_alt = sp.diag(0, sp.Rational(1, 4), sp.Rational(3, 4), 2)
    assert interacting == sp.diag(*(sp.exp(-ell*h_alt[j, j]) for j in range(4)))
    assert sp.diag(*(sp.exp(2*sp.pi*sp.I*h_alt[j, j]) for j in range(4))) == quarter
    sea_scalar = sp.Pow(2, sp.Rational(1, 12))
    scalar_body, scalar_action = action_split(sea_scalar * interacting)
    assert scalar_body == sea_scalar and scalar_action == iaction
    assert scale(kernel(sea_scalar * interacting), 1 / sea_scalar) == kernel(interacting)
    h_total = h_alt - sp.eye(4) / 96
    expected_loop = sp.exp(-sp.I * sp.pi / 48) * quarter
    actual_loop = sp.diag(*(sp.exp(2*sp.pi*sp.I*h_total[j, j]) for j in range(4)))
    assert all(sp.simplify(actual_loop[j, j] - expected_loop[j, j]) == 0 for j in range(4))
    assert kernel(gaussian_fock) != kernel(interacting)
    assert gaussian_fock[:, :3] == interacting[:, :3]
    return {
        "status": "EXACT_CONDITIONAL_EVEN_CAR_KERNEL_RECONSTRUCTION",
        "identity_kernel_is_exp_bars_fields": True,
        "generic_even_two_mode_roundtrip": True,
        "noncommuting_chronological_sewing": True,
        "gaussian_action_is_bilinear": True,
        "interacting_quartic_action_nonzero": [str(iaction[m]) for m in quartic_masks],
        "quarter_time_equation_exact": True,
        "sea_scalar_retained_in_body": str(scalar_body),
        "full_hamiltonian_loop_phase_exact": True,
        "interacting_ordinary_trace": str(coherent_trace(interacting)),
        "quarter_inserted_trace": str(coherent_trace(interacting, quarter)),
        "same_low_columns_preserve_distinct_extension": True,
        "P1_P2_source_selection_proved": False,
        "scope": "Exact finite coherent-kernel chart for parity-even operators, one to five modes; no physical source selection or continuum claim.",
    }


_RATIONAL = re.compile(r"[+-]?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?\Z")


def _fraction_string(value):
    if not isinstance(value, str) or not _RATIONAL.fullmatch(value):
        raise ValueError("Gaussian rational components must be integer/fraction strings")
    return sp.Rational(Fraction(value).numerator, Fraction(value).denominator)


def _parse_sparse_diagonal(payload):
    if payload.get("schema") != "tfpt.positive-word-kernel.finite-completion.v1":
        raise ValueError("unexpected completion certificate schema")
    if payload.get("ok") is not True or payload.get("complete_finite_kernel_certificate") is not True:
        raise ValueError("completion certificate flags do not authorize this adapter")
    if payload.get("m") != 5 or payload.get("dimension") != 32:
        raise ValueError("adapter accepts only the certified five-mode source-clustering files")
    entries = payload.get("T", {}).get("entries")
    if not isinstance(entries, list):
        raise ValueError("missing sparse Gaussian-rational T entries")
    matrix = sp.zeros(32)
    seen = set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("malformed sparse matrix entry")
        row, column = entry.get("row"), entry.get("column")
        if type(row) is not int or type(column) is not int or not (0 <= row < 32 and 0 <= column < 32):
            raise ValueError("matrix indices must be integers in [0,31]")
        if (row, column) in seen:
            raise ValueError("duplicate sparse matrix entry")
        seen.add((row, column))
        matrix[row, column] = (_fraction_string(entry.get("re"))
                               + sp.I * _fraction_string(entry.get("im")))
    # Restrict the adapter to the actual exact diagonal certificates, where
    # Hermiticity and positivity can be checked without trusting a stored flag.
    if any(matrix[i, j] != 0 for i in range(32) for j in range(32) if i != j):
        raise ValueError("adapter currently accepts diagonal exact T certificates only")
    if any(sp.im(matrix[i, i]) != 0 or sp.re(matrix[i, i]) < 0 for i in range(32)):
        raise ValueError("diagonal T must be exactly Hermitian positive semidefinite")
    dense = payload.get("exact_matrices", {}).get("T")
    if not isinstance(dense, list) or len(dense) != 32 or any(not isinstance(row, list) or len(row) != 32 for row in dense):
        raise ValueError("missing full exact_matrices.T matrix")
    for i in range(32):
        for j in range(32):
            cell = dense[i][j]
            if not isinstance(cell, dict):
                raise ValueError("malformed exact dense T matrix")
            def part(pair):
                if (not isinstance(pair, list) or len(pair) != 2
                        or any(type(x) is not int for x in pair) or pair[1] <= 0):
                    raise ValueError("dense rational entries must be [numerator, positive_denominator]")
                return sp.Rational(pair[0], pair[1])
            value = part(cell.get("re")) + sp.I * part(cell.get("im"))
            if value != matrix[i, j]:
                raise ValueError("sparse T disagrees with exact dense matrix")
    return matrix


def completion_report(path):
    source = Path(path)
    raw = source.read_bytes()
    payload = json.loads(raw.decode("utf-8"))
    matrix = _parse_sparse_diagonal(payload)
    n = 5
    k = kernel(matrix)
    reconstructed = matrix_from_kernel(k, n)
    omega = payload.get("Omega_basis_index")
    if type(omega) is not int or not 0 <= omega < 32:
        raise ValueError("invalid inferred reference basis index")
    r_data = payload.get("R", {})
    exponents = r_data.get("i_exponents")
    if r_data.get("format") != "diagonal roots of unity" or not isinstance(exponents, list) or len(exponents) != 32:
        raise ValueError("missing normalized quarter implementer data")
    if any(type(q) is not int or q not in range(4) for q in exponents):
        raise ValueError("quarter implementer exponents must lie in Z/4")
    normalized_quarter = sp.diag(*(sp.I**q for q in exponents))
    original_data = payload.get("original_internal_W", {})
    original_exponents = original_data.get("i_exponents")
    if (original_data.get("format") != "diagonal roots of unity"
            or not isinstance(original_exponents, list) or len(original_exponents) != 32):
        raise ValueError("missing original internal quarter action")
    if any(type(q) is not int or q not in range(4) for q in original_exponents):
        raise ValueError("original quarter exponents must lie in Z/4")
    expected_original = [i.bit_count() % 4 for i in range(32)]
    if original_exponents != expected_original:
        raise ValueError("original_internal_W is not the recorded all-positive i^N field action")
    omega_phase = original_exponents[omega]
    expected_normalized = [(q - omega_phase) % 4 for q in original_exponents]
    if exponents != expected_normalized:
        raise ValueError("normalized R does not equal original i^N divided by its Omega phase")
    original_quarter = sp.diag(*(sp.I**q for q in original_exponents))
    if matrix * normalized_quarter != normalized_quarter * matrix:
        raise ValueError("T fails exact covariance with the certificate's normalized quarter operator")
    if matrix * original_quarter != original_quarter * matrix:
        raise ValueError("T fails exact covariance with the certificate's original i^N operator")
    if exponents[omega] != 0:
        raise ValueError("normalized quarter operator does not fix inferred reference vector")
    trace = coherent_trace(matrix)
    original_quarter_trace = coherent_trace(matrix, original_quarter)
    normalized_quarter_trace = coherent_trace(matrix, normalized_quarter)
    eigenvalues = [sp.re(matrix[i, i]) for i in range(32)]
    strictly_positive = all(value > 0 for value in eigenvalues)
    if matrix[0, 0] == 0:
        body, action = sp.Integer(0), {}
        action_status = "unavailable: empty-vacuum body T[0,0] is zero"
        body_text = "0"
        action_coefficients = []
        action_degree = None
    else:
        body, action = action_split(matrix)
        action_status = "exact nilpotent logarithm; scalar retained separately"
        body_text = str(body)
        action_coefficients = [{"mask": mask, "coefficient": str(value)} for mask, value in sorted(action.items())]
        action_degree = max((mask.bit_count() for mask in action if mask), default=0)
    closed_form = None
    if payload.get("objective") == "min-trace":
        vacuum_projector = sp.zeros(32)
        vacuum_projector[omega, omega] = 1
        minimum_candidate = (sp.eye(32) + 3 * vacuum_projector) / 4
        bilinear_part = {}
        for mode in range(n):
            bilinear_part = add(bilinear_part,
                                multiply(variable(mode), variable(n + mode)))
        projector_kernel = kernel(vacuum_projector)
        expected_action = add(scale(bilinear_part, -1), scale(projector_kernel, -3))
        closed_form = {
            "label": "min-trace search-objective artifact; not physical vertex selection",
            "matrix_formula": "T=(I+3|31><31|)/4",
            "kernel_formula": "K=exp(B)/4+3M/4",
            "action_formula_when_T00_nonzero": "S_nil=-B-3M",
            "reference_index_is_31": omega == 31,
            "matrix_matches": omega == 31 and matrix == minimum_candidate,
            "kernel_matches": k == add(scale(kernel(sp.eye(32)), sp.Rational(1, 4)),
                                        scale(projector_kernel, sp.Rational(3, 4))),
            "nilpotent_action_matches": matrix[0, 0] != 0 and action == expected_action,
        }
    return {
        "schema": "tfpt.coherent-kernel.functional-reconstruction.v1",
        "input": str(source),
        "input_file_sha256": hashlib.sha256(raw).hexdigest(),
        "completion_input_sha256": payload.get("input_sha256"),
        "source_certificate_schema": payload["schema"],
        "source_verdict": payload.get("verdict"),
        "source_regularization_and_assumptions": {
            "objective_scope": payload.get("objective_scope"),
            "inferred_occupations": payload.get("inferred_occupations"),
            "Omega_basis_index": omega,
            "R": r_data,
            "original_internal_W": payload.get("original_internal_W"),
            "time_scope": payload.get("time_scope"),
            "extra_realization_ansatz": payload.get("extra_realization_ansatz"),
        },
        "finite_matrix_checks": {
            "diagonal_hermitian_psd_checked_exactly": True,
            "sparse_dense_matrix_reconstruction": reconstructed == matrix,
            "normalized_quarter_covariance": True,
            "original_internal_quarter_covariance": True,
            "normalized_R_matches_original_divided_by_Omega_phase": True,
            "T_singular": not strictly_positive,
            "finite_selfadjoint_H_minus_log_T_exists": strictly_positive,
        },
        "kernel_coefficients": [{"mask": mask, "coefficient": str(value)} for mask, value in sorted(k.items())],
        "empty_vacuum_body_T_00": str(matrix[0, 0]),
        "inferred_reference_eigenvalue": str(matrix[omega, omega]),
        "scalar_body_lambda": body_text,
        "nilpotent_action_status": action_status,
        "nilpotent_action_degree": action_degree,
        "nilpotent_action_coefficients": action_coefficients,
        "min_trace_closed_form": closed_form,
        "roundtrip_exact": reconstructed == matrix and (matrix[0, 0] == 0 or
             scale(exponential(scale(action, -1), 2*n), body) == k),
        "ordinary_trace": str(trace),
        "original_internal_quarter_trace": str(original_quarter_trace),
        "normalized_GNS_quarter_trace": str(normalized_quarter_trace),
        "quarter_time_equation_derived": False,
        "quarter_inserted_trace": str(original_quarter_trace),
        "quarter_inserted_trace_implementer": "original_internal_W = i^N",
        "physical_source_selection_proved": False,
        "scope": "Conditional five-mode finite CAR completion; source-clustering certificate, with no P1/P2 selection or continuum claim.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--completion", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.completion is None:
        raise SystemExit("--completion PATH is required")
    result = completion_report(args.completion)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
