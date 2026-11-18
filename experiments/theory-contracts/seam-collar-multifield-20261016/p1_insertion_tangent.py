#!/usr/bin/env python3
"""Test the existing covariance reader and E8 insertion action; no fitted vertex.

Covariance matrices below are derivative controls, not physical P1 sources.
E8 charges and correlators come from the existing stipulated current dictionary.
The analytic scope and the missing source jet are stated in P1-INSERTION-TANGENT.txt.
"""
from __future__ import annotations

from collections import Counter
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

import e8_boundary_response as e8


def existing_induction():
    repo = Path(__file__).resolve().parents[3]
    sys.path.insert(0, str(repo / "verification"))
    spec = importlib.util.spec_from_file_location(
        "existing_v258", repo / "verification/v258_dirac_covariance_induction.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.induce


def logit_tangent(covariance, variation):
    """Frechet derivative of h(C)=log(1-C)-log(C), including repeated eigenvalues."""
    c, u = np.linalg.eigh(covariance)
    if np.any(c <= 0) or np.any(c >= 1):
        raise ValueError("Requires a faithful finite covariance, 0<C<I")
    h = np.log1p(-c) - np.log(c)
    coeff = np.empty((len(c), len(c)))
    for i in range(len(c)):
        for j in range(len(c)):
            if abs(c[i] - c[j]) < 1e-12:
                coeff[i, j] = -1 / (c[i] * (1 - c[i]))
            else:
                coeff[i, j] = (h[i] - h[j]) / (c[i] - c[j])
    v = u.conj().T @ variation @ u
    return u @ (coeff * v) @ u.conj().T, coeff


def covariance_checks():
    induce = existing_induction()
    # Distinct and repeated eigenvalues, then a noncommuting complex frame.
    unitary = np.array([[1, 1, 1, 1], [1, 1, -1, -1],
                        [1, -1, 1, -1], [1, -1, -1, 1]], complex) / 2
    unitary = np.diag([1, 1j, -1, -1j]) @ unitary
    variation = np.array([[1, 2j, 3, 1j], [-2j, -1, 1+1j, 2],
                          [3, 1-1j, 2, -1j], [-1j, 2, 1j, -2]], complex) / 20
    errors = []
    all_negative = True
    step = 1e-5
    for eigenvalues in ([.2, .3, .6, .8], [.2, .2, .6, .6]):
        cov = unitary @ np.diag(eigenvalues) @ unitary.conj().T
        analytic, coeff = logit_tangent(cov, variation)
        numerical = (induce(cov + step*variation) - induce(cov - step*variation))/(2*step)
        errors.append(float(np.max(np.abs(analytic-numerical))))
        all_negative &= bool(np.all(coeff < 0))
    # All 16 Hermitian basis directions: exact zero blocks remain zero.
    cov = np.diag([.2, .3, .6, .8])
    directions = []
    for i in range(4):
        a = np.zeros((4, 4), complex)
        a[i, i] = 1
        directions.append(a)
        for j in range(i):
            for phase in (1, 1j):
                a = np.zeros((4, 4), complex)
                a[i, j], a[j, i] = phase, np.conj(phase)
                directions.append(a)
    for a in directions:
        da, _ = logit_tangent(cov, a)
        assert np.array_equal(np.abs(a) > 0, np.abs(da) > 0)
    # Moving carrier embedding: even fixed ambient C can have nonzero compressed jet.
    cross = np.diag([1, -1, 2, -2]) / 40
    ambient = np.block([[.5*np.eye(4), cross], [cross, .5*np.eye(4)]])
    q0 = np.vstack([np.eye(4), np.zeros((4, 4))])
    q1 = np.vstack([np.zeros((4, 4)), np.eye(4)])
    def compressed(t):
        q = np.cos(t)*q0 + np.sin(t)*q1
        return q.T @ ambient @ q
    exact = q1.T@ambient@q0 + q0.T@ambient@q1
    moving_error = float(np.max(np.abs((compressed(step)-compressed(-step))/(2*step)-exact)))
    assert max(errors) < 1e-8 and all_negative and moving_error < 1e-8
    return {"existing_reader": "verification/v258_dirac_covariance_induction.py:induce",
            "noncommuting_derivative_max_errors": errors,
            "all_divided_difference_coefficients_negative": all_negative,
            "hermitian_support_directions": len(directions),
            "moving_embedding_derivative_max_error": moving_error,
            "scope": "Derivative controls; no P1 covariance, embedding or Higgs vertex selected"}


def retained_collar_checks():
    """Actual previously retained particle/hole transfer eigenvalues, no new coupling."""
    sp = e8.sp
    r = sp.diag(sp.Rational(1, 4), sp.Rational(1, 64))
    eye = sp.eye(2)
    a = (eye+r*r).inv()
    graph = eye.col_join(r)
    projection = graph*a*graph.T
    assert projection*projection == projection
    assert projection*graph == graph
    # Restriction to its own range is the identity, hence not faithful mixed C.
    q = graph*sp.diag(sp.sqrt(sp.Rational(16, 17)), sp.sqrt(sp.Rational(4096, 4097)))
    assert q.T*q == eye and sp.simplify(q.T*projection*q) == eye
    occupation = eye-a
    h = sp.diag(*(sp.log((1-x)/x) for x in occupation.diagonal()))
    ell = 8*sp.log(2)
    b = sp.diag(sp.Rational(1, 4), sp.Rational(3, 4))
    assert sp.simplify(h-2*ell*b) == sp.zeros(2)
    return {"transfer_eigenvalues": ["1/4", "1/64"],
            "calderon_projection_exact": True,
            "own_range_covariance": "I (not in the domain of finite logit induction)",
            "mixed_one_face_occupations": list(map(str, occupation.diagonal())),
            "one_face_modular_generator": list(map(str, h.diagonal())),
            "exact_recovery": "h(1-A)=2*ell*B, ell=8*log(2), B=diag(1/4,3/4)",
            "scope": "Retained conditional collar, not a physical P1 field family"}


def current_action_checks():
    roots = e8.roots()
    rootset = set(roots)
    spinors = tuple(p for p in roots if e8.grade(p) == 1)
    targets = set(p for p in roots if e8.grade(p) == 3)
    mediators = tuple(r for r in roots if e8.grade(r) == 2)
    outgoing, ranks, family_edges = Counter(), Counter(), Counter()
    values = Counter()
    sample = None
    for r in mediators:
        images = set()
        for p in spinors:
            s = e8.add(r, p)
            allowed = s in rootset
            assert allowed == (e8.dot(r, p) == -1)
            if not allowed:
                continue
            assert s in targets
            # Actual existing three-current reader, not an inserted Yukawa matrix.
            response = e8.current_response((r, p, e8.neg(s)), (3, 2, 1))
            norm = e8.current_response((s, e8.neg(s)), (2, 1))
            coefficient = 2*response/norm
            assert coefficient == e8.epsilon(r, p)
            assert response != 0
            source_family = p[5:]
            target_conjugate_family = e8.neg(s)[5:]
            assert source_family != target_conjugate_family
            family_edges[(source_family, target_conjugate_family)] += 1
            outgoing[p] += 1
            images.add(s)
            values[str(response)] += 1
            if sample is None:
                sample = {"r": list(map(str, r)), "p": list(map(str, p)),
                          "target": list(map(str, s)), "response_3_2_1": str(response),
                          "zero_mode_coefficient": str(coefficient)}
        ranks[len(images)] += 1
    # Check every possible same-family pairing against all existing grade-two currents.
    same_family_pairs = [(p, q) for p in spinors for q in spinors if p[5:] == q[5:]]
    for p, q in same_family_pairs:
        assert e8.neg(e8.add(p, q)) not in set(mediators)
    # One direct existing-reader zero per same-family pair, with each mediator sampled cyclically.
    for i, (p, q) in enumerate(same_family_pairs):
        assert e8.current_response((mediators[i % len(mediators)], p, q), (3, 2, 1)) == 0
    assert len(spinors) == 64 and len(mediators) == 60
    assert set(outgoing.values()) == {15} and sum(outgoing.values()) == 960
    assert ranks == {16: 60} and len(family_edges) == 12
    assert set(family_edges.values()) == {80}
    return {"spinor_current_dimension_each": 64, "grade_two_currents": 60,
            "candidate_root_pairs": len(spinors)*len(mediators),
            "nonzero_root_actions": sum(outgoing.values()),
            "outgoing_per_spinor": 15, "rank_each_root_action": 16,
            "ordered_distinct_family_channels": len(family_edges),
            "root_actions_per_family_channel": 80,
            "same_family_pairs_without_any_grade_two_root": len(same_family_pairs),
            "three_point_values": dict(values), "sample": sample,
            "scope": "Internal E8 current action; bosonic currents, not physical fermions or selected Phi"}


def main():
    result = {"verdict": "SOURCE_JET_NOT_DERIVED_CONDITIONAL_ACTION_EXPLICIT",
              "covariance": covariance_checks(), "retained_collar": retained_collar_checks(),
              "e8_action": current_action_checks(),
              "raw_P1_Phi_map_derived": False, "physical_gates_closed": []}
    if len(sys.argv) > 1:
        Path(sys.argv[1]).write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
