"""Evaluate the common exterior probe on the two conjugate NS modes."""
import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import block_diag, eigh, solve


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with contextlib.redirect_stdout(io.StringIO()):
    probe = runpy.run_path(str(base / "check_original_interference_cut_20261109.py"))

old = probe["old"]
identity, swap, gamma = [probe[name] for name in ("identity", "swap", "gamma")]
chains, block, read = [probe[name] for name in ("chains", "block", "phase_read")]
trace, center = [probe[name] for name in ("exact_from_chains", "center")]


def conjugate(a):
    return swap @ a.conj() @ swap


def half(radius, sign, reference=False, reverse=False):
    """Recompute both generators using the same recorded Higgs and link history."""
    k = identity.copy()
    sea = 0.
    segments = list(enumerate(old["durations"][:3]))
    for i, width in reversed(segments) if reverse else segments:
        _, phi, p = old["existing"]["generator"](i, "S" if reference else "A", radius)
        h = sign * p * gamma - 1j * gamma @ phi
        eigenvalues, u = eigh(h)
        k = ((u * np.exp(-width * eigenvalues)) @ u.conj().T) @ k
        sea += width * eigenvalues[:len(eigenvalues)//2].sum()
    return k, sea


def evaluate(original, sea, reference, ref_sea, factors, marks):
    """One simultaneous color selection on the entire exterior algebra."""
    det_globals = probe["det_mark"].__globals__
    det_globals["identity"] = np.eye(len(original))
    try:
        denominator = trace([], [(reference, False)], factors, ref_sea)
        assert abs(denominator.imag) < 1e-13 * denominator.real and denominator.real > 0
        shift = 2 * sea - np.log(denominator.real)
        gram = block(chains(original, *marks), factors, shift)
        scale = np.max(abs(gram))
        assert np.max(abs(gram - gram.conj().T)) < 1e-11 * scale
        assert min(np.linalg.eigvalsh(gram)) > 0
        empty = trace([(original, False)], [(original, False)], factors, shift)
        return gram, denominator.real, empty.real
    finally:
        det_globals["identity"] = identity


def cauchy_graph(k):
    graph = np.vstack((identity, k))
    return graph @ solve(graph.conj().T @ graph, graph.conj().T, assume_a="pos")


marks = (probe["yukawa"](probe["cut"], "e"),
         probe["yukawa"](probe["first"], "nu"),
         probe["mark"](probe["weak"]))
assert np.max(abs(conjugate(gamma) + gamma)) == 0
for a in (*marks, *center, old["source"]["majorana"]):
    assert np.max(abs(conjugate(a) - a)) < 1e-13
for a, parity in zip(marks, (-1, -1, 1)):
    assert np.max(abs(gamma @ conjugate(a) @ gamma - parity * a)) < 1e-13

rows = {}
for radius in (-.001, .001):
    signs = (1, -1)
    halves = [half(radius, sign) for sign in signs]
    refs = [half(radius, sign, reference=True) for sign in signs]
    k, s = zip(*halves)
    ref_half, ref_half_sea = zip(*refs)
    refs_whole = [r.conj().T @ r for r in ref_half]
    assert np.max(abs(k[1] - conjugate(k[0]))) < 1e-12
    assert np.max(abs(refs_whole[1] - conjugate(refs_whole[0]))) < 1e-12
    assert abs(s[0] - s[1]) < 1e-12
    face_swap = block_diag(swap, swap)
    graph_defect = np.max(abs(
        face_swap @ cauchy_graph(k[0]).conj() @ face_swap - cauchy_graph(k[1])))
    assert graph_defect < 1e-12
    # A literal face exchange turns the graph of a transfer into the graph
    # of its inverse. Conjugation on the same two faces did not do that.
    literal_face_exchange = np.linalg.norm(conjugate(k[0]) - solve(k[1], identity))
    assert literal_face_exchange > 9
    # With gamma included and the ordered history reversed, the graph can
    # be intertwined at a fixed p. This is a distinct operation.
    k_reversed, _ = half(radius, 1, reverse=True)
    swap_faces = np.roll(np.eye(2 * len(identity)), len(identity), axis=0)
    gamma_real = block_diag(gamma @ swap, gamma @ swap)
    action = swap_faces @ gamma_real
    reversed_graph_defect = np.max(abs(
        action @ cauchy_graph(k[0]).conj() @ action.conj().T - cauchy_graph(k_reversed)))
    assert reversed_graph_defect < 1e-12
    ref_sea = 2 * np.array(ref_half_sea)
    singles = [evaluate(k[i], s[i], refs_whole[i], ref_sea[i], center, marks)[0]
               for i in range(2)]
    assert np.max(abs(singles[1] - singles[0].conj())) < 1e-11 * np.max(abs(singles[0]))
    assert abs(singles[0][0, 1].imag) > 4e-10
    assert np.max(abs(singles[0] - singles[1])) > 8e-10

    # For a spatially constant Higgs/mark in this circle chart, the two
    # momentum modes are reducing. The same source parameter acts on both.
    # The exterior of their direct sum is used before the single Z3 average.
    gram, denominator, empty = evaluate(block_diag(*k), sum(s),
        block_diag(*refs_whole), sum(ref_sea),
        [block_diag(z, z) for z in center],
        [block_diag(a, a) for a in marks])
    assert np.max(abs(gram.imag)) < 1e-11 * np.max(abs(gram))
    assert gram[0, 1].real > 0
    # This is a check against simply averaging the two normalized single-mode Grams.
    assert np.max(abs(gram - (singles[0] + singles[1])/2)) > 1e-7
    results = read(gram)
    assert abs(results["I90"] - results["arm_sum"]) < 1e-11 * results["arm_sum"]
    rows[str(radius)] = {
        "positive_p": read(singles[0]),
        "negative_p": read(singles[1]),
        "paired_joint": results,
        "paired_reference": denominator,
        "paired_empty_over_reference": empty,
        "pair_antiunitary_defect_relative": float(np.max(abs(gram.imag)) / np.max(abs(gram))),
        "conditional_cauchy_graph_intertwining_error": float(graph_defect),
        "literal_face_swap_at_same_history_frobenius_error": float(literal_face_exchange),
        "gamma_J_face_swap_on_reversed_history_graph_error": float(reversed_graph_defect),
        "full_physical_R113_I0": None,
        "full_physical_R113_I90": None,
    }

print(json.dumps({
    "scope": "Fixed spatially constant original R032 history. First conjugate NS pair only; 192 one-body modes, full exterior, one common original color-center insertion. Not the original dynamical R113 measure or Padm.",
    "rows": rows,
    "native_gate_status": probe["gate"]["status"],
    "native_gate_missing": probe["gate"]["missing_operations"] + probe["gate"]["fiber_missing_operations"],
}, indent=2))
