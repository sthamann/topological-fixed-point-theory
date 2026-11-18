"""Exact controls for two source-to-carrier implications; no physical model.

The mathematical proofs and source scope are in PROOF.txt. Integer arithmetic
checks tensor invariance on a Lie-algebra generating set, not a numerical fit.
Run from this directory or any working directory; the report is printed as JSON.
"""

import hashlib
import importlib.util
import itertools
import json
from pathlib import Path


def wedge_action(indices, row, col):
    """Action of the matrix unit E_row,col on an ordered exterior basis."""
    terms = []
    for position, item in enumerate(indices):
        if item != col:
            continue
        changed = list(indices)
        changed[position] = row
        if len(set(changed)) != len(changed):
            continue
        inversions = sum(changed[i] > changed[j]
                         for i in range(len(changed))
                         for j in range(i + 1, len(changed)))
        terms.append((tuple(sorted(changed)), (-1) ** inversions))
    return terms


def volume(first, rest, dimension):
    word = (first,) + rest
    if sorted(word) != list(range(dimension)):
        return 0
    inversions = sum(word[i] > word[j]
                     for i in range(dimension)
                     for j in range(i + 1, dimension))
    return (-1) ** inversions


def epsilon(first, second):
    return 0 if first == second else (1 if (first, second) == (0, 1) else -1)


def invariant_control(b):
    """(V tensor W) tensor Lambda^(b-1)V tensor W -> det V tensor det W."""
    exterior = list(itertools.combinations(range(b), b - 1))
    coordinates = list(itertools.product(range(b), range(2), exterior, range(2)))

    def tensor(i, j, subset, h):
        return volume(i, subset, b) * epsilon(j, h)

    residuals = []
    # All matrix units of sl(V), plus diagonal differences.
    for row in range(b):
        for col in range(b):
            if row == col:
                continue
            for i, j, subset, h in coordinates:
                value = tensor(row, j, subset, h) if i == col else 0
                value += sum(sign * tensor(i, j, changed, h)
                             for changed, sign in wedge_action(subset, row, col))
                residuals.append(value)
    for k in range(b - 1):
        for i, j, subset, h in coordinates:
            weight = ((i == k) - (i == k + 1)
                      + (k in subset) - (k + 1 in subset))
            residuals.append(weight * tensor(i, j, subset, h))
    # sl(W): the two elementary roots and its diagonal generator.
    for row, col in [(0, 1), (1, 0)]:
        for i, j, subset, h in coordinates:
            value = tensor(i, row, subset, h) if j == col else 0
            value += tensor(i, j, subset, row) if h == col else 0
            residuals.append(value)
    for i, j, subset, h in coordinates:
        weight = (j == 0) - (j == 1) + (h == 0) - (h == 1)
        residuals.append(weight * tensor(i, j, subset, h))
    # Primitive central cocharacter: (-2 on V, b on W), since b is odd.
    charges = [b - 2, -2 * (b - 1), b]
    central_residual = sum(charges)
    assert central_residual == 0
    assert all(value == 0 for value in residuals)
    nonzero = sum(tensor(*args) != 0 for args in coordinates)
    assert nonzero == 2 * b
    lower_neutral = []
    for degree in (1, 2):
        for powers in itertools.product(range(degree + 1), repeat=3):
            if sum(powers) == degree and sum(x * q for x, q in zip(powers, charges)) == 0:
                lower_neutral.append(powers)
    assert not lower_neutral
    assert (b - 1) % 2 == 0  # both fermion legs are in the even exterior packet
    # At total V-degree d<b, the SU(b) center exp(2pi i/b) acts nontrivially.
    # At degree b, the determinant is allowed; central neutrality then needs W-degree 2.
    smaller_positive_carrier_degrees = []
    for degree_v in range(b + 2):
        for degree_w in range(b + 2 - degree_v):
            total = degree_v + degree_w
            if (0 < total < b + 2 and degree_v % b == 0 and degree_w % 2 == 0
                    and -2 * degree_v + b * degree_w == 0):
                smaller_positive_carrier_degrees.append((degree_v, degree_w))
    assert not smaller_positive_carrier_degrees
    return {
        "rank_minus": b,
        "rank_plus": 2,
        "fermion_exterior_degrees": [2, b - 1],
        "field_dimensions": [2 * b, b, 2],
        "tensor_entries": len(coordinates),
        "nonzero_entries": nonzero,
        "exact_lie_residuals": len(residuals),
        "max_abs_lie_residual": max(map(abs, residuals)),
        "central_charges": charges,
        "central_residual": central_residual,
        "lower_field_degree_neutral_monomials": lower_neutral,
        "lower_positive_carrier_degree_candidates": smaller_positive_carrier_degrees,
        "ordinary_polynomial_indecomposable": True,
        "minimum_positive_carrier_degree": b + 2,
        "full_even_packet_dimension": 2 ** (b + 1),
        "claim_scope": "representation inference only; not an admissible P1/P2 source",
    }


def cp1_h0(degree):
    return max(degree + 1, 0)


def cp1_h1(degree):
    return max(-degree - 1, 0)


def source_vertex_control(root):
    """Reuse actual Berezin integration; supplied deformations are diagnostics."""
    import sympy as s
    path = root / "experiments/theory-contracts/seam-collar-multifield-20261016/joint_source.py"
    spec = importlib.util.spec_from_file_location("tfpt_existing_joint_source", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    phi = s.Symbol("phi", real=True)
    # Actual lowest charged particle/hole closure from the retained quarter source.
    k0 = s.diag(s.Rational(17, 16), s.Rational(4097, 4096))
    probes = [("constant_kinetic", s.zeros(2)),
              ("supplied_offdiagonal_deformation", s.Matrix([[0, 1], [1, 0]]))]
    records = []
    background_coefficients = []
    for name, deformation in probes:
        # This is the existing literal Grassmann integration, not an inverse formula.
        polynomial = module.berezin_source(k0 + phi * deformation)
        z0 = polynomial[0]
        response = s.Matrix(2, 2, lambda i, j:
                            s.cancel(polynomial.get((1 << i) | (1 << (2 + j)), 0) / z0))
        assert s.simplify(response - (k0 + phi * deformation).inv()) == s.zeros(2)
        background_coefficients.append({mask: s.simplify(value.subs(phi, 0))
                                        for mask, value in polynomial.items()})
        vertex = response.diff(phi).subs(phi, 0)
        assert s.simplify(vertex + k0.inv() * deformation * k0.inv()) == s.zeros(2)
        amputated = s.simplify(-k0 * vertex * k0)
        assert amputated == deformation
        records.append({"diagnostic": name,
                        "partition_at_zero": str(z0.subs(phi, 0)),
                        "mixed_source_derivative": [[str(v) for v in row] for row in vertex.tolist()],
                        "amputated_vertex": [[str(v) for v in row] for row in amputated.tolist()]})
    keys = set(background_coefficients[0]) | set(background_coefficients[1])
    assert all(background_coefficients[0].get(key, 0) == background_coefficients[1].get(key, 0)
               for key in keys)
    assert records[0]["mixed_source_derivative"] != records[1]["mixed_source_derivative"]
    return {"same_complete_zero_phi_source": True,
            "same_absolute_partition_at_zero": True,
            "distinct_mixed_source_derivatives": True,
            "controls": records,
            "source_coupling_derived": False,
            "scope": "Exact identifiability control on the existing free source, not two P1/P2 worlds."}


def main():
    root = Path(__file__).resolve().parents[3]
    sources = [
        "_archive/tfpt-45/source_extracts/01_boundary_kernel_source.tex",
        "_archive/tfpt-45/source_extracts/02_carrier_source.tex",
        "experiments/lean4-carrier-rigidity/TfptCarrier/HiggsIndexShadow.lean",
        "experiments/lean4-carrier-rigidity/TfptCarrier/YukawaPrimitive.lean",
        "experiments/theory-contracts/seam-collar-multifield-20261016/joint_source.py",
    ]
    hashes = {path: hashlib.sha256((root / path).read_bytes()).hexdigest() for path in sources}
    bundle_controls = []
    for degrees in ([1], [1, 0], [1, 0, 0], [3, -2]):
        h0 = sum(map(cp1_h0, degrees))
        h1 = sum(map(cp1_h1, degrees))
        assert sum(degrees) == 1
        assert h0 - h1 == len(degrees) + 1
        bundle_controls.append({"splitting_degrees": degrees, "determinant_degree": 1,
                                "fiber_rank": len(degrees), "h0": h0, "h1": h1,
                                "dolbeault_index": h0 - h1})
    report = {
        "contract": "UR.SEAM.JOINT_SOURCE_ORIGIN_AUDIT.01",
        "verdict": "PARTIAL_SOURCE_DERIVATION_NOT_CLOSED",
        "source_sha256": hashes,
        "higgs_determinant_controls": bundle_controls,
        "yukawa_controls": [invariant_control(3), invariant_control(5)],
        "existing_source_vertex_control": source_vertex_control(root),
        "excluded_implication": "invariant-ring primitiveness alone selects b=3",
        "not_excluded": "a source-derived stronger field-degree or global-minimality selection",
        "physical_source_derived": False,
        "p1_p2_countermodel_constructed": False,
        "lean_rebuilt": False,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
