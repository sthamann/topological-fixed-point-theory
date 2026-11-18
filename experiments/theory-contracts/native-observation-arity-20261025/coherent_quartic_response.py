#!/usr/bin/env python3
"""Original quartic source and its exact three-register coherent response.

The six directions are the marked E8/J ovoid rays from the October 24
certificate.  This checker independently rebuilds the tensor polynomials
from the original word supports and checks complete polynomial identities.
It does not infer a physical preparation or a state update rule.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools as it
import json
from pathlib import Path
import sys

import sympy as sp


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PRIOR = HERE.parent / "atlas-defect-synthesis-20261024"
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tfpt_explorer.process import _code_basis, _hamming_words  # noqa: E402


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, "cannot load source")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial_json(value, variables):
    terms = sp.Poly(sp.expand(value), *variables).terms()
    return [[list(exponent), str(coefficient)] for exponent, coefficient in terms]


def main(output: Path):
    _, supports = _code_basis()
    sizes = [len(support) for support in supports]
    require(sizes == [4, 12, 12, 12, 24], "original code supports changed")
    bridge_path = PRIOR / "doily_quartic_bridge.py"
    bridge = load_module("doily_for_coherent_response", bridge_path).build_report()
    rays = [sp.Matrix(q) for q in bridge["selection_data"]["native_ovoid_rays"]]
    require(len(rays) == 6, "need the actual six marked rays")
    source_ray_equations = []
    for q in rays:
        a, b, c, d, e = q
        residuals = [a*a - d*d - 4*e*e, b*b - d*d, c*c - d*d,
                     a*b + c*d, a*c + b*d, a*d + b*c,
                     b*e, c*e, d*e]
        require(residuals == [0] * 9, "marked ovoid is not a source product ray")
        source_ray_equations.append(list(map(str, residuals)))
    x = sp.symbols("x0:4")
    y = sp.symbols("y0:4")  # Algebraically independent conjugate variables.
    word_polynomials = [
        sp.expand(sum(sp.prod(x[j] for j in word) for word in support))
        for support in supports
    ]
    phi = sp.Matrix([
        word_polynomials[i] / sp.sqrt(sizes[i]) for i in range(5)
    ])
    diagonal = sp.diag(*[sp.sqrt(sp.Rational(size, 4)) for size in sizes])
    original_W = load_module(
        "erasure_for_coherent_response", PRIOR / "native_erasure_readout.py"
    ).decoder_isometry_gauge_first(supports)
    oldest_erasure_path = (
        HERE.parent / "joint-native-algebra-selection-20261016/operator_checker.py"
    )
    original_logical_first = load_module(
        "original_erasure_for_coherent_response", oldest_erasure_path
    ).erasure_isometry(supports)
    original_column_checks = [
        original_W[:, 5 * gauge + logical]
        == original_logical_first[:, 4 * logical + gauge]
        for gauge, logical in it.product(range(4), range(5))
    ]
    require(all(original_column_checks), "W differs from the original erasure isometry")
    require(original_W.H * original_W == sp.eye(20), "W is not isometric")
    cubic_input = sp.Matrix([
        sp.prod(x[j] for j in word) for word in it.product(range(4), repeat=3)
    ])
    decoded = original_W.H * cubic_input
    derivative_checks = []
    for gauge, logical in it.product(range(4), range(5)):
        residual = sp.expand(
            decoded[5 * gauge + logical] - sp.diff(phi[logical], x[gauge]) / 2
        )
        derivative_checks.append(residual == 0)
    require(all(derivative_checks), "original W is not half the quartic derivative")

    six_quartics = [sp.expand((q.T * diagonal * phi)[0]) for q in rays]
    require(sp.expand(sum(six_quartics)) == 0, "six coordinates do not sum to zero")
    igusa = sp.expand(
        sum(f**2 for f in six_quartics)**2
        - 4 * sum(f**4 for f in six_quartics)
    )
    require(igusa == 0, "marked six-coordinate Igusa relation failed")
    # The scalar in the actual same-source product is also the genus-two
    # complete weight enumerator of the SAME original binary Hamming code.
    # This is an exact polynomial equality, not a numerical coefficient match.
    fock_phi = phi / sp.sqrt(sp.factorial(4))
    scalar_source = sp.expand(sum(f*f for f in fock_phi))
    hamming_words = _hamming_words()
    require(len(hamming_words) == 16, "original Hamming code size changed")
    hamming_tetrads = {
        tuple(i for i, bit in enumerate(word) if bit)
        for word in hamming_words if sum(word) == 4
    }
    ray_tetrads = set()
    for row in bridge["quartic_bridge"]["per_reflection"]:
        # A Gaussian ray is a mu4 orbit, not a single real root. Recover all
        # four original real roots before comparing their binary supports.
        # Taking only the 60 canonical representatives would miss 4 tetrads.
        for quarter_turn in range(4):
            coordinates = []
            for real, imaginary in row["gaussian_ray"]:
                phase_orbit = [(real, imaginary), (-imaginary, real),
                               (-real, -imaginary), (imaginary, -real)]
                coordinates.extend(phase_orbit[quarter_turn])
            support = tuple(i for i, value in enumerate(coordinates) if value)
            if len(support) == 4:
                ray_tetrads.add(support)
    require(len(hamming_tetrads) == 14 and hamming_tetrads == ray_tetrads,
            "enumerator code and actual Construction-A root tetrads differ")
    complete_enumerator = sp.expand(sum(
        sp.prod(x[2*left[k] + right[k]] for k in range(8))
        for left, right in it.product(hamming_words, repeat=2)
    ))
    octic = sp.expand(
        sum(z**8 for z in x)
        + 14*sum(x[i]**4*x[j]**4 for i, j in it.combinations(range(4), 2))
        + 168*sp.prod(z*z for z in x)
    )
    require(complete_enumerator == octic, "actual Hamming pair enumerator differs from octic")
    require(sp.expand(octic - 96*scalar_source) == 0,
            "normalized original source product differs from Hamming enumerator")
    require(sp.expand(octic - sum(f*f for f in six_quartics)/3) == 0,
            "marked six-quartic sum does not equal original Hamming enumerator")
    fock_norm_squared = sp.simplify(sum(
        sp.conjugate(coefficient)*coefficient
        * sp.prod(sp.factorial(exponent) for exponent in powers)
        for powers, coefficient in sp.Poly(scalar_source, *x).terms()
    ))
    require(fock_norm_squared == 140, "actual scalar source Fock normalization changed")
    gradients = [sp.Matrix([sp.diff(f, z) for z in x]) for f in six_quartics]
    gradient_probability_polynomials = []
    direct_probability_polynomials = []
    substitution = dict(zip(x, y))
    for q, gradient in zip(rays, gradients):
        gradient_probability = sp.expand(
            sum(gradient[j] * gradient[j].subs(substitution) for j in range(4)) / 48
        )
        q_orthonormal = diagonal * q
        effect5 = q_orthonormal * q_orthonormal.T / 12
        direct = sp.expand(
            (decoded.subs(substitution).T
             * sp.kronecker_product(sp.eye(4), effect5) * decoded)[0]
        )
        require(direct == gradient_probability, "gradient and Born rule disagree")
        gradient_probability_polynomials.append(gradient_probability)
        direct_probability_polynomials.append(direct)
    norm = sum(x[j] * y[j] for j in range(4))
    require(
        sp.expand(sum(gradient_probability_polynomials) - norm**3) == 0,
        "six probabilities are not normalized on the entire coherent source",
    )

    e0_substitution = {x[0]: 1, y[0]: 1}
    e0_substitution.update({z: 0 for z in (*x[1:], *y[1:])})
    e0_probabilities = [sp.simplify(p.subs(e0_substitution))
                        for p in gradient_probability_polynomials]
    require(
        e0_probabilities == [sp.Rational(1, 12), sp.Rational(1, 12),
                             sp.Rational(1, 12), sp.Rational(1, 3),
                             sp.Rational(1, 12), sp.Rational(1, 3)],
        "coherent witness changed",
    )
    # Also test a phase-sensitive source, rather than only a real basis input.
    complex_vector = [sp.Rational(1, 2), sp.I / 2, -sp.Rational(1, 2), sp.I / 2]
    complex_substitution = dict(zip(x, complex_vector))
    complex_substitution.update(dict(zip(y, map(sp.conjugate, complex_vector))))
    complex_probabilities = [sp.simplify(p.subs(complex_substitution))
                             for p in gradient_probability_polynomials]
    require(all(p.is_nonnegative for p in complex_probabilities),
            "complex source produced a negative Born probability")
    require(sum(complex_probabilities) == 1, "complex witness not normalized")

    report = {
        "claim": "MATH.NATIVE.COHERENT-QUARTIC-DERIVATIVE.01",
        "additional_claims": ["MATH.NATIVE.HAMMING-OCTIC-SOURCE.01"],
        "verdict": "PASS_IN_DECLARED_NATIVE_TENSOR_CATEGORY",
        "original_support_sizes": sizes,
        "original_W_shape": list(original_W.shape),
        "original_erasure_column_comparisons": original_column_checks,
        "coefficientwise_derivative_identities": derivative_checks,
        "marked_native_ovoid_rays": [list(q) for q in rays],
        "original_product_nine_quadric_residuals": source_ray_equations,
        "six_quartics": [str(f) for f in six_quartics],
        "six_sum_zero": True,
        "marked_igusa_polynomial_identity": True,
        "original_hamming_codeword_pairs": len(hamming_words)**2,
        "same_original_construction_A_hamming_tetrads": sorted(map(list, hamming_tetrads)),
        "hamming_genus_two_complete_enumerator": str(complete_enumerator),
        "same_source_octic_identities": "I8=96 sum_i f_i^2=(1/3) sum_A R_A^2",
        "scalar_source_fock_norm_squared": str(fock_norm_squared),
        "normalized_scalar_source": "F8hat=I8/(96 sqrt(140))",
        "born_equals_gradient_for_all_complex_inputs": True,
        "probability_sum_equals_norm_power_six": True,
        "probability_polynomials_x_and_conjugate_x": [
            polynomial_json(p, (*x, *y)) for p in gradient_probability_polynomials
        ],
        "e0_triple_probabilities": list(map(str, e0_probabilities)),
        "complex_source": list(map(str, complex_vector)),
        "complex_source_probabilities": list(map(str, complex_probabilities)),
        "source_sha256": {
            str(path.relative_to(REPO)): sha256(path) for path in (
                REPO / "tfpt_explorer/process.py", bridge_path,
                PRIOR / "native_erasure_readout.py", oldest_erasure_path,
                HERE.parent / "joint-native-algebra-selection-20261016"
                / "INTRINSIC-EVENT-GEOMETRY-AND-COMPOSITION.txt",
            )
        },
        "not_inferred": [
            "Preparation of identical registers from an unknown input",
            "Physical selection of the coherent input or of its tensor degree",
            "Physical source, time, state update, or spacetime dimensions",
            "Independence of this native construction from the E8/J input",
            "A new discovery of the already recorded Igusa quartic",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    print(json.dumps({"output": str(output), "claim": report["claim"],
                      "e0_probabilities": report["e0_triple_probabilities"],
                      "checks": "exact polynomial identities passed"}))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "coherent_quartic_response_results.json")
    args = parser.parse_args()
    main(args.output)
