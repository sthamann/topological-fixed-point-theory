"""Exact necessary type test on the existing current Weyl relation.

This does not construct a physical P1 source. It checks the conditional
OS-word/current-Weyl formula in the preceding reflection report before
attempting to reconstruct states or optimize any source coefficients.
Only Python's exact integer/Fraction arithmetic is used.
"""

from fractions import Fraction
from itertools import product
import json
from pathlib import Path


def transpose(a):
    return [list(x) for x in zip(*a)]


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)]
            for row in a]


def scale(c, a):
    return [[c * x for x in row] for row in a]


def kron(a, b):
    return [[x * y for x in row_a for y in row_b]
            for row_a in a for row_b in b]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def evaluate_weyl_word(labels, symplectic):
    """Direct multiplication, phase in pi/4 units and sigma in pi/2 units."""
    total = [0] * len(symplectic)
    phase = 0
    for label in labels:
        pairing = sum(total[i] * symplectic[i][j] * label[j]
                      for i in range(len(total)) for j in range(len(total)))
        phase -= pairing
        total = [x + y for x, y in zip(total, label)]
    return phase % 8, total


def apply_matrix(matrix, vector):
    return [sum(x * y for x, y in zip(row, vector)) for row in matrix]


def run():
    # One real Fourier pair. Reflection fixes cosine and reverses sine.
    # For phi_1=1, psi_1=i the ORIGINAL formula gives [Jphi,Jpsi]=2i.
    phi1, psi1 = 1 + 0j, 1j
    raw = phi1.conjugate() * psi1 - phi1 * psi1.conjugate()
    require(raw == 2j, "Original first-mode current convention changed")
    require((phi1 * psi1.conjugate() - phi1.conjugate() * psi1) == -raw,
            "Orientation reversal must flip the actual current cocycle")
    omega = [[0, 1], [-1, 0]]
    identity = [[1, 0], [0, 1]]
    reflection = [[1, 0], [0, -1]]
    swap = [[0, 1], [1, 0]]
    # The first two are actual level conventions already present in the
    # source. The last two are ONLY candidate target algebra types.
    channel_cases = [
        ("original_cover_level1", [[1]], [[1]], False),
        ("folded_total_level8", [[8]], [[1]], False),
        ("oppositely_oriented_face_pair_no_swap", [[8, 0], [0, -8]],
         identity, True),
        ("oppositely_oriented_face_pair_with_swap", [[8, 0], [0, -8]],
         swap, True),
    ]
    rows = []
    for (name, metric, face_map, conditional_target), c, o, r, field_sign in product(
            channel_cases, (1, -1), (1, -1), (1, -1), (1, -1)):
        symplectic = kron(metric, omega)
        mode_map = identity if r == 1 else reflection
        operator_map = scale(field_sign, kron(face_map, mode_map))
        transformed = multiply(multiply(transpose(operator_map), symplectic),
                               operator_map)
        # c=+1: linear; c=-1: conjugate-linear.
        # o=+1: multiplicative; o=-1: order reversing.
        # Weyl relation W(f)W(g)=exp(-i sigma(f,g)/2)W(f+g)
        # requires sigma(Lf,Lg)=c*o*sigma(f,g).
        passes = transformed == scale(c * o, symplectic)
        # Independently multiply actual words. This does not use the
        # c*o identity to predict the observed phases.
        word_passes = True
        for i, j, fscale, gscale in product(range(len(symplectic)),
                                           range(len(symplectic)),
                                           (Fraction(1), Fraction(1, 3)),
                                           (Fraction(1), Fraction(1, 5))):
            f = [fscale if x == i else 0 for x in range(len(symplectic))]
            g = [gscale if x == j else 0 for x in range(len(symplectic))]
            word = [f, g, [-x for x in f], [-x for x in g]]
            source_phase, source_total = evaluate_weyl_word(word, symplectic)
            require(not any(source_total), "Commutator word did not close")
            image_word = [apply_matrix(operator_map, x) for x in word]
            if o == -1:
                image_word.reverse()
            image_phase, image_total = evaluate_weyl_word(image_word, symplectic)
            require(not any(image_total), "Reflected commutator did not close")
            word_passes &= image_phase == (c * source_phase) % 8
        require(word_passes == passes, "Direct Weyl multiplication disagrees with type theorem")
        actual_sign = next((s for s in (1, -1)
                            if transformed == scale(s, symplectic)), None)
        rows.append(dict(
            channel_case=name, conditional_target=conditional_target,
            scalar_type="linear" if c == 1 else "conjugate_linear",
            product_type="multiplicative" if o == 1 else "order_reversing",
            spatial_map="identity" if r == 1 else "reflection",
            field_sign=field_sign, required_cocycle_sign=c * o,
            actual_cocycle_sign=actual_sign, weyl_relation_preserved=passes,
            independently_multiplied_weyl_words_match=word_passes))

    single_os = [x for x in rows if not x["conditional_target"]
                 and x["scalar_type"] == "conjugate_linear"
                 and x["product_type"] == "order_reversing"
                 and x["spatial_map"] == "reflection"]
    paired_os = [x for x in rows
                 if x["channel_case"] == "oppositely_oriented_face_pair_with_swap"
                 and x["scalar_type"] == "conjugate_linear"
                 and x["product_type"] == "order_reversing"
                 and x["spatial_map"] == "reflection"]
    require(all(not x["weyl_relation_preserved"] for x in single_os),
            "Claimed single-chiral exclusion failed")
    require(all(x["weyl_relation_preserved"] for x in paired_os),
            "Conditional opposite-face algebraic repair failed")
    require(len(rows) == 64 and sum(x["weyl_relation_preserved"] for x in rows) == 32,
            "Incomplete finite type table")

    # A Weyl character / additive scalar shift cannot repair the cocycle:
    # its phases cancel identically in the four-unitary group commutator.
    character_exponent = Fraction(2, 7) + Fraction(3, 11) - Fraction(2, 7) - Fraction(3, 11)
    require(character_exponent == 0, "Central character did not cancel")

    return {
        "verdict": "SINGLE_CHIRAL_OS_ANTI_ANTI_FORMULA_EXCLUDED_FOR_ORDINARY_WEYL_PRODUCT",
        "scope": "Necessary algebraic type test; no original P1 source, state, OS positivity, time or spacetime reconstruction proved.",
        "enumerated_variants": len(rows),
        "passed_algebraic_type_test": sum(x["weyl_relation_preserved"] for x in rows),
        "formula": "sigma(Lf,Lg) = scalar_conjugation_sign * product_order_sign * sigma(f,g)",
        "single_chiral_os_reflection_variants": single_os,
        "conditional_opposite_face_os_reflection_variants": paired_os,
        "central_phase_witness": {
            "choose_original_group_commutator": "exp(-i*pi/2) = -i",
            "conjugate_linear_image_required": "+i",
            "order_reversal_plus_single_spatial_reflection_gives": "-i",
            "contradiction": True,
        },
        "constant_field_sign_and_weyl_character_do_not_repair": True,
        "all_variants": rows,
        "complete_TFPT_solution": False,
        "ledger_changed": False,
    }


if __name__ == "__main__":
    result = run()
    destination = Path(__file__).with_name("weyl-reflection-type-results.json")
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("all_variants", "single_chiral_os_reflection_variants",
                                   "conditional_opposite_face_os_reflection_variants")},
                     indent=2, sort_keys=True))
