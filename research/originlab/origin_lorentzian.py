"""Exact conditional II_{9,1}=E8+U parent with marked null reduction.

The charge lattice and explicit norm-preserving transitions are constructed.
The null transvections act identically on the E8 quotient and cannot be
unitary in this finite coordinate representation for any positive metric.
An additionally chosen root reflection does act on the quotient, but does
not select a physical process, time, Hamiltonian, state, positive seam kernel,
or c3. No original TFPT source is replaced.

Lattice classification context: Plamadeala, Mulligan and Nayak,
Phys. Rev. B 88, 045131 (2013), section II.C, equations (16), (18), (21),
https://arxiv.org/abs/1304.0772 . The finite checks below do not prove that
classification theorem or physically select its hypotheses.
"""

from fractions import Fraction


def require(condition, why):
    if not condition:
        raise ArithmeticError(why)


def transpose(matrix):
    return list(map(list, zip(*matrix)))


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, column))
             for column in zip(*b)] for row in a]


def determinant(matrix):
    require(all(len(row) == len(matrix) for row in matrix),
            "determinant requires a square matrix")
    a = [[Fraction(value) for value in row] for row in matrix]
    result = Fraction(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result = -result
        value = a[i][i]
        result *= value
        for j in range(i + 1, len(a)):
            factor = a[j][i] / value
            a[j][i] = Fraction(0)
            for k in range(i + 1, len(a)):
                a[j][k] -= factor * a[i][k]
    return result


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def e8_cartan():
    # Diagram arms have lengths 1,2,4; this is an E8 simple-root basis.
    matrix = [[2 * int(i == j) for j in range(8)] for i in range(8)]
    for i, j in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (1, 6), (6, 7)]:
        matrix[i][j] = matrix[j][i] = -1
    return matrix


def transvection(cartan, vector):
    n = len(cartan)
    require(len(vector) == n and all(isinstance(x, int) for x in vector),
            "transvection requires an integral vector of the lattice rank")
    cv = [sum(cartan[i][j] * vector[j] for j in range(n)) for i in range(n)]
    norm = sum(vector[i] * cv[i] for i in range(n))
    require(norm % 2 == 0, "null transvection requires even vector norm")
    matrix = identity(n + 2)
    for i in range(n):
        matrix[i][n + 1] = vector[i]
        matrix[n][i] = -cv[i]
    matrix[n][n + 1] = -norm // 2
    return matrix


def root_reflection(gram, root):
    """Integral norm-two reflection; choosing this transition is extra input."""
    n = len(gram)
    require(all(len(row) == n for row in gram)
            and len(root) == n and all(isinstance(x, int) for x in root),
            "root reflection requires a square Gram and integral root")
    dual = [sum(gram[i][j] * root[j] for j in range(n))
            for i in range(n)]
    require(sum(root[i] * dual[i] for i in range(n)) == 2,
            "root reflection requires a norm-two root")
    return [[int(i == j) - root[i] * dual[j] for j in range(n)]
            for i in range(n)]


def controls():
    cartan = e8_cartan()
    minors = [int(determinant([row[:n] for row in cartan[:n]]))
              for n in range(1, 9)]
    require(minors == [2, 3, 4, 5, 6, 7, 4, 1],
            "E8 positive-definite unimodular Gram certificate changed")
    require(cartan == transpose(cartan)
            and all(cartan[i][i] == 2 for i in range(8)),
            "E8 Gram is not symmetric even integral")
    gram = [row + [0, 0] for row in cartan]
    gram += [[0] * 8 + [0, 1], [0] * 8 + [1, 0]]
    require(determinant(gram) == -1,
            "E8 plus hyperbolic plane lost unimodularity")
    # Sylvester's criterion proves the E8 block positive; U has eigenvalues +/-1.
    require(all(value > 0 for value in minors)
            and gram[8][8] == gram[9][9] == 0
            and gram[8][9] == gram[9][8] == 1,
            "parent signature or marked primitive null pair changed")
    vectors = [[int(i == j) for i in range(8)] for j in range(8)]
    vectors += [[1, -2, 3, -1, 0, 2, -3, 1], [2, 0, -1, 4, -2, 0, 1, -1]]
    null_e = [[0] for _ in range(8)] + [[1], [0]]
    for vector in vectors:
        matrix = transvection(cartan, vector)
        require(multiply(transpose(matrix), multiply(gram, matrix)) == gram,
                "null transvection failed to preserve the exact Gram form")
        require(determinant(matrix) == 1,
                "null transvection is not unimodular")
        require(multiply(matrix, transvection(cartan, [-x for x in vector]))
                == identity(10), "null transvection inverse failed")
        require(multiply(matrix, null_e) == null_e,
                "null transvection failed to fix the marked null direction")
        # For b=0, T_v(x,a,0)=(x,a-(x,v),0).  The changed a coordinate
        # is discarded by e^perp/Z e, so the induced E8 action is identity.
        require([row[:8] for row in matrix[:8]] == identity(8)
                and matrix[9][:9] == [0] * 9,
                "null transvection is not identity on the E8 quotient")
    for vector in vectors:
        for other in vectors:
            require(multiply(transvection(cartan, vector),
                             transvection(cartan, other)) ==
                    transvection(cartan, [a + b for a, b in zip(vector, other)]),
                    "additive null transvection composition failed")
    # e=(0,...,0,1,0). (e,z)=b, so e^perp={(x,a,0)}.
    # Quotienting Ze removes a and leaves precisely the E8 Gram matrix.
    quotient_gram = [row[:8] for row in gram[:8]]
    require(quotient_gram == cartan,
            "marked primitive null quotient is not the E8 Gram matrix")

    # A nontrivial unipotent operator cannot be unitary for any positive
    # definite inner product: positive-metric unitary operators are
    # diagonalizable, whereas the following minimal polynomial is (X-1)^3.
    root_v = vectors[0]
    root_transvection = transvection(cartan, root_v)
    nilpotent = [[root_transvection[i][j] - int(i == j)
                  for j in range(10)] for i in range(10)]
    nilpotent_squared = multiply(nilpotent, nilpotent)
    nilpotent_cubed = multiply(nilpotent_squared, nilpotent)
    require(any(any(row) for row in nilpotent),
            "unipotent negative witness became identity")
    require(any(any(row) for row in nilpotent_squared),
            "unipotent negative witness lost its degree-three Jordan block")
    require(not any(any(row) for row in nilpotent_cubed),
            "unipotent negative witness is not nilpotent of degree three")

    # Contrast with a transition that is visible after the same quotient.
    # alpha=(v,1,0) has norm two and (e,alpha)=0.  Its reflection fixes e
    # and induces x -> x-(x,v)v below.  Selecting it remains extra input.
    alpha = root_v + [1, 0]
    reflection = root_reflection(gram, alpha)
    require(multiply(transpose(reflection), multiply(gram, reflection)) == gram,
            "additional root reflection failed Gram preservation")
    require(multiply(reflection, reflection) == identity(10)
            and determinant(reflection) == -1,
            "additional root reflection is not a unimodular involution")
    require(multiply(reflection, null_e) == null_e,
            "additional root reflection failed to preserve the marked null vector")
    cv = [sum(cartan[i][j] * root_v[j] for j in range(8))
          for i in range(8)]
    expected_reflection = [[int(i == j) - root_v[i] * cv[j]
                            for j in range(8)] for i in range(8)]
    reflection_quotient = [row[:8] for row in reflection[:8]]
    require(reflection_quotient == expected_reflection
            and reflection_quotient != identity(8),
            "additional root reflection has wrong or trivial E8 quotient action")
    require(multiply(reflection_quotient, [[x] for x in root_v])
            == [[-x] for x in root_v],
            "additional quotient reflection does not reverse its root")
    return {
        "status": "EXACT_CONDITIONAL_LORENTZIAN_PARENT_AND_NULL_QUOTIENT",
        "parent_lattice": "II_{9,1}=E8 direct_sum U",
        "parent_gram": gram,
        "parent_determinant": int(determinant(gram)),
        "parent_signature": [9, 1],
        "e8_leading_principal_minors": minors,
        "primitive_null_vectors": [[0] * 8 + [1, 0], [0] * 8 + [0, 1]],
        "primitive_null_pairing": 1,
        "null_quotient_gram": quotient_gram,
        "null_quotient": "e^perp/Z e is E8",
        "transvections_checked": len(vectors),
        "transvection_compositions_checked": len(vectors) ** 2,
        "transvection_quotient_actions_checked": len(vectors),
        "transvections_induce_identity_on_E8_null_quotient": True,
        "nontrivial_E8_quotient_dynamics_from_transvections_derived": False,
        "transvection_formula": "T_v(x,a,b)=(x+bv,a-(x,v)-b(v,v)/2,b)",
        "general_norm_identity": "Q(T_v z)=Q(z), Q(x,a,b)=(x,x)+2ab",
        "general_composition_identity": "T_v T_w=T_(v+w)",
        "transvection_unipotent_witness": {
            "root_vector": root_v,
            "N_equals_T_minus_identity_nonzero": True,
            "N_squared_nonzero": True,
            "N_cubed_zero": True,
            "minimal_polynomial": "(X-1)^3",
            "positive_definite_unitary_coordinate_realization_possible": False,
            "reason": "Nontrivial unipotence obstructs diagonalizability; positive-metric unitary matrices are diagonalizable.",
            "alternative_unitary_lift": "Permutation action on l2(L) is possible but adds a Hilbert representation, not a selected physical state or clock.",
        },
        "additional_quotient_visible_transition": {
            "selection_is_extra_input": True,
            "parent_root_alpha": alpha,
            "formula": "r_alpha(z)=z-(z,alpha)*alpha",
            "order": 2,
            "determinant": -1,
            "preserves_marked_null_vector": True,
            "quotient_action": "x -> x-(x,v)*v",
            "quotient_matrix": reflection_quotient,
            "quotient_action_is_nontrivial": True,
            "E8_or_native_four_clock_selected": False,
            "physical_process_or_time_selected": False,
        },
        "external_classification_reference": "https://arxiv.org/abs/1304.0772",
        "classification_theorem_proved_by_finite_matrix_checks": False,
        "assumptions": [
            "even unimodular charge lattice",
            "one negative direction",
            "minimal nonzero positive transverse quotient",
            "marked primitive null direction",
        ],
        "physical_selection_of_these_assumptions_derived": False,
        "physical_time_from_negative_direction_derived": False,
        "distinguished_transfer_or_Hamiltonian_selected": False,
        "positive_Hilbert_reconstruction_derived": False,
        "original_one_sided_spherical_P1_source_identified": False,
        "c3_one_over_8pi_derived_from_parent": False,
        "physical_3plus1D_spacetime_derived": False,
        "F1_solved": False,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(controls(), sort_keys=True, indent=2))
