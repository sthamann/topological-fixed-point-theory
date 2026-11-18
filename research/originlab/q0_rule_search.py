"""Search positive two-window forms on a finite affine word space."""

from itertools import combinations_with_replacement

import sympy as sp

from affine_operators import LatticeFock, require
from source import MARKS, compression_defect, evaluate_original


def affine_words(depth):
    """Generate current words in the level-one root subalgebra of the E8 chart."""
    fock = LatticeFock([[2]])
    cartan = (1,)
    currents = {
        "H": lambda n, v: fock.heisenberg(cartan, n, v),
        "E": lambda n, v: fock.vertex(cartan, n, v),
        "F": lambda n, v: fock.combine((-1, fock.vertex((-1,), n, v))),
    }
    by_grade = [[((), fock.vacuum())]]
    for grade in range(1, depth + 1):
        words = []
        for mode in range(1, grade + 1):
            for name, current in currents.items():
                for word, state in by_grade[grade - mode]:
                    image = current(-mode, state)
                    if image:
                        words.append((((name, -mode), *word), image))
        by_grade.append(sorted(words, key=lambda entry: (len(entry[0]), entry[0])))
    words = [(grade, word, vector) for grade, entries in enumerate(by_grade)
             for word, vector in entries]
    keys = sorted({key for _, _, vector in words for key in vector})
    columns = sp.Matrix([[vector.get(key, 0) for _, _, vector in words] for key in keys])
    _, pivots = columns.rref()
    basis = [words[i][2] for i in pivots]
    gram = sp.Matrix([[fock.inner(u, v) for v in basis] for u in basis])
    inv = gram.inv()

    def operator(name, mode):
        forms = sp.Matrix([[fock.inner(u, currents[name](mode, v)) for v in basis]
                           for u in basis])
        return inv * forms

    zero = {name: operator(name, 0) for name in currents}
    require(zero["E"].T * gram == gram * zero["F"], "zero-mode adjoint")
    require(zero["E"] * zero["F"] - zero["F"] * zero["E"] == zero["H"], "root Ward bracket")
    require(zero["H"] * zero["E"] - zero["E"] * zero["H"] == 2 * zero["E"], "Cartan Ward bracket")
    require(gram.is_positive_definite, "affine reference Gram")
    require(columns.rank() == sum(len(fock.basis((q,), grade - q*q))
                                  for grade in range(depth + 1) for q in range(-2, 3)
                                  if q*q <= grade), "truncated lattice sector not spanned")
    coordinates = columns[:, list(pivots)].gauss_jordan_solve(columns)[0]
    return words, pivots, coordinates, gram, zero, operator


def invariants(forms, test, dimension):
    if not forms:
        return []
    images = [test(form) for form in forms]
    rows = sp.Matrix.hstack(*(image.reshape(len(image), 1) for image in images))
    return [sum((coefficient * form for coefficient, form in zip(vector, forms)),
                sp.zeros(dimension)) for vector in rows.nullspace()]


def independent_probes(words, coordinates, directions):
    """Find independent one-window observations by increasing word cost."""
    probes = sorted(combinations_with_replacement(range(len(words)), 2),
                    key=lambda pair: (sum(len(words[i][1]) for i in pair),
                                      max(words[i][0] for i in pair), pair))
    values = [coordinates.T * direction * coordinates for direction in directions]
    tested = sp.zeros(0, len(directions))
    chosen = []
    for a, b in probes:
        row = sp.Matrix([[value[a, b] for value in values]])
        expanded = tested.col_join(row)
        if expanded.rank() > tested.rows:
            chosen.append({"left": words[a][1], "right": words[b][1],
                           "linear_response": [str(value) for value in row]})
            tested = expanded
            if tested.rows == len(directions):
                return chosen
    raise ArithmeticError("affine words must separate every finite form direction")


def search(depth=3):
    if type(depth) is not int or not 1 <= depth <= 4:
        raise ValueError("affine depth must be between one and four")
    words, pivots, coordinates, gram, zero, operator = affine_words(depth)
    size = gram.rows
    identity = sp.eye(size)
    inverse = gram.inv()
    grades = [words[i][0] for i in pivots]
    charges = zero["H"].diagonal()
    require(zero["H"].is_diagonal(), "native affine word basis must be Cartan homogeneous")
    metric_pairs = list(combinations_with_replacement(range(size), 2))
    graded_pairs = [(i,j) for i,j in metric_pairs if grades[i] == grades[j]]
    neutral_pairs = [(i,j) for i,j in graded_pairs if charges[i] == charges[j]]
    forms = [(identity[:, i] * identity[:, j].T if i == j else
              identity[:, i] * identity[:, j].T + identity[:, j] * identity[:, i].T)
             for i, j in neutral_pairs]
    stages = [{"rule": rule, "linear_dimension": len(pairs)} for rule,pairs in (
        ("reflected_symmetric_form",metric_pairs),("affine_grading",graded_pairs),("cartan_ward",neutral_pairs))]
    # G-selfadjointness of T is automatic for T=G^-1 Q and Q symmetric.
    # The E constraint also gives the F constraint by adjointness.
    root = zero["E"]
    adjoint_in_form_basis = gram*root*inverse
    forms = invariants(forms, lambda q: q*root-adjoint_in_form_basis*q, size)
    stages.append({"rule":"root_ward","linear_dimension":len(forms)})
    require(all(q*zero["F"] == gram*zero["F"]*inverse*q for q in forms), "lowering Ward")

    contacts = (0, next(i for i, pivot in enumerate(pivots) if words[pivot][0] == 1))
    baseline = gram / 2
    fixed = {}
    for name, i in zip(("vacuum_contact", "first_current_contact"), contacts):
        fixed[i] = baseline[i, i]
        previous = len(forms)
        forms = invariants(forms, lambda q: sp.Matrix([[q[i, i]]]), size)
        stages.append({"rule": name, "linear_dimension": len(forms),
                       "independent": len(forms) < previous, "control_value": str(fixed[i])})

    result = {
        "status": "EXACT_FINITE_UNIQUENESS" if not forms else "EXACT_FINITE_UNDERDETERMINED",
        "input_role": "level_one_affine_comparison_only",
        "basis_dimension": size,
        "grade_dimensions": [sum(words[i][0] == grade for i in pivots) for grade in range(depth + 1)],
        "word_count": len(words),
        "stages": stages,
        "original_values": dict.fromkeys(MARKS),
        "absolute_original_reference": None,
        "physical_solution": False,
    }
    if not forms:
        return result

    # Choose an exact interior pair. Candidate synthesis and positivity are
    # separate: the nullspace alone is not evidence of positive solutions.
    direction = forms[0]
    epsilon = sp.S.One
    while True:
        left, right = baseline - epsilon*direction, baseline + epsilon*direction
        if all(q.is_positive_definite and (gram - q).is_positive_definite for q in (left, right)):
            break
        epsilon /= 2
    candidates = []
    for q in (left, right):
        require(all(q[i, i] == target for i, target in fixed.items()), "contact preserved")
        transfer = inverse*q
        two = q*transfer
        require(compression_defect(gram, q, two) == sp.zeros(size), "two-window sewing")
        # Distinct operator insertions, with physical adjoints represented
        # on this finite space. These are affine probes, not Higgs marks.
        u = operator("H", -1) + operator("H", 1)
        y = operator("E", -1) + operator("F", 1)
        require(u.T*gram == gram*u and y.T*gram == gram*y, "probe adjoints")
        candidates.append({
            "empty": two,
            "cartan": q*u*transfer,
            "root": q*y*transfer,
            "ordered": q*u*y*transfer,
        })
    winner = None
    for name, first in candidates[0].items():
        second = candidates[1][name]
        difference = coordinates.T*(second - first)*coordinates
        for a, b in difference.todok():
            cost = len(words[a][1]) + len(words[b][1]) + {"empty": 0, "cartan": 1, "root": 1, "ordered": 2}[name]
            key = (cost, max(words[a][0], words[b][0]), name, a, b)
            if winner is None or key < winner:
                winner = key
    require(winner is not None, "candidate pair is separated by a declared word")
    cost, _, name, a, b = winner
    v1, v2 = coordinates[:, a], coordinates[:, b]
    va, vb = v1.T*candidates[0][name]*v2, v1.T*candidates[1][name]*v2
    unit = sp.sqrt((v1.T*gram*v1)[0] * (v2.T*gram*v2)[0])
    result.update(
        independent_one_window_probes=independent_probes(words, coordinates, forms),
        independent_words=[[list(letter) for letter in words[i][1]] for i in pivots],
        reference_gram=[[str(value) for value in row] for row in gram.tolist()],
        one_window_pair=[[[str(value) for value in row] for row in form.tolist()]
                         for form in (left, right)],
        separating_word={"left": [list(letter) for letter in words[a][1]],
                         "right": [list(letter) for letter in words[b][1]], "insertion": name,
                         "cost": cost, "left_candidate": str(va[0]), "right_candidate": str(vb[0]),
                         "unit_word_candidates": [str(va[0]/unit), str(vb[0]/unit)]},
        two_window_candidate_pair=[
            {name: [[str(x) for x in row] for row in q.tolist()] for name, q in forms.items()}
            for forms in candidates
        ],
        first_original_test="Read the separating two-window form only after a native word and cap binding.",
    )
    return result


def native_comparison(result, original):
    separator = result.get("separating_word")
    probe = original.get("q0_rule_probe")
    if separator is None:
        return {"status": "NO_SEPARATOR_AT_THIS_DEPTH"}
    if probe is None:
        return {"status": "SOURCE_PROBE_REQUIRED", "original_status": original["status"]}
    if (separator["insertion"] != "empty" or
            any(probe[side] != separator[side] for side in ("left", "right"))):
        return {"status": "DIFFERENT_WORD_PROBE", "probe": probe}
    # The source reader emitted this exact expression. By retaining the
    # source role, an arithmetic fixture cannot become physical evidence.
    expected = separator["unit_word_candidates"]
    observed = probe["unit_word_response"]
    survivors = [i for i, value in enumerate(expected) if sp.simplify(
        sp.sympify(value) - sp.sympify(observed)) == 0]
    return {"status": "SCOPED_PROBE_MATCH" if survivors else "PAIR_EXCLUDED_BY_PROBE",
            "surviving_pair_indices": survivors, "probe": probe,
            "physical_source_identity_proved": False}


def run(root, domain=None, depth=3):
    result = search(depth)
    if root is not None:
        result["native_comparison"] = native_comparison(result, evaluate_original(root, domain))
    return result
