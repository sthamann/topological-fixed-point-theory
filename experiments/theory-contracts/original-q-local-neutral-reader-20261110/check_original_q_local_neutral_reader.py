"""Compare the neutral shifted-stress reader with the actual old Fold vertices."""
from contextlib import redirect_stdout
from functools import cache
import io
import json
from math import factorial, prod
from pathlib import Path
import runpy
import sys

repo = Path(sys.argv[1]).resolve()
with redirect_stdout(io.StringIO()):
    source = runpy.run_path(str(repo / "experiments/theory-contracts/chat-triangulation-20261108/"
                               "check_tfpt_family_middle_ordered_cut_20261110.py"))

sp, field, q, r, b = (source[key] for key in ("sp", "Q", "q", "r", "d1"))
old, vertex, inner = (source[key] for key in ("old", "vertex", "inner"))
assert source["axes"][:, 4] == q and q.dot(q) == 2
q_component = [field.from_sympy(sp.simplify(q.dot(axis)/2))
               for axis in (source["axes"][:, i] for i in range(6))]
r_component = [field.from_sympy(sp.simplify(r.dot(axis)/2))
               for axis in (source["axes"][:, i] for i in range(6))]


@cache
def partitions(grade, minimum=1):
    if grade == 0:
        return ((),)
    return tuple((n,) + rest for n in range(minimum, grade + 1)
                 for rest in partitions(grade - n, n))


@cache
def gram(grade):
    return sp.diag(*(prod((2*n)**word.count(n)*factorial(word.count(n))
                          for n in set(word)) for word in partitions(grade)))


def affine_stress(grade, step, charge):
    """Shifted L_-step at actual Cartan charge a_0=charge."""
    src, dst = partitions(grade), partitions(grade + step)
    indices = {word: i for i, word in enumerate(dst)}
    out = sp.zeros(len(dst), len(src))
    for j, word in enumerate(src):
        out[indices[tuple(sorted(word + (step,)))], j] += sp.Rational(charge-1, 2)
        for k in range(1, step):
            out[indices[tuple(sorted(word + (k, step-k)))], j] += sp.Rational(1, 4)
        for i, n in enumerate(word):
            out[indices[tuple(sorted(word[:i] + word[i+1:] + (n+step,)))], j] += n
    return out


def vertex_square(grade, charge, root_sign):
    momentum = r if charge == 0 else q
    target_grade = grade + root_sign - 1 - root_sign*charge
    if target_grade < 0:
        return sp.zeros(len(partitions(grade)))
    target = [tuple((n, 4) for n in word) for word in partitions(target_grade)]
    actual = [vertex({tuple((n, 4) for n in word): field.one},
                     root_sign*q, momentum, root_sign) for word in partitions(grade)]
    assert all(word in target for column in actual for word in column)
    matrix = sp.Matrix(len(target), len(actual),
                       lambda i, j: field.to_sympy(actual[j].get(target[i], field.zero)))
    return gram(grade).inv() * matrix.T * gram(target_grade) * matrix


def affine_blocks(charge, max_grade):
    sectors, squares, signed_roots, sizes = {}, {}, {}, {}
    for grade in range(max_grade + 1):
        dim = len(partitions(grade))
        metric = gram(grade)
        descendants = {}
        for step in range(1, grade + 1):
            shift = affine_stress(grade-step, step, charge)
            for n, states in sectors[grade-step].items():
                descendants.setdefault(n, []).extend(shift * state for state in states)
        levels = {n: sp.Matrix.hstack(*states).columnspace()
                  for n, states in descendants.items() if states}
        rows = [state for states in levels.values() for state in states]
        highest = (sp.Matrix.hstack(*rows).T * metric).nullspace() if rows else list(sp.eye(dim).columnspace())
        new_labels = [n for n in range(grade + 1) if n*(n+1) == grade]
        assert len(highest) == len(new_labels)
        if highest:
            levels[new_labels[0]] = highest
        raise_square = vertex_square(grade, charge, 1)
        if charge == 0:
            assert raise_square - vertex_square(grade, charge, -1) == sp.eye(dim)
        neutral_reader = sp.zeros(dim)
        signed_root = sp.zeros(dim)
        for n, states in levels.items():
            block = sp.Matrix.hstack(*states)
            assert len(states) == len(partitions(grade-n*(n+1))) - len(partitions(grade-(n+1)*(n+2)))
            eig = (n+1)**2 - charge//2
            assert raise_square * block == eig * block
            projector = block * (block.T * metric * block).inv() * block.T * metric
            neutral_reader += eig * projector
            signed_root += (-1)**n * (n+1) * projector
        assert neutral_reader == raise_square and sum(map(len, levels.values())) == dim
        squares[grade], signed_roots[grade], sectors[grade] = neutral_reader, signed_root, levels
        sizes[grade] = {str(n): len(states) for n, states in levels.items()}
    return squares, signed_roots, sizes


def change_basis(state, sign, split_r=False):
    """Separate the named orthogonal Cartan axes, or invert the change."""
    terms = []
    for word, coefficient in state.items():
        expanded = [((), coefficient)]
        for n, axis in word:
            extended = [(w + ((n, axis),), a) for w, a in expanded]
            if axis != 4 and q_component[axis]:
                extended.extend((w + ((n, 4),), sign*q_component[axis]*a) for w, a in expanded)
            if split_r and axis != 5 and r_component[axis]:
                extended.extend((w + ((n, 5),), sign*r_component[axis]*a) for w, a in expanded)
            expanded = extended
        terms.extend((tuple(sorted(w)), a) for w, a in expanded)
    return source["combine"](terms)


def apply_reader(reader, state, parity=False):
    state = change_basis(state, 1)
    grouped = {}
    for word, coefficient in state.items():
        qword = tuple(n for n, axis in word if axis == 4)
        passive = tuple((n, axis) for n, axis in word if axis != 4)
        grouped.setdefault(passive, {})[qword] = coefficient
    out = {}
    for passive, block in grouped.items():
        grade = sum(next(iter(block)))
        words = partitions(grade)
        vector = sp.Matrix([field.to_sympy(block.get(word, field.zero)) for word in words])
        image = reader[grade] * vector
        for word, coefficient in zip(words, image, strict=True):
            if coefficient:
                if parity:
                    coefficient *= (-1)**len(word)
                out[tuple(sorted(passive + tuple((n, 4) for n in word)))] = field.from_sympy(coefficient)
    return change_basis(out, -1)


def residual_inertia(late, early, split_r):
    """Signs of the exact contrast after tracing over the indicated Cartan oscillators."""
    active = {4, 5} if split_r else {4}
    split = []
    for vector in (late, early):
        groups = {}
        for word, coefficient in change_basis(vector, 1, split_r).items():
            left = tuple(token for token in word if token[1] in active)
            right = tuple(token for token in word if token[1] not in active)
            groups.setdefault(left, {})[right] = field.to_sympy(coefficient)
        split.append(groups)
    passive = sorted({word for groups in split for terms in groups.values() for word in terms})
    columns, weights = [], []
    for sign, groups in zip((1, -1), split, strict=True):
        for active_word, terms in groups.items():
            weight = prod((2*n)**active_word.count((n, axis))*factorial(active_word.count((n, axis)))
                          for n, axis in set(active_word))
            columns.append(sp.Matrix([terms.get(word, 0) for word in passive]))
            weights.append(sign * weight)
    vectors = sp.Matrix.hstack(*columns)
    _, basis_columns = vectors.rref()
    basis = vectors[:, list(basis_columns)]
    _, basis_rows = basis.T.rref()
    coordinates = (basis[list(basis_rows), :].inv() * vectors[list(basis_rows), :]).applyfunc(sp.expand)
    assert (basis * coordinates - vectors).applyfunc(sp.expand) == sp.zeros(*vectors.shape)
    # The passive monomials form a basis. Congruence to the independent
    # columns preserves inertia and removes the zero-only directions.
    contrast = (coordinates * sp.diag(*weights) * coordinates.H).applyfunc(sp.expand)
    coeff = [sp.simplify(c) for c in contrast.charpoly().all_coeffs()]
    assert all(c.is_real for c in coeff)
    zeros = len(passive) - max(i for i, c in enumerate(coeff) if c)
    degree = len(coeff)-1
    def sign_variations(coefficients):
        signs = [sp.sign(c) for c in coefficients if c]
        return sum(a != b for a, b in zip(signs, signs[1:]))
    positive = sign_variations(coeff)
    negative = sign_variations([(-1)**(degree-i)*c for i, c in enumerate(coeff)])
    assert positive + negative + zeros == len(passive)
    return [positive, negative, zeros]


neutral, neutral_signed_root, neutral_sizes = affine_blocks(0, 7)
charged, _, charged_sizes = affine_blocks(2, 7)
assert source["axes"].rank() == 6
output = {}
for alpha in (sp.Rational(1, 4), sp.Rational(3, 4)):
    drift = sp.Matrix([j + alpha - 4 for j in range(8)])
    connection = [field.from_sympy(z) for z in old["cl"](source["V"].T * drift / 8)]
    packets = [old["vertex_zero"](old["series_response"](f, connection, 4))
               for f in source["polys"]]
    rows = []
    for packet in packets:
        b_packet = vertex(packet, b, q)
        pre = vertex(packet, q, q, 1)
        early = vertex(pre, b, 2*q) if pre else {}
        post = vertex(b_packet, q, r, 1)
        reverse = vertex(b_packet, -q, r, -1)
        assert inner(packet, apply_reader(charged, packet)) == inner(pre, pre)
        late = inner(b_packet, apply_reader(neutral, b_packet))
        assert late == inner(post, post)
        assert apply_reader(neutral_signed_root, b_packet, parity=True) == post
        assert late - inner(reverse, reverse) == inner(b_packet, b_packet)
        early_norm = inner(early, early)
        overlap = inner(post, early)
        det = sp.simplify(late*early_norm - overlap*sp.conjugate(overlap))
        contrast = sp.simplify(late-early_norm)
        trace_norm = sp.sqrt(contrast**2+4*det)
        rows.append({"B_norm2": str(inner(b_packet, b_packet)),
                     "early_q_leg_norm2_before_B": str(inner(pre, pre)),
                     "early_BJ_norm2": str(early_norm),
                     "late_JB_norm2": str(late),
                     "overlap_late_early": str(overlap),
                     "output_gram_determinant": str(det),
                     "unweighted_late_minus_early": str(contrast),
                     "sharp_relative_form_bound_for_late_dominance": float(contrast/trace_norm)
                     if trace_norm else None,
                     "trace_q_contrast_inertia_pos_neg_zero": residual_inertia(post, early, False),
                     "trace_qr_contrast_inertia_pos_neg_zero": residual_inertia(post, early, True),
                     "reverse_q_leg_norm2_after_B": str(inner(reverse, reverse))})
    output[str(alpha)] = rows

print(json.dumps({"charge_zero_multiplicities": neutral_sizes,
                  "charge_two_multiplicities": charged_sizes,
                  "original_packets": output,
                  "physical_weights_computed": False}, indent=2))
