"""Binary context and quaternionic two-window controls for the P1/P2 bridge."""

from itertools import combinations, product

import sympy as sp

from reconstruction import PAULI, require
from selection import BIT_VECTORS, full_face_refinements, index, mod2


def controls():
    checks = []

    def check(name, condition):
        require(condition, name)
        checks.append(name)

    def same(left, right):
        return (left-right).applyfunc(sp.simplify) == sp.zeros(*left.shape)

    identity, x, y, z = PAULI
    tensor = sp.kronecker_product
    eye4 = sp.eye(4)
    t = sp.Rational(3, 5)
    a = sp.diag(z, z)
    b_plus, b_minus = t*z+sp.Rational(4, 5)*x, -t*z+sp.Rational(4, 5)*x
    b = sp.diag(b_plus, b_minus)
    central_angle = (a*b+b*a)/2
    check("both projections are sharp and have independent binary marginals", a*a == b*b == eye4
          and all(sp.trace(value) == 0 for value in (a, b, a*b)))
    check("the central angles in the two blocks do not vanish", central_angle == sp.diag(t, t, -t, -t)
          and central_angle*a == a*central_angle and central_angle*b == b*central_angle)
    check("one four-letter word detects cancellation of the two angles", sp.trace((a*b)**2)/4 == -sp.Rational(7, 25)
          and same(central_angle**2, ((a*b)**2+(b*a)**2+2*eye4)/4))
    marks = [(identity+sign*axis)/2 for axis, sign in product((x, z), (1, -1))]
    check("four rays are two antipodal pairs with orthogonal Bloch axes", len({tuple(mark) for mark in marks}) == 4
          and marks[0]+marks[1] == marks[2]+marks[3] == identity
          and all(sp.trace(u*v) == sp.Rational(1, 2) for u, v in product(marks[:2], marks[2:])))

    words = {(a, b, c, d): (sp.I**(a*b+c*d)*tensor(x**a*z**b, x**c*z**d)).applyfunc(sp.simplify)
             for a, b, c, d in BIT_VECTORS}
    linear_j = x*z
    thetas = (eye4, tensor(linear_j, identity), tensor(linear_j, linear_j))
    check("real and quaternionic choices have squares +1, -1, +1", all(
        lift*lift.conjugate() == sign*eye4 for lift, sign in zip(thetas, (1, -1, 1))
    ))
    sizes = [sum(lift*word.conjugate()*lift.inv() == word for word in words.values())-1 for lift in thetas]
    check("two real, mixed and two quaternionic windows give dimensions nine, five, nine", sizes == [9, 5, 9])
    theta = thetas[1]

    def quadratic(v):
        a, b, c, d = v
        return (a+b+a*b+c*d) % 2

    def polar(u, v):
        a, b, c, d = u
        e, f, g, h = v
        return (a*f+b*e+c*h+d*g) % 2

    labels = tuple(v for v in words if any(v) and not quadratic(v))
    check("antiunitary signs on all sixteen Hermitian Pauli words", all(
        theta*word.conjugate()*theta.inv() == (-1)**quadratic(v)*word for v, word in words.items()
    ))
    check("the sign refinement polarizes to the full two-qubit commutator", all(
        (quadratic(tuple(a^b for a, b in zip(u, v)))+quadratic(u)+quadratic(v)) % 2 == polar(u, v)
        and words[u]*words[v] == (-1)**polar(u, v)*words[v]*words[u]
        for u, v in product(BIT_VECTORS, repeat=2)
    ))
    # Ordering puts the two one-window axes first and the common-quarter fixed axis last.
    gamma = [tensor(identity, x), tensor(identity, z), tensor(x, y), tensor(z, y), tensor(y, y)]
    label_order = [next(v for v in labels if words[v] == operator) for operator in gamma]
    check("five axes carry Clifford products and have common volume +I or -I", len(labels) == 5
          and all(g*g == eye4 for g in gamma) and all(g*h == -h*g for g, h in combinations(gamma, 2))
          and sp.prod(gamma) in (eye4, -eye4))
    check("two axes are supported on the second window and three on both", [sum(v[i:i+2] == (0, 0)
          for v in labels) for i in (0, 2)] == [2, 0])

    native = full_face_refinements()["candidates"][0]["values"]
    native_points = [sp.Matrix(v) for v in BIT_VECTORS if any(v) and native[index(v)] == 0]
    change = mod2(sp.Matrix.hstack(*native_points[:4])*sp.Matrix.hstack(*map(sp.Matrix, label_order[:4])).inv())
    check("explicit map to the existing unmarked four-face code", len(native_points) == 5
          and all(native[index(mod2(change*sp.Matrix(v)))] == quadratic(v) for v in BIT_VECTORS))

    quarter = (identity+linear_j)/sp.sqrt(2)
    paired_quarter = tensor(quarter, quarter)
    permutation = [next(i for i, target in enumerate(gamma) if any(
        same(paired_quarter*g*paired_quarter.H, sign*target) for sign in (1, -1)
    )) for g in gamma]
    check("common local quarter has one fixed projective axis and two swaps", permutation == [1, 0, 3, 2, 4]
          and same(paired_quarter**4, eye4) and not same(paired_quarter**2, eye4))

    # Reflection normals lift the three adjacent transpositions of (0 1 2 3).
    permutation_lift = sp.prod([(gamma[k]-gamma[k+1])/sp.sqrt(2) for k in range(3)])
    cycle = [1, 2, 3, 0, 4]
    result = {"fourth_power_minus_identity": 0, "fourth_power_noncentral": 0}
    example = None
    for signs in product((1, -1), repeat=4):
        parity = sp.prod(signs)
        signed = (*signs, -parity)
        output_flips = [gamma[cycle[k]] for k, sign in enumerate(signed) if sign == -1]
        lift = (sp.prod(output_flips)*permutation_lift).applyfunc(sp.simplify)
        check(f"unitary quaternionic lift with cycle signs {signs}", same(lift*lift.H, eye4)
              and same(theta*lift.conjugate()*theta.inv(), lift) and all(
            same(lift*g*lift.H, sign*gamma[target]) for g, target, sign in zip(gamma, cycle, signed)
        ))
        fourth = (lift**4).applyfunc(sp.simplify)
        if parity == 1:
            check(f"central spin sign with cycle signs {signs}", fourth == -eye4)
            result["fourth_power_minus_identity"] += 1
            example = lift
        else:
            check(f"noncentral spin fourth power with cycle signs {signs}",
                  not same(fourth*gamma[0]*fourth.H, gamma[0]))
            result["fourth_power_noncentral"] += 1

    # A real multiplicity doublet cancels the central case without
    # altering the action on the five Clifford axes or quaternionic type.
    auxiliary = sp.Matrix([[1, -1], [1, 1]])/sp.sqrt(2)
    total = tensor(example, auxiliary)
    total_theta = tensor(theta, identity)
    check("one real multiplicity doublet cancels the spin sign", same(total**4, sp.eye(8))
          and not same(total**2, sp.eye(8)) and same(total_theta*total.conjugate()*total_theta.inv(), total))
    check("the extra doublet leaves the projective four-cycle unchanged", all(any(
        same(total*tensor(g, identity)*total.H, sign*tensor(gamma[target], identity)) for sign in (1, -1)
    ) for g, target in zip(gamma, cycle)))

    # Keep the 4 P1 rays and the 4 P2 directions in their actual operator
    # spaces: the former have a linear relation that the latter do not.
    spatial_cycle = (marks[0], marks[2], marks[1], marks[3])
    fourgram = sp.Matrix(4, 4, lambda i, j: sp.trace(spatial_cycle[i]*spatial_cycle[j]))
    check("a linear single-copy map cannot transport the four independent code directions",
          fourgram.rank() == 3 and spatial_cycle[0]+spatial_cycle[2] == spatial_cycle[1]+spatial_cycle[3]
          and sp.Matrix.hstack(*(g.reshape(16, 1) for g in gamma[:4])).rank() == 4)
    squared_gram = fourgram.applyfunc(lambda value: value**2)
    check("quadratic mark tensors supply four independent directions",
          squared_gram.eigenvals() == {sp.Rational(3, 2): 1, 1: 2, sp.Rational(1, 2): 1})
    check("local rotation cycles the four marked rays", all(
        same(quarter*mark*quarter.H, spatial_cycle[(i-1) % 4]) for i, mark in enumerate(spatial_cycle)
    ))
    target_effects = [(eye4+operator)/2 for operator in gamma[:4]]
    check("opposite mark effects obey an order bound that the target sharp pair violates",
          spatial_cycle[0]*spatial_cycle[2] == sp.zeros(2)
          and (target_effects[0]+target_effects[2]-eye4).eigenvals()
          == {sp.sqrt(2)/2: 2, -sp.sqrt(2)/2: 2})
    check("the product of adjacent sharp effects also has incompatible corner constants",
          target_effects[0]*target_effects[1]*target_effects[0] == target_effects[0]/2
          and squared_gram[0, 1] == sp.Rational(1, 4))

    return {
        "status": "EXACT_SCOPED_CONTROLS_PASSED_F1_OPEN",
        "checks": checks,
        "four_letter_witness": {"anticommuting_target": -1, "global_trace_control": "-7/25"},
        "fixed_traceless_dimensions": sizes,
        "same_window_and_bilocal_dimensions": [2, 3],
        "native_unmarked_quadratic_isometry": [list(map(int, row)) for row in change.tolist()],
        "local_quarter_permutation": permutation,
        "single_block_signed_four_cycle": result,
        "explicit_doubled_real_lift_dimension": 8,
        "single_and_quadratic_mark_span_ranks": [3, 4],
        "physical_solution": False,
    }
