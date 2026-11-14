"""Time and charge tests of the marked E8 word on its original affine vacuum."""

from itertools import combinations, product
from math import comb
from pathlib import Path

import sympy as sp

from original_paths import original_word
from reconstruction import require, scalar
from source_return import REGIONAL
from verify import digest


CARTAN = Path("docs/research/fundamental-plan/TFPT-Cartan-Prozess-2026-11-12.md")
PAIR = Path("docs/research/fundamental-plan/TFPT-R113-GELADENE-PAAR-RUECKKEHR-20261112.txt")
COLLAR = Path("docs/research/fundamental-plan/check_tfpt_charged_collar_sewing_20261111.py")
MODULAR = Path("docs/research/fundamental-plan/TFPT-MODULARE-VIERSCHNITTE-E8-20261111.txt")
SOURCES = {
    CARTAN: "1f05e02b144af965776f7c6fd1047707ad13d73f91730e550b5b3cd7ddd5fc18",
    PAIR: "9798e2bf3de6c9042fcad177b9032c973fe977f1cb1b57c435197a91391e8d8d",
    COLLAR: "1ef8234a6fceb9033c89ca6c34dc721c40fe39bc18bdc76a3516067fb098feb9",
    REGIONAL: "6888980fcb80f292c4da087a5bc3f486395097be914794cb22b6ce876f40c64d",
    MODULAR: "d3d5496c527a8a4446ee113e969fa59400a752937a6007336ec8eade2de36571",
}
LEGS = ("J^(-p1)", "J^gamma", "J^p2")
GLUE = (1, 2, 1)
TRIALITY = (0, 1, 2)
HYPERCHARGES = ("0", "-1/3", "1/3")


def marked_charges():
    # The fixed physical D8 chart, in the left-to-right order of insertions.
    p1 = sp.ones(8, 1)/2
    p2 = sp.Matrix([-1, 1, 1, 1, 1, -1, 1, 1])/2
    weights = (-p1, p1-p2, p2)
    require(sum(weights, sp.zeros(8, 1)) == sp.zeros(8, 1), "marked word must be neutral")
    rows = []
    for leg, weight, glue, triality, hypercharge in zip(LEGS, weights, GLUE, TRIALITY, HYPERCHARGES):
        y = weight.dot(sp.Matrix([-2, -2, -2, 3, 3, 0, 0, 0]))/6
        require(weight.dot(weight) == 2 and (2*sum(weight[5:])) % 4 == glue and
                (weight[0]+weight[1]-2*weight[2]) % 3 == triality and str(y) == hypercharge,
                "wrong marked charge in the physical D8 chart")
        rows.append({"field": leg, "root": [str(w) for w in weight], "Y": hypercharge,
                     "color_triality": triality, "internal_glue_class": glue, "D_phase": str(sp.I**glue)})
    return rows


def radial_word(x, y):
    return -sp.I*x*y/((1-sp.I*x)*(sp.I+y)*(1+x*y))


def radial_grid(left, right):
    return sp.Matrix([[sp.cancel(radial_word(x, y), extension=sp.I) for y in right] for x in left])


def spectral_check(size):
    x, y = sp.symbols("x y")
    kernel = sp.Matrix(size, size, lambda m, n: min(m+1, n+1))
    triangle = sp.Matrix(size, size, lambda m, n: 1 if n <= m else 0)
    series = -sum(kernel[m, n]*sp.I**(m+n)*x**(m+1)*y**(n+1)
                  for m in range(size) for n in range(size))
    residue = sp.Poly(sp.expand(series*(1-sp.I*x)*(1-sp.I*y)*(1+x*y)+x*y), x, y)
    require(kernel == triangle*triangle.T and
            all(m > size or n > size for (m, n), _ in residue.terms()),
            "double radial coefficients do not reproduce the original kernel")
    points = [sp.Rational(1, j+2) for j in range(size)]
    # Strip only the nonzero row and column factors of the original word.
    grid = sp.Matrix([[1/(1+x_i*y_j) for y_j in points] for x_i in points])
    vandermonde = sp.prod(b-a for a, b in combinations(points, 2))
    determinant = (-1)**comb(size, 2)*vandermonde**2/sp.prod(1+a*b for a in points for b in points)
    require(grid.det() == determinant != 0, "Cauchy time grid loses rank")
    return {"coefficient_law": "-i**(m+n-2)*min(m,n), m,n >= 1",
            "row_column_phase_stripped_rank": size, "coefficient_block_determinant": str(kernel.det()),
            "distinct_positive_time_grid_rank": size, "analytic_separation_rank": "infinite"}


def parabolic_action(point, time):
    conjugate = sp.conjugate(point)
    return sp.Matrix([[1-time/2, point*time/2], [-conjugate*time/2, 1+time/2]])


def balanced_clock():
    a = sp.symbols("a", real=True)
    b = 2*a/(1-a*a)
    first, second = parabolic_action(1, a), parabolic_action(-1, b)
    rotation = sp.diag((1-a)/(1+a), (1+a)/(1-a))
    require((first*second*first-rotation).applyfunc(sp.cancel) == sp.zeros(2),
            "positive three-step composition does not yield the radial transfer")
    require((parabolic_action(-1, a)*parabolic_action(1, b)*parabolic_action(-1, a)-rotation)
            .applyfunc(sp.cancel) == sp.zeros(2), "opposite orientations must give the same transfer")
    c, d = sp.symbols("c d", real=True)
    two = parabolic_action(1, c)*parabolic_action(-1, d)
    require(sp.expand(two[0, 1]+two[1, 0]) == c*d and
            sp.expand(two[0, 1]-two[1, 0]) == c-d,
            "two-step minimum test changed")
    m = sp.symbols("m", real=True)
    sandwich = parabolic_action(1, c)*parabolic_action(-1, m)*first
    require(sp.expand(sandwich[0, 1]+sandwich[1, 0]) == m*c-m*a and
            sp.expand(sandwich[0, 1]-sandwich[1, 0]) == a*m*c+a+c-m,
            "the diagonal composition must balance the outer steps")
    geometric_quarter = sp.diag(sp.I, 1)
    require(geometric_quarter*parabolic_action(sp.I, a)*geometric_quarter.inv() ==
            parabolic_action(-1, a), "quarter transport of the marked inclusion changed")
    # After sending the shared endpoints to infinity and zero, then
    # changing coordinate w=i*u, the adjacent generators are real matrices.
    coordinate = sp.diag(-sp.I, 1)*sp.Matrix([[1-sp.I, -1-sp.I], [1, -1]])
    require((coordinate*parabolic_action(1, a)*coordinate.inv()).applyfunc(sp.simplify) ==
            sp.Matrix([[1, a], [0, 1]]) and
            (coordinate*parabolic_action(sp.I, a)*coordinate.inv()).applyfunc(sp.simplify) ==
            sp.Matrix([[1, 0], [a/2, 1]]), "adjacent positive generators lost their common real form")
    scale = sp.symbols("q", real=True, positive=True)
    target = coordinate*sp.diag(scale, 1/scale)*coordinate.inv()
    require(sp.simplify(sp.im(target[0, 0])-(scale-1/scale)/2) == 0,
            "radial transfer must lie outside the adjacent generators' real form")
    scale = sp.S.Half
    # In the response at x=1/2 the marked width is log(2), so the
    # two positive parabolic step lengths are rational.
    outer, middle = (1-scale)/(1+scale), (1/scale-scale)/2
    composed = parabolic_action(1, outer)*parabolic_action(-1, middle)*parabolic_action(1, outer)
    require(composed == sp.diag(scale, 1/scale), "dyadic marked width composition changed")
    return {"transfer": "exp(-a*P_+) exp(-b*P_-) exp(-a*P_+) = exp(-t*L0)",
            "steps_in_radial_width_t": {"outer": "tanh(t/4)", "middle": "sinh(t/2)"},
            "reversed_signs_also_work": True,
            "opposite_pair": ["1", "-1"], "other_quarter_pair": ["I", "-I"],
            "original_inclusion_pair": ["1", "I"],
            "path_to_antipodal_pair": "conjugate the positive translation at I by the geometric quarter",
            "adjacent_pair_alone_matches_radial_mobius_element": False,
            "minimum_number_of_positive_segments": 3,
            "radial_factor_1_over_4_steps": [str(outer), str(middle), str(outer)],
            "independent_spatial_directions_added": 0,
            "requires_same_original_mobius_vacuum_and_positive_translation_pair": True}


def clock_checks():
    s, t = sp.symbols("s t", real=True)
    value = sp.cancel(radial_word(1, 1), extension=sp.I)
    flows = {}
    fixed_mark_words = {-sp.S.One: -1/(4*(s-sp.I)), sp.S.One: -1/(4*(t-sp.I)),
                        sp.I: 1/(2*(s+t+2*sp.I))}
    for point, expected in fixed_mark_words.items():
        # The inclusion pair shares 1 and i as respective endpoints. The
        # antipodal partner at -1 is available by geometric quarter covariance.
        def move(mark, time):
            top, den = parabolic_action(point, time)*sp.Matrix([mark, 1])
            return top/den, 1/den**2

        mid, mid_jacobian = move(sp.I, s)
        end, end_jacobian = move(-sp.S.One, s+t)
        word = sp.cancel(-sp.I*mid_jacobian*end_jacobian/((1-mid)*(mid-end)*(1-end)), extension=sp.I)
        require(sp.simplify(word-expected) == 0, "translated original word changed")
        derivatives = [sp.simplify(sp.diff(word, time).subs({s: 0, t: 0})) for time in (s, t)]
        flows[str(point)] = {"word": str(expected), "first_derivatives": list(map(str, derivatives))}

    # On a weight-one vector, P = L0-(L1+L_-1)/2 has Jacobi moments.
    # They agree with 4 E exp(-2 E) dE on the positive real axis.
    v = [sp.S.One]
    for k in range(1, 7):
        moved = [sp.S.Zero]*(len(v)+1)
        for n, amplitude in enumerate(v):
            moved[n] += (n+1)*amplitude
            moved[n+1] -= amplitude/2
            if n:
                moved[n-1] -= n*(n+1)*amplitude/2
        v = moved
        require(v[0] == sp.factorial(k+1)/2**k, "positive translation moments changed")

    x, y = sp.symbols("x y")
    word = radial_word(x, y)
    radial_derivatives = [sp.simplify((-q*sp.diff(word, q)).subs({x: 1, y: 1})) for q in (x, y)]
    require(radial_derivatives == [-sp.Rational(1, 8)]*2 and value == -sp.I/4,
            "the original point word or derivatives changed")
    steps = [(1+time/2)**-2 for time in (sp.S.One, sp.S(2))]
    return {"radial": {"zero_width": str(value), "first_derivatives": list(map(str, radial_derivatives)),
                       "on_normalized_current": "exp(-t)", "rank_one_on_grade_one_current": True},
            "parabolic_on_original_marks": flows,
            "parabolic_on_normalized_current": {"response": "(1+t/2)**(-2)",
                "spectral_measure_density": "4*E*exp(-2*E)", "first_generator_moment": "1",
                "second_generator_moment": "3/2", "steps_1_and_2": list(map(str, steps)),
                "held_out_rank_one_defect": str(steps[1]-steps[0]**2)},
            "balanced_positive_radial_composition": balanced_clock(),
            "single_fixed_mark_parabolic_matches_full_radial_word_under_common_rescaling": False}


def quadratic_floor():
    # Recompute the minimum using the D8+ roots and Weyl's dimension formula.
    positives = []
    for j, k in combinations(range(8), 2):
        for sign in (1, -1):
            root = [0]*8
            root[j], root[k] = 2, 2*sign
            positives.append(tuple(root))
    positives.extend((1, *tail) for tail in product((-1, 1), repeat=7) if tail.count(-1) % 2 == 0)
    positive_set = set(positives)
    simple = [r for r in positives if not any(tuple(x-y for x, y in zip(r, a)) in positive_set
                                              for a in positives)]
    require(len(simple) == 8 and len(positive_set) == 120, "wrong E8 positive or simple system")
    # Components in the simple-root basis are the pairings with the
    # fundamental weights since the root system is simply laced.
    simple_basis = sp.Matrix(simple).T/2
    root_basis = simple_basis.inv()
    coefficients = [tuple(root_basis*sp.Matrix(p)/2) for p in positives]
    require(all(all(c.is_Integer and c >= 0 for c in row) for row in coefficients), "wrong positive system")
    dimensions = [sp.prod((sum(p)+p[i])/sum(p) for p in coefficients) for i in range(8)]
    require(sorted(dimensions) == [248, 3875, 30380, 147250, 2450240, 6696000, 146325270, 6899079264],
            "unexpected E8 fundamental dimensions")
    least_dimension = int(min(dimensions))
    weights = root_basis.T
    weyl_vector = weights*sp.ones(8, 1)
    require(min(weights.T*weights) > 0, "Casimir monotonicity needs the positive inverse Cartan matrix")
    indices = {}
    for dim, weight in zip(dimensions, weights.columnspace()):
        index = dim*weight.dot(weight+2*weyl_vector)/248
        require(index.is_even, "quadratic trace index must be an even integer")
        indices[int(dim)] = int(index)
    require(indices[248] == min(indices.values()) == 60, "wrong fundamental quadratic trace indices")
    return {"fundamental_representation_dimensions": sorted(map(int, dimensions)),
            "fundamental_trace_index_by_dimension": indices,
            "least_nontrivial_complex_dimension": least_dimension,
            "lower_bound_complex_car_modes_for_scalar_quadratic_embedding": least_dimension//2,
            "minimum_level_for_linear_majorana_action_and_quadratic_free_currents": min(indices.values())//2,
            "required_original_affine_level": 1,
            "existing_fiber_modes": 96, "existing_fiber_majoranas": 192,
            "full_e8_by_number_preserving_forms_on_the_192_dimensional_doubled_fibre": False,
            "full_e8_on_single_96_mode_fiber_by_scalar_quadratics": False}


def original_time(root):
    original = original_word(root)
    refs = {str(p): digest(root/p) for p in SOURCES}
    require(all(refs[str(p)] == checksum for p, checksum in SOURCES.items()),
            "review changed affine time or charged source principles")
    charges = marked_charges()
    require(charges[1]["Y"] == original["r113_binding"]["charge_of_original_middle"]["Y"],
            "marked middle no longer agrees with the original path")
    return {"status": "EXACT_ORIGINAL_AFFINE_TIME_TESTS", "source_refs": {**original["source_refs"], **refs},
            "marked_charges": charges, "right_to_left_color_path": [0, 2, 0, 0],
            "radial_spectrum": spectral_check(4), "clocks": clock_checks(),
            "finite_degree_two_restriction": quadratic_floor(), "weighted_r113_values": None}


def checked_marked_time(data, source_id, time_id, contract, clock_root=None):
    """Compare independently supplied two-time forms through the original endpoint."""
    require(data.get("source_id") == source_id and data.get("time_id") == time_id and
            all(isinstance(data.get(k), str) and data[k].strip() for k in (
                "native_field_and_charge_binding", "physical_to_radial_time_binding", "absolute_reference_binding"
            )), "marked time grid must bind the native charge, physical time and absolute reference")
    require(data.get("native_circle_marks") == ["1", "I", "-1"] and
            data.get("ordered_fields") == list(LEGS) and data.get("internal_glue_classes") == list(GLUE) and
            data.get("color_trialities") == list(TRIALITY) and data.get("hypercharges") == list(HYPERCHARGES),
            "marked time grid must use the signed original three legs")
    axes = [data[key] for key in ("left_axis", "right_axis")]
    coordinates, time_grid = [], []
    for axis in axes:
        times = [scalar(row["physical_time"]) for row in axis]
        ratios = [scalar(row["radial_factor"]) for row in axis]
        require(len(times) >= 2 and all(t.is_Rational and t > 0 for t in times) and
                all(a < b for a, b in zip(times, times[1:])) and
                all(x.is_Rational and 0 < x < 1 for x in ratios) and
                all(a > b for a, b in zip(ratios, ratios[1:])),
                "need at least two strictly ordered time pairs with exact decreasing radial factors")
        coordinates.append(ratios)
        time_grid.append(times)
    forms = data["ordered_word_forms"]
    require(len(forms) == len(axes[0]) and all(len(row) == len(axes[1]) for row in forms),
            "the ordered forms must cover the entire marked time grid")
    observed = sp.Matrix([[contract(form) for form in row] for row in forms])
    expected = radial_grid(*coordinates)
    mismatches = [[m, n] for m in range(expected.rows) for n in range(expected.cols)
                  if sp.cancel(observed[m, n]-expected[m, n], extension=sp.I) != 0]
    result = {"status": "REFUTED_SUPPLIED_RADIAL_TIME_GRID" if mismatches else "EXACT_SUPPLIED_RADIAL_TIME_GRID",
            "mismatch_indices": mismatches, "supplied_rank": observed.rank(), "affine_rank": expected.rank(),
            "observed": [[str(x) for x in row] for row in observed.tolist()],
            "affine": [[str(x) for x in row] for row in expected.tolist()],
            "physical_source_identity_proved": False}
    from folded_time import observed_two_time_grid
    result.update(observed_two_time_grid(observed,expected,coordinates,time_grid,clock_root))
    return result
