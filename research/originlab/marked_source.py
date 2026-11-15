"""Finite affine tests for a proposed identification of the P1 source."""

from collections import Counter
from itertools import product
from math import isqrt

import sympy as sp

from affine_operators import MARKED_GAMMA, MARKED_P1, MARKED_P2, ROOTS, charge_signature, components, dot, neg
from reconstruction import require


def _glue_degree(vector):
    glue, colour, _ = charge_signature(vector)
    # In the marked chart p2 has (glue, colour) = (1, 2).
    return (9 * glue + 8 * colour) % 12


def _shells(limit):
    """Count E8 lattice charges by norm and degree in doubled coordinates."""
    bound = isqrt(8 * limit)

    def extend(prefix, square, remaining, parity):
        if not remaining:
            if sum(prefix) % 4 == 0:
                yield square // 8, _glue_degree(prefix)
            return
        for entry in range(-bound, bound + 1):
            if entry % 2 == parity and square + entry**2 <= 8 * limit:
                yield from extend((*prefix, entry), square + entry**2, remaining - 1, parity)

    return Counter(pair for parity in (0, 1) for pair in extend((), 0, 8, parity))


def _character():
    shells = _shells(3)
    # Eight independent Heisenberg oscillators at each positive integer mode.
    oscillators = [1, 0, 0, 0]
    for mode in range(1, 4):
        for _ in range(8):
            for grade in range(mode, 4):
                oscillators[grade] += oscillators[grade - mode]
    counts = [[sum(shells[charge_grade, degree] * oscillators[grade - charge_grade]
                   for charge_grade in range(grade + 1)) for degree in range(12)]
              for grade in range(4)]
    degrees = [_glue_degree(root) for root in (MARKED_P2, MARKED_GAMMA, neg(MARKED_P1))]
    require(degrees == [1, 2, 9] and
            [sum(row) for row in counts] == [1, 248, 4124, 34752] and
            Counter(_glue_degree(root) for root in ROOTS) ==
            Counter({k: n - 8 * (k == 0) for k, n in enumerate(counts[1])}),
            "marked lattice character failed the affine grade or root census")
    folded = [[row[k] + row[k + 6] for k in range(6)] for row in counts]
    require(folded[1] == [54, 48, 33, 32, 33, 48] and
            folded[2] == [762, 720, 657, 608, 657, 720],
            "the joint parity and colour character changed")
    neutral_components = components(root for root in ROOTS if _glue_degree(root) == 0)
    order_two_weights = Counter(
        (str(charge_signature(root)[2]), (root[3] - root[4]) // 2)
        for root in ROOTS if _glue_degree(root) == 6
    )
    require(neutral_components == ["A1", "A1", "A2", "A3"] and order_two_weights ==
            {(charge, weak): 6 for charge in ("-1/2", "1/2") for weak in (-1, 1)},
            "the coarse-neutral class contains the wrong root actions")
    return {"marked_glue_path_mod_12": [0, 1, 3, 0],
            "marked_step_degrees_mod_12": degrees,
            "neutral_algebra_mod_12": [*neutral_components, "u1"],
            "class_six_root_weights": [
                {"hypercharge": charge, "weak_doubled_weight": weak, "multiplicity": multiplicity}
                for (charge, weak), multiplicity in sorted(order_two_weights.items())
            ],
            "affine_grade_counts_mod_12": counts,
            "affine_grade_counts_mod_6": folded}


def _metric_constraints(neutral_orbits, colour_simple, bracket):
    names = ("E6", "A2", "colour_1", "colour_2")

    def sector(root):
        _, colour, _ = charge_signature(root)
        if colour:
            return colour + 1
        return int(any(dot(root, alpha) for alpha in colour_simple))

    require(neutral_orbits == (81, 81) and
            Counter(map(sector, ROOTS)) == {0: 72, 1: 6, 2: 81, 3: 81},
            "neutral covariance must have four multiplicity-free summands")
    links = {}
    for root in sorted(ROOTS):
        target, coefficient = bracket(MARKED_P2, root)
        if coefficient:
            source_sector, target_sector = sector(root), sector(target)
            if source_sector != target_sector:
                links.setdefault((source_sector, target_sector), (root, target))
    unit = sp.eye(len(names))
    constraints = sp.Matrix.vstack(*(unit.row(j) - unit.row(i) for i, j in sorted(links)))
    kernel = constraints.nullspace()
    require(kernel == [sp.ones(4, 1)], "the charged Ward adjoint must fix the relative grade-one metric")
    nonuniform_metric = sp.Matrix([1, 1, 1, 2])
    require(any(constraints * nonuniform_metric), "the charged constraint must reject a neutral-only metric")
    return {"summands": dict(zip(names, (78, 8, 81, 81))),
            "neutral_covariant_hermitian_parameters": len(names),
            "with_marked_charged_adjoint_hermitian_parameters": len(kernel),
            "charged_ward_constraint_rank": constraints.rank(),
            "nonuniform_neutral_metric_rejected_by_charged_adjoint": True,
            "charged_root_witnesses": [
                {"from": names[i], "to": names[j], "input_doubled_root": root,
                 "output_doubled_root": target}
                for (i, j), (root, target) in sorted(links.items())
            ],
            "vacuum_normalization_requires": "level-one creation/annihilation and the same bound empty reference"}


def _clock_lift(character):
    """Compare the polarization square-root construction to the native glue action."""
    reflection = sp.Matrix([[0, 1], [1, 0]])
    polarization = sp.diag(1, -1)
    identity = sp.eye(2)
    clock = reflection * polarization
    half_step = (identity + clock) / sp.sqrt(2)
    require(half_step**2 == clock and half_step * polarization * half_step.H == reflection,
            "the half-step must rotate the two polarizations")
    # Rephase so that the lift, rather than its conjugation alone, has order four.
    lift = ((1-sp.I) * half_step / sp.sqrt(2)).applyfunc(sp.expand)
    tensor_lift = sp.kronecker_product(lift, lift)

    def adjoint_dimensions(eigenvalues):
        return [sum(m*n for (left, m), (right, n) in product(eigenvalues.items(), repeat=2)
                    if sp.simplify(left - sp.I**grade * right) == 0) for grade in range(4)]

    lift_eigenvalues = tensor_lift.eigenvals()
    one_cell = adjoint_dimensions(lift.eigenvals())
    two_cells = adjoint_dimensions(lift_eigenvalues)
    grade_one = Counter(charge_signature(root)[0] for root in ROOTS)
    grade_one[0] += 8
    native = [grade_one[g] for g in range(4)]
    counts = character["affine_grade_counts_mod_12"]
    require(native == [sum(counts[1][g::4]) for g in range(4)] and
            one_cell == [2, 1, 0, 1] and two_cells == [6, 4, 2, 4],
            "the native or reconstructed quarter action changed")
    traces = [sum((-1)**degree * value for degree, value in enumerate(row)) for row in counts]
    require(traces == [1, -8, 28, -64], "the native glue parity trace changed")
    # A spectral embedding can leave the identified grade space invariant under
    # the clock without leaving it invariant under both elementary operations.
    needed_per_parity = max(sum(native[::2]), sum(native[1::2]))
    multiplicity = needed_per_parity // tensor_lift.rows
    balanced_cover = [multiplicity * (lift_eigenvalues.get(sp.I**g,0) +
                        lift_eigenvalues.get(-sp.I**g,0)) for g in range(4)]
    require(sum(balanced_cover) == 2*needed_per_parity and
            all(cover >= target for cover,target in zip(balanced_cover,native)),
            "the phase-shifted tensor copies must attain the spectral lower bound")
    xx = sp.kronecker_product(reflection, reflection)
    zz = sp.kronecker_product(polarization, polarization)
    double_identity = sp.eye(4)
    vacuum = (double_identity + xx) * (double_identity + zz) / 4
    require(vacuum.rank() == 1 and xx*zz == zz*xx and
            (tensor_lift*xx*tensor_lift.H - zz).applyfunc(sp.expand).is_zero_matrix and
            (tensor_lift*zz*tensor_lift.H - xx).applyfunc(sp.expand).is_zero_matrix,
            "the tensor clock must preserve the Bell vacuum and exchange the stabilizers")
    bell_energies = double_identity - (xx + zz)/2
    require(bell_energies.eigenvals() == {0:1,1:2,2:1},
            "the normalized additive Bell generator has the wrong spectrum")
    return {
        "status": "EXCLUDED_SCOPED_POLARIZATION_LIFTS",
        "one_cell_adjoint_dimensions": one_cell,
        "two_cell_adjoint_dimensions": two_cells,
        "native_grade_one_glue_dimensions": native,
        "native_affine_grade_parity_traces": traces,
        "minimum_grade_one_nullity_of_parity_reversing_operator":abs(traces[1]),
        "one_cell_obstruction": "The half-step lacks adjoint phase -1, present on native grade-one currents.",
        "two_cell_obstruction": "The squared tensor lift changes sign by conjugation with the first cell polarization. Its parity trace is zero on any common invariant grade space, unlike native grade one.",
        "two_cell_scope": "Commuting cells with a first-cell polarization preserving the identified native grade-one space. Relative vacuum phase can change the sign of the parity trace, not make it zero.",
        "minimal_balanced_grade_one_spectral_cover_dimension":sum(balanced_cover),
        "balanced_grade_one_spectral_cover_dimensions":balanced_cover,
        "spectral_cover_scope":"Copies with a central sign shift match multiplicities only. They carry no supplied current, reflection or source intertwiner for the embedded subspace.",
        "bell_clock_control":{
            "vacuum_rank":vacuum.rank(), "normalized_additive_energies":[0,1,1,2],
            "general_stabilizer_diagonal_clock_invariant_energies":["0","a","a","b"],
            "scope":"Clock covariance exchanges the two single stabilizer excitations. It requires b=2a only for an additive stabilizer Hamiltonian; the general diagonal family retains b/a."},
    }


def controls(neutral_orbits, colour_simple, bracket):
    """Check an identified affine grade, before the physical time insertion."""
    character = _character()
    return {"scope": "fixed level-one affine vacuum module, before time and selection",
            "character": character,
            "clock_lift": _clock_lift(character),
            "equal_time_grade_one_metric": _metric_constraints(neutral_orbits, colour_simple, bracket),
            "physical_source_identity_proved": False}
