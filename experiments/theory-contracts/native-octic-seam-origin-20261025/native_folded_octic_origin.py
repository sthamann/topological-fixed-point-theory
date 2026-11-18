#!/usr/bin/env python3
"""Exact folded-source overlap with the native normalized octic.

This checker rebuilds the original O/M/U+ transport.  It does not fit a
Cartan direction to the octic.  The primary result is the level-eight Ward
matrix element in both original quarter sectors.  A bounded secondary check
records the first full-root-OPE overlap without replacing its full state by
the holomorphic polynomial component.
"""

from __future__ import annotations

import itertools as it
import json
from collections import Counter
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]

import sys

if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tfpt_explorer import cartan_source as cs  # noqa: E402


def need(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def clean(value):
    if isinstance(value, sp.MatrixBase):
        return value.applyfunc(sp.expand)
    return sp.expand(value)


def fock_norm_squared(polynomial: sp.Expr, variables: tuple[sp.Symbol, ...]):
    """Positive Heisenberg-Fock norm: ||x^multi||^2=product multi!."""
    poly = sp.Poly(sp.expand(polynomial), *variables)
    return sp.simplify(sum(
        sp.conjugate(coefficient) * coefficient
        * sp.prod(sp.factorial(power) for power in powers)
        for powers, coefficient in poly.terms()
    ))


def mode_partitions(weight: int, mode: int = 8):
    """All occupation vectors with sum n*m_n=weight, without mode truncation."""
    if mode == 0:
        if weight == 0:
            yield ()
        return
    for count in range(weight // mode + 1):
        for lower in mode_partitions(weight - mode * count, mode - 1):
            yield lower + (count,)


def full_translation_descendants(d_squared: sp.Expr):
    """Keep all oscillator modes of (L_-1+d.h_-1)^k Omega, k<=8.

    In y_n=d.h_-n, L_-1 is the derivation y_n -> n*y_(n+1).
    The norm uses [d.h_n,d.h_-n]=n*d_squared.  This is a second,
    full-word check of the sl2 norm, independent of its recursion.
    """
    y = sp.symbols("y1:9")
    state = sp.Integer(1)
    all_levels = {}
    for level in range(9):
        partitions = list(mode_partitions(level))
        bell_state = sp.expand(sum(
            sp.factorial(level) * sp.prod(
                y[mode - 1] ** occupation
                / (sp.Integer(mode) ** occupation * sp.factorial(occupation))
                for mode, occupation in enumerate(powers, 1)
            )
            for powers in partitions
        ))
        need(sp.expand(state - bell_state) == 0,
             f"actual translation word and complete Bell coefficient differ at {level}")
        poly = sp.Poly(state, *y)
        norm = sp.simplify(sum(
            coefficient**2 * sp.prod(
                sp.factorial(occupation) * (mode * d_squared) ** occupation
                for mode, occupation in enumerate(powers, 1)
            )
            for powers, coefficient in poly.terms()
        ))
        need(norm == sp.factorial(level) * sp.rf(d_squared, level),
             f"full oscillator norm differs from sl2 at level {level}")
        need(len(poly.terms()) == sp.partition(level),
             f"an oscillator partition was discarded at level {level}")
        all_levels[str(level)] = {
            "oscillator_partitions": len(partitions),
            "full_fock_norm_squared": str(norm),
            "pure_mode_minus1_coefficient": str(poly.coeff_monomial(y[0] ** level)),
        }
        if level < 8:
            state = sp.expand(
                y[0] * state
                + sum(mode * y[mode] * sp.diff(state, y[mode - 1])
                      for mode in range(1, 8))
            )
    # At zero shift y_n=0: the unshifted vacuum is translation invariant.
    need(state.subs(dict.fromkeys(y, 0)) == 0,
         "zero-shift vacuum control failed")
    return all_levels


def original_transport():
    lattice = cs._lattice_data()
    target = sp.Matrix.hstack(*map(sp.Matrix, cs._TARGET_SIMPLE))
    orthogonal = clean(target * lattice["basis"].inv() / sp.sqrt(2))
    positive = sp.zeros(8, 4)
    for column in range(4):
        positive[2 * column, column] = 1 / sp.sqrt(2)
        positive[2 * column + 1, column] = -sp.I / sp.sqrt(2)
    fold_path = (
        REPO
        / "experiments/theory-contracts/v16-operator-pullback-20261014"
        / "fold-results.json"
    )
    fold = sp.Matrix(json.loads(fold_path.read_text())["explicit_Weyl_matrix"])
    need(clean(fold.T * fold) == sp.eye(8), "original M is not orthogonal")
    need(clean(orthogonal.T * orthogonal) == sp.eye(8), "original O is not orthogonal")
    need(clean(positive.H * positive) == sp.eye(4), "original U+ is not Hermitian-isometric")
    need(clean(positive.T * positive) == sp.zeros(4), "original U+ is not isotropic")
    return lattice, orthogonal, positive, fold


def main() -> None:
    lattice, orthogonal, positive, fold = original_transport()
    source_to_positive = clean(positive.T * orthogonal.T * fold)
    need(clean(source_to_positive * source_to_positive.H) == sp.eye(4),
         "original CAR current map does not retain the positive metric")
    need(clean(source_to_positive * source_to_positive.T) == sp.zeros(4),
         "original CAR current map does not retain isotropic locality")
    x = sp.symbols("x0:4")
    i8 = sp.expand(
        sum(value**8 for value in x)
        + 14 * sum(x[a]**4 * x[b]**4 for a, b in it.combinations(range(4), 2))
        + 168 * sp.prod(value**2 for value in x)
    )
    f8 = i8 / (96 * sp.sqrt(140))
    need(fock_norm_squared(f8, x) == 1, "native F8 normalization changed")
    need(sp.expand(sum(value * sp.diff(f8, value) for value in x) - 8 * f8) == 0,
         "octic relative oscillator energy changed")

    expected_i8 = {
        sp.Rational(1, 4): sp.Rational(5416993, 2**32) - sp.I * sp.Rational(584073, 2**27),
        sp.Rational(3, 4): -sp.Rational(3356639, 2**32) - sp.I * sp.Rational(495495, 2**27),
    }
    expected_projection_fraction = {
        sp.Rational(1, 4): sp.Rational(2611534296241, 59248596692359104307200),
        sp.Rational(3, 4): sp.Rational(23879517043811, 781004229126551829504000),
    }
    complete_words = full_translation_descendants(sp.Rational(85, 128))
    cases = {}
    for alpha in (sp.Rational(1, 4), sp.Rational(3, 4)):
        delta = sp.Matrix([(index + alpha) / 8 - sp.Rational(1, 2) for index in range(8)])
        d = clean(fold * delta)
        b = clean(positive.T * orthogonal.T * d)
        need(clean(source_to_positive * delta - b) == sp.zeros(4, 1),
             "direct original CAR current map differs from the M/O route")
        h = sp.simplify((d.T * d)[0] / 2)
        need(h == sp.Rational(85, 256), "original folded conformal weight changed")
        need(sp.simplify((b.H * b)[0]) == h, "positive Cartan norm is not half of d squared")

        value_i8 = sp.simplify(i8.subs(dict(zip(x, b))))
        need(value_i8 == expected_i8[alpha], "original O/M octic evaluation changed")
        onepoint = sp.simplify(f8.subs(dict(zip(x, b))))
        # The Hilbert product is linear in its second slot.  The original
        # fold kernel therefore uses the conjugate Riesz coefficient.
        ward_amplitude = sp.simplify(sp.factorial(8) * sp.conjugate(onepoint))
        need(ward_amplitude != 0, "folded Ward octic amplitude vanished")

        descendant_norm_squared = sp.simplify(sp.factorial(8) * sp.rf(2 * h, 8))
        # This is also the eight-step sl2 recursion product
        # product_{k=0}^7 (k+1)(2h+k).
        recurrence_norm = sp.prod((k + 1) * (2 * h + k) for k in range(8))
        need(descendant_norm_squared == recurrence_norm, "Virasoro descendant norm recursion failed")
        projection_fraction = sp.factor(sp.simplify(sp.expand_complex(
            ward_amplitude * sp.conjugate(ward_amplitude) / descendant_norm_squared
        )))
        need(
            projection_fraction == expected_projection_fraction[alpha],
            "normalized octic projection fraction changed",
        )

        # L_1^d F8 is the directional derivative in the supplied shift.
        # Its x0^7 coefficient already proves that the shifted state is not
        # an L^d-primary.  No large Fock-space projection is needed.
        directional = sp.expand(sum(b[index] * sp.diff(f8, x[index]) for index in range(4)))
        witness = sp.Poly(directional, *x).coeff_monomial(x[0] ** 7)
        need(sp.simplify(witness - b[0] / (12 * sp.sqrt(140))) == 0 and witness != 0,
             "shifted octic unexpectedly remained primary")

        # Actual source relation: Lhat_-1=Lcal^alpha_-8/8.  Hence the same
        # matrix element for eight original source modes carries 8^8.
        source_lcal_amplitude = sp.simplify(8**8 * ward_amplitude)
        need(source_lcal_amplitude != 0, "original Lcal_-8 amplitude vanished")

        cases[str(alpha)] = {
            "delta": [str(value) for value in delta],
            "d": [str(value) for value in d],
            "d_squared": str(2 * h),
            "h": str(h),
            "b_Uplus_T_O_T_d": [str(value) for value in b],
            "b_norm_squared": str((b.H * b)[0]),
            "I8_of_b": str(value_i8),
            "F8hat_onepoint": str(onepoint),
            "bra_ward_amplitude_8factorial_conjugate_F8hat_b": str(ward_amplitude),
            "descendant_norm_squared": str(descendant_norm_squared),
            "normalized_projection_fraction": str(projection_fraction),
            "L1d_nonprimary_witness_x0_power7": str(witness),
            "bra_Lcal_minus8_power8_amplitude": str(source_lcal_amplitude),
            "mobius_full_state_level8_amplitude": (
                "(1-|w|^2)^(85/256) * w^8 * conjugate(F8hat(b))"
            ),
        }

    # Bounded root-OPE result retained from the preceding premise.  The
    # physical roots are alpha=r/sqrt(2), and U+ is the same original
    # polarization used above.  This tests the actual coefficient rather
    # than a fitted embedding.
    physical_roots = [sp.Matrix(root) / sp.sqrt(2) for root in lattice["roots"]]
    root_weights = [clean(positive.T * root) for root in physical_roots]
    need(len(root_weights) == 240, "actual root count changed")
    need(sp.partition(8) == 22, "weight eight no longer has 22 oscillator partitions")
    moment_polynomials = {}
    for degree in (2, 4, 6, 8):
        moment_polynomials[degree] = sp.expand(sum(
            sum(weight[index] * x[index] for index in range(4)) ** degree
            for weight in root_weights
        ))
    need(all(moment_polynomials[degree] == 0 for degree in (2, 4, 6)),
         "a lower pure-holomorphic root moment survived")
    need(sp.expand(moment_polynomials[8] - sp.Rational(15, 2) * i8) == 0,
         "actual root octic restriction coefficient changed")
    root_ope_overlap = sp.simplify(
        sp.Rational(15, 2) * 96 * sp.sqrt(140) / sp.factorial(8)
    )
    need(root_ope_overlap == sp.sqrt(140) / 56, "root-OPE F8 overlap changed")

    integer_roots = [tuple(map(int, root)) for root in lattice["roots"]]
    pair_inner_products = Counter(
        sum(left[index] * right[index] for index in range(8)) // 2
        for left in integer_roots
        for right in integer_roots
    )
    need(pair_inner_products == Counter({0: 30240, 1: 13440, -1: 13440, 2: 240, -2: 240}),
         "E8 ordered root-pair census changed")
    full_root_ope_norm_squared = sp.simplify(sum(
        count * sp.binomial(inner + 7, 8)
        for inner, count in pair_inner_products.items()
    ))
    need(full_root_ope_norm_squared == 15600, "full weight-eight root-OPE norm changed")

    report = {
        "status": "PASS_EXACT_ORIGINAL_FOLDED_OCTIC_OVERLAP",
        "claim": "MATH.NATIVE.FOLDED-OCTIC-WARD.01",
        "direct_original_car_current_dictionary": {
            "A": [[str(value) for value in source_to_positive.row(row)] for row in range(4)],
            "A_times_A_adjoint": "I4",
            "A_times_A_transpose": "zero4",
            "field": ":I8(A*(j_sea(z)+delta/z))/(96*sqrt(140)):",
            "normal_ordering": (
                "retained bosonic current normal order, equivalently the original "
                "point-split CAR limit with its derivative counterterms"
            ),
            "charge": 0,
            "relative_fold_energy": 8,
            "original_source_excitation_energy": 64,
            "root_or_spinor_vertex_extension_required": False,
            "support": "the same eight original preimage arcs",
        },
        "complete_translation_word_check": {
            "definition": "y_n=d.h_-n; L_-1^d acts as y1+sum_n n*y_(n+1)*partial_y_n",
            "levels": complete_words,
            "zero_shift_unshifted_vacuum_control": "PASS",
        },
        "cases": cases,
        "root_ope_bounded_result": {
            "phase_correct_bilocal": (
                "sum_beta epsilon(beta,-beta)^-1 Y(e^beta,z)e^-beta"
            ),
            "weight8_state": (
                "S8=sum_beta [z^8] exp(sum_n>=1 beta(-n) z^n/n) Omega"
            ),
            "all_22_oscillator_partitions_retained": True,
            "pure_holomorphic_root_moments_2_4_6": "zero",
            "pure_holomorphic_root_moment_8": "(15/2) I8",
            "F8hat_overlap": str(root_ope_overlap),
            "full_S8_norm_squared": str(full_root_ope_norm_squared),
            "exclusivity": "not inferred from the overlap",
        },
        "scope": (
            "Exact supplied E8_1 shifted-current dictionary and original O/M/U+ fold. "
            "The result proves a nonzero source matrix element and a full Mobius-orbit "
            "amplitude. It does not identify a passive dictionary transport with an "
            "active physical preparation. A profile is an intervention input; the P1 "
            "net/state/reflection/time embedding is not inferred from this matrix element."
        ),
    }
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
