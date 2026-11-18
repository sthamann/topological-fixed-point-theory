#!/usr/bin/env python3
"""Exact charge, cocycle, sector-lift and quarter-time checks for W_ext."""

from __future__ import annotations

from fractions import Fraction as Q
from itertools import product
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
I = 1j


def need(value, label):
    if not bool(value):
        raise RuntimeError(label)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def cyclic(vector):
    return (vector[-1],) + tuple(vector[:-1])


def affine(vector):
    out = list(cyclic(vector))
    out[0] += 1
    return tuple(out)


def eps_z(left, right):
    return -1 if sum(left[i] * right[j] for i in range(8) for j in range(i)) % 2 else 1


def eta_r(vector):
    return -1 if vector[7] * sum(vector[:7]) % 2 else 1


ONES = (1,) * 8
ZERO = (0,) * 8
E0 = (1,) + (0,) * 7


def normalize(base, half):
    """Unique n+c*s section, c in {0,1}, with 2s=(1,...,1)."""
    if half < 2:
        return tuple(base), half
    return add(base, ONES), half - 2


def add_labels(left, right):
    return normalize(add(left[0], right[0]), left[1] + right[1])


def eps_ext(left, right):
    """Associative section cocycle for E_(n+c*s)=E_n X^c.

    Relations: X E_n X^-1=i^Q(n) E_n and X^2=-E_(1,...,1).
    It restricts to the ordered CAR cocycle on Z8 and to the previously
    fixed even E8 cocycle on D8 union (D8+s).
    """
    n, c = left
    m, d = right
    value = (I ** (c * sum(m))) * eps_z(n, m)
    if c and d:
        value *= -eps_z(add(n, m), ONES)
    return complex(round(value.real), round(value.imag))


def prior_eps_e8(left, right):
    """The fixed even E8 cocycle from SEAM-JOINT-WARD-STATE, section 20.

    This formula is defined on D8 union (D8+s), hence both integer bases
    have even coordinate sum in the common n+c*s section.
    """
    n, c = left
    m, d = right
    need(sum(n) % 2 == 0 and sum(m) % 2 == 0,
         "prior E8 cocycle domain")
    sign_exponent = c * (sum(m) // 2) + c * d
    value = ((-1) ** sign_exponent) * eps_z(n, m)
    if c and d:
        value *= eps_z(add(n, m), ONES)
    return value


def sector(label):
    n, c = label
    return (("0", "v") if c == 0 else ("s", "c"))[sum(n) % 2]


def r_label(label):
    n, c = label
    return cyclic(n), c


def t_label(label):
    n, c = label
    return affine(n), c


def eta_ext(label, kappa=I):
    n, c = label
    return eta_r(n) * (kappa ** c)


def tau_gauge(label):
    """Intertwiner-basis gauge: identity on NS and the E8 vacuum VOA."""
    n, c = label
    return (-1) ** (c * sum(n))


def x_action(label):
    """Left multiplication by the marked spinor intertwiner X=E_s."""
    return eps_ext((ZERO, 1), label), add_labels((ZERO, 1), label)


def inverse_t_label(label):
    n, c = label
    shifted = list(n)
    shifted[0] -= 1
    return tuple(shifted[1:] + shifted[:1]), c


def w_action(label, kappa):
    return eta_ext(label, kappa), t_label(label)


def w_inverse_action(label, kappa):
    source = inverse_t_label(label)
    return eta_ext(source, kappa).conjugate(), source


def x_inverse_action(label, action=x_action):
    n, c = label
    source = (n, 0) if c == 1 else (tuple(a - 1 for a in n), 1)
    phase, target = action(source)
    need(target == label, "inverse marked spinor intertwiner")
    return phase.conjugate(), source


def mixed_word(label, kappa, action=x_action):
    """Apply W X W^* X^* from right to left."""
    phase = 1
    current = label
    operators = (
        lambda y: x_inverse_action(y, action),
        lambda y: w_inverse_action(y, kappa),
        action,
        lambda y: w_action(y, kappa),
    )
    for operator in operators:
        factor, current = operator(current)
        phase *= factor
    return phase, current


def main():
    # The old ordered CAR lift is exact, including its relative phases.
    bits = tuple(product((0, 1), repeat=8))
    for n in bits:
        for m in bits:
            need(eta_r(add(n, m)) * eps_z(n, m)
                 == eta_r(n) * eta_r(m) * eps_z(cyclic(n), cyclic(m)),
                 "ordered CAR cyclic-lift cocycle")

    # The extension is associative and restricts to both required sublattices.
    generators = [(ZERO, 0)] + [
        (tuple(1 if i == j else 0 for i in range(8)), 0) for j in range(8)
    ] + [(ZERO, 1), (ONES, 0)]
    samples = list(generators)
    samples += [add_labels(x, y) for x in generators for y in generators]
    samples = list(dict.fromkeys(samples))
    for x in samples:
        for y in samples:
            xy = add_labels(x, y)
            for z in generators:
                need(eps_ext(x, y) * eps_ext(xy, z)
                     == eps_ext(y, z) * eps_ext(x, add_labels(y, z)),
                     "extended cocycle associativity")
    for n in bits:
        for m in bits:
            need(eps_ext((n, 0), (m, 0)) == eps_z(n, m),
                 "exact ordered CAR restriction")
            if sum(n) % 2 == 0 and sum(m) % 2 == 0:
                for c in (0, 1):
                    for d in (0, 1):
                        need(eps_ext((n, c), (m, d))
                             == prior_eps_e8((n, c), (m, d)),
                             "exact fixed E8 cocycle restriction")

    # Both kappa choices are the forced square roots on the Ramond modules.
    for kappa in (I, -I):
        need(abs(kappa * kappa + 1) == 0, "Ramond lift kappa squared")
        for x in samples:
            for y in generators:
                need(eta_ext(add_labels(x, y), kappa) * eps_ext(x, y)
                     == eta_ext(x, kappa) * eta_ext(y, kappa)
                     * eps_ext(r_label(x), r_label(y)),
                     "extended cyclic lift")

    reps = {
        "0": (ZERO, 0),
        "v": (E0, 0),
        "s": (ZERO, 1),
        "c": (E0, 1),
    }
    sector_map = {name: sector(t_label(label)) for name, label in reps.items()}
    need(sector_map == {"0": "v", "v": "0", "s": "c", "c": "s"},
         "four-sector affine action")

    # The two kappa representatives are related by a basis gauge which is
    # pointwise the identity on both supplied local subtheories 0+v and 0+s.
    for label in samples:
        n, c = label
        if c == 0 or sum(n) % 2 == 0:
            need(tau_gauge(label) == 1,
                 "tau is identity on NS and E8 vacuum subtheories")
        for kappa in (I, -I):
            conjugated = (tau_gauge(label) * eta_ext(label, kappa)
                           * tau_gauge(t_label(label)))
            need(conjugated == eta_ext(label, -kappa),
                 "tau conjugates W_kappa to W_minus_kappa")

        x_phase, x_target = x_action(label)
        transported_x_phase = (tau_gauge(label) * x_phase
                               * tau_gauge(x_target))
        expected_ratio = -1 if sector(label) in ("v", "c") else 1
        need(transported_x_phase == expected_ratio * x_phase,
             "tau changes X only on the v-to-c intertwiner channel")

    # A simultaneous basis transport leaves every mixed word invariant.  If
    # X is instead held fixed as an extra marking, W X W^* X^*=-i*kappa I
    # distinguishes the two representatives relative to that marking.
    for kappa in (I, -I):
        def transported_x(label):
            phase, target = x_action(label)
            return (tau_gauge(label) * phase * tau_gauge(target), target)

        for label in samples:
            fixed_phase, fixed_target = mixed_word(label, kappa)
            moved_phase, moved_target = mixed_word(label, -kappa, transported_x)
            need(fixed_target == label and fixed_phase == -I * kappa,
                 "fixed-X mixed word")
            need((moved_phase, moved_target) == (fixed_phase, fixed_target),
                 "fully transported mixed word is gauge invariant")

    # W_ext^8 is the genuine integer carry E_q on every sector, including
    # its ordered cocycle phase. It is never collapsed to a finite register.
    for kappa in (I, -I):
        for label in samples:
            start = label
            current = label
            phase = 1
            for _ in range(8):
                phase *= eta_ext(current, kappa)
                current = t_label(current)
            expected_label = add_labels(start, (ONES, 0))
            expected_phase = eps_ext((ONES, 0), start)
            need(current == expected_label and phase == expected_phase,
                 "eighth power is exact integer charge translation")

    # The two lifts of the same E8 g-character differ only on the twisted
    # module. The distinction is observable once W1 is fixed.
    def g_cartan(label):
        return (-1) ** label[1]                    # NS+, R-

    def g_standard(label):
        n, c = label
        return (-1) ** (sum(n) + c)               # 0+,v-,s-,c+

    cartan_signs = {name: g_cartan(label) for name, label in reps.items()}
    standard_signs = {name: g_standard(label) for name, label in reps.items()}
    need(cartan_signs == {"0": 1, "v": 1, "s": -1, "c": -1},
         "Cartan lift signs")
    need(standard_signs == {"0": 1, "v": -1, "s": -1, "c": 1},
         "standard spin-sector lift signs")
    need(all(g_cartan(t_label(x)) == g_cartan(x) for x in reps.values()),
         "Cartan lift commutes with W1")
    need(all(g_standard(t_label(x)) == -g_standard(x) for x in reps.values()),
         "standard spin-sector lift anticommutes with W1")
    need(all(g_cartan(x) == 1 for name, x in reps.items() if name in ("0", "v")),
         "retained original Cartan holonomy is identity on NS")

    # Charge and quarter-time covariance are polynomial identities, hence
    # they hold for integer and half-integer charges alike.
    p = sp.Matrix(sp.symbols("p0:8", real=True))
    R = sp.zeros(8)
    for k in range(8):
        R[k, (k - 1) % 8] = 1
    e0 = sp.eye(8)[:, 0]
    one = sp.ones(8, 1)
    quarter_states = {}
    pstar = sp.Matrix([sp.Rational(1, 2)] * 4 + [-sp.Rational(1, 2)] * 4)
    for alpha in (sp.Rational(1, 4), sp.Rational(3, 4)):
        delta = sp.Matrix([(k + alpha) / 8 - sp.Rational(1, 2) for k in range(8)])
        K = lambda x: sp.expand(x.dot(x) / 2 + delta.dot(x))
        Tp = R * p + e0
        need(sp.expand(K(Tp) - K(p) - (one.dot(p) + alpha) / 8) == 0,
             "quarter Ward identity for arbitrary charge")
        need(K(pstar) == 0 and K(R * pstar + e0) == alpha / 8,
             "half-charge quarter-state transport")
        need(K(sp.zeros(8, 1)) == 0 and K(e0) == alpha / 8,
             "vacuum quarter-state transport")
        quarter_states[str(alpha)] = {
            "K_vacuum": "0", "K_T_vacuum": str(alpha / 8),
            "K_pstar": "0", "K_T_pstar": str(alpha / 8),
        }

    # Conjugate the entire construction by the already fixed M. This is the
    # actual D8_F frame, not the canonical compiler D8 audit frame.
    fold = json.loads((ROOT / "experiments/theory-contracts/"
        "v16-operator-pullback-20261014/fold-results.json").read_text())
    M = sp.Matrix([[sp.Rational(x) for x in row] for row in fold["explicit_Weyl_matrix"]])
    RF = M * R * M.T
    hF = M * e0
    need(M.T * M == sp.eye(8) and RF.T * RF == sp.eye(8), "fixed M-frame orthogonality")
    x = sp.Matrix(sp.symbols("x0:8", real=True))
    for alpha in (sp.Rational(1, 4), sp.Rational(3, 4)):
        delta = sp.Matrix([(k + alpha) / 8 - sp.Rational(1, 2) for k in range(8)])
        d = M * delta
        KF = lambda y: sp.expand(y.dot(y) / 2 + d.dot(y))
        QF = (M * one).dot(x)
        need(sp.expand(KF(RF * x + hF) - KF(x) - (QF + alpha) / 8) == 0,
             "actual M-frame Ward identity")

    result = {
        "verdict": "FULL_FOUR_MODULE_AFFINE_DEFECT_UNIQUE_UP_TO_INTERTWINER_BASIS",
        "sector_order": ["0", "v", "s", "c"],
        "T1_sector_map": sector_map,
        "E8_untwisted_to_g_twist": {"0+s": "v+c"},
        "fermion_NS_preserved": "0+v",
        "Ramond_sum_preserved": "s+c",
        "ordered_CAR_phase": "eta_R(n)=(-1)^(n_7 sum_(j<7)n_j)",
        "extended_cocycle": {
            "relations": ["X E_n X^-1=i^Q(n) E_n", "X^2=-E_(1,...,1)"],
            "restriction_to_Z8": "exact ordered CAR cocycle",
            "restriction_to_E8": "existing even E8 cocycle",
            "ordinary_cross_locality": False,
        },
        "Ramond_R8_lift_kappa_choices": ["i", "-i"],
        "kappa_relation": "intertwiner-basis gauge U_tau; not an independent physical choice",
        "tau_gauge": "U_tau|n+c*s>=(-1)^(c*Q(n)); identity on 0+v and 0+s",
        "fixed_X_word": "W_kappa X W_kappa^* X^*=-i*kappa I",
        "fully_transported_mixed_words": "unchanged under U_tau",
        "W_ext_eighth_power": "E_(1,...,1) with its cocycle phase",
        "Cartan_lift_signs": cartan_signs,
        "standard_spin_sector_lift_signs": standard_signs,
        "Cartan_lift_W1_relation": "commutes; retained NS Cartan four-word is +I",
        "standard_lift_W1_relation": "anticommutes; corresponding four-word is -I",
        "twisted_module_lift_phase_fixed_by_E8_character": False,
        "quarter_states": quarter_states,
        "source_Ward": "W_ext^* Q W_ext=Q+1; W_ext^* (8K_delta) W_ext=8K_delta+Q+alpha",
        "actual_target_frame": {
            "D8": "D8_F=M(D8_CAR)",
            "linear_part": "R_F=M R8 M^-1",
            "carry": [str(x) for x in hF],
            "canonical_compiler_D8_identified": False,
        },
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
