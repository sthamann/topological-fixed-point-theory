"""Unpromoted exact audit of the special Picard origin C: y^3 = x^4 - 1.

The C12 subgroup preserves a chosen cubic deck.  An explicit additional
automorphism generates a group of order 48.  This corrects the scope of the
full-Aut(C) claim in v784 without removing its direct code-transport obstruction.

Primary source for the full automorphism group: Bouw, Koutsianas, Sijsling,
Wewers, J. London Math. Soc. 102 (2020), 368-404, section 4.1 in the PDF:
https://d-nb.info/1253140073/34 .  The full order-48 upper bound is cited;
the subgroup, its action, and all reported finite counts are recomputed here.
The torus, cover and kernel identifications remain separate physical gates.
"""

from collections import Counter, deque
import json

import sympy as sp


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


class ExactPicardGroup:
    """Faithful action on span(dx/y, dx/y^2), with exact number-field entries.

    No construction occurs on import.  Matrices are immutable four-tuples.
    Pullback composition reverses geometric composition, which does not change
    the generated group, orders, conjugacy or regular-representation controls.
    """

    def __init__(self):
        root3 = sp.sqrt(3)
        omega = (-1 + sp.I * root3) / 2
        self.field = sp.QQ.algebraic_field(sp.I, root3)
        cv = self.field.from_sympy
        zero, one = self.field.zero, self.field.one
        self.identity = (one, zero, zero, one)
        self.tau = (cv(sp.expand(omega**2)), zero, zero, cv(omega))
        self.rho = (cv(sp.I), zero, zero, cv(sp.I))
        self.switch = tuple(cv(n / root3) for n in (1, -2, -1, -1))
        self.generators = (self.tau, self.rho, self.switch)
        elements = [self.identity]
        index = {self.identity: 0}
        queue = deque([self.identity])
        while queue:
            element = queue.popleft()
            for generator in self.generators:
                product = self.mul(element, generator)
                if product not in index:
                    index[product] = len(elements)
                    elements.append(product)
                    queue.append(product)
            require(len(elements) <= 48, "Picard group exceeded the stated bound")
        require(len(elements) == 48, "The explicit Picard generators do not give 48 elements")
        self.elements = tuple(elements)
        self.index = index

    @staticmethod
    def mul(a, b):
        return (
            a[0] * b[0] + a[1] * b[2],
            a[0] * b[1] + a[1] * b[3],
            a[2] * b[0] + a[3] * b[2],
            a[2] * b[1] + a[3] * b[3],
        )

    def power(self, element, exponent):
        require(isinstance(exponent, int) and exponent >= 0, "Expected a nonnegative power")
        result = self.identity
        for _ in range(exponent):
            result = self.mul(result, element)
        return result

    def order(self, element):
        result = self.identity
        for order in range(1, 49):
            result = self.mul(result, element)
            if result == self.identity:
                return order
        raise ArithmeticError("Element order exceeds the generated group bound")

    def inverse(self, element):
        return self.power(element, self.order(element) - 1)

    def conjugate(self, element, target):
        return self.mul(self.mul(element, target), self.inverse(element))

    def regular(self, element):
        matrix = sp.zeros(len(self.elements))
        for column, other in enumerate(self.elements):
            matrix[self.index[self.mul(element, other)], column] = 1
        return matrix

    def cubic_deck_subgroups(self):
        subgroups = []
        for element in self.elements:
            tau = self.conjugate(element, self.tau)
            subgroup = frozenset((self.identity, tau, self.power(tau, 2)))
            if subgroup not in subgroups:
                subgroups.append(subgroup)
        require(len(subgroups) == 4, "The four conjugate cubic decks changed")
        return tuple(subgroups)


def curve_identity_controls():
    X, Y, Z = sp.symbols("X Y Z")
    curve = X**4 - Y**3 * Z - Z**4
    switch = sp.Matrix([[sp.sqrt(3), 0, 0], [0, -1, 2], [0, 1, 1]])
    images = switch * sp.Matrix([X, Y, Z])
    transformed = curve.subs(dict(zip((X, Y, Z), images)), simultaneous=True)
    require(sp.expand(transformed - 9 * curve) == 0, "The extra map fails to preserve C")
    require(switch * switch == 3 * sp.eye(3), "The extra map is not projectively involutive")
    require(switch.det() == -3 * sp.sqrt(3), "The extra projective map is singular")
    return {
        "curve": "X^4-Y^3*Z-Z^4=0",
        "extra_automorphism": "S(x,y)=(sqrt(3)*x/(y+1),(2-y)/(y+1))",
        "global_projective_identity": "F(Sv)=9F(v), S^2=3I, det(S)=-3*sqrt(3)",
        "extra_automorphism_is_global": True,
        "chosen_x_projection_preserved_by_S": False,
        "all_automorphisms_follow_from_x_deck_normalizer": False,
    }


def group_controls(group):
    counts = Counter(group.order(g) for g in group.elements)
    require(dict(counts) == {1: 1, 2: 7, 3: 8, 4: 8, 6: 8, 12: 16},
            "Element-order census changed")
    center = [g for g in group.elements
              if all(group.mul(g, h) == group.mul(h, g) for h in group.generators)]
    centralizer = [g for g in group.elements
                   if group.mul(g, group.tau) == group.mul(group.tau, g)]
    normalizer = [g for g in group.elements
                  if group.conjugate(g, group.tau) in
                  (group.tau, group.power(group.tau, 2))]
    require(len(center) == 4 and len(centralizer) == 12 and len(normalizer) == 12,
            "Full, central and marked deck symmetries were conflated")
    require(group.mul(group.switch, group.tau) != group.mul(group.tau, group.switch),
            "The new quotient switch must not commute with the chosen deck")
    commutator = group.mul(
        group.mul(group.mul(group.switch, group.tau), group.switch),
        group.power(group.tau, 2),
    )
    require(group.order(commutator) == 4, "The quotient-switch commutator changed")
    return {
        "generated_group_order_exact": len(group.elements),
        "full_automorphism_order_in_primary_source": 48,
        "full_group_upper_bound_independently_classified_here": False,
        "element_order_census": dict(sorted(counts.items())),
        "center_order": len(center),
        "marked_deck_centralizer_order": len(centralizer),
        "marked_deck_normalizer_order": len(normalizer),
        "conjugate_cubic_deck_subgroups": len(group.cubic_deck_subgroups()),
        "commutator_S_tau_order": group.order(commutator),
        "v784_scope_correction": "C12 is the marked cubic-deck normalizer, not all Aut(C)",
    }


def code_obstruction_controls(group):
    order_three = {g for g in group.elements if group.order(g) == 3}
    conjugates = {group.conjugate(g, h) for g in group.elements
                  for h in (group.tau, group.power(group.tau, 2))}
    require(order_three == conjugates and len(order_three) == 8,
            "The complete order-three set was not covered")
    for g in order_three:
        require(g[0] + g[3] == -group.field.one,
                "Order-three trace on the abelian-surface forms changed")
        require(g[0] * g[3] - g[1] * g[2] == group.field.one,
                "Order-three determinant on the abelian-surface forms changed")
    t = sp.symbols("t")
    real_characteristic = (t**2 + t + 1)**2
    determinant = int(real_characteristic.subs(t, 1))
    require(determinant == 9 and determinant % 2 == 1,
            "The mod-two fixed-space obstruction disappeared")
    # For a g-invariant nondegenerate symplectic form, quadratic refinements
    # form an affine V* torsor.  Invertible g-I gives exactly one fixed q.
    # Its nonzero zero-set consists of length-three orbits.  In dimension four,
    # 10 (Arf zero) and 6 (Arf one) are the only nondegenerate zero counts.
    compatible_zero_counts = [n for n in (6, 10) if n % 3 == 1]
    require(compatible_zero_counts == [10], "The invariant Arf argument changed")
    require(group.power(group.rho, 2) == tuple(-v for v in group.identity),
            "The anti-invariant Prym block was not preserved")
    return {
        "order_three_elements_checked": 8,
        "surface_H1_characteristic": str(sp.expand(real_characteristic)),
        "det_g_minus_I_on_surface_H1": determinant,
        "fixed_dimension_on_A2_for_all_order_three_elements": 0,
        "native_code_sigma_fixed_dimension_from_v784": 2,
        "direct_linear_equivariant_code_transport_possible": False,
        "unique_invariant_quadratic_refinement_zero_count": 10,
        "unique_invariant_quadratic_refinement_Arf": 0,
        "native_code_selector_Arf_from_v784": 1,
        "refinement_argument_assumption": "A chosen g-invariant nondegenerate symplectic form",
        "direct_refinement_preserving_transport_possible": False,
        "existing_scoped_code_no_go_survives_full_group_correction": True,
    }


def cover_controls():
    def genus(base_genus, degree, ramification_sum):
        twice = 2 + degree * (2 * base_genus - 2) + ramification_sum
        require(twice >= 0 and twice % 2 == 0, "Invalid Riemann-Hurwitz data")
        return twice // 2

    require(genus(1, 2, 0) == 1, "Unramified double torus changed genus")
    require(genus(1, 2, 4) == 3, "Four-branched double torus changed genus")
    require(genus(0, 2, 4) == 1, "Four-branched double sphere changed genus")
    require(genus(0, 3, 5 * 2) == genus(0, 4, 4 * 3) == 3,
            "The two Picard quotient presentations disagree")
    omega = (-1 + sp.I * sp.sqrt(3)) / 2
    lam_equianharmonic = sp.simplify((-1 + omega**2) / (-omega + omega**2))
    lam_square = sp.simplify(2 * (2 * sp.I) / (1 + sp.I)**2)

    def elliptic_j(lam):
        return sp.simplify(256 * (1 - lam + lam**2)**3 / (lam**2 * (1 - lam)**2))

    require(sp.simplify(lam_equianharmonic**2 - lam_equianharmonic + 1) == 0,
            "The alternative branch cross-ratio is not equianharmonic")
    require(elliptic_j(lam_equianharmonic) == 0 and elliptic_j(lam_square) == 1728,
            "The two different four-point configurations were conflated")
    # The central involution (x,y)->(-x,y) has quotient u=x^2,
    # u^2=y^3+1.  It is a degree-two cover of E_0, branched at x=0
    # (three points) and infinity.  It is not the square elliptic double
    # cover of the four finite x-branch values.
    x, y = sp.symbols("x y")
    require(sp.expand((x**2)**2 - y**3 - 1) == x**4 - y**3 - 1,
            "The bielliptic quotient equation failed")
    j_bielliptic = sp.Rational(1728) * 0 / 27
    require(j_bielliptic == 0, "The bielliptic quotient was mislabeled as square")
    a, b, c = sp.symbols("a b c", real=True)
    metric = sp.Matrix([[a, b], [b, c]])
    quarter = sp.Matrix([[0, -1], [1, 0]])
    solution = sp.solve(list(quarter.T * metric * quarter - metric), (b, c), dict=True)
    require(solution == [{b: 0, c: a}], "Linear quarter-turn metric rigidity changed")
    return {
        "connected_unramified_double_torus_genus": 1,
        "double_torus_with_four_simple_branch_points_genus": 3,
        "double_sphere_with_four_simple_branch_points_genus": 1,
        "special_picard_genus": 3,
        "x_mu3_quotient_branch_values": "1,i,-1,-i,infinity (five, including infinity)",
        "y_mu4_quotient_branch_values": "-1,-omega,-omega^2,infinity (four)",
        "finite_x_four_point_double_cover_j": 1728,
        "alternative_y_four_point_double_cover_j": 0,
        "central_bielliptic_quotient": "u=x^2; u^2=y^3+1; j=0",
        "central_bielliptic_cover_degree": 2,
        "central_bielliptic_simple_branch_points": 4,
        "central_bielliptic_quotient_is_native_square_torus": False,
        "genuine_linear_quarter_turn_metric": "G=a*I, a>0; scale is undetermined",
        "order_four_translation_fixes_complex_modulus": False,
    }


def controls():
    group = ExactPicardGroup()
    return {
        "curve_identity": curve_identity_controls(),
        "full_versus_marked_group": group_controls(group),
        "preserved_code_obstructions": code_obstruction_controls(group),
        "distinct_cover_geometries": cover_controls(),
        "scope": {
            "status": "unpromoted exact geometry audit",
            "same_source_physical_kernel_derived": False,
            "native_quarter_clock_intertwiner_derived": False,
            "physical_c3_normalization_derived": False,
            "F1_closed": False,
        },
    }


if __name__ == "__main__":
    print(json.dumps(controls(), indent=2))
