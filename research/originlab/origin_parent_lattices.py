"""Exact negative controls for identifying a larger parent from torus data.

E8 + E8 and D16+ have the same genus-one theta series but different
root/current algebras. This module checks finite root data and a finite
theta expansion. It does not infer an all-order theorem from that finite
agreement, identify the original TFPT source, or solve F1.

Full theta/CFT statement: Gaberdiel and Volpato, JHEP 06 (2009) 048,
https://arxiv.org/abs/0903.4107 . Their physical genus <= 4 vacuum
amplitudes also agree. The graph trace(A**4) here is NOT such an amplitude.
"""

from collections import Counter
from itertools import combinations, product


def require(condition, why):
    if not condition:
        raise ArithmeticError(why)


def dn_roots_doubled(n):
    """Return D_n roots in twice Euclidean coordinates (norm squared 8)."""
    require(isinstance(n, int) and n >= 4, "D_n root rank must be >= 4")
    roots = []
    for i, j in combinations(range(n), 2):
        for a, b in product((-2, 2), repeat=2):
            root = [0] * n
            root[i], root[j] = a, b
            roots.append(tuple(root))
    return roots


def e8_roots_doubled():
    roots = dn_roots_doubled(8)
    roots += [signs for signs in product((-1, 1), repeat=8)
              if sum(value < 0 for value in signs) % 2 == 0]
    require(len(roots) == len(set(roots)) == 240,
            "E8 root enumeration did not produce 240 distinct roots")
    return roots


def root_graph_audit(roots):
    """Edges mean root inner product +1 or -1; opposite roots are excluded."""
    count = len(roots)
    require(count > 0 and len(set(roots)) == count,
            "root graph needs nonempty distinct roots")
    require(all(len(root) == len(roots[0]) for root in roots),
            "root coordinates have inconsistent ranks")
    adjacency = [set() for _ in roots]
    histograms = [Counter() for _ in roots]
    for i in range(count):
        for j in range(i, count):
            dot4 = sum(a * b for a, b in zip(roots[i], roots[j]))
            require(dot4 % 4 == 0, "root inner product is not integral")
            dot = dot4 // 4
            histograms[i][dot] += 1
            if i != j:
                histograms[j][dot] += 1
            if abs(dot) == 1:
                adjacency[i].add(j)
                adjacency[j].add(i)
    require(all(hist == histograms[0] for hist in histograms),
            "root pair distributions differ within this root system")
    remaining, components = set(range(count)), []
    while remaining:
        seed = next(iter(remaining))
        remaining.remove(seed)
        todo, size = [seed], 0
        while todo:
            i = todo.pop()
            size += 1
            newly_reached = adjacency[i] & remaining
            remaining.difference_update(newly_reached)
            todo.extend(newly_reached)
        components.append(size)
    trace_a3 = sum(len(adjacency[i] & adjacency[j])
                   for i in range(count) for j in adjacency[i])
    trace_a4 = sum(len(a & b) ** 2 for a in adjacency for b in adjacency)
    require(all(sum(value * value for value in root) == 8 for root in roots),
            "root norm squared differs from two")
    return {
        "root_count": count,
        "norm_squared": 2,
        "root_inner_products_per_root": {
            str(dot): value for dot, value in sorted(histograms[0].items())
        },
        "graph_degree": sorted({len(neighbors) for neighbors in adjacency}),
        "graph_component_sizes": sorted(components),
        "trace_adjacency_cubed": trace_a3,
        "trace_adjacency_fourth": trace_a4,
        "root_pair_dot_products_checked": count * (count + 1) // 2,
    }


def polynomial_multiply(a, b, cutoff):
    out = [0] * (cutoff + 1)
    for i, coefficient in enumerate(a[:cutoff + 1]):
        if not coefficient:
            continue
        for j, other in enumerate(b[:cutoff + 1 - i]):
            if other:
                out[i + j] += coefficient * other
    return out


def polynomial_power(a, n, cutoff):
    out = [1] + [0] * cutoff
    for _ in range(n):
        out = polynomial_multiply(out, a, cutoff)
    return out


def theta_dn_plus(n, q_order=6):
    """Finite theta(D_n+) expansion via theta3, theta4, theta2; t=q^(1/8)."""
    require(isinstance(n, int) and n > 0 and n % 8 == 0,
            "even unimodular D_n+ control requires rank divisible by eight")
    require(isinstance(q_order, int) and q_order >= 0,
            "theta cutoff must be a nonnegative integer")
    cutoff = 8 * q_order
    theta3 = [0] * (cutoff + 1)
    theta4 = [0] * (cutoff + 1)
    theta2 = [0] * (cutoff + 1)
    for m in range(-cutoff - 1, cutoff + 2):
        exp_integer, exp_half = 4 * m * m, (2 * m + 1) ** 2
        if exp_integer <= cutoff:
            theta3[exp_integer] += 1
            theta4[exp_integer] += 1 if m % 2 == 0 else -1
        if exp_half <= cutoff:
            theta2[exp_half] += 1
    terms = [polynomial_power(poly, n, cutoff)
             for poly in (theta3, theta4, theta2)]
    coefficients = [sum(poly[i] for poly in terms)
                    for i in range(cutoff + 1)]
    require(all(value % 2 == 0 for value in coefficients),
            "D_n+ theta formula produced a half-integral coefficient")
    coefficients = [value // 2 for value in coefficients]
    require(all(value == 0 for i, value in enumerate(coefficients) if i % 8),
            "D_n+ theta formula retained a nonintegral q power")
    return [coefficients[8 * i] for i in range(q_order + 1)]


def controls():
    e8 = e8_roots_doubled()
    e8_pair = [root + (0,) * 8 for root in e8]
    e8_pair += [(0,) * 8 + root for root in e8]
    d16 = dn_roots_doubled(16)
    lhs, rhs = root_graph_audit(e8_pair), root_graph_audit(d16)
    require(lhs["root_count"] == rhs["root_count"] == 480,
            "rank-sixteen parents must both have 480 roots")
    expected_histogram = {"-2": 1, "-1": 56, "0": 366, "1": 56, "2": 1}
    require(lhs["root_inner_products_per_root"] ==
            rhs["root_inner_products_per_root"] == expected_histogram,
            "larger-parent pair-distribution negative control changed")
    require(lhs["graph_degree"] == rhs["graph_degree"] == [112],
            "larger-parent root-product graph degree changed")
    require(lhs["trace_adjacency_cubed"] ==
            rhs["trace_adjacency_cubed"] == 3010560,
            "larger-parent three-step degeneracy changed")
    require(lhs["graph_component_sizes"] == [240, 240]
            and rhs["graph_component_sizes"] == [480],
            "root-product graph no longer distinguishes the two parents")
    require(lhs["trace_adjacency_fourth"] == 319979520
            and rhs["trace_adjacency_fourth"] == 237404160,
            "larger-parent fourth graph moment discriminator changed")
    order = 6
    theta_e8 = theta_dn_plus(8, order)
    theta_pair = polynomial_multiply(theta_e8, theta_e8, order)
    theta_d16 = theta_dn_plus(16, order)
    expected_theta = [1, 480, 61920, 1050240, 7926240, 37500480, 135480960]
    require(theta_pair == theta_d16 == expected_theta,
            "finite genus-one theta degeneracy changed")
    return {
        "status": "EXACT_FINITE_PARENT_SPECTRUM_NONUNIQUENESS_CONTROL",
        "E8_plus_E8": lhs,
        "D16_plus_root_system": rhs,
        "q_convention": "q exponent = norm squared / 2",
        "checked_theta_order": order,
        "checked_theta_coefficients": order + 1,
        "theta_E8_plus_E8": theta_pair,
        "theta_D16_plus": theta_d16,
        "finite_agreement_proves_all_order_theta_equality": False,
        "all_order_theta_equality_uses_external_modular_form_theorem": True,
        "full_theta_theorem_reference": "https://arxiv.org/abs/0903.4107",
        "graph_trace_four_is_genus_four_vacuum_amplitude": False,
        "spectra_and_pair_counts_identify_parent_operator_algebra": False,
        "operator_product_graph_distinguishes_this_pair": True,
        "original_P1_source_operator_identified": False,
        "F1_solved": False,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(controls(), sort_keys=True, indent=2))
