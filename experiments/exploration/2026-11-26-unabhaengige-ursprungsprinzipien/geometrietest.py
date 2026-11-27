#!/usr/bin/env python3
"""Scharfer Strukturtest fuer zwei unabhaengige Geometrie-Urspruenge.

Der Haupttest beginnt mit einer glatten komplexen Flaeche als kleinstem
vollstaendigen Schnitt zweier Quadriken. Er setzt weder einen TFPT-Randkern
noch vier Marken, D5, 5, 3+2, 16 oder E8 voraus.

Der zweite Test prueft nur ein moegliches physikalisches Woerterbuch fuer
isolierte 4d-N=2-Hypersurflaechensingularitaeten. Sein Fehlschlag ist kein
Ausschluss von Singularitaetsgeometrie im Allgemeinen.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations, product
from math import factorial, prod

import sympy as sp


def d5_structure() -> dict[str, object]:
    """Berechne Wurzelzahl, Cartanmatrix und Halbspinorgewichte von D5."""

    roots: set[tuple[int, ...]] = set()
    for i, j in combinations(range(5), 2):
        for si, sj in product((-1, 1), repeat=2):
            root = [0] * 5
            root[i], root[j] = si, sj
            roots.add(tuple(root))

    simple = (
        (1, -1, 0, 0, 0),
        (0, 1, -1, 0, 0),
        (0, 0, 1, -1, 0),
        (0, 0, 0, 1, -1),
        (0, 0, 0, 1, 1),
    )
    cartan = sp.Matrix([[sum(x * y for x, y in zip(a, b)) for b in simple] for a in simple])

    # Doppelte Halbspinorgewichte: (±1,...,±1) mit gerader Minuszahl.
    half_spin = {
        signs
        for signs in product((-1, 1), repeat=5)
        if sum(sign < 0 for sign in signs) % 2 == 0
    }
    s5_orbit_profile: dict[int, int] = {}
    for weight in half_spin:
        minus_count = sum(sign < 0 for sign in weight)
        s5_orbit_profile[minus_count] = s5_orbit_profile.get(minus_count, 0) + 1

    assert len(roots) == 40
    assert {sum(x * x for x in root) for root in roots} == {2}
    assert cartan.det() == 4
    assert len(half_spin) == 16
    assert s5_orbit_profile == {0: 1, 2: 10, 4: 5}
    return {
        "roots": roots,
        "cartan": cartan,
        "half_spin": half_spin,
        "s5_orbit_profile": s5_orbit_profile,
    }


def del_pezzo_pencil() -> dict[str, object]:
    """Pruefe den Fuenferdiskriminanten eines glatten Quadrikenbueschels."""

    u, v = sp.symbols("u v")
    lambdas = (0, 1, 2, 3, 4)
    pencil = sp.diag(*(u + lam * v for lam in lambdas))
    discriminant = sp.factor(pencil.det())
    expected = sp.prod(u + lam * v for lam in lambdas)

    ranks = []
    for lam in lambdas:
        # Der Punkt [u:v]=[-lambda:1] ist der zugehoerige singulaere Bueschelpunkt.
        ranks.append(pencil.subs({u: -lam, v: 1}).rank())

    assert sp.expand(discriminant - expected) == 0
    assert sp.Poly(discriminant, u, v).total_degree() == 5
    assert ranks == [4] * 5
    assert len(set(lambdas)) == 5  # Das diagonale Modell hat glatten Basisschnitt.

    return {
        "discriminant": discriminant,
        "singular_members": len(lambdas),
        "ranks": ranks,
    }


def inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Zaehle positive, negative und verschwindende Eigenwerte exakt."""

    positive = negative = zero = 0
    for eigenvalue, multiplicity in matrix.eigenvals().items():
        if eigenvalue.is_positive:
            positive += multiplicity
        elif eigenvalue.is_negative:
            negative += multiplicity
        elif eigenvalue.is_zero:
            zero += multiplicity
        else:
            raise AssertionError(f"Inertia fuer {eigenvalue} nicht exakt entscheidbar")
    return positive, negative, zero


def pencil_is_smooth_over_c(
    metric: sp.Matrix,
    form: sp.Matrix,
    roots: tuple[sp.Expr, ...],
) -> bool:
    """Pruefe: einfache Rang-4-Quadrik und ihr Vertex liegt nicht im Basisschnitt."""

    for root in roots:
        singular_quadric = form - root * metric
        if singular_quadric.rank() != 4:
            return False
        kernel = singular_quadric.nullspace()
        if len(kernel) != 1:
            return False
        vertex = kernel[0]
        if sp.simplify((vertex.T * metric * vertex)[0]) == 0:
            return False
    return True


def real_witt_pencils() -> dict[str, object]:
    """Zwei Signatur-(4,1)-Bueschel: einmal 3+2, einmal 5+0."""

    lam = sp.symbols("lambda")
    metric = sp.diag(1, 1, 1, 1, -1)

    # Auf der letzten positiven und der negativen Richtung entsteht der
    # J-selbstadjungierte Block [[0,1],[-1,0]] mit Eigenwerten ±i.
    paired = sp.zeros(5)
    paired[0, 0], paired[1, 1], paired[2, 2] = 1, 2, 3
    paired[3, 4] = paired[4, 3] = 1
    paired_operator = metric * paired
    paired_polynomial = sp.factor((paired - lam * metric).det())
    paired_roots = (sp.Integer(1), sp.Integer(2), sp.Integer(3), sp.I, -sp.I)

    # Dieselbe Signatur beider Formen, aber nun sind alle fuenf Marken reell.
    real = sp.diag(1, 2, 3, 4, -5)
    real_operator = metric * real
    real_polynomial = sp.factor((real - lam * metric).det())
    real_roots = tuple(sp.Integer(value) for value in range(1, 6))

    assert inertia(metric) == (4, 1, 0)
    assert inertia(paired) == (4, 1, 0)
    assert inertia(real) == (4, 1, 0)
    assert paired_operator.T * metric == metric * paired_operator
    assert real_operator.T * metric == metric * real_operator
    assert paired_polynomial == (lam - 1) * (lam - 2) * (lam - 3) * (lam**2 + 1)
    assert real_polynomial == sp.prod(lam - value for value in range(1, 6))
    assert pencil_is_smooth_over_c(metric, paired, paired_roots)
    assert pencil_is_smooth_over_c(metric, real, real_roots)

    neutral_plane = metric.extract((3, 4), (3, 4))
    complex_block = paired_operator.extract((3, 4), (3, 4))
    assert inertia(neutral_plane) == (1, 1, 0)
    assert complex_block**2 == -sp.eye(2)

    # Die drei reellen Eigenlinien plus die reelle Ebene des konjugierten
    # Paars sind eine operatorinvariante 3+2-Zerlegung von R^5.
    carrier_projector = sp.diag(1, 1, 1, 0, 0)
    assert carrier_projector * paired_operator == paired_operator * carrier_projector

    return {
        "metric_inertia": inertia(metric),
        "paired_form_inertia": inertia(paired),
        "real_form_inertia": inertia(real),
        "paired_polynomial": paired_polynomial,
        "real_polynomial": real_polynomial,
        "paired_smooth": True,
        "real_smooth": True,
        "neutral_plane_inertia": inertia(neutral_plane),
        "carrier_block_dimensions": (3, 2),
        "complex_block_square": complex_block**2,
    }


def split_obstruction() -> dict[str, object]:
    """Pruefe, ob die volle S5-Wirkung eine 3+2-Teilung auszeichnet."""

    points = tuple(range(5))
    two_subsets = tuple(frozenset(pair) for pair in combinations(points, 2))
    group = tuple(permutations(points))

    def image(subset: frozenset[int], perm: tuple[int, ...]) -> frozenset[int]:
        return frozenset(perm[i] for i in subset)

    globally_fixed = [subset for subset in two_subsets if all(image(subset, g) == subset for g in group)]
    chosen = frozenset((0, 1))
    stabilizer = [g for g in group if image(chosen, g) == chosen]
    orbit = {image(chosen, g) for g in group}

    assert len(group) == factorial(5) == 120
    assert len(two_subsets) == 10
    assert globally_fixed == []
    assert len(stabilizer) == factorial(2) * factorial(3) == 12
    assert len(orbit) == 10

    # Negativkontrolle: Eine zusaetzlich faktorisierte binaere Quintik traegt 3+2.
    u, v = sp.symbols("u v")
    cubic = u**3 - 2 * v**3
    quadratic = u**2 - 3 * v**2
    controlled = sp.expand(cubic * quadratic)
    assert sp.Poly(cubic, u, v).total_degree() == 3
    assert sp.Poly(quadratic, u, v).total_degree() == 2
    assert sp.discriminant(controlled.subs(v, 1), u) != 0

    return {
        "partitions": len(two_subsets),
        "fixed_partitions": len(globally_fixed),
        "stabilizer_size": len(stabilizer),
        "orbit_size": len(orbit),
        "control_factor_degrees": (3, 2),
        "control_discriminant": sp.factor(controlled),
    }


def brieskorn_pham_spectrum(exponents: tuple[int, int, int, int]) -> dict[str, object]:
    """Exaktes 4d-N=2-Deformationsspektrum einer diagonalen IHS in C^4."""

    weights = tuple(Fraction(1, a) for a in exponents)
    delta = sum(weights, Fraction(0)) - 1
    if delta <= 0:
        raise ValueError("Die kanonische Dreiform hat keine positive Skalendimension.")

    dimensions: list[Fraction] = []
    for alpha in product(*(range(a - 1) for a in exponents)):
        monomial_weight = sum((Fraction(n, a) for n, a in zip(alpha, exponents)), Fraction(0))
        dimensions.append((1 - monomial_weight) / delta)

    dimensions.sort(reverse=True)
    milnor = prod(a - 1 for a in exponents)
    coulomb_rank = sum(value > 1 for value in dimensions)
    flavor_rank = sum(value == 1 for value in dimensions)

    assert len(dimensions) == milnor
    assert milnor == 2 * coulomb_rank + flavor_rank
    assert all(a + b == 2 for a, b in zip(dimensions, reversed(dimensions)))
    return {
        "exponents": exponents,
        "milnor": milnor,
        "coulomb_rank": coulomb_rank,
        "flavor_rank": flavor_rank,
        "dimensions": dimensions,
    }


def n2_dictionary_test() -> dict[str, object]:
    """Teste die ausdrueckliche Lesart 5=Coulombrang und 3=Flavourrang."""

    e8_suspension = brieskorn_pham_spectrum((2, 2, 3, 5))
    required_milnor = 2 * 5 + 3

    # Fuer diagonale Brieskorn-Pham-Formen erzwingt mu=13 wegen der Primzahl 13
    # bis auf Permutation genau (2,2,2,14). Der Census macht das explizit.
    mu13 = []
    for exponents in combinations_with_replacement(range(2, 15), 4):
        if prod(a - 1 for a in exponents) == required_milnor:
            mu13.append(brieskorn_pham_spectrum(exponents))

    assert e8_suspension["milnor"] == 8
    assert (e8_suspension["coulomb_rank"], e8_suspension["flavor_rank"]) == (4, 0)
    assert [entry["exponents"] for entry in mu13] == [(2, 2, 2, 14)]
    assert [(entry["coulomb_rank"], entry["flavor_rank"]) for entry in mu13] == [(6, 1)]
    assert not any(
        (entry["coulomb_rank"], entry["flavor_rank"]) == (5, 3)
        for entry in mu13
    )

    return {
        "required_milnor": required_milnor,
        "e8_suspension": e8_suspension,
        "diagonal_mu13": mu13,
    }


def fmt_spectrum(values: list[Fraction]) -> str:
    return ", ".join(str(value) for value in values)


def main() -> None:
    pencil = del_pezzo_pencil()
    d5 = d5_structure()
    real_pencils = real_witt_pencils()
    split = split_obstruction()
    n2 = n2_dictionary_test()

    print("DEL-PEZZO-TEST: BESTANDEN")
    print(f"  Diskriminante: {pencil['discriminant']}")
    print(f"  singulaere Quadriken: {pencil['singular_members']}, Raenge: {pencil['ranks']}")
    print(f"  D5-Wurzeln: {len(d5['roots'])}, det(Cartan): {d5['cartan'].det()}")
    print(f"  Halbspinorgewichte / Geradenzaehler: {len(d5['half_spin'])}")
    print(f"  S5-Orbits der Halbspinorgewichte: {d5['s5_orbit_profile']}")
    print("REELLER WITT-TEST: BESTANDEN")
    print(f"  Signaturen J, A_Paar, A_reell: {real_pencils['metric_inertia']}, "
          f"{real_pencils['paired_form_inertia']}, {real_pencils['real_form_inertia']}")
    print(f"  mit Paar: {real_pencils['paired_polynomial']}")
    print(f"  ohne Paar: {real_pencils['real_polynomial']}")
    print(f"  beide komplex glatt: {real_pencils['paired_smooth'] and real_pencils['real_smooth']}")
    print(f"  neutrale Paarebene: {real_pencils['neutral_plane_inertia']}, "
          f"Traegerbloecke: {real_pencils['carrier_block_dimensions']}")
    print("UNMARKIERTER 3+2-MONODROMIETEST: NEGATIVKONTROLLE")
    print(f"  3+2-Partitionen: {split['partitions']}, S5-fixiert: {split['fixed_partitions']}")
    print(f"  Wahl einer Partition: Stabilisator {split['stabilizer_size']}, Orbit {split['orbit_size']}")
    print(f"  Negativkontrolle mit gesetzter Faktorisierung: {split['control_discriminant']}")
    print("4D-N=2-WOERTERBUCH: NEGATIV")
    e8 = n2["e8_suspension"]
    print(f"  E8-Suspension: mu={e8['milnor']}, (r,f)=({e8['coulomb_rank']},{e8['flavor_rank']})")
    print(f"  Spektrum: {fmt_spectrum(e8['dimensions'])}")
    for control in n2["diagonal_mu13"]:
        print(
            "  mu=13-Kontrolle "
            f"{control['exponents']}: (r,f)=({control['coulomb_rank']},{control['flavor_rank']}), "
            f"Spektrum={fmt_spectrum(control['dimensions'])}"
        )
    print("GESAMTURTEIL: Reelles Grad-4-del-Pezzo-Bueschel ist Kandidat Stufe 1 fuer 5, 3+2, D5 und 16.")
    print("OFFEN: Signatur (4,1) erlaubt das Paar, erzwingt es aber nicht; die physische Eichwirkung fehlt.")


if __name__ == "__main__":
    main()
