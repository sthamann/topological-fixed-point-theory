"""Exact response bounds for four marks under tetrahedral symmetry."""

import sympy as sp

from reconstruction import matrix, require


def response_bounds(moments):
    """Third moment interval and attaining spectra on the closed interval [0, 1]."""
    mass, first, second = moments[:3]
    if not (0 <= second <= first <= mass and first**2 <= mass*second):
        return {"status": "INFEASIBLE_MOMENTS"}

    if first == 0 or first == mass:
        spectrum = [(sp.S.One if first else sp.S.Zero, mass)] if mass else []
        spectra = [spectrum, spectrum]
    else:
        above_zero = second/first
        below_one = (first-second)/(mass-first)
        lower_weight = first/above_zero
        upper_weight = (mass-first)/(1-below_one)
        spectra = [[(sp.S.Zero, mass-lower_weight), (above_zero, lower_weight)],
                   [(below_one, upper_weight), (sp.S.One, mass-upper_weight)]]

    bounds = [sum(weight*value**3 for value, weight in spectrum) for spectrum in spectra]
    if len(moments) == 4:
        status = "OBSERVED_THIRD_COMPATIBLE" if bounds[0] <= moments[3] <= bounds[1] else "REFUTED_THIRD"
    else:
        status = "THIRD_FIXED" if bounds[0] == bounds[1] else "THIRD_SEPARATES"
    return {
        "status": status,
        "third_interval": list(map(str, bounds)),
        "spectral_witnesses": [
            [{"transfer_eigenvalue": str(value), "weight": str(weight)}
             for value, weight in spectrum if weight]
            for spectrum in spectra
        ],
        "extra_copies_required": int(bool(first**2 < mass*second)),
    }


def analyze(data):
    """Check exact moment forms in the original order of a declared four-mark frame."""
    require(isinstance(data, dict) and set(data) <= {"moments", "conserved_distinct_weights"},
            "unknown four-mark transfer inputs")
    moments = data.get("moments")
    require(isinstance(moments, list) and len(moments) in (3, 4),
            "supply moment forms at steps zero, one, two and optionally three")
    conserved = data.get("conserved_distinct_weights", False)
    require(type(conserved) is bool, "conserved_distinct_weights must be Boolean")
    scope = {"input_role": "declared_finite_oracle", "physical_source_derived": False,
             "transfer_spectrum_domain": "[0, 1]", "conserved_distinct_weights": conserved}
    diagonals, edges = [], []
    for step, rows in enumerate(moments):
        form = matrix(rows)
        require(form.shape == (4, 4), "each moment must use the four-mark basis")
        diagonal, edge = form[0, 0], form[0, 1]
        if form != form.H or form != (diagonal-edge)*sp.eye(4)+edge*sp.ones(4):
            return {**scope, "status": "REFUTED_TETRAHEDRAL_FORM", "first_defect_step": step}
        if conserved and edge:
            return {**scope, "status": "REFUTED_CONSERVED_CARTAN_FORM", "first_defect_step": step,
                    "mixed_entry": str(edge)}
        diagonals.append(diagonal)
        edges.append(edge)

    if conserved:
        sectors = {"four_weight_lines": (4, diagonals)}
    else:
        sectors = {"singlet": (1, [d+3*e for d, e in zip(diagonals, edges)]),
                   "triplet": (3, [d-e for d, e in zip(diagonals, edges)])}
    reports = {name: {"multiplicity": multiplicity, **response_bounds(values)}
               for name, (multiplicity, values) in sectors.items()}
    refuted = any(report["status"] in ("REFUTED_THIRD", "INFEASIBLE_MOMENTS") for report in reports.values())
    result = {**scope, "status": "REFUTED_TRANSFER_MOMENTS" if refuted else "EXACT_FOUR_MARK_TRANSFER_BOUNDS",
              "sectors": reports}
    if not refuted:
        result["minimum_additional_directions_from_first_three"] = sum(
            report["multiplicity"]*report["extra_copies_required"] for report in reports.values())
        result["next_probe"] = next(({"step": 3, "sector": name, "interval": report["third_interval"]}
                                     for name, report in reports.items() if report["status"] == "THIRD_SEPARATES"), None)
    return result


def control():
    """A stationary singlet and a triplet with two distinct compatible spectra."""
    singlet = sp.ones(4)/4
    triplet = sp.eye(4)-singlet
    forms = [singlet + value*triplet for value in (1, sp.Rational(1, 2), sp.Rational(5, 16))]
    report = analyze({"moments": [[list(map(str, row)) for row in form.tolist()] for form in forms]})
    return {**report, "input_role": "finite_control"}
