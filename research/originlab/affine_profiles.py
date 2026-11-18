"""Check level-one profiled root returns on the common reflected Gram form."""

from itertools import combinations_with_replacement

import sympy as sp

from language import positive
from reconstruction import quotient, require, scalar
from return_moments import rows


def check_returns(gram, neutral, charged, products, neutral_transfer=None):
    """Read K(a,b)=E(a)*E(b) on neutral grade zero before outer selection.

    Profiles are self-adjoint. Index zero is the unit profile. The declared
    products give real coefficients in the same ordered list of profiles.
    """
    count, width = len(charged), len(neutral)
    require(width and count >= 2 and all(len(indices) == width for indices in charged),
            "supply neutral words and matching lists for the charged returns")
    indices = [index for group in charged for index in group]
    require(all(type(i) is int and 0 <= i < gram.rows for i in (*neutral, *indices)),
            "profile return index is outside the original word basis")
    require(all(len(indices) == len(set(indices)) for indices in (neutral, *charged)),
            "repeated word index within one profile")
    require(not any(gram.extract(neutral, indices)), "neutral and charged words must be orthogonal")
    pairs = list(combinations_with_replacement(range(count), 2))
    require(len(products) == len(pairs) and all(len(row) == count for row in products),
            "supply coefficients for every product in lexicographic upper-triangle order")
    products = [list(map(scalar, row)) for row in products]
    require(all(value.is_Rational for row in products for value in row) and
            sp.Matrix(products[:count]) == sp.eye(count), "invalid real profile products or unit")
    g = gram.extract(neutral, neutral)
    forms = {(a,b): gram.extract(charged[a], charged[b]) for a,b in pairs}
    if neutral_transfer is not None:
        require(positive(neutral_transfer), "neutral transfer form must be Hermitian PSD")
        forms["transfer"] = neutral_transfer
    _, metric, ops = quotient(g, forms)
    unit = sp.eye(metric.rows)
    if ops[0,0] != unit:
        return {"status":"UNIT_RETURN_MISMATCH", "neutral_rank":metric.rows}

    central = [ops[0,b] for b in range(count)]
    adjoints = [metric.inv()*operator.H*metric for operator in central]
    ward = []
    closure = []
    for (a,b), product in zip(pairs, products):
        expected = sum((value*operator for value,operator in zip(product, central)), sp.zeros(metric.rows))
        if ops[a,b] != expected or (a == 0 and metric*central[b] != central[b].H*metric):
            ward.append([a,b])
        if ops[a,b] != adjoints[a]*central[b]:
            closure.append([a,b])

    outside = (gram.extract(indices, indices) -
               gram.extract(indices, charged[0])*forms[0,0].pinv()*gram.extract(charged[0], indices))
    require(positive(outside), "negative charged Gram residual outside the unit return")
    status = "COMPATIBLE_ON_SUPPLIED_PROFILES"
    if ward:
        status = "PROFILE_WARD_MISMATCH"
    elif closure:
        status = "ADDITIONAL_NEUTRAL_RETURNS_REQUIRED"
    result = {
        "status": status,
        "neutral_rank": metric.rows,
        "charged_rank_outside_unit_return": outside.rank(),
        "profile_product_mismatches": ward,
        "profile_return_mismatches": closure,
        "central_action": [rows(operator) for operator in central],
        "scope": "supplied ground and charged grade-one forms, common level-one root Ward relation",
        "physical_source_identity_proved": False,
    }
    if neutral_transfer is not None and not ward and not closure:
        transfer = ops["transfer"]
        residuals = [metric*(transfer*operator-operator*transfer) for operator in central]
        mixing = [index for index,residual in enumerate(residuals) if residual != sp.zeros(metric.rows)]
        result["neutral_transfer"] = {
            "status":"PROFILE_MIXING_DETECTED" if mixing else "COMMUTES_ON_SUPPLIED_FORMS",
            "noncommuting_profile_indices":mixing,
            "commutator_forms":[rows(residual) for residual in residuals],
        }
    return result


def controls():
    # Two evaluations of the profile s with s^2=1. Both ground directions
    # are retained, with a nonorthonormal input metric.
    e1,e2,e3,e4,e5 = (sp.eye(5)[:,j] for j in range(5))
    unit = [e3, 2*e4]
    words = sp.Matrix.hstack(e1, 2*e2, *unit, e3, -2*e4)
    products = [[1,0], [0,1], [1,0]]
    g = words.H*words
    direct = check_returns(g, [0,1], [[2,3],[4,5]], products, sp.diag(1,2)/2)
    mixed = check_returns(g, [0,1], [[2,3],[4,5]], products, sp.Matrix([[1,1],[1,2]])/2)

    # The same scalar mixed return and s^2 Ward moment as a symmetric
    # average, with the missing variance moved to an extra charged state.
    words = sp.Matrix.hstack(e1, e3, e5)
    extra = check_returns(words.H*words, [0], [[1],[2]], products)
    require(direct["status"] == "COMPATIBLE_ON_SUPPLIED_PROFILES" and
            direct["central_action"][1] == [["1","0"],["0","-1"]] and
            direct["charged_rank_outside_unit_return"] == 0 and
            direct["neutral_transfer"]["noncommuting_profile_indices"] == [] and
            mixed["neutral_transfer"]["noncommuting_profile_indices"] == [1] and
            extra["status"] == "ADDITIONAL_NEUTRAL_RETURNS_REQUIRED" and
            extra["profile_return_mismatches"] == [[1,1]] and
            extra["charged_rank_outside_unit_return"] == 1,
            "profile closure or its negative control did not replay")
    return {"two_evaluation_sectors": direct, "transfer_mixing_evaluation_sectors": mixed["neutral_transfer"],
            "scalar_average_with_extra_charged_path": extra}
