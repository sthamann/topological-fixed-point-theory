"""Exact bounded searches on the harmonic seam and its transfer alternatives."""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations_with_replacement, permutations
from math import gcd

from reconstruction import require, scalar


PERMUTATIONS = tuple(permutations(range(4)))


def projective_stabilizer(numerator, denominator):
    """Permutations preserving the cross-ratio of (infinity, 0, 1, lambda)."""
    points = ((1, 0), (0, 1), (1, 1), (numerator, denominator))
    brackets = tuple(tuple(x*w-y*z for z, w in points) for x, y in points)
    return tuple((a, b, c, d) for a, b, c, d in PERMUTATIONS
                 if denominator*brackets[a][c]*brackets[b][d] ==
                 numerator*brackets[a][d]*brackets[b][c])


def is_four_cycle(perm):
    first = perm[0]
    return first != 0 and perm[first] != 0 and perm[perm[first]] != 0


def rational_observation(data, name):
    if name not in data:
        return None
    value = scalar(data[name])
    require(value.is_Rational, f"{name} must be rational")
    return Fraction(int(value.p), int(value.q))


def modulus_sweep(bound, ratio):
    digest = sha256()
    counts = {"V4": 0, "D4": 0}
    quarter_survivors = []
    ratio_survivors = []
    for q in range(1, bound+1):
        for p in range(1, q+1):
            if gcd(p, q) != 1:
                continue
            stabilizer = projective_stabilizer(-p*p, q*q)
            order = len(stabilizer)
            require(order == (8 if p == q else 4), "unexpected projective stabilizer")
            counts["D4" if order == 8 else "V4"] += 1
            digest.update(f"{p}/{q} {order}\n".encode())
            if any(is_four_cycle(perm) for perm in stabilizer):
                quarter_survivors.append(str(Fraction(p, q)))
            if ratio is not None and ratio*p == q:
                ratio_survivors.append(str(Fraction(p, q)))
    return {"evaluated": sum(counts.values()), "projective_stabilizers": counts,
            "quarter_survivors_t": quarter_survivors, "stream_sha256": digest.hexdigest(),
            "odd_green_ratio_observed": str(ratio) if ratio is not None else None,
            "odd_green_ratio_survivors_t": ratio_survivors if ratio is not None else None,
            "continuous_t_if_odd_green_ratio_bound": str(1/ratio) if ratio is not None else None}


def transfer_sweep(grid, adjacent, opposite):
    """Sweep (alternating, doublet) eigenvalues for the stationary singlet model."""
    digest = sha256()
    matched = 0
    witnesses = []
    for alternating in range(grid+1):
        for doublet in range(grid+1):
            # A common denominator 4*grid keeps the two orbital tests exact.
            edge_num = grid-alternating
            opposite_num = grid+alternating-2*doublet
            digest.update(f"{alternating} {doublet} {edge_num} {opposite_num}\n".encode())
            if adjacent is not None and edge_num != 4*grid*adjacent:
                continue
            if opposite is not None and opposite_num != 4*grid*opposite:
                continue
            matched += 1
            if len(witnesses) < 2:
                witnesses.append([str(Fraction(v, grid)) for v in (alternating, doublet)])
    return {"evaluated": (grid+1)**2, "accidentally_a4_symmetric": grid+1,
            "matches_supplied_orbit_observations": matched, "matching_example_eigenvalues": witnesses,
            "stream_sha256": digest.hexdigest(), "eigenvalue_order": ["alternating", "doublet"]}


def bridge_observations(rows):
    require(isinstance(rows, dict) and rows, "supply named rational modulus/ratio observations")
    observations = {}
    for t, ratio in rows.items():
        t = rational_observation({"t": t}, "t")
        ratio = rational_observation({"ratio": ratio}, "ratio")
        require(0 < t < 1 and 0 < ratio <= 1 and t not in observations,
                "moduli must be distinct off-square points with ratios in (0, 1]")
        observations[t] = ratio
    return observations


def reciprocal_bridge_sweep(data):
    allowed = {"factor_limit", "max_factors", "max_power", "fit", "withheld"}
    require(isinstance(data, dict) and set(data) <= allowed, "unknown reciprocal bridge input")
    limit = data.get("factor_limit", 64)
    factors = data.get("max_factors", 3)
    power = data.get("max_power", 3)
    require(type(limit) is int and 2 <= limit <= 128, "factor_limit must be in 2..128")
    require(type(factors) is int and 0 <= factors <= 4, "max_factors must be in 0..4")
    require(type(power) is int and 0 <= power <= 8, "max_power must be in 0..8")
    fit = bridge_observations(data["fit"])
    withheld = bridge_observations(data["withheld"])
    require(not (fit.keys() & withheld.keys()), "fit and withheld moduli must be disjoint")
    parameters = tuple(sorted({Fraction(j) for j in range(2, limit+1)} | {Fraction(34, 7)}))
    moments = {t: {a: (1+a*t)/(a+t) for a in parameters} for t in (*fit, *withheld)}

    def response(t, exponent, word):
        value = t**exponent
        for a in word:
            value *= moments[t][a]
        return value

    evaluated = fitted = retained = 0
    examples = []
    digest = sha256()
    for length in range(factors+1):
        for word in combinations_with_replacement(parameters, length):
            for exponent in range(power+1):
                evaluated += 1
                if not all(response(t, exponent, word) == value for t, value in fit.items()):
                    continue
                fitted += 1
                prediction = {str(t): str(response(t, exponent, word)) for t in withheld}
                agrees = all(prediction[str(t)] == str(value) for t, value in withheld.items())
                retained += int(agrees)
                model = {"power": exponent, "factors": list(map(str, word)),
                         "withheld_prediction": prediction, "passes_withheld": agrees}
                digest.update(f"{model}\n".encode())
                if len(examples) < 8:
                    examples.append(model)
    return {"scope": "declared reciprocal rational transfer-ratio family",
            "evaluated": evaluated, "fit_survivors": fitted, "withheld_survivors": retained,
            "degree_bound": power+factors, "distinct_points": len(fit)+len(withheld),
            "fit_survivor_stream_sha256": digest.hexdigest(), "fit_survivor_examples": examples}


def run(data):
    allowed = {"max_denominator", "transfer_grid", "odd_green_ratio", "adjacent", "opposite", "rational_bridge"}
    require(isinstance(data, dict) and set(data) <= allowed, "unknown harmonic seam search inputs")
    bound, grid = data.get("max_denominator", 256), data.get("transfer_grid", 128)
    require(type(bound) is int and 1 <= bound <= 4096, "max_denominator must be in 1..4096")
    require(type(grid) is int and 1 <= grid <= 4096, "transfer_grid must be in 1..4096")
    ratio = rational_observation(data, "odd_green_ratio")
    require(ratio is None or ratio >= 1, "odd_green_ratio must be at least one in this chamber")
    adjacent = rational_observation(data, "adjacent")
    opposite = rational_observation(data, "opposite")
    moduli = modulus_sweep(bound, ratio)
    transfers = transfer_sweep(grid, adjacent, opposite)
    result = {"status": "EXACT_BOUNDED_HARMONIC_SEAM_SEARCH", "physical_source_derived": False,
              "max_denominator": bound, "transfer_grid": grid,
              "moduli": moduli, "d4_transfers": transfers,
              "total_evaluated": moduli["evaluated"]+transfers["evaluated"],
              "missing_physical_bridge": "same-source Ward and odd Green forms, independently transported off square"}
    if "rational_bridge" in data:
        bridge = reciprocal_bridge_sweep(data["rational_bridge"])
        result["rational_bridge"] = bridge
        result["total_evaluated"] += bridge["evaluated"]
    return result
