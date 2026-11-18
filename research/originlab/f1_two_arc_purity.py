"""Search over source rules on the four seam marks, read as purity of two arcs.

The four marks cut the seam circle into four arcs; A is the union of two
opposite arcs. A rule assigns non-negative integer multiplicities to the three
pairings of the marks and yields the ratio

    X(x) = Tr rho_A^2 / (Tr rho_A1^2 Tr rho_A2^2),   x = cross ratio of the marks.

The three recorded pairing weights (session of 18 November 2026, rounds 17 to
25) contribute 1/(1-x), 1-x and x^2/(1-x). The guards are exact. Nothing here
is an original producer: that the P1 form is this purity is the hypothesis
H-R2, and `original_r113_values` stays None.
"""
from itertools import product

import sympy as sp

x = sp.symbols("x", positive=True)
# opposite pairing (alternating signs), the two side pairings, each times |P_12|
PAIRINGS = {"opposite": 1/(1-x), "side_far": 1-x, "side_near": x**2/(1-x)}
OCTETS = 1        # recorded weight 1/2 at each mark = c/16 with c = 8 (round 17)
CURRENTS = 248


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def ratio(rule):
    total = sum(n*PAIRINGS[name]**OCTETS for n, name in zip(rule, PAIRINGS))
    return sp.simplify(total/(rule[0]+rule[1]))


def guards(rule):
    """Exact conditions a purity of a pure seam state has to meet."""
    n_opp, n_far, n_near = rule
    out = {}
    # distant arcs decouple in both channels (the arcs, and their complement)
    out["distant_arcs_decouple"] = n_far == n_near and n_opp+n_far > 0
    if not out["distant_arcs_decouple"]:
        return out
    value = ratio(rule)
    out["pure_state"] = sp.simplify(value/x-(value/x).subs(x, 1-x)) == 0
    series = sp.series(value, x, 0, 3).removeO()
    first, second = series.coeff(x, 1), series.coeff(x, 2)
    out["fermion_pairs_exchanged"] = sp.nsimplify(16*first)
    out["integral_exchange"] = (16*first).is_integer and first >= 0
    # integer spin (recorded axiom A2): start and end marks are not distinguished
    out["integer_spin"] = n_opp == n_far == n_near
    out["weight_two_coefficient"] = second
    return out


def search(limit=4):
    """Enumerate multiplicities up to `limit` and keep primitive survivors."""
    rows = []
    for rule in product(range(limit+1), repeat=3):
        if not any(rule) or sp.gcd(sp.gcd(rule[0], rule[1]), rule[2]) != 1:
            continue
        g = guards(rule)
        if g.get("distant_arcs_decouple") and g["pure_state"] and g["integral_exchange"]:
            rows.append({"rule": rule, "ratio": str(ratio(rule)),
                         "fermion_pairs_exchanged": str(g["fermion_pairs_exchanged"]),
                         "integer_spin": g["integer_spin"],
                         "square": str(ratio(rule).subs(x, sp.Rational(1, 2))),
                         "counterwitness": str(ratio(rule).subs(x, sp.Rational(1, 5)))})
    return rows


def witnesses():
    """Exact facts of the session the rule has to reproduce."""
    turn = sp.Matrix([[0, 1], [-1, 0]])
    require(turn**2 == -sp.eye(2) and set(turn.eigenvals()) == {sp.I, -sp.I},
            "two fermion copies pass a mark by a quarter turn")
    twists = sorted(sp.Rational(1, 2)+k/2 for k in (sp.Rational(1, 2), -sp.Rational(1, 2)))
    require(twists == [sp.Rational(1, 4), sp.Rational(3, 4)], "recorded quarter twists")
    for alpha in twists:
        for charge in range(-5, 6):
            energy = (sum(j+alpha for j in range(charge)) if charge >= 0
                      else sum(j-alpha for j in range(1, -charge+1)))
            require(energy == sp.Rational(charge*charge, 2)+(alpha-sp.Rational(1, 2))*charge,
                    "recorded charge function F_alpha")
    e8 = ratio((1, 1, 1))
    free = ratio((1, 0, 0))
    require(sp.simplify(e8-(1-x+x**2)/(1-x)) == 0 and sp.simplify(free-1/(1-x)) == 0, "closed forms")
    require(sp.Rational(CURRENTS, 256)+sp.Rational(1, 2)**2*2/16 == 1
            and sp.series(e8, x, 0, 3).removeO() == 1+x**2, "currents and stress at weight two")
    j = 256*(1-x+x**2)**3/(x**2*(1-x)**2)
    require(j.subs(x, sp.Rational(1, 2)) == 1728, "square torus")
    return {"quarter_turn_order": 4, "twists": [str(t) for t in twists],
            "inverse_temperature_over_pi": 4,
            "e8_ratio": str(e8), "free_ratio": str(free),
            "square": {"e8": str(e8.subs(x, sp.Rational(1, 2))), "free": str(free.subs(x, sp.Rational(1, 2)))},
            "counterwitness": {"e8": str(e8.subs(x, sp.Rational(1, 5))), "free": str(free.subs(x, sp.Rational(1, 5)))}}


def reflections():
    """The quarter turn at a mark as product of two reflections, and the mark types."""
    swap, sign = sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
    turn = sign*swap
    require(turn == sp.Matrix([[0, 1], [-1, 0]]) and swap*sign*swap*sign == -sp.eye(2),
            "quarter turn is sign times exchange")

    def weight(matrix, fermions=8):
        return fermions*sum(mult*sp.nsimplify(sp.arg(value)/(2*sp.pi))**2/2
                            for value, mult in matrix.eigenvals().items())
    weights = {"1": weight(sp.eye(2)), "S": weight(swap), "P": weight(sign), "PS": weight(turn)}
    require(weights == {"1": 0, "S": 1, "P": 1, "PS": sp.Rational(1, 2)}, "weights of the mark types")
    return {
        "weights": {name: str(value) for name, value in weights.items()},
        "fock_relation": "S P S P = (-1)^N on two copies; anticommuting only at odd total parity",
        # verdicts of the session script check_marks_u_y_dictionary.py (round 26)
        "U_is_S_and_Y_is_P_on_equal_copies": "EXCLUDED_FOR_GROUP_A1: Z_U = Z_Y and mixed answer >= 0",
        "U_is_S_and_Y_is_P_on_one_shared_line": "FORM_OF_GROUP_A9: Z_U = cos(phi), Z_Y = 0, Z_UY = i sin(phi)",
        # round 27, check_shared_line_phase.py: coherent_origin.marked_operators read back
        "recorded_cuts": "swap = S, parity = P, phase = i S P; the cap is an eigenstate of P (line in one copy)",
        "recorded_marks_need": "no shared line and no phase; equal copies give rate q = 1",
        "zero_modes_at_a_sign_mark": "GROUP_A1_TABLE_REPRODUCED: same algebra as the recorded D8 chart",
    }


def twist_rule():
    """One recorded charge function for all mark types: F_alpha(Q) = Q^2/2 + (alpha - 1/2) Q."""
    energy = lambda alpha, charge: sp.Rational(charge*charge, 2)+(alpha-sp.Rational(1, 2))*charge
    swap, sign = sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
    rows = {}
    for name, matrix in (("PS", sign*swap), ("P", sign), ("S", swap)):
        twists = sorted((sp.Rational(1, 2)+sp.nsimplify(sp.arg(value)/(2*sp.pi))) % 1
                        for value, mult in matrix.eigenvals().items() for _ in range(mult))
        ground = []
        for alpha in twists:
            levels = {charge: energy(alpha, charge) for charge in range(-4, 5)}
            low = min(levels.values())
            ground.append(sum(value == low for value in levels.values()))
        rows[name] = {"twists": [str(t) for t in twists], "ground_states_per_fermion": ground}
    require(rows["PS"] == {"twists": ["1/4", "3/4"], "ground_states_per_fermion": [1, 1]},
            "quarter marks have a unique empty state (R066 law)")
    require(rows["P"] == rows["S"] == {"twists": ["0", "1/2"], "ground_states_per_fermion": [2, 1]},
            "sign and exchange marks carry one doublet per fermion (flat sector)")
    return {"marks": rows, "flat_sector_dimension_for_eight_fermions": 2**8,
            "matter_half": 2**7,
            "native_check": "flat without potential, not flat in any of the 27 native configurations"}


def charged_seam(max_order=8, max_fermions=64):
    """Two real copies as one charged field; search over rotation angle and number of fermions."""
    turn, swap, sign = sp.Matrix([[0, 1], [-1, 0]]), sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
    plus, minus = sp.Matrix([1, sp.I]), sp.Matrix([1, -sp.I])
    require(turn.T*plus == -sp.I*plus and sign.T*plus == minus and swap.T*plus == sp.I*minus,
            "quarter turn rotates the charge, sign and exchange conjugate it")
    marks = [1, sp.I, -1, -sp.I]
    copies = [sp.Matrix(v) for v in ([1, 0], [0, 1], [-1, 0], [0, -1])]
    on_marks = [tuple(marks.index(sp.simplify(f(m))) for m in marks)
                for f in (lambda z: sp.I*z, lambda z: sp.I*sp.conjugate(z), sp.conjugate)]
    on_copies = [tuple(copies.index(mat.T*v) for v in copies) for mat in (turn, swap, sign)]
    require(on_marks == on_copies, "square of marks and two copies carry the same dihedral action")
    # a rotation of the charge by 2 pi k/n on N complex fermions has weight N (k/n)^2 / 2 at each mark
    hits = [{"fermions": n_f, "angle_over_2pi": str(sp.Rational(k, n))}
            for n_f in range(1, max_fermions+1) for n in range(2, max_order+1) for k in range(1, n//2+1)
            if sp.gcd(k, n) == 1 and sp.Rational(n_f*k*k, 2*n*n) == sp.Rational(1, 2)]
    quarter = [row for row in hits if row["angle_over_2pi"] == "1/4"]
    require(quarter == [{"fermions": 16, "angle_over_2pi": "1/4"}],
            "quarter angle with mark weight 1/2 needs sixteen complex fermions")
    return {"dihedral_order": 8, "recorded_mirror_is": "exchange of the two copies",
            "weight_half_solutions": hits, "quarter_angle_solution": quarter[0],
            "fermions_per_copy_complex": 8,
            "purity_as_charge_count": "Tr rho_A^2 of one real copy = |<exp(i pi Q_A / 2)>| of the charged field",
            "native_check": "real structure maps shift 1/4 to shift 3/4 in all 27 configurations",
            # round 30, check_native_quarter_count.py
            # corrected in round 33 (check_native_dirac_cap.py): the native potential term is -i gamma D_F
            "native_cuts": {"U": "grading (charge sign)", "Y": "direction of the native potential -i gamma D_F",
                            "cap": "unique real eigenline of D_F; eta = sign of its value",
                            "marks": ["1", "0", "0", "i*eta"], "rate": "1", "positive_lift": "1/2",
                            "status": "HYPOTHESIS_NATIVE_REALISATION_OF_RECORDED_CANONICAL_MODEL"},
            # round 31, check_native_charge_marks.py
            "native_charge_marks": {"scalar_marker": "real", "charge_marker": "imaginary", "mixed": "imaginary",
                                    "electric_single": "0", "flat_density_hypercharge_mixed": "imaginary, not zero",
                                    "matches_recorded_E8_contract_sign": False},
            # round 32, check_native_real_lines.py
            "native_real_lines": {"quarter_shift": "none in 27 configurations",
                                  "shift_zero": "12 unique real eigenlines, Majorana neutrino lines",
                                  "marks_on_them": ["1", "0", "+-1", "0"], "eta_fixed": False},
            # round 34, check_native_candidate_marks.py: H-D evaluated on all 27 native configurations
            "native_candidate": {"cap": "mode 0 times lightest unique eigenline of D_F", "eta": "sign of the mass value",
                                 "baseline_marks": ["1", "0", "0", "i*eta"], "rate_range": ["0.98842", "1"],
                                 "rate_below_one_iff": "profile not constant", "original_producer": False,
                                 # round 35, check_candidate_line_and_step.py
                                 "line_rule": "lowest level of |h_N| at the baseline", "zero_step_marks_exact": True,
                                 "step_fixed_by_record": False},
            # round 36, check_commutator_rule_and_note16.py
            "commutator_rule": {"statement": "<HXHX> - <H^2><X^2> = -(1/2)<|[H,X]|^2> on the flat state",
                                "group_A1": ["-1/32", "-1/72"], "group_A9": "ABAB = -1 for anticommuting cuts",
                                "native_hypercharge": "-0.129976", "common_quarter_rotation_fourth_power_on_matter": "+1"}}


def controls(limit=4):
    rows = search(limit)
    integer = [row for row in rows if row["integer_spin"]]
    ordered = [row for row in rows if row["fermion_pairs_exchanged"] == "16"]
    require(len(integer) == 1 and integer[0]["rule"] == (1, 1, 1), "integer spin selects one rule")
    require(len(ordered) == 1 and ordered[0]["rule"] == (1, 0, 0), "sixteen free fermions select one rule")
    return {
        "status": "EXACT_FINITE_RULE_SEARCH",
        "hypothesis": "H-R2: the weighted reflected form on the four marks is the purity of two opposite arcs",
        "witnesses": witnesses(),
        "reflections": reflections(),
        "twist_rule": twist_rule(),
        "charged_seam": charged_seam(),
        "survivors_before_spin_condition": rows,
        "integer_spin_survivor": integer[0],
        "free_seam_survivor": ordered[0],
        "assumptions": ["the P1 boundary form is the exchange of two copies on the arcs",
                        "the seam state is pure", "the recorded mark weight 1/2 fixes one octet"],
        "first_native_test": ("per octet in the continuum the ratio at the square is 3/2 or 2; the native operator has "
                              "four circle modes, hence four places: check_native_sea_purity.py of the session"),
        "original_r113_values": None,
        "physical_source_identity_proved": False,
    }
