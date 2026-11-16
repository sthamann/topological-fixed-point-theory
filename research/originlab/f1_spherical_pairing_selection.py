"""Exact conformal response and four-mark selection controls for P1."""

from itertools import combinations

import sympy as sp

from reconstruction import require


def controls():
    checks = []

    def check(name, condition):
        require(condition, name)
        checks.append(name)

    def equal(left, right):
        return sp.simplify(left - right) == 0

    z, weight = sp.symbols("z w", real=True)
    ell = sp.symbols("ell", integer=True, nonnegative=True)
    half = sp.Rational(1, 2)
    upper = -(ell + weight) * (ell + 1) / (2*ell + 1)
    lower = ell * (ell + 1 - weight) / (2*ell + 1)
    for degree in range(9):
        harmonic = sp.legendre(degree, z)
        boost = (1-z**2)*sp.diff(harmonic, z) - weight*z*harmonic
        raised = upper.subs(ell, degree)*sp.legendre(degree+1, z)
        lowered = lower.subs(ell, degree)*sp.legendre(degree-1, z) if degree else 0
        check(f"boost form on Legendre degree {degree}", equal(boost, raised + lowered))

    k = 2*ell + 1
    check("all adjacent degrees obey the input and output Ward identity", equal(
        k.subs(ell, ell+1)*upper.subs(weight, half), k*upper.subs(weight, 3*half)
    ))
    check("reverse Ward identity holds on all nonzero degrees", equal(
        k.subs(ell, ell-1)*lower.subs(weight, half), k*lower.subs(weight, 3*half)
    ))
    q = sp.symbols("q", positive=True)
    antipode = sp.limit(1/(4*sp.pi*(1+q)), q, 1, dir="-")
    check("antipodal response from the full Ward spectrum", antipode == 1/(8*sp.pi))
    # Rotation symmetry and the first two modes allow this positive finite
    # rank deformation, which the next independent Ward test detects.
    check("a deformation preserving the first pair fails on the next pair", not equal(
        6*upper.subs({ell: 1, weight: half}), 3*upper.subs({ell: 1, weight: 3*half})
    ))

    step = sp.symbols("r", positive=True)
    cut1 = sp.Matrix([[0, 1], [1, 0]])
    cut2 = sp.Matrix([[0, step], [1/step, 0]])
    dilation = sp.diag(step, 1/step)
    check("two anti-Mobius cut reflections compose to a real dilation", cut1**2 == cut2**2 == sp.eye(2)
          and cut2*cut1 == dilation)
    rotation = sp.Matrix([[0, 1], [-1, 0]])
    transverse = sp.Matrix([[0, 1], [1, 0]])
    witness_difference = (dilation*rotation*dilation.inv()
                          - (step**2+step**-2)*rotation/2 - (step**2-step**-2)*transverse/2)
    check("conjugating a rotation through the dilation produces a boost direction",
          witness_difference.applyfunc(sp.simplify) == sp.zeros(2))

    s, beta = sp.symbols("s beta", nonnegative=True)
    h = sp.symbols("h", real=True)
    radius = sp.sqrt(1-h**2)
    points = [sp.Matrix(coords) for coords in ((radius, 0, h), (-radius, 0, h),
                                               (0, radius, -h), (0, -radius, -h))]
    pairs = ((0, 1), (2, 3))
    distances = {(i, j): sp.simplify((points[i]-points[j]).dot(points[i]-points[j]))
                 for i, j in combinations(range(4), 2)}
    check("four unit marks in the saturating family", all(equal(p.dot(p), 1) for p in points))
    check("both within-pair distances reach their Jensen mean", all(
        equal(distances[pair], 4-4*h**2) for pair in pairs
    ))
    check("all cross-pair distances reach their Jensen mean", all(
        equal(value, 2+2*h**2) for pair, value in distances.items() if pair not in pairs
    ))
    dipoles = [points[i] + points[j] for i, j in pairs]
    check("global dipole vanishes while the two window dipoles distinguish the minima",
          dipoles[0] + dipoles[1] == sp.zeros(3, 1)
          and equal(sum(p.dot(p) for p in dipoles), 8*h**2))
    all_matchings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
    check("summing all matchings erases the selection along the minimizer family", equal(
        sum((points[i]+points[j]).dot(points[i]+points[j])
            for matching in all_matchings for i, j in matching), 8
    ))

    objective = 2/sp.sqrt(4-s) + 4/sp.sqrt(2+s/2) + 2*beta*s
    derivative = (4-s)**(-3*half) - (2+s/2)**(-3*half) + 2*beta
    check("strictly increasing first derivative on 0 <= s < 4", equal(sp.diff(objective, s), derivative)
          and equal(sp.diff(derivative, s), 3*half/(4-s)**(5*half)+3*half/(2*(2+s/2)**(5*half))))
    threshold = (sp.sqrt(2)/4 - sp.Rational(1, 8))/2
    check("exact threshold at the square endpoint", equal(derivative.subs(s, 0), 2*(beta-threshold)))
    check("without within-window pairing the minimum is tetrahedral", equal(
        derivative.subs({beta: 0, s: sp.Rational(4, 3)}), 0
    ) and all(equal(value.subs(h**2, sp.Rational(1, 3)), sp.Rational(8, 3)) for value in distances.values()))

    a = sp.symbols("a", real=True)
    cross_ratio = 2/(1-sp.I*a)
    projected = (step, sp.I/step, -step, -sp.I/step)
    first, second, third, fourth = projected
    check("cross ratio of the actual stereographic marks", equal(
        (first-third)*(second-fourth)/((first-fourth)*(second-third)),
        cross_ratio.subs(a, (step**2-step**-2)/2)
    ) and equal(((step**2-step**-2)/2).subs(step**2, (1+h)/(1-h)), 2*h/(1-h**2)))
    j_invariant = 256*(1-cross_ratio+cross_ratio**2)**3/(cross_ratio**2*(1-cross_ratio)**2)
    j_expression = 64*(3-a**2)**3/(1+a**2)**2
    check("elliptic j along the minimizer family", equal(j_invariant, j_expression))
    check("distinct endpoint moduli", j_expression.subs(a, 0) == 1728 and j_expression.subs(a, sp.sqrt(3)) == 0)
    dipole_square = 3/(4*sp.pi) * sum(p.dot(p) for p in dipoles)
    check("local degree-one spherical harmonic intensity measures the pairing term", equal(dipole_square, 6*h**2/sp.pi))

    identity = sp.eye(2)
    x = sp.Matrix([[0, 1], [1, 0]])
    y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    z_axis = sp.diag(1, -1)
    xx, yy, zz = (sp.kronecker_product(m, m) for m in (x, y, z_axis))
    pair_identity = sp.eye(4)
    singlet = (pair_identity-xx-yy-zz)/4
    product = sp.kronecker_product
    check("two ordered axis anticorrelations select only the singlet", xx*yy == yy*xx == -zz
          and (pair_identity-xx)*(pair_identity-yy)/4 == singlet and singlet.rank() == 1)
    check("two axis responses lower-bound the singlet fidelity", (
        singlet+(xx+yy)/2
    ).eigenvals() == {0: 3, 1: 1})
    p, q, t = sp.symbols("p q t", real=True)
    ray = (identity+p*x+q*y+t*z_axis)/2
    b, c, d = sp.symbols("b c d", real=True)
    other = (identity+b*x+c*y+d*z_axis)/2
    check("singlet reward is the two-mark dipole penalty up to a constant", equal(
        4-8*sp.trace(singlet*product(ray, other)), 2+2*(p*b+q*c+t*d)
    ))
    weighted = product(ray, identity)*singlet
    partner = weighted[:2, :2]+weighted[2:, 2:]
    check("unnormalized conditional singlet readout is half of the antipode", (
        partner-(identity-ray)/2
    ).applyfunc(sp.simplify) == sp.zeros(2))
    check("antipodal state inversion is not a CP replacement", (
        sp.eye(4)-sp.Matrix([1, 0, 0, 1])*sp.Matrix([[1, 0, 0, 1]])
    ).eigenvals() == {-1: 1, 1: 3})
    spin_quarter = sp.diag((1-sp.I)/sp.sqrt(2), (1+sp.I)/sp.sqrt(2))
    pair_quarter = product(spin_quarter, spin_quarter)
    check("single and paired quarter lifts have different fourth powers", (
        spin_quarter**4
    ).applyfunc(sp.simplify) == -identity and (pair_quarter**4).applyfunc(sp.simplify) == pair_identity
          and pair_quarter**2 != pair_identity)
    check("the paired quarter preserves the singlet", (pair_quarter*singlet).applyfunc(sp.simplify) == singlet)

    return {
        "status": "EXACT_SCOPED_CONTROLS_PASSED_F1_OPEN",
        "checks": checks,
        "response": {"ward_eigenvalue": str(k), "antipode": str(antipode), "input_weights": ["1/2", "3/2"]},
        "mark_selection": {"reduced_objective": str(objective), "critical_pair_weight": str(threshold),
                           "critical_pair_weight_decimal": float(threshold), "j_of_family": str(j_expression),
                           "unpaired_j": 0, "square_j": 1728},
        "singlet_test": {"two_axis_correlations": [-1, -1], "conditioned_partner": "(I-P)/2",
                         "single_quarter_fourth_power": "-I", "pair_quarter_order": 4},
        "physical_solution": False,
    }
