"""Classify the simultaneous three-cycle and quarter action on one marked sphere."""

from itertools import permutations

import sympy as sp

from reconstruction import require


THREE_CYCLE = (1,2,0,3)
QUARTER = (1,2,3,0)
PERMUTATIONS = tuple(permutations(range(4)))
MODULUS = sp.Symbol("lambda")


def _cross_ratios():
    points = ((1,0),(0,1),(1,1),(MODULUS,1))

    def det(a,b):
        return a[0]*b[1]-a[1]*b[0]

    def cross(p):
        a,b,c,d = (points[i] for i in p)
        return sp.cancel(det(a,c)*det(b,d)/(det(a,d)*det(b,c)))

    return {perm:cross(perm) for perm in PERMUTATIONS}


def _symmetries(ratios, value):
    holomorphic, antiholomorphic = [], []
    for perm, ratio in ratios.items():
        change = ratio.subs(MODULUS,value)
        if sp.simplify(change-value)==0:
            holomorphic.append(perm)
        if sp.simplify(change-sp.conjugate(value))==0:
            antiholomorphic.append(perm)
    return holomorphic, antiholomorphic


def _conjugate(action, permutation):
    return tuple(action[permutation[action.index(i)]] for i in range(4))


def controls():
    ratios = _cross_ratios()
    family, quarter = ratios[THREE_CYCLE], ratios[QUARTER]
    family_polynomial = sp.fraction(sp.cancel(family-MODULUS))[0]
    quarter_polynomial = sp.fraction(sp.cancel(quarter-MODULUS))[0]
    require(sp.gcd(family_polynomial,quarter_polynomial)==1 and
            len(set(ratios.values()))==6,"the two holomorphic permutations must have disjoint fixed moduli")
    candidates = sorted(sp.solve(family_polynomial,MODULUS),key=sp.default_sort_key)
    require(len(candidates)==2 and all(sp.simplify(family.subs(MODULUS,v)-v)==0 and
                                      sp.simplify(quarter.subs(MODULUS,v)-sp.conjugate(v))==0
                                      for v in candidates), "the semilinear branch must be the conjugate pair")
    square_holo,square_anti = _symmetries(ratios,sp.Integer(2))
    generic_holo,generic_anti = _symmetries(ratios,sp.Integer(3))
    modulus = (1+sp.I*sp.sqrt(3))/2
    tetra_holo,tetra_anti = _symmetries(ratios,modulus)
    require([len(x) for x in (square_holo,square_anti,generic_holo,generic_anti,tetra_holo,tetra_anti)]==
            [8,8,4,4,12,12] and set(tetra_holo).isdisjoint(tetra_anti) and
            QUARTER in tetra_anti and THREE_CYCLE in tetra_holo, "incorrect marked-sphere symmetry action")

    z,radius_squared = sp.symbols("z t")

    def apply(action, point):
        a,b,c,d = action
        n,m = point
        return sp.Matrix([a*n+b*m,c*n+d*m])

    points = (sp.Matrix([1,0]),sp.Matrix([0,1]),sp.Matrix([1,1]),sp.Matrix([modulus,1]))
    family_map = (0,1,-1,1)             # 1/(1-z)
    quarter_map = (0,1,-modulus,1)       # 1/(1-lambda*conj(z))
    reflected = (0,2,1,3)
    reflection_map = (-1,1,0,1)        # 1-conj(z)
    for action, permutation, anti in ((family_map,THREE_CYCLE,False),(quarter_map,QUARTER,True),
                                      (reflection_map,reflected,True)):
        for index,point in enumerate(points):
            image = apply(action,point.conjugate() if anti else point)
            require(sp.simplify(image.row_join(points[permutation[index]]).det())==0,
                    "the projective map does not act on the four original labels as declared")
    # A fixed point of the quarter would solve z=1+lambda*|z|^2.
    # Then |z|^2=t requires t = |1+lambda*t|^2, which has no real root.
    modulus_norm = sp.expand((1+modulus*radius_squared)*(1+sp.conjugate(modulus)*radius_squared)-radius_squared)
    require(modulus_norm==1+radius_squared**2,"the semilinear quarter should have no real fixed-point radius")
    elliptic_j = 256*(1-MODULUS+MODULUS**2)**3/(MODULUS**2*(1-MODULUS)**2)
    require(sp.simplify(elliptic_j.subs(MODULUS,modulus))==0 and elliptic_j.subs(MODULUS,2)==1728,
            "tetrahedral and harmonic elliptic moduli coincide")
    # On C/(Z+lambda Z), the four ramification points are the half-periods.
    # Encode their doubled lattice coordinates; the maps preserve the
    # lattice because lambda^2=lambda-1 and conjugate(lambda)=1-lambda.
    require(sp.simplify(modulus**2-modulus+1)==0 and modulus+sp.conjugate(modulus)==1,
            "invalid hexagonal lattice basis")
    periods = ((0,1),(1,0),(1,1),(0,0))
    for k,(a,b) in enumerate(periods):
        require(periods[THREE_CYCLE[k]]==((a+b)%2,a%2) and periods[QUARTER[k]]==((a+b)%2,(1-b)%2),
                "elliptic lifts do not reproduce the four branch permutations")
    # Q(z)=conjugate(z)+lambda/2 squares to translation by 1/2.
    # The other lift differs by the deck involution and has the same square
    # modulo the lattice. On differentials both anti lifts square to +1.
    lift_square = sp.conjugate(sp.conjugate(z)+modulus/2)+modulus/2
    require(sp.simplify(lift_square-z)==sp.S.Half,"the tetrahedral quarter lift lost its half-period translation")

    # If all even permutations preserved one circle, A4 would act on
    # it faithfully. Every 3-cycle preserves the circle orientation and
    # generates A4 with its conjugates. A finite oriented circle group
    # is cyclic (its action on any point orbit preserves the cyclic order
    # and has trivial stabilizers). The two exhibited 3-cycles do not commute.
    conjugate_cycle = _conjugate(QUARTER,THREE_CYCLE)
    left = tuple(THREE_CYCLE[conjugate_cycle[i]] for i in range(4))
    right = tuple(conjugate_cycle[THREE_CYCLE[i]] for i in range(4))
    require(conjugate_cycle in tetra_holo and left!=right,"the circle obstruction needs two noncommuting even rotations")
    reflection_orbit = sorted({_conjugate(action,reflected) for action in tetra_holo})
    require(len(reflection_orbit)==6 and all(reflection in tetra_anti and
            sum(i==j for i,j in enumerate(reflection))==2 for reflection in reflection_orbit),
            "the marked cut orbit should consist of the six reflective transpositions")

    # Even allowing all conjugated cuts, effective Möbius generator fields
    # still lie in sl2(C). Its bracket is invertible on the complex wedge
    # square, so commuting fields are complex-collinear.
    cartan = sp.diag(1,-1)
    rising, falling = sp.Matrix([[0,1],[0,0]]),sp.Matrix([[0,0],[1,0]])
    brackets = [a*b-b*a for a,b in ((cartan,rising),(cartan,falling),(rising,falling))]
    bracket = sp.Matrix([[value[0,0],value[0,1],value[1,0]] for value in brackets]).T
    require(bracket.det()==-4,"unexpected abelian rank of the complex Möbius generators")
    return {"status":"EXACT_MARKED_SPHERE_ACTION_CLASSIFICATION",
            "marked_points":["infinity","0","1","lambda"],
            "three_cycle":list(THREE_CYCLE),"quarter":list(QUARTER),
            "both_holomorphic_moduli":[],"semilinear_moduli":list(map(str,candidates)),
            "holomorphic_group_orders":{"generic":len(generic_holo),"harmonic":len(square_holo),
                                        "equianharmonic":len(tetra_holo)},
            "semilinear_tetrahedral_group_order":len(tetra_holo)+len(tetra_anti),
            "harmonic_elliptic_j":1728,"tetrahedral_elliptic_j":0,
            "example_maps":{"family":str(1/(1-z)),"quarter":str(1/(1-modulus*sp.conjugate(z))),
                            "reflection":str(1-sp.conjugate(z))},
            "quarter_fixed_points":0,"reflection_fixed_circle":"Re(z)=1/2, with infinity",
            "elliptic_quarter_square_translation":"1/2 in the lattice Z+lambda*Z",
            "quarter_square_on_elliptic_differential":"+1",
            "primitive_complex_linear_quarter_on_same_elliptic_line":False,
            "common_invariant_reflection_circle_for_even_mark_action":False,
            "noncommuting_order_three_circle_witness":[list(THREE_CYCLE),list(conjugate_cycle)],
            "two_mark_reflection_orbit":[list(reflection) for reflection in reflection_orbit],
            "effective_mobius_translation_rank_bound_over_R":2,
            "three_spatial_translations_from_mobius_only":False,
            "next_original_observation":"Bind how the native family action and the quarter move the actual reflected cuts. To realize three commuting spatial directions, the physical regional action must extend the effective Möbius action on this sphere.",
            "physical_source_identity_proved":False}
