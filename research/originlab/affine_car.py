"""Free CAR source for the D8 pair currents and the neutral weak doublet word.

Eight complex free fermions in the Neveu-Schwarz sector carry so(16) at
level one. Their bilinears are the 112 D8 pair currents of the fixed D8
chart. The neutral doublet word J^(-delta) J^gamma J^alpha uses only such
currents, so a positive quadratic CAR source can evaluate it. The marked
word J^(-p1) J^gamma J^p2 has spinor ends and is outside this source.

Scope. This is a declared free model on the D8 vector fermions. These are
not the 96 native matter fibres, and no weighted R113 value is computed.
"""
from fractions import Fraction
from itertools import product

from affine_operators import (DOUBLET_ALPHA, DOUBLET_DELTA, MARKED_GAMMA, MARKED_P1, MARKED_P2, ROOTS,
                              ROOT_SET, WEAK_ROOT, LatticeFock, add, charge_signature, dot, integral,
                              lattice_frame, marked_word, neg, require)

VACUUM = {(): Fraction(1)}


def fermion(dagger, species, mode, vector):
    """psi_species,mode/2 or its conjugate field mode; the mode is an odd integer, doubled."""
    creates = mode < 0
    # psi creates particles (kind 0) and removes antiparticles (kind 1); psi-dagger the reverse.
    kind = int(dagger) if creates else 1-int(dagger)
    letter = (kind, species, abs(mode))
    out = {}
    for occupied, value in vector.items():
        if creates:
            if letter in occupied:
                continue
            position = sum(1 for other in occupied if other < letter)
            key = tuple(sorted((*occupied, letter)))
        else:
            if letter not in occupied:
                continue
            position = occupied.index(letter)
            key = occupied[:position]+occupied[position+1:]
        out[key] = out.get(key, 0)+(-value if position % 2 else value)
    return {key: value for key, value in out.items() if value}


def bilinear(root):
    """The ordered fermion pair of a D8 root: (dagger, species) twice."""
    places = [(i, x) for i, x in enumerate(root) if x]
    require(len(places) == 2 and all(abs(x) == 2 for _, x in places), "a D8 pair current is required")
    (i, a), (j, b) = places
    if a > 0 and b > 0:
        return (False, i), (False, j)
    if a < 0 and b < 0:
        return (True, j), (True, i)                   # adjoint of psi_i psi_j
    return ((False, i), (True, j)) if a > 0 else ((False, j), (True, i))


def current(root, mode, vector):
    """J^root_mode = sum over r+s=mode of the two fermion modes. No ordering term for distinct species."""
    first, second = bilinear(root)
    top = max((letter[2] for occupied in vector for letter in occupied), default=1)
    bound = top+2*abs(mode)+2
    out = {}
    for s in range(-bound, bound+1):
        if s % 2 == 0:
            continue
        r = 2*mode-s
        image = fermion(*first, r, fermion(*second, s, vector))
        for key, value in image.items():
            out[key] = out.get(key, 0)+value
    return {key: value for key, value in out.items() if value}


def number(species, vector):
    """Zero mode of the Cartan current of one species."""
    out = {}
    for occupied, value in vector.items():
        charge = sum(1 if kind == 0 else -1 for kind, index, _ in occupied if index == species)
        if charge:
            out[occupied] = value*charge
    return out


def inner(left, right):
    return sum(value*right.get(key, 0) for key, value in left.items())


def scale(factor, vector):
    return {key: factor*value for key, value in vector.items() if factor*value}


def plus(*vectors):
    out = {}
    for vector in vectors:
        for key, value in vector.items():
            out[key] = out.get(key, 0)+value
    return {key: value for key, value in out.items() if value}


def weak_casimir(vector):
    """Casimir of the weak su(2) zero modes, E = psi_4 psi_5*, F = E*, H = N_4 - N_5."""
    up, down = WEAK_ROOT, neg(WEAK_ROOT)
    cartan = plus(number(3, vector), scale(-1, number(4, vector)))
    square = plus(number(3, cartan), scale(-1, number(4, cartan)))
    return plus(scale(Fraction(1, 4), square), scale(Fraction(1, 2), current(up, 0, current(down, 0, vector))),
                scale(Fraction(1, 2), current(down, 0, current(up, 0, vector))))


def triangles_through_gamma():
    """All neutral three-leg words through the marked middle current, by class."""
    rows = {}
    for alpha in ROOTS:
        if dot(MARKED_GAMMA, alpha) != -1:
            continue
        third = neg(add(MARKED_GAMMA, alpha))
        require(third in ROOT_SET, "pairing -1 must close a triangle")
        key = frozenset((alpha, third))               # one triangle per unordered pair of outer legs
        pair = all(x % 2 == 0 for x in alpha) and all(x % 2 == 0 for x in third)
        odd = dot(alpha, WEAK_ROOT) % 2 == 1
        rows[key] = (odd, pair)
    counts = {}
    for odd, pair in rows.values():
        name = ("weak_doublet" if odd else "weak_even")+("_pair_currents" if pair else "_spinor_ends")
        counts[name] = counts.get(name, 0)+1
    return counts


def controls(size=5):
    alpha, gamma, delta = DOUBLET_ALPHA, MARKED_GAMMA, DOUBLET_DELTA
    require(all(all(x % 2 == 0 for x in root) for root in (alpha, gamma, delta)),
            "the doublet word must consist of D8 pair currents")
    require(any(x % 2 for x in MARKED_P1) and any(x % 2 for x in MARKED_P2),
            "the marked ends must be spinor currents outside the free source")

    # Level one and positivity on the first charged vectors.
    for root in (alpha, gamma, delta, WEAK_ROOT):
        for m in range(1, size+1):
            vector = current(root, -m, VACUUM)
            require(inner(vector, vector) == m, "free pair current must have level-one norm m")
            back = current(neg(root), m, vector)
            require(back == {(): Fraction(m)}, "the adjoint pair current must return m times the vacuum")

    # The neutral doublet word and the weak Casimir inside it.
    word, inserted = {}, {}
    for m, n in product(range(1, size+1), repeat=2):
        first = current(alpha, -n, VACUUM)
        word[m, n] = current(neg(delta), m, current(gamma, n-m, first)).get((), Fraction(0))
        inserted[m, n] = current(neg(delta), m, current(gamma, n-m, weak_casimir(first))).get((), Fraction(0))
    orientation = word[1, 1]
    require(orientation in (1, -1) and all(orientation*value == min(m, n) for (m, n), value in word.items()),
            "free source must reproduce the kernel min(m,n)")
    require(all(inserted[key] == Fraction(3, 4)*word[key] for key in word),
            "weak Casimir must give 3/4 inside the free doublet word")
    singlet = current(gamma, -2, VACUUM)
    triplet = current(WEAK_ROOT, -1, VACUUM)
    require(weak_casimir(singlet) == {} and weak_casimir(triplet) == scale(2, triplet),
            "weak Casimir controls 0 and 2 failed")

    # Agreement with the lattice operators, which also carry the spinor ends.
    frame, _ = lattice_frame()
    fock = LatticeFock([[dot(a, b) for b in frame] for a in frame])
    a, g, d = (integral(frame, r) for r in (alpha, gamma, delta))
    lattice = {}
    for m, n in word:
        state = fock.vertex(tuple(-x for x in d), m, fock.vertex(g, n-m, fock.vertex(a, -n, fock.vacuum())))
        lattice[m, n] = state.get(((0,)*8, ()), Fraction(0))
    require(all(abs(lattice[key]) == abs(word[key]) for key in word), "lattice and free sources must agree")
    marked = marked_word(fock, frame, 2)
    require(abs(marked[2, 2]) == 2, "the marked word exists on the lattice source")

    pairs = [r for r in ROOTS if all(x % 2 == 0 for x in r)]
    counts = triangles_through_gamma()
    require(len(pairs) == 112 and sum(counts.values()) == 28, "triangle census through gamma changed")
    return {
        "status": "EXACT_FREE_CAR_DOUBLET_SOURCE",
        "complex_fermions": 8, "sector": "Neveu-Schwarz", "pair_currents": len(pairs), "level": 1,
        "doublet_word": "J^(-delta) J^gamma J^alpha with all three legs fermion bilinears",
        "doublet_word_coefficients": f"min(m,n) for m,n <= {size}",
        "weak_casimir_inside_word": "3/4",
        "weak_casimir_controls": {"singlet": "0", "triplet": "2"},
        "agrees_with_lattice_operators": True,
        "marked_word_in_free_source": False,
        "marked_word_reason": "its ends p1 and p2 are spinor currents, which are not fermion bilinears",
        "triangles_through_marked_gamma": counts,
        "doublet_signatures": {name: [charge_signature(r)[1], str(charge_signature(r)[2])]
                               for name, r in (("alpha", alpha), ("gamma", gamma), ("delta", delta))},
        "physical_source_identity_proved": False, "weighted_r113_values": None,
        "scope": "Declared free model on the sixteen real vector fermions of D8. The native matter "
                 "fibres are spinors of D5 and are a different carrier. A quadratic native source "
                 "could evaluate this word only after its own field map to these currents is given.",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(controls(), indent=2))
