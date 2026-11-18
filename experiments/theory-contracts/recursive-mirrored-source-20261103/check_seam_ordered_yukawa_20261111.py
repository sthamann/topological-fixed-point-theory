"""Original E8 vacuum Ward continuations in the source's marked charged slots."""
from itertools import combinations, permutations, product
import json
from pathlib import Path
import sys

import sympy as sp

repo = Path(sys.argv[1]).resolve()
contracts = repo / "experiments/theory-contracts"
sys.path.insert(0, str(contracts / "seam-collar-multifield-20261016"))
from e8_boundary_response import current_response, epsilon, vertex_response

source = contracts / "compiler-vacuum-current-cubic-20260918/certificate.json"
cubic = json.loads(source.read_text())["cubic"]
base = [tuple(sp.Rational(x, 2) for x in p) for p in cubic["base_roots_doubled"]]
labels = {}
for i, p in enumerate(base):
    if all(abs(2*x) == 1 for x in p[:5]):
        occupation = [x + sp.Rational(1, 2) for x in p[:5]]
        names = {(0, 0): "nu", (1, 1): "Q", (2, 0): "u", (0, 2): "e",
                 (2, 2): "d", (3, 1): "L"}
        labels[i] = names[sum(occupation[:3]), sum(occupation[3:])]
higgs_up = next(i for i, p in enumerate(base) if p[:5] == (0, 0, 0, 0, 1))
higgs_down = next(i for i, p in enumerate(base) if p[:5] == (0, 0, 0, 0, -1))
slots = {"e": [], "u": [], "d": []}
for entry in cubic["monomials"]:
    for tag, higgs, left in (("e", higgs_down, "L"), ("u", higgs_up, "Q"), ("d", higgs_down, "Q")):
        if higgs not in entry["indices"]:
            continue
        matter = [i for i in entry["indices"] if i != higgs]
        if {labels.get(i) for i in matter} == {left, tag}:
            li = next(i for i in matter if labels[i] == left)
            ri = next(i for i in matter if labels[i] == tag)
            slots[tag].append((li, ri, higgs))
assert [len(slots[tag]) for tag in ("e", "u", "d")] == [1, 3, 3]

def vector(*tail):
    return (0,) * 5 + tail

def charge(i, j):
    return sum(a*b for a, b in zip(i, j))

def matter_vertex(slot, family):
    p = base[slot]
    if family == 0:
        return p, 1
    tail = [-1, 0, 0]
    tail[family] = 1
    shift = vector(*[-x for x in tail])
    return tuple(a+b for a, b in zip(p, shift)), epsilon(shift, p)

em = (sp.Rational(-1, 3),) * 3 + (1, 0, 0, 0, 0)
family_cartan = vector(1, 0, -1)
alpha = vector(1, -1, 0)
beta = vector(0, 1, -1)

def hermitian_roots(root):
    dual = tuple(-x for x in root)
    return ((root, 1), (dual, epsilon(dual, root)))

def one_matrix(slot, higgs_family, extra=(), cartan=()):
    position = (3, 2, 1)
    out = sp.zeros(3)
    for i, j in product(range(3), repeat=2):
        p, p_phase = matter_vertex(slot[0], i)
        q, q_phase = matter_vertex(slot[1], j)
        h, h_phase = matter_vertex(slot[2], higgs_family)
        phase = p_phase*q_phase*h_phase
        direct = current_response(tuple(extra)+(p,q,h), tuple(range(3+len(extra),3,-1))+position, cartan)
        swapped = current_response(tuple(extra)+(q,p,h), tuple(range(3+len(extra),3,-1))+position, cartan)
        out[i,j] = phase * (direct + swapped) / 2
    return out

def raw_coefficient(slot):
    roots, signs = zip(*(matter_vertex(index, family) for family, index in enumerate(slot)))
    v = vertex_response(roots, (3,2,1))*sp.prod(signs)
    assert abs(v) == sp.Rational(1,2)
    return v

def normalized_matrices(slot):
    raw = raw_coefficient(slot)
    pieces = {}
    pieces["bare"] = [one_matrix(slot, h) / raw for h in range(3)]
    for name, d in (("em", em), ("family", family_cartan)):
        pieces[name] = [one_matrix(slot, h, cartan=((d,4),)) / raw for h in range(3)]
    for key, cartan in (("paired_roots", ()), ("paired_roots_em", ((em,6),))):
        pieces[key] = [
            sum((s*t*one_matrix(slot, h, extra=(p,q),cartan=cartan)
                 for (p,s),(q,t) in product(hermitian_roots(alpha),hermitian_roots(beta))),
                sp.zeros(3)) / raw
            for h in range(3)
        ]
    return pieces

coefficients = {tag: [normalized_matrices(slot) for slot in seq] for tag, seq in slots.items()}
for family in range(3):
    for tag, seq in coefficients.items():
        for block in seq:
            assert block["bare"][family] == sp.zeros(3)
    for tag in ("e","u","d"):
        for block in coefficients[tag]:
            assert block["family"][family] == coefficients["e"][0]["family"][family]
            assert block["paired_roots"][family] == coefficients["e"][0]["paired_roots"][family]

h = sp.Matrix(sp.symbols("h1:4", real=True))
def on_higgs(mats):
    return sum((x*m for x,m in zip(h,mats)), sp.zeros(3))

e = coefficients["e"][0]
u = coefficients["u"][0]
d = coefficients["d"][0]
assert [charge(em, base[i]) for i in slots["e"][0]] == [-1, 1, 0]
assert [charge(em, base[i]) for i in slots["u"][0]] == [sp.Rational(2,3),sp.Rational(-2,3),0]
assert on_higgs(e["em"]) == -sp.Rational(3,2)*on_higgs(u["em"])
assert on_higgs(u["em"]) == -2*on_higgs(d["em"])
assert on_higgs(e["em"]).rank() == 2
print("Native one-current electric Ward insertion: e/u=3/2 and u/d=2 on nonzero slots; rank two.")

b = on_higgs(e["family"])
assert b.T == b and b.diagonal() == sp.zeros(1,3)
gram = b.T*b
assert sp.expand(sp.trace(gram)**2 - 2*sp.trace(gram**2)) == 0
assert b.det() != 0
print("Native one-current family Ward insertion has determinant",sp.factor(b.det()),
      "and a degenerate singular-value triangle.")

ordered = on_higgs(e["paired_roots"])
print("Two ordered Hermitian family root currents:", ordered)
print("Two-current determinant:", sp.factor(ordered.det()))
print("Two-current diagonal:", ordered.diagonal())
assert ordered.T == ordered and ordered.det() != 0
higgs_test = {h[0]: 1, h[1]: 1, h[2]: 1}
numeric_gram = ordered.subs(higgs_test).T*ordered.subs(higgs_test)
assert sp.discriminant(numeric_gram.charpoly().as_expr(),numeric_gram.charpoly().gen) != 0
print("Normalized two-current matrices agree in e, u and d, including all three colors.")

sys.path.insert(0, str(repo / "verification"))
from v18_quark_yukawa import lY, phi0, yuk
from mpmath import mp
p_ratio = mp.mpf(55)/117 * (mp.mpf(34)/47) * (mp.mpf(3)/26)
reduced_ud = p_ratio * (1-phi0)**2 / phi0
reduced_eu = mp.fprod(yuk(x) for x in ("e","mu","tau")) / (
    lY**6 * mp.fprod(yuk(x) for x in ("u","c","t")))
assert mp.mpf("0.65") < reduced_ud < mp.mpf("0.68")
assert reduced_eu > mp.mpf("1.5")
for fields in (("e","mu","tau"),("u","c","t"),("d","s","b")):
    singulars = sorted([yuk(f) for f in fields])
    assert singulars[-1] > sum(singulars[:-1])
print("Largest source singular value exceeds the sum of the other two in every sector.")
print("Determinant after removal of all length factors: u/d (exact source formula) =",mp.nstr(reduced_ud,12),
      "; e/u (old table) =",mp.nstr(reduced_eu,12))

# Add the original neutral EM Ward current at the same ordered points in all sectors.
# Its coefficient and the Higgs-family vector remain variables.
electric = {}
for tag, rows in slots.items():
    for idx, slot in enumerate(rows):
        q = charge(em,base[slot[0]])
        assert charge(em,base[slot[1]]) == -q and charge(em,base[slot[2]]) == 0
        shape = on_higgs(coefficients[tag][idx]["paired_roots_em"]) / q
        electric[tag] = shape
        assert shape.T == -shape and shape != sp.zeros(3)
        assert shape == electric["e"]
q = sp.symbols("q", real=True)
mixed = sp.Poly(sp.expand((ordered + q*electric["e"]).det()),q)
assert mixed.degree() == 2 and mixed.nth(1) == 0
det_e,det_u,det_d = (mixed.eval(x) for x in (1,sp.Rational(2,3),sp.Rational(1,3)))
assert sp.expand(3*det_e-8*det_u+5*det_d) == 0
print("Adding one neutral EM current yields B(q)=B0+q B1, B1 skew; exact signed identity 3 De-8 Du+5 Dd=0.")

# Compare against a generous envelope from the printed e/u amplitudes.
from v18_quark_yukawa import TAB
half_printing = mp.mpf("0.0005")
min_eu = mp.fprod(mp.mpf(TAB[x][1])-half_printing for x in ("e","mu","tau")) / (
         mp.fprod(mp.mpf(TAB[x][1])+half_printing for x in ("u","c","t")))
max_eu = mp.fprod(mp.mpf(TAB[x][1])+half_printing for x in ("e","mu","tau")) / (
         mp.fprod(mp.mpf(TAB[x][1])-half_printing for x in ("u","c","t")))
target_e_over_d = (min_eu*reduced_ud, max_eu*reduced_ud)
envelopes = []
for sign_e, sign_u in product((-1,1),repeat=2):
    envelope = sorted(3*sign_e*x - 8*sign_u*reduced_ud + 5 for x in target_e_over_d)
    assert envelope[0] * envelope[1] > 0
    envelopes.append(min(abs(x) for x in envelope))
print("All four real determinant sign choices fail even with table rounding. Closest residual:",
      mp.nstr(min(envelopes),12))


# Coincide the two matter currents at 2, keeping family currents and Higgs separated.
# This extracts only the simple pole and regular term of the original lattice OPE.
def local_coefficients(roots):
    if any(sum(column) != 0 for column in zip(*roots)):
        return 0, 0
    left, right = len(roots)-3, len(roots)-2
    positions = tuple(range(len(roots),3,-1))+(2,2,1)
    power = charge(roots[left], roots[right])
    assert power >= -1
    value = sp.S.One
    derivative = sp.S.Zero
    for i, j in combinations(range(len(roots)), 2):
        value *= epsilon(roots[i], roots[j])
        if (i,j) == (left,right):
            continue
        distance = positions[i] - positions[j]
        exponent = charge(roots[i], roots[j])
        value *= sp.Integer(distance)**exponent
        if left in (i,j):
            derivative += exponent * (1 if i == left else -1) / distance
    return { -1: (value, value*derivative), 0: (0, value) }.get(power, (0, 0))

def local_matrices(slot, second_terms=None, first_terms=None, conjugate=False):
    if second_terms is None:
        second_terms = hermitian_roots(beta)
    if first_terms is None:
        first_terms = hermitian_roots(alpha)
    even, pole = sp.zeros(3), sp.zeros(3)
    raw = raw_coefficient(slot)
    for i,j,k in product(range(3), repeat=3):
        p, sp_ = matter_vertex(slot[0], i)
        q, sq = matter_vertex(slot[1], j)
        higgs, sh = conjugated_down_higgs(k) if conjugate else matter_vertex(slot[2], k)
        for (r,s),(t,v) in product(first_terms,second_terms):
            pm, fm = local_coefficients((r,t,p,q,higgs))
            pn, fn = local_coefficients((r,t,q,p,higgs))
            assert pm + pn == 0
            factor = h[k]*sp_*sq*sh*s*v/(2*raw)
            even[i,j] += factor*(fm+fn)
            pole[i,j] += factor*(pm-pn)
    return even,pole

local_even,local_pole = local_matrices(slots["e"][0])
for slots_in_sector in slots.values():
    for slot in slots_in_sector:
        assert local_matrices(slot) == (local_even,local_pole)
assert local_even.T == local_even and local_pole.T == -local_pole
assert local_even.det() != 0 and local_pole != sp.zeros(3)
local_gram = local_even.subs(higgs_test).T * local_even.subs(higgs_test)
char = local_gram.charpoly()
assert sp.discriminant(char.as_expr(), char.gen) != 0
local_det = sp.Poly(sp.expand((local_even+q*local_pole).det()),q)
assert local_det.degree() == 2 and local_det.nth(1) == 0
print("Coincident matter: regular symmetric tensor", local_even)
print("Coincident matter: original simple-pole tensor", local_pole)
print("Coincident matter: full rank for a real Higgs with h2 != 0; det =",sp.factor(local_even.det()))
print("Smooth vacuum-normalized EM dressing retains only the linear charge term at coincidence.")

# An invariant quadratic insertion contracts both currents with the native metric.
# It provides the comparison when no ordered pair of family directions is selected.
family_roots = [vector(*tail) for tail in set(permutations((1,-1,0)))]
for family in range(3):
    slot = slots["e"][0]
    casimir = sp.zeros(3)
    for p in family_roots:
        dual = tuple(-x for x in p)
        casimir += epsilon(dual,p)*one_matrix(slot,family,extra=(p,dual))
    for d in (alpha,vector(1,1,-2)):
        casimir += one_matrix(slot,family,cartan=((d,5),(d,4)))/charge(d,d)
    assert casimir == sp.zeros(3)
print("The full invariant quadratic SU(3)_F current insertion gives zero in the symmetric channel.")

# Independent check of the regular term by the unchanged separated-point engine.
x = sp.Rational(1,10**9)
for p,r in product(range(3),repeat=2):
    left, s = matter_vertex(slots["e"][0][0],p)
    right, t = matter_vertex(slots["e"][0][1],r)
    root_h, u = matter_vertex(slots["e"][0][2],1)
    nearby = 0
    for (a,sa),(b,sb) in product(hermitian_roots(alpha),hermitian_roots(beta)):
        nearby += sa*sb*(vertex_response((a,b,left,right,root_h),(5,4,2+x,2,1))
                       + vertex_response((a,b,right,left,root_h),(5,4,2+x,2,1)))
    nearby *= s*t*u/(2*raw_coefficient(slots["e"][0]))
    target = local_even[p,r].subs({h[0]:0,h[1]:1,h[2]:0})
    assert abs(nearby-target) < sp.Rational(1,10**7)
print("The unchanged separated-point evaluator independently recovers the finite OPE limit.")

# Actual single-Higgs conjugation in the real lattice module:
# H_d is the adjoint of the physical up component, not the second (27,3) vertex.
def conjugated_down_higgs(family):
    p, s = matter_vertex(higgs_down,family)
    dual = tuple(-x for x in p)
    return dual, s*epsilon(dual,p)

for left,right,_ in slots["u"]:
    for i,j,k in product(range(3),repeat=3):
        p,_ = matter_vertex(left,i)
        q_,_ = matter_vertex(right,j)
        higgs,_ = conjugated_down_higgs(k)
        assert sum(p[5:])+sum(q_[5:])+sum(higgs[5:]) == 2
        assert vertex_response((p,q_,higgs),(3,2,1)) == 0

complement = [vector(*tail) for tail in sorted(set(permutations((0,-1,-1))))]

def complement_local(slot, root, second_terms=None):
    raw = raw_coefficient(slot)
    out = sp.zeros(3)
    for i,j,k in product(range(3),repeat=3):
        p,ps = matter_vertex(slot[0],i)
        r,rs = matter_vertex(slot[1],j)
        hc,hs = conjugated_down_higgs(k)
        inserts = [((root,),1)] if second_terms is None else [
            ((root,q),sgn) for q,sgn in second_terms]
        for prefix,sgn in inserts:
            pole1,finite1 = local_coefficients(prefix+(p,r,hc))
            pole2,finite2 = local_coefficients(prefix+(r,p,hc))
            assert pole1+pole2 == 0
            out[i,j] += sgn*h[k]*ps*rs*hs*(finite1+finite2)/(2*raw)
    return out

for root in complement:
    assert charge(root,root) == 2 and sum(root[5:]) == -2 and charge(root,em) == 0
first_color = [complement_local(slots["u"][0],root) for root in complement]
for slot in slots["u"]:
    assert [complement_local(slot,root) for root in complement] == first_color
print("True conjugate Higgs in the up slot has residual family weight +2; all SU(3)_F currents and EM carry 0.")
for root,matrix in zip(complement,first_color):
    print("Family complement",root[5:], "symmetric local matrix:",matrix,"det =",sp.factor(matrix.det()))

v = sp.Matrix(sp.symbols("v1:4",real=True))
assert (v*h.T+h*v.T).det() == 0
for i,component in enumerate((-sp.eye(3)[:,2],sp.eye(3)[:,1],sp.eye(3)[:,0])):
    assert first_color[i] == (component*h.T+h*component.T)/4
print("Any single complement current gives a symmetric outer-product sum of rank at most two.")
for root in complement:
    mixed_up = complement_local(slots["u"][0],root,hermitian_roots(beta))
    assert mixed_up == mixed_up.T
    for slot in slots["u"]:
        assert complement_local(slot,root,hermitian_roots(beta)) == mixed_up
    index = complement.index(root)
    expected_det = (-h[0]**2*h[1]/6144,h[0]**2*h[2]/6144,0)[index]
    assert sp.expand(mixed_up.det()-expected_det) == 0
    print("Complement plus native beta", root[5:], "gives Up", mixed_up, "det =",sp.factor(mixed_up.det()))

# Transport the already used old q marking as a lattice vector, retaining its orientation.
fold_file = contracts / "v16-operator-pullback-20261014/fold-results.json"
fold = sp.Matrix(json.loads(fold_file.read_text())["explicit_Weyl_matrix"])
assert fold*fold.T == sp.eye(8)
old_q = fold.T*(sp.eye(8)[:,7]-sp.eye(8)[:,6])
marked = tuple(fold*old_q)
assert marked == tuple(-x for x in beta)
oriented = ((marked,1),)
marked_e, _ = local_matrices(slots["e"][0],oriented)
assert marked_e.det() == 0
assert all(local_matrices(slot,oriented)[0] == marked_e for slot in slots["d"])
print("Old oriented family mark after the source fold:", marked[5:],
      "; electron/down determinant:",sp.factor(marked_e.det()))
for root in complement:
    marked_up = complement_local(slots["u"][0],root,oriented)
    assert marked_up.det() == 0
    assert all(complement_local(slot,root,oriented)==marked_up for slot in slots["u"])
    print("Old oriented mark with complement",root[5:],"Up determinant:",sp.factor(marked_up.det()))

# With unit coefficients for the candidate first current and its complement,
# the first real Hermitian branch fails the exact u/d determinant test for all real h.
primary_up = complement_local(slots["u"][0],complement[0],hermitian_roots(beta))
exact_signed_ratio = sp.cancel(-primary_up.det()/local_even.det())
assert sp.cancel(exact_signed_ratio - 6561*h[0]**2/
                 (321489*h[0]**2 + 20736*h[1]**2 + 409600*h[2]**2)) == 0
assert reduced_ud > 1/mp.mpf(49)
print("Equal-weight first Hermitian pair gives |Du/Dd| <= 1/49, below 0.661554.")

reduced_ed = mp.fprod(mp.mpf(TAB[f][1]) for f in ("e","mu","tau"))/mp.fprod(
    mp.mpf(TAB[f][1]) for f in ("d","s","b"))
assert mp.fprod(mp.mpf(TAB[f][1])-half_printing for f in ("e","mu","tau")) > mp.fprod(
    mp.mpf(TAB[f][1])+half_printing for f in ("d","s","b"))
print("Pure family dressing leaves E=D at fixed neutral Higgs; the old table needs |De/Dd| =",
      mp.nstr(reduced_ed,12), "after pure length factors.")

positive = ((beta,1),)
print("Opposite direction alone: electron/down det =",sp.factor(local_matrices(slots["e"][0],positive)[0].det()),
      "; Up det =",sp.factor(complement_local(slots["u"][0],complement[0],positive).det()))
# One repeated Hermitian source direction at the two marked positions.
x_beta,x_complement = sp.symbols("b g",real=True)
for root in complement:
    combined = tuple((p,x_beta*c) for p,c in hermitian_roots(beta)) + tuple(
        (p,x_complement*c) for p,c in hermitian_roots(root))
    single_e, _ = local_matrices(slots["e"][0],combined,combined)
    single_d, _ = local_matrices(slots["d"][0],combined,combined)
    single_up, _ = local_matrices(slots["u"][0],combined,combined,conjugate=True)
    assert single_e == single_d and sp.expand(single_up.det()) == 0
    assert (sp.expand(single_e.det()) == 0) == (root == complement[-1])
    print("Repeated direction with complement",root[5:],": det(E) =",sp.factor(single_e.det()),
          "; det(U) =",sp.factor(single_up.det()))

# Type the actual folded spin character separately from the affine quarter grade.
def old_spin(root):
    momentum = fold.T*sp.Matrix(root)
    return sp.simplify(sp.exp(2*sp.pi*sp.I*momentum[0]))

assert old_spin(alpha) == 1 and old_spin(beta) == -1
assert [old_spin(root) for root in complement] == [1,-1,-1]
relative_threshold = (49*reduced_ud)**(mp.mpf(1)/3)
assert relative_threshold > mp.mpf("3.18")
print("Old fold spin character: alpha and only the first complement are even; beta is odd.")
print("If the first mark is spin-even and positions stay fixed, |g/a| must exceed",
      mp.nstr(relative_threshold,12), "to pass the exact residual determinant.")
