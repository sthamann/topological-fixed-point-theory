"""Independent audit of pasted folding proposal; preserves source holonomy."""
import argparse
import itertools
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import sympy as S

parser = argparse.ArgumentParser()
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
checks = 0


def ck(ok):
    global checks
    if not ok:
        raise AssertionError(checks)
    checks += 1


roots = []
for i, j in itertools.combinations(range(8), 2):
    for a, b in itertools.product((-1, 1), repeat=2):
        r = [F(0)]*8
        r[i], r[j] = F(a), F(b)
        roots.append(tuple(r))
roots += [tuple(F(x, 2) for x in r) for r in itertools.product((-1, 1), repeat=8)
          if sum(x < 0 for x in r) % 2 == 0]
rootset = set(roots)
dot = lambda x, y: sum(a*b for a, b in zip(x, y))


def perm(order):
    return S.Matrix([[int(j == order[i]) for j in range(8)] for i in range(8)])


hh = S.Matrix([[1, 1, 1, 1], [1, 1, -1, -1],
               [1, -1, 1, -1], [1, -1, -1, 1]])/2
pp = perm([0, 4, 2, 6, 1, 5, 3, 7])
qq = perm([2, 3, 4, 6, 7, 0, 1, 5])
mm = qq*S.diag(hh, hh)*pp
ck(mm.T*mm == S.eye(8) and mm.det() == 1)
converted = {tuple(F(x) for x in mm*S.Matrix(r)) for r in roots}
ck(converted == rootset)
nu = S.Matrix([S.Rational(k, 4) for k in range(8)])
target = S.Matrix([S.Rational(1, 2)]*5+[0]*3)
lattice_shift = mm*nu-target
ck(all(F(x).denominator == 2 for x in lattice_shift) and sum(lattice_shift) % 2 == 0)
mixed_cartan_root = S.Matrix([0, 0, 1, 0, 0, 1, 0, 0])
transported_charge = mm*S.ones(8, 1)
ck(transported_charge == 2*mixed_cartan_root)
ck(tuple(map(F, mixed_cartan_root)) in rootset)
for r in roots:
    mr = tuple(F(x) for x in mm*S.Matrix(r))
    ck((dot(tuple(map(F, nu)), r)-dot(tuple(map(F, target)), mr)).denominator == 1)


def fixed_components(torus):
    fixed = {r for r in roots if dot(torus, r).denominator == 1}
    components = []
    while fixed:
        todo = [fixed.pop()]
        comp = set(todo)
        while todo:
            r = todo.pop()
            neighbors = {t for t in fixed if dot(r, t)}
            fixed -= neighbors
            comp |= neighbors
            todo.extend(neighbors)
        rank = S.Matrix(list(comp)).rank()
        integer = sum(all(x.denominator == 1 for x in r) for r in comp)
        components.append({'roots':len(comp), 'rank':rank,
                           'integer_roots':integer, 'spinor_roots':len(comp)-integer})
    return sorted(components, key=lambda c:c['roots'])


def root_order(torus):
    result = 1
    for r in roots:
        result = S.ilcm(result, dot(torus, r).denominator)
    return int(result)


clocks = {}
# Independent exact check of the inherited Virasoro central term, including
# its required ground shift. c'=8 follows from the original c=1 bracket.
for m in range(-9,10):
    central_from_cover = F((8*m)**3-8*m,12*64)-2*m*F(21,64)
    ck(central_from_cover == F(8*(m**3-m),12))
for a in (F(0), F(1, 2), F(1, 4), F(3, 4)):
    mu = tuple((k+a)/8 for k in range(8))
    clock = tuple(2*x for x in mu)
    deck = tuple(4*x for x in mu)
    transported_clock = mm*S.Matrix(clock)
    source_charge_term = S.Rational(a.numerator, a.denominator)/4*transported_charge
    ck(transported_clock == target+lattice_shift+source_charge_term)
    for r in roots:
        mr = tuple(F(x) for x in mm*S.Matrix(r))
        expected = tuple(map(F, target+source_charge_term))
        ck((dot(clock, r)-dot(expected, mr)).denominator == 1)
    comps = fixed_components(clock)
    grades = Counter(str(dot(clock, r) % 1) for r in roots)
    clocks[str(a)] = {'clock_order_on_E8':root_order(clock),
                     'clock_fixed_components':comps,
                     'deck_fixed_components':fixed_components(deck),
                     'clock_phase_counts':dict(sorted(grades.items())),
                     'folded_ground_weight':str(sum((x-F(1,2))**2/2 for x in mu))}
    ck(sum((x-F(1,2))**2/2 for x in mu)
       == F(21,64)+(a-F(1,2))**2/16)
    if a == 0:
        ck([c['roots'] for c in comps] == [12,40])
        ck([c['rank'] for c in comps] == [3,5])
        ck([(c['integer_roots'],c['spinor_roots']) for c in comps] == [(4,8),(16,24)])
        ck(root_order(clock) == 4)
        ck([c['roots'] for c in fixed_components(deck)] == [112])
        ck([grades.get(str(F(k,4)),0) for k in range(4)] == [52,64,60,64])
    if a in (F(1,4), F(3,4)):
        ck(root_order(clock) != 4)

# Exact single-channel lattice fold: bijective reindexing, not extra copies.
# U shifts x by1 and the last link carries the original boundary phase.
# Fold x=j*L+r into (r,j). Its block-cut transition has g^8=phase I.
for length in (2,3,5):
    n = 8*length
    for boundary in (1, -1, 1j, -1j):
        def shift(x):
            return (x+1) % n, boundary if x == n-1 else 1
        def folded(x):
            j,r = divmod(x,length)
            return 8*r+j
        ck(len({folded(x) for x in range(n)}) == n)
        for x in range(n):
            y,p = shift(x)
            j,r = divmod(x,length)
            if r+1 < length:
                yy,pph = 8*(r+1)+j,1
            else:
                yy,pph = (j+1)%8, boundary if j==7 else 1
            ck(folded(y)==yy and p==pph)
        for j in range(8):
            result, ph = j*length, 1
            for _ in range(n):
                result, p = shift(result)
                ph *= p
            ck(result==j*length and ph==boundary)

output = {'research_id':'UR.SEAM.V16_OPERATOR_PULLBACK.01','check_verdict':'PASS',
          'checks':checks,'explicit_Weyl_matrix':[[str(x) for x in row] for row in mm.tolist()],
          'E8_lattice_shift':[str(x) for x in lattice_shift], 'spin_and_holonomy_cases':clocks,
          'transported_source_U1_charge':[str(x) for x in transported_charge],
          'source_clock_deformation':'M nu_alpha = t + lattice_shift + (alpha/2)(e3+e6)',
          'charge_scope':'Cartan charge transported within the stipulated folded E8 extension; not a raw-source field map or selected compensation',
          'fold_source_copies_added':0,
          'operator_relation':'F h_cover F^-1 = 8(-i d/dx) with cyclic boundary holonomy',
          'state_relation':'C_fold=F C_cover F*, generally not an eight-copy product vacuum',
          'quotient_Virasoro':'Lhat_n=L_(8n)/8+21/64 delta_n0, central charge8',
          'locality_relation':'A base interval lifts to eight separated cover arcs',
          'complete_TFPT_solution':False,'raw_Moebius_quotient_selected':False,
          'E8_extension_from_raw_source':False,
          'scope':'Exact periodic root-lattice conjugacy and finite-cover operator fold; raw source, locality, Pin and extension assumptions remain explicit.'}
args.out.parent.mkdir(parents=True,exist_ok=True)
args.out.write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({'checks':checks,'check_verdict':'PASS',
                  'clock_orders':{a:x['clock_order_on_E8'] for a,x in clocks.items()}}))
