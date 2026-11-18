"""Calculate positive charged Higgs transitions using the original E8 root OPE."""
import json
from itertools import combinations
from pathlib import Path
import sys

import sympy as s


REPO = Path(sys.argv[1])
CONTRACTS = REPO / 'experiments/theory-contracts'
sys.path.insert(0, str(CONTRACTS / 'seam-collar-multifield-20261016'))
from e8_boundary_response import add, adjoint_sign, dot, epsilon, neg, vertex_response

HALF = s.Rational(1, 2)
FAMILIES = [tuple(HALF*i for i in row) for row in [(-1,1,1),(1,-1,1),(1,1,-1),(-1,-1,-1)]]


def family_vertex(root, index):
    shift = (0,)*5 + tuple(i-j for i,j in zip(FAMILIES[index], FAMILIES[0]))
    return add(root, shift), epsilon(shift, root)


def charged_slots():
    original = CONTRACTS / 'compiler-vacuum-current-cubic-20260918/certificate.json'
    cubic = json.loads(original.read_text())['cubic']
    base = [tuple(HALF*v for v in row) for row in cubic['base_roots_doubled']]
    def particle(root):
        if not all(abs(v) == HALF for v in root[:5]):
            return None
        color = sum(v+HALF for v in root[:3])
        weak = sum(v+HALF for v in root[3:5])
        return {(1,1):'Q',(2,0):'u',(0,2):'e',(2,2):'d',(3,1):'L'}.get((color,weak))
    result = {tag:[] for tag in ('e','d','u')}
    for row in cubic['monomials']:
        for tag, left, sign in [('e','L',-1),('d','Q',-1),('u','Q',1)]:
            roots = [base[i] for i in row['indices']]
            if not any(p[:5] == (0,0,0,0,sign) for p in roots):
                continue
            l = [p for p in roots if particle(p) == left]
            r = [p for p in roots if particle(p) == tag]
            if l and r:
                result[tag].append((l[0],r[0]))
    assert [len(result[t]) for t in ('e','d','u')] == [1,3,3]
    return result


def transition(slot, terms):
    left, right = slot
    terms = dict(terms)
    result = s.zeros(4)
    checked = 0
    for i in range(4):
        p, p_sign = family_vertex(left,i)
        for j in range(4):
            q, q_sign = family_vertex(right,j)
            h = neg(add(p,q))
            if h not in terms:
                continue
            assert dot(h,q) == -1
            # bra corresponding to the adjoint of the marked left state
            sign = p_sign*q_sign*adjoint_sign(p)
            amplitude = sign*epsilon(h,q)
            direct = sign*2*vertex_response((h,q,p),(3,2,1))/vertex_response((neg(p),p),(2,1))
            assert s.simplify(direct-amplitude) == 0
            result[i,j] = terms[h]*amplitude
            checked += 1
    assert checked == 12
    return result


def hodge_conjugate(form):
    star = s.zeros(4)
    for i in range(4):
        for j in range(4):
            if i != j:
                k,l = sorted(set(range(4))-{i,j})
                star[i,j] = s.LeviCivita(i,j,k,l)*s.conjugate(form[k,l])
    return star


def norm(matrix):
    return s.simplify(matrix.H*matrix)


def positive_source_test(form, h):
    coordinates = [form.subs(dict(zip(h,[int(i==j) for i in range(6)]))) for j in range(6)]
    supports = {tuple(i for i in range(4) if norm(a)[i,i]): norm(a) for a in coordinates}
    assert len(supports) == 6
    assert sum(supports.values(),s.zeros(4))/6 == s.eye(4)/2
    pure_a2 = sum((p for ij,p in supports.items() if 3 not in ij),s.zeros(4))/3
    assert pure_a2 == s.diag(s.Rational(2,3),s.Rational(2,3),s.Rational(2,3),0)

    # Candidate original minimal 2+2 clutch grade, plus the existing 3+1 mark.
    # Invariance under the remaining family transposition leaves two equal
    # diagonal responses on the retained triplet in a grade-diagonal model.
    weights = {}
    for positive_pair in combinations(range(4),2):
        beta = [HALF if i in positive_pair else -HALF for i in range(4)]
        for omitted in range(4):
            weights.clear()
            k = s.zeros(4)
            for i,j in supports:
                labels = (beta[i]+beta[j],int(omitted in (i,j)))
                if labels not in weights:
                    weights[labels] = s.Symbol(f'w{len(weights)}')
                k += weights[labels]*supports[i,j]
            assert len({k[i,i] for i in range(4) if i!=omitted}) <= 2

    x = s.Symbol('x',positive=True)
    target = [s.Rational(55,117)**2, (s.Rational(34,47)*(1-x))**2,
              (s.Rational(3,26)*(1-x)/x)**2]
    # Each pair product differs from 1 for the transcendental physical x.
    assert all(s.factor(target[i]*target[j]-1)!=0 for i in range(3) for j in range(i+1,3))
    seed = 1/(6*s.pi)+3/(256*s.pi**4)
    strength = 1/target[0]-1
    seed_lo = s.Rational(5317,100000)
    seed_hi = s.Rational(2659,50000)
    seed_formula = lambda pi: 1/(6*pi)+3/(256*pi**4)
    assert seed_lo < seed_formula(s.Rational(314160,100000))
    assert seed_formula(s.Rational(314159,100000)) < seed_hi
    assert strength > (target[2]-1).subs(x,seed_lo)
    assert strength > s.Rational(1,4)

    # At minimum common strength kappa_1=0. The trace-two positive
    # extension can be written using the three complement root planes.
    kappa = [s.factor(((1+strength)*r-1)/(strength*(1+r))) for r in target]
    assert kappa[0]==0
    k2,k3 = kappa[1:]
    mixture = {(1,2):k2+k3-1,(1,3):1-k3,(2,3):1-k2}
    assert s.simplify(sum(mixture.values())-1)==0
    assert (k2+k3-1).subs(x,seed_hi) > 0
    assert (1-k2).subs(x,seed_lo)>0 and (1-k3).subs(x,seed_lo)>0
    mixed = sum((value*supports[key] for key,value in mixture.items()),s.zeros(4)).applyfunc(s.factor)
    assert s.simplify(mixed-s.diag(*kappa, 2-sum(kappa))) == s.zeros(4)
    down = s.eye(3)+strength*mixed[:3,:3]
    up = s.eye(3)+strength*(s.eye(3)-mixed[:3,:3])
    assert s.simplify(down*up.inv()-s.diag(*target)) == s.zeros(3)
    # The unobserved fourth diagonal is bounded for every allowed strength.
    limiting_observed_sum = sum(r/(1+r) for r in target)
    assert limiting_observed_sum.subs(x,seed_lo) < s.Rational(1309,1000)
    hidden_floor = 2-limiting_observed_sum
    # At minimum strength K_11=0 forces every contributing form into 234.
    # In that 3D space K=I-density_of_Hodge_normals, whose rank is 3 here.
    covariance_normals = s.eye(3)-mixed[1:,1:]
    assert all(float(v.subs(x,seed))>0 for v in covariance_normals.diagonal())
    # Bound the determinant for any positive mixture, even if the bare
    # Higgs vertex is not diagonal. Increasing kappa to sum 2 increases
    # the ratio. Write m_i=1-kappa_i and eliminate m_3 using sum(m)=1.
    t,m,n=s.symbols('t m n',real=True)
    holes=(m,n,1-m-n)
    denominator=s.prod(1+t*val for val in holes)
    numerator=s.prod(1+t*(1-val) for val in holes)
    assert s.expand((1+t)*denominator-numerator-(t+2)*t**3*s.prod(holes))==0
    ordered_equal_weights_strength_floor = 49**2*s.prod(target)-1
    # Uniform cubic-vacuum root Gram and the known equal root average carry
    # no such norm contrast; this mixture is a conditional feasibility test.
    return {'target': [float(r.subs(x,seed)) for r in target],
            'target_products': [float((target[i]*target[j]).subs(x,seed))
                                for i,j in [(0,1),(0,2),(1,2)]],
            'clutch_grade_diagonal_checks': 24,
            'minimum_common_strength': float(strength),
            'fourth_direction_strict_floor': float(hidden_floor.subs(x,seed)),
            'minimum_covariance_rank_at_threshold': 3,
            'ordered_equal_weights_strength_floor': float(ordered_equal_weights_strength_floor.subs(x,seed)),
            'positive_mixture': {str(tuple(i+1 for i in key)):float(v.subs(x,seed)) for key,v in mixture.items()},
            'complement_orbit_occupancies': [float(a.subs(x,seed)) for a in mixed.diagonal()]}


def main():
    h = s.symbols('h1:7')
    down = []
    for j in range(3):
        for sign in (-1,1):
            down.append(tuple([0,0,0,0,-1]+[sign if k==j else 0 for k in range(3)]))
    terms = tuple(zip(down,h))
    adjoint = [(neg(a),adjoint_sign(a)*s.conjugate(value)) for a,value in terms]
    slots = charged_slots()
    maps = {tag:[transition(slot,adjoint if tag=='u' else terms) for slot in rows] for tag,rows in slots.items()}
    reference = maps['e'][0]
    dual = hodge_conjugate(reference)
    assert reference.T == -reference and dual.T == -dual
    assert all(a in (reference,-reference) for a in maps['d'])
    assert all(a in (dual,-dual) for a in maps['u'])
    mass = sum(a*s.conjugate(a) for a in h)
    assert s.simplify(norm(reference)+norm(dual)-mass*s.eye(4)) == s.zeros(4)
    assert s.Matrix(6,6,lambda i,j: adjoint_sign(down[i])*
                    vertex_response((neg(down[i]),down[j]),(2,1))) == s.eye(6)
    for rows in slots.values():
        for _,right in rows:
            states = [family_vertex(right,k)[0] for k in range(4)]
            assert s.Matrix(4,4,lambda i,j: adjoint_sign(states[i])*
                            vertex_response((neg(states[i]),states[j]),(2,1))) == s.eye(4)
    result = {'charged_channels':sum(map(len,slots.values())), 'nonzero_original_OPE_checks':84,
              'representative_down_and_e':str(reference), 'representative_up':str(dual),
              'vacuum_right_gram': 'I4 for all 7 original slots',
              'complement': 'A^H A + (*conj A)^H (*conj A) = sum(abs(h)^2) I4'}
    result.update(positive_source_test(reference,h))
    # Try the existing neutral charge-weighted current as the sole carrier
    # correction to the fixed-Higgs complementary norm. Its zero mode on
    # a root primary is p^T R_Y p/2, with no oscillator contribution.
    sys.path.insert(0,str(REPO))
    from tfpt_explorer.source_charge_response import exact_charge_response
    charge_operator=exact_charge_response()['R_Y']
    shifts={}
    for tag, rows in slots.items():
        leg_values=[]
        for leg in (0,1):
            values={ (p.T*charge_operator*p)[0]/2 for slot in rows for i in range(3)
                     for p in [s.Matrix(family_vertex(slot[leg],i)[0])] }
            assert len(values)==1
            leg_values.append(values.pop())
        shifts[tag]=leg_values
    assert {tag: vals[1] for tag,vals in shifts.items()} == dict(e=s.Rational(31,12),u=s.Rational(11,12),d=-s.Rational(1,12))
    rounded=s.Rational(1,2000)
    leptons=list(map(s.Rational,['0.466','1.085','0.899']))
    downs=list(map(s.Rational,['0.936','0.943','0.480']))
    archived_interval=[s.prod(v-rounded for v in leptons)/s.prod(v+rounded for v in downs),
                       s.prod(v+rounded for v in leptons)/s.prod(v-rounded for v in downs)]
    ratios={}
    seed_bounds=(s.Rational(5317,100000),s.Rational(2659,50000))
    def target_at(x):
        return (s.Rational(55,117)**2,(s.Rational(34,47)*(1-x))**2,
                (s.Rational(3,26)*(1-x)/x)**2)
    for legs in [(1,),(0,1)]:
        charges={tag:sum(vals[i] for i in legs) for tag,vals in shifts.items()}
        assert charges['d']-charges['u']==-1
        power=3*(charges['e']-charges['d'])/4
        ratios[str(legs)]=[pair_product**float(power) for pair_product in result['target_products']]
        for i,j in combinations(range(3),2):
            upper,lower=[(r[i]*r[j])**int(2*power) for r in map(target_at,seed_bounds)]
            assert upper < archived_interval[0]**2 or lower > archived_interval[1]**2
    result['charge_weighted_original_zero_mode']={tag:[str(v) for v in vals] for tag,vals in shifts.items()}
    result['single_charge_exponential_archived_ed_tests']=ratios
    result['archived_ed_interval']=[float(v) for v in archived_interval]
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
