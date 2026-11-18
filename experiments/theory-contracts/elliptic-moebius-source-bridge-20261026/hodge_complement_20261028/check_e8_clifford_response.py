"""Compare the actual E8 current bracket to the complete CAR response.

No stochastic time law or source insertion is selected by this check.
All source D8 root actions, including their Chevalley signs, are compared.
"""
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0, str(REPO))
from tfpt_explorer.current_block_geometry import _actual_current_source


def need(condition, name):
    if not condition:
        raise RuntimeError(name)


def ladder(mask, i, create):
    occupied = (mask >> i) & 1
    if occupied == create:
        return None
    sign = -1 if (mask & ((1 << i)-1)).bit_count() % 2 else 1
    return mask ^ (1 << i), sign


def bilinear(mask, root):
    axes = [i for i, x in enumerate(root) if x]
    i, j = axes
    first = ladder(mask, j, root[j] > 0)
    if first is None:
        return None
    second = ladder(first[0], i, root[i] > 0)
    if second is None:
        return None
    return second[0], first[1] * second[1]


def source_hash(path):
    return {'path': path, 'sha256': hashlib.sha256((REPO/path).read_bytes()).hexdigest()}


def main():
    source = _actual_current_source()
    algebra = source['chevalley']
    roots = list(algebra.roots)
    halves = sorted(r for r in roots if all(abs(x) == 1 for x in r))
    d8 = sorted(r for r in roots if sum(x != 0 for x in r) == 2)
    cross = [r for r in d8 if sum(x != 0 for x in r[:5]) == 1]
    masks = {r: sum((x > 0) << i for i, x in enumerate(r)) for r in halves}
    need(len(halves) == 128 and len(d8) == 112 and len(cross) == 60,
         'root sectors 128,112,60')
    need(all(mask.bit_count() % 2 == 0 for mask in masks.values()), 'even eight-Fock parent')

    # A diagonal weight-space phase is the only allowed state-basis freedom.
    # Root phases are harmless orthonormal Kraus-coordinate choices.
    hindex = {r: i for i, r in enumerate(halves)}
    bindex = {r: len(halves)+i for i, r in enumerate(d8)}
    nvars = len(halves)+len(d8)
    equations = []
    actions = {}
    for beta in d8:
        action = {}
        for alpha in halves:
            raw = algebra.bracket(algebra.ridx[beta], algebra.ridx[alpha])
            candidate = tuple(a+b for a, b in zip(alpha, beta))
            expected = bilinear(masks[alpha], beta)
            if candidate not in hindex:
                need(not raw and expected is None, 'zero action matches actual source')
                continue
            need(set(raw) == {algebra.ridx[candidate]}, 'actual root action has no omitted component')
            coefficient = raw[algebra.ridx[candidate]]
            need(coefficient in (1, -1) and expected is not None, 'unit coefficient and CAR support')
            need(expected[0] == masks[candidate], 'all Cartan charges match')
            rhs = int(coefficient != expected[1])
            row = (1 << hindex[alpha]) ^ (1 << hindex[candidate]) ^ (1 << bindex[beta])
            equations.append((row, rhs))
            action[alpha] = (candidate, int(coefficient))
        actions[beta] = action

    pivots = {}
    for row, rhs in equations:
        while row:
            pivot = row.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = (row, rhs)
                break
            other, val = pivots[pivot]
            row ^= other
            rhs ^= val
        need(bool(row) or rhs == 0, 'source/CAR phase equations consistent')
    solution = 0
    for pivot in sorted(pivots):
        row, rhs = pivots[pivot]
        value = rhs ^ ((row & solution).bit_count() % 2)
        if value:
            solution |= 1 << pivot
    need(all((row & solution).bit_count() % 2 == rhs for row, rhs in equations),
         'complete source-root intertwiner including all cocycle signs')

    carrier_parity = lambda r: (masks[r] & 31).bit_count() % 2
    g1 = [r for r in halves if carrier_parity(r) == 0]
    g3 = [r for r in halves if carrier_parity(r) == 1]
    need(len(g1) == len(g3) == 64, 'two full source modules64+64')
    input_gain = Counter()
    output_gain = Counter()
    for beta in cross:
        for alpha, (target, coef) in actions[beta].items():
            if alpha in g1:
                need(target in g3, 'g2 maps entire g1 into g3')
                input_gain[alpha] += coef*coef
                output_gain[target] += coef*coef
    need(set(input_gain.values()) == {15} and len(input_gain) == 64, 'sum B*B=15I64')
    need(set(output_gain.values()) == {15} and len(output_gain) == 64, 'sum BB*=15I64')

    # Compare every source Kraus column, with its solved phases, to the
    # tensor product of the actual carrier and family CAR ladders.
    factor_columns = 0
    for beta in cross:
        i = next(i for i in range(5) if beta[i])
        j = next(j for j in range(5, 8) if beta[j])
        for alpha in g1:
            cm = masks[alpha] & 31
            fm = masks[alpha] >> 5
            a = ladder(cm, i, beta[i] > 0)
            b = ladder(fm, j-5, beta[j] > 0)
            actual = actions[beta].get(alpha)
            if a is None or b is None:
                need(actual is None, 'tensor ladder has exact zero support')
                continue
            target, coef = actual
            need(masks[target] == a[0]+(b[0] << 5), 'tensor target charges')
            phase = (-1) ** (((solution >> hindex[alpha]) ^
                             (solution >> hindex[target]) ^
                             (solution >> bindex[beta])) & 1)
            need(coef == phase*a[1]*b[1], 'entire signed Kraus factorization')
            factor_columns += 1
    need(factor_columns == 960, '60 Kraus maps with16 source columns each')

    # The complete carrier response keeps both chiralities and all sectors.
    even = [m for m in range(32) if m.bit_count() % 2 == 0]
    odd = [m for m in range(32) if m.bit_count() % 2]
    population = {}
    for m in even+odd:
        targets = Counter()
        for i in range(5):
            for create in (False, True):
                out = ladder(m, i, create)
                if out:
                    targets[out[0]] += F(out[1]*out[1], 5)
        need(sum(targets.values()) == 1, 'complete Clifford carrier channel is TP')
        population[m] = targets
    shell_matrix = [[F(0) for _ in range(3)] for _ in range(3)]
    for m in even:
        shell_in = {0:0,4:1,2:2}[m.bit_count()]
        folded = Counter({31 ^ out: val for out, val in population[m].items()})
        bins = Counter()
        for out, val in folded.items():
            bins[{0:0,4:1,2:2}[out.bit_count()]] += val
        need(sum(bins.values()) == 1, 'folded Clebsch step TP')
        if not any(shell_matrix[shell_in]):
            shell_matrix[shell_in] = [bins[j] for j in range(3)]
        need(shell_matrix[shell_in] == [bins[j] for j in range(3)], 'exact shell lumping')
        if m.bit_count() == 2:
            jump = {out:val for out,val in folded.items() if out.bit_count() == 2}
            carry = {out:val for out,val in folded.items() if out.bit_count() == 4}
            need(len(jump) == 3 and sum(jump.values()) == F(3,5), 'Petersen jump3/5')
            need(all(not (m & out) and val == F(1,5) for out,val in jump.items()),
                 'Petersen exact disjoint matrix A/5')
            need(len(carry) == 2 and sum(carry.values()) == F(2,5), 'genuine carry2/5 retained')

    result = {
        'scope':'Actual weight-one E8 current response under the given D5+A3 frame; no physical event/time selection.',
        'source_modules': {'g1':len(g1),'g2_root_kraus':len(cross),'g3':len(g3)},
        'all_D8_root_intertwiner_equations':len(equations),
        'phase_rank':len(pivots),'phase_variables':nvars,
        'signed_tensor_columns_checked':factor_columns,
        'source_gain':'15 I64','source_response_denominator':15,
        'carrier_channel':'1/5 sum_i (a_i rho a_i* + a_i* rho a_i)',
        'family_channel':'1/3 sum_j (b_j rho b_j* + b_j* rho b_j)',
        'factorization':'E8_g2_response = carrier_channel tensor family_channel in the full signed D8 source intertwiner',
        'carrier_population':'Q5 on32 parity-alternating words; Clebsch16 after antipodal complement labels',
        'Clebsch_shell_order':['0','weight4','weight2'],
        'Clebsch_shell_row_matrix':[[str(x) for x in row] for row in shell_matrix],
        'Petersen_jump':'A/5 with probability3/5',
        'retained_carry':'Lambda1 before complement label; five states, weight2/5; not reset2I/5',
        'state_phases':{','.join(map(str,r)):str((-1)**((solution>>hindex[r])&1)) for r in halves},
        'root_phases':{','.join(map(str,r)):str((-1)**((solution>>bindex[r])&1)) for r in d8},
        'sources':[source_hash('tfpt_explorer/current_block_geometry.py'),
                   source_hash('verification/v498_celestial_wp5b_singular_vector.py'),
                   source_hash('verification/v810_cutcode_moment.py')],
    }
    Path(__file__).with_name('e8-clifford-response-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('state_phases','root_phases')},indent=2))


if __name__ == '__main__':
    main()
