"""Exact E8 Chevalley test of the A4 x A4 Petersen response.

This checks the actual current source, its cocycle signs, both A4 actions,
the 50 current insertions, and the normalized inclusive channel.  It does
not turn that inclusive algebraic response into a physical P1 instrument.
"""
from collections import Counter
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


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def positive(root):
    return next(x for x in root if x) > 0


def perm_sign(word):
    return -1 if sum(word[i] > word[j] for i in range(len(word))
                     for j in range(i+1, len(word))) % 2 else 1


def wedge_action(subset, u, v):
    """E_uv on an ordered exterior word."""
    if v not in subset or u in subset:
        return None
    word = list(subset)
    word[word.index(v)] = u
    return tuple(sorted(word)), perm_sign(word)


def carrier_weight(root):
    total = sum(root[:5])
    return tuple(5*x-total for x in root[:5])


def carrier_type(root):
    values = sorted(carrier_weight(root))
    table = {
        (-2, -2, -2, -2, 8): '5',
        (-8, 2, 2, 2, 2): 'bar5',
        (-4, -4, -4, 6, 6): '10',
        (-6, -6, 4, 4, 4): 'bar10',
        (0, 0, 0, 0, 0): 'zero',
    }
    return table.get(tuple(values), 'carrier_root')


def carrier_label(root, representation):
    weight = carrier_weight(root)
    if representation == '5':
        return weight.index(8)
    if representation == 'bar5':
        return weight.index(-8)
    if representation == '10':
        return tuple(i for i, x in enumerate(weight) if x == 6)
    if representation == 'bar10':
        return tuple(i for i, x in enumerate(weight) if x == -6)
    raise RuntimeError('unknown carrier representation')


def source_hash(path):
    return {
        'path': path,
        'sha256': hashlib.sha256((REPO/path).read_bytes()).hexdigest(),
    }


def solve_f2(equations, variables):
    pivots = {}
    for row, rhs in equations:
        while row:
            pivot = row.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = (row, rhs)
                break
            other, value = pivots[pivot]
            row ^= other
            rhs ^= value
        need(bool(row) or rhs == 0, 'joint cocycle equations consistent')
    solution = 0
    for pivot in sorted(pivots):
        row, rhs = pivots[pivot]
        value = rhs ^ ((row & solution).bit_count() % 2)
        if value:
            solution |= 1 << pivot
    need(all((row & solution).bit_count() % 2 == rhs
             for row, rhs in equations), 'joint cocycle solution exact')
    need(all(pivot < variables for pivot in pivots), 'phase pivots in range')
    return solution, pivots


def main():
    source = _actual_current_source()
    algebra = source['chevalley']
    roots = sorted(algebra.roots)
    root_set = set(roots)
    need(len(roots) == 240, 'E8 has240 roots')

    # The first A4 is the given SU5 inside D5.  Its root centralizer is
    # intrinsic: all first five doubled coordinates are equal.
    carrier_roots = [r for r in roots
                     if sum(x != 0 for x in r[:5]) == 2
                     and all(x == 0 for x in r[5:])
                     and sum(r[:5]) == 0]
    family_roots = [r for r in roots if len(set(r[:5])) == 1]
    need(len(carrier_roots) == len(family_roots) == 20,
         'given carrier A4 has an exact20-root A4 centralizer')

    # Derive the family simple roots from the actual centralizer.  Lexical
    # positivity fixes a chamber.  The chain orientation is then named so
    # that the input five is bar5, a representation convention rather than
    # a new embedding.
    family_positive = [r for r in family_roots if positive(r)]
    positive_set = set(family_positive)
    family_simple_unordered = [
        r for r in family_positive
        if not any(sub(r, a) in positive_set for a in family_positive)
    ]
    need(len(family_simple_unordered) == 4, 'centralizer has four simple roots')
    adjacency = {
        i: [j for j, b in enumerate(family_simple_unordered)
            if dot(a, b) == -4]
        for i, a in enumerate(family_simple_unordered)
    }
    need(sorted(len(v) for v in adjacency.values()) == [1, 1, 2, 2],
         'centralizer Dynkin graph is A4')

    endpoint_orders = []
    for start in sorted(i for i, v in adjacency.items() if len(v) == 1):
        order = []
        previous = None
        current = start
        while True:
            order.append(current)
            following = [j for j in adjacency[current] if j != previous]
            if not following:
                break
            previous, current = current, following[0]
        endpoint_orders.append([family_simple_unordered[i] for i in order])

    standard5 = (
        (1, 0, 0, 0),
        (-1, 1, 0, 0),
        (0, -1, 1, 0),
        (0, 0, -1, 1),
        (0, 0, 0, -1),
    )
    standard5_set = set(standard5)
    bar5_set = {tuple(-x for x in w) for w in standard5}

    sectors = {name: [r for r in roots if carrier_type(r) == name]
               for name in ('5', 'bar5', '10', 'bar10')}
    need(all(len(v) == 50 for v in sectors.values()),
         'four nonneutral A4xA4 sectors have50 roots each')

    def old_d5a3_sector(root):
        if sum(x != 0 for x in root) == 2:
            carrier_support = sum(x != 0 for x in root[:5])
            return 'g0' if carrier_support in (0, 2) else 'g2'
        need(all(abs(x) == 1 for x in root), 'remaining E8 root is halfroot')
        return 'g1' if sum(x > 0 for x in root[:5]) % 2 == 0 else 'g3'

    old_sector_overlap = {
        name: dict(sorted(Counter(old_d5a3_sector(r) for r in values).items()))
        for name, values in sectors.items()
    }
    need(old_sector_overlap == {
        '5': {'g2': 30, 'g3': 20},
        'bar5': {'g1': 20, 'g2': 30},
        '10': {'g0': 10, 'g1': 40},
        'bar10': {'g0': 10, 'g3': 40},
    }, 'new A4 sectors cross the old D5xA3 grading exactly')

    def family_charge(root, simple):
        return tuple(dot(root, a)//4 for a in simple)

    orientations = [
        simple for simple in endpoint_orders
        if {family_charge(r, simple) for r in sectors['10']} == bar5_set
    ]
    need(len(orientations) == 1, 'family chain orientation names input as bar5 uniquely')
    family_simple = orientations[0]
    cartan = [[dot(a, b)//4 for b in family_simple] for a in family_simple]
    need(cartan == [[2, -1, 0, 0], [-1, 2, -1, 0],
                    [0, -1, 2, -1], [0, 0, -1, 2]],
         'derived family Cartan is standard A4')

    std_index = {w: i for i, w in enumerate(standard5)}
    bar_index = {tuple(-x for x in w): i for i, w in enumerate(standard5)}
    pair_index = {
        tuple(standard5[a][k]+standard5[b][k] for k in range(4)): (a, b)
        for a in range(5) for b in range(a+1, 5)
    }

    def family_label(root, representation):
        charge = family_charge(root, family_simple)
        if representation == '5':
            return std_index[charge]
        if representation == 'bar5':
            return bar_index[charge]
        if representation == '10':
            return pair_index[charge]
        if representation == 'bar10':
            return pair_index[tuple(-x for x in charge)]
        raise RuntimeError('unknown family representation')

    inputs = sectors['10']
    insertions = sectors['5']
    outputs = sectors['bar10']
    modules = {
        'input': inputs,
        'insertion': insertions,
        'output': outputs,
    }
    state_module = {r: name for name, values in modules.items() for r in values}

    labels = {}
    subsets = {}
    all_five = set(range(5))
    for root in inputs:
        ij = carrier_label(root, '10')
        q = family_label(root, 'bar5')
        labels[root] = (ij, q)
        subsets[root] = (ij, tuple(sorted(all_five-{q})))
    for root in insertions:
        m = carrier_label(root, '5')
        ab = family_label(root, '10')
        labels[root] = (m, ab)
        subsets[root] = ((m,), ab)
    for root in outputs:
        kl = carrier_label(root, 'bar10')
        r = family_label(root, '5')
        labels[root] = (kl, r)
        subsets[root] = (tuple(sorted(all_five-set(kl))), (r,))
    need(all(len({labels[r] for r in values}) == 50 for values in modules.values()),
         'all three tensor weight charts are bijective')

    # Joint cocycle variables are phases on the 150 sector roots and on the
    # 40 neutral A4 roots.  No factorized state gauge is assumed.
    phase_roots = inputs + insertions + outputs + carrier_roots + family_roots
    need(len(set(phase_roots)) == 190, 'phase roots are disjoint')
    phase_index = {r: i for i, r in enumerate(phase_roots)}
    equations = []
    bracket_actions = {}

    for beta in insertions:
        m, ab = labels[beta]
        action = {}
        for alpha in inputs:
            ij, q = labels[alpha]
            candidate = add(alpha, beta)
            raw = algebra.bracket(algebra.ridx[beta], algebra.ridx[alpha])
            expected_nonzero = m not in ij and q in ab
            if not expected_nonzero:
                need(not raw and candidate not in outputs,
                     'zero bracket equals tensor support zero')
                continue
            need(candidate in outputs and set(raw) == {algebra.ridx[candidate]},
                 'nonzero bracket has the expected output root only')
            kl, family_out = labels[candidate]
            need(set(kl) == all_five-set(ij)-{m}, 'carrier complement target')
            need(family_out == (ab[1] if q == ab[0] else ab[0]),
                 'family antisymmetric target')
            coefficient = int(raw[algebra.ridx[candidate]])
            # In the covariance chart bar10 is represented directly as
            # Lambda^3, so the carrier sign is the wedge sign.  The input
            # bar5 is represented as Lambda^4; converting its complement
            # word to a covector contributes the displayed Hodge sign.
            carrier_sign = perm_sign((m,)+ij)
            family_word = tuple(sorted(all_five-{q}))
            hodge_sign = perm_sign(family_word+(q,))
            family_sign = hodge_sign*(1 if q == ab[0] else -1)
            expected_sign = carrier_sign*family_sign
            rhs = int(coefficient != expected_sign)
            row = ((1 << phase_index[alpha]) ^ (1 << phase_index[beta])
                   ^ (1 << phase_index[candidate]))
            equations.append((row, rhs))
            action[alpha] = (candidate, coefficient, expected_sign)
        bracket_actions[beta] = action

    bracket_equations = len(equations)
    need(bracket_equations == 600, '50 rank12 brackets give600 signed edges')

    # Add both neutral A4 actions on all three modules.  This proves that the
    # same phase chart transports the full representation action, including
    # Cartan through the root commutators.
    carrier_root_label = {}
    for beta in carrier_roots:
        plus = [i for i, x in enumerate(beta[:5]) if x == 2]
        minus = [i for i, x in enumerate(beta[:5]) if x == -2]
        need(len(plus) == len(minus) == 1, 'carrier root is E_uv')
        carrier_root_label[beta] = (plus[0], minus[0])

    standard_root_label = {
        tuple(standard5[u][k]-standard5[v][k] for k in range(4)): (u, v)
        for u in range(5) for v in range(5) if u != v
    }
    family_root_label = {
        beta: standard_root_label[family_charge(beta, family_simple)]
        for beta in family_roots
    }

    covariance_edges = 0
    for beta in carrier_roots + family_roots:
        is_carrier = beta in carrier_root_label
        u, v = (carrier_root_label if is_carrier else family_root_label)[beta]
        for alpha in inputs + insertions + outputs:
            carrier_subset, family_subset = subsets[alpha]
            moved = wedge_action(carrier_subset if is_carrier else family_subset, u, v)
            candidate = add(alpha, beta)
            raw = algebra.bracket(algebra.ridx[beta], algebra.ridx[alpha])
            if moved is None:
                need(not raw and candidate not in state_module,
                     'neutral A4 zero action equals exterior support')
                continue
            moved_subset, expected_sign = moved
            expected_subsets = ((moved_subset, family_subset) if is_carrier
                                else (carrier_subset, moved_subset))
            need(candidate in state_module
                 and state_module[candidate] == state_module[alpha],
                 'neutral A4 action stays in its module')
            need(subsets[candidate] == expected_subsets,
                 'neutral A4 action has exact tensor target')
            need(set(raw) == {algebra.ridx[candidate]},
                 'neutral A4 action has one actual target')
            coefficient = int(raw[algebra.ridx[candidate]])
            rhs = int(coefficient != expected_sign)
            row = ((1 << phase_index[alpha]) ^ (1 << phase_index[candidate])
                   ^ (1 << phase_index[beta]))
            equations.append((row, rhs))
            covariance_edges += 1
    need(covariance_edges == 1500, 'both A4 actions contribute1500 signed edges')

    solve_f2(equations[:bracket_equations], len(phase_roots))
    solve_f2(equations[bracket_equations:], len(phase_roots))
    solution, pivots = solve_f2(equations, len(phase_roots))

    input_gain = Counter()
    output_gain = Counter()
    ranks = {}
    signed_factor_edges = 0
    dual_coordinate_edges = 0
    xcharge_edges = Counter()
    X = (-2, -2, -2, 3, 3)
    for beta, action in bracket_actions.items():
        m, unused_ab = labels[beta]
        ranks[beta] = len(action)
        for alpha, (target, coefficient, expected_sign) in action.items():
            phase = (-1) ** (((solution >> phase_index[alpha])
                             ^ (solution >> phase_index[beta])
                             ^ (solution >> phase_index[target])) & 1)
            need(coefficient*phase == expected_sign,
                 'actual bracket equals C_m tensor A_ab after cocycle transport')
            # The covariance test uses Lambda^4 for the input bar5 and
            # Lambda^3 for the output bar10.  Convert both to the fixed dual
            # coordinates.  This must turn the wedge coefficient literally
            # into epsilon(m,i,j,k,l) times the antisymmetric A_ab entry.
            ij, q = labels[alpha]
            kl, unused_r = labels[target]
            input_family_word = tuple(sorted(all_five-{q}))
            h_q = perm_sign(input_family_word+(q,))
            output_carrier_word = tuple(sorted(all_five-set(kl)))
            h_kl = perm_sign(output_carrier_word+kl)
            converted = h_kl*coefficient*phase*h_q
            epsilon_sign = perm_sign((m,)+ij+kl)
            family_dual_sign = 1 if q == unused_ab[0] else -1
            need(converted == epsilon_sign*family_dual_sign,
                 'dual coordinate amplitude is literal epsilon tensor A_ab')
            dual_coordinate_edges += 1
            input_gain[alpha] += coefficient*coefficient
            output_gain[target] += coefficient*coefficient
            signed_factor_edges += 1
            qin = sum(X[i] for i in ij)
            qinsert = X[m]
            qout = -sum(X[k] for k in kl)
            need(qout == qin+qinsert, 'original P2 X charge is conserved')
            xcharge_edges[(qin, qinsert, qout)] += 1

    need(set(ranks.values()) == {12} and len(ranks) == 50,
         'all50 actual insertion maps have rank12')
    need(set(input_gain.values()) == {12} and len(input_gain) == 50,
         'sum Bdagger B is12 I50')
    need(set(output_gain.values()) == {12} and len(output_gain) == 50,
         'sum B Bdagger is12 I50')
    need(signed_factor_edges == 600, 'all600 signed factors verified')
    need(dual_coordinate_edges == 600,
         'all600 fixed dual-coordinate amplitudes verified')

    # Independent finite channel counts: C_m/sqrt3 is Petersen, while
    # A_ab/2 is the complete-graph wedge/contraction channel on five states.
    carrier_population = {}
    for ij in sorted({carrier_label(r, '10') for r in inputs}):
        targets = Counter()
        for m in range(5):
            if m not in ij:
                kl = tuple(sorted(all_five-set(ij)-{m}))
                targets[kl] += 1
        need(len(targets) == 3 and set(targets.values()) == {1},
             'Hodge channel has three Petersen neighbours')
        need(all(set(ij).isdisjoint(kl) for kl in targets),
             'Hodge population is Petersen disjointness')
        carrier_population[','.join(map(str, ij))] = [','.join(map(str, kl))
                                                     for kl in sorted(targets)]
    family_population = {}
    for q in range(5):
        targets = Counter()
        for a in range(5):
            for b in range(a+1, 5):
                if q in (a, b):
                    targets[b if q == a else a] += 1
        need(len(targets) == 4 and set(targets.values()) == {1},
             'family wedge channel has four neighbours')
        family_population[str(q)] = sorted(targets)

    result = {
        'scope': ('Actual weight-one E8 current response in the given carrier '
                  'SU5 frame; no physical source selection or time law.'),
        'root_partition': {
            'carrier_A4_roots': len(carrier_roots),
            'centralizer_A4_roots': len(family_roots),
            'input_(10,bar5)': len(inputs),
            'insertion_(5,10)': len(insertions),
            'output_(bar10,5)': len(outputs),
            'conjugate_sector_(bar5,bar10)': len(sectors['bar5']),
        },
        'overlap_with_old_D5xA3_grading': old_sector_overlap,
        'family_simple_roots': [list(r) for r in family_simple],
        'family_cartan': cartan,
        'bracket_signed_edges': bracket_equations,
        'literal_epsilon_dual_coordinate_edges': dual_coordinate_edges,
        'neutral_A4xA4_signed_edges': covariance_edges,
        'phase_variables': len(phase_roots),
        'phase_rank': len(pivots),
        'kraus_count': len(insertions),
        'kraus_rank': 12,
        'source_gain': '12 I50',
        'source_response_denominator': 12,
        'factorization': ('E8 response = Hodge Petersen channel tensor '
                          'family wedge channel in the joint cocycle chart'),
        'carrier_channel': '1/3 sum_m C_m rho C_m*; population A_Petersen/3',
        'family_channel': '1/4 sum_ab A_ab rho A_ab*; population (J-I)/4 on five',
        'p2_X': list(X),
        'p2_X_convention': ('q_input=X_i+X_j, q_insertion=X_m, '
                            'q_output=-(X_k+X_l); equivalently '
                            'q_output+q_environment=q_input with '
                            'q_environment=-X_m'),
        'p2_X_charge_edge_counts': {
            ','.join(map(str, key)): value
            for key, value in sorted(xcharge_edges.items())
        },
        'carrier_population': carrier_population,
        'family_population': family_population,
        'chart_caveat': ('The diagonal cocycle transport is a joint tensor '
                         'chart; no product phase gauge or raw-coordinate '
                         'partial trace is asserted.'),
        'dual_coordinate_conversion': {
            'input_bar5': 'h_q = sign(complement(q),q)',
            'output_bar10': 'h_kl = sign(complement(kl),k,l)',
            'verified_amplitude': 'h_kl * B_wedge * h_q = epsilon(m,i,j,k,l) A_ab',
        },
        'state_phases': {
            ','.join(map(str, r)): str((-1)**((solution >> phase_index[r]) & 1))
            for r in inputs + insertions + outputs
        },
        'input_bar5_phase_chart': {
            ','.join(map(str, r)): str(
                (-1)**((solution >> phase_index[r]) & 1)
                * perm_sign(tuple(sorted(all_five-{labels[r][1]}))+(labels[r][1],))
            )
            for r in inputs
        },
        'insertion_phase_chart': {
            ','.join(map(str, r)): str((-1)**((solution >> phase_index[r]) & 1))
            for r in insertions
        },
        'output_bar10_phase_chart': {
            ','.join(map(str, r)): str(
                (-1)**((solution >> phase_index[r]) & 1)
                * perm_sign(tuple(sorted(all_five-set(labels[r][0])))+labels[r][0])
            )
            for r in outputs
        },
        'neutral_root_phases': {
            ','.join(map(str, r)): str((-1)**((solution >> phase_index[r]) & 1))
            for r in carrier_roots + family_roots
        },
        'sources': [
            source_hash('tfpt_explorer/current_block_geometry.py'),
            source_hash('verification/v498_celestial_wp5b_singular_vector.py'),
        ],
    }
    output_path = Path(__file__).with_name('e8-a4-results.json')
    output_path.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('state_phases', 'input_bar5_phase_chart',
                                   'insertion_phase_chart',
                                   'output_bar10_phase_chart',
                                   'neutral_root_phases',
                                   'carrier_population', 'family_population')},
                     indent=2))


if __name__ == '__main__':
    main()
