"""Exact OriginLab-style campaign for coherent OPH histories and charge transport.

Run with Python and SymPy: python3 -B originlab_oph_f1_search.py [--json]
The actual carrier is OPH model/polynomial_sv/carrier.py (CELL_SEAMS); the local
phase is the one-load seam Hamiltonian of SPEC.md at kappa=8, duration pi/2.
The directed chart, coherent accept/reject projection, and private route bit
are model inputs. No native TFPT source packet is populated by this campaign.
"""

import argparse
from itertools import combinations, product
import json

import sympy as sp


EDGES = {
    (0, 1), (0, 5), (0, 7), (0, 10), (0, 11), (1, 5), (1, 7), (1, 8), (1, 9), (2, 3),
    (2, 4), (2, 6), (2, 10), (2, 11), (3, 4), (3, 6), (3, 8), (3, 9), (4, 5), (4, 9),
    (4, 11), (5, 9), (5, 11), (6, 7), (6, 8), (6, 10), (7, 8), (7, 10), (8, 9), (10, 11),
}
I2 = sp.eye(2)
HOLD = sp.diag(-sp.I, 1)


def exact(matrix):
    return matrix.applyfunc(sp.expand)


def directed_chart_operators():
    neighbors = {i: {b if i == a else a for a, b in EDGES if i in (a, b)} for i in range(12)}
    face_queue = [(0, 1, 7)]
    flags = {}
    for a, b, c in face_queue:
        if (a, b) in flags:
            assert flags[a, b] == c
            continue
        for i, j, opposite in ((a, b, c), (b, c, a), (c, a, b)):
            flags[i, j] = opposite
            other, = (neighbors[i] & neighbors[j]) - {opposite}
            face_queue.append((j, i, other))
    assert len(flags) == 60
    eye = sp.eye(12)

    def frame(edge):
        a, b = edge
        c = flags[a, b]
        d = flags[c, b]
        assert len({a, b, c, d}) == 4
        return sp.Matrix.hstack(eye[:, a]-eye[:, b], eye[:, c]-eye[:, d])/sp.sqrt(2)

    def turn(edge):
        a, b = edge
        return b, flags[a, b]

    def half_pulse(a, b):
        assert tuple(sorted((a, b))) in EDGES
        difference = eye[:, a]-eye[:, b]
        return eye - (1+sp.I)*difference*difference.T/2

    transfers = set()
    for edge in flags:
        target = turn(edge)
        assert turn(turn(target)) == edge
        basis = frame(edge)
        target_basis = frame(target)
        assert basis.H*basis == I2
        a, b = edge
        # On the local support the eigenvalues are 8 and 9. Off the addressed
        # edge they are zero. At duration pi/2 only the odd line gains -i.
        pulse = half_pulse(a, b)
        assert exact(pulse*basis) == exact(basis*HOLD)
        c, d = target
        e = flags[c, d]
        f = flags[e, d]
        parity = half_pulse(c, d)**2*half_pulse(e, f)**2
        projection = (eye-parity)/2
        assert exact(projection) == target_basis*target_basis.H
        accepted = projection*basis
        rejected = (eye-projection)*basis
        assert exact(accepted.H*accepted+rejected.H*rejected) == I2
        transfers.add(tuple(exact(target_basis.H*basis)))
    assert len(transfers) == 1
    advance = sp.Matrix(2, 2, next(iter(transfers)))
    assert advance == sp.Matrix([[-1, -1], [1, 0]])/2
    return len(flags), advance


def closed_words(advance, depth):
    for letters in product('HT', repeat=depth):
        if letters.count('T') % 3:
            continue
        result = I2
        for gate in letters:
            result = exact((HOLD if gate == 'H' else advance)*result)
        yield ''.join(letters), result


def is_diagonal(matrix):
    return matrix[0, 1] == matrix[1, 0] == 0


def paired_return_search(advance, max_steps):
    """Find equal clock length, diagonal, complementary projective inverses."""
    depths = []
    for depth in range(1, max_steps + 1):
        histories = list(closed_words(advance, depth))
        diagonal_paths = [(word, value) for word, value in histories if is_diagonal(value)]
        witnesses = []
        for (left, a), (right, b) in combinations(diagonal_paths, 2):
            for power in range(4):
                upper = exact(HOLD**power*(a+b)/2)
                lower = exact(HOLD**power*(a-b)/2)
                if not all(x.is_positive for x in (upper[0, 0], upper[1, 1], lower[0, 0], lower[1, 1])):
                    continue
                if upper[0, 0] == upper[1, 1]:
                    continue
                reverse = exact(upper*lower)
                if reverse == reverse[0, 0]*I2:
                    witnesses.append((left, right, power, upper, lower))
        depths.append({'length': depth, 'closed_words': len(histories), 'diagonal_pairs': len(diagonal_paths)*(len(diagonal_paths)-1)//2,
                       'witnesses': len(witnesses)})
        if witnesses:
            return depths, {item[:3]: item[3:] for item in witnesses}
    return depths, {}


def loop_and_charge_checks(advance):
    assert advance**3 == I2/8 and HOLD**4 == I2
    assert advance.inv().applyfunc(sp.expand) == 8*advance**2
    # Determinants of the doubled turn and of the hold are Gaussian units;
    # all projective histories with no coherent addition stay in PGL(2,Z[i]).
    assert (2*advance).det() == 1 and HOLD.det() == -sp.I
    upper = exact(sp.I*HOLD*advance*HOLD*advance.inv())
    lower = exact(sp.I*HOLD*advance.inv()*HOLD*advance)
    assert upper == sp.Matrix([[1, 1+sp.I], [0, 1]])
    assert lower == sp.Matrix([[1, 0], [1-sp.I, 1]])
    assert exact(HOLD*upper*HOLD.inv())[0, 1] == 1-sp.I
    assert exact(HOLD*lower*HOLD.inv())[1, 0] == 1+sp.I

    a, b = HOLD**6, advance**6
    k_plus, k_minus = (a+b)/2, (a-b)/2
    accepted = k_plus.H*k_plus + k_minus.H*k_minus
    missing = I2-accepted
    assert accepted == sp.Rational(4097, 8192)*I2
    assert missing == sp.Rational(4095, 8192)*I2
    interference = k_plus.H*k_plus - k_minus.H*k_minus
    assert interference == sp.diag(-1, 1)/64
    ratio = sp.cancel(k_plus[0, 0]*(-1)/k_plus[1, 1])
    assert ratio == sp.Rational(63, 65)

    # Exhaust all grade-character insertions i^(code*grade) into the two
    # single-event operators. A shared phase cancels out projectively.
    same, opposite = 0, 0
    for hold_grade, turn_grade in product(range(4), repeat=2):
        compared = [sp.I**((6*(turn_grade-hold_grade)*grade) % 4) for grade in range(4)]
        if compared == [1]*4:
            same += 1
            continue
        assert compared == [1, -1, 1, -1]
        odd_branch = exact(HOLD**2*(a-b)/2)
        assert odd_branch == sp.diag(65, 63)/128
        opposite += 1
    assert (same, opposite) == (8, 8)
    # Integer powers of matched upper (and lower) loops can separate the
    # factors for opposite dilations. Further conjugation sends the step to zero.
    p, q = ratio.p, ratio.q
    assert p/ratio-q == q*ratio-p == 0
    assert p*ratio-q != 0 and q/ratio-p != 0
    return ratio, interference, same, opposite


def history_coherence_check(advance):
    # Keep the two route dilations fixed but allow unobserved environment
    # states of overlap z. Exactly |z| <= 1 is allowed by their Gram matrix.
    x, y = sp.symbols('x y', real=True)
    z = x+sp.I*y
    a, b = HOLD**6, advance**6
    assert a.H*a == I2 and b.H*b == I2/4096
    difference = exact((z*b.H*a+sp.conjugate(z)*a.H*b)/2)
    assert difference == x*sp.diag(-1, 1)/64
    vector = sp.Matrix([1, 1])/sp.sqrt(2)
    columns = sp.Matrix.hstack(a*vector, b*vector)
    joint = columns*sp.Matrix([[1, z], [sp.conjugate(z), 1]])*columns.H/4
    assert sp.factor(joint.det()-(1-x*x-y*y)/65536) == 0
    return {'range_on_fixed_basis_input': ['-1/64', '1/64'],
            'passive_committed_records': '0',
            'maximal_allowed_interference_fixing_z': '1/64',
            'finite_density_requires_full_coherence_for_single_filter': True}


def two_site_contact():
    swap = sp.Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    full_swap = sp.kronecker_product(swap, swap)
    echo = ((1-sp.I)*sp.eye(16)+(1+sp.I)*full_swap)/2
    assert exact(echo.H*echo) == sp.eye(16)
    # A site gauge charge of 0 or 1 is moved by the swap. For arbitrary
    # on-site Hilbert spaces the same swap commutes with M tensor M.
    charge = sp.diag(0, 1)
    charge_a = sp.kronecker_product(sp.eye(4), charge, I2)
    charge_b = sp.kronecker_product(sp.eye(4), I2, charge)
    assert exact(echo*(charge_a+charge_b)-(charge_a+charge_b)*echo) == sp.zeros(16)
    assert exact(echo*charge_a-charge_a*echo) != sp.zeros(16)
    a, b, c, d = sp.symbols('a b c d')
    geometry = sp.kronecker_product(sp.Matrix([[a, b], [c, d]]), sp.Matrix([[a, b], [c, d]]))
    simultaneous = sp.kronecker_product(geometry, sp.eye(4))
    assert exact(echo*simultaneous-simultaneous*echo) == sp.zeros(16)
    # Starting in |01>spin|01>charge, the output has two equal Schmidt weights.
    state = echo[:, 5].reshape(4, 4)
    assert exact(sp.Matrix([sp.trace((state*state.H)**2)])) == sp.Matrix([sp.Rational(1, 2)])
    # The very same exchange generator gives a positive Euclidean transfer.
    # A mixed readout with the two site projectors has common normalization.
    transfer = (3*sp.eye(16)+full_swap)/4
    filtered = transfer[:, 5]
    spin_a = sp.kronecker_product(charge, I2, sp.eye(4))
    outcomes = [(filtered.H*o*filtered)[0] for o in (sp.eye(16), spin_a, charge_a, spin_a*charge_a)]
    assert outcomes == [sp.Rational(5, 8), sp.Rational(1, 16), sp.Rational(1, 16), sp.Rational(1, 16)]
    norm, spin, matter, mixed = outcomes
    connected = mixed/norm - spin*matter/norm**2
    assert connected == sp.Rational(9, 100)
    return {'reduced_purity_after_unitary': '1/2', 'euclidean_connected_contact': str(connected)}


def run(max_steps=6):
    charts, advance = directed_chart_operators()
    depths, witnesses = paired_return_search(advance, max_steps)
    ratio, interference, same, opposite = loop_and_charge_checks(advance)
    witness = ('HHHHHH', 'TTTTTT', 2)
    if max_steps >= 6:
        assert depths[-1]['length'] == 6 and witness in witnesses
        assert witnesses[witness] == (sp.diag(63, 65)/128, sp.diag(65, 63)/128)
    return {
        'status': 'EXACT_DECLARED_HISTORY_SEARCH',
        'physical_source_derived': False,
        'physical_solution': False,
        'charts': charts,
        'depths': depths,
        'witness': {'word_a': witness[0], 'word_b': witness[1], 'hold_correction': 2,
                    'positive_operators': [['63/128', '65/128'], ['65/128', '63/128']]}
        if witness in witnesses else None,
        'withheld_interference_eigenvalues': [str(interference[i, i]) for i in range(2)],
        'history_coherence': history_coherence_check(advance),
        'positive_return_ratio': str(ratio),
        'grade_twist_menu': {'projectively_factorizing_on_selected_return': same, 'opposite_spinor_dilations': opposite},
        'two_site_contact': two_site_contact(),
        'open': ['native coherent history producer', 'native charged and regional source binding', 'interacting continuum'],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-steps', type=int, default=6)
    parser.add_argument('--json', action='store_true', help='Print a machine-readable campaign receipt.')
    args = parser.parse_args()
    if args.max_steps < 1:
        parser.error('--max-steps must be positive')
    report = run(args.max_steps)
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        for key, value in report.items():
            print(f'{key}: {value}')
