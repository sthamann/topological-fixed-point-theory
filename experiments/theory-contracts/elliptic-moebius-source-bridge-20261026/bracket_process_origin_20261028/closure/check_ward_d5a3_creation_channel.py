"""Exact Ward/U5 selection of the creation half of the actual E8 g2 channel.

This check uses the pinned Ward root, the original P2 X-mod-5 filter, and the
actual E8 Chevalley source.  It proves algebraic current support and the
resulting normalized bracket channel.  It does not assert that P1 physically
applies the current projector or forgets its event record.
"""
from collections import Counter, deque
import hashlib
import json
from pathlib import Path
import sys


BASE = Path(__file__).resolve().parent
REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0, str(REPO))
from tfpt_explorer.current_block_geometry import _actual_current_source

CANON = (REPO / 'experiments/theory-contracts/'
         'elliptic-moebius-source-bridge-20261026/'
         'hodge_complement_20261028')
CANON_CHECKER = CANON / 'check_e8_clifford_response.py'
CANON_RESULT = CANON / 'e8-clifford-response-results.json'
WARD_SOURCE = (REPO / 'experiments/theory-contracts/'
               'seam-joint-ward-state-20261014/PROOF.txt')
P2_SOURCE = (REPO / 'documents/'
             'archive--tfpt-45--source_extracts--02_carrier_source/'
             'sections/003-subsection-character-trace-identities-and-the-abelian-in.tex')
CURRENT_SOURCE = REPO / 'tfpt_explorer/current_block_geometry.py'
CHEVALLEY_SOURCE = REPO / 'verification/v498_celestial_wp5b_singular_vector.py'
QUARTER_REPORT = BASE.parent / 'original/SOURCE-CHARGE-GLUE-QUARTER-IDENTITY.txt'
QUARTER_RESULT = BASE.parent / 'source-charge-glue-quarter-results.json'
OUTPUT = BASE / 'ward-d5a3-creation-channel-results.json'

EXPECTED = {
    CANON_CHECKER: '75ce80dbec9d943573ba96e64eaf85088c0cbc47577673cbe99673ba34957aa9',
    CANON_RESULT: '67fe9e3f06140c5c23452e9a7bfaefc02441ff1e7b55708e55409cd0a398e7ee',
    WARD_SOURCE: '7c0cd92dfc5eca88260baae0af9edad66f958a27ce0d15d7965414b3971f8da0',
    P2_SOURCE: 'dbd7f1020e414aeaa35ab3a7191fc5863d6b0afda270a2e347e9a0d6af2b0284',
    CURRENT_SOURCE: 'c9efcac4ba33d4ee3dc20f8c3c69a0f6ff6b2084893a5edf4aaae6c204d9056b',
    CHEVALLEY_SOURCE: '9f5a2d62523283595b754fea6f70ddc6501a251daf86752cbc67ded252b96180',
    QUARTER_REPORT: 'adbf414bf1105de1087740b3bd1707e9f0dcfe7504eb0da4467a0e01a2289b9c',
    QUARTER_RESULT: '19c5761a525a5b3be0265acc8f33598a106a9fe3aca73974fba7ccd0519fede6',
}


def need(condition, name):
    if not condition:
        raise RuntimeError(name)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_label(path):
    # Keep the certificate identical after copying the whole proof folder.
    if path == QUARTER_REPORT:
        return '../original/SOURCE-CHARGE-GLUE-QUARTER-IDENTITY.txt'
    if path == QUARTER_RESULT:
        return '../source-charge-glue-quarter-results.json'
    return str(path.relative_to(REPO))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def key(root):
    return ','.join(map(str, root))


def ladder(mask, i, create):
    occupied = (mask >> i) & 1
    if occupied == create:
        return None
    sign = -1 if (mask & ((1 << i) - 1)).bit_count() % 2 else 1
    return mask ^ (1 << i), sign


def orbit(seed, generators, root_set):
    reached = {seed}
    pending = deque([seed])
    edges = 0
    while pending:
        current = pending.popleft()
        for generator in generators:
            candidate = add(current, generator)
            if candidate in root_set:
                edges += 1
                if candidate not in reached:
                    reached.add(candidate)
                    pending.append(candidate)
    return reached, edges


def main():
    for path, expected in EXPECTED.items():
        need(sha256(path) == expected, f'pinned source changed: {path}')
    ward_text = WARD_SOURCE.read_text()
    p2_text = P2_SOURCE.read_text()
    need("q -> q'=2(e3+e6)" in ward_text,
         'pinned Ward dictionary contains q-prime=2(e3+e6)')
    need('not a derivation of the raw source map Xi' in ward_text,
         'Ward source retains the raw-Xi scope boundary')
    need('X:=6Y=3N_2-2N_3' in p2_text,
         'P2 source defines the marked 3+2 charge X')
    need('U_5:=e^{2\\pi iX/5}' in p2_text,
         'P2 source defines the five-fold phase filter U5')

    source = _actual_current_source()
    algebra = source['chevalley']
    roots = tuple(sorted(algebra.roots))
    root_set = set(roots)
    halves = tuple(sorted(r for r in roots if all(abs(x) == 1 for x in r)))
    d8 = tuple(sorted(r for r in roots if sum(x != 0 for x in r) == 2))
    mixed = tuple(r for r in d8
                  if sum(x != 0 for x in r[:5]) == 1
                  and sum(x != 0 for x in r[5:]) == 1)
    need((len(roots), len(halves), len(d8), len(mixed)) == (240, 128, 112, 60),
         'actual E8 source sector counts')

    ward = (0, 0, 2, 0, 0, 2, 0, 0)
    need(ward in root_set, 'Ward vector is an actual E8 root')

    # The first-five carrier sign is the exterior creation/annihilation sign
    # in the already verified D5+A3 CAR intertwiner.
    creation = {
        root for root in mixed
        if next(root[i] for i in range(5) if root[i]) == 2
    }
    annihilation = set(mixed) - creation
    need(len(creation) == len(annihilation) == 30,
         'mixed g2 splits into 30 creation plus 30 annihilation roots')
    need(ward in creation and tuple(-x for x in ward) in annihilation,
         'Ward orientation selects the creation half and reflection its dual')

    # Original neutral SU5_C x D3=A3_F root algebra.  This is an internal
    # current algebra statement, not a claim that full SU5 is physical gauge.
    su5 = {
        root for root in d8
        if all(x == 0 for x in root[5:]) and sum(root[:5]) == 0
    }
    d3 = {
        root for root in d8
        if all(x == 0 for x in root[:5])
    }
    need((len(su5), len(d3)) == (20, 12),
         'neutral SU5 and D3 root counts are 20 and 12')
    ward_orbit, neutral_edges = orbit(ward, su5 | d3, root_set)
    need(ward_orbit == creation,
         'neutral SU5 x D3 closure of Ward root is exactly creation30')

    # P2 X=(-2,-2,-2,3,3) supplies a stronger spectral description.  In the
    # adjoint E8 current representation the root current E_beta has X weight
    # X(beta).  All creation roots have residue 3 and all annihilation roots
    # residue 2 modulo 5.  Thus the original Fourier filter separates them.
    x5 = (-2, -2, -2, 3, 3)

    def x_charge(root):
        value = sum(x5[i] * root[i] for i in range(5))
        need(value % 2 == 0, 'doubled root gives integral P2 X charge')
        return value // 2

    creation_x = Counter(x_charge(root) for root in creation)
    annihilation_x = Counter(x_charge(root) for root in annihilation)
    need(creation_x == Counter({-2: 18, 3: 12}),
         'creation30 retains the marked 3+2 multiplicities 18+12')
    need(annihilation_x == Counter({2: 18, -3: 12}),
         'annihilation30 has the conjugate 3+2 multiplicities 18+12')
    need({x_charge(root) % 5 for root in creation} == {3},
         'U5 residue3 is exactly the creation orientation')
    need({x_charge(root) % 5 for root in annihilation} == {2},
         'U5 residue2 is exactly the annihilation orientation')
    need({root for root in mixed if x_charge(root) % 5 == 3} == creation,
         'the U5 Fourier projector selects exactly creation30 inside g2')

    # U5 alone has a 50-dimensional residue-3 sector on the whole E8 root
    # system.  The 30 is the joint eigenspace with the existing internal Z4
    # glue grade2.  This prevents promoting an in-g2 statement to a false
    # global projector identity.
    def old_internal_z4_grade(root):
        if sum(x != 0 for x in root) == 2:
            carrier_support = sum(x != 0 for x in root[:5])
            return 0 if carrier_support in (0, 2) else 2
        need(all(abs(x) == 1 for x in root),
             'non-coordinate E8 root is a doubled halfroot')
        return 1 if sum(x > 0 for x in root[:5]) % 2 == 0 else 3

    joint_census = Counter(
        (old_internal_z4_grade(root), x_charge(root) % 5)
        for root in roots
    )
    expected_joint_census = Counter({
        (0, 0): 32, (0, 1): 10, (0, 4): 10,
        (1, 0): 4, (1, 1): 40, (1, 2): 20,
        (2, 2): 30, (2, 3): 30,
        (3, 0): 4, (3, 3): 20, (3, 4): 40,
    })
    need(joint_census == expected_joint_census,
         'complete 240-root old-Z4 by U5-residue census')
    full_residue3 = {root for root in roots if x_charge(root) % 5 == 3}
    joint_creation = {
        root for root in full_residue3 if old_internal_z4_grade(root) == 2
    }
    residue3_remainder = full_residue3 - joint_creation
    need(len(full_residue3) == 50,
         'U5 residue3 has dimension50 on the full E8 root system')
    need(joint_creation == creation and len(residue3_remainder) == 20,
         'old internal grade2 intersect U5 residue3 is creation30 plus remainder20')
    need({old_internal_z4_grade(root) for root in residue3_remainder} == {3},
         'the additional residue3 roots have old internal grade3')

    quarter_result = json.loads(QUARTER_RESULT.read_text())
    need(quarter_result['formula'] ==
         'g_alpha^2 = exp(2 pi i alpha Q_source/4) C_internal^(-1)',
         'charged fold square recovers inverse internal glue after charge removal')
    need(quarter_result['full_U5_residue3_dimension'] == 50 and
         quarter_result['creation_joint_dimension'] == 30,
         'independent charged-fold proof has the same 50 to30 joint filter')

    # Physical fixed-mark SU3 x SU2 plus family D3 does not by itself close
    # the 30.  It splits the Ward orbit into the color18 and weak12 blocks.
    su3 = {r for r in su5 if all(r[i] == 0 for i in (3, 4))}
    su2 = {r for r in su5 if all(r[i] == 0 for i in (0, 1, 2))}
    physical_neutral = su3 | su2 | d3
    color_orbit, _ = orbit(ward, physical_neutral, root_set)
    weak_seed = (0, 0, 0, 2, 0, 2, 0, 0)
    weak_orbit, _ = orbit(weak_seed, physical_neutral, root_set)
    need(len(color_orbit) == 18 and len(weak_orbit) == 12,
         'fixed 3+2 physical gauge splits creation30 into color18 plus weak12')
    need(color_orbit | weak_orbit == creation and not color_orbit & weak_orbit,
         'the two fixed-mark gauge orbits partition creation30')

    # Reuse the complete signed D8/CAR intertwiner, but test the selected
    # 30-root restriction and its actual Chevalley action independently.
    canonical = json.loads(CANON_RESULT.read_text())
    state_phase = canonical['state_phases']
    root_phase = canonical['root_phases']
    masks = {r: sum((x > 0) << i for i, x in enumerate(r)) for r in halves}
    input_states = {
        r for r in halves if (masks[r] & 31).bit_count() == 2
    }
    output_states = {
        r for r in halves if (masks[r] & 31).bit_count() == 3
    }
    need(len(input_states) == len(output_states) == 40,
         'typed source and target are Lambda2 C5 tensor4 and Lambda3 C5 tensor4bar')

    input_gain = Counter()
    output_gain = Counter()
    signed_columns = 0
    color_gain = Counter()
    weak_gain = Counter()
    for beta in creation:
        carrier_i = next(i for i in range(5) if beta[i])
        family_j = next(j for j in range(5, 8) if beta[j]) - 5
        for alpha in input_states:
            raw = algebra.bracket(algebra.ridx[beta], algebra.ridx[alpha])
            carrier = ladder(masks[alpha] & 31, carrier_i, True)
            family = ladder(masks[alpha] >> 5, family_j, beta[family_j + 5] > 0)
            if carrier is None or family is None:
                need(not raw, 'selected tensor ladder and actual bracket have same zeros')
                continue
            candidate = add(alpha, beta)
            need(candidate in output_states, 'creation bracket has typed Lambda3 output')
            need(set(raw) == {algebra.ridx[candidate]},
                 'actual selected bracket has one root-vector output')
            coefficient = raw[algebra.ridx[candidate]]
            need(coefficient in (-1, 1), 'actual selected bracket coefficient is signed unit')
            need(masks[candidate] == carrier[0] + (family[0] << 5),
                 'selected bracket has exact carrier and family target masks')
            phase = (int(state_phase[key(alpha)]) *
                     int(state_phase[key(candidate)]) *
                     int(root_phase[key(beta)]))
            need(coefficient == phase * carrier[1] * family[1],
                 'selected current has the pinned signed tensor factorization')
            input_gain[alpha] += coefficient * coefficient
            output_gain[candidate] += coefficient * coefficient
            (color_gain if carrier_i < 3 else weak_gain)[alpha] += 1
            signed_columns += 1

    need(signed_columns == 360, '30 selected currents have 360 signed nonzero columns')
    need(set(input_gain.values()) == {9} and len(input_gain) == 40,
         'selected sum B*B is 9 I40')
    need(set(output_gain.values()) == {9} and len(output_gain) == 40,
         'selected sum BB* is 9 I40')
    need(len(set(color_gain.values())) > 1 and len(set(weak_gain.values())) > 1,
         'neither fixed-mark gauge block alone has scalar source gain')
    need(all(color_gain[s] + weak_gain[s] == 9 for s in input_states),
         'color and weak blocks combine to the scalar gain9')

    # The carrier marginal is already the pure Petersen channel: creation on
    # Lambda2 followed by the canonical volume-form complement.  There is no
    # annihilation/carry branch in the residue-3 current sector.
    petersen_edges = 0
    for mask in range(32):
        if mask.bit_count() != 2:
            continue
        targets = set()
        for i in range(5):
            out = ladder(mask, i, True)
            if out is None:
                continue
            folded = 31 ^ out[0]
            need(folded.bit_count() == 2 and not (mask & folded),
                 'creation plus Hodge complement is a Petersen edge')
            targets.add(folded)
            petersen_edges += 1
        need(len(targets) == 3, 'each duad has exactly three Petersen outputs')
    need(petersen_edges == 30, 'directed Petersen support has 30 edges')

    # Full root-current closure is much larger and therefore is not the
    # permitted event algebra.  This records the need for the U5 sector rule.
    full_orbit, _ = orbit(ward, set(roots), root_set)
    need(full_orbit == root_set, 'unrestricted E8 root closure leaves the selected sector')

    result = {
        'verdict': 'PASS',
        'scope': (
            'Exact algebraic support theorem in the actual E8 current source. '
            'The joint projectors designate a charged field component, not '
            'an executed measurement. The check does not derive raw Xi, the '
            'local oriented P1 field/process map, an event record, or time.'
        ),
        'atlas_revision_read': '63fab31610450cc88e168232',
        'ward_root_doubled': list(ward),
        'mixed_g2': 60,
        'selected_creation_currents': len(creation),
        'dual_annihilation_currents': len(annihilation),
        'neutral_SU5xD3_orbit': len(ward_orbit),
        'neutral_orbit_directed_edges': neutral_edges,
        'P2_X_creation_spectrum': {str(k): v for k, v in sorted(creation_x.items())},
        'P2_X_annihilation_spectrum': {str(k): v for k, v in sorted(annihilation_x.items())},
        'U5_current_residue': {'creation': 3, 'annihilation': 2},
        'full_240_root_joint_grade_census': {
            f'z4_{grade}_u5_{residue}': count
            for (grade, residue), count in sorted(joint_census.items())
        },
        'full_E8_U5_residue3': len(full_residue3),
        'joint_old_Z4_grade2_U5_residue3': len(joint_creation),
        'residue3_old_Z4_grade3_remainder': len(residue3_remainder),
        'charged_fold_glue_identity': quarter_result['formula'],
        'fixed_mark_gauge_orbits': {'color': len(color_orbit), 'weak': len(weak_orbit)},
        'typed_input_dimension': len(input_states),
        'typed_output_dimension': len(output_states),
        'signed_columns_checked': signed_columns,
        'selected_source_gain': '9 I40',
        'normalized_selected_channel': (
            '(1/3 sum_m a_m^dagger rho a_m) tensor '
            '(1/3 sum_j(b_j rho b_j^dagger+b_j^dagger rho b_j))'
        ),
        'carrier_marginal': 'full Hodge channel Phi; Petersen A/3 is its population matrix',
        'carry_branch': 'absent in U5 residue3; it lies in the conjugate residue2 sector',
        'normalization_identity': '1/9 = (5/3)*(1/15), conditioning the old 60-current response on its 3/5 creation branch',
        'unrestricted_E8_root_orbit': len(full_orbit),
        'sources': [
            {'path': source_label(path), 'sha256': sha256(path)}
            for path in EXPECTED
        ],
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + '\n'
    OUTPUT.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
