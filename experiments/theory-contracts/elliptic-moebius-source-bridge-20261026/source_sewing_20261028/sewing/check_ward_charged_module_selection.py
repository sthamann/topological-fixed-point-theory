"""Exact Ward-vector closure test in the existing E8 A4 x A4 chart.

The script imports the canonical A4 response checker and its actual E8 current
source.  It adds no embedding, representation matrix, or fitted coefficient.
It verifies only the source-pinned Ward vector, its neutral-current orbit, the
old D5 x A3 grades, the P2 X charges, and the rank of its single bracket map.
"""
from collections import Counter, deque
import hashlib
import importlib.util
import json
from pathlib import Path


BASE = Path(__file__).resolve().parent
REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
CANON = (REPO / 'experiments/theory-contracts/'
         'elliptic-moebius-source-bridge-20261026/'
         'hodge_complement_20261028')
CANON_CHECKER = CANON / 'check_e8_a4_response.py'
CANON_RESULTS = CANON / 'e8-a4-results.json'
WARD_SOURCE = (REPO / 'experiments/theory-contracts/'
               'seam-joint-ward-state-20261014/PROOF.txt')
OUTPUT = BASE / 'ward-charged-module-selection-results.json'

EXPECTED_CANON_CHECKER_SHA256 = (
    '9659806dfc92ec5f7d322a066d8e6f5b7b7c80e9a5cc285257d5d44bc57ea1f0'
)
EXPECTED_CANON_RESULTS_SHA256 = (
    '2816fbc59ad82469571831c4ae852753a4a147d07fd247d835b0b8adfb49c334'
)
EXPECTED_WARD_SOURCE_SHA256 = (
    '7c0cd92dfc5eca88260baae0af9edad66f958a27ce0d15d7965414b3971f8da0'
)


def require(condition, name):
    if not condition:
        raise RuntimeError(name)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_canonical_checker():
    spec = importlib.util.spec_from_file_location(
        'canonical_e8_a4_response', CANON_CHECKER
    )
    require(spec is not None and spec.loader is not None,
            'canonical A4 checker import specification exists')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def old_d5a3_sector(root):
    """The exact classification used by the canonical A4 response checker."""
    if sum(x != 0 for x in root) == 2:
        carrier_support = sum(x != 0 for x in root[:5])
        return 'g0' if carrier_support in (0, 2) else 'g2'
    require(all(abs(x) == 1 for x in root),
            'non-coordinate E8 root is a doubled halfroot')
    return 'g1' if sum(x > 0 for x in root[:5]) % 2 == 0 else 'g3'


def main():
    require(sha256(CANON_CHECKER) == EXPECTED_CANON_CHECKER_SHA256,
            'canonical A4 checker source is pinned')
    require(sha256(CANON_RESULTS) == EXPECTED_CANON_RESULTS_SHA256,
            'canonical A4 result is pinned')
    require(sha256(WARD_SOURCE) == EXPECTED_WARD_SOURCE_SHA256,
            'original Ward source is pinned')

    ward_text = WARD_SOURCE.read_text()
    require("q -> q'=2(e3+e6)" in ward_text,
            'Ward source contains the marked q to q-prime dictionary')
    require('not a derivation of the raw source map Xi' in ward_text,
            'Ward source preserves the conditional raw-Xi boundary')

    canonical = load_canonical_checker()
    canonical_result = json.loads(CANON_RESULTS.read_text())
    family_simple = tuple(
        tuple(root) for root in canonical_result['family_simple_roots']
    )
    require(len(family_simple) == 4,
            'canonical family A4 has four simple roots')

    source = canonical._actual_current_source()
    algebra = source['chevalley']
    roots = tuple(sorted(algebra.roots))
    root_set = set(roots)
    require(len(roots) == 240, 'actual E8 current source has 240 roots')

    ward = (0, 0, 2, 0, 0, 2, 0, 0)
    require(ward in root_set, 'marked Ward vector is an actual E8 root')
    require(canonical.carrier_type(ward) == '5',
            'marked Ward vector is in the carrier 5')
    carrier_weight = canonical.carrier_weight(ward)
    require(carrier_weight == (-2, -2, 8, -2, -2),
            'marked Ward vector has the exact carrier weight')
    family_weight = tuple(canonical.dot(ward, root)//4
                          for root in family_simple)
    require(family_weight == (0, 1, 0, 0),
            'marked Ward vector has the exact family Lambda2 weight')

    target = {root for root in roots
              if canonical.carrier_type(root) == '5'}
    require(len(target) == 50, 'target (5,10) sector has 50 roots')

    carrier_neutral = {
        root for root in roots
        if sum(x != 0 for x in root[:5]) == 2
        and all(x == 0 for x in root[5:])
        and sum(root[:5]) == 0
    }
    family_neutral = {
        root for root in roots if len(set(root[:5])) == 1
    }
    neutral = carrier_neutral | family_neutral
    require(len(carrier_neutral) == 20 and len(family_neutral) == 20,
            'the two neutral A4 factors each have 20 roots')
    require(len(neutral) == 40,
            'neutral A4 x A4 current set has 40 roots')

    reached = {ward}
    pending = deque([ward])
    directed_edges = 0
    while pending:
        current = pending.popleft()
        for neutral_root in neutral:
            candidate = canonical.add(current, neutral_root)
            if candidate in root_set:
                require(candidate in target,
                        'neutral current does not leave the target module')
                directed_edges += 1
                if candidate not in reached:
                    reached.add(candidate)
                    pending.append(candidate)
    require(reached == target,
            'neutral current closure of the Ward vector is the full (5,10) sector')
    require(directed_edges == 500,
            'full neutral closure has exactly 500 directed root edges')

    grading = dict(sorted(Counter(old_d5a3_sector(root)
                                  for root in target).items()))
    require(grading == {'g2': 30, 'g3': 20},
            'target sector splits as 30 old-g2 plus 20 old-g3 roots')

    p2_x = (-2, -2, -2, 3, 3)
    x_spectrum = dict(sorted(Counter(
        p2_x[canonical.carrier_label(root, '5')] for root in target
    ).items()))
    require(x_spectrum == {-2: 30, 3: 20},
            'P2 X spectrum on the insertion sector is 30 plus 20')
    require(p2_x[canonical.carrier_label(ward, '5')] == -2,
            'marked Ward vector lies in the X=-2 part')

    reflected = tuple(-x for x in ward)
    require(reflected in root_set, 'compact reflection is an E8 root')
    require(canonical.carrier_type(reflected) == 'bar5',
            'compact reflection lies in the conjugate carrier sector')
    reflected_family_weight = tuple(canonical.dot(reflected, root)//4
                                    for root in family_simple)
    require(reflected_family_weight == (0, -1, 0, 0),
            'compact reflection has conjugate family weight')

    input_sector = [root for root in roots
                    if canonical.carrier_type(root) == '10']
    output_sector = {root for root in roots
                     if canonical.carrier_type(root) == 'bar10'}
    require(len(input_sector) == len(output_sector) == 50,
            'single bracket has the canonical 50 by 50 source and target sectors')

    nonzero_columns = []
    output_roots = set()
    ward_index = algebra.ridx[ward]
    for input_root in input_sector:
        raw = algebra.bracket(ward_index, algebra.ridx[input_root])
        if not raw:
            continue
        require(len(raw) == 1,
                'root-vector bracket has one root-vector output')
        output_index, coefficient = next(iter(raw.items()))
        candidate = canonical.add(ward, input_root)
        require(candidate in output_sector,
                'nonzero marked bracket lands in the bar10 output sector')
        require(output_index == algebra.ridx[candidate],
                'actual Chevalley bracket has the root-sum output')
        require(coefficient in (-1, 1),
                'actual Chevalley bracket coefficient is a cocycle sign')
        nonzero_columns.append(input_root)
        output_roots.add(candidate)
    require(len(nonzero_columns) == 12,
            'single marked bracket has 12 nonzero input columns')
    require(len(output_roots) == 12,
            'single marked bracket has rank 12')

    require(canonical_result['kraus_rank'] == 12,
            'canonical inclusive response records rank 12 per bracket')
    require(canonical_result['source_gain'] == '12 I50',
            'canonical full 50-current sum has gain 12 I50')

    result = {
        'verdict': 'PASS',
        'scope': (
            'Conditional algebraic selection from the pinned Ward dictionary '
            'and the actual local E8 current closure. This does not derive the '
            'raw source map Xi or the inclusive reflected sewing operation.'
        ),
        'sources': {
            'ward_source': {
                'path': str(WARD_SOURCE),
                'sha256': sha256(WARD_SOURCE),
                'pinned_statement': "q -> q'=2(e3+e6)",
            },
            'canonical_a4_checker': {
                'path': str(CANON_CHECKER),
                'sha256': sha256(CANON_CHECKER),
            },
            'canonical_a4_results': {
                'path': str(CANON_RESULTS),
                'sha256': sha256(CANON_RESULTS),
            },
        },
        'ward_root_doubled': list(ward),
        'carrier_weight': list(carrier_weight),
        'family_dynkin_weight': list(family_weight),
        'selected_module': '(5,10)',
        'neutral_roots': len(neutral),
        'neutral_orbit_size': len(reached),
        'neutral_directed_edges': directed_edges,
        'orbit_equals_full_target_sector': reached == target,
        'old_D5xA3_grading': grading,
        'p2_X': list(p2_x),
        'p2_X_insertion_charge_spectrum': {
            str(charge): multiplicity
            for charge, multiplicity in x_spectrum.items()
        },
        'ward_p2_X_charge': -2,
        'reflected_module': '(bar5,bar10)',
        'single_ward_bracket_nonzero_columns': len(nonzero_columns),
        'single_ward_bracket_rank': len(output_roots),
        'inclusive_current_count': canonical_result['kraus_count'],
        'inclusive_source_gain': canonical_result['source_gain'],
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + '\n'
    OUTPUT.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
