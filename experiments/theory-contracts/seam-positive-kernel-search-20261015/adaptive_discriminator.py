"""Rank exact, already certified response questions by candidate separation.

This is a finite-candidate diagnostic, not a physical source selector.  The
segment entry point verifies the stored original certificates first.  Model
proposals and source metadata are never silently turned into constraints.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

from certify_response_segment import certify


def import_hypotheses(proposals):
    """Validate structure only. No claim, source or premise is verified here."""
    if not isinstance(proposals, list):
        raise ValueError('AI hypothesis input must be a list')
    required = {'hypothesis', 'source_ids', 'premises', 'additional_assumptions', 'tests'}
    for proposal in proposals:
        if not isinstance(proposal, dict) or set(proposal) != required:
            raise ValueError('Hypothesis requires exactly the declared descriptive fields')
        if not isinstance(proposal['hypothesis'], str) or not proposal['hypothesis'].strip():
            raise ValueError('Hypothesis text is required')
        for field in ('source_ids', 'premises', 'additional_assumptions', 'tests'):
            if not isinstance(proposal[field], list) or any(not isinstance(x, str) or not x.strip() for x in proposal[field]):
                raise ValueError('Hypothesis lists must contain nonempty strings')
        if not proposal['source_ids'] or not proposal['premises'] or not proposal['tests']:
            raise ValueError('Source IDs, premises and tests must be supplied')
    return {'items': proposals, 'source_ids_verified': False,
            'premises_verified': False, 'added_as_constraints': False,
            'scope': 'Unverified external proposals; review against current atlas and original sources'}


def rank_responses(responses, candidate_count):
    """Rank exact real vectors; affine duplicates share one representative.

    Pair separation is counted only among the supplied candidates.  Affine
    rank is likewise restricted to these vectors, not all possible processes.
    Costs are explicit positive rationals, not inferred probabilities.
    """
    if candidate_count < 2:
        raise ValueError('At least two certified candidates are required')
    groups = {}
    for response in responses:
        word = response['word']
        values = tuple(Fraction(x) for x in response['values'])
        cost = Fraction(response.get('cost', len(word)))
        if len(values) != candidate_count or cost <= 0:
            raise ValueError('Wrong candidate count or nonpositive question cost')
        delta = tuple(v-values[0] for v in values[1:])
        pivot = next((v for v in delta if v), None)
        if pivot is None:
            continue
        key = tuple(v/pivot for v in delta)
        groups.setdefault(key, []).append({'word': word, 'values': values, 'cost': cost})

    remaining = {(i, j) for i in range(candidate_count) for j in range(i+1, candidate_count)}
    output = []
    for key, members in groups.items():
        members.sort(key=lambda x: (x['cost'], len(x['word']), x['word']))
    while groups:
        choices = []
        for key, members in groups.items():
            best = members[0]
            separated = {p for p in remaining if best['values'][p[0]] != best['values'][p[1]]}
            choices.append((Fraction(len(separated), 1)/best['cost'], len(separated),
                            best['cost'], best['word'], key, separated))
        # Exact rational score; deterministic cheapest/lexical tie breaking.
        choices.sort(key=lambda x: (-x[0], -x[1], x[2], x[3]))
        score, gain, cost, word, key, separated = choices[0]
        members = groups.pop(key)
        best = members[0]
        output.append({'word': word, 'candidate_values': [str(x) for x in best['values']],
                       'cost': str(cost), 'new_candidate_pairs_separated': gain,
                       'pair_separation_per_cost': str(score),
                       'affine_equivalent_words': [m['word'] for m in members]})
        remaining.difference_update(separated)
    differences = [[Fraction(v)-Fraction(r['candidate_values'][0])
                    for v in r['candidate_values'][1:]] for r in output]
    matrix = sp.Matrix([[sp.Rational(v.numerator, v.denominator) for v in row]
                        for row in differences]) if output else sp.zeros(0, candidate_count-1)
    return {'ranked_questions': output, 'affine_signal_rank': int(matrix.rank()),
            'affine_equivalence_groups': len(output), 'candidate_pairs_unseparated': len(remaining),
            'score_scope': 'exact pair separation divided by declared cost on supplied candidates; not Shannon information or a global bound'}


def analyze_segment(run):
    # certify replays both full stored candidates and all supplied word moments.
    result = certify(run)
    s = sp.Symbol('s', real=True)
    entries = []
    for row in result['responses']:
        expression = sp.sympify(row['polynomial'], locals={'s': s})
        ends = [expression.subs(s, k) for k in (0, 1)]
        if any(not value.is_Rational for value in ends):
            raise ValueError('Response must have exact real rational values')
        entries.append({'word': row['word'], 'values': [str(x) for x in ends]})
    ranking = rank_responses(entries, 2)
    return {'schema': 'tfpt.adaptive-discriminator.v1',
            'certificate_status': result['status'], 'candidate_count': 2,
            'segment_parameter': result['parameter'],
            'segment_replay_moments': result['all_replayed_moment_polynomials_constant'],
            'candidate_carrier': 'same declared 32-dimensional finite CAR representation, state and internal quarter',
            'physical_comparison': 'original P1 carrier and field map are not identified',
            'ranking': ranking, 'source_sha256': result['source_sha256'],
            'physical_source_identified': False, 'fundamental_solution': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--hypotheses', type=Path,
                        help='Optional descriptive JSON list; never inserted as constraints')
    args = parser.parse_args()
    result = analyze_segment(args.run)
    if args.hypotheses:
        result['unverified_ai_proposals'] = import_hypotheses(json.loads(args.hypotheses.read_text()))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
