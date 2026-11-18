"""Reproduce the bounded original charged-pair control, not a TFPT selection."""
import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
SEARCH = HERE.parent
COLLAR = SEARCH.parent / 'seam-collar-multifield-20261016'
sys.path[:0] = [str(SEARCH), str(COLLAR)]

from charged_pair_control import moment, source_lift, RESPONSE
from finite_completion import complete
from kernel_search import load_source_control
from run_search import campaign


def write(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')


def matrix(value):
    def number(entry):
        return sp.Rational(*entry['re']) + sp.I*sp.Rational(*entry['im'])
    return sp.Matrix([[number(entry) for entry in row] for row in value])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True, help='fresh output directory')
    out = parser.parse_args().out
    out.mkdir(parents=True, exist_ok=False)
    original = json.loads((COLLAR/'results/charged-pair-control.json').read_text())
    source_values = {tuple(row['word']): (Fraction(row['re']), Fraction(row['im']))
                     for row in original['fixed_moments']}
    control = load_source_control(HERE/'source-control.json')
    for row in control['fixed_moments']:
        if source_values[tuple(row['word'])] != (Fraction(row['re']), Fraction(row['im'])):
            raise ValueError('Declared initial moment differs from original source')
    a, b = source_values[(0,-1,0)][0], source_values[(2,-1,2)][0]
    if a != Fraction(1,4) or b != Fraction(1,64) or a-16*b:
        raise ValueError('Original source does not satisfy declared relation')

    summary = campaign(out/'campaign',m=2,degree=2,**control)
    rows = [{'word': [], 're': 1.0, 'im': 0.0}]
    rows += [{'word': row['word'], 're': float(Fraction(row['re'])),
              'im': float(Fraction(row['im']))}
             for row in original['fixed_moments'] if len(row['word']) <= 4]
    payload = {'m': 2, 'degree': 2, 'transfer_floor': 0, 'moments': rows,
               'linear_relations': control['linear_relations']}
    write(out/'exact-control-input.json', payload)
    certificate = complete(payload,'min-trace')
    write(out/'exact-completion.json', certificate)
    if not certificate['ok']:
        raise RuntimeError('Exact completion did not pass: '+str(certificate))

    word = [0,2,-1,2,0]
    original_calculated = moment(word,source_lift(RESPONSE))
    original_stored = source_values[tuple(word)][0]
    candidate = next(row for row in certificate['candidate_predictions'] if row['word'] == word)
    candidate_value = Fraction(candidate['value']['re'])
    exact = certificate['exact_matrices']
    omega = matrix(exact['Omega'])
    gammas = [matrix(gamma) for gamma in exact['gammas']]
    transfer = matrix(exact['T'])
    operation = sp.eye(transfer.rows)
    for token in word:
        operation *= transfer if token == -1 else gammas[token]
    independently_calculated = sp.expand((omega.adjoint()*operation*omega)[0])
    if (original_calculated != original_stored or original_stored != Fraction(1,256)
            or candidate_value != 0 or independently_calculated != candidate_value):
        raise AssertionError('Independent source/candidate word comparison failed')
    write(out/'result.json', {
        'scope': 'Conditional given alpha=1/4 and declared ell=4log4 single-chiral particle/hole source; no P1/P2 physical-source selection.',
        'solver_status': summary['solver_status'],
        'campaign_verdict': summary['verdict'],
        'campaign_completions': summary['completions'],
        'exact_control_result': certificate['verdict'],
        'exact_control_input_original_moments': len(rows)-1,
        'normalization_moments': 1,
        'all_degree_four_words_asserted_equal': False,
        'distinguishing_word': word,
        'original_word_computed_from_source_operator': str(original_calculated),
        'original_word_stored_source': str(original_stored),
        'exact_completion_word': str(candidate_value),
        'exact_completion_word_recomputed_from_certificate_matrices': str(independently_calculated),
        'notes': 'Original single-field values are supplied for exact control; unequal degree-five value proves those inputs do not determine the whole original kernel. Internal quarter is distinct from charged physical quarter.',
    })
    print(json.dumps(json.loads((out/'result.json').read_text()),indent=2))


if __name__ == '__main__':
    main()
