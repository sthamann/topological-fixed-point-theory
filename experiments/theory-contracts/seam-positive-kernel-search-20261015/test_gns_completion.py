"""Source-linked tests for faithful mixed-state completion and query feedback."""
import copy
import json
from pathlib import Path
import sys
import unittest

import sympy as sp

from finite_completion import complete
from run_search import discriminate
from word_algebra import WordAlgebra, T
from gns_completion import first_missing_probe


def control():
    sibling = Path(__file__).resolve().parents[1]/'seam-collar-multifield-20261016'
    if str(sibling) not in sys.path:
        sys.path.insert(0, str(sibling))
    from collar_control import thermal_gns_transfer
    _, _, transfer, omega = thermal_gns_transfer()
    gammas = [sp.kronecker_product(g, sp.eye(2)) for g in
              (sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]))]

    def moment(word):
        vector = omega
        for token in reversed(word):
            vector = (transfer if token == T else gammas[token])*vector
        return sp.expand((omega.adjoint()*vector)[0])

    payload = {'m': 1, 'degree': 2, 'transfer_kind': 'thermal_gns',
               'reference_eigenvalue': 1, 'transfer_floor': '1/4', 'transfer_ceiling': 4,
               'source_observables': {'source': str(sibling/'collar_control.py'),
                 'scope': 'Existing given quarter thermal GNS control, not P1-derived input'},
               'moments': [{'word': list(w), 're': str(sp.re(moment(w))),
                            'im': str(sp.im(moment(w)))} for w in WordAlgebra(1).words(4)]}
    return payload, moment


def campaign():
    payload, oracle = control()
    # The probe family is generated from the CAR cyclic basis. Enumerate its
    # shortest new Hermitian word, rather than supply the expected source value.
    probe = first_missing_probe(payload)
    first = complete(payload, 'min-probe', representation='faithful_car_gns', probe_word=probe)
    second = complete(payload, 'max-probe', representation='faithful_car_gns', probe_word=probe)
    if not first['ok'] or not second['ok']:
        raise RuntimeError(str((first, second)))
    difference = discriminate(first, second, payload)
    question = tuple(difference['differences'][0]['word'])
    answer = oracle(question)
    updated = copy.deepcopy(payload)
    updated['moments'].append({'word': list(question), 're': str(sp.re(answer)), 'im': str(sp.im(answer))})
    after_min = complete(updated, 'min-probe', representation='faithful_car_gns', probe_word=question)
    after_max = complete(updated, 'max-probe', representation='faithful_car_gns', probe_word=question)
    if not after_min['ok'] or not after_max['ok']:
        raise RuntimeError(str((after_min, after_max)))
    return {'input': payload, 'first': first, 'second': second, 'discriminator': difference,
            'query': list(question), 'source_answer': str(answer), 'after_input': updated,
            'after_min': after_min, 'after_max': after_max,
            'same_operator_after_feedback': after_min['transfer'] == after_max['transfer'],
            'origin': 'Supplied existing thermal quarter control; no P1 derivation',
            'full_tfpt_solution': False}


class FaithfulGnsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = campaign()

    def test_mixed_source_recovered_without_pure_fock_fallback(self):
        r = self.result
        self.assertTrue(r['first']['ok'])
        self.assertEqual(r['first']['dimension'], 4)
        self.assertFalse(r['first']['pure_state_assumed'])
        self.assertFalse(complete(r['input'])['ok'])

    def test_missing_response_is_generated_then_queried(self):
        r = self.result
        self.assertEqual(r['discriminator']['shortest_found_degree'], 5)
        self.assertTrue(r['discriminator']['shortest_certified'])
        # The canonical Clifford word carries the reversed-pair minus sign.
        self.assertEqual(r['source_answer'], '-1')
        self.assertTrue(r['same_operator_after_feedback'])
        self.assertFalse(r['first']['unique_in_fixed_faithful_representation'])
        self.assertEqual(r['after_min']['one_T_metric_entries_determined'], 16)
        self.assertTrue(r['after_min']['unique_in_fixed_faithful_representation'])

    def test_recovered_operator_has_actual_source_spectrum(self):
        transfer = sp.Matrix([[sp.Rational(x['re'])+sp.I*sp.Rational(x['im']) for x in row]
                              for row in self.result['after_min']['transfer']])
        self.assertEqual(transfer.eigenvals(), {sp.Integer(1): 2, sp.Rational(1,4): 1, sp.Integer(4): 1})

    def test_bad_multi_transfer_data_never_certified(self):
        bad = copy.deepcopy(self.result['after_input'])
        for row in bad['moments']:
            if row['word'] == [0, T, T, 0]:
                row['re'] = '0'
        receipt = complete(bad, 'feasibility', representation='faithful_car_gns')
        self.assertFalse(receipt['ok'])
        self.assertIn('Exact replay failed', receipt['reason'])

    def test_nonfaithful_state_requires_quotient(self):
        bad = copy.deepcopy(self.result['input'])
        for row in bad['moments']:
            if row['word'] == [0, 1]:
                row['im'] = '1'
        self.assertFalse(complete(bad, 'feasibility', representation='faithful_car_gns')['ok'])

    def test_new_exact_control_is_not_dropped_outside_old_prefix(self):
        payload = copy.deepcopy(self.result['input'])
        payload['fixed_moments'] = [{'word': self.result['query'], 're': '-1', 'im': '0'}]
        receipt = complete(payload, 'feasibility', representation='faithful_car_gns')
        self.assertTrue(receipt['ok'], receipt.get('reason'))
        self.assertTrue(receipt['unique_in_fixed_faithful_representation'])

    def test_impossible_external_control_is_replayed(self):
        payload = copy.deepcopy(self.result['input'])
        payload['fixed_moments'] = [{'word':[T]*5,'re':'2','im':'0'}]
        receipt = complete(payload,'feasibility',representation='faithful_car_gns')
        self.assertFalse(receipt['ok'])
        self.assertIn('Exact replay failed',receipt['reason'])

    def test_two_mode_neutral_transfers_in_existing_source(self):
        from collar_control import build_source_choi_control
        source = build_source_choi_control()
        q = source.transfer
        transfer = sp.kronecker_product(q, q.inv().T)
        omega = sp.Matrix([q[i,j]/sp.sqrt(sp.trace(q*q)) for i in range(4) for j in range(4)])
        gamma = [sp.kronecker_product(g,sp.eye(4)) for g in source.gammas]
        algebra = WordAlgebra(2)
        basis = [w for w in algebra.words(4) if T not in w]
        def moment(word):
            vector = omega
            for token in reversed(word):
                vector = (transfer if token == T else gamma[token])*vector
            return sp.expand((omega.adjoint()*vector)[0])
        words = set(basis)
        words.update(left+(T,)+right for left in basis for right in basis)
        payload = {'m':2,'transfer_kind':'thermal_gns','transfer_floor':'1/16',
                   'transfer_ceiling':16,'reference_eigenvalue':1,
                   'moments':[{'word':list(w),'re':str(sp.re(moment(w))),
                               'im':str(sp.im(moment(w)))} for w in sorted(words,key=lambda x:(len(x),x))]}
        receipt = complete(payload,'feasibility',representation='faithful_car_gns')
        self.assertTrue(receipt['ok'],receipt.get('reason'))
        self.assertEqual(receipt['dimension'],16)
        self.assertEqual(receipt['one_T_metric_entries_determined'],256)
        # Actual neutral transfers between two orthogonal source arc modes.
        c0=(gamma[0]+sp.I*gamma[1])/2
        c1=(gamma[2]+sp.I*gamma[3])/2
        n0,n1=c0.adjoint()*c0,c1.adjoint()*c1
        e01,e10=c0.adjoint()*c1,c1.adjoint()*c0
        self.assertEqual(e01*e10,n0*(sp.eye(16)-n1))
        self.assertEqual(e01*e10-e10*e01,n0-n1)


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--output':
        Path(sys.argv[2]).write_text(json.dumps(campaign(), indent=2, sort_keys=True)+'\n')
    else:
        unittest.main()
