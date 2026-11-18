"""Finite completion must replay every original word, not just fit linear data."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
import sympy as sp
from word_algebra import words, T
from finite_completion import complete, Representation, matrix_json


def fixture(m=1):
    dim=1<<m
    diagonal=[sp.Rational(1,2) if m==1 else sp.Rational(1,4)]*(dim-1)+[1]
    transfer=sp.SparseMatrix(sp.diag(*diagonal))
    representation=Representation(m,dim-1,transfer)
    moments=[]
    for word in words(4,m=m):
        value=representation.moment(word)
        moments.append({'word':list(word),'re':float(sp.re(value)),'im':float(sp.im(value))})
    return {'m':m,'degree':2,'transfer_floor':.25,'moments':moments}


class FiniteCompletionTests(unittest.TestCase):
    def test_direct_completion_rejects_unsupported_source_metadata(self):
        for metadata in [
            {'transfer_kind':'bounded_positive'},
            {'reference_eigenvalue':1/16},
            {'transfer_ceiling':4},
            {'transfer_kind':'bounded_positive','reference_eigenvalue':1/16,'transfer_ceiling':1/16},
            {'reference_eigenvalue':float('nan')},
        ]:
            with self.subTest(metadata=metadata):
                payload=fixture();payload.update(metadata)
                result=complete(payload)
                self.assertFalse(result['ok'])
                self.assertFalse(result['complete_finite_kernel_certificate'])
                self.assertEqual(result['verdict'],'INCONCLUSIVE_FINITE_FOCK_COMPLETION')
                self.assertIn('supports only vacuum_effect',result['reason'])

    def test_positive_exact_completion_and_JSON_matrices(self):
        payload=fixture()
        result=complete(payload,'min-trace')
        self.assertTrue(result['ok'],result)
        self.assertEqual(result['dimension'],2)
        self.assertEqual(result['inferred_occupations'],[1])
        self.assertEqual(result['exact_matrices']['T'],matrix_json(sp.diag(sp.Rational(1,2),1)))
        self.assertEqual(result['exact_checks']['all_input_moments_replayed'],len(payload['moments']))
        self.assertGreater(result['exact_checks']['multi_T_moments_replayed'],0)
        self.assertTrue(result['finite_selfadjoint_H_minus_log_T_exists'])
        self.assertTrue(result['finite_unitary_field_time_is_UCP_automorphism'])
        self.assertTrue(result['exact_checks']['quarter_fixes_Omega'])
        self.assertEqual(result['R']['i_exponents'][result['Omega_basis_index']],0)
        self.assertNotEqual(result['original_internal_W']['i_exponents'][result['Omega_basis_index']],0)
        self.assertFalse(result['flatness_claimed'])
        self.assertFalse(result['physical_source_selection_proved'])
        json.dumps(result)

    def test_minmax_both_exact_when_inputs_determine_transfer(self):
        low=complete(fixture(),'min-trace');high=complete(fixture(),'max-trace')
        self.assertTrue(low['ok'] and high['ok'])
        self.assertEqual(low['exact_matrices']['T'],high['exact_matrices']['T'])

    def test_many_T_mutation_is_rejected_despite_same_linear_constraints(self):
        payload=fixture()
        target=next(row for row in payload['moments'] if row['word']==[0,T,T,0])
        target['re']=.5
        result=complete(payload,'min-trace')
        self.assertFalse(result['ok'])
        self.assertEqual(result['verdict'],'INCONCLUSIVE_MULTI_T_REPLAY_FAILED')
        self.assertEqual(result['replay_failure']['word'],[0,T,T,0])

    def test_nonpure_occupations_are_inconclusive(self):
        payload=fixture()
        next(row for row in payload['moments'] if row['word']==[0,1])['im']=0
        result=complete(payload)
        self.assertFalse(result['ok'])
        self.assertIn('not pure',result['reason'])

    def test_input_rational_recovery_is_required(self):
        payload=fixture()
        next(row for row in payload['moments'] if row['word']==[0,T,T,0])['re']=.120715417
        result=complete(payload)
        self.assertFalse(result['ok'])
        self.assertIn('rational recovery failed',result['reason'])

    def test_CAR_moment_not_replaced_by_target_state(self):
        payload=fixture()
        next(row for row in payload['moments'] if row['word']==[0])['re']=.25
        result=complete(payload)
        self.assertFalse(result['ok'])
        self.assertIn('Pure CAR state disagrees',result['reason'])

    def test_clustering_predicate_detects_27_stationary_directions(self):
        result=complete(fixture(5),'max-trace')
        self.assertTrue(result['ok'],result)
        self.assertFalse(result['finite_clustering'])
        self.assertEqual(result['stationary_eigenspace_dimension'],27)
        self.assertEqual(result['clustering']['declared_gap_bound'],'0')

    def test_gap_constrained_minmax_are_both_clustering(self):
        payload=fixture(5);payload['cluster_gap']=.25
        for objective,expected in [('min-trace','1/4'),('max-trace','3/4')]:
            result=complete(payload,objective)
            self.assertTrue(result['ok'],result)
            self.assertTrue(result['finite_clustering'])
            self.assertEqual(result['stationary_eigenspace_dimension'],1)
            self.assertEqual(result['clustering']['declared_gap_bound'],'1/4')
            self.assertEqual(result['candidate_predictions'][0]['value']['re'],expected)
            self.assertTrue(result['finite_selfadjoint_H_minus_log_T_exists'])

    def test_floor_above_cluster_cap_is_honestly_infeasible(self):
        payload=fixture();payload['cluster_gap']=.5;payload['transfer_floor']=.75
        result=complete(payload)
        self.assertFalse(result['ok'])
        self.assertIn('Infeasible bounds',result['reason'])


if __name__=='__main__': unittest.main()
