"""Independent finite-matrix controls for the unknown-moment search.

The matrices here are test witnesses only.  MomentSearch never receives them.
Raw words are evaluated by matrix action, without using the word reducer.
"""
from __future__ import annotations

import unittest
from fractions import Fraction
import json
from pathlib import Path
import tempfile

import numpy as np

from exact_certificate import recover
from kernel_search import MomentSearch, normalize_fixed_moments, load_source_control
from word_algebra import T, WordAlgebra


def majoranas(m):
    """Jordan-Wigner representation, independent of the search's word algebra."""
    identity = np.eye(2, dtype=complex)
    z = np.diag([1, -1]).astype(complex)
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    result = {}
    for j in range(m):
        for offset, pauli in enumerate((x, y)):
            matrix = np.ones((1, 1), dtype=complex)
            for k in range(m):
                matrix = np.kron(matrix, z if k < j else pauli if k == j else identity)
            result[2*j + offset] = matrix
    return result


def vacuum(m):
    result = np.zeros(2**m, dtype=complex)
    result[0] = 1
    return result


def word_vector(word, matrices, reference):
    result = reference.copy()
    for token in reversed(word):
        result = matrices[token] @ result
    return result


def evaluate(word, matrices, reference):
    return np.vdot(reference, word_vector(word, matrices, reference))


def assign(search, matrices, reference):
    search.y.value = np.array([evaluate(w, matrices, reference) for w in search.words])


def violation(search):
    return max(float(np.max(c.violation())) for c in search.constraints)


def payload(m, degree, matrices, reference, max_degree=None):
    rows = []
    for word in WordAlgebra(m).words(2 * degree if max_degree is None else max_degree):
        value = evaluate(word, matrices, reference)
        rows.append({'word': list(word), 're': float(value.real), 'im': float(value.imag)})
    return {'m': m, 'degree': degree, 'transfer_floor': 0, 'moments': rows}


class KernelSearchTests(unittest.TestCase):
    def test_scaled_source_reference_preserves_all_moments_without_claiming_absolute_Z(self):
        # Declared integer vacuum shift at ell=4log2 gives lambda=1/16;
        # the integer shift leaves the full-period operator phase unchanged.
        lam=Fraction(1,16)
        matrices=majoranas(1);matrices[T]=np.diag([float(lam),float(lam/4)])
        reference=vacuum(1)
        controls=[]
        for word in WordAlgebra(1).words(4):
            value=evaluate(word,matrices,reference)
            controls.append({'word':list(word),'re':str(Fraction(float(value.real))),
                             'im':str(Fraction(float(value.imag)))})
        scaled=MomentSearch(1,2,transfer_kind='bounded_positive',transfer_ceiling=float(lam),
                            reference_eigenvalue=float(lam),fixed_moments=controls)
        assign(scaled,matrices,reference)
        self.assertLess(violation(scaled),1e-12)
        self.assertEqual(scaled.moment((T,T)).value,float(lam**2))
        self.assertEqual(scaled.moment((T,0,T,0)).value,float(lam**2/4))
        self.assertEqual(scaled.moment((0,T,0,T)).value,float(lam**2/4))
        old=MomentSearch(1,2,transfer_kind='bounded_positive',fixed_moments=controls)
        assign(old,matrices,reference)
        self.assertGreater(violation(old),.9)
        result=scaled.solve(solver='SCS',tolerance=1e-8)
        self.assertIn(result['solver_status'],('optimal','optimal_inaccurate'))
        self.assertEqual(old.solve(solver='SCS')['solver_status'],'infeasible')
        self.assertEqual(result['reference_eigenvalue'],float(lam))
        self.assertTrue(result['conditional_normalized_state'])
        self.assertFalse(result['absolute_partition_function_recovered'])
        self.assertIn('not a recovered Z(0)',result['state_normalization_scope'])
        self.assertIn('not a value derived from P1',result['reference_eigenvalue_scope'])

    def test_reference_eigenvalue_and_general_bounds_validation(self):
        for eigenvalue in (0,-1,float('nan'),float('inf')):
            with self.assertRaises(ValueError):MomentSearch(1,2,reference_eigenvalue=eigenvalue)
        for lower,upper,eigenvalue in ((.1,.2,.05),(0,.1,.2),(float('nan'),1,1),(0,float('inf'),1)):
            with self.assertRaises(ValueError):MomentSearch(1,2,lower,transfer_ceiling=upper,
                reference_eigenvalue=eigenvalue,transfer_kind='bounded_positive')
        with self.assertRaises(ValueError):MomentSearch(1,2,reference_eigenvalue=.5)
        # Neither bound is tied to 1 in the declared general bounded class.
        search=MomentSearch(1,2,1.5,transfer_ceiling=3,reference_eigenvalue=2,transfer_kind='bounded_positive')
        self.assertEqual(search.reference_eigenvalue,2)

    def test_optional_pure_graph_effect_allows_product_below_single_response(self):
        # OPTIONAL pure-graph gapped-parent/effect control, not thermal GNS time.
        # Exterior transfer of K on the declared two-mode Fock realization.
        # The reference is the occupied graph orbital, not a guessed basis state.
        t=1/4
        graph=np.array([[16,4],[4,1]],dtype=float)/17
        one_particle=t*(np.eye(2)-graph)+(1/t)*graph
        transfer=np.zeros((4,4),dtype=complex)
        transfer[0,0]=t;transfer[1:3,1:3]=t*one_particle
        transfer[3,3]=t*np.linalg.det(one_particle)
        matrices=majoranas(2);matrices[T]=transfer
        reference=np.array([0,4,1,0],dtype=complex)/np.sqrt(17)
        self.assertAlmostEqual(evaluate((0,T,0),matrices,reference).real,1/4)
        pair_word=(0,2,T,2,0)
        self.assertAlmostEqual(evaluate(pair_word,matrices,reference).real,1/16)
        for floor in (0,1/4):
            search=MomentSearch(m=2,degree=3,transfer_floor=floor)
            assign(search,matrices,reference)
            if floor==0:self.assertLess(violation(search),1e-9)
            else:self.assertGreater(violation(search),.18)
        controls=[]
        for word in [(j,T,j) for j in range(4)]+[pair_word]:
            value=evaluate(word,matrices,reference)
            controls.append({'word':list(word),'re':str(Fraction(float(value.real)).limit_denominator(1024)),
                             'im':str(Fraction(float(value.imag)).limit_denominator(1024))})
        accepted=MomentSearch(2,3,0,fixed_moments=controls,source_observables={'label':'optional pure-graph gapped-parent/effect control; not universal collar/KMS source acceptance'}).solve(solver='SCS',tolerance=1e-7)
        rejected=MomentSearch(2,3,.25,fixed_moments=controls).solve(solver='SCS',tolerance=1e-7)
        self.assertIn(accepted['solver_status'],('optimal','optimal_inaccurate'))
        self.assertEqual(rejected['solver_status'],'infeasible')

    def test_thermal_GNS_transfer_needs_source_ceiling4_not_contraction(self):
        # Left CAR action on Hilbert-Schmidt space and sqrt(rho), rho=diag(16,1)/17.
        matrices={j:np.kron(g,np.eye(2)) for j,g in majoranas(1).items()}
        matrices[T]=np.diag([1,4,1/4,1])
        reference=np.array([4,0,0,1],dtype=complex)/np.sqrt(17)
        np.testing.assert_allclose(matrices[T]@reference,reference)
        self.assertEqual(sorted(np.linalg.eigvalsh(matrices[T])),[1/4,1,1,4])
        self.assertAlmostEqual(evaluate((0,T,0),matrices,reference).real,8/17)
        self.assertAlmostEqual(evaluate((0,T,T,0),matrices,reference).real,1)
        controls=[]
        for word in WordAlgebra(1).words(4):
            value=evaluate(word,matrices,reference)
            controls.append({'word':list(word),'re':str(Fraction(float(value.real)).limit_denominator(65536)),
                             'im':str(Fraction(float(value.imag)).limit_denominator(65536))})
        for ceiling in (1,4):
            search=MomentSearch(1,2,transfer_kind='thermal_gns',transfer_ceiling=ceiling,fixed_moments=controls)
            assign(search,matrices,reference)
            if ceiling==4:self.assertLess(violation(search),1e-9)
            else:self.assertGreater(violation(search),.5)
            result=search.solve(solver='SCS',tolerance=1e-7)
            self.assertEqual(result['transfer_ceiling'],ceiling)
            self.assertIn('class mismatch',result['transfer_scope'])
            if ceiling==4:self.assertIn(result['solver_status'],('optimal','optimal_inaccurate'))
            else:self.assertEqual(result['solver_status'],'infeasible')

    def test_transfer_type_and_bound_must_be_explicit_and_finite(self):
        for ceiling in (float('nan'),float('inf'),.9):
            with self.assertRaises(ValueError):MomentSearch(1,2,transfer_ceiling=ceiling)
        with self.assertRaises(ValueError):MomentSearch(1,2,transfer_ceiling=4)
        with self.assertRaises(ValueError):MomentSearch(1,2,transfer_kind='invented')

    def test_feasibility_default_does_not_force_binary_occupations(self):
        result=MomentSearch(1,2).solve(solver='SCS')
        self.assertEqual(result['objective'],'feasibility')
        self.assertEqual(result['transfer_floor'],0)
        self.assertAlmostEqual(result['occupations'][0],.5,places=5)
        self.assertIn('feasibility_only',result['objective_role'])

    def test_safe_fixed_moment_API_preserves_CAR_sign_and_provenance(self):
        rows=[{'word':[1,0],'re':'0','im':'1/3'}]
        normalized=normalize_fixed_moments(1,2,rows)
        self.assertEqual(normalized,[{'word':[0,1],'re':'0','im':'-1/3'}])
        with self.assertRaises(ValueError):normalize_fixed_moments(1,2,rows+[{'word':[0,1],'re':0,'im':'1/3'}])
        with self.assertRaises(ValueError):normalize_fixed_moments(1,2,[{'word':[],'re':'__import__("os")'}])
        with self.assertRaises(ValueError):normalize_fixed_moments(1,2,[{'word':[-1]*5,'re':0}])
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'source.json';path.write_text(json.dumps({'fixed_moments':rows,'source_observables':{'label':'source control'}}))
            control=load_source_control(path)
            self.assertEqual(control['source_observables']['label'],'source control')
            self.assertEqual(len(control['source_observables']['input_sha256']),64)

    def test_five_mode_supplied_kernel_matches_full_complex_gram(self):
        search = MomentSearch(m=5, degree=2)
        matrices = majoranas(5)
        # Rational diagonal contraction: one fixed vacuum and distinct excitations.
        matrices[T] = np.diag([1] + [k / 32 for k in range(1, 32)])
        reference = vacuum(5)
        assign(search, matrices, reference)
        self.assertLess(violation(search), 1e-9)
        vectors = np.column_stack([word_vector(w, matrices, reference) for w in search.basis])
        expected = vectors.conj().T @ vectors
        np.testing.assert_allclose(search.gram.value, expected, atol=1e-12)
        np.testing.assert_allclose(search.gram.value, search.gram.value.conj().T, atol=1e-12)
        self.assertGreaterEqual(np.linalg.eigvalsh(expected)[0], -1e-12)

    def test_imaginary_cross_moment_survives_quarter_invariance(self):
        search = MomentSearch(m=1, degree=2)
        matrices = majoranas(1)
        matrices[T] = np.diag([1, 0.5])
        assign(search, matrices, vacuum(1))
        self.assertAlmostEqual(search.moment((0, 1)).value, 1j)
        self.assertAlmostEqual(search.moment((1, 0)).value, -1j)
        self.assertAlmostEqual(search.moment((0, T, 1)).value, 0.5j)
        self.assertLess(violation(search), 1e-10)

    def test_stationarity_does_not_erase_internal_noncommuting_transfer(self):
        search = MomentSearch(m=1, degree=2)
        matrices = majoranas(1)
        matrices[T] = np.diag([1, 0.25])
        assign(search, matrices, vacuum(1))
        self.assertGreater(np.linalg.norm(matrices[T] @ matrices[0] - matrices[0] @ matrices[T]), 0.5)
        self.assertAlmostEqual(search.moment((T, 0, 0)).value, 1)
        self.assertAlmostEqual(search.moment((0, T, 0)).value, 0.25)
        self.assertAlmostEqual(search.moment((0, T, T, 0)).value, 0.0625)
        self.assertLess(violation(search), 1e-10)

    def test_positive_contraction_without_fixed_reference_is_rejected(self):
        search = MomentSearch(m=1, degree=2)
        matrices = majoranas(1)
        matrices[T] = np.diag([0.75, 0.5])
        assign(search, matrices, vacuum(1))
        # All positivity conditions hold; this must fail the fixed-reference condition.
        self.assertGreaterEqual(np.linalg.eigvalsh(search.gram.value)[0], -1e-12)
        self.assertGreaterEqual(np.linalg.eigvalsh(search.transfer.value)[0], -1e-12)
        self.assertGreaterEqual(np.linalg.eigvalsh(search.lower_gram.value-search.transfer.value)[0], -1e-12)
        self.assertGreater(violation(search), 0.2)

    def test_transfer_outside_each_bound_is_rejected_despite_positive_gram(self):
        for excited_eigenvalue in (-0.25, 1.25):
            with self.subTest(excited_eigenvalue=excited_eigenvalue):
                search = MomentSearch(m=1, degree=2)
                matrices = majoranas(1)
                matrices[T] = np.diag([1, excited_eigenvalue])
                assign(search, matrices, vacuum(1))
                self.assertGreaterEqual(np.linalg.eigvalsh(search.gram.value)[0], -1e-12)
                self.assertGreater(violation(search), 0.2)

    def test_positive_transfer_floor_is_an_actual_constraint(self):
        search = MomentSearch(m=1, degree=2, transfer_floor=0.5)
        matrices = majoranas(1)
        matrices[T] = np.diag([1, 0.25])
        assign(search, matrices, vacuum(1))
        self.assertGreater(violation(search), 0.2)

    def test_contraction_localizer_detects_transfer_leaking_out_of_lower_span(self):
        search = MomentSearch(m=1, degree=2)
        # Two copies of the CAR module.  The previous level sees only e0,e1;
        # T couples e1 to e3.  Its compression is a contraction, but T is not.
        matrices = {j: np.kron(np.eye(2), matrix) for j, matrix in majoranas(1).items()}
        matrices[T] = np.diag([1, 0.5, 0.5, 0.5]).astype(complex)
        matrices[T][1, 3] = matrices[T][3, 1] = 0.75
        reference = np.array([1, 0, 0, 0], dtype=complex)
        assign(search, matrices, reference)
        self.assertGreaterEqual(np.linalg.eigvalsh(search.gram.value)[0], -1e-12)
        self.assertGreaterEqual(np.linalg.eigvalsh(search.transfer.value)[0], -1e-12)
        self.assertGreaterEqual(np.linalg.eigvalsh(search.lower_gram.value-search.transfer.value)[0], -1e-12)
        self.assertAlmostEqual(search.moment((0, T, 0)).value, 0.5)
        self.assertAlmostEqual(search.moment((0, T, T, 0)).value, 0.8125)
        self.assertLess(np.linalg.eigvalsh(search.transfer.value-search.transfer_squared.value)[0], -0.3)
        self.assertGreater(violation(search), 0.3)

    def test_noninvariant_reference_is_rejected_even_with_identity_transfer(self):
        search = MomentSearch(m=1, degree=2)
        matrices = majoranas(1)
        matrices[T] = np.eye(2)
        reference = np.array([1, 1], dtype=complex) / np.sqrt(2)
        assign(search, matrices, reference)
        self.assertGreaterEqual(np.linalg.eigvalsh(search.gram.value)[0], -1e-12)
        self.assertAlmostEqual(search.moment((T,)).value, 1)
        self.assertGreater(violation(search), 0.9)

    def test_nonhermitian_moments_are_not_hidden_by_psd_symmetrization(self):
        search = MomentSearch(m=1, degree=2)
        matrices = majoranas(1)
        matrices[T] = np.diag([1, 0.5])
        assign(search, matrices, vacuum(1))
        values = search.y.value.copy()
        values[search.index[(0,)]] = 0.2j
        search.y.value = values
        self.assertGreater(violation(search), 0.1)

    def test_flat_truncated_model_requires_independent_extension_data(self):
        matrices = majoranas(1)
        matrices[T] = np.diag([1, 0.5])
        search = MomentSearch(m=1, degree=2)
        assign(search, matrices, vacuum(1))
        self.assertEqual(np.linalg.matrix_rank(search.lower_gram.value), 2)
        self.assertEqual(np.linalg.matrix_rank(search.gram.value), 2)
        result = recover(payload(1, 2, matrices, vacuum(1)))
        self.assertFalse(result['ok'])
        self.assertEqual(result['verdict'], 'EXACT_RECOVERY_FAILED')
        self.assertIn('missing recovered moment', result['error'])

    def test_positive_but_nonflat_kernel_is_not_accepted(self):
        # A faithful one-mode state has a four-dimensional GNS space.  At d=2
        # the previous level spans only {I,X,Y}, so a flat claim is premature.
        matrices = majoranas(1)
        matrices[T] = np.eye(2)
        rows = []
        for word in WordAlgebra(1).words(6):
            operator = np.eye(2, dtype=complex)
            for token in word:
                operator = operator @ matrices[token]
            value = np.trace(operator) / 2
            rows.append({'word': list(word), 're': float(value.real), 'im': float(value.imag)})
        result = recover({'m': 1, 'degree': 2, 'transfer_floor': 0, 'moments': rows})
        self.assertFalse(result['ok'])
        self.assertIn('not represented by the recovered flat basis', result['error'])

    def test_actual_solve_returns_only_a_numerical_candidate(self):
        search = MomentSearch(m=1, degree=2)
        result = search.solve(objective='diversity', seed=3, solver='CLARABEL')
        self.assertIn(result['solver_status'], ('optimal', 'optimal_inaccurate'))
        self.assertEqual(result['verdict'], 'NUMERICAL_PSEUDOMOMENT_CANDIDATE_REQUIRES_EXACT_RECOVERY')
        self.assertLess(result['max_constraint_violation'], 1e-6)
        self.assertGreaterEqual(min(result['min_eigenvalues'].values()), -1e-6)
        self.assertIn('not P1/P2 consequences', result['selection_scope'])
        self.assertIn('physical chiral holonomy', result['not_assumed'])
        self.assertFalse(recover(result)['ok'])


if __name__ == '__main__':
    unittest.main()
