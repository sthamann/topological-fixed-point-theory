"""Independent evaluation of constructed sources, separate from their claims.

Construction candidates may include hypotheses. Neither a model response,
catalogue role string nor a source hash establishes original TFPT identity.
Only trusted executable checkers can close the declared target obligations.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import sympy as s
from source_identity import original_identity_results, providers_for_repository

Q_IDS = (
    'original_joint_weighted_Z0_ZU_ZY_ZUY',
    'original_OS_history_RP',
    'actual_quantum_Ward_SD_contacts',
    'whole_physical_P_Ref_Eta_domains',
    'original_regulator_limit_order',
)
F1_IDS = ('constructed_common_values', *Q_IDS,
          'primitive_original_source_identity', 'original_marked_time_word_match')
TOE_IDS = tuple('TFPT.TOE.COMPLETE.01.T'+str(i) for i in range(1, 9))


def equality(a, b):
    d = s.Matrix(a)-s.Matrix(b) if isinstance(a, s.MatrixBase) else s.sympify(a)-s.sympify(b)
    values = list(d) if isinstance(d, s.MatrixBase) else [d]
    reduced = [s.simplify(x) for x in values]
    if all(x == 0 for x in reduced):
        return True
    if any(x.is_zero is False for x in reduced):
        return False
    return None


def common_values(candidate):
    """Recalculate all four jets and contractions from the constructed kernel."""
    needed = ['matrix', 'source_symbols', 'dictionary', 'ordered_source_forms']
    missing = [x for x in needed if x not in candidate]
    if missing:
        return {'status': 'MISSING_INPUT', 'missing_inputs': missing}
    K = s.Matrix(candidate['matrix'])
    ju, jy = candidate['source_symbols']
    d = candidate['dictionary']
    if set(('P', 'cap', 'reference_scalar', 'outer_scalar'))-set(d):
        return {'status': 'MISSING_INPUT', 'missing_inputs': ['actual_outer_dictionary']}
    P, cap = s.Matrix(d['P']), s.Matrix(d['cap'])
    if P.shape != K.shape or cap.shape != K.shape:
        return {'status': 'FAIL', 'counterexample': 'outer dictionary dimension mismatch'}
    ref, scalar = s.sympify(d['reference_scalar']), s.sympify(d['outer_scalar'])
    if ref == 0:
        return {'status': 'FAIL', 'counterexample': 'zero reference denominator'}
    at0 = {ju: 0, jy: 0}
    forms = {
        'Z0': K.subs(at0),
        'ZU': K.diff(ju).subs(at0),
        'ZY': K.diff(jy).subs(at0),
        'ZUY': K.diff(ju, jy).subs(at0),
    }
    tests, values = {}, {}
    for name, form in forms.items():
        tests[name+'_form'] = equality(form, candidate['ordered_source_forms'][name])
        values[name] = s.simplify(scalar*s.trace(cap*P*form*P)/ref)
        tests[name+'_value'] = equality(values[name], candidate[name])
    if values['Z0'] == 0:
        return {'status': 'FAIL', 'counterexample': 'zero whole-source denominator',
                'unnormalized_values': {k: str(v) for k, v in values.items()}}
    tests['nonzero_normalizer'] = (None if values['Z0'].is_zero is None else not values['Z0'].is_zero)
    connected = s.simplify(values['ZUY']/values['Z0']-
                          values['ZU']*values['ZY']/values['Z0']**2)
    tests['connected'] = equality(connected, candidate['normalized']['connected_UY'])
    state = 'FAIL' if False in tests.values() else 'NOT_PROVED' if None in tests.values() else 'PASS'
    return {'status': state, 'tests': tests,
            'unnormalized_values': {k: str(v) for k, v in values.items()},
            'connected_UY': str(connected),
            'scope': 'exact calculation of this constructed marked kernel, including its stated outer data; original identity is a separate requirement'}


class OriginalCompletion:
    """An injectable, code-bound verifier, with no fixed completion verdict.

    Extensions are trusted Python callables registered by the implementation,
    never module paths, commands or result labels supplied by a candidate.
    """
    def __init__(self, repo, target='joint_source', trusted_checkers=None, original_providers=None):
        if target not in ('joint_source', 'tfpt_complete', 'tfpt_validated'):
            raise ValueError('unknown original target')
        self.repo, self.target = Path(repo), target
        self.trusted = dict(trusted_checkers or {})
        self.original_providers = (providers_for_repository(self.repo) if original_providers is None
                                   else dict(original_providers))
        self.required_obligation_ids = F1_IDS + (TOE_IDS if target != 'joint_source' else ())
        if target == 'tfpt_validated':
            self.required_obligation_ids += ('original_independent_holdouts',)
        q = self.repo/'experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/runtime/source_adapter.py'
        self.q_path = q
        program = self.repo/'semantic_explorer/research/program.json'
        self.evidence_revision = json.loads(program.read_text())['revision']
        self.version = hashlib.sha256(Path(__file__).read_bytes()+
            Path(__file__).with_name('source_identity.py').read_bytes()+q.read_bytes()+program.read_bytes()).hexdigest()
        self.target_scope = 'original_TFPT_'+target

    def __call__(self, candidate, context=None):
        result = {'constructed_common_values': common_values(candidate)}
        binding = candidate.get('runtime_binding')
        result.update(original_identity_results(candidate, self.original_providers,
                      {'runtime_binding': binding, 'validator_context': context}))
        if binding is not None:
            # Only actual arrays/domain/entries are consumed by Q's reviewed
            # callbacks; supplied PASS labels are not used.
            spec = importlib.util.spec_from_file_location('tfpt_R_Q_checks', self.q_path.parents[1]/'checks/original_checks.py')
            q = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(q)
            registry = q.ensemble_test_registry(binding['domain'])
            # A correct independent packet must not be attached to an
            # unrelated proposed kernel. Verify the actual word embedding,
            # all four forms, cap, P and absolute scalars first.
            relation = self.runtime_relation(candidate, binding, q)
            for name in Q_IDS:
                checked = registry[name](binding['result'], binding['entries'], binding['domain'])
                if relation['status'] != 'PASS':
                    result[name] = {'status': relation['status'],
                                    'candidate_to_runtime_relation': relation,
                                    'separate_runtime_check': checked}
                else:
                    result[name] = checked
                if result[name]['status'] not in ('PASS', 'FAIL', 'MISSING_INPUT', 'NOT_PROVED', 'ERROR'):
                    original = result[name]['status']
                    result[name]['native_check_status'] = original
                    result[name]['status'] = ('NOT_PROVED' if original.startswith('UNDECIDED') else 'ERROR')
        else:
            for name in Q_IDS:
                result[name] = {'status': 'MISSING_INPUT',
                                'missing_inputs': ['actual_original_runtime_binding'],
                                'scope': 'existing original source requirement; candidate trace alone is not this check'}
        for name in self.required_obligation_ids:
            if name in self.trusted:
                result[name] = self.trusted[name](candidate)
            elif name not in result:
                result[name] = {'status': 'MISSING_INPUT',
                                'missing_inputs': ['original_executable_evidence_for_'+name]}
        for name, row in result.items():
            row['id'] = name
        # Hypothesis strings and candidate-supplied status flags have no
        # authority. Conversely, a trusted proof can close a hypothesis:
        # the presence of a hypothesis string is not a permanent False gate.
        return {'obligations': list(result.values()), 'scope': self.target_scope,
                'evidence': {
                    'computed_values': result['constructed_common_values'].get('unnormalized_values'),
                    'connected_UY': result['constructed_common_values'].get('connected_UY'),
                    'hypotheses': candidate.get('hypotheses', []),
                    'atlas_revision_used': self.evidence_revision}}

    @staticmethod
    def runtime_relation(candidate, binding, q):
        if 'candidate_to_word_map' not in binding:
            return {'status': 'MISSING_INPUT', 'missing_inputs': ['actual_candidate_to_original_word_embedding']}
        domain = binding['domain']
        payload = domain.get('original_check_inputs', {}).get(Q_IDS[0], {})
        try:
            d = q.bind_current_caller(payload, binding['result'], binding['entries'], domain)
        except q.MissingInput as exc:
            return {'status': 'MISSING_INPUT', 'missing_inputs': [str(exc)]}
        except (ValueError, KeyError, TypeError) as exc:
            return {'status': 'FAIL', 'counterexample': str(exc)}
        B = s.Matrix(binding['candidate_to_word_map'])
        K = s.Matrix(candidate['matrix'])
        if B.rows != K.rows:
            return {'status': 'FAIL', 'counterexample': 'candidate/word embedding dimensions differ'}
        dic = candidate['dictionary']
        tests = {'whole_candidate_span_retained': B.rank() == K.rows,
                 'same_word_gram': equality(B.H*B, q.matrix(d['gram_form'])),
                 'same_P': equality(B.H*s.Matrix(dic['P'])*B, q.matrix(d['physical_selector_form'])),
                 'same_cap': equality(B.H*s.Matrix(dic['cap'])*B, q.matrix(d['boundary_pairing_form'])),
                 'same_reference': equality(dic['reference_scalar'], q.scalar(d['reference_scalar'])),
                 'same_outer_scalar': equality(dic['outer_scalar'], q.scalar(d['scalar_factor']))}
        for name in ('Z0', 'ZU', 'ZY', 'ZUY'):
            tests[name+'_same_whole_form'] = equality(
                B.H*s.Matrix(candidate['ordered_source_forms'][name])*B,
                q.matrix(d['observed_retained_source_forms'][name]))
            tests[name+'_same_value'] = equality(candidate[name], q.scalar(binding['result'][name]))
        return {'status': 'FAIL' if False in tests.values() else 'NOT_PROVED' if None in tests.values() else 'PASS',
                'tests': tests, 'scope': 'actual entire finite retained word and outer-data matching; no continuum identity inferred'}
