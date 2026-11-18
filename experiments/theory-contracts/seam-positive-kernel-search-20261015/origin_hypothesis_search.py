"""Source-indexed local hypothesis proposals with only preapproved checks.

This driver does not promote proposed premises to constraints and does not
claim that a static corpus inventory rechecked the underlying proofs.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
ENGINE = ROOT / 'experiments/theory-contracts/seam-equivariant-active-search-20261024'
POSITIVE = ROOT / 'experiments/theory-contracts/seam-positive-kernel-search-20261015'
PROGRAM = ROOT / 'semantic_explorer/research/program.json'
TARGET = 'problem:original-paper-source-chain-reconciliation'
TOPIC_TERMS = ('p1', 'seam', 'naht', 'quillen', 'c3', 'horizon', 'clock',
               'genesis', 'ursprung', 'information', 'entrop', 'reflection',
               'double', 'cover', 'kernel', 'determinant', 'e8')
CHECK_SCRIPT = 'experiments/theory-contracts/seam-positive-kernel-search-20261015/genesis_checks.py'
INTEGRAL_CHECK_SCRIPT = 'experiments/theory-contracts/seam-positive-kernel-search-20261015/integral_origin_checks.py'
SOURCE_CHARGE_CHECK_SCRIPT = 'experiments/theory-contracts/seam-positive-kernel-search-20261015/source_charge_checks.py'
ALLOWED = {
    'quillen-relative-normalization': 'normalized_zero_source_has_shape_response',
    'theta-heat-shape-response': 'charge_curvature_fixes_theta_shape_derivative',
    'doubled-source-counting': 'two_quarter_blocks_are_two_faces',
    'variable-mark-completion': 'primitive_family_projection_allows_higher_marks',
    'primitive-integral-completion': 'primitive_integral_projection_allows_n_other_than_4',
    'function-differential-typing': 'function_quotient_equals_three_differentials',
    'matter-word-spinor-current': 'even_matter_words_contain_bosonic_spinor_current',
    'charged-one-point': 'invariant_state_has_nonzero_charge_one_onepoint',
    'clifford-vector-spinor-typing': 'clifford_vector_is_spinor_current_under_same_spin_action',
}
CHECK_SCRIPTS = {key: (INTEGRAL_CHECK_SCRIPT if key in
                 {'primitive-integral-completion', 'function-differential-typing'} else CHECK_SCRIPT)
                 for key in ALLOWED}
for _source_check in ('matter-word-spinor-current', 'charged-one-point',
                      'clifford-vector-spinor-typing'):
    CHECK_SCRIPTS[_source_check] = SOURCE_CHARGE_CHECK_SCRIPT


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def language(value):
    if isinstance(value, dict):
        return value.get('de') or value.get('en') or ''
    return str(value or '')


def index_corpus():
    """Static full inventory; filenames/status are not validation receipts."""
    program = json.loads(PROGRAM.read_text())
    with (ROOT/'verification/script_registry.csv').open(newline='') as stream:
        registry = list(csv.DictReader(stream))
    with (ROOT/'verification/status_ledger.csv').open(newline='') as stream:
        ledger = list(csv.DictReader(stream))
    findings = []
    for problem in program['problems']:
        for finding in problem.get('findings', []):
            findings.append({'id': finding['id'], 'problem_id': problem['id'],
                             'title': language(finding.get('title')),
                             'statement': language(finding.get('statement')),
                             'status': finding.get('status'),
                             'premise_ids': finding.get('premise_ids', []),
                             'scope': language(finding.get('scope')),
                             'sources': finding.get('sources', [])})
    from registry_ast_inventory import inventory
    static_sources = inventory(ROOT)
    return {'atlas_revision': program['revision'],
            'program': program,
            'registry': registry, 'ledger': ledger, 'findings': findings,
            'registry_original_sources': static_sources,
            'canonical_research_documents': canonical_documents(),
            'scope': 'Static full index; source proofs and existing statuses are not independently revalidated.'}


def canonical_documents():
    records = []
    for name in ('PLAN.md', 'progress.json', 'RESULTS.md', 'SEARCH-METHODS.md', 'SOURCE-DICTIONARY.md', 'AUDIT.md'):
        path = ROOT/'docs/research/fundamental-plan'/name
        records.append({'path': str(path.relative_to(ROOT)), 'sha256': digest(path),
                        'content': path.read_text(), 'scope': 'Historical records preserve their original scope, not new proof receipts'})
    return records


def reuse_inventory(path):
    """Reuse parsed source only after checking every inventoried input."""
    idx = json.loads(path.read_text())
    if json.loads(PROGRAM.read_text()) != idx['program']:
        raise ValueError('Atlas changed since inventory; rebuild index')
    for name, key in [('script_registry.csv', 'registry'), ('status_ledger.csv', 'ledger')]:
        with (ROOT/'verification'/name).open(newline='') as stream:
            if list(csv.DictReader(stream)) != idx[key]:
                raise ValueError(name+' changed since inventory')
    for row in idx['registry_original_sources']['records']:
        if row['status'] != 'parsed' or digest(ROOT/row['source_path']) != row['sha256']:
            raise ValueError('Original verification source changed or was not parsed')
    # Continuation instructions can change without changing scientific sources.
    idx['canonical_research_documents'] = canonical_documents()
    return idx


def source_ids(index):
    ids = {item['id'] for item in index['findings']}
    ids.update(p['id'] for p in index['program']['problems'])
    ids.update(p['id'] for p in index['program']['premises'])
    ids.update(row['claim_id'] for row in index['ledger'] if row.get('claim_id'))
    ids.update('claim:'+row['claim_id'] for row in index['ledger'] if row.get('claim_id'))
    ids.update(row['script'] for row in index['registry'] if row.get('script'))
    ids.update('script:'+row['script'] for row in index['registry'] if row.get('script'))
    return ids


def premise_ids(index):
    # Being a known finding or a filename does not make it a premise.
    ids = {p['id'] for p in index['program']['premises']}
    ids.update(row['claim_id'] for row in index['ledger'] if row.get('claim_id'))
    ids.update('claim:'+row['claim_id'] for row in index['ledger'] if row.get('claim_id'))
    return ids


def evidence(index, test_catalog):
    target = next(p for p in index['program']['problems'] if p['id'] == TARGET)
    def score(row):
        text = (row['title'] + ' ' + row['statement'] + ' ' + row['scope']).lower()
        return sum(term in text for term in TOPIC_TERMS)
    related = sorted((x for x in index['findings'] if x['problem_id'] != TARGET and score(x)),
                     key=lambda x: (-score(x), x['id']))[:4]
    needed_ids = {sid for check in test_catalog for sid in check.get('source_ids', [])}
    needed_ids.add('finding:original-quarter-pair-fixes-square-charge-law-half-response-and-real-quillen-norm')
    needed_ids.add('finding:original-elliptic-fixed-holonomy-insertion-gives-local-operation-on-unchanged-source')
    related_ids = {x['id'] for x in related}
    related.extend(x for x in index['findings'] if x['id'] in needed_ids-related_ids)
    premises = {pid for x in related for pid in x['premise_ids']}
    # Target findings and scopes are always first. Catalogue tests are data,
    # not new source assumptions.
    selected_scripts = {'v58_seam_horizon_chain', 'v813_p1_index_kms',
                        'v189_riemann_roch_carrier', 'v471_seam_horizon_replica'}
    selected_scripts.update(sid.split(':', 1)[1] for sid in needed_ids
                            if sid.startswith('script:'))
    source_excerpts = []
    for row in index['registry_original_sources']['records']:
        if row['registry']['script'] not in selected_scripts or row['status'] != 'parsed':
            continue
        source_excerpts.append({'path': row['source_path'], 'sha256': row['sha256'],
             'docstring_excerpt': row['docstring'][:2500],
             'excerpt_complete': len(row['docstring']) <= 2500,
             'checks': row['checks'][:8], 'checks_total': len(row['checks']),
             'scope': 'Static original source excerpt, not an executed or validated proof'})
    return {'target_problem': target, 'related_findings': related,
            'original_verification_excerpts': source_excerpts,
            'premise_definitions': [p for p in index['program']['premises'] if p['id'] in premises],
            'available_tests': test_catalog,
            'counts': {'registry_entries': len(index['registry']),
                       'ledger_entries': len(index['ledger']),
                       'atlas_findings': len(index['findings']),
                       'atlas_problems': len(index['program']['problems'])},
            'atlas_revision': index['atlas_revision'], 'scope': index['scope']}


def run_checks(proposal, catalog):
    """Only built-in exact commands; even catalogue text cannot add code."""
    results = []
    for test_id in proposal['tests']:
        item = catalog[test_id]
        expected = ['python3', '-B', '-OO', CHECK_SCRIPTS[test_id], '--check', test_id]
        if item.get('command') != expected:
            raise ValueError('Curated command differs from fixed executable allowlist')
        proc = subprocess.run([sys.executable, *expected[1:]], cwd=ROOT, capture_output=True,
                              text=True, timeout=item.get('timeout_seconds', 60), check=False)
        try:
            output = json.loads(proc.stdout)
        except json.JSONDecodeError:
            raise ValueError('Curated check did not return JSON')
        if proc.returncode != 0 or output.get('id') != test_id or output.get('status') != 'pass':
            raise ValueError('Curated check failed or returned a mismatched ID')
        results.append({'test_id': test_id, 'exit_code': proc.returncode,
                        'output': output, 'stderr': proc.stderr[-3000:],
                        'meaning': item.get('question', ''),
                        'scope': item.get('scope', ''),
                        'not_physical_proof': True})
    return results


def validate_proposal_ids(proposal, known, known_premises, catalog):
    """Check typed IDs and test scopes; known IDs are not proof validation."""
    unknown_sources = [x for x in proposal['source_ids'] if x not in known]
    unknown_premises = [x for x in proposal['premises'] if x not in known_premises]
    unknown_tests = [x for x in proposal['tests'] if x not in catalog]
    if unknown_sources or unknown_premises or unknown_tests:
        raise ValueError('Unknown source/premise/test IDs: '+str(
            (unknown_sources, unknown_premises, unknown_tests)))
    declared = set(proposal['premises'])
    missing = {test: sorted(set(catalog[test].get('premise_ids', []))-declared)
               for test in proposal['tests']}
    missing = {key: value for key, value in missing.items() if value}
    if missing:
        raise ValueError('Test scope prerequisites not explicitly declared: '+str(missing))


def proposal_schema(catalog):
    strings = {'type': 'array', 'items': {'type': 'string'}}
    props = {'hypothesis': {'type': 'string'}, 'source_ids': strings,
             'premises': strings, 'additional_assumptions': strings,
             'tests': {'type': 'array', 'items': {'type': 'string', 'enum': list(catalog)}},
             'proposed_predictions': {'type': 'array', 'minItems': 1, 'items': {
                 'type': 'object', 'additionalProperties': False,
                 'required': ['check_id', 'claim_id', 'expects'],
                 'properties': {'check_id': {'type': 'string', 'enum': list(catalog)},
                    'claim_id': {'type': 'string', 'enum': [ALLOWED[k] for k in catalog]},
                    'expects': {'type': 'boolean'}}}}}
    if 'matter-word-spinor-current' in catalog:
        props['source_candidate'] = {'type': 'string', 'maxLength': 12000}
    return {'type': 'object', 'additionalProperties': False,
            'required': list(props), 'properties': props}


def scoped_predictions_match(entry):
    return (entry.get('source_candidate_assessment', {}).get('verdict') != 'typing_error'
            and all(p['matches'] for p in entry.get('prediction_checks', [])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', default='qwen3.6:35b')
    parser.add_argument('--catalog', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--rounds', type=int, default=2)
    parser.add_argument('--index-only', action='store_true')
    parser.add_argument('--reuse-inventory', type=Path)
    parser.add_argument('--question', default='Propose a new minimal source principle joining P1, double cover, information and the compiler, with an explicit mathematical rule and an independently checkable consequence. Do not merely restate the source audit or methodological corrections. Name which declared input the principle would replace, and keep unproved bridges explicit.')
    args = parser.parse_args()
    if args.rounds < 1:
        parser.error('rounds must be positive')
    if args.output.exists():
        parser.error('output already exists; use a new run path')
    catalog_raw = json.loads(args.catalog.read_text())
    # This file must be supplied and reviewed by the project owner.  Never
    # accept commands, paths or executable text from the model response.
    if catalog_raw.get('schema') != 'tfpt.genesis-checks.v1':
        raise ValueError('Unexpected curated catalogue schema')
    catalog = {entry['id']: entry for entry in catalog_raw['checks']}
    if len(catalog) != len(catalog_raw['checks']) or set(catalog) - set(ALLOWED):
        raise ValueError('Duplicate curated test ID')
    idx = reuse_inventory(args.reuse_inventory) if args.reuse_inventory else index_corpus()
    pack = evidence(idx, [{'id': key, 'meaning': item.get('question', ''),
                           'scope': item.get('scope', ''),
                           'premise_ids': item.get('premise_ids', []),
                           'source_ids': item.get('source_ids', []),
                           'allowed_claim_id': ALLOWED[key]}
                          for key, item in catalog.items()])
    result = {'schema': 'tfpt.origin-hypothesis-search.v1',
              'started_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
              'inventory': pack['counts'], 'atlas_revision': idx['atlas_revision'],
              'sources_sha256': {'program': digest(PROGRAM),
                                 'registry': digest(ROOT/'verification/script_registry.csv'),
                                 'ledger': digest(ROOT/'verification/status_ledger.csv'),
                                 'curated_catalog': digest(args.catalog)},
              'evidence': pack, 'rounds': [], 'physical_solution': False,
              'scope': 'Inventory and proposals do not constitute renewed proof verification.',
              'driver_sha256': digest(Path(__file__)),
              'checker_sha256': digest(ROOT/CHECK_SCRIPT),
              'checker_sha256_by_id': {key: digest(ROOT/CHECK_SCRIPTS[key]) for key in catalog}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    inventory_path = args.output.with_name(args.output.stem+'.inventory.json')
    with inventory_path.open('x') as stream:
        json.dump(idx, stream, indent=2, ensure_ascii=False)
        stream.write('\n')
    result['full_static_inventory'] = str(inventory_path)
    def save():
        args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n')
    if args.index_only:
        result['status'] = 'STATIC_INDEX_ONLY'
        save(); return
    sys.path.insert(0, str(ENGINE))
    sys.path.insert(0, str(POSITIVE))
    from tfpt_local.engine import Ollama, parse_response
    from adaptive_discriminator import import_hypotheses
    known = source_ids(idx)
    known_premises = premise_ids(idx)
    system = ('Work on the original TFPT source chain, preserving its complete declared P1 kernel. '
              'Use only listed known source IDs and curated test IDs. Return JSON object with fields '
              'hypothesis (string), source_ids (list of IDs), premises (list of known premise/claim IDs), '
              'additional_assumptions (list of strings), tests (list of curated test IDs), '
              'proposed_predictions (list of {check_id, claim_id, expects}, where expects is boolean). '
              'Each prediction must use the specified check ID and its listed allowed_claim_id. '
              'Do not assume an optional CAR clock condition is a universal P1 requirement. '
              'Distinguish conjecture, existing original result, scoped counterexample, and actual physical proof.')
    if 'matter-word-spinor-current' in catalog:
        system += (' Also return source_candidate: a JSON-encoded string describing the actual operator expression. '
                   'Allowed AST: {"op":"primitive","name":NAME} with NAME matter_spinor, tensor_boson, '
                   'scalar_quarter_current or clifford_vector; {"op":"product","args":[AST,...]}, '
                   '{"op":"sum","args":[AST,...]}, {"op":"adjoint","arg":AST}, '
                   '{"op":"derivative","arg":AST}, or '
                   '{"op":"proposed_sector_insertion","source_ref":"original source handle or explicit unproved proposal"}. '
                   'The last form is never an admitted original field merely because you name it. '
                   'The program checks the actual expression independently of your predictions. '
                   'Do not hide an extra sector generator in a scalar-current expression.')
    messages = [{'role': 'system', 'content': system},
                {'role': 'user', 'content': args.question+'\nEVIDENCE:\n'+json.dumps(pack, ensure_ascii=False)}]
    result['question'] = args.question
    try:
        backend = Ollama(model=args.model, schema=False, num_predict=1600,
                         response_schema=proposal_schema(catalog))
        result['backend_identity'] = backend.identity
        for iteration in range(args.rounds):
            raw, metrics = backend(messages)
            entry = {'round': iteration+1, 'raw': raw, 'metrics': metrics}
            try:
                proposal = parse_response(raw)
                if not isinstance(proposal, dict) or 'proposed_predictions' not in proposal:
                    raise ValueError('Missing proposed_predictions')
                base = {key: value for key, value in proposal.items()
                        if key not in {'proposed_predictions', 'source_candidate'}}
                import_hypotheses([base])
                predictions = proposal['proposed_predictions']
                if not isinstance(predictions, list) or not predictions:
                    raise ValueError('At least one testable prediction is required')
                for prediction in predictions:
                    if (not isinstance(prediction, dict) or
                        set(prediction) != {'check_id', 'claim_id', 'expects'} or
                        type(prediction['expects']) is not bool or
                        prediction['check_id'] not in base['tests'] or
                        ALLOWED.get(prediction['check_id']) != prediction['claim_id']):
                        raise ValueError('Prediction must use an allowed check/claim and a boolean')
                validate_proposal_ids(proposal, known, known_premises, catalog)
                assessment = None
                if 'matter-word-spinor-current' in catalog:
                    from source_charge_checks import assess_expression
                    encoded = proposal.get('source_candidate')
                    if not isinstance(encoded, str) or len(encoded) > 12000:
                        raise ValueError('A bounded JSON source_candidate expression is required')
                    assessment = assess_expression(json.loads(encoded))
                test_results = run_checks(proposal, catalog)
                checks = {x['test_id']: x['output'] for x in test_results}
                verdicts = []
                for prediction in predictions:
                    output = checks[prediction['check_id']]
                    facts = output.get('inference_claims', [])
                    claim = next((x for x in facts if x.get('id') == prediction['claim_id']), None)
                    if claim is None or type(claim.get('holds')) is not bool:
                        raise ValueError('Check lacks specified boolean inference claim')
                    verdicts.append({**prediction, 'observed': claim['holds'],
                                     'matches': claim['holds'] == prediction['expects']})
                entry.update(proposal=proposal, schema_and_ids='accepted',
                             proof_status='unverified_hypothesis',
                             exclusion_scope='manual atlas premise/strategy review required; not a lexical no-go',
                             tests=test_results, prediction_checks=verdicts)
                if assessment is not None:
                    entry['source_candidate_assessment'] = assessment
                    if assessment['verdict'] == 'excluded_scoped':
                        entry['proof_status'] = 'excluded_source_candidate_in_declared_inventory'
                    elif assessment['verdict'] == 'typing_error':
                        entry['proof_status'] = 'rejected_source_candidate_type'
                feedback = ('Curated scoped checks and prediction comparisons: '+json.dumps(verdicts)+
                            '. Evidence: '+json.dumps(test_results, ensure_ascii=False)+
                            '. Revise the proposal and predictions; a matching check does not prove the entire hypothesis.')
                if assessment is not None:
                    feedback += (' Actual source-expression decision: '+json.dumps(assessment)+
                                 '. A failed operator cannot be repaired by unchanged explanatory prose. '
                                 'A proposed sector insertion needs its actual original field binding and full locality/time/norm proof.')
            except (ValueError, TypeError, KeyError, RecursionError, subprocess.TimeoutExpired) as error:
                entry.update(schema_and_ids='rejected', error=str(error))
                feedback = 'Validation/check failed: '+str(error)+'. Return revised JSON using only IDs from the evidence packet.'
            result['rounds'].append(entry)
            messages.extend([{'role': 'assistant', 'content': raw}, {'role': 'user', 'content': feedback}])
            save()
        accepted = [x for x in result['rounds'] if x.get('schema_and_ids') == 'accepted']
        result['accepted_proposal_rounds'] = len(accepted)
        result['matching_scoped_prediction_rounds'] = sum(
            scoped_predictions_match(x) for x in accepted)
        if 'matter-word-spinor-current' in catalog:
            result['source_candidate_verdict_counts'] = {
                verdict: sum(x.get('source_candidate_assessment', {}).get('verdict') == verdict
                             for x in accepted)
                for verdict in ('excluded_scoped', 'needs_original_source', 'typing_error')}
            result['original_source_binding_proved'] = False
        result['status'] = ('PROPOSALS_AND_CURATED_CHECKS_RECORDED' if accepted
                            else 'NO_VALID_PROPOSAL_RECORDED')
    except Exception as error:
        result['status'] = 'BACKEND_OR_DRIVER_ERROR'
        result['error'] = type(error).__name__+': '+str(error)
    result['finished_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
    save()


if __name__ == '__main__':
    main()
