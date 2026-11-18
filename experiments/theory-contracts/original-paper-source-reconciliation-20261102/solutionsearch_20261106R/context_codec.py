"""Readable, reversible model-context packing; no scientific text truncation.

Full provenance remains in the locally saved original-evidence.json. References
are data, not proof attestations. History dictionaries use explicit JSON pools,
not binary compression, summaries, gzip or executable deserialization.
"""
from collections import Counter
from copy import deepcopy
import hashlib
import json


PROVENANCE_FIELDS = ('source_files', 'source_pool', 'local_attachments')


def canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def _handle_manifest(value):
    if not isinstance(value, dict):
        return None
    ids = list(value)
    # Exact range notation for generated consecutive IDs, else explicit IDs.
    for prefix in ('s', 'f'):
        if ids == [prefix + str(i + 1) for i in range(len(ids))]:
            return {'prefix': prefix, 'first': 1, 'last': len(ids),
                    'rule': 'Every integer in the inclusive interval, no omissions.'}
    return {'ids': ids}


def model_knowledge_packet(fullpack, archive_path='original-evidence.json'):
    """Move only raw provenance indices outside the model prompt.

    The caller must persist fullpack at archive_path before transmitting the
    model view. Existing finding source handles still point to that archive.
    Every other packet value is unchanged.
    """
    if 'local_provenance_archive' in fullpack:
        raise ValueError('Knowledge packet already contains codec metadata')
    result = {k: deepcopy(v) for k, v in fullpack.items() if k not in PROVENANCE_FIELDS}
    fields = {}
    for key in PROVENANCE_FIELDS:
        if key not in fullpack:
            continue
        value = fullpack[key]
        fields[key] = {'json_pointer': '/' + key, 'sha256_canonical_json': digest(value),
                       'count': len(value), 'handle_manifest': _handle_manifest(value)}
    result['local_provenance_archive'] = {
        'schema': 'tfpt.local-provenance-index.v1', 'path': str(archive_path),
        'full_packet_sha256_canonical_json': digest(fullpack), 'fields': fields,
        'interpretation': ('Die vollständigen Quellenindizes bleiben unverändert im lokalen Archiv. '
                           'source_handle verweist dort auf source_pool und anschließend source_files. '
                           'Wissenschaftliche Befunde weiterhin über finding_ids zitieren. '
                           'Ein Quellenhash belegt Dateifrische und keinen mathematischen Beweis.'),
    }
    return result


def expand_model_knowledge_packet(modelpack, archive_fullpack):
    archive = modelpack['local_provenance_archive']
    if digest(archive_fullpack) != archive['full_packet_sha256_canonical_json']:
        raise ValueError('Local provenance archive does not match this model packet')
    result = {k: deepcopy(v) for k, v in modelpack.items() if k != 'local_provenance_archive'}
    for key, metadata in archive['fields'].items():
        value = archive_fullpack[key]
        if digest(value) != metadata['sha256_canonical_json'] or len(value) != metadata['count']:
            raise ValueError('Local provenance field mismatch: ' + key)
        result[key] = deepcopy(value)
    if result != archive_fullpack:
        raise ValueError('Model scientific content or archive content was modified')
    return result


def _pack_program(program, nodes, node_index):
    local_ids, overrides = [], {}
    for i, instruction in enumerate(program):
        if not isinstance(instruction, dict) or set(instruction) != {'op', 'id', 'args'}:
            raise ValueError('Program is not a native JSON instruction list')
        args = instruction['args']
        if any(type(a) is not int or not 0 <= a < i for a in args):
            raise ValueError('Program arguments must reference earlier instructions')
        node = {'op': instruction['op'], 'id': instruction['id'],
                'args': [local_ids[a] for a in args]}
        key = canonical(node)
        if key not in node_index:
            nid = 'n' + str(len(nodes) + 1)
            node_index[key] = nid
            nodes[nid] = node
        nid = node_index[key]
        # The global DAG shares identical subexpressions. Overrides preserve
        # exact original indices when a flat program repeats a shared node.
        first_args = [local_ids.index(local_ids[a]) for a in args]
        if first_args != args:
            overrides[str(i)] = list(args)
        local_ids.append(nid)
    result = {'$program_nodes': local_ids}
    if overrides:
        result['argument_index_overrides'] = overrides
    return result


def _unpack_program(program, nodes):
    ids, out = program['$program_nodes'], []
    overrides = program.get('argument_index_overrides', {})
    for i, nid in enumerate(ids):
        node = nodes[nid]
        args = overrides.get(str(i))
        if args is None:
            args = [ids[:i].index(n) for n in node['args']]
        if any(type(a) is not int or not 0 <= a < i for a in args):
            raise ValueError('Invalid restored program argument')
        if [ids[a] for a in args] != node['args']:
            raise ValueError('Argument override does not match its actual DAG node')
        out.append({'op': node['op'], 'id': node['id'], 'args': list(args)})
    return out


def compact_trial_history(history):
    """Losslessly share program DAGs and repeated check/version values.

    Model reading rule: {$shared_value: rN} means the exact JSON value stored
    under shared_values[rN]. Nodes retain operation IDs and child node IDs.
    Every trial ID, computed value, failure, open obligation and version stays.
    """
    if not isinstance(history, dict) or 'whole_programs' not in history:
        raise ValueError('Expected the complete prior-trial history dictionary')
    nodes, node_index = {}, {}
    transformed = deepcopy(history)
    for trial in transformed['whole_programs']:
        if 'program' in trial:
            trial['program'] = _pack_program(trial['program'], nodes, node_index)
    trials = transformed['whole_programs']
    if trials and all(set(t) == set(trials[0]) for t in trials):
        columns = list(trials[0])
        transformed['whole_programs'] = {'$record_table': {
            'columns': columns, 'rows': [[t[c] for c in columns] for t in trials]}}
    operation_pool, operation_ids, compact_nodes = {}, {}, {}
    for nid, node in nodes.items():
        operation = {'op': node['op'], 'id': node['id']}
        key = canonical(operation)
        if key not in operation_ids:
            oid = 'o' + str(len(operation_pool) + 1)
            operation_ids[key] = oid
            operation_pool[oid] = operation
        compact_nodes[nid] = [operation_ids[key], node['args']]
    counts = Counter()

    def count(value):
        if isinstance(value, str):
            if len(value) >= 64:
                counts[canonical(value)] += 1
        elif isinstance(value, (dict, list)):
            key = canonical(value)
            if len(key) >= 80:
                counts[key] += 1
            for child in (value.values() if isinstance(value, dict) else value):
                count(child)
    count(transformed)
    shared, refs = {}, {}

    def encode(value, allow_ref=True):
        if not isinstance(value, (dict, list, str)):
            return value
        key = canonical(value)
        if allow_ref and counts[key] > 1:
            if key not in refs:
                rid = 'r' + str(len(refs) + 1)
                refs[key] = rid
                shared[rid] = encode(value, allow_ref=False)
            return {'$shared_value': refs[key]}
        if isinstance(value, str):
            return value
        if isinstance(value, list):
            return [encode(v) for v in value]
        return {k: encode(v) for k, v in value.items()}

    content = encode(transformed, allow_ref=False)
    return {
        'schema': 'tfpt.trial-history-readable-codec.v1',
        'reading_rule': ('{$shared_value:rN} is exactly shared_values[rN], recursively. '
                         'program.$program_nodes lists the complete flat instruction sequence; '
                         'program_node_pool rows follow program_node_columns and give the operation and child node IDs. '
                         'operation_ref resolves through program_operation_pool. '
                         '{$record_table:{columns,rows}} is a table: each row has one value per column, in order. '
                         'An argument normally references the first earlier occurrence of its child node; '
                         'argument_index_overrides preserves different original indices. '
                         'All trial IDs, checks, open obligations and recorded evaluation versions are retained. '
                         'These finite hypothesis results are not original physical identity proofs.'),
        'original_history_sha256_canonical_json': digest(history),
        'original_trial_count': len(history['whole_programs']),
        'program_node_columns': ['operation_ref', 'args'],
        'program_node_pool': compact_nodes, 'program_operation_pool': operation_pool,
        'shared_values': shared, 'history': content,
    }


def expand_trial_history(compact):
    if compact.get('schema') != 'tfpt.trial-history-readable-codec.v1':
        raise ValueError('Unknown history codec schema')
    shared = compact['shared_values']

    def decode(value, active=()):
        if isinstance(value, list):
            return [decode(v, active) for v in value]
        if not isinstance(value, dict):
            return value
        if set(value) == {'$shared_value'}:
            ref = value['$shared_value']
            if ref in active:
                raise ValueError('Cyclic shared history value')
            return decode(shared[ref], active + (ref,))
        return {k: decode(v, active) for k, v in value.items()}

    history = decode(compact['history'])
    trials = history['whole_programs']
    if isinstance(trials, dict) and set(trials) == {'$record_table'}:
        table = trials['$record_table']
        if any(len(r) != len(table['columns']) for r in table['rows']):
            raise ValueError('Trial record table width mismatch')
        history['whole_programs'] = [dict(zip(table['columns'], row)) for row in table['rows']]
    if compact['program_node_columns'] != ['operation_ref', 'args']:
        raise ValueError('Unknown program-node columns')
    nodes = {nid: {**compact['program_operation_pool'][node[0]],
                   'args': node[1]} for nid, node in compact['program_node_pool'].items()}
    for trial in history['whole_programs']:
        if 'program' in trial:
            trial['program'] = _unpack_program(trial['program'], nodes)
    if len(history['whole_programs']) != compact['original_trial_count'] or digest(history) != compact['original_history_sha256_canonical_json']:
        raise ValueError('Restored history does not match original content')
    return history
