"""Deterministic frontiers of whole, executable source-program hypotheses.

Structure-family counts describe AST construction heuristics, not independent
physical theories. Constants do not derive the proposed coefficient relations.
Only the current registered constructors are emitted; no endpoint U/Y alias,
new Hamiltonian, physical selector or state is supplied here.
"""
from hashlib import sha256
from itertools import combinations, zip_longest
import json


def leaf(name):
    return {'op': 'input', 'id': name, 'args': []}


def call(name, *args):
    return {'op': 'apply', 'id': name, 'args': list(args)}


def ast_cost(program):
    return 1 + sum(ast_cost(x) for x in program['args'])


def program_id(program):
    raw = json.dumps(program, sort_keys=True, separators=(',', ':'), allow_nan=False)
    return sha256(raw.encode()).hexdigest()


def _coefficient(which):
    c3, ratio = call('native.c3'), call('native.c6_ratio')
    if which == 0:
        return ratio
    if which == 1:
        return c3
    if which == 2:
        return call('native.scalar_multiply', call('native.scalar_i'), ratio)
    raise ValueError('unknown declared coefficient program')


def _source_coefficient(which, mark):
    coefficient = _coefficient(which)
    symbols = {'U': leaf('jU'), 'Y': leaf('jY'),
               'UY': call('native.scalar_multiply', leaf('jU'), leaf('jY'))}
    return call('native.scalar_multiply', coefficient, symbols[mark])


def _base(n, seed=0):
    body = call('native.c3') if seed == 0 else leaf('one')
    q = call('native.c6_ratio')
    if seed:
        q = call('construct.negate', q)
    action = call('native.quartic_source_program', leaf('chart_'+str(n)),
                  body, q, leaf('source_symbols'))
    operator = call('native.nil_action', action)
    goal = _finish(operator, n)
    return action, operator, body, goal


def _finish(operator, n):
    return call('native.source_jets', operator, leaf('dictionary_'+str(n)))


def _term(n, mask, coefficient):
    return call('construct.kernel_term', leaf('monomial_'+str(n)+'_'+str(mask)), coefficient)


def _mask_frontier(n):
    """Interleave off-diagonal bilinears, quartics and genuine higher vertices."""
    diagonal = {(1 << j) | (1 << (n+j)) for j in range(n)}
    old_quartic = (1 << 0) | (1 << 1) | (1 << n) | (1 << (n+1))
    groups = []
    for degree in (2, 4, 6, 8):
        masks = [sum(1 << j for j in indices)
                 for indices in combinations(range(2*n), degree)]
        masks = [m for m in masks if m not in diagonal and m != old_quartic]
        if masks:
            groups.append(masks)
    return [m for batch in zip_longest(*groups) for m in batch if m is not None]


def _altered_action(n, mask, coefficient, seed=0, second=None):
    action, _, body, parent = _base(n, seed)
    polynomial = call('construct.action_polynomial', action)
    polynomial = call('construct.kernel_add', polynomial, _term(n, mask, coefficient))
    if second is not None:
        polynomial = call('construct.kernel_add', polynomial, _term(n, second[0], second[1]))
    changed = call('native.action_from_polynomial', polynomial, body)
    return call('native.nil_action', changed), parent


def _record(program, n, parent, change_type, family, changes):
    return {
        'id': program_id(program), 'program': program, 'ast_cost': ast_cost(program),
        'change_type': change_type, 'structure_family': family,
        'parent_id': program_id(parent), 'what_changes': changes,
        'what_is_held_fixed': (
            'The same stated '+str(n)+'-mode hypothesis chart, the original Grassmann ordering, '
            'formal jU/jY coordinates, the whole ordered contraction, and the explicitly hypothetical '
            'P/cap/reference/time dictionary. P remains outside the whole marked word. '
            'No endpoint occupation is identified with the physical Higgs or Yukawa mark.'),
        'source_relation_status': 'hypothesis',
        'family_count_interpretation': 'heuristic AST structure; no mathematical independence claim',
    }


def _action_vertices():
    for index in range(128):
        for n in (2, 4):
            masks = _mask_frontier(n)
            if index >= len(masks):
                continue
            mask = masks[index]
            for coefficient in (0, 1, 2):
                operator, parent = _altered_action(n, mask, _coefficient(coefficient))
                yield _record(_finish(operator, n), n, parent, 'action_vertex',
                              'action_even_degree_'+str(mask.bit_count()),
                              'Add the specified even action monomial mask '+str(mask)+
                              ' with a coefficient program built from the existing scalar seeds; retain the entire finite exponential.')


def _marked_action_vertices():
    for index in range(48):
        for n in (2, 4):
            masks = _mask_frontier(n)
            if index >= len(masks):
                continue
            mask = masks[index]
            for mark in ('U', 'Y', 'UY'):
                operator, parent = _altered_action(n, mask, _source_coefficient(index % 3, mark))
                yield _record(_finish(operator, n), n, parent, 'marked_action_vertex',
                              'action_'+mark+'_even_degree_'+str(mask.bit_count()),
                              'Change the proposed '+mark+' source insertion through action monomial '+str(mask)+
                              '; differentiate the whole resulting kernel, including the mixed contact when both sources occur.')


def _kernel_corrections():
    for index in range(64):
        for n in (2, 4):
            masks = _mask_frontier(n)
            if index >= len(masks):
                continue
            mask = masks[index]
            _, operator, _, parent = _base(n)
            original = call('construct.kernel_of_operator', operator)
            for mark in (None, 'U', 'UY'):
                coefficient = _coefficient(index % 3) if mark is None else _source_coefficient(index % 3, mark)
                changed = call('construct.kernel_add', original, _term(n, mask, coefficient))
                program = _finish(call('native.even_kernel', changed), n)
                yield _record(program, n, parent, 'whole_kernel_correction',
                              'kernel_'+(mark or 'unmarked')+'_degree_'+str(mask.bit_count()),
                              'Add a whole-kernel coefficient at mask '+str(mask)+
                              (' with formal '+mark+' dependence.' if mark else ', rather than changing only a scalar quartic strength.'))


def _coupled_vertices():
    for index in range(32):
        for n in (2, 4):
            masks = _mask_frontier(n)
            if len(masks) < 2:
                continue
            first, second = masks[index % len(masks)], masks[(index+1) % len(masks)]
            for mark in (None, 'UY'):
                coefficient = _coefficient(index % 3)
                other = (_coefficient((index+1) % 3) if mark is None
                         else _source_coefficient((index+1) % 3, mark))
                operator, parent = _altered_action(n, first, coefficient, second=(second, other))
                yield _record(_finish(operator, n), n, parent, 'coupled_action_vertices',
                              'two_even_vertices_'+str(first.bit_count())+'_'+str(second.bit_count())+
                              ('_mixed_source' if mark else ''),
                              'Retain both action vertices '+str(first)+' and '+str(second)+
                              ' inside the same untruncated Grassmann exponential; no independent normalization of the terms.')


def _whole_sewing():
    for index in range(32):
        for n in (2, 4):
            masks = _mask_frontier(n)
            if not masks:
                continue
            A, parent = _altered_action(n, masks[index % len(masks)], _coefficient(index % 3))
            B, _ = _altered_action(n, masks[(index+1) % len(masks)], _coefficient((index+1) % 3))
            _, C, _, _ = _base(n, seed=1)
            operations = (
                ('ordered_AB', call('native.sew', A, B)),
                ('ordered_BA', call('native.sew', B, A)),
                ('reflected_A', call('native.sew', call('native.adjoint', A), A)),
                ('ordered_ABC', call('native.sew', call('native.sew', A, B), C)),
                ('ordered_ACB', call('native.sew', call('native.sew', A, C), B)),
                ('sandwiched_B', call('native.sew', call('native.sew', call('native.adjoint', A), B), A)),
                ('ordered_ABAB', call('native.sew', call('native.sew', A, B), call('native.sew', A, B))),
                ('reflected_AB', call('native.sew', call('native.adjoint', call('native.sew', A, B)), call('native.sew', A, B))),
            )
            for family, operator in operations:
                yield _record(_finish(operator, n), n, parent, 'whole_ordered_sewing', family,
                              'Compose complete weighted operators in the stated '+family+
                              ' order using the original Berezin sewing. Every intermediate state remains in the whole word; no semigroup or physical time is inferred.')


def iter_constructions(limit=256, max_cost=64, offset=0):
    """Yield a bounded, deduplicated, family-interleaved whole-program frontier."""
    if type(limit) is not int or limit < 0 or type(offset) is not int or offset < 0:
        raise ValueError('finite nonnegative limit and offset required')
    if type(max_cost) is not int or max_cost < 1:
        raise ValueError('positive finite AST cost required')
    streams = [_action_vertices(), _marked_action_vertices(), _kernel_corrections(),
               _coupled_vertices(), _whole_sewing()]
    seen = set()
    skipped = emitted = 0
    for rows in zip_longest(*streams):
        for row in rows:
            if row is None or row['ast_cost'] > max_cost or row['id'] in seen:
                continue
            seen.add(row['id'])
            if skipped < offset:
                skipped += 1
                continue
            if emitted >= limit:
                return
            yield row
            emitted += 1


def make_batch(size=256, max_cost=64, offset=0):
    return list(iter_constructions(size, max_cost, offset))


def batch_summary(batch):
    families = {}
    for row in batch:
        families[row['structure_family']] = families.get(row['structure_family'], 0)+1
    return {'whole_candidate_programs': len(batch),
            'structure_family_counts': families,
            'counting_scope': 'heuristic AST structure classes, not independent theories',
            'maximum_AST_cost': max((r['ast_cost'] for r in batch), default=0),
            'source_relation_status': 'hypotheses; original whole-source acceptance is separate'}


def main(argv=None):
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', type=int, default=256)
    parser.add_argument('--max-cost', type=int, default=64)
    parser.add_argument('--offset', type=int, default=0)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args(argv)
    batch = make_batch(args.size, args.max_cost, args.offset)
    summary = batch_summary(batch)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps({'summary': summary, 'candidates': batch}, indent=2)+'\n')
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
