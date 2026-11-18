"""Forgetful image in the previously classified quarter-clock moduli.

This consumes an accepted model presentation when --proposal is supplied.
Direct compilation from original relations is an explicitly labelled baseline.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

from .engine import save_json
from .mathcore import digest
from .research import Transport, Task, act, linear


def span(basis):
    values = {0}
    for v in basis:
        values |= {x ^ v for x in values}
    return values


def forgetful(task, rows):
    solved = linear(rows,task.n)
    original = linear(task.rows,task.n)
    assert solved['rank']==original['rank']
    assert all(not any(task.residual(v)) for v in solved['null_basis'])
    image_basis = [v & ((1 << 32)-1) for v in solved['null_basis']]
    dim = linear(image_basis,32)['rank']
    gauge = [sum(((1 << i)^act(task.group.A,1 << i)) << (8*j) for j in range(4)) for i in range(8)]
    gauge_rank = linear(gauge,32)['rank']

    def alpha(v):
        return sum(act(task.group.A,v >> (8*j)&255) << (8*j) for j in range(4))

    differences = [v ^ alpha(v) for v in image_basis]
    induced_rank = linear(gauge+differences,32)['rank']-gauge_rank
    q = dim-gauge_rank
    fixed = q-induced_rank
    burnside = ((1 << q)+(1 << fixed))//2
    # Separate finite orbit calculation, independent of fixed-space formula.
    values, translations = span(image_basis), span(gauge)
    assert len(values)==1 << dim
    remaining = set(values)
    orbit_sizes = Counter()
    while remaining:
        v = next(iter(remaining))
        orbit = {v ^ h for h in translations} | {alpha(v)^h for h in translations}
        assert orbit <= values
        remaining.difference_update(orbit)
        orbit_sizes[len(orbit)] += 1
    assert sum(orbit_sizes.values())==burnside
    return {'image_dimension_in_ordered_meridian_characters':dim,
            'raw_extension_characters_per_meridian_tuple':1 << (len(solved['null_basis'])-dim),
            'character_gauge_rank':gauge_rank,'quotient_dimension_before_clock_conjugation':q,
            'clock_fixed_quotient_dimension':fixed,'underlying_Gamma_classes':burnside,
            'underlying_Gamma_orbits_by_size':{str(k):v for k,v in orbit_sizes.items()},
            'proof':'Linear projection, character translations, involution Burnside formula; independently enumerated all image orbits',
            'scope':'Image under forgetting the extra equivariant lifts in the established ordered quarter-clock moduli; not mapping class quotient'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--proposal',help='Saved model run with a verified last program')
    p.add_argument('--out',default='results/research-v3/forgetful.json')
    args = p.parse_args()
    group = Transport()
    result = {'source_id':group.source.identity,'findings':{},'physical_completion':False}
    if args.proposal:
        run = json.loads(Path(args.proposal).read_text())
        name = run['identity']['task']; kind,s = name.rsplit('-',1)
        task = Task(group,kind,int(s)); program = run['attempts'][-1]['program']
        assert task.check(program)['status']=='verified'
        result['input'] = {'origin':'model_derived_presentation','path':args.proposal,'sha256':digest(Path(args.proposal).read_text())}
        result['findings'][name] = forgetful(task,task.compile(program))
    else:
        result['input'] = {'origin':'direct_original_group_word_compiler_baseline','model_contribution':False}
        for kind in ['reflection','rotation','dihedral']:
            task = Task(group,kind)
            result['findings'][task.id] = forgetful(task,task.rows)
    save_json(args.out,result); print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
