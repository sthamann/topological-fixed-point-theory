#!/usr/bin/env python3
"""Feed the operational test classes to the EXISTING genesis_rules.search.

Inside the repository no argument is needed. --originlab-root supports an
explicit, pinned read-only source directory. No new search engine is created.
These are necessary-constraint censuses, not claims that every abstract tuple
is realized by a physical model or a quantum code.
"""
import argparse
import json
from pathlib import Path
import sys
from sector_origin import sector_universe
from causal_origin import bracket_invariants


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--originlab-root',type=Path,
                        default=Path(__file__).resolve().parents[3]/'research'/'originlab')
    args=parser.parse_args()
    if not (args.originlab_root/'genesis_rules.py').is_file():
        parser.error('genesis_rules.py missing; supply the existing repository originlab directory')
    sys.path.insert(0,str(args.originlab_root.resolve()))
    from genesis_rules import search
    sectors=search(sector_universe(),['both_noncommutative','nonisomorphic'],
                   lambda o:o['carrier'],lambda o:(o['p'],o['q']))
    cones=[dict(d=d,clauses={'nonzero_equivariant_bracket'}
                if bracket_invariants(d)['invariant_dimension'] else set()) for d in range(2,9)]
    cone_result=search(cones,['nonzero_equivariant_bracket'],lambda o:o['d'],lambda o:o['d'])
    packing=[]
    for n in range(1,9):
        for q in range(2,7):
            clauses=set()
            if n >= 5: clauses.add('one_error_singleton_bound')
            if q**(n-1) == 1+n*(q*q-1): clauses.add('perfect_nondegenerate_packing')
            packing.append(dict(n=n,q=q,clauses=clauses))
    packing_result=search(packing,['one_error_singleton_bound','perfect_nondegenerate_packing'],
                          lambda o:o['n'],lambda o:(o['n'],o['q']))
    return {'uses_existing_search':'research/originlab/genesis_rules.py::search',
            'sector_rules':sectors,'cone_rules':cone_result,'packing_constraints':packing_result,
            'scope':'finite rule ablation; conditional minimality, not physical selection'}


if __name__ == '__main__':
    print(json.dumps(main(),indent=2,ensure_ascii=False))
