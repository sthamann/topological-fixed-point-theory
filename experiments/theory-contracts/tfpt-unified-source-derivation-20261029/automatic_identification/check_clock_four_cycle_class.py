"""Exact full finite class gate for Claude9's SM-commuting family C4 lift.

The unchanged native helper performs its own old root/flag checks on import.
Those checks are reused, not a new scientific result of this application.
"""
from pathlib import Path
from collections import Counter
from itertools import product
import importlib.util
import json

ROOT=Path(__file__).resolve().parents[4]
path=ROOT/'experiments/theory-contracts/e-series-joint-source-20261023/compiler_closure.py'
spec=importlib.util.spec_from_file_location('existing_native_compiler',path)
old=importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def require(ok,message):
    if not ok:
        raise RuntimeError(message)


def run():
    rows=[]
    for a in (1,3,5,7):
        for b in range(4):
            # Four times the actual torus parameter. Family eigenphases are
            # (1,3,-3,-1)/8 on the four tetrahedral Spin6 weights.
            av=(a,a,a,b,b,2,-1,0)
            require(sum(av)%2==0, '4v belongs to actual E8 lattice')
            c=[8,0,0,0]
            for r in old.ROOTS:
                d=sum(x*y for x,y in zip(av,r))
                require(d%2==0,'quarter phase integrality')
                c[(d//2)%4]+=1
            want=[64,64,56,64] if b%2==0 else [70,56,66,56]
            require(c==want,'full exact E8 weight-one character')
            rows.append({'a':a,'b':b,'char':c,
                         'matches_native_J':c==[60,64,60,64]})
    require(not any(x['matches_native_J'] for x in rows),
            'none in the native full field conjugacy class')
    out={'scope':'Given fixed Spin10/SU4 embedding, SM-commuting Spin10 torus with fourth power minus-one and one SU4 four-cycle lift. Complete covering representatives; no preferred field conjugator or physical Yukawa dictionary.',
         'family_four_cycle_torus_four_v':[2,-1,0],
         'candidates':rows,'native_matches':0,
         'root_arithmetic':'exact doubled native E8 roots',
         'complete_class_argument':'CLOCK-CLASS-PROOF.txt'}
    Path(__file__).with_name('clock-four-cycle-class-test.json').write_text(
        json.dumps(out,indent=2)+'\n')
    histogram=Counter()
    matches=[]
    d5=[tuple(x//2 for x in r[:5]) for r in old.INTEGER_ROOTS
        if all(x==0 for x in r[5:])]
    for n0 in product(range(4),repeat=5):
        if sum(n0)%2==0:
            continue
        for central in (0,1):
            n=list(n0)
            n[0]+=4*central
            av=tuple(n)+(2,-1,0)
            c=[8,0,0,0]
            for r in old.ROOTS:
                d=sum(x*y for x,y in zip(av,r))
                require(d%2==0,'general lift phase integral')
                c[(d//2)%4]+=1
            histogram[tuple(c)]+=1
            if c==[60,64,60,64]:
                require(sum(x%2 for x in n)==3 and
                        sorted(x%4 for x in n if x%2==0)==[0,2],
                        'native class matching condition')
                fixed=[r for r in d5 if sum(x*y for x,y in zip(n,r))%4==0]
                require(len(fixed)==6,'Spin10 centralizer has six roots')
                # All fixed roots use exactly the three odd planes, and have
                # the A2 sum/difference pattern. Together with Cartan rank5
                # they give A2 plus three abelian directions.
                require(all(all(r[j]==0 for j in range(5) if n[j]%2==0)
                            for r in fixed),'fixed roots confined to three planes')
                matches.append(n)
    require(sum(histogram.values())==1024 and len(matches)==320,
            'whole general Spin10 lift class')
    general={'scope':'All Spin10 torus lifts with g5^4=-1 and fixed SU4 four-cycle conjugacy class, before fixed SM commutation. Complete covering representatives modulo4D5; original root helper.',
             'representatives':1024,
             'character_histogram':{str(k):v for k,v in sorted(histogram.items())},
             'native_character_matches':320,
             'matching_condition':'Exactly3 odd Spin vector planes and2 even planes with residues0,2',
             'matching_Spin10_centralizer_root_count':6,
             'Spin10_fixed_semisimple_algebra_for_matches':'A2 only, plus3 abelian directions',
             'not_a_new_original_source_selection':True}
    Path(__file__).with_name('clock-general-spin10-class.json').write_text(
        json.dumps(general,indent=2)+'\n')
    print(json.dumps({'covering_representatives':len(rows),
                      'characters':sorted(set(tuple(x['char']) for x in rows)),
                      'native_matches':0,'general_representatives':1024,
                      'general_native_matches':320}))


if __name__=='__main__':
    run()
