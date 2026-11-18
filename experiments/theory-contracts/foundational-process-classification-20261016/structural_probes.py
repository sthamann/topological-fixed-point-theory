"""Independent finite causal, information and composition classes.

No TFPT, quantum, complex, metric or fixed carrier data are imported.
The finite cutoffs are computational limits, not physical dimensions.
"""
from itertools import product, permutations
from pathlib import Path
import json


def canonical_relation(n,r):
    return min(tuple(sorted((p[a],p[b]) for a,b in r)) for p in permutations(range(n)))


def poset_census(n):
    pairs=[(a,b) for a in range(n) for b in range(a+1,n)]
    classes={}
    for bits in product((0,1),repeat=len(pairs)):
        r={p for p,v in zip(pairs,bits) if v}
        if any((a,c) not in r for a,b in r for bb,c in r if b==bb): continue
        classes[canonical_relation(n,r)]=r
    beta_counts={}
    for r in classes.values():
        ideals=[frozenset(i for i in range(n) if (mask>>i)&1) for mask in range(1<<n)
                if all(not ((mask>>b)&1) or ((mask>>a)&1) for a,b in r)]
        aut=[p for p in permutations(range(n)) if {(p[a],p[b]) for a,b in r}==r]
        assert all(i|j in ideals and i&j in ideals for i in ideals for j in ideals)
        irreducible=[]
        for i in ideals:
            if i and not any(j<i and k<i and j|k==i for j in ideals for k in ideals):irreducible.append(i)
        principals={frozenset({x}|{a for a,b in r if b==x}) for x in range(n)}
        assert set(irreducible)==principals
        seen=set();beta=0
        for i in ideals:
            if len(i)==n or i in seen:continue
            seen.update(frozenset(p[a] for a in i) for p in aut)
            enabled=set(range(n))-i
            enabled={x for x in enabled if all(a in i for a,b in r if b==x)}
            stab=[p for p in aut if frozenset(p[a] for a in i)==i]
            visited=set();orbits=0
            for x in enabled:
                if x not in visited:
                    orbits+=1;visited.update(p[x] for p in stab)
            beta+=orbits-1
        beta_counts[beta]=beta_counts.get(beta,0)+1
    return {'unlabeled_posets':len(classes),'invariant_Markov_weight_dimensions':beta_counts}


def canonical_function(f):
    n=len(f);best=None
    for p in permutations(range(n)):
        target=[None]*n
        for x in range(n):target[p[x]]=p[f[x]]
        val=tuple(target)
        if best is None or val<best:best=val
    return best


def sink_components(n,edges):
    reach=[[i==j or (i,j) in edges for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):reach[i][j]|=reach[i][k] and reach[k][j]
    remaining=set(range(n));result=[]
    while remaining:
        x=min(remaining);c={y for y in remaining if reach[x][y] and reach[y][x]};remaining-=c
        if all(b in c for a,b in edges if a in c):result.append(sorted(c))
    return result


def future_partition(f,o):
    partition=tuple(o)
    while True:
        signatures=[(partition[i],partition[f[i]]) for i in range(len(f))]
        unique={sig:j for j,sig in enumerate(dict.fromkeys(signatures))}
        new=tuple(unique[s] for s in signatures)
        if new==partition:return new
        partition=new


def monoid(gens):
    n=len(gens[0]);found={tuple(range(n))};changed=True
    while changed:
        changed=False
        for f in tuple(found):
            for g in gens:
                x=tuple(g[f[i]] for i in range(n))
                if x not in found:found.add(x);changed=True
    return found


def run():
    posets={n:poset_census(n) for n in range(1,5)}
    assert [posets[n]['unlabeled_posets'] for n in range(1,5)]==[1,2,5,16]
    funcs={n:len({canonical_function(f) for f in product(range(n),repeat=n)}) for n in range(1,5)}
    assert list(funcs.values())==[1,3,7,19]
    f=(1,0,2);edges={(x,f[x]) for x in range(3)}
    records=sink_components(3,edges);noisy=sink_components(3,edges|{(0,2)})
    assert records==[[0,1],[2]] and noisy==[[2]]
    assert len(set(future_partition(f,(0,0,0))))==1
    assert len(set(future_partition(f,(0,1,2))))==3
    flip=(1,0,2);reset=(0,0,2);m=monoid((flip,reset))
    rs=tuple(reset[flip[x]] for x in range(3));sr=tuple(flip[reset[x]] for x in range(3))
    assert rs!=sr and rs[0]==0 and sr[0]==1
    assert all(g[2]==2 and g[0] in (0,1) and g[1] in (0,1) for g in m)
    # All generated maps preserve two record sectors, one supports period2.
    # There is no conscious/spatial observer claim.
    return {'verdict':'FINITE_ALTERNATIVE_CLASSES_ONLY',
            'causal_posets':posets,'unlabeled_deterministic_step_worlds':funcs,
            'information_example':{'step':f,'persistent_components':records,'after_one_extra_error_edge':noisy,
                                   'constant_observation_quotient_size':1,'full_observation_quotient_size':3},
            'composition_example':{'flip':flip,'sector_preserving_reset':reset,'generated_monoid_size':len(m),
                                   'reset_after_flip':rs,'flip_after_reset':sr,'stable_record_sectors':[[0,1],[2]],
                                   'period_in_first_sector':2},
            'physical_time_metric_quantum_or_TFPT_structure_derived':False,
            'models_selected_by_matching_physical_numbers':False}


if __name__=='__main__':
    result=run();Path(__file__).with_name('structural_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
