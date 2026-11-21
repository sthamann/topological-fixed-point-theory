#!/usr/bin/env python3
"""Smallest nontrivial associativity-coherence loop, no imposed arity four or five."""
from functools import lru_cache
from itertools import combinations
import json

@lru_cache(None)
def bracketings(leaves):
    if len(leaves)==1:return (leaves[0],)
    out=[]
    for k in range(1,len(leaves)):
        out.extend((a,b) for a in bracketings(leaves[:k]) for b in bracketings(leaves[k:]))
    return tuple(out)


def rotations(node):
    if isinstance(node,str):return set()
    a,b=node
    out=set()
    if isinstance(a,tuple):
        x,y=a
        out.add((x,(y,b)))
    if isinstance(b,tuple):
        y,z=b
        out.add(((a,y),z))
    out |= {(x,b) for x in rotations(a)}
    out |= {(a,y) for y in rotations(b)}
    return out


def strtree(t):
    return t if isinstance(t,str) else '('+strtree(t[0])+strtree(t[1])+')'


def simple_paths(graph,source,target):
    if source==target:return [[source]]
    results=[]
    def rec(path):
        at=path[-1]
        if at==target:
            results.append(path)
            return
        for y in graph[at]:
            if y not in path:rec(path+[y])
    rec([source])
    return results


def test():
    all_counts={}
    for n in range(2,7):
        trees=bracketings(tuple('abcdef'[:n]))
        all_counts[n]=len(trees)
    assert [all_counts[i] for i in (2,3,4,5,6)]==[1,2,5,14,42]
    first_cycle=None
    for n in range(2,7):
        T=bracketings(tuple('abcdef'[:n]))
        G={t:rotations(t) for t in T}
        assert all(y in T for nbrs in G.values() for y in nbrs)
        vertices=len(T); edges=sum(len(nbrs) for nbrs in G.values())//2
        # Connected associahedron; a finite connected graph has a cycle iff E>=V.
        if edges>=vertices:
            first_cycle=(n,vertices,edges,G)
            break
    n,v,e,G=first_cycle
    assert (n,v,e)==(4,5,5)
    assert all(len(nbrs)==2 for nbrs in G.values())
    a=bracketings(tuple('abcd'))[0]
    b=bracketings(tuple('abcd'))[-1]
    paths=simple_paths(G,a,b)
    assert sorted(len(p)-1 for p in paths)==[2,3]
    # Exhaustive same-grammar baseline for ordered binary bracketings of arity 2..6.
    same_complexity_hits=sum(cnt==5 for cnt in all_counts.values())
    assert same_complexity_hits==1
    return {'bracketing_counts':all_counts,
            'first_cycle_arity':n,'first_cycle_vertices':v,'first_cycle_edges':e,
            'first_cycle_degree':2,
            'path_lengths_between_extreme_parenthesizations':[2,3],
            'five_bracketings':[strtree(t) for t in G],
            'null_arity_range':'2..6',
            'null_equal_five_bracketings':f'{same_complexity_hits}/{len(all_counts)}',
            'interpretation_limit':'2 and 3 count associator steps, not physical SU3/SU2 sites'}


if __name__=='__main__':print(json.dumps(test(),indent=2,ensure_ascii=False))
