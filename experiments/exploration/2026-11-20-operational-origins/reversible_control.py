#!/usr/bin/env python3
"""Exact negative control for the smallest reversible binary interactions.

This is about gate arity, not particle families or spatial dimension. No
lattice or clock geometry is inferred. Enumerates ALL permutations at n<=3.
"""
from itertools import permutations
from math import prod
import json


def is_affine(perm,n):
    offset=perm[0]
    columns=[perm[1<<i]^offset for i in range(n)]
    for x in range(1<<n):
        predicted=offset
        for i,column in enumerate(columns):
            if (x>>i)&1: predicted ^= column
        if predicted != perm[x]: return False
    return True


def main():
    rows=[]
    for n in (1,2,3):
        total=affine=0
        for p in permutations(range(1<<n)):
            total += 1
            affine += is_affine(p,n)
        formula=(1<<n)*prod((1<<n)-(1<<i) for i in range(n))
        assert affine == formula
        rows.append({'bits':n,'reversible_gates':total,'affine_gates':affine,
                     'nonlinear_gates':total-affine})
    return {'census':rows,'minimum_non_affine_gate_arity':3,
            'stage':0,'verdict':'OFFEN for causal dynamics; no TFPT structure identified',
            'scope':'binary reversible local gates, unconditioned count; not a no-go for quantum or larger rules'}


if __name__ == '__main__':
    print(json.dumps(main(),indent=2))
