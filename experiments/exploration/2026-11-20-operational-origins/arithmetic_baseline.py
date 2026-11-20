#!/usr/bin/env python3
"""Transparent expression census, independent of structural verdicts.

Ordered binary expression trees; atoms 1,q,n evaluated at q=2,n=5; exactly
0,1,2 binary operations from +,-,*,/,**. Division by zero is excluded.
Power requires an integer exponent in 0..8; 0**0 is excluded. Finite numeric
resource bound abs(numerator),denominator <=10**6 applies after EVERY node.
Counts are syntactic, so commuted/identity expressions count separately.
No empirical number or tolerance is used. This is NOT a p-value or an
independent rediscovery of n,q. It shows why mere integer hits are cheap.
"""
from fractions import Fraction as F
import json


def apply(op,a,b):
    if op == '+': value=a+b
    elif op == '-': value=a-b
    elif op == '*': value=a*b
    elif op == '/':
        if b == 0: return None
        value=a/b
    else:
        if b.denominator != 1 or not 0 <= b <= 8 or (a == 0 and b == 0): return None
        value=a**int(b)
    return value if abs(value.numerator) <= 10**6 and value.denominator <= 10**6 else None


def main():
    levels=[[(F(1),'1'),(F(2),'q'),(F(5),'n')]]
    for size in (1,2):
        expressions=[]
        for left_size in range(size):
            right_size=size-1-left_size
            for a,at in levels[left_size]:
                for b,bt in levels[right_size]:
                    for op in ('+','-','*','/','**'):
                        value=apply(op,a,b)
                        if value is not None: expressions.append((value,f'({at}{op}{bt})'))
        levels.append(expressions)
    return {'grammar':__doc__,'q':2,'n':5,
            'by_exact_binary_operation_count':[
                {'operations':k,'valid_expressions':len(expressions),
                 'target_hits':{str(target):{'count':sum(v==target for v,t in expressions),
                                           'examples':[t for v,t in expressions if v==target][:8]}
                                for target in (3,4,5,16,60,63,68,120,240)}}
                for k,expressions in enumerate(levels)]}


if __name__ == '__main__':
    print(json.dumps(main(),indent=2,ensure_ascii=False))
