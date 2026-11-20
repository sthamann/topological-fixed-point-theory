#!/usr/bin/env python3
"""EXHAUSTIVE matched-complexity null for an AME(5,2) 2->3 tensor.
Fix the same 16 even-parity computational basis strings as the stabilizer state,
with equal amplitudes +/-1/4. Exhaustively enumerate all 2^15 relative sign
assignments and test the ten exact 2-site reduced density matrices.
This is deliberately CONDITIONAL on parity support and equal amplitudes.
"""
from itertools import combinations
import json
from sympy import Rational
from five_qubit_qec import GENS, exact_stabilizer_state


def bit_index(idx, positions):
    bits=tuple((idx>>(4-i))&1 for i in range(5))
    return sum(bits[p]<<(len(positions)-1-j) for j,p in enumerate(positions))


def constraints():
    support=[b for b in range(32) if b.bit_count()%2==0]
    pos={b:i for i,b in enumerate(support)}
    cuts=[]
    for A in combinations(range(5),2):
        C=tuple(i for i in range(5) if i not in A)
        m={(bit_index(b,A),bit_index(b,C)):pos[b] for b in support}
        checks=[]
        for a1,a2 in combinations(range(4),2):
            pairs=[(m[(a1,c)],m[(a2,c)]) for c in range(8)
                   if (a1,c) in m and (a2,c) in m]
            if pairs:checks.append(tuple(pairs))
        assert len(checks)==2 and all(len(check)==4 for check in checks)
        cuts.append(tuple(checks))
    return support,tuple(cuts)


def check_constraint(signbits, pairs):
    balance=0
    for i,j in pairs:
        balance += 1 if (((signbits>>i)^(signbits>>j))&1)==0 else -1
    return balance==0


def test():
    support,cuts=constraints()
    assert len(support)==16 and len(cuts)==10
    psi=exact_stabilizer_state(GENS)
    signbits=sum((1<<i) for i,b in enumerate(support) if psi[b]==-Rational(1,4))
    if signbits&1:
        signbits ^= 0xFFFF # quotient global sign
    assert all(check_constraint(signbits,pairs) for cut in cuts for pairs in cut)
    passing=0
    examples=[]
    # position 0 sets global sign +1; count 32768 relative choices.
    for bits in range(0,1<<16,2):
        ncuts=0
        for cut in cuts:
            if all(check_constraint(bits,pairs) for pairs in cut):
                ncuts+=1
            else:
                break
        if ncuts==10:
            passing+=1
            if len(examples)<6: examples.append(f'{bits:016b}')
    return {'null_description':'all real flat +/-1/4 states on even-parity support of 16 out of 32 qubit strings; global sign fixed',
            'null_sign_assignments':1<<15,
            'all_ten_cuts_exactly_AME':passing,
            'conditional_success_fraction':f'{passing}/{1<<15}',
            'TFPT_free_parameters_used':0,
            'example_sign_patterns':examples,
            'conditioned_structure':'choice of 5 binary qubits, fixed even-parity support, equal-magnitude amplitudes'}


if __name__=='__main__':print(json.dumps(test(),indent=2,ensure_ascii=False))
