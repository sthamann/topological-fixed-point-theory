#!/usr/bin/env python3
"""Adversarial control: a completely non-correcting four-check qubit system
has the SAME abstract F2^4 syndrome space and therefore yields the SAME
Hamming-to-E8 construction after choosing a hyperplane and parity extension.
This rules out treating that E8 step as *specific* evidence for the perfect
five-qubit code; it remains a valid but non-discriminating mathematical bridge.
"""
from collections import Counter
import json
from five_qubit_qec import code_data, xz, syndrome
from syndrome_geometry import hyperplane, hamming_extension


def test():
    controls=['ZIIII','IZIII','IIZII','IIIZI']
    gs=code_data(controls)[0]
    one_errors=[syndrome(xz('I'*i+p+'I'*(4-i)),gs)
                for i in range(5) for p in 'XYZ']
    assert len(one_errors)==15
    assert len(set(one_errors))==5
    assert 0 in one_errors
    counts=dict(sorted(Counter(one_errors).items()))
    # A generic linear 4-bit space automatically offers the same 15
    # hyperplanes, irrespective of whether its measurement scheme can correct.
    num=0
    for a in range(1,16):
        _, w7,w8=hamming_extension(hyperplane(a))
        assert w8=={0:1,4:14,8:1}
        assert 16+16*w8[4]==240
        num+=1
    assert num==15
    return {'controls':controls,
            'syndrome_rank':4,
            'distinct_single_error_syndromes_including_zero':len(set(one_errors)),
            'syndrome_multiplicities':counts,
            'corrects_single_unknown_qubit_error':False,
            'hyperplanes_available_in_same_F2_4':15,
            'E8_isomorphism_from_any_hyperplane':True,
            'matched_structural_null':'15/15 rank-3 hyperplanes yield the same Hamming/E8 output even for a non-correcting code',
            'verdict':'The syndrome-to-E8 chain is structurally valid but generic, not evidence of unique physical selection.'}

if __name__=='__main__':print(json.dumps(test(),indent=2))
