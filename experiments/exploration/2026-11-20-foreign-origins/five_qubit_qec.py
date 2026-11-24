#!/usr/bin/env python3
"""Exact GF(2) and SymPy verification of [[5,1,3]] and an AME(5,2) 2->3 tensor.
Physical assumptions (qubits, correction of one arbitrary Pauli error, minimality)
are declared, not inferred from TFPT.
"""
from itertools import combinations, product
from collections import Counter
import json
from sympy import Matrix, eye, kronecker_product, sqrt, simplify, I as imaginary

N=5
ALPH='IXYZ'
SEED='XZZXI'
GENS=[SEED[-i:]+SEED[:-i] if i else SEED for i in range(4)]


def xz(w):
    x=z=0
    for i,ch in enumerate(w):
        if ch in 'XY': x |= (1<<i)
        if ch in 'ZY': z |= (1<<i)
    return (x,z)


def mul(a,b): return (a[0]^b[0], a[1]^b[1])  # disregards global phase only


def commute(a,b):
    return ((a[0]&b[1]).bit_count() + (a[1]&b[0]).bit_count())%2==0


def weight(a): return (a[0]|a[1]).bit_count()


def span(gs):
    s={(0,0)}
    for g in gs:
        s |= {mul(h,g) for h in tuple(s)}
    return s


def all_paulis(n=N):
    for seq in product(ALPH,repeat=n):
        yield xz(''.join(seq))


def code_data(words):
    gs=list(map(xz,words)); st=span(gs)
    if len(st)!=16 or any(not commute(a,b) for a,b in combinations(gs,2)):
        return None
    cen={w for w in all_paulis() if all(commute(w,g) for g in gs)}
    logical=cen-st
    d=min(map(weight,logical))
    return gs,st,cen,logical,d


def syndrome(err,gs):
    return sum((not commute(err,g))<<i for i,g in enumerate(gs))


def exact_stabilizer_state(words):
    I=eye(2)
    X=Matrix([[0,1],[1,0]])
    Z=Matrix([[1,0],[0,-1]])
    Y=Matrix([[0,-imaginary],[imaginary,0]])
    paulis={'I':I,'X':X,'Y':Y,'Z':Z}
    psi=Matrix([1]+[0]*31)
    for w in words:
        op=kronecker_product(*[paulis[c] for c in w])
        psi=(psi+op*psi)/2
    # Projected |00000> is a nonzero logical-Z eigenstate.
    norm2=simplify((psi.conjugate().T*psi)[0])
    assert norm2>0
    psi=psi/sqrt(norm2)
    return psi


def reduced_gram(psi, positions):
    positions=tuple(positions)
    outpos=tuple(i for i in range(5) if i not in positions)
    M=Matrix.zeros(2**len(outpos),2**len(positions))
    for idx in range(32):
        bits=tuple((idx>>(4-i))&1 for i in range(5))
        a=sum(bits[p]<<(len(positions)-1-j) for j,p in enumerate(positions))
        b=sum(bits[p]<<(len(outpos)-1-j) for j,p in enumerate(outpos))
        M[b,a]=psi[idx]
    return M.conjugate().T*M


def cyclic_seed_baseline():
    total=0; commuting_independent=0; perfect=0; examples=[]
    for blank in range(5):
        for letters in product('XYZ',repeat=4):
            chars=iter(letters)
            word=''.join('I' if j==blank else next(chars) for j in range(5))
            words=[word[-i:]+word[:-i] if i else word for i in range(4)]
            total+=1
            c=code_data(words)
            if c is None:continue
            commuting_independent+=1
            if c[-1]>=3:
                perfect+=1
                if len(examples)<8:examples.append(word)
    return {'total_cyclic_seeds_one_identity':total,
            'commuting_independent':commuting_independent,
            'distance_at_least_3':perfect,
            'examples':examples}


def test():
    # Minimal length follows from Singleton: n-k >= 2(d-1), with k=1, d>=3.
    min_n=1+2*(3-1)
    assert min_n==5
    gs,st,cen,logical,d=code_data(GENS)
    assert d==3 and len(st)==16 and len(cen)==64 and len(logical)==48
    errors=[xz('IIIII')]
    for i in range(5):
        for ch in 'XYZ':errors.append(xz('I'*i+ch+'I'*(4-i)))
    syndromes={syndrome(e,gs) for e in errors}
    assert len(errors)==16 and len(syndromes)==16 and syndromes==set(range(16))
    # For every 3-site authorized subset, find two anti-commuting logical Paulis there.
    authorized=0
    for triple in combinations(range(5),3):
        mask=sum(1<<i for i in triple)
        ops=[a for a in logical if ((a[0]|a[1]) & ~mask)==0]
        assert any(not commute(a,b) for a,b in combinations(ops,2))
        authorized+=1
    # For every 2-site subset, there is no logical operator supported there.
    assert all(not any(((a[0]|a[1]) & ~sum(1<<i for i in duo))==0
                       for a in logical) for duo in combinations(range(5),2))
    # Cyclic relabelling is an exact code symmetry; no 3-site color set invariant.
    gs_rot=[xz(w[-1:]+w[:-1]) for w in GENS]
    assert span(gs_rot)==st
    # The cyclic fifth generator equals the signed matrix product, not only modulo Pauli phase.
    from sympy import eye as _eye
    ops={"I":_eye(2),"X":Matrix([[0,1],[1,0]]),"Z":Matrix([[1,0],[0,-1]])}
    def mat(w):return kronecker_product(*[ops[c] for c in w])
    product4=_eye(32)
    for w in GENS:product4=product4*mat(w)
    assert product4==mat(SEED[-4:]+SEED[:-4])
    assert not any({(i+1)%5 for i in t}==set(t) for t in combinations(range(5),3))
    psi=exact_stabilizer_state(GENS)
    marginals=[reduced_gram(psi,duo) for duo in combinations(range(5),2)]
    assert all(g==eye(4)/4 for g in marginals)
    assert simplify((psi.conjugate().T*psi)[0])==1
    amplitude_counts=Counter(str(z) for z in psi if z!=0)
    return {'singleton_min_physical_qubits':min_n,'code':'[[5,1,3]]',
            'generators':GENS,'stabilizer_size':len(st),
            'normalizer_size':len(cen),'distance':d,
            'one_error_syndromes':len(syndromes),
            'all_3_site_subsets_reconstruct':authorized,
            'all_2_site_subsets_hide':10,
            'AME_5_2_two_site_marginals_I_over_4':len(marginals),
            '2_to_3_isometry_all_splits':len(marginals),
            'state_nonzero_amplitudes':sum(amplitude_counts.values()),
            'state_amplitude_counts':dict(amplitude_counts),
            'cyclic_transitivity_forbids_invariant_3_site_split':True,
            'null_cyclic_seeds':cyclic_seed_baseline()}


if __name__=='__main__':
    print(json.dumps(test(),indent=2,ensure_ascii=False))
