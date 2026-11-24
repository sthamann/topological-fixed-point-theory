"""Smallest protected binary quantum memory and its potential real marking.

Run with Python 3.10+. No packages, approximation, or randomization.
A signed stabilizer sector is counted once for each underlying Pauli space.
"""
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, prod

def span(vectors):
    result = {0}
    for v in vectors:
        result |= {a ^ v for a in result}
    return result

def pauli(word):
    x = sum(1 << i for i,c in enumerate(word) if c in 'XY')
    z = sum(1 << i for i,c in enumerate(word) if c in 'ZY')
    return x | (z << len(word))

def anti(a,b,n):
    mask = (1 << n)-1
    return ((((a & mask) & (b >> n)) ^ ((b & mask) & (a >> n))).bit_count() & 1)

def syndrome(p,checks,n):
    return sum(anti(p,c,n) << i for i,c in enumerate(checks))

def all_spreads(unused,lines):
    if not unused:
        yield frozenset()
        return
    anchor = min(unused)
    for line in lines:
        if anchor in line and line <= unused:
            for tail in all_spreads(unused-line,lines):
                yield tail | {line}

def linear_map(columns):
    images = [0]
    for v in columns:
        images += [a ^ v for a in images]
    return tuple(images)

def compose(a,b):
    return tuple(a[x] for x in b)

def test():
    # The quantum Singleton bound n-k >= 2(d-1) applies to general codes.
    logical_qubits = 1
    errors_to_correct = 1
    distance_needed = 2*errors_to_correct + 1
    n = logical_qubits + 2*(distance_needed-1)
    seed = 'XZZXI'
    checks = tuple(pauli(seed[i:]+seed[:i]) for i in range(n-1))
    stable = span(checks)
    assert len(stable) == 1 << (n-logical_qubits)
    assert all(not anti(a,b,n) for a,b in combinations(checks,2))
    syndrome_lines = []
    for i in range(n):
        sx = syndrome(1 << i, checks,n)
        sz = syndrome(1 << (i+n), checks,n)
        syndrome_lines.append(frozenset([sx,sz,sx^sz]))
    universe = frozenset(range(1,1<<len(checks)))
    assert frozenset().union(*syndrome_lines) == universe
    assert all(len(line)==3 for line in syndrome_lines)
    assert all(not (a & b) for a,b in combinations(syndrome_lines,2))
    logical_support = {
        (p | (p >> n)) & ((1<<n)-1)
        for p in range(1<<(2*n))
        if p not in stable and syndrome(p,checks,n)==0
    }
    d = min(s.bit_count() for s in logical_support)
    correctable = {}
    for erased in range(n+1):
        recover = 0
        for sites in combinations(range(n),erased):
            mask = sum(1 << i for i in sites)
            recover += all(s & mask != s for s in logical_support)
        correctable[erased] = recover
    assert correctable == {0:1,1:5,2:10,3:0,4:0,5:0}
    print('minimum n, attained d:',n,d)
    print('correctable erasure sets by lost-site count:',correctable)

    # Count the complete same-size class of labelled check spaces.
    lines = {frozenset([a,b,a^b]) for a,b in combinations(universe,2)}
    spreads = set(all_spreads(universe, lines))
    assert len(spreads) == 56
    # In any basis of each local plane, AB^T+BA^T is unchanged.
    # Check commutation of the four rows for every spread.
    for spread in spreads:
        columns = [sorted(line)[:2] for line in spread]
        for i,j in combinations(range(4),2):
            assert sum((((a>>i)&(b>>j)) ^ ((a>>j)&(b>>i))) & 1
                       for a,b in columns) % 2 == 0
    ordered_bases = prod(16-2**i for i in range(4))
    # Order the five syndrome lines, choose an ordered local basis in each,
    # then remove the arbitrary choice of the ordered global syndrome basis.
    good,rest = divmod(len(spreads)*factorial(n)*6**n,ordered_bases)
    assert not rest
    # Successive checks lie in an orthogonal complement and outside their span.
    total,rest = divmod(prod(2**(2*n-i)-2**i for i in range(4)),ordered_bases)
    assert not rest
    print('complete distinct syndrome table:',good,'/',total,'=',Fraction(good,total))

    # GL(4,2) acts on the concrete spread. The projective line over F4
    # emerges via its native order-three scalar, recovered by the kernel below.
    indices = {l:i for i,l in enumerate(syndrome_lines)}
    perms = {}
    for basis in permutations(universe,4):
        trans = linear_map(basis)
        if len(set(trans)) != 16:
            continue
        images = [frozenset(trans[x] for x in l) for l in syndrome_lines]
        if all(l in indices for l in images):
            perms[trans] = tuple(indices[l] for l in images)
    site_identity = tuple(range(n))
    identity = tuple(range(16))
    kernel = [t for t,p in perms.items() if p==site_identity]
    assert len(kernel)==3
    scalar = next(t for t in kernel if t!=identity)
    scalar_inv = compose(scalar,scalar)
    assert scalar_inv != identity and compose(scalar_inv,scalar) == identity
    # v -> scalar(v) satisfies x^2+x+1=0, defining F4 multiplication.
    assert all((scalar_inv[x] ^ scalar[x] ^ x) == 0 for x in range(16))
    all_actions = set(perms.values())
    involutions = {p for p in all_actions if p!=site_identity and compose(p,p)==site_identity}
    frob_lifts = [t for t in perms if compose(t,t)==identity and
                  compose(t,scalar)==compose(scalar_inv,t)]
    frob_actions = {perms[t] for t in frob_lifts}
    marks = {frozenset(i for i in range(n) if p[i]==i) for p in frob_actions}
    all_triplets = {frozenset(s) for s in combinations(range(n),3)}
    invariant = {s for s in all_triplets
                 if all(frozenset(p[i] for i in s)==s for p in all_actions)}
    assert len(perms)==360 and len(all_actions)==120
    assert len(involutions)==25 and len(frob_lifts)==30 and len(frob_actions)==10
    assert marks==all_triplets and not invariant
    print('incidence automorphisms / site permutations:',len(perms),'/',len(all_actions))
    print('Frobenius site involutions / site involutions:',len(frob_actions),'/',len(involutions))
    print('Frobenius fixed triples / unconditionally invariant triples:',len(marks),'/',len(invariant))
    # Standard physical complex conjugation K sends Y to -Y. Pauli labels
    # lose signs, so K fixes every local label and induces the identity.
    assert all(sum(c=='Y' for c in seed[i:]+seed[:i])%2==0 for i in range(n-1))

if __name__ == '__main__':
    test()
