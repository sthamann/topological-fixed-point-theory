"""Which permutations of the five return labels preserve full flag multiplication?"""
from itertools import permutations
from math import comb
import sympy as sp
from flag_process import FLAGS, I, N, change, tr


def main():
    a, b = change(0), change(1)
    m = [I, a, b, a*b, b*a, a*b*a]
    sizes = [sum(x != 0 for x in z[0, :]) for z in m]
    assert sizes == [1, 2, 2, 4, 4, 8]
    gram = [tr(z.T*z) for z in m]
    # Their disjoint masks partition matrix entries, so one entry per mask
    # reads all six basis coefficients in any product.
    sample = [next(k for k in range(N*N) if z[k] != 0) for z in m]
    constants = {}
    for i, u in enumerate(m):
        for j, v in enumerate(m):
            uv = u*v
            c = tuple(uv[k]/z[k] for k, z in zip(sample, m))
            assert uv == sum((w*z for w, z in zip(c, m)), sp.zeros(N))
            constants[i, j] = c
    def multiplicative(p):
        return all(constants[p[i],p[j]][p[k]] == constants[i,j][k]
                   for i in range(6) for j in range(6) for k in range(6))
    mappings = []
    for odd in permutations((1,2,5)):
        for even in permutations((3,4)):
            mappings.append((0,odd[0],odd[1],even[0],even[1],odd[2]))
    preserving = [p for p in mappings if multiplicative(p)]
    assert len(preserving) == 2
    group = list(permutations(range(3)))
    def mul(p, q):
        return tuple(p[i] for i in q)
    e, s, t = (0,1,2), (1,0,2), (0,2,1)
    labels = [e,s,t,mul(s,t),mul(t,s),mul(s,mul(t,s))]
    def action(p):
        pinv = tuple(p.index(i) for i in e)
        return tuple(labels.index(mul(p,mul(w,pinv))) for w in labels)
    assert set(preserving) == {action(e),action(labels[5])}
    print("Parity-preserving label permutations preserving multiplication:",
          len(preserving), "/",len(mappings))
    print("Triangle relabelings preserving all flag products:",
          sum(multiplicative(action(p)) for p in group), "/",len(group))
    print("Squared normalized HS norms of raw odd galleries:",[gram[i] for i in (1,2,5)])
    # A chart names the transposition at distance three. Its stabilizer
    # fixes the chart or reverses its two elementary types by point-line duality.
    odd = (1,2,5)
    representatives = {i:next(p for p in group if action(p)[5]==i) for i in odd}
    dual = sp.Matrix([[int(f==(l,p)) for f in FLAGS] for p,l in FLAGS])
    assert dual.T == dual and dual*a*dual == b
    def inv(p):
        return tuple(p.index(i) for i in e)
    def block(c, i):
        return m[action(inv(representatives[c]))[i]]
    for g in group:
        for c in odd:
            dest = action(g)[c]
            residual = mul(inv(representatives[dest]),mul(g,representatives[c]))
            assert residual in (e,labels[5])
            transport = I if residual == e else dual
            for i in range(6):
                assert transport*block(c,i)*transport.T == block(dest,action(g)[i])
    # A single pair of named generators need not braid in another chart;
    # there the two elementary generators have different names.
    braid_per_chart = []
    for c in odd:
        u, v = block(c,1),block(c,2)
        braid_per_chart.append(u*v*u == v*u*v)
    assert sum(braid_per_chart)==1
    print("Three-chart covariance: 6 x 3 x 6 blocks; fixed named pair braids in",
          sum(braid_per_chart),"/ 3 charts")
    q=sp.symbols("q", positive=True)
    # If chart i has weight pi, its long transposition has q^3 reachable
    # flags; that label is elementary in each other chart. Both 3-cycles
    # always have q^2 reachable flags. Exact equal visibility fixes pi.
    pi=sp.factor((q**2-q)/(q**3-q))
    assert pi == 1/(q+1) and sp.solve(3*pi-1,q)==[2]
    mismatch=sp.factor((2*q+q**3)/3-q**2)
    assert mismatch==q*(q-2)*(q-1)/3
    mean_sizes=[sum(sp.Integer(sizes[action(inv(p))[i]]) for p in representatives.values())/3
                for i in range(6)]
    assert mean_sizes==[1,4,4,4,4,4]
    print("Equal raw visibility in every active class: chart weights",pi,
          "; unique thick order q=2; raw counts",mean_sizes)
    print("Uniform chart raw excess of each transposition over each 3-cycle:",mismatch)
    inverse = (0,1,2,4,3,5)
    good = 0
    for p in permutations(range(6)):
        # Supports of distinct Weyl classes are disjoint.
        masses = [int(p[i]==inverse[i]) for i in range(6)]
        good += all(masses)
    assert good == 1
    print("Permutations of reverse Weyl types giving all six successes:",good,"/ 720")
    # Six weights in units of 1/6. Choi orthogonality makes the law
    # injective on the identity and the five distinct pinching channels.
    print("Uniform recovered channel on the 1/6 probability grid: 1 /",comb(11,5))


if __name__ == "__main__":
    main()
